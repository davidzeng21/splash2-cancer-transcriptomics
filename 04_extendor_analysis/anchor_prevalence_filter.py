#!/usr/bin/env python3
"""
Anchor-level tumor (and normal) prevalence from SATC dumps, then filter
extendor_proportions.parquet by prevalent anchors and log10FC thresholds.

Does not modify extendor_analysis.py or regenerate the parquet.
"""

from __future__ import annotations

import argparse
import math
from datetime import datetime, timezone
from pathlib import Path

import polars as pl


def parse_sample_metadata(input_file: Path) -> tuple[set[str], set[str]]:
    """Same rules as extendor_analysis.parse_sample_metadata (suffix -tumor / -normal)."""
    tumor_samples: set[str] = set()
    normal_samples: set[str] = set()
    with input_file.open("r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 1:
                continue
            sample_name = parts[0]
            if sample_name.endswith("-tumor"):
                tumor_samples.add(sample_name)
            elif sample_name.endswith("-normal"):
                normal_samples.add(sample_name)
    print(
        f"Found {len(tumor_samples)} tumor samples and {len(normal_samples)} normal samples"
    )
    return tumor_samples, normal_samples


def frac_slug(min_frac: float) -> str:
    """Stable tag for filenames, e.g. 0.10 -> 10pct."""
    pct = min_frac * 100.0
    if abs(pct - round(pct)) < 1e-9:
        return f"{int(round(pct))}pct"
    return str(min_frac).replace(".", "p")


def tumor_prevalence_threshold(
    n_tumor_samples: int,
    min_frac: float,
    rounding: str,
) -> float | int:
    """Return numeric threshold for n_tumor_nonzero (int for ceil/floor, float for strict_gt)."""
    raw = n_tumor_samples * min_frac
    if rounding == "ceil":
        return max(1, int(math.ceil(raw - 1e-12)))
    if rounding == "floor":
        return max(0, int(math.floor(raw + 1e-12)))
    if rounding == "strict_gt":
        return raw
    raise ValueError(f"Unknown --minTumorFracRounding: {rounding}")


def resolve_dump_glob(dump_path: Path, dump_glob: str) -> str:
    """Single .dump file or directory + glob for polars scan_csv."""
    if dump_path.is_file():
        return str(dump_path.resolve())
    return str((dump_path / dump_glob).resolve())


def collect_streaming(lf: pl.LazyFrame) -> pl.DataFrame:
    return lf.collect(engine="streaming")


def anchor_prevalence_from_dumps(
    dump_glob: str,
    tumor_samples: set[str],
    normal_samples: set[str],
) -> pl.DataFrame:
    tumor_list = list(tumor_samples)
    normal_list = list(normal_samples)
    lf = pl.scan_csv(
        dump_glob,
        separator="\t",
        has_header=False,
        new_columns=["sample", "anchor", "target", "count"],
        schema_overrides={"count": pl.Int32},
    )
    lf = lf.select(pl.col("sample"), pl.col("anchor")).with_columns(
        pl.when(pl.col("sample").is_in(tumor_list))
        .then(pl.lit("tumor"))
        .when(pl.col("sample").is_in(normal_list))
        .then(pl.lit("normal"))
        .otherwise(pl.lit("unknown"))
        .alias("group")
    )
    lf = lf.filter(pl.col("group") != "unknown")
    lf = lf.group_by(["anchor", "group"]).agg(
        pl.col("sample").n_unique().alias("n_samples")
    )
    long_df = collect_streaming(lf)
    wide = long_df.pivot(on="group", index="anchor", values="n_samples")
    # Ensure columns exist
    for col, default in ("tumor", 0), ("normal", 0):
        if col not in wide.columns:
            wide = wide.with_columns(pl.lit(default).alias(col))
    wide = wide.rename({"tumor": "n_tumor_nonzero", "normal": "n_normal_nonzero"})
    wide = wide.with_columns(
        pl.col("n_tumor_nonzero").fill_null(0).cast(pl.Int32),
        pl.col("n_normal_nonzero").fill_null(0).cast(pl.Int32),
    )
    n_tumor = len(tumor_samples)
    n_normal = len(normal_samples)
    wide = wide.with_columns(
        (pl.col("n_tumor_nonzero") / float(n_tumor)).alias("tumor_frac"),
        (pl.col("n_normal_nonzero") / float(n_normal)).alias("normal_frac"),
    )
    return wide.sort("anchor")


def build_whitelist(
    prev: pl.DataFrame,
    rounding: str,
    min_frac: float,
    n_tumor_samples: int,
) -> pl.DataFrame:
    thr = tumor_prevalence_threshold(n_tumor_samples, min_frac, rounding)
    if rounding == "strict_gt":
        return prev.filter(pl.col("n_tumor_nonzero") > thr)
    return prev.filter(pl.col("n_tumor_nonzero") >= int(thr))


def count_tsv_data_rows(path: Path) -> int:
    with path.open("rb") as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)


def filter_parquet_write_tsv(
    prop_parquet: Path,
    whitelist: list[str],
    fc: float,
    out_tsv: Path,
) -> int:
    wl = pl.DataFrame({"anchor": whitelist})
    lf = (
        pl.scan_parquet(str(prop_parquet))
        .join(wl.lazy(), on="anchor", how="semi")
        .filter(pl.col("log10FC") >= float(fc))
    )
    df = collect_streaming(lf)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(out_tsv, separator="\t")
    return len(df)


def append_summary_stats(
    summary_path: Path,
    *,
    min_frac: float,
    rounding: str,
    thr_display: str,
    n_tumor: int,
    n_normal: int,
    n_anchors_total: int,
    n_anchors_prevalent: int,
    fc_counts: dict[float, int],
    legacy_lines: dict[float, int | None],
    frac_tag: str,
) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "",
        f"--- anchor prevalence filter ({ts}) ---",
        f"min_tumor_frac: {min_frac}",
        f"min_tumor_frac_rounding: {rounding}",
        f"min_tumor_samples_required: {thr_display}",
        f"n_tumor_samples: {n_tumor}",
        f"n_normal_samples: {n_normal}",
        f"n_anchors_total: {n_anchors_total:,}",
        f"n_anchors_tumor_prevalent_ge_{frac_tag}: {n_anchors_prevalent:,}",
    ]
    for fc in sorted(fc_counts):
        lines.append(
            f"n_extendors_log10FC_ge{fc:g}_and_anchor_prevalent_ge_{frac_tag}: "
            f"{fc_counts[fc]:,}"
        )
    for fc in sorted(legacy_lines):
        leg = legacy_lines[fc]
        if leg is not None:
            lines.append(
                f"legacy_n_data_rows_log10FC_ge{fc:g}_count_normal_ge10_tsv: {leg:,}"
            )
    lines.append("")
    with summary_path.open("a") as f:
        f.write("\n".join(lines))


def parse_fc_thresholds(s: str) -> list[float]:
    out: list[float] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(float(part))
    if not out:
        raise ValueError("fcThresholds is empty")
    return out


def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Anchor-level tumor prevalence from SATC dumps; filter extendor parquet."
    )
    p.add_argument(
        "--inputFile",
        type=str,
        default="input_CPTAC_107.txt",
        help="Sample list (tumor/normal by suffix), same as extendor_analysis",
    )
    p.add_argument(
        "--dumpFolder",
        type=str,
        required=True,
        help="Directory containing *.dump files, or path to a single .dump file",
    )
    p.add_argument(
        "--dumpGlob",
        type=str,
        default="*.dump",
        help="Glob under dumpFolder (ignored if dumpFolder is a file)",
    )
    p.add_argument(
        "--propParquet",
        type=str,
        required=True,
        help="extendor_proportions.parquet path",
    )
    p.add_argument(
        "--outFolder",
        type=str,
        required=True,
        help="Output directory (anchor_prevalence table, whitelist, filtered TSVs)",
    )
    p.add_argument(
        "--minTumorFrac",
        type=float,
        default=0.10,
        help="Minimum fraction of tumor samples with count>0 for an anchor",
    )
    p.add_argument(
        "--minTumorFracRounding",
        type=str,
        choices=("ceil", "floor", "strict_gt"),
        default="ceil",
        help="How to compare n_tumor_nonzero to n_tumor*minTumorFrac",
    )
    p.add_argument(
        "--fcThresholds",
        type=str,
        default="2,3",
        help="Comma-separated log10FC thresholds for filtered extendor TSVs",
    )
    p.add_argument(
        "--summaryStats",
        type=str,
        default=None,
        help="Append summary block here (default: extendor_analysis/summary_stats.txt next to outFolder)",
    )
    p.add_argument(
        "--noSummaryAppend",
        action="store_true",
        help="Do not append summary_stats.txt",
    )
    p.add_argument(
        "--legacyTsvParent",
        type=str,
        default=None,
        help="Directory to find extendor_proportions_log10FC_ge*_count_normal_ge10.tsv for comparison lines",
    )
    return p.parse_args()


def main() -> None:
    args = get_args()
    cwd = Path.cwd()
    input_file = Path(args.inputFile).expanduser()
    if not input_file.is_absolute():
        input_file = (cwd / input_file).resolve()

    dump_path = Path(args.dumpFolder).expanduser()
    if not dump_path.is_absolute():
        dump_path = (cwd / dump_path).resolve()

    prop_parquet = Path(args.propParquet).expanduser()
    if not prop_parquet.is_absolute():
        prop_parquet = (cwd / prop_parquet).resolve()

    out_folder = Path(args.outFolder).expanduser()
    if not out_folder.is_absolute():
        out_folder = (cwd / out_folder).resolve()

    out_folder.mkdir(parents=True, exist_ok=True)

    tumor_samples, normal_samples = parse_sample_metadata(input_file)
    dump_glob = resolve_dump_glob(dump_path, args.dumpGlob)
    print(f"Scanning dumps: {dump_glob}")

    prev = anchor_prevalence_from_dumps(dump_glob, tumor_samples, normal_samples)
    n_anchors_total = len(prev)
    thr = tumor_prevalence_threshold(
        len(tumor_samples), args.minTumorFrac, args.minTumorFracRounding
    )
    prevalent_tbl = build_whitelist(
        prev,
        args.minTumorFracRounding,
        args.minTumorFrac,
        len(tumor_samples),
    )
    whitelist = prevalent_tbl["anchor"].to_list()
    n_prev = len(whitelist)
    frac_tag = frac_slug(args.minTumorFrac)

    prev_tsv = out_folder / "anchor_prevalence.tsv"
    prev.write_csv(prev_tsv, separator="\t")
    print(f"Wrote {prev_tsv} ({n_anchors_total:,} anchors)")

    wl_path = out_folder / f"anchor_tumor_frac_ge_{frac_tag}_whitelist.txt"
    with wl_path.open("w") as f:
        for a in whitelist:
            f.write(f"{a}\n")
    print(f"Wrote {wl_path} ({n_prev:,} anchors)")

    fc_list = parse_fc_thresholds(args.fcThresholds)
    fc_counts: dict[float, int] = {}
    for fc in fc_list:
        out_tsv = out_folder / (
            f"extendor_proportions_log10FC_ge{fc:g}_anchor_tumor_frac_ge_{frac_tag}.tsv"
        )
        n_rows = filter_parquet_write_tsv(prop_parquet, whitelist, fc, out_tsv)
        fc_counts[fc] = n_rows
        print(f"Wrote {out_tsv} ({n_rows:,} rows)")

    if args.noSummaryAppend:
        return

    if args.summaryStats:
        summary_path = Path(args.summaryStats).expanduser()
        if not summary_path.is_absolute():
            summary_path = (cwd / summary_path).resolve()
    else:
        summary_path = out_folder.parent / "summary_stats.txt"

    legacy_parent = args.legacyTsvParent
    if legacy_parent:
        leg_root = Path(legacy_parent).expanduser()
        if not leg_root.is_absolute():
            leg_root = (cwd / leg_root).resolve()
    else:
        leg_root = out_folder.parent

    legacy_lines: dict[float, int | None] = {}
    for fc in fc_list:
        leg_path = leg_root / (
            f"extendor_proportions_log10FC_ge{fc:g}_count_normal_ge10.tsv"
        )
        if leg_path.is_file():
            legacy_lines[fc] = count_tsv_data_rows(leg_path)
        else:
            legacy_lines[fc] = None

    thr_display = (
        f">{thr:.12g}"
        if args.minTumorFracRounding == "strict_gt"
        else str(int(thr))
    )
    append_summary_stats(
        summary_path,
        min_frac=args.minTumorFrac,
        rounding=args.minTumorFracRounding,
        thr_display=thr_display,
        n_tumor=len(tumor_samples),
        n_normal=len(normal_samples),
        n_anchors_total=n_anchors_total,
        n_anchors_prevalent=n_prev,
        fc_counts=fc_counts,
        legacy_lines=legacy_lines,
        frac_tag=frac_tag,
    )
    print(f"Appended summary block to {summary_path}")


if __name__ == "__main__":
    main()
