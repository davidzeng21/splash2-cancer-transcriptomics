#!/usr/bin/env python3
"""
Filter ge2 BLAST RNA / GRCh38 outputs to extendors present in
extendor_proportions_log10FC_ge2_anchor_tumor_frac_ge_10pct.tsv.
Writes mirrored outputs under anchor_prevalence/blast_ge2_filtered/.
"""

from __future__ import annotations

import argparse
import re
from collections.abc import Callable
from pathlib import Path

import polars as pl


def load_prevalent_ge2_extendors(tsv: Path) -> set[str]:
    df = pl.read_csv(tsv, separator="\t", columns=["extendor"])
    return set(df["extendor"].to_list())


def filter_tsv_by_extendor_col(
    src: Path,
    dst: Path,
    ok: set[str],
    col: str,
) -> int:
    df = pl.read_csv(src, separator="\t")
    if col not in df.columns:
        raise ValueError(f"{src}: missing column {col!r}")
    filt = df.filter(pl.col(col).is_in(list(ok)))
    dst.parent.mkdir(parents=True, exist_ok=True)
    filt.write_csv(dst, separator="\t")
    return len(filt)


def extendor_from_blast_qseqid_cell(cell: str) -> str:
    return cell.split("|", 1)[0]


def filter_blast_tab_no_header(
    src: Path,
    dst: Path,
    ok: set[str],
    extendor_from_line: Callable[[str], str],
) -> int:
    dst.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with src.open() as inf, dst.open("w") as outf:
        for line in inf:
            if not line.strip():
                continue
            e = extendor_from_line(line)
            if e in ok:
                outf.write(line)
                kept += 1
    return kept


def filter_large_blast_tabular(
    src: Path,
    dst: Path,
    ok: set[str],
    *,
    progress_every: int = 2_000_000,
) -> tuple[int, int]:
    """
    Stream-filter BLAST -outfmt 6 (tabular) lines: qseqid is column 1 as extendor|metadata.
    Uses binary I/O for ~76M+ line files (e.g. ge2_human_genome.txt).
    Returns (n_lines_kept, n_lines_read).
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    total = 0
    buf = 8 * 1024 * 1024
    with src.open("rb", buffering=buf) as inf, dst.open("wb", buffering=buf) as outf:
        for raw in inf:
            total += 1
            tab = raw.find(b"\t")
            if tab <= 0:
                continue
            pipe = raw.find(b"|", 0, tab)
            key = raw[: pipe if pipe != -1 else tab]
            try:
                ext = key.decode("utf-8")
            except UnicodeDecodeError:
                ext = key.decode("utf-8", errors="replace")
            if ext in ok:
                outf.write(raw)
                kept += 1
            if progress_every and total % progress_every == 0:
                print(
                    f"  ge2_human_genome.txt: {total:,} lines read, {kept:,} kept",
                    flush=True,
                )
    return kept, total


def filter_rna_no_hit_txt(src: Path, dst: Path, ok: set[str]) -> int:
    dst.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with src.open() as inf, dst.open("w") as outf:
        header = inf.readline()
        outf.write(header)
        for line in inf:
            if not line.strip():
                continue
            ext = line.split("\t", 1)[0].strip()
            if ext in ok:
                outf.write(line)
                kept += 1
    return kept


def build_gene_summary_from_best_hits(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.group_by("gene_name")
        .agg(
            pl.len().alias("n_extendors"),
            pl.col("log10FC").median().alias("median_log10FC"),
            pl.col("log10FC").max().alias("max_log10FC"),
            pl.col("pident").median().alias("median_pident"),
        )
        .sort("n_extendors", descending=True)
    )


def write_gene_summary_tsv(df: pl.DataFrame, path: Path) -> None:
    """Match loose spacing style of original ge2_gene_summary.tsv (tab-separated)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "gene_name\tn_extendors\tmedian_log10FC\tmax_log10FC\tmedian_pident",
    ]
    for row in df.iter_rows(named=True):
        lines.append(
            f"{row['gene_name']}\t{row['n_extendors']}\t{row['median_log10FC']}\t"
            f"{row['max_log10FC']}\t{row['median_pident']}"
        )
    path.write_text("\n".join(lines) + "\n")


def parse_section_genes(header_line: str) -> frozenset[str]:
    """First line of a ### section (includes leading ###)."""
    inner = header_line.strip().removeprefix("###").strip()
    left = inner.split("—", 1)[0].strip()
    if " / " in left:
        parts = [p.strip().split()[0] for p in left.split(" / ")]
        return frozenset(parts)
    return frozenset({left.split()[0]})


def update_section_header_counts(
    header_line: str,
    genes: frozenset[str],
    gene_counts: dict[str, int],
) -> str | None:
    """Return updated header or None if section should be dropped."""
    counts = [gene_counts.get(g, 0) for g in sorted(genes)]
    if sum(counts) == 0:
        return None
    if len(genes) == 2:
        a, b = sorted(genes)
        ca, cb = gene_counts.get(a, 0), gene_counts.get(b, 0)
        return re.sub(
            r"—\s*\d+\s*\+\s*\d+\s*extendors",
            f"— {ca} + {cb} extendors",
            header_line,
        )
    g = next(iter(genes))
    n = gene_counts[g]
    line = re.sub(r"—\s*\d+\s*extendors", f"— {n} extendors", header_line)
    line = re.sub(r"—\s*\d+\s*extendor", f"— {n} extendor", line)
    return line


def rebuild_key_findings_table(md_preamble: str, gene_df: pl.DataFrame) -> str:
    """Replace the markdown table under ## Key Findings with recomputed rows."""
    top = gene_df.head(21)
    rows = [
        "| Gene | N extendors | Median log10FC | Role in LUAD |",
        "|------|------------|---------------|-------------|",
    ]
    role_stub = "See gene-by-gene section below."
    for row in top.iter_rows(named=True):
        rows.append(
            f"| {row['gene_name']} | {row['n_extendors']} | "
            f"{row['median_log10FC']:.4g} | {role_stub} |"
        )
    new_table = "\n".join(rows) + "\n"
    pattern = r"(## Key Findings at a Glance\n\n)(?:\|[^\n]+\n)+"
    if not re.search(pattern, md_preamble):
        return md_preamble
    return re.sub(pattern, r"\1" + new_table, md_preamble, count=1)


def write_filtered_biology_summary(
    src: Path,
    dst: Path,
    gene_df: pl.DataFrame,
    gene_counts: dict[str, int],
    n_ge2_no_human_filtered: int,
) -> None:
    text = src.read_text()
    parts = re.split(r"\n(?=### )", text, maxsplit=1)
    preamble = parts[0]
    body = parts[1] if len(parts) > 1 else ""

    note = (
        "\n\n*This version retains only extendors whose anchor passes "
        "≥10% tumor-sample prevalence (same subset as "
        "`extendor_proportions_log10FC_ge2_anchor_tumor_frac_ge_10pct.tsv`).*\n"
    )
    if "≥10% tumor-sample prevalence" not in preamble:
        preamble = preamble.rstrip() + note + "\n"

    preamble = rebuild_key_findings_table(preamble, gene_df)

    section_blocks = re.split(r"\n(?=### )", body)
    out_sections: list[str] = []
    for block in section_blocks:
        if not block.strip():
            continue
        lines = block.split("\n")
        header = lines[0]
        genes = parse_section_genes(header)
        new_header = update_section_header_counts(header, genes, gene_counts)
        if new_header is None:
            continue
        lines[0] = new_header
        out_sections.append("\n".join(lines))

    filtered_body = "\n\n".join(out_sections)
    filtered_body = re.sub(
        r"\(40 in ge3, 71 in ge2\)",
        f"(40 in ge3, {n_ge2_no_human_filtered} in ge2 after anchor prevalence filter)",
        filtered_body,
        count=1,
    )

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(preamble + "\n" + filtered_body)


def write_filtering_notes(path: Path) -> None:
    path.write_text(
        """# blast_ge2_filtered — filtering notes

## What is filtered

- **TSV summaries** (`*_best_hits.tsv`, `ge2_no_human_hits.tsv`): rows whose extendor
  is not in `extendor_proportions_log10FC_ge2_anchor_tumor_frac_ge_10pct.tsv`.
- **BLAST tabular** (`refseq_rna/*_human.txt`, `refseq_rna_select/*_human.txt`,
  `human_database_GRCh38.p13/ge2_human_genome.txt`): each output line is one HSP;
  lines are kept iff qseqid (column 1) begins with `extendor` before `|`, and that
  extendor is in the same prevalent-ge2 set.

## Why some files barely change

Only **one** extendor was dropped from the original ge2 RNA best-hits relative to
the anchor-prevalence ge2 table (DEFA5 example in your cohort). Full tabular files
**do** drop **all** HSP lines for removed extendors (e.g. 8 lines in
`ge2_refseq_rna_human.txt`, 4 in select, 10 in `ge2_human_genome.txt` for that hit).

## Why `ge2_rna_no_hit.txt` can stay the same

That file lists extendors with **no** RefSeq RNA BLAST hit. An extendor removed for
low anchor prevalence but that **had** RNA hits never appears there, so its removal
does not change the no-hit list.

"""
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--prevalentGe2Tsv",
        type=Path,
        required=True,
        help="extendor_proportions_log10FC_ge2_anchor_tumor_frac_ge_10pct.tsv",
    )
    ap.add_argument(
        "--blastRoot",
        type=Path,
        required=True,
        help=".../blast_tumor_specific/BLAST_results",
    )
    ap.add_argument(
        "--outRoot",
        type=Path,
        required=True,
        help="anchor_prevalence/blast_ge2_filtered",
    )
    ap.add_argument(
        "--skipLargeGenome",
        action="store_true",
        dest="skip_large_genome",
        help="Skip streaming filter of ge2_human_genome.txt (very large)",
    )
    ap.add_argument(
        "--onlyLargeGenome",
        action="store_true",
        dest="only_large_genome",
        help="Only run ge2_human_genome.txt streaming filter (requires existing outRoot)",
    )
    args = ap.parse_args()

    ok = load_prevalent_ge2_extendors(args.prevalentGe2Tsv)
    print(f"Prevalent ge2 extendors: {len(ok):,}")

    rna = args.blastRoot / "RNA"
    hg = args.blastRoot / "human_database_GRCh38.p13"
    out_rna = args.outRoot / "RNA"
    out_hg = args.outRoot / "human_database_GRCh38.p13"

    if args.only_large_genome:
        src_g = hg / "ge2_human_genome.txt"
        dst_g = out_hg / "ge2_human_genome.txt"
        print(f"Streaming filter: {src_g} -> {dst_g}")
        k, t = filter_large_blast_tabular(src_g, dst_g, ok)
        print(f"ge2_human_genome.txt: kept {k:,} / {t:,} lines")
        write_filtering_notes(args.outRoot / "FILTERING_NOTES.md")
        return

    # RNA TSV hits
    filter_tsv_by_extendor_col(
        rna / "ge2_rna_best_hits.tsv",
        out_rna / "ge2_rna_best_hits.tsv",
        ok,
        "extendor_id",
    )
    filter_tsv_by_extendor_col(
        rna / "refseq_rna" / "ge2_refseq_rna_best_hits.tsv",
        out_rna / "refseq_rna" / "ge2_refseq_rna_best_hits.tsv",
        ok,
        "extendor_id",
    )
    filter_tsv_by_extendor_col(
        rna / "refseq_rna_select" / "ge2_refseq_select_rna_best_hits.tsv",
        out_rna / "refseq_rna_select" / "ge2_refseq_select_rna_best_hits.tsv",
        ok,
        "extendor_id",
    )

    # Gene summary + md from filtered main RNA best hits
    rna_hits = pl.read_csv(out_rna / "ge2_rna_best_hits.tsv", separator="\t")
    gsum = build_gene_summary_from_best_hits(rna_hits)
    write_gene_summary_tsv(gsum, out_rna / "ge2_gene_summary.tsv")
    gene_counts = {r["gene_name"]: r["n_extendors"] for r in gsum.iter_rows(named=True)}

    # Human genome TSV
    filter_tsv_by_extendor_col(
        hg / "ge2_human_genome_best_hits.tsv",
        out_hg / "ge2_human_genome_best_hits.tsv",
        ok,
        "extendor_id",
    )
    n_no_h = filter_tsv_by_extendor_col(
        hg / "ge2_no_human_hits.tsv",
        out_hg / "ge2_no_human_hits.tsv",
        ok,
        "extendor",
    )
    print(f"Filtered ge2_no_human_hits rows: {n_no_h}")

    # RNA no-hit list + BLAST tabular (qseqid has extendor|meta)
    filter_rna_no_hit_txt(rna / "ge2_rna_no_hit.txt", out_rna / "ge2_rna_no_hit.txt", ok)
    filter_blast_tab_no_header(
        rna / "refseq_rna" / "ge2_refseq_rna_human.txt",
        out_rna / "refseq_rna" / "ge2_refseq_rna_human.txt",
        ok,
        lambda ln: extendor_from_blast_qseqid_cell(ln.split("\t", 1)[0]),
    )
    filter_blast_tab_no_header(
        rna / "refseq_rna_select" / "ge2_refseq_select_rna_human.txt",
        out_rna / "refseq_rna_select" / "ge2_refseq_select_rna_human.txt",
        ok,
        lambda ln: extendor_from_blast_qseqid_cell(ln.split("\t", 1)[0]),
    )

    write_filtered_biology_summary(
        rna / "top_gene_biology_summary.md",
        out_rna / "top_gene_biology_summary.md",
        gsum,
        gene_counts,
        n_ge2_no_human_filtered=n_no_h,
    )

    if not args.skip_large_genome:
        src_g = hg / "ge2_human_genome.txt"
        dst_g = out_hg / "ge2_human_genome.txt"
        print(f"Streaming filter (large): {src_g} -> {dst_g}")
        k, t = filter_large_blast_tabular(src_g, dst_g, ok)
        print(f"ge2_human_genome.txt: kept {k:,} / {t:,} lines")
    else:
        print("Skipped ge2_human_genome.txt (--skipLargeGenome)")

    write_filtering_notes(args.outRoot / "FILTERING_NOTES.md")
    print(f"Wrote filtered BLAST bundle under {args.outRoot}")


if __name__ == "__main__":
    main()
