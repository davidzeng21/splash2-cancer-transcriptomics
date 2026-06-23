#!/usr/bin/env python3
import os
import random
import sys

DRY_RUN = "--dry-run" in sys.argv

WORKDIR = "/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched"
FASTQ_LIST = os.path.join(WORKDIR, "fastq_list.txt")
FASTQ_DIR = os.path.join(WORKDIR, "fastqs")
REMOVED_LIST = os.path.join(WORKDIR, "removed_samples.txt")
SEED = 42
N_KEEP = 107

with open(FASTQ_LIST) as f:
    lines = [l.strip() for l in f if l.strip()]

sample_ids = sorted(set(
    os.path.basename(l).replace("_normal.fastq.gz", "").replace("_tumor.fastq.gz", "")
    for l in lines
))

print(f"Total samples found: {len(sample_ids)}")
assert len(sample_ids) == 215, f"Expected 215 samples, got {len(sample_ids)}"

random.seed(SEED)
keep = set(random.sample(sample_ids, N_KEEP))
remove = set(sample_ids) - keep

if DRY_RUN:
    print("=== DRY RUN (no files will be modified) ===\n")

print(f"Keeping {len(keep)} samples ({len(keep)*2} fastq files)")
print(f"Removing {len(remove)} samples ({len(remove)*2} fastq files)")

print(f"\n--- KEEP ({len(keep)} samples) ---")
for sid in sorted(keep):
    print(f"  {sid}")

print(f"\n--- REMOVE ({len(remove)} samples) ---")
for sid in sorted(remove):
    print(f"  {sid}")

if DRY_RUN:
    print("\n=== DRY RUN complete. Re-run without --dry-run to execute. ===")
    sys.exit(0)

# Save removed samples list (with their fastq paths) before deleting
removed_lines = []
for sid in sorted(remove):
    removed_lines.append(f"fastqs/{sid}_normal.fastq.gz")
    removed_lines.append(f"fastqs/{sid}_tumor.fastq.gz")

with open(REMOVED_LIST, "w") as f:
    f.write("\n".join(removed_lines) + "\n")
print(f"\nSaved removed sample list to {REMOVED_LIST} ({len(removed_lines)} entries)")

# Delete the files
deleted = 0
for sid in sorted(remove):
    for suffix in ["_normal.fastq.gz", "_tumor.fastq.gz"]:
        fpath = os.path.join(FASTQ_DIR, sid + suffix)
        if os.path.exists(fpath):
            os.remove(fpath)
            deleted += 1
        else:
            print(f"  WARNING: {fpath} not found")

print(f"Deleted {deleted} files")

# Update fastq_list.txt with only kept samples
kept_lines = []
for sid in sorted(keep):
    kept_lines.append(f"fastqs/{sid}_normal.fastq.gz")
    kept_lines.append(f"fastqs/{sid}_tumor.fastq.gz")

with open(FASTQ_LIST, "w") as f:
    f.write("\n".join(kept_lines) + "\n")

print(f"Updated {FASTQ_LIST} with {len(kept_lines)} entries")

# Verify
remaining = os.listdir(FASTQ_DIR)
print(f"Files remaining in fastqs/: {len(remaining)}")
