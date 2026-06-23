library(stringdist)
library(Biostrings)
library(stringr)
library(GenomicAlignments)
library(data.table)

riffle <- function(a, b) {
  seqmlab <- seq(length = length(a))
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab])
}

##################################################
############ INPUT ARGUMENTS #####################
##################################################
args <- commandArgs(trailingOnly = TRUE)
output_directory       = args[1]   # directory for writing output files (must end with /)
compactor_file         = args[2]   # path to SPLASH2 compactors output (after_correction.scores.top_effect_size_bin.tsv)
STAR_reference         = args[3]   # path to STAR index directory
bowtie2_reference      = args[4]   # path to bowtie2 index for the reference genome (prefix)
annotated_splice_juncs_file    = args[5]   # annotated splice junctions file
annotated_exon_boundaries_file = args[6]   # annotated exon boundaries file
gene_coords_file               = args[7]   # gene coordinates file
anchor_length          = as.numeric(args[8])  # anchor k-mer length (25 for SPLASH2)

if (is.na(anchor_length)) anchor_length <- 25

cat("Output directory:", output_directory, "\n")
cat("Compactor file:", compactor_file, "\n")
cat("Anchor length:", anchor_length, "\n")

setwd(output_directory)
system(paste("mkdir -p ", output_directory, "STAR_alignment", sep = ""))
system(paste("mkdir -p ", output_directory, "Bowtie_alignment", sep = ""))

##################################################
##### STEP 1: Read and reshape SPLASH2 compactors
##################################################

compactors_dt = fread(compactor_file)
cat("Read", nrow(compactors_dt), "compactor rows for", length(unique(compactors_dt$anchor)), "unique anchors\n")

setnames(compactors_dt, "compactor", "compactor_seq")

compactors_dt = compactors_dt[!duplicated(paste(anchor, compactor_seq, sep = "--"))]

compactors_dt[, anchor_total_support := sum(exact_support), by = anchor]
setorder(compactors_dt, anchor, -exact_support)
compactors_dt[, num_compactor_per_anchor := .N, by = anchor]
compactors_dt = compactors_dt[num_compactor_per_anchor > 1]
compactors_dt[, compactor_order := 1:.N, by = anchor]

anchor_rank_dt = data.table(unique(compactors_dt$anchor), 1:length(unique(compactors_dt$anchor)))
names(anchor_rank_dt) = c("anchor", "anchor_index")
compactors_dt = merge(compactors_dt, anchor_rank_dt, all.x = TRUE, all.y = FALSE, by = "anchor")

cat("After filtering (>1 compactor per anchor):", nrow(compactors_dt), "rows,",
    length(unique(compactors_dt$anchor)), "anchors\n")

##################################################
##### STEP 2: Compute distances between top 2 compactors
##################################################

compactors_dt_high_rank = compactors_dt[compactor_order < 3]

compactors_dt_high_rank[, length_compactor := nchar(compactor_seq)]
compactors_dt_high_rank[, min_length_compactor := min(length_compactor), by = anchor]
compactors_dt_high_rank[, compactor_trimmed := substr(compactor_seq, 1, min_length_compactor)]
compactors_dt_high_rank[, ex_anchor_part := substr(compactor_trimmed, anchor_length + 1, nchar(compactor_trimmed))]

compactors_dt_high_rank_shifted = compactors_dt_high_rank[, data.table::shift(.SD, 1, NA, "lead", TRUE), .SDcols = 1:ncol(compactors_dt_high_rank)]
compactors_dt_high_rank = cbind(compactors_dt_high_rank, compactors_dt_high_rank_shifted[, list(anchor_lead_1, ex_anchor_part_lead_1)])
compactors_dt_high_rank = compactors_dt_high_rank[anchor_lead_1 == anchor]

cat("Computing Levenshtein and Hamming distances for", nrow(compactors_dt_high_rank), "anchor pairs...\n")

compactors_dt_high_rank[, lev_dist := stringdist(ex_anchor_part, ex_anchor_part_lead_1, method = "lv"), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[, ham_dist := stringdist(ex_anchor_part, ex_anchor_part_lead_1, method = "hamming"), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[, lev_operations := as.character(attributes(adist(ex_anchor_part, ex_anchor_part_lead_1, count = T))$trafos), by = 1:nrow(compactors_dt_high_rank)]

safe_max_run <- function(ops_string, char) {
  chars <- strsplit(ops_string, "")[[1]]
  r <- rle(chars == char)
  true_runs <- r$lengths[r$values == TRUE]
  if (length(true_runs) == 0) return(0)
  return(max(true_runs))
}

compactors_dt_high_rank[, run_length_D := safe_max_run(lev_operations, "D"), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[, run_length_I := safe_max_run(lev_operations, "I"), by = 1:nrow(compactors_dt_high_rank)]

compactors_dt = merge(compactors_dt,
                      unique(compactors_dt_high_rank[, list(anchor, ham_dist, lev_dist, lev_operations, run_length_D, run_length_I)]),
                      all.x = TRUE, all.y = FALSE, by = "anchor")

compactors_dt[, anchor_event := ""]
compactors_dt[ham_dist == lev_dist, anchor_event := paste("Base_pair_change_", ham_dist, sep = "")]

cat("Distance computation complete\n")

##################################################
##### STEP 3: STAR alignment of compactors to genome
##################################################

compactors_dt[, compactor_index := paste(">", anchor_index, "_", compactor_order, sep = "")]

compactors_fasta = riffle(compactors_dt$compactor_index, compactors_dt$compactor_seq)
compactors_fasta = data.table(compactors_fasta)
write.table(compactors_fasta, paste(output_directory, "compactors_fasta.fa", sep = ""),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

cat("Running STAR alignment...\n")
system(paste("STAR --runThreadN 4 --genomeDir ", STAR_reference,
             " --readFilesIn ", output_directory, "compactors_fasta.fa",
             " --outFileNamePrefix ", output_directory, "STAR_alignment/Compactors",
             " --twopassMode Basic --alignIntronMax 1000000",
             " --chimJunctionOverhangMin 10 --chimSegmentReadGapMax 0",
             " --chimOutJunctionFormat 1 --chimSegmentMin 12",
             " --chimScoreJunctionNonGTAG -4 --chimNonchimScoreDropMin 10",
             " --outSAMtype SAM --chimOutType SeparateSAMold",
             " --outSAMunmapped None --clip3pAdapterSeq AAAAAAAAA",
             " --outSAMattributes NH HI AS nM NM", sep = ""))

alignment_info_compactors = fread(paste(output_directory, "STAR_alignment/CompactorsAligned.out.sam", sep = ""),
                                  header = FALSE, skip = "NH:")
alignment_info_compactors[, V1 := paste(">", V1, sep = "")]
alignment_info_compactors[, num_alignments := .N, by = V1]
alignment_info_compactors[, STAR_num_mismatches := as.numeric(strsplit(V16, split = ":")[[1]][3]), by = V16]

compactors_dt = merge(compactors_dt,
                      alignment_info_compactors[!duplicated(V1), list(V1, V2, V3, V4, V6, num_alignments, STAR_num_mismatches)],
                      all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "V1")
setnames(compactors_dt, c("V2", "V3", "V4", "V6", "num_alignments"),
         c("STAR_flag", "STAR_chr", "STAR_coord", "STAR_CIGAR", "STAR_num_alignments"))

compactors_dt[, compactor_index := gsub(">", "", compactor_index)]
compactors_dt[, is.aligned_STAR := 0]
compactors_dt[, is.STAR_chimeric := 0]
compactors_dt[, is.STAR_SJ := 0]

num_chimeric_alignments = fread(paste(output_directory, "STAR_alignment/CompactorsLog.final.out", sep = ""),
                                sep = "|", skip = 35)
num_chimeric_alignments[, V2 := gsub("\t", "", V2)]
if (num_chimeric_alignments$V2[1] != 0) {
  chimeric_alignment_info = fread(paste(output_directory, "STAR_alignment/CompactorsChimeric.out.sam", sep = ""),
                                  header = FALSE, skip = "NH:")
  compactors_dt[compactor_index %in% chimeric_alignment_info$V1, is.STAR_chimeric := 1]
}
compactors_dt[STAR_CIGAR %like% "N", is.STAR_SJ := 1]
compactors_dt[!is.na(STAR_chr) & (is.STAR_chimeric == 0), is.aligned_STAR := 1]

cat("STAR alignment complete\n")

##################################################
##### STEP 4: Gene name assignment
##################################################

system(paste("samtools view -S -b ", output_directory, "STAR_alignment/CompactorsAligned.out.sam > ",
             output_directory, "STAR_alignment/CompactorsAligned.out.bam", sep = ""))
system(paste("bedtools bamtobed -split -i ", output_directory, "STAR_alignment/CompactorsAligned.out.bam",
             " | sed '/^chr/!d' | sort -k1,1 -k2,2n > ",
             output_directory, "STAR_alignment/called_exons.bed", sep = ""))
system(paste("bedtools intersect -a ", output_directory, "STAR_alignment/called_exons.bed",
             " -b ", gene_coords_file,
             " -wb -loj | cut -f 4,10 | bedtools groupby -g 1 -c 2 -o distinct > ",
             output_directory, "STAR_alignment/compactor_genes.txt", sep = ""))

compactor_genes = fread(paste(output_directory, "STAR_alignment/compactor_genes.txt", sep = ""),
                        sep = "\t", header = FALSE)
names(compactor_genes) = c("compactor_index", "compactor_gene")
compactors_dt = merge(compactors_dt, compactor_genes[!duplicated(compactor_index)],
                      all.x = TRUE, all.y = FALSE, by = "compactor_index")
compactors_dt[compactor_gene == ".", compactor_gene := NA]
compactors_dt[, num_compactor_gene_anchor := length(unique(compactor_gene)), by = anchor]

cat("Gene assignment complete\n")

##################################################
##### STEP 5: Extract and annotate splice junctions
##################################################

alignment_info_compactors = fread(paste(output_directory, "STAR_alignment/CompactorsAligned.out.sam", sep = ""),
                                  header = FALSE, skip = "NH:")
alignment_info_compactors = alignment_info_compactors[!duplicated(V1)][V6 %like% "N"]

if (nrow(alignment_info_compactors) > 0) {
  write.table(alignment_info_compactors,
              paste(output_directory, "STAR_alignment/top_splice_alignments.out.sam", sep = ""),
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  system(paste("samtools view -H ", output_directory, "STAR_alignment/CompactorsAligned.out.bam > ",
               output_directory, "STAR_alignment/sam_header.txt", sep = ""))
  system(paste("cat ", output_directory, "STAR_alignment/sam_header.txt ",
               output_directory, "STAR_alignment/top_splice_alignments.out.sam > ",
               output_directory, "STAR_alignment/top_splice_alignments_with_header.out.sam", sep = ""))
  system(paste("samtools view -S -b ", output_directory, "STAR_alignment/top_splice_alignments_with_header.out.sam > ",
               output_directory, "STAR_alignment/top_splice_alignments_with_header.out.bam", sep = ""))
  system(paste("bedtools bamtobed -split -i ", output_directory, "STAR_alignment/top_splice_alignments_with_header.out.bam > ",
               output_directory, "STAR_alignment/extracted_splice_junction.bed", sep = ""))

  splice_junctions = fread(paste(output_directory, "STAR_alignment/extracted_splice_junction.bed", sep = ""), header = FALSE)
  splice_junctions_shifted = splice_junctions[, data.table::shift(.SD, 1, NA, "lead", TRUE), .SDcols = 1:6]
  splice_junctions = cbind(splice_junctions, splice_junctions_shifted)
  splice_junctions = splice_junctions[V4 == V4_lead_1]
  splice_junctions[, splice_junc := paste(V1, ":", V3, ":", V2_lead_1 + 1, sep = "")]
  splice_junctions[, all_splice_juncs := paste(splice_junc, collapse = "--"), by = V4]

  ##### annotating AS splice sites
  known_splice_sites = fread(annotated_splice_juncs_file)
  known_splice_sites[, chr_V2 := paste(V1, V2, sep = ":")]
  known_splice_sites[, chr_V3 := paste(V1, V3, sep = ":")]
  known_splice_sites[, num_uniq_V2_for_V3 := length(unique(chr_V2)), by = chr_V3]
  known_splice_sites[, num_uniq_V3_for_V2 := length(unique(chr_V3)), by = chr_V2]
  alt_v2 = known_splice_sites[num_uniq_V3_for_V2 > 1]$V2
  alt_v3 = known_splice_sites[num_uniq_V2_for_V3 > 1]$V3
  total = c(alt_v2, alt_v3, alt_v2 - 1, alt_v3 - 1, alt_v2 + 1, alt_v3 + 1)
  splice_junctions[, SSA_AS_annot := 0]
  splice_junctions[, SSB_AS_annot := 0]
  splice_junctions[V3 %in% total, SSA_AS_annot := 1]
  splice_junctions[V2_lead_1 %in% total, SSB_AS_annot := 1]
  splice_junctions[, SS_AS_annot := paste(SSA_AS_annot, ":", SSB_AS_annot, sep = "")]
  splice_junctions[, all_SS_AS_annot := paste(SS_AS_annot, collapse = "--"), by = V4]

  ##### annotating exon boundaries
  known_exon_boundaries = fread(annotated_exon_boundaries_file, header = TRUE, sep = "\t")
  known_exon_boundaries = known_exon_boundaries[!duplicated(paste(V1, V2, V3))]
  total = c(known_exon_boundaries$chr_V2, known_exon_boundaries$chr_V2_1, known_exon_boundaries$chr_V2_2,
            known_exon_boundaries$chr_V3, known_exon_boundaries$chr_V3_1, known_exon_boundaries$chr_V3_2)
  splice_junctions[, SSA_annot := 0]
  splice_junctions[, SSB_annot := 0]
  splice_junctions[paste(V1, V3, sep = "") %in% total, SSA_annot := 1]
  splice_junctions[paste(V1, V2_lead_1, sep = "") %in% total, SSB_annot := 1]
  splice_junctions[, SS_annot := paste(SSA_annot, ":", SSB_annot, sep = "")]
  splice_junctions[, all_SS_annot := paste(SS_annot, collapse = "--"), by = V4]

  splice_junctions = splice_junctions[!duplicated(V4)]
  compactors_dt = merge(compactors_dt,
                        splice_junctions[, list(V4, all_splice_juncs, all_SS_AS_annot, all_SS_annot)],
                        all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "V4")
  cat("Splice junction annotation complete\n")
} else {
  compactors_dt[, all_splice_juncs := NA_character_]
  compactors_dt[, all_SS_AS_annot := NA_character_]
  compactors_dt[, all_SS_annot := NA_character_]
  cat("No splice junctions found in STAR alignment\n")
}

##################################################
##### STEP 6: Bowtie2 duplicate anchor removal
##################################################

anchors_fasta = riffle(paste(">", compactors_dt[!duplicated(anchor)]$anchor_index, sep = ""),
                       compactors_dt[!duplicated(anchor)]$anchor)
write.table(anchors_fasta, paste(output_directory, "anchors_fasta.fa", sep = ""),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

cat("Building Bowtie2 index of compactors and aligning anchors...\n")
system(paste("bowtie2-build ", output_directory, "compactors_fasta.fa ",
             output_directory, "Bowtie_alignment/compactors_bowtie2_index", sep = ""))
system(paste("bowtie2 -f -k 10 -x ", output_directory, "Bowtie_alignment/compactors_bowtie2_index",
             " -U ", output_directory, "anchors_fasta.fa",
             " -S ", output_directory, "Bowtie_alignment/anchors_to_compactors_bowtie.sam", sep = ""))

Bowtie_alignment_info = fread(cmd = paste0("grep -v '^@' ",
                                            output_directory, "Bowtie_alignment/anchors_to_compactors_bowtie.sam"),
                              header = FALSE, fill = TRUE)
Bowtie_alignment_info = Bowtie_alignment_info[V12 == "AS:i:0"]
Bowtie_alignment_info[, V1 := as.numeric(V1)]  # anchor index is numeric for comparison
Bowtie_alignment_info[, aligned_compactor_anchor_index := as.numeric(strsplit(V3, split = "_")[[1]][1]), by = V3]
Bowtie_alignment_info$aligned_compactor_anchor_index = as.numeric(Bowtie_alignment_info$aligned_compactor_anchor_index)
Bowtie_alignment_info[, min_aligned_compactor_anchor_index := min(aligned_compactor_anchor_index), by = V1]
Bowtie_alignment_info[, duplicate_anchor := 0]
Bowtie_alignment_info[min_aligned_compactor_anchor_index < V1, duplicate_anchor := 1]
Bowtie_alignment_info[duplicate_anchor == 1 & aligned_compactor_anchor_index == min_aligned_compactor_anchor_index,
                      high_ranked_compactors := paste(V3, collapse = ":"), by = V1]
Bowtie_alignment_info = Bowtie_alignment_info[duplicate_anchor == 1 & !is.na(high_ranked_compactors)]
Bowtie_alignment_info = Bowtie_alignment_info[!duplicated(V1)]
compactors_dt = merge(compactors_dt,
                      Bowtie_alignment_info[, list(V1, duplicate_anchor, high_ranked_compactors)],
                      all.x = TRUE, all.y = FALSE, by.x = "anchor_index", by.y = "V1")
compactors_dt[is.na(duplicate_anchor), duplicate_anchor := 0]

cat("Duplicate anchor removal complete\n")

##################################################
##### STEP 7: Final classification
##################################################

compactors_dt[, is_STAR_SJ_in_top_two := 0]
compactors_dt[(compactor_order < 3), is_STAR_SJ_in_top_two := sum(is.STAR_SJ), by = anchor]
compactors_dt[, is_STAR_SJ_in_top_two := max(is_STAR_SJ_in_top_two), by = anchor]

compactors_dt[, is_aligned_STAR_in_top_two := 0]
compactors_dt[(compactor_order < 3), is_aligned_STAR_in_top_two := sum(is.aligned_STAR), by = anchor]
compactors_dt[, is_aligned_STAR_in_top_two := max(is_aligned_STAR_in_top_two), by = anchor]

compactors_dt[is_STAR_SJ_in_top_two > 0 & (ham_dist != lev_dist | ham_dist > 5) & num_compactor_gene_anchor == 1,
              anchor_event := "Splicing"]

compactors_dt[(lev_dist < run_length_D + run_length_I + 1) & (run_length_D > 1 | run_length_I > 1),
              anchor_event := "Internal_splicing"]
compactors_dt[, help := 2 * pmax(run_length_D, run_length_I) + 1]
compactors_dt[lev_dist < help & (run_length_D > 1 | run_length_I > 1),
              anchor_event := "Internal_splicing"]
compactors_dt[, help := NULL]

compactors_dt[!is.na(STAR_CIGAR),
              num_STAR_match := sum(explodeCigarOpLengths(STAR_CIGAR, ops = c("M"))[[1]]),
              by = STAR_CIGAR]

cat("\nClassification summary:\n")
print(compactors_dt[!duplicated(anchor), .N, by = anchor_event][order(-N)])

##################################################
##### STEP 8: Write output
##################################################

write.table(compactors_dt, paste(output_directory, "classified_compactors.tsv", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nOutput written to:", paste(output_directory, "classified_compactors.tsv", sep = ""), "\n")
cat("Done.\n")
