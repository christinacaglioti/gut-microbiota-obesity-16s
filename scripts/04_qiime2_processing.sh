#!/usr/bin/env bash
# ------------------------------------------------------------------
# 04_qiime2_processing.sh
# Secondary analysis in QIIME 2: quality control, primer removal,
# denoising (DADA2), taxonomy (SILVA) and phylogeny.
#
# Usage (inside the QIIME 2 conda environment):
#   conda activate rachis-qiime2-2026.10
#   bash scripts/04_qiime2_processing.sh
# ------------------------------------------------------------------
set -euo pipefail

mkdir -p results

# 1. Quality overview of the raw reads ------------------------------
#    Result: 74 samples, 40,221-357,170 read pairs per sample
#    (median ~129,800), reads 2 x 300 bp. Quality drops after
#    ~250 bp in R1 and ~210 bp in R2.
qiime demux summarize \
  --i-data data/reads_paired.qza \
  --o-visualization results/demux_summary.qzv

# 2. Primer removal (cutadapt) --------------------------------------
#    Primers confirmed with scripts/03_check_primers.py.
#    341F: CCTACGGGNGGCWGCAG / 805R: GACTACHVGGGTATCTAATCC
#    --p-match-read-wildcards: R1 reads often start with "N".
#    --p-discard-untrimmed: drop reads without the expected primer.
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences data/reads_paired.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-match-read-wildcards \
  --p-discard-untrimmed \
  --p-cores 2 \
  --o-trimmed-sequences data/reads_trimmed.qza \
  --o-stats data/cutadapt_stats.qza

qiime demux summarize \
  --i-data data/reads_trimmed.qza \
  --o-visualization results/demux_trimmed_summary.qzv

# 3. Denoising with DADA2 -------------------------------------------
#    TODO: truncation lengths chosen from demux_trimmed_summary.qzv

# 4. Taxonomic classification (SILVA 138.2) -------------------------
#    TODO

# 5. Phylogenetic tree (MAFFT + FastTree) ---------------------------
#    TODO