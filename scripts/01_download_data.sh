#!/usr/bin/env bash
# ------------------------------------------------------------------
# 01_download_data.sh
# Download 16S rRNA amplicon reads and SRA metadata for
# BioProject PRJNA962831 (Ahmad et al., 2023) using q2-fondue.
#
# Usage (inside the QIIME 2 conda environment):
#   conda activate rachis-qiime2-2026.10
#   bash scripts/01_download_data.sh your.email@example.com
# ------------------------------------------------------------------
set -euo pipefail

EMAIL="${1:?Please provide your e-mail address (required by NCBI)}"
BIOPROJECT="PRJNA962831"

mkdir -p data metadata

# 1. Accession list -> QIIME 2 artifact
printf "ID\n%s\n" "$BIOPROJECT" > metadata/ids.tsv
qiime tools import \
  --type NCBIAccessionIDs \
  --input-path metadata/ids.tsv \
  --output-path data/ids.qza

# 2. SRA metadata (one row per sequencing run)
qiime fondue get-metadata \
  --i-accession-ids data/ids.qza \
  --p-email "$EMAIL" \
  --o-metadata data/sra_metadata.qza \
  --o-failed-runs data/failed_metadata.qza \
  --verbose

qiime tools export \
  --input-path data/sra_metadata.qza \
  --output-path data/sra_metadata

# 3. Raw paired-end reads (FASTQ) for all runs
qiime fondue get-sequences \
  --i-accession-ids data/ids.qza \
  --p-email "$EMAIL" \
  --p-retries 3 \
  --o-single-reads data/reads_single.qza \
  --o-paired-reads data/reads_paired.qza \
  --o-failed-runs data/failed_sequences.qza \
  --verbose