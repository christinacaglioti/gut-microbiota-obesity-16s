#!/usr/bin/env bash
# ------------------------------------------------------------------
# 04_qiime2_processing.sh
# Secondary analysis in QIIME 2 (rachis 2026.10): quality control,
# primer removal, denoising (DADA2), taxonomy (SILVA 144), phylogeny,
# and export of processed data for the R analysis.
#
# Requirements:
#   - Run scripts/01_download_data.sh first (creates data/reads_paired.qza)
#   - Run scripts/02_make_metadata.py first (creates metadata/sample-metadata.tsv)
#   - Machine with >= 16 GB RAM (DADA2 failed on 8 GB, see step 3)
#
# Usage (inside the QIIME 2 conda environment):
#   conda activate rachis-qiime2-2026.10
#   bash scripts/04_qiime2_processing.sh
#
# Note: the full run takes several hours (DADA2 is the slowest step).
# ------------------------------------------------------------------
set -euo pipefail

mkdir -p results/tables results/processed

# 1. Quality overview of the raw reads ------------------------------
#    Result: 74 samples, 40,221-357,170 read pairs per sample
#    (median ~129,800; total 10.4 M), reads 2 x 300 bp.
#    Quality drops after ~250 bp in R1 and ~210 bp in R2.
qiime demux summarize \
  --i-data data/reads_paired.qza \
  --o-visualization results/demux_summary.qzv

# 2. Primer removal (cutadapt) --------------------------------------
#    Primers confirmed with scripts/03_check_primers.py.
#    341F: CCTACGGGNGGCWGCAG / 805R: GACTACHVGGGTATCTAATCC
#    --p-match-read-wildcards: most R1 reads start with "N".
#    --p-discard-untrimmed: drop reads without the expected primer.
#    Result: 97.7 % of read pairs retained (10,196,394);
#    read length after trimming: R1 = 284 nt, R2 = 279 nt.
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
#    Truncation chosen where median quality starts to drop:
#      R1 at 260 nt, R2 at 200 nt.
#    Overlap check: amplicon between primers ~426 nt
#      260 + 200 = 460 -> ~34 nt overlap (DADA2 needs >= 12).
#    maxEE = 2 (default) removes reads with too many expected errors.
#    Run on a 4-core / 16 GB machine with 2 threads: with 8 GB RAM,
#    DADA2 failed with "'names' attribute must be the same length as
#    the vector" (a memory error).
#    Result (mean over 74 samples): 80.7 % passed filter, 77.6 % merged,
#    73.8 % non-chimeric; 30,812-241,617 reads per sample; 4,662 ASVs.
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs data/reads_trimmed.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 200 \
  --p-n-threads 2 \
  --o-table data/table.qza \
  --o-representative-sequences data/rep-seqs.qza \
  --o-denoising-stats data/denoising-stats.qza \
  --o-base-transition-stats data/dada2-base-transition-stats.qza \
  --verbose

# Keep the per-sample denoising statistics as a small table
qiime tools export \
  --input-path data/denoising-stats.qza \
  --output-path data/denoising-stats
cp data/denoising-stats/stats.tsv results/tables/dada2_denoising_stats.tsv

# Feature table summary (reads per sample / per ASV)
qiime feature-table summarize \
  --i-table data/table.qza \
  --m-metadata-file metadata/sample-metadata.tsv \
  --output-dir results/table-summary

# 4. Taxonomic classification ---------------------------------------
#    Pre-trained SILVA 144 SSURef NR99 uniform classifier, V3V4
#    (341f-806r), provided by the SILVA team (trained with QIIME 2
#    2026.7; compatible with 2026.4 and later).
#    805R and 806R bind at nearly the same position -> same region.
#    Result: 88.1 % of ASVs classified to genus level.
SILVA_URL="https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V3V4-341f-806r/uniform"
SILVA_FILE="SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza"

wget -O data/silva-144-v3v4-classifier.qza "$SILVA_URL/$SILVA_FILE"
wget -O data/silva-144-v3v4-classifier.qza.md5 "$SILVA_URL/$SILVA_FILE.md5"

# Integrity check: the MD5 of the download must match SILVA's .md5 file
EXPECTED_MD5=$(cut -d " " -f1 data/silva-144-v3v4-classifier.qza.md5)
ACTUAL_MD5=$(md5sum data/silva-144-v3v4-classifier.qza | cut -d " " -f1)
if [ "$EXPECTED_MD5" != "$ACTUAL_MD5" ]; then
  echo "ERROR: MD5 mismatch for the SILVA classifier" >&2
  exit 1
fi

qiime feature-classifier classify-sklearn \
  --i-classifier data/silva-144-v3v4-classifier.qza \
  --i-reads data/rep-seqs.qza \
  --p-n-jobs 1 \
  --o-classification data/taxonomy.qza

# 5. Phylogenetic tree (MAFFT + FastTree) ---------------------------
#    Needed for phylogenetic diversity metrics (Faith PD, UniFrac).
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences data/rep-seqs.qza \
  --p-n-threads 2 \
  --output-dir data/phylogeny

# 6. Export processed data for the R analysis -----------------------
#    Small files tracked in the repository, so the downstream analysis
#    can be reproduced without re-running DADA2:
#      feature-table.biom, taxonomy.tsv, tree.nwk
qiime tools export --input-path data/table.qza --output-path results/processed
qiime tools export --input-path data/taxonomy.qza --output-path results/processed
qiime tools export --input-path data/phylogeny/rooted_tree.qza --output-path results/processed

echo "Done. Processed data in results/processed/"