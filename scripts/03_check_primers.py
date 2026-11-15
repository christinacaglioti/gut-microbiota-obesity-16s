"""
03_check_primers.py
Check whether the Illumina V3-V4 primers (341F / 805R) are still present at
the start of the reads, using the first reads of one sample.

Result for PRJNA962831: primers are present in both R1 and R2.
Note: the first base of most R1 reads is an "N" (uncalled base in the first
sequencing cycle), so the forward primer only matches when read wildcards are
allowed -> use --p-match-read-wildcards in cutadapt.

Usage:
  python scripts/03_check_primers.py
"""
import gzip
import re
import zipfile

READS = "data/reads_paired.qza"
N_READS = 5000

# Primers with IUPAC wildcards written as regular expressions
FWD = re.compile(r"[CN]CTACGG[ACGTN]GGC[ATN]GCAG")       # 341F: CCTACGGGNGGCWGCAG
REV = re.compile(r"GACTAC[ACTN][ACGN]GGGTATCTAATCC")     # 805R: GACTACHVGGGTATCTAATCC

z = zipfile.ZipFile(READS)
files = sorted(f for f in z.namelist() if f.endswith(".fastq.gz"))
r1 = next(f for f in files if "_R1_" in f)
r2 = r1.replace("_R1_", "_R2_")


def check(name, primer):
    seqs = []
    with z.open(name) as fh, gzip.open(fh, "rt") as g:
        for i, line in enumerate(g):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= N_READS:
                    break
    hits = sum(bool(primer.match(s)) for s in seqs)
    print(f"{name.split('/')[-1]}: primer at read start in "
          f"{hits}/{len(seqs)} reads ({100 * hits / len(seqs):.1f} %)")


check(r1, FWD)
check(r2, REV)