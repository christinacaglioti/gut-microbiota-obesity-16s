"""
02_make_metadata.py
Build the QIIME 2 sample metadata file from the SRA metadata.

The public SRA metadata of PRJNA962831 contain no explicit obesity status.
Groups are inferred from the sample codes ("Isolate" field):
  - codes ending in "A" (HMG..A, QH..A) -> obese  (n = 43)
  - codes without suffix (H.., ZU..)    -> lean   (n = 31)
This 43/31 split matches the original publication (Ahmad et al., 2023).

Usage:
  python scripts/02_make_metadata.py
"""
import pandas as pd

SRA_METADATA = "data/sra_metadata/sra-metadata.tsv"
OUTPUT = "metadata/sample-metadata.tsv"

df = pd.read_csv(SRA_METADATA, sep="\t")
isolate = df["Isolate [sample]"]

meta = pd.DataFrame({
    "sample-id": df["ID"],
    "isolate": isolate,
    "site_code": isolate.str.extract(r"^([A-Za-z]+)")[0],
    "group": isolate.str.endswith("A").map({True: "obese", False: "lean"}),
}).sort_values("sample-id")

# Sanity check: group sizes must match the publication
counts = meta["group"].value_counts()
assert counts.get("obese") == 43 and counts.get("lean") == 31, counts

meta.to_csv(OUTPUT, sep="\t", index=False)
print(f"Wrote {OUTPUT}")
print(counts)