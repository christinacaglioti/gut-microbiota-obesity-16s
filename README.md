# gut-microbiota-obesity-16s
Reanalysis of public 16S rRNA amplicon data (NCBI BioProject PRJNA962831): gut microbiota in adults with obesity vs. normal weight
## Data

This project reanalyses publicly available 16S rRNA gene amplicon data from:

> Ahmad, M. A., Karavetian, M., Moubareck, C. A., Wazz, G., Mahdy, T., & Venema, K. (2023).
> Association of the gut microbiota with clinical variables in obese and lean Emirati subjects.
> *Frontiers in Microbiology*, 14, 1182460. https://doi.org/10.3389/fmicb.2023.1182460

| | |
|---|---|
| Repository | NCBI Sequence Read Archive (SRA) |
| BioProject | [PRJNA962831](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA962831) |
| Samples | 74 stool samples (adults, United Arab Emirates) |
| Target | 16S rRNA gene, V3–V4 region |
| Platform | Illumina MiSeq, 2 × 300 bp paired-end |

Raw reads and SRA metadata were retrieved with the QIIME 2 plugin
[q2-fondue](https://github.com/bokulich-lab/q2-fondue) (Ziemski et al., 2022).

### Group labels (obese vs. lean)

The public SRA metadata do **not** include an explicit obesity status or BMI column.
Group labels were therefore **inferred from the sample codes** (`Isolate` field):

| Sample code pattern | n | Assigned group |
|---|---|---|
| `HMG…A`, `QH…A` (suffix "A") | 43 | obese |
| `H…`, `ZU…` (no suffix) | 31 | lean |

This 43/31 split matches exactly the group sizes reported in the original publication
(43 obese, 31 lean). The suffix "A" is consistent with the pre-surgery time point of the
bariatric pre–post study the samples come from. However, this assignment is an
**inference, not an official label**, and has not been confirmed by the original authors.

### Limitations of the public metadata

- No BMI, waist circumference, diet or medication data are available in SRA, so
  analyses adjusting for these covariates (or correlating taxa with BMI) are not possible.
- Large files (raw reads, QIIME 2 artifacts) are not stored in this repository
  (see `.gitignore`); they can be regenerated with the scripts in `scripts/`.