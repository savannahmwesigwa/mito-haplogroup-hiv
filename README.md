# Mitochondrial haplogroups and HIV progression in children (CAfGEN)

Analysis code for:
Mwesigwa et al. Mitochondrial L2 Haplogroup and Disease Progression among Ugandan Children with Rapidly Progressive HIV Infection (in revision, Scientific Reports).

## Structure

```
R/00_setup.R                  packages, paths, data import, derived variables, helpers
R/01_manuscript_analysis.R    manuscript tables, figures, statistics; manuscript vs re-run check table
R/02_revision_analyses.R      additional analyses for revision
R/03_tables.R                 Tables 1-2 and Supplementary Tables S1-S7 as Word files
data/                         input data (not included, see Data availability)
output/                       figures and results (created on run)
```

## Requirements

- R >= 4.3
- Packages:

```r
install.packages(c("here", "readxl", "dplyr", "tidyr", "ggplot2",
                   "survival", "survminer", "waffle", "coxphf",
                   "flextable", "officer"))
```

If `waffle` is not available on CRAN for your R version:

```r
install.packages("remotes")
remotes::install_github("hrbrmstr/waffle")
```

## Input data

`data/Mito_haplo_data.xlsx`, sheet `Combined`, one row per participant.
Columns used:

| Column | Description |
|---|---|
| `Haplogroup_2_subclades` | major haplogroup (e.g. L2) |
| `Haplogroup_3_subclades` | subclade (e.g. L2b) |
| `RP_or_LTNP` | phenotype: RP, LTNP |
| `Gender` | F, M |
| `Country` | UGR (Uganda), BWR (Botswana) |
| `Sequencing` | WES, WGS |
| `Time_to_progression` | months from birth to date of progression |
| `Birth_date` | date of birth |
| `ART_inititation_date` | date of ART initiation |
| `Age_at_sampleCollection_Months` | age at enrollment (months) |
| `ART_duration_at_Collection_recalculated` | HAART duration at enrollment (months) |

## Run

From the project root:

```bash
Rscript R/01_manuscript_analysis.R > output/01_results.txt 2>&1
Rscript R/02_revision_analyses.R   > output/02_results.txt 2>&1
Rscript R/03_tables.R              > output/03_results.txt 2>&1
```

## Data availability

Participant-level data are not included in this repository.
Sequencing data: European Genome-phenome Archive (EGA), accession [EGAS...].
Clinical data: available from the CAfGEN consortium on request, subject to ethics approval.

## Citation

If you use this code, please cite the archived release (see `CITATION.cff`).

## License

MIT (see `LICENSE`)
