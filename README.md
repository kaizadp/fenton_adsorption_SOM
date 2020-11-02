## Reactive oxygen species alter chemical composition and adsorptive fractionation of soil-derived organic matter

Kaizad F. Patel, Václav Tejnecký, Tsutomu Ohno, Vanessa L. Bailey, Rachel L. Sleighter, and Patrick G. Hatcher

This repository contains data files and processing/visualization scripts for the published article, Patel et al. 2020 Geoderma doi: 
[10.1016/j.geoderma.2020.114805](https://doi.org/10.1016/j.geoderma.2020.114805)

Data are archived at the Environmental Data Initiative, doi: XXX-XXX

---

## Directory Structure

```
home
|---- code/
|        |---- 0.packages.R
|        |---- 1-preprocessing.R
|        |---- 2-fticr_abundances.R
|        |---- 3-doc_conc.R
|        
|---- data/
|        |---- data_for_EDI/
|        |       |---- fenton_adsorption_weoc.csv
|        |       |---- PatelKF_fenton_adsorption_metadata_20201022.docx
|        |       |---- postfenton_hardwood.csv
|        |       |---- postfenton_softwood.csv
|        |       |---- prefenton_hardwood.csv
|        |       |---- prefenton_softwood.csv
|        |       |---- soil_key.csv
|        |---- fticr/
|        |       |---- PostFenHWAdsorp.csv
|        |       |---- PostFentonSWAdsorp
|        |       |---- PreFentonHWAdsorp
|        |       |---- PreFentonSWAdsorp
|        |---- doc_concentrations.xlsx
|        |---- soil_key.csv
|
|---- fticr/ (processed output FTICR files, used for further data analysis/visualization)
|        |---- fticr_fenton_loss.csv
|        |---- fticr_fenton.csv
|        |---- fticr_goethite_relabund.csv
|        |---- fticr_goethite_sorbed.csv
|        |---- fticr_goethite.csv
|        |---- fticr_master_hcoc.csv
|        |---- fticr_master_long.csv
|        |---- fticr_meta.csv
|        |---- fticr_raw_percentile.csv
|        |---- fticr_rawmaster_long.csv
|        |---- fticr_relative_abundance.csv
|        |---- fticr_summary_elements.csv
|        |---- meta_class.csv
|        |---- meta_elements.csv
|        |---- meta_hcoc.csv
|
|---- reports-markdown/
|        |---- images/
|        |---- markdown_report.md
|        |---- markdown_report.Rmd
|
|---- adsorptive_fractionation_som.Rproj
|---- README.md
```