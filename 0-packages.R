## Functions
# Kaizad F. Patel
# September 2019

## packages ####
library(readxl)
library(ggplot2)       # 2.1.0
library(dplyr)         # 0.5.0
library(readr)         # 1.0.0
library(lubridate)     # 1.6.0
library(stringr)       # 1.1.0
library(luzlogr)       # 0.2.0
library(tidyr)
library(readr)
library(tidyverse)
library(dplyr)
library(Rmisc)
library(ggplot2)
library(data.table)
library(cowplot)
library(qwraps2)
library(knitr)
library(reshape2)
library(ggalt)
library(ggExtra)
library(stringi)
library(nlme)
library(car)
library(agricolae)
library(googlesheets)
library(gsheet)

# DATA_DIR               <- "data/"
# OUTPUT_DIR		         <- "outputs/"

## SET OUTPUT FILES
# from SCRIPT 1 
FTICR_META <- "fticr/fticr_meta.csv" # all metadata about formula, etc. assignment for each m/z value
HCOC <- "fticr/meta_hcoc.csv"
ELEMENTS <-  "fticr/meta_elements.csv"
CLASS <- "fticr/meta_class.csv"


FTICR_MASTER_LONG <- "fticr/fticr_master_long.csv" #
FTICR_RAWMASTER_LONG <- "fticr/fticr_rawmaster_long.csv"
FTICR_FENTON <- "fticr/fticr_fenton.csv" # pre- and post-Fenton data, intensities only
FTICR_GOETHITE <- "fticr/fticr_goethite.csv" # pre- and post-Goethite data, intensities only

# from SCRIPT 2
PERCENTILE <- "fticr/fticr_raw_percentile.csv" # relative intensities and percentile classification for native and Fenton extracts
RELATIVE_ABUND <- "fticr/fticr_relative_abundance.csv"
SUMMARY_ELEMENTS <- "fticr/fticr_summary_elements.csv"
GOETHITE_ADSORPTION <- "fticr/fticr_goethite_sorbed.csv"
