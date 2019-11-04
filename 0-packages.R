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

# create a custom ggplot theme
theme_kp <- function() {  # this for all the elements common across plots
  theme_bw() %+replace%
    theme(legend.position = "top",
          legend.key=element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'lines'),
          panel.border = element_rect(color="black",size=1.5, fill = NA),
          
          plot.title = element_text(hjust = 0.05, size = 14),
          axis.text = element_text(size = 14, face = "bold", color = "black"),
          axis.title = element_text(size = 14, face = "bold", color = "black"),
          
          # formatting for facets
          panel.background = element_blank(),
          strip.background = element_rect(colour="white", fill="white"), #facet formatting
          panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
          panel.spacing.y = unit(1.5, "lines"), #facet spacing for x axis
          strip.text.x = element_text(size=12, face="bold"), #facet labels
          strip.text.y = element_text(size=12, face="bold", angle = 270) #facet labels
    )
}
# create a custom ggplot function for Van Krevelen plots
gg_vankrev <- function(data,mapping){
  ggplot(data,mapping)+
# plot points
    geom_point(size=2, alpha = 0.5)+ # set size and transparency
# axis labels
    ylab("H/C")+
    xlab("O/C")+
# axis limits
    xlim(0,1.25)+
    ylim(0,2.5)+
# add boundary lines for Van Krevelen regions
    geom_segment(x = 0.0, y = 1.5, xend = 1.2, yend = 1.5,color="black",linetype="longdash")+
    geom_segment(x = 0.0, y = 2, xend = 1.2, yend = 2,color="black",linetype="longdash")+
    geom_segment(x = 0.0, y = 1, xend = 1.2, yend = 0.75,color="black",linetype="longdash")+
    geom_segment(x = 0.0, y = 0.7, xend = 1.2, yend = 0.5,color="black",linetype="longdash")
}

## to make the Van Krevelen plot:
# replace the initial `ggplot` function with `gg_vankrev` and use as normal

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
FENTON_LOSS <- "fticr/fticr_fenton_loss.csv"
GOETHITE_ADSORPTION <- "fticr/fticr_goethite_sorbed.csv"
GOETHITE_RELABUND <- "fticr/fticr_goethite_relabund.csv"
