### SOM adsorptive fractionation
### script for analyzing processed FTICR data

### Kaizad F. Patel
### September 2019

# ------------------------------------------
# ------------------------------------------

source("0-packages.R")


## Step 1: load file from Google Drive

fticr_data <- function() {
  #  https://drive.google.com/file/d/1rcm6i8J6KQTeXFgcBqapfbQMKjkoHHNz/view?usp=sharing    
  key <- gs_key("1rcm6i8J6KQTeXFgcBqapfbQMKjkoHHNz")
    old <- getwd()
    setwd("data/")
    on.exit(setwd(old))
    gs_download(key, overwrite = TRUE)
}

fticr_data = gs_download(gs_key("1rcm6i8J6KQTeXFgcBqapfbQMKjkoHHNz"))

