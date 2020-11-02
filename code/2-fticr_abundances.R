#   Reactive oxygen species alter chemical composition and adsorptive fractionation of soil-derived organic matter
#   Kaizad Patel

#   FTICR-RELATIVE ABUNDANCES
#   use this script to calculate relative abundances and peak counts


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

source("code/0-packages.R")

# 0. input files ----------------------------------------------------------
meta = read.csv(FTICR_META)# <- "fticr/fticr_meta.csv" # all metadata about formula, etc. assignment for each m/z value
hcoc = read.csv(HCOC)
class = read.csv(CLASS)
elements = read.csv(ELEMENTS)
master = read.csv(FTICR_MASTER_LONG)# <- "fticr/fticr_master_long.csv" #
rawmaster = read.csv(FTICR_RAWMASTER_LONG)# <- "fticr/fticr_rawmaster_long.csv"
fenton = read.csv(FTICR_FENTON)# <- "fticr/fticr_fenton.csv" # pre- and post-Fenton data, intensities only
goethite = read.csv(FTICR_GOETHITE)# <- "fticr/fticr_goethite.csv" # pre- and post-Goethite data, intensities only

#
# ---------------------------------------------------------------------------- ----

# 1. create hcoc file for van krevelens ----
master_hcoc = 
  master %>% left_join(hcoc, by = "Mass")

### OUTPUT
write.csv(master_hcoc, MASTER_HCOC, row.names = FALSE)
# this file will be used for the Van Krevelen plots (preFenton and postFenton)

#
# ---------------------------------------------------------------------------- ----

# 2. RELATIVE ABUNDANCE FOR PRE- AND POST-FENTON GROUPS ----
# `rawmaster` is the longform master file. 

# summarizing by groups
rawmaster %>% 
  filter(Goethite == "PreGoethite") %>% # keep only pre-Goethite data. we don't want post-adsorption data
  group_by(Forest,Fenton,soil,Class) %>% 
  dplyr::summarize(compounds = n())%>% # this gives total COUNTS for each group
# in the same command, we will also create a column for total intensity
  ungroup() %>% # remove the previous grouping
  group_by(soil) %>% # create a new grouping 
  dplyr::mutate(total = sum(compounds)) %>%  # add a column for total intensity
# we can also calculate the relative abundance in the same command
  mutate(relabund = (compounds/total)*100) %>% # calculate relative abundance as a %
  mutate(relabund = round(relabund,2))-> # round to two decimal places
  raw_coregroups

# ^^^ this file has relative abundance of each group for each core

# now we need to summarize this for each treatment. combine all cores
raw_coregroups %>% 
  group_by(Forest, Fenton, Class) %>% 
  dplyr::summarise(rel_abund = mean(relabund)) %>% 
  dplyr::mutate(rel_abund = round(rel_abund,2))->
  raw_groups

### OUTPUT
write.csv(raw_groups,RELATIVE_ABUND, row.names = FALSE)

#

## Elements ----

master %>% 
  left_join(elements, by = "Mass") %>% 
  filter(Treatment=="PreFenton"| Treatment=="PostFenton") %>%  # choose only the preGoethite samples
  gather(element, el_count, C:P) %>% 
# replace all 0 by NA and then remove NA to help with calculations  
  na_if(.,"0") %>% # replace NA with 0
  na.omit() %>% 
  group_by(Forest, Treatment, element) %>% 
  dplyr::summarise(avg = mean(el_count)) %>% 
  ungroup %>% 
  dplyr::mutate(avg = round(avg,0))->
  master_el

### OUTPUT
write.csv(master_el, SUMMARY_ELEMENTS, row.names = FALSE)  
  


#
# ---------------------------------------------------------------------------- ----

# 2. PEAK COUNTS ----
## 2.1 INITIAL PEAK COUNTS ----
# we want to determine the total peaks in HW vs. SW
# as well as the number of peaks in each group type

# rawmaster

rawmaster %>% 
  filter(Treatment=="PreFenton") %>% 
  filter(intensity>0) %>% 
  group_by(Mass, Forest,Class) %>% 
  dplyr::summarise(intensity = mean(intensity)) %>% 
  ungroup %>% 
  group_by(Forest,Class) %>% 
  dplyr::summarize(peaks = n()) %>% 
# get totals
  ungroup %>% 
  group_by(Forest) %>% 
  dplyr::mutate(total = sum(peaks))->
  counts
  

### 2.1.1 peaks common to both forests ----
master %>% 
  filter(Treatment=="PreFenton") %>% 
  filter(presence>0) %>% 
  group_by(Mass) %>% 
  dplyr::summarise(counts = n())->peakcounts_forest

total_peaks = nrow(peakcounts_forest)

peakcounts_forest %>% filter(counts==2) %>% dplyr::summarise(common_peaks = n())
common_peaks = nrow(peakcounts_forest$counts==2)

#
## 2.2 FENTON PEAK COUNTS  ----
# to determine peak counts in pre- vs. post-Fenton extracts
rawmaster %>% 
  filter(Goethite=="PreGoethite") %>% 
  filter(intensity>0) %>% 
  group_by(Mass,Fenton,Class) %>% 
  dplyr::summarise(intensity = mean(intensity)) %>% 
  ungroup %>% 
  group_by(Fenton,Class) %>% 
  dplyr::summarize(peaks = n()) %>% 
  ungroup %>% 
  group_by(Fenton) %>% 
  dplyr::mutate(total = sum(peaks))->
  fenton_counts
  
## 2.3 GOETHITE PEAK COUNTS  ----
rawmaster %>% 
  filter(intensity>0) %>% 
  group_by(Mass,Goethite,Fenton,Class) %>% 
  dplyr::summarise(intensity = mean(intensity)) %>% 
  ungroup %>% 
  group_by(Fenton,Goethite,Class) %>% 
  dplyr::summarize(peaks = n()) %>% 
  ungroup %>% 
  group_by(Fenton,Goethite) %>% 
  dplyr::mutate(total = sum(peaks))->
  goethite_counts


#
# ---------------------------------------------------------------------------- ----
# FENTON relative abundance lost vs. gained ----
# merge `fenton` file with `relative_intensity_percentile`

fenton_loss = 
  fenton %>% 
  left_join(hcoc, by = "Mass") %>% 
  left_join(select(elements, Mass, O), by = "Mass")

### OUTPUT
write.csv(fenton_loss, FENTON_LOSS)

#
# ---------------------------------------------------------------------------- ----
# GOETHITE adsorbed vs. non-adsorbed  ----
## sorbed vs. unbound ----

goethite %>% 
  mutate(fenton = factor(fenton, levels = c("PreFenton","PostFenton"))) %>% # order the levels
  dplyr::group_by(Forest,fenton) %>% 
# replace all NA with 0
  replace(., is.na(.),0) %>%
# create a binary `adsorbed` column for adsorbed/not adsorbed
  dplyr::mutate(adsorbed = case_when(PreGoethite>0 &PostGoethite==0~"sorbed",
                                     PreGoethite>0&PostGoethite>0~"unbound",
                                     PreGoethite==0&PostGoethite>0~"new")) %>% 
  dplyr::select(Mass,Forest,fenton,adsorbed) %>% 
  left_join(hcoc, by="Mass") %>% 
  left_join(class, by="Mass")->  goethite_sorption


### OUTPUT
write.csv(goethite_sorption, GOETHITE_ADSORPTION, row.names = FALSE)


#
# ---------------------------------------------------------------------------- ----
# GOETHITE relative abundance of sorbed vs. unsorbed groups ----

# first, subset the goethite_relabund file

goethite_sorption %>% 
  group_by(Forest,fenton,adsorbed,Class) %>% 
  dplyr::summarize(compounds = n()) %>% # this gives total intensity for each group
  # in the same command, we will also create a column for total intensity
  ungroup() %>% # remove the previous grouping
  group_by(Forest, fenton, adsorbed) %>% # create a new grouping 
  dplyr::mutate(total = sum(compounds)) %>%  # add a column for total intensity
  # we can also calculate the relative abundance in the same command
  mutate(relabund = (compounds/total)*100) %>% # calculate relative abundance as a %
  mutate(relabund = round(relabund,2))-> # round to two decimal places
  goethite_relabund


## use this in the graph for relative distribution

### OUTPUT
write.csv(goethite_relabund,GOETHITE_RELABUND, row.names = FALSE)

#
