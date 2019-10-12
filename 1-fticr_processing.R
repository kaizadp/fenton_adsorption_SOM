#   adsorptive fractionation of SOM
#   Kaizad Patel
---------------------- 
----------------------

source("0-packages.R")
# ---------------------------------------------------------------------------- ---- 

# 1. load files ---------------------- ####
## 1.1 create metadata file ----
# importing full data file
fticr_meta = read_excel("stomfiles/Master-Formularity-reprocessing V2.xlsx", sheet = "RAW")

# remove sample data columns. All sample names are Pre-xxx or Post-xxx
# remove unnamed columns 39-40-41
fticr_meta %>% 
  select(-starts_with("Pre")) %>% 
  select(-starts_with("Post")) %>% 
  select(-c("...39","...40","...41")) %>% 
  select(-c("C13", "Candidates","outlier")) %>% 
  mutate(Mass = round(Mass,4))-> # round this to 4 decimal places
  fticr_meta

#fticr_meta$Mass = round(fticr_meta$Mass,4)

# create a new column for NOSC
fticr_meta %>% 
  mutate(NOSC = 4 - (((4*C) + H - (3*N) - (2*O) - (2*S))/C)) %>% 
  dplyr::rename(HC = H.C) %>% 
  dplyr:: rename(OC = O.C) %>% 
  mutate(OC = round(HC,2),
         HC = round(OC,2))->
  fticr_meta

# create a smaller meta file for only HC OC
fticr_meta %>% 
  select(Mass, HC, OC) ->
  fticr_meta_hcoc

# create a subset of META for only relevant columns
fticr_meta %>% 
  select(Mass, HC, OC,
         El_comp,
         NOSC, 
         Class,
         KM_CH2, KMD_CH2,
         DBE, AI_0,
         C:Na)->
  fticr_meta_subset

# create a subset with only Mass and Class
fticr_meta %>% 
  select(Mass, Class)->
  fticr_meta_class

# create a subset with only Mass and elements
fticr_meta %>% 
  select(Mass:Na)->
  fticr_meta_elements

### OUTPUT
# SAVE META FILES
write_csv(fticr_meta, FTICR_META)
write_csv(fticr_meta_hcoc, HCOC)
write_csv(fticr_meta_subset, path = "fticr/fticr_meta_subset.csv")
write_csv(fticr_meta_elements, path = "fticr/fticr_meta_elements.csv")


## 1.2 import data file ----
### 1.2.1 final merged datafile ----

fticr_data = read_excel("Master-Formularity-reprocessing V2.xlsx", sheet = "RAW")

# subset only "Mass" and sample columns
fticr_data %>% 
  select("Mass",starts_with("Pre"), starts_with("Post")) %>% 
  dplyr::rename(PostFentonGoethiteSW = `Post FentonGoethiteSW`) %>% 
  dplyr::rename(PostFentonSW = `Post FentonSW`) %>% 
  mutate(Mass = round(Mass,4))->
  fticr_data

# collapse columns by forest type.
# first, separate the data file into columns for hw or sw

fticr_data %>% 
  select(Mass,ends_with("HW"))->
  fticr_data_hw

fticr_data %>% 
  select(Mass,ends_with("SW"))->
  fticr_data_sw  

# then, rename to remove the "HW"/"SW" qualifier in the column name. Instead, replace it with a new column "Forest"

fticr_data_hw %>% 
  mutate(Forest = "HW") %>% 
  dplyr::rename(PreFenton = PreFentonHW) %>% 
  dplyr::rename(PreFentonGoethite = PreFentonGoethiteHW) %>% 
  dplyr::rename(PostFenton = PostFentonHW) %>% 
  dplyr::rename(PostFentonGoethite = PostFentonGoethiteHW) ->
  fticr_data_hw
  
fticr_data_sw %>% 
  mutate(Forest = "SW") %>% 
  dplyr::rename(PreFenton = PreFentonSW) %>% 
  dplyr::rename(PreFentonGoethite = PreFentonGoethiteSW) %>% 
  dplyr::rename(PostFenton = PostFentonSW) %>% 
  dplyr::rename(PostFentonGoethite = PostFentonGoethiteSW) ->
  fticr_data_sw  

# then, combine both files  
## This is now the MASTER DATA file
fticr_data_2 = rbind(fticr_data_hw,fticr_data_sw)

## OUTPUT
# write_csv(fticr_data_2, path = "fticr/fticr_data_master.csv")

# gather, to make a longform file

fticr_data_3 = fticr_data_2 %>% 
  gather(treatment, intensity, PreFenton:PostFentonGoethite)

## OUTPUT
write_csv(fticr_data_3, FTICR_MASTER_LONG)


### 1.2.2 raw files ----
soil_key = read_csv("data/soil_key.csv")
soil_key$Fenton = factor(soil_key$Fenton, levels = c("PreFenton","PostFenton"))
soil_key$Goethite = factor(soil_key$Goethite, levels = c("PreGoethite","PostGoethite"))
soil_key$Treatment = factor(soil_key$Treatment, levels = c("PreFenton","PostFenton","PreFentonGoethite","PostFentonGoethite"))

# data files for each forest/treatment
hw_pref = read_csv("data/hw_prefenton.csv")
hw_postf = read_csv("data/hw_postfenton.csv")
sw_pref = read_csv("data/sw_prefenton.csv")
sw_postf = read_csv("data/sw_postfenton.csv")

# combine all raw data files
fticr_data_raw = merge(hw_pref, hw_postf, by = "mass")
fticr_data_raw = merge(fticr_data_raw, sw_pref, by = "mass")
fticr_data_raw = merge(fticr_data_raw, sw_postf, by = "mass")

names(fticr_data_raw)

# rename all the f-ing columns to be consistent with the KEY
fticr_data_raw %>% 
  dplyr::rename(Mass = mass) %>% 
  mutate(Mass = round(Mass,4)) %>% # round the mass to four decimals
  
  dplyr::rename(Soil_1 = PreHW1.csv) %>% 
  dplyr::rename(Soil_2 = PreHW2.csv) %>%
  dplyr::rename(Soil_3 = PreHW3.csv) %>%
  
  dplyr::rename(Soil_4 = PreSW4.csv) %>%
  dplyr::rename(Soil_5 = PreSW5.csv) %>%
  dplyr::rename(Soil_6 = PreSW6.csv) %>%
  
  dplyr::rename(Soil_7 = PosFHW1.csv) %>%
  dplyr::rename(Soil_8 = PosFHW2.csv) %>%
  dplyr::rename(Soil_9 = PosFHW3.csv) %>%
  
  dplyr::rename(Soil_10 = PosSW10.csv) %>%
  dplyr::rename(Soil_11 = PosSW11.csv) %>%
  dplyr::rename(Soil_12 = PosSW12.csv) %>%
  
  dplyr::rename(Soil_13 = PreHWGt13.csv) %>%
  dplyr::rename(Soil_14 = PreHWGt14.csv) %>%
  dplyr::rename(Soil_15 = PreHWGt15.csv) %>%
  
  dplyr::rename(Soil_16 = PreSWGt16.csv) %>%
  dplyr::rename(Soil_17 = PreSWGt17.csv) %>%
  dplyr::rename(Soil_18 = PreSWGt18.csv) %>%
  
  dplyr::rename(Soil_19 = `PosFHW-Gt1.csv`) %>%
  dplyr::rename(Soil_20 = `PosFHW-Gt2.csv`) %>%
  dplyr::rename(Soil_21 = `PosFHW-Gt3.csv`) %>%
  
  dplyr::rename(Soil_22 = PosSW22.csv) %>%
  dplyr::rename(Soil_23 = PosSW23.csv) %>%
  dplyr::rename(Soil_24 = PosSW24.csv) ->
  fticr_data_raw

# merge with meta_class
fticr_data_raw2 = merge(fticr_meta_class,fticr_data_raw,by = "Mass", all.y = T)
## some rows have unassigned groups, coded as NA. WHY???
## remove NA

fticr_data_raw2 = fticr_data_raw2[complete.cases(fticr_data_raw2),]

# write_csv(fticr_data_raw, path = "fticr/fticr_data_raw.csv")


# convert to long-form and then merge with the soil KEY
fticr_data_raw2 %>% 
  gather(soil, intensity, Soil_1:Soil_24)->
  fticr_data_raw_long

fticr_data_raw_long = merge(fticr_data_raw_long,soil_key, by = "soil", all.x = T)

### OUTPUT 
write_csv(fticr_data_raw_long, FTICR_RAWMASTER_LONG)

#

# ---------------------------------------------------------------------------- ---- 

# 2 process data ---------------------- ####
## create a new file for Fenton ----

# select only the relevant columns. don't include the Goethite columns
fticr_data_2 %>% 
  select(Mass, Forest, PreFenton, PostFenton)->
  fticr_data_fenton

#

## create a new file for goethite ----
# split the MASTER file into two files, for pre-Goethite vs. post-Goethite. gather and add columns indicating whether pre or post F. and then combine. 
fticr_data_2 %>% 
  select(Mass, Forest, PreFenton,PostFenton)->
  fticr_data_preg

fticr_data_2 %>% 
  select(Mass, Forest, PreFentonGoethite,PostFentonGoethite)->
  fticr_data_postg

fticr_data_preg %>% 
  gather(Fenton, PreGoethite, PreFenton:PostFenton)->
  fticr_data_preg2

fticr_data_postg %>% 
  dplyr::rename(PreFenton = PreFentonGoethite) %>% 
  dplyr::rename(PostFenton = PostFentonGoethite) %>% 
  gather(Fenton, PostGoethite, PreFenton:PostFenton)->
  fticr_data_postg2

#fticr_data_goethite2 = cbind(fticr_data_preg2,fticr_data_postg2)
fticr_data_goethite = merge(fticr_data_preg2,fticr_data_postg2)
# `cbind` keeps duplicate columns, `merge` deletes duplicate columns

### OUTPUT
write.csv(fticr_data_goethite,FTICR_GOETHITE,na="")
#
# ---------------------------------------------------------------------------- ---- 

## 2.1 PROCESSING FENTON: determining molecules lost and gained ----
# dataframe name is fticr_data_fenton

# create new column for molecules gained/lost during the Fenton reaction
setDT(fticr_data_fenton)[PreFenton >0 & PostFenton == 0, loss := "lost"]
fticr_data_fenton[PreFenton == 0 & PostFenton > 0, loss := "gained"]
fticr_data_fenton[PreFenton > 0 & PostFenton > 0, loss := "conserved"]
fticr_data_fenton$loss = ordered(fticr_data_fenton$loss, levels = c("lost", "gained", "conserved"))

### OUTPUT
write_csv(fticr_data_fenton, FTICR_FENTON, na = "")

# ---------------------------------------------------------------------- -

# STOP SCRIPT HERE
# EVERYTHING ELSE GOES TO THE NEXT SCRIPT FOR ABUNDANCE CALCULATIONS 

# ---------------------------------------------------------------------- -
# ---------------------------------------------------------------------- -
# ---------------------------------------------------------------------- -



### 2.1.1 fenton relative abundance lost vs. gained MOVE TO OTHER SCRIPT? ----
      # subset the fenton2 to get just the lost/gained column
      
      fticr_data_fenton %>% 
        select(Mass, Forest, loss)->
        fticr_data_fenton_loss
      
      # merge loss file with data_4_hcoc
      
      fticr_fenton_loss_relabund = merge(fticr_data_4_hcoc, fticr_data_fenton_loss, by = c("Mass", "Forest"))


## 2.2 PROCESSING GOETHITE: determining adsorbed vs. not adsorbed molecules ----
### 2.2.1. adsorbed/not adsorbed (binary classification) ----

# Using S/N method of Avneri-Katz 2017
# find minimum intensity
# divide all by minimum. if >2, SN = 10
# but first, convert all zero to NA

fticr_data_goethite[fticr_data_goethite==0]<-NA
minimum = min(c(fticr_data_goethite$PreGoethite, fticr_data_goethite$PostGoethite), na.rm = TRUE)

# then convert NA back to 0
fticr_data_goethite[is.na(fticr_data_goethite)]<-0


## don't need this any more
      # fticr_data_goethite %>% 
      #  mutate(PreG_sn = PreGoethite/minimum) %>% 
      #  mutate(PostG_sn = PostGoethite/minimum)->
      #  fticr_data_goethite

# create a column for adsorbed/not adsorbed
setDT(fticr_data_goethite)[PreGoethite/minimum >2 & PostGoethite/minimum < 1, adsorbed := "adsorbed"]
fticr_data_goethite[PreGoethite/minimum >2 & PostGoethite/minimum > 1, adsorbed := "not adsorbed"]

# create a column for new molecules created post-adsorption
setDT(fticr_data_goethite)[PreGoethite/minimum == 0 & PostGoethite/minimum > 1, new := "new molecules"]

### OUTPUT
write_csv(fticr_data_goethite,path = "fticr/fticr_data_goethite.csv")

#

### 2.2.2. relative strength of sorption ----
# technique from Williams, Borch et al. 2018. Soil Systems
# Calculate relative abundance of each formula in the PreG and PostG samples. Substract PostG-PreG to calculate delta-abundance.
# Use delta-abundance to group the molecules into seven classes:
# < -0.00015 = most sorbed | -0.00010 = more sorbed | -0.00005 = sorbed | 0.00005 = minimal change | 0.00010 = unbound | 0.00015 = more unbound | > 0.00015 = most unbound,
# ALL the calculations are done in one step.

fticr_data_goethite %>% 
  mutate(Fenton = factor(Fenton, levels = c("PreFenton","PostFenton"))) %>% # order the levels
  dplyr::group_by(Forest,Fenton) %>% 
  dplyr::mutate(preg_total = sum(PreGoethite)) %>% 
  dplyr::mutate(postg_total = sum(PostGoethite)) %>% 
  mutate(preg_rel_abund = PreGoethite/preg_total) %>% 
  mutate(postg_rel_abund = PostGoethite/postg_total) %>%
  mutate(delta_abund = postg_rel_abund - preg_rel_abund) %>% 
  mutate(sorption_frac = cut(delta_abund, 
                             breaks = c(-Inf,-0.00015, -0.00010, -0.00005, 0.00005,0.00010,0.00015,Inf),
                             labels = c("most sorbed", "more sorbed", "sorbed", "minimal change","unbound","more unbound","most unbound")))->
  fticr_data_goethite_relabund

fticr_data_goethite_relabund = merge(fticr_data_goethite_relabund,fticr_meta_hcoc, by = "Mass", all.x = T)  

### OUTPUT
write_csv(fticr_data_goethite_relabund, path = "fticr/fticr_data_goethite_relabund.csv")

