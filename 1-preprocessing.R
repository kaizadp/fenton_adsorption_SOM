# 1-preprocessing

source("0-packages.R")

### all input files are in `stomfiles` folder
# we do not have separate files for smaple data vs. meta data, so first we need to create the separate files


## INPUT FILES ----
HW_PREFENTONGOETHITE = read.csv("stomfiles/PreFentonHWAdsorp-Master.csv") 
HW_POSTFENTONGOETHITE = read.csv("stomfiles/PostFenHWAdsorp-Master.csv") #needs cleaning
SW_PREFENTONGOETHITE = read.csv("stomfiles/PreFentonSWAdsorp.csv") #ok
SW_POSTFENTONGOETHITE = read.csv("stomfiles/PostFentonSWAdsorp.csv")

## INPUT -- META ----

# because different files have potentially different sets of peaks, we want to import all four files, get the relevant columns, combine, and then remove duplicates.
# this will ensure we have captured all the necessary peaks for the meta-data file
# tiring, yes

# remove the sample data columns. these are coded as xxx.csv
HW_PREFENTONGOETHITE %>% 
  dplyr::select(-ends_with(".csv")) -> 
  meta_HW_PREFENTONGOETHITE
  
HW_POSTFENTONGOETHITE %>% 
  select(-ends_with(".csv"))->
  meta_HW_POSTFENTONGOETHITE

SW_PREFENTONGOETHITE %>% 
  select(-ends_with(".csv"))->
  meta_SW_PREFENTONGOETHITE

SW_POSTFENTONGOETHITE %>% 
  select(-ends_with(".csv")) -> 
  meta_SW_POSTFENTONGOETHITE

# set all column names to be consistent across all
# manually checked before doing so -- the relevant columns have the same names and are in the same positions.
# only the irrelevant columns (AvgInitial vs. InitialAvg) are mislabelled, so it is ok

names(meta_HW_POSTFENTONGOETHITE) = names(meta_HW_PREFENTONGOETHITE)
names(meta_SW_POSTFENTONGOETHITE) = names(meta_HW_PREFENTONGOETHITE)
names(meta_SW_PREFENTONGOETHITE) = names(meta_HW_PREFENTONGOETHITE)

# confirm that the column names are identical across all sheets
identical(names(meta_HW_PREFENTONGOETHITE),names(meta_HW_POSTFENTONGOETHITE))
identical(names(meta_HW_POSTFENTONGOETHITE),names(meta_SW_PREFENTONGOETHITE))
identical(names(meta_SW_PREFENTONGOETHITE),names(meta_SW_POSTFENTONGOETHITE))

# combine all four files using rbind, which stacks the sheets one below the other
# we need to ensure all sheets have the same column names, hence the exercise above
meta_RAW = rbind(meta_HW_PREFENTONGOETHITE,meta_HW_POSTFENTONGOETHITE, meta_SW_PREFENTONGOETHITE,meta_SW_POSTFENTONGOETHITE)

# this is f-ed up
# there isn't a single column with Class groupings, but one column for each grouping
# NA is coded as ---- here,
# so we need to code it back to NA and then create a single column that has all the Class assignments

# there are six Class assignments
  # 1. CondAr
  # 2. Aromatic
  # 3. HighUnsatLignin
  # 4. Aliph_noN
  # 5. SatFACarb
  # 6. Aliph_N

meta_RAW %>% 
  na_if(.,"----") %>% 
  mutate(Class = case_when(!is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&is.na(UnsatAliph.N)&is.na(SatFatAcCarb)~"CondAr",
                           is.na(PolyCyArom)&(Aromatic=="Aromatic")&is.na(HighUnsatLign)&is.na(UnsatAliph.N)&is.na(SatFatAcCarb)~"Aromatic",
                           is.na(PolyCyArom)&is.na(Aromatic)&(HighUnsatLign=="HUnSatLig")&is.na(UnsatAliph.N)&is.na(SatFatAcCarb)~"HighUnsatLign",
                           is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&(UnsatAliph.N=="AlipatNoN")&is.na(SatFatAcCarb)~"UnSatAliph_noN",
                           is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&is.na(UnsatAliph.N)&(SatFatAcCarb=="SatFACarb")~"SatFatAcCarb",
                           is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&is.na(UnsatAliph.N)&is.na(SatFatAcCarb)&(UnSatAlip.N=="Alipat+N")~"Aliphat+N"))->
  meta_RAW

# now select only the relevant columns
# rename the columns as needed
# remove duplicates
meta_RAW %>%
  select(C,H,N,O,S,P,
         mass,`H.C`,`O.C`,
         DBE, CRAM, `AI.mod`,
         Class,
         KM, NKM, KMD) %>% 
  dplyr::rename(
                HC = `H.C`,
                OC = `O.C`,
                Mass = mass,
                AI_mod = `AI.mod`) %>% 
  mutate(Mass = round(Mass,4), # round Mass to 4 decimal places. do this for all files so it is easy to merge later 
         HC = round(HC,2),
         OC = round(OC,2))%>% 
  distinct()-> # this removes duplicates
  meta_RAW_distinct

### OUTPUT
write.csv(meta_RAW_distinct,"stomfiles/meta_RAW.csv")


#
## INPUT -- DATA ----
# these contain peaks seen in all 3 replicates. data have been pre-filtered.

## select specific columns within each file and then merge all
HW_PREFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  data_HW_PREFENTONGOETHITE

HW_POSTFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  data_HW_POSTFENTONGOETHITE

SW_PREFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  data_SW_PREFENTONGOETHITE

SW_POSTFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  data_SW_POSTFENTONGOETHITE

RAW_MERGED = merge(data_HW_PREFENTONGOETHITE, data_HW_POSTFENTONGOETHITE, by = "mass", all.x = T, all.y = T)
RAW_MERGED = merge(RAW_MERGED, data_SW_PREFENTONGOETHITE, by = "mass", all.x = T, all.y = T)
RAW_MERGED = merge(RAW_MERGED, data_SW_POSTFENTONGOETHITE, by = "mass", all.x = T, all.y = T)

RAW_MERGED %>% 
  dplyr::rename(Mass = mass) %>% 
  mutate(Mass = round(Mass,4)) %>% # round the mass to four decimals
  gather(code, intensity, ends_with(".csv")) %>% # each sample is a different column, combine/gather them
  na_if(.,"NA") %>% # replace "NA" lettering and 0 with NA
  na_if(.,0) %>% 
  na.omit()-> # remove all rows with NA values
  RAW_DATA

# merge this with the meta file
RAW_DATA2 = merge(RAW_DATA,meta_RAW_distinct,by = "Mass", all.y = T)

# now merge this with soil_key
SOIL_KEY = read.csv("data/soil_key.csv")

RAW_DATA_LONG = merge(SOIL_KEY,RAW_DATA2, by = "code")

write.csv(RAW_DATA_LONG,FTICR_RAWMASTER_LONG)



## compare meta files from raw vs. master formularity NOT DOING THIS NOW
    # meta_processed = read.csv(FTICR_META)
    # 
    # meta_processed2 = meta_processed %>% 
    #   select(Mass, El_comp, Class, HC, OC)
    # meta_RAW_distinct2 = meta_RAW_distinct %>% 
    #   select(Mass, HC, OC, Class)
    # 
    # meta_combined = merge(meta_processed2, meta_RAW_distinct2, by = "Mass")

#

## PROCESSING DATA FILES ----

