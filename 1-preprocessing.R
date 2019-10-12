# 1-preprocessing

source("0-packages.R")

### all input files are in `stomfiles` folder
# we do not have separate files for smaple data vs. meta data, so first we need to create the separate files


## INPUT FILES -- META ----
# because different files have potentially different sets of peaks, we want to import all four files, get the relevant columns, combine, and then remove duplicates.
# this will ensure we have captured all the necessary peaks for the meta-data file
# tiring, yes


meta_HW_PREFENTONGOETHITE = read.csv("stomfiles/PreFentonHWAdsorp-Master.csv") 
meta_HW_POSTFENTONGOETHITE = read.csv("stomfiles/PostFenHWAdsorp-Master.csv") #needs cleaning
meta_SW_PREFENTONGOETHITE = read.csv("stomfiles/PreFentonSWAdsorp.csv") #ok
meta_SW_POSTFENTONGOETHITE = read.csv("stomfiles/PostFentonSWAdsorp.csv")

# remove the sample data columns. these are coded as xxx.csv
meta_HW_PREFENTONGOETHITE %>% 
  dplyr::select(-ends_with(".csv")) -> 
  meta_HW_PREFENTONGOETHITE
  
meta_HW_POSTFENTONGOETHITE %>% 
  select(-ends_with(".csv"))->
  meta_HW_POSTFENTONGOETHITE

meta_SW_PREFENTONGOETHITE %>% 
  select(-ends_with(".csv"))->
  meta_SW_PREFENTONGOETHITE

meta_SW_POSTFENTONGOETHITE %>% 
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
         `m.z`,DBE, CRAM, `AI.mod`,
         Class,
         KM, NKM, KMD) %>% 
  dplyr::rename(
                HC = `H.C`,
                OC = `O.C`,
                Mass = mass,
                AI_mod = `AI.mod`) %>% 
  mutate(Mass = round(Mass,4))%>% # round Mass to 4 decimal places. do this for all files so it is easy to merge later 
  distinct()-> # this removes duplicates
  meta_RAW_distinct

### OUTPUT
write.csv(meta_RAW_distinct,"stomfiles/meta_RAW.csv")


#
## INPUT FILES -- DATA ----
# these contain peaks seen in all 3 replicates. data have been pre-filtered.


HW_FENTON = read_xlsx("stomfiles/FentonHWEffects.xlsx", sheet = "merged") #ok
SW_FENTON = read_xlsx("stomfiles/FentonSWEffects.xlsx", sheet = "merged") #ok

HW_PREFENTONGOETHITE = read_xlsx("stomfiles/PreFentonHWAdsorp-Master.xlsx", sheet = "merged") #ok
HW_POSTFENTONGOETHITE = read_xlsx("stomfiles/PostFenHWAdsorp-Master.xlsx", sheet = "RAW") #needs cleaning

SW_PREFENTONGOETHITE = read_xlsx("stomfiles/PreFentonSWAdsorp.xlsx", sheet = "merged") #ok
SW_POSTFENTONGOETHITE = read_xlsx("stomfiles/PostFentonSWAdsorp.xlsx", sheet = "merged")

# `HW_POSTFENTONGOETHITE` has extra crap at the end, so remove by subsetting
HW_POSTFENTONGOETHITE = HW_POSTFENTONGOETHITE[1:2231,]

## select specific columns within each file and then merge all
HW_FENTON %>% 
  select(mass, ends_with(".csv"))->
  HW_FENTON2
  
SW_FENTON %>% 
  select(mass, ends_with(".csv"))->
  SW_FENTON2

HW_PREFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  HW_PREFENTONGOETHITE2


HW_POSTFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  HW_POSTFENTONGOETHITE2


SW_PREFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  SW_PREFENTONGOETHITE2


SW_POSTFENTONGOETHITE %>% 
  select(mass, ends_with(".csv"))->
  SW_POSTFENTONGOETHITE2


RAW_MERGED = merge(HW_PREFENTONGOETHITE2, HW_POSTFENTONGOETHITE2, by = "mass", all.x = T, all.y = T)
RAW_MERGED = merge(RAW_MERGED, SW_PREFENTONGOETHITE2, by = "mass", all.x = T, all.y = T)
RAW_MERGED = merge(RAW_MERGED, SW_POSTFENTONGOETHITE2, by = "mass", all.x = T, all.y = T)


RAW_MERGED %>% 
  dplyr::rename(Mass = mass) %>% 
  mutate(Mass = round(Mass,4)) %>% # round the mass to four decimals
  
  dplyr::rename (
    Soil_1 = PreHW1.csv,
    Soil_2 = PreHW2.csv,
    Soil_3 = PreHW3.csv,
    Soil_4 = PreSW4.csv,
    Soil_5 = PreSW5.csv,
    Soil_6 = PreSW6.csv,
    Soil_7 = PosFHW1.csv,
    Soil_8 = PosFHW2.csv,
    Soil_9 = PosFHW3.csv,
    Soil_10 = PosSW10.csv,
    Soil_11 = PosSW11.csv,
    Soil_12 = PosSW12.csv,
    Soil_13 = PreHWGt13.csv,
    Soil_14 = PreHWGt14.csv,
    Soil_15 = PreHWGt15.csv,
    Soil_16 = PreSWGt16.csv,
    Soil_17 = PreSWGt17.csv,
    Soil_18 = PreSWGt18.csv,
    Soil_19 = `PosFHW-Gt1.csv`,
    Soil_20 = `PosFHW-Gt2.csv`,
    Soil_21 = `PosFHW-Gt3.csv`,
    Soil_22 = PosSW22.csv,
    Soil_23 = PosSW23.csv,
    Soil_24 = PosSW24.csv) %>% 
  na_if(.,"NA") %>% 
  na_if(.,0)  ->
  RAW_DATA

RAW_DATA2 = merge(meta_RAW_distinct,RAW_DATA,by = "Mass", all.y = T)

RAW_DATA2 %>% 
  gather(soil, intensity, starts_with("Soil_")) %>% 
  na.omit()->
  RAW_DATA_LONG

RAW_DATA_LONG = merge(RAW_DATA_LONG, soil_key, by = "soil")

    # RAW_DATA_LONG %>% 
    #   na_if(.,"----") %>% 
    #   mutate(Class = case_when(!is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&is.na(UnSatAliph_noN)&is.na(SatFatAcCarb)~"CondAr",
    #          is.na(PolyCyArom)&(Aromatic=="Aromatic")&is.na(HighUnsatLign)&is.na(UnSatAliph_noN)&is.na(SatFatAcCarb)~"Aromatic",
    #          is.na(PolyCyArom)&is.na(Aromatic)&(HighUnsatLign=="HUnSatLig")&is.na(UnSatAliph_noN)&is.na(SatFatAcCarb)~"HighUnsatLign",
    #          is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&(UnSatAliph_noN=="AlipatNoN")&is.na(SatFatAcCarb)~"UnSatAliph_noN",
    #          is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&is.na(UnSatAliph_noN)&(SatFatAcCarb=="SatFACarb")~"SatFatAcCarb",
    #          is.na(PolyCyArom)&is.na(Aromatic)&is.na(HighUnsatLign)&is.na(UnSatAliph_noN)&is.na(SatFatAcCarb)&(UnSatAlip.N=="Alipat+N")~"Aliphat+N"))->
    #   RAW_DATA_LONG
  
RAW_DATA_LONG %>% 
  select(-PolyCyArom,-Aromatic,-HighUnsatLign,-UnSatAliph_noN,-SatFatAcCarb,-UnSatAlip.N)->
  RAW_DATA_LONG


write.csv(RAW_DATA_LONG,FTICR_RAWMASTER_LONG)



## compare meta files from raw vs. master formularity ----

meta_processed = read.csv(FTICR_META)

meta_processed2 = meta_processed %>% 
  select(Mass, El_comp, Class, HC, OC)
meta_RAW_distinct2 = meta_RAW_distinct %>% 
  select(Mass, HC, OC, Class)

meta_combined = merge(meta_processed2, meta_RAW_distinct2, by = "Mass")


