#   adsorptive fractionation of SOM
#   Kaizad Patel
---------------------- 
----------------------

source("code/0-packages.R")


# INPUT FILES
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
############################
############################
# STOP HERE
############################
############################

#
# ---------------------------------------------------------------------------- ---- 


#


# ---------------------------------------------------------------------------- ---- 


# 4 summary of groups ----



#
## 4.2 summary of peaks for fenton and goethite ----
### PEAKS, not INTENSITIES

# create a single file that has fenton anf goethite information
# merge the fenton and goethite files
data_fg_merged = merge(fticr_data_fenton,fticr_data_goethite, by = c("Mass","Forest"), all = T)
names(data_fg_merged)
# select only the relevant columns for fenton loss, goethite adsorbed, and goethite new
data_fg_merged %>% 
  select(Mass, Forest, loss, adsorbed, new) %>% 
  dplyr::rename(fenton_loss = loss) %>% 
  dplyr::rename(goethite_new = new) %>% 
  mutate(Mass = round(Mass,4))->
  data_fg_merged2

# merge with meta_class
data_fg_merged3 = merge(data_fg_merged2,fticr_meta_class, by = "Mass", all.x = T)

# get counts by group and treatment
data_fentoncounts = summarySE(data_fg_merged3, measurevar = "Mass", 
                            groupvars = c("Forest","fenton_loss","Class"), 
                            na.rm = TRUE)
data_fentoncounts %>% 
  select(1:4)->
  data_fentoncounts2
  
  
data_adsorbedcounts = summarySE(data_fg_merged3, measurevar = "Mass", 
                                groupvars = c("Forest","adsorbed","Class"), 
                                na.rm = TRUE)
data_adsorbedcounts %>% 
  select(1:4)->
  data_adsorbedcounts2

data_goethite_newcounts = summarySE(data_fg_merged3, measurevar = "Mass", 
                                groupvars = c("Forest","goethite_new","Class"), 
                                na.rm = TRUE)
data_goethite_newcounts %>% 
  select(1:4)->
  data_goethite_newcounts2

### OUTPUT
write_csv(data_fentoncounts2, path = "fticr/data_fentoncounts2.csv")
write_csv(data_adsorbedcounts2, path = "fticr/data_adsorbedcounts2.csv")
write_csv(data_goethite_newcounts2, path = "fticr/data_goethite_newcounts2.csv")


## create summary tables of the counts
data_fentoncounts_summarytable = dcast(data_fentoncounts2[!is.na(data_fentoncounts2$fenton_loss),],
                                                Forest+fenton_loss~Class,value.var = "N") 

data_adsorbedcounts_summarytable = dcast(data_adsorbedcounts2[!is.na(data_adsorbedcounts2$adsorbed),],
                                       Forest+adsorbed~Class,value.var = "N") 

data_goethite_newcounts_summarytable = dcast(data_goethite_newcounts2[!is.na(data_goethite_newcounts2$goethite_new),],
                                       Forest+goethite_new~Class,value.var = "N") 

write_csv(data_fentoncounts_summarytable, path = "fticr/data_fentoncounts_summarytable.csv")
write_csv(data_adsorbedcounts_summarytable, path = "fticr/data_adsorbedcounts_summarytable.csv")
write_csv(data_goethite_newcounts_summarytable, path = "fticr/data_goethite_newcounts_summarytable.csv")

#

## 4.3 summary of elements for initial ####

# merge data_raw_long with meta_class
fticr_data_2_el = merge(fticr_data_2, fticr_meta_elements, by = "Mass", all.x = T)

# remove unnecessary columns and gather the elements
fticr_data_2_el %>% 
  select(-PreFentonGoethite, -PostFentonGoethite) %>% 
  gather(fenton, intensity, PreFenton:PostFenton) %>% 
  gather(element, el_ratio, C:Na)->
  fticr_data_2_el2

# replace 0 counts with NA and then remove all NA
fticr_data_2_el2[fticr_data_2_el2==0]<-NA
fticr_data_2_el2 = fticr_data_2_el2[complete.cases(fticr_data_2_el2),]
# get element counts
fticr_data_el_count = summarySE(fticr_data_2_el2, measurevar = "el_ratio",
                                groupvars = c("Forest","fenton","element"),
                                na.rm = TRUE)
fticr_data_el_count$summary = paste(round(fticr_data_el_count$el_ratio,1), "\u00B1", round(fticr_data_el_count$se,1))

data_prefenton_el_summarytable = dcast(fticr_data_el_count,
                                    element~Forest+fenton, value.var = "summary")

#### OUTPUT
write_csv(data_prefenton_el_summarytable, path = "fticr/fticr_prefenton_el_counts.csv")

#
## 4.4 summary of elements for fenton and goethite ----

# merge data_fg_merged2 with meta_class
data_fg_merged_el = merge(data_fg_merged2,fticr_meta_elements, by = "Mass", all.x = T)

# gather the elements into a single column
data_fg_merged_el %>% 
  gather(element, el_ratio, C:Na)->
  data_fg_merged_el2

# get element counts for fenton
data_fenton_el = summarySE(data_fg_merged_el2, measurevar = "el_ratio", 
                              groupvars = c("Forest","fenton_loss","element"), 
                              na.rm = TRUE)

# remove NA
data_fenton_el = data_fenton_el[complete.cases(data_fenton_el),]
data_fenton_el$el_ratio = round(data_fenton_el$el_ratio,1)

data_fenton_el_summarytable = dcast(data_fenton_el,
                                    element~Forest+fenton_loss, value.var = "el_ratio")

#
# for goethite, use file fticr_data_goethite
# merge with fticr_meta_elements and then repeat the same process

fticr_data_goethite_el = merge(fticr_data_goethite, fticr_meta_elements, by = "Mass", all.x = T)

# gather
fticr_data_goethite_el %>% 
  gather(element, el_ratio, C:Na)->
  fticr_data_goethite_el2

data_goethite_adsorbed_el = summarySE(fticr_data_goethite_el2, measurevar = "el_ratio", 
                           groupvars = c("Forest","Fenton","adsorbed","element"), 
                           na.rm = TRUE)
# remove NA
data_goethite_adsorbed_el = data_goethite_adsorbed_el[complete.cases(data_goethite_adsorbed_el),]
data_goethite_adsorbed_el$el_ratio = round(data_goethite_adsorbed_el$el_ratio,1)

data_goethite_adsorbed_el_summarytable = dcast(data_goethite_adsorbed_el,
                                    adsorbed+element~Forest+Fenton, value.var = "el_ratio")  

#
# now repeat the goethite stuff for new molecules
data_goethite_new_el = summarySE(fticr_data_goethite_el2, measurevar = "el_ratio", 
                                      groupvars = c("Forest","Fenton","new","element"), 
                                      na.rm = TRUE)
# remove NA
data_goethite_new_el = data_goethite_new_el[complete.cases(data_goethite_new_el),]
data_goethite_new_el$el_ratio = round(data_goethite_new_el$el_ratio,1)

data_goethite_new_el_summarytable = dcast(data_goethite_new_el,
                                               element~Forest+Fenton, value.var = "el_ratio")  


### OUTPUT
write_csv(data_fenton_el_summarytable, path = "fticr/data_fenton_el_counts.csv")
write_csv(data_goethite_adsorbed_el_summarytable, path = "fticr/data_goethite_ads_el_counts.csv")
write_csv(data_goethite_new_el_summarytable, path = "fticr/data_goethite_new_el_counts.csv")



#
# ---------------------------------------------------------------------------- ---- 
# 4. COUNTS for Venn diagram ----
fticr_data_3 %>% 
  filter(intensity>0) %>% 
  group_by(Forest, treatment) %>% 
  dplyr::summarise(counts = length(Mass))->
  fticr_data_counts

fticr_data_counts2 = summarySE(fticr_data_3[!fticr_data_3$intensity==0,], 
                               measurevar = "Mass", groupvars = c("treatment"), na.rm = TRUE)
fticr_data_counts3 = summarySE(fticr_data_3[!fticr_data_3$intensity==0,], 
                               measurevar = "Mass", groupvars = c("Forest","treatment"), na.rm = TRUE)

# counts of HW vs. SW in native (initial) SOM 
fticr_data_3 %>% 
  filter(treatment=="PreFenton") %>% # select only the initial samples 
  dcast(Mass~Forest+treatment, value.var = "intensity") %>% # wide-form, so forests are in separate columns
  na_if(.,"0") %>% # replace NA with 0
  #code for unique vs. both forests
  dplyr::mutate(counts = case_when(HW_PreFenton>0 & SW_PreFenton>0 ~"both",
                                   HW_PreFenton>0 & is.na(SW_PreFenton) ~"HW",
                                   is.na(HW_PreFenton) & SW_PreFenton>0 ~"SW")) %>% 
  group_by(counts)%>% 
  dplyr::summarise(initial = length(Mass))-> # create a final dataframe with just the counts
  fticr_initial_counts
  
# counts of PreFenton vs. PostFenton. ignoring forest type
fticr_data_3 %>% 
  filter(treatment=="PreFenton"| treatment == "PostFenton") %>% # select only the Pre-Goethite samples 
  na_if(.,"0") %>%
  group_by(Mass,treatment) %>% 
  dplyr::summarise(intensity= mean(intensity)) %>% 
  dcast(Mass~treatment, value.var = "intensity") %>% # wide-form, so treatments are in separate columns
  na_if(.,"0") %>% # replace NA with 0
  #code for unique vs. both forests
  dplyr::mutate(counts = case_when(PreFenton>0 & PostFenton>0 ~"both",
                                   PreFenton>0 & is.na(PostFenton) ~"PostFenton",
                                   is.na(PreFenton) & PostFenton>0 ~"PreFenton")) %>% 
  group_by(counts)%>% 
  dplyr::summarise(Fenton = length(Mass))-> # create a final dataframe with just the counts
  fticr_fenton_counts  
  

# ---------------------------------------------------------------------------- ---- 
## ANOVA for fenton effects
# test if Fenton reaction altered the relative abundances of the different groups in the native soil.
# we don't need the postGoethite data for this
# use dataframe `fticr_data_relabundance`
# one single dplyr command will do everything

fticr_data_relabundance %>% 
  filter(Goethite=="PreGoethite") %>% # keep only PreGoethite data
  select(-total) %>% # the total column, which was used to calculate relative abundance, is useless here
  gather(group, relabund, AminoSugar:Tannin) %>% # convert wide into long-form
  group_by(Forest,group) %>% 
  # this is tricky
  # in one line, calculate the lm -> anova and then retrieve the p-value for that and round it to 4 decimal places
  dplyr::summarise(p_value = round(anova(lm(log(relabund)~Fenton))$`Pr(>F)`[1],4))->
  fticr_relabund_fenton_p_value

# effect of forest type on relative abundance
# same as above, but Fenton and Forest are swapped
fticr_data_relabundance %>% 
  filter(Goethite=="PreGoethite") %>% # keep only PreGoethite data
  select(-total) %>% # the total column, which was used to calculate relative abundance, is useless here
  gather(group, relabund, AminoSugar:Tannin) %>% # convert wide into long-form
  group_by(Fenton,group) %>% 
  # this is tricky
  # in one line, calculate the lm -> anova and then retrieve the p-value for that and round it to 4 decimal places
  dplyr::summarise(p_value = round(anova(lm(log(relabund)~Forest))$`Pr(>F)`[1],4))->
  fticr_relabund_fenton_p_value2

#  
# ---------------------------------------------------------------------------- ---- 
# ---------------------------------------------------------------------------- ---- 
# ---------------------------------------------------------------------------- ---- 

# 4. NOSC processing  ---------------------- ####
fticr_data_nosc = merge (fticr_data_3,fticr_meta_subset, by = "Mass", all = T)

# remove all rows without a value
# the absences are coded as 0, which will be included in the histogram. 
# so convert all zeroes to NA and then remove all NA

fticr_data_nosc %>%
  mutate_all(~replace(., . == 0, NA))->
  fticr_data_nosc

fticr_data_nosc = fticr_data_nosc[complete.cases(fticr_data_nosc),]

# now make a subset with only relevant columns

fticr_data_nosc %>% 
  select(Mass, Forest, treatment, intensity, NOSC, Class)->
  fticr_data_nosc2

# split the treatment column into fenton and goethite
# this makes faceting easier for the graphs

setDT(fticr_data_nosc2)[treatment == "PreFenton", fenton := "initial"]
fticr_data_nosc2[treatment == "PreFentonGoethite", fenton := "initial"]
fticr_data_nosc2[treatment == "PostFenton", fenton := "post-Fenton"]
fticr_data_nosc2[treatment == "PostFentonGoethite", fenton := "post-Fenton"]

setDT(fticr_data_nosc2)[treatment == "PreFenton", goethite := "pre-Goethite"]
fticr_data_nosc2[treatment == "PostFenton", goethite := "pre-Goethite"]
fticr_data_nosc2[treatment == "PreFentonGoethite", goethite := "post-Goethite"]
fticr_data_nosc2[treatment == "PostFentonGoethite", goethite := "post-Goethite"]

### OUTPUT
write_csv(fticr_data_nosc2, path = "fticr/fticr_data_nosc.csv")


#
#  ----------------------  ---------------------- #### ####
#  ----------------------  ---------------------- #### ####
### ### ### ### ### ###  ------
# begin calculating for summary
data_fg_merged2_long %>% 
  group_by(Forest,Forest,fenton,Class) %>% 
  dplyr::summarize(compounds = sum(intensity)) ->
  fticr_fg_soil_groups

# fticr_data_groups$compounds = as.numeric(fticr_data_groups$compounds)

fticr_data_soil_groups_wide = spread(fticr_data_soil_groups,Class,compounds)
fticr_data_soil_groups_wide = as.data.frame(fticr_data_soil_groups_wide)

# create a `total` column adding counts across all "group" columns (columns 4-10)
fticr_data_soil_groups_wide %>%
  mutate(total = rowSums(.[4:10])) ->
  fticr_data_soil_groups_wide

## relative abundance:
# split the dataset into (a) just the abundance values for easy calculations, and (b) the core key. Then combine again.
fticr_data_soil_groups_wide[,-c(1:3)] %>% 
  sapply('/', fticr_data_soil_groups_wide$total)->
  fticr_data_abundance

fticr_data_abundance = data.frame(fticr_data_abundance)
soilnames = data.frame(fticr_data_soil_groups_wide[,c(1:3)])

fticr_data_relabundance = cbind(soilnames,fticr_data_abundance)

### OUTPUT
write_csv(fticr_data_relabundance,path = "fticr_pore_relabund_soils.csv")

## relative abundance by treatment/site
# convert to long form and then do summary
fticr_data_relabundance_long = fticr_data_relabundance %>% 
  gather(Class, relabund, AminoSugar:total)


fticr_relabundance_summary = summarySE(fticr_data_relabundance_long, measurevar = "relabund", 
                                       groupvars = c("Forest","Treatment","Class"),na.rm = TRUE)
fticr_relabundance_summary$relativeabundance = paste((round(fticr_relabundance_summary$relabund,3)),
                                                     "\u00B1",
                                                     round(fticr_relabundance_summary$se,3))

fticr_relabundance_summary_summarytable = dcast(fticr_relabundance_summary,
                                                Forest+Treatment~Class,value.var = "relativeabundance") 

# move Unnamed and total columns to the end
# "Unnamed" is "Other" for pores
# fticr_pore_relabundance_summarytable %>% 
#  select(-Other,Other) %>% 
#  select(-total,total) ->
#  fticr_pore_relabundance_summarytable

# remove +/- SE values for the total column
fticr_relabundance_summary_summarytable$total="1"
## ## some cells have +/- 0. probably because n=1 for those. (??!!) double-check. 

### OUTPUT
# write.csv(fticr_soil_relabundance_summarytable,"fticr_soil_relabundance_groups.csv")
write_csv(fticr_relabundance_summary_summarytable,path = "fticr_relabundance_groups.csv")

#


