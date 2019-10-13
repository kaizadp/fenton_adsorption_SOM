#   adsorptive fractionation of SOM
#   Kaizad Patel
---------------------- 
----------------------


# DO NOT SOURCE SCRIPT #1.
# Run Script #1 separately, then start a new session (ctrl+shift+F10) and run this script.

source("0-packages.R")


# INPUT FILES
# use this file for meta for now. this may move to a different folder/name later
meta = read.csv("stomfiles/meta_RAW.csv")
hcoc = meta %>% select(Mass,HC,OC)
#meta = read.csv(FTICR_META)# <- "fticr/fticr_meta.csv" # all metadata about formula, etc. assignment for each m/z value
#hcoc = read.csv(HCOC)
master = read.csv(FTICR_MASTER_LONG)# <- "fticr/fticr_master_long.csv" #
rawmaster = read.csv(FTICR_RAWMASTER_LONG)# <- "fticr/fticr_rawmaster_long.csv"
fenton = read.csv(FTICR_FENTON)# <- "fticr/fticr_fenton.csv" # pre- and post-Fenton data, intensities only
goethite = read.csv(FTICR_GOETHITE)# <- "fticr/fticr_goethite.csv" # pre- and post-Goethite data, intensities only

# ---------------------------------------------------------------------------- ----
# this was just to compare raw counts vs. processed counts in the final. 
# not relevant any more because the processed came from the raw.

rawmaster %>% 
  filter(Treatment=="PreFenton") %>% 
  filter(intensity>0) %>% 
  group_by(Mass, Forest) %>% 
  dplyr::summarise(intensity = mean(intensity))  %>% 
  spread(Forest, intensity)->
  raw
write.csv(raw, "raw.csv",na="")

master %>% 
  filter(treatment=="PreFenton") %>% 
  filter(intensity>0) %>% 
  group_by(Mass, Forest) %>% 
  dplyr::summarise(intensity = mean(intensity)) %>% 
  spread(Forest, intensity)->
  processed
write.csv(processed, "master.csv", na="")
# ---------------------------------------------------------------------------- ----

# 1. relative intensity of each formula. and percentile ----
# this portion of the script will assign the molecules into quartiles based on relative abundance.
# this classification will be used in Van Krevelen diagrams

# `rawmaster` is the longform master file. calculate relative abundance of each molecule

master %>% 
  mutate(intensity = as.numeric(intensity)) %>% # set intensity as a numeric variable
  filter(Treatment == "PreFenton"| Treatment=="PostFenton") %>% # keep only pre-Goethite data. we don't want post-adsorption data
  na.omit %>% # remove all NA, or it won't calculate
  dplyr::group_by(Forest, Treatment) %>% 
  dplyr::mutate(total = sum(intensity)) %>% # add a new column calculating total intensity
  dplyr::mutate(rel_abund = (intensity/total)*100)-> # calculate relative intensity as a %
  relative_intensity

# then create a column for quartiles
relative_intensity %>% 
  dplyr::group_by(Forest, Treatment) %>% 
  dplyr::mutate(percentile = ntile(rel_abund, 100)) %>% 
  mutate(perc = cut(percentile, 
                    breaks = c(-Inf,25, 50, 75, Inf),
                    labels = c("lowest 25%", "third 25 %", "second 25 %", "top 25 %")))->
  relative_intensity_percentile

# remove unnecessary columns
relative_intensity_percentile %>% 
  select(-intensity, -total, -percentile)->
  relative_intensity_percentile
  
# merge with the hcoc file
relative_intensity_percentile = merge(relative_intensity_percentile,hcoc, by = "Mass", all.x=T)
# relative_intensity_percentile = merge(relative_intensity_percentile,hcoc, by = "Mass", all.x = T)

ggplot(relative_intensity_percentile, aes(x = OC,y = HC, color = perc))+
  geom_point(alpha = 0.5)+
  scale_color_brewer(palette = "Reds")+
  facet_grid(Forest~Treatment)

### OUTPUT
write_csv(relative_intensity_percentile, PERCENTILE)
# this file will be used for the Van Krevelen plots (preFenton and postFenton)

#
# ---------------------------------------------------------------------------- ----

# 2. RELATIVE ABUNDANCE FOR PRE- AND POST-FENTON GROUPS ----
# `rawmaster` is the longform master file. 

# summarizing by groups
rawmaster %>% 
  mutate(intensity = as.numeric(intensity)) %>% # set intensity as a numeric variable
  filter(Goethite == "PreGoethite") %>% # keep only pre-Goethite data. we don't want post-adsorption data
  group_by(Forest,Fenton,soil,Class) %>% 
  dplyr::summarize(compounds = sum(as.numeric(intensity), na.rm = TRUE)) %>% # this gives total intensity for each group
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
  dplyr::summarise(rel_abund = mean(relabund), 
                   se = sd(relabund)/sqrt(n())) %>% 
  dplyr::mutate(rel_abund = round(rel_abund,2),
                se = round(se,2),
                relativeabund = paste(rel_abund,"\u00B1",se))->
  raw_groups

# now do Tukey HSD

fit_hsd <- function(dat) {
  a <-aov(relabund ~ Fenton, data = dat)
  h <-HSD.test(a,"Fenton")
  #create a tibble with one column for each treatment
  #the hsd results are row1 = drought, row2 = saturation, row3 = time zero saturation, row4 = field moist. hsd letters are in column 2
  tibble(`PreFenton` = h$groups["PreFenton",2], 
         `PostFenton` = h$groups["PostFenton",2])
}

raw_coregroups %>% 
  group_by(Forest, Class) %>% 
  do(fit_hsd(.)) %>% 
# ^ the script above creates a data.frame with columns `Forest`, `Class`, `PreFenton`,`PostFenton`
# in the same command, we are gathering the PreFenton and PostFenton columns into a single column, `Fenton`. hashtag efficiency
  gather(Fenton, hsd, 3:4)->
  hsd

# now combine `raw_groups` with `hsd`
raw_groups_hsd = merge(raw_groups,hsd, by = c("Forest", "Fenton","Class"))

# now combine the `relativeabund` and `hsd` columns and then remove all unnecessary columns
raw_groups_hsd %>% 
  mutate(relabund_hsd = paste(relativeabund, hsd)) %>% 
  select(-se, -relativeabund,-hsd)->
  raw_groups_hsd

### OUTPUT
write.csv(raw_groups_hsd,RELATIVE_ABUND)
# write_csv(fticr_relabundance_summary_summarytable,path = "output/table1_relabundance_groups_bytrt.csv")


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

fenton_loss = merge(fenton, relative_intensity_percentile, by = c("Mass", "Forest"))

ggplot(fenton_loss[!fenton_loss$loss=="conserved",], aes(x = OC,y = HC, color = loss))+
  geom_point(alpha = 0.5)+
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(~Forest)

# ---------------------------------------------------------------------------- ----
# GOETHITE adsorbed vs. non-adsorbed  ----
## .1 determining adsorbed vs. not adsorbed molecules ----
      ## NOT DOING THIS FOR NOW
      ##    # (binary classification) 
      ##    # Using S/N method of Avneri-Katz 2017
      ##    # find minimum intensity
      ##    # divide all by minimum. if >2, SN = 10
      ##    # but first, convert all zero to NA
      ##    
      ##    fticr_data_goethite[fticr_data_goethite==0]<-NA
      ##    minimum = min(c(fticr_data_goethite$PreGoethite, fticr_data_goethite$PostGoethite), na.rm = TRUE)
      ##    
      ##    # then convert NA back to 0
      ##    fticr_data_goethite[is.na(fticr_data_goethite)]<-0
      ##    
      ##    # create a column for adsorbed/not adsorbed
      ##    setDT(fticr_data_goethite)[PreGoethite/minimum >2 & PostGoethite/minimum < 1, adsorbed := "adsorbed"]
      ##    fticr_data_goethite[PreGoethite/minimum >2 & PostGoethite/minimum > 1, adsorbed := "not adsorbed"]
      ##    
      ##    # create a column for new molecules created post-adsorption
      ##    setDT(fticr_data_goethite)[PreGoethite/minimum == 0 & PostGoethite/minimum > 1, new := "new molecules"]
      ##    
      ##    ### OUTPUT
      ##    write_csv(fticr_data_goethite,path = "fticr/fticr_data_goethite.csv")
      ## 
#

## .2 relative strength of sorption ----

# technique from Williams, Borch et al. 2018. Soil Systems
# Calculate relative abundance of each formula in the PreG and PostG samples. 
# Subtract PostG-PreG to calculate delta-abundance.
# Use delta-abundance to group the molecules into seven classes:
# < -0.00015 = most sorbed | -0.00010 = more sorbed | -0.00005 = sorbed | 0.00005 = minimal change | 0.00010 = unbound | 0.00015 = more unbound | > 0.00015 = most unbound,
# ALL the calculations are done in one step.

data_goethite %>% 
  mutate(fenton = factor(fenton, levels = c("PreFenton","PostFenton"))) %>% # order the levels
  dplyr::group_by(Forest,fenton) %>% 
# replace all NA with 0
  replace(., is.na(.),0) %>% 
# new columns for pre-goethite and post-goethite total intensities  
  dplyr::mutate(preg_total = sum(PreGoethite, na.rm = TRUE),
                postg_total = sum(PostGoethite, na.rm = TRUE)) %>% 
# new columns for relative abundance as fraction
  mutate(preg_rel_abund = PreGoethite/preg_total) %>% 
  mutate(postg_rel_abund = PostGoethite/postg_total) %>%
# subtract post-pre relative abundance
  mutate(delta_abund = postg_rel_abund - preg_rel_abund) %>% 
# create a column for binning 
  mutate(sorption_frac = cut(delta_abund, 
                             breaks = c(-Inf,-0.00015, -0.00010, -0.00005, 0.00005,0.00010,0.00015,Inf),
                             labels = c("most sorbed", "more sorbed", "sorbed", "minimal change","unbound","more unbound","most unbound"))) ->
# cleaning up: remove unnecessary columnns  
#  select(-(preg_total:delta_abund))->
  data_goethite_relabund

data_goethite_relabund = merge(data_goethite_relabund,hcoc, by = "Mass", all.x = T)  

### OUTPUT
write_csv(data_goethite_relabund, GOETHITE_ADSORPTION)

ggplot(data_goethite_relabund, aes(x = OC,y = HC, color = sorption_frac))+
  geom_point(alpha = 0.2)+
  scale_color_brewer(palette = "PuOr")+
  facet_wrap(~fenton)+
  theme(legend.position = "top")

#
# ---------------------------------------------------------------------------- ----
# GOETHITE relative abundance of sorbed vs. unsorbed groups ----

# first, subset the goethite_relabund file

data_goethite_relabund %>% 
  select(Mass, Forest, fenton, PreGoethite, sorption_frac)->
  data_goethite_adsorbed

# the adsorbed_frac column has multiple levels. choose only the "most sorbed" and "most unbound"
data_goethite_adsorbed %>% 
  filter(sorption_frac=="most sorbed"| sorption_frac=="most unbound")->
  data_goethite_adsorbed  

# merge with the class meta file
data_goethite_adsorbed = merge(data_goethite_adsorbed, meta_CLASS, by = "Mass", all.x = T)

    ## # remove the "Other" class
    ## data_goethite_adsorbed %>% 
    ##   filter(!Class=="Other")->
    ##   data_goethite_adsorbed

    # fticr_data_goethite_relabund_adsorbed = fticr_data_goethite_relabund_adsorbed[
    #  !fticr_data_goethite_relabund_adsorbed$Class=="Other",]


# now follow steps for relative abundance of groups
# 

data_goethite_adsorbed %>% 
  group_by(Forest,fenton,sorption_frac,Class) %>% 
  dplyr::summarize(compounds = sum(as.numeric(PreGoethite), na.rm = TRUE)) %>% # this gives total intensity for each group
  # in the same command, we will also create a column for total intensity
  ungroup() %>% # remove the previous grouping
  group_by(Forest, fenton, sorption_frac) %>% # create a new grouping 
  dplyr::mutate(total = sum(compounds)) %>%  # add a column for total intensity
  # we can also calculate the relative abundance in the same command
  mutate(relabund = (compounds/total)*100) %>% # calculate relative abundance as a %
  mutate(relabund = round(relabund,2))-> # round to two decimal places
  data_goethite_adsorbed_relabund


ggplot(data_goethite_adsorbed_relabund, aes(x = Class, y = relabund, fill = sorption_frac))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_grid(fenton~Forest)
## use this in the graph for relative distribution


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


### UNNECESSARY ---------------------
### UNNECESSARY SCRIPT. DELETE LATER. 
###

hw_fenton = read_excel("Master-Formularity-reprocessing V2.xlsx", 
                        sheet = "HWPrePost") %>%
                         select (H.C,
                                O.C,
                                PreFentonHW,
                                PostFentonHW,
                                Group,
                                O)

sw_fenton = read_excel("Master-Formularity-reprocessing V2.xlsx", 
                      sheet = "SWPrePost")%>%
                      select (H.C,
                             O.C,
                             PreFentonSW,
                             PostFentonSW)

hw_pref_ads = read_excel("Master-Formularity-reprocessing V2.xlsx", 
                      sheet = "HW-PreF-Ads") %>%
                      select (H.C,
                             O.C,
                             Adsorbed,
                             New)

hw_postf_ads = read_excel("Master-Formularity-reprocessing V2.xlsx", 
                        sheet = "HW-PostF-Ads") %>%
                          select (H.C,
                                 O.C,
                                 Adsorbed,
                                 New)

sw_pref_ads = read_excel("Master-Formularity-reprocessing V2.xlsx", 
                        sheet = "SW-PreF-Ads") %>%
                          select (H.C,
                                 O.C,
                                 Adsorbed,
                                 New)

sw_postf_ads = read_excel("Master-Formularity-reprocessing V2.xlsx", 
                        sheet = "SW-PostF-Ads") %>%
                          select (H.C,
                                 O.C,
                                 Adsorbed,
                                 New)

## create new columns

hw_fenton$OC = hw_fenton$O.C
hw_fenton$HC = hw_fenton$H.C

hw_fenton$O = factor(hw_fenton$O, levels=c("0","1","2","3","4","5","6","7","8","9","10",
                                           "11","12","13","14","15","16","17","18","19","20",
                                           "21","22","23","24","25"))

names(hw_fenton)

# creating loss/gain column 
setDT(hw_fenton)[PreFentonHW >0 & PostFentonHW == 0, loss := "lost"]
setDT(hw_fenton)[PreFentonHW == 0 & PostFentonHW > 0, loss := "gained"]
setDT(hw_fenton)[PreFentonHW > 0 & PostFentonHW > 0, loss := "conserved"]
hw_fenton$loss = ordered(hw_fenton$loss, levels = c("lost", "gained", "conserved"))

setDT(sw_fenton)[PreFentonSW >0 & PostFentonSW == 0, loss := "lost"]
setDT(sw_fenton)[PreFentonSW == 0 & PostFentonSW > 0, loss := "gained"]
setDT(sw_fenton)[PreFentonSW > 0 & PostFentonSW > 0, loss := "conserved"]
sw_fenton$loss = ordered(sw_fenton$loss, levels = c("lost", "gained", "conserved"))


hw_lost=summarySE(hw_fenton[loss=="lost",],measurevar = "PreFentonHW", groupvars=c("loss","Group"),na.rm=TRUE);hw_lost
hw_gained=summarySE(hw_fenton[loss=="gained",],measurevar = "PostFentonHW", groupvars=c("Group"),na.rm=TRUE);hw_gained

attach(hw_fenton)
attach(sw_fenton)
#
#### -------- convex hull -- ??? ----------
library(alphahull)
alpha = 100
alphashape = ashape(hw_fenton[which(hw_fenton$conserved=="conserved")],c("O.C","H.C"),alpha = alpha)

find_hull <- function(hw_fenton) hw_fenton[chull(hw_fenton$`O.C`, hw_fenton$`H.C`), ]
hulls <- ddply(hw_fenton, "conserved",find_hull)

chulls <- ddply(hw_fenton, .(conserved), function(hw_fenton) hw_fenton[chull(hw_fenton$OC, hw_fenton$HC), ])

cHullPoints = chull(hw_fenton)
polygon(hw_fenton[cHullPoints,])

hull = hw_fenton %>%
  slice(chull(OC, HC))
##
## 1. Fenton van krevelen ---------------------- ## -----
vankrev_hw_fentonloss=ggplot (hw_fenton,
                             aes (x=`O.C`,y=`H.C`,
                                 color=loss, 
                                 shape=loss,na.rm=TRUE))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,4,1),limits=c("lost","gained","conserved"))+
  scale_color_manual(values=c("darkred","darkgreen","grey"),limits=c("lost","gained","conserved"))+
  
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  #annotate("text", label = "HW", x = 0.1, y = 2.4)+ 
  ggtitle("HW fenton")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hw_fentonloss


vankrev_sw_fentonloss=ggplot(sw_fenton,aes(x=`O.C`,y=`H.C`,color=loss, shape=loss,na.rm=TRUE))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,4,1),limits=c("lost","gained","conserved"))+
  scale_color_manual(values=c("darkred","darkgreen","grey"),limits=c("lost","gained","conserved"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  #annotate("text", label = "SW fenton", x = 0.1, y = 2.4)+ 
  ggtitle("SW fenton")+

  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_sw_fentonloss

#
## intensity vs. O ---------------------- ## -----
hw_pref_o=summarySE(hw_fenton,measurevar = "PreFentonHW", groupvars=c("O"),na.rm=TRUE);hw_pref_o

ggplot(hw_pref_o,
       aes(x = O, y = PreFentonHW), geom_bar())


## 2. adsorbed vs. not adsorbed ----
names(hw_pref_ads)
hw_pref_ads$OC = hw_pref_ads$O.C
hw_pref_ads$HC = hw_pref_ads$H.C
hw_pref_ads$New = as.factor(hw_pref_ads$New)
hw_pref_ads$Adsorbed = factor(hw_pref_ads$Adsorbed, 
                              levels = c("adsorbed","not adsorbed"))

attach(hw_pref_ads)

names(hw_postf_ads)
hw_postf_ads$OC = hw_postf_ads$O.C
hw_postf_ads$HC = hw_postf_ads$H.C
hw_postf_ads$New = as.factor(hw_postf_ads$New)
hw_postf_ads$Adsorbed = as.factor(hw_postf_ads$Adsorbed)
attach(hw_postf_ads)

names(sw_pref_ads)
sw_pref_ads$OC = sw_pref_ads$O.C
sw_pref_ads$HC = sw_pref_ads$H.C
sw_pref_ads$New = as.factor(sw_pref_ads$New)
sw_pref_ads$Adsorbed = as.factor(sw_pref_ads$Adsorbed)
attach(sw_pref_ads)

names(sw_postf_ads)
sw_postf_ads$OC = sw_postf_ads$O.C
sw_postf_ads$HC = sw_postf_ads$H.C
sw_postf_ads$New = as.factor(sw_postf_ads$New)
sw_postf_ads$Adsorbed = as.factor(sw_postf_ads$Adsorbed)
attach(sw_postf_ads)


## van krevelen for adsorbed vs. not adsorbed ----
## hw ----
# hw preFenton adsorbed vs. not adsorbed
vankrev_hwpref=ggplot(hw_pref_ads,
                      aes (x=OC,y=HC,
                           color=Adsorbed,
                           shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1,4))+ #,breaks = c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #((breaks=c("adsorbed","not adsorbed")))+
  xlab("O/C")+
  ylab("H/C")+
   
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  ggtitle("HW preF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
    
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="top")+
    
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpref

# hw preFenton new compounds
vankrev_hwpref_new = ggplot(hw_pref_ads,
                            aes(x=OC,y=HC,
                                color=New,shape=New))+
  geom_point()+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
 # scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
 # scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="none")+
  
 # annotate("text", label = "HW, preF Goethite NEW", x = 0.1, y = 2.4)+ 
  ggtitle("New compounds: HW preF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpref_new


# hw postFenton adsorbed vs. not adsorbed
vankrev_hwpostf=ggplot(hw_postf_ads,
                       aes(x=OC,y=HC,
                           color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="top")+
  
  #annotate("text", label = "HW, postF", x = 0.1, y = 2.4)+ 
  ggtitle("HW postF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpostf

# hw postFenton new compounds
vankrev_hwpostf_new = ggplot(hw_postf_ads,
                             aes(x=OC,y=HC,
                                 color=New,shape=New))+
  geom_point()+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
  # scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="none")+
  
  #annotate("text", label = "HW, postF Goethite NEW", x = 0.1, y = 2.4)+ 
  ggtitle("HW postF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpostf_new

## sw ----

# sw preFenton adsorbed vs. not adsorbed
vankrev_swpref=ggplot(sw_pref_ads,
                      aes (x=OC,y=HC,
                           color=Adsorbed,
                           shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="top")+
  
  #annotate("text", label = "sw, preF", x = 0.1, y = 2.4)+ 
  ggtitle("SW preF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpref

# sw preFenton new compounds
vankrev_swpref_new = ggplot(sw_pref_ads,
                            aes(x=OC,y=HC,
                                color=New,shape=New))+
  geom_point()+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
  # scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="none")+
  
  # annotate("text", label = "sw, preF Goethite NEW", x = 0.1, y = 2.4)+ 
  ggtitle("SW preF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpref_new


# sw postFenton adsorbed vs. not adsorbed
vankrev_swpostf=ggplot(sw_postf_ads,
                       aes(x=OC,y=HC,
                           color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="top")+
  
  #annotate("text", label = "sw, postF", x = 0.1, y = 2.4)+ 
  ggtitle("SW postF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpostf

# sw postFenton new compounds
vankrev_swpostf_new = ggplot(sw_postf_ads,
                             aes(x=OC,y=HC,
                                 color=New,shape=New))+
  geom_point()+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual(values=c(1,4))+ #,breaks=c("adsorbed","not adsorbed"))+
  # scale_color_discrete()+ #breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="none")+
  
  #annotate("text", label = "sw, postF Goethite NEW", x = 0.1, y = 2.4)+ 
  ggtitle("SW postF")+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpostf_new

plot_grid (vankrev_hwpref,
          vankrev_hwpostf,
          vankrev_swpref,
          vankrev_swpostf,
          vankrev_hwpref_new,
          vankrev_hwpostf_new,
          vankrev_swpref_new,
          vankrev_swpostf_new,
          nrow = 2,ncol = 4,
          align="hv")

#
#### ------------------------------

### old ----
FTICRMS_by_treatment <- read_excel("FTICRMS_by treatment.xlsx", sheet = "Summary")
View(FTICRMS_by_treatment)
attach(FTICRMS_by_treatment)
names(FTICRMS_by_treatment)

## testing for normality
shapiro.test(sqrt(C))
qqnorm(sqrt(C))
qqline(log(C))


#### creating summary tables for each variable ####

library(dplyr)
library(tidyr)

df1=data.frame(FTICRMS_by_treatment[Forest=='HW'&Fenton=='PreF'&Goethite=="PreG",])
df2=data.frame(FTICRMS_by_treatment[Forest=='HW'&Fenton=='PostF'&Goethite=="PreG",])
df3=data.frame(FTICRMS_by_treatment[Forest=='SW'&Fenton=='PreF'&Goethite=="PreG",])
df4=data.frame(FTICRMS_by_treatment[Forest=='SW'&Fenton=='PostF'&Goethite=="PreG",])
df5=data.frame(FTICRMS_by_treatment[Forest=='HW'&Fenton=='PreF'&Goethite=="PostG",])
df6=data.frame(FTICRMS_by_treatment[Forest=='HW'&Fenton=='PostF'&Goethite=="PostG",])
df7=data.frame(FTICRMS_by_treatment[Forest=='SW'&Fenton=='PreF'&Goethite=="PostG",])
df8=data.frame(FTICRMS_by_treatment[Forest=='SW'&Fenton=='PostF'&Goethite=="PostG",])


##elements

cols=c('C','H','O','N','S','P','mz','HC','OC','AImod','DBE')
stargazer(df1[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df2[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df3[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df4[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df5[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df6[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df7[,cols],type="text",summary.stat=c("mean","sd"))
stargazer(df8[,cols],type="text",summary.stat=c("mean","sd"))


##groups
cols2=c('PolyCyPer','AromPer','UnsatLigPer','UnsatNoNPer','SatFarAcPer','UnsatNPer')
stargazer(df1[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df2[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df3[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df4[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df5[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df6[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df7[,cols2],type="text",summary.stat=c("mean","sd"))
stargazer(df8[,cols2],type="text",summary.stat=c("mean","sd"))

#### creating summary tables for adsorbed vs. not adsorbed -- using stargazer ####
hw_preF = read_excel("Adsorption compounds.xlsx", sheet = "hw_preF")
sw_preF = read_excel("Adsorption compounds.xlsx", sheet = "sw_preF")
hw_postF = read_excel("Adsorption compounds.xlsx", sheet = "hw_postF")
sw_postF = read_excel("Adsorption compounds.xlsx", sheet = "sw_postF")

library(dplyr)
library(tidyr)

cols=c('C','H','O','N','S','P','mz','HC','OC','AImod','DBE')

  attach(hw_preF)

df_hwpref_ads=data.frame(hw_preF[adsorbed=="adsorbed",])
df_hwpref_notads=data.frame(hw_preF[adsorbed=="not adsorbed",])

hwpref_ads_summary=stargazer(df_hwpref_ads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(hwpref_ads_summary,file="hwpref_ads_summary.csv")
hwpref_notads_summary=stargazer(df_hwpref_notads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(hwpref_notads_summary,file="hwpref_notads_summary.txt")

  attach(hw_postF)

df_hwpostf_ads=data.frame(hw_postF[adsorbed=="adsorbed",])
df_hwpostf_notads=data.frame(hw_postF[adsorbed=="not adsorbed",])

hwpostf_ads_summary=stargazer(df_hwpostf_ads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(hwpostf_ads_summary,file="hwpostf_ads_summary.csv")
hwpostf_notads_summary=stargazer(df_hwpostf_notads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(hwpostf_notads_summary,file="hwpostf_notads_summary.csv")

  attach(w_preF)

df_swpref_ads=data.frame(sw_preF[adsorbed=="adsorbed",])
df_swpref_notads=data.frame(sw_preF[adsorbed=="not adsorbed",])

swpref_ads_summary=stargazer(df_swpref_ads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(swpref_ads_summary,file="swpref_ads_summary.csv")
swpref_notads_summary=stargazer(df_swpref_notads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(swpref_notads_summary,file="swpref_notads_summary.csv")


  attach(sw_postF)

df_swpostf_ads=data.frame(sw_postF[adsorbed=="adsorbed",])
df_swpostf_notads=data.frame(sw_postF[adsorbed=="not adsorbed",])

swpostf_ads_summary=stargazer(df_swpostf_ads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(swpostf_ads_summary,file="swpostf_ads_summary.csv")
swpostf_notads_summary=stargazer(df_swpostf_notads[,cols],type="text",summary.stat=c("mean","sd"))
write.csv(swpostf_notads_summary,file="swpostf_notads_summary.csv")

#

## testing ANOVA for elements -- untransformed ####
## 1. for forest type 
anova1c= anova(lm(C~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1c
anova1h= anova(lm(H~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1h
anova1n= anova(lm(N~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1n
anova1o= anova(lm(O~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1o
anova1s= anova(lm(S~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1s
anova1p= anova(lm(P~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1p

anova1mz= anova(lm(mz~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1mz
anova1hc= anova(lm(HC~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1hc
anova1oc= anova(lm(OC~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1oc
anova1ai= anova(lm(AImod~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1ai
anova1dbe= anova(lm(DBE~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1dbe


## 2. for Fenton: HW
anova2c= anova(lm(C~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2c
anova2h= anova(lm(H~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2h
anova2n= anova(lm(N~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2n
anova2o= anova(lm(O~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2o
anova2s= anova(lm(S~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2s
anova2p= anova(lm(P~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2p
anova2mz= anova(lm(mz~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2mz
anova2hc= anova(lm(HC~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2hc
anova2oc= anova(lm(OC~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2oc
anova2ai= anova(lm(AImod~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2ai
anova2dbe= anova(lm(DBE~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2dbe

## 3. for Fenton: SW
anova3c= anova(lm(C~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3c
anova3h= anova(lm(H~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3h
anova3n= anova(lm(N~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3n
anova3o= anova(lm(O~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3o
anova3s= anova(lm(S~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3s
anova3p= anova(lm(P~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3p
anova3mz= anova(lm(mz~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3mz
anova3hc= anova(lm(HC~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3hc
anova3oc= anova(lm(OC~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3oc
anova3ai= anova(lm(AImod~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3ai
anova3dbe= anova(lm(DBE~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3dbe


## testing ANOVA for groups -- untransformed ####
## 1. for forest type 
anova1grp1= anova(lm(PolyCyPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp1
anova1grp2= anova(lm(AromPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp2
anova1grp3= anova(lm(UnsatLigPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp3
anova1grp4= anova(lm(UnsatNoNPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp4
anova1grp5= anova(lm(SatFarAcPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp5
anova1grp6= anova(lm(UnsatNPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp6

## 2. for Fenton: HW
anova2grp1= anova(lm(PolyCyPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp1
anova2grp2= anova(lm(AromPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp2
anova2grp3= anova(lm(UnsatLigPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp3
anova2grp4= anova(lm(UnsatNoNPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp4
anova2grp5= anova(lm(SatFarAcPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp5
anova2grp6= anova(lm(UnsatNPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp6

## 3. for Fenton: SW
anova3grp1= anova(lm(PolyCyPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp1
anova3grp2= anova(lm(AromPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp2
anova3grp3= anova(lm(UnsatLigPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp3
anova3grp4= anova(lm(UnsatNoNPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp4
anova3grp5= anova(lm(SatFarAcPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp5
anova3grp6= anova(lm(UnsatNPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp6




## testing ANOVA for elements -- log transformed ####
## 1. for forest type 
anova1c= anova(lm(log(C)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1c
anova1h= anova(lm(log(H)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1h
anova1n= anova(lm(log(N)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1n
anova1o= anova(lm(log(O)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1o
anova1s= anova(lm(log(S)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1s
anova1p= anova(lm(log(P)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1p

anova1mz= anova(lm(log(mz)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1mz
anova1hc= anova(lm(log(HC)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1hc
anova1oc= anova(lm(log(OC)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1oc
anova1ai= anova(lm(log(AImod)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1ai
anova1dbe= anova(lm(log(DBE)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1dbe


## 2. for Fenton: HW
anova2c= anova(lm(log(C)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2c
anova2h= anova(lm(log(H)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2h
anova2n= anova(lm(log(N)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2n
anova2o= anova(lm(log(O)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2o
anova2s= anova(lm(log(S)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2s
anova2p= anova(lm(log(P)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2p
anova2mz= anova(lm(log(mz)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2mz
anova2hc= anova(lm(log(HC)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2hc
anova2oc= anova(lm(log(OC)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2oc
anova2ai= anova(lm(log(AImod)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2ai
anova2dbe= anova(lm(log(DBE)~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2dbe

## 3. for Fenton: SW
anova3c= anova(lm(log(C)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3c
anova3h= anova(lm(log(H)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3h
anova3n= anova(lm(log(N)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3n
anova3o= anova(lm(log(O)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3o
anova3s= anova(lm(log(S)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3s
anova3p= anova(lm(log(P)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3p
anova3mz= anova(lm(log(mz)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3mz
anova3hc= anova(lm(log(HC)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3hc
anova3oc= anova(lm(log(OC)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3oc
anova3ai= anova(lm(log(AImod)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3ai
anova3dbe= anova(lm(log(DBE)~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3dbe


## testing ANOVA for groups -- log transformed ####
## 1. for forest type 
anova1grp1= anova(lm(PolyCyPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp1
anova1grp2= anova(lm(AromPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp2
anova1grp3= anova(lm(UnsatLigPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp3
anova1grp4= anova(lm(UnsatNoNPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp4
anova1grp5= anova(lm(SatFarAcPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp5
anova1grp6= anova(lm(UnsatNPer~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',])); anova1grp6

## 2. for Fenton: HW
anova2grp1= anova(lm(PolyCyPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp1
anova2grp2= anova(lm(AromPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp2
anova2grp3= anova(lm(UnsatLigPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp3
anova2grp4= anova(lm(UnsatNoNPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp4
anova2grp5= anova(lm(SatFarAcPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp5
anova2grp6= anova(lm(UnsatNPer~Fenton,data=FTICRMS_by_treatment[Forest=='HW'&Goethite=='PreG',])); anova2grp6

## 3. for Fenton: SW
anova3grp1= anova(lm(PolyCyPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp1
anova3grp2= anova(lm(AromPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp2
anova3grp3= anova(lm(UnsatLigPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp3
anova3grp4= anova(lm(UnsatNoNPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp4
anova3grp5= anova(lm(SatFarAcPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp5
anova3grp6= anova(lm(UnsatNPer~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',])); anova3grp6







## testing Wilcoxon for elements ####
wilcox.test(Forest,C)
wilcox.test(C~Forest)

wilcox1c= wilcox.test(C~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1c
wilcox1h= wilcox.test(H~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1h
wilcox1n= wilcox.test(N~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1n
wilcox1o= wilcox.test(O~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1o
wilcox1s= wilcox.test(S~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1s
wilcox1p= wilcox.test(P~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1p


wilcox1mz= kruskal.test(mz~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1mz
wilcox1hc= wilcox.test(log(HC)~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1hc
wilcox1oc= wilcox.test(OC~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1oc
wilcox1ai= wilcox.test(AImod~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1ai
wilcox1dbe= wilcox.test(DBE~Forest,data=FTICRMS_by_treatment[Fenton=='PreF'&Goethite=='PreG',]); wilcox1dbe

wilcox3mz= wilcox.test(mz~Fenton,data=FTICRMS_by_treatment[Forest=='SW'&Goethite=='PreG',]); wilcox3mz



## new compounds post-Fenton####
library(readxl)
newcompounds <- read_excel("FTICRMS_by treatment.xlsx",sheet = "New compounds")
names(newcompounds)
newcompounds$`HC`=newcompounds$`H/C`; newcompounds$`OC`=newcompounds$`O/C`
attach(newcompounds)
library(ggplot2)

vankrev_new=ggplot(newcompounds,aes(x=OC,y=HC,color=Vegetation,shape=Vegetation))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("hardwood","softwood"))+
  scale_color_discrete(breaks=c("hardwood","softwood"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "postF new compounds", x = 0.1, y = 2.5)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_new

#
## adsorption van krevelen####
library(readxl)
ads_hwpref = read_excel("Adsorption compounds.xlsx",sheet = "hw_preF")
ads_swpref = read_excel("Adsorption compounds.xlsx",sheet = "sw_preF")
ads_hwpostf = read_excel("Adsorption compounds.xlsx",sheet = "hw_postF")
ads_swpostf = read_excel("Adsorption compounds.xlsx",sheet = "sw_postF")

ads_hwpref2=data.frame(ads_hwpref[adsorbed2=="adsorbed",])

  attach(TrialAdsHWPreF)
  names(TrialAdsHWPreF)

ads_hwpref$`HC`=ads_hwpref$`H/C`; ads_hwpref$`OC`=ads_hwpref$`O/C`
ads_swpref$`HC`=ads_swpref$`H/C`; ads_swpref$`OC`=ads_swpref$`O/C`
ads_hwpostf$`HC`=ads_hwpostf$`H/C`; ads_hwpostf$`OC`=ads_hwpostf$`O/C`
ads_swpostf$`HC`=ads_swpostf$`H/C`; ads_swpostf$`OC`=ads_swpostf$`O/C`

names(ads_hwpref)
library(ggplot2)

vankrev_hwpref=ggplot(ads_hwpref,aes(x=OC,y=HC,color=adsorbed2,shape=adsorbed2))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "HW, preF", x = 0.05, y = 2.4)+ 
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpref

vankrev_swpref=ggplot(ads_swpref,aes(x=OC,y=HC,color=adsorbed2,shape=adsorbed2))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "SW, preF", x = 0.05, y = 2.4)+ 
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpref

vankrev_hwpostf=ggplot(ads_hwpostf,aes(x=OC,y=HC,color=adsorbed2,shape=adsorbed2))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "HW, postF", x = 0.05, y = 2.4)+ 
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpostf

vankrev_swpostf=ggplot(ads_swpostf,aes(x=OC,y=HC,color=adsorbed2,shape=adsorbed2))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "SW, postF", x = 0.05, y = 2.4)+ 
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpostf

library(cowplot)
plot_grid(vankrev_hwpref,vankrev_swpref,vankrev_hwpostf,vankrev_swpostf,labels=c("A","B","C","D"),align="hv")


## adsorption van Krevelen clean ####
library(readxl)
ads_hwpref = read_excel("Adsorption compounds.xlsx",sheet = "hw_preF")
ads_swpref = read_excel("Adsorption compounds.xlsx",sheet = "sw_preF")
ads_hwpostf = read_excel("Adsorption compounds.xlsx",sheet = "hw_postF")
ads_swpostf = read_excel("Adsorption compounds.xlsx",sheet = "sw_postF")
ads_hwpref$`HC`=ads_hwpref$`H/C`; ads_hwpref$`OC`=ads_hwpref$`O/C`
ads_swpref$`HC`=ads_swpref$`H/C`; ads_swpref$`OC`=ads_swpref$`O/C`
ads_hwpostf$`HC`=ads_hwpostf$`H/C`; ads_hwpostf$`OC`=ads_hwpostf$`O/C`
ads_swpostf$`HC`=ads_swpostf$`H/C`; ads_swpostf$`OC`=ads_swpostf$`O/C`

library(Rmisc)
ads_hwpreF_summary=summarySE(ads_hwpref,measurevar = "adsnumeric",groupvars=c("HC","OC"))
ads_hwpostF_summary=summarySE(ads_hwpostf,measurevar = "adsnumeric",groupvars=c("HC","OC"))
ads_swpreF_summary=summarySE(ads_swpref,measurevar = "adsnumeric",groupvars=c("HC","OC"))
ads_swpostF_summary=summarySE(ads_swpostf,measurevar = "adsnumeric",groupvars=c("HC","OC"))

write.csv(ads_hwpreF_summary,file = "ads_hwpreFsummary.csv")
write.csv(ads_hwpostF_summary,file = "ads_hwpostFsummary.csv")
write.csv(ads_swpreF_summary,file = "ads_swpreFsummary.csv")
write.csv(ads_swpostF_summary,file = "ads_swpostFsummary.csv")

## compile the summaries into Excel, code 1 = adsorbed, 100 = not adsorbed

#ads_hwpref2=data.frame(ads_hwpref[adsorbed2=="adsorbed",])
ads_clean = read_excel("Adsorption compounds.xlsx",sheet = "ads_summary")
names(ads_clean)
attach(ads_clean)

library(ggplot2)

vankrev_hwpref=ggplot(ads_clean[Treatment=="hw_preF",],aes(x=OC,y=HC,color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpref

vankrev_swpref=ggplot(ads_clean[Treatment=="sw_preF",],aes(x=OC,y=HC,color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "SW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpref


vankrev_hwpostf=ggplot(ads_clean[Treatment=="hw_postF",],aes(x=OC,y=HC,color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "HW, postF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_hwpostf

vankrev_swpostf=ggplot(ads_clean[Treatment=="sw_postF",],aes(x=OC,y=HC,color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "SW, postF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_swpostf

library(cowplot)
plot_grid(vankrev_hwpref,vankrev_swpref,vankrev_hwpostf,vankrev_swpostf,labels=c("A","B","C","D"),align="hv")

## van Krevelen for HW-SW combined
#preFenton

vankrev_pref=ggplot(ads_clean[Fenton=="preF",],aes(x=OC,y=HC,color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "pre-Fenton", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_pref


vankrev_postf=ggplot(ads_clean[Fenton=="postF",],aes(x=OC,y=HC,color=Adsorbed,shape=Adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("adsorbed","not adsorbed"))+
  scale_color_discrete(breaks=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  annotate("text", label = "post-Fenton", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  #annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  #annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  #annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  #annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  #annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_postf

library(cowplot)
plot_grid(vankrev_pref,vankrev_postf,labels=c("A","B"),align="hv")


#
## van krev for vegetation ####
library(readxl)
veg_vankrev = read_excel("Adsorption compounds.xlsx",sheet = "Veg_pref")

veg_vankrev$`HC`=veg_vankrev$`H/C`; veg_vankrev$`OC`=veg_vankrev$`O/C`


names(veg_vankrev)
library(ggplot2)

vankrev_veg=ggplot(veg_vankrev,aes(x=OC,y=HC,color=Vegetation,shape=Vegetation))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(1,4),breaks=c("Hardwood","Softwood"))+
  scale_color_discrete(breaks=c("Hardwood","Softwood"))+
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));vankrev_veg

#
## adsorption by group ####
library(readxl)
Adsorption_compounds = read_excel("Adsorption compounds.xlsx",sheet = "ads_grps_AVERAGE")
attach(Adsorption_compounds)
names(Adsorption_compounds)
str(Adsorption_compounds)
Adsorption_compounds$Group=factor(Adsorption_compounds$Group,levels=c("Cond.Ar","Aromatic","Lignin-like","Carb-like","Aliph-N","Aliph-noN"))

library(Rmisc)
grp=summarySEwithin(Adsorption_compounds, measurevar = "average", withinvars=c("Group","Adsorbed"))
grp_hw=summarySEwithin(Adsorption_compounds[Forest=="HW",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)
grp_sw=summarySEwithin(Adsorption_compounds[Forest=="SW",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)
grp_preF=summarySEwithin(Adsorption_compounds[Fenton=="preF",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)
grp_postF=summarySEwithin(Adsorption_compounds[Fenton=="postF",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)

grp
grp_hw
grp_sw

# bar plots for O horizon contents
library(ggplot2)
g1=ggplot(grp,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "All", x = 6, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,45)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g1

g2=ggplot(grp_hw,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  annotate("text", label = "*", x = 3, y = 55,size=10)+  # add significance letters
  annotate("text", label = "Hardwoods", x = 5, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,60)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g2

g3=ggplot(grp_sw,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  annotate("text", label = "*", x = 3, y = 50,size=10)+  # add significance letters
  annotate("text", label = "Softwoods", x = 5, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,60)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g3

library(cowplot)
plot_grid(g2,g3,g1,labels = c("A","B","C","D"),align="hv")

g4=ggplot(grp_preF,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  annotate("text", label = "pre-Fenton", x = 5, y = 44)+
  annotate("text", label = "*", x = 3, y = 55,size=10)+  # add significance letters
  annotate("text", label = "*", x = 4, y = 10,size=10)+
  annotate("text", label = "*", x = 5, y = 25,size=10)+
  annotate("text", label = "*", x = 6, y = 25,size=10)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,60)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g4

g5=ggplot(grp_postF,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  annotate("text", label = "post-Fenton", x = 5, y = 44)+
  annotate("text", label = "*", x = 2, y = 10,size=10)+
  annotate("text", label = "*", x = 3, y = 50,size=10)+
  annotate("text", label = "*", x = 4, y = 25,size=10)+
  
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,60)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g5

library(cowplot)
plot_grid(g4,g5,labels = c("A","B"),align="hv")


## anova for adsorption groups ####
#don't do Wilcoxon, use ANOVA

#prefenton
anova_ads1= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Fenton=='preF'&Group=='Cond.Ar',]); anova_ads1
anova_ads2= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Fenton=='preF'&Group=='Aromatic',]); anova_ads2
anova_ads3= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='preF'&Group=='Lignin-like',])); anova_ads3
anova_ads4= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='preF'&Group=='Carb-like',])); anova_ads4
anova_ads5= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='preF'&Group=='Aliph-N',])); anova_ads5
anova_ads6= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='preF'&Group=='Aliph-noN',])); anova_ads6

#postfenton
anova_ads7= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='postF'&Group=='Cond.Ar',])); anova_ads7
anova_ads8= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='postF'&Group=='Aromatic',])); anova_ads8
anova_ads9= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='postF'&Group=='Lignin-like',])); anova_ads9
anova_ads10= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='postF'&Group=='Carb-like',])); anova_ads10
anova_ads11= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='postF'&Group=='Aliph-N',])); anova_ads11
anova_ads12= anova(lm(average~Adsorbed,data=Adsorption_compounds[Fenton=='postF'&Group=='Aliph-noN',])); anova_ads12

#hw
anova_ads13= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='HW'&Group=='Cond.Ar',],exact=FALSE); anova_ads13
anova_ads14= anova(lm(average~Adsorbed,data=Adsorption_compounds[Vegetation=='HW'&Group=='Aromatic',])); anova_ads14
anova_ads15= anova(lm(average~Adsorbed,data=Adsorption_compounds[Vegetation=='HW'&Group=='Lignin-like',])); anova_ads15
anova_ads16= anova(lm(average~Adsorbed,data=Adsorption_compounds[Vegetation=='HW'&Group=='Carb-like',])); anova_ads16
anova_ads17= anova(lm(average~Adsorbed,data=Adsorption_compounds[Vegetation=='HW'&Group=='Aliph-N',])); anova_ads17
anova_ads18= anova(lm(average~Adsorbed,data=Adsorption_compounds[Vegetation=='HW'&Group=='Aliph-noN',])); anova_ads18

#sw
wilcox_ads19= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='SW'&Group=='Cond.Ar',],exact=FALSE); wilcox_ads19
wilcox_ads20= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='SW'&Group=='Aromatic',],exact=FALSE); wilcox_ads20
wilcox_ads21= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='SW'&Group=='Lignin-like',],exact=FALSE); wilcox_ads21
wilcox_ads22= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='SW'&Group=='Carb-like',],exact=FALSE); wilcox_ads22
wilcox_ads23= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='SW'&Group=='Aliph-N',]); wilcox_ads23
wilcox_ads24= wilcox.test(average~Adsorbed,data=Adsorption_compounds[Vegetation=='SW'&Group=='Aliph-noN',]); wilcox_ads24

Adsorption_compounds[Vegetation=='SW'&Group=='Lignin-like',]
#
## adsorption by group NEW ####
# for all replicates, by forest and Fenton
library(readxl)
Adsorption_compounds = read_excel("Adsorption compounds.xlsx",sheet = "ads_grps_NEW")
attach(Adsorption_compounds)
names(Adsorption_compounds)
str(Adsorption_compounds)
Adsorption_compounds$Groups=factor(Adsorption_compounds$Groups,levels=c("Cond.Ar","Aromatic","Lignin-like","Carb-like","Aliph-N","Aliph-noN"))

library(Rmisc)
grp=summarySEwithin(Adsorption_compounds, measurevar = "average", withinvars=c("Group","Adsorbed"))
grp_hw=summarySEwithin(Adsorption_compounds[Forest=="HW",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)
grp_sw=summarySEwithin(Adsorption_compounds[Forest=="SW",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)
grp_preF=summarySEwithin(Adsorption_compounds[Fenton=="preF",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)
grp_postF=summarySEwithin(Adsorption_compounds[Fenton=="postF",], measurevar = "average", withinvars=c("Group","Adsorbed"),na.rm=TRUE)

grp_hwpreF=summarySEwithin(Adsorption_compounds[Forest=="HW"&Fenton=="preF",], measurevar = "Percentage", withinvars=c("Groups","Adsorbed"),na.rm=TRUE)
grp_swpreF=summarySEwithin(Adsorption_compounds[Forest=="SW"&Fenton=="preF",], measurevar = "Percentage", withinvars=c("Groups","Adsorbed"),na.rm=TRUE)
grp_hwpostF=summarySEwithin(Adsorption_compounds[Forest=="HW"&Fenton=="postF",], measurevar = "Percentage", withinvars=c("Groups","Adsorbed"),na.rm=TRUE)
grp_swpostF=summarySEwithin(Adsorption_compounds[Forest=="SW"&Fenton=="postF",], measurevar = "Percentage", withinvars=c("Groups","Adsorbed"),na.rm=TRUE)


grp
grp_hw
grp_sw

# bar plots 
library(ggplot2)
g1=ggplot(grp,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "All", x = 6, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,45)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g1

g2=ggplot(grp_hw,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "Hardwoods", x = 5, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,45)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g2

g3=ggplot(grp_sw,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "Softwoods", x = 5, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,45)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g3

library(cowplot)
plot_grid(g2,g3,g1,labels = c("A","B","C","D"),align="hv")

g4=ggplot(grp_preF,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "pre-Fenton", x = 5, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,45)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g4

g5=ggplot(grp_postF,aes(y=average,x=Group,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`average`-se, ymax=`average`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "post-Fenton", x = 5, y = 44)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,45)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g5

library(cowplot)
plot_grid(g4,g5,labels = c("A","B"),align="hv")

library(ggplot2)

g_hw_pref=ggplot(grp_hwpreF,aes(y=Percentage,x=Groups,fill=Adsorbed)) +
  geom_bar(stat="summary",width=0.5,position=position_dodge(0.6),color="black",size=1)+
  geom_errorbar(aes(ymin=`Percentage`-se, ymax=`Percentage`+se),width=0.2,position=position_dodge(0.6),color="black",size=1) +
  #annotate("text", label = "B", x = 0.85, y = 53)+  # add significance letters
  annotate("text", label = "HW preF", x = 6, y = 25)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x=expression(""),y=expression(bold("% of total intensity")))+
  ylim(0,30)+
  #ylab("TC (Mg ha"^-1* ")")+
  #xlab("")+
  theme(legend.position=c(0.2,0.8))+
  theme(legend.key = element_blank())+  #remove white box around legend symbols
  #scale_fill_manual(values=c("grey80","darkred"))+
  guides(fill=guide_legend(title=""))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,face="bold",color="black"),axis.title=element_text(size=14,face="bold",color="black"))+
  theme(axis.text.x = element_text(angle=45,hjust=1));g_hw_pref
