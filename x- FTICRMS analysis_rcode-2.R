#   adsorptive fractionation of SOM
#   Kaizad Patel
---------------------- 
----------------------


## i. packages ---------------------- ####
library(readxl)
library(tibble)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(Rmisc)
#

source("0-packages.R")


## ii. load files ---------------------- ####

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

## 3. rename columns ---------------------- ####

hw_fenton %>% 
  dplyr::rename(OC = O.C) %>% 
  dplyr::rename(HC = H.C) %>% 
  mutate(O = as.factor(as.integer(O)))->
  hw_fenton

sw_fenton %>% 
  dplyr::rename(OC = O.C) %>% 
  dplyr::rename(HC = H.C) ->
  sw_fenton

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
                             aes (x=`OC`,y=`HC`,
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


vankrev_sw_fentonloss=ggplot(sw_fenton,aes(x=`OC`,y=`HC`,color=loss, shape=loss,na.rm=TRUE))+
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
       aes(x = O, y = PreFentonHW))+
  geom_point()


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
