## Adsorptive fractionation of SOM
##### output script for FTICR data

# Kaizad F. Patel
# September 2019

source("0-packages.R")

# 1. import files ----

fticr_data_fenton = read_csv("fticr/fticr_data_fenton.csv",
                             col_types = cols(
                               Mass = col_number(),
                               Forest = col_factor(),
                               PreFenton = col_number(),
                               PostFenton = col_number(),
                               loss = col_factor()))

fticr_data_goethite = read.csv("fticr/fticr_data_goethite.csv")
fticr_data_2 = read.csv("fticr/fticr_data_master.csv")
fticr_meta_subset = read_csv("fticr/fticr_meta_subset.csv",
                             col_types = cols(
                               Mass = col_number(),
                               HC = col_number(),
                               OC = col_number(),
                               NOSC = col_number(),
                               Class = col_factor(),
                               KM_CH2 = col_number(),
                               KMD_CH2 = col_number(),
                               AI_0 = col_number()
                             ))

fticr_meta_subset$Mass = round(fticr_meta_subset$Mass,4)

# merge files with the meta for van krevelen plots
fticr_data_fenton2 = merge(fticr_data_fenton, fticr_meta_subset,by = "Mass", all.x = T)
fticr_data_fenton2 = cbind(fticr_data_fenton, fticr_meta_subset)
fticr_data_fenton3 = left_join(fticr_data_fenton,fticr_meta_subset, by = "Mass")

fticr_data_goethite = merge(fticr_meta_subset,fticr_data_goethite, all = TRUE)
fticr_data_2 = merge(fticr_meta_subset,fticr_data_2)


# 2. baseline van Krevelen ----

gg_vk_native=ggplot (fticr_data_2[fticr_data_2$PreFenton>0,],
                     aes (x=`OC`,y=`HC`, color = Forest, shape = Forest, na.rm=TRUE))+
  geom_point(size=1,stroke=1, alpha = 0.5)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("blue","grey3"))+
  
  stat_ellipse()+
  
  xlab("O/C")+
  ylab("H/C")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position=c(0.8,0.2))+
  
  #annotate("text", label = "HW", x = 0.1, y = 2.4)+ 
  ggtitle("Native DOC")+
  
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
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_native




# 3. Fenton van krevelen ---------------------- ## -----

gg_vk_hw_fentonloss = ggplot (fticr_data_fenton[fticr_data_fenton$Forest=="HW",],
                            aes (x=`OC`,y=`HC`,
                                 color=loss, shape=loss,na.rm=TRUE))+
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
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_hw_fentonloss


gg_vk_sw_fentonloss = ggplot (fticr_data_fenton[fticr_data_fenton$Forest=="SW",],
                              aes(x=`OC`,y=`HC`,
                                         color=loss, shape=loss,na.rm=TRUE))+
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
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_sw_fentonloss








#
# 4. van krevelen for adsorbed vs. not adsorbed ----
## hw prefenton

gg_vk_adsorbed_hw_pref = ggplot(fticr_data_goethite[fticr_data_goethite$Forest == "HW" &
                                                 fticr_data_goethite$Fenton == "PreFenton",],
                      aes (x=OC,y=HC,
                           color=adsorbed,
                           shape=adsorbed))+
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
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_adsorbed_hw_pref

## sw prefenton

gg_vk_adsorbed_sw_pref = ggplot(fticr_data_goethite[fticr_data_goethite$Forest == "SW" &
                                                 fticr_data_goethite$Fenton == "PreFenton",],
                           aes (x=OC,y=HC,
                                color=adsorbed,
                                shape=adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1,4))+ #,breaks = c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #((breaks=c("adsorbed","not adsorbed")))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
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
  
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position="top")+
  
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.border=element_rect(color="black",size=1.5))+
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_adsorbed_sw_pref

## hw postfenton

gg_vk_adsorbed_hw_postf = ggplot(fticr_data_goethite[fticr_data_goethite$Forest == "HW" &
                                                 fticr_data_goethite$Fenton == "PostFenton",],
                           aes (x=OC,y=HC,
                                color=adsorbed,
                                shape=adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1,4))+ #,breaks = c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #((breaks=c("adsorbed","not adsorbed")))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  ggtitle("HW postF")+
  
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
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_adsorbed_hw_postf

## sw postfenton

gg_vk_adsorbed_sw_postf = ggplot(fticr_data_goethite[fticr_data_goethite$Forest == "SW" &
                                                 fticr_data_goethite$Fenton == "PostFenton",],
                           aes (x=OC,y=HC,
                                color=adsorbed,
                                shape=adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1,4))+ #,breaks = c("adsorbed","not adsorbed"))+
  scale_color_discrete()+ #((breaks=c("adsorbed","not adsorbed")))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  ggtitle("SW postF")+
  
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
  theme(axis.text=element_text(size=12,color="black"),axis.title=element_text(size=14,color="black",face="bold"));gg_vk_adsorbed_sw_postf


###




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
# 5. NOSC graph ----
fticr_data_3 = read_csv("fticr_data_master_longform.csv")
# create nosc file by importing only Mass and nosc columns

fticr_meta_nosc = read_csv("fticr_meta_subset.csv",
                           col_types = cols_only(
                             Mass = col_number(), 
                             NOSC = col_number()))

# before merging, make sure "Mass" column has the same number of decimals

fticr_data_3 %>% 
  mutate(Mass = paste(round(Mass,4)))->
  fticr_data_3
fticr_meta_nosc %>% 
  mutate(Mass = paste(round(Mass,4)))->
  fticr_meta_nosc

fticr_data_4 = merge (fticr_meta_nosc, fticr_data_3, by = "Mass", all = T)

# convert all zeroes to NA
fticr_data_4 %>%
  mutate_all(~replace(., . == 0, NA))->
  fticr_data_4

fticr_data_4 = fticr_data_4[complete.cases(fticr_data_4),]

# graph
gg_nosc_initial

ggplot(fticr_data_4[fticr_data_4$treatment=="PreFenton",], 
         aes(x = NOSC, fill = Forest))+
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.5, color = "black")+
  xlim(-2.5, 2)+
  ylim(0,200)+
  theme(legend.position = "top")+
  ggtitle("PreFenton")

ggplot(fticr_data_4[fticr_data_4$treatment=="PostFenton",], 
       aes(x = NOSC, fill = Forest))+
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.5, color = "black")+
  xlim(-2.5, 2)+
  ylim(0,200)+
  theme(legend.position = "top")+
  ggtitle("PostFenton")

ggplot(fticr_data_4[fticr_data_4$treatment=="PreFentonGoethite",], 
       aes(x = NOSC, fill = Forest))+
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.5, color = "black")+
  xlim(-2.5, 2)+
  ylim(0,200)+
  theme(legend.position = "top")+
  ggtitle("PreFentonGoethite")

ggplot(fticr_data_4[fticr_data_4$treatment=="PostFentonGoethite",], 
       aes(x = NOSC, fill = Forest))+
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.5, color = "black")+
  xlim(-2.5, 2)+
  ylim(0,200)+
  theme(legend.position = "top")+
  ggtitle("PostFentonGoethite")

#
### ### ### ----
## intensity vs. O ---------------------- ## -----
hw_pref_o=summarySE(hw_fenton,measurevar = "PreFentonHW", groupvars=c("O"),na.rm=TRUE);hw_pref_o

ggplot(hw_pref_o,
       aes(x = O, y = PreFentonHW), geom_bar())


# 4. adsorbed vs. not adsorbed ----
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


#### ------------------------------

