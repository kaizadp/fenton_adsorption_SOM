## Adsorptive fractionation of SOM
##### output script for FTICR data

# Kaizad F. Patel
# September 2019

## note to self: use FACETS for multi-panel graphs where possible, since it greatly shortens the code. 

source("0-packages.R")

# 1. import files ----

fticr_data_fenton = read.csv("fticr/fticr_data_fenton.csv")
fticr_data_goethite = read.csv("fticr/fticr_data_goethite.csv")
fticr_data_2 = read.csv("fticr/fticr_data_master.csv")
fticr_data_3 = read.csv("fticr/fticr_data_master_longform.csv")
fticr_data_4_hcoc = read.csv("fticr/fticr_data_hcoc_percentile.csv")
fticr_data_nosc = read.csv("fticr/fticr_data_nosc.csv")
fticr_meta_subset = read.csv("fticr/fticr_meta_subset.csv")
fticr_meta_hcoc=read.csv("fticr/fticr_meta_hcoc.csv")
# fticr_meta_nosc=read.csv("fticr/fticr_meta_nosc.csv")


# merge files with the meta for van krevelen plots
fticr_data_fenton2 = merge(fticr_data_fenton, fticr_meta_subset,by = "Mass", all.x = T)
fticr_data_goethite2 = merge(fticr_data_goethite,fticr_meta_subset, by="Mass", all.x = T)
fticr_data_2_hcoc = merge(fticr_data_2,fticr_meta_hcoc, by="Mass", all.x = T)

fticr_data_goethite2$Fenton = factor(fticr_data_goethite2$Fenton, levels = c("PreFenton","PostFenton"))

# 2. baseline van Krevelen ----

## 2.1 
gg_vk_native=ggplot (fticr_data_2_hcoc[fticr_data_2$PreFenton>0,],
                     aes (x=`OC`,y=`HC`, color = Forest, shape = Forest, na.rm=TRUE))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,1))+
# scale_color_manual(values=c("blue","grey20"))+
  
#  stat_ellipse()+
  
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

gg_vk_native_marginal=ggMarginal(gg_vk_native,groupColour = TRUE,groupFill = TRUE)

save_plot("output/vankrev_nativesom.tiff", gg_vk_native_marginal, 
          base_width = 7, base_height = 7)
#

## 2.2 native van krevelen with percentile abundance 

gg_vk_native_perc = 
  ggplot(fticr_data_4_hcoc[fticr_data_4_hcoc$treatment=="PreFenton",], #[fticr_data_goethite2_elcomp$Fenton=="PreFenton"&fticr_data_goethite2_elcomp$Forest=="SW",],
       aes (x=OC,y=HC, color=perc))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual (values=c(1,4),limits = c("adsorbed","not adsorbed"))+
  scale_color_brewer(palette = "Reds")+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  #geom_hline(yintercept = 2.0,color="black",linetype="solid")+ #aliph
  #geom_hline(yintercept = 1.5,color="black",linetype="solid")+ #aliph
  geom_vline(xintercept = 0.3,color="black",linetype="solid")+
  geom_vline(xintercept = 0.6,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 1.5, xend = 1.2, yend = 1.5,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 2, xend = 1.2, yend = 2,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 1, xend = 1.2, yend = 0.75,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 0.7, xend = 1.2, yend = 0.5,color="black",linetype="solid")+
  
  annotate("text", label = "saturated aliphatic", x = 1, y = 2.25, size = 3)+ 
  annotate("text", label = "unsaturated aliphatic", x = 1, y = 1.75, size = 3)+ 
  annotate("text", label = "highly unsaturated", x = 1, y = 1.25, size = 3)+ 
  annotate("text", label = "aromatic", x = 1, y = 0.65, size = 3)+ 
  annotate("text", label = "condensed aromatic", x = 1, y = 0.25, size = 3)+ 
  
  
  annotate("text", label = "low-O", x = 0.125, y = 0.0)+ 
  annotate("text", label = "high-O", x = 1.0, y = 0.0)+ 
  
  # facet_wrap(~El_comp)+ #facet with two variables
  facet_grid(~Forest) + 
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

save_plot("output/fig1_vankrev_nativesom_perc.tiff", gg_vk_native_perc, 
          base_width = 10, base_height = 6)


## 2.1 van krevelen with relative abundances ----

# order the treatments
fticr_data_4_hcoc %>% 
  mutate(treatment=factor(treatment,
                          levels = c("PreFenton","PostFenton","PreFentonGoethite","PostFentonGoethite")))->
  fticr_data_4_hcoc


gg_vk_treatments_perc = 
  ggplot(fticr_data_4_hcoc, 
       aes (x=OC,y=HC, color=perc))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual (values=c(1,4),limits = c("adsorbed","not adsorbed"))+
  scale_color_brewer(palette = "Reds")+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_vline(xintercept = 0.3,color="black",linetype="dashed")+
  geom_vline(xintercept = 0.6,color="black",linetype="dashed")+
  
# annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
#  annotate("text", label = "aromatic", x = 1.0, y = 0.5)+ 
#  annotate("text", label = "low-O", x = 0.125, y = 0.0)+ 
#  annotate("text", label = "high-O", x = 1.0, y = 0.0)+ 
  
  # facet_wrap(~El_comp)+ #facet with two variables
  facet_grid(treatment~Forest) + 
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

save_plot("output/fig2_vankrev_treatments_perc.tiff", gg_vk_treatments_perc, 
          base_width = 7, base_height = 15)


#
# 3. Fenton van krevelen ---------------------- ## -----
## VK with facet
gg_vk_fentonloss =
  ggplot (na.omit(fticr_data_fenton2[!fticr_data_fenton2$loss=="conserved",]),
                              aes (x=`OC`,y=`HC`, color=loss, shape=loss,na.rm=TRUE))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,4),limits=c("gained","lost"))+
  #scale_color_discrete(limits=c("lost","gained"))+
  scale_color_brewer(palette = "Set1",limits=c("gained","lost"))+
  facet_wrap(~Forest)+
  #  stat_ellipse(aes(fticr_data_fenton2[fticr_data_fenton2$loss=="conserved",], color = loss))+
  
  xlab("O/C")+
  ylab("H/C")+

  #annotate("text", label = "HW", x = 0.1, y = 2.4)+ 
  #ggtitle("HW fenton")+
  
  geom_vline(xintercept = 0.3,color="black",linetype="solid")+
  geom_vline(xintercept = 0.6,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 1.5, xend = 1.2, yend = 1.5,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 2, xend = 1.2, yend = 2,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 1, xend = 1.2, yend = 0.75,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 0.7, xend = 1.2, yend = 0.5,color="black",linetype="solid")+
  
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

  save_plot("output/fig3_vankrev_fenton_loss.tiff", gg_vk_fentonloss, 
            base_width = 10, base_height = 6)
  

## individual
gg_vk_hw_fentonloss = ggplot (fticr_data_fenton2[fticr_data_fenton2$Forest=="HW" &
                                                   !fticr_data_fenton2$loss=="conserved",],
                            aes (x=`OC`,y=`HC`,
                                 color=loss, shape=loss,na.rm=TRUE))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,4),limits=c("lost","gained"))+
  scale_color_discrete(limits=c("lost","gained"))+
  
#  stat_ellipse(aes(fticr_data_fenton2[fticr_data_fenton2$loss=="conserved",], color = loss))+
 
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


gg_vk_sw_fentonloss = ggplot (fticr_data_fenton2[fticr_data_fenton2$Forest=="SW" &
                                                   !fticr_data_fenton2$loss=="conserved",],
                              aes(x=`OC`,y=`HC`,
                                         color=loss, shape=loss,na.rm=TRUE))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual(values=c(19,4),limits=c("lost","gained"))+
  scale_color_discrete(limits=c("lost","gained"))+

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

gg_vk_fentonloss = plot_grid(gg_vk_hw_fentonloss,gg_vk_sw_fentonloss, ncol = 2, align = "hv", axis = "bt")

save_plot("output/vankrev_fentonloss.tiff", gg_vk_fentonloss, 
          base_width = 15, base_height = 7)


## fenton vk with percentiles 


gg_vk_fenton_perc = 
  ggplot(fticr_data_4_hcoc[fticr_data_4_hcoc$treatment=="PostFenton",], 
         aes (x=OC,y=HC, color=perc))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual (values=c(1,4),limits = c("adsorbed","not adsorbed"))+
  scale_color_brewer(palette = "Reds")+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_vline(xintercept = 0.3,color="black",linetype="dashed")+
  geom_vline(xintercept = 0.6,color="black",linetype="dashed")+
  
  # annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  #  annotate("text", label = "aromatic", x = 1.0, y = 0.5)+ 
  #  annotate("text", label = "low-O", x = 0.125, y = 0.0)+ 
  #  annotate("text", label = "high-O", x = 1.0, y = 0.0)+ 
  
  # facet_wrap(~El_comp)+ #facet with two variables
  facet_grid(~Forest) + 
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

save_plot("output/fig2_vankrev_fenton_perc.tiff", gg_vk_fenton_perc, 
          base_width = 10, base_height = 6)


#

# 3.1 fenton characterization -----------
#fticr_data_fenton2
fticr_data_fenton_el = merge(fticr_data_fenton2, fticr_meta_elements, by = "Mass")
ggplot(fticr_data_fenton_el, aes(x = loss, y = O))+
  geom_boxplot()+
  facet_wrap(~Class)
  
  

#
# 4. van krevelen for adsorbed vs. not adsorbed ----

# 4.1 adsorbed/non-adsorbed using facet ----
gg_vk_adsorbed_facet = 
  ggplot(fticr_data_goethite2,
                                aes (x=OC,y=HC, color=adsorbed, shape=adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1,4),limits = c("adsorbed","not adsorbed"))+
  scale_color_discrete(limits=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 

  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  facet_grid(Fenton~Forest+El_comp)+ #facet with two variables
  
#  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
#  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
#  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
#  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
#  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme_bw()+
  theme(panel.grid=element_blank(),
  legend.position="top",
  strip.background = element_rect(colour="white", fill="white"), #facet formatting
  panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
  strip.text.x = element_text(size=12, face="bold"), #facet labels
  strip.text.y = element_text(size=12, face="bold"), #facet labels
  legend.title=element_blank(),
  legend.text=element_text(size=12),
  panel.border=element_rect(color="black",size=1.5),
  axis.text=element_text(size=12,color="black"),
  axis.title=element_text(size=14,color="black",face="bold"));gg_vk_adsorbed_facet

save_plot("output/vankrev_adsorbed.tiff", gg_vk_adsorbed_facet, 
          base_width = 10, base_height = 10)

#

# 4.2 post-adsorption new using facet ----
gg_vk_ads_new_facet = ggplot(fticr_data_goethite2,
                              aes (x=OC,y=HC, color=new, shape=new))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1),limits = c("new molecules"))+
  scale_color_discrete(limits=c("new molecules"))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
  facet_grid(Fenton~Forest)+ #facet with two variables
  
  #  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  #  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  #  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  #  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  #  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"));gg_vk_ads_new_facet

save_plot("output/vankrev_ads_newmolecules.tiff", gg_vk_ads_new_facet, 
          base_width = 10, base_height = 10)


#
# 4.2 post-adsorption new single panel ----
gg_vk_ads_new = 
  ggplot(fticr_data_goethite2[fticr_data_goethite2$new=="new molecules",],
                             aes (x=OC,y=HC, color=Fenton))+
  geom_point(size=1,stroke=1, alpha = 0.5)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
#  scale_shape_manual (values=c(1),limits = c("new molecules"))+
#  scale_color_discrete(limits=c("new molecules"))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph

  geom_vline(xintercept = 0.3,color="black",linetype="dashed")+ #aliph
  geom_vline(xintercept = 0.6,color="black",linetype="dashed")+ #aliph
  
#  facet_grid(Fenton~Forest)+ #facet with two variables
  
  #  annotate("text", label = "lignin", x = 0.25, y = 0.8)+ 
  #  annotate("text", label = "condensed aromatic", x = 0.3, y = 0.3)+ 
  #  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
  #  annotate("text", label = "aromatic", x = 1.0, y = 1.3)+ 
  #  annotate("text", label = "carbohydrate", x = 1.0, y = 2.3)+ 
  ##group boundaries from Ohno et al. 2014| doi: 10.1021/es405570c
  
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

gg_vk_ads_new_marginal=ggMarginal(gg_vk_ads_new,groupColour = TRUE,groupFill = TRUE)

save_plot("output/fig4_vankrev_ads_newmolecules.tiff", gg_vk_ads_new, 
          base_width = 6, base_height = 6)


#
# 4.3 adsorption fractionation ----

gg_vk_adsorp_frac = 
  ggplot(fticr_data_goethite_relabund[!fticr_data_goethite_relabund$sorption_frac=="minimal change",], 
       aes(x = OC, y = HC, color = sorption_frac))+
  geom_point(size=1,stroke=1, alpha = 0.5)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  # scale_shape_manual (values=c(1,4),limits = c("adsorbed","not adsorbed"))+
  scale_color_brewer(palette = "PuOr")+
  #scale_color_viridis_d(option = "inferno")+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
 
  geom_vline(xintercept = 0.3,color="black",linetype="solid")+
  geom_vline(xintercept = 0.6,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 1.5, xend = 1.2, yend = 1.5,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 2, xend = 1.2, yend = 2,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 1, xend = 1.2, yend = 0.75,color="black",linetype="solid")+
  geom_segment(x = 0.0, y = 0.7, xend = 1.2, yend = 0.5,color="black",linetype="solid")+
  
  
#  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
#  geom_vline(xintercept = 0.3,color="black",linetype="dashed")+
#  geom_vline(xintercept = 0.6,color="black",linetype="dashed")+
  
#  annotate("text", label = "aliphatic", x = 0.4, y = 2.4)+ 
#  annotate("text", label = "aromatic", x = 1.0, y = 0.5)+ 
#  annotate("text", label = "low-O", x = 0.125, y = 0.0)+ 
#  annotate("text", label = "high-O", x = 1.0, y = 0.0)+ 
  
  # facet_wrap(~El_comp)+ #facet with two variables
  facet_grid(Fenton~Forest) + 
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="right",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

save_plot("output/fig4a_vankrev_adsorbed_frac.tiff", gg_vk_adsorp_frac, 
          base_width = 10, base_height = 8)

#
# 4.3 adsorbed/non-adsorbed by El_comp ----

#trial
# filter only CHO and CHON compounds. because the others have a low occurrence

fticr_data_goethite2 %>% 
  filter(El_comp=="CHO"|El_comp=="CHON")->
  fticr_data_goethite2_elcomp

ggplot(fticr_data_goethite2_elcomp, #[fticr_data_goethite2_elcomp$Fenton=="PreFenton"&fticr_data_goethite2_elcomp$Forest=="SW",],
       aes (x=OC,y=HC, color=adsorbed, shape=adsorbed))+
  geom_point(size=1,stroke=1)+
  xlim(0.0,1.2)+
  ylim(0.0,2.5)+
  scale_shape_manual (values=c(1,4),limits = c("adsorbed","not adsorbed"))+
  scale_color_discrete(limits=c("adsorbed","not adsorbed"))+
  xlab("O/C")+
  ylab("H/C")+
  
  #annotate("text", label = "HW, preF", x = 0.1, y = 2.4)+ 
  
  geom_hline(yintercept = 1.5,color="black",linetype="dashed")+ #aliph
  geom_rect(xmin=0.1,xmax=0.67,ymin=0.7,ymax=1.5,fill=NA,color="black",linetype="dashed")+ #lignin
  geom_rect(xmin=0.0,xmax=0.67,ymin=0.2,ymax=0.7,fill=NA,color="black",linetype="dashed")+ #cond.ar
  geom_rect(xmin=0.67,xmax=1.2,ymin=1.5,ymax=2.4,fill=NA,color="black",linetype="dashed")+ #carb
  
 # facet_wrap(~El_comp)+ #facet with two variables
facet_grid(Fenton~Forest+El_comp) + 
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

#
#
# 4.3 adsorbed/non-adsorbed by O, N ----
fticr_data_goethite3 = merge(fticr_data_goethite_relabund_subset, fticr_meta_elements, by = "Mass", all.x = T)
fticr_data_goethite3_N = fticr_data_goethite3[fticr_data_goethite3$N>0,]

ggplot(fticr_data_goethite3_N[fticr_data_goethite3_N$sorption_frac=="most sorbed"|
                              fticr_data_goethite3_N$sorption_frac=="most unbound",],
        aes(y = PreGoethite, x = as.factor(O)))+
  geom_bar(stat = "identity")+
  facet_wrap(~sorption_frac)
  geom_boxplot()

fticr_data_goethite4 = merge(fticr_data_goethite_relabund, fticr_meta_elements, by = "Mass", all.x = T)
fticr_data_goethite4 = merge(fticr_data_goethite4, fticr_meta_subset, by = "Mass", all.x = T)

ggplot(fticr_data_goethite4[fticr_data_goethite4$sorption_frac=="most sorbed"|
                                  fticr_data_goethite4$sorption_frac=="most unbound",],
         aes(y = delta_abund, x = as.factor(O),  fill=Class))+
    geom_bar(stat = "identity")+
    facet_wrap(~Fenton)+
  scale_fill_brewer(palette = "Set1")+
  geom_hline(yintercept = 0)
geom_boxplot()
  
ggplot(fticr_data_goethite4[fticr_data_goethite4$Class=="Lignin" &
                            fticr_data_goethite4$sorption_frac=="most sorbed" |
                              fticr_data_goethite4$sorption_frac=="most unbound",],
       aes(y = delta_abund, x = as.factor(O), fill=sorption_frac))+
  geom_bar(stat = "identity")+
  facet_wrap(~Fenton)

ggplot(fticr_data_goethite4[fticr_data_goethite4$Class=="Lignin"|
                              fticr_data_goethite4$Class=="Protein"|
                              fticr_data_goethite4$Class=="Tannin"|
                              fticr_data_goethite4$Class=="Carb",],
       aes(y = preg_rel_abund, x = as.factor(O), fill = Fenton))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(~Class)+
  scale_x_discrete(breaks = c("0", "5","10","15","20","25","30"))


# create a summary dataframe of total relative intensity vs. number of O
fticr_data_goethite4 %>% 
  group_by(Fenton, O, Class) %>% 
  dplyr::summarize(intensity = sum(preg_rel_abund))->
  fticr_data_goethite4_summary

# use this file for graph of rel_abund vs. O
# rel_abund is as a fraction, so multipy by 100 for %
ggplot(fticr_data_goethite4_summary[fticr_data_goethite4_summary$Class=="Lignin"|
                                      fticr_data_goethite4_summary$Class=="Protein"|
                                      fticr_data_goethite4_summary$Class=="Tannin"|
                                      fticr_data_goethite4_summary$Class=="Carb",],
       aes(y = intensity*100, x = as.factor(O), fill = Fenton))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(~Class)+
  ylab("% intensity")+
  scale_x_discrete(breaks = c("0", "5","10","15","20","25","30"))


# to get median:
fticr_data_goethite4 %>% 
  group_by(sorption_frac) %>% 
  dplyr::summarise(median = median(as.numeric(O)))->
  median

median(as.numeric(fticr_fenton_loss_relabund_o$O))


ggplot(fticr_data_goethite4[fticr_data_goethite4$Class=="Lignin",], 
       aes(x = O, y = preg_rel_abund, fill = Fenton))+
  geom_bar(stat = "identity")+
  facet_wrap(~Fenton)

#
# 4.3 individual plots -- don't do ----


## hw prefenton

gg_vk_adsorbed_hw_pref = ggplot(fticr_data_goethite2[fticr_data_goethite2$Forest == "HW" &
                                                 fticr_data_goethite2$Fenton == "PreFenton",],
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

gg_vk_adsorbed_sw_pref = ggplot(fticr_data_goethite2[fticr_data_goethite2$Forest == "SW" &
                                                 fticr_data_goethite2$Fenton == "PreFenton",],
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

gg_vk_adsorbed_hw_postf = ggplot(fticr_data_goethite2[fticr_data_goethite2$Forest == "HW" &
                                                 fticr_data_goethite2$Fenton == "PostFenton",],
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

gg_vk_adsorbed_sw_postf = ggplot(fticr_data_goethite2[fticr_data_goethite2$Forest == "SW" &
                                                 fticr_data_goethite2$Fenton == "PostFenton",],
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

## sw 

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
# molecules sorbed vs. unbound ----

fticr_goethite_data_relabundance_long2 = fticr_goethite_data_relabundance_long[!fticr_goethite_data_relabundance_long$Class=="total",]

gg_relabund_adsorbed = 
  ggplot(fticr_goethite_data_relabundance_long2,
       aes(x = Class, y = relabund, fill = sorption_frac))+
  geom_bar(stat = "summary", position = position_dodge(width = 0.7), color = "black", width = 0.5)+
  #scale_fill_brewer(palette = "PuOr")+
  scale_fill_manual(values = c("chocolate3","mediumpurple3"))+
  facet_grid(Fenton~Forest)+
  ylab("% of total intensity")+
  xlab("")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title=element_text(size=14,color="black",face="bold"))

save_plot("output/fig5_adsorbed_relabund.tiff", gg_relabund_adsorbed, 
          base_width = 8, base_height = 7)
  
#

## fenton: loss/gain vs. O ----
fticr_fenton_loss_relabund_o = merge(fticr_fenton_loss_relabund, fticr_meta_elements, by = "Mass", all.x = T)
# select only lost and gained values
fticr_fenton_loss_relabund_o = fticr_fenton_loss_relabund_o[fticr_fenton_loss_relabund_o$loss=="lost"|
                                                              fticr_fenton_loss_relabund_o$loss=="gained",]
fticr_fenton_loss_relabund_o$O = as.factor(fticr_fenton_loss_relabund_o$O)


fticr_fenton_loss_relabund_o <- fticr_fenton_loss_relabund_o %>%
  mutate(rel_abund_new = ifelse(loss == "gained",
                                 rel_abund,
                                 -1*rel_abund))

gg_fenton_loss_relabund=
  ggplot(na.omit(fticr_fenton_loss_relabund_o), aes(x = O, y = rel_abund_new, fill = loss))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~Forest)+
  scale_x_discrete(breaks = c("0", "5","10","15","20","25","30"))+
  ylab("% of total intensity")+
  theme_bw()+
  theme(
    #panel.grid=element_blank(),
        legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold"))

save_plot("output/fig_fentonloss_relabund.tiff", gg_fenton_loss_relabund, 
          base_width = 8, base_height = 5)
#

  ggplot(na.omit(fticr_fenton_loss_relabund_o), aes(y = as.numeric(O), x = loss, fill = loss))+
  geom_boxplot()+
    stat()

fticr_fenton_loss_relabund_o %>% 
  group_by(loss, Forest) %>% 
  dplyr::summarise(median = median(as.numeric(O)))
  
  median(as.numeric(fticr_fenton_loss_relabund_o$O))

#

 # 5. NOSC graph ----
# graph with facets

fticr_data_nosc$goethite = factor(fticr_data_nosc$goethite, levels = c("pre-Goethite", "post-Goethite"))

gg_nosc =
  ggplot(fticr_data_nosc[fticr_data_nosc$goethite=="pre-Goethite",], 
         aes(x = NOSC, fill = Forest))+
  geom_histogram(binwidth = 0.10, position = "identity", alpha = 0.4, color = "black")+
  xlim(-2, 2)+
  ylim(0,75)+
  facet_grid(fenton~Forest)+
  
  theme_bw()+
  theme(legend.position="top",
        strip.background = element_rect(colour="white", fill="white"), #facet formatting
        panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
        strip.text.x = element_text(size=12, face="bold"), #facet labels
        strip.text.y = element_text(size=12, face="bold"), #facet labels
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.border=element_rect(color="black",size=1.5),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black",face="bold")); gg_nosc

save_plot("output/fig5_nosc_0.10.tiff", gg_nosc, 
          base_width = 10, base_height = 8)

#

### individual plots -- don't do ----


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
# 6. Number of O ----
fticr_data_goethite_oxygen = merge(fticr_data_goethite_relabund, fticr_meta_elements, by = "Mass", all.x = T)

# select only relevant columns
fticr_data_goethite_oxygen %>% 
  select(Mass, Forest, Fenton, PreGoethite, preg_rel_abund, adsorbed, new, O)->
  fticr_data_goethite_oxygen2

ggplot(fticr_data_goethite_oxygen2[!is.na(fticr_data_goethite_oxygen2$adsorbed),], aes(x = factor(O), y = preg_rel_abund, color = adsorbed))+
  geom_bar(stat = "identity")

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

