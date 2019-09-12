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
