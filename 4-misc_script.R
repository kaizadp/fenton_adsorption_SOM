master %>% 
  group_by(Treatment,Mass) %>% 
  dplyr::summarise(intensity = mean(intensity)) %>% 
  
  
  left_join(dplyr::select(meta_elements, O, Mass), by = "Mass") %>% 
  dplyr::filter(Treatment=="PreFenton"| Treatment=="PostFenton")->master_el


      #   Class             Treatment  median
      #   <fct>             <fct>       <dbl>
      # 1 Carbohydrate-like PostFenton     14
      # 2 Carbohydrate-like PreFenton      13
      # 3 Lignin-like       PostFenton     11
      # 4 Lignin-like       PreFenton      10


master_el %>% 
  left_join(class, by = "Mass") %>%
  filter(Class=="Carbohydrate-like"|Class=="Lignin-like") %>% 
  group_by(Treatment,Class,O) %>% 
  dplyr::summarise(count = n()) %>% 
 
  ungroup %>% 
  group_by(Treatment,Class) ->temp_master

ggplot(temp_master[temp_master$Class=="Carbohydrate-like"| temp_master$Class=="Lignin-like",], aes(x = O, y =count, fill = Treatment))+
  geom_bar(stat = "identity")+
  facet_grid(Treatment~ Class)+
  theme_kp()





ggplot(master_el, aes(x = O, fill = Treatment))+
  geom_histogram(binwidth = 1)+
  facet_wrap(Forest~Treatment)


ggplot(temp_master, aes(x = O, y = count, fill = Treatment))+
  geom_bar(stat = "identity")+
  facet_wrap(Treatment~.)


master_el %>% 
  group_by(Treatment) %>% 
  dplyr::mutate(median = median(O))->master_el

master_el %>% 
  group_by(Treatment) %>% 
  dplyr::summarise(n = n())


sorbed %>% 
  left_join(select(meta_elements,Mass,O), by = "Mass")->sorbed_temp

sorbed_temp %>% 
  group_by(adsorbed) %>% 
  dplyr::summarise(median = median(O))




ggplot(sorbed_temp[sorbed_temp$adsorbed=="sorbed"&(temp$Class=="Carbohydrate-like"|temp$Class=="Lignin-like"),], aes(x = O, fill = adsorbed))+
  geom_histogram(binwidth = 1)+
  facet_grid(Forest+fenton~adsorbed+Class)


fenton_el %>% 
  left_join(class, by = "Mass") %>% 
  group_by(Class, Forest, Treatment, loss) %>% 
  dplyr::summarize(median = median(O))->temp


ggplot(temp, aes(x = O, fill = loss))+
  geom_histogram()+
  facet_grid(Forest+loss~Class)




#### molecular analysis of fenton loss-gained molecules of carb and lignin


fenton_el %>%
  group_by(Mass,loss) %>% 
  dplyr::summarise(HC = mean(HC),
                   OC = mean(OC),
                   O = mean(O)) %>% 
  left_join(class, by = "Mass") %>% 
  filter(Class=="Carbohydrate-like"|Class=="Lignin-like") ->fenton_temp1

fenton_temp1 %>% 
  group_by(loss) %>% 
  dplyr::summarize(median = median(O),
                   n = n())


    #   Class             loss      median
    #   <fct>             <fct>      <dbl>
    # 1 Carbohydrate-like conserved     14
    # 2 Carbohydrate-like gained        14
    # 3 Carbohydrate-like lost          11
    # 4 Lignin-like       conserved     10
    # 5 Lignin-like       gained        14
    # 6 Lignin-like       lost           9


fenton_temp1%>% 
  group_by(Class,loss, O) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(count2 = if_else(loss=="lost",-1*count,as.numeric(count)))->fenton_temp2


ggplot(fenton_temp[!fenton_temp$loss=="conserved",], aes(x = O, fill = loss))+
  geom_histogram(binwidth = 1)+
  facet_grid(Class~loss)
  


ggplot(fenton_temp2[!fenton_temp2$loss=="conserved",], aes(x=as.factor(O), y=count2, fill=loss))+
  geom_bar(stat = "identity", alpha=0.8)+
  scale_fill_brewer(palette="Dark2")+
  scale_x_discrete(breaks = c("0", "5","10","15","20","25","30"))+
  ylab("number of peaks")+
  xlab("number of O")+
  facet_wrap(~Class)+
  theme_kp()


#### molecular analysis of sorbed
sorbed %>%
  group_by(Mass,fenton,adsorbed) %>% 
  dplyr::summarise(HC = mean(HC),
                   OC = mean(OC))%>% 
  left_join(class, by = "Mass") %>%
  left_join(select(meta_elements,Mass,O),  by = "Mass") %>% 
  filter(Class=="Carbohydrate-like"|Class=="Lignin-like") ->sorbed_temp1


sorbed_temp1%>% 
  group_by(Class,adsorbed, fenton, O) %>% 
  dplyr::summarise(count = n()) ->sorbed_temp2


ggplot(sorbed_temp2[sorbed_temp2$adsorbed=="sorbed",], aes(x=as.factor(OC), y=count, fill=fenton))+
  geom_bar(stat = "identity", alpha=0.8)+
  scale_fill_brewer(palette="Dark2")+
#  scale_x_discrete(breaks = c("0", "5","10","15","20","25","30"))+
 # scale_x_continuous(limits = c(0,30))+
  ylab("number of peaks")+
  xlab("number of O")+
  facet_grid(fenton~Class)+
  theme_kp()


sorbed_temp1 %>% 
  group_by(adsorbed, fenton,Class) %>% 
  dplyr::summarize(median = median(O)) 


### sorbed -- comparing only peaks seen in the original extracts

sorbed %>% 
  spread(fenton, adsorbed)%>% 
  filter(!is.na(PreFenton)) %>% 
  gather(key = "fenton", value ="sorbed", PreFenton:PostFenton) %>% 
  filter(Class=="Carbohydrate-like"|Class=="Lignin-like") ->sorbed_initial_temp


gg_vankrev(na.omit(sorbed_initial_temp), aes(x = OC, y=HC, color = sorbed))+
  facet_grid(fenton~Class)

###

## more meta 

fticr_meta = read.csv(FTICR_META)

fticr_meta %>% 
  mutate(NOSC = 4-(((4*C)+H-(3*N)-(2*O)-(2*S))/C))->
  fticr_meta


master2  = 
  master %>% 
  left_join(fticr_meta, by = "Mass")

ggplot(master2[master2$Treatment=="PreFenton"|master2$Treatment=="PostFenton",], 
       aes(x = NOSC, color = Treatment, fill = Treatment))+
  geom_histogram(binwidth = 0.1, alpha = 0.2,position = "identity")+
  facet_wrap(~Class)


ggplot(master2[master2$Treatment=="PreFenton"|master2$Treatment=="PostFenton",], 
       aes(x = KM, y = KMD, color = Treatment, fill = Treatment, shape = Treatment))+
  geom_point(alpha = 0.5, size = 2)+
  scale_shape_manual(values = c(1,16))+
  facet_wrap(~Class)

fentonloss2 = 
  fentonloss %>% left_join(fticr_meta, by  ="Mass")

ggplot(fentonloss2[!fentonloss2$loss=="conserved",],
       aes(x = KM, y = KMD, color = loss, fill = loss, shape = loss))+
  geom_point(alpha = 0.5, size = 2)+
  scale_shape_manual(values = c(1,16))+
  facet_wrap(~Class)

  
 ## peaks lost or gained
fentonloss %>% 
  group_by(Mass,loss) %>% 
  dplyr::summarize(n = mean(Mass)) %>% 
  left_join(class, by = "Mass")->loss

loss %>% 
  group_by(Class, loss) %>% 
  dplyr::summarise(n = n()) %>% 
  spread(loss, n)
  
sorbed2 = 
  sorbed %>% 
  left_join(fticr_meta, by = "Mass")

ggplot(sorbed2[sorbed2$adsorbed=="sorbed",], aes(x = NOSC, color = fenton, fill = fenton))+
  geom_histogram(binwidth = 0.1, alpha = 0.2,position = "identity")+
  facet_wrap(~Class.x)
  
ggplot(sorbed2[sorbed2$adsorbed=="sorbed",], aes(x = O, color = fenton, fill = fenton))+
  geom_histogram(binwidth = 0.1, alpha = 0.2,position = position_dodge(width=0.7) )+
  facet_wrap(~Class.x)


fenton_el %>% left_join(class, by = "Mass")->fenton_el

ggplot(fenton_el, aes(x = O, color = loss, fill = loss))+
  geom_histogram(alpha  = 0.2, position = "identity", binwidth = 1)+
  facet_grid(.~Class)
       

gg_vankrev(fenton_el, aes(x = OC, y = HC, color = loss, shape = loss))+
  geom_point(alpha = 0.1)




ggplot(fenton_el, aes(x = Mass, color = loss, fill = loss))+
  geom_histogram(alpha  = 0.2, position = "identity")+
  facet_grid(loss~Forest)




       



