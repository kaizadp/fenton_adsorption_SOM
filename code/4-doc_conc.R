#   Reactive oxygen species alter chemical composition and adsorptive fractionation of soil-derived organic matter
#   Kaizad Patel

#   DOC PROCESSING
#   use this script to process DOC concentration data

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


source("code/0-packages.R")

## step 1. extraction parameters -------------------------------------------------------------------------
g = 1 #weight of goethite used, units = g
mL = 30 # volume of solution used, mL

## step 2. calculating adsorbed C -------------------------------------------------------------------------
conc = read_excel("data/doc_concentrations.xlsx", sheet = "DOC conc") %>% 
  dplyr::rename(doc_mg_L = `DOC blank corrected, mg/L`) %>% 
  select(Sample_no, Forest, Fenton, Goethite, doc_mg_L) %>% 
  na.omit() %>% 
  spread(Goethite,doc_mg_L) %>% 
  dplyr::mutate(adsorbed_mg_L = PreGoethite-PostGoethite,
                adsorbed_mg_g = adsorbed_mg_L*(1/1000)*(mL/g))

ggplot(conc, aes(x = Forest, y = adsorbed_mg_g, color = Fenton))+
geom_point()

## step 3. summarizing -------------------------------------------------------------------------
conc_summary = 
  conc %>% 
  group_by(Forest, Fenton) %>% 
  dplyr::summarise(doc_mg_L = mean(PreGoethite),
                   doc_se = sd(PreGoethite)/sqrt(n()),
                   adsorbed_mgg = mean(adsorbed_mg_g),
                   adsorbed_se = sd(adsorbed_mg_g)/sqrt(n()))

conc %>% 
  dplyr::summarize(adsorbed_mean = mean(adsorbed_mg_g),
                   adsorbed_se = sd(adsorbed_mg_g)/sqrt(n()))

## step 4. statistics -------------------------------------------------------------------------
## 4.1. effect of Fenton reaction on DOC concentration
car::Anova(lm(log10(PreGoethite) ~ Fenton, data = conc %>% filter(Forest=="HW")), type="III")
car::Anova(lm(log10(PreGoethite) ~ Fenton, data = conc %>% filter(Forest=="SW")), type="III")
