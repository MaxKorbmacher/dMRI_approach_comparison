# check loadings (partly visually)
#
# read loading matrices
ukbl = read.csv("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/corrected_loadings_ukb.csv")
abcdl = read.csv("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/corrected_loadings_abcd.csv")
#
# highest loadings of fourth component
ukbl%>%filter(Comp==4) %>% filter(B==max(B))
test = ukbl%>%filter(Comp==4) # CG appears at multiple top positions, hence:
ukbl%>%filter(Comp==4) %>% filter(grepl("CG", Variable)) %>%
  summarise(mean(abs(B)), sd(abs(B)), median(abs(B)),mad(abs(B)))
test = ukbl%>%filter(Comp==5) # UF appears at multiple top positions, hence:
ukbl%>%filter(Comp==5) %>% filter(grepl("UF", Variable)) %>%
  summarise(mean(abs(B)), sd(abs(B)), median(abs(B)),mad(abs(B)))


test = abcdl%>%filter(Comp==4) #UF
abcdl%>%filter(Comp==4) %>% filter(grepl("UF", Variable)) %>%
  summarise(mean(abs(B)), sd(abs(B)), median(abs(B)),mad(abs(B)))

test = abcdl%>%filter(Comp==5) # CG_hippocampus
abcdl%>%filter(Comp==5) %>% filter(grepl("CG", Variable)) %>%
  summarise(mean(abs(B)), sd(abs(B)), median(abs(B)),mad(abs(B)))



# SMT
## UKB
test = ukbl%>%filter(Comp==3 & Approach=="SMT") #%>% filter(B==max(B))
test = ukbl%>%filter(Comp==4 & Approach=="SMT") #%>% filter(B==max(B))
## ABCD
test = abcdl%>%filter(Comp==3 & Approach=="SMT") #%>% filter(B==max(B))
test = abcdl%>%filter(Comp==4 & Approach=="SMT") #%>% filter(B==max(B))
#
# BRIA
test = ukbl%>%filter(Comp==4 & Approach=="BRIA") #%>% filter(B==max(B))
test = abcdl%>%filter(Comp==4 & Approach=="BRIA") #%>% filter(B==max(B))
#
# DKI
# ukb
ukbl%>%filter(Comp==3 & Approach=="DKI") %>% filter(grepl("CG", Variable) |
                                                             grepl("CST", Variable)) %>%
  summarise(mean(abs(B)), sd(abs(B)), median(abs(B)),mad(abs(B)))
# abcd
test = abcdl%>%filter(Comp==3 & Approach=="DKI")
abcdl%>%filter(Comp==3 & Approach=="DKI") %>% filter(grepl("SLF", Variable))%>%
  summarise(mean(abs(B)), sd(abs(B)), median(abs(B)),mad(abs(B)))
