# check loadings (partly visually)
# Max Korbmacher, 17 Jan 2024
#
#
# PREP ####
# read loading matrices
ukbl = read.csv("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/corrected_loadings_ukb.csv")
abcdl = read.csv("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/corrected_loadings_abcd.csv")

ukbl_test = read.csv("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/corrected_loadings_ukb_TEST.csv")
abcdl_test = read.csv("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/corrected_loadings_abcd_TEST.csv")
#
#
# Context: SEX CLASSIFICATION ####
# ########### #
#
# BRIA ####
# highest loadings of fourth BRIA component
#
# in UKB
test = ukbl%>%filter(Comp==4) %>% filter(Approach =="BRIA") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==4) %>% filter(Approach =="BRIA") %>% mutate(B=abs(B))
# largest in CG, CST, hippocampus
#
test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]
test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")

test[grepl("CST", test$Variable),]
test3[grepl("CST", test2$Variable),]
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")

paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")
#
#
#
# in ABCD
test2 = abcdl%>%filter(Comp==4) %>% filter(Approach =="BRIA")%>% mutate(B=abs(B))
# FMIN, FMAJ, CG_hippocampus, UF identified as highly correlated
test2[grepl("FMIN", test2$Variable),]
test2[grepl("FMAJ", test2$Variable),] # not so impressive after all
test2[grepl("hippocampus", test2$Variable),] 
test2[grepl("UF", test2$Variable),] 
#

test4 = abcdl_test%>%filter(Comp==4) %>% filter(Approach =="BRIA")%>% mutate(B=abs(B))
test4[grepl("FMIN", test2$Variable),]
test4[grepl("FMAJ", test2$Variable),] # not so impressive after all
test4[grepl("hippocampus", test2$Variable),]
test4[grepl("UF", test2$Variable),]

# there is no particular one metric standing out
# how about regions?
paste("FMIN=",round(mean(test2[grepl("FMIN", test2$Variable),]$B),2),"±",round(sd(test2[grepl("FMIN", test2$Variable),]$B),2),sep = "")
paste("FMIN=",round(mean(test4[grepl("FMIN", test2$Variable),]$B),2),"±",round(sd(test4[grepl("FMIN", test2$Variable),]$B),2),sep = "")

paste("FMAJ=",round(mean(test2[grepl("FMAJ", test2$Variable),]$B),2),"±",round(sd(test2[grepl("FMAJ", test2$Variable),]$B),2),sep = "")
paste("FMAJ=",round(mean(test4[grepl("FMAJ", test2$Variable),]$B),2),"±",round(sd(test4[grepl("FMAJ", test2$Variable),]$B),2),sep = "")

paste("hippocampus=",round(mean(test2[grepl("hippocampus", test2$Variable),]$B),2),"±",round(sd(test2[grepl("hippocampus", test2$Variable),]$B),2),sep = "")
paste("hippocampus=",round(mean(test4[grepl("hippocampus", test4$Variable),]$B),2),"±",round(sd(test4[grepl("hippocampus", test2$Variable),]$B),2),sep = "")

mean(test2[grepl("UF", test2$Variable),]$B)
mean(test4[grepl("UF", test2$Variable),]$B)
paste("UF=",round(mean(test2[grepl("UF", test2$Variable),]$B),2),"±",round(sd(test2[grepl("UF", test2$Variable),]$B),2),sep = "")
paste("UF=",round(mean(test4[grepl("UF", test2$Variable),]$B),2),"±",round(sd(test4[grepl("UF", test2$Variable),]$B),2),sep = "")


########### #

# DKI ####
# in UKB
test = ukbl%>%filter(Comp==2) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==2) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)

paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
# !
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#
#
#
# in ABCD
test = ukbl%>%filter(Comp==3) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==3) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
# !
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#

# DTI ####
# in UKB
test = ukbl%>%filter(Comp==4) %>% filter(Approach =="DTI") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==4) %>% filter(Approach =="DTI") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

#
# ABCD
test = abcdl%>%filter(Comp==5) %>% filter(Approach =="DTI") %>% mutate(B=abs(B))
test3 = abcdl_test%>%filter(Comp==5) %>% filter(Approach =="DTI") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")




# SMT ####
# in UKB
test = ukbl%>%filter(Comp==5) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==5) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]
#
# ABCD
test = abcdl%>%filter(Comp==2) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test3 = abcdl_test%>%filter(Comp==2) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]







# SMTmc ####
# in UKB
test = ukbl%>%filter(Comp==4) %>% filter(Approach =="SMTmc") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==4) %>% filter(Approach =="SMTmc") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]
#
# ABCD
test = abcdl%>%filter(Comp==5) %>% filter(Approach =="SMTmc") %>% mutate(B=abs(B))
test3 = abcdl_test%>%filter(Comp==5) %>% filter(Approach =="SMTmc") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]




# WMTI ####
# in UKB
test = ukbl%>%filter(Comp==5) %>% filter(Approach =="WMTI") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==5) %>% filter(Approach =="WMTI") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]
#
# ABCD
test = abcdl%>%filter(Comp==5) %>% filter(Approach =="WMTI") %>% mutate(B=abs(B))
test3 = abcdl_test%>%filter(Comp==5) %>% filter(Approach =="WMTI") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]

#
#

# OLD ####

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


# Context: AGE ASSOCIATIONS ####
# Note: besides DKI and SMT, all other PC associations follow a similar patterns as for sex classifications
# Hence, we examine only DKI and SMT
# DKI ####
# in UKB
test = ukbl%>%filter(Comp==5) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==5) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)

paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
# !
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#
#
#
# in ABCD
test = abcdl%>%filter(Comp==5) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test3 = abcdl_test%>%filter(Comp==5) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
# !
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

# SMT ####
# in UKB
test = ukbl%>%filter(Comp==5) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==5) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]
#
# ABCD
test = abcdl%>%filter(Comp==2) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test3 = abcdl_test%>%filter(Comp==2) %>% filter(Approach =="SMT") %>% mutate(B=abs(B))
test$B = ifelse(test$B > 1, 1, test$B)
test3$B = ifelse(test3$B > 1, 1, test3$B)
paste("UF=",round(mean(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("UF", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("UF=",round(mean(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("UF", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
#!
paste("ATR=",round(mean(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("ATR", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("ATR=",round(mean(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("ATR", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CST", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CST=",round(mean(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CST", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),"±",round(sd(test[grepl("CG", gsub("CG_hippocampus","not_relevant",test$Variable)),]$B),2),sep = "")
paste("CG=",round(mean(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("CG", gsub("CG_hippocampus","not_relevant",test3$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),"±",round(sd(test[grepl("hippocampus", gsub("0","0",test$Variable)),]$B),2),sep = "")
paste("hippocampus=",round(mean(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),"±",round(sd(test3[grepl("hippocampus", gsub("0","0",test3$Variable)),]$B),2),sep = "")

test[grepl("FMAJ", test$Variable),]
test3[grepl("FMAJ", test3$Variable),]

test[grepl("FMIN", test$Variable),]
test3[grepl("FMIN", test3$Variable),]



# first dki 
test = ukbl%>%filter(Comp==1) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test3 = ukbl_test%>%filter(Comp==1) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test = abcd%>%filter(Comp==1) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))
test3 = abc_test%>%filter(Comp==1) %>% filter(Approach =="DKI") %>% mutate(B=abs(B))

