# data prep for UKB and ABCD tract-level White Matter Assessments 
# Max Korbmacher, 16 August 2024
#
# R version: 4.2.1
#
# PREP ####
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, fastICA, reshape2, ggpubr, 
               factoextra, dataPreparation,dplyr, 
               data.table, update = F)
#
# load data
## UKB
BRIA = read.csv("/cluster/projects/p33/groups/imaging/ukbio/TBSS/stats/diffusion_data/batch_2023/models/all-participants/BRIA.csv")
DKI = read.csv("/cluster/projects/p33/groups/imaging/ukbio/TBSS/stats/diffusion_data/batch_2023/models/all-participants/DKI.csv")
DTI = read.csv("/cluster/projects/p33/groups/imaging/ukbio/TBSS/stats/diffusion_data/batch_2023/models/all-participants/DTI.csv")
SMT =  read.csv("/cluster/projects/p33/groups/imaging/ukbio/TBSS/stats/diffusion_data/batch_2023/models/all-participants/SMT.csv")
SMTmc =  read.csv("/cluster/projects/p33/groups/imaging/ukbio/TBSS/stats/diffusion_data/batch_2023/models/all-participants/SMT_mc.csv")
WMTI  =  read.csv("/cluster/projects/p33/groups/imaging/ukbio/TBSS/stats/diffusion_data/batch_2023/models/all-participants/WMTI.csv")
## ABCD
abcd = read.csv("/cluster/projects/p33/users/maxk/abcd/tabular/data/dMRI.csv")
#
# rename crucial columns
DKI = DKI %>% rename_all(~stringr::str_replace_all(.,"DKI.","")) %>% select(-c(sex,age,site_t3,site_t4,X))
DTI = DTI %>% rename_all(~stringr::str_replace_all(.,"DTI.","")) %>% select(-c(sex,age,site_t3,site_t4,X))
SMT = SMT %>% rename_all(~stringr::str_replace_all(.,"SMT.","")) %>% select(-c(sex,age,site_t3,site_t4,X))
SMTmc = SMTmc %>% rename_all(~stringr::str_replace_all(.,"SMT_mc.","")) %>% select(-c(sex,age,site_t3,site_t4,X))
WMTI = WMTI %>% rename_all(~stringr::str_replace_all(.,"WMTI.","")) %>% select(-c(sex,age,site_t3,site_t4,X))

# merge UK
ukb=merge(BRIA, DKI, by = "eid")
ukb=merge(ukb, DTI, by = "eid")
ukb=merge(ukb, SMT, by = "eid")
ukb=merge(ukb, SMTmc, by = "eid")
ukb=merge(ukb, WMTI, by = "eid")

#
# select skeleton average metrics
abcd_mean = abcd %>% select(eid, sex, scanner, age, contains("Mean")) %>% select(!contains("RSI")) %>% select(!contains("axIAD"))
ukb_mean = ukb %>% select(eid, sex, site_t3, age, contains("Mean")) %>% select(!contains("RSI")) %>% select(!contains("axIAD"))%>% 
  select(!contains("smt_fa")) %>% select(!contains("smt_trans"))
#names(ukb_mean) = c("eid", "sex", "scanner", "age", "v_intra_skeleton_Mean","v_extra_skeleton_Mean", "v_csf_skeleton_Mean", "micro")
colnames(abcd_mean) = gsub("skeleton_", "", colnames(abcd_mean))
## remove doubled diffusion approaches
colnames(abcd_mean) = gsub("WMTI_", "", colnames(abcd_mean))
colnames(abcd_mean) = gsub("Bayes_", "", colnames(abcd_mean))
colnames(ukb_mean) = gsub("SMT_mc.", "", colnames(ukb_mean))
colnames(ukb_mean) = gsub("SMT.", "", colnames(ukb_mean))
colnames(ukb_mean) = gsub("DTI.", "", colnames(ukb_mean))
colnames(ukb_mean) = gsub("DKI.", "", colnames(ukb_mean))
colnames(ukb_mean) = gsub("BRIA.", "", colnames(ukb_mean))
colnames(ukb_mean) = gsub("WMTI.", "", colnames(ukb_mean))
colnames(ukb_mean) = gsub("site", "scanner", colnames(ukb_mean))
colnames(ukb_mean) = gsub("micro_", "micro", colnames(ukb_mean))
colnames(ukb_mean) = gsub("MK", "mk", colnames(ukb_mean))
colnames(ukb_mean) = gsub("AK", "ak", colnames(ukb_mean))
colnames(ukb_mean) = gsub("RK", "rk", colnames(ukb_mean))
colnames(ukb_mean) = gsub("MD", "md", colnames(ukb_mean))
colnames(ukb_mean) = gsub("AD", "ad", colnames(ukb_mean))
colnames(ukb_mean) = gsub("RD", "rd", colnames(ukb_mean))
colnames(ukb_mean) = gsub("microadC", "microADC", colnames(ukb_mean))
colnames(ukb_mean) = gsub('\\.', "", colnames(ukb_mean))
colnames(ukb_mean) = gsub('axEad', "axEAD", colnames(ukb_mean))
colnames(ukb_mean) = gsub('radEad', "radEAD", colnames(ukb_mean))
colnames(abcd_mean) = gsub('axIAD', "axEAD", colnames(abcd_mean))
ukb_mean = rename(ukb_mean, scanner = scanner_t3)
dat_mean = rbind(ukb_mean, abcd_mean)
#
# second, select tracts for the diffusion approaches of interest
# these include BRIA, DKI, DTI, SMT, SMTmc, WMTI
abcd.dat = abcd %>% select(eid, sex, scanner, age, contains("_ATR"),contains("_CST"),
                           contains("_FMIN"),contains("_FMAJ"),contains("IFOF"),contains("_ILF"),
                           contains("_SLFR"),contains("_SLFL"),
                           contains("_UFL"),contains("_UFR"),contains("SLTFR"),
                           contains("_SLFT"), contains("CG_hippo"), contains("CCG_hippo"), contains("CGR"), 
                           contains("CGL")) %>% select(!contains("RSI")) %>% select(!contains("Cingulum")) %>% select(!contains("axIAD"))
ukb$site = ukb$site_t3
ukb.dat = ukb %>% select(eid, sex, site, age, contains("ATRL"),contains("ATRR"),contains("CST"),contains("CCG"),contains("CING"),
                         contains("FMIN"),contains("FMAJ"),contains("IFOF"),contains("ILF"),contains("SLFL"),contains("SLFR"),contains("UFL"),
                         contains("UFR"),contains("SLFTL"),contains("SLTFR")) %>% select(!contains("smt_fa"),!contains("smt_trans"),!ends_with("_SD"), !contains("Cingulum"),
                         )
# note: for abcd data, SMT FA and trans are missing.
# also: SLTF (right) tract name has a typo "SLTFR"
# standardize tract names
## remove "skeleton_" from middle of column names
colnames(abcd.dat) = gsub("skeleton_", "", colnames(abcd.dat))
## remove doubled diffusion approaches
colnames(abcd.dat) = gsub("WMTI_", "", colnames(abcd.dat))
colnames(abcd.dat) = gsub("Bayes_", "", colnames(abcd.dat))
colnames(ukb.dat) = gsub("SMT_mc.", "", colnames(ukb.dat))
colnames(ukb.dat) = gsub("SMT.", "", colnames(ukb.dat))
colnames(ukb.dat) = gsub("DTI.", "", colnames(ukb.dat))
colnames(ukb.dat) = gsub("DKI.", "", colnames(ukb.dat))
colnames(ukb.dat) = gsub("BRIA.", "", colnames(ukb.dat))
colnames(ukb.dat) = gsub("WMTI.", "", colnames(ukb.dat))
colnames(ukb.dat) = gsub("site", "scanner", colnames(ukb.dat))
colnames(ukb.dat) = gsub("micro_", "micro", colnames(ukb.dat))
colnames(ukb.dat) = gsub("MK", "mk", colnames(ukb.dat))
colnames(ukb.dat) = gsub("AK", "ak", colnames(ukb.dat))
colnames(ukb.dat) = gsub("RK", "rk", colnames(ukb.dat))
colnames(ukb.dat) = gsub("MD", "md", colnames(ukb.dat))
colnames(ukb.dat) = gsub("AD", "ad", colnames(ukb.dat))
colnames(ukb.dat) = gsub("RD", "rd", colnames(ukb.dat))
colnames(ukb.dat) = gsub("microadC", "microADC", colnames(ukb.dat))
colnames(ukb.dat) = gsub("CING.1", "CG_hippocampus", colnames(ukb.dat))
colnames(ukb.dat) = gsub("CING", "CG_hippocampus", colnames(ukb.dat))
colnames(ukb.dat) = gsub("Cingulumcingulategyrus", "CG", colnames(ukb.dat))
colnames(ukb.dat) = gsub('\\.', "", colnames(ukb.dat))
colnames(ukb.dat) = gsub('axEad', "axEAD", colnames(ukb.dat))
colnames(ukb.dat) = gsub('radEad', "radEAD", colnames(ukb.dat))
colnames(abcd.dat) = gsub('axIAD', "axEAD", colnames(abcd.dat))
colnames(abcd.dat) = gsub('CCG_hippocampus', "CG_hippocampus", colnames(abcd.dat))
colnames(abcd.dat) = gsub('CG_hippocampus', "CG_hippocampus", colnames(abcd.dat))
#
# final check whether there are any unexpected mismatches
ukb.dat[!colnames(ukb.dat) %in% colnames(abcd.dat)] %>% names() # yes, SMT_fa and trans are missing, as outlined above
abcd.dat[!colnames(abcd.dat) %in% colnames(ukb.dat)] %>% names()
#
# remove outliers 5 standard deviations from mean
# removeOutliers <- function(data, column_name) {
#   SD = sd(data[[column_name]],na.rm=T)
#   MEAN = mean(data[[column_name]],na.rm=T)
#   lower_bound <- MEAN-5*SD
#   upper_bound <- MEAN+5*SD
#   return(data[data[[column_name]] >= lower_bound & as.numeric(data[[column_name]]) <= upper_bound, ])
# }
ukb.dat[5:ncol(ukb.dat)]=data.frame(lapply(ukb.dat[5:ncol(ukb.dat)], as.numeric)) # make all numeric columns numeric
ukb.dat = (na.omit(ukb.dat, cols = names(ukb.dat[5:ncol(ukb.dat)])))
# for (i in 5:ncol(ukb.dat)){
#   ukb.dat = ukb.dat[ukb.dat$eid %in% (removeOutliers(ukb.dat, names(ukb.dat[5:ncol(ukb.dat)])[i]))$eid,]
# }
ukb.dat = remove_sd_outlier(ukb.dat, cols = "auto", n_sigmas = 5, verbose = TRUE)
print(paste("Excluded UKB participants:",nrow(ukb)-nrow(ukb.dat),sep = ""))
before = nrow(abcd.dat)
abcd.dat = remove_sd_outlier(abcd.dat, cols = "auto", n_sigmas = 5, verbose = TRUE)
# for (i in 5:ncol(abcd.dat[5:ncol(abcd.dat)])){
#   abcd.dat = removeOutliers(abcd.dat, names(ukb.dat[5:ncol(abcd.dat)])[i])
# }

print(paste("Excluded ABCD participants:",before-nrow(abcd.dat),sep = ""))
# check which participants were excluded
ukb[!(ukb$eid) %in% (ukb.dat$eid),]$eid

# match data frames
ukb.dat = data.frame(ukb.dat)
abcd.dat = data.frame(abcd.dat)
ukb.dat = (ukb.dat[c(colnames(ukb.dat) %in% colnames(abcd.dat))])
dat = rbind(abcd.dat,ukb.dat)
dat$data = c(replicate(nrow(abcd.dat), "ABCD"),replicate(nrow(ukb.dat), "UKB"))
write.csv(dat,"/cluster/projects/p33/users/maxk/wm_ica/dat.csv", row.names = F)
write.csv(ukb.dat,"/cluster/projects/p33/users/maxk/wm_ica/ukb_dat.csv", row.names = F)
write.csv(abcd.dat,"/cluster/projects/p33/users/maxk/wm_ica/abcd_dat.csv", row.names = F)
#
# do the same for the skeleton mean values
dat_mean = dat_mean[dat_mean$eid %in% dat$eid,]
abcd_mean = abcd_mean[abcd_mean$eid %in% abcd.dat$eid,]
ukb_mean = ukb_mean[ukb_mean$eid %in% ukb.dat$eid,]
write.csv(dat_mean,"/cluster/projects/p33/users/maxk/wm_ica/dat_mean.csv", row.names = F)
write.csv(ukb_mean,"/cluster/projects/p33/users/maxk/wm_ica/ukb_mean.csv", row.names = F)
write.csv(abcd_mean,"/cluster/projects/p33/users/maxk/wm_ica/abcd_mean.csv", row.names = F)
#