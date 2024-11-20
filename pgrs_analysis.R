# analyse new pgrs estimates
# Max Korbmacher, 20 Nov 2024
# DEFINE PATH FOR DATA AND OUTPUT
PATH = "/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/"
#
#
# PREP ####
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, fastICA, reshape2, ggpubr, 
               factoextra, dataPreparation,dplyr, 
               data.table, lme4, pheatmap, 
               gridExtra, grid, lattice,ggplotify,
               patchwork, factoextra, stats, 
               lmerTest, graphics,rlist,pwr,
               viridis,cowplot,
               update = F)
# load data
## tract
ukb = read.csv(paste(PATH,"ukb_dat.csv",sep=""))
abcd = read.csv(paste(PATH,"abcd_dat.csv",sep=""))
data=list(ukb,abcd)
## pgrs
ukb.pgrs = read.csv(paste(PATH,"UKB_PGRS.csv",sep=""))
abcd.pgrs = read.csv(paste(PATH,"ABCD_PGRS.csv",sep=""))
abcd.pgrs$eid=paste("NDAR_",abcd.pgrs$eid,sep="")# unified eid
#names(abcd.pgrs)=names(ukb.pgrs)=c("eid","ADHD","ASD","BIP","MDD","OCD","SCZ","AD")
pgrs=list(ukb.pgrs,abcd.pgrs)
## demo
ukb.demo = read.csv(paste(PATH,"ukb_anthro.csv",sep=""))
abcd.demo = read.csv(paste(PATH,"abcd_anthro.csv",sep=""))
demo=list(ukb.demo,abcd.demo)
#
# estimate principal components of white matter tract metrics
PC=list()
for(i in 1:length(data)){
  dkipc = data[[i]] %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  dtipc =  data[[i]] %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  briapc =  data[[i]] %>% select(starts_with("v_"), starts_with("micro_"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
  smtpc =  data[[i]] %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  smtmcpc =  data[[i]] %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
  wmtipc =  data[[i]] %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
  PC[[i]] = data.frame(DKI = as.numeric(dkipc$x[,1]), DTI = as.numeric(dtipc$x[,1]), BRIA = as.numeric(briapc$x[,1]), 
                  SMT = as.numeric(smtpc$x[,1]), SMTmc = as.numeric(smtmcpc$x[,1]), WMTI = as.numeric(wmtipc$x[,1]))
  PC[[i]]$eid = data[[i]]$eid
  PC[[i]] = merge(PC[[i]],demo[[i]],by="eid")
  PC[[i]] = merge(PC[[i]],pgrs[[i]],by="eid")
}
