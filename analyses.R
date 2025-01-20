# analyses to examine the extend to which DTI explains white matter microstructure
# Max Korbmacher, 27 August 2024
#
# R version: 4.2.1
#
#
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
#
# DEFINE PATH FOR DATA AND OUTPUT
PATH = "/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/WMM_Approaches_comparison/"
# create subdir called Figures [prompts warning and skips, if folder exists already]
dir.create(file.path(PATH, "Figures"))
#
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
               viridis,cowplot,neuroCombat,caret,
               PCAtest,pscl,rstatix,
               update = F)
#
# load data
## skeleton averages
#dat_mean = read.csv(paste(PATH,"dat_mean.csv",sep=""))
ukb_mean = read.csv(paste(PATH,"ukb_mean.csv",sep=""))
abcd_mean = read.csv(paste(PATH,"abcd_mean.csv",sep=""))
## tract averages
#dat = read.csv(paste(PATH,"dat.csv",sep=""))
ukb = read.csv(paste(PATH,"ukb_dat.csv",sep=""))
abcd = read.csv(paste(PATH,"abcd_dat.csv",sep=""))
#
# volumetrics
vols = read.csv(paste(PATH,"vols.csv",sep=""))
# transform abcd age from months to years
abcd_mean$age = abcd_mean$age/12
abcd$age = abcd$age/12
ukb$sex = ukb_mean$sex=factor(ukb_mean$sex)
levels(ukb$sex) = levels(ukb_mean$sex) = c("female", "male")
abcd$sex = abcd_mean$sex = factor(abcd$sex)
dat_mean = rbind(ukb_mean,abcd_mean)
dat = rbind(ukb,abcd)
#
#
#
# 00) Harmonisation: COMBAT ####
# we correct for scanner effects and keep variance patterns in tact for sex and age
## starting with mean values
covars = dat_mean %>% dplyr::select(eid,sex,scanner,age)
dat_mean = neuroCombat(t(dat_mean%>%dplyr::select(-eid,-sex,-scanner,-age)),batch=dat_mean$scanner,mod=model.matrix(~covars$age+covars$sex))
dat_mean = data.frame(t(dat_mean$dat.combat))
dat_mean = cbind(covars,dat_mean)
# specify data type
dat_mean$data = c(replicate(nrow(ukb_mean),"UKB"),replicate(nrow(abcd_mean),"ABCD"))
#
## continuing with tract values
dat = neuroCombat(t(dat%>%dplyr::select(-eid,-sex,-scanner,-age)),batch=dat$scanner,mod=model.matrix(~covars$age+covars$sex))
dat = data.frame(t(dat$dat.combat))
dat = cbind(covars,dat)
# specify data type
dat$data = c(replicate(nrow(ukb_mean),"UKB"),replicate(nrow(abcd_mean),"ABCD"))
#
# merge volumetric data for later correction of vol effects
dat_mean = left_join(dat_mean,vols,by=c("eid","data"))
dat = left_join(dat,vols,by=c("eid","data"))
#
# sep train-test
my.ids <- createDataPartition(dat_mean$data, p = 0.2) #
dat_mean_val = dat_mean[as.numeric(my.ids[[1]]), ]
dat_val = dat[as.numeric(my.ids[[1]]), ]
dat_mean = dat_mean[!dat_mean$eid %in% dat_mean_val$eid,]
dat = dat[!dat$eid %in% dat_val$eid,]
# get unique sets for training data
abcd = dat %>% filter(data == "ABCD")
ukb = dat %>% filter(data == "UKB")
abcd_mean = dat_mean %>% filter(data == "ABCD")
ukb_mean = dat_mean %>% filter(data == "UKB")
# validation sets
abcd_val = dat_val %>% filter(data == "ABCD")
ukb_val = dat_val %>% filter(data == "UKB")
abcd_mean_val = dat_mean_val %>% filter(data == "ABCD")
ukb_mean_val = dat_mean_val %>% filter(data == "UKB")
#
# write the harmonized mri training frames
write.csv(x = ukb,file = paste(PATH,"ukb_dat1.csv",sep=""))
write.csv(x = abcd,file = paste(PATH,"abcd_dat1.csv",sep=""))
write.csv(x = ukb_mean,file = paste(PATH,"ukb_mean1.csv",sep=""))
write.csv(x = abcd_mean,file = paste(PATH,"abcd_mean1.csv",sep=""))
# as well as test frames
write.csv(x = ukb_val,file = paste(PATH,"ukb_dat_val.csv",sep=""))
write.csv(x = abcd_val,file = paste(PATH,"abcd_dat_val.csv",sep=""))
write.csv(x = ukb_mean_val,file = paste(PATH,"ukb_mean_val.csv",sep=""))
write.csv(x = abcd_mean_val,file = paste(PATH,"abcd_mean_val.csv",sep=""))
#
# 0) Power #### 
# quick and dirty power analysis for univariate associations with the lowest power
##ABCD data only
pwr_out = pwr.f2.test(u = 4,v = 3691,sig.level = 0.05/708,power = .95)
print(paste("lowest observable f2 and R2 = ",(pwr_out$f2),"(ABCD)",sep = ""))
print(paste("lowest observable r = ",sqrt(pwr_out$f2),"(ABCD)",sep = ""))
## UKB data only
pwr_out = pwr.f2.test(u = 4,v = 29653,sig.level = 0.05/708,power = .95)
print(paste("lowest observable f2 and R2 = ",(pwr_out$f2),"(UKB)",sep = ""))
print(paste("lowest observable r = ",sqrt(pwr_out$f2),"(UKB)",sep = ""))
# ## all data
# pwr_out = pwr.f2.test(u = 4,v = 47371,sig.level = 0.05/708,power = .95)
# print(paste("lowest observable f2 and R2 = ",(pwr_out$f2),"(all data)",sep = ""))
# print(paste("lowest observable r = ",sqrt(pwr_out$f2),"(all data)",sep = ""))
# u = 4 for metric, sex, site, age
# v = 5,147 (ABCD datasets) - 5
# 1) Check skeleton averages ####
# 
# make list of dependent variables
dependent = dat_mean %>% select(-c(md_Mean,rd_Mean, ad_Mean, FA_Mean, eid, sex, age, scanner, data)) %>% names
# standardize independent and dependent variables to obtain standardized (and comparable coefficients)
stand = function(dataframe){
  newdataframe = dataframe %>% select(-c(eid, sex, age, scanner, data)) %>% scale %>% data.frame()
  newdataframe = cbind(newdataframe, dataframe %>% select(c(eid, sex, age, scanner)))
  return(newdataframe)
}
# create list of dfs
mean_list = list(stand(ukb_mean), stand(abcd_mean), stand(dat_mean))
# 1.1) Correlations controlled for age, sex, scanner, and cortical volume ####
# create a list for the results
# res = list()
# for (df in 1:length(mean_list)){
#   AD = AD.SE = FA = FA.SE = RD = RD.SE = MD = MD.SE = c()
#   for (i in 1:length(dependent)){
#     f1 = formula(paste(dependent[i],"~ad_Mean+age+sex+(scanner)",sep="")) # regular linear models used due to singular fit 
#     f2 = formula(paste(dependent[i],"~FA_Mean+age+sex+(scanner)",sep="")) # otherwise lmer with site as random intercept
#     f3 = formula(paste(dependent[i],"~rd_Mean+age+sex+(scanner)",sep=""))
#     f4 = formula(paste(dependent[i],"~md_Mean+age+sex+(scanner)",sep=""))
#     m1 = lm(f1, data=mean_list[[df]])
#     m2 = lm(f2, data=mean_list[[df]])
#     m3 = lm(f3, data=mean_list[[df]])
#     m4 = lm(f4, data=mean_list[[df]])
#     AD[i] = summary(m1)$coefficients[2]
#     AD.SE[i] = summary(m1)$coefficients[2,2]
#     FA[i] = summary(m2)$coefficients[2]
#     FA.SE[i] = summary(m2)$coefficients[2,2]    
#     RD[i] = summary(m3)$coefficients[2]
#     RD.SE[i] = summary(m3)$coefficients[2,2]
#     MD[i] = summary(m4)$coefficients[2]
#     MD.SE[i] = summary(m4)$coefficients[2,2]
#     }
#   res[[df]] = data.frame(dependent, AD, AD.SE, FA, FA.SE, RD, RD.SE, MD, MD.SE)
# }
deps = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
# for (i in 1:length(res)){
#   res[[i]]$dependent = deps
#   res[[i]] = melt(res[[i]]%>%select(!contains("SE"))) %>% rename(Beta = value)
# }
# p1 = ggplot(res[[1]], aes(x=variable, y=dependent)) + 
#   geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
#   scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
#                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
#   theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
# p2 = ggplot(res[[2]], aes(x=variable, y=dependent)) + 
#   geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
#   scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
#                        breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
#   theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
# p3 = ggplot(res[[3]], aes(x=variable, y=dependent)) + 
#   geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
#   scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
#                        breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
#   theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
# plot1 = ggarrange(p1,p2,p3,nrow=1,labels = c("a", "b", "c"), common.legend = T, legend = "right") #c("UKB", "ABCD", "Both")
# ggsave(paste(PATH,"Figures/plot1.pdf",sep=""),plot1,width=9,height=7)
# 1.2) Crude correlations ####
# create a list for the results
# res2 = list()
# for (df in 1:length(mean_list)){
#   AD = AD.SE = FA = FA.SE = RD = RD.SE = MD = MD.SE = c()
#   for (i in 1:length(dependent)){
#     f1 = formula(paste(dependent[i],"~ad_Mean",sep=""))
#     f2 = formula(paste(dependent[i],"~FA_Mean",sep=""))
#     f3 = formula(paste(dependent[i],"~rd_Mean",sep=""))
#     f4 = formula(paste(dependent[i],"~md_Mean",sep=""))
#     m1 = lm(f1, data=mean_list[[df]])
#     m2 = lm(f2, data=mean_list[[df]])
#     m3 = lm(f3, data=mean_list[[df]])
#     m4 = lm(f4, data=mean_list[[df]])
#     AD[i] = summary(m1)$coefficients[2] 
#     AD.SE[i] = summary(m1)$coefficients[2,2]
#     FA[i] = summary(m2)$coefficients[2]
#     FA.SE[i] = summary(m2)$coefficients[2,2]    
#     RD[i] = summary(m3)$coefficients[2]
#     RD.SE[i] = summary(m3)$coefficients[2,2]
#     MD[i] = summary(m4)$coefficients[2]
#     MD.SE[i] = summary(m4)$coefficients[2,2]
#   }
#   res2[[df]] = data.frame(dependent, AD, AD.SE, FA, FA.SE, RD, RD.SE, MD, MD.SE)
# }
deps = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
# for (i in 1:length(res2)){
#   res2[[i]]$dependent = deps
#   res2[[i]] = melt(res2[[i]]%>%select(!contains("SE"))) %>% rename(Beta = value)
# }
# p1 = ggplot(res2[[1]], aes(x=variable, y=dependent)) + 
#   geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
#   scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
#                        breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
#   theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
# p2 = ggplot(res2[[2]], aes(x=variable, y=dependent)) + 
#   geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
#   scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
#                        breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
#   theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
# p3 = ggplot(res2[[3]], aes(x=variable, y=dependent)) + 
#   geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
#   scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
#                        breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
#   theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
# plot2 = ggarrange(p1,p2,p3,nrow=1,labels = c("a", "b", "c"), common.legend = T, legend = "right") #c("UKB", "ABCD", "Both")
# ggsave(paste(PATH,"Figures/plot2.pdf",sep=""),plot2,width=9,height=7)
#
#
#
#
#
# Summary stats for the skeleton associations
# Note, here we assess only the corrected associations, not the crude ones.
# res[[1]] %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
# a = res[[1]] %>% group_by(variable) %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
# res[[2]] %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
# b = res[[2]] %>% group_by(variable) %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
# res[[3]] %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
# c = res[[3]] %>% group_by(variable) %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
# cbind(a,b,c)
# rm(a,b,c)
#
# 2) Crude correlations between all skeleton averages ####
deps2 = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK","DTI-FA", "DTI-MD","DTI-RD","DTI-AD", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
cormat = function(df){
  cor1 = df%>%select(-c(eid,age,sex,scanner,CortexVol))
  #roworder = names(cor1)
  names(cor1) = deps2
  cor1 = cor(cor1)
  return(cor1)
}
a = as.ggplot(pheatmap(cormat(mean_list[[1]]),cluster_rows = T, cluster_cols = T))
b = as.ggplot(pheatmap(cormat(mean_list[[2]]),cluster_rows = T, cluster_cols = T))
c = as.ggplot(pheatmap(cormat(mean_list[[3]]),cluster_rows = T, cluster_cols = T))
# d = as.ggplot(pheatmap(cormat(mean_list[[1]]),cluster_rows = F, cluster_cols = T))
# e = as.ggplot(pheatmap(cormat(mean_list[[2]]),cluster_rows = F, cluster_cols = T))
# f = as.ggplot(pheatmap(cormat(mean_list[[3]]),cluster_rows = F, cluster_cols = T))
plot3 = a+b+c#+d+e+f
# note: top row = unclustered ukb,abcd,both; bottom row = clustered
ggsave(paste(PATH,"Figures/skeleton_correlations.pdf",sep=""),plot3,width=20,height=7)
rm(mean_list,plot3)
# 3) Principal component analysis (PCA) ####
# 3.1) Estimation and visualization of the PCs ####
# first, test:
#   (1) the hypothesis that there is more correlational structure among the observed variables than expected by random chance, 
#   (2) the statistical significance of each PC, and 
#   (3) the contribution of each observed variable to each significant PC. (although this last point is not used)
# This was done one-by-one:
#
# UNTICK FOR CHECKS
#
## UKB
result = PCAtest(ukb %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")),100,100,0.05)
# result = PCAtest(ukb %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")),100,100,0.05)
# result = PCAtest(ukb %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")),10,10,0.05)
# result = PCAtest(ukb %>% select(starts_with("smt_md"), starts_with("smt_long")),100,100,0.05)
# result = PCAtest(ukb %>%select(starts_with("smt_mc")),100,100,0.05)
# result = PCAtest(ukb %>%select(starts_with("axEAD"), starts_with("radEAD")),100,100,0.05)
## ABCD
result = PCAtest(abcd %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")),100,100,0.05)
# result = PCAtest(abcd %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")),100,100,0.05)
# result = PCAtest(abcd %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")),10,10,0.05)
# result = PCAtest(abcd %>% select(starts_with("smt_md"), starts_with("smt_long")),100,100,0.05) # only 3 sig
# result = PCAtest(abcd %>%select(starts_with("smt_mc")),100,100,0.05)
# result = PCAtest(abcd %>%select(starts_with("axEAD"), starts_with("radEAD")),100,100,0.05)
#
#
# estimate principal components and visualise the first 10 components
pcplots = function(data){
dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
a = fviz_eig(dkipc, main = "DKI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
b = fviz_eig(dtipc, main = "DTI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
c = fviz_eig(briapc, main = "BRIA", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
d = fviz_eig(smtpc, main = "SMT", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
e = fviz_eig(smtmcpc, main = "SMTmc", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
f = fviz_eig(wmtipc, main = "WMTI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
plot4 = ggarrange(a,b,c,d,e,f)
return(plot4)
}
plot4.1 = pcplots(dat%>%select(-c(eid,age,sex,scanner)))
plot4.2 = pcplots(ukb%>%select(-c(eid,age,sex,scanner)))
plot4.3 = pcplots(abcd%>%select(-c(eid,age,sex,scanner)))
ggsave(paste(PATH,"Figures/scree.pdf",sep=""),plot4.1,width=14,height=10)
ggsave(paste(PATH,"Figures/UKB_scree.pdf",sep=""),plot4.2,width=14,height=10)
ggsave(paste(PATH,"Figures/ABCD_scree.pdf",sep=""),plot4.3,width=14,height=10)
rm(plot4.1,plot4.2,plot4.3)
#
# 3.2) Associations of the PCs with each other ####
# extract the first component
pccor = function(data){
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
  PC = data.frame(DKI = as.numeric(dkipc$x[,1]), DTI = as.numeric(dtipc$x[,1]), BRIA = as.numeric(briapc$x[,1]), 
                  SMT = as.numeric(smtpc$x[,1]), SMTmc = as.numeric(smtmcpc$x[,1]), WMTI = as.numeric(wmtipc$x[,1]))
  plot5 = as.ggplot(pheatmap(cor(PC),cluster_rows = F, cluster_cols = F, display_numbers = round(cor(PC),2)))
  return(plot5)
}
a = pccor(dat%>%select(-c(eid,age,sex,scanner)))
b = pccor(ukb%>%select(-c(eid,age,sex,scanner)))
c = pccor(abcd%>%select(-c(eid,age,sex,scanner)))
plot5 = ggarrange(b,c,a,ncol=1,common.legend = T)
ggsave(paste(PATH,"Figures/plot5.pdf",sep=""),plot5,width=4,height=10)
#
# estimate corrected correlations, standard errors and p-values
cors = function(data){
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
  PC = data.frame(DKI = as.numeric(dkipc$x[,1]), DTI = as.numeric(dtipc$x[,1]), BRIA = as.numeric(briapc$x[,1]), 
                  SMT = as.numeric(smtpc$x[,1]), SMTmc = as.numeric(smtmcpc$x[,1]), WMTI = as.numeric(wmtipc$x[,1]))
  PC = data.frame(scale(PC))
  PC$age = as.numeric(scale(data$age))
  PC$sex = factor(data$sex)
  PC$scanner = factor(data$scanner)
  WMMapproach = c("DKI", "BRIA", "SMT", "SMTmc", "WMTI")
  B = SE = p = c()
  for (i in 1:length(WMMapproach)){
    #f1 = formula(paste(WMMapproach[i],"~DTI+age+sex+(1|scanner)",sep="")) # singulary issues when 
    f1 = formula(paste(WMMapproach[i],"~DTI+age+sex+(scanner)",sep=""))
    mod = lm(f1,PC)
    B[i] = summary(mod)$coefficients[2,1]
    SE[i] = summary(mod)$coefficients[2,2]
    p[i] = summary(mod)$coefficients[2,4]
  }
  export = data.frame(B,SE,p)
  names(export) = c("B", "SE","p")
  rownames(export) = c("DKI", "BRIA", "SMT", "SMTmc", "WMTI")
  return(export)
}
tmp = rbind(cors(ukb), cors(abcd),cors(dat))
tmp$Approach = row.names(cors(ukb))
tmp$Data = c(replicate(nrow(tmp)/3,"UKB"),replicate(nrow(tmp)/3,"ABCD"),replicate(nrow(tmp)/3,"Both"))
#base::order(tmp$Data)=1:length((tmp$Data))
tmp$B = abs(tmp$B)
plot6 = ggplot(tmp, aes(x=Approach, y=B, color = Data)) +
  geom_errorbar(aes(ymin=B-SE, ymax=B+SE,group=Data), width=.2,position=position_dodge(.6)) +
  geom_point(position = position_dodge(width = 0.6))+
  scale_color_manual(values=c("#880700", "#E69F00", "#3a81b5"))+ ##999999','#E69F00', '#56B4E9'
  #geom_bar( aes(x=row.names(tmp), y=B), stat="identity", fill="skyblue", alpha=0.7) +
  #geom_pointrange( aes(x=Approach, y=B, ymin=B-SE, ymax=B+SE, fill = Data), colour="black", alpha=0.06, size=0.6)+
  coord_flip() + theme_bw() + xlab("Diffusion Approach") + ylab("Adjusted absolute standardized beta coefficients") + theme(text = element_text(size=15))
ggsave(paste(PATH,"Figures/plot6.pdf",sep=""),plot6,width=8,height=7)
#
# 3.3) Associations of the PCs with global/skeleton metrics ####
# estimate corrected correlations, standard errors and p-values
cors = function(data){ # mean_data ... data argument is for the data frame containing tract-level data, mean_data for skeleton-level data
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = F, center = F)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = F, center = F)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = F, center = F)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = F, center = F)
  smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = F, center = F)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD"),starts_with("awf")) %>% prcomp(scale. = F, center = F)
  pcdatalist = list(dkipc,dtipc,briapc,smtpc,smtmcpc,wmtipc)
  PC=list()
  for (i in 1:5){
    PC[[i]] = data.frame(DKI = as.numeric(dkipc$x[,i]), DTI = as.numeric(dtipc$x[,i]), BRIA = as.numeric(briapc$x[,i]), 
                    SMT = as.numeric(smtpc$x[,i]), SMTmc = as.numeric(smtmcpc$x[,i]), WMTI = as.numeric(wmtipc$x[,i]))
    PC[[i]] = data.frame(scale(PC[[i]]))
    PC[[i]]$eid = data$eid
    PC[[i]] = merge(data,PC[[i]])
  }
  # for (n in 1:length(predictor_vars)){
  #   mean_data[predictor_vars[n]] = c(scale(mean_data[predictor_vars[n]]))
  # }
  WMMapproach = c("DKI", "DTI","BRIA", "SMT", "SMTmc", "WMTI")
  ex=list()
  for (i in 1:length(WMMapproach)){
    predictor_vars = rownames(pcdatalist[[i]]$rotation)
    export = list()
    for (j in 1:length(predictor_vars)){
      Comp=c()
      B = SE = p = c()
      for (l in 1:length(PC)){
        f1 = formula(paste(WMMapproach[i],"~",predictor_vars[j],"+age+sex+scanner+CortexVol",sep=""))
        mod = lm(f1,PC[[l]])
        B[l] = summary(mod)$coefficients[2,1]
        SE[l] = summary(mod)$coefficients[2,2]
        p[l] = summary(mod)$coefficients[2,4]
        Comp[l] = l
      }
      export[[j]] = data.frame(B,SE,p,Comp)
      export[[j]]$Variable = predictor_vars[j]
      export[[j]]$Approach = WMMapproach[i]
      names(export[[j]]) = c("B", "SE","p","Comp","Variable","Approach")
    }
    ex[[i]] = rbindlist(export)
    #rownames(export[[j]]) = WMMapproach
  }
  names(ex)=WMMapproach
  return(ex)
}
ukblist = cors(ukb %>% mutate(across(where(is.numeric), scale)))
abcdlist = cors(abcd %>% mutate(across(where(is.numeric), scale)))
datlist = cors(dat %>% mutate(across(where(is.numeric), scale)))
corrected_cors_ukb = data.frame(rbindlist(ukblist))
corrected_cors_abcd = data.frame(rbindlist(abcdlist))
corrected_cors_both = data.frame(rbindlist(datlist))
write.csv(x=corrected_cors_ukb, file=paste(PATH,"corrected_loadings_ukb.csv",sep=""))
write.csv(x=corrected_cors_abcd, file=paste(PATH,"corrected_loadings_abcd.csv",sep=""))
write.csv(x=corrected_cors_both, file=paste(PATH,"corrected_loadings_both.csv",sep=""))
#
#
#
# 4) Age-associations ####
# 4.1) Lifespan plots ####
# start with visualising age-associations (skeleton level)
preds=names(dat_mean[5:30]) # create vector of variables to be predicted
lifespan=function(predictor_var, ylab){
  ggplot(stand(dat_mean), aes_string(x="age", y=predictor_var)) + geom_point(alpha = 0.05) + #geom_density_2d(col="black")+
    stat_smooth(method = "gam",formula = y ~ s(x, k = 4), col = "blue", se = TRUE) +theme_bw()+
    xlab("Age")+ylab(ylab)
}
plotlist=list()
for (i in 1:length(preds)){
  plotlist[[i]]=lifespan(preds[i],deps2[i])
}
plot9 = ggarrange(plotlist = plotlist)
ggsave(paste(PATH,"Figures/plot9.jpg",sep=""),plot9,width=15,height=10)
# 4.2) Linear models for comparison of beta coefficients ASSESSING AGE ASSOCIATIONS ####
# prep models for skeleton-level associations
res = list()
mean_list=list(ukb_mean,abcd_mean)
Beta=SE=c()
for (df in 1:length(mean_list)){
  for (i in 1:length(preds)){
    tmp_dat = stand(mean_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+age+sex+scanner+CortexVol",sep="")) # regular linear models used due to singular fit 
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = summary(m1)$coefficients[2]
    SE[i] = summary(m1)$coefficients[2,2] # not necessary for the plotting
  }
  res[[df]] = data.frame(preds, Beta, SE)
}
res[[1]]$Data = "UKB"
res[[2]]$Data = "ABCD"
plot_dat = list.rbind(res) # merge list/data frames
#
plot_dat$Metric=deps2 # add plotable metric names
# split data up into 3 parts for better plotting
BRIA.df=dplyr::filter(plot_dat, grepl("BRIA",Metric))
DTI.DKI.df = rbind(plot_dat %>% filter(grepl("DKI",Metric)),plot_dat %>% filter(grepl("DTI",Metric)))
SMT.WMTI.df = rbind(plot_dat %>% filter(grepl("SMT",Metric)),plot_dat %>% filter(grepl("WMTI",Metric)))
res=list(BRIA.df,DTI.DKI.df,SMT.WMTI.df)
for (i in 1:length(res)){
  res[[i]] = melt(res[[i]]%>%select(-SE,-preds))%>% dplyr::rename(Beta = value)
}
p1 = ggplot(res[[1]], aes(x=Data, y=Metric)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2 = ggplot(res[[2]], aes(x=Data, y=Metric)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p3 = ggplot(res[[3]], aes(x=Data, y=Metric)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
#plot10 = ggarrange(p1,p2,p3,nrow=1, common.legend = T, legend = "bottom") #c("UKB", "ABCD", "Both")
plot10=plot_grid(p1,p2,p3,nrow=1, align="v")
ggsave(paste(PATH,"Figures/age_plot_skeleton.pdf",sep=""),plot10,width=11,height=5)
#
#
#
# 5.) SEX CLASSIFICATIONS ####
pcest = function(data){
  dkipc = data.frame((data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T))$x[,1:5])
  dtipc =  data.frame((data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T))$x[,1:5])
  briapc =  data.frame((data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T))$x[,1:5])
  smtpc =  data.frame((data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T))$x[,1:5])
  smtmcpc =  data.frame((data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T))$x[,1:5])
  wmtipc =  data.frame((data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T))$x[,1:5])
  d = list(briapc,dkipc,dtipc,smtpc,smtmcpc,wmtipc)
  return(d)
}
# first five components for each: BRIA, DKI, DTI, SMT, SMTmc, WMTI
ukb_pc = pcest(ukb%>%select(-c(eid,age,sex,scanner,CortexVol))) # re-estimate PCs
abcd_pc = pcest(abcd%>%select(-c(eid,age,sex,scanner,CortexVol)))
#
#
#
###### SCALE THE TEST DATA BASED ON THE TRAINING DATA
#
# scale the test data based on the training data
ukb_scaled = scale(ukb%>%select(-c(eid,age,sex,scanner,CortexVol,data)))
ukb_val_scaled = data.frame(scale(ukb_val%>%select(-c(eid,age,sex,scanner,CortexVol,data)), center=attr(ukb_scaled, "scaled:center"), scale=attr(ukb_scaled, "scaled:scale")))
abcd_scaled = scale(abcd%>%select(-c(eid,age,sex,scanner,CortexVol,data)))
abcd_val_scaled = data.frame(scale(abcd_val%>%select(-c(eid,age,sex,scanner,CortexVol,data)), center=attr(abcd_scaled, "scaled:center"), scale=attr(abcd_scaled, "scaled:scale")))
#
# get the list into the same format as for the training data
sel = function(data){ 
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) 
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) 
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) 
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long"))
  smtmcpc =  data %>% select(starts_with("smt_mc")) 
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) 
  d = list(briapc,dkipc,dtipc,smtpc,smtmcpc,wmtipc)
  return(d)
}
abcd_pc_val = sel(abcd_val_scaled)
ukb_pc_val = sel(ukb_val_scaled)
# get the respective PCA models from the training data
pcmf = function(data){ 
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk"))  %>% prcomp(scale. = T, center = T)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa"))  %>% prcomp(scale. = T, center = T)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax"))  %>% prcomp(scale. = T, center = T)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  smtmcpc =  data %>% select(starts_with("smt_mc"))  %>% prcomp(scale. = T, center = T)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD"))  %>% prcomp(scale. = T, center = T)
  d = list(briapc,dkipc,dtipc,smtpc,smtmcpc,wmtipc)
  return(d)
}
abcd_pcm = pcmf(abcd)
ukb_pcm = pcmf(ukb)
###### THEN PREDICT IN THE TEST SET (get loadings)
# these individual level component scores or loadings can then be used as to predict things, such as sex and age
#
# Note that the predict.pcrcomp() or simply predict() call has the scaling based on the original model included
# Hence, the scaling done above resulting in the lists abcd_pc_val, and ukb_pc_val is not necessary
#
abcd_test_loadings = ukb_test_loadings = list()
for (i in 1:length(abcd_pcm)){
  abcd_test_loadings[[i]] = data.frame(predict(abcd_pcm[[i]], newdata = abcd_val))[1:5]
  ukb_test_loadings[[i]] = data.frame(predict(ukb_pcm[[i]], newdata = ukb_val))[1:5]
}
#
# We have already estimated the corrected correlations between single variables and PCs for the training set.
# Now, we do it for the test set to validate the associations
cors = function(data, data_pc_val, data_test_loadings){ # mean_data ... data argument is for the data frame containing tract-level data, mean_data for skeleton-level data
  PC=list()
  for (i in 1:5){
    PC[[i]] = cbind(data %>% select(age,sex,scanner,CortexVol), data_test_loadings[[i]],data_pc_val)
  }
  WMMapproach = c("BRIA","DKI", "DTI", "SMT", "SMTmc", "WMTI")
  ex=list()
  for (i in 1:length(WMMapproach)){
    predictor_vars = names(data_pc_val[[i]])
    export = list()
    for (j in 1:length(predictor_vars)){
      Comp=c()
      B = SE = p = c()
      for (l in 1:length(PC)){
        f1 = formula(paste("PC",l,"~",predictor_vars[j],"+age+sex+scanner+CortexVol",sep=""))
        mod = lm(f1,PC[[l]])
        B[l] = summary(mod)$coefficients[2,1]
        SE[l] = summary(mod)$coefficients[2,2]
        p[l] = summary(mod)$coefficients[2,4]
        Comp[l] = l
      }
      export[[j]] = data.frame(B,SE,p,Comp)
      export[[j]]$Variable = predictor_vars[j]
      export[[j]]$Approach = WMMapproach[i]
      names(export[[j]]) = c("B", "SE","p","Comp","Variable","Approach")
    }
    ex[[i]] = rbindlist(export)
    #rownames(export[[j]]) = WMMapproach
  }
  names(ex)=WMMapproach
  return(ex)
}
abcd_val_list = cors(abcd_val, abcd_pc_val, abcd_test_loadings) # get test result lists
ukb_val_list = cors(abcd_val, abcd_pc_val, abcd_test_loadings)
ukb_val_list = data.frame(rbindlist(ukb_val_list)) # put them into data frame format
abcd_val_list = data.frame(rbindlist(abcd_val_list))
write.csv(x=ukb_val_list, file=paste(PATH,"corrected_loadings_ukb_TEST.csv",sep="")) # save
write.csv(x=abcd_val_list, file=paste(PATH,"corrected_loadings_abcd_TEST.csv",sep=""))
#
#
### NOW, FINALLY CLASSIFY
#
# we focus first model performance expressed as accuracy
# Other options which can be considered: LHR, pseudo R2, variable importance
#
out=list()
mods = c("BRIA","DKI","DTI","SMT","SMTmc","WMTI") # dMRI models
sex_pred = function(orig_train,orig_test, train_list, test_list, datasource){
  for (i in 1:length(train_list)){
    df = cbind(orig_train%>%select(c(eid,age,sex,scanner,CortexVol)),train_list[[i]])
    m = glm(sex~PC1+PC2+PC3+PC4+PC5+CortexVol+age+scanner, family = "binomial", data = df)
    testing = cbind(orig_test%>%select(c(eid,age,sex,scanner,CortexVol)),test_list[[i]])
    pred = predict(m, newdata=testing,type="response")
    pred = factor(round(pred))
    levels(pred) = c("female","male")
    accuracy = table(pred, testing[,"sex"])
    Test = sum(diag(accuracy))/sum(accuracy)
    #confusionMatrix(data=(pred), testing$sex)
    pred = predict(m, newdata=df,type="response")
    pred = factor(round(pred))
    levels(pred) = c("female","male")
    accuracy = table(pred, df$sex)
    Train = sum(diag(accuracy))/sum(accuracy)
    #confusionMatrix(data=(pred), df$sex)
    out[[i]] = data.frame(Model = mods[i], round(Train,4)*100, round(Test,4)*100)
    names(out[[i]]) = c("Model", paste(datasource, "Train"),paste(datasource, "Test"))
  }
  return(out)
}
sex_accuracy = merge(list.rbind(sex_pred(ukb, ukb_val,ukb_pc,ukb_test_loadings,"UKB")),
      list.rbind(sex_pred(abcd, abcd_val,abcd_pc,abcd_test_loadings,"ABCD")),
      by = "Model"
      )
write.csv(x = sex_accuracy,file = paste(PATH,"Tables/sex_accuracy.csv",sep=""))
#
# check also the contribution of cortical volume to the prediction expressed
# these are absolute z values for PC1-5, age and cortex vol
out = out2 = list()
for (i in 1:length(ukb_pc)){
    df = cbind(ukb%>%select(c(eid,age,sex,scanner,CortexVol)),ukb_pc[[i]])
    m = glm(sex~PC1+PC2+PC3+PC4+PC5+CortexVol+age+scanner, family = "binomial", data = df)
    out[[i]] = varImp(m)[1:7,]
    df = cbind(abcd%>%select(c(eid,age,sex,scanner,CortexVol)),abcd_pc[[i]])
    m = glm(sex~PC1+PC2+PC3+PC4+PC5+CortexVol+age+scanner, family = "binomial", data = df)
    out2[[i]] = varImp(m)[1:7,]
}
cbind(colMeans(list.rbind(out)),colMeans(list.rbind(out2)))


# m0 = glm(sex~PC1+PC2+PC3+PC4+PC5+age+scanner, family = "binomial", data = df)
# m = glm(sex~PC1+PC2+PC3+PC4+PC5+age+scanner+CortexVol, family = "binomial", data = df)
# m1 = glm(sex~age+scanner+CortexVol, family = "binomial", data = df)
#
# # Likelihood ratio tests comparing Cortical Volume only models vs models with added dMRI PCs
# anova(m,m1)
# #
# # McFadden pseudo R2
# pR2(m0)[4] 
# pR2(m)[4]
# pR2(m1)[4]
# #
# # variable importance
# varImp(m0)
# varImp(m)
# varImp(m1)
#
#
# in a second step, we can focus on the predictors, which are the peincipal components
# we use a simple function, with three main elements:
# orig_dat is the original data, which is only used here to extract covariates
# pc_dat is the data containing the pcs, meaning each diffusion model
# description is the description of the data, e.g. "abcd training data"
apply_glm = function(orig_dat,pc_dat,description){
  df = cbind(orig_dat%>%select(c(eid,age,sex,scanner,CortexVol)),pc_dat)
  m = glm(sex~PC1+PC2+PC3+PC4+PC5+age+scanner+CortexVol, family = "binomial", data = df)
  d = data.frame(Description = description, t(m$coefficients[2:6]))
  return(d)
}
# all of this can be easily put together:
sex_preds = list() # empty list
for (i in 1:length(ukb_pc)){
  sex_preds[[i]] = data.frame(Model = mods[i],rbind(apply_glm(ukb,ukb_pc[[i]],"UKB Training"),
      apply_glm(abcd,abcd_pc[[i]],"ABCD Training"),
      apply_glm(ukb_mean_val,ukb_test_loadings[[i]],"UKB Test"),
      apply_glm(abcd_mean_val,abcd_test_loadings[[i]],"ABCD Test")))
}
# and now, magically transform it all to a cluster of odds ratios
sex_plots = sex_plots_val = list()
for (i in 1:length(ukb_pc)){
sex_plots[[i]] = ggplot(melt(sex_preds[[i]]) %>% mutate(OR = exp(value)), aes(x=variable, y=Description)) +
  geom_tile(colour="black", size=0.25, aes(fill=OR)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 1, mid = "white") +
  #                     breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75)) +
  geom_text(aes(label=round(OR,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")+
  ggtitle(mods[i])
sex_plots_val[[i]] = ggplot(melt(sex_preds[[i]]) %>% mutate(OR = exp(value)), aes(x=variable, y=Description)) +
  geom_tile(colour="black", size=0.25, aes(fill=OR)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 1, mid = "white") +
  #                     breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75)) +
  geom_text(aes(label=round(OR,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")+
  ggtitle(mods[i])
}
sex_plot = ggarrange(plotlist = sex_plots, common.legend = T)
ggsave(paste(PATH,"Figures/sex_plot.pdf",sep=""),sex_plot,width=11,height=5)
#
#
#
# Sex by skeleton average classifications
metric = names(ukb_mean)[5:30]
datalist = list(ukb_mean, ukb_mean_val, abcd_mean, abcd_mean_val)
CoefList = list()
for (j in 1:length(datalist)){
  data = datalist[[j]]
  Coefs = list()
  for (i in 1:length(metric)){
    m = lm(formula(paste(metric[i],"~sex+age+scanner+CortexVol")),data = data)
    Name = metric[i]
    B = (lm.beta::lm.beta(m)$coefficients[2])
    Coefs[[i]] = data.frame(Name,B)
  }
  CoefList[[j]] = list.rbind(Coefs)
}
d = data.frame(Name = CoefList[[1]]$Name,UKB = CoefList[[1]]$B, UKB_val = CoefList[[2]]$B, ABCD = CoefList[[3]]$B, ABCD_val = CoefList[[4]]$B)
d = melt(d)
ggplot(d, aes(x=variable, y=Name)) +
  geom_tile(colour="black", size=0.25, aes(fill=value)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(value,3)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
#
# do the same on the tract level
metric = names(ukb %>% dplyr::select(-c(eid,sex,age,scanner,CortexVol,data)))
datalist = list(ukb, ukb_val, abcd, abcd_val)
CoefList = list()
for (j in 1:length(datalist)){
  data = datalist[[j]]
  Coefs = list()
  for (i in 1:length(metric)){
    m = lm(formula(paste(metric[i],"~sex+age+scanner+CortexVol")),data = data)
    Name = metric[i]
    B = (lm.beta::lm.beta(m)$coefficients[2])
    Coefs[[i]] = data.frame(Name,B)
  }
  CoefList[[j]] = list.rbind(Coefs)
}
d = data.frame(Name = CoefList[[1]]$Name,UKB = CoefList[[1]]$B, UKB_val = CoefList[[2]]$B, ABCD = CoefList[[3]]$B, ABCD_val = CoefList[[4]]$B)
#
# Check which effects replicated in the internal validation sets
# effects replicated in UKB
sum(ifelse(d$UKB>0,1,0) == ifelse(d$UKB_val>0,1,0)) / nrow(d)
# effects replicated in ABCD
sum(ifelse(d$ABCD>0,1,0) == ifelse(d$ABCD_val>0,1,0)) / nrow(d)
# effects replicated in both datasets
sum((ifelse(d$UKB>0,1,0) == ifelse(d$UKB_val>0,1,0)) == (ifelse(d$ABCD>0,1,0) == ifelse(d$ABCD_val>0,1,0)))/nrow(d)
#
# check for the strength and replicability of all, hippo and UF effects
### all
d %>% summarize(M_UKB = mean(abs(UKB)),M_UKB_val = mean(abs(UKB_val)),
                                     M_ABCD = mean(abs(ABCD)),M_ABCD_val = mean(abs(ABCD_val)))

### Uncinate fasciculus
d[grepl("UF",d$Name),] %>% summarize(M_UKB = mean(abs(UKB)),M_UKB_val = mean(abs(UKB_val)),
                                   M_ABCD = mean(abs(ABCD)),M_ABCD_val = mean(abs(ABCD_val)))
### hippocampus
d[grepl("hippocampus",d$Name),] %>% summarize(M_UKB = mean(abs(UKB)),M_UKB_val = mean(abs(UKB_val)),
                                     M_ABCD = mean(abs(ABCD)),M_ABCD_val = mean(abs(ABCD_val)))




# 6. AGE ASSOCIATIONS ####
age_pred = function(orig_train,orig_test, train_list, test_list, datasource){
  for (i in 1:length(train_list)){
    
    df = cbind(orig_train%>%select(c(eid,age,sex,scanner,CortexVol)),train_list[[i]])
    m = lm(age~PC1+PC2+PC3+PC4+PC5+CortexVol+sex, data = df)
    testing = cbind(orig_test%>%select(c(eid,age,sex,scanner,CortexVol)),test_list[[i]])
    Test = cor(predict(m, newdata=testing),testing$age,use = "pairwise.complete.obs")
    Train = cor(predict(m, newdata=df),df$age,use = "pairwise.complete.obs")
    out[[i]] = data.frame(Model = mods[i], round(Train,4), round(Test,4))
    names(out[[i]]) = c("Model", paste(datasource, "Train"),paste(datasource, "Test"))
  }
  return(out)
}
age_corrs = merge(list.rbind(age_pred(ukb, ukb_val,ukb_pc,ukb_test_loadings,"UKB")),
                     list.rbind(age_pred(abcd, abcd_val,abcd_pc,abcd_test_loadings,"ABCD")),
                     by = "Model"
)
write.csv(x = age_corrs,file = paste(PATH,"Tables/age_corrs.csv",sep=""))
#
# beta coefficients of principal components
#
apply_lm = function(orig_dat,pc_dat,description){
  df = cbind(orig_dat%>%select(c(eid,age,sex,scanner,CortexVol)),pc_dat)
  m = lm(age~PC1+PC2+PC3+PC4+PC5+sex+scanner+CortexVol, data = df)
  d = data.frame(Description = description, Component = names(lm.beta::lm.beta(m)$coefficients[2:6]), 
                 lm.beta::lm.beta(m)$coefficients[2:6])
  return(d)
}
# all of this can be easily put together:
age_preds = list() # empty list
for (i in 1:length(ukb_pc)){
  age_preds[[i]] = data.frame(Model = mods[i],rbind(apply_lm(ukb,ukb_pc[[i]],"UKB Training"),
                                                    apply_lm(abcd,abcd_pc[[i]],"ABCD Training"),
                                                    apply_lm(ukb_mean_val,ukb_test_loadings[[i]],"UKB Test"),
                                                    apply_lm(abcd_mean_val,abcd_test_loadings[[i]],"ABCD Test")))
}
age_plots = list()
for (i in 1:length(ukb_pc)){
  age_plots[[i]] = ggplot(melt(age_preds[[i]]) %>% mutate(Beta = value), aes(x=Component, y=Description)) +
    geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
    scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white") +
    #                     breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75)) +
    geom_text(aes(label=round(Beta,2)),size=4)+
    theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")+
    ggtitle(mods[i])
}
age_plot = ggarrange(plotlist = age_plots, common.legend = F)
ggsave(paste(PATH,"Figures/age_plot.pdf",sep=""),age_plot,width=13,height=5)

# small correlations with age might be due to the cortical volume being more age-associated
apply_lm = function(orig_dat,pc_dat,description){
  df = cbind(orig_dat%>%select(c(eid,age,sex,scanner,CortexVol)),pc_dat)
  m = lm(age~CortexVol+PC1+PC2+PC3+PC4+PC5+sex+scanner, data = df)
  d = data.frame(Description = description, Component = names(lm.beta::lm.beta(m)$coefficients[2]), 
                 lm.beta::lm.beta(m)$coefficients[2])
  return(d)
}
# all of this can be easily put together:
age_preds = list() # empty list
for (i in 1:length(ukb_pc)){
  age_preds[[i]] = data.frame(Model = mods[i],rbind(apply_lm(ukb,ukb_pc[[i]],"UKB Training"),
                                                    apply_lm(abcd,abcd_pc[[i]],"ABCD Training"),
                                                    apply_lm(ukb_mean_val,ukb_test_loadings[[i]],"UKB Test"),
                                                    apply_lm(abcd_mean_val,abcd_test_loadings[[i]],"ABCD Test")))
}
list.rbind(age_preds)
# Well ...
# the small correlations in abcd data are not due to confounding cortical volumes
#
#
#
age_pred = function(orig_train,orig_test, train_list, test_list, datasource){
  for (i in 1:length(train_list)){
    df = cbind(orig_train%>%select(c(eid,age,sex,scanner,CortexVol)),train_list[[i]])
    m = lm(age~PC1+PC2+PC3+PC4+PC5+CortexVol+sex, data = df)
    testing = cbind(orig_test%>%select(c(eid,age,sex,scanner,CortexVol)),test_list[[i]])
    Test = MAE(predict(m, newdata=testing),testing$age,na.rm = T)
    Train = MAE(predict(m, newdata=df),df$age,na.rm = T)
    out[[i]] = data.frame(Model = mods[i], round(Train,4), round(Test,4))
    names(out[[i]]) = c("Model", paste(datasource, "Train"),paste(datasource, "Test"))
  }
  return(out)
}
age_rmse = merge(list.rbind(age_pred(ukb, ukb_val,ukb_pc,ukb_test_loadings,"UKB")),
                  list.rbind(age_pred(abcd, abcd_val,abcd_pc,abcd_test_loadings,"ABCD")),
                  by = "Model"
)
age_rmse # added as a supplement

# Skeleton level age associations
res = list()
mean_list=list(ukb_mean,abcd_mean)
Beta=SE=p=c()
for (df in 1:length(mean_list)){
  for (i in 1:length(preds)){
    tmp_dat = stand(mean_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+age+sex+scanner+CortexVol",sep="")) # regular linear models used due to singular fit
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = lm.beta::lm.beta(m1)$coefficients[2]
    #SE[i] = summary(m1)$coefficients[3,2] # not necessary for the plotting
    p[i] = summary(m1)$coefficients[2,4]
  }
  res[[df]] = data.frame(preds, Beta, p) #SE
}
res[[1]]$Data = "UKB"
res[[2]]$Data = "ABCD"
plot_dat = list.rbind(res) # merge list/data frames
plot_dat$Metric=deps2 # add plotable metric names
# split data up into 3 parts for better plotting
BRIA.df=dplyr::filter(plot_dat, grepl("BRIA",Metric))
DTI.DKI.df = rbind(plot_dat %>% filter(grepl("DKI",Metric)),plot_dat %>% filter(grepl("DTI",Metric)))
SMT.WMTI.df = rbind(plot_dat %>% filter(grepl("SMT",Metric)),plot_dat %>% filter(grepl("WMTI",Metric)))
res=list(BRIA.df,DTI.DKI.df,SMT.WMTI.df)
for (i in 1:length(res)){
  res[[i]] = melt(res[[i]]%>%select(-SE,-p,-preds))%>% rename(Beta = value)
}
p1 = ggplot(res[[1]], aes(x=Data, y=Metric)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2 = ggplot(res[[2]], aes(x=Data, y=Metric)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p3 = ggplot(res[[3]], aes(x=Data, y=Metric)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
#plot10 = ggarrange(p1,p2,p3,nrow=1, common.legend = T, legend = "bottom") #c("UKB", "ABCD", "Both")
plot11=plot_grid(p1,p2,p3,nrow=1, align="v")
ggsave(paste(PATH,"Figures/age_plot_skeleton.pdf",sep=""),plot11,width=11,height=5)
#
# HIPPOCAMPUS tract level associations
# define hippo tracts
preds = abcd %>% select(contains("hippo")) %>% names
preds = preds[!grepl("awf",preds)]
res = res2 = list()
# predict for training data
tract_list = list(abcd,ukb)
Beta=SE=p=c()
for (df in 1:length(tract_list)){
  for (i in 1:length(preds)){
    tmp_dat = stand(tract_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+age+sex+scanner+CortexVol",sep="")) # regular linear models used due to singular fit
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = lm.beta::lm.beta(m1)$coefficients[2]
    #SE[i] = summary(m1)$coefficients[3,2] # not necessary for the plotting
    p[i] = summary(m1)$coefficients[2,4]
  }
  res[[df]] = data.frame(preds, Beta, p) #SE
}
# now, do the same for test data (standardisation already done based on training data)
# Note that these values have been standardized based on the training samples
tract_list=list(data.frame(abcd_pc_val), data.frame(ukb_pc_val))
orig_list=list(abcd_val,ukb_val)
for (df in 1:length(tract_list)){
  for (i in 1:length(preds)){
    tmp_dat = cbind(orig_list[[df]] %>% select(age,sex,scanner,CortexVol), tract_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+age+sex+scanner+CortexVol",sep="")) # regular linear models used due to singular fit
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = lm.beta::lm.beta(m1)$coefficients[2]
    #SE[i] = summary(m1)$coefficients[3,2] # not necessary for the plotting
    p[i] = summary(m1)$coefficients[2,4]
  }
  res2[[df]] = data.frame(preds, Beta, p) #SE
}
sum((res[[1]]$Beta > 0) == (res2[[1]]$Beta > 0))/nrow(res[[1]])
sum((res[[2]]$Beta > 0) == (res2[[2]]$Beta > 0))/nrow(res[[2]])

mean(abs(res[[1]]$Beta))
mean(abs(res2[[1]]$Beta))

mean(abs(res[[2]]$Beta))
mean(abs(res2[[2]]$Beta))
#
#
# 7. ASYMMETRIES ######
left = ukb %>% select(-CortexVol) %>% select(ends_with("L")) %>% names
right = ukb %>% select(-scanner) %>% select(ends_with("R")) %>% names
dflist = list(ukb,ukb_val,abcd,abcd_val)
the_result = list()
for (j in 1:length(dflist)){
  resu = list()
  data = dflist[[j]]
  for (i in 1:length(left)){
    tmp = data.frame(value = c(unlist(data[left[i]]), unlist(data[right[i]])),
                     hemi = c(replicate(nrow(data),"left"),replicate(nrow(data),"right")))
    d = cohens_d(formula = value~hemi,data = tmp)$effsize # NOTE: THIS IS LEFT MINUS RIGHT!
    t = t.test(value~hemi,tmp)$p.value # This is the p-val only
    resu[[i]] = data.frame(Name = left[i], Cohens_d = d,p = t)
  }
  the_result[[j]] = resu
}
deps = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", "DTI-FA", "DTI-MD", "DTI-RD", "DTI-AD", 
         "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
ukb_hemi = list.rbind(the_result[[1]])
ukb_val_hemi = list.rbind(the_result[[2]])
abcd_hemi = list.rbind(the_result[[3]])
abcd_val_hemi = list.rbind(the_result[[4]])
#hemilist = list(ukb_hemi,ukb_val_hemi,abcd_hemi,abcd_val_hemi) # put them all together into a list
# define metric names in a plot-friendly way
ukb_hemi$Metric = ukb_val_hemi$Metric = abcd_hemi$Metric = abcd_val_hemi$Metric = deps
# define tract names in a plot-friendly way
## this defines the tracts
gsub("L","",gsub("microFA_","",ukb_hemi$Name[grepl("microFA_",ukb_hemi$Name)]))
## yet, the formatting is not perfect, hence done by hand:
ukb_hemi$Tract = ukb_val_hemi$Tract = abcd_hemi$Tract = abcd_val_hemi$Tract = 
  c(replicate(length(deps),"ATR"),
    replicate(length(deps),"CST"),
    replicate(length(deps),"CG"),
    replicate(length(deps),"CG-hippocampus"),
    replicate(length(deps),"IFOF"),
    replicate(length(deps),"ILF"),
    replicate(length(deps),"SLF"),
    replicate(length(deps),"UF"),
    replicate(length(deps),"SLFT"))
# Make a plotting function
plot_hemi = function(title,data_hemi){
    p=ggplot(data_hemi, aes(x=Metric, y=Tract)) +
    geom_tile(colour="black", size=0.25, aes(fill=Cohens_d)) +
    scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white") +
    #breaks = c(-1, -0.5, 0, 0.25, 0.5, 1)) +
    geom_text(aes(label=round(Cohens_d,2)),size=4)+
    theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
  p = p+theme(axis.text.x = element_text(size = 8, angle = 90))+ggtitle(title)
  return(p)
}
hemiplots = list(plot_hemi("UKB  training",ukb_hemi),plot_hemi("UKB validation",ukb_val_hemi),
                 plot_hemi("ABCD training",abcd_hemi),plot_hemi("ABCD validation",abcd_val_hemi))
ggarrange(plotlist = hemiplots)
