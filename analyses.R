# analyses to examine the extend to which DTI explains white matter microstructure
# Max Korbmacher, 27 August 2024
#
# R version: 4.2.1
#
#
#
# DEFINE PATH FOR DATA AND OUTPUT
PATH = "your/data/path/
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
               viridis,cowplot,
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
# transform abcd age from months to years
abcd_mean$age = abcd_mean$age/12
abcd$age = abcd$age/12
ukb$sex = ukb_mean$sex=factor(ukb_mean$sex)
levels(ukb$sex) = levels(ukb_mean$sex) = c("female", "male")
abcd$sex = abcd_mean$sex = factor(abcd$sex)
dat_mean = rbind(ukb_mean,abcd_mean)
dat = rbind(ukb,abcd)
#
# 0) Power #### 
# quick and dirty power analysis for univariate associations with the lowest power
##ABCD data only
pwr_out = pwr.f2.test(u = 4,v = 5142,sig.level = 0.05,power = .95)
print(paste("lowest observable f2 and R2 = ",(pwr_out$f2),"(ABCD)",sep = ""))
print(paste("lowest observable r = ",sqrt(pwr_out$f2),"(ABCD)",sep = ""))
## UKB data only
pwr_out = pwr.f2.test(u = 4,v = 42224,sig.level = 0.05,power = .95)
print(paste("lowest observable f2 and R2 = ",(pwr_out$f2),"(UKB)",sep = ""))
print(paste("lowest observable r = ",sqrt(pwr_out$f2),"(UKB)",sep = ""))
## all data
pwr_out = pwr.f2.test(u = 4,v = 47371,sig.level = 0.05,power = .95)
print(paste("lowest observable f2 and R2 = ",(pwr_out$f2),"(all data)",sep = ""))
print(paste("lowest observable r = ",sqrt(pwr_out$f2),"(all data)",sep = ""))
# u = 4 for metric, sex, site, age
# v = 5,147 (ABCD datasets) - 5
# 1) Check skeleton averages ####
# 
# make list of dependent variables
dependent = dat_mean %>% select(-c(md_Mean,rd_Mean, ad_Mean, FA_Mean, eid, sex, age, scanner)) %>% names
# standardize independent and dependent variables to obtain standardized (and comparable coefficients)
stand = function(dataframe){
  newdataframe = dataframe %>% select(-c(eid, sex, age, scanner)) %>% scale %>% data.frame()
  newdataframe = cbind(newdataframe, dataframe %>% select(c(eid, sex, age, scanner)))
  return(newdataframe)
}
# create list of dfs
mean_list = list(stand(ukb_mean), stand(abcd_mean), stand(dat_mean))
# 1.1) Correlations controlled for age, sex and scanner ####
# create a list for the results
res = list()
for (df in 1:length(mean_list)){
  AD = AD.SE = FA = FA.SE = RD = RD.SE = MD = MD.SE = c()
  for (i in 1:length(dependent)){
    f1 = formula(paste(dependent[i],"~ad_Mean+age+sex+(scanner)",sep="")) # regular linear models used due to singular fit 
    f2 = formula(paste(dependent[i],"~FA_Mean+age+sex+(scanner)",sep="")) # otherwise lmer with site as random intercept
    f3 = formula(paste(dependent[i],"~rd_Mean+age+sex+(scanner)",sep=""))
    f4 = formula(paste(dependent[i],"~md_Mean+age+sex+(scanner)",sep=""))
    m1 = lm(f1, data=mean_list[[df]])
    m2 = lm(f2, data=mean_list[[df]])
    m3 = lm(f3, data=mean_list[[df]])
    m4 = lm(f4, data=mean_list[[df]])
    AD[i] = summary(m1)$coefficients[2]
    AD.SE[i] = summary(m1)$coefficients[2,2]
    FA[i] = summary(m2)$coefficients[2]
    FA.SE[i] = summary(m2)$coefficients[2,2]    
    RD[i] = summary(m3)$coefficients[2]
    RD.SE[i] = summary(m3)$coefficients[2,2]
    MD[i] = summary(m4)$coefficients[2]
    MD.SE[i] = summary(m4)$coefficients[2,2]
    }
  res[[df]] = data.frame(dependent, AD, AD.SE, FA, FA.SE, RD, RD.SE, MD, MD.SE)
}
deps = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
for (i in 1:length(res)){
  res[[i]]$dependent = deps
  res[[i]] = melt(res[[i]]%>%select(!contains("SE"))) %>% rename(Beta = value)
}
p1 = ggplot(res[[1]], aes(x=variable, y=dependent)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                      breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2 = ggplot(res[[2]], aes(x=variable, y=dependent)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p3 = ggplot(res[[3]], aes(x=variable, y=dependent)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
plot1 = ggarrange(p1,p2,p3,nrow=1,labels = c("a", "b", "c"), common.legend = T, legend = "right") #c("UKB", "ABCD", "Both")
ggsave(paste(PATH,"Figures/plot1.pdf",sep=""),plot1,width=9,height=7)

# 1.2) Crude correlations ####
# create a list for the results
res2 = list()
for (df in 1:length(mean_list)){
  AD = AD.SE = FA = FA.SE = RD = RD.SE = MD = MD.SE = c()
  for (i in 1:length(dependent)){
    f1 = formula(paste(dependent[i],"~ad_Mean",sep=""))
    f2 = formula(paste(dependent[i],"~FA_Mean",sep=""))
    f3 = formula(paste(dependent[i],"~rd_Mean",sep=""))
    f4 = formula(paste(dependent[i],"~md_Mean",sep=""))
    m1 = lm(f1, data=mean_list[[df]])
    m2 = lm(f2, data=mean_list[[df]])
    m3 = lm(f3, data=mean_list[[df]])
    m4 = lm(f4, data=mean_list[[df]])
    AD[i] = summary(m1)$coefficients[2] 
    AD.SE[i] = summary(m1)$coefficients[2,2]
    FA[i] = summary(m2)$coefficients[2]
    FA.SE[i] = summary(m2)$coefficients[2,2]    
    RD[i] = summary(m3)$coefficients[2]
    RD.SE[i] = summary(m3)$coefficients[2,2]
    MD[i] = summary(m4)$coefficients[2]
    MD.SE[i] = summary(m4)$coefficients[2,2]
  }
  res2[[df]] = data.frame(dependent, AD, AD.SE, FA, FA.SE, RD, RD.SE, MD, MD.SE)
}
deps = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
for (i in 1:length(res2)){
  res2[[i]]$dependent = deps
  res2[[i]] = melt(res2[[i]]%>%select(!contains("SE"))) %>% rename(Beta = value)
}
p1 = ggplot(res2[[1]], aes(x=variable, y=dependent)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2 = ggplot(res2[[2]], aes(x=variable, y=dependent)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p3 = ggplot(res2[[3]], aes(x=variable, y=dependent)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
plot2 = ggarrange(p1,p2,p3,nrow=1,labels = c("a", "b", "c"), common.legend = T, legend = "right") #c("UKB", "ABCD", "Both")
ggsave(paste(PATH,"Figures/plot2.pdf",sep=""),plot2,width=9,height=7)

#
#
#
#
#
# Summary stats for the skeleton associations
# Note, here we assess only the corrected associations, not the crude ones.
res[[1]] %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
a = res[[1]] %>% group_by(variable) %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
res[[2]] %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
b = res[[2]] %>% group_by(variable) %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
res[[3]] %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
c = res[[3]] %>% group_by(variable) %>% summarize(M = mean(abs(Beta)),SD = sd(abs(Beta)))
cbind(a,b,c)
rm(a,b,c)
#
# 2) Crude correlations between all skeleton averages ####
deps2 = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK","DTI-FA", "DTI-MD","DTI-RD","DTI-AD", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
cormat = function(df){
  cor1 = df%>%select(-c(eid,age,sex,scanner))
  #roworder = names(cor1)
  names(cor1) = deps2
  cor1 = cor(cor1)
  return(cor1)
}
a = as.ggplot(pheatmap(cormat(mean_list[[1]]),cluster_rows = F, cluster_cols = F))
b = as.ggplot(pheatmap(cormat(mean_list[[2]]),cluster_rows = F, cluster_cols = F))
c = as.ggplot(pheatmap(cormat(mean_list[[3]]),cluster_rows = F, cluster_cols = F))
d = as.ggplot(pheatmap(cormat(mean_list[[1]]),cluster_rows = F, cluster_cols = T))
e = as.ggplot(pheatmap(cormat(mean_list[[2]]),cluster_rows = F, cluster_cols = T))
f = as.ggplot(pheatmap(cormat(mean_list[[3]]),cluster_rows = F, cluster_cols = T))
plot3 = a+b+c+d+e+f
# note: top row = unclustered ukb,abcd,both; bottom row = clustered
ggsave(paste(PATH,"Figures/plot3.pdf",sep=""),plot3,width=20,height=14)
rm(mean_list,plot3)
# 3) Principal component analysis (PCA) ####
# 3.1) Estimation and visualization of the PCs ####
# estimate principal components and visualise the first 10 components
pcplots = function(data){
dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
a = fviz_eig(dkipc, main = "DKI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
b = fviz_eig(dtipc, main = "DTI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
briapc =  data %>% select(starts_with("v_"), starts_with("micro_"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
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
ggsave(paste(PATH,"Figures/plot4_1.pdf",sep=""),plot4.1,width=14,height=10)
ggsave(paste(PATH,"Figures/plot4_2.pdf",sep=""),plot4.2,width=14,height=10)
ggsave(paste(PATH,"Figures/plot4_3.pdf",sep=""),plot4.3,width=14,height=10)
rm(plot4.1,plot4.2,plot4.3)
#
# 3.2) Associations of the PCs with each other ####
# extract the first component
pccor = function(data){
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro_"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
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
  briapc =  data %>% select(starts_with("v_"), starts_with("micro_"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
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
# 3.3) Associations of the PCs with global/skeleton DTI metrics ####
# estimate corrected correlations, standard errors and p-values
cors = function(data,mean_data){ #data argument is for the data frame containing tract-level data, mean_data for skeleton-level data
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro_"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
  PC = data.frame(DKI = as.numeric(dkipc$x[,1]), DTI = as.numeric(dtipc$x[,1]), BRIA = as.numeric(briapc$x[,1]), 
                  SMT = as.numeric(smtpc$x[,1]), SMTmc = as.numeric(smtmcpc$x[,1]), WMTI = as.numeric(wmtipc$x[,1]))
  PC = data.frame(scale(PC))
  PC$eid = data$eid
  predictor_vars = c("ad_Mean","rd_Mean","md_Mean", "FA_Mean")
  for (n in 1:length(predictor_vars)){
    mean_data[predictor_vars[n]] = c(scale(mean_data[predictor_vars[n]]))
  }
  PC = merge(mean_data,PC)
  WMMapproach = c("BRIA","DKI", "DTI", "SMT", "SMTmc", "WMTI")
  B = SE = p = c()
  export = list()
  for (j in 1:length(predictor_vars)){
    for (i in 1:length(WMMapproach)){
      #f1 = formula(paste(WMMapproach[i],"~DTI+age+sex+(1|scanner)",sep="")) # singulary issues using lmer 
      f1 = formula(paste(WMMapproach[i],"~",predictor_vars[j],"+age+sex+(scanner)",sep=""))
      mod = lm(f1,PC)
      B[i] = summary(mod)$coefficients[2,1]
      SE[i] = summary(mod)$coefficients[2,2]
      p[i] = summary(mod)$coefficients[2,4]
    }
    export[[j]] = data.frame(B,SE,p)
    names(export[[j]]) = c("B", "SE","p")
    rownames(export[[j]]) = WMMapproach
  }
  names(export)=predictor_vars
  return(export)
}
ukblist = cors(ukb,ukb_mean)
abcdlist = cors(abcd,abcd_mean)
datlist = cors(dat,dat_mean)
plt_list = list()
for (i in 1:length(ukblist)){
  tmp = rbind(ukblist[[i]], abcdlist[[i]], datlist[[i]])
  tmp$Approach = row.names(ukblist[[i]])
  tmp$Data = c(replicate(nrow(tmp)/3,"UKB"),replicate(nrow(tmp)/3,"ABCD"),replicate(nrow(tmp)/3,"Both"))
  tmp$B = abs(tmp$B)
  plt_list[[i]] = ggplot(tmp, aes(x=Approach, y=B, color = Data)) +
    geom_errorbar(aes(ymin=B-SE, ymax=B+SE,group=Data), width=.2,position=position_dodge(.6)) +
    geom_point(position = position_dodge(width = 0.6))+
    scale_color_manual(values=c("#880700", "#E69F00", "#3a81b5"))+ ##999999','#E69F00', '#56B4E9'
    coord_flip() + theme_bw() + theme(text = element_text(size=15)) + ylab("")+xlab("")+ylim(0.25,1.3)
}
plot7 = ggarrange(plotlist = plt_list, common.legend = T, legend = "right", labels = c("AD", "RD", "MD", "FA"))
plot7 = annotate_figure(plot7, left = text_grob("Diffusion Approach", rot = 90), bottom = text_grob("Adjusted Absolute Standardized Beta Coefficients"))
ggsave(paste(PATH,"Figures/plot7.pdf",sep=""),plot7,width=8,height=7)
# have a quick look at the stats for each data set
statcheck = function(dflist){
  for (i in 1:length(dflist)){
    print(paste(names(dflist)[i]," Mean = ",mean(abs(dflist[[names(dflist)[i]]]$B))," ± ",sd(abs(dflist[[names(dflist)[i]]]$B)),sep=""))
  }
}
print("Checking the associations between DTI skeleton-level metrics and PCs of the different approaches")
print("In UKb:")
statcheck(ukblist)
print("In ABCD")
statcheck(abcdlist)
print("In combined data:")
statcheck(datlist)
rm(ukblist,abcdlist,datlist,plt_list)
#
# 4) Tract-level associations (within-tract associations) ####
# to automatically detect all tracts there are in the data, we take one of the three data frames (containing identical tracts and tract labels!)
tract_names = abcd%>%select(contains("smt_mc_intra_"))%>%names
# we then remove the specific metric name and get the clean tract names
tract_names = gsub('smt_mc_intra_', "", tract_names)
# we use these names to loop over the data frames to associate DTI metrics with the advanced diffusion metrics
# to get the lst of the advanced diffusion metrics, we use the same logic, but this time sorting first by tract (except DTI metrics)
metric_names = abcd_mean %>% select(contains("Mean")) %>% select(-ad_Mean,-rd_Mean,-md_Mean,-FA_Mean) %>% names
metric_names = gsub('_Mean', "", metric_names)
dti_names = c("ad", "rd", "md", "FA") # similar: abcd%>%select(contains("CGL")) %>% select(ad_CGL, FA_CGL, md_CGL, rd_CGL) %>% names
# now loop over this
tmp1=tmp2=tmp3=data.frame(matrix(nrow = length(metric_names), ncol = 4*3)) # ncol = #DTImetricts*3(beta,sd,p)
tmp1_lst=tmp2_lst=tmp3_lst=list()
stand_ukb = stand(ukb)
stand_abcd = stand(abcd)
#dat=subset(dat, select = -data)
stand_dat = stand(dat)
for (tract in 1:length(tract_names)){
  # each data source by tract (standardize numerics as well on the way there)
  ukb_tmp = stand_ukb %>% select(contains(tract_names[tract]), c(sex,age,scanner)) %>% mutate(age=scale(age))
  abcd_tmp = stand_abcd %>% select(contains(tract_names[tract]), c(sex,age,scanner)) %>% mutate(age=scale(age))
  dat_tmp = stand_dat %>% select(contains(tract_names[tract]), c(sex,age,scanner)) %>% mutate(age=scale(age))
  # now, we loop through the metrics
  for (metric in 1:length(metric_names)){
    # make a global formula applied to these tracts (remember DTI predicts all else)
    ad_formula = formula(paste(metric_names[metric],"_",tract_names[tract],"~ad_",tract_names[tract],"+age+sex+scanner", sep=""))
    rd_formula = formula(paste(metric_names[metric],"_",tract_names[tract],"~rd_",tract_names[tract],"+age+sex+scanner", sep=""))
    md_formula = formula(paste(metric_names[metric],"_",tract_names[tract],"~md_",tract_names[tract],"+age+sex+scanner", sep=""))
    FA_formula = formula(paste(metric_names[metric],"_",tract_names[tract],"~FA_",tract_names[tract],"+age+sex+scanner", sep=""))
    # run linear model for each metric
    run_lm= function(df){
      # define models
      mod1=lm(ad_formula,df)
      mod2=lm(rd_formula,df)
      mod3=lm(md_formula,df)
      mod4=lm(FA_formula,df)
      # extract coefficients
      B.ad = summary(mod1)$coefficients[2,1]
      SE.ad = summary(mod1)$coefficients[2,2]
      p.ad = summary(mod1)$coefficients[2,4]
      B.rd = summary(mod2)$coefficients[2,1]
      SE.rd = summary(mod2)$coefficients[2,2]
      p.rd = summary(mod2)$coefficients[2,4]
      B.md = summary(mod3)$coefficients[2,1]
      SE.md = summary(mod3)$coefficients[2,2]
      p.md = summary(mod3)$coefficients[2,4]
      B.FA = summary(mod4)$coefficients[2,1]
      SE.FA = summary(mod4)$coefficients[2,2]
      p.FA = summary(mod4)$coefficients[2,4]
      tmp=c(B.ad,SE.ad,p.ad,B.rd,SE.rd,p.rd,B.md,SE.md,p.md,B.FA,SE.FA,p.FA)
      return(tmp)
      }
    tmp1[metric,] = run_lm(ukb_tmp)
    tmp2[metric,] = run_lm(abcd_tmp)
    tmp3[metric,] = run_lm(dat_tmp)
  }
  tmp1_lst[[tract]]=tmp1
  tmp2_lst[[tract]]=tmp2
  tmp3_lst[[tract]]=tmp3
}
# add names to the data frames (including nicer metric names, order is correct)
metric_names = c("SMTmc-Diff", "SMT-Diff","BRIA-DAXextra", "BRIA-DAXintra","BRIA-DRADextra","BRIA-microADC",
  "BRIA-microAX", "BRIA-microFA","BRIA-microRD", "BRIA-vCSF","BRIA-Vextra","BRIA-Vintra","DKI-AK","DKI-MK","DKI-RK",
  "SMT-long","SMTmc-extraMD","SMTmc-extratrans", "SMTmc_intra","WMTI-AWF","WMTI-axEAD", "WMTI-radEAD")
for (i in 1:length(tmp1_lst)){
  names(tmp1_lst[[i]])=names(tmp2_lst[[i]])=names(tmp3_lst[[i]])=c("B.ad","SE.ad","p.ad","B.rd","SE.rd","p.rd","B.md","SE.md","p.md","B.FA","SE.FA","p.FA")
  tmp1_lst[[i]]$Metric=tmp2_lst[[i]]$Metric=tmp3_lst[[i]]$Metric=metric_names
  tmp1_lst[[i]]$Data = "UKB"
  tmp2_lst[[i]]$Data = "ABCD"
  tmp3_lst[[i]]$Data = "Both"
}
# bring tables into right shape to be outputted as figure
tract_names= gsub('CG_hippocampus', "CG-Hip", tract_names) # more representable tract names
counter=c(replicate(4,LETTERS[1:(length(metric_names))]))#make a counter for the metrics instead of displaying all
bigplot=function(tmplst){
  plotlist=list()
  for (i in 1:length(tmplst)){
    plotlist[[i]]=tmplst[[i]] %>% select(contains("B."),Metric)%>%`colnames<-`(c("AD","RD","MD","FA","Metric"))%>% melt%>%mutate(Tract=tract_names[i],Metric=counter,Beta=value)%>%
      ggplot(aes(x=Metric, y=variable)) + 
      geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
      scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                           breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
      theme(axis.text.y = element_text(size = 8),base_size = 8) + 
      theme_bw() + xlab("") + ylab("") + labs(title=tract_names[i])
    plotlist[[i]]=plotlist[[i]]+theme(plot.title = element_text(size = 10,face="italic"))
  }
  plot=ggarrange(plotlist = plotlist,common.legend = T, legend = "right")
  plot=plot+theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
  return(plot)
}
# plot it
plot8=ggarrange(annotate_figure(bigplot(tmp1_lst), left = text_grob("UKB", rot = 90)),
                annotate_figure(bigplot(tmp2_lst), left = text_grob("ABCD", rot = 90)),
                annotate_figure(bigplot(tmp3_lst), left = text_grob("Both", rot = 90)),ncol=1)
ggsave(paste(PATH,"Figures/plot8.pdf",sep=""),plot8,width=17,height=15)
# plot option for long horizontal fig
# plot8=ggarrange(bigplot(tmp1_lst), bigplot(tmp2_lst),bigplot(tmp3_lst),nrow=1)
# ggsave(paste(PATH,"Figures/plot8.pdf",sep=""),plot8,width=40,height=8)
# one can also go for single plots, which hasnt been done here
# plot8.1=bigplot(tmp1_lst)
# plot8.2=bigplot(tmp2_lst)
# plot8.3=bigplot(tmp3_lst)
# ggsave(paste(PATH,"Figures/plot8_1.pdf",sep=""),plot8.1,width=20,height=8)
# ggsave(paste(PATH,"Figures/plot8_2.pdf",sep=""),plot8.2,width=20,height=8)
# ggsave(paste(PATH,"Figures/plot8_3.pdf",sep=""),plot8.3,width=20,height=8)
#
# write full table with SE and p-vals as csv
for (i in 1:length(tmp1_lst)){
  tmp1_lst[[i]]$Tract=tract_names[i]
  tmp2_lst[[i]]$Tract=tract_names[i]
  tmp3_lst[[i]]$Tract=tract_names[i]
}
tab1 = rbind(data.frame(list.rbind(tmp1_lst)),data.frame(list.rbind(tmp2_lst)),data.frame(list.rbind(tmp3_lst)))
write.csv(x=tab1, file=paste(PATH,"tract_level_associations.csv",sep=""))
# check some stats
# have a quick look at the stats for each data set
tract_avg = tab1 %>%group_by(Data)%>%summarize(AD.mean = mean(abs(B.ad)),AD.sd = sd(abs(B.ad)),
                                   RD.mean = mean(abs(B.rd)),rd.sd = sd(abs(B.rd)),
                                   md.mean = mean(abs(B.md)),md.sd = sd(abs(B.md)),
                                   FA.mean = mean(abs(B.FA)),FA.sd = sd(abs(B.FA)))
write.csv(x=tract_avg, file=paste(PATH,"tract_level_associations_avg.csv",sep=""))
# check stats by metric
tract_avg = tab1 %>%group_by(Data,Metric)%>%summarize(AD.mean = mean(abs(B.ad)),AD.sd = sd(abs(B.ad)),
                                               RD.mean = mean(abs(B.rd)),rd.sd = sd(abs(B.rd)),
                                               md.mean = mean(abs(B.md)),md.sd = sd(abs(B.md)),
                                               FA.mean = mean(abs(B.FA)),FA.sd = sd(abs(B.FA)))
write.csv(x=tract_avg, file=paste(PATH,"tract_level_associations_by_Data_and_Metric.csv",sep=""))

#remove obsolete dfs
rm(tmp1_lst,tmp2_lst,tmp3_lst,tab1)

# 5) Age-associations ####
# 5.1) Lifespan plots ####
# start with visualising age-associations (skeleton level)
preds=names(dat_mean[5:length(dat_mean)]) # create vector of variables to be predicted
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
# 5.1) Linear models for comparison of beta coefficients ASSESSING AGE ASSOCIATIONS ####
# prep models for skeleton-level associations
res = list()
mean_list=list(ukb_mean,abcd_mean)
Beta=SE=c()
for (df in 1:length(mean_list)){
  for (i in 1:length(preds)){
    tmp_dat = stand(mean_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+age+sex+scanner",sep="")) # regular linear models used due to singular fit 
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = summary(m1)$coefficients[2]
    SE[i] = summary(m1)$coefficients[2,2] # not necessary for the plotting
  }
  res[[df]] = data.frame(preds, Beta, SE)
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
  res[[i]] = melt(res[[i]]%>%select(-SE,-preds))%>% rename(Beta = value)
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
ggsave(paste(PATH,"Figures/plot10.pdf",sep=""),plot10,width=11,height=5)
#
#
#
# 5.2) Linear models for comparison of beta coefficients ASSESSING SEX DIFFERENCES ####
# prep models for skeleton-level associations
res = list()
mean_list=list(ukb_mean,abcd_mean)
Beta=SE=p=c()
for (df in 1:length(mean_list)){
  for (i in 1:length(preds)){
    tmp_dat = stand(mean_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+age+sex+scanner",sep="")) # regular linear models used due to singular fit 
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = summary(m1)$coefficients[3]
    SE[i] = summary(m1)$coefficients[3,2] # not necessary for the plotting
    p[i] = summary(m1)$coefficients[3,4]
  }
  res[[df]] = data.frame(preds, Beta, SE, p)
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
ggsave(paste(PATH,"Figures/plot11.pdf",sep=""),plot11,width=11,height=5)
