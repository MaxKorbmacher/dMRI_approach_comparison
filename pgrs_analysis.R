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
               viridis,cowplot,effsize,
               update = F)
# load data
## tract
ukb = read.csv(paste(PATH,"ukb_dat.csv",sep=""))
abcd = read.csv(paste(PATH,"abcd_dat.csv",sep=""))
tracts=list(ukb,abcd)
## mean values
ukb = read.csv(paste(PATH,"ukb_mean.csv",sep=""))
abcd = read.csv(paste(PATH,"abcd_mean.csv",sep=""))
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
# estimate principal components of white matter tract metrics and put data together
PC=list()
pgrs_names=names(abcd.pgrs)[2:7]
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
  #d=data[[i]]%>%select(eid,age,sex,scanner) # in case the full data frame (containign all teh tract metrics) is used
  PC[[i]] = merge(PC[[i]],data[[i]],by="eid")
  # add pgrs pc
  psypc = PC[[i]][pgrs_names] %>% prcomp(scale. = T, center = T) # psychiatric disorder pc
  PC[[i]]$psypc1 = as.numeric(psypc$x[,1])
  #PC[[i]]$psypc2 = as.numeric(psypc$x[,2])
}
#
# Analyse ####
#
#
#
# PGRS sample differences (population versus imaging sample)
pgrs_names = (pgrs[[1]]%>%names)[2:8]
#ds=data.frame(matrix(ncol=length(PC),nrow = length(pgrs_names)))
ds=list()
for(l in 1:length(PC)){
  d=low=up=p=c()
  for (i in 1:length(pgrs_names)){
    d[i]=cohen.d(unlist(pgrs[[l]][pgrs_names][i]),unlist(PC[[l]][pgrs_names][i]))$estimate
    low[i]=as.numeric(cohen.d(unlist(pgrs[[l]][pgrs_names][i]),unlist(PC[[l]][pgrs_names][i]))$conf.int[1])
    up[i]=as.numeric(cohen.d(unlist(pgrs[[l]][pgrs_names][i]),unlist(PC[[l]][pgrs_names][i]))$conf.int[2])
    p[i]=t.test(unlist(pgrs[[l]][pgrs_names][i]),unlist(PC[[l]][pgrs_names][i]))$p.value
  }
  ds[[l]]=data.frame(pgrs_names,d,low,up,p)
}
ds=list.rbind(ds)
ds$data = c(replicate((nrow(ds))/2,"UKB"),replicate((nrow(ds))/2,"ABCD"))
write.csv(x = ds,file = paste(PATH,"Tables/diff_pgrs.csv",sep=""))
# mean(unlist(pgrs[[l]][pgrs_names][i])) # for ref of the direction of effects. This is the population.
# mean(unlist(PC[[l]][pgrs_names][i])) # This is the imaging sample
#
#
#
#
#
# Now, we start with looking at the associations between the PC of the psychiatric disorder | AD PGRS (y) ~ skeleton average/mean scores
pgrs_names=c("psypc1","AD")
skeleton_metrics=PC[[1]]%>%select(contains("_Mean"))%>%names
res=list() # empty list to be filled by loop
for (i in 1:length(PC)){
  res1=data.frame(matrix(ncol=length(skeleton_metrics), nrow=length(pgrs_names)))
  for (j in 1:length(skeleton_metrics)){
    res11 = c()
    for (k in 1:length(pgrs_names)){
      f1=formula(paste(pgrs_names[k],"~scale(",skeleton_metrics[j],")+age+sex+scanner+ethnicity",sep="")) # formula for skeleton mean associations
      m=lm(f1,data=PC[[i]])
      res11[k] = summary(m)$coefficients[2]
    }
    res1[,j] = res11
  }
  res[[i]] = res1
  names(res[[i]])=skeleton_metrics
  row.names(res[[i]])=pgrs_names
}
res = list.rbind(res)
names(res) = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", 
         "DTI-FA","DTI-MD","DTI-RD","DTI-AD",
         "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
res$Group = c("UKB-PsyPC","UKB-AD","ABCD-PsyPC","ABCD-AD")
res=melt(res)
names(res)=c("Group","Metric","Beta")

p1 = ggplot(res, aes(x=Group, y=Metric)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.04, -0.02, 0, 0.02, 0.04)) +
  geom_text(aes(label = round(Beta,3)))+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p1
ggsave(paste(PATH,"Figures/pgrs_skeleton.pdf",sep=""),p1,width=10,height=5)





# Second, we look at how the first five PCs of WM tracts associate with PGRS PCs
tracts = names(tracts[[1]]%>%select(-eid,-sex,-scanner,-age)) # use these names to then loop over data frames to get the first 5 PCs
#diffusion_components=names(PC[[1]])[2:7] # now, not needed anymore

res=list() # empty list to be filled by loop
for (i in 1:length(PC)){
  res1=data.frame(matrix(ncol=length(skeleton_metrics), nrow=length(pgrs_names)))
  for (j in 1:length(skeleton_metrics)){
    res11 = c()
    for (k in 1:length(pgrs_names)){
      f1=formula(paste(pgrs_names[k],"~scale(",skeleton_metrics[j],")+age+sex+scanner+ethnicity",sep="")) # formula for skeleton mean associations
      m=lm(f1,data=PC[[i]])
      res11[k] = summary(m)$coefficients[2]
    }
    res1[,j] = res11
  }
  res[[i]] = res1
}
res





# Plot #### (WORK IN PROGRESS)
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
