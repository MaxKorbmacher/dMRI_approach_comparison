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
# make some functions
simple_range_extracter <- function(p, scale) {
  d <- ggplot_build(p)
  d$plot$scales$get_scales(scale$aesthetics)$range$range
}
get_shared_scale <- function(..., scale) {
  plots <- list(...)
  ranges <- purrr::map(plots, ~simple_range_extracter(., scale))
  single_range <- range(unlist(ranges))
  scale$limits <- single_range
  scale
}
# Main function
set_scale_union <- function(..., scale) {
  exprs <- rlang::enexprs(...)
  scale <- get_shared_scale(..., scale = scale)
  var_nms <- purrr::map_chr(exprs, rlang::as_name)
  edit_plots_in_place(var_nms, env = parent.frame(),
                      scale = scale)
  # Invisibly return the scale, in case you need it later
  invisible(scale)
}

# Sub-function
edit_plots_in_place <- function(names, env, scale) {
  vars <- rlang::env_has(env = env, nms = names)
  if (!all(vars))
    stop("Environment does not have variables for ",
         paste(names(vars[!vars]), collapse=", "))
  purrr:::walk(names, function(nm) {
    og_plot <- rlang::env_get(env, nm = nm)
    message("Changing plot `", nm, "`")
    # Muffles messages about already having scales
    withCallingHandlers(
      assign(x = nm, envir = env,
             value = og_plot + scale),
      message = function(err) {
        if (grepl("already present", err$message))
          invokeRestart("muffleMessage")
      })
  })
}
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
# PGRS PC
pgrs_names = (pgrs[[1]]%>%names)[2:7]
pgrs_comp = pgrs[[1]][pgrs_names]%>% prcomp(scale. = T, center = T)
fviz_eig(pgrs_comp, main = "UKB Psych PC", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
pgrs_comp = pgrs[[2]][pgrs_names]%>% prcomp(scale. = T, center = T)
fviz_eig(pgrs_comp, main = "ABCD Psych PC", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
#
#
#
#
# Now, we start with looking at the associations between the PGRS (y) ~ skeleton average/mean scores
pgrs_names=(pgrs[[1]]%>%names)[2:length(pgrs[[1]]%>%names)]
skeleton_metrics=PC[[1]]%>%select(contains("_Mean"))%>%names
res=res02=list() # empty list to be filled by loop
for (i in 1:length(PC)){
  res1=res2=data.frame(matrix(ncol=length(skeleton_metrics), nrow=length(pgrs_names)))
  for (j in 1:length(skeleton_metrics)){
    res22 = res11 = c()
    for (k in 1:length(pgrs_names)){
      f1=formula(paste(pgrs_names[k],"~scale(",skeleton_metrics[j],")+age+sex+scanner+ethnicity",sep="")) # formula for skeleton mean associations
      m=lm(f1,data=PC[[i]])
      res11[k] = summary(m)$coefficients[2]
      res22[k] =summary(m)$coefficients[2,4]
    }
    res1[,j] = res11
    res2[,j] = res22
  }
  res[[i]] = res1 # beta coefficients
  res02[[i]] = res2 # p-values
  names(res02[[i]])=names(res[[i]])=skeleton_metrics
  row.names(res02[[i]])=row.names(res[[i]])=pgrs_names
}
res = list.rbind(res) # beta coefficients
res02 = list.rbind(res02) # p-values
names(res02) = names(res) = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
         "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
         "DKI-MK", "DKI-RK", "DKI-AK", 
         "DTI-FA","DTI-MD","DTI-RD","DTI-AD",
         "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
         "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
res02$Group = res$Group = c(paste("UKB-",rownames(res02)[1:(length(pgrs_names))],sep=""),paste("ABCD-",rownames(res02)[1:(length(pgrs_names))],sep=""))# c("UKB-PsyPC","UKB-AD","ABCD-PsyPC","ABCD-AD")
res=melt(res)
res02=melt(res02)
res02$cor_p = res02$value*708
res02%>%filter(cor_p<.05) # show the metrics that survived multiple comparison
names(res)=c("Group","Metric","Beta")
p1 = ggplot(res, aes(x=Group, y=Metric)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.02,-0.01, 0,0.01, 0.02)) +
  geom_text(aes(label = round(Beta,3)))+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
ggsave(paste(PATH,"Figures/pgrs_skeleton.pdf",sep=""),p1,width=14,height=8)
#
#
#
#
#
# Second, we look at how the first five PCs of WM tracts associate with PGRS PCs
#tracts = names(tracts[[1]]%>%select(-eid,-sex,-scanner,-age)) # use these names to then loop over data frames to get the first 5 PCs
diffusion_components=names(PC[[1]])[2:7]
res02=res=list() # empty list to be filled by loop
for (i in 1:length(PC)){
  res2=res1=data.frame(matrix(ncol=length(diffusion_components), nrow=length(pgrs_names)))
  for (j in 1:length(diffusion_components)){
    res11 = res22 = c()
    for (k in 1:length(pgrs_names)){
      f1=formula(paste(pgrs_names[k],"~scale(",diffusion_components[j],")+age+sex+scanner+ethnicity",sep="")) # formula for skeleton mean associations
      m=lm(f1,data=PC[[i]])
      res11[k] = summary(m)$coefficients[2]
      res22[k] =summary(m)$coefficients[2,4]
    }
    res1[,j] = res11
    res2[,j] = res22
  }
  names(res2)=names(res1)=diffusion_components
  res2$PGRS = res1$PGRS = pgrs_names
  res[[i]] = res1
  res02[[i]] = res2
}
#res # 1) UKB, 2) ABCD
r1=melt(res[[1]])
r2=melt(res[[2]])
names(r1) = names(r2) = c("PGRS","variable","Beta")
p1=ggplot((r1), aes(x=variable, y=PGRS)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.02, -0.015, -0.01,-0.005, 0)) +
  geom_text(aes(label = round(Beta,3)))+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p1 = annotate_figure(p1, top = text_grob("UKB", rot = 0))
p2=ggplot((r2), aes(x=variable, y=PGRS)) + 
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
                       breaks = c(-0.0075,-0.005,-0.0025, 0, 0.0025, 0.005,0.0075)) +
  geom_text(aes(label = round(Beta,3)))+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2 = annotate_figure(p2, top = text_grob("ABCD", rot = 0))
p3 = ggarrange(p1,p2,common.legend = F)
ggsave(paste(PATH,"Figures/PGRS_WM_PC.pdf",sep=""),p3,width=12,height=4)
#
# check sig
melt(res02[[1]]) %>% filter(value<.05/708)
melt(res02[[2]]) %>% filter(value<.05/708)
