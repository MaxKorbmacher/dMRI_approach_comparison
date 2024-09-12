# White matter microstructure characteristics in ABCD and UKB
# Here, the focus is to examine age, sex, and asymmetry
# Max Korbmacher, 12 September 2024
#
# R version: 4.2.1
#
#
#
# DEFINE PATH FOR DATA AND OUTPUT
PATH = "/your/path/"
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
               viridis,cowplot,regclass,
               MuMIn,
               pls,plsVarSel, # pls 
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
#
# Functions (if not here, those appear below)
#
# standardize independent and dependent variables to obtain standardized (and comparable coefficients)
stand = function(dataframe){
  newdataframe = dataframe %>% select(-c(eid, sex, age, scanner)) %>% scale %>% data.frame()
  newdataframe = cbind(newdataframe, dataframe %>% select(c(eid, sex, age, scanner)))
  return(newdataframe)
}
#
# for common colouring scale
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
#
#
#
#
# 1) Descriptives and other results not directly answering the research questions ####
# 1.1) Power #### 
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
# u = 4 for metric, sex, scanner, age
# v = 5,147 (ABCD datasets) - 5
#
#
# create list of dfs
mean_list = list(stand(ukb_mean), stand(abcd_mean), stand(dat_mean))
#
# 1.2) Crude correlations between all skeleton averages ####
deps1 = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
          "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
          "DKI-MK", "DKI-RK", "DKI-AK","DTI-FA", "DTI-MD","DTI-RD","DTI-AD", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
          "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
cormat = function(df){
  cor1 = df%>%select(-c(eid,age,sex,scanner))
  #roworder = names(cor1)
  names(cor1) = deps1
  cor1 = cor(cor1)
  return(cor1)
}
a = as.ggplot(pheatmap(cormat(mean_list[[1]]),cluster_rows = F, cluster_cols = T))
b = as.ggplot(pheatmap(cormat(mean_list[[2]]),cluster_rows = F, cluster_cols = T))
c = as.ggplot(pheatmap(cormat(mean_list[[3]]),cluster_rows = F, cluster_cols = T))
d = as.ggplot(pheatmap(cormat(mean_list[[1]]),cluster_rows = T, cluster_cols = T))
e = as.ggplot(pheatmap(cormat(mean_list[[2]]),cluster_rows = T, cluster_cols = T))
f = as.ggplot(pheatmap(cormat(mean_list[[3]]),cluster_rows = T, cluster_cols = T))
plot3 = a+b+c+d+e+f
# note: top row = unclustered ukb,abcd,both; bottom row = clustered
ggsave(paste(PATH,"Figures/skeleton_correlations.pdf",sep=""),plot3,width=20,height=14)
#
#
#
#
#
#
#
#
#
#
# 2) Sex classification ####
# use logistic regression for sex classification and then compare model performance across diffusion approaches
# 
# first, make separate data frames for each approach
sepdat = function(data){
  data = stand(data)
  # WMTI
  wmti = data %>% select(starts_with("axEAD"), starts_with("radEAD"), c(age, sex, scanner))
  # BRIA
  bria = data %>% select(starts_with("v_"), starts_with("micro_"), starts_with("Drad"), starts_with("Dax"), c(age, sex, scanner))
  # DKI
  dki = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk"), c(age, sex, scanner))
  # DTI
  dti = data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa"), c(age, sex, scanner)) 
  # SMT
  smt =  data %>% select(starts_with("smt_md"), starts_with("smt_long"), c(age, sex, scanner))
  # SMTmc
  smtmc =  data %>% select(starts_with("smt_mc"), c(age, sex, scanner))
  # all
  All = data%>% select(-eid)
  dflist = list(bria, dki,dti,smt,smtmc,wmti,All)
  return(dflist)
}
# write data and empty list
abcd_list = sepdat(abcd)
ukb_list = sepdat(ukb)
pcplots=pcplots2=mod=mod2=res1=res2=tmp_dat1=tmp_dat2=list()
tmp=tmp2=data.frame(matrix(ncol=3, nrow = length(abcd_list)))
# make accuracy function
acc = function(fitted.model,data){
  f1 = factor(round(predict(fitted.model, type="response")))
  levels(f1)= levels=c("female","male")
  cm=as.matrix(table(Actual = data$sex, Predicted=f1))
  Accuracy = data.frame(t(diag(cm)/apply(cm, 2, sum))) # Per sex/class accuracy
  Accuracy$Overall = sum(diag(cm))/sum(cm) # diag(cm) = number of correctly classified instances per class
  return(Accuracy)
}
app = c("BRIA", "DKI","DTI","SMT","SMTmc","WMTI", "All")
for (i in 1:length(abcd_list)){
  test = abcd_list[[i]] %>% select(-c(sex,age,scanner)) %>%prcomp(scale. = T, center = T) # scaling numerical data
  tmp_dat1[[i]] = data.frame(test$x[,1:5], abcd_list[[i]] %>%select(c(sex,age,scanner))) # putting data together
  pcplots[[i]]=fviz_eig(test, main = app[i], addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100) # scree plots
  tmp_dat1[[i]]$age=scale(tmp_dat1[[i]]$age)
  mod[[i]] = glm(factor(sex) ~ PC1+PC2+PC3+PC4+PC5+age+scanner, data = tmp_dat1[[i]], family = "binomial") # log regression
  res1[[i]] = data.frame(summary(mod[[i]])$coef) # save coefficients, z, p, etc.
  tmp[i,]=acc(mod[[i]],abcd_list[[i]]) # make accuracy table
  test = ukb_list[[i]] %>% select(-c(sex,age,scanner)) %>%prcomp(scale. = T, center = T) # repeat steps from above
  pcplots2[[i]]=fviz_eig(test, main = app[i], addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  tmp_dat2[[i]] = data.frame(test$x[,1:5], ukb_list[[i]] %>%select(c(sex,age,scanner)))
  tmp_dat2[[i]]$age=scale(tmp_dat2[[i]]$age)
  mod2[[i]] = glm(factor(sex) ~ PC1+PC2+PC3+PC4+PC5+age+scanner, data = tmp_dat2[[i]], family = "binomial")
  res2[[i]] = data.frame(summary(mod2[[i]])$coef)
  tmp2[i,]=acc(mod2[[i]],ukb_list[[i]])
}
abcdscree=ggarrange(plotlist = pcplots)
ukbscree=ggarrange(plotlist = pcplots2)
ggsave(paste(PATH,"Figures/ABCD_scree.pdf",sep=""),abcdscree,width=17,height=9)
ggsave(paste(PATH,"Figures/UKB_scree.pdf",sep=""),ukbscree,width=17,height=9)

# 2.1 Model performance plotting ####
names(tmp)=names(tmp2)=c("Female","Male","Overall")
tmp$Approach = tmp2$Approch = app
tmp_list=list(melt(tmp),melt(tmp2))
names(tmp_list[[1]])=names(tmp_list[[2]])=c("Approach","Group", "Specificity")
# plot
p1 = ggplot(tmp_list[[1]], aes(y=Approach, x=Group)) +
  geom_tile(colour="black", size=0.25, aes(fill=Specificity)) +
  #scale_y_continuous(limits = c(0,1),breaks = seq(.6,.9,.1))+
  #scale_fill_gradient(high="#880700", low="white")+
  #scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
  #                     breaks = c(0, 0.7,.75,.8)) +
  geom_text(aes(label=round(Specificity,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("") +labs(title="ABCD")
p2 = ggplot(tmp_list[[2]], aes(y=Approach, x=Group)) +
  geom_tile(colour="black", size=0.25, aes(fill=Specificity)) +
  #scale_y_continuous(limits = c(0,1),breaks = seq(.6,.9,.1))+
  #scale_fill_gradient(high="#880700", low="white")+
  #scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white",
  #                     breaks = c(0, 0.7,.75,.8)) +
  geom_text(aes(label=round(Specificity,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("") +labs(title="UKB")
# unify the colouring scale and finish off
set_scale_union(p1,p2,scale=scale_fill_gradient(high="#880700", low="white"))
sex_class = ggarrange(p1,p2,common.legend = T, legend = "right")
ggsave(paste(PATH,"Figures/sex_plot_tracts.pdf",sep=""),sex_class,width=5,height=4)
rm(p1,p2,tmp,tmp2,mod,mod2)
# 2.2 Predictors ####
# Note that the obtained coefficients are odds ratios.
# It represents the ratio of the odds that an event will occur (event = 1) given 
# the presence of the predictor x (x = 1), compared to the odds of the event occurring in the absence of that predictor (x = 0).
# In other words, 1 = even, smaller than 1 is female, bigger than 1 male.
for (i in 1:length(app)){
  res1[[i]]$Approach = app[i]
  res2[[i]]$Approach = app[i]
}
res1 = list.rbind(res1)
res2 = list.rbind(res2)
res2$Predictor=rownames(res2)
res1$Predictor = rownames(res1)
res1=res1%>%filter(grepl("PC",Predictor))
res1$Predictor=c("PC1","PC2","PC3","PC4","PC5")
res2=res2%>%filter(grepl("PC",Predictor))
res2$Predictor=c("PC1","PC2","PC3","PC4","PC5")
res1$Approach=res2$Approach=c(replicate(5,app[1]),replicate(5,app[2]),
                              replicate(5,app[3]),replicate(5,app[4]),
                              replicate(5,app[5]), replicate(5,app[6]),
                              replicate(5,app[7]))
res1$Beta=res1$Estimate
res2$Beta=res2$Estimate
p1=ggplot(res1, aes(x=Approach, y=Predictor)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white") +
  geom_text(aes(label=round(Beta,2)),size=4)+ labs(title="ABCD")+ 
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2=ggplot(res2, aes(x=Approach, y=Predictor)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white") +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")+
  labs(title="UKB")
set_scale_union(p1,p2,scale=scale_fill_gradient2(high="#880700", low="#3a81b5", 
                                                 midpoint = 0, mid = "white"))
sex_pc=ggarrange(p1,p2,common.legend = T,legend="right")
ggsave(paste(PATH,"Figures/sex_pc_tract.pdf",sep=""),sex_pc,width=10,height=4)
#
# combine specificity and PC sex figures
sex_class=annotate_figure(sex_class,top = "Prediction Accuracy")
sex_pc=annotate_figure(sex_pc,top = "Principal Components' Contributions to Predictions")
sex_fig=ggarrange(sex_class,sex_pc)
ggsave(paste(PATH,"Figures/sex_accuracy_and_pc_tract.pdf",sep=""),sex_fig,width=17,height=4)
#write.csv(res1, paste(PATH,"Tables/sex_predictors_abcd.csv",sep=""))
#write.csv(res2, paste(PATH,"Tables/sex_predictors_ukb.csv",sep=""))
# correlate each predictor with each of the first 5 principal components
PCnames=c("PC1","PC2","PC3","PC4","PC5")
dflist=mod=mod2=list()
for (i in 1:length(tmp_dat1)){
  variable_names = names(abcd_list[[i]]%>% select(-c(age,sex,scanner)))
  df=as.data.frame(matrix(ncol=10,nrow=length(variable_names)))
  for (components in 1:length(PCnames)){
    for (predictors in 1:length(variable_names)){
      f1=formula(paste(PCnames[components],"~",variable_names[predictors],"+age+sex+(scanner)",sep=""))
      dat1=cbind(tmp_dat1[[i]],abcd_list[[i]]%>% select(-c(age,sex,scanner)))
      dat2=cbind(tmp_dat2[[i]],ukb_list[[i]]%>% select(-c(age,sex,scanner)))
      mod = lm(f1,data=dat1)
      mod2 = lm(f1,data=dat2)
      df[predictors,]=c(summary(mod)$coefficients[2],summary(mod)$coefficients[3],
                        summary(mod)$coefficients[4],summary(mod)$coefficients[5],
                        summary(mod)$coefficients[6],summary(mod2)$coefficients[2],
                        summary(mod2)$coefficients[3],summary(mod2)$coefficients[4],
                        summary(mod2)$coefficients[5],summary(mod2)$coefficients[6])
    }
  }
  dflist[[i]]=df
}
# make 6 plots with one abcd panel and one ukb panel each
for (i in 1:length(dflist)){ # first melt the frames and name for ggplot readibility
  names(dflist[[i]]) = c("ABCD_PC1","ABCD_PC2","ABCD_PC3","ABCD_PC4","ABCD_PC5", "UKB_PC1", "UKB_PC2","UKB_PC3", "UKB_PC4", "UKB_PC5")
  dflist[[i]] = melt(dflist[[i]])
  dflist[[i]]$Metric = names(abcd_list[[i]]%>% select(-c(age,sex,scanner)))
  names(dflist[[i]]) = c("var","Beta", "Metric")
  dflist[[i]]$var = sub("_"," ",dflist[[i]]$var)
}
# then plot into a plotlist
plots=list()
for (i in 1:length(dflist)){
  plots[[i]]=ggplot(dflist[[i]], aes(x=var, y=Metric)) +
    geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
    scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white")+
    #breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
    geom_text(aes(label=round(Beta,2)),size=2)+
    theme(axis.text.y = element_text(size = 4)) + theme_bw() + xlab("") + ylab("")+
    theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=6))
}
# save each of the six individual diffusion apprach plots + 1 plot for the model combining all the approaches
for (i in 1:length(plots)){
  ggsave(file = paste(PATH,"Figures/Suppl_plot_",2+i,".pdf",sep=""),plots[[i]],width=12,height=14)
}
# redo the last (massive) plot
# split up in several smaller figures
rework=dflist[[i]] %>% filter(var=="ABCD PC1")
rework2=dflist[[i]] %>% filter(var=="UKB PC1")
# make an own function
plotfunc = function(data){
  p=ggplot(data, aes(x=var, y=Metric)) +
    geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
    scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white")+
    #breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
    geom_text(aes(label=round(Beta,2)),size=2)+
    theme(axis.text.y = element_text(size = 4)) + theme_bw() + xlab("") + ylab("")+
    theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=6))
  return(p)
}
p1=plotfunc(head(rework[order(rework$Beta, decreasing= T),], n = 50))
p2=plotfunc(head(rework[order(rework$Beta, decreasing= F),], n = 50))
p3=plotfunc(head(rework2[order(rework2$Beta, decreasing= T),], n = 50))
p4=plotfunc(head(rework2[order(rework2$Beta, decreasing= F),], n = 50))
p5=ggarrange(p1,p2,p3,p4,nrow=1)
ggsave(file = paste(PATH,"Figures/Suppl_plot_",9,".pdf",sep=""),p5,width=14,height=5)
#
#
#
#
#
#
#
# 3) Age prediction ####
# 3.1) Global ####
#
# we use simple linear models to then being able to compare the obtained coefficients and variance explained
# This is done for each white matter approach
preds=names(dat_mean[5:length(dat_mean)]) # create vector of variables to be predicted
res = list()
mean_list=list(ukb_mean,abcd_mean)
Beta=SE=c()
for (df in 1:length(mean_list)){
  for (i in 1:length(preds)){
    tmp_dat = stand(mean_list[[df]])
    tmp_dat$age = scale(tmp_dat$age)
    f1 = formula(paste("age~",preds[i],"+sex+scanner",sep="")) # regular linear models used due to singular fit 
    m1 = lm(f1, data=tmp_dat)
    Beta[i] = summary(m1)$coefficients[2]
    SE[i] = summary(m1)$coefficients[2,2] # not necessary for the plotting
  }
  res[[df]] = data.frame(preds, Beta, SE)
}
res[[1]]$Data = "UKB"
res[[2]]$Data = "ABCD"
plot_dat = list.rbind(res) # merge list/data frames
deps2 = c("BRIA-Vintra", "BRIA-vextra", "BRIA-vCSF", "BRIA-microRD", "BRIA-microFA", 
          "BRIA-microAX", "BRIA-microADC", "BRIA-DRADextra", "BRIA-DAXintra", "BRIA-DAXextra",
          "DKI-MK", "DKI-RK", "DKI-AK","DTI-FA", "DTI-MD","DTI-RD","DTI-AD", "SMT-long", "SMT-MD", "SMTmc-intra", "SMTmc-extraMD", "SMTmc-extratrans",
          "SMTmc-Diff", "WMTI-axEAD", "WMTI-AWF", "WMTI-radEAD")
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
set_scale_union(p1,p2,p3,
                scale=scale_fill_gradient2(
                  high="#880700", low="#3a81b5", 
                  midpoint = 0,mid = "white",
                  breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
age_asso=ggarrange(p1,p2,p3,nrow=1, align="v", common.legend = T, legend = "bottom")
ggsave(paste(PATH,"Figures/age_plot_skeleton.pdf",sep=""),age_asso,width=8,height=4)
#plot_grid(p1,p2,p3,nrow=1, align="v")
# 3.2 Regional ####
# Due to multicollinearity, we first reduce the dimensions for each approach and then all together using PCA
# the two lists to be used (each contained data frame represents a WMM approach): abcd_list, ukb_list
abcd.r2=ukb.r2=c()
model1=model2=list()
for (i in 1:length(abcd_list)){
  #test = abcd_list[[i]] %>% select(-c(sex,age,scanner)) %>%prcomp(scale. = T, center = T)
  #tmp_dat1[[i]] = data.frame(test$x[,1:5], abcd_list[[i]] %>%select(c(sex,age,scanner)))
  model1[[i]]=lm(age~PC1+PC2+PC3+PC4+PC5+sex+scanner, tmp_dat1[[i]])
  abcd.r2[i]=summary(model1[[i]])$r.squared
  #test = ukb_list[[i]] %>% select(-c(sex,age,scanner)) %>%prcomp(scale. = T, center = T)
  #tmp_dat2[[i]] = data.frame(test$x[,1:5], ukb_list[[i]] %>%select(c(sex,age,scanner)))
  model2[[i]]=lm(age~PC1+PC2+PC3+PC4+PC5+sex+scanner, tmp_dat2[[i]])
  ukb.r2[i]=summary(model2[[i]])$r.squared
}
#
#
# 3.2.1 Model performance ####
R2 = data.frame(app,abcd.r2,ukb.r2)
names(R2)=c("Approach","ABCD","UKB")
R2=melt(R2)
names(R2)=c("Approach","Data","R2")
pca_age=ggplot(R2, aes(x=Data, y=Approach)) +
  geom_tile(colour="black", size=0.25, aes(fill=R2)) +
  scale_fill_gradient(high="#880700", low="white") +
  geom_text(aes(label=round(R2,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + 
  xlab("") + ylab("")
pca_age$labels$fill = "Variance Explained"
ggsave(paste(PATH,"Figures/age_plot_R2_tract.pdf",sep=""),pca_age,width=5,height=4)
# 3.2.2 Predictors ####
# we can use the same inference as for the sex predictions, plotting the contr. of the five PCs
coefs1=coefs2=list()
for (i in 1:length(abcd_list)){
  coefs1[[i]] = data.frame(summary(model1[[i]])$coefficients)
  coefs2[[i]] = data.frame(summary(model2[[i]])$coefficients)
  coefs1[[i]]$Approach=app[i]
  coefs2[[i]]$Approach=app[i]
}

res1 = list.rbind(coefs1)
res2 = list.rbind(coefs2)
res2$Predictor=rownames(res2)
res1$Predictor = rownames(res1)
res1=res1%>%filter(grepl("PC",Predictor))
res1$Predictor=c("PC1","PC2","PC3","PC4","PC5")
res2=res2%>%filter(grepl("PC",Predictor))
res2$Predictor=c("PC1","PC2","PC3","PC4","PC5")
res1$Approach=res2$Approach=c(replicate(5,app[1]),replicate(5,app[2]),
                              replicate(5,app[3]),replicate(5,app[4]),
                              replicate(5,app[5]), replicate(5,app[6]),
                              replicate(5,app[7]))
res1$Beta=res1$Estimate
res2$Beta=res2$Estimate
p1=ggplot(res1, aes(x=Approach, y=Predictor)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white") +
  geom_text(aes(label=round(Beta,2)),size=4)+ labs(title="ABCD")+ 
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")
p2=ggplot(res2, aes(x=Approach, y=Predictor)) +
  geom_tile(colour="black", size=0.25, aes(fill=Beta)) +
  scale_fill_gradient2(high="#880700", low="#3a81b5", midpoint = 0, mid = "white") +
  geom_text(aes(label=round(Beta,2)),size=4)+
  theme(axis.text.y = element_text(size = 8)) + theme_bw() + xlab("") + ylab("")+
  labs(title="UKB")
set_scale_union(p1,p2,scale=scale_fill_gradient2(high="#880700", low="#3a81b5", 
                                                 midpoint = 0, mid = "white"))
age_pc=ggarrange(p1,p2,common.legend = T,legend="right")
ggsave(paste(PATH,"Figures/age_pc_tract.pdf",sep=""),age_pc,width=10,height=4)
# combined age plot
pca_age=annotate_figure(pca_age,top = "Variance in age explained by diffusion approaches")
age_pc=annotate_figure(age_pc,top = "White matter tract scalars' principal components' contributions to age preditions")
age_plot=ggarrange(pca_age, age_pc)
ggsave(paste(PATH,"Figures/age_plot_tract.pdf",sep=""),age_plot,width=17,height=4)
#
#
#
# 4) Asymmetry assessment ####
# Estimate tract-level asymmetry & plot them
# 5) Anthropometrics and blood pressure ####
# 6) Polygenic risk ####
