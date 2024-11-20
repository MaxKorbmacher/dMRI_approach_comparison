# wrangle abcd and ukb body anthropometrics
# 15 September 2024
#
# get height, weight, systolic, diastolic, and pulse pressure
# load pkgs
library(dplyr)
#
######################### ABCD
# load in data
anthro=read.csv("/cluster/projects/p33/groups/imaging/abcd/pheno/ABCDStudyNDA_5.1/core/physical-health/ph_y_anthro.csv")
eth=read.csv("/cluster/projects/p33/groups/imaging/abcd/pheno/ABCDStudyNDA_5.1/core/abcd-general/abcd_p_demo.csv")

# we are interested in anthroweightcalc = weight, anthroheightcalc = height
weihwi = anthro %>% filter(eventname=="baseline_year_1_arm_1") %>% select(src_subject_id, anthroweightcalc, anthroheightcalc)
#
# now, we get the blood pressure variables
anthro=read.csv("/cluster/projects/p33/groups/imaging/abcd/pheno/ABCDStudyNDA_5.1/core/physical-health/ph_y_bp.csv")
bp = anthro %>% select(src_subject_id, blood_pressure_sys_mean, blood_pressure_dia_mean, blood_pressure_pulse_mean)
anthro = merge(weihwi,bp,by="src_subject_id")
names(anthro) = c("eid", "weight","height","systolic","diastolic","PP")
# transform to metric system
anthro$height=anthro$height*2.54
anthro$weight = anthro$weight*0.45359237
anthro$BMI=anthro$weight/((anthro$height/100)^2)
anthro=anthro%>% select(-c(weight, height))
eth=eth %>% select(src_subject_id,race_ethnicity)
names(eth)=c("eid","ethnicity")
antrho=merge(eth,anthro,by="eid")
write.csv(anthro, "/cluster/projects/p33/users/maxk/wm_ica/abcd_anthro.csv", row.names = F)
#
#
#
######################### UKB
# load raw data
ukb=read.csv("/cluster/projects/p33/users/maxk/UKB/BP_BAG/data/BP.csv")
eth=read.csv("/cluster/projects/p33/users/maxk/UKB/data/ethnicity.txt.csv")
# note, ethnicity coding: https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=1001
eth=eth%>%select(eid,X21000.0.0) # least na in the first ethnicity field, hence selected
names(eth)=c("eid","ethnicity")
ukb$systolic1 = ukb$X4080.2.0
ukb$systolic2 = ukb$X4080.2.1
ukb$diastolic1 = ukb$X4079.2.0
ukb$diastolic2 = ukb$X4079.2.1
# estimate averages
ukb$systolic = (ukb$systolic1+ukb$systolic2)/2
ukb$diastolic = (ukb$diastolic1+ukb$diastolic2)/2
ukb$PP = ((ukb$systolic1-ukb$diastolic1)+(ukb$systolic2-ukb$diastolic2))/2
# BMI
ukb$BMI = ukb$X21001.2.0
# select only relevant vars
ukb = ukb %>% select(eid, systolic, diastolic, PP, BMI)
ukb=merge(ukb,eth,by="eid")
write.csv(ukb, "/cluster/projects/p33/users/maxk/wm_ica/ukb_anthro.csv", row.names = F)
