library(dplyr)
library(lsr)
library(psych)
library(cluster)
library(rstudioapi)
library(fpc)
library(factoextra)
library(psych)
library(clValid)
library(cowplot)
library(reshape2)
library(ggpubr)
library(forcats)
library(ggplot2)
library(arules)

#set working directory to location of script--not needed when called 
#through Rmarkdown doc. Uncomment below if running locally/manually
#this.dir = getActiveDocumentContext()$path
#setwd(dirname(this.dir))

source('functions.R')
source('ARM_ORconf.R')
LOC_ndup = read.csv('Data/LOC_noDuplicates.csv')

#### Demographics ####
#get crosstab table for LOC and count number of NA overall
#and for differences between LOC and no LOC
#note: sd.function.na and range.function.na all come from
#the functions.R file so that must be sourced above to work

## LOC
#get crosstab table for LOC and count number of NA
loc_tab = xtabs(~loc1, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$loc1), ])

## Age
age_mean = mean(LOC_ndup$cAge_mo, na.rm = TRUE)/12
age_sd = sd(LOC_ndup$cAge_mo, na.rm = TRUE)/12
age_range = range(LOC_ndup$cAge_mo, na.rm = TRUE)/12
nrow(LOC_ndup[is.na(LOC_ndup$cAge_mo), ])

age_loc_t = t.test((cAge_mo/12)~loc1, data = LOC_ndup)
age_loc_sd = sd.function.na(LOC_ndup, LOC_ndup$cAge_mo, LOC_ndup$loc1)/12
age_loc_range = range.function.na(LOC_ndup, (LOC_ndup$cAge_mo/12), LOC_ndup$loc1)
age_loc_d = cohensD(LOC_ndup[LOC_ndup$loc1 == 'Yes', ]$cAge_mo, LOC_ndup[LOC_ndup$loc1 == 'No', ]$cAge_mo)

nrow(LOC_ndup[is.na(LOC_ndup$cAge_mo) && LOC_ndup$loc1 == 'Yes', ])
nrow(LOC_ndup[is.na(LOC_ndup$cAge_mo) && LOC_ndup$loc1 == 'No', ])
## Gender
gender_tab = xtabs(~sex, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$sex), ])

gender_loc_tab = xtabs(~sex + loc1, data=LOC_ndup)
gender_loc_chi = chisq.test(gender_loc_tab)
nrow(LOC_ndup[which(is.na(LOC_ndup$sex) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$sex) && LOC_ndup$loc1 == 'No'), ])

## cBMIp
cBMIp_mean = mean(LOC_ndup$cBodyMass_p, na.rm = TRUE)
cBMIp_sd = sd(LOC_ndup$cBodyMass_p, na.rm = TRUE)
cBMIp_range = range(LOC_ndup$cBodyMass_p, na.rm = TRUE)
nrow(LOC_ndup[is.na(LOC_ndup$cBodyMass_p), ])

cBMIp_loc_t = t.test(cBodyMass_p~loc1, data = LOC_ndup)
cBMIp_loc_sd = sd.function.na(LOC_ndup, LOC_ndup$cBodyMass_p, LOC_ndup$loc1)
cBMIp_loc_range = range.function.na(LOC_ndup, (LOC_ndup$cBodyMass_p), LOC_ndup$loc1)
cBMIp_loc_d = cohensD(LOC_ndup[LOC_ndup$loc1 == 'Yes', ]$cBodyMass_p, LOC_ndup[LOC_ndup$loc1 == 'No', ]$cBodyMass_p)

nrow(LOC_ndup[is.na(LOC_ndup$cBodyMass_p) && LOC_ndup$loc1 == 'Yes', ])
nrow(LOC_ndup[is.na(LOC_ndup$cBodyMass_p) && LOC_ndup$loc1 == 'No', ])

## BMI class
BMI_tab = xtabs(~cBodyMass_class, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$cBodyMass_class), ])

BMI_loc_tab = xtabs(~cBodyMass_class + loc1, data = LOC_ndup)
BMI_loc_fisher = fisher.test(BMI_loc_tab)
BMInoUW_chi = chisq.test(matrix(c(93, 22, 13, 8, 14, 7), byrow = TRUE, ncol = 2))
BMInoUW_loc_fisher = fisher.test(matrix(c(93, 22, 13, 8, 14, 7), byrow = TRUE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$cBodyMass_class) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$cBodyMass_class) & LOC_ndup$loc1 == 'No'), ])

##objective measurement of parent height/weight in FBS, EBS, DMK, and Test-Retest studies
HWmeasure = nrow(LOC_ndup[grep("DMK*|Test*|FBS*|EBS*", LOC_ndup$StudyID), ])
parent_HWmeasure = xtabs(~relationship, data = LOC_ndup[grep("DMK*|Test*|FBS*|EBS*", LOC_ndup$StudyID), ])

## Race
race_tab = xtabs(~cRace, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$cRace), ])

race_loc_tab = xtabs(~cRace + loc1, data = LOC_ndup)
race_loc_fisher = fisher.test(matrix(c(4,31,1,2,114,4), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$cRace) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$cRace) & LOC_ndup$loc1 == 'No'), ])

## Ethnicity
ethnicity_tab = xtabs(~cEthnicity, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$cEthnicity), ])

ethnicity_loc_tab = xtabs(~cEthnicity + loc1, data = LOC_ndup)
ethnicity_loc_fisher = fisher.test(matrix(c(2,33,5,76), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$cEthnicity) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$cEthnicity) & LOC_ndup$loc1 == 'No'), ])

## Income
income_tab = xtabs(~income, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$income), ])

income_loc_tab = xtabs(~income + loc1, data = LOC_ndup)
income_loc_fisher = fisher.test(matrix(c(6,21,10,40,54,24), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$income) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$income) & LOC_ndup$loc1 == 'No'), ])

## mEducation
mEducation_tab = xtabs(~mEducation, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$mEducation), ])

mEducation_loc_tab = xtabs(~mEducation + loc1, data = LOC_ndup)
mEducation_loc_fisher = fisher.test(matrix(c(10,26,11,105), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$mEducation) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$mEducation) & LOC_ndup$loc1 == 'No'), ])

## dEducation
dEducation_tab = xtabs(~dEducation, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$dEducation), ])

dEducation_loc_tab = xtabs(~dEducation + loc1, data = LOC_ndup)
dEducation_loc_fisher = fisher.test(matrix(c(9, 25, 23, 92), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$dEducation) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$dEducation) & LOC_ndup$loc1 == 'No'), ])

## breast fed total
BreastFed_tab = xtabs(~BreastFed, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$BreastFed), ])

BreastFed_loc_tab = xtabs(~BreastFed + loc1, data = LOC_ndup)
BreastFed_loc_chi = chisq.test(BreastFed_loc_tab)
nrow(LOC_ndup[which(is.na(LOC_ndup$BreastFed) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$BreastFed) & LOC_ndup$loc1 == 'No'), ])

## breast fed mo
BreastFed_mo_tab = xtabs(~BreastFed_mo, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$BreastFed_mo), ])

BreastFed_mo_loc_tab = xtabs(~BreastFed_mo + loc1, data = LOC_ndup)
BreastFed_mo_loc_fisher = fisher.test(matrix(c(10,3,13,5,5,30,12,35,31,11), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$BreastFed_mo) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$BreastFed_mo) & LOC_ndup$loc1 == 'No'), ])

## mBMI
mBMI_tab = xtabs(~mBodyMass_class, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$mBodyMass_class), ])

mBMI_loc_tab = xtabs(~mBodyMass_class + loc1, data = LOC_ndup)
mBMI_loc_fisher = fisher.test(matrix(c(8,4,5,15,2,11,13,19,57,2), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$mBodyMass_class) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$mBodyMass_class) & LOC_ndup$loc1 == 'No'), ])

## dBMI
dBMI_tab = xtabs(~dBodyMass_class, data = LOC_ndup)
nrow(LOC_ndup[is.na(LOC_ndup$dBodyMass_class), ])

dBMI_loc_tab = xtabs(~dBodyMass_class + loc1, data = LOC_ndup)
dBMI_loc_fisher = fisher.test(matrix(c(5,7,15,9,0,12,22,46,34,0), byrow = FALSE, ncol = 2))
dBMI_loc_fisher_noUW = fisher.test(matrix(c(5,7,15,9,12,22,46,34), byrow = FALSE, ncol = 2))
nrow(LOC_ndup[which(is.na(LOC_ndup$dBodyMass_class) & LOC_ndup$loc1 == 'Yes'), ])
nrow(LOC_ndup[which(is.na(LOC_ndup$dBodyMass_class) & LOC_ndup$loc1 == 'No'), ])

##### CEBQ ####

#recode CEBQ data as numeric
#create dataset to add numeric CEBQ values
cebq_LOC_ndup = LOC_ndup[c(2, 52:86)]
#for the 35 cebq questions
for(c in 2:36){
  #create new column
  cebq_LOC_ndup[[c+35]] = NA
  
  #get original varname
  varname = names(cebq_LOC_ndup)[[c]]
  
  #add "_raw" to varname and set as the name for the new column
  names(cebq_LOC_ndup)[[c+35]] = paste(varname, "raw", sep = "_")
  
  #set that needs to be reverse coded
  if(c == 4|c == 5|c == 11|c == 17|c == 33){
    cebq_LOC_ndup[c+35] = ifelse(is.na(cebq_LOC_ndup[[c]]), NA, ifelse(
      cebq_LOC_ndup[[c]] == "Never", 5, ifelse(
        cebq_LOC_ndup[[c]] == "Rarely", 4, ifelse(
          cebq_LOC_ndup[[c]] == "Sometimes", 3, ifelse(
            cebq_LOC_ndup[[c]] == "Often", 2, 1
          )
        )
      )
    )) 
  }
  #non-reverse coded
  else {
    cebq_LOC_ndup[c+35] = ifelse(is.na(cebq_LOC_ndup[[c]]), NA, ifelse(
      cebq_LOC_ndup[[c]] == "Never", 1, ifelse(
        cebq_LOC_ndup[[c]] == "Rarely", 2, ifelse(
          cebq_LOC_ndup[[c]] == "Sometimes", 3, ifelse(
            cebq_LOC_ndup[[c]] == "Often", 4, 5
          )
        )
      )
    )) 
  }  
}
#add loc
cebq_LOC_ndup = merge(cebq_LOC_ndup, LOC_ndup[c(2, 33)], by = 'StudyID')

#get means and descriptives for CEBQ questionnaire
## FR
#numeric mean for subscale for each person as new var
cebq_LOC_ndup$cebqFR = rowMeans(data.frame(cebq_LOC_ndup$cebq12_raw, cebq_LOC_ndup$cebq14_raw, cebq_LOC_ndup$cebq19_raw, cebq_LOC_ndup$cebq28_raw, cebq_LOC_ndup$cebq34_raw), na.rm = T)

#internal reliability
alpha_cebqFR = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq12_raw, cebq_LOC_ndup$cebq14_raw, cebq_LOC_ndup$cebq19_raw, cebq_LOC_ndup$cebq28_raw, cebq_LOC_ndup$cebq34_raw)))

cebqFR_mean = mean(cebq_LOC_ndup$cebqFR, na.rm = TRUE)
cebqFR_sd = sd(cebq_LOC_ndup$cebqFR, na.rm = TRUE)
cebqFR_range = range(cebq_LOC_ndup$cebqFR, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqFR), ])

cebqFR_loc_t = t.test(cebqFR~loc1, data = cebq_LOC_ndup)
cebqFR_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqFR, cebq_LOC_ndup$loc1)
cebqFR_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqFR, cebq_LOC_ndup$loc1)
cebqFR_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqFR, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqFR)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqFR) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqFR) && cebq_LOC_ndup$loc1 == 'No', ])

## EOE

cebq_LOC_ndup$cebqEOE = rowMeans(data.frame(cebq_LOC_ndup$cebq2_raw, cebq_LOC_ndup$cebq13_raw, cebq_LOC_ndup$cebq15_raw, 
  cebq_LOC_ndup$cebq23_raw, cebq_LOC_ndup$cebq27_raw), na.rm = T)
alpha_cebqEOE = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq2_raw, cebq_LOC_ndup$cebq13_raw, cebq_LOC_ndup$cebq15_raw, 
  cebq_LOC_ndup$cebq23_raw, cebq_LOC_ndup$cebq27_raw)))

cebqEOE_mean = mean(cebq_LOC_ndup$cebqEOE, na.rm = TRUE)
cebqEOE_sd = sd(cebq_LOC_ndup$cebqEOE, na.rm = TRUE)
cebqEOE_range = range(cebq_LOC_ndup$cebqEOE, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEOE), ])

cebqEOE_loc_t = t.test(cebqEOE~loc1, data = cebq_LOC_ndup)
cebqEOE_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqEOE, cebq_LOC_ndup$loc1)
cebqEOE_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqEOE, cebq_LOC_ndup$loc1)
cebqEOE_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqEOE, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqEOE)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEOE) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEOE) && cebq_LOC_ndup$loc1 == 'No', ])

## EF
cebq_LOC_ndup$cebqEF = rowMeans(data.frame(cebq_LOC_ndup$cebq1_raw, cebq_LOC_ndup$cebq5_raw, cebq_LOC_ndup$cebq20_raw, 
  cebq_LOC_ndup$cebq22_raw), na.rm = T)
alpha_cebqEF = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq1_raw, cebq_LOC_ndup$cebq5_raw, cebq_LOC_ndup$cebq20_raw, 
  cebq_LOC_ndup$cebq22_raw)))

cebqEF_mean = mean(cebq_LOC_ndup$cebqEF, na.rm = TRUE)
cebqEF_sd = sd(cebq_LOC_ndup$cebqEF, na.rm = TRUE)
cebqEF_range = range(cebq_LOC_ndup$cebqEF, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEF), ])

cebqEF_loc_t = t.test(cebqEF~loc1, data = cebq_LOC_ndup)
cebqEF_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqEF, cebq_LOC_ndup$loc1)
cebqEF_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqEF, cebq_LOC_ndup$loc1)
cebqEF_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqEF, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqEF)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEF) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEF) && cebq_LOC_ndup$loc1 == 'No', ])

## DD
cebq_LOC_ndup$cebqDD = rowMeans(data.frame(cebq_LOC_ndup$cebq6_raw, cebq_LOC_ndup$cebq29_raw, cebq_LOC_ndup$cebq31_raw), na.rm = T)
alpha_cebqDD = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq6_raw, cebq_LOC_ndup$cebq29_raw, cebq_LOC_ndup$cebq31_raw)))

cebqDD_mean = mean(cebq_LOC_ndup$cebqDD, na.rm = TRUE)
cebqDD_sd = sd(cebq_LOC_ndup$cebqDD, na.rm = TRUE)
cebqDD_range = range(cebq_LOC_ndup$cebqDD, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqDD), ])

cebqDD_loc_t = t.test(cebqDD~loc1, data = cebq_LOC_ndup)
cebqDD_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqDD, cebq_LOC_ndup$loc1)
cebqDD_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqDD, cebq_LOC_ndup$loc1)
cebqDD_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqDD, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqDD)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqDD) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqDD) && cebq_LOC_ndup$loc1 == 'No', ])

## SR
cebq_LOC_ndup$cebqSR = rowMeans(data.frame(cebq_LOC_ndup$cebq3_raw, cebq_LOC_ndup$cebq17_raw, cebq_LOC_ndup$cebq21_raw, 
  cebq_LOC_ndup$cebq26_raw, cebq_LOC_ndup$cebq30_raw), na.rm = T)
alpha_cebqSR = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq3_raw, cebq_LOC_ndup$cebq17_raw, cebq_LOC_ndup$cebq21_raw, 
  cebq_LOC_ndup$cebq26_raw, cebq_LOC_ndup$cebq30_raw)))

cebqSR_mean = mean(cebq_LOC_ndup$cebqSR, na.rm = TRUE)
cebqSR_sd = sd(cebq_LOC_ndup$cebqSR, na.rm = TRUE)
cebqSR_range = range(cebq_LOC_ndup$cebqSR, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqSR), ])

cebqSR_loc_t = t.test(cebqSR~loc1, data = cebq_LOC_ndup)
cebqSR_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqSR, cebq_LOC_ndup$loc1)
cebqSR_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqSR, cebq_LOC_ndup$loc1)
cebqSR_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqSR, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqSR)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqSR) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqSR) && cebq_LOC_ndup$loc1 == 'No', ])

## SE
cebq_LOC_ndup$cebqSE = rowMeans(data.frame(cebq_LOC_ndup$cebq4_raw, cebq_LOC_ndup$cebq8_raw, cebq_LOC_ndup$cebq18_raw, 
  cebq_LOC_ndup$cebq35_raw), na.rm = T)
alpha_cebqSE = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq4_raw, cebq_LOC_ndup$cebq8_raw, cebq_LOC_ndup$cebq18_raw, 
  cebq_LOC_ndup$cebq35_raw)))

cebqSE_mean = mean(cebq_LOC_ndup$cebqSE, na.rm = TRUE)
cebqSE_sd = sd(cebq_LOC_ndup$cebqSE, na.rm = TRUE)
cebqSE_range = range(cebq_LOC_ndup$cebqSE, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqSE), ])

cebqSE_loc_t = t.test(cebqSE~loc1, data = cebq_LOC_ndup)
cebqSE_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqSE, cebq_LOC_ndup$loc1)
cebqSE_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqSE, cebq_LOC_ndup$loc1)
cebqSE_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqSE, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqSE)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqSE) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqSE) && cebq_LOC_ndup$loc1 == 'No', ])

## EUE
cebq_LOC_ndup$cebqEUE = rowMeans(data.frame(cebq_LOC_ndup$cebq9_raw, cebq_LOC_ndup$cebq11_raw, cebq_LOC_ndup$cebq25_raw), na.rm = T)
alpha_cebqEUE = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq9_raw, cebq_LOC_ndup$cebq11_raw, cebq_LOC_ndup$cebq25_raw)))

cebqEUE_mean = mean(cebq_LOC_ndup$cebqEUE, na.rm = TRUE)
cebqEUE_sd = sd(cebq_LOC_ndup$cebqEUE, na.rm = TRUE)
cebqEUE_range = range(cebq_LOC_ndup$cebqEUE, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEUE), ])

cebqEUE_loc_t = t.test(cebqEUE~loc1, data = cebq_LOC_ndup)
cebqEUE_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqEUE, cebq_LOC_ndup$loc1)
cebqEUE_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqEUE, cebq_LOC_ndup$loc1)
cebqEUE_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqEUE, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqEUE)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEUE) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqEUE) && cebq_LOC_ndup$loc1 == 'No', ])

## FF
cebq_LOC_ndup$cebqFF = rowMeans(data.frame(cebq_LOC_ndup$cebq7_raw, cebq_LOC_ndup$cebq10_raw, cebq_LOC_ndup$cebq16_raw, 
  cebq_LOC_ndup$cebq24_raw, cebq_LOC_ndup$cebq32_raw, cebq_LOC_ndup$cebq33_raw), na.rm = T)
alpha_cebqFF = summary(psych::alpha(data.frame(cebq_LOC_ndup$cebq7_raw, cebq_LOC_ndup$cebq10_raw, cebq_LOC_ndup$cebq16_raw, 
  cebq_LOC_ndup$cebq24_raw, cebq_LOC_ndup$cebq32_raw, cebq_LOC_ndup$cebq33_raw)))

cebqFF_mean = mean(cebq_LOC_ndup$cebqFF, na.rm = TRUE)
cebqFF_sd = sd(cebq_LOC_ndup$cebqFF, na.rm = TRUE)
cebqFF_range = range(cebq_LOC_ndup$cebqFF, na.rm = TRUE)
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqFF), ])

cebqFF_loc_t = t.test(cebqFF~loc1, data = cebq_LOC_ndup)
cebqFF_loc_sd = sd.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqFF, cebq_LOC_ndup$loc1)
cebqFF_loc_range = range.function.na(cebq_LOC_ndup, cebq_LOC_ndup$cebqFF, cebq_LOC_ndup$loc1)
cebqFF_loc_d = cohensD(cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'Yes', ]$cebqFF, cebq_LOC_ndup[cebq_LOC_ndup$loc1 == 'No', ]$cebqFF)

nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqFF) && cebq_LOC_ndup$loc1 == 'Yes', ])
nrow(cebq_LOC_ndup[is.na(cebq_LOC_ndup$cebqFF) && cebq_LOC_ndup$loc1 == 'No', ])

##### CFQ ####
#as as for cebq above
cfq_LOC_ndup = LOC_ndup[c(2, 87:117)]


for(c in 2:32){
  cfq_LOC_ndup[[c+31]] = NA
  varname = names(cfq_LOC_ndup)[[c]]
  names(cfq_LOC_ndup)[[c+31]] = paste(varname, "raw", sep = "_")
  
  if(c < 5){
    cfq_LOC_ndup[c+31] = ifelse(is.na(cfq_LOC_ndup[[c]]), NA, ifelse(
      cfq_LOC_ndup[[c]] == "Never", 1, ifelse(
        cfq_LOC_ndup[[c]] == "Seldom", 2, ifelse(
          cfq_LOC_ndup[[c]] == "HalfTimes", 3, ifelse(
            cfq_LOC_ndup[[c]] == "MostTimes", 4, 5
          )
        )
      )
    )) 
  }
  if(c > 4 & c < 15){
    cfq_LOC_ndup[c+31] = ifelse(is.na(cfq_LOC_ndup[[c]]), NA, ifelse(
      cfq_LOC_ndup[[c]] == "MarkedlyUnderwt", 1, ifelse(
        cfq_LOC_ndup[[c]] == "Underwt", 2, ifelse(
          cfq_LOC_ndup[[c]] == "Avg", 3, ifelse(
            cfq_LOC_ndup[[c]] == "Overwt", 4, 5
          )
        )
      )
    )) 
  }
  if( c > 14 & c < 18){
    cfq_LOC_ndup[c+31] = ifelse(is.na(cfq_LOC_ndup[[c]]), NA, ifelse(
      cfq_LOC_ndup[[c]] == "Unconcerned", 1, ifelse(
        cfq_LOC_ndup[[c]] == "SlightlyUncon", 2, ifelse(
          cfq_LOC_ndup[[c]] == "Concerned", 3, ifelse(
            cfq_LOC_ndup[[c]] == "Slightlycon", 4, 5
          )
        )
      )
    )) 
  }
  if(c > 17 & c < 30){
    cfq_LOC_ndup[c+31] = ifelse(is.na(cfq_LOC_ndup[[c]]), NA, ifelse(
      cfq_LOC_ndup[[c]] == "Disagree", 1, ifelse(
        cfq_LOC_ndup[[c]] == "SlightlyDis", 2, ifelse(
          cfq_LOC_ndup[[c]] == "Neutral", 3, ifelse(
            cfq_LOC_ndup[[c]] == "SlightlyAgree", 4, 5
          )
        )
      )
    )) 
  }
  if(c > 29){
    cfq_LOC_ndup[c+31] = ifelse(is.na(cfq_LOC_ndup[[c]]), NA, ifelse(
      cfq_LOC_ndup[[c]] == "Never", 1, ifelse(
        cfq_LOC_ndup[[c]] == "Rarely", 2, ifelse(
          cfq_LOC_ndup[[c]] == "Sometimes", 3, ifelse(
            cfq_LOC_ndup[[c]] == "Mostly", 4, 5
          )
        )
      )
    )) 
  }
}

#add loc
cfq_LOC_ndup = merge(cfq_LOC_ndup, LOC_ndup[c(2, 33)], by = 'StudyID')

#get means and descriptives for CFQ questionnaire
## PR
cfq_LOC_ndup$cfqPR = rowMeans(cfq_LOC_ndup[33:35], na.rm = T)
alpha_cfqPR = summary(psych::alpha(cfq_LOC_ndup[33:35]))

cfqPR_mean = mean(cfq_LOC_ndup$cfqPR, na.rm = TRUE)
cfqPR_sd = sd(cfq_LOC_ndup$cfqPR, na.rm = TRUE)
cfqPR_range = range(cfq_LOC_ndup$cfqPR, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPR), ])

cebqPR_loc_t = t.test(cfqPR~loc1, data = cfq_LOC_ndup)
cebqPR_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPR, cfq_LOC_ndup$loc1)
cebqPR_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPR, cfq_LOC_ndup$loc1)
cebqPR_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqPR, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqPR)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPR) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPR) && cfq_LOC_ndup$loc1 == 'No', ])

## PPW
cfq_LOC_ndup$cfqPPW = rowMeans(cfq_LOC_ndup[36:39], na.rm = T)
alpha_cfqPPW = summary(psych::alpha(cfq_LOC_ndup[36:39]))

cfqPPW_mean = mean(cfq_LOC_ndup$cfqPPW, na.rm = TRUE)
cfqPPW_sd = sd(cfq_LOC_ndup$cfqPPW, na.rm = TRUE)
cfqPPW_range = range(cfq_LOC_ndup$cfqPPW, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPPW), ])

cebqPPW_loc_t = t.test(cfqPPW~loc1, data = cfq_LOC_ndup)
cebqPPW_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPPW, cfq_LOC_ndup$loc1)
cebqPPW_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPPW, cfq_LOC_ndup$loc1)
cebqPPW_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqPPW, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqPPW)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPPW) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPPW) && cfq_LOC_ndup$loc1 == 'No', ])

## PCW
cfq_LOC_ndup$cfqPCW = rowMeans(cfq_LOC_ndup[40:45], na.rm = T)
alpha_cfqPCW = summary(psych::alpha(cfq_LOC_ndup[40:45]))

cfqPCW_mean = mean(cfq_LOC_ndup$cfqPCW, na.rm = TRUE)
cfqPCW_sd = sd(cfq_LOC_ndup$cfqPCW, na.rm = TRUE)
cfqPCW_range = range(cfq_LOC_ndup$cfqPCW, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPCW), ])

cebqPCW_loc_t = t.test(cfqPCW~loc1, data = cfq_LOC_ndup)
cebqPCW_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPCW, cfq_LOC_ndup$loc1)
cebqPCW_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPCW, cfq_LOC_ndup$loc1)
cebqPCW_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqPCW, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqPCW)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPCW) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPCW) && cfq_LOC_ndup$loc1 == 'No', ])

## CONC
cfq_LOC_ndup$cfqCONC = rowMeans(cfq_LOC_ndup[46:48], na.rm = T)
alpha_cfqCONC = summary(psych::alpha(cfq_LOC_ndup[46:48]))

cfqCONC_mean = mean(cfq_LOC_ndup$cfqCONC, na.rm = TRUE)
cfqCONC_sd = sd(cfq_LOC_ndup$cfqCONC, na.rm = TRUE)
cfqCONC_range = range(cfq_LOC_ndup$cfqCONC, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqCONC), ])

cebqCONC_loc_t = t.test(cfqCONC~loc1, data = cfq_LOC_ndup)
cebqCONC_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqCONC, cfq_LOC_ndup$loc1)
cebqCONC_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqCONC, cfq_LOC_ndup$loc1)
cebqCONC_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqCONC, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqCONC)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqCONC) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqCONC) && cfq_LOC_ndup$loc1 == 'No', ])

## REST
cfq_LOC_ndup$cfqREST = rowMeans(cfq_LOC_ndup[49:56], na.rm = T)
alpha_cfqREST = summary(psych::alpha(cfq_LOC_ndup[49:56]))

cfqREST_mean = mean(cfq_LOC_ndup$cfqREST, na.rm = TRUE)
cfqREST_sd = sd(cfq_LOC_ndup$cfqREST, na.rm = TRUE)
cfqREST_range = range(cfq_LOC_ndup$cfqREST, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqREST), ])

cebqREST_loc_t = t.test(cfqREST~loc1, data = cfq_LOC_ndup)
cebqREST_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqREST, cfq_LOC_ndup$loc1)
cebqREST_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqREST, cfq_LOC_ndup$loc1)
cebqREST_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqREST, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqREST)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqREST) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqREST) && cfq_LOC_ndup$loc1 == 'No', ])

## PE
cfq_LOC_ndup$cfqPE = rowMeans(cfq_LOC_ndup[57:60], na.rm = T)
alpha_cfqPE = summary(psych::alpha(cfq_LOC_ndup[57:60]))

cfqPE_mean = mean(cfq_LOC_ndup$cfqPE, na.rm = TRUE)
cfqPE_sd = sd(cfq_LOC_ndup$cfqPE, na.rm = TRUE)
cfqPE_range = range(cfq_LOC_ndup$cfqPE, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPE), ])

cebqPE_loc_t = t.test(cfqPE~loc1, data = cfq_LOC_ndup)
cebqPE_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPE, cfq_LOC_ndup$loc1)
cebqPE_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqPE, cfq_LOC_ndup$loc1)
cebqPE_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqPE, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqPE)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPE) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqPE) && cfq_LOC_ndup$loc1 == 'No', ])

## MON
cfq_LOC_ndup$cfqMON = rowMeans(cfq_LOC_ndup[61:63], na.rm = T)
alpha_cfqMON = summary(psych::alpha(cfq_LOC_ndup[61:63]))

cfqMON_mean = mean(cfq_LOC_ndup$cfqMON, na.rm = TRUE)
cfqMON_sd = sd(cfq_LOC_ndup$cfqMON, na.rm = TRUE)
cfqMON_range = range(cfq_LOC_ndup$cfqMON, na.rm = TRUE)
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqMON), ])

cebqMON_loc_t = t.test(cfqMON~loc1, data = cfq_LOC_ndup)
cebqMON_loc_sd = sd.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqMON, cfq_LOC_ndup$loc1)
cebqMON_loc_range = range.function.na(cfq_LOC_ndup, cfq_LOC_ndup$cfqMON, cfq_LOC_ndup$loc1)
cebqMON_loc_d = cohensD(cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'Yes', ]$cfqMON, cfq_LOC_ndup[cfq_LOC_ndup$loc1 == 'No', ]$cfqMON)

nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqMON) && cfq_LOC_ndup$loc1 == 'Yes', ])
nrow(cfq_LOC_ndup[is.na(cfq_LOC_ndup$cfqMON) && cfq_LOC_ndup$loc1 == 'No', ])


###############################################
####            ASSOCIATION RULES          ####
###############################################
#mother completed
ParentComplete_noNA = xtabs(~relationship, data = LOC_ndup[!is.na(LOC_ndup$loc1), ])
ParentComplete = xtabs(~relationship, data = LOC_ndup)

#### convert to transactions ####
#all questions expanded below
LOC_arules = LOC_ndup[!is.na(LOC_ndup$loc1), c(4, 119:127, 129:134, 135:168, 33, 52:117, 169:335)]

LOC_arules_trans = as(LOC_arules, "transactions")

#get all transactions (going to be the LHS items)--only have to do this once then comment out
#transDataFrame = as.matrix(LOC_arules_trans@itemInfo$labels)
#export to manually add the category labes you want
#write.csv(transDataFrame, file = "ResultsOutput/indQ_TransactionsLabels.csv", row.names = FALSE)

#read in transaction labels
transDataLabels = read.csv("ResultsOutput/indQ_TransactionsLabels.csv", header = TRUE)

#### item frequency ####
LOC_arules_LOC = LOC_arules[LOC_arules$loc1 == 'Yes', ]
LOC_arules_noLOC = LOC_arules[LOC_arules$loc1 == 'No', ]

LOC_arules_trans_LOConly = as(LOC_arules_LOC, "transactions")
LOC_arules_trans_noLOConly = as(LOC_arules_noLOC, "transactions")

LOC_freq = data.frame(itemFrequency(LOC_arules_trans_LOConly,type="absolute"))
names(LOC_freq) = 'Freq'
LOC_freq$item = row.names(LOC_freq)
LOC_freq = LOC_freq[order(LOC_freq$Freq, decreasing = TRUE), ]

LOC_freq_lolipop = ggplot(head(LOC_freq, 20), aes(x=fct_reorder(item, Freq, .desc =TRUE),y=Freq)) +
  geom_point(size = 3, colour = "cornflowerblue") + 
  geom_segment(aes(xend = item, yend = 0), size = 1.2, color = 'cornflowerblue')+
  geom_text(aes(label=round(Freq/121, digits = 2)),hjust=-0.25, vjust = 1.5)+
  #geom_label(aes(item, Freq+1.5, label = Freq), colour = "darkred", nudge_x = 0.35, size = 4)+
  labs(y= "Frequency", x="Item")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90), axis.ticks = element_blank())

noLOC_freq = data.frame(itemFrequency(LOC_arules_trans_noLOConly,type="absolute"))
names(noLOC_freq) = 'Freq'
noLOC_freq$item = row.names(noLOC_freq)
noLOC_freq = noLOC_freq[order(noLOC_freq$Freq, decreasing = TRUE), ]

noLOC_freq_lolipop = ggplot(head(noLOC_freq, 20), aes(x=fct_reorder(item, Freq, .desc =TRUE),y=Freq)) +
  geom_point(size = 3, colour = "cornflowerblue") + 
  geom_segment(aes(xend = item, yend = 0), size = 1.2, color = 'cornflowerblue')+
  geom_text(aes(label=round(Freq/121, digits = 2)),hjust=-0.25, vjust = 1.5)+
  #geom_label(aes(item, Freq+1.5, label = Freq), colour = "darkred", nudge_x = 0.35, size = 4)+
  labs(y= "Frequency", x="Item")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90), axis.ticks = element_blank())

#### find Association rules####
#LOC
#25% of 37 (n participants indicating LOC) = 9.25; 9/177 = 0.051; 
#remove NAs and there are 158 total so 9/158 = 0.06 
#25% -- support will be set between .05 (n=9) and .056 (n=10)
#33% -- support of 6.7%, n=12, almost 1/3 of LOC
#no NAs 33% -- support of 0.076

#No LOC
#25% of 121 (n participants indicating no LOC) = 30; 30/177 = 0.17; 
#remove NAs and there are 158 total so 30/158 = 0.19 
#33% 0f 121 (n participants indicating no LOC) = 39; 40/177 = 0.23; 
#remove NAs and there are 158 total so 40/158 = 0.255 

###############################################
####      LOC-yes; Single category Questionnaire Variables Included      ####
###############################################

### loc-yes 1, conf=.33, sup = 0.055 ####
#get rules
rules_lhs1_conf33 = apriori(LOC_arules_trans, parameter = list(maxlen = 2, supp = 0.06, conf = 0.33), appearance = list(default="lhs", rhs="loc1=Yes"))

#prune rules that are redundant-more important for multiple predictors as may get two rules with 
#different predictor order, thus are the same
pruned.rules_lhs1_conf33 = rules_lhs1_conf33[!is.redundant(rules_lhs1_conf33)]

#add odds ratios, chi-squared and its p, and fisher p to rules information
quality(pruned.rules_lhs1_conf33) = cbind(quality(pruned.rules_lhs1_conf33), 
  oddsRatio = interestMeasure(pruned.rules_lhs1_conf33, 
    measure = "oddsRatio", 
    transactions = LOC_arules_trans),
  addVal = interestMeasure(pruned.rules_lhs1_conf33, 
    measure = "addedValue", 
    transactions = LOC_arules_trans),
  kappa = interestMeasure(pruned.rules_lhs1_conf33, 
    measure = "kappa", 
    transactions = LOC_arules_trans),
  fisher.p = interestMeasure(pruned.rules_lhs1_conf33, 
    measure = "fishersExactTest",
    transactions = LOC_arules_trans))

#remove empty rules
pruned.rules_lhs1_conf33 = subset(pruned.rules_lhs1_conf33, subset = size(lhs(pruned.rules_lhs1_conf33))!=0)

#check added value and correlation
pruned.rules_lhs1_conf33_qdat = data.frame(quality(pruned.rules_lhs1_conf33))
nrow(pruned.rules_lhs1_conf33_qdat[pruned.rules_lhs1_conf33_qdat$kappa>=0.21, ])
nrow(pruned.rules_lhs1_conf33_qdat[pruned.rules_lhs1_conf33_qdat$addVal>=0.05, ])

pruned.rules_lhs1_conf33 = subset(pruned.rules_lhs1_conf33, subset = quality(pruned.rules_lhs1_conf33)[8]>=0.21)
pruned.rules_lhs1_conf33 = subset(pruned.rules_lhs1_conf33, subset = quality(pruned.rules_lhs1_conf33)[7]>=0.05)

#add addjustment for multiple comparisons to data
quality(pruned.rules_lhs1_conf33) = cbind(quality(pruned.rules_lhs1_conf33), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs1_conf33)$fisher.p, method = 'holm'), 4))

#sort rules
pruned.rules_lhs1_conf33 = sort(pruned.rules_lhs1_conf33, by=c("fisher.p","fisher.padj_holm"), decreasing = FALSE)

#HERE
#make list of significant predictors
lhs1_conf33_vars = c("cfq26=Disagree", "cfq28=Disagree", "cebq34_gRarely", "cfq23_lsNeutralAD", "cebq12_gRarely",
  "cebq7_gSometimes", "cfq20=Agree", "cebq33_gSometimes", "mEducation_HS", "cfq29=Always")

#get odds ratio confidence intervals for significant uncorrected fisher pvalue rules
#(all rules in this case)
OR_CI_pruned.rules_lhs1_conf33 = ARM_ORconf(pruned.rules_lhs1_conf33, LOC_arules, 1, "Yes")

#add odds ratio CI to rules information
quality(pruned.rules_lhs1_conf33) = cbind(quality(pruned.rules_lhs1_conf33),
  OR_CI_pruned.rules_lhs1_conf33[3:4])

#organize columns 
quality(pruned.rules_lhs1_conf33) = quality(pruned.rules_lhs1_conf33)[c(5, 1:4, 7:8, 6, 11:12, 9:10)]

#convert rules to a dataframe
DataFrame_pruned.rules_lhs1_conf33 = DATAFRAME(pruned.rules_lhs1_conf33, separate = TRUE, setStart = '', itemSep = '+', setEnd = '')

#merge with transaction labels dataframe to get category labels for each rule
DataFrame_pruned.rules_lhs1_conf33 = merge(DataFrame_pruned.rules_lhs1_conf33, transDataLabels, by.x = "LHS", by.y = "LHS_label")

#drop extra/unused levels of category variable
DataFrame_pruned.rules_lhs1_conf33$Cat = factor(DataFrame_pruned.rules_lhs1_conf33$Cat)

#make crosstab of rule categories
xtabs_pruned.rules_lhs1_conf33 = xtabs(~Cat, data = DataFrame_pruned.rules_lhs1_conf33)

#write out if want
write.csv(DataFrame_pruned.rules_lhs1_conf33, file = "ResultsOutput/indQ_LOC-YES_lhs1_conf33.csv", row.names = FALSE)

### loc-yes 2, conf=.50, sup = 0.08 (1/3) ####
#get rules
rules_lhs2_conf50 = apriori(LOC_arules_trans, parameter = list(maxlen = 3, supp = 0.08, conf = 0.50), appearance = list(default="lhs", rhs="loc1=Yes"))

#prune rules that are redundant-more important for multiple predictors as may get two rules with 
#different predictor order, thus are the same
pruned.rules_lhs2_conf50 = rules_lhs2_conf50[!is.redundant(rules_lhs2_conf50)]

#add odds ratio, chi-square and its p, and fisher p to rules information
quality(pruned.rules_lhs2_conf50) = cbind(quality(pruned.rules_lhs2_conf50), 
  oddsRatio = interestMeasure(pruned.rules_lhs2_conf50, 
    measure = "oddsRatio", 
    transactions = LOC_arules_trans),
  addVal = interestMeasure(pruned.rules_lhs2_conf50, 
    measure = "addedValue", 
    transactions = LOC_arules_trans),
  kappa = interestMeasure(pruned.rules_lhs2_conf50, 
    measure = "kappa", 
    transactions = LOC_arules_trans),
  fisher.p = interestMeasure(pruned.rules_lhs2_conf50, 
    measure = "fishersExactTest", 
    transactions = LOC_arules_trans))

#remove empty rules
pruned.rules_lhs2_conf50 = subset(pruned.rules_lhs2_conf50, subset = size(lhs(pruned.rules_lhs2_conf50))!=0)

#prune rules for correlation and added value
pruned.rules_lhs2_conf50_qdat = data.frame(quality(pruned.rules_lhs2_conf50))
nrow(pruned.rules_lhs2_conf50_qdat[pruned.rules_lhs2_conf50_qdat$kappa >=0.21, ])
nrow(pruned.rules_lhs2_conf50_qdat[pruned.rules_lhs2_conf50_qdat$addVal >=0.05, ])

pruned.rules_lhs2_conf50 = subset(pruned.rules_lhs2_conf50, subset = quality(pruned.rules_lhs2_conf50)[8]>=0.20)
pruned.rules_lhs2_conf50 = subset(pruned.rules_lhs2_conf50, subset = quality(pruned.rules_lhs2_conf50)[7]>=0.05)

### no subset--all rules  ####
#rename so don't overwrite original
pruned.rules_lhs2_conf50_all = pruned.rules_lhs2_conf50

#add correction for multiple comparisons
quality(pruned.rules_lhs2_conf50_all) = cbind(quality(pruned.rules_lhs2_conf50_all), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs2_conf50_all)$fisher.p, method = 'holm'), 4))

#add correction for multiple comparisons
pruned.rules_lhs2_conf50_all = sort(pruned.rules_lhs2_conf50_all, by=c("fisher.p","fisher.padj_holm"), decreasing = FALSE)

#get odds ratio confidence intervals for significant uncorrected fisher pvalue rules
#(all rules in this case)
OR_CI_pruned.rules_lhs2_conf50_all = ARM_ORconf(pruned.rules_lhs2_conf50_all, LOC_arules, 2, "Yes")

#add odds ration CI to rules information
quality(pruned.rules_lhs2_conf50_all) = cbind(quality(pruned.rules_lhs2_conf50_all),
  OR_CI_pruned.rules_lhs2_conf50_all[4:5])

#organize columns
quality(pruned.rules_lhs2_conf50_all) = quality(pruned.rules_lhs2_conf50_all)[c(5, 1:4, 7:8, 6, 11:12, 9:10)]

### clusters loc-yes 2, conf=.50, sup = 0.23 (1/3) to ####
#get matrix information for final set (holm sig/trend)
im_lhs2_conf50_nosubset = quality(pruned.rules_lhs2_conf50_all)

#clustering Rules lhs2 conf50 lhs1 sig chi
dist_pruned.rules_lhs2_conf50 = dissimilarity(pruned.rules_lhs2_conf50,args = list(transactions = LOC_arules_trans), method = "gupta")

## mediods partitioning--determine number of clusters ####
#within cluster sums of squares--look for elbow in graph
pruned.rules_lhs2_conf50_nclustSSW = fviz_nbclust(as.matrix(dist_pruned.rules_lhs2_conf50), cluster::pam, method = "wss", k.max=15)+ geom_vline(xintercept = 4, linetype = 2)+labs(subtitle = "Elbow method")

#silhouette width
pruned.rules_lhs2_conf50_nclustSil = fviz_nbclust(as.matrix(dist_pruned.rules_lhs2_conf50), cluster::pam, method = "silhouette", k.max=15)+ labs(subtitle = "Silhouette method")

#get combind plot
pruned.rules_lhs2_conf50_nclustPlotGrid = plot_grid(pruned.rules_lhs2_conf50_nclustSSW, pruned.rules_lhs2_conf50_nclustSil, labels = c("", ""))

#test internal validity metrics
#pruned.rules_lhs2_conf50_nclust_internV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50), 2:10, clMethods=c("pam"), validation="internal")

#test stability metrics (remove one collumn at a time sequentially and look at stability)
#pruned.rules_lhs2_conf50_nclust_stabV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50), 2:10, clMethods=c("pam"), validation="stability")

#cluster validity table
pruned.rules_lhs2_conf50_nclust_Vsum = data.frame(matrix(c("Connectivity", 0.000,    "pam", 2,
  "Dunn",         0.593,    "pam", 10,
  "Silhouette",   0.498,    "pam", 4,
  "APN",          0.0000,    "pam", 4,
  "AD",           0.5054,    "pam", 10,
  "ADM",          0.0000,    "pam", 4,
  "FOM",          0.0753,    "pam", 10), 
  byrow = TRUE, nrow = 7, ncol = 4))
names(pruned.rules_lhs2_conf50_nclust_Vsum) = c("measure", "score_at5", "method", "clusters")     

## run pam to get cluster membership ####
pruned.rules_lhs2_conf50_4clust = pam(dist_pruned.rules_lhs2_conf50, diss = TRUE, k=4)

#silhouette widths
pruned.rules_lhs2_conf50_4clust_asw = data.frame(c(1, 2, 3, 4), pruned.rules_lhs2_conf50_4clust$silinfo$clus.avg.widths)
names(pruned.rules_lhs2_conf50_4clust_asw) = c("cluster", "asw")

#add cluster belonging to  quality matrix
im_lhs2_conf50_nosubset$clust_gupta4 = pruned.rules_lhs2_conf50_4clust$clustering

#add to rules quality
quality(pruned.rules_lhs2_conf50_all) = cbind(quality(pruned.rules_lhs2_conf50_all), Cluster_gupta4 = im_lhs2_conf50_nosubset$clust_gupta4)

#convert rules to dataframe
DataFrame_pruned.rules_lhs2_conf50_nosubset = DATAFRAME(pruned.rules_lhs2_conf50_all, separate = TRUE, setStart = '', itemSep = '+', setEnd = '')

#need the following two split the LHS into 2 columns
DataFrame_pruned.rules_lhs2_conf50_nosubset = with(DataFrame_pruned.rules_lhs2_conf50_nosubset,
  cbind(colsplit(DataFrame_pruned.rules_lhs2_conf50_nosubset$LHS,"\\+", c('LHS1', 'LHS2')),
    DataFrame_pruned.rules_lhs2_conf50_nosubset[2:15]))


#merge with transaction labels dataframe to get category labels for each rule in LHS1
DataFrame_pruned.rules_lhs2_conf50_nosubset = merge(DataFrame_pruned.rules_lhs2_conf50_nosubset, transDataLabels, by.x = "LHS1", by.y = "LHS_label")

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_nosubset)[17] = "Cat1"

#remove/drop unused levels of Cat1
DataFrame_pruned.rules_lhs2_conf50_nosubset$Cat1 = factor(DataFrame_pruned.rules_lhs2_conf50_nosubset$Cat1)

#merge with transaction labels dataframe to get category labels for each rule in LHS2
DataFrame_pruned.rules_lhs2_conf50_nosubset = merge(DataFrame_pruned.rules_lhs2_conf50_nosubset, transDataLabels, by.x = "LHS2", by.y = "LHS_label")

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_nosubset)[18] = "Cat2"

#remove/drop unused levels of Cat2
DataFrame_pruned.rules_lhs2_conf50_nosubset$Cat2 = factor(DataFrame_pruned.rules_lhs2_conf50_nosubset$Cat2)

#make crosstab of rule categories
xtabs_pruned.rules_lhs2_conf50_nosubset = xtabs(~Cat1 + Cat2, data = DataFrame_pruned.rules_lhs2_conf50_nosubset)

#write out if want
write.csv(DataFrame_pruned.rules_lhs2_conf50_nosubset, file = "ResultsOutput/indQ_LOC-YES_lhs2_conf50_nosubset.csv", row.names = TRUE)

## get crosstables for each cluster--5 cluster solution ####
#--need to first subset and drop unsued levels then make table for each cluster
#cluster 4-1
#subset to just cluster 1
pruned.rules_lhs2_conf50_nosubset_clust4.1 = DataFrame_pruned.rules_lhs2_conf50_nosubset[DataFrame_pruned.rules_lhs2_conf50_nosubset$Cluster_gupta4 == 1, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_nosubset_clust4.1$LHS1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.1$LHS1)
pruned.rules_lhs2_conf50_nosubset_clust4.1$LHS2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.1$LHS2)
pruned.rules_lhs2_conf50_nosubset_clust4.1$Cat1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.1$Cat1)
pruned.rules_lhs2_conf50_nosubset_clust4.1$Cat2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.1$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1 = as.matrix(xtabs(~Cat1 + Cat2, data = pruned.rules_lhs2_conf50_nosubset_clust4.1))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1 = rbind( xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1 = cbind(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.1))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.1_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_nosubset_clust4.1))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.1_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_nosubset_clust4.1))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.1 = merge(qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.1_LHS1, qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.1_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 4-2
#subset to just cluster 2
pruned.rules_lhs2_conf50_nosubset_clust4.2 = DataFrame_pruned.rules_lhs2_conf50_nosubset[DataFrame_pruned.rules_lhs2_conf50_nosubset$Cluster_gupta4 == 2, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_nosubset_clust4.2$LHS1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.2$LHS1)
pruned.rules_lhs2_conf50_nosubset_clust4.2$LHS2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.2$LHS2)
pruned.rules_lhs2_conf50_nosubset_clust4.2$Cat1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.2$Cat1)
pruned.rules_lhs2_conf50_nosubset_clust4.2$Cat2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.2$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2 = as.matrix(xtabs(~Cat1 + Cat2, data = pruned.rules_lhs2_conf50_nosubset_clust4.2))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2 = rbind( xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2 = cbind(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.2))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.2_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_nosubset_clust4.2))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.2_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_nosubset_clust4.2))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.2 = merge(qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.2_LHS1, qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.2_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 4-3
#subset to just cluster 3
pruned.rules_lhs2_conf50_nosubset_clust4.3 = DataFrame_pruned.rules_lhs2_conf50_nosubset[DataFrame_pruned.rules_lhs2_conf50_nosubset$Cluster_gupta4 == 3, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_nosubset_clust4.3$LHS1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.3$LHS1)
pruned.rules_lhs2_conf50_nosubset_clust4.3$LHS2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.3$LHS2)
pruned.rules_lhs2_conf50_nosubset_clust4.3$Cat1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.3$Cat1)
pruned.rules_lhs2_conf50_nosubset_clust4.3$Cat2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.3$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3 = as.matrix(xtabs(~Cat1 + Cat2, data = pruned.rules_lhs2_conf50_nosubset_clust4.3))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3 = rbind( xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3 = cbind(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.3))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.3_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_nosubset_clust4.3))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.3_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_nosubset_clust4.3))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.3 = merge(qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.3_LHS1, qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.3_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 4-4
#subset to just cluster 2
pruned.rules_lhs2_conf50_nosubset_clust4.4 = DataFrame_pruned.rules_lhs2_conf50_nosubset[DataFrame_pruned.rules_lhs2_conf50_nosubset$Cluster_gupta4 == 4, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_nosubset_clust4.4$LHS1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.4$LHS1)
pruned.rules_lhs2_conf50_nosubset_clust4.4$LHS2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.4$LHS2)
pruned.rules_lhs2_conf50_nosubset_clust4.4$Cat1 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.4$Cat1)
pruned.rules_lhs2_conf50_nosubset_clust4.4$Cat2 = factor(pruned.rules_lhs2_conf50_nosubset_clust4.4$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4 = as.matrix(xtabs(~Cat1 + Cat2, data = pruned.rules_lhs2_conf50_nosubset_clust4.4))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4 = rbind( xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4))
xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4 = cbind(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_nosubset_clust4.4))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.4_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_nosubset_clust4.4))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.4_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_nosubset_clust4.4))
qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.4 = merge(qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.4_LHS1, qfreq_pruned.rules_lhs2_conf50_nosubset_clust4.4_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

### subset to holm sigtrend from lhs1 ####
#subset to rules that contain one of the significat single predictors
pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = subset(pruned.rules_lhs2_conf50_all, subset=lhs %in% lhs1_conf33_vars)

#sort
pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = sort(pruned.rules_lhs2_conf50_lhs1_holm_sigtrend, by=c("fisher.p","fisher.padj_holm"), decreasing = FALSE)

#convert rules to dataframe
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = DATAFRAME(pruned.rules_lhs2_conf50_lhs1_holm_sigtrend, separate = TRUE, setStart = '', itemSep = '+', setEnd = '')

#need the following two split the LHS into 2 columns
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = with(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend, 
  cbind(colsplit(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS, "\\+", c('LHS1', 'LHS2')), 
    DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend[2:15]))


#merge with transaction labels dataframe to get category labels for each rule in LHS1
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = merge(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend, transDataLabels, by.x = "LHS1", by.y = "LHS_label")

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend)[17] = "Cat1"

#remove/drop unused levels of Cat1
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$Cat1 = factor(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$Cat1)

#merge with transaction labels dataframe to get category labels for each rule in LHS2
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = merge(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend, transDataLabels, by.x = "LHS2", by.y = "LHS_label")

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend)[18] = "Cat2"

#remove/drop unused levels of Cat2
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$Cat2 = factor(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$Cat2)

#make crosstab of rule categories
xtabs_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend = xtabs(~Cat1 + Cat2, data = DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend)


#--need to limit data to just those with one of the significant predictors first
DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$OR_ExceedCI = 
  ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cfq26=Disagree" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cfq26=Disagree", 
    ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[1], "Y"," N"), 
    ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cfq28=Disagree" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cfq28=Disagree", 
      ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[2], "Y", "N"),
      ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cebq34_gRarely" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cebq34_gRarely", 
        ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[3], "Y", "N"),
        ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cfq23_lsNeutralAD" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cfq23_lsNeutralAD", 
          ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[4], "Y", "N"),
          ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cebq12_gRarely" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cebq12_gRarely", 
            ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[5], "Y", "N"),
            ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cebq7_gSometimes" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cebq7_gSometimes", 
              ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[6], "Y", "N"),
              ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cfq20=Agree" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cfq20=Agree", 
                ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[7], "Y", "N"),
                ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cebq33_gSometimes" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cebq33_gSometimes", 
                  ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[8], "Y", "N"),
                  ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "mEducation_HS" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "mEducation_HS", 
                    ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[9], "Y", "N"),
                    ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS1 == "cfq29=Always" | DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$LHS2 == "cfq29=Always", 
                      ifelse(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend$oddsRatio > OR_CI_pruned.rules_lhs1_conf33$OR_upperCI[10], "Y", "N"), NA))))))))))


write.csv(DataFrame_pruned.rules_lhs2_conf50_lhs1_holm_sigtrend, file = "ResultsOutput/indQ_LOC-YES_lhs2_conf50_additive.csv", row.names = FALSE)
      
### Logits to look for interaction for rules with 2 items and OR outside of 1 item OR 95% bound ####

## CEBQ7_gSometimesy ####
LOC_arules$loc1 = factor(LOC_arules$loc1)

CEBQ7_logit = glm(loc1~cebq7_gSometimes, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ7_sum = summary(CEBQ7_logit)
CEBQ7_coef = coef(CEBQ7_sum)
CEBQ7_OR = exp(coef(CEBQ7_sum))
CEBQ7_tab = data.frame(CEBQ7_coef, CEBQ7_OR)
colnames(CEBQ7_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ7_gSometimesy+CEB1_gSometimes
LOC_arules$cfq26_Disagree = ifelse(LOC_arules$cebq1_gRarely == "Disagree", TRUE, FALSE)

CEBQ7.CEBQ1_logit = glm(loc1~cebq7_gSometimes*cebq1_gSometimes, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ7.CEBQ1_sum = summary(CEBQ7.CEBQ1_logit)
CEBQ7.CEBQ1_coef = coef(CEBQ7.CEBQ1_sum)[2:4, ]
CEBQ7.CEBQ1_OR = exp(coef(CEBQ7.CEBQ1_sum)[2:4,1:2])
CEBQ7.CEBQ1_tab = data.frame(CEBQ7.CEBQ1_coef, CEBQ7.CEBQ1_OR)
colnames(CEBQ7.CEBQ1_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ7_gSometimes+CFQ28_Disagree - no interaction 
LOC_arules$cfq28_lsNeutralAD = ifelse(LOC_arules$cfq28 == "Disagree" | LOC_arules$cfq28 == "SlightlyDis", TRUE, FALSE)

CEBQ7.CFQ28_logit = glm(loc1~cebq7_gSometimes*cfq28_lsNeutralAD, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ7.CFQ28_sum = summary(CEBQ7.CFQ28_logit)
CEBQ7.CFQ28_coef = coef(CEBQ7.CFQ28_sum)[2:4, ]
CEBQ7.CFQ28_OR = exp(coef(CEBQ7.CFQ28_sum)[2:4,1:2])
CEBQ7.CFQ28_tab = data.frame(CEBQ7.CFQ28_coef, CEBQ7.CFQ28_OR)
colnames(CEBQ7.CFQ28_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ33_gSometimes ####
CEBQ33_logit = glm(loc1~cebq33_gSometimes, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ33_sum = summary(CEBQ33_logit)
CEBQ33_coef = coef(CEBQ33_sum)
CEBQ33_OR = exp(coef(CEBQ33_sum))
CEBQ33_tab = data.frame(CEBQ33_coef, CEBQ33_OR)
colnames(CEBQ33_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ33_gSometimes+BF_ls7mo
CEBQ33.BFls7mo_logit = glm(loc1~cebq33_gSometimes*BreastFed_ls7mo, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ33.BFls7mo_sum = summary(CEBQ33.BFls7mo_logit)
CEBQ33.BFls7mo_coef = coef(CEBQ33.BFls7mo_sum)[2:4, ]
CEBQ33.BFls7mo_OR = exp(coef(CEBQ33.BFls7mo_sum)[2:4,1:2])
CEBQ33.BFls7mo_tab = data.frame(CEBQ33.BFls7mo_coef, CEBQ33.BFls7mo_OR)
colnames(CEBQ33.BFls7mo_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ33_gSometimes+CEB1_gSometimes
CEBQ33.CEBQ1_logit = glm(loc1~cebq33_gSometimes*cebq1_gSometimes, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ33.CEBQ1_sum = summary(CEBQ33.CEBQ1_logit)
CEBQ33.CEBQ1_coef = coef(CEBQ33.CEBQ1_sum)[2:4, ]
CEBQ33.CEBQ1_OR = exp(coef(CEBQ33.CEBQ1_sum)[2:4,1:2])
CEBQ33.CEBQ1_tab = data.frame(CEBQ33.CEBQ1_coef, CEBQ33.CEBQ1_OR)
colnames(CEBQ33.CEBQ1_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ33_gSometimes + CFQ28_Disagree
LOC_arules$cfq28_Disagree = ifelse(LOC_arules$cfq28 == "Disagree", TRUE, FALSE)
CFQ28.CEBQ33_logit = glm(loc1~cebq33_gSometimes*cfq28_Disagree, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CFQ28.CEBQ33_sum = summary(CFQ28.CEBQ33_logit)
CFQ28.CEBQ33_coef = coef(CFQ28.CEBQ33_sum)[2:4, ]
CFQ28.CEBQ33_OR = exp(coef(CFQ28.CEBQ33_sum)[2:4,1:2])
CFQ28.CEBQ33_tab = data.frame(CFQ28.CEBQ33_coef, CFQ28.CEBQ33_OR)
colnames(CFQ28.CEBQ33_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ34_gRarely ####
CEBQ34_logit = glm(loc1~cebq34_gRarely, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ34_sum = summary(CEBQ34_logit)
CEBQ34_coef = coef(CEBQ34_sum)
CEBQ34_OR = exp(coef(CEBQ34_sum))
CEBQ34_tab = data.frame(CEBQ34_coef, CEBQ34_OR)
colnames(CEBQ34_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ34_gRarely+cfq27_lsNeutralAD
LOC_arules$cfq27_lsNeutralAD = ifelse(LOC_arules$cfq27 == "Disagree" | LOC_arules$cfq27 == "SlightlyDis", TRUE, FALSE)
CEBQ34.CFQ27_logit = glm(loc1~cebq34_gRarely*cfq27_lsNeutralAD, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ34.CFQ27_sum = summary(CEBQ34.CFQ27_logit)
CEBQ34.CFQ27_coef = coef(CEBQ34.CFQ27_sum)[2:4, ]
CEBQ34.CFQ27_OR = exp(coef(CEBQ34.CFQ27_sum)[2:4,1:2])
CEBQ34.CFQ27_tab = data.frame(CEBQ34.CFQ27_coef, CEBQ34.CFQ27_OR)
colnames(CEBQ34.CFQ27_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ34_gRarely+CFQ28_Disagree
CEBQ34.CFQ28_logit = glm(loc1~cebq34_gRarely*cfq28_Disagree, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ34.CFQ28_sum = summary(CEBQ34.CFQ28_logit)
CEBQ34.CFQ28_coef = coef(CEBQ34.CFQ28_sum)[2:4, ]
CEBQ34.CFQ28_OR = exp(coef(CEBQ34.CFQ28_sum)[2:4,1:2])
CEBQ34.CFQ28_tab = data.frame(CEBQ34.CFQ28_coef, CEBQ34.CFQ28_OR)
colnames(CEBQ34.CFQ28_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## cfq20_Agree ####
LOC_arules$cfq20_Agree = ifelse(LOC_arules$cfq20 == "Agree", TRUE, FALSE)

CFQ20_logit = glm(loc1~cfq20_Agree, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CFQ20_sum = summary(CFQ20_logit)
CFQ20_coef = coef(CFQ20_sum)
CFQ20_OR = exp(coef(CFQ20_sum))
CFQ20_tab = data.frame(CFQ20_coef, CFQ20_OR)
colnames(CFQ20_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## cfq20_Agree+cfq21_lsNeutralAD 
LOC_arules$cfq21_lsNeutralAD = ifelse(LOC_arules$cfq21 == "Disagree" | LOC_arules$cfq21 == "SlightlyDis", TRUE, FALSE)

CFQ20.CFQ21_logit = glm(loc1~cfq20_Agree*cfq21_lsNeutralAD, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CFQ20.CFQ21_sum = summary(CFQ20.CFQ21_logit)
CFQ20.CFQ21_coef = coef(CFQ20.CFQ21_sum)[2:4, ]
CFQ20.CFQ21_OR = exp(coef(CFQ20.CFQ21_sum)[2:4,1:2])
CFQ20.CFQ21_tab = data.frame(CFQ20.CFQ21_coef, CFQ20.CFQ21_OR)
colnames(CFQ20.CFQ21_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CFQ26_Disagre ####
LOC_arules$cfq26_Disagree = ifelse(LOC_arules$cfq26 == "Disagree", TRUE, FALSE)

CFQ26_logit = glm(loc1~cfq26_Disagree, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CFQ26_sum = summary(CFQ26_logit)
CFQ26_coef = coef(CFQ26_sum)
CFQ26_OR = exp(coef(CFQ26_sum))
CFQ26_tab = data.frame(CFQ26_coef, CFQ26_OR)
colnames(CFQ26_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")

## CEBQ12_gRarely+CFQ26_Disagree - no Interaction 
LOC_arules$cfq26_Disagree = ifelse(LOC_arules$cfq26 == "Disagree", TRUE, FALSE)

CEBQ12.CFQ26_logit = glm(loc1~cebq12_gRarely*cfq26_Disagree, family=binomial(link="logit"), data = LOC_arules[!is.na(LOC_arules$loc1), ])
CEBQ12.CFQ26_sum = summary(CEBQ12.CFQ26_logit)
CEBQ12.CFQ26_coef = coef(CEBQ12.CFQ26_sum)[2:4, ]
CEBQ12.CFQ26_OR = exp(coef(CEBQ12.CFQ26_sum)[2:4,1:2])
CEBQ12.CFQ26_tab = data.frame(CEBQ12.CFQ26_coef, CEBQ12.CFQ26_OR)
colnames(CEBQ12.CFQ26_tab) = c("Beta", "SE", "z", "P", "e^beta", "e^se")


### trace transactions matching rules to participants ####
st_pruned.rules_lhs2_conf50 = supportingTransactions(pruned.rules_lhs2_conf50_all, LOC_arules_trans)
indQ_IDmatch_trans = data.frame(LOC_ndup[!is.na(LOC_ndup$loc1), ]$StudyID, rownames(LOC_arules), LOC_ndup[!is.na(LOC_ndup$loc1), ]$loc1)
indQ_IDmatch_rules = data.frame(indQ_IDmatch_trans, t(as(st_pruned.rules_lhs2_conf50, 'matrix')))
names(indQ_IDmatch_rules)[1:3] = c('StudyID', "TransID", "LOC")
indQ_IDmatch_rules = indQ_IDmatch_rules[indQ_IDmatch_rules$LOC == 'Yes' & !is.na(indQ_IDmatch_rules$LOC), ]

for(c in 4:51){
  indQ_IDmatch_rules[[c]] = ifelse(indQ_IDmatch_rules[[c]] == "TRUE", 1, 0)
}

indQ_IDmatch_rules$nRulesTrue = rowSums(indQ_IDmatch_rules[,4:51])
indQ_n_noRule = nrow(indQ_IDmatch_rules[indQ_IDmatch_rules$nRulesTrue == 0, ])

indQ_endorse_LOC_conf50 = ggplot(indQ_IDmatch_rules, aes(x = nRulesTrue))+
  geom_histogram(aes(y = ..density..), binwidth = 2, color = "grey30", fill = "white") +
  geom_density(alpha = .2, fill = "antiquewhite3")+
  labs(x="Number of Rules Matched")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

###############################################
####      LOC-no; Single category Questionnaire Variables Included      ####
###############################################
###loc-no 1, sup = 0.19 (n = 30), conf=33####
#get rules with only 1 predictor
rules_lhs1_conf33_no = apriori(LOC_arules_trans, parameter = list(maxlen = 2, supp = 0.19, conf = 0.33), 
  appearance = list(default="lhs", rhs="loc1=No"))

#prune rules that are redundant-more important for multiple predictors as may get two rules with 
#different predictor order, thus are the same
pruned.rules_lhs1_conf33_no = rules_lhs1_conf33_no[!is.redundant(rules_lhs1_conf33_no)]

#add odds ratio,chi-square (need to specify also want pvalue) and fisher exact p test to the rules information
quality(pruned.rules_lhs1_conf33_no) = cbind(quality(pruned.rules_lhs1_conf33_no), 
  oddsRatio = interestMeasure(pruned.rules_lhs1_conf33_no, 
    measure = "oddsRatio", 
    transactions = LOC_arules_trans),
  addVal = interestMeasure(pruned.rules_lhs1_conf33_no, 
    measure = "addedValue", 
    transactions = LOC_arules_trans),
  kappa = interestMeasure(pruned.rules_lhs1_conf33_no, 
    measure = "kappa", 
    transactions = LOC_arules_trans),
  fisher.p = interestMeasure(pruned.rules_lhs1_conf33_no, 
    measure = "fishersExactTest",
    transactions = LOC_arules_trans))

#remove empty rules
pruned.rules_lhs1_conf33_no = subset(pruned.rules_lhs1_conf33_no, subset = size(lhs(pruned.rules_lhs1_conf33_no))!=0)

#check added value and correlation
pruned.rules_lhs1_conf33_no_qdat = data.frame(quality(pruned.rules_lhs1_conf33_no))
nrow(pruned.rules_lhs1_conf33_no_qdat[pruned.rules_lhs1_conf33_no_qdat$kappa >=0.21, ])
nrow(pruned.rules_lhs1_conf33_no_qdat[pruned.rules_lhs1_conf33_no_qdat$addVal >=0.05, ])

pruned.rules_lhs1_conf33_no_red = subset(pruned.rules_lhs1_conf33_no, quality(pruned.rules_lhs1_conf33_no)[8]>=0.21)
pruned.rules_lhs1_conf33_no_red = subset(pruned.rules_lhs1_conf33_no_red, quality(pruned.rules_lhs1_conf33_no_red)[7]>=0.05)

#add holm correction to fisher p-values
quality(pruned.rules_lhs1_conf33_no_red) = cbind(quality(pruned.rules_lhs1_conf33_no_red), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs1_conf33_no_red)$fisher.p, method = 'holm'), 4))

#sort rules by pvalues
pruned.rules_lhs1_conf33_no_red = sort(pruned.rules_lhs1_conf33_no_red, by=c("fisher.p", "fisher.padj_holm"), decreasing = FALSE)

#make list of significant variables to reduce multiple predictor rules later on (when we test for additive effects)
lhs1_conf33_no_vars_holm_sig = c("cebq34_lsSometimes")

#get confidence intervals for odds ratio of rules
OR_CI_pruned.rules_lhs1_conf33_no = ARM_ORconf(pruned.rules_lhs1_conf33_no_red, LOC_arules, 1, "No")

#get odds ratio upper confidence bound to use for comparision with multiple predictors odds ratios
OR_upperCI_pruned.rules_lhs1_conf33_no = OR_CI_pruned.rules_lhs1_conf33_no$OR_upperCI[1]

#add CIs to rules qualities and re-order columns nicely
quality(pruned.rules_lhs1_conf33_no_red) = cbind(quality(pruned.rules_lhs1_conf33_no_red),
  OR_CI_pruned.rules_lhs1_conf33_no[3:4])
quality(pruned.rules_lhs1_conf33_no_red) = quality(pruned.rules_lhs1_conf33_no_red)[c(5, 1:4, 7:8, 6, 11:12, 9:10)]

#convert rules into a data frame to look at/use xtable with in markdown document
DataFrame_pruned.rules_lhs1_conf33_no_sig = DATAFRAME(pruned.rules_lhs1_conf33_no_red, separate = TRUE, setStart = '', itemSep = ' + ', setEnd = '')

#merge with transaction labels so each rule gets a category assignment
DataFrame_pruned.rules_lhs1_conf33_no_sig = merge(DataFrame_pruned.rules_lhs1_conf33_no_sig, transDataLabels, by.x = "LHS", by.y = "LHS_label")
DataFrame_pruned.rules_lhs1_conf33_no_sig$Cat = factor(DataFrame_pruned.rules_lhs1_conf33_no_sig$Cat)

#write out data
write.csv(DataFrame_pruned.rules_lhs1_conf33_no_sig, file = "ResultsOutput/indQ_LOC-NO_lhs1_conf33.csv", row.names = FALSE)

###loc-no 2, sup = 0.255 (1/3 no, n = 40), conf=50####
#get the rules ####
rules_lhs2_conf50_no = apriori(LOC_arules_trans, parameter = list(maxlen = 3, supp = 0.255, conf = 0.50), appearance = list(default="lhs", rhs="loc1=No"))

#prune rules that are redundant-more important for multiple predictors as may get two rules with 
#different predictor order, thus are the same
pruned.rules_lhs2_conf50_no = rules_lhs2_conf50_no[!is.redundant(rules_lhs2_conf50_no)]

#add odds ratio, chi-square and its pvalue, and fisher to the rules information
quality(pruned.rules_lhs2_conf50_no) = cbind(quality(pruned.rules_lhs2_conf50_no), 
  oddsRatio = interestMeasure(pruned.rules_lhs2_conf50_no, 
    measure = "oddsRatio", 
    transactions = LOC_arules_trans),
  addVal = interestMeasure(pruned.rules_lhs2_conf50_no, 
    measure = "addedValue", 
    transactions = LOC_arules_trans),
  kappa = interestMeasure(pruned.rules_lhs2_conf50_no, 
    measure = "kappa", 
    transactions = LOC_arules_trans),
  fisher.p = interestMeasure(pruned.rules_lhs2_conf50_no, 
    measure = "fishersExactTest",
    transactions = LOC_arules_trans))

#no subset--all rules ####
#remove empty rules
pruned.rules_lhs2_conf50_no_all = subset(pruned.rules_lhs2_conf50_no, subset = size(lhs(pruned.rules_lhs2_conf50_no))!=0)

#check added value and correlation
pruned.rules_lhs2_conf50_no_all_qdat = data.frame(quality(pruned.rules_lhs2_conf50_no_all))
nrow(pruned.rules_lhs2_conf50_no_all_qdat[pruned.rules_lhs2_conf50_no_all_qdat$kappa >=0.21, ])
nrow(pruned.rules_lhs2_conf50_no_all_qdat[pruned.rules_lhs2_conf50_no_all_qdat$addVal >=0.05, ])

pruned.rules_lhs2_conf50_no_red = subset(pruned.rules_lhs2_conf50_no_all, quality(pruned.rules_lhs2_conf50_no_all)[8]>=0.21)
pruned.rules_lhs2_conf50_no_red = subset(pruned.rules_lhs2_conf50_no_red, quality(pruned.rules_lhs2_conf50_no_red)[7]>=0.05)

#add adjusted pvalues to rules information--adjusting for all rules here, not just subset
quality(pruned.rules_lhs2_conf50_no_red) = cbind(quality(pruned.rules_lhs2_conf50_no_red), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs2_conf50_no_red)$fisher.p, method = 'holm'), 4))

#sort
pruned.rules_lhs2_conf50_no_red = sort(pruned.rules_lhs2_conf50_no_red, by=c("fisher.p", "fisher.padj_holm"), decreasing = FALSE)

#get odds ration confidence intervals for those with significant fisher 
OR_CI_pruned.rules_lhs2_conf50_no = ARM_ORconf(pruned.rules_lhs2_conf50_no_red, LOC_arules, 2, "No")

quality(pruned.rules_lhs2_conf50_no_red) = cbind(quality(pruned.rules_lhs2_conf50_no_red),
  OR_CI_pruned.rules_lhs2_conf50_no[4:5])
quality(pruned.rules_lhs2_conf50_no_red) = quality(pruned.rules_lhs2_conf50_no_red)[c(5, 1:4, 7:8, 6, 11:12, 9:10)]

### clusters loc-no 2, conf=.50, sup = 0.23 (1/3) to ####
#get matrix information for final set
im_lhs2_conf50_no_nosubset = quality(pruned.rules_lhs2_conf50_no_red)

#clustering Rules lhs2 conf50_no lhs1 sig chi
dist_pruned.rules_lhs2_conf50_no = dissimilarity(pruned.rules_lhs2_conf50_no_red,args = list(transactions = LOC_arules_trans), method = "gupta")

## mediods partitioning--determine number of clusters ####
#within cluster sums of squares--look for elbow in graph
pruned.rules_lhs2_conf50_no_nclustSSW = fviz_nbclust(as.matrix(dist_pruned.rules_lhs2_conf50_no), cluster::pam, method = "wss", k.max=10)+ geom_vline(xintercept = 2, linetype = 3)+labs(subtitle = "Elbow method")

#silhouette width
pruned.rules_lhs2_conf50_no_nclustSil = fviz_nbclust(as.matrix(dist_pruned.rules_lhs2_conf50_no), cluster::pam, method = "silhouette", k.max=10)+ labs(subtitle = "Silhouette method")

#get combind plot
pruned.rules_lhs2_conf50_no_nclustPlotGrid = plot_grid(pruned.rules_lhs2_conf50_no_nclustSSW, pruned.rules_lhs2_conf50_no_nclustSil, labels = c("", ""))

#test internal validity metrics
pruned.rules_lhs2_conf50_no_nclust_internV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50_no), 2:10, clMethods=c("pam"), validation="internal")

#test stability metrics (remove one collumn at a time sequentially and look at stability)
pruned.rules_lhs2_conf50_no_nclust_stabV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50_no), 2:10, clMethods=c("pam"), validation="stability")

#cluster validity table

pruned.rules_lhs2_conf50_no_nclust_Vsum = data.frame(matrix(c("Connectivity",  5.7579,    "pam", 2,
  "Dunn",         'inf',    "pam", 10,
  "Silhouette",   0.6516,    "pam", 2,
  "APN",          0.0000,    "pam", 5,
  "AD",           0.0000,    "pam", 10,
  "ADM",          0.0000,    "pam", 5,
  "FOM",          0.0000,    "pam", 10), 
  byrow = TRUE, nrow = 7, ncol = 4))
names(pruned.rules_lhs2_conf50_no_nclust_Vsum) = c("measure", "score_at3", "method", "clusters") 

## run pam to get cluster membership ####
pruned.rules_lhs2_conf50_no_2clust = pam(dist_pruned.rules_lhs2_conf50_no, diss = TRUE, k=2)

#silhouette widths
pruned.rules_lhs2_conf50_no_2clust_asw = data.frame(c(1, 2), pruned.rules_lhs2_conf50_no_2clust$silinfo$clus.avg.widths)
names(pruned.rules_lhs2_conf50_no_2clust_asw) = c("cluster", "asw")

#add cluster belonging to  quality matrix
im_lhs2_conf50_no_nosubset$clust_gupta2 = pruned.rules_lhs2_conf50_no_2clust$clustering

#add to rules quality
quality(pruned.rules_lhs2_conf50_no_red) = cbind(quality(pruned.rules_lhs2_conf50_no_red), Cluster_gupta2 = im_lhs2_conf50_no_nosubset$clust_gupta2)

#convert rules to dataframe
DataFrame_pruned.rules_lhs2_conf50_no_nosubset = DATAFRAME(pruned.rules_lhs2_conf50_no_red, separate = TRUE, setStart = '', itemSep = '+', setEnd = '')

#need the following two split the LHS into 2 columns
DataFrame_pruned.rules_lhs2_conf50_no_nosubset = with(DataFrame_pruned.rules_lhs2_conf50_no_nosubset, 
  cbind(colsplit(DataFrame_pruned.rules_lhs2_conf50_no_nosubset$LHS, "\\+", c('LHS1', 'LHS2')), 
    DataFrame_pruned.rules_lhs2_conf50_no_nosubset[2:15]))


#merge with transaction labels dataframe to get category labels for each rule in LHS1
DataFrame_pruned.rules_lhs2_conf50_no_nosubset = merge(DataFrame_pruned.rules_lhs2_conf50_no_nosubset, transDataLabels, by.x = "LHS1", by.y = "LHS_label")

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_no_nosubset)[17] = "Cat1"

#remove/drop unused levels of Cat1
DataFrame_pruned.rules_lhs2_conf50_no_nosubset$Cat1 = factor(DataFrame_pruned.rules_lhs2_conf50_no_nosubset$Cat1)

#merge with transaction labels dataframe to get category labels for each rule in LHS2
DataFrame_pruned.rules_lhs2_conf50_no_nosubset = merge(DataFrame_pruned.rules_lhs2_conf50_no_nosubset, transDataLabels, by.x = "LHS2", by.y = "LHS_label", all.x = TRUE)

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_no_nosubset)[18] = "Cat2"

#remove/drop unused levels of Cat2
DataFrame_pruned.rules_lhs2_conf50_no_nosubset$Cat2 = factor(DataFrame_pruned.rules_lhs2_conf50_no_nosubset$Cat2)

#make crosstab of rule categories
xtabs_pruned.rules_lhs2_conf50_no_nosubset = xtabs(~Cat1 + Cat2, data = DataFrame_pruned.rules_lhs2_conf50_no_nosubset)

#look at cluster overlap
xtabsClust_pruned.rules_lhs2_conf50_no_nosubset = xtabs(~Cluster_gupta2, data = DataFrame_pruned.rules_lhs2_conf50_no_nosubset)

#write out if want
write.csv(DataFrame_pruned.rules_lhs2_conf50_no_nosubset, file = "ResultsOutput/indQ_LOC-NO_lhs2_conf50_nosubset.csv", row.names = FALSE)

## get crosstables for each cluster--4 cluster solution ####
#--need to first subset and drop unsued levels then make table for each cluster
#cluster 2-1
#subset to just cluster 1
pruned.rules_lhs2_conf50_no_nosubset_clust2.1 = DataFrame_pruned.rules_lhs2_conf50_no_nosubset[DataFrame_pruned.rules_lhs2_conf50_no_nosubset$Cluster_gupta2 == 1, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_no_nosubset_clust2.1$LHS1 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.1$LHS1)
pruned.rules_lhs2_conf50_no_nosubset_clust2.1$LHS2 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.1$LHS2)
pruned.rules_lhs2_conf50_no_nosubset_clust2.1$Cat1 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.1$Cat1)
pruned.rules_lhs2_conf50_no_nosubset_clust2.1$Cat2 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.1$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_no_nosubset_clust2.1))
xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1 = rbind( xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1))
xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1 = cbind(xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.1))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.1_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_no_nosubset_clust2.1))
qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.1_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_no_nosubset_clust2.1))
qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.1 = merge(qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.1_LHS1, qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.1_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 2-2
#subset to just cluster 2
pruned.rules_lhs2_conf50_no_nosubset_clust2.2 = DataFrame_pruned.rules_lhs2_conf50_no_nosubset[DataFrame_pruned.rules_lhs2_conf50_no_nosubset$Cluster_gupta2 == 2, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_no_nosubset_clust2.2$LHS1 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.2$LHS1)
pruned.rules_lhs2_conf50_no_nosubset_clust2.2$LHS2 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.2$LHS2)
pruned.rules_lhs2_conf50_no_nosubset_clust2.2$Cat1 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.2$Cat1)
pruned.rules_lhs2_conf50_no_nosubset_clust2.2$Cat2 = factor(pruned.rules_lhs2_conf50_no_nosubset_clust2.2$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2 = as.matrix(xtabs(~Cat1 + Cat2, data = pruned.rules_lhs2_conf50_no_nosubset_clust2.2))
xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2 = rbind( xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2))
xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2 = cbind(xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_no_nosubset_clust2.2))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.2_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_no_nosubset_clust2.2))
qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.2_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_no_nosubset_clust2.2))
qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.2 = merge(qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.2_LHS1, qfreq_pruned.rules_lhs2_conf50_no_nosubset_clust2.2_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)


write.csv(DataFrame_pruned.rules_lhs2_conf50_no_nosubset, file = "ResultsOutput/indQ_LOC-NO_lhs2_conf50_nosubset_clust.csv", row.names = FALSE)



#***subset to holm sig from lhs1 (based on holm sig trend decision) ####
#get the subset of rules
pruned.rules_lhs2_conf50_no_lhs1_holm_sig = subset(pruned.rules_lhs2_conf50_no_red, subset=lhs %in% lhs1_conf33_no_vars_holm_sig)

#get Dataframe
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig = DATAFRAME(pruned.rules_lhs2_conf50_no_lhs1_holm_sig, separate = TRUE, setStart = '', itemSep = '+', setEnd = '')

#need the following two split the LHS into 2 columns
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig = with(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig, 
  cbind(colsplit(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$LHS, "\\+", c('LHS1', 'LHS2')), 
    DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig[2:15]))

#this is because a single predictor was repeated in LHS2 during the split and it shouldn't have
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$LHS2[5] = NA

#merge
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig = merge(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig, transDataLabels, by.x = "LHS1", by.y = "LHS_label")
names(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig)[17] = "Cat1"
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$Cat1 = factor(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$Cat1)
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig = merge(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig, transDataLabels, by.x = "LHS2", by.y = "LHS_label")
names(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig)[18] = "Cat2"
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$Cat2 = factor(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$Cat2)

#xtabs
xtab_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig = xtabs(~Cat1 + Cat2, data = DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig)

#compare Odds ratio to upper CI of single predictor 
DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$OR_ExceedCI = ifelse(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig$oddsRatio > OR_upperCI_pruned.rules_lhs1_conf33_no, "Y", "N")


#reduce to final ones
write.csv(DataFrame_pruned.rules_lhs1_lhs2_conf50_no_lhs1_holm_sig, file = "ResultsOutput/indQ_LOC-NO_lhs2_conf50_additive.csv", row.names = FALSE)

### trace transactions matching rules to participants ####
st_pruned.rules_lhs2_conf50_no = supportingTransactions(pruned.rules_lhs2_conf50_no_red, LOC_arules_trans)
indQ_IDmatch_trans_no = data.frame(LOC_ndup[!is.na(LOC_ndup$loc1), ]$StudyID, rownames(LOC_arules), LOC_ndup[!is.na(LOC_ndup$loc1), ]$loc1)
indQ_IDmatch_rules_no = data.frame(indQ_IDmatch_trans_no, t(as(st_pruned.rules_lhs2_conf50_no, 'matrix')))
names(indQ_IDmatch_rules_no)[1:3] = c('StudyID', "TransID", "LOC")
indQ_IDmatch_rules_no = indQ_IDmatch_rules_no[indQ_IDmatch_rules_no$LOC == 'No' & !is.na(indQ_IDmatch_rules_no$LOC), ]

for(c in 4:14){
  indQ_IDmatch_rules_no[[c]] = ifelse(indQ_IDmatch_rules_no[[c]] == "TRUE", 1, 0)
}

indQ_IDmatch_rules_no$nRulesTrue = rowSums(indQ_IDmatch_rules_no[,4:14])
indQ_n_noRule = nrow(indQ_IDmatch_rules_no[indQ_IDmatch_rules_no$nRulesTrue == 0, ])

indQ_endorse_LOC_conf50_no = ggplot(indQ_IDmatch_rules_no, aes(x = nRulesTrue))+
  geom_histogram(aes(y = ..density..), binwidth = 2, color = "grey30", fill = "white") +
  geom_density(alpha = .2, fill = "antiquewhite3")+
  labs(x="Number of Rules Matched")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

###############################################
####      High Risk Sample Single category Questionnaire Variables Included      ####
####      ALL LOC STEP 1 INCLUDED (CFQ 26 = ‘Disagree’, CFQ 20 = ‘Agree’, and CEBQ 34 > ‘Rarely’)      ####
###############################################

##HighRisk/Protective (n=70), loc-no (n = 42, 35% of no LOC sample), sup = 0.15 (25%=11), conf=33####
#subset to those that had opposites predicting no LOC: cebq34_gRarely, mEducation_HS, cfq23_lsNeutralAD
LOC_ndup$HR = ifelse(!is.na(LOC_ndup$cebq34_gRarely) & LOC_ndup$cebq34_gRarely == TRUE, 'Yes',
  ifelse(!is.na(LOC_ndup$mEducation_HS) & LOC_ndup$mEducation_HS == TRUE, 'Yes',
    ifelse(!is.na(LOC_ndup$cfq23_lsNeutralAD) & LOC_ndup$cfq23_lsNeutralAD == TRUE, 'Yes', 'No')))
LOC_ndup_HRall = LOC_ndup[LOC_ndup$HR == 'Yes', ]


LOC_arules$HR = ifelse(!is.na(LOC_arules$cebq34_gRarely) & LOC_arules$cebq34_gRarely == TRUE, 'Yes',
  ifelse(!is.na(LOC_arules$mEducation_HS) & LOC_arules$mEducation_HS == TRUE, 'Yes',
    ifelse(!is.na(LOC_arules$cfq23_lsNeutralAD) & LOC_arules$cfq23_lsNeutralAD == TRUE, 'Yes', 'No')))

LOC_arules_HRall = LOC_arules[LOC_arules$HR == 'Yes', ]


#convert data to transactions (dichotomize)
LOC_arules_trans_HRall = as(LOC_arules_HRall[1:287], "transactions")

#get rules
rules_lhs1_conf33_HRallno = apriori(LOC_arules_trans_HRall, parameter = list(maxlen = 2, supp = 0.15, conf = 0.33), 
  appearance = list(default="lhs", rhs="loc1=No"))

#prune rules that are redundant-more important for multiple predictors as may get two rules with 
#different predictor order, thus are the same
pruned.rules_lhs1_conf33_HRallno = rules_lhs1_conf33_HRallno[!is.redundant(rules_lhs1_conf33_HRallno)]

#add odds ratio, chi-squre and its p, and fisher to rules information
quality(pruned.rules_lhs1_conf33_HRallno) = cbind(quality(pruned.rules_lhs1_conf33_HRallno), 
  oddsRatio = interestMeasure(pruned.rules_lhs1_conf33_HRallno, 
    measure = "oddsRatio", 
    transactions = LOC_arules_trans_HRall),
  addVal = interestMeasure(pruned.rules_lhs1_conf33_HRallno, 
    measure = "addedValue", 
    transactions = LOC_arules_trans_HRall),
  kappa = interestMeasure(pruned.rules_lhs1_conf33_HRallno, 
    measure = "kappa", 
    transactions = LOC_arules_trans_HRall),
  fisher.p = interestMeasure(pruned.rules_lhs1_conf33_HRallno, 
    measure = "fishersExactTest",
    transactions = LOC_arules_trans_HRall))

#remove empty rules
pruned.rules_lhs1_conf33_HRallno = subset(pruned.rules_lhs1_conf33_HRallno, subset = size(lhs(pruned.rules_lhs1_conf33_HRallno))!=0)

#check added value and correlation
pruned.rules_lhs1_conf33_HRallno_qdat = data.frame(quality(pruned.rules_lhs1_conf33_HRallno))
nrow(pruned.rules_lhs1_conf33_HRallno_qdat[pruned.rules_lhs1_conf33_HRallno_qdat$kappa >=0.21, ])
nrow(pruned.rules_lhs1_conf33_HRallno_qdat[pruned.rules_lhs1_conf33_HRallno_qdat$addVal >=0.05, ])

pruned.rules_lhs1_conf33_HRallno_red = subset(pruned.rules_lhs1_conf33_HRallno, quality(pruned.rules_lhs1_conf33_HRallno)[8]>=0.21)
pruned.rules_lhs1_conf33_HRallno_red = subset(pruned.rules_lhs1_conf33_HRallno_red, quality(pruned.rules_lhs1_conf33_HRallno_red)[7]>=0.05)


#add multiple comparison correction
quality(pruned.rules_lhs1_conf33_HRallno_red) = cbind(quality(pruned.rules_lhs1_conf33_HRallno_red), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs1_conf33_HRallno_red)$fisher.p, method = 'holm'), 4))

#sort rules
pruned.rules_lhs1_conf33_HRallno_red = sort(pruned.rules_lhs1_conf33_HRallno_red, by=c("fisher.padj_holm", "fisher.p"), decreasing = FALSE)

#get confidence intervals for odds ratios
OR_CI_pruned.rules_lhs1_conf33_HRallno = ARM_ORconf(pruned.rules_lhs1_conf33_HRallno_red, LOC_arules_HRall, 1, "No")

#add confidence intervals to rules information
quality(pruned.rules_lhs1_conf33_HRallno_red) = cbind(quality(pruned.rules_lhs1_conf33_HRallno_red),
  OR_CI_pruned.rules_lhs1_conf33_HRallno[3:4])

#organize columns
quality(pruned.rules_lhs1_conf33_HRallno_red) = quality(pruned.rules_lhs1_conf33_HRallno_red)[c(5, 1:4, 7:8, 6, 11:12, 9:10)]

#convert rules into a data frame to look at/use xtable with in markdown document
DataFrame_pruned.rules_lhs1_conf33_HRallno_sig = DATAFRAME(pruned.rules_lhs1_conf33_HRallno_red, separate = TRUE, setStart = '', itemSep = ' + ', setEnd = '')

#merge with transaction labels so each rule gets a category assignment
DataFrame_pruned.rules_lhs1_conf33_HRallno_sig = merge(DataFrame_pruned.rules_lhs1_conf33_HRallno_sig, transDataLabels, by.x = "LHS", by.y = "LHS_label")
DataFrame_pruned.rules_lhs1_conf33_HRallno_sig$Cat = factor(DataFrame_pruned.rules_lhs1_conf33_HRallno_sig$Cat)

#get crosstab table for rule categories to summarize rules
xtab_pruned.rules_lhs1_conf33_no_sig = xtabs(~Cat, data = DataFrame_pruned.rules_lhs1_conf33_HRallno_sig)
write.csv(DataFrame_pruned.rules_lhs1_conf33_HRallno_sig, file = "ResultsOutput/indQ_LOC-HR_NO_lhs1_conf33.csv", row.names = TRUE)


##HighRisk/Protective (n=70), has at least 2; loc-no2 (42), sup = 0.20 (33%=14), conf=50####
rules_lhs2_conf50_HRallno = apriori(LOC_arules_trans_HRall, parameter = list(maxlen = 3, supp = 0.20, conf = 0.50), 
  appearance = list(default="lhs", rhs="loc1=No"))

pruned.rules_lhs2_conf50_HRallno = rules_lhs2_conf50_HRallno[!is.redundant(rules_lhs2_conf50_HRallno)]

quality(pruned.rules_lhs2_conf50_HRallno) = cbind(quality(pruned.rules_lhs2_conf50_HRallno), 
  oddsRatio = interestMeasure(pruned.rules_lhs2_conf50_HRallno, 
    measure = "oddsRatio", 
    transactions = LOC_arules_trans_HRall),
  addVal = interestMeasure(pruned.rules_lhs2_conf50_HRallno, 
    measure = "addedValue", 
    transactions = LOC_arules_trans_HRall),
  kappa = interestMeasure(pruned.rules_lhs2_conf50_HRallno, 
    measure = "kappa", 
    transactions = LOC_arules_trans_HRall),
  fisher.p = interestMeasure(pruned.rules_lhs2_conf50_HRallno, 
    measure = "fishersExactTest",
    transactions = LOC_arules_trans_HRall))

#remove empty rules
pruned.rules_lhs2_conf50_HRallno = subset(pruned.rules_lhs2_conf50_HRallno, subset = size(lhs(pruned.rules_lhs2_conf50_HRallno))!=0)

#check added value and correlation
pruned.rules_lhs2_conf50_HRallno_qdat = data.frame(quality(pruned.rules_lhs2_conf50_HRallno))
nrow(pruned.rules_lhs2_conf50_HRallno_qdat[pruned.rules_lhs2_conf50_HRallno_qdat$kappa >=0.21, ])
nrow(pruned.rules_lhs2_conf50_HRallno_qdat[pruned.rules_lhs2_conf50_HRallno_qdat$addVal >=0.05, ])

pruned.rules_lhs2_conf50_HRallno_red = subset(pruned.rules_lhs2_conf50_HRallno, quality(pruned.rules_lhs2_conf50_HRallno)[8]>=0.21)
pruned.rules_lhs2_conf50_HRallno_red = subset(pruned.rules_lhs2_conf50_HRallno_red, quality(pruned.rules_lhs2_conf50_HRallno_red)[7]>=0.05)


#add multiple comparison to rules information
quality(pruned.rules_lhs2_conf50_HRallno_red) = cbind(quality(pruned.rules_lhs2_conf50_HRallno_red), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs2_conf50_HRallno_red)$fisher.p, method = 'holm'), 4))

#sort rules
pruned.rules_lhs2_conf50_HRallno_red = sort(pruned.rules_lhs2_conf50_HRallno_red, by=c("fisher.padj_holm", "fisher.p"), decreasing = FALSE)

#reduce to significant
pruned.rules_lhs2_conf50_HRallno_redsig = pruned.rules_lhs2_conf50_HRallno_red[1:226]

#get confidence intervals for odds ratios
OR_CI_pruned.rules_lhs2_conf50_HRallno = ARM_ORconf(pruned.rules_lhs2_conf50_HRallno_redsig, LOC_arules_HRall, 2, "No")

#add confidence interval information to rules information
quality(pruned.rules_lhs2_conf50_HRallno_redsig) = cbind(quality(pruned.rules_lhs2_conf50_HRallno_redsig),
  OR_CI_pruned.rules_lhs2_conf50_HRallno[4:5])

#organize columns
quality(pruned.rules_lhs2_conf50_HRallno_redsig) = quality(pruned.rules_lhs2_conf50_HRallno_redsig)[c(5, 1:4, 7:8, 6, 11:12, 9:10)]

### High Risk clusters loc-no 2, conf=.50, sup = 0.23 (1/3) to ####
#get matrix information for final set (holm sig/trend)
im_matrix_lhs2_conf50_HRallno = quality(pruned.rules_lhs2_conf50_HRallno_redsig)
row.names(im_matrix_lhs2_conf50_HRallno) = seq(1, nrow(im_matrix_lhs2_conf50_HRallno), by=1)

#clustering Rules lhs2 conf50 lhs1 sig chi
dist_pruned.rules_lhs2_conf50_HRallno_gupta = dissimilarity(pruned.rules_lhs2_conf50_HRallno_redsig, args = list(transactions = LOC_arules_trans_HRall), method = "gupta")

## mediods partitioning--determine number of clusters ####
#within cluster sums of squares--look for elbow in graph
pruned.rules_lhs2_conf50_noHR_nclustSSW = fviz_nbclust(as.matrix(dist_pruned.rules_lhs2_conf50_HRallno_gupta), cluster::pam, method = "wss", k.max=15)+ geom_vline(xintercept =6, linetype = 3)+ geom_vline(xintercept = 8, linetype = 3)+labs(subtitle = "Elbow method")

#silhouette width
pruned.rules_lhs2_conf50_noHR_nclustSil = fviz_nbclust(as.matrix(dist_pruned.rules_lhs2_conf50_HRallno_gupta), cluster::pam, method = "silhouette", k.max=15)+ labs(subtitle = "Silhouette method")

#get combind plot
pruned.rules_lhs2_conf50_noHR_nclustPlotGrid = plot_grid(pruned.rules_lhs2_conf50_noHR_nclustSSW, pruned.rules_lhs2_conf50_noHR_nclustSil, labels = c("", ""))

#test internal validity metrics
#pruned.rules_lhs2_conf50_noHR_nclust_internV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50_HRallno_gupta), 2:10, clMethods=c("pam"), validation="internal")

#test stability metrics (remove one collumn at a time sequentially and look at stability)
#pruned.rules_lhs2_conf50_noHR_nclust_stabV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50_HRallno_gupta), 2:10, clMethods=c("pam"), validation="stability")

#cluster validity table
pruned.rules_lhs2_conf50_noHR_nclust_Vsum = data.frame(matrix(c("Connectivity", 18.1659,    "pam",      2,
  "Dunn",         0.2245,    "pam",       6,
  "Silhouette",   0.3262,    "pam",       8,
  "APN",          0.0090,    "pam",       2,
  "AD",           1.5123,    "pam",       10,
  "ADM",          0.0261,    "pam",       9,
  "FOM",          0.0785,    "pam",       10), 
  byrow = TRUE, nrow = 7, ncol = 4))
names(pruned.rules_lhs2_conf50_noHR_nclust_Vsum) = c("measure", "score", "method", "clusters")

## run pam to get cluster membership ####
pruned.rules_lhs2_conf50_noHR_8clust = pam(dist_pruned.rules_lhs2_conf50_HRallno_gupta, diss = TRUE, k=8)

#silhouette widths
pruned.rules_lhs2_conf50_noHR_8clust_asw = data.frame(c(1, 2, 3, 4, 5, 6, 7, 8), pruned.rules_lhs2_conf50_noHR_8clust$silinfo$clus.avg.widths)
names(pruned.rules_lhs2_conf50_noHR_8clust_asw) = c("cluster", "asw")

#add cluster belonging to  quality matrix
im_matrix_lhs2_conf50_HRallno$clust_gupta8 = pruned.rules_lhs2_conf50_noHR_8clust$clustering

#add to rules quality
quality(pruned.rules_lhs2_conf50_HRallno_redsig) = cbind(quality(pruned.rules_lhs2_conf50_HRallno_redsig), Cluster_gupta8 = im_matrix_lhs2_conf50_HRallno$clust_gupta8)

#convert rules to dataframe
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset = DATAFRAME(pruned.rules_lhs2_conf50_HRallno_redsig, separate = TRUE, setStart = '', itemSep = '+', setEnd = '')

#need the following two split the LHS into 2 columns
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset = with(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset, 
  cbind(colsplit(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$LHS, "\\+", c('LHS1', 'LHS2')), 
    DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[2:15]))


#merge with transaction labels dataframe to get category labels for each rule in LHS1
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset = merge(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset, transDataLabels, by.x = "LHS1", by.y = "LHS_label", all.x = TRUE)

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset)[17] = "Cat1"

#remove/drop unused levels of Cat1
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cat1 = factor(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cat1)

#merge with transaction labels dataframe to get category labels for each rule in LHS2
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset = merge(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset, transDataLabels, by.x = "LHS2", by.y = "LHS_label", all.x = TRUE)

#rename variable so can add two category variabls
names(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset)[18] = "Cat2"

#remove/drop unused levels of Cat2
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cat2 = factor(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cat2)

#reorder so can call first row in output - 1 sig after correction
DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[order(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$fisher.padj_holm),]

#make crosstab of rule categories
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset = xtabs(~Cat1 + Cat2, data = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset)

#look at cluster overlap
xtabsClust_pruned.rules_lhs2_conf50_noHR_nosubset = xtabs(~Cluster_gupta8, data = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset)

#write out if want
write.csv(DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset, file = "ResultsOutput/indQ_LOC-NO.HR_lhs2_conf50_nosubset.csv", row.names = FALSE)

## get crosstables for each cluster--5 cluster solution ####
#--need to first subset and drop unsued levels then make table for each cluster
#cluster 8-1
#subset to just cluster 1
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 1, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.1_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-2
#subset to just cluster 2

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 2, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-3
#subset to just cluster 3

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 3, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.3_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-4
#subset to just cluster 4

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 4, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.4_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-5
#subset to just cluster 5

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 5, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.5_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-6
#subset to just cluster 6

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 6, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.6_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-7
#subset to just cluster 7

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 7, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.7_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

#cluster 8-8
#subset to just cluster 8

pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8 = DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset[DataFrame_pruned.rules_lhs2_conf50_noHR_nosubset$Cluster_gupta8 == 8, ]

#drop extra levels not in this cluster
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$LHS1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$LHS1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$LHS2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$LHS2)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$Cat1 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$Cat1)
pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$Cat2 = factor(pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8$Cat2)

#make cross tabs
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8 = as.matrix(xtabs(~Cat2 + Cat1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8 = rbind( xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8, 
  coltotals = colSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8))
xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8 = cbind(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8,
  rowtotals = rowSums(xtabs_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8))

#get frequency for each question
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8_LHS1 = as.data.frame(xtabs(~LHS1, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8_LHS2 = as.data.frame(xtabs(~LHS2, data = pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8))
qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8 = merge(qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.2_LHS1, qfreq_pruned.rules_lhs2_conf50_noHR_nosubset_clust8.8_LHS2, by.x = "LHS1", by.y = "LHS2", all.x = TRUE, all.y = TRUE)

### trace transactions matching rules to participants ####
st_pruned.rules_lhs2_conf50_noHR = supportingTransactions(pruned.rules_lhs2_conf50_HRallno_redsig, LOC_arules_trans_HRall)
indQ_IDmatch_trans_noHR = data.frame(LOC_ndup_HRall[!is.na(LOC_ndup_HRall$loc1), ]$StudyID, rownames(LOC_arules_HRall), LOC_ndup_HRall[!is.na(LOC_ndup_HRall$loc1), ]$loc1)
indQ_IDmatch_rules_noHR = data.frame(indQ_IDmatch_trans_noHR, t(as(st_pruned.rules_lhs2_conf50_noHR, 'matrix')))
names(indQ_IDmatch_rules_noHR)[1:3] = c('StudyID', "TransID", "LOC")
indQ_IDmatch_rules_noHR = indQ_IDmatch_rules_noHR[indQ_IDmatch_rules_noHR$LOC == 'No' & !is.na(indQ_IDmatch_rules_noHR$LOC), ]

for(c in 4:229){
  indQ_IDmatch_rules_noHR[[c]] = ifelse(indQ_IDmatch_rules_noHR[[c]] == "TRUE", 1, 0)
}

indQ_IDmatch_rules_noHR$nRulesTrue = rowSums(indQ_IDmatch_rules_noHR[,4:229])
indQ_n_noRuleHR = nrow(indQ_IDmatch_rules_noHR[indQ_IDmatch_rules_noHR$nRulesTrue == 0, ])

indQ_endorse_LOC_conf50_noHR = ggplot(indQ_IDmatch_rules_noHR, aes(x = nRulesTrue))+
  geom_histogram(aes(y = ..density..), binwidth = 2, color = "grey30", fill = "white") +
  geom_density(alpha = .2, fill = "antiquewhite3")+
  labs(x="Number of Rules Matched - High Risk No LOC")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

###############################################

###############################################

###############################################
