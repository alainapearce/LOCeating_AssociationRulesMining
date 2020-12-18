# This script was written by Alaina Pearce in 2020 
# to binarize data to use in Association Rules Mining 
# Orignial publication in Apetite (2021) entitled 
# 'Using Association Rules Mining to Characterize of 
# Loss of Control Eating in Childhood'
# 
#     Copyright (C) 2021 Alaina L Pearce
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(lsr)

source('functions.R')
LOC_dat = read.csv('Data/LOC_binarized_data.csv')

#### Demographics ####
#get crosstab table for LOC and count number of NA overall
#and for differences between LOC and no LOC
#note: sd.function.na and range.function.na all come from
#the functions.R file so that must be sourced above to work

## LOC
#get crosstab table for LOC and count number of NA
loc_tab = xtabs(~loc1, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$loc1), ])

## Age
age_mean = mean(LOC_dat$cAge_mo, na.rm = TRUE)/12
age_sd = sd(LOC_dat$cAge_mo, na.rm = TRUE)/12
age_range = range(LOC_dat$cAge_mo, na.rm = TRUE)/12
nrow(LOC_dat[is.na(LOC_dat$cAge_mo), ])

age_loc_t = t.test((cAge_mo/12)~loc1, data = LOC_dat)
age_loc_sd = sd.function.na(LOC_dat, LOC_dat$cAge_mo, LOC_dat$loc1)/12
age_loc_range = range.function.na(LOC_dat, (LOC_dat$cAge_mo/12), LOC_dat$loc1)
age_loc_d = cohensD(LOC_dat[LOC_dat$loc1 == 'Yes', ]$cAge_mo, LOC_dat[LOC_dat$loc1 == 'No', ]$cAge_mo)

nrow(LOC_dat[is.na(LOC_dat$cAge_mo) && LOC_dat$loc1 == 'Yes', ])
nrow(LOC_dat[is.na(LOC_dat$cAge_mo) && LOC_dat$loc1 == 'No', ])

## Gender
gender_tab = xtabs(~sex, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$sex), ])

gender_loc_tab = xtabs(~sex + loc1, data=LOC_dat)
gender_loc_chi = chisq.test(gender_loc_tab)
nrow(LOC_dat[which(is.na(LOC_dat$sex) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$sex) && LOC_dat$loc1 == 'No'), ])

## cBMIp
cBMIp_mean = mean(LOC_dat$cBodyMass_p, na.rm = TRUE)
cBMIp_sd = sd(LOC_dat$cBodyMass_p, na.rm = TRUE)
cBMIp_range = range(LOC_dat$cBodyMass_p, na.rm = TRUE)
nrow(LOC_dat[is.na(LOC_dat$cBodyMass_p), ])

cBMIp_loc_t = t.test(cBodyMass_p~loc1, data = LOC_dat)
cBMIp_loc_sd = sd.function.na(LOC_dat, LOC_dat$cBodyMass_p, LOC_dat$loc1)
cBMIp_loc_range = range.function.na(LOC_dat, (LOC_dat$cBodyMass_p), LOC_dat$loc1)
cBMIp_loc_d = cohensD(LOC_dat[LOC_dat$loc1 == 'Yes', ]$cBodyMass_p, LOC_dat[LOC_dat$loc1 == 'No', ]$cBodyMass_p)

nrow(LOC_dat[is.na(LOC_dat$cBodyMass_p) && LOC_dat$loc1 == 'Yes', ])
nrow(LOC_dat[is.na(LOC_dat$cBodyMass_p) && LOC_dat$loc1 == 'No', ])

## BMI class
BMI_tab = xtabs(~cBodyMass_class, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$cBodyMass_class), ])

BMI_loc_tab = xtabs(~cBodyMass_class + loc1, data = LOC_dat)
BMI_loc_fisher = fisher.test(BMI_loc_tab)
BMInoUW_chi = chisq.test(matrix(c(93, 22, 13, 8, 14, 7), byrow = TRUE, ncol = 2))
BMInoUW_loc_fisher = fisher.test(matrix(c(93, 22, 13, 8, 14, 7), byrow = TRUE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$cBodyMass_class) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$cBodyMass_class) & LOC_dat$loc1 == 'No'), ])

##objective measurement of parent height/weight in FBS, EBS, DMK, and Test-Retest studies
HWmeasure = nrow(LOC_dat[grep("DMK*|Test*|FBS*|EBS*", LOC_dat$StudyID), ])
parent_HWmeasure = xtabs(~relationship, data = LOC_dat[grep("DMK*|Test*|FBS*|EBS*", LOC_dat$StudyID), ])

## Race
race_tab = xtabs(~cRace, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$cRace), ])

race_loc_tab = xtabs(~cRace + loc1, data = LOC_dat)
race_loc_fisher = fisher.test(matrix(c(4,31,1,2,114,4), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$cRace) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$cRace) & LOC_dat$loc1 == 'No'), ])

## Ethnicity
ethnicity_tab = xtabs(~cEthnicity, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$cEthnicity), ])

ethnicity_loc_tab = xtabs(~cEthnicity + loc1, data = LOC_dat)
ethnicity_loc_fisher = fisher.test(matrix(c(2,33,5,76), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$cEthnicity) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$cEthnicity) & LOC_dat$loc1 == 'No'), ])

## Income
income_tab = xtabs(~income, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$income), ])

income_loc_tab = xtabs(~income + loc1, data = LOC_dat)
income_loc_fisher = fisher.test(matrix(c(6,21,10,40,54,24), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$income) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$income) & LOC_dat$loc1 == 'No'), ])

## mEducation
mEducation_tab = xtabs(~mEducation, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$mEducation), ])

mEducation_loc_tab = xtabs(~mEducation + loc1, data = LOC_dat)
mEducation_loc_fisher = fisher.test(matrix(c(10,26,11,105), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$mEducation) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$mEducation) & LOC_dat$loc1 == 'No'), ])

## dEducation
dEducation_tab = xtabs(~dEducation, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$dEducation), ])

dEducation_loc_tab = xtabs(~dEducation + loc1, data = LOC_dat)
dEducation_loc_fisher = fisher.test(matrix(c(9, 25, 23, 92), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$dEducation) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$dEducation) & LOC_dat$loc1 == 'No'), ])

## breast fed total
BreastFed_tab = xtabs(~BreastFed, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$BreastFed), ])

BreastFed_loc_tab = xtabs(~BreastFed + loc1, data = LOC_dat)
BreastFed_loc_chi = chisq.test(BreastFed_loc_tab)
nrow(LOC_dat[which(is.na(LOC_dat$BreastFed) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$BreastFed) & LOC_dat$loc1 == 'No'), ])

## breast fed mo
BreastFed_mo_tab = xtabs(~BreastFed_mo, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$BreastFed_mo), ])

BreastFed_mo_loc_tab = xtabs(~BreastFed_mo + loc1, data = LOC_dat)
BreastFed_mo_loc_fisher = fisher.test(matrix(c(10,3,13,5,5,30,12,35,31,11), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$BreastFed_mo) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$BreastFed_mo) & LOC_dat$loc1 == 'No'), ])

## mBMI
mBMI_tab = xtabs(~mBodyMass_class, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$mBodyMass_class), ])

mBMI_loc_tab = xtabs(~mBodyMass_class + loc1, data = LOC_dat)
mBMI_loc_fisher = fisher.test(matrix(c(8,4,5,15,2,11,13,19,57,2), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$mBodyMass_class) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$mBodyMass_class) & LOC_dat$loc1 == 'No'), ])

## dBMI
dBMI_tab = xtabs(~dBodyMass_class, data = LOC_dat)
nrow(LOC_dat[is.na(LOC_dat$dBodyMass_class), ])

dBMI_loc_tab = xtabs(~dBodyMass_class + loc1, data = LOC_dat)
dBMI_loc_fisher = fisher.test(matrix(c(5,7,15,9,0,12,22,46,34,0), byrow = FALSE, ncol = 2))
dBMI_loc_fisher_noUW = fisher.test(matrix(c(5,7,15,9,12,22,46,34), byrow = FALSE, ncol = 2))
nrow(LOC_dat[which(is.na(LOC_dat$dBodyMass_class) & LOC_dat$loc1 == 'Yes'), ])
nrow(LOC_dat[which(is.na(LOC_dat$dBodyMass_class) & LOC_dat$loc1 == 'No'), ])

##### CEBQ ####

#recode CEBQ data as numeric
#create dataset to add numeric CEBQ values
cebq_LOC_dat = LOC_dat[c(2, 52:86)]
#for the 35 cebq questions
for(c in 2:36){
  #create new column
  cebq_LOC_dat[[c+35]] = NA
  
  #get original varname
  varname = names(cebq_LOC_dat)[[c]]
  
  #add "_raw" to varname and set as the name for the new column
  names(cebq_LOC_dat)[[c+35]] = paste(varname, "raw", sep = "_")
  
  #set that needs to be reverse coded
  if(c == 4|c == 5|c == 11|c == 17|c == 33){
    cebq_LOC_dat[c+35] = ifelse(is.na(cebq_LOC_dat[[c]]), NA, ifelse(
      cebq_LOC_dat[[c]] == "Never", 5, ifelse(
        cebq_LOC_dat[[c]] == "Rarely", 4, ifelse(
          cebq_LOC_dat[[c]] == "Sometimes", 3, ifelse(
            cebq_LOC_dat[[c]] == "Often", 2, 1
          )
        )
      )
    )) 
  }
  #non-reverse coded
  else {
    cebq_LOC_dat[c+35] = ifelse(is.na(cebq_LOC_dat[[c]]), NA, ifelse(
      cebq_LOC_dat[[c]] == "Never", 1, ifelse(
        cebq_LOC_dat[[c]] == "Rarely", 2, ifelse(
          cebq_LOC_dat[[c]] == "Sometimes", 3, ifelse(
            cebq_LOC_dat[[c]] == "Often", 4, 5
          )
        )
      )
    )) 
  }  
}
#add loc
cebq_LOC_dat = merge(cebq_LOC_dat, LOC_dat[c(2, 33)], by = 'StudyID')

#get means and descriptives for CEBQ questionnaire
## FR
#numeric mean for subscale for each person as new var
cebq_LOC_dat$cebqFR = rowMeans(data.frame(cebq_LOC_dat$cebq12_raw, cebq_LOC_dat$cebq14_raw, cebq_LOC_dat$cebq19_raw, cebq_LOC_dat$cebq28_raw, cebq_LOC_dat$cebq34_raw), na.rm = T)

#internal reliability
alpha_cebqFR = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq12_raw, cebq_LOC_dat$cebq14_raw, cebq_LOC_dat$cebq19_raw, cebq_LOC_dat$cebq28_raw, cebq_LOC_dat$cebq34_raw)))

cebqFR_mean = mean(cebq_LOC_dat$cebqFR, na.rm = TRUE)
cebqFR_sd = sd(cebq_LOC_dat$cebqFR, na.rm = TRUE)
cebqFR_range = range(cebq_LOC_dat$cebqFR, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqFR), ])

cebqFR_loc_t = t.test(cebqFR~loc1, data = cebq_LOC_dat)
cebqFR_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqFR, cebq_LOC_dat$loc1)
cebqFR_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqFR, cebq_LOC_dat$loc1)
cebqFR_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqFR, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqFR)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqFR) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqFR) && cebq_LOC_dat$loc1 == 'No', ])

## EOE

cebq_LOC_dat$cebqEOE = rowMeans(data.frame(cebq_LOC_dat$cebq2_raw, cebq_LOC_dat$cebq13_raw, cebq_LOC_dat$cebq15_raw, 
  cebq_LOC_dat$cebq23_raw, cebq_LOC_dat$cebq27_raw), na.rm = T)
alpha_cebqEOE = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq2_raw, cebq_LOC_dat$cebq13_raw, cebq_LOC_dat$cebq15_raw, 
  cebq_LOC_dat$cebq23_raw, cebq_LOC_dat$cebq27_raw)))

cebqEOE_mean = mean(cebq_LOC_dat$cebqEOE, na.rm = TRUE)
cebqEOE_sd = sd(cebq_LOC_dat$cebqEOE, na.rm = TRUE)
cebqEOE_range = range(cebq_LOC_dat$cebqEOE, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEOE), ])

cebqEOE_loc_t = t.test(cebqEOE~loc1, data = cebq_LOC_dat)
cebqEOE_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqEOE, cebq_LOC_dat$loc1)
cebqEOE_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqEOE, cebq_LOC_dat$loc1)
cebqEOE_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqEOE, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqEOE)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEOE) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEOE) && cebq_LOC_dat$loc1 == 'No', ])

## EF
cebq_LOC_dat$cebqEF = rowMeans(data.frame(cebq_LOC_dat$cebq1_raw, cebq_LOC_dat$cebq5_raw, cebq_LOC_dat$cebq20_raw, 
  cebq_LOC_dat$cebq22_raw), na.rm = T)
alpha_cebqEF = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq1_raw, cebq_LOC_dat$cebq5_raw, cebq_LOC_dat$cebq20_raw, 
  cebq_LOC_dat$cebq22_raw)))

cebqEF_mean = mean(cebq_LOC_dat$cebqEF, na.rm = TRUE)
cebqEF_sd = sd(cebq_LOC_dat$cebqEF, na.rm = TRUE)
cebqEF_range = range(cebq_LOC_dat$cebqEF, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEF), ])

cebqEF_loc_t = t.test(cebqEF~loc1, data = cebq_LOC_dat)
cebqEF_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqEF, cebq_LOC_dat$loc1)
cebqEF_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqEF, cebq_LOC_dat$loc1)
cebqEF_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqEF, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqEF)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEF) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEF) && cebq_LOC_dat$loc1 == 'No', ])

## DD
cebq_LOC_dat$cebqDD = rowMeans(data.frame(cebq_LOC_dat$cebq6_raw, cebq_LOC_dat$cebq29_raw, cebq_LOC_dat$cebq31_raw), na.rm = T)
alpha_cebqDD = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq6_raw, cebq_LOC_dat$cebq29_raw, cebq_LOC_dat$cebq31_raw)))

cebqDD_mean = mean(cebq_LOC_dat$cebqDD, na.rm = TRUE)
cebqDD_sd = sd(cebq_LOC_dat$cebqDD, na.rm = TRUE)
cebqDD_range = range(cebq_LOC_dat$cebqDD, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqDD), ])

cebqDD_loc_t = t.test(cebqDD~loc1, data = cebq_LOC_dat)
cebqDD_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqDD, cebq_LOC_dat$loc1)
cebqDD_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqDD, cebq_LOC_dat$loc1)
cebqDD_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqDD, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqDD)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqDD) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqDD) && cebq_LOC_dat$loc1 == 'No', ])

## SR
cebq_LOC_dat$cebqSR = rowMeans(data.frame(cebq_LOC_dat$cebq3_raw, cebq_LOC_dat$cebq17_raw, cebq_LOC_dat$cebq21_raw, 
  cebq_LOC_dat$cebq26_raw, cebq_LOC_dat$cebq30_raw), na.rm = T)
alpha_cebqSR = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq3_raw, cebq_LOC_dat$cebq17_raw, cebq_LOC_dat$cebq21_raw, 
  cebq_LOC_dat$cebq26_raw, cebq_LOC_dat$cebq30_raw)))

cebqSR_mean = mean(cebq_LOC_dat$cebqSR, na.rm = TRUE)
cebqSR_sd = sd(cebq_LOC_dat$cebqSR, na.rm = TRUE)
cebqSR_range = range(cebq_LOC_dat$cebqSR, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqSR), ])

cebqSR_loc_t = t.test(cebqSR~loc1, data = cebq_LOC_dat)
cebqSR_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqSR, cebq_LOC_dat$loc1)
cebqSR_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqSR, cebq_LOC_dat$loc1)
cebqSR_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqSR, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqSR)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqSR) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqSR) && cebq_LOC_dat$loc1 == 'No', ])

## SE
cebq_LOC_dat$cebqSE = rowMeans(data.frame(cebq_LOC_dat$cebq4_raw, cebq_LOC_dat$cebq8_raw, cebq_LOC_dat$cebq18_raw, 
  cebq_LOC_dat$cebq35_raw), na.rm = T)
alpha_cebqSE = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq4_raw, cebq_LOC_dat$cebq8_raw, cebq_LOC_dat$cebq18_raw, 
  cebq_LOC_dat$cebq35_raw)))

cebqSE_mean = mean(cebq_LOC_dat$cebqSE, na.rm = TRUE)
cebqSE_sd = sd(cebq_LOC_dat$cebqSE, na.rm = TRUE)
cebqSE_range = range(cebq_LOC_dat$cebqSE, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqSE), ])

cebqSE_loc_t = t.test(cebqSE~loc1, data = cebq_LOC_dat)
cebqSE_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqSE, cebq_LOC_dat$loc1)
cebqSE_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqSE, cebq_LOC_dat$loc1)
cebqSE_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqSE, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqSE)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqSE) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqSE) && cebq_LOC_dat$loc1 == 'No', ])

## EUE
cebq_LOC_dat$cebqEUE = rowMeans(data.frame(cebq_LOC_dat$cebq9_raw, cebq_LOC_dat$cebq11_raw, cebq_LOC_dat$cebq25_raw), na.rm = T)
alpha_cebqEUE = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq9_raw, cebq_LOC_dat$cebq11_raw, cebq_LOC_dat$cebq25_raw)))

cebqEUE_mean = mean(cebq_LOC_dat$cebqEUE, na.rm = TRUE)
cebqEUE_sd = sd(cebq_LOC_dat$cebqEUE, na.rm = TRUE)
cebqEUE_range = range(cebq_LOC_dat$cebqEUE, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEUE), ])

cebqEUE_loc_t = t.test(cebqEUE~loc1, data = cebq_LOC_dat)
cebqEUE_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqEUE, cebq_LOC_dat$loc1)
cebqEUE_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqEUE, cebq_LOC_dat$loc1)
cebqEUE_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqEUE, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqEUE)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEUE) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqEUE) && cebq_LOC_dat$loc1 == 'No', ])

## FF
cebq_LOC_dat$cebqFF = rowMeans(data.frame(cebq_LOC_dat$cebq7_raw, cebq_LOC_dat$cebq10_raw, cebq_LOC_dat$cebq16_raw, 
  cebq_LOC_dat$cebq24_raw, cebq_LOC_dat$cebq32_raw, cebq_LOC_dat$cebq33_raw), na.rm = T)
alpha_cebqFF = summary(psych::alpha(data.frame(cebq_LOC_dat$cebq7_raw, cebq_LOC_dat$cebq10_raw, cebq_LOC_dat$cebq16_raw, 
  cebq_LOC_dat$cebq24_raw, cebq_LOC_dat$cebq32_raw, cebq_LOC_dat$cebq33_raw)))

cebqFF_mean = mean(cebq_LOC_dat$cebqFF, na.rm = TRUE)
cebqFF_sd = sd(cebq_LOC_dat$cebqFF, na.rm = TRUE)
cebqFF_range = range(cebq_LOC_dat$cebqFF, na.rm = TRUE)
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqFF), ])

cebqFF_loc_t = t.test(cebqFF~loc1, data = cebq_LOC_dat)
cebqFF_loc_sd = sd.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqFF, cebq_LOC_dat$loc1)
cebqFF_loc_range = range.function.na(cebq_LOC_dat, cebq_LOC_dat$cebqFF, cebq_LOC_dat$loc1)
cebqFF_loc_d = cohensD(cebq_LOC_dat[cebq_LOC_dat$loc1 == 'Yes', ]$cebqFF, cebq_LOC_dat[cebq_LOC_dat$loc1 == 'No', ]$cebqFF)

nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqFF) && cebq_LOC_dat$loc1 == 'Yes', ])
nrow(cebq_LOC_dat[is.na(cebq_LOC_dat$cebqFF) && cebq_LOC_dat$loc1 == 'No', ])

##### CFQ ####
#as as for cebq above
cfq_LOC_dat = LOC_dat[c(2, 87:117)]


for(c in 2:32){
  cfq_LOC_dat[[c+31]] = NA
  varname = names(cfq_LOC_dat)[[c]]
  names(cfq_LOC_dat)[[c+31]] = paste(varname, "raw", sep = "_")
  
  if(c < 5){
    cfq_LOC_dat[c+31] = ifelse(is.na(cfq_LOC_dat[[c]]), NA, ifelse(
      cfq_LOC_dat[[c]] == "Never", 1, ifelse(
        cfq_LOC_dat[[c]] == "Seldom", 2, ifelse(
          cfq_LOC_dat[[c]] == "HalfTimes", 3, ifelse(
            cfq_LOC_dat[[c]] == "MostTimes", 4, 5
          )
        )
      )
    )) 
  }
  if(c > 4 & c < 15){
    cfq_LOC_dat[c+31] = ifelse(is.na(cfq_LOC_dat[[c]]), NA, ifelse(
      cfq_LOC_dat[[c]] == "MarkedlyUnderwt", 1, ifelse(
        cfq_LOC_dat[[c]] == "Underwt", 2, ifelse(
          cfq_LOC_dat[[c]] == "Avg", 3, ifelse(
            cfq_LOC_dat[[c]] == "Overwt", 4, 5
          )
        )
      )
    )) 
  }
  if( c > 14 & c < 18){
    cfq_LOC_dat[c+31] = ifelse(is.na(cfq_LOC_dat[[c]]), NA, ifelse(
      cfq_LOC_dat[[c]] == "Unconcerned", 1, ifelse(
        cfq_LOC_dat[[c]] == "SlightlyUncon", 2, ifelse(
          cfq_LOC_dat[[c]] == "Concerned", 3, ifelse(
            cfq_LOC_dat[[c]] == "Slightlycon", 4, 5
          )
        )
      )
    )) 
  }
  if(c > 17 & c < 30){
    cfq_LOC_dat[c+31] = ifelse(is.na(cfq_LOC_dat[[c]]), NA, ifelse(
      cfq_LOC_dat[[c]] == "Disagree", 1, ifelse(
        cfq_LOC_dat[[c]] == "SlightlyDis", 2, ifelse(
          cfq_LOC_dat[[c]] == "Neutral", 3, ifelse(
            cfq_LOC_dat[[c]] == "SlightlyAgree", 4, 5
          )
        )
      )
    )) 
  }
  if(c > 29){
    cfq_LOC_dat[c+31] = ifelse(is.na(cfq_LOC_dat[[c]]), NA, ifelse(
      cfq_LOC_dat[[c]] == "Never", 1, ifelse(
        cfq_LOC_dat[[c]] == "Rarely", 2, ifelse(
          cfq_LOC_dat[[c]] == "Sometimes", 3, ifelse(
            cfq_LOC_dat[[c]] == "Mostly", 4, 5
          )
        )
      )
    )) 
  }
}

#add loc
cfq_LOC_dat = merge(cfq_LOC_dat, LOC_dat[c(2, 33)], by = 'StudyID')

#get means and descriptives for CFQ questionnaire
## PR
cfq_LOC_dat$cfqPR = rowMeans(cfq_LOC_dat[33:35], na.rm = T)
alpha_cfqPR = summary(psych::alpha(cfq_LOC_dat[33:35]))

cfqPR_mean = mean(cfq_LOC_dat$cfqPR, na.rm = TRUE)
cfqPR_sd = sd(cfq_LOC_dat$cfqPR, na.rm = TRUE)
cfqPR_range = range(cfq_LOC_dat$cfqPR, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPR), ])

cebqPR_loc_t = t.test(cfqPR~loc1, data = cfq_LOC_dat)
cebqPR_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPR, cfq_LOC_dat$loc1)
cebqPR_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPR, cfq_LOC_dat$loc1)
cebqPR_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqPR, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqPR)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPR) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPR) && cfq_LOC_dat$loc1 == 'No', ])

## PPW
cfq_LOC_dat$cfqPPW = rowMeans(cfq_LOC_dat[36:39], na.rm = T)
alpha_cfqPPW = summary(psych::alpha(cfq_LOC_dat[36:39]))

cfqPPW_mean = mean(cfq_LOC_dat$cfqPPW, na.rm = TRUE)
cfqPPW_sd = sd(cfq_LOC_dat$cfqPPW, na.rm = TRUE)
cfqPPW_range = range(cfq_LOC_dat$cfqPPW, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPPW), ])

cebqPPW_loc_t = t.test(cfqPPW~loc1, data = cfq_LOC_dat)
cebqPPW_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPPW, cfq_LOC_dat$loc1)
cebqPPW_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPPW, cfq_LOC_dat$loc1)
cebqPPW_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqPPW, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqPPW)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPPW) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPPW) && cfq_LOC_dat$loc1 == 'No', ])

## PCW
cfq_LOC_dat$cfqPCW = rowMeans(cfq_LOC_dat[40:45], na.rm = T)
alpha_cfqPCW = summary(psych::alpha(cfq_LOC_dat[40:45]))

cfqPCW_mean = mean(cfq_LOC_dat$cfqPCW, na.rm = TRUE)
cfqPCW_sd = sd(cfq_LOC_dat$cfqPCW, na.rm = TRUE)
cfqPCW_range = range(cfq_LOC_dat$cfqPCW, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPCW), ])

cebqPCW_loc_t = t.test(cfqPCW~loc1, data = cfq_LOC_dat)
cebqPCW_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPCW, cfq_LOC_dat$loc1)
cebqPCW_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPCW, cfq_LOC_dat$loc1)
cebqPCW_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqPCW, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqPCW)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPCW) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPCW) && cfq_LOC_dat$loc1 == 'No', ])

## CONC
cfq_LOC_dat$cfqCONC = rowMeans(cfq_LOC_dat[46:48], na.rm = T)
alpha_cfqCONC = summary(psych::alpha(cfq_LOC_dat[46:48]))

cfqCONC_mean = mean(cfq_LOC_dat$cfqCONC, na.rm = TRUE)
cfqCONC_sd = sd(cfq_LOC_dat$cfqCONC, na.rm = TRUE)
cfqCONC_range = range(cfq_LOC_dat$cfqCONC, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqCONC), ])

cebqCONC_loc_t = t.test(cfqCONC~loc1, data = cfq_LOC_dat)
cebqCONC_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqCONC, cfq_LOC_dat$loc1)
cebqCONC_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqCONC, cfq_LOC_dat$loc1)
cebqCONC_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqCONC, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqCONC)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqCONC) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqCONC) && cfq_LOC_dat$loc1 == 'No', ])

## REST
cfq_LOC_dat$cfqREST = rowMeans(cfq_LOC_dat[49:56], na.rm = T)
alpha_cfqREST = summary(psych::alpha(cfq_LOC_dat[49:56]))

cfqREST_mean = mean(cfq_LOC_dat$cfqREST, na.rm = TRUE)
cfqREST_sd = sd(cfq_LOC_dat$cfqREST, na.rm = TRUE)
cfqREST_range = range(cfq_LOC_dat$cfqREST, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqREST), ])

cebqREST_loc_t = t.test(cfqREST~loc1, data = cfq_LOC_dat)
cebqREST_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqREST, cfq_LOC_dat$loc1)
cebqREST_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqREST, cfq_LOC_dat$loc1)
cebqREST_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqREST, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqREST)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqREST) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqREST) && cfq_LOC_dat$loc1 == 'No', ])

## PE
cfq_LOC_dat$cfqPE = rowMeans(cfq_LOC_dat[57:60], na.rm = T)
alpha_cfqPE = summary(psych::alpha(cfq_LOC_dat[57:60]))

cfqPE_mean = mean(cfq_LOC_dat$cfqPE, na.rm = TRUE)
cfqPE_sd = sd(cfq_LOC_dat$cfqPE, na.rm = TRUE)
cfqPE_range = range(cfq_LOC_dat$cfqPE, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPE), ])

cebqPE_loc_t = t.test(cfqPE~loc1, data = cfq_LOC_dat)
cebqPE_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPE, cfq_LOC_dat$loc1)
cebqPE_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqPE, cfq_LOC_dat$loc1)
cebqPE_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqPE, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqPE)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPE) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqPE) && cfq_LOC_dat$loc1 == 'No', ])

## MON
cfq_LOC_dat$cfqMON = rowMeans(cfq_LOC_dat[61:63], na.rm = T)
alpha_cfqMON = summary(psych::alpha(cfq_LOC_dat[61:63]))

cfqMON_mean = mean(cfq_LOC_dat$cfqMON, na.rm = TRUE)
cfqMON_sd = sd(cfq_LOC_dat$cfqMON, na.rm = TRUE)
cfqMON_range = range(cfq_LOC_dat$cfqMON, na.rm = TRUE)
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqMON), ])

cebqMON_loc_t = t.test(cfqMON~loc1, data = cfq_LOC_dat)
cebqMON_loc_sd = sd.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqMON, cfq_LOC_dat$loc1)
cebqMON_loc_range = range.function.na(cfq_LOC_dat, cfq_LOC_dat$cfqMON, cfq_LOC_dat$loc1)
cebqMON_loc_d = cohensD(cfq_LOC_dat[cfq_LOC_dat$loc1 == 'Yes', ]$cfqMON, cfq_LOC_dat[cfq_LOC_dat$loc1 == 'No', ]$cfqMON)

nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqMON) && cfq_LOC_dat$loc1 == 'Yes', ])
nrow(cfq_LOC_dat[is.na(cfq_LOC_dat$cfqMON) && cfq_LOC_dat$loc1 == 'No', ])