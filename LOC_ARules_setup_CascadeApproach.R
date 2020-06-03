library(reporttools)
library(xtable)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(MASS)
library(emmeans)
library(reshape)
library(lsr)
library(ordinal)
#library(LDdiag)
library(stats)
library(lubridate)
library(mediation)
library(lavaan)
library(foreign)
library(eeptools)
library(plyr)
library(arules)
library(arulesViz)
source('functions.R')

LOC = read.csv('Data/LOC_redDat.csv')

#########################################
#### Expand Basch Coding--logical ####
#Age ####
LOC$cAge_yr = LOC$cAge_mo/12

LOC$cAge_7.8yr = ifelse(LOC$cAge_yr<9, "T", "F")
LOC$cAge_7.8yr = as.logical(LOC$cAge_7.8yr)

LOC$cAge_7.9yr = ifelse(LOC$cAge_yr<10, "T", "F")
LOC$cAge_7.9yr = as.logical(LOC$cAge_7.9yr)

LOC$cAge_8.9yr = ifelse(LOC$cAge_yr >=8 & LOC$cAge_yr<10, "T", "F")
LOC$cAge_8.9yr = as.logical(LOC$cAge_8.9yr)

LOC$cAge_8.10yr = ifelse(LOC$cAge_yr >=8 & LOC$cAge_yr<11,"T", "F")
LOC$cAge_8.10yr = as.logical(LOC$cAge_8.10yr)

LOC$cAge_9.10yr = ifelse(LOC$cAge_yr>=9 & LOC$cAge_yr<11, "T", "F")
LOC$cAge_9.10yr = as.logical(LOC$cAge_9.10yr)

LOC$cAge_9.12yr = ifelse(LOC$cAge_yr>=9 & LOC$cAge_yr<13,"T", "F")
LOC$cAge_9.12yr = as.logical(LOC$cAge_9.12yr)

LOC$cAge_10.11yr = ifelse(LOC$cAge_yr>=10 & LOC$cAge_yr<12, "T", "F")
LOC$cAge_10.11yr = as.logical(LOC$cAge_10.11yr)

LOC$cAge_10.12yr = ifelse(LOC$cAge_yr>=10 & LOC$cAge_yr<13, "T", "F")
LOC$cAge_10.12yr = as.logical(LOC$cAge_10.12yr)

LOC$cAge_11.12yr = ifelse(LOC$cAge_yr>=11 & LOC$cAge_yr<13, "T", "F")
LOC$cAge_11.12yr = as.logical(LOC$cAge_11.12yr)

LOC$cAge_round = round(LOC$cAge_yr)
LOC$cAge_round = as.logical(LOC$cAge_round)

#Anthropometrics ####
#child BMI
LOC$cBMI = LOC$cBodyMass_class

LOC$cBMI_OWOB = ifelse(is.na(LOC$cBodyMass_class), NA, ifelse(
  LOC$cBodyMass_class == "OW" | LOC$cBodyMass_class == "OB", "T", "F"))
LOC$cBMI_OWOB = as.logical(LOC$cBMI_OWOB)

#birth weight 2 length percentile
LOC$BirthW2L_p_g90 = ifelse(is.na(LOC$birthW2L_p), NA, ifelse(
  LOC$birthW2L_p == "<P95" | LOC$birthW2L_p == "<97" | LOC$birthW2L_p == ">97" , "T", "F"))
LOC$BirthW2L_p_g90 = as.logical(LOC$BirthW2L_p_g90)

LOC$BirthW2L_p_g75 = ifelse(is.na(LOC$birthW2L_p), NA, ifelse(
  LOC$birthW2L_p == "<90" | LOC$birthW2L_p == "<P95" | 
    LOC$birthW2L_p == "<97" | LOC$birthW2L_p == ">97" , "T", "F"))
LOC$BirthW2L_p_g75 = as.logical(LOC$BirthW2L_p_g75)

LOC$BirthW2L_p_ls25 = ifelse(is.na(LOC$birthW2L_p), NA, ifelse(
  LOC$birthW2L_p == "<P3" | LOC$birthW2L_p == "<P5" | 
    LOC$birthW2L_p == "<P10" | LOC$birthW2L_p == "<P25", "T", "F"))
LOC$BirthW2L_p_ls25 = as.logical(LOC$BirthW2L_p_ls25)

LOC$BirthW2L_p_ls10 = ifelse(is.na(LOC$birthW2L_p), NA, ifelse(
  LOC$birthW2L_p == "<P3" | LOC$birthW2L_p == "<P5" | LOC$birthW2L_p == "<P10", "T", "F"))
LOC$BirthW2L_p_ls10 = as.logical(LOC$BirthW2L_p_ls10)

#dad BMI
LOC$dBMI = LOC$dBodyMass_class
LOC$dBMI_OB2.3 = ifelse(is.na(LOC$dBodyMass_class), NA, ifelse(
  LOC$dBodyMass_class == "Obese3" | LOC$dBodyMass_class == "Obese2", "T", "F"))
LOC$dBMI_OB2.3 = as.logical(LOC$dBMI_OB2.3)

LOC$dBMI_OB1.3 = ifelse(is.na(LOC$dBodyMass_class), NA, ifelse(
  LOC$dBodyMass_class == "Obese3" | LOC$dBodyMass_class == "Obese2" |
    LOC$dBodyMass_class == "Obese1", "T", "F"))
LOC$dBMI_OB1.3 = as.logical(LOC$dBMI_OB1.3)

LOC$dBMI_OBOW = ifelse(is.na(LOC$dBodyMass_class), NA, ifelse(
  LOC$dBodyMass_class == "Obese3" | LOC$dBodyMass_class == "Obese2" |
    LOC$dBodyMass_class == "Obese1" | LOC$dBodyMass_class == "Overweight", "T", "F"))
LOC$dBMI_OBOW = as.logical(LOC$dBMI_OBOW)

#mom BMI
LOC$mBMI = LOC$mBodyMass_class
LOC$mBMI_OB2.3 = ifelse(is.na(LOC$mBodyMass_class), NA, ifelse(
  LOC$mBodyMass_class == "Obese3" | LOC$mBodyMass_class == "Obese2", "T", "F"))
LOC$mBMI_OB2.3 = as.logical(LOC$mBMI_OB2.3)

LOC$mBMI_OB1.3 = ifelse(is.na(LOC$mBodyMass_class), NA, ifelse(
  LOC$mBodyMass_class == "Obese3" | LOC$mBodyMass_class == "Obese2" |
    LOC$mBodyMass_class == "Obese1", "T", "F"))
LOC$mBMI_OB1.3 = as.logical(LOC$mBMI_OB1.3)

LOC$mBMI_OBOW = ifelse(is.na(LOC$mBodyMass_class), NA, ifelse(
  LOC$mBodyMass_class == "Obese3" | LOC$mBodyMass_class == "Obese2" |
    LOC$mBodyMass_class == "Obese1" | LOC$mBodyMass_class == "Overweight", "T", "F"))
LOC$mBMI_OBOW = as.logical(LOC$mBMI_OBOW)

#Breast Fed ####
#at least
LOC$BreastFed_least10mo = ifelse(is.na(LOC$BreastFed_mo), NA, ifelse(
  LOC$BreastFed_mo == "<12" | LOC$BreastFed_mo == "10-12", "T", "F"))
LOC$BreastFed_least10mo = as.logical(LOC$BreastFed_least10mo)

LOC$BreastFed_least7mo = ifelse(is.na(LOC$BreastFed_mo), NA, ifelse(
  LOC$BreastFed_mo == "<12" | LOC$BreastFed_mo == "10-12" | 
    LOC$BreastFed_mo == "7-9", "T", "F"))
LOC$BreastFed_least7mo = as.logical(LOC$BreastFed_least7mo)

LOC$BreastFed_least4mo = ifelse(is.na(LOC$BreastFed_mo), NA, ifelse(
  LOC$BreastFed_mo == "<12" | LOC$BreastFed_mo == "10-12" | 
    LOC$BreastFed_mo == "7-9" | LOC$BreastFed_mo == "4-6", "T", "F"))
LOC$BreastFed_least4mo = as.logical(LOC$BreastFed_least4mo)

#less than
LOC$BreastFed_ls10mo = ifelse(is.na(LOC$BreastFed_mo), NA, ifelse(
  LOC$BreastFed_mo == "7-9" | LOC$BreastFed_mo == "4-6"| 
    LOC$BreastFed_mo == "1-3" | LOC$BreastFed_mo == "Never", "T", "F"))
LOC$BreastFed_ls10mo = as.logical(LOC$BreastFed_ls10mo)

LOC$BreastFed_ls7mo = ifelse(is.na(LOC$BreastFed_mo), NA, ifelse(
  LOC$BreastFed_mo == "4-6" | LOC$BreastFed_mo == "1-3"| 
    LOC$BreastFed_mo == "Never", "T", "F"))
LOC$BreastFed_ls7mo = as.logical(LOC$BreastFed_ls7mo)

LOC$BreastFed_ls4mo = ifelse(is.na(LOC$BreastFed_mo), NA, ifelse(
  LOC$BreastFed_mo == "1-3" | LOC$BreastFed_mo == "Never", "T", "F"))
LOC$BreastFed_ls4mo = as.logical(LOC$BreastFed_ls4mo)


#income ####
LOC$income_g75K = ifelse(is.na(LOC$income), NA, ifelse(
  LOC$income == ">$100K" | LOC$income == "$76-100K", "T", "F"))
LOC$income_g75K = as.logical(LOC$income_g75K)

LOC$income_g50K = ifelse(is.na(LOC$income), NA, ifelse(
  LOC$income == ">$100K" | LOC$income == "$76-100K" | 
    LOC$income == "$51-75K" , "T", "F"))
LOC$income_g50K = as.logical(LOC$income_g75K)

LOC$income_ls50K = ifelse(is.na(LOC$income), NA, ifelse(
  LOC$income == "<$20K" | LOC$income == "$21-35K" |
    LOC$income == "$36-50K" , "T", "F"))
LOC$income_ls50K = as.logical(LOC$income_ls50K)

LOC$income_ls35K = ifelse(is.na(LOC$income), NA, ifelse(
  LOC$income == "<$20K" | LOC$income == "$21-35K", "T", "F"))
LOC$income_ls35K = as.logical(LOC$income_ls35K)


#education ####
#mom
LOC$mEducation_gBS = ifelse(is.na(LOC$mEducation), NA, ifelse(
  LOC$mEducation == "MastersDegree" | LOC$mEducation == "DoctoralDegree" |
    LOC$mEducation == "Completed graduate school" | 
    LOC$mEducation == "SomeGraduateSchool", "T", "F"))
LOC$mEducation_gBS = as.logical(LOC$mEducation_gBS)

LOC$mEducation_gAssoc = ifelse(is.na(LOC$mEducation), NA, ifelse(
  LOC$mEducation == "MastersDegree" | LOC$mEducation == "DoctoralDegree" |
    LOC$mEducation == "Completed graduate school" | 
    LOC$mEducation == "SomeGraduateSchool" | 
    LOC$mEducation == "BachelorDegree", "T", "F"))
LOC$mEducation_gAssoc = as.logical(LOC$mEducation_gAssoc)

LOC$mEducation_gHS = ifelse(is.na(LOC$mEducation), NA, ifelse(
  LOC$mEducation == "MastersDegree" | LOC$mEducation == "DoctoralDegree" |
    LOC$mEducation == "Completed graduate school" | 
    LOC$mEducation == "SomeGraduateSchool" | 
    LOC$mEducation == "BachelorDegree" | LOC$mEducation == "AssociatesDegree" |
    LOC$mEducation == "SomeCollege" | 
    LOC$mEducation == "TechnicalDegree", "T", "F"))
LOC$mEducation_gHS = as.logical(LOC$mEducation_gHS)

LOC$mEducation_HS = ifelse(is.na(LOC$mEducation), NA, ifelse(
  LOC$mEducation == "HighSchool", "T", "F"))
LOC$mEducation_HS = as.logical(LOC$mEducation_HS)

LOC$mEducation_lsHS = ifelse(is.na(LOC$mEducation), NA, ifelse(
  LOC$mEducation == "SomeHighSchool", "T", "F"))
LOC$mEducation_lsHS = as.logical(LOC$mEducation_lsHS)

#dad
LOC$dEducation_gBS = ifelse(is.na(LOC$dEducation), NA, ifelse(
  LOC$dEducation == "MastersDegree" | LOC$dEducation == "DoctoralDegree" |
    LOC$dEducation == "Completed graduate school" | 
    LOC$dEducation == "SomeGraduateSchool", "T", "F"))
LOC$dEducation_gBS = as.logical(LOC$dEducation_gBS)

LOC$dEducation_gAssoc = ifelse(is.na(LOC$dEducation), NA, ifelse(
  LOC$dEducation == "MastersDegree" | LOC$dEducation == "DoctoralDegree" |
    LOC$dEducation == "Completed graduate school" | 
    LOC$dEducation == "SomeGraduateSchool" | 
    LOC$dEducation == "BachelorDegree", "T", "F"))
LOC$dEducation_gAssoc = as.logical(LOC$dEducation_gAssoc)

LOC$dEducation_gHS = ifelse(is.na(LOC$dEducation), NA, ifelse(
  LOC$dEducation == "MastersDegree" | LOC$dEducation == "DoctoralDegree" |
    LOC$dEducation == "Completed graduate school" | 
    LOC$dEducation == "SomeGraduateSchool" | 
    LOC$dEducation == "BachelorDegree" | LOC$dEducation == "AssociatesDegree" |
    LOC$dEducation == "SomeCollege" | 
    LOC$dEducation == "TechnicalDegree", "T", "F"))
LOC$dEducation_gHS = as.logical(LOC$dEducation_gHS)

LOC$dEducation_HS = ifelse(is.na(LOC$dEducation), NA, ifelse(
  LOC$dEducation == "HighSchool", "T", "F"))
LOC$dEducation_HS = as.logical(LOC$dEducation_HS)

LOC$dEducation_lsHS = ifelse(is.na(LOC$dEducation), NA, ifelse(
  LOC$dEducation == "SomeHighSchool", "T", "F"))
LOC$dEducation_lsHS = as.logical(LOC$dEducation_lsHS)

#family eating ####
#eating out

LOC$EatOut_g1wk = ifelse(is.na(LOC$EatOut), NA, ifelse(
  LOC$EatOut == "2/wk" | LOC$EatOut == "3/wk" | 
    LOC$EatOut == "4/wk" | LOC$EatOut == "4+/wk", "T", "F"))
LOC$EatOut_g1wk = as.logical(LOC$EatOut_g1wk)

LOC$EatOut_g2wk = ifelse(is.na(LOC$EatOut), NA, ifelse(
  LOC$EatOut == "3/wk" | LOC$EatOut == "4/wk" | LOC$EatOut == "4+/wk", "T", "F"))
LOC$EatOut_g2wk = as.logical(LOC$EatOut_g2wk)

LOC$EatOut_1.4mo = ifelse(is.na(LOC$EatOut), NA, ifelse(
  LOC$EatOut == "1/mo" | LOC$EatOut == "2/mo" | LOC$EatOut == "1/wk", "T", "F"))
LOC$EatOut_1.4mo = as.logical(LOC$EatOut_1.4mo)

#FamilyEat together
LOC$EatFamily_g1wk = ifelse(is.na(LOC$EatFamily), NA, ifelse(
  LOC$EatFamily == "2/wk" | LOC$EatFamily == "3/wk" | 
    LOC$EatFamily == ">4/wk", "T", "F"))
LOC$EatFamily_g1wk = as.logical(LOC$EatFamily_g1wk)

LOC$EatFamily_g2wk = ifelse(is.na(LOC$EatFamily), NA, ifelse(
  LOC$EatFamily == "3/wk" | LOC$EatFamily == ">4/wk", "T", "F"))
LOC$EatFamily_g2wk = as.logical(LOC$EatFamily_g2wk)

LOC$EatFamily_1.4mo = ifelse(is.na(LOC$EatFamily), NA, ifelse(
  LOC$EatFamily == "1/mo" | LOC$EatFamily == "1/wk", "T", "F"))
LOC$EatFamily_1.4mo = as.logical(LOC$EatFamily_1.4mo)

#CEBQ ####
for(c in 52:86){
  if(c == 54|c == 55|c == 61|c == 67|c == 83){
    #often and always
    name_gSometimesRev = paste(names(LOC)[[c]], "gSometimesRev", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gSometimesRev[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Never" | LOC[[c]] == "Rarely", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    # #sometimes, often, always
    name_gRarelyRev = paste(names(LOC)[[c]], "gRarelyRev", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gRarelyRev[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Never" | LOC[[c]] == "Rarely" |
        LOC[[c]] == "Sometimes", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #Rarely, Never
    name_lsSometimesRev = paste(names(LOC)[[c]], "lsSometimesRev", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_lsSometimesRev[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Often" | LOC[[c]] == "Always" |
        LOC[[c]] == "Sometimes", "F", "T"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
  else {
    #often and always
    name_gSometimes = paste(names(LOC)[[c]], "gSometimes", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gSometimes[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "Often", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #sometimes, often, always
    name_gRarely = paste(names(LOC)[[c]], "gRarely", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gRarely[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "Often" |
        LOC[[c]] == "Sometimes", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #Never, Rarely
    name_lsSometimes = paste(names(LOC)[[c]], "lsSometimes", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_lsSometimes[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "Often" |
        LOC[[c]] == "Sometimes", "F", "T"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
}

#CFQ ####
for(c in 87:117){
  if(c < 90){
    #most and always
    name_gHalf = paste(names(LOC)[[c]], "gHalf", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gHalf[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "MostTimes", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #half, most, always
    name_gSeldom = paste(names(LOC)[[c]], "gSeldom", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gSeldom[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "MostTimes" |
        LOC[[c]] == "HalfTimes", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
  if(c > 89 & c < 100){
    #overweight and Markedly overweight
    name_gOW = paste(names(LOC)[[c]], "gOW", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gOW[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "MarkedlyOverwt" | LOC[[c]] == "Overwt", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #underweight and markedly underweight
    name_lsUW = paste(names(LOC)[[c]], "lsUW", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_lsUW[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "MarkedlyUnderwt" | LOC[[c]] == "Underwt", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
  if(c > 99 & c < 103){
    #slightly and very concerned
    name_gNeutralCon = paste(names(LOC)[[c]], "gNeutralCon", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gNeutralCon[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "SlightlyCon" | LOC[[c]] == "VeryCon", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #slightly unconcerned and  unconcerned
    name_lsNeutralCon = paste(names(LOC)[[c]], "lsNeutralCon", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_lsNeutralCon[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "SlightlyUncon" | LOC[[c]] == "Unconcerned", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
  if( c > 102 & c < 115){
    #slightly agree and agree
    name_gNeutralAD = paste(names(LOC)[[c]], "gNeutralAD", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gNeutralAD[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "SlightlyAgree" | LOC[[c]] == "Agree", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #slightly disagree and  disagree
    name_lsNeutralAD = paste(names(LOC)[[c]], "lsNeutralAD", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_lsNeutralAD[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "SlightlyDis" | LOC[[c]] == "Disagree", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
  if(c > 114){
    #most and always
    name_gSometimes = paste(names(LOC)[[c]], "gSometimes", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gSometimes[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "Mostly", "T", "F")) 
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
    
    #Sometimes, most, always
    name_gRarely = paste(names(LOC)[[c]], "gRarely", sep = "_")
    LOC[[(length(variable.names(LOC))+1)]] = NA
    names(LOC)[length(variable.names(LOC))] = name_gRarely[1]
    
    LOC[[length(variable.names(LOC))]] = ifelse(is.na(LOC[[c]]), NA, ifelse(
      LOC[[c]] == "Always" | LOC[[c]] == "Mostly" |
        LOC[[c]] == "Sometimes", "T", "F"))
    LOC[[length(variable.names(LOC))]] = as.logical(LOC[[length(variable.names(LOC))]])
  }
}

#########################################
##remove duplicates
LOC_1V = LOC[LOC$DuplicateEx == "N", ]
write.csv(LOC, file = 'Data/LOC_all.csv', row.names = FALSE)
write.csv(LOC_1V, file = 'Data/LOC_noDuplicates.csv', row.names = FALSE)

