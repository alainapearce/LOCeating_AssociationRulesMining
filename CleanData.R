############ Basic Data Load/Setup ############
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
source('/Users/azp271/Library/Mobile Documents/com~apple~CloudDocs/Graduate/Research/functions.R')

################################################
LOC_dat = read.csv('Data/All_data.csv', na.strings = c('<NA>', '#NULL!', '', ' '))
LOC_dat = LOC_dat[!grepl("cceb", LOC_dat$StudyID) & 
    !grepl("Dairy", LOC_dat$StudyID) & 
    !grepl("CTSI", LOC_dat$StudyID) &
    !grepl("VegBrand", LOC_dat$StudyID) &
    !grepl("Spice", LOC_dat$StudyID), ]

################################################

#### Combined Loss of Control Databases into 1 ####
LOC_dat$DuplicateEx = ifelse(LOC_dat$StudyID == "FBS_1371" | LOC_dat$StudyID == "DMK_111" |
                               LOC_dat$StudyID == "EBS_119" | LOC_dat$StudyID == "DMK_122" |
                               LOC_dat$StudyID == "FBS_1021" | LOC_dat$StudyID == "DMK_102" |
                               LOC_dat$StudyID == "TestRetest_107" | LOC_dat$StudyID == "EBS_108" |
                               LOC_dat$StudyID == "FBS_1281" | LOC_dat$StudyID == "DMK_115" | 
                               LOC_dat$StudyID == "FBS_1091" | LOC_dat$StudyID == "DMK_116" |
                               LOC_dat$StudyID == "FBS_1481" | LOC_dat$StudyID == "DMK_125" |
                               LOC_dat$StudyID == "EBS_115" | LOC_dat$StudyID == "DMK_157" |
                               LOC_dat$StudyID == "EBS_118" | LOC_dat$StudyID == "DMK_135" |
                               LOC_dat$StudyID == "EBS_112" | LOC_dat$StudyID == "EBS_111" |
                               LOC_dat$StudyID == "DMK_119" | LOC_dat$StudyID == "EBS_101" |
                               LOC_dat$StudyID == "DMK_120" | LOC_dat$StudyID == "DMK_101" |
                               LOC_dat$StudyID == "DMK_121" | LOC_dat$StudyID == "DMK_108" |
                               LOC_dat$StudyID == "DMK_107" | LOC_dat$StudyID == "DMK_129" |
                               LOC_dat$StudyID == "TestRetest_101" | LOC_dat$StudyID == "TestRetest_111" |
                               LOC_dat$StudyID == "TestRetest_121" | LOC_dat$StudyID == "TestRetest_122" |
                               LOC_dat$StudyID == "TestRetest_125", "Y", "N")

LOC_dat = LOC_dat[c(1:2, 6:226)]
################################################

#Birth Weight 2 Length ####
wtleninf_fp = 'Data/wtleninf.csv'

wtleninf = read.csv(wtleninf_fp)
wtleninf_boy = wtleninf[wtleninf$Sex == 1, ]
wtleninf_girl = wtleninf[wtleninf$Sex == 2, ]

LOC_dat$cBirthLength_cm = as.numeric(as.character(LOC_dat$cBirthLength_cm))
LOC_dat$cBirthWeight_kg = as.numeric(as.character(LOC_dat$cBirthWeight_kg))

LOC_dat$w2l_row = NA

for (p in 1:length(LOC_dat$StudyID)){
  if (!is.na(LOC_dat[p, 'cBirthLength_cm'])){
    if (as.character(LOC_dat[p, 'sex']) == "Male"){
      round_cm = round_any(LOC_dat[p, 'cBirthLength_cm'], 0.5, floor)
      op = round_cm %% 1
      if (op == 0){
        round_cm = round_cm + 0.5
      }
      row = which(wtleninf_boy$Length == round_cm)
      LOC_dat[p, 'w2l_row'] = row[1]
    }
    else {
      if (as.character(LOC_dat[p, 'sex']) == "Female"){
        round_cm = round_any(LOC_dat[p, 'cBirthLength_cm'], 0.5, floor)
        op = round_cm %% 1
        if (op == 0){
          round_cm = round_cm + 0.5
        }
        row = which(wtleninf_girl$Length == round_cm)
        LOC_dat[p, 'w2l_row'] = row[1]
      }
    }
  }  
}

LOC_dat$birthW2L_p = ifelse(LOC_dat$sex == "Boy", ifelse(
  !is.na(LOC_dat$w2l_row), ifelse(
    LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 6], "<P3", ifelse(
      LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 7], "<P5", ifelse(
        LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 8], "<P10", ifelse(
          LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 9], "<P25", ifelse(
            LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 10], "<P50", ifelse(
              LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 11], "<P75", ifelse(
                LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 12], "<P90", ifelse(
                  LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 13], "<P95", ifelse(
                    LOC_dat$cBirthWeight_kg < wtleninf_boy[LOC_dat$w2l_row, 14], "<P97", ">97"
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
   NA),                         
ifelse(!is.na(LOC_dat$w2l_row), ifelse(
  LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 6], "<P3", ifelse(
    LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 7], "<P5", ifelse(
      LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 8], "<P10", ifelse(
        LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 9], "<P25", ifelse(
          LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 10], "<P50", ifelse(
            LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 11], "<P75", ifelse(
              LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 12], "<P90", ifelse(
                LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 13], "<P95", ifelse(
                  LOC_dat$cBirthWeight_kg < wtleninf_girl[LOC_dat$w2l_row, 14], "<P97", ">97"
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
  NA)
)

LOC_dat$birthW2L_p = factor(LOC_dat$birthW2L_p)

#BreastFeeding ####
LOC_dat$BreastFed = factor(LOC_dat$BreastFed)

LOC_dat$BreastFed_numeric = as.numeric(as.character(LOC_dat$BreastFed_wk))
LOC_dat$BreastFed_numeric = LOC_dat$BreastFed_numeric/4

LOC_dat$BreastFed_mo = ifelse(LOC_dat$BreastFed == "No", "Never", ifelse(
  LOC_dat$BreastFed_wk == "4-12 wks", "1-3", ifelse(
    LOC_dat$BreastFed_wk == "16-24 wks", "4-6", ifelse(
      LOC_dat$BreastFed_wk == ">24 wks", "7-9", ifelse(
        LOC_dat$BreastFed_numeric == 0, "Never", ifelse(
          LOC_dat$BreastFed_numeric < 1, "<1", ifelse(
            LOC_dat$BreastFed_numeric < 4, "1-3", ifelse(
              LOC_dat$BreastFed_numeric < 7, "4-6", ifelse(
                LOC_dat$BreastFed_numeric < 10, "7-9", ifelse(
                  LOC_dat$BreastFed_numeric <= 12, "10-12", ifelse(
                    LOC_dat$BreastFed_numeric > 12, ">12", NA)
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

LOC_dat$BreastFed_mo = factor(LOC_dat$BreastFed_mo)

#EatOut EatFam
LOC_dat$EatOut = ifelse(LOC_dat$EatOut == "Four or more times a week", "4+/wk", as.character(LOC_dat$EatOut))

LOC_dat$EatOut = factor(LOC_dat$EatOut)

LOC_dat$EatFamily = ifelse(LOC_dat$EatFamily == "1", "1/wk", ifelse(
  LOC_dat$EatFamily == "2", "2/wk", ifelse(
    LOC_dat$EatFamily == "3", "3/wk", ifelse(
      LOC_dat$EatFamily == "4" | LOC_dat$EatFamily == "5" | LOC_dat$EatFamily == "6" | LOC_dat$EatFamily == "7" | LOC_dat$EatFamily == "8", ">4/wk", ifelse(
        LOC_dat$EatFamily == " ", NA, as.character(LOC_dat$EatFamily))
      )
    )
  )
)

LOC_dat$EatFamily = factor(LOC_dat$EatFamily)

#CFQ ####
for(c in 136:166){
  LOC_dat[[c]] = ifelse(LOC_dat[[c]] == "" | LOC_dat[[c]] == " " | LOC_dat[[c]] == "#NULL!", NA, ifelse(
    LOC_dat[[c]] == "Half of the time" | LOC_dat[[c]] == "HalfoftheTime", "HalfTimes", ifelse(
      LOC_dat[[c]] == "Most of the time", "MostTimes", as.character(LOC_dat[[c]])
      )
    )
  )
  
  LOC_dat[[c]] = ifelse(LOC_dat[[c]] == "Average", "Avg", ifelse(
    LOC_dat[[c]] == "Markedly Underweight", "MarkedlyUnderwt", ifelse(
      LOC_dat[[c]] == "Underweight", "Underwt", ifelse(
        LOC_dat[[c]] == "Overweight", "Overwt", ifelse(
          LOC_dat[[c]] == "Markedly overweight", "MarkedlyOverwt", ifelse(
            LOC_dat[[c]] == "Not Applicable", NA, ifelse(
              LOC_dat[[c]] == "N/A", NA, as.character(LOC_dat[[c]])
              )
            )
          )
        )
      )  
    )
  )
  
  
  
  LOC_dat[[c]] = ifelse(LOC_dat[[c]] == "Slightly unconcerned", "SlightlyUncon", ifelse(
    LOC_dat[[c]] == "Slightly concerned", "SlightlyCon", ifelse(
      LOC_dat[[c]] == "Very concerned", "VeryCon", as.character(LOC_dat[[c]])
      )  
    )
  )
  
  LOC_dat[[c]] = ifelse(LOC_dat[[c]] == "Slightly agree", "SlightlyAgree", ifelse(
    LOC_dat[[c]] == "Slightly disagree", "SlightlyDis", as.character(LOC_dat[[c]])
    )  
  )

  LOC_dat[[c]] = factor(LOC_dat[[c]])
}
################################################
#### Make Final Variable Selection and Export ####
LOC_redDat = LOC_dat[c(1:20, 226, 22:25, 27, 35:40, 205:223, 56:90, 136:166)]
write.csv(LOC_redDat, 'Data/LOC_redDat.csv', row.names = FALSE)
