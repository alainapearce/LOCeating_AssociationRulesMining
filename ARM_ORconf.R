#function to calculate odds ratio confidence intervals
#for association rules

ARM_ORconf = function(rules, dataset, nlhs, outcome){
  rules_data = DATAFRAME(rules)
  
  library(dplyr)
  library(lsr)
  if(nlhs==1){
    ORdat = data.frame(matrix(, nrow = nrow(rules_data), ncol = 5))
    colnames(ORdat) = c("LHS", "OR", "OR_lowerCI", "OR_upperCI", "CramersV")
    
    for(r in 1:nrow(rules_data)){
      lhs_item = substr(as.character(rules_data[r, 1]), 2, (nchar(as.character(rules_data[r, 1]))-1))
      
      ORdat[r, 1] = lhs_item
      
      True_str = TRUE
      
      True_str = ifelse(grepl("=", lhs_item), strsplit(lhs_item, "=")[[1]][2], TRUE)
      lhs_item = ifelse(grepl("=", lhs_item), strsplit(lhs_item, "=")[[1]][1], as.character(lhs_item))

      dat_new = dplyr::select(dataset, one_of(c(lhs_item, "loc1")))
      
      if(outcome == "Yes"){
        a = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "Yes" & !is.na(dat_new[2]), ])
        b = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "No", ])
        c = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "Yes" & !is.na(dat_new[2]), ])
        d = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "No", ])
        
        OR = (a*d)/(b*c)
        SE_logOR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
        logOR_CIlower = exp(log(OR)-(1.96*SE_logOR))
        logOR_CIupper = exp(log(OR)+(1.96*SE_logOR))
      } else if (outcome == "No"){
        a = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "Yes", ])
        b = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "No" & !is.na(dat_new[2]), ])
        c = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "Yes" , ])
        d = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "No" & !is.na(dat_new[2]), ])
        
        OR = (b*c)/(a*d)
        SE_logOR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
        logOR_CIlower = exp(log(OR)-(1.96*SE_logOR))
        logOR_CIupper = exp(log(OR)+(1.96*SE_logOR))
      }

      ORdat[r, 2] = OR
      ORdat[r, 3] = logOR_CIlower
      ORdat[r, 4] = logOR_CIupper
      
      #cramers d
      cramers_tab = data.frame(matrix(c(a, b, c, d), ncol = 2, byrow = TRUE))
      V = cramersV(cramers_tab)
      ORdat[r, 5]=V
    } 
  } else {
    ORdat = data.frame(matrix(, nrow = nrow(rules_data), ncol = 6))
    colnames(ORdat) = c("LHS1", "LHS2", "OR", "OR_lowerCI", "OR_upperCI", "CramersV")
    for(r in 1:nrow(rules_data)){
      lhs_split = as.data.frame.character(strsplit(as.character(rules_data[r, 1]), ',')[[1]])
      
      if (nrow(lhs_split) == 1){
        lhs_item = substr(as.character(rules_data[r, 1]), 2, (nchar(as.character(rules_data[r, 1]))-1))
        
        ORdat[r, 1] = lhs_item
        ORdat[r, 2] = NA
        
        True_str = TRUE
        
        True_str = ifelse(grepl("=", lhs_item), strsplit(lhs_item, "=")[[1]][2], TRUE)
        lhs_item = ifelse(grepl("=", lhs_item), strsplit(lhs_item, "=")[[1]][1], as.character(lhs_item))
        
        dat_new = dplyr::select(dataset, one_of(c(lhs_item, "loc1")))
        
        if(outcome == "Yes"){
          a = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "Yes" & !is.na(dat_new[2]), ])
          b = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "No", ])
          c = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "Yes" & !is.na(dat_new[2]), ])
          d = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "No", ])
          
          OR = (a*d)/(b*c)
          SE_logOR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
          logOR_CIlower = exp(log(OR)-(1.96*SE_logOR))
          logOR_CIupper = exp(log(OR)+(1.96*SE_logOR))
        } else if (outcome == "No"){
          a = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "Yes", ])
          b = nrow(dat_new[dat_new[1] == True_str & !is.na(dat_new[1]) & dat_new[2] == "No" & !is.na(dat_new[2]), ])
          c = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "Yes" , ])
          d = nrow(dat_new[dat_new[1] != True_str & dat_new[2] == "No" & !is.na(dat_new[2]), ])
          
          OR = (b*c)/(a*d)
          SE_logOR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
          logOR_CIlower = exp(log(OR)-(1.96*SE_logOR))
          logOR_CIupper = exp(log(OR)+(1.96*SE_logOR))
        }

        ORdat[r, 3] = OR
        ORdat[r, 4] = logOR_CIlower
        ORdat[r, 5] = logOR_CIupper
        
        #cramers d
        cramers_tab = data.frame(matrix(c(a, b, c, d), ncol = 2, byrow = TRUE))
        V = cramersV(cramers_tab)
        ORdat[r, 6]=V
        
      } else if (nrow(lhs_split) == 2){
        lhs_item1 = substr(as.character(lhs_split[1, 1]), 2, nchar(as.character(lhs_split[1, 1])))
        lhs_item2 = substr(as.character(lhs_split[2, 1]), 1, (nchar(as.character(lhs_split[2, 1]))-1))
      
        ORdat[r, 1] = lhs_item1
        ORdat[r, 2] = lhs_item2
      
        True_str1 = TRUE
        True_str2 = TRUE

        True_str1 = ifelse(grepl("=", lhs_item1), strsplit(lhs_item1, "=")[[1]][2], TRUE)
        True_str2 = ifelse(grepl("=", lhs_item2), strsplit(lhs_item2, "=")[[1]][2], TRUE)
        
        lhs_item1 = ifelse(grepl("=", lhs_item1), strsplit(lhs_item1, "=")[[1]][1], as.character(lhs_item1))
        lhs_item2 = ifelse(grepl("=", lhs_item2), strsplit(lhs_item2, "=")[[1]][1], as.character(lhs_item2))
      
        dat_new = dplyr::select(dataset, one_of(c(lhs_item1, lhs_item2, "loc1")))
        
        if(outcome == "Yes"){
          a = nrow(dat_new[dat_new[1] == True_str1 & dat_new[2] == True_str2 &
                             !is.na(dat_new[1]) & !is.na(dat_new[2]) & 
                             !is.na(dat_new[3]) & dat_new[3] == "Yes", ])
          b = nrow(dat_new[dat_new[1] == True_str1 & dat_new[2] == True_str2 & 
                             !is.na(dat_new[1]) & !is.na(dat_new[2]) & 
                             dat_new[3] == "No", ])
          
          c_1 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] == True_str2 & !is.na(dat_new[2]) &
                               dat_new[3] == "Yes" & !is.na(dat_new[3]), ])
          c_2 = nrow(dat_new[dat_new[1] == True_str1 & !is.na(dat_new[1]) & dat_new[2] != True_str2 & 
                               dat_new[3] == "Yes" & !is.na(dat_new[3]), ])
          c_3 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] != True_str2 & 
                               dat_new[3] == "Yes" & !is.na(dat_new[3]), ])
          c = c_1 + c_2 + c_3
          
          d_1 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] == True_str2 & !is.na(dat_new[2]) &
                               dat_new[3] == "No", ])
          d_2 = nrow(dat_new[dat_new[1] == True_str1 & !is.na(dat_new[1]) & dat_new[2] != True_str2 & 
                               dat_new[3] == "No", ])
          d_3 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] != True_str2 & 
                               dat_new[3] == "No", ])
          d = d_1 + d_2 + d_3
          
          OR = (a*d)/(b*c)
          SE_logOR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
          logOR_CIlower = exp(log(OR)-(1.96*SE_logOR))
          logOR_CIupper = exp(log(OR)+(1.96*SE_logOR))
        } else if (outcome == "No"){
          a = nrow(dat_new[dat_new[1] == True_str1 & dat_new[2] == True_str2 &
                             !is.na(dat_new[1]) & !is.na(dat_new[2]) & dat_new[3] == "Yes", ])
          b = nrow(dat_new[dat_new[1] == True_str1 & dat_new[2] == True_str2 & 
                             !is.na(dat_new[1]) & !is.na(dat_new[2]) & 
                             dat_new[3] == "No" & !is.na(dat_new[3]), ])
          
          c_1 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] == True_str2 & !is.na(dat_new[2]) &
                               dat_new[3] == "Yes", ])
          c_2 = nrow(dat_new[dat_new[1] == True_str1 & !is.na(dat_new[1]) & dat_new[2] != True_str2 & 
                               dat_new[3] == "Yes", ])
          c_3 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] != True_str2 & 
                               dat_new[3] == "Yes", ])
          c = c_1 + c_2 + c_3
          d_1 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] == True_str2 & !is.na(dat_new[2]) &
                               dat_new[3] == "No" & !is.na(dat_new[3]), ])
          d_2 = nrow(dat_new[dat_new[1] == True_str1 & !is.na(dat_new[1]) & dat_new[2] != True_str2 & 
                               dat_new[3] == "No" & !is.na(dat_new[3]), ])
          d_3 = nrow(dat_new[dat_new[1] != True_str1 & dat_new[2] != True_str2 & 
                               dat_new[3] == "No" & !is.na(dat_new[3]), ])
          d = d_1 + d_2 + d_3
          
          OR = (b*c)/(a*d)
          SE_logOR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
          logOR_CIlower = exp(log(OR)-(1.96*SE_logOR))
          logOR_CIupper = exp(log(OR)+(1.96*SE_logOR))
        }
        
        
        ORdat[r, 3] = OR
        ORdat[r, 4] = logOR_CIlower
        ORdat[r, 5] = logOR_CIupper
        
        #cramers d
        cramers_tab = data.frame(matrix(c(a, b, c, d), ncol = 2, byrow = TRUE))
        V = cramersV(cramers_tab)
        ORdat[r, 6]=V
      }    
    }
  }
  return(ORdat)
}
