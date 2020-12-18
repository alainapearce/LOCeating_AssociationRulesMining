
library(arules)
library(ggplot2)
library(factoextra)
library(cluster)
library(cowplot)
library(reshape2)
library(clValid)




source('functions.R')
source('ARM_ORconf.R')
LOC_dat = read.csv('Data/LOC_binarized_data.csv')

###############################################
####            ASSOCIATION RULES          ####
###############################################
#mother completed
ParentComplete_noNA = xtabs(~relationship, data = LOC_dat[!is.na(LOC_dat$loc1), ])
ParentComplete = xtabs(~relationship, data = LOC_dat)

#### convert to transactions ####
#all questions expanded below
LOC_arules = LOC_dat[!is.na(LOC_dat$loc1), c(4, 119:127, 129:134, 135:168, 33, 52:117, 169:335)]

LOC_arules_trans = as(LOC_arules, "transactions")

#get all transactions (going to be the LHS items)--only have to do this once then comment out
# transDataFrame = as.matrix(LOC_arules_trans@itemInfo$labels)

#export to manually add the category labes you want
#write.csv(transDataFrame, file = "ResultsOutput/indQ_TransactionsLabels.csv", row.names = FALSE)

#read in transaction labels
transDataLabels = read.csv("ResultsOutput/indQ_TransactionsLabels.csv", header = TRUE)

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

pruned.rules_lhs1_conf33 = subset(pruned.rules_lhs1_conf33, subset = quality(pruned.rules_lhs1_conf33)[8]>=0.21)

#add addjustment for multiple comparisons to data
quality(pruned.rules_lhs1_conf33) = cbind(quality(pruned.rules_lhs1_conf33), 
  fisher.padj_holm = round(p.adjust(quality(pruned.rules_lhs1_conf33)$fisher.p, method = 'holm'), 4))

#sort rules
pruned.rules_lhs1_conf33 = sort(pruned.rules_lhs1_conf33, by=c("fisher.p","fisher.padj_holm"), decreasing = FALSE)

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

pruned.rules_lhs2_conf50 = subset(pruned.rules_lhs2_conf50, subset = quality(pruned.rules_lhs2_conf50)[8]>=0.20)

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

#test internal validity metrics - these take a long time to run
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

pruned.rules_lhs1_conf33_no_red = subset(pruned.rules_lhs1_conf33_no, quality(pruned.rules_lhs1_conf33_no)[8]>=0.21)

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

pruned.rules_lhs2_conf50_no_red = subset(pruned.rules_lhs2_conf50_no_all, quality(pruned.rules_lhs2_conf50_no_all)[8]>=0.21)

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

#test internal validity metrics - can take a while to run
# pruned.rules_lhs2_conf50_no_nclust_internV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50_no), 2:10, clMethods=c("pam"), validation="internal")

#test stability metrics (remove one collumn at a time sequentially and look at stability)
# pruned.rules_lhs2_conf50_no_nclust_stabV <- clValid(as.matrix(dist_pruned.rules_lhs2_conf50_no), 2:10, clMethods=c("pam"), validation="stability")

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
