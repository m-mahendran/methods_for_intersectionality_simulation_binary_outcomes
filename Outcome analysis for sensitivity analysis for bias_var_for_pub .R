# SENSITIVITY RESULTS ANALYSIS: for the bias/variance sensitivity analysis,run these 
# analyses to assess data. Run with each of the three following sets of data separately: 


#rare prevalence
results2000<- read.csv("C:/Bin_rare_cat_bias_var.2000.csv")
results5000<- read.csv("C:/Bin_rare_cat_bias_var.5000.csv")

#common prevalence
results2000<- read.csv("C:/Bin_common_cat_bias_var.2000.csv")
results5000<- read.csv("C:/Bin_common_cat_bias_var.5000.csv")

#50% prevalence
results2000<- read.csv("C:/Bin_50_preval_cat_bias_var.2000.csv")
results5000<- read.csv("C:/Bin_50_preval_cat_bias_var.5000.csv")

#####################################################################

#run for each of the three sets: 

library(tidyverse)

create_long_bias_var <- function(n, results) {
  
  a <- cbind.data.frame(rep(n, 192), rep("CC", 192), 
                        apply((results[c(which(colnames(results)=="X1.CC.bias")
                                         :which(colnames(results)=="X192.CC.bias"))]*100), 
                              2, mean, na.rm = TRUE), 
                        apply((results[c(which(colnames(results)=="X1.CC.pred")
                                         :which(colnames(results)=="X192.CC.pred"))]*100), 
                              2, var, na.rm = TRUE))
  colnames(a) <- c ("n", "Method", "BIAS", "VAR")
  
  b <- cbind.data.frame(rep(n, 192), rep("MLM", 192), 
                        apply((results[c(which(colnames(results)=="X1.MLM.bias")
                                         :which(colnames(results)=="X192.MLM.bias"))]*100), 
                              2, mean, na.rm = TRUE), 
                        apply((results[c(which(colnames(results)=="X1.MLM.pred")
                                         :which(colnames(results)=="X192.MLM.pred"))]*100), 
                              2, var, na.rm = TRUE))
  colnames(b) <- c ("n", "Method", "BIAS", "VAR")
  
  c <- cbind.data.frame(rep(n, 192), rep("MF", 192), 
                        apply((results[c(which(colnames(results)=="X1.MF.bias")
                                         :which(colnames(results)=="X192.MF.bias"))]*100), 
                              2, mean, na.rm = TRUE), 
                        apply((results[c(which(colnames(results)=="X1.MF.pred")
                                         :which(colnames(results)=="X192.MF.pred"))]*100), 
                              2, var, na.rm = TRUE))
  colnames(c) <- c ("n", "Method", "BIAS", "VAR")
  
  d <- cbind.data.frame(rep(n, 192), rep("BF", 192), 
                        apply((results[c(which(colnames(results)=="X1.BF.bias")
                                         :which(colnames(results)=="X192.BF.bias"))]*100), 
                              2, mean, na.rm = TRUE), 
                        apply((results[c(which(colnames(results)=="X1.BF.pred")
                                         :which(colnames(results)=="X192.BF.pred"))]*100), 
                              2, var, na.rm = TRUE))
  colnames(d) <- c ("n", "Method", "BIAS", "VAR")
  
  final <- rbind.data.frame(a, b, c, d)
  
}


bias_var_2000   <- create_long_bias_var(2000, results2000)
bias_var_5000   <- create_long_bias_var(5000, results5000)

final_bias_var <- rbind.data.frame(bias_var_2000, bias_var_5000)
final_bias_var <- na.omit(final_bias_var)



final_bias_var %>%
  group_by(Method, n) %>%
  summarise(
    median_bias = median(BIAS),
    min_bias = min(BIAS),
    max_bias = max(BIAS),
    median_var = median(VAR),
    min_var = min(VAR),
    max_var = max(VAR)
  )


