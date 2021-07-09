# RESULTS ANALYSIS: for binary outcome with mixed inputs models,run these 
# analyses to assess data. Run with each of the two following sets of data separately: 

#Common binary outcome, mixed inputs
  results2000<- read.csv ("T:/Bin_common_mixed.2000.csv")
  results5000<- read.csv ("T:/Bin_common_mixed.5000.csv")
  results50000<- read.csv ("T:/Bin_common_mixed.50000.csv")
  results200000<- read.csv ("T:/Bin_common_mixed.2e+05.csv")
  
  a_list <- c( "T:/Bin_common_mixed.2000.csv",
               "T:/Bin_common_mixed.5000.csv",
               "T:/Bin_common_mixed.50000.csv",
               "T:/Bin_common_mixed.2e+05.csv") 

#Rare binary outcome, mixed inputs
  results2000<- read.csv ("T:/Bin_rare_mixed.2000.csv")
  results5000<- read.csv ("T:/Bin_rare_mixed.5000.csv")
  results50000<- read.csv ("T:/Bin_rare_mixed.50000.csv")
  results200000<- read.csv ("T:/Bin_rare_mixed.2e+05.csv")
  
  a_list <- c( "T:/Bin_rare_mixed.2000.csv",
               "T:/Bin_rare_mixed.5000.csv",
               "T:/Bin_rare_mixed.50000.csv",
               "T:/Bin_rare_mixed.2e+05.csv")

#####################################################################

#run for each of the two sets: 

###########primary outcome box plots

  
  create_long_MAD <- function(n, results) {
    
    a <- cbind.data.frame(rep(n, 1000), rep("CART", 1000), results$MAD_mean.CART)
    colnames(a) <- c ("n", "Method", "MAD") 
    
    b <- cbind.data.frame(rep(n, 1000), rep("CTree", 1000), results$MAD_mean.CTree)
    colnames(b) <- c ("n", "Method", "MAD") 
    
    c <- cbind.data.frame(rep(n, 1000), rep("RF", 1000), results$MAD_mean.RF)
    colnames(c) <- c ("n", "Method", "MAD") 
    
    d <- cbind.data.frame(rep(n, 1000), rep("BF", 1000), results$MAD_mean.BF)
    colnames(d) <- c ("n", "Method", "MAD") 
    
    e <- cbind.data.frame(rep(n, 1000), rep("OF", 1000), results[978])
    colnames(e) <- c ("n", "Method", "MAD") 
    
    f <- cbind.data.frame(rep(n, 1000), rep("CC", 1000), results$MAD_mean.CC)
    colnames(f) <- c ("n", "Method", "MAD") 
    
    g <- cbind.data.frame(rep(n, 1000), rep("MLM", 1000), results$MAD_mean.MLM)
    colnames(g) <- c ("n", "Method", "MAD") 
    
    h <- cbind.data.frame(rep(n, 1000), rep("MF", 1000), results$MAD_mean.MF)
    colnames(h) <- c ("n", "Method", "MAD") 
    
    final <- rbind.data.frame(a, b, c, d, e, f, g, h)
    
  }
  
  MAD_2000   <- create_long_MAD(2000, results2000)
  MAD_5000   <- create_long_MAD(5000, results5000)
  MAD_50000  <- create_long_MAD(50000, results50000)
  MAD_200000 <- create_long_MAD(200000, results200000)
  
  final_primary_outcome_MAD <- rbind.data.frame(MAD_2000, MAD_5000, MAD_50000, MAD_200000)
  
  
  library(ggplot2)
  
  p <- ggplot(final_primary_outcome_MAD, aes(x=n, y=MAD, group = interaction(n, Method), fill=Method)) + geom_boxplot(outlier.shape=NA, width=1) +
    scale_y_continuous(limits = c(0, 1.9)) + scale_x_continuous(trans='log10', breaks=c(2000, 5000, 50000, 200000), labels=scales::comma)
  
  p_colour <- p + scale_fill_manual( values=c("#56B4E9", "#F0E442","#E69F00",
                                               "#009E73", "#999999",
                                               "#D55E00", "#CC79A7", "#0072B2"),
                                      labels= c("Regression (correctly-specified)", "CART", "Cross-classification", 
                                                "CTree", "Regression (main effects - non-intersectional)", 
                                                "MAIHDA", "Regression (saturated)", "Random Forest"),
                                      guide = guide_legend(nrow=2, byrow=TRUE))
  
  p_final <- p_colour + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"),
                              legend.direction = "horizontal", legend.position = c(0.5,0.95))
  

  p_final
  
###########secondary outcomes for machine learning methods
  
  
  output_ML_final <- NULL
  
  for (a in a_list)  {
    
    results <- read.csv(a)
    
    #CART
    
    #1557: x6CART
    #1558: splitvarCART
      CARTx6  <- as.vector(apply(results[1557], 2, mean, na.rm=TRUE)) 
      
      results$splitx1 <- ifelse(grepl("x1", results[,1558])==TRUE, 1, 0)
      results$splitx2 <- ifelse(grepl("x2", results[,1558])==TRUE, 1, 0)
      results$splitx3 <- ifelse(grepl("x3", results[,1558])==TRUE, 1, 0)
      results$splitx4 <- ifelse(grepl("x4", results[,1558])==TRUE, 1, 0)
      results$splitx5 <- ifelse(grepl("x5", results[,1558])==TRUE, 1, 0)
      
      CARTx1_to_x5  <- as.vector(apply(results[,c("splitx1", "splitx2", "splitx3", "splitx4", "splitx5")], 2, mean, na.rm=TRUE))
      
    #CTree
    
    #1559: x6CTree
    #1560: splitvarCTree
      CTreex6  <- as.vector(apply(results[1559], 2, mean, na.rm=TRUE))
      
      results$splitx1 <- ifelse(grepl("x1", results[,1560])==TRUE, 1, 0)
      results$splitx2 <- ifelse(grepl("x2", results[,1560])==TRUE, 1, 0)
      results$splitx3 <- ifelse(grepl("x3", results[,1560])==TRUE, 1, 0)
      results$splitx4 <- ifelse(grepl("x4", results[,1560])==TRUE, 1, 0)
      results$splitx5 <- ifelse(grepl("x5", results[,1560])==TRUE, 1, 0)
      
      CTreex1_to_x5 <- as.vector(apply(results[,c("splitx1", "splitx2", 
                                                      "splitx3", "splitx4", "splitx5")], 2, mean, na.rm=TRUE))
      
      
    #Random Forest
      
    VIM.impurity<- as.vector(apply(results[c(1561:1566)], 2, mean, na.rm=TRUE))

    VIM.pvalue<- results[c(1573:1578)]
    VIM.pvalue.sig2.05 <- ifelse(VIM.pvalue < 0.05, 1, 0) 
    VIM.pvalue.sig.05 <- as.vector(apply((VIM.pvalue.sig2.05), 2, mean, na.rm=TRUE))
    
    
    #save results
    output_ML <- c(CARTx1_to_x5, CARTx6, CTreex1_to_x5, CTreex6, 
                    VIM.impurity, VIM.pvalue.sig.05)
    output_ML_final <- rbind(output_ML_final, output_ML) 
  }
  
  colnames(output_ML_final) <- c("CART_x1", "CART_x2", "CART_x3", 
                                 "CART_x4", "CART_x5","CART_x6",
                                 "CTree_x1", "CTree_x2", "CTree_x3", 
                                 "CTree_x4", "CTree_x5", "CTree_x6",
                                 "VIM.imp_x1", "VIM.imp_x2", "VIM.imp_x3", 
                                 "VIM.imp_x4", "VIM.imp_x5", "VIM.imp_x6",
                                 "VIM.perm_x1", "VIM.perm_x2", "VIM.perm_x3", 
                                 "VIM.perm_x4", "VIM.perm_x5", "VIM.perm_x6") 
  
  output_ML_final

  
  ###########convergence results for saturated regression
  
  output_covergence <- NULL
  
  for (a in a_list)  {
    results <- read.csv(a) 
    c <- table(results[1579])
    saturated.converged <- c
    output_covergence <- rbind(output_covergence,saturated.converged ) 
  }
  
  output_covergence