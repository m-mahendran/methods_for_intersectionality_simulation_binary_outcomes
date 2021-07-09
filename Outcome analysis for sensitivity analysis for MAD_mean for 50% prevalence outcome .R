
################################################################################
########### Results analysis, for sensitivity analysis of MAD_mean of ##########
###################### binary outcome with 50% prevalence ######################
################################################################################
#run for each of the 4 sets to create MAD_mean boxplots: 

#Common binary outcome, categorical inputs, mixed level of interaction and main effects
results2000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval.2000.csv")
results5000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval.5000.csv")

#Common binary outcome, categorical inputs, interaction greater than main effects
results2000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval_v2.2000.csv")
results5000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval_v2.5000.csv")

#Common binary outcome, mixed inputs, mixed level of interaction and main effects
results2000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval_mixed.2000.csv")
results5000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval_mixed.5000.csv")

#Common binary outcome, mixed inputs, interaction greater than main effects
results2000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval_mixed_v2.2000.csv")
results5000<- read.csv ("T:/Bauer_grp/Intersectionality/Project- Simulation/publication 2 - categorical outcomes/new analyses - higher prevalence/Bin_50_preval_mixed_v2.5000.csv")



#####################################################################

# run analysis below for each of the 4 sets:

###########

create_long_MAD <- function(n, results) {
  
  a <- cbind.data.frame(rep(n, 1000), rep("BF", 1000), results$MAD_mean.BF)
  colnames(a) <- c ("n", "Method", "MAD") 
  
  b <- cbind.data.frame(rep(n, 1000), rep("MLM", 1000), results$MAD_mean.MLM)
  colnames(b) <- c ("n", "Method", "MAD") 
  
  c <- cbind.data.frame(rep(n, 1000), rep("MF", 1000), results$MAD_mean.MF)
  colnames(c) <- c ("n", "Method", "MAD") 
  
  d <- cbind.data.frame(rep(n, 1000), rep("CC", 1000), results$MAD_mean.CC)
  colnames(d) <- c ("n", "Method", "MAD") 
  
  final <- rbind.data.frame( a, b, c, d)
  
}

MAD_2000   <- create_long_MAD(2000, results2000)
MAD_5000   <- create_long_MAD(5000, results5000)

final_primary_outcome_MAD <- rbind.data.frame(MAD_2000, MAD_5000)
final_primary_outcome_MAD <- na.omit(final_primary_outcome_MAD)


#plot 

library(ggplot2)

p <- ggplot(final_primary_outcome_MAD, aes(x=n, y=MAD, group = interaction(n, Method), fill=Method)) + geom_boxplot(outlier.shape=NA, width=1) +
  scale_y_continuous(limits = c(0, 0.5)) + scale_x_continuous(trans='log10', breaks=c(2000, 5000), labels=scales::comma)

p_colour <- p + scale_fill_manual( values=c("#56B4E9", "#E69F00",
                                            "#999999",
                                            "#D55E00"),
                                   labels= c("Regression (correctly-specified)",  "Cross-classification", 
                                             "Regression (main effects - non-intersectional)", 
                                             "MAIHDA"),
                                   guide = guide_legend(nrow=2, byrow=TRUE))


p_final <- p_colour + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"),
                            legend.direction = "horizontal", legend.position = c(0.5,0.95))


p_final

