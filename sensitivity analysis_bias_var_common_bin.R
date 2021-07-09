################################################################################
#################### Sensitivity analysis, bias and variance ###################
########### Common binary outcome, all categorical inputs simulation ###########
################################################################################


Packages <- c("lsr", "foreach", "doParallel", "lme4","rpart", "partykit", 
              "tuneRanger", "mlr", "truncnorm", "merTools", "CHAID", "plyr",
              "parallel", "MASS")
lapply(Packages, library, character.only = TRUE)


#############functions for simulation

#set known Truth (needs to reset for each round)
createTruth <- function(Int, B1.1, B1.2, B1.3, B2, B3, B4, B5 ,Bint12, Bint13, Bint2) {
  
  x1<- c(0, 1, 2, 3)
  x2 <- c(0,1)
  x3 <- c(0,1)
  x4 <- c(0, 1)
  x5 <- c(0, 1)
  x6 <- c(0, 1, 2)
  
  #create every combo
  Truth <- expand.grid(x1=x1,x2=x2,x3=x3, x4=x4, x5=x5, x6=x6)
  
  
  Truth$cluster <- factor(100000*(Truth$x1+1) + 10000*(Truth$x2+1) + 1000*(Truth$x3+1) + 100*(Truth$x4+1) + 10*(Truth$x5 + 1) + 1*(Truth$x6+1))
  
  Truth$int12 <- 0
  Truth$int12 [Truth$x1==2 & Truth$x2==1] <- 1  
  
  Truth$int13 <- 0
  Truth$int13 [Truth$x1==3 & Truth$x2==1] <- 1  
  
  Truth$x1.1 <- ifelse(Truth$x1==1, 1, 0)
  Truth$x1.2 <- ifelse(Truth$x1==2, 1, 0)
  Truth$x1.3 <- ifelse(Truth$x1==3, 1, 0)
  
  Truth$z <- (Int + B1.1*Truth$x1.1 + B1.2*Truth$x1.2 + B1.3*Truth$x1.3 + B2*Truth$x2 + B3*Truth$x3 + B4*Truth$x4 + B5*Truth$x5 +  
                Bint12*Truth$int12 + Bint13*Truth$int13+ Bint2*Truth$x3*Truth$x4*Truth$x5)
  Truth$pr = exp(Truth$z)  # pass through log function
  
  Truth$x1<- as.factor(Truth$x1) 
  Truth$x2<- as.factor(Truth$x2) 
  Truth$x3<- as.factor(Truth$x3)
  Truth$x4<- as.factor(Truth$x4) 
  Truth$x5<- as.factor(Truth$x5) 
  Truth$x6 <- as.factor(Truth$x6)
  
  Truth
}



#############looping
n_list = c(2000, 5000)

registerDoParallel(2)


for (n in n_list){
  
  MADtable = NULL
  
  results = NULL
  
  results <- foreach (i=1:1000, .combine=rbind, .packages=c('rpart', 'partykit', 'tuneRanger', 'lme4', 'lsr', 'truncnorm', 'mlr', 'merTools', 'CHAID', 'plyr', 'sandwich', 'lmtest')) %dopar% {
    
    set.seed (i)
    
    D <- data.frame (x1 = sample(0:3, n, replace=TRUE),
                     x2 = sample(c(0,1),n, prob=c(0.80, 0.20), replace = TRUE),
                     x3 = sample(c(0,1),n, replace = TRUE),
                     x5 = sample(c(0,1), n, prob=c(0.75, 0.25), replace=TRUE),
                     x6 = sample(0:2, n, replace=TRUE))
    
    
    #create x3 as x4 mediator
    D$prX4 <- ifelse(D$x3<1,0.40, 0.7)
    D$x4=rbinom(n,1,D$prX4)
    
    
    D$x1 <- as.numeric(D$x1)
    D$x2 <- as.numeric(D$x2)
    D$x3 <- as.numeric(D$x3)
    D$x4 <- as.numeric(D$x4)
    D$x5 <- as.numeric(D$x5)
    D$x6 <- as.numeric(D$x6)
    
    
    #create intersection labels
    D$cluster <- factor(100000*(D$x1+1) + 10000*(D$x2+1) + 1000*(D$x3+1) + 100*(D$x4+1) + 10*(D$x5 + 1) + 1*(D$x6+1))
    
    
    D$int12 <- 0
    D$int12 [D$x1==2 & D$x2==1] <- 1  
    
    D$int13 <- 0
    D$int13 [D$x1==3 & D$x2==1] <- 1  
    
    D$x1.1 <- ifelse(D$x1==1, 1, 0)
    D$x1.2 <- ifelse(D$x1==2, 1, 0)
    D$x1.3 <- ifelse(D$x1==3, 1, 0)
    
    
    # natural log the bimodal function, so can create RR once exponentiated
    
    int <- -1.50
    b1.1 <- 0.4963259
    b1.2 <- 0.3327284
    b1.3 <- 0.3925075
    b2 <-  0.2672962
    b3 <- -0.5797873
    b4 <- -0.88203
    b5 <- -0.5378758
    b6 <- 0  
    bint12 <- -0.6366067
    bint13 <- -0.9815759
    bint2 <- 0.2179539
    
    D$z <- int + b1.1*D$x1.1 + b1.2*D$x1.2 + b1.3*D$x1.3+  b2*D$x2 + b3*D$x3 + b4*D$x4 + b5*D$x5 + bint12*D$int12 + bint13*D$int13+ bint2*D$x3*D$x4*D$x5
    D$pr = exp(D$z)  # pass through log function
    D$y = rbinom(n,1,D$pr) 
    preval <- mean(D$y)
    
    
    D$x1<- as.factor(D$x1) 
    D$x2<- as.factor(D$x2) 
    D$x3<- as.factor(D$x3)
    D$x4<- as.factor(D$x4) 
    D$x5<- as.factor(D$x5) 
    D$x6<- as.factor(D$x6) 
    
    Truth<- createTruth(int, b1.1, b1.2, b1.3, b2, b3, b4, b5, bint12,bint13, bint2)
    
    
    ##############Methods
    
    #main effects regression
    glmmain1<- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = D,
                   family = poisson(link = "log"))
    
    
    
    #best-specified regression
    glmbest1<- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x1*x2 +x3*x4*x5, data = D,
                   family = poisson(link = "log"))
    
   
    
    #cross classified
    crossclass <- aggregate(D$y==1, list(by=D$cluster), mean)
    colnames(crossclass) <- c("cluster","y")
    Dcrossclass <- join(Truth["cluster"], crossclass, by="cluster", type="left", match="first")
    #only merges for clusters that exist in D(:. won't have rows for clusters missing)
    
    
    #MAIHDA
    MLM<- glmer(y ~ x1 + x2 +x3 + x4 + x5  + x6 + (1|cluster), data=D, family=binomial,  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
    
    
    ##############intersection predictions function
    create_binom_predictions <- function(model, abbrev) {
      
      Bias <- (model- Truth$pr)
      predy <- (model)
      
      names2 <- c(1:192)
      names(Bias) <- c(paste(names2, rep(abbrev, 192), rep("bias", 192), sep="."))
      names(predy) <- c(paste(names2, rep(abbrev, 192), rep("pred", 192), sep="."))
      
      MAD_mean <- mean(abs(Bias)/abs(preval), na.rm=TRUE)
      
      names(MAD_mean) <- c(paste("MAD_mean", abbrev, sep="."))
      
      return(c(MAD_mean, Bias, predy))
      
    }
    
    ##############intersection predictions
    
    # to predict outcome, exponentiate
    BF.predict <- predict(glmbest1,Truth, se.fit=FALSE, type="response", interval = NULL)
    BF.predict.results <- create_binom_predictions(BF.predict, "BF")
    
    # to predict outcome, exponentiate
    MF.predict <- predict(glmmain1,Truth, se.fit=FALSE, type="response", interval = NULL)
    MF.predict.results <- create_binom_predictions(MF.predict, "MF")
    
    CC.predict <- Dcrossclass[[2]]
    CC.predict.results <- create_binom_predictions(CC.predict, "CC")
    
    MLM.predict1<- (predict (MLM, newdata=Truth, type="response", allow.new.levels = TRUE))
    # allow.new.levels statement for any clusters not included. Therefore remove any predictions made with clusters we don't have info on 
    MLM.predict2 <- ifelse(is.na(Dcrossclass$y)=="TRUE", NA, MLM.predict1)
    MLM.predict.results <- create_binom_predictions(MLM.predict2, "MLM")
    
    #create a sum of NA clusters... these clusters aren't calculated for cross.class variable, and MLM
    NAcluster <- sum(is.na(Dcrossclass[2]))
    
    
    MADtable = c(i,preval, int,  b1.1, b1.2, b1.3, b2, b3, b4, b5, bint12, bint13, bint2, NAcluster, 
                 BF.predict.results,  MF.predict.results,
                 CC.predict.results, MLM.predict.results
                 
    )
    
    
  }
  
  
  #then save the table to working directory
  write.csv(results,paste("Bin_common_cat_bias_var", n, "csv", sep = "."))   
  
  
}

