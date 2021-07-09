################################################################################
################ Common binary outcome, mixed inputs simulation ################
################################################################################

Packages <- c("lsr", "foreach", "doParallel", "lme4","rpart", "partykit", 
              "tuneRanger", "mlr", "truncnorm", "merTools", "CHAID", "plyr",
              "parallel", "MASS")
lapply(Packages, library, character.only = TRUE)



#############functions for simulation

  #bimodal dist for coeff
    Sample.negative.coeff <- function (n) {
      y0 <- rtruncnorm(n,a=0.20, b=0.89, mean=((0.89+0.2)/2), sd = 0.30)
    }
    
    Sample.positive.coeff <- function (n) {
      y0 <- rtruncnorm(n,a=1.11, b=1.80, mean=((1.11+1.80)/2), sd = 0.30)
    }

  
  #SET KNOWN TRUTH (NEEDS TO RESET FOR EACH ROUND)
    createTruth <- function(Int, B1, B2, B3, B4, B5 ,Bint1, Bint2) {
      
      
      x1<- c(-1.273, -0.325, 0.325, 1.273)
      x2 <- c(0,1)
      x3 <- c(0,1)
      x4 <- c(0, 1)
      x5 <- c(0, 1)
      x6 <- c(-1.090, 0, 1.090)
      
      #create every combo
      Truth <- expand.grid(x1=x1,x2=x2,x3=x3, x4=x4, x5=x5, x6=x6)
      
      Truth$x1cat <- quantileCut(Truth$x1, 4, labels(1:4)) 
      Truth$x6cat <- quantileCut(Truth$x6, 3, labels(1:3))
      Truth$x1cat <- as.numeric(Truth$x1cat) 
      Truth$x6cat <- as.numeric(Truth$x6cat) 
      
      Truth$cluster <- factor(100000*(Truth$x1cat+1) + 10000*(Truth$x2+1) + 1000*(Truth$x3+1) + 100*(Truth$x4+1) + 10*(Truth$x5 + 1) + 1*(Truth$x6cat+1))
      
      Truth$int1 <- 0
      Truth$int1 [Truth$x1>1 & Truth$x2==1] <- 1
      Truth$int2 <- 0  
      Truth$int2 [Truth$x3==1 & Truth$x4==1 & Truth$x5==1] <- 1#INTERACTION AFTER CUTOFF
      
      
      Truth$z <- (Int + B1*Truth$x1  + B2*Truth$x2 + B3*Truth$x3 + B4*Truth$x4 + B5*Truth$x5 +  
                    Bint1*Truth$int1*Truth$x1*Truth$x2 + Bint2*Truth$int2*Truth$x3*Truth$x4*Truth$x5)
      Truth$pr = exp(Truth$z)  # pass through log function
      
      Truth$x1cat<- as.factor(Truth$x1cat) 
      Truth$x2<- as.factor(Truth$x2) 
      Truth$x3<- as.factor(Truth$x3)
      Truth$x4<- as.factor(Truth$x4) 
      Truth$x5<- as.factor(Truth$x5) 
      Truth$x6cat <- as.factor(Truth$x6cat)
      
      Truth
    }



#############looping
n_list = c(2000, 5000, 50000, 200000)

registerDoParallel(5)


for (n in n_list){
  
  MADtable = NULL
  MADtable <- matrix(nrow=1000, ncol=1580)
  
  results = NULL
  results <- matrix(nrow=1000, ncol=1580)
  
  results <- foreach (i=1:1000, .combine=rbind, .packages=c('rpart', 'partykit', 'tuneRanger', 'lme4', 'lsr', 'truncnorm', 'mlr', 'merTools', 'CHAID', 'plyr')) %dopar% {
    
    set.seed (i)
    
    D <- data.frame (x1 = rnorm(n,0,1),
                     x2 = sample(c(0,1),n, prob=c(0.80, 0.20), replace = TRUE),
                     x3 = sample(c(0,1),n, replace = TRUE),
                     x5 = sample(c(0,1), n, prob=c(0.75, 0.25), replace=TRUE),
                     x6 = rnorm(n, 0, 1))
    
    
    #create x3 as x4 mediator
    D$prX4 <- ifelse(D$x3<1,0.40, 0.7)
    D$x4=rbinom(n,1,D$prX4)
    
    #create categories for cts variables
    D$x1cat <- quantileCut(D$x1, 4, labels(1:4))
    D$x6cat <- quantileCut(D$x6, 3, labels(1:3))
    
    
    D$x1cat <- as.numeric(D$x1cat)
    D$x2 <- as.numeric(D$x2)
    D$x3 <- as.numeric(D$x3)
    D$x4 <- as.numeric(D$x4)
    D$x5 <- as.numeric(D$x5)
    D$x6cat <- as.numeric(D$x6cat)
    
    
    #create intersection labels
    D$cluster <- factor(100000*(D$x1cat+1) + 10000*(D$x2+1) + 1000*(D$x3+1) + 100*(D$x4+1) + 10*(D$x5 + 1) + 1*(D$x6cat+1))
    
    D$int1 <- 0
    D$int1 [D$x1>1 & D$x2==1] <- 1
    D$int2 <- 0  
    D$int2 [D$x3==1 & D$x4==1 & D$x5==1] <- 1 #INTERACTION b/w binaries
    
    
    # natural log the bimodal function, so can create RR once exponentiated
    int <- -1.5
    b1 <- log(Sample.positive.coeff(1))
    b2 <-  log(Sample.positive.coeff(1))
    b3 <- log(Sample.negative.coeff(1))
    b4 <- log(Sample.negative.coeff(1))
    b5 <- log(Sample.negative.coeff(1))
    b6 <- 0  
    bint1 <- log(Sample.negative.coeff(1))
    bint2 <- log(Sample.positive.coeff(1))
    
    D$z <- int + b1*D$x1 + b2*D$x2 + b3*D$x3 + b4*D$x4 + b5*D$x5 + + bint1*D$int1*D$x1*D$x2 + bint2*D$int2*D$x3*D$x4*D$x5
    D$pr = exp(D$z)  # pass through log function
    D$y = rbinom(n,1,D$pr) 
    preval <- mean(D$y)
    
    #resample if needed so probability will not exceed 1
    while (is.na(preval)==TRUE) {
      b1 <- log(Sample.positive.coeff(1))
      b2 <-  log(Sample.positive.coeff(1))
      b3 <- log(Sample.negative.coeff(1))
      b4 <- log(Sample.negative.coeff(1))
      b5 <- log(Sample.negative.coeff(1))
      b6 <- 0  
      bint1 <- log(Sample.negative.coeff(1))
      bint2 <- log(Sample.positive.coeff(1))
      
      D$z <- int + b1*D$x1 + b2*D$x2 + b3*D$x3 + b4*D$x4 + b5*D$x5 + + bint1*D$int1*D$x1*D$x2 + bint2*D$int2*D$x3*D$x4*D$x5
      D$pr = exp(D$z)  # pass through log function
      D$y = rbinom(n,1,D$pr) 
      preval <- mean(D$y)}
    
    D$x1cat<- as.factor(D$x1cat) 
    D$x2<- as.factor(D$x2) 
    D$x3<- as.factor(D$x3)
    D$x4<- as.factor(D$x4) 
    D$x5<- as.factor(D$x5) 
    D$x6cat<- as.factor(D$x6cat) 
    
    Truth<- createTruth(int,b1, b2, b3, b4, b5, bint1, bint2)
    
    
    ##############Methods
    
    #main effects regression
    glmmain1<- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = D,
                   family = poisson(link = "log"))
    
      #collect fit statistics
      mainfit.converged <- glmmain1$converged
    
    #best-specified regression
    glmbest1<- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x1*x2 +x3*x4*x5, data = D,
                   family = poisson(link = "log"))
    
      #collect fit statistics
      bestfit.converged <- glmbest1$converged
    
    
    #saturated regression
    Overfit <- try( glm(y ~ (x1 + x2 + x3 + x4 + x5 + x6)^6, data = D,
                        family = poisson(link = "log")), silent=TRUE)
    if (class(Overfit)=="try-error") {
      Overfit.converged <- NA
    } else {
      
      Overfit.converged <- Overfit$converged
      
    }
    
    #cross classified
    crossclass <- aggregate(D$y==1, list(by=D$cluster), mean)
    colnames(crossclass) <- c("cluster","y")
    Dcrossclass <- join(Truth["cluster"], crossclass, by="cluster", type="left", match="first")
    #only merges for clusters that exist in D(:. won't have rows for clusters missing)
    
    
    #MAIHDA
    MLM<- glmer(y ~ x1 + x2 +x3 + x4 + x5  + x6 + (1|cluster), data=D, family=binomial,  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
    
    #For decision tree methods, convert y to factor 
    D$y <- as.factor (D$y)
    
    #CART
    CART <-rpart(y~x1+x2+x3+x4+x5+x6,D, method="class")
    CARTpruned<- prune(CART, cp=CART$cptable[which.min(CART$cptable[,"xerror"]),"CP"])
    splitvarCART <- paste(unique((CARTpruned$frame[,1])), collapse=".")
    x6CART <- ifelse(grepl("x6", splitvarCART)=="TRUE", 1, 0)
    
    #CTree 
    CTree <- ctree (y~x1+x2+x3+x4+x5+x6,data = D, control = ctree_control(mincriterion = 0.95))
    splitvarCTree <- paste(names(varimp(CTree)), collapse=".") 
    x6CTree <- ifelse(grepl("x6", splitvarCTree)=="TRUE", 1, 0) 
    
    #Random forest
    RFmodel <- tuneMtryFast(formula = y~x1+x2+x3+x4+x5+x6, data = D, stepFactor = 1, improve = 0.05,
                            trace = TRUE, plot = TRUE, doBest = T, importance = "impurity", seed=i, probability=TRUE)
    newRFimp <- importance(RFmodel)
    
      #altmann method for VIM only recorded for first 200 iterations
      if (i < 201){
        RFmodel2 <- tuneMtryFast(formula = y~x1+x2+x3+x4+x5+x6, data = D, stepFactor = 1, improve = 0.05,
                                 trace = TRUE, plot = TRUE, doBest = T, importance = "permutation", seed=i, probability=TRUE)
        RFimp2<- importance_pvalues(RFmodel2 , method = 'altmann',    formula = y~x1+x2+x3+x4+x5+x6, data = D, seed=i)
        newRFimp2 <- as.vector(RFimp2) 
      } else {
        newRFimp2 <-  rep(NA, 12)
      }
    
    
    
    ##############intersection predictions function
    create_binom_predictions <- function(model, abbrev) {
      
      Ipredy <- (model- Truth$pr)
      
      names2 <- c(1:192)
      names(Ipredy) <- c(paste(names2, rep(abbrev, 192), sep="."))
      
      MAD_mean <- mean(abs(Ipredy)/abs(preval), na.rm=TRUE)

      names(MAD_mean) <- c(paste("MAD_mean", abbrev, sep="."))
      
      return(c(MAD_mean, Ipredy))
      
    }
    
    ##############intersection predictions
    CART.predict <- (predict(CARTpruned, Truth))[,2]
    CART.predict.results <- create_binom_predictions(CART.predict, "CART")
    
    CTree.predict <- predict(CTree, Truth, type="prob")[,2]
    CTree.predict.results <- create_binom_predictions(CTree.predict, "CTree")
    
    RF.predict <- predict(RFmodel, Truth)$predictions[,2]
    RF.predict.results <- create_binom_predictions(RF.predict, "RF")
    
    # to predict outcome, exponentiate
    BF.predict <- predict(glmbest1,Truth, se.fit=FALSE, type="response", interval = NULL)
    BF.predict.results <- create_binom_predictions(BF.predict, "BF")
    
    # to predict outcome, exponentiate
    MF.predict <- predict(glmmain1,Truth, se.fit=FALSE, type="response", interval = NULL)
    MF.predict.results <- create_binom_predictions(MF.predict, "MF")
    
    if (class(Overfit)=="try-error") { 
      OF.predict.results <- rep(NA, 193)} else {
        OF.predict <- predict(Overfit,Truth, se.fit=FALSE, type="response", interval = NULL)
        OF.predict.results <- create_binom_predictions(OF.predict, "OF") }
    
    CC.predict <- Dcrossclass[[2]]
    CC.predict.results <- create_binom_predictions(CC.predict, "CC")
    
    MLM.predict1<- (predict (MLM, newdata=Truth, type="response", allow.new.levels = TRUE))
    # allow.new.levels statement for any clusters not included. Therefore remove any predictions made with clusters we don't have info on 
    MLM.predict2 <- ifelse(is.na(Dcrossclass$y)=="TRUE", NA, MLM.predict1)
    MLM.predict.results <- create_binom_predictions(MLM.predict2, "MLM")
    
    #create a sum of NA clusters... these clusters aren't calculated for cross.class variable, and MLM
    NAcluster <- sum(is.na(Dcrossclass[2]))
    
    
    MADtable = c(i,preval, int,  b1, b2, b3, b4, b5, bint1, bint2, NAcluster, 
                 CART.predict.results, CTree.predict.results, RF.predict.results,
                 BF.predict.results, MF.predict.results, OF.predict.results, 
                 CC.predict.results, MLM.predict.results,
                 x6CART, splitvarCART, 
                 x6CTree , splitvarCTree,
                 newRFimp, newRFimp2, 
                 Overfit.converged, bestfit.converged, mainfit.converged   
    )
    
    
  }
  
  
  
  #then save the table to working directory
  write.csv(results,paste("Bin_common_mixed", n, "csv", sep = "."))   
  
  #MADTABLE: 1580 in width
  
}
