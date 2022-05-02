setwd("D:\\PUBLIKACJE\\____DepDegree\\CLEAN-DATA")
library(Rnalytica)
library(corrplot)
library(caret)


#############################
# Load and prepare the data #
#############################

d <- read.csv("github-UBD.csv")

#dependent variable 'bug' as factor
d$bug <- as.factor(d$bug)
levels(d$bug) <- c("bug","clean")

#log(x+1) except for density metrics
d[,-c(1, 2, 4, 8, 9, 21, 22, 27, 66, 67)] <- log(d[,-c(1, 2, 4, 8, 9, 21, 22, 27, 66, 67)]+1)

# define sets of variables
indVars <- names(d)[c(3:63)]   # names of independent variables
indDFvars <- names(d)[c(64:65)] # names of independent data-flow variables
depVar <- names(d)[67]      # name of dependent variable


###############################
# Functions used in the study #
###############################

splitdata <- function(d, pr) {
  trainRowNumbers <- createDataPartition(d$bug, p=pr, list=FALSE)
  trainData <- d[trainRowNumbers,]
  testData <- d[-trainRowNumbers,]
  rv <- list("trainData"=trainData, "testData"=testData)
  return(rv)
}

autoSp <- function(d, features, thr, vif) {
  res <- AutoSpearman(dataset = d, metrics = features, 
                             spearman.threshold = thr, 
                             vif.threshold = vif, 
                             groups = FALSE, 
                             verbose = F)
  return(res)
}


########################
# Correlation analysis #
########################

# Identify uncorrelated metrics
aSp <- autoSp(d, c(indVars, indDFvars), thr=0.7, vif=5)

#correlation plots
res <- cor(as.matrix(d[, aSp]))
corrplot(res, order = 'hclust', addCoef.col = 'black', tl.pos = 'd',
         cl.pos = 'n', col = COL1('Blues'), method='number', 
         tl.cex = 0.7, number.cex = 0.7, type="upper")

res <- cor(as.matrix(d[, c(indVars, "DD", "DDD")]))
plot(res[1:(nrow(res)-2),c("DD", "DDD")], xlim=c(0, 1), ylim=c(0,1))
abline(a=0, b=1)


###############################################
# Within-Project Defect Prediction Experiment #
###############################################

# control parameters for caret train() function
fitControl <- trainControl(
  method = 'boot632',              # .632 bootstrap
  savePredictions = 'final',       # saves predictions for optimal tuning parameter
  classProbs = TRUE,               # should class probabilities be returned
  summaryFunction=twoClassSummary  # results summary function
) 

# choose models and projects from the Unified Bug Dataset
MODELS <- c('naive_bayes', 'glm','rf')
PROJECTS <- c('broadleaf-3.0.10', 'elasticsearch-0.90.11', 'hazelcast-3.3',
              'mcMMO-1.4.06', 'netty')

filename <- 'results.csv'         # output file with experiment results
N <- 100                          # number of iterations

cat("Iter;",sep="",file=filename,append=TRUE)
cat(indVars,sep=";",file=filename,append=TRUE)
cat(";",file=filename,append=TRUE)
cat(indDFvars,sep=";",file=filename,append=TRUE)
cat(";Project;Model;F1X;F1D;PrecX;PrecD;RecallX;RecallD\n",sep="",file=filename,append=TRUE)

for (p in PROJECTS) {
  
  # filter by project and remove any non-variable column
  d2 <- d[d$ProjectName==p,c(indVars, indDFvars, depVar)]
  
  nbugs <- nrow(d2[d2$bug=="bug",])
  nclean <- nrow(d2[d2$bug=="clean",])
  
  for (i in 1:N) {
  cat("proj=",p," clean=",nclean," bugs=",nbugs, " iter=",i, "\n",sep="")
  
  # generate training and testing samples
  spl <- splitdata(d2, 0.8)
  trainData <- spl[["trainData"]]
  testData <- spl[["testData"]]
  
  # let AutoSpearman select uncorrelated metrics for this project data
  VARS2 <- autoSp(trainData, c(indVars, indDFvars), thr=0.7, vif=5)
  # VARS will have the same set of variables with data-flow metrics removed
  VARS <- setdiff(VARS2, indDFvars)
  cat("#VARS=",length(VARS)," #VARS2=",length(VARS2),"\n")
  
  # it may happen that neither DD nor DDD will be selected
  # in this case ignore the iteration [in our experiments it never happened]
  if (length(VARS2) > length(VARS))  
  {
  
  #construct defect models for selected variables with and without a given data flow metric
  for (m in MODELS) {
    mdl = train(x = as.data.frame(trainData[,VARS]), y = trainData$bug,  
                method=m, tuneLength = 5, metric='ROC', trControl = fitControl)
    mdl2 = train(x = as.data.frame(trainData[,VARS2]), y = trainData$bug,  
                 method=m, tuneLength = 5, metric='ROC', trControl = fitControl)
    
    predicted <- predict(mdl, testData)
    predicted2 <- predict(mdl2, testData)
    
    F1X <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
    F1D <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
    PrecX <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
    PrecD <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
    RecallX <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]
    RecallD <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]

    cat(i,";",sep="",file=filename,append=TRUE)
    cat(1*(indVars %in% VARS),sep=";",file=filename,append=TRUE)
    cat(";",sep="",file=filename,append=TRUE)
    
    # write down which metrics were selected by AutoSpearman
    cat(1*(indDFvars %in% setdiff(VARS2, VARS)),sep=";",file=filename,append=TRUE)
    
    cat(";",p,sep="",file=filename,append=TRUE)
    cat(";",m,";",sep="",file=filename,append=TRUE)
    cat(F1X,";",F1D,";",PrecX,";",PrecD,";",RecallX,";",RecallD,"\n",sep="",file=filename,append=TRUE)
    
  }
  } else {
    cat("AutoSpearman did not select any data-flow metric\n")
  }
}
}

# summarize the results
df <- read.csv("results.csv", sep=";")
aggregate(cbind(F1X, F1D, PrecX, PrecD, RecallX, RecallD) ~ Model+Project, data=df, FUN=mean)



# AutoSpearman used always DDD and never DD; perform similar computations but, for DD

filename <- 'resultsDD.csv'       # output file with experiment results
N <- 100                          # number of iterations

cat("Iter;",sep="",file=filename,append=TRUE)
cat(indVars,sep=";",file=filename,append=TRUE)
cat(";",file=filename,append=TRUE)
cat(indDFvars,sep=";",file=filename,append=TRUE)
cat(";Project;Model;F1X;F1D;PrecX;PrecD;RecallX;RecallD\n",sep="",file=filename,append=TRUE)

for (p in PROJECTS) {
  
  # filter by project and remove any non-variable column
  d2 <- d[d$ProjectName==p,c(indVars, indDFvars, depVar)]
  
  nbugs <- nrow(d2[d2$bug=="bug",])
  nclean <- nrow(d2[d2$bug=="clean",])
  
  for (i in 1:N) {
  cat("proj=",p," clean=",nclean," bugs=",nbugs, " iter=",i, "\n",sep="")
  
  # generate training and testing samples
  spl <- splitdata(d2, 0.8)
  trainData <- spl[["trainData"]]
  testData <- spl[["testData"]]
  
  # let AutoSpearman select uncorrelated metrics for this project data, excluding data flow metrics
  VARS <- autoSp(trainData, indVars, thr=0.7, vif=5)
  # VARS2 will have the same set of variables with DD as additional metric
  VARS2 <- c(VARS, "DD")
  cat("#VARS=",length(VARS)," #VARS2=",length(VARS2),"\n")
  
  
  #construct defect models for selected variables with and without a given data flow metric
  for (m in MODELS) {
    mdl = train(x = as.data.frame(trainData[,VARS]), y = trainData$bug,  
                method=m, tuneLength = 5, metric='ROC', trControl = fitControl)
    mdl2 = train(x = as.data.frame(trainData[,VARS2]), y = trainData$bug,  
                 method=m, tuneLength = 5, metric='ROC', trControl = fitControl)
    
    predicted <- predict(mdl, testData)
    predicted2 <- predict(mdl2, testData)
    
    F1X <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
    F1C <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
    PrecX <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
    PrecC <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
    RecallX <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]
    RecallD <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]

    cat(i,";",sep="",file=filename,append=TRUE)
    cat(1*(indVars %in% VARS),sep=";",file=filename,append=TRUE)
    cat(";",sep="",file=filename,append=TRUE)
    
    # write down which metrics were selected by AutoSpearman
    cat(1*(indDFvars %in% setdiff(VARS2, VARS)),sep=";",file=filename,append=TRUE)
    
    cat(";",p,sep="",file=filename,append=TRUE)
    cat(";",m,";",sep="",file=filename,append=TRUE)
    cat(F1X,";",F1D,";",PrecX,";",PrecD,";",RecallX,";",RecallD,"\n",sep="",file=filename,append=TRUE)
    
  }
}
}

# summarize the results
df <- read.csv("results.csv", sep=";")
aggregate(cbind(F1X, F1D, PrecX, PrecD, RecallX, RecallD) ~ Model+Project, data=df, FUN=mean)



##############################################
# Cross-Project Defect Prediction Experiment #
##############################################

MODELS <- c('naive_bayes', 'glm','rf')
PROJECTS <- c('broadleaf-3.0.10', 'elasticsearch-0.90.11', 'hazelcast-3.3',
              'mcMMO-1.4.06', 'netty')

filename <- 'results-cross.csv'

cat("Iter;",sep="",file=filename,append=TRUE)
cat(indVars,sep=";",file=filename,append=TRUE)

cat(";Project;Model;crF1X;crF1DD;crF1DDD;crPrecX;crPrecDD;crPrecDDD;crRecX;crRecDD;crRecDDD\n",sep="",file=filename,append=TRUE)

for (p in PROJECTS) {

  cat("proj=",p," iter=",i, "\n",sep="")

  # generate training and testing samples
  trainData <- d[d$ProjectName==p,]
  testData <- d[d$ProjectName!=p,]
  
  # let AutoSpearman select uncorrelated metrics for this project data
  # excluding DDD
  VARS <- autoSp(trainData, indVars, thr=0.7, vif=5)
  # VARS2 = VARS with DDD added
  VARS2 <- c(VARS, "DD")
  VARS3 <- c(VARS, "DDD")

  #construct defect models for selected variables with and without a given data flow metric
  for (m in MODELS) {
    
    # train on a given project
    mdl = train(x = as.data.frame(trainData[,VARS]), y = trainData$bug,  
                method=m, tuneLength = 5, metric='ROC', trControl = fitControl)
    mdl2 = train(x = as.data.frame(trainData[,VARS2]), y = trainData$bug,  
                 method=m, tuneLength = 5, metric='ROC', trControl = fitControl)
    mdl3 = train(x = as.data.frame(trainData[,VARS3]), y = trainData$bug,  
                 method=m, tuneLength = 5, metric='ROC', trControl = fitControl)    
    
	# predict on other projects    
    predicted <- predict(mdl, testData)
    predicted2 <- predict(mdl2, testData)
    predicted3 <- predict(mdl2, testData)
	
    # calculate cross-project performance for both models (without and with DDD)     
    F1X <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
    F1DD <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
	  F1DDD <- confusionMatrix(predicted3, testData$bug, positive="bug", mode="everything")[[4]]["F1"]
    PrecX <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
    PrecDD <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
	  PrecDDD <- confusionMatrix(predicted3, testData$bug, positive="bug", mode="everything")[[4]]["Precision"]
    RecallX <- confusionMatrix(predicted, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]
    RecallDD <- confusionMatrix(predicted2, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]
	  RecallDDD <- confusionMatrix(predicted3, testData$bug, positive="bug", mode="everything")[[4]]["Recall"]
      
    cat(i,";",sep="",file=filename,append=TRUE)
    cat(1*(indVars %in% VARS),sep=";",file=filename,append=TRUE)

    cat(";",p,sep="",file=filename,append=TRUE)
    cat(";",m,";",sep="",file=filename,append=TRUE)
    cat(F1X,";",F1DD,";",F1DDD,";",PrecX,";",PrecDD,";",PrecDDD,";",RecallX,";",RecallDD,";",RecallDDD,
        "\n",sep="",file=filename,append=TRUE)
        
      }
}

# summarize the results
df <- read.csv("results-cross.csv", sep=";")
aggregate(cbind(crF1X, crF1DD, crF1DDD, crPrecX, crPrecDD, crPrecDDD, crRecX, crRecDD, crRecDDD) ~ Model+Project, data=df, FUN=mean)
