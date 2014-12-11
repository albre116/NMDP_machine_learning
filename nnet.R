####################################################
####neural network approach
####this approach will also do the L1-L2 regularized sparse vector machine
####based upon kernel methods (should be roughly equivalent to SVM)
####################################################
rm(list=ls())
gc()
graphics.off()
op<-par()
if(!require("MASS"))
  (install.packages("MASS"))
if(!require("ggplot2"))
  (install.packages("ggplot2"))
if(!require("kernlab"))
  (install.packages("kernlab"))
if(!require("e1071"))
  (install.packages("e1071"))
if(!require("gbm"))
  (install.packages("gbm"))
if(!require("dplyr"))
  (install.packages("dplyr"))
if(!require("reshape2"))
  (install.packages("reshape2"))
if(!require("nnet"))
  (install.packages("nnet"))
if(!require("RCurl"))
  (install.packages("RCurl"))
if(!require("devtools"))
  (install.packages("devtools"))
if(!require("reshape"))
  (install.packages("reshape"))
if(!require("RSNNS"))
  (install.packages("RSNNS"))
if(!require("Rcpp"))
  (install.packages("Rcpp"))
if(!require("neuralnet"))
  (install.packages("neuralnet"))

###Import Data
RAW_DATA<-read.csv("Raw_Data.csv")
RAW_DATA$Frequency.1*RAW_DATA$Frequency.2
permute<-rbinom(nrow(RAW_DATA),1,prob=0.5)==1
RAW_DATA[permute,]<-RAW_DATA[permute,c(1:6,8,7,10,9,11:13)]
DATA<-RAW_DATA[c("RID","Race","Frequency.1",
                 "Frequency.2","Rank.1","Rank.2",
                 "Total.Genotype.Frequency","Productivity.Group")]

DATA$Productivity.Group<-ordered(DATA$Productivity.Group,levels=c("D","C","B","A"))
colnames(DATA)<-c("RID","Race","H1","H2","Rank_H1","Rank_H2","GF","Productivity")
DATA<-DATA[complete.cases(DATA),]###remove NA's
DATA<-DATA[DATA$GF!=0 & DATA$H1!=0 & DATA$H2!=0,]
rm(RAW_DATA)

###Format Predictors to be ln-scale
DATA$GF<-log(DATA$GF)
DATA$H1<-log(DATA$H1)
DATA$H2<-log(DATA$H2)




####Split Train and Test
set.seed(1103)
train_idx<-sample(1:nrow(DATA),floor(nrow(DATA)*0.8))
logical<-rep(FALSE,nrow(DATA))
logical[train_idx]<-TRUE
train_idx<-logical
test_idx<-!logical
rm(logical)
TRAIN<-DATA[train_idx,]
TEST<-DATA[test_idx,]


x<-TRAIN[c("H1","H2")]
x<-cbind(x,decodeClassLabels(TRAIN[,c("Race")]))
#x$intercept=1
y<-TRAIN[c("Productivity")]
y<-decodeClassLabels(y[,1])
x_tst<-TEST[c("H1","H2")]
x_tst<-cbind(x_tst,decodeClassLabels(TEST[,c("Race")]))
#x_tst$intercept=1
y_tst<-TEST[c("Productivity")]
y_tst<-decodeClassLabels(y_tst[,1])

bestmod <- mlp(x,y, size=c(3,5), learnFuncParams=c(0.1),
             maxit=50)
summary(bestmod)
bestmod$snnsObject$getUnitDefinitions()
weightMatrix(bestmod)
extractNetInfo(bestmod)
predictions <- predict(bestmod,x_tst)
plotRegressionError(predictions[,2], y_tst[,2])
confusionMatrix(y_tst,predictions)
plotROC(predictions[,2],y_tst[,2])



#import plot function from Github
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')

plot.nnet(bestmod)
plot.nnet(bestmod,nid=F)###don't show thickness
plot.nnet(bestmod,wts.only=T)

#import sensitivity analysis function
source_url('https://gist.githubusercontent.com/fawda123/6860630/raw/b8bf4a6c88d6b392b1bfa6ef24759ae98f31877c/lek_fun.r')

#sensitivity analsyis, note 'exp.in' argument
lek.fun(bestmod,exp.in=x,var.sens=c("X1","X2"),
        split.vals=c(0,.50,1))









