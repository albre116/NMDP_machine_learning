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
if(!require("glmnet"))
  (install.packages("glmnet"))


###Import Data
RAW_DATA<-read.csv("Raw_Data.csv")
RAW_DATA$Productivity.Group<-as.character(RAW_DATA$Productivity.Group)
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="B"]="A"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="C"]="B"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="D"]="C"


min_thresh<-min(RAW_DATA$Total.Genotype.Frequency[ RAW_DATA$Total.Genotype.Frequency!=0],na.rm=T)
RAW_DATA$Total.Genotype.Frequency[ RAW_DATA$Total.Genotype.Frequency==0]<-min_thresh
RAW_DATA$Total.Genotype.Frequency[ is.na(RAW_DATA$Total.Genotype.Frequency)]=min_thresh

min_thresh<-min(RAW_DATA$Frequency.1[ RAW_DATA$Frequency.1!=0],na.rm=T)
RAW_DATA$Frequency.1[ RAW_DATA$Frequency.1==0]<-min_thresh
RAW_DATA$Frequency.1[ is.na(RAW_DATA$Frequency.1)]=min_thresh

min_thresh<-min(RAW_DATA$Frequency.2[ RAW_DATA$Frequency.2!=0],na.rm=T)
RAW_DATA$Frequency.2[ RAW_DATA$Frequency.2==0]<-min_thresh
RAW_DATA$Frequency.2[ is.na(RAW_DATA$Frequency.2)]=min_thresh


permute<-rbinom(nrow(RAW_DATA),1,prob=0.5)==1
RAW_DATA[permute,]<-RAW_DATA[permute,c(1:6,8,7,10,9,11:13)]
DATA<-RAW_DATA[c("RID","Race","Frequency.1",
                 "Frequency.2","Rank.1","Rank.2",
                 "Total...of.potential.10.10s.......50.",
                 "Total...of.potential.9.10s......50.",
                 "Total.Genotype.Frequency","Productivity.Group")]

DATA$Productivity.Group<-ordered(DATA$Productivity.Group,levels=c("C","B","A"))
colnames(DATA)<-c("RID","Race","H1","H2","Rank_H1","Rank_H2","ten_of_tens","nine_of_tens","GF","Productivity")
DATA<-DATA[complete.cases(DATA),]###remove NA's
DATA<-DATA[DATA$GF!=0 & DATA$H1!=0 & DATA$H2!=0,]
rm(RAW_DATA)

###Format Predictors to be ln-scale
DATA$GF<-log(DATA$GF)
DATA$H1<-log(DATA$H1)
DATA$H2<-log(DATA$H2)
DATA$ten_of_tens<-log(DATA$ten_of_tens+1)
DATA$nine_of_tens<-log(DATA$nine_of_tens+1)

DATA[c("H1","H2","GF","ten_of_tens","nine_of_tens")]<-normalizeData(DATA[c("H1","H2","GF","ten_of_tens","nine_of_tens")],type="0_1")
summary(DATA)

####Split Train and Test
set.seed(1103)
train_idx<-sample(1:nrow(DATA),floor(nrow(DATA)*0.5))
logical<-rep(FALSE,nrow(DATA))
logical[train_idx]<-TRUE
train_idx<-logical
test_idx<-!logical
rm(logical)
TRAIN<-DATA[train_idx,]
TEST<-DATA[test_idx,]
