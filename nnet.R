####################################################
####neural network approach
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
RAW_DATA$Productivity.Group<-as.character(RAW_DATA$Productivity.Group)
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="B"]="A"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="C"]="B"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="D"]="C"


permute<-rbinom(nrow(RAW_DATA),1,prob=0.5)==1
RAW_DATA[permute,]<-RAW_DATA[permute,c(1:6,8,7,10,9,11:13)]
DATA<-RAW_DATA[c("RID","Race","Frequency.1",
                 "Frequency.2","Rank.1","Rank.2",
                 "Total.Genotype.Frequency","Productivity.Group")]

DATA$Productivity.Group<-ordered(DATA$Productivity.Group,levels=c("C","B","A"))
colnames(DATA)<-c("RID","Race","H1","H2","Rank_H1","Rank_H2","GF","Productivity")
DATA<-DATA[complete.cases(DATA),]###remove NA's
DATA<-DATA[DATA$GF!=0 & DATA$H1!=0 & DATA$H2!=0,]
rm(RAW_DATA)

###Format Predictors to be ln-scale
DATA$GF<-log(DATA$GF)
DATA$H1<-log(DATA$H1)
DATA$H2<-log(DATA$H2)
DATA[c("H1","H2","GF")]<-normalizeData(DATA[c("H1","H2","GF")],type="0_1")
summary(DATA)

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
y<-TRAIN[c("Productivity")]
y<-decodeClassLabels(y[,1])
x_tst<-TEST[c("H1","H2")]
x_tst<-cbind(x_tst,decodeClassLabels(TEST[,c("Race")]))
y_tst<-TEST[c("Productivity")]
y_tst<-decodeClassLabels(y_tst[,1])

tmp_data<-splitForTrainingAndTest(x,y,ratio=0.3)


parameterGrid <-  expand.grid(c(0,3),c(3,5,9,15), c(0.00316, 0.0147,0.1))
colnames(parameterGrid) <-  c("nHidden_l1","nHidden_l2", "learnRate") 
rownames(parameterGrid) <- paste("nnet-", apply(parameterGrid, 1, function(x) {paste(x,sep="", collapse="-")}), sep="") 
models<-apply(parameterGrid, 1, function(p) { 
  s<-p[c(1,2)]
  if(s[1]==0){s=s[2]}
  l<-p[3]
  mlp(tmp_data$inputsTrain, tmp_data$targetsTrain, size=s, learnFunc="Std_Backpropagation",
      learnFuncParams=c(l, 0.1), maxit=200, inputsTest=tmp_data$inputsTest,
      targetsTest=tmp_data$targetsTest)
  })


pp<-ceiling(sqrt(nrow(parameterGrid)))
par(mfrow=c(pp,pp)) 
for(modInd in 1:length(models)) { 
      plotIterativeError(models[[modInd]], main=names(models)[modInd]) 
    }

par(op)

###measure error using multinomial deviance as loss function
trainErrors <-  data.frame(lapply(models, function(mod) {
  correct<-apply(tmp_data$targetsTrain,1,function(b){
    which(b==1)
  })
  correct<-as.matrix(data.frame(1:length(correct),correct))
  error<-sqrt(sum(-log(mod$fitted.values[correct]))) 
            error 
          })) 


testErrors <-  data.frame(lapply(models, function(mod) { 
        pred <-  predict(mod,tmp_data$inputsTest) 
        correct<-apply(tmp_data$targetsTest,1,function(b){
          which(b==1)
        })
        error<-sqrt(sum(-log(pred[correct]))) 
        error 
      })) 
t(trainErrors)
t(testErrors)


trainErrors[which(min(trainErrors) == trainErrors)]
testErrors[which(min(testErrors) == testErrors)]
param_chosen<-which(min(testErrors) == testErrors)###you can change this to pick other models
#param_chosen<-11
p<-as.numeric(parameterGrid[param_chosen,])
s<-p[c(1,2)]
if(s[1]==0){s=s[2]}
l<-p[3]
bestmod<-mlp(x,y, size=s, learnFunc="Std_Backpropagation",
    learnFuncParams=c(l, 0.05), maxit=200)


####now that we have fixed a model we can cheat and look at the test error
bestmod$snnsObject$getUnitDefinitions()
weightMatrix(bestmod)
extractNetInfo(bestmod)
plotIterativeError(bestmod)
predictions <- predict(bestmod,x_tst)
plotRegressionError(predictions[,4], y_tst[,4])
confusionMatrix(y_tst,predictions)
plotROC(predictions[,4],y_tst[,4])



#import plot function from Github
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')

plot.nnet(bestmod)
#plot.nnet(bestmod,nid=F)###don't show thickness
plot.nnet(bestmod,wts.only=T)

#import sensitivity analysis function
source_url('https://gist.githubusercontent.com/fawda123/6860630/raw/b8bf4a6c88d6b392b1bfa6ef24759ae98f31877c/lek_fun.r')

#sensitivity analsyis, note 'exp.in' argument
lek.fun(bestmod,exp.in=x,var.sens=c("X1","X2"),
        split.vals=c(0,.50,1))



###Plot the linear fit
######Plot 2-D grid of Classification Rule
create_grid<-function(x1,x2,n=1000){
  min_x1<-min(x1)
  max_x1<-max(x1)
  x1_seq<-seq(min_x1,max_x1,length.out=floor(sqrt(n)))
  min_x2<-min(x2)
  max_x2<-max(x2)
  x2_seq<-seq(min_x2,max_x2,length.out=floor(sqrt(n)))
  grid<-expand.grid(H1=x1_seq,H2=x2_seq)
  return(grid)
}

s_grid<-create_grid(DATA$H1,DATA$H2,n=10000)

par_disp<-unique(TRAIN$Race)
race=par_disp[1]

###look at the response surface after applying the softmax function
  
  for(race in par_disp){
    idx<-as.logical(x_tst[,colnames(x_tst)==race])
    CUT_TEST<-x_tst[idx, ]
    grid=data.frame(s_grid,CUT_TEST[1,colnames(CUT_TEST) %in% par_disp,drop=F])
    func=predict(bestmod,grid)
    colnames(func)<-colnames(y_tst)
    func<-data.frame(grid,func)
    tmp<-predict(bestmod,CUT_TEST)
    colnames(tmp)<-colnames(y_tst)
    class<-apply(tmp,1,function(b){
      which(b==max(b))
    })
    class<-colnames(tmp)[class]
    
    response<-apply(y_tst[idx,],1,function(b){
      which(b==max(b))
    })
    response<-colnames(y_tst)[response]
    concordance<-numeric(length(response))
    
    for (i in 1:nrow(CUT_TEST)){
      concordance[i]=as.numeric(response[i] %in% class[i])                                   
    }
    correct<-round(mean(concordance),3)*100
    
    
    plot_dat<-melt(func,id=c("H1","H2","AFA","API","CAU","HIS","UNK"))
    plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
      geom_tile()+
      stat_contour(data=plot_dat,aes(z=value))+
      ggtitle(paste("Probability Distribution for",race,"Percent Correct=",correct,"%"))+
      facet_wrap(~variable,ncol=2)
    print(plot)  
  }
  
  







