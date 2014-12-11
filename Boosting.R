####################################################
####More General Kernel Method Approaches
####################################################
rm(list=ls())
gc()
graphics.off()
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


###fit a boosted model
gbm1 <-
  gbm(Productivity~H1+H2+Race,         # formula
      data=TRAIN,                   # dataset
      var.monotone=c(0,0,0), # -1: monotone decrease,
      # +1: monotone increase,
      #  0: no monotone restrictions
      distribution="multinomial",     # see the help for other choices
      n.trees=1000,                # number of trees
      shrinkage=0.05,              # shrinkage or learning rate,
      # 0.001 to 0.1 usually work
      interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
      bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
      train.fraction = 0.5,        # fraction of data for training,
      # first train.fraction*N used for training
      n.minobsinnode = 10,         # minimum total weight needed in each node
      cv.folds = 3,                # do 3-fold cross-validation
      keep.data=TRUE,              # keep a copy of the dataset with the object
      verbose=FALSE,               # don't print out progress
      n.cores=1)                   # use only a single core (detecting #cores is
# # error-prone, so avoided here)
# # check performance using an out-of-bag estimator
# # OOB underestimates the optimal number of iterations
# best.iter <- gbm.perf(gbm1,method="OOB")
# print(best.iter)
# 
# # check performance using a 50% heldout test set
# best.iter <- gbm.perf(gbm1,method="test")
# print(best.iter)

# check performance using 5-fold cross-validation
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)

# plot the performance # plot variable influence
summary(gbm1,n.trees=1)         # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees

# compactly print the first and last trees for curiosity
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot(gbm1,1,best.iter)
plot(gbm1,2,best.iter)
plot(gbm1,3,best.iter)
par(mfrow=c(1,1))
# contour plot of variables 1 and 2 after "best" iterations
#plot(gbm1,1:2,best.iter)
# lattice plot of variables 2 and 3
###3 way plot
#plot(gbm1,1:3,best.iter)
# do another 100 iterations

#best.iter=1000
bestmod<-gbm1


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

grid<-create_grid(DATA$H1,DATA$H2,n=10000)


par_disp<-unique(TRAIN$Race)
race=par_disp[1]


for(race in par_disp){
  CUT_TEST<-TEST[TEST$Race==race, ]
  grid$Race=CUT_TEST$Race[1]
  func=predict(bestmod,newdata = grid,best.iter,type ="response" )
  func<-data.frame(grid,func)
  tmp<-predict(bestmod,newdata=CUT_TEST,best.iter,type="response")
  class<-apply(tmp,1,function(b){
    which(b==max(b))
  })
  class<-colnames(tmp)[class]
  concordance<-numeric(length(class))
  for (i in 1:nrow(CUT_TEST)){
    concordance[i]=as.numeric(CUT_TEST$Productivity[i] %in% class[i])                                   
  }
  correct<-round(mean(concordance),3)*100
  
  #   base_layer<-ggplot(data=grid,aes(x=H1,y=H2,colour=class))+
  #     ggtitle(paste("Decision Boundaries, Race=",race))+
  #     geom_point(size=3.5,alpha=1,shape=15)
  #   
  #   train_plot<-base_layer+geom_text(data=CUT_TEST,aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  #     ggtitle(paste("Decision Boundaries for",race,"Percent Correct=",correct,"%"))+facet_wrap(~Productivity,ncol=2)
  #   print(train_plot)
  
  plot_dat<-melt(func,id=c("H1","H2","Race"))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",race,"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)
  print(plot)  
}

####if more iterations are desired
gbm1 <- gbm.more(gbm1,1000,
                 verbose=FALSE) # stop printing detailed progress



