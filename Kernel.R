####################################################
####More General Kernel Method Approaches
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


#####fit a ksvm kernel to output class probabilities (versus call)
######Fit a simple SVM using the e1071 package
ww<-table(TRAIN$Productivity)/nrow(TRAIN)
# ###this shows the effect of loss function
# weights["A"]=0.97
# weights["B"]=0.01
# weights["C"]=0.01
# weights["D"]=0.01

#ww=c("A"=0.25,"B"=0.25,"C"=0.25,"D"=0.25)
#TRAIN$Productivity<-as.character(TRAIN$Productivity)
cost=c(0.001,0.01,1,5,10,100)
error<-cost
for (i in 1:length(cost)){
mod<-ksvm(Productivity~H1+H2+Race,data=TRAIN,kernel="rbfdot",class.weights =NULL,
              kpar="automatic",cross=3,C=cost[i],prob.model=TRUE)
print(paste("Loop",i,",C=",cost[i]))
error[i]<-mod@cross
}

plot(error,type="b")
C=cost[which(error==min(error))]

bestmod<-ksvm(Productivity~H1+H2+Race,data=TRAIN,kernel="rbfdot", class.weights=NULL,
          kpar="automatic",cross=3,C=C,prob.model=TRUE)



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
  grid$class<-predict(bestmod,newdata = grid)
  func=predict(bestmod,newdata = grid,type ="probabilities" )
  func<-data.frame(grid,func)
  class<-predict(bestmod,newdata=CUT_TEST)
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
    
  plot_dat<-melt(func,id=c("H1","H2","Race","class"))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",race,"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)
  print(plot)  
}



##########################################################################
#####Now Apply a gaussian Process to the data using RBF kernel
#####This has a natural probabilisitc interpretation that SVM doesn't have
##########################################################################

bestmod<-gausspr(Productivity~H1+H2+Race,data=TRAIN,kernel="rbfdot",
              kpar="automatic")


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
  func=predict(bestmod,newdata = grid,type ="probabilities" )
  func<-data.frame(grid,func)
  class<-predict(bestmod,newdata=CUT_TEST)
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




