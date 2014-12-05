####################################################
####SVM and Kernel Method Approach
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
train_idx<-sample(1:nrow(DATA),floor(nrow(DATA)*0.8))
logical<-rep(FALSE,nrow(DATA))
logical[train_idx]<-TRUE
train_idx<-logical
test_idx<-!logical
rm(logical)
TRAIN<-DATA[train_idx,]
TEST<-DATA[test_idx,]


######Fit a simple SVM using the e1071 package
weights<-table(TRAIN$Productivity)/nrow(TRAIN)
# ###this shows the effect of loss function
# weights["A"]=0.97
# weights["B"]=0.01
# weights["C"]=0.01
# weights["D"]=0.01
svmfit=svm(Productivity~H1+H2,data=TRAIN,kernel="polynomial",
           cost=3,class.weights=weights,degree=1,scale=T)
plot(svmfit,TRAIN,H1~H2)

tune.out<-tune(svm,Productivity~H1+H2,data=TRAIN,kernel="polynomial",
         class.weights=weights,degree=1,scale=T,
     ranges=list(cost=c(0.001,0.01,1,5,10,100)))
summary(tune.out)

bestmod=tune.out$best.model
summary(bestmod)

plot(bestmod,TRAIN,H1~H2)

tune.out<-tune(svm,Productivity~H1+H2,data=TRAIN,kernel="radial",
               class.weights=weights,scale=T,
               ranges=list(cost=c(0.001,0.01,1,5,10,100)))
summary(tune.out)

bestmod=tune.out$best.model
summary(bestmod)

plot(bestmod,TRAIN,H1~H2)


#####Get training error
class<-predict(bestmod)
concordance<-numeric(length(class))
for (i in 1:nrow(TRAIN)){
  concordance[i]=as.numeric(TRAIN$Productivity[i] %in% class[i])                                   
}
print(round(mean(concordance),3))
correct<-round(mean(concordance),3)*100




#####Get test error
class<-predict(bestmod,newdata=TEST)
concordance<-numeric(length(class))
for (i in 1:nrow(TEST)){
  concordance[i]=as.numeric(TEST$Productivity[i] %in% class[i])                                   
}
print(round(mean(concordance),3))

correct<-round(mean(concordance),3)*100


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
grid$class<-predict(bestmod,newdata = grid)

base_layer<-ggplot(data=grid,aes(x=H1,y=H2,colour=class))+
  ggtitle("Decision Boundaries")+
  geom_point(size=3.5,alpha=1,shape=15)
print(base_layer)
train_plot<-base_layer+geom_text(data=TEST,aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle("Decision Boundaries and All Classes")
print(train_plot)


###look at each independently
type="D"
train_plot<-base_layer+geom_text(data=TEST[TEST$Productivity==type,],aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle(paste("Decision Boundaries and Class",type))
print(train_plot)
type="C"
train_plot<-base_layer+geom_text(data=TEST[TEST$Productivity==type,],aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle(paste("Decision Boundaries and Class",type))
print(train_plot)
type="B"
train_plot<-base_layer+geom_text(data=TEST[TEST$Productivity==type,],aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle(paste("Decision Boundaries and Class",type))
print(train_plot)
type="A"
train_plot<-base_layer+geom_text(data=TEST[TEST$Productivity==type,],aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle(paste("Decision Boundaries and Class",type))
print(train_plot)






