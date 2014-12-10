####################################################
####Cumulative (proportional odds) logit analysis
####################################################
rm(list=ls())
gc()
graphics.off()
op<-par()
if(!require("MASS"))
  (install.packages("MASS"))
if(!require("ggplot2"))
  (install.packages("ggplot2"))


###Import Data
RAW_DATA<-read.csv("Raw_Data.csv")
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



###Try a basic Proportional Odds Model on the genotype Frequencies
fit<-polr(Productivity~GF+Race,data=TRAIN)
summary(fit)
probs<-fitted(fit)
display<-round(probs,2)
geno_values<-seq(min(DATA$GF),max(DATA$GF),length.out=10000)

####make charts for different races
par_disp<-unique(TRAIN$Race)
par(mfrow=c(3,2))
race=par_disp[1]


for(race in par_disp){
pi_preds<-as.data.frame(predict(fit,newdata=data.frame("GF"=geno_values,"Race"=race),
                                type="probs"))

#calculate cutoff values
Cutoff<-array(3)
for(i in 1:3){
  CUT<-array(1:10000)
  for(j in 1:10000){
    CUT[j]=pi_preds[j,i+1]>pi_preds[j,i]}
  intersect=10000-sum(CUT)
  Cutoff[i]<-(geno_values[intersect])
}
print(round(Cutoff,2))
t=signif(exp(Cutoff),3)

###########################################
#####Plot of the Test Data With Cutoffs
####Concordance of the TEST DATA
CUT_TEST<-TEST[TEST$Race==race, ]
class<-predict(fit,newdata=CUT_TEST)
concordance<-numeric(length(class))
for (i in 1:nrow(CUT_TEST)){
  concordance[i]=as.numeric(CUT_TEST$Productivity[i] %in% class[i])                                   
}
print(round(mean(concordance),3))

correct<-round(mean(concordance),3)*100

plot(pi_preds$D~geno_values,ylim=c(0,1),type="l",lwd=2,
     ylab="Membership Probability",xlab="log(Genotype Frequency)")
title(main=paste0("Test Data Model Fit: ",correct,"% concordant ",race),cex.main=.8)
lines(pi_preds$C~geno_values,lwd=2,lty=2)
lines(pi_preds$B~geno_values,lwd=2,lty=3)
lines(pi_preds$A~geno_values,lwd=2,lty=4)
legend("topleft",legend=c("Group D","Group C","Group B", "Group A"),lty=c(1:4),lwd=2,cex=.7)
abline(v=Cutoff[1])
abline(v=Cutoff[2])
abline(v=Cutoff[3])

y<-as.numeric(CUT_TEST$Productivity)/4
x<-CUT_TEST$GF
text(x=x,y=y,labels=as.character(CUT_TEST$Productivity))
}

par(op)###reset graphics



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

plot(H1~H2,data=DATA)###check data support
plot(GF~I(H1+H2),data=DATA)###check that GF=H1*H2
float<-lm(GF~I(H1+H2),data=DATA)
abline(float)
grid<-create_grid(DATA$H1,DATA$H2,n=10000)
grid$GF<-(grid$H1+grid$H2)*float$coefficients[2]+float$coefficients[1]


race=par_disp[1]
for(race in par_disp){
  
CUT_TEST<-TEST[TEST$Race==race, ]
grid$Race=race
grid$class<-predict(fit,newdata = grid)
class<-predict(fit,newdata=CUT_TEST)
concordance<-numeric(length(class))
for (i in 1:nrow(CUT_TEST)){
  concordance[i]=as.numeric(CUT_TEST$Productivity[i] %in% class[i])                                   
}
correct<-round(mean(concordance),3)*100
base_layer<-ggplot(data=grid,aes(x=H1,y=H2,colour=class))+
  ggtitle(paste("Decision Boundaries, Race=",race))+
  geom_point(size=3.5,alpha=1,shape=15)

train_plot<-base_layer+geom_text(data=CUT_TEST,aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle(paste("Decision Boundaries for",race,"Percent Correct=",correct,"%"))+facet_wrap(~Productivity,ncol=2)
print(train_plot)

}


















