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
train_idx<-sample(1:nrow(DATA),floor(nrow(DATA)*0.5))
logical<-rep(FALSE,nrow(DATA))
logical[train_idx]<-TRUE
train_idx<-logical
test_idx<-!logical
rm(logical)
TRAIN<-DATA[train_idx,]
TEST<-DATA[test_idx,]



###Try a basic Proportional Odds Model on the genotype Frequencies
fit<-polr(Productivity~GF*Race,data=TRAIN)
summary(fit)
probs<-fitted(fit)
display<-round(probs,2)
geno_values<-seq(min(DATA$GF),max(DATA$GF),length.out=10000)


###################################
####TRAIN Data Fit
####make charts for different races
###################################
par_disp<-unique(TRAIN$Race)
#par(op)
#par(mfrow=c(3,2))
race=par_disp[1]


correct<-data.frame()
for(race in par_disp){ 
  CUT_TRAIN<-TRAIN[TRAIN$Race==race, ]
  class<-predict(fit,newdata=CUT_TRAIN)
  concordance<-numeric(length(class))
  for (i in 1:nrow(CUT_TRAIN)){
    concordance[i]=as.numeric(CUT_TRAIN$Productivity[i] %in% class[i])                                   
  }
  
  correct<-rbind(correct,data.frame(pct_correct=mean(concordance),Race=race))
  
}


####set up race grids for prediction values
grid_values<-data.frame()

for(race in par_disp){ 
  grid<-data.frame("GF"=geno_values,"Race"=race)  
  preds<-as.data.frame(predict(fit,newdata=grid,type="probs"))
  grid$Race=paste(race,paste0(round(correct$pct_correct[correct$Race==race],3)*100,"% ","Correct"))
  grid_values<-rbind(grid_values,cbind(grid,preds)) 
}

grid<-melt(grid_values,id=c("GF","Race"))

colnames(grid)[3]<-"Prognosis"


####calculate cutoff values from grid
cutoff_values<-data.frame()
par_disp=unique(grid_values$Race)
for(race in par_disp){ 
  dat<-subset(grid_values,Race==race)
  dat<-dat[c("C","B","A")]
  Cutoff<-data.frame(ncol(dat)-1)
  for(i in 1:2){
    CUT<-array(1:10000)
    for(j in 1:10000){
      CUT[j]=dat[j,i+1]>dat[j,i]}
    intersect=10000-sum(CUT)
    Cutoff[i]<-(geno_values[intersect])
  }
  
  cutoff_values<-rbind(cutoff_values,data.frame(cut=t(Cutoff),Race=race))   
}



class(grid$Race)

TRAIN$Race<-as.character(TRAIN$Race)
mapping<-unique(grid$Race)
for ( b in unique(TRAIN$Race)){
  TRAIN$Race[TRAIN$Race==b]=mapping[grep(b,mapping)]
  
}

TRAIN$value<-0
TRAIN$value[ TRAIN$Productivity=="C"]<-1/6
TRAIN$value[ TRAIN$Productivity=="B"]<-1/6+1/3
TRAIN$value[ TRAIN$Productivity=="A"]<-1/6+2/3



plot_boundaries_train<-ggplot(data=grid)+geom_line(aes(x=GF,y=value,colour=Prognosis,lty=Prognosis),lwd=1.25)+
  facet_wrap(~Race)+geom_vline(data=cutoff_values,aes(xintercept=cut),lwd=1,lty=1)+
  geom_text(data=TRAIN,aes(x=GF,y=value,label=Productivity),size=5,colour="black")+
  ylab("Prognosis Probability")+xlab("Log(Genotype Frequency)")+
  ggtitle(paste0("Training Data: n=",nrow(TRAIN)," Observations Used for Model Fit"))

print(plot_boundaries_train)

print(cutoff_values)

###################################
####TEST Data Fit
####make charts for different races
###################################
par_disp<-unique(TEST$Race)
#par(op)
#par(mfrow=c(3,2))
race=par_disp[1]


correct<-data.frame()
for(race in par_disp){ 
  CUT_TEST<-TEST[TEST$Race==race, ]
  class<-predict(fit,newdata=CUT_TEST)
  concordance<-numeric(length(class))
  for (i in 1:nrow(CUT_TEST)){
    concordance[i]=as.numeric(CUT_TEST$Productivity[i] %in% class[i])                                   
  }
  
  correct<-rbind(correct,data.frame(pct_correct=mean(concordance),Race=race))
  
}


####set up race grids for prediction values
grid_values<-data.frame()

for(race in par_disp){ 
  grid<-data.frame("GF"=geno_values,"Race"=race)  
  preds<-as.data.frame(predict(fit,newdata=grid,type="probs"))
  grid$Race=paste(race,paste0(round(correct$pct_correct[correct$Race==race],3)*100,"% ","Correct"))
  grid_values<-rbind(grid_values,cbind(grid,preds)) 
}

grid<-melt(grid_values,id=c("GF","Race"))

colnames(grid)[3]<-"Prognosis"


####calculate cutoff values from grid
cutoff_values<-data.frame()
par_disp=unique(grid_values$Race)
for(race in par_disp){ 
  dat<-subset(grid_values,Race==race)
  dat<-dat[c("C","B","A")]
  Cutoff<-data.frame(ncol(dat)-1)
  for(i in 1:2){
    CUT<-array(1:10000)
    for(j in 1:10000){
      CUT[j]=dat[j,i+1]>dat[j,i]}
    intersect=10000-sum(CUT)
    Cutoff[i]<-(geno_values[intersect])
  }
  
  cutoff_values<-rbind(cutoff_values,data.frame(cut=t(Cutoff),Race=race))   
}



class(grid$Race)

TEST$Race<-as.character(TEST$Race)
mapping<-unique(grid$Race)
for ( b in unique(TEST$Race)){
  TEST$Race[TEST$Race==b]=mapping[grep(b,mapping)]
  
}

TEST$value<-0
TEST$value[ TEST$Productivity=="C"]<-1/6
TEST$value[ TEST$Productivity=="B"]<-1/6+1/3
TEST$value[ TEST$Productivity=="A"]<-1/6+2/3




plot_boundaries_test<-ggplot(data=grid)+geom_line(aes(x=GF,y=value,colour=Prognosis,lty=Prognosis),lwd=1.25)+
  facet_wrap(~Race)+geom_vline(data=cutoff_values,aes(xintercept=cut),lwd=1,lty=1)+
  geom_text(data=TEST,aes(x=GF,y=value,label=Productivity),size=5,colour="black")+
  ylab("Prognosis Probability")+xlab("Log(Genotype Frequency)")+
  ggtitle(paste0("testing Data: n=",nrow(TEST)," Observations Held Out from Model Fit"))

print(plot_boundaries_test)

print(cutoff_values)


#########################################################
#########################################################
####Look at 2D Support
#########################################################
#########################################################

###Import Data
RAW_DATA<-read.csv("Raw_Data.csv")
RAW_DATA$Productivity.Group<-as.character(RAW_DATA$Productivity.Group)
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="B"]="A"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="C"]="B"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="D"]="C"


min_thresh_gf<-min(RAW_DATA$Total.Genotype.Frequency[ RAW_DATA$Total.Genotype.Frequency!=0],na.rm=T)
RAW_DATA$Total.Genotype.Frequency[ RAW_DATA$Total.Genotype.Frequency==0]<-min_thresh_gf
RAW_DATA$Total.Genotype.Frequency[ is.na(RAW_DATA$Total.Genotype.Frequency)]=min_thresh_gf

min_thresh_h1<-min(RAW_DATA$Frequency.1[ RAW_DATA$Frequency.1!=0],na.rm=T)
RAW_DATA$Frequency.1[ RAW_DATA$Frequency.1==0]<-min_thresh_h1
RAW_DATA$Frequency.1[ is.na(RAW_DATA$Frequency.1)]=min_thresh_h1

min_thresh_h2<-min(RAW_DATA$Frequency.2[ RAW_DATA$Frequency.2!=0],na.rm=T)
RAW_DATA$Frequency.2[ RAW_DATA$Frequency.2==0]<-min_thresh_h2
RAW_DATA$Frequency.2[ is.na(RAW_DATA$Frequency.2)]=min_thresh_h2


permute<-rbinom(nrow(RAW_DATA),1,prob=0.5)==1
RAW_DATA[permute,]<-RAW_DATA[permute,c(1:6,8,7,10,9,11:13)]
DATA<-RAW_DATA[c("RID","Race","Frequency.1",
                 "Frequency.2","Rank.1","Rank.2",
                 "Total.Genotype.Frequency","Productivity.Group")]

DATA$Productivity.Group<-ordered(DATA$Productivity.Group,levels=c("C","B","A"))
colnames(DATA)<-c("RID","Race","H1","H2","Rank_H1","Rank_H2","GF","Productivity")
DATA<-DATA[complete.cases(DATA),]###remove NA's
use_in_fit<-DATA$GF>min_thresh_gf & DATA$H1>min_thresh_h1 & DATA$H2>min_thresh_h2
rm(RAW_DATA)

###Format Predictors to be ln-scale
DATA$GF<-log(DATA$GF)
DATA$H1<-log(DATA$H1)
DATA$H2<-log(DATA$H2)
DATA[c("H1","H2","GF")]<-normalizeData(DATA[c("H1","H2","GF")],type="0_1")
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

plot(H1~H2,data=DATA[use_in_fit,])###check data support
plot(GF~I(H1+H2),data=DATA[use_in_fit,])###check that GF=H1*H2
float<-lm(GF~I(H1+H2),data=DATA[use_in_fit,])
abline(float)
grid<-create_grid(DATA$H1,DATA$H2,n=10000)
grid$GF<-(grid$H1+grid$H2)*float$coefficients[2]+float$coefficients[1]

par_disp<-unique(TEST$Race)
race=par_disp[1]
for(race in par_disp){
CUT_TEST<-TEST[TEST$Race==race, ]
grid$Race=CUT_TEST$Race[1]
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


















