####################################################
####Cumulative (proportional odds) logit analysis
####################################################
source("data_prep.R")
####summarize deomgraphics data
table(TRAIN$Race)
table(TEST$Race)

ftable(TRAIN$Race,TRAIN$Productivity)
ftable(TEST$Race,TEST$Productivity)


###Try a basic Proportional Odds Model on the genotype Frequencies
###first fit weights to the different classes
ww<-table(TRAIN$Productivity,TRAIN$Race)
ww<-as.data.frame(ww)
for(race in unique(ww$Var2)){
  idx<-ww$Var2==race
  ww$Freq[idx]<-ww$Freq[idx]/sum(ww$Freq[idx])
  ww$Freq[idx]<-1/ww$Freq[idx]
}
TRAIN$weights_fit<-0
for(i in 1:nrow(ww)){
  idx<-TRAIN$Race==ww$Var2[i] & as.character(TRAIN$Productivity)==ww$Var1[i]
  TRAIN$weights_fit[idx]<-ww$Freq[i]
}

fit<-polr(Productivity~GF*Race,weights=weights_fit,data=TRAIN)
#fit<-polr(Productivity~GF*Race,data=TRAIN)
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
  class_pred<-predict(fit,newdata=CUT_TRAIN)
  out<-data.frame(table(class_pred,CUT_TRAIN$Productivity))
  for(id in unique(out$Var2)){
    out$Freq[out$Var2==id]<-out$Freq[out$Var2==id]/sum(out$Freq[out$Var2==id])
  }
  
  out<-out[ out$class_pred==out$Var2,]
  c<-paste0(out$Var2,"=",round(out$Freq,2)*100,"%",collapse="; ")
  correct<-rbind(correct,data.frame(pct_correct=c,Race=race))
  
}


####set up race grids for prediction values
grid_values<-data.frame()

for(race in par_disp){ 
  grid<-data.frame("GF"=geno_values,"Race"=race)  
  preds<-as.data.frame(predict(fit,newdata=grid,type="probs"))
  grid$Race=paste0(race,": (",correct$pct_correct[correct$Race==race],")")
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

cutoff_natural_scale<-cutoff_values
cutoff_natural_scale$cut<-exp(cutoff_natural_scale$cut)
cutoff_natural_scale$cut<-signif(cutoff_natural_scale$cut,2)
print(cutoff_natural_scale)


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
  theme(strip.text.x = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))

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


####get the confusion matrix for all combinations
all<-data.frame(prediction=predict(fit,newdata=TEST))
all$prediction<-all$prediction
all$truth<-TEST$Productivity
all$race<-TEST$Race

result<-ftable(all$race,all$truth,all$prediction)
print(result)
result<-data.frame(result)
colnames(result)<-c("Race","Truth","Prediction","Freq")

for(u in unique(result$Race)){
  for(b in unique(result$Truth)){
    idx<-result$Race==u & result$Truth==b
    result$Freq[idx]<-result$Freq[idx]/sum(result$Freq[idx])
  }
}

plot_dat<-arrange(result,Race,Truth,Prediction)
plot_dat<-ddply(plot_dat,c("Race","Truth"),transform,label_y=cumsum(Freq)-0.5*Freq)


plot_bars<-ggplot(data=plot_dat,aes(x=Truth,y=Freq,fill=Prediction))+
  geom_bar(stat="identity")+facet_wrap(~Race)+
  geom_text(aes(y=label_y,label=round(Freq*100,1)))+
  ggtitle("Percent Concordance on Validation Data by Class")+
  xlab("True Search Productivity Classification")+ylab("Percent")
print(plot_bars)



correct<-data.frame()
for(race in par_disp){ 
  CUT_TEST<-TEST[TEST$Race==race, ]
  class_pred<-predict(fit,newdata=CUT_TEST)
  out<-data.frame(table(class_pred,CUT_TEST$Productivity))
  for(id in unique(out$Var2)){
    out$Freq[out$Var2==id]<-out$Freq[out$Var2==id]/sum(out$Freq[out$Var2==id])
  }
  
  out<-out[ out$class_pred==out$Var2,]
  c<-paste0(out$Var2,"=",round(out$Freq,2)*100,"%",collapse="; ")
  correct<-rbind(correct,data.frame(pct_correct=c,Race=race))
}


####set up race grids for prediction values
grid_values<-data.frame()

for(race in par_disp){ 
  grid<-data.frame("GF"=geno_values,"Race"=race)  
  preds<-as.data.frame(predict(fit,newdata=grid,type="probs"))
  grid$Race=paste0(race,": (",correct$pct_correct[correct$Race==race],")")
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
  theme(strip.text.x = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))

print(plot_boundaries_test)

print(cutoff_values)

#########################################################
#########################################################
####Look at 2D Support
#########################################################
#########################################################

source("data_prep.R")


###Try a basic Proportional Odds Model on the genotype Frequencies
###first fit weights to the different classes
ww<-table(TRAIN$Productivity,TRAIN$Race)
ww<-as.data.frame(ww)
for(race in unique(ww$Var2)){
  idx<-ww$Var2==race
  ww$Freq[idx]<-ww$Freq[idx]/sum(ww$Freq[idx])
  ww$Freq[idx]<-1/ww$Freq[idx]
}
TRAIN$weights_fit<-0
for(i in 1:nrow(ww)){
  idx<-TRAIN$Race==ww$Var2[i] & as.character(TRAIN$Productivity)==ww$Var1[i]
  TRAIN$weights_fit[idx]<-ww$Freq[i]
}

fit<-polr(Productivity~GF*Race,weights=weights_fit,data=TRAIN)



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
out<-data.frame(table(class,CUT_TEST$Productivity))
  for(id in unique(out$Var2)){
    out$Freq[out$Var2==id]<-out$Freq[out$Var2==id]/sum(out$Freq[out$Var2==id])
  }
out<-out[ out$class==out$Var2,]
c<-paste0(out$Var2,"=",round(out$Freq,2)*100,"%",collapse="; ")

base_layer<-ggplot(data=grid,aes(x=H1,y=H2,colour=class))+
  ggtitle(paste("Decision Boundaries, Race=",race))+
  geom_point(size=3.5,alpha=1,shape=15)

train_plot<-base_layer+geom_text(data=CUT_TEST,aes(x=H1,y=H2,label=Productivity),size=5,colour="black")+
  ggtitle(paste("Decision Boundaries for",race,"Percent Correct=",c))+facet_wrap(~Productivity,ncol=2)
print(train_plot)

}


















