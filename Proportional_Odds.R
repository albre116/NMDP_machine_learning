####################################################
####Cumulative (proportional odds) logit analysis
####################################################
rm(list=ls())
gc()
graphics.off()
library(MASS)

###Import Data
RAW_DATA<-read.csv("Raw_Data.csv")
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


###Try a basic Proportional Odds Model on the genotype Frequencies
fit<-polr(Productivity~GF,data=DATA)
summary(fit)
probs<-fitted(fit)
display<-round(probs,2)
print(display[1:30,])
geno_values<-seq(min(Geno_freq),max(Geno_freq),length.out=10000)
pi_preds<-as.data.frame(predict(fit,newdata=data.frame("Geno_freq"=geno_values),
                                type="probs"))
class<-predict(fit)
concordance<-numeric(length(class))
for (i in 1:nrow(DATA)){
  concordance[i]=as.numeric(DATA$Productivity.Group[i] %in% class[i])                                   
}
print(round(mean(concordance),3))



#calculate cutoff values
Cutoff<-array(3)
for(i in 1:3){
  CUT<-array(1:10000)
  for(j in 1:10000){
    CUT[j]=pi_preds[j,i+1]>pi_preds[j,i]}
  intersect=10000-sum(CUT)
  Cutoff[i]<-(geno_values[intersect])
}
print(Geno_freq)
print(round(Cutoff,2))
t=signif(exp(Cutoff),3)
print(t)




plot(pi_preds$D~geno_values,xlim=c(-30,15),ylim=c(0,1),type="l",lwd=2,
     ylab="Membership Probability",xlab="log(Genotype Frequency)")
title(main="Probability of Search Productivity Group by Genotype Freq",cex.main=.8)
lines(pi_preds$C~geno_values,lwd=2,lty=2)
lines(pi_preds$B~geno_values,lwd=2,lty=3)
lines(pi_preds$A~geno_values,lwd=2,lty=4)
legend("topleft",legend=c("Group D","Group C","Group B", "Group A"),lty=c(1:4),lwd=2,cex=.7)
legend("center",legend=c(as.expression(paste("Race: ",races[m], sep="")),as.expression(paste("Worksheet #",wksheet,sep=""))),bty="n")
abline(v=Cutoff[1])
abline(v=Cutoff[2])
abline(v=Cutoff[3])

val<-exp(Cutoff[3])+.1*exp(Cutoff[3])
plot(pi_preds$D~exp(geno_values),xlim=c(0,val),ylim=c(0,1),type="l",lwd=2,
     main="Probability of Search Productivity Group by Genotype Freq",cex.main=.8,ylab="Membership Probability",
     xlab="Genotype Frequency")
lines(pi_preds$C~exp(geno_values),lwd=2,lty=2,)
lines(pi_preds$B~exp(geno_values),lwd=2,lty=3)
lines(pi_preds$A~exp(geno_values),lwd=2,lty=4)
legend("topleft",legend=c("Group D","Group C","Group B", "Group A"),lty=c(1:4),lwd=2,cex=.7)
legend("center",legend=c(as.expression(paste("Race: ",races[m], sep="")),as.expression(paste("Worksheet #",wksheet,sep=""))),bty="n")
abline(v=exp(Cutoff[1]))
abline(v=exp(Cutoff[2]))
abline(v=exp(Cutoff[3]))


plot.new()
legend("center",legend=c(paste("D to C: ",round(Cutoff[1],2),t[1],sep=" "),
                         paste("C to B: ",round(Cutoff[2],2),t[2],sep="  "),paste("B to A: ",round(Cutoff[3],2),t[3],sep="  "),
                         paste("Concordance: ",round(mean(concordance),3),sep="")))
