####################################################
####Cumulative (proportional odds) logit analysis ##
####Mark Albrecht, edited by Rachel Fonstad
####11/14/2013, 12/06/2013
####################################################
rm(list=ls())
gc()
graphics.off()
op<-par(mfrow=c(2,2))


library(MASS)
setwd("K:/xover/cibmtr/Rachel/Search Strategy")


##change dataset here
DATA<-read.csv("Data for Bioinformatics all UNK are CAU.csv")


summary(DATA)

response2<-DATA$Productivity.Group

response<-ifelse(response2=='A','A',ifelse(response2=='B','A',ifelse(response2=='C','C','D')))



response<-ordered(response,levels=c("D","C","A"))
Geno_freq<-DATA$Total.Genotype.Frequency
Geno_freq[which(Geno_freq==0)]=NA
summary(Geno_freq)
pull<-which(is.na(Geno_freq))
Geno_freq<-Geno_freq[-pull]
response<-response[-pull]

summary(Geno_freq)
Geno_freq<-log(Geno_freq)
summary(Geno_freq)


numvec<-as.numeric(response)


fit<-polr(response~Geno_freq)
summary(fit)
probs<-fitted(fit)
display<-round(probs,2)
print(display[1:30,])
geno_values<-seq(min(Geno_freq),max(Geno_freq),length.out=10000)
pi_preds<-as.data.frame(predict(fit,newdata=data.frame("Geno_freq"=geno_values),
                  type="probs"))
class<-predict(fit)
concordance<-numeric(length(class))
i=1
for (i in 1:length(response)){
  concordance[i]=as.numeric(response[i] %in% class[i])                                   
}
print(round(mean(concordance),3))



#calculate cutoff values
Cutoff<-array(2)
for(i in 1:2){
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


plot(numvec~Geno_freq,type="n")
text(Geno_freq,numvec,labels=as.character(response))
abline(v=Cutoff[1])
abline(v=Cutoff[2])



plot(pi_preds$D~geno_values,xlim=c(min(geno_values),max(geno_values)),ylim=c(0,1),type="l",lwd=2,
     ylab="Membership Probability",xlab="log(Genotype Frequency)")
title(main="Probability of Search Productivity Group by Genotype Freq",cex.main=.8)
lines(pi_preds$C~geno_values,lwd=2,lty=3)
lines(pi_preds$A~geno_values,lwd=2,lty=4)
legend("topleft",legend=c("Group D","Group C", "Group A/B"),lty=c(1:4),lwd=2,cex=.7)
abline(v=Cutoff[1])
abline(v=Cutoff[2])


val<-exp(Cutoff[3])+.1*exp(Cutoff[3])
plot(pi_preds$D~exp(geno_values),xlim=c(0,max(exp(geno_values))/10),ylim=c(0,1),type="l",lwd=2,
     main="Probability of Search Productivity Group by Genotype Freq",cex.main=.8,ylab="Membership Probability",
     xlab="Genotype Frequency")
lines(pi_preds$C~exp(geno_values),lwd=2,lty=3)
lines(pi_preds$A~exp(geno_values),lwd=2,lty=4)
legend("topleft",legend=c("Group D","Group C", "Group A/B"),lty=c(1:4),lwd=2,cex=.7)
abline(v=exp(Cutoff[1]))
abline(v=exp(Cutoff[2]))


response<-ordered(response,levels=c("A","C","D"))
class<-ordered(class,levels=c("A","C","D"))
counts<-table(response,class)
tvals<-t(prop.table(counts,1))*100
round(tvals,1)
barplot(tvals,main="GF precision",xlab="Productivity group",legend=rownames(tvals))

table(response)

Cutoff[1]
Cutoff[2]
exp(Cutoff[1])
exp(Cutoff[2])

