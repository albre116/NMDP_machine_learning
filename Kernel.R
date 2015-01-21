####################################################
####More General Kernel Method Approaches
####################################################
source("data_prep.R")

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
race=par_disp[3]


for(race in par_disp){
  CUT_TRAIN<-TRAIN[TRAIN$Race==race, ]
  weights<-100/table(CUT_TRAIN$Productivity)/nrow(CUT_TRAIN)
  
  ###now, in this example the class.weights function doesn't work properly
  ###what is an equivalent solution from the "Loss" function perspective???
  ###simply sample data in equal proportions to create balanced data
  nsample<-table(CUT_TRAIN$Productivity)
  CUT_ADJUST<-data.frame()
  for(i in names(nsample)){
    idx<-rep_len(1:nsample[i],max(nsample))
    TMP<-CUT_TRAIN[CUT_TRAIN$Productivity==i,]
    CUT_ADJUST<-rbind(CUT_ADJUST,TMP[idx,])
  }
    
  table(CUT_ADJUST$Productivity)
  
  
  cost=c(0.001,0.01,1,5,10,100)
  error<-cost
  i=1
  for (i in 1:length(cost)){
    mod<-ksvm(Productivity~H1+H2,data=CUT_ADJUST,kernel="rbfdot", class.weights =NULL,
              kpar="automatic",cross=3,C=cost[i],prob.model=TRUE)
    print(paste("Loop",i,",C=",cost[i]))
    error[i]<-mod@cross
  }
  
  plot(error,type="b")
  C=cost[which(error==min(error))]
  
  bestmod<-ksvm(Productivity~H1+H2+Race,data=CUT_ADJUST,kernel="rbfdot", class.weights=NULL,
                kpar="automatic",cross=3,C=C,prob.model=TRUE)
  
    
  ###do the holdout validation
  CUT_TEST<-TEST[TEST$Race==race, ]
  grid$Race=CUT_TEST$Race[1]
  grid$class<-predict(bestmod,newdata = grid)
  class<-predict(bestmod,newdata=CUT_TEST)
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
  
  ###check out probabilities instead of class calls

  func=predict(bestmod,newdata = grid,type ="probabilities" )
  func<-data.frame(grid,func)
  class<-predict(bestmod,newdata=CUT_TEST)

  
  plot_dat<-melt(func,id=c("H1","H2","Race","class"))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",race,"Percent Correct=",c))+
    facet_wrap(~variable,ncol=2)
  print(plot)  
  
}




##########################################################################
###Now apply a sparse vector machine (l1 & l2 regularized) 
###that should perform similar to
###the kernel based methods above
###Use a penalized GLM fit to the data to pull out relevant kernel elements
###Plus model everything together
##########################################################################
ww<-table(TRAIN$Productivity,TRAIN$Race)
ww<-as.data.frame(ww)
for(race in unique(ww$Var2)){
  idx<-ww$Var2==race
  ###get the count for this race
  
  ###adjust within category proportions
  ww$Freq[idx]<-ww$Freq[idx]/sum(ww$Freq[idx])
  ww$Freq[idx]<-1/ww$Freq[idx]

}

TRAIN$weights_fit<-1
for(i in 1:nrow(ww)){
  idx<-TRAIN$Race==ww$Var2[i] & as.character(TRAIN$Productivity)==ww$Var1[i]
  TRAIN$weights_fit[idx]<-ww$Freq[i]
}

par_disp<-unique(TRAIN$Race)
race=par_disp[2]

####only do API otherwise takes forever
for(race in race){

###take subsample or algorithm will take forever
TRAIN_SUB<-TRAIN[ TRAIN$Race==race,]
TEST_SUB<-TEST[ TEST$Race==race,]

y=as.factor(as.character(TRAIN_SUB$Productivity))
y_tst=as.factor(as.character(TEST_SUB$Productivity))
x=as.matrix(TRAIN_SUB[c("H1","H2")])
#x<-cbind(x,decodeClassLabels(TRAIN_SUB[,c("Race")]))
x_tst=as.matrix(TEST_SUB[c("H1","H2")])
#x_tst<-cbind(x_tst,decodeClassLabels(TEST_SUB[,c("Race")]))
weights<-TRAIN_SUB$weights_fit


###fit kernel for the x matrix
models<-list()
param<-c(0.25,1,1.5)
i=1
for(i in 1:length(param)){
  print(paste("on loop",i))
  rbf<-rbfdot(sigma=param[i])
  kx<-kernelMatrix(rbf,x)
  kx_tst<-kernelMatrix(rbf,x_tst,x)
  cvfit=cv.glmnet(kx,y,family="multinomial",weights=weights)
  models[[i]]<-cvfit
}

win<-ceiling(sqrt(length(param)))
par(mfrow=c(win,win))
lapply(models,function(b){plot(b)})
cv_result<-lapply(models,function(b){min(b$cvm)})
cv_result<-data.frame("lambda"=param,"Min Deviance"=unlist(cv_result))
print(cv_result)
cvfit<-models[[which(cv_result$Min.Deviance==min(cv_result$Min.Deviance))]]
par(op)

#preds<-predict(cvfit,newx=kx_tst,s="lambda.min",type="response")
coef(cvfit,s="lambda.min")

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


CUT_TEST<-x_tst
grid_p<-as.matrix(grid)
grid_p<-kernelMatrix(rbf,grid_p,x)
func=predict(cvfit,newx=grid_p,s="lambda.min",type="response")
names<-colnames(func)
func<-matrix(func,ncol=ncol(func))
colnames(func)<-names
func<-data.frame(grid,func)
kx_tst<-kernelMatrix(rbf,CUT_TEST,x)
tmp=predict(cvfit,newx=kx_tst,s="lambda.min",type="response")
names<-colnames(tmp)
tmp<-matrix(tmp,ncol=ncol(tmp))
colnames(tmp)<-names
class<-apply(tmp,1,function(b){
    which(b==max(b))
})
class<-colnames(tmp)[class]
  
response<-y_tst
concordance<-numeric(length(response))
  
for (i in 1:nrow(CUT_TEST)){
    concordance[i]=as.numeric(response[i] %in% class[i])                                   
}
correct<-round(mean(concordance),3)*100
  
  
plot_dat<-melt(func,id=c("H1","H2"))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",race,"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)
  print(plot)  
}




##########################################################################
###Apply an elastic net penalty for multiple regression to the
###To the continous response
###So that we can do multiple regression on Y
##########################################################################
####needs to be blanaced not correct as of now####
y=TRAIN[c("ten_of_tens","nine_of_tens")]
y_tst=TEST[c("ten_of_tens","nine_of_tens")]
x=as.matrix(TRAIN[c("H1","H2")])
x<-cbind(x,decodeClassLabels(TRAIN[,c("Race")]))
x_tst=as.matrix(TEST[c("H1","H2")])
x_tst<-cbind(x_tst,decodeClassLabels(TEST[,c("Race")]))

###look at structure of Y's to confirm association
plot(ten_of_tens~nine_of_tens,data=y)


###fit kernel for the x matrix
models<-list()
param<-c(0.25,1,1.5)
for(i in 1:length(param)){
  print(paste("on loop",i))
  rbf<-rbfdot(sigma=param[i])
  kx<-kernelMatrix(rbf,x)
  kx_tst<-kernelMatrix(rbf,x_tst,x)
  cvfit=glmnet(kx,y,family="mgaussian")###x-fold doesn't work, but no big deal
  #cvfit=cv.glmnet(kx,y,family="mgaussian")
  models[[i]]<-cvfit
}

win<-ceiling(sqrt(length(param)))
par(mfrow=c(win,win))
lapply(models,function(b){plot(b)})
cv_result<-lapply(models,function(b){min(b$cvm)})
cv_result<-data.frame("lambda"=param,"Min Deviance"=unlist(cv_result))
print(cv_result)
cvfit<-models[[which(cv_result$Min.Deviance==min(cv_result$Min.Deviance))]]
par(op)

#preds<-predict(cvfit,newx=kx_tst,s="lambda.min",type="response")
coef(cvfit,s="lambda.min")


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

for(race in par_disp){
  idx<-as.logical(x_tst[,colnames(x_tst)==race])
  CUT_TEST<-x_tst[idx, ]
  grid=data.frame(s_grid,CUT_TEST[1,colnames(CUT_TEST) %in% par_disp,drop=F])
  grid<-as.matrix(grid)
  grid<-kernelMatrix(rbf,grid,x)
  func=predict(cvfit,newx=grid,s="lambda.min",type="response")
  names<-colnames(func)
  func<-matrix(func,ncol=ncol(func))
  colnames(func)<-names
  func<-data.frame(s_grid,func)
  kx_tst<-kernelMatrix(rbf,CUT_TEST,x)
  tmp=predict(cvfit,newx=kx_tst,s="lambda.min",type="response")
  names<-colnames(tmp)
  tmp<-matrix(tmp,ncol=ncol(tmp))
  colnames(tmp)<-names
  class<-apply(tmp,1,function(b){
    which(b==max(b))
  })
  class<-colnames(tmp)[class]
  
  response<-y_tst[idx]
  concordance<-numeric(length(response))
  
  for (i in 1:nrow(CUT_TEST)){
    concordance[i]=as.numeric(response[i] %in% class[i])                                   
  }
  correct<-round(mean(concordance),3)*100
  
  
  plot_dat<-melt(func,id=c("H1","H2"))
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


####This needs to be class balanced####

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


##########################################################################
#####Now Apply a custom kernel to the data
#####Using a SNP variant type data Set
#####Only for example
##########################################################################
###this needs to be balanced###

source("gausspr_mod.R")

knew<-function (scale=0.001) 
{
  rval <- function(x, y = NULL) {
    if (!is(x, "vector")) 
      stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y)) 
      stop("y must a vector")
    if (is(x, "vector") && is.null(y)) {
      return(1)
    }
    if (is(x, "vector") && is(y, "vector")) {
      if (!length(x) == length(y)) 
        stop("number of dimension must be the same on both data points")
      return((sum(x*y) +1)*exp(-scale*sum((x-y)^2)))
    }
  }
  return(new("kernel", .Data = rval, kpar = list(scale = scale)))
}


k<-knew()


bestmod<-gausspr_mod(Productivity~H1+H2+Race,data=TRAIN,kernel=k)




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



