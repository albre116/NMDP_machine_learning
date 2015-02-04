#####################################################
####Boosting Tutorial
#####################################################
source("data_prep.R")
mu=c(1,2)#mean of simulated data
sig=matrix(c(1,0.5,0.5,2),nrow=2,byrow=T)
x<-mvrnorm(n=1000,mu=mu,Sigma=sig)
y<-1.5+3*x[,1]+rnorm(nrow(x),mean=0,sd=1) ###make simple association with x1

data=data.frame(x=x[,1],y=y)
ggplot(data,aes(x,y))+geom_point()+ggtitle("Scatterplot of Y vs X")

L2_Loss<-function(Beta,y,x){
  residual <- (1/2)*(y-(Beta[1]+Beta[2]*x))^2  
  SSE <- sum(residual)
  return(-SSE)
}

Gradient_Residual<-function(Beta,y,x){
  residual <- -(y-(Beta[1]+Beta[2]*x))  
  return(residual)
}

predict_fit <- function(Beta,x){
  f_x <- Beta[1]+Beta[2]*x
  return(f_x)
}


beta_grid<-expand.grid(seq(0,3,length.out = 50),
                       seq(1.5,4.5,length.out=50))
colnames(beta_grid) <- c("B0","B1")
beta_grid$Loss<-apply(beta_grid,1,L2_Loss,y=data$y,x=data$x)
####look at the true loss function 
ggplot(beta_grid,aes(x=B0,y=B1,z=Loss))+geom_contour()



beta=c(0,1.5)###  starting point for all values
learn_rate=0.1
nboost=100
beta_history <- matrix(nrow=nboost+1,ncol=2)
colnames(beta_history) <- c("B0","B1")
beta_history <- as.data.frame(beta_history)
beta_history[1,] <- beta

####look at starting point
f_x <- predict_fit(beta,data$x)
plot_f_x<-data.frame(x=data$x,y=f_x) %>% arrange(x)
ggplot(data,aes(x,y))+geom_point()+ggtitle("Starting Model Fit")+
  geom_line(data=plot_f_x,aes(x=x,y=y))

###run boosting
for (i in 1:nboost){

###compute the negative gradient
neg_grad<--(Gradient_Residual(beta,y=data$y,x=data$x))

###fit the simple base learner
grad_plot <- data.frame(x=data$x,neg_grad=neg_grad)
fit <- lm(neg_grad~x,data=grad_plot)
beta[1] <- beta[1]+learn_rate*fit$coefficients[1]
beta[2] <- beta[2]+learn_rate*fit$coefficients[2]
beta_history[1+i,] <- beta
}

####look at fit after M boosting iterations
f_x <- predict_fit(beta,data$x)
plot_f_x<-data.frame(x=data$x,y=f_x) %>% arrange(x)
ggplot(data,aes(x,y))+geom_point()+
  geom_line(data=plot_f_x,aes(x=x,y=y))+
  ggtitle(paste("Model Fit (f_x) after M=",nboost,"Boosting Iterations"))

####look at the gradient residuals of the last fit to which the base
####learner is being fit
ggplot(data=grad_plot,aes(x=x,y=neg_grad))+geom_point()+
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2])+
  ggtitle(paste("Gradient Residuals After M=",nboost,"Boosting Iterations and Base Learner Fit"))

####look at the descent path of the beta values
####as the simple base learners are combined
ggplot()+geom_contour(data=beta_grid,aes(x=B0,y=B1,z=Loss))+
  geom_line(data=beta_history,aes(x=B0,y=B1),lwd=1)+
  geom_point(data=beta_history,aes(x=B0,y=B1),size=3)+
  ggtitle(paste("Gradient Descent Path Across Loss Function After M=",nboost,"Boosting Iterations"))


###Print final model fit
print(paste0("Y=",beta[1]," + ",beta[2],"*X"))

####################################################
####Now lets make this more realistic
####And use the mboost package
####to specify other base learners
####################################################






####################################################
####Boosting using trees with GBM package
####################################################
source("data_prep.R")

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
race=par_disp[2]

###look at the latent funciton F(X) that is squished by the sigmoid softmax function

for(race in par_disp){
  CUT_TEST<-TEST[TEST$Race==race, ]
  grid$Race=CUT_TEST$Race[1]
  func=predict(bestmod,newdata = grid,best.iter,type ="link" )
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
  
  
  plot_dat<-melt(func,id=c("H1","H2","Race"))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",levels_race[race],"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)
  print(plot)  
}

###look at the response surface after applying the softmax function

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
  
  
  plot_dat<-melt(func,id=c("H1","H2","Race"))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",levels_race[race],"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)
  print(plot)  
}




