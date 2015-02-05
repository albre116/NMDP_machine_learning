#####################################################
####Boosting Tutorial
#####################################################
source("data_prep.R")
mu=c(1,2)#mean of simulated data
sig=matrix(c(1,0,0,2),nrow=2,byrow=T)
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

###first lets do the prior example with the mboost package
###start with the linear regression model using least squares
lm <- lm(y~x,data=data)
coef(lm)

###boosted version
glm1 <- glmboost(y~x,data=data,control=boost_control(mstop=100))
coef(glm1,off2int = T)
###look at coefficient paths as a function of boosting iterations
plot(glm1,off2int = T)

###Now lets say we didn't know that X2 was inactive
data<-data.frame(y=y,x)
glm1 <- glmboost(y~.,data=data,control=boost_control(mstop=100))
coef(glm1,off2int = T)
###look at coefficient paths as a function of boosting iterations
plot(glm1,off2int = T)
###this did a nice job of regularization


#######Less of a toy problem
mu=c(1,2)#mean of simulated data
sig=matrix(c(1,0,0,2),nrow=2,byrow=T)
x<-mvrnorm(n=1000,mu=mu,Sigma=sig)
y<-1.5+3*sin(2*pi*x[,1])+rnorm(nrow(x),mean=0,sd=1) ###make simple association with x1

data=data.frame(x,y)
ggplot(data,aes(x=X1,y=y))+geom_point()+ggtitle("Scatterplot of Y vs X")

gam1 <- gamboost(y~bols(X1)+bols(X2,intercept=F),data=data)
coef(gam1,off2int = T)
###not very good
preds <- data.frame(f_x=predict(gam1),data.frame(model.frame(gam1))["X1"]) %>% arrange(X1)
ggplot(data,aes(x=X1,y=y))+geom_point()+
  ggtitle("Fit to Data of Linear Bols") +
  geom_line(data=preds,aes(x=X1,y=f_x))

####lets do a spline basis for X1 & X2
gam1 <- gamboost(y~bols(X1)+bols(X2)+bbs(X1,df=20)+bbs(X2,df=20),data=data)
plot(gam1)
###better
preds <- data.frame(f_x=predict(gam1),data.frame(model.frame(gam1))["X1"]) %>% arrange(X1)
ggplot(data,aes(x=X1,y=y))+geom_point()+
  ggtitle("Fit to Data of Linear Bols") +
  geom_line(data=preds,aes(x=X1,y=f_x))


###we can also do this using trees
###the gbm package is better for this
###fit a boosted model
gbm1 <-
  gbm(y~X1+X2,         # formula
      data=data,                   # dataset
      var.monotone=c(0,0), # -1: monotone decrease,
      # +1: monotone increase,
      #  0: no monotone restrictions
      distribution="gaussian",     # see the help for other choices
      n.trees=1000,                # number of trees
      shrinkage=0.025,              # shrinkage or learning rate,
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


# check performance using 3-fold cross-validation
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)


# plot the performance # plot variable influence
summary(gbm1,n.trees=1)         # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees

# compactly print the first and last trees for curiosity
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))

# create marginal plots
# plot variable X1,X2 after "best" iterations
par(mfrow=c(1,2))
plot(gbm1,1,best.iter)
plot(gbm1,2,best.iter)
par(mfrow=c(1,1))


preds <- data.frame(f_x=predict(gbm1),X1=data$X1) %>% arrange(X1)
ggplot(data,aes(x=X1,y=y))+geom_point()+
  ggtitle("Fit to Data of Tree Basis Bols") +
  geom_line(data=preds,aes(x=X1,y=f_x))



####################################################
####Boosting using trees with GBM package
####For our exact problem
####################################################
source("data_prep.R")

###construct weights function for full model
###to be fit at once
weights <-as.data.frame(table(TRAIN$Race,TRAIN$Productivity))
weights$weight <- unlist(lapply(weights$Freq,function(x){
  1/x
}))

TRAIN$Weight <- 0
for (r in weights$Var1){
  for(p in weights$Var2){
    w <- weights$weight[ weights$Var1==r & weights$Var2==p]
    TRAIN$Weight[ TRAIN$Race==r & TRAIN$Productivity==p] <- w
  }
}


###fit a boosted model
gbm1 <-
  gbm(Productivity~H1+H2+Race,         # formula
      data=TRAIN,                   # dataset
      var.monotone=c(0,0,0), # -1: monotone decrease,
      # +1: monotone increase,
      #  0: no monotone restrictions
      weights=Weight,
      distribution="multinomial",     # see the help for other choices
      n.trees=1000,                # number of trees
      shrinkage=0.01,              # shrinkage or learning rate,
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
  cols <- rev(brewer.pal(11, 'RdYlBu'))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Latent Logit for",race,"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)+ scale_fill_gradientn(colours = cols)
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
  cols <- rev(brewer.pal(11, 'RdYlBu'))
  plot<-ggplot(data=plot_dat,aes(x=H1,y=H2,fill=value))+
    geom_tile()+
    stat_contour(data=plot_dat,aes(z=value))+
    ggtitle(paste("Probability Distribution for",race,"Percent Correct=",correct,"%"))+
    facet_wrap(~variable,ncol=2)+ scale_fill_gradientn(colours = cols)
  print(plot)  
}




####################################
####Work with BART Package
####Lets look at some genomic data
####################################
##simulate data (example from Friedman MARS paper)
f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
set.seed(99)
bartFit = bart(x,y,ndpost=200) #default is ndpost=1000, this is to run example fast.
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))

##simulate data (example from Friedman MARS paper)
##simulate data 
f = function(x) { return(.5*x[,1] + 2*x[,2]*x[,3]) }
sigma=.2 # y = f(x) + sigma*z
n=100 #number of observations
set.seed(27)
x = matrix(2*runif(n*3)-1,ncol=3) ; colnames(x) = c('rob','hugh','ed')
Ey = f(x)
y = Ey +  sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
par(mfrow=c(1,3)) #first two for pdbart, third for pd2bart
##pdbart: one dimensional partial dependence plot
set.seed(99)
pdb1 = pdbart(x,y,xind=c(1,2),
              levs=list(seq(-1,1,.2),seq(-1,1,.2)),pl=FALSE,
              keepevery=10,ntree=100,nskip=100,ndpost=200) #should run longer!
plot(pdb1,ylim=c(-.6,.6))
##pd2bart: two dimensional partial dependence plot
set.seed(99)
pdb2 = pd2bart(x,y,xind=c(2,3),
               levquants=c(.05,.1,.25,.5,.75,.9,.95),pl=FALSE,
               ntree=100,keepevery=10,verbose=FALSE,nskip=100,ndpost=200) #should run longer!
plot(pdb2)
##compare BART fit to linear model and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,pdb1$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))


###pull in a balanced data set of genomic expression
###for non-balanced you should be using weights
###unless you expect your margin to come into play
data(promotergene)
summary(promotergene)

x=promotergene %>% select(-Class)
y=promotergene %>% select(Class)
y=y$Class  ###we have to make the response a factor instead of a data.frame


bfit <- bart(x,y)












