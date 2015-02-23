#########makes a shiny demo of the basic parameters
#########for boosting (at least the idea of gradient descent)
###http://rstudio.github.io/shinydashboard/get_started.html

header <- dashboardHeader(title = "Boosting Demo")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Simulation Settings",tabName="Simulation",icon = icon("dashboard")),
    menuItem("Gradient Descent",tabName="Gradient",icon = icon("th"))
    )
  )

body <- dashboardBody(
  tabItems(
    #first tab content
    tabItem(tabName="Simulation",
            box(status = "primary",title = "Settings", solidHeader = TRUE, width = NULL,collapsible = T,
                fluidRow(
                  column(width=2,
                         h6("Mean Values"),
                         htable("means",colHeaders = "provided"),
                         h6("Variance-Covariance Matrix"),
                         htable("varmatrix",colHeaders = "provided",rowNames = "provided"),
                         numericInput("nsim","Number of Simulations",value=1000,min=10,max=10000,step=1),
                         h6("Y=intercept+B1*X1+0*X2+N(0,error_var)"),
                         numericInput("intercept","intercept",value=1.5),
                         numericInput("B1","B1",value=3),
                         numericInput("error_var","error variance",value=1),
                         actionButton("update","Update Parameters",icon = icon("refresh"))
                         ),
                  column(width=10,
                         plotOutput("xyplot")
                         )
                   )
                ),###end box content
            box(status = "primary",title = "Loss Surface", solidHeader = TRUE, width = NULL,collapsible = T,
                fluidRow(
                  column(width=2,
                         numericInput("B0min","Min B0",0),
                         numericInput("B0max","Max B0",3),
                         numericInput("B1min","Min B1",1.5),
                         numericInput("B1max","Max B1",4.5)
                  ),
                  column(width=10,
                         plotOutput("loss_surface")
                  )
                )
            )###end box content
            ),
    #second tab content
    tabItem(tabName="Gradient",
            box(status = "primary",title = "Settings", solidHeader = TRUE, width = NULL,collapsible = T,
                fluidRow(
                  column(width=2,
                         h6("Starting Values"),
                         htable("initial",colHeaders = "provided"),
                         numericInput("learnrate","Learning Rate",0.1,min=0),
                         numericInput("boostiter","Boosting Iterations",5,min=1),
                         actionButton("update_init","Run Simulation",icon = icon("refresh"))
                         
                ),
                column(width=10,
                       plotOutput("boost_fit")
                       )
                )
    ),###end box
    box(status = "primary",title = "Descent Path", solidHeader = TRUE, width = NULL,collapsible = T,
        plotOutput("descent_path")
        
        
        )###end box
    )
    
  )###end all tab items
)






server <-   function(input, output, session) {
  ###generic functions
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
  
  
  ####end generic functions
  

  output$means <- renderHtable({
    if (is.null(input$means)){
      means <- as.data.frame(matrix(c(1,2),nrow=1,byrow=T))
      colnames(means)<-c("X1","X2")
      return(means)
    } else{
      return(input$means)
    }
  }) 
  
  output$varmatrix <- renderHtable({
    if (is.null(input$varmatrix)){
      varmatrix <- as.data.frame(matrix(c(1,0,0,2),nrow=2,byrow=T))
      colnames(varmatrix)<-c("X1","X2")
      rownames(varmatrix)<-c("X1","X2")
      return(varmatrix)
    } else{
      return(input$varmatrix)
    }
  }) 
  
  
  PARAMS <- reactive({
    input$update
    isolate(if (is.null(input$varmatrix)){
      varmatrix <- as.data.frame(matrix(c(1,0,0,2),nrow=2,byrow=T))
      colnames(varmatrix)<-c("X1","X2")
      rownames(varmatrix)<-c("X1","X2")}else{
        varmatrix <- input$varmatrix
      })
    
    isolate(if (is.null(input$means)){
      means <- as.data.frame(matrix(c(1,2),nrow=1,byrow=T))
      colnames(means)<-c("X1","X2")}else{
        means <- input$means
      })
      
    out <- list(mu=means,sigma=varmatrix)
    return(out)
  })
  
  SIM <- reactive({
    mu <- as.matrix(PARAMS()[["mu"]])
    r <- nrow(mu)
    c <- ncol(mu)
    mu <- matrix(as.numeric(unlist(mu)),nrow=r,ncol=c,byrow = F)
    sigma <- as.matrix(PARAMS()[["sigma"]])
    r <- nrow(sigma)
    c <- ncol(sigma)
    sigma <- matrix(as.numeric(unlist(sigma)),nrow=r,ncol=c,byrow = F)
    isolate(x<-mvrnorm(n=input$nsim,mu=mu,Sigma=sigma))
    isolate(y<-input$intercept+input$B1*x[,1]+rnorm(nrow(x),mean=0,sd=input$error_var)) ###make simple association with x
    return(list(x=x,y=y))
  })
  
  output$xyplot <- renderPlot({
    x <- SIM()[["x"]]
    y <- SIM()[["y"]]
    data=data.frame(x=x[,1],y=y)
    plot <- ggplot(data,aes(x,y))+geom_point()+ggtitle("Scatterplot of Y vs X1")
    print(plot)
  })


  BETA_GRID <- reactive({
    x <- SIM()[["x"]]
    y <- SIM()[["y"]]
    data=data.frame(x=x[,1],y=y)
    beta_grid<-expand.grid(seq(input$B0min,input$B0max,length.out = 50),
                           seq(input$B1min,input$B1max,length.out=50))
    colnames(beta_grid) <- c("B0","B1")
    beta_grid$Loss<-apply(beta_grid,1,L2_Loss,y=data$y,x=data$x)
    return(beta_grid)
  })

  
  
  output$loss_surface <- renderPlot({
    beta_grid <- BETA_GRID()
    plot <- ggplot(beta_grid,aes(x=B0,y=B1,z=Loss))+geom_contour()
    print(plot)
  })

  
  ####code up the gradient descent part
  
  output$initial <- renderHtable({
    if (is.null(input$initial)){
      initial <- as.data.frame(matrix(c(0,1.5),nrow=1,byrow=T))
      colnames(initial)<-c("B0","B1")
      return(initial)
    } else{
      return(input$initial)
    }
  }) 
  
  INIT <- reactive({
    input$update_init
    isolate(if (is.null(input$initial)){
      means <- as.data.frame(matrix(c(0,1.5),nrow=1,byrow=T))
      colnames(means)<-c("B0","B1")}else{
        means <- input$initial
      })
    out <- means
    return(out)
  })
  
  
  GRADIENT <- reactive({
    beta <- as.matrix(INIT())
    x <- SIM()[["x"]]
    y <- SIM()[["y"]]
    data=data.frame(x=x[,1],y=y)
    r <- nrow(beta)
    c <- ncol(beta)
    beta <- matrix(as.numeric(unlist(beta)),nrow=r,ncol=c,byrow = F)
    isolate(learn_rate <- input$learnrate)
    isolate(nboost <- input$boostiter)
    beta_history <- matrix(nrow=nboost+1,ncol=2)
    colnames(beta_history) <- c("B0","B1")
    beta_history <- as.data.frame(beta_history)
    beta_history[1,] <- beta
    
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
    
    return(list(beta=beta,nboost=nboost,beta_history=beta_history))
    
  })
  
  output$boost_fit <- renderPlot({
    beta <- GRADIENT()[["beta"]]
    nboost <- GRADIENT()[["nboost"]]
    x <- SIM()[["x"]]
    y <- SIM()[["y"]]
    data=data.frame(x=x[,1],y=y)
    fit_slope <- paste0("Y=",round(beta[1],3)," + ",round(beta[2],3),"*X")
    ####look at fit after M boosting iterations
    f_x <- predict_fit(beta,data$x)
    plot_f_x<-data.frame(x=data$x,y=f_x) %>% arrange(x)
    plot <- ggplot(data,aes(x,y))+geom_point()+
      geom_line(data=plot_f_x,aes(x=x,y=y))+
      ggtitle(paste("Model Fit: (f_x)=",fit_slope,"  After M=",nboost,"Boosting Iterations"))
    print(plot)
    
  })
  
  
  output$descent_path <- renderPlot({
    beta <- GRADIENT()[["beta"]]
    nboost <- GRADIENT()[["nboost"]]
    beta_history <- GRADIENT()[["beta_history"]]
    x <- SIM()[["x"]]
    y <- SIM()[["y"]]
    data=data.frame(x=x[,1],y=y)
    beta_grid <- BETA_GRID()
    ####look at the descent path of the beta values
    ####as the simple base learners are combined
    plot <- ggplot()+geom_contour(data=beta_grid,aes(x=B0,y=B1,z=Loss))+
      geom_line(data=beta_history,aes(x=B0,y=B1),lwd=1)+
      geom_point(data=beta_history,aes(x=B0,y=B1),size=3)+
      ggtitle(paste("Gradient Descent Path Across Loss Function After M=",nboost,"Boosting Iterations"))
    print(plot)
  })
}

