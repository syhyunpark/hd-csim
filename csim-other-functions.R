library(SAM)

###############################################################
# Modified covariate method with efficiency augmentationl;
# mc is a function for fitting the modified covariate (linear) model of Tian et al. (2015)
# for estimating the treatment-by-covariates interactions,
# with efficiency augmentation.
###############################################################
mc <- function(y, Tr, X, eff.aug = TRUE, X.aug= NULL, use.lasso=TRUE, n.folds=20, boot.CI=FALSE, n.boot=500)
{
  
  ################
  ## estimation ##
  ################
  mc.obj <- fit.mc(y,Tr, X, use.lasso=use.lasso, n.folds=n.folds, eff.aug=eff.aug, X.aug=X.aug)
  
  ##################################
  ## Bootstrap CIs for alpha.coef ##
  ##################################
  boot.results <- NULL
  if(boot.CI)
  {
    datax <- data.frame(y=y, Tr=Tr, X=X)
    # define the "bs" function for bootstraping
    bs <- function(datax, indices,
                   use.lasso=FALSE, n.folds=20, eff.aug=FALSE)
    {
      mc.obj.b <- fit.mc(y=datax[indices,][,1], Tr= datax[indices,][,2], X=datax[indices,][,-c(1,2)],
                         use.lasso=use.lasso, n.folds=n.folds, eff.aug=eff.aug)
      return(mc.obj.b$alpha.coef)
    }
    boot.results <- boot(data = datax, statistic = bs, R = n.boot,
                         use.lasso=use.lasso, n.folds=n.folds, eff.aug=eff.aug)
    #boot.ci(boot.results, type = "norm", index = 1)
  }
  
  mc.obj$boot.results <- boot.results
  
  return(mc.obj)
}




###############################################################
# Modified covariate method with efficiency augmentationl;
# mc is a function for fitting the modified covariate (linear) model of Tian et al. (2015)
# for estimating the treatment-by-covariates interactions,
# with efficiency augmentation.
###############################################################
fit.mc <- function(y, Tr, X, eff.aug = TRUE, X.aug = NULL, use.lasso=TRUE, n.folds=20)
{
  
  K <- length(unique(Tr))
  if(K > 2)  stop("The number of treatment groups should be 2.")
  
  # center and scale X
  Xc <- scale(X, center = TRUE, scale = TRUE)
  scale.X <-  attr(Xc, "scaled:scale")
  center.X <- attr(Xc, "scaled:center")
  
  # center y within each treatment, and store the treatment-specific mean at intercept.y
  dat <- data.frame(y=y, Tr=Tr, Xc= Xc)
  dat.list <- dlply(dat, .(dat$Tr), function(x) as.matrix(x));
  intercept.y <- vector("list", K); 
  yc = Xc = Tr = NULL; 
  for(t in 1:K)
  {
    tmp <- scale(dat.list[[t]][,1], center = TRUE, scale = FALSE)
    intercept.y[[t]] <- attr(tmp, "scaled:center")
    yc <- c(yc, tmp);  rm(tmp)
    Tr <- c(Tr, dat.list[[t]][,2]);
    Xc <- rbind(Xc, dat.list[[t]][,-c(1,2)])
  }
  rm(dat); rm(dat.list)
  
  # efficiency augmentation
  eta.coef <- NULL
  if(eff.aug)
  {
    if(is.null(X.aug)) X.aug <- Xc; 
    X.aug <- as.matrix(X.aug) 
    
    if(use.lasso)
    {
      lasso.temp <- cv.glmnet(x=X.aug, y=yc, nfolds = n.folds)
      lamda.opt <- lasso.temp$lambda[which.min(lasso.temp$cvm)]
      lasso.fit <- glmnet(x=X.aug, y=yc, lambda =lamda.opt, alpha=1)
      eta.coef <- as.numeric(lasso.fit$beta)
      yc  <- yc - predict(lasso.fit, X.aug)
    }else{
      lm.fit <- lm(yc ~ X.aug)
      eta.coef <- coef(lm.fit)[-1]
      yc <- yc - predict(lm.fit)
    }
  }
  
  # modify the covariates
  Xc.modified <- as.matrix(Xc * (Tr-1.5) )
  
  # estimate the interaction term
  if(use.lasso)
  {
    lasso.temp <- cv.glmnet(x= Xc.modified, y= yc, nfolds = n.folds)
    lamda.opt <- lasso.temp$lambda[which.min(lasso.temp$cvm)]
    lasso.fit <- glmnet(x= Xc.modified, y= yc, lambda = lamda.opt, alpha = 1)
    alpha.coef <- as.numeric(lasso.fit$beta)
    MSE <- mean((yc - predict(lasso.fit, Xc.modified) )^2)
  }else{
    lm.fit2 <- lm(yc  ~ Xc.modified)
    alpha.coef <- coef(lm.fit2)[-1]
    MSE <- mean((yc - predict(lm.fit2))^2)
  }
  names(alpha.coef) <- colnames(X)
  
  results <- list(alpha.coef = alpha.coef, eta.coef= eta.coef, MSE = MSE,
                  intercept.y = intercept.y, scale.X = scale.X, center.X = center.X)
  
  class(results) <- c("mc", "list" )
  
  return(results)
}



###############################################################
# mc prediction function;
# predict.mc makes predictions from the modfiied covariate model
# given a \code{mc} object and pretreatment covariates.
# The function returns predicted outcomes for each treatment and treatment selection rules.
###############################################################
predict.mc <- function(mc.obj, newx)
{
  if(!inherits(mc.obj, "mc"))
    stop("Object must be of class `mc' ")
  
  newx.scaled <- scale(newx, center = mc.obj$center.X, scale = mc.obj$scale.X);
  
  # compute treatment-specific predicted value
  pred.new <- matrix(0, nrow(newx), 2)
  for(t in 1:2)
  {
    if(is.null(mc.obj$eta.coef))
    {
      pred.new[, t] <- mc.obj$intercept.y[[t]] +  (newx.scaled *(t-1.5)) %*% mc.obj$alpha.coef
    }else{
      pred.new[, t] <- mc.obj$intercept.y[[t]] +  (newx.scaled *(t-1.5)) %*% mc.obj$alpha.coef + newx.scaled %*% mc.obj$eta.coef
    }
  }
  
  # compute treatment assignment
  trt.rule <- apply(pred.new, 1, which.max)  # WLLG, we assume a higher y is prefered.
  colnames(pred.new) <- c("Tr1", "Tr2")
  results <- list(trt.rule = trt.rule, pred.new = pred.new)
  
  return(results)
}




#####################################################################
# K-separate linear regression models
# K.LR is a function for fitting a system of K separate linear regression models,
# one for each treatment group, for estimating the treatment-by-covariates interaction effects.
#####################################################################
K.LR <- function(y, Tr, X, use.lasso=TRUE)
{
  K <- length(unique(Tr))
  dat <- data.frame(Tr=Tr, y=y, X=X)
  dat_list <- dlply(dat,.(dat$Tr),function(x) as.matrix(x[,-1]) )  # list
  
  K.LR.obj <- vector("list", K)
  for(t in 1:K)
  {
    if(use.lasso)
    {
      lm.temp <- cv.glmnet(x= dat_list[[t]][,-1], y=dat_list[[t]][,1], nfolds=20)
      K.LR.obj[[t]] <- glmnet(x=  dat_list[[t]][,-1], y=dat_list[[t]][,1], lambda = lm.temp$lambda[which.min(lm.temp$cvm)], alpha = 1)
    }else{
      K.LR.obj[[t]] <- lm(dat_list[[t]][,1] ~ dat_list[[t]][,-1])
    }
  }
  
  K.LR.obj$use.lasso <- use.lasso;
  K.LR.obj$K <- K;
  
  class(K.LR.obj) <- c("K.LR", "list" )
  return(K.LR.obj)
}



#####################################################################
# L.LR prediction function;
# it makes predictions from the fitted K separate linear models,
# given a K.LR object and pretreatment covariates.
# The function returns predicted outcomes for each treatment and treatment selection rules.
#####################################################################
predict.K.LR <- function(K.LR.obj, newx)
{
  if(!inherits(K.LR.obj, "K.LR"))
    stop("Object must be of class `K.LR' ")
  
  K <- K.LR.obj$K
  n <- nrow(newx)
  pred_K.LR_list <-  vector("list", K)
  newx <- as.matrix(newx)
  for(t in 1:K)
  {
    if(K.LR.obj$use.lasso)
    {
      pred_K.LR_list[[t]] <- predict(K.LR.obj[[t]], newx)
    }else{
      pred_K.LR_list[[t]] <-  as.matrix(cbind(1, newx)) %*% K.LR.obj[[t]]$coefficients 
    }
  }
  
  pred.new  <- NULL
  for(t in 1:K)  pred.new  <- cbind(pred.new, pred_K.LR_list[[t]] )
  
  regime <- apply(pred.new, 1, which.max)  # WLLG, we assume a higher y is prefered.
  results <- list(regime = regime, pred.new = pred.new)
  
  return(results)
}



#####################################################################
# K-separate sparse additive models
# K.SAM is a function for fitting a system of K separate (sparse) additve models (Ravikumar et al. 2009),
# one for each treatment group, for estimating the treatment-by-covariates (possibly nonlinear) interaction effects.
#####################################################################
K.SAM <- function(y, Tr, X, lambda.opt = NULL, n.folds = 5)
{
  
  p <- ncol(X)
  K <- length(unique(Tr))
  
  # center and scale X
  Xc <- scale(X, center = TRUE, scale = TRUE)
  scale.X <-  attr(Xc, "scaled:scale")
  center.X <- attr(Xc, "scaled:center")
  
  # create a list, dat.list, grouped by the treatment indicator
  dat <- data.frame(Tr=Tr, y=y, Xc = Xc);
  dat.list <- dlply(dat, .(dat$Tr), function(x) {as.matrix(x[,-1])})
  K.SAM.fit =  range.info <- vector("list", K)
  
  for(t in 1:K)
  {
    datax.t <- data.frame(t, dat.list[[t]])  # data from the tth treatment group
    range.info[[t]] <-  apply(datax.t[,-c(1,2)], 2, range)
    
    if(is.null(lambda.opt))
    {
      # perform n.folds cross validation to chooose lambda
      #n.folds <- 5;
      dataxx <- datax.t[sample(nrow(datax.t)), ]  # randomly shuffle the datax
      folds <- cut(seq(1,nrow(dataxx)), breaks= n.folds, labels=FALSE)
      PE.storage <- matrix(NA, n.folds, 30 )
      for(i in 1:n.folds)
      {
        testIndexes <- which(folds==i, arr.ind=TRUE)
        testData <- dataxx[testIndexes, ]
        trainData <- dataxx[-testIndexes, ]
        range.train <- apply(trainData[,-c(1,2)], 2, range)
        for(j in 1:p)
        {
          testData[,-c(1,2)][ ,j] [ testData[,-c(1,2)][ ,j] <= range.train[1,j] ] <- range.train[1,j]
          testData[,-c(1,2)][ ,j] [ testData[,-c(1,2)][ ,j] >= range.train[2,j] ] <- range.train[2,j]
        }
        K.SAM.t.temp <- samQL(trainData[,-c(1,2)] , trainData$y)
        K.SAM.t.temp$knots <- K.SAM.t.temp$nkots
        pred.K.SAM.temp <- predict(K.SAM.t.temp, testData[,-c(1,2)])$values
        for(l in 1:30)
        {
          PE.storage[i, l] <- sum((testData$y - pred.K.SAM.temp[,l])^2)
        }
      }
      lambdas <- K.SAM.t.temp$lambda
      mean.PE <- apply(PE.storage, 2, mean)
      lambda.opt <- lambdas[which.min(mean.PE)]
    }
    
    K.SAM.t <- samQL(datax.t[,-c(1,2)] , datax.t[,2], lambda = lambda.opt)
    K.SAM.fit[[t]] <- K.SAM.t
    
    class(K.SAM.fit[[t]]) <- c("samQL", "list")
  }
  results <- list(scale.X =scale.X, center.X= center.X,
                  range.info= range.info, K.SAM.fit = K.SAM.fit,
                  p=p, K=K)
  class(results) <- c("K.SAM", "list" )
  
  return(results)
}


#####################################################################
# K.SAM prediction function;
# itmakes predictions from the fitted K separate sparse additive models,
# given a K.SAM object and pretreatment covariates.
# The function returns predicted outcomes for each treatment and treatment selection rules.
#####################################################################
predict.K.SAM <- function(K.SAM.obj, newx)
{
  if(!inherits(K.SAM.obj, "K.SAM"))
    stop("Object must be of class `K.SAM' ")
  
  K.SAM.fit  <- K.SAM.obj$K.SAM.fit
  K <- K.SAM.obj$K
  p <- K.SAM.obj$p
  newx.scaled <- scale(newx, center = K.SAM.obj$center.X, scale = K.SAM.obj$scale.X);
  
  dat_K.SAM_list <- vector("list", K)
  pred_K.SAM_list <-  vector("list", K)
  for(t in 1:K)
  {
    for(j in 1:p)
    {
      newx.scaled[ ,j] [ newx.scaled[ ,j] <= K.SAM.obj$range.info[[t]][1,j] ] <- K.SAM.obj$range.info[[t]][1,j]
      newx.scaled[ ,j] [ newx.scaled[ ,j] >= K.SAM.obj$range.info[[t]][2,j] ] <- K.SAM.obj$range.info[[t]][2,j]
    }
    newx.scaled <- as.data.frame(newx.scaled)
    pred_K.SAM_list[[t]] <- predict(K.SAM.fit[[t]], newx.scaled)$values
  }
  pred.new  <- NULL
  for(t in 1:K)   pred.new  <- cbind(pred.new, as.vector(pred_K.SAM_list[[t]]) )
  
  regime <- apply(pred.new, 1, which.max)  # WLLG, we assume a higher y is prefered.
  results <- list(regime = regime, pred.new = pred.new)
  return(results)
}



#####################################################################
# Linear generated effect modifiers method (workhorse);
# lgem is a function for fitting the generated treatment effect modifiers (GEM) model
# by the "Numerator" criterion of Petkova et al. (2016).
# The GEM variable is defined as a linear combination of pretreatment covariates X,
# optimized to exhibit a strong (linear) interaction effect with the treatment indicator,
# under a linear model framework.
#####################################################################
fit.lgem <- function(y, Tr, X, use.lasso = FALSE, eff.aug= FALSE, n.folds=20)
{
  
  K <- length(unique(Tr));
  p <- ncol(X)
  
  # center and scale X
  Xc <- scale(X, center = TRUE, scale = TRUE)
  scale.X <-  attr(Xc, "scaled:scale")
  center.X <- attr(Xc, "scaled:center")
  
  # center y within each treatment, and store the treatment-specific mean at intercept.y
  dat <- data.frame(Tr=Tr, y=y, Xc = Xc)
  dat.list <- dlply(dat, .(dat$Tr), function(x) {as.matrix(x[,-1])});
  intercept.y <- vector("list", K);
  yc = Xc <- NULL
  for(t in 1:K)
  {
    dat.list[[t]][,1] = temp <- scale(dat.list[[t]][,1], center = TRUE, scale = FALSE)
    intercept.y[[t]] <- attr(temp, "scaled:center");
    rm(temp);
    yc <- c(yc, dat.list[[t]][,1])
    Xc <- rbind(Xc, dat.list[[t]][,-1])
  }
  
  # efficiency augmentation
  eta.coef <- NULL
  if(eff.aug)
  {
    if(use.lasso)
    {
      lasso.temp <- cv.glmnet(x=Xc, y=yc, nfolds = n.folds)
      lamda.opt <- lasso.temp$lambda[which.min(lasso.temp$cvm)]
      lasso.fit <- glmnet(x=Xc, y=yc, lambda =lamda.opt, alpha=1)
      eta.coef <- drop(lasso.fit$beta)
      yc  <- yc - predict(lasso.fit, Xc)
    }else{
      lm.fit <- lm(yc ~ Xc)
      eta.coef <- coef(lm.fit)[-1]
      yc <- yc - predict(lm.fit)
    }
    
    temp <- data.frame(Tr=Tr, yc=yc)
    temp.list <- dlply(temp, .(temp$Tr), function(x) {as.matrix(x[,-1])});
    for(t in 1:K)
    {
      dat.list[[t]][,1] <- temp.list[[t]][,1];
    }
    rm(temp);
    rm(temp.list)
  }
  
  Beta <- matrix(data = NA, nrow = p, ncol = K)
  if(use.lasso)
  {
    for(t in 1:K)
    {
      lasso.fit <- cv.glmnet(x=  dat.list[[t]][,-1], y=dat.list[[t]][,1], nfolds=n.folds)
      lasso.lambda <- lasso.fit$lambda
      lambda.opt <- which.min(lasso.fit$cvm)
      beta.t <- round(coef(glmnet(x= dat.list[[t]][,-1], y=dat.list[[t]][,1], lambda = lasso.lambda[lambda.opt], alpha = 1)), 1)[-1]
      Beta[,t] = drop(beta.t)
    }
  }else{
    for(t in 1:K)
    {
      beta.t <- lm.ridge( dat.list[[t]][,1]  ~ dat.list[[t]][,-1] -1, lambda = 0.01)$coef
      Beta[,t] = as.vector(beta.t)
    }
  }
  
  Beta.mean <- apply(Beta, 1, mean)
  Beta.c  <- t(scale(t(Beta), center = Beta.mean, scale = FALSE))
  B <- Beta.c %*% t(Beta.c)
  alphax <- eigen(B)$vectors[,1]  # take the first eigenvector
  if(alphax[1]< 0) alphax <- -1*alphax
  
  coef.list = gamma.list <- vector("list")
  for(t in 1:K)
  {
    gamma.list[[t]] <- drop(alphax %*% Beta.c[,t])
    coef.list[[t]] <- alphax *  gamma.list[[t]]
  }
  
  
  results <- list(alpha.coef = alphax,
                  gamma.list = gamma.list,
                  coef.list = coef.list,
                  eta.coef = eta.coef,
                  intercept.y = intercept.y,
                  scale.X = scale.X,
                  center.X = center.X,
                  K=K, p=p)
  
  class(results) <- c("lgem", "list" )
  
  return(results)
}


#####################################################################
# Linear generated effect modifiers method (main);
# lgem is a wrapper function for fitting the generated treatment effect modifiers (GEM) model
# by the "Numerator" criterion of Petkova et al. (2016).
# The GEM variable is defined as a linear combination of pretreatment covariates X,
# optimized to exhibit a strong (linear) interaction effect with the treatment indicator,
# under a linear model framework.
#####################################################################
lgem <- function(y, Tr, X, use.lasso = FALSE,
                 n.folds =20, eff.aug= FALSE, boot.CI=FALSE, n.boot = 200, plots=TRUE)
{
  
  ################
  ## estimation ##
  ################
  lgem.obj <- fit.lgem(y,Tr, X, use.lasso=use.lasso, n.folds=n.folds, eff.aug=eff.aug)
  
  ###################
  ## visualization ##
  ###################
  # plots
  #gt.plot = alpha.coef.plot <- NULL
  if(plots &  (sum(lgem.obj$alpha.coef!=0)!=0) )
  {
    # 1) link function (gt) plots
    dat  <- data.frame(y = y, x = scale(X) %*% lgem.obj$alpha.coef, Treatment = factor(Tr))  # factor(Tr, labels=c("Placebo","Active drug"))
    lgem.obj$gt.plot  <- ggplot(dat, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
      geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
      scale_colour_brewer(palette = "Set1", direction = -1) +
      geom_smooth(method=lm, se=TRUE, fullrange=TRUE, alpha = 0.35) +
      theme( axis.title.x=element_text(size=15,face="bold"))  +
      theme(title =element_text(size=12)) +
      xlab(expression(paste(alpha^T, "x"))) +  #ylab("Improvement") +
      ylab("y") + theme_bw()
    
    # 2) single-index coefficient (alpha.coef) plot
    dotplot_identity = function(frame, xvar, yvar, colour, colorvar=NULL)
    {
      if(is.null(colorvar))  gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar))
      else  gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar, color=colorvar))
      return(gplot + geom_point(colour = colour) + geom_linerange(aes(ymin=0)))
    }
    
    lgem.obj$alpha.plot <- dotplot_identity(frame = data.frame(Z = factor(colnames(X), level = colnames(X)), coef = lgem.obj$alpha.coef), xvar="Z", yvar="coef", colour = "blue") + scale_y_continuous(name = "coef.") + ggtitle(expression(alpha))  + xlab(" ")
    
  }
  
  
  ##################################
  ## Bootstrap CIs for alpha.coef ##
  ##################################
  if(boot.CI)
  {
    datax <- data.frame(y=y, Tr=Tr, X=X)
    # define the "bs" function for bootstraping
    bs <- function(datax, indices,
                   use.lasso=FALSE, n.folds=20, eff.aug=FALSE)
    {
      lgem.obj.b <- fit.lgem(y=datax[indices,][,1], Tr= datax[indices,][,2], X=datax[indices,][,-c(1,2)],
                             use.lasso=use.lasso, n.folds=n.folds, eff.aug=eff.aug)
      return(lgem.obj.b$alpha.coef)
    }
    lgem.obj$boot.results <- boot(data = datax, statistic = bs, R = n.boot,
                                  use.lasso=use.lasso, n.folds=n.folds, eff.aug=eff.aug)
    #boot.ci(lgem.obj$boot.results, type = "norm", index = 1)
  }
  
  
  return(lgem.obj)
}


#####################################################################
# lgem prediction function;
# it makes predictions from the linear GEM model,
# given a lgem object and pretreatment covariates.
# The function returns predicted outcomes for each treatment and treatment selection rules.
#####################################################################
predict.lgem <- function(lgem.obj, newx)
{
  if(!inherits(lgem.obj, "lgem"))
    stop("Object must be of class `lgem' ")
  
  K <- lgem.obj$K
  pred_lgem_list <-  vector("list", K)
  pred.new  <- NULL;
  newx.scaled <- scale(newx, center = lgem.obj$center.X, scale = lgem.obj$scale.X)
  for(t in 1:K)
  {
    temp <- lgem.obj$intercept.y[[t]] + newx.scaled %*% lgem.obj$coef.list[[t]]
    pred.new  <- cbind(pred.new, temp)
    rm(temp)
  }
  
  trt.rule <- apply(pred.new, 1, which.max)  # WLLG, we assume a higher y is prefered.
  
  results <- list(trt.rule = trt.rule, pred.new = pred.new)
  return(results)
}




#####################################################################
# A dataset generation function
# dataGenerationFn generates an example dataset under a model that contains a main effect component,
# a treatment-by-covariates interaction effect component, and a random noise component.
#####################################################################
dataGenerationFn <- function(n = 200, p=50, w = 1, delta = 1,  true.alpha = NULL, true.eta = NULL,
                             sigma = 0.4, correlationX= 0, sigmaX = 1, sim.seed = NULL, obs = FALSE)
{
  if(is.null(sim.seed))   sim.seed <- 10000*runif(1);
  set.seed(sim.seed)
  
  if(is.null(true.alpha))
  {
    true.alpha <- c(c(1, 0.5, 0.25, 0.125), rep(0, p-4));  # only the first 4 components are nonzero.
    #true.alpha <- true.alpha/sqrt(sum(true.alpha^2))
  }
  
  if(length(true.alpha)!= p)   stop("true.alpha must be of length p");
  if(p < 10)   stop("p should be at least 10");
  
  if(is.null(true.eta))  # generate true.eta randomly.
  {
    eta.hold <- rnorm(10, 0, 1);
    eta.hold  <- eta.hold /sqrt(sum(eta.hold^2) )
    eta.hold
    true.eta <- c(eta.hold, rep(0, p-10))   # only the first 10 components are nonzero.
  }
  
  # the link function (that defines the interaction effect component);
  # w is the nonlinearity parameter (w=1: linear; w=2: almost linear; w=3: hihgly nonlinear)
  g <- function(u, w)
  {
    if(w==1) return(0.5* u)
    if(w==2) return(sin(u) - u)
    if(w==3) return(cos(u) - 0.5)
  }
  
  # the main effect function;
  # delta is the intensity parametr (delta = 1: moderate main effevct; delta=2: big main effect)
  m <- function(u, delta= 1)   2+ delta*cos(u*0.5*pi)
  
  
  Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
  ePsix <- eigen(Psix)
  X <- sapply(1:p, function(x) rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)

  if(obs)
  {
    mu.Tr  <- -1 + 0.8*X[,1]^2 + 0.8*X[,2]^2
    p.Tr <- 1/(1+ exp(-mu.Tr) )
    Tr <-  rbinom(n, 1, p.Tr) +1 
  }else{
    Tr <- drop(rbinom(n, 1, 0.5) + 1)
  }
  
  # X main effect
  main.effect  <-  m(drop(X %*% true.eta), delta)
  
  # Tr-by-X interaction effect
  TIE <- g(drop(X %*% true.alpha), w)
  interaction.effect <- (-1)^Tr* TIE
  
  # noise
  noise <- sigma * rnorm(n)
  
  # response
  Ey <- main.effect + interaction.effect
  y  <- Ey   + noise
  
  ## some information about the dataset
  Y1 <- -TIE;
  Y2 <- TIE;
  # the signal to noise ratio  (SNR)
  SNR <- var(interaction.effect)/(var(main.effect) + var(noise))
  SNR
  var.IE <- var(interaction.effect)
  var.ME <- var(main.effect)
  var.noise <- var(noise)
  optTr <- as.numeric(Y2 > Y1) + 1  # this takes 1 or 2
  value.opt <- mean(Ey[Tr == optTr ])
  value.opt
  
  results <- list(y=y, Tr=Tr, X=X, main.effect=main.effect, SNR=SNR, true.alpha = true.alpha, true.eta = true.eta,
                  var.IE=var.IE, var.ME=var.ME, var.noise=var.noise,
                  optTr = optTr, value.opt = value.opt, n=n, p=p, delta = delta, w=w)
  
  return(results)
}



#####################################################################
# A performance mesure function;
# performance_measure assesses the performance of any given treatment selection rule,
# in terms of the "Value" of the treatment decision rule
# and the proportion of correct decisions (PCD).
#####################################################################
performance_measure <- function(pred.test, # predicted value for the testing data
                                y, Tr, X, # these are the testing data
                                value.opt= NULL, optTr =NULL)
{
  
  Tr.coded <- as.numeric(as.factor(Tr))
  
  if(is.null(value.opt))  value.opt <- 1
  
  #n <- length(y)
  ## Prediction error estimation
  #predicted <- numeric(n)
  #for(i in 1:n)   predicted[i] <- pred.test[i, Tr[i]]
  #MSE <-  mean((y - predicted)^2)
  #MAD <- mean(abs(y - predicted))
  #MSE.sd <- sd((y - predicted)^2) /sqrt(n)
  
  # Value estimation
  regime <- apply(pred.test, 1, which.max)
  right.trt.index <- regime == Tr.coded
  
  value <- sum(y[right.trt.index])/sum(right.trt.index)
  value.s <- value/value.opt
  
  #percent correct decision
  if(is.null(optTr)) optTr <-  numeric(length(y))
  pcd <- mean(regime == optTr)
  
  results <- list(#MSE=MSE, MAD=MAD, 
    value =value, value.s= value.s, pcd = pcd)
  return(results)
}



######################################################################
## END OF THE FILE
######################################################################
