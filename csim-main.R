# High dimensional constrained single index model (HD-CSIM) for interactions
# the main functions 
source('csim-other-functions.R');
source('csim-block-cd.R');

library(ggplot2)
library(plyr)
library(splines)
library(boot)
library(glmnet)
library(MASS)
library(magic)  # for adiag()
library(mgcv)

######################################################################
## the wrapper function for fitting the constrained single index model (with variable selection)
## y              -- treatment outcomes, n-by-1 vector
## Tr             -- treatment indicators, n-by-1 vector; each element represents one of the K available treatment options
## X              -- pre-treatment covarate matrix, n-by-p matrix
## ortho.constr   -- the constraint that separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect.
## sparse         -- if TRUE, apply L1 regularization when estimating the single index coefficients.
## type           -- can choose bewteen "CV", "AIC", and "BIC" for the sparsity tuninig parameter selection.
## seed           -- when type="CV", randomization seed for cross validation.
## nbasis.t       -- a length K+1 vector of the numbers of B-spline basis funtions to be used for each treatment group; the last element is for the "main effect" model.
## rho.grid       -- a grid vector of (ridge-type) smoothing parameters for the link functions
## lam.by         -- a value specifying the grid of the sparsity tuning parameters [1, 1+lam.by, ... 1 + n.lam*lam.by].
## n.lam          -- a value specifying the grid of the sparsity tuning parameters [1, 1+lam.by, ... 1 + n.lam*lam.by].
## n.max          -- a maximum number of nonzero (active) coeffients in C-SIM.
## eps            -- a value specifying the converge criterion of algorithm.
## it.max         -- an integer value specifying the maximum number of iterations for each coordinate.
## trace           -- if TRUE, return tracee of the fitting procedure; the default is FALSE.
## eff.aug        -- if TRUE, perform efficiency augmentation (using a L1 regularized linear regression for the main effects of X); the default is FALSE.
## X.aug          -- a design matrix to be used for efficinecy augmentation; the default is NULL.
## coef.ini       -- an initial solution for alpha.coef; the default is NULL.
## i.fx           -- an index to be fixed throughout estimation for model identifiability; the default is NULL.
## mc.ini         -- if TRUE, use an estimate from the modified covariate (MC) approach as an initial solution for alpha.coef; only applicable for K=2 case; the default is FALSE.
## linear.link.fn -- if TRUE, restrict the link functions to be linear; the default is FALSE.
## boot.ci        -- if TRUE, return a boot object that can be used to construct a bootstrap  confidence intervals for the single index coefficients; the default is FALSE.
## n.boot         -- when boot.ci=TRUE, a value specifying the number of bootstrap replications.
######################################################################
csim  <- function(y, Tr, X,
                  ortho.constr = TRUE,
                  sparse = TRUE,
                  type = "AIC", seed = 1234,
                  lam.by = 0.03, n.lam=100,
                  n.max=10, eps = 10^-4, it.max= 70, trace = FALSE,
                  nbasis.t = NULL, rho.grid = c(0, 0.25, 0.5),
                  eff.aug = FALSE, X.aug = NULL,
                  linear.link.fn = FALSE,
                  coef.ini = NULL,
                  i.fx= NULL, mc.ini = FALSE,
                  plots = TRUE, boot.CI = FALSE, n.boot=500)
{


  ################
  ## estimation ##
  ################
  # estimate the single-index coefficients of the model; the sparsity tuning parameter of the L1 regularization is chosen by a 5 fold cross-validation for prediction errors
  coef.obj   <- fit.csim.cv(y, Tr, X,
                           ortho.constr=ortho.constr, sparse=sparse,
                           type = type, seed=seed, lam.by = lam.by, n.lam=n.lam,
                           n.max = n.max, eps = eps, it.max=it.max, trace=trace,
                           nbasis.t = nbasis.t, rho.grid = rho.grid, eff.aug = eff.aug, X.aug=X.aug,
                           linear.link.fn = linear.link.fn, coef.ini= coef.ini, i.fx=i.fx, mc.ini=mc.ini)

  # estimate the link functions using B-splines, given the single-index coefficients.
  link.fn.obj <- fit.link.fn.gcv(coef.obj$yc, coef.obj$Tr, u= drop(coef.obj$Xc %*% coef.obj$coef.opt),
                       nbasis.t = coef.obj$nbasis.t, rho.grid = rho.grid, linear.link.fn = linear.link.fn, ortho.constr = ortho.constr)

  K <- length(unique(Tr))
  t.ind  <- unlist(lapply(1:K, function(x) rep(x, link.fn.obj$smoother$nbasis.t[-(K+1)][x])))
  beta.t.coef <- split(link.fn.obj$beta.t.coef, t.ind)
  beta.0.coef <- link.fn.obj$beta.0.coef


  ###################
  ## visualization ##
  ###################
  # plots
  link.fn.plot = coef.plot <- NULL
  if(plots)
  {
    # 1) link function (gt) plots
    y.means <- rep(0, length(y));
    for(t in 1:K)  y.means  <- y.means +  coef.obj$intercept.y[[t]]*(coef.obj$Tr==t);

    dat  <- data.frame(y =  y.means + coef.obj$yc- link.fn.obj$smoother$B0 %*% beta.0.coef,
                       x = link.fn.obj$smoother$u.t[[K+1]],
                       #Treatment = factor(Tr, labels=c("Placebo","Active drug")) )
                       Treatment = factor(coef.obj$Tr));

    link.fn.plot  <- ggplot(dat, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
      geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
      scale_colour_brewer(palette = "Set1", direction = -1)  +
      theme( axis.title.x=element_text(size=15,face="bold"))  +
      theme(title =element_text(size=12)) +
      xlab(expression(paste(alpha*minute, "x"))) +  #ylab("Improvement") +
      ylab("(Adjusted) y") + theme_bw(base_size = 14);
    if(linear.link.fn)
    {
      link.fn.plot = link.fn.plot + geom_smooth(method=lm, se=TRUE, fullrange=FALSE, alpha = 0.35)
    }else{
      tmp1 <- 0;  for(t in 1:K) tmp1 <- tmp1 + link.fn.obj$smoother$nbasis.t[t];
      link.fn.plot <- link.fn.plot + geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= floor(tmp1/K), sp= link.fn.obj$rho.opt), se=TRUE, fullrange=F, alpha = 0.35)
      link.fn.plot
    }

    # 2) single-index coefficient (alpha.coef) plot
    dotplot_identity = function(frame, xvar, yvar, colour, colorvar=NULL)
    {
      if(is.null(colorvar))  gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar))
      else  gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar, color=colorvar))
      return(gplot + geom_point(colour = colour) + geom_linerange(aes(ymin=0)))
    }
    coef.plot <- dotplot_identity(frame = data.frame(Z = factor(1:ncol(X)), coef = coef.obj$coef.opt), xvar="Z", yvar="coef", colour = "blue") + scale_y_continuous(name = "coef.")  + xlab(expression(alpha))

  }


  ##################################
  ## Bootstrap CIs for alpha.coef ##
  ##################################
  boot.results <- NULL
  if(boot.CI)
  {
    datax <- data.frame(y=y, Tr=Tr, X=X)

    # define the "bs" function for bootstraping
    bs <- function(datax, indices,
                   ortho.constr=ortho.constr, sparse= sparse,
                   type = type, seed = seed, lam.by = lam.by, n.lam=n.lam,
                   n.max = n.max, eps = eps, it.max=it.max,
                   nbasis.t = nbasis.t, rho.grid = rho.grid, eff.aug = eff.aug, linear.link.fn=linear.link.fn, i.fx=i.fx)
    {
      coef.obj.b <- fit.csim.cv(y= datax[indices,][,1], Tr= datax[indices,][,2], X = datax[indices,][,-c(1,2)],
                               ortho.constr = ortho.constr, sparse=sparse,
                               type = type, seed = seed, lam.by = lam.by, n.lam=n.lam,
                               n.max = n.max, eps = eps, it.max=it.max,
                               nbasis.t = nbasis.t, rho.grid = rho.grid, eff.aug = eff.aug, linear.link.fn=linear.link.fn, i.fx=i.fx)
      return(coef.obj.b$coef.opt)
    }

    boot.results <- boot(data = datax, statistic = bs, R = n.boot,
                         ortho.constr=ortho.constr, sparse = sparse, type = type, seed = seed, lam.by = lam.by, n.lam=n.lam,
                         n.max = n.max, eps = eps, it.max=it.max,
                         nbasis.t = nbasis.t, rho.grid = rho.grid, eff.aug = eff.aug, linear.link.fn=linear.link.fn, i.fx = coef.obj$i.fx )
    #boot.ci(boot.results, type = "norm", index = 1)
  }


  results <-  list(coef.obj = coef.obj, link.fn.obj=link.fn.obj,
                   alpha.coef = coef.obj$coef.opt,
                   alpha.coef.pst = coef.obj$coef.pst.opt,
                   sol.path = coef.obj$coef,
                   beta.t.coef = beta.t.coef,
                   beta.0.coef = beta.0.coef,
                   eta.coef = coef.obj$eta.coef,
                   y.hat = link.fn.obj$y.hat,
                   working.y.hat = link.fn.obj$working.y.hat,
                   resid = link.fn.obj$resid,
                   working.resid = link.fn.obj$working.resid,
                   smoother = link.fn.obj$smoother,
                   intercept.y = coef.obj$intercept.y,
                   center.X = coef.obj$center.X,
                   scale.X = coef.obj$scale.X,
                   link.fn.plot= link.fn.plot,
                   coef.plot = coef.plot,
                   boot.results = boot.results,
                   p= ncol(X), n= length(y), K=K)

  class(results) <- c("csim", "list")
  return(results)
}




########################################################
# this is a wrapper function for cross-validation;
# fit.csim.cv performs a 5 fold cross validation to choose the sparsity tuning parameter,
# and returns the model with the optimized tuning parameter.
# If type="BIC", then the function returns the model with the lowest value of BIC.
# the core functions adopted from the code from Professor Peter Radchenko (Univ. of Sydney)
########################################################
fit.csim.cv <- function(y, Tr, X,
                        ortho.constr = TRUE,
                        sparse = TRUE,
                        type = "AIC", seed = 1234, lam.by = 0.03, n.lam=100,
                        n.max=10, eps = 10^-4, it.max=50, trace=FALSE,
                        nbasis.t = NULL, rho.grid = c(0, 0.25, 0.5),
                        eff.aug = FALSE, X.aug = NULL,
                        linear.link.fn = FALSE, coef.ini = NULL, i.fx=NULL, mc.ini=FALSE)
{


  if(type=="CV" & sparse)
  {
    err0 = err1 <- rep(0, n.lam);
    n <- length(y)
    set.seed(seed);
    sam <- sample(1:n);
    dat  <- data.frame(y, Tr, X)[sam, ]
    for(k in 1:5)
    {
      ind.test <- sam[(((k-1)*floor(n/5)+1):(k*floor(n/5)))];

      datax.test <- data.preprocess(dat[ind.test,][,1], dat[ind.test,][,2], dat[ind.test,][,-c(1,2)]);
      y.test <- datax.test$y;
      Tr.test <- datax.test$Tr;
      X.test <- datax.test$Xsc

      y.train <- dat[-ind.test,][,1];
      Tr.train <- dat[-ind.test,][,2];
      X.train <- dat[-ind.test,][,-c(1,2)]

      # fit csim (based on the training set) and obtain a seq of the single-index coefficients over the increasing grids of the tuning parameters
      coef.obj.k <- fit.csim(y.train, Tr.train, X.train,
                            ortho.constr=ortho.constr, sparse=sparse,
                            type = type, lam.by = lam.by, n.lam=n.lam,
                            n.max = n.max, eps = eps, it.max=it.max,
                            nbasis.t = nbasis.t, rho.grid = rho.grid,
                            eff.aug = eff.aug, X.aug=X.aug, linear.link.fn=linear.link.fn,
                            trace= trace, coef.ini=coef.ini, i.fx=i.fx,  mc.ini=mc.ini)

      # appropriately center/scale the testing set
      y.means <- rep(0, length(y.test));
      for(t in 1:length(unique(Tr)))   y.means  <- y.means +  coef.obj.k$intercept.y[[t]] *(Tr.test ==t);
      yc.test <- y.test - y.means

      Xc.test <- scale(X.test, center = coef.obj.k$center.X, scale = coef.obj.k$scale.X)

      # apply the fitted models (i.e., the fitted single-index coefficients) to the validation sets and estimate the mean squares errors
      err.k0 <- apply(coef.obj.k$coef.pst, 1, function(h)  csim.mse.test(coef=h,
                                                                       yc= coef.obj.k$yc, Tr= coef.obj.k$Tr, Xc= coef.obj.k$Xc,
                                                                       yc.test = yc.test, Tr.test = Tr.test, Xc.test = Xc.test,
                                                                       nbasis.t =nbasis.t, rho.grid=rho.grid) )
      err.k1 <- apply(coef.obj.k$coef, 1, function(h)  csim.mse.test(coef=h,
                                                                   yc= coef.obj.k$yc, Tr= coef.obj.k$Tr, Xc= coef.obj.k$Xc,
                                                                   yc.test = yc.test, Tr.test = Tr.test, Xc.test = Xc.test,
                                                                   nbasis.t =nbasis.t, rho.grid=rho.grid) )

      if(length(err.k0) < length(err0))
      {
        err0 <- err0[1:length(err.k0)] + err.k0
      }else{
        err0 <- err0 + err.k0[1:length(err0)]
      }

      if(length(err.k1) < length(err1))
      {
        err1 <- err1[1:length(err.k1)] + err.k1
      }else{
        err1 <- err1+err.k1[1:length(err1)]
      }
    }

    i0 <- floor(median(which(err0==min(err0))));
    i1 <- floor(median(which(err1==min(err1))))
    coef.obj  <- fit.csim(y, Tr, X,
                         ortho.constr = ortho.constr, sparse=sparse,
                         type = type, lam.by = lam.by, n.lam= max(i0, i1),
                         n.max = 1000, eps = eps, it.max=it.max,
                         nbasis.t = nbasis.t, rho.grid = rho.grid,
                         eff.aug = eff.aug, X.aug=X.aug, linear.link.fn=linear.link.fn,
                         trace= trace, coef.ini= coef.ini, i.fx=i.fx,  mc.ini=mc.ini)

    coef.opt <- coef.obj$coef[max(nrow(coef.obj$coef),i1), ]
    coef.pst.opt <- coef.obj$coef.pst[max(nrow(coef.obj$coef.pst),i0), ]

  }else{
    coef.obj <- fit.csim(y, Tr, X,
                        ortho.constr = ortho.constr, sparse=sparse,
                        type = type, lam.by = lam.by, n.lam=n.lam,
                        n.max = n.max, eps = eps, it.max=it.max,
                        nbasis.t = nbasis.t, rho.grid = rho.grid,
                        eff.aug = eff.aug, X.aug=X.aug, linear.link.fn=linear.link.fn,
                        trace= trace, coef.ini= coef.ini, i.fx=i.fx,  mc.ini=mc.ini)
    if(sparse)
    {
      plot(coef.obj$lam.seq, coef.obj$MSE, xlab= expression(lambda), ylab = "MSE")
      if(type=="BIC")
      {
        plot(coef.obj$lam.seq, coef.obj$BIC, xlab= expression(lambda), ylab = "BIC")
        i1 <- floor(median(which(coef.obj$BIC==min(coef.obj$BIC))));
        i0 <- floor(median(which(coef.obj$BIC==min(coef.obj$BIC))))
      }else{
        plot(coef.obj$lam.seq, coef.obj$AIC, xlab= expression(lambda), ylab = "AIC")
        i1 <- floor(median(which(coef.obj$AIC==min(coef.obj$AIC))));
        i0 <- floor(median(which(coef.obj$AIC==min(coef.obj$AIC))))
      }
    }else{
      i1 = i0 <- 1
    }
    coef.opt <- coef.obj$coef[i1,];
    coef.pst.opt <- coef.obj$coef.pst[i0,]
  }

  coef.obj$coef.opt <- coef.opt
  coef.obj$coef.pst.opt <- coef.pst.opt

  return(coef.obj)
}


########################################################
# a subfunction used in cross validation
########################################################
csim.mse.test <- function(coef,
                          yc, Tr, Xc, yc.test, Tr.test, Xc.test,
                          nbasis.t= c(6,6,8), rho.grid=0, linear.link.fn=FALSE)
{

  u <- Xc %*% coef;
  link.fn.obj <- fit.link.fn.gcv(yc, Tr, u=u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link.fn = linear.link.fn, ortho.constr = FALSE)

  t2 <- Xc.test %*% coef;
  if( sum(t2< link.fn.obj$smoother$u.min | t2 > link.fn.obj$smoother$u.max) >0)
  {
    rm.ind <- c((1:length(yc.test))[t2< link.fn.obj$smoother$u.min], (1:length(yc.test))[t2> link.fn.obj$smoother$u.max])
    t2 <- t2[-rm.ind];
    yc.test <- yc.test[-rm.ind];
    Tr.test <- Tr.test[-rm.ind];
    Xc.test <- Xc.test[-rm.ind,];
  }
  smoother.test  <-  smoother.fn(Tr.test, u=t2, nbasis.t= link.fn.obj$smoother$nbasis.t, rho.grid = link.fn.obj$rho.opt, linear.link.fn =linear.link.fn)
  mse <- mean((yc.test - smoother.test$Bt %*% link.fn.obj$beta.t.coef)^2)

  return(mse)
}



########################################################
# this is the workhorse function for fitting the C-SIM.
# It returns the sequence of the model coefficients implied by the tuning paramters lam
# is fit by block coordinate descent algorithm.
########################################################
fit.csim <- function(y, Tr, X,
                     ortho.constr = TRUE,
                     sparse = TRUE,
                     type = "BIC", lam.by = 0.03, n.lam=100,
                     n.max=10, eps = 10^-4, it.max=50, trace=FALSE,
                     nbasis.t = NULL, rho.grid = c(0, 0.25, 0.5),
                     eff.aug = FALSE,
                     X.aug = NULL,   # one can provide a design matrix for the effiency augmentation; if NULL, then the scaled X matrix is taken as X.aug.
                     linear.link.fn = FALSE,
                     coef.ini= NULL,  # an initial solution for alpha.coef.
                     i.fx = NULL,  # i.fx is the index to be fixed throughout the estimation; the default is NULL, hence it will be estimated.
                     mc.ini =FALSE)
{

  ###################
  ## pre-procssing ##
  ###################
  # order the observations by the order of the treatment indicators
  datax <- data.preprocess(y, Tr, Xsc = X);
  y <- datax$y;
  Tr <- datax$Tr;
  X <- datax$Xsc

  p <- ncol(X);
  n <- length(y);
  K <- length(unique(Tr))

  # center and scale X
  Xc <- scale(X, center = TRUE, scale = TRUE)
  scale.X <-  attr(Xc, "scaled:scale")
  center.X <- attr(Xc, "scaled:center")
  # center y within each treatment, and store the treatment-specific mean at intercept.y
  dat <- data.frame(Tr=Tr, y=y)
  dat.list <- dlply(dat, .(dat$Tr), function(x) {as.matrix(x[,-1])});
  intercept.y <- vector("list", K);
  yc <- NULL
  for(t in 1:K)
  {
    dat.list[[t]] <- scale(dat.list[[t]], center = TRUE, scale = FALSE)
    intercept.y[[t]] <- attr(dat.list[[t]], "scaled:center")
    yc <- c(yc, dat.list[[t]])
  }
  rm(dat);
  rm(dat.list)

  if(eff.aug)  # perform efficiency augmentation
  {
    if(is.null(X.aug))  X.aug <- as.matrix(Xc)

    # choose an optimal sparsity tuning parameter
    lasso.temp <- cv.glmnet(x= X.aug, y=yc, nfolds=20)
    lamda.opt <- lasso.temp$lambda[which.min(lasso.temp$cvm)]

    # fit the linear model with lasso
    lasso.fit <- glmnet(x= X.aug, y=yc, lambda =lamda.opt, alpha=1)

    # residualize the outcomes
    yc <- drop( yc - predict(lasso.fit, X.aug) )
    eta.coef <- lasso.fit$beta
  }else{
    eta.coef <- NULL
  }

  if(is.null(nbasis.t))
  {
    nt <- summary(as.factor(Tr))
    for(t in 1:K)
    {
      nbasis.t[t] <- floor(nt[t]^{1/5.5})  + 4
    }
    nbasis.t[K+1] <- floor(n^{1/5.5}) + 6
  }


  # initialize the single-index coefficient vector, coef
  if(is.null(coef.ini))
  {
    coef.ini <- rep(0,p);
    if(mc.ini)
    {
      tmp <-  as.numeric(mc(y, Tr, X, eff.aug = FALSE, use.lasso = sparse)$alpha.coef)
      if(sum(tmp!=0) ==0)
      {
        if(is.null(i.fx))  i.fx <- 1;
        coef.ini[i.fx] <- 1
      }else{
        if(is.null(i.fx))  i.fx <- which.max(abs(tmp));
        coef.ini  <- tmp / tmp[i.fx ]
      }
      rm(tmp)
    }else{
      if(is.null(i.fx))
      {
        tmp <- rep(10^4,p);
        for(j in  1:p)
        {
          tmp[j] <- fit.link.fn.gcv(yc, Tr, u= Xc[,j], nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link.fn = linear.link.fn, ortho.constr =ortho.constr, ini=TRUE)$MSE
        }
        i.fx  <-  order(tmp, decreasing= FALSE)[1]
      }
      coef.ini[i.fx] <- 1;
    }
  }else{
    if(is.null(i.fx))  i.fx  <- which.max(abs(coef.ini));
    coef.ini  <- coef.ini /coef.ini[i.fx]
  }

  coef = coef.pst <- rbind(coef.ini)
  MSE = BIC = AIC= lam.seq <- rep(10^4, n.lam)  # some large number

  if(sparse)
  {
    link.fn.obj <- fit.link.fn.gcv(yc, Tr, u= Xc %*% coef.ini, nbasis.t = nbasis.t, rho.grid = rho.grid,
                         linear.link.fn = linear.link.fn, ortho.constr =ortho.constr)
    d.link.fn.obj <- deriv.link.fn(link.fn.obj)

    thr1 <- abs(apply(Xc, 2, function(v){ crossprod(d.link.fn.obj$d.link.fn1, v)} ));
    thr1[i.fx] <- -10^10;
    # find the one that has the largest absolute gradient
    i.m <- which.max(thr1)

    # increase the L1 norm of solution
    lam <- sum(abs(coef.ini)) + lam.by;
    w1 <- coef.ini
    # i.m is the reference index
    w1[i.m] <- lam.by*sign(-crossprod(d.link.fn.obj$d.link.fn1, Xc[,i.m]));
    # compute the new single-index
    u <- Xc %*% w1

    link.fn.obj <- fit.link.fn.gcv(yc, Tr, u=u, nbasis.t = nbasis.t, rho.grid =  link.fn.obj$rho.opt,
                         linear.link.fn = linear.link.fn, ortho.constr = ortho.constr)
    d.link.fn.obj <- deriv.link.fn(link.fn.obj)

    # thr1 is the absolute value of the gradient  (here, thr1[i.fx] is excluded )
    act1 <- (thr1 >= thr1[i.m]);
    act.old <- act1;
    new.act <- act1;
    fl.new <- TRUE;
    dcoef <- 1;
    it2 <- 0

    lam.seq[dcoef] <- lam;
    MSE[dcoef] <- link.fn.obj$MSE;
    BIC[dcoef] <- log(MSE[dcoef]) + sum(w1 != 0) /n * log(n);
    AIC[dcoef] <- log(MSE[dcoef]) + sum(w1 != 0) /n * 2 + 2*sum(w1 != 0)*(sum(w1 != 0)+1)/(n-sum(w1 != 0)-1) /n
    MSE.old <- MSE[1]

    if(trace) cat("\n\n","act1 = ",(1:p)[act1],"  non-zeros1  = ",(1:p)[w1!=0],"\n\n")

    while((sum(w1!=0) < min(n.max,p)) & (dcoef < n.lam) )
    {
      if(trace) cat("\n","^^^^^^^^^^^",lam,"^^^^^^^^^^^^^^^^^","\n")
      # active indices
      ind1.a <- (1:p)[act1];
      flag <- FALSE
      # coefficient associated with the active indices
      w.a <- w1[ind1.a];
      ind.m <- (1:length(ind1.a))[ind1.a==i.m];
      ind.old <- (1:length(ind1.a))[is.element(ind1.a, (1:p)[act.old])];
      # X associated with the active indices
      X.a <- as.matrix(Xc[,ind1.a])

      obj.cd <- block.cd(w.a=w.a, ind.m=ind.m, ind.old=ind.old, u=u, y=yc, Tr= Tr, X.a=X.a,
                         d.link.fn.obj = d.link.fn.obj, nbasis.t = nbasis.t, rho.grid = link.fn.obj$rho.opt,
                         linear.link.fn= linear.link.fn, ortho.constr = ortho.constr, it.max=it.max, eps=eps, trace= trace)

      w1[ind1.a] <- obj.cd$w.a;
      u <- obj.cd$u;
      it <- obj.cd$it
      w1[abs(w1)< eps] <- 0;
      act1[w1==0] <- FALSE;
      # if all coefficients in the new active set are zero, the flag is on.
      if(sum(w1[new.act]!=0)==0)  flag <- TRUE
      tmp <- rep(0, p);
      tmp[ind1.a] <- abs(w1[ind1.a]);
      i.m <- which.max(tmp);
      rm(tmp)

      link.fn.obj <- fit.link.fn.gcv(yc, Tr, u=u, nbasis.t = nbasis.t, rho.grid = link.fn.obj$rho.opt,
                           linear.link.fn = linear.link.fn, ortho.constr = ortho.constr)
      d.link.fn.obj <- deriv.link.fn(link.fn.obj)

      tmp1 <- apply(Xc, 2, function(v){crossprod(d.link.fn.obj$d.link.fn1, v)});
      tmp2 <- -sign(w1);
      tmp2[tmp2==0] <- sign(tmp1[tmp2==0])
      # updated gradient
      thr1 <- tmp1*tmp2;
      thr1[i.fx] <- -10^10;
      thr.entr <- max(thr1[act1]);
      rm(tmp1);
      rm(tmp2)
      if(sum(!act1)>0)
      {
        if((max(thr1[!act1]) <= thr.entr)|flag)
        {
          if(fl.new & (!flag))
          {
            act1.pst <- act1;
            act1.pst[i.fx] <- TRUE;
            w.pst <- rep(0,p);
            w.pst[act1.pst] <- forw.cd(w1[act1.pst], u, y =yc, Tr= Tr, Xc[,act1.pst],
                                       nbasis.t= nbasis.t, rho.grid = rho.grid, it.max=it.max, eps=eps)$w.a
            coef.pst <- rbind(coef.pst, w.pst);
            #coef.pst <- rbind(coef.pst,w.pst/sqrt(sum(w.pst^2)));
            dcoef <- nrow(coef.pst)
          }
          else{
            coef.pst <- rbind(coef.pst, coef.pst[dim(coef.pst)[1],]);
            dcoef <- nrow(coef.pst)
          }
          new.act <- act1;
          fl.new <- FALSE
          if(trace) cat("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n")
          if(trace) cat("lam = ",lam,"   MSE=", link.fn.obj$MSE,"\n",it,"it"," i.fx=",i.fx," i.m=",i.m," l1=", sum(abs(w1)),"\n","w1=",w1[w1!=0],"\n")
          lam <- lam + lam.by;
          coef <- rbind(coef, w1)
          w1[i.m] <- w1[i.m] + lam.by*sign(w1[i.m])
          u <- u + lam.by*sign(w1[i.m])*Xc[,i.m]

          link.fn.obj <- fit.link.fn.gcv(yc, Tr, u=u, nbasis.t = nbasis.t, rho.grid = rho.grid,
                               linear.link.fn = linear.link.fn, ortho.constr = ortho.constr)
          d.link.fn.obj <- deriv.link.fn(link.fn.obj)

          lam.seq[dcoef] <- lam;
          MSE[dcoef] <- link.fn.obj$MSE;
          BIC[dcoef] <- log(MSE[dcoef]) + sum(w1 != 0) /n * log(n)
          AIC[dcoef] <- log(MSE[dcoef]) + sum(w1 != 0) /n * 2 + 2*sum(w1 != 0)*(sum(w1 != 0)+1)/(n-sum(w1 != 0)-1) /n
          if(MSE.old*(1+ log(p)*eps) < MSE[dcoef])  break
          MSE.old <- MSE[dcoef]
        }
        else{
          if(trace) cat(it,"it"," i.fx=",i.fx," i.m=",i.m," l1=",sum(abs(w1)),"\n","w1=",w1[w1!=0],"\n")
          act.old <- act1;
          new.act <- (1:p)[thr1 > thr.entr];
          fl.new <- TRUE;
          act1[thr1 > thr.entr] <- TRUE;
          if(trace) cat("new act =   ",new.act,"\n")
          it2 <- it2 + 1;
          if(it2 > it.max) break
        }
      }
      if(trace) cat("non-zeros= ",(1:p)[w1!=0],"\n")
      if(trace) cat("act= ",(1:p)[act1],"\n")
      if(trace) cat("thresh= ",thr1,"\n")
      if(sum(act1)==0) break
    }
  }else{
    # u is an initial single index variable
    w.pst <- forw.cd(coef, u= Xc %*% coef.ini, y=yc, Tr= Tr, X.a =Xc, linear.link.fn=linear.link.fn, ortho.constr = ortho.constr,
                     nbasis.t= nbasis.t, rho.grid = rho.grid, it.max=it.max, eps=eps, trace= trace, i.fx = i.fx)$w.a
    coef = coef.pst <- w.pst
  }


  return(list(coef=coef, coef.pst=coef.pst, eta.coef=eta.coef,
              intercept.y = intercept.y, scale.X = scale.X, center.X = center.X,
              yc =yc, Tr= Tr, Xc = Xc, nbasis.t = nbasis.t, rho.grid=rho.grid,
              BIC =BIC[BIC<10^4], AIC = AIC[AIC<10^4], MSE=MSE[MSE<10^4], lam.seq = lam.seq[lam.seq<10^4], i.fx=i.fx) )
}


########################################################
# This function makes predictions from the constrained single index model,
# given a csim object and pretreatment covariates.
# The function returns predicted outcomes for each treatment and treatment selection rules.
########################################################
predict.csim  <-  function(csim.obj, newx)
{
  if(!inherits(csim.obj, "csim"))   # checks input
    stop("obj must be of class `csim'")

  n <- nrow(newx)
  K <- csim.obj$K
  if(ncol(newx) != csim.obj$p)
    stop("newx needs to be of p columns ")

  alpha.coef <- csim.obj$alpha.coef
  beta.t.coef <- csim.obj$beta.t.coef
  link.fn.obj <- csim.obj$link.fn.obj

  newx.scaled <- scale(newx, center = csim.obj$center.X, scale = csim.obj$scale.X);
  Index  <- newx.scaled %*% alpha.coef;
  Index[Index < link.fn.obj$smoother$u.min] <- link.fn.obj$smoother$u.min;
  Index[Index > link.fn.obj$smoother$u.max] <- link.fn.obj$smoother$u.max

  # compute treatment-specific predicted value
  pred.new <- matrix(0, n, K)
  for(t in 1:K)
  {
    if(link.fn.obj$smoother$linear.link.fn)
    {
      pred.new[, t]  <- cbind(1, Index) %*% beta.t.coef[[t]]
    }else{
      temp <- splineDesign(link.fn.obj$smoother$knots.t[[t]], x= Index, outer.ok = TRUE);
      pred.new[, t] <- csim.obj$intercept.y[[t]] +  temp  %*% beta.t.coef[[t]];
      rm(temp)
    }
  }

  # compute treatment assignment
  trt.rule <- apply(pred.new, 1, which.max)  # WLLG, we assume a higher y is prefered.
  if(K==2)  colnames(pred.new) <- c("Tr1", "Tr2")
  results <- list(trt.rule = trt.rule, pred.new = pred.new)

  return(results)
}



########################################################
# A pre-processing  function that can be  used to create a data frame that consists of the observations
# in the order of the treatment indicators;
# i.e., (y, Tr, X) with Tr==1 come in the first rows of the data frame,
# (y, Tr, X) with Tr==2 second, ..., and (y, Tr, X) with Tr==K come last.
########################################################
data.preprocess <- function(y=NULL, Tr,  Xsc=NULL, X1D=NULL, optTr=NULL, main.effect= NULL)
{

  K <- length(unique(Tr))
  n <- length(y)
  p <- length(X1D)

  if(!is.null(Xsc))
  {
    q <- dim(Xsc)[2]
    dat.temp <- data.frame(Tr=Tr, X = Xsc)   # organize your data in a dataframe
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
    Xsc.new <- NULL
    for(t in 1:K) Xsc.new <- rbind(Xsc.new, dat.temp.list[[t]])
    Xsc <- Xsc.new
  }

  if(!is.null(X1D))
  {
    for(j in 1:p){
      X1D.j <- X1D[[j]]
      dat.temp <- data.frame(Tr=Tr, X = X1D.j)   # organize your data in a dataframe
      dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
      X1D.new <- NULL
      for(t in 1:K) X1D.new <- rbind(X1D.new, dat.temp.list[[t]])
      X1D[[j]] <- X1D.new
    }
  }

  if(!is.null(main.effect))
  {
    dat.temp <- data.frame(Tr=Tr, main.effect=main.effect)   # organize your data in a dataframe
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x)  as.matrix(x[,-1]) )
    main.effect.new <- NULL
    for(t in 1:K) main.effect.new <- rbind(main.effect.new, dat.temp.list[[t]])
    main.effect <- main.effect.new
  }

  if(!is.null(y))
  {
    dat.temp <- data.frame(Tr=Tr, y = y)   # organize your data in a dataframe
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x)  as.matrix(x[,-1]) )
    y.new <- NULL
    for(t in 1:K) y.new <- rbind(y.new, dat.temp.list[[t]])
    y <- y.new
  }

  if(!is.null(optTr))  # this is only avaialbe for a simulated dataset where we know the true model
  {
    optTr <- optTr
    dat.temp <- data.frame(Tr=Tr, optTr)
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
    optTr.new <- NULL
    for(t in 1:K) optTr.new <- rbind(optTr.new, dat.temp.list[[t]])
    optTr <- optTr.new  # this is in (-1, 1) values
    optTr <- 0.5*optTr +1.5  # this is in (1, 2) values
  }

  dat.temp <- data.frame(Tr=Tr, Tr )   # organize your data in a dataframe
  dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
  Tr.new <- NULL
  for(t in 1:K) Tr.new <- rbind(Tr.new, dat.temp.list[[t]])
  Tr <- as.vector(Tr.new)

  return(list(y=y, Tr=Tr, Xsc=Xsc, X1D=X1D, optTr=optTr, main.effect=main.effect))
}



######################################################################
## END OF THE FILE
######################################################################
