# CSIM paper example Section 5.1 (K=2 treatment case)
# variable selection (VS) performance

source('csim-main.R');

# higher dimensional setting
HD.Simulation.K2  = function(n.rep = 1,
                             sim.seed = 1234,
                             n = 200, p=50, sigma = 0.4,
                             correlationX= 0, sigmaX = 1,
                             n.test = 10000,
                             w = 1, delta = 1,
                             n.max=10, n.lam=100, trace = FALSE,
                             nbasis.t = c(6,6,8), rho.grid = c(0,0.25,0.5), type = "AIC")
{

  true.alpha <- c(c(1, 0.5, 0.25, 0.125), rep(0, p-4))

  g <- function(u, w)
  {
    if(w==1) return(0.5* u)
    if(w==2) return(sin(u) - u)
    if(w==3) return(cos(u) - 0.5)
  }

  m <- function(u, delta= 1)   2+ delta*cos(u*0.5*pi)

  method.names <- c("C-SIM", "MCA", "K-LR")
  n.col <- length(method.names)

  # model size recorder
  Model.Size.Recorder <- matrix(NA, nrow=n.rep, ncol= n.col )
  colnames(Model.Size.Recorder) <- method.names
  # C and IC recorder
  # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
  Model.Sel.Recorder.alpha.C <- matrix(NA, nrow=n.rep, ncol=n.col )
  colnames(Model.Sel.Recorder.alpha.C) <- method.names
  # the number of noise predictors incorrectly identified as signal.
  Model.Sel.Recorder.IC <- matrix(NA, nrow=n.rep, ncol=n.col )
  colnames(Model.Sel.Recorder.IC) <- method.names

  Value.Recorder  <- matrix(NA, nrow=n.rep, ncol=n.col)
  colnames(Value.Recorder) <- method.names
  Per.correct.Recorder  <- matrix(NA, nrow=n.rep, ncol=n.col)
  colnames(Per.correct.Recorder) <- method.names


  results.aggregated <- NULL
  for(sim_run in 1:n.rep)
  {
    print(sim_run)
    set.seed(sim.seed + sim_run)

    ##################################
    #### testing data generation  ####
    ##################################

    Tr.test <- drop(rbinom(n.test, 1, 0.5) + 1)
    Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
    ePsix <- eigen(Psix)
    X.test <- sapply(1:p, function(x) rnorm(n.test)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)

    eta.hold <- rnorm(12, 0, 1);
    eta.hold  <- eta.hold /sqrt(sum(eta.hold^2) ) #* delta
    eta.hold
    eta <- c(eta.hold, rep(0, p-12))

    main.effect  <-  m(drop(X.test %*% eta), delta)

    TIE <- g(drop(X.test %*% true.alpha), w)
    interaction.effect <- (-1)^Tr.test* TIE
    noise <- sigma * rnorm(n.test)

    Ey <- main.effect + interaction.effect
    y.test <- Ey + noise
    Y1 <- -TIE
    Y2 <- TIE
    SNR <- var(interaction.effect)/(var(main.effect) + var(noise))
    SNR
    var(interaction.effect)
    var(main.effect)
    var(noise)
    optTr <- as.numeric(Y2 > Y1) + 1  # this takes 1 or 2
    value.opt <- mean(Ey[Tr.test  == optTr ])
    value.opt


    ##################################
    #### training data generation ####
    ##################################
    Tr <- drop(rbinom(n, 1, 0.5) + 1)
    Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
    ePsix <- eigen(Psix)
    X <- sapply(1:p, function(x) rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)

    main.effect  <-  m(drop(X %*% eta), delta)

    TIE <- g(drop(X %*% true.alpha), w)

    interaction.effect <- (-1)^Tr* TIE
    noise <- sigma * rnorm(n)
    y <-  main.effect + interaction.effect  + noise


    ########
    ########

    #1. C-SIM
    csim.object  <- csim(y, Tr, X,
                         type = type,  # tuning parameter selection
                         lam.by = 0.03,
                         eps = 10^-4, n.max=n.max, n.lam=n.lam, it.max=50, trace=trace,
                         nbasis.t = nbasis.t,  # a vector of the numbers of basis funtions to be used for each treatment group
                         rho.grid = rho.grid, # grid of (ridge-type) smoothing parameters for the links
                         eff.aug = FALSE,  # If TRUE, use a linear model as a main effect working model
                         linear.link = FALSE,
                         plots = FALSE)

    csim.object$alpha.coef
    #  y.test, Tr.test, X.test are testing data
    pred.test <- predict.csim(csim.object, X.test)$pred.new
    performance.object <- performance_measure(pred.test, y.test, Tr.test, X.test, value.opt=value.opt, optTr=optTr)
    performance.object$pcd
    performance.object$value.s

    csim.res <- c(csim.pcd = performance.object$pcd, csim.value = performance.object$value.s)
    csim.res

    contrast <- csim.object$alpha.coef
    Model.Size.Recorder[sim_run,1] <- sum( contrast !=0)
    #C and IC recorder
    # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
    Model.Sel.Recorder.alpha.C[sim_run,1] <- sum(contrast[true.alpha!=0] != 0)
    # the number of noise predictors incorrectly identified as signal.
    Model.Sel.Recorder.IC[sim_run,1] <- sum(contrast[ ! (true.alpha!=0) ] != 0)

    Value.Recorder[sim_run,1] <- performance.object$value.s
    Per.correct.Recorder[sim_run,1] <- performance.object$pcd


    # 2. mc
    mc.object <- mc(y, Tr, X, eff.aug  = TRUE, use.lasso=TRUE)
    pred.test <- predict.mc(mc.object, X.test)$pred.new
    performance.object <- performance_measure(pred.test, y.test, Tr.test, X.test, value.opt=value.opt, optTr=optTr)
    performance.object$pcd
    performance.object$value.s
    mc.res <- c( mc.pcd = performance.object$pcd,  mc.value = performance.object$value.s)
    mc.res

    contrast <- mc.object$alpha.coef
    Model.Size.Recorder[sim_run,2] <- sum( contrast !=0)
    #C and IC recorder
    # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
    Model.Sel.Recorder.alpha.C[sim_run,2] <- sum(contrast[true.alpha!=0] != 0)
    # the number of noise predictors incorrectly identified as signal.
    Model.Sel.Recorder.IC[sim_run,2] <- sum(contrast[ ! (true.alpha!=0) ] != 0)

    Value.Recorder[sim_run,2] <- performance.object$value.s
    Per.correct.Recorder[sim_run,2] <- performance.object$pcd



    # 3. system of K linear regression with lasso
    K.LR.object <- K.LR(y, Tr, X)
    pred.test <- predict(K.LR.object, X.test)$pred.new
    performance.object <- performance_measure(pred.test, y.test, Tr.test, X.test, value.opt=value.opt, optTr=optTr)
    performance.object$pcd
    performance.object$value.s
    K.LR.res <- c( K.LR.pcd = performance.object$pcd,  K.LR.value = performance.object$value.s)
    K.LR.res

    contrast = as.numeric(coef(K.LR.object[[2]] ))[-1] -  as.numeric(coef(K.LR.object[[1]] ))[-1]
    Model.Size.Recorder[sim_run,3] <- sum( contrast !=0)
    #C and IC recorder
    # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
    Model.Sel.Recorder.alpha.C[sim_run,3] <- sum(contrast[true.alpha!=0] != 0)
    # the number of noise predictors incorrectly identified as signal.
    Model.Sel.Recorder.IC[sim_run,3] <- sum(contrast[ ! (true.alpha!=0) ] != 0)
    Value.Recorder[sim_run,3] <- performance.object$value.s
    Per.correct.Recorder[sim_run,3] <- performance.object$pcd


    results.ii <- c(csim.res,
                    mc.res,
                    K.LR.res)

    results.aggregated <- rbind(results.aggregated, results.ii)
    # print(sim_run)
  } # the end of the for loop


  results.ii

  results.mean <- apply(results.aggregated, 2, mean)
  results.sd <- apply(results.aggregated, 2, sd)

  mean.Model.Size.Recorder <- colMeans(Model.Size.Recorder)
  mean.Model.Sel.Recorder.IC <- colMeans(Model.Sel.Recorder.IC)
  mean.Model.Sel.Recorder.alpha.C <- colMeans(Model.Sel.Recorder.alpha.C)
  mean.Per.correct.Recorder <-colMeans(Per.correct.Recorder)
  mean.Value.Recorder <-colMeans(Value.Recorder)

  var.Model.Size.Recorder <- apply(Model.Size.Recorder, 2, var)
  var.Value.Recorder <- apply(Value.Recorder, 2,var)
  var.Model.Sel.Recorder.IC <- apply(Model.Sel.Recorder.IC, 2,var)
  var.Model.Sel.Recorder.alpha.C <- apply(Model.Sel.Recorder.alpha.C , 2,var)
  var.Per.correct.Recorder <- apply(Per.correct.Recorder , 2,var)

  final<- NULL
  final <- rbind(mean.Value.Recorder, var.Value.Recorder,
                 mean.Per.correct.Recorder, var.Per.correct.Recorder,
                 mean.Model.Size.Recorder, var.Model.Size.Recorder,
                 mean.Model.Sel.Recorder.alpha.C, var.Model.Sel.Recorder.alpha.C,
                 mean.Model.Sel.Recorder.IC, var.Model.Sel.Recorder.IC)
  print(final)

  return(list(final = final, SNR = SNR,
              results.aggregated = results.aggregated, results.mean = results.mean, results.sd =results.sd,
              Value.Recorder=Value.Recorder, Per.correct.Recorder=Per.correct.Recorder,
              Model.Size.Recorder=Model.Size.Recorder,
              Model.Sel.Recorder.alpha.C=Model.Sel.Recorder.alpha.C,
              Model.Sel.Recorder.IC= Model.Sel.Recorder.IC))
}


# variable selection performance

scenarios <- expand.grid(n = c(100,
                               200,
                               300,
                               400,
                               500,
                               600,
                               700,
                               800,
                               900,
                               1000),
                         p = c(50, 500),
                         delta = c(1, 2),
                         w = c(1, 3))
scenarios
n.rep <- 2 # 200
sim.seed <- 1234

results.HD.K2  <- list()
for(r in 1:nrow(scenarios))
{
  n = scenarios[r, 1]   # the number of observations
  p = scenarios[r, 2]   # the number of covariates
  delta = scenarios[r, 3]   # intensity of main effect; delta=1 moderate, delta=2 large
  w = scenarios[r, 4]   # degree of nonlinearity; w=1 linear, w=3 highly nonlinear (cosine)

  results.HD.K2.r  <- HD.Simulation.K2(n.rep = n.rep, sim.seed= sim.seed,
                                     n = n, p=p, w=w, delta = delta,
                                     nbasis.t = NULL, rho.grid = 0,
                                     n.max = 10, n.lam =100, trace = F,
                                     sigma = 0.4, correlationX= 0, n.test = 10000,
                                     sigmaX = 1, type = "AIC")
  print(results.HD.K2.r$final)
  print(r)

  results.HD.K2[[r]] <- results.HD.K2.r$final
}

results.HD.K2
 
#save.file <- paste("csim-vs-K2-", r, ".RData", sep="")
#save.image(save.file)




######################################################################
## END OF THE FILE
######################################################################
