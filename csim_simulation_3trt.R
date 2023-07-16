# CSIM paper example Section 5.2 (K=3 treatment case)
source('csim-main.R')


# higher dimensional setting
HD.Simulation.K3  = function(n.rep = 1, sim.seed = 1234,
                             n = 200, p=50, sigma = 0.4,
                             correlationX= 0, sigmaX = 1,
                             n.test = 10000,
                             delta = 1, trace = FALSE,
                             nbasis.t = c(6,6,6,8),
                             rho.grid = c(0,0.25,0.5), type = "AIC")
{

  true.alpha <- c(c(1, 0.5, 0.25, 0.125), rep(0, p-4))

  # define the treamtent-specific link functions
  g1.tmp = function(x) x^(1)*(1-x)^(4) /(gamma(2)*gamma(5)/gamma(7))
  g2.tmp = function(x) x^(1)*(1-x)^(1) /(gamma(2)*gamma(2)/gamma(4))
  g3.tmp = function(x) x^(4)*(1-x)^(0) /(gamma(5)*gamma(1)/gamma(6))
  g0 <- function(x) (g1.tmp(x) + g2.tmp(x) + g3.tmp(x))/3

  g1 <- function(x)  g1.tmp(x) - g0(x)
  g2 <- function(x)  g2.tmp(x) - g0(x)
  g3 <- function(x)  g3.tmp(x) - g0(x)

  # define the main effect functions
  m <- function(u, delta= 1)   delta*cos(u*0.5*pi)


  method.names <- c("C-SIM", "K-LR", "K-SAM")
  n.col <- length(method.names)

  # model size recorder
  Model.Size.Recorder <- matrix(NA, nrow=n.rep, ncol= n.col )
  colnames(Model.Size.Recorder)<- method.names
  # C and IC recorder
  # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
  Model.Sel.Recorder.alpha.C <- matrix(NA, nrow=n.rep, ncol=n.col )
  colnames(Model.Sel.Recorder.alpha.C)<- method.names
  # the number of noise predictors incorrectly identified as signal.
  Model.Sel.Recorder.IC <- matrix(NA, nrow=n.rep, ncol=n.col )
  colnames(Model.Sel.Recorder.IC)<- method.names
  Value.Recorder  <- matrix(NA, nrow=n.rep, ncol=n.col )
  colnames(Value.Recorder )<- method.names
  Per.correct.Recorder  <- matrix(NA, nrow=n.rep, ncol=n.col )
  colnames(Per.correct.Recorder)<- method.names


  results.aggregated <- NULL
  for(sim_run in 1:n.rep)
  {
    print(sim_run)
    set.seed(sim.seed + sim_run)

    ##################################
    #### testing data generation  ####
    ##################################
    Tr.test <- sample(c(1,2,3), n.test, replace=TRUE)
    Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
    ePsix <- eigen(Psix)
    X.test <- sapply(1:p, function(x) rnorm(n.test)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)

    eta.hold <- rnorm(12, 0, 1);
    eta.hold  <- eta.hold /sqrt(sum(eta.hold^2) ) #* delta
    eta.hold
    eta <- c(eta.hold, rep(0, p-12))

    main.effect  <-  m(drop(X.test %*% eta), delta)

    # interaction effect
    z.test <- X.test %*% true.alpha
    index.range <- max(range(abs(z.test)))
    IntegralTransform  <- function(index.range){
      w <- 5;  integrand <- function(t)  gamma(w+1)/(gamma((w+1)/2)^2 * 2^w) *(1-t^2)^((w-1)/2)
      function(v) stats::integrate(integrand, lower = -1.01, upper = v/index.range)$value
    }
    u.test <- sapply(z.test, IntegralTransform(index.range) )
    interaction.effect <- g1(u.test)*(Tr.test ==1) + g2(u.test)*(Tr.test ==2) + g3(u.test)*(Tr.test ==3)

    # noise
    noise <- sigma * rnorm(n.test)

    Ey <- main.effect + interaction.effect
    y.test <- Ey + noise

    # obtain an optimal ITR
    Y1 <-  g1(u.test)
    Y2 <-  g2(u.test)
    Y3 <-  g3(u.test)
    optTr <- apply(cbind(Y1, Y2, Y3), 1, which.max)


    SNR <- var(interaction.effect)/(var(main.effect) + var(noise))
    SNR
    value.opt <- mean(Ey[Tr.test  == optTr])
    value.opt


    ##################################
    #### training data generation ####
    ##################################
    Tr <- sample(c(1,2,3), n, replace=TRUE)
    Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
    ePsix <- eigen(Psix)
    X <- sapply(1:p, function(x) rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)

    # main effect
    main.effect  <-  m(drop(X %*% eta), delta)

    # interaction effect
    z <- X %*% true.alpha
    index.range <- max(range(abs(z)))
    u <- sapply(z, IntegralTransform(index.range) )
    interaction.effect <- g1(u)*(Tr ==1) + g2(u)*(Tr ==2) + g3(u)*(Tr ==3)

    # noise
    noise <- sigma*rnorm(n)

    Ey <- main.effect + interaction.effect
    y <- Ey + noise


    #1. C-SIM
    csim.object  <- csim(y, Tr, X,
                         type = type,  # tuning parameter selection
                         lam.by = 0.03,
                         eps = 10^-4, it.max=50, trace=trace,
                         nbasis.t = nbasis.t,  # a vector of the numbers of basis funtions to be used for each treatment group
                         rho.grid = rho.grid,   # grid of (ridge-type) smoothing parameters for the links
                         eff.aug = FALSE,  # If TRUE, use a linear model as a main effect working model
                         linear.link = FALSE,
                         plots = FALSE)

    csim.object$alpha.coef
    # say, y.test, Tr.test, X.test are testing data
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



    # 2. system of K linear regression with lasso
    K.LR.object <- K.LR(y, Tr, X)
    pred.test <- prediction.object <- predict(K.LR.object, X.test)$pred.new
    performance.object <- performance_measure(pred.test, y.test, Tr.test, X.test, value.opt=value.opt, optTr=optTr)
    performance.object$pcd
    performance.object$value.s
    K.LR.res <- c( K.LR.pcd = performance.object$pcd,  K.LR.value = performance.object$value.s)
    K.LR.res


    contrast = as.numeric(coef(K.LR.object[[2]] ))[-1] -  as.numeric(coef(K.LR.object[[1]] ))[-1]
    Model.Size.Recorder[sim_run,2] <- sum( contrast !=0)
    #C and IC recorder
    # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
    Model.Sel.Recorder.alpha.C[sim_run,2] <- sum(contrast[true.alpha!=0] != 0)
    # the number of noise predictors incorrectly identified as signal.
    Model.Sel.Recorder.IC[sim_run,2] <- sum(contrast[ ! (true.alpha!=0) ] != 0)
    Value.Recorder[sim_run,2] <- performance.object$value.s
    Per.correct.Recorder[sim_run,2] <- performance.object$pcd


    # 3. K separate SAM
    K.SAM.object <- K.SAM(y, Tr, X)
    pred.test <- predict(K.SAM.object, X.test)$pred.new
    performance.object <- performance_measure(pred.test, y.test, Tr.test, X.test, value.opt=value.opt, optTr=optTr)
    performance.object$pcd
    performance.object$value.s
    K.SAM.res <- c( K.SAM.pcd = performance.object$pcd,  K.SAM.value = performance.object$value.s)
    K.SAM.res

    contrast  =  K.SAM.object$K.SAM[[2]]$func_norm -  K.SAM.object$K.SAM[[1]]$func_norm
    Model.Size.Recorder[sim_run,3] <- sum( contrast !=0)
    #C and IC recorder
    # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
    Model.Sel.Recorder.alpha.C[sim_run,3] <- sum(contrast[true.alpha!=0] != 0)
    # the number of noise predictors incorrectly identified as signal.
    Model.Sel.Recorder.IC[sim_run,3] <- sum(contrast[ ! (true.alpha!=0) ] != 0)

    Value.Recorder[sim_run,3] <- performance.object$value.s
    Per.correct.Recorder[sim_run,3] <- performance.object$pcd



    results.ii <- c(csim.res,
                    K.LR.res,
                    K.SAM.res)
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




scenarios <- expand.grid(n = c(250, 500),  # the number of observations
                         p = c(50, 500),  # the number of covariates
                         delta = c(1, 2))   # intensity of main effect; delta=1 moderate, delta=2 large
scenarios
n.rep <- 2 # 200  # the number of simulation replications for each scenario
sim.seed <- 1234

results.HD.K3 = results.aggregated.HD.K3  <- list()
for(r in 1:nrow(scenarios))
{

  n = scenarios[r, 1]
  p = scenarios[r, 2]
  delta = scenarios[r, 3]

  # note that most of the computing time is from fitting the "K.SAM" (particularly for large p(=500) scenarios)
  results.HD.K3.r  <- HD.Simulation.K3(n.rep = n.rep, sim.seed = sim.seed,
                                     n = n, p=p, delta = delta,
                                     trace = FALSE,
                                     sigma = 0.4, correlationX= 0, n.test = 10000,
                                     sigmaX = 1, type = "AIC")
  print(results.HD.K3.r$final)
  results.HD.K3[[r]] <- results.HD.K3.r$final
  results.aggregated.HD.K3[[r]] <- results.HD.K3.r$results.aggregated
}

results.HD.K3
print(r)
print(scenarios)

#save.file <- paste("csim-K3-", r, ".RData", sep="")
#save.image(save.file)



###########
## plots ##
###########
library(ggplot2)


# 1)
# for n =200 and delta = 1
Value  <- c(results.aggregated.HD.K3[[1]][1:n.rep,2],
            results.aggregated.HD.K3[[1]][1:n.rep,4],
            results.aggregated.HD.K3[[1]][1:n.rep,6],
            results.aggregated.HD.K3[[3]][1:n.rep,2],
            results.aggregated.HD.K3[[3]][1:n.rep,4],
            results.aggregated.HD.K3[[3]][1:n.rep,6])
Method <- rep(factor(c(rep("CSIM", n.rep),  rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM",  "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 3*n.rep), rep(500, 3*n.rep)))
data1 <- data.frame(Value, Method, p)



# 2)
# for n =400 and delta = 1
Value  <- c(results.aggregated.HD.K3[[2]][1:n.rep,2],
            results.aggregated.HD.K3[[2]][1:n.rep,4],
            results.aggregated.HD.K3[[2]][1:n.rep,6],
            results.aggregated.HD.K3[[4]][1:n.rep,2],
            results.aggregated.HD.K3[[4]][1:n.rep,4],
            results.aggregated.HD.K3[[4]][1:n.rep,6])
Method <- rep(factor(c(rep("CSIM", n.rep),  rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 3*n.rep), rep(500, 3*n.rep)))
data2 <- data.frame(Value, Method, p)


# 3)
# for n =200 and delta = 2
Value  <- c(results.aggregated.HD.K3[[5]][1:n.rep,2],
            results.aggregated.HD.K3[[5]][1:n.rep,4],
            results.aggregated.HD.K3[[5]][1:n.rep,6],
            results.aggregated.HD.K3[[7]][1:n.rep,2],
            results.aggregated.HD.K3[[7]][1:n.rep,4],
            results.aggregated.HD.K3[[7]][1:n.rep,6])
Method <- rep(factor(c(rep("CSIM", n.rep),  rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 3*n.rep), rep(500, 3*n.rep)))
data3 <- data.frame(Value, Method, p)


# 4)
# for n =400 and delta = 2
Value  <- c(results.aggregated.HD.K3[[6]][1:n.rep,2],
            results.aggregated.HD.K3[[6]][1:n.rep,4],
            results.aggregated.HD.K3[[6]][1:n.rep,6],
            results.aggregated.HD.K3[[8]][1:n.rep,2],
            results.aggregated.HD.K3[[8]][1:n.rep,4],
            results.aggregated.HD.K3[[8]][1:n.rep,6])
Method <- rep(factor(c(rep("CSIM", n.rep),  rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM", "MCA", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 3*n.rep), rep(500, 3*n.rep)))
data4 <- data.frame(Value, Method, p)




dodge <- position_dodge(width = 0.6)


P1  <- ggplot(data1, aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.6, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Moder. Main Eff. & "~n== 250)) + xlab(expression(p)) + xlab(expression(p)) + ylab(expression(Value / Value[opt]) )
P1



P2  <- ggplot(data2,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.6, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Moder. Main Eff. & "~n== 500))  + xlab(expression(p)) + xlab(expression(p))
P2



P3  <- ggplot(data3,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.6, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Big Main Eff. & "~n== 250)) + xlab(expression(p))
P3


P4 <- ggplot(data4,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.6, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Big Main Eff. & "~n== 500))  + xlab(expression(p))
P4




######################################################################
## END OF THE FILE
######################################################################
