# CSIM paper example in Section 5.1 (K=2 treatment case)
source('csim-main.R')


# higher dimensional setting
HD.Simulation.K2  = function(n.rep = 1,
                             sim.seed = 1234,
                             n = 200, p=50, sigma = 0.4,
                             correlationX= 0, sigmaX = 1,
                             n.test = 10000,
                             w = 1, delta = 1,
                             trace = FALSE,
                             nbasis.t = c(6,6,8), rho.grid = 0, type = "AIC")
{

  true.alpha <- c(c(1, 0.5, 0.25, 0.125), rep(0, p-4))

  g <- function(u, w)
  {
    if(w==1) return(0.5* u)
    if(w==2) return(sin(u) - u)
    if(w==3) return(cos(u) - 0.5)
  }

  m <- function(u, delta= 1)   2+ delta*cos(u*0.5*pi)

  method.names <- c("CSIM", "MC", "K-LR", "K-SAM")
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
    eta.hold  <- eta.hold /sqrt(sum(eta.hold^2) )
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

    #1. CSIM
    csim.object  <- csim(y, Tr, X,
                         type = type,  # tuning parameter selection
                         lam.by = 0.03,
                         eps = 10^-4, it.max=50, trace=trace,
                         nbasis.t = nbasis.t,  # a vector of the numbers of basis funtions to be used for each treatment group
                         rho.grid = rho.grid,  # grid of (ridge-type) smoothing parameters for the links
                         eff.aug = FALSE,  # If TRUE, use a linear model as a main effect working model
                         linear.link = FALSE,
                         plots = FALSE)

    csim.object$alpha.coef
    # y.test, Tr.test, X.test are testing data
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



    # 4. K separate SAM
    K.SAM.object = K.SAM.obj <- K.SAM(y, Tr, X)
    pred.test <- predict(K.SAM.object, X.test)$pred.new
    performance.object <- performance_measure(pred.test, y.test, Tr.test, X.test, value.opt=value.opt, optTr=optTr)
    performance.object$pcd
    performance.object$value.s
    K.SAM.res <- c( K.SAM.pcd = performance.object$pcd,  K.SAM.value = performance.object$value.s)
    K.SAM.res

    contrast  =  K.SAM.object$K.SAM[[2]]$func_norm -  K.SAM.object$K.SAM[[1]]$func_norm
    Model.Size.Recorder[sim_run,4] <- sum( contrast !=0)
    #C and IC recorder
    # the number of "effect modifying (interaction)" signal predictors correctly identified as signal.
    Model.Sel.Recorder.alpha.C[sim_run,4] <- sum(contrast[true.alpha!=0] != 0)
    # the number of noise predictors incorrectly identified as signal.
    Model.Sel.Recorder.IC[sim_run,4] <- sum(contrast[ ! (true.alpha!=0) ] != 0)

    Value.Recorder[sim_run,4] <- performance.object$value.s
    Per.correct.Recorder[sim_run,4] <- performance.object$pcd



    results.ii <- c(csim.res,
                    mc.res,
                    K.LR.res,
                    K.SAM.res)
    results.aggregated <- rbind(results.aggregated, results.ii)
    # print(sim_run)
  } # the end of the for loop


  results.ii

  results.mean <- apply(results.aggregated, 2, mean)
  results.median <- apply(results.aggregated, 2, median)
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
              results.aggregated = results.aggregated,
              results.mean = results.mean, results.median = results.median, results.sd =results.sd,
              Value.Recorder=Value.Recorder, Per.correct.Recorder=Per.correct.Recorder,
              Model.Size.Recorder=Model.Size.Recorder,
              Model.Sel.Recorder.alpha.C=Model.Sel.Recorder.alpha.C,
              Model.Sel.Recorder.IC= Model.Sel.Recorder.IC))
}





scenarios <- expand.grid(n = c(250, 500),   # the number of observations
                         p = c(50, 500),  # the number of covariates
                         delta = c(1, 2),  # intensity of main effect; delta=1 moderate, delta=2 large
                         w = c(1,3))  # degree of nonlinearity; w=1 linear, w=3 highly nonlinear (cosine)
scenarios
n.rep <- 2 # 200  # the number of simulation replications for each scenario
sim.seed <- 1234

results.HD.K2 = results.aggregated.HD.K2 <- list()
for(r in 1:nrow(scenarios))
{
  n = scenarios[r, 1]
  p = scenarios[r, 2]
  delta = scenarios[r, 3]
  w = scenarios[r, 4]

  # note that most of the computing time is from fitting the "K.SAM" (particularly for large p(=500) scenarios)
  results.HD.K2.r  <- HD.Simulation.K2(n.rep = n.rep, sim.seed = sim.seed,
                                     n = n, p=p, w=w, delta = delta,
                                     trace = FALSE,
                                     sigma = 0.4, correlationX= 0, n.test = 10000,
                                     sigmaX = 1, type = "AIC")
  results.HD.K2[[r]] <- results.HD.K2.r$final
  results.aggregated.HD.K2[[r]] <- results.HD.K2.r$results.aggregated
}

results.HD.K2
print(r)
print(scenarios)

#save.file <- paste("csim-K2-", r, ".RData", sep="")
#save.image(save.file)






###########
## plots ##
###########
library(ggplot2)


dodge <- position_dodge(width = 0.6)
scenarios

# 1) for w= 1
# for n =200 and delta = 1
Value  <- c(results.aggregated.HD.K2[[1]][1:n.rep,2],
            results.aggregated.HD.K2[[1]][1:n.rep,4],
            results.aggregated.HD.K2[[1]][1:n.rep,6],
            results.aggregated.HD.K2[[1]][1:n.rep,8],
            results.aggregated.HD.K2[[3]][1:n.rep,2],
            results.aggregated.HD.K2[[3]][1:n.rep,4],
            results.aggregated.HD.K2[[3]][1:n.rep,6],
            results.aggregated.HD.K2[[3]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data1 <- data.frame(Value, Method, p)



# 2) for w= 1
# for n =400 and  delta = 1
Value  <- c(results.aggregated.HD.K2[[2]][1:n.rep,2],
            results.aggregated.HD.K2[[2]][1:n.rep,4],
            results.aggregated.HD.K2[[2]][1:n.rep,6],
            results.aggregated.HD.K2[[2]][1:n.rep,8],
            results.aggregated.HD.K2[[4]][1:n.rep,2],
            results.aggregated.HD.K2[[4]][1:n.rep,4],
            results.aggregated.HD.K2[[4]][1:n.rep,6],
            results.aggregated.HD.K2[[4]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data2 <- data.frame(Value, Method, p)


# 3) for w= 1
# for n =200 and  delta = 2
Value  <- c(results.aggregated.HD.K2[[5]][1:n.rep,2],
            results.aggregated.HD.K2[[5]][1:n.rep,4],
            results.aggregated.HD.K2[[5]][1:n.rep,6],
            results.aggregated.HD.K2[[5]][1:n.rep,8],
            results.aggregated.HD.K2[[7]][1:n.rep,2],
            results.aggregated.HD.K2[[7]][1:n.rep,4],
            results.aggregated.HD.K2[[7]][1:n.rep,6],
            results.aggregated.HD.K2[[7]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data3 <- data.frame(Value, Method, p)




# 4) for w= 1
# for n =400 and  delta = 2
Value  <- c(results.aggregated.HD.K2[[6]][1:n.rep,2],
            results.aggregated.HD.K2[[6]][1:n.rep,4],
            results.aggregated.HD.K2[[6]][1:n.rep,6],
            results.aggregated.HD.K2[[6]][1:n.rep,8],
            results.aggregated.HD.K2[[8]][1:n.rep,2],
            results.aggregated.HD.K2[[8]][1:n.rep,4],
            results.aggregated.HD.K2[[8]][1:n.rep,6],
            results.aggregated.HD.K2[[8]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data4 <- data.frame(Value, Method, p)




# 5) for w= 3
# for n =200 and  delta = 1
Value  <- c(results.aggregated.HD.K2[[9]][1:n.rep,2],
            results.aggregated.HD.K2[[9]][1:n.rep,4],
            results.aggregated.HD.K2[[9]][1:n.rep,6],
            results.aggregated.HD.K2[[9]][1:n.rep,8],
            results.aggregated.HD.K2[[11]][1:n.rep,2],
            results.aggregated.HD.K2[[11]][1:n.rep,4],
            results.aggregated.HD.K2[[11]][1:n.rep,6],
            results.aggregated.HD.K2[[11]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data5 <- data.frame(Value, Method, p)




# 6) for w= 3
# for n =400 and  delta = 1
Value  <- c(results.aggregated.HD.K2[[10]][1:n.rep,2],
            results.aggregated.HD.K2[[10]][1:n.rep,4],
            results.aggregated.HD.K2[[10]][1:n.rep,6],
            results.aggregated.HD.K2[[10]][1:n.rep,8],
            results.aggregated.HD.K2[[12]][1:n.rep,2],
            results.aggregated.HD.K2[[12]][1:n.rep,4],
            results.aggregated.HD.K2[[12]][1:n.rep,6],
            results.aggregated.HD.K2[[12]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data6 <- data.frame(Value, Method, p)




# 7) for w= 3
# for n =200 and  delta = 2
Value  <- c(results.aggregated.HD.K2[[13]][1:n.rep,2],
            results.aggregated.HD.K2[[13]][1:n.rep,4],
            results.aggregated.HD.K2[[13]][1:n.rep,6],
            results.aggregated.HD.K2[[13]][1:n.rep,8],
            results.aggregated.HD.K2[[15]][1:n.rep,2],
            results.aggregated.HD.K2[[15]][1:n.rep,4],
            results.aggregated.HD.K2[[15]][1:n.rep,6],
            results.aggregated.HD.K2[[15]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data7 <- data.frame(Value, Method, p)



# 8) for w= 3
# for n =400 and  delta = 2
Value  <- c(results.aggregated.HD.K2[[14]][1:n.rep,2],
            results.aggregated.HD.K2[[14]][1:n.rep,4],
            results.aggregated.HD.K2[[14]][1:n.rep,6],
            results.aggregated.HD.K2[[14]][1:n.rep,8],
            results.aggregated.HD.K2[[16]][1:n.rep,2],
            results.aggregated.HD.K2[[16]][1:n.rep,4],
            results.aggregated.HD.K2[[16]][1:n.rep,6],
            results.aggregated.HD.K2[[16]][1:n.rep,8])
Method <- rep(factor(c(rep("CSIM", n.rep), rep("MC", n.rep), rep("K-LR", n.rep), rep("K-SAM", n.rep)),
                     levels = c("CSIM","MC", "K-LR", "K-SAM")), 2)
p  <- factor(c(rep(50, 4*n.rep), rep(500, 4*n.rep)))
data8 <- data.frame(Value, Method, p)




P1  <- ggplot(data1, aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Moder. Main Eff. & "~n== 250))  + xlab(expression(p)) + xlab(expression(p)) + ylab(expression(Value / Value[opt]) )
P1



P2  <- ggplot(data2,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Moder. Main Eff. & "~n== 500))  + xlab(expression(p)) +
  xlab(expression(p))
P2



P3  <- ggplot(data3,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Big Main Eff. & "~n== 250)) + xlab(expression(p))
P3


P4 <- ggplot(data4,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Big Main Eff. & "~n== 500)) + xlab(expression(p))
P4


P5  <- ggplot(data5,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ylab(expression(Value / Value[opt]) ) +
  ggtitle(bquote("Moder. Main Eff. & "~n== 250))  + xlab(expression(p))
P5



P6  <- ggplot(data6, aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Moder. Main Eff. & "~n== 500))  + xlab(expression(p))
P6



P7  <- ggplot(data7,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right")+
  ggtitle(bquote("Big Main Eff. & "~n== 250))  + xlab(expression(p))
P7


P8  <- ggplot(data8,  aes(x=p, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value") + ylim(0.8, 1) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=p, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="PuRd") + theme_minimal() + theme(legend.position="right") +
  ggtitle(bquote("Big Main Eff. & "~n== 500)) + xlab(expression(p))
P8



######################################################################
## END OF THE FILE
######################################################################
