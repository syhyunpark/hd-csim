# Simulation illustration for consistency

source('csim-main.R')
library(ggplot2)


######################################################

alphax <- c(1,1)
mux <- c(1,-1)
grid.length <- 100
n.rep  <- 200
theta = seq(0, pi, length.out=grid.length)

link.fn <- function(u, w=1)
{
 if(w==1)  return(cos(u) - 0.5)
 if(w==2)  return(u/2.15)
}

scenarios <- expand.grid(delta = c(1, 3, 5), w = c(1, 2))
scenarios

n <- 200
sigma <- 0.2
nbasis.t <- c(5,5,8)
rho.grid <- 0

results <- list()
sim.seed <- 1234

for(r in 1:nrow(scenarios))
{
  delta <- scenarios[r, 1]
  w <- scenarios[r, 2]

  results.csim = results.simml = results.mc  <- matrix(rep(0, n.rep*grid.length), ncol =n.rep)

  for(rep.number in 1:n.rep)
  {
    set.seed(sim.seed + rep.number)
    Tr <- rbinom(n, 1, 0.5) + 1
    x1 <- runif(n, -1, 1);
    x2 <- runif(n, -1, 1)
    X <- cbind(x1, x2)

    main.effect <- delta * cos(X %*% mux)
    interaction.effect <- (-1)^Tr *link.fn(X %*% alphax, w=w)
    var(main.effect)/var(interaction.effect)

    noise <- sigma * rnorm(n)
    y <- main.effect + interaction.effect + noise

    # some pre-procssing for data to be used
    datax <- data.preprocess(y, Tr, Xsc = X)
    y <- datax$y;
    Tr <- datax$Tr;
    X <- datax$Xsc

    K <- length(unique(Tr))
    # center y within each treatment, and store the treatment-specific mean at intercept.y
    dat <- data.frame(Tr=Tr, y=y)
    dat.list <- dlply(dat, .(dat$Tr), function(x) {as.matrix(x[,-1])});
    yc <- NULL
    for(t in 1:K)
    {
      dat.list[[t]] <- scale(dat.list[[t]], center = TRUE, scale = FALSE)
      yc <- c(yc, dat.list[[t]])
    }
    rm(dat); rm(dat.list);

    # look at the criterion values for a full range of possible linear combinations
    for (i in 1:length(theta))
    {
      itheta <- theta[i]
      alphai = sqrt(2) *c(cos(itheta), sin(itheta))
      ui <- drop(X %*% alphai)

      gt.obj1 <- fit.link.fn.gcv(yc, Tr, ui, nbasis.t = nbasis.t, rho.grid = rho.grid, ortho.constr=TRUE, ini=TRUE)
      gt.obj2 <- fit.link.fn.gcv(yc, Tr, ui, nbasis.t = nbasis.t, rho.grid = rho.grid, ortho.constr=FALSE, ini=TRUE)

      X.modified <- as.matrix(ui * (Tr-1.5) )
      mc.obj <- lm(yc  ~ X.modified )

      results.csim[i, rep.number] <- gt.obj1$MSE
      results.simml[i, rep.number] <- gt.obj2$MSE
      results.mc[i, rep.number] <- mean(mc.obj$residuals^2)
    }
    print(rep.number)
  }

  results[[r]] <- list(results.csim = results.csim, results.simml=results.simml, results.mc= results.mc)
}


plots <- list()

for(r in 1:nrow(scenarios))
{
  results1.mean <- apply(results[[r]]$results.csim, 1, mean)
  results2.mean <- apply(results[[r]]$results.simml, 1, mean)
  results3.mean <- apply(results[[r]]$results.mc, 1, mean)

  results1.mean <- as.data.frame(cbind(theta, results1.mean))
  colnames(results1.mean) <- c("Angle", "Objective")
  range.1 <- range(results1.mean[,2]);
  results1.mean$Objective.s <-  (results1.mean[,2] - range.1[1]) /(range.1[2] - range.1[1])

  results2.mean <- as.data.frame(cbind(theta, results2.mean))
  colnames(results2.mean) <- c("Angle", "Objective")
  range.2 <- range(results2.mean[,2]);
  results2.mean$Objective.s <-  (results2.mean[,2] - range.2[1]) /(range.2[2] - range.2[1])

  results3.mean <- as.data.frame(cbind(theta, results3.mean))
  colnames(results3.mean) <- c("Angle", "Objective")
  range.3 <- range(results3.mean[,2]);
  results3.mean$Objective.s <-  (results3.mean[,2] - range.3[1]) /(range.3[2] - range.3[1])

  Method = c(rep("Constrained SIM", grid.length), rep("Naive SIM", grid.length), rep("Modified covariates", grid.length))
  Method <- factor(Method, levels= c("Constrained SIM", "Naive SIM","Modified covariates"))
  dat1 <- data.frame(rbind(results1.mean, results2.mean, results3.mean), Method = Method)

  plot1 <- ggplot(data = dat1, aes(x = Angle, y=Objective.s, col = Method, linetype = Method)) +
    geom_line(lwd=0.7) + theme_bw() +
    scale_y_continuous("Criterion value", limits = c(0,1)) +
    geom_vline(aes(xintercept= atan(alphax[2]/alphax[1])), colour="grey", linetype="dashed", size=0.75)+
    geom_vline(aes(xintercept= pi+atan(mux[2]/mux[1])), colour="grey", linetype="dotted", size=0.75) +
    xlab(expression(theta)) + scale_x_continuous(breaks = c(0, 0.78, 1.57, 2.35, 3.14),
                                                 labels = c("0", expression(paste(pi, "/4")),
                                                            expression(paste(pi, "/2")), expression(paste("3",pi, "/4")), expression(pi)))+
    theme(axis.text.x = element_text(color = c("black", "purple", "black", "black", "black", "black")))

  plots[[r]] <- plot1

}

plots

######################################################################
## END OF THE FILE
######################################################################
