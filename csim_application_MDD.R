# Constrained SIM application to a depression (MDD) dataset (Section 7 of the paper)
source('csim-main.R')


## load the dataset
load("datasetChangeScore8w.rda")
y1 <- datasetChangeScore8w$y1  # outcome vector for Trt group 1
y2 <- datasetChangeScore8w$y2  # Trt group 2
X1 <- datasetChangeScore8w$X1  # predictor matrix for Trt group 1
X2 <- datasetChangeScore8w$X2  # Trt group 2
n1 <- length(y1)
n2 <- length(y2)
n <- n1 + n2

Tr <- c(rep(1,n1), rep(2,n2))  # a vector of treatment indicators
y <- c(y1,y2)    # outcome vector y
X.temp <- rbind(X1, X2)   # X matrix
X <- X.temp[,-c(2, 6, 7, 8, 9, 10, 11)]  # disregard discrete predictors
X[,3] <- log(X[,3]+1)   # log transformation for the variable "Dur. of MDE"


rho.grid = c(0, 0.25, 0.5)
nbasis.t <- c(6,6,8)

## estimate the constrained single index model
set.seed(1234)
# takes about 20 seconds, including tuning parameter selection 
csim.obj   <- csim(y, Tr, X, type = "AIC", rho.grid = rho.grid, nbasis.t=nbasis.t, eff.aug = TRUE, n.max = 7, trace = FALSE,  it.max = 100, plots=FALSE)
i.fx <- which(csim.obj$alpha.coef==1)
csim.obj$alpha.coef


## Value comparison
results.aggregated2  <- NULL   # will store results here
n.rep <- 500  # the number of replications for computing the "Value"

for(rep.number in 1:n.rep)
{
  print(rep.number)

  set.seed(1234 + rep.number)
  folds <- sample(1:2, n, replace = TRUE, prob = c(0.83, 0.17))
  y.train <- y[folds==1]
  y.test <- y[folds==2]
  X.train <- X[folds==1,]
  X.test <- X[folds==2,]
  Tr.train <- Tr[folds==1]
  Tr.test <- Tr[folds==2]


  ## 1. estimate the constrained single index model, based on the training set, with variable selection
  csim.obj   <- csim(y.train, Tr.train, X.train, type = "AIC",
                     nbasis.t=nbasis.t, rho.grid=rho.grid, eff.aug = TRUE, i.fx=i.fx, n.max = 7, it.max = 100, plots=F)
  csim.obj$alpha.coef
  ## performance assessment:
  pred.test <- predict(csim.obj, X.test)$pred.new
  csim.value  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  csim.value


  ## 2. estimate the c-sim without variable selection
  csim.obj2   <- csim(y.train, Tr.train, X.train, sparse=FALSE,
                      nbasis.t=nbasis.t, rho.grid=rho.grid, eff.aug = TRUE,  i.fx=i.fx, plots=F)
  csim.obj2$alpha.coef
  ## performance assessment:
  pred.test <- predict(csim.obj2, X.test)$pred.new
  csim.value2  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  csim.value2


  ## 3. estimate the c-sim with variable selection, without efficiency augmentation
  csim.obj3   <- csim(y.train, Tr.train, X.train, sparse=TRUE,
                      nbasis.t=nbasis.t, rho.grid=rho.grid, eff.aug =FALSE, i.fx=i.fx, n.max = 7, it.max = 100, plots=F)
  csim.obj3$alpha.coef
  pred.test <- predict(csim.obj3, X.test)$pred.new
  csim.value3  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  csim.value3


  ## 4. estimate the c-sim without variable selection, without efficiency augmentation
  csim.obj4   <- csim(y.train, Tr.train, X.train, sparse=FALSE,
                      nbasis.t=nbasis.t, rho.grid=rho.grid, eff.aug =FALSE, i.fx=i.fx, plots=F)
  csim.obj4$alpha.coef
  pred.test <- predict(csim.obj4, X.test)$pred.new
  csim.value4  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  csim.value4



  ## 5. fit the modificed covariate model with effiency augmentation
  mc.obj <- mc(y.train, Tr.train, X.train, eff.aug = TRUE, use.lasso=TRUE)
  mc.obj$alpha.coef
  pred.test <- predict.mc(mc.obj, X.test)$pred.new
  mc.value  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  mc.value

  ## 6. fit the modificed covariate model with effiency augmentationm, without variable selection
  mc.obj2 <- mc(y.train, Tr.train, X.train, eff.aug = TRUE, use.lasso=FALSE)
  mc.obj2$alpha.coef
  pred.test <- predict.mc(mc.obj2, X.test)$pred.new
  mc.value2  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  mc.value2



  ## 7. fit a system of K linear regression with lasso
  K.LR.obj <- K.LR(y.train, Tr.train, X.train, use.lasso=TRUE)
  pred.test <- predict(K.LR.obj, X.test)$pred.new
  K.LR.value <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  K.LR.value


  ## 8. fit a system of K separate sparse addtive models, one for each treatment group
  K.SAM.obj <- K.SAM(y.train, Tr.train, X.train)
  pred.test <- predict(K.SAM.obj, X.test)$pred.new
  K.SAM.value <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  K.SAM.value


  all.placebo <- mean(y.test[Tr.test == 1])
  all.placebo

  all.drug <- mean(y.test[Tr.test == 2])
  all.drug

  results.ii <- c(csim.value, csim.value2, csim.value3, csim.value4, mc.value, mc.value2, K.LR.value, K.SAM.value,  all.placebo, all.drug)

  print(results.ii)
  results.aggregated2 <- rbind(results.aggregated2, results.ii)
}


results.aggregated <- results.aggregated2[complete.cases(results.aggregated2), ]
colnames(results.aggregated) <-  c("CSIM(VS)", "CSIM(no VS)", "CSIM(VS, no eff.)", "CSIM(no VS, no eff.)", "MC(VS)", "MC(no VS)",
                                   "K-LR", "K-SAM", "all PBO", "all DRUG")
apply(results.aggregated, 2, mean)
apply(results.aggregated, 2, sd)
apply(results.aggregated, 2, median)

#save.file <- "csim-EMBARC.RData"
#save.image(save.file)






###########
## plots ##
###########

library(gridExtra)
library(grid)

### value plot
n.rep <- nrow(results.aggregated)
apply(results.aggregated[1:n.rep,], 2, mean)
apply(results.aggregated[1:n.rep,], 2, median)


Value  <- c(
  results.aggregated[1:n.rep,1], # C-SIM
  results.aggregated[1:n.rep,5], # MC
  results.aggregated[1:n.rep,2], # C-SIM(no reg.)
  results.aggregated[1:n.rep,6], # MC(no reg.)
  results.aggregated[1:n.rep,7], # K-LR
  results.aggregated[1:n.rep,8], # K-SAM
  results.aggregated[1:n.rep,9],
  results.aggregated[1:n.rep,10])
Method <- factor(c(
  rep("CSIM", n.rep), rep("MC", n.rep),
  rep("CSIM (no VS)", n.rep), rep("MC (no VS)", n.rep),
  rep("K-LR", n.rep), rep("K-SAM", n.rep),
  rep("All PBO", n.rep), rep("All DRUG", n.rep)),
  levels = c("CSIM", "MC",
             "K-LR", "K-SAM",
             "All PBO", "All DRUG",
             "CSIM (no VS)", "MC (no VS)") )
value.data1 <- data.frame(Value, Method)



dodge <- position_dodge(width = 0.6)
P1 <- ggplot(value.data1, aes(x=Method, y=Value, fill=Method)) +
  geom_boxplot(width=.2, outlier.colour=NA, #position = dodge,
               aes(x=Method, y=Value, fill=Method, shape = Method, color = Method) ) +
  scale_fill_manual(values=c("lightcyan", "lightcyan", "lightcyan",
                             "lightcyan", #"lightcyan", "lightcyan",
                             "royalblue", "firebrick3", "lightcyan", "lightcyan")) +   #geom_violin(trim=TRUE, fill="white", position = dodge) +
  labs(x=" ", y = "Value") + ylim(2.5, 13) + #ylim(0.5, 14.5) +
  geom_hline(yintercept = 7.32,  colour="firebrick3", linetype="dashed", size = 0.7) +
  geom_hline(yintercept = 6.33,  colour="royalblue", linetype="dashed", size = 0.7) +
  stat_summary(fun.y=mean, geom="point", size=1, color="red") +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5))

P1 + theme_light(base_size=14)+ theme(legend.position="none")





#############################################
## analysis results appeared in the paper  ##
#############################################
rho.grid <- c(0, 0.25, 0.5)
nbasis.t <- c(6,6,8)
set.seed(1234)
csim.obj  <- csim(y, Tr, X, type = "AIC",  rho.grid = rho.grid, nbasis.t=nbasis.t,
                  eff.aug = TRUE, trace = FALSE, n.max = 7, it.max = 100, plots=TRUE)

csim.obj$alpha.coef
csim.obj$link.fn.plot
csim.obj$eta.coef




###
cfs.obj   <- fit.csim.cv(y, Tr, X,
                         type = "AIC",
                         n.max = 7, it.max=100,
                         nbasis.t = nbasis.t, rho.grid = rho.grid, eff.aug = TRUE)

# estimate the link functions using B-splines, given the single-index coefficients.
gt.obj <- fit.link.fn.gcv(cfs.obj$yc, cfs.obj$Tr, u= drop(cfs.obj$Xc %*% cfs.obj$coef.opt),
                     nbasis.t = cfs.obj$nbasis.t, rho.grid = rho.grid)
rho.opt = 1  #gt.obj$rho.opt 

K <- length(unique(Tr))
y.means <- rep(0, length(y));
for(t in 1:K)  y.means  <- y.means +  cfs.obj$intercept.y[[t]]*(cfs.obj$Tr==t);

dat  <- data.frame(y =  y.means + cfs.obj$yc- gt.obj$smoother$B0 %*% gt.obj$beta.0.coef,
                   x = gt.obj$smoother$u.t[[K+1]],
                   Treatment = factor(Tr, labels=c("Placebo","Active drug")) );

gt.plot  <- ggplot(dat, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
  geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
  scale_colour_brewer(palette = "Set1", direction = -1)  +
  theme( axis.title.x=element_text(size=15,face="bold"))  +
  theme(title =element_text(size=12)) +
  xlab(expression(paste(beta*minute, "x"))) +  ylab("(Adjusted) Improvement") + theme_light(base_size = 13);

tmp1 <- 0;  for(t in 1:K) tmp1 <- tmp1 + gt.obj$smoother$nbasis.t[t];

gt.plot <- gt.plot + geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= floor(tmp1/K), sp= rho.opt), se=TRUE, fullrange=TRUE, alpha = 0.35)
gt.plot 

plot1  <- gt.plot + ylab("(Adjusted) Response") + ggtitle("Single-index plot") # + ylim(c(-10,25))

plot1



rho.opt <- 1 # 0.05
###
  dat1  <- data.frame(y= y.means + cfs.obj$yc- gt.obj$smoother$B0 %*% gt.obj$beta.0.coef, x = X[,1],Treatment = factor(Tr, labels=c("Placebo","Active drug")) )
  gt.plot1  <- ggplot(dat1, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
   geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
   scale_colour_brewer(palette = "Set1", direction = -1)  +
   theme( axis.title.x=element_text(size=15,face="bold"))  +
   theme(title =element_text(size=12)) + theme_light(base_size = 14);

  gt.plot1 <- gt.plot1 + geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= 6, sp= rho.opt), se=TRUE, fullrange=TRUE, alpha = 0.35)
  P1 <- gt.plot1 + ylab("(Adjusted) Response") + xlab("Age at evaluation")
  P1


  dat2  <- data.frame(y=y.means + cfs.obj$yc- gt.obj$smoother$B0 %*% gt.obj$beta.0.coef, x = X[,2],Treatment = factor(Tr, labels=c("Placebo","Active drug")) )
  gt.plot2  <- ggplot(dat2, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
    geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
    scale_colour_brewer(palette = "Set1", direction = -1)  +
    theme( axis.title.x=element_text(size=15,face="bold"))  +
    theme(title =element_text(size=12)) + theme_light(base_size = 14);

  gt.plot2 <- gt.plot2 + geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= 6, sp= rho.opt+1), se=TRUE, fullrange=TRUE, alpha = 0.35)
  P2 <- gt.plot2 + ylab("(Adjusted) Response") + xlab("Symptom severity")
  P2


  dat3  <- data.frame(y=y.means + cfs.obj$yc- gt.obj$smoother$B0 %*% gt.obj$beta.0.coef, x = X[, 3],Treatment = factor(Tr, labels=c("Placebo","Active drug")) )
  gt.plot3  <- ggplot(dat3, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
    geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
    scale_colour_brewer(palette = "Set1", direction = -1)  +
    theme( axis.title.x=element_text(size=15,face="bold"))  +
    theme(title =element_text(size=12)) + theme_light(base_size = 14);

  gt.plot3 <- gt.plot3 + geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= 6, sp= rho.opt), se=TRUE, fullrange=TRUE, alpha = 0.35)
  P3 <- gt.plot3 + ylab("(Adjusted) Response") + xlab("(log) Dur. of MDE")
  P3

  
  dat4  <- data.frame(y=y.means + cfs.obj$yc- gt.obj$smoother$B0 %*% gt.obj$beta.0.coef, x = X[,11],Treatment = factor(Tr, labels=c("Placebo","Active drug")) )
  gt.plot4  <- ggplot(dat4, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
    geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
    scale_colour_brewer(palette = "Set1", direction = -1)  +
    theme( axis.title.x=element_text(size=15,face="bold"))  +
    theme(title =element_text(size=12)) + theme_light(base_size = 14);

  gt.plot4 <- gt.plot4 + geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= 6, sp= rho.opt), se=TRUE, fullrange=TRUE, alpha = 0.35)
  P4 <- gt.plot4 + ylab("(Adjusted) Response") + xlab("Flanker Accuracy")
  P4

 
#grid.arrange(P1 + theme(legend.position="none"), P2 + theme(legend.position="none"), 
#             P3 + theme(legend.position="none"), P4 + theme(legend.position="none"),
#             nrow = 1)

library(lemon)
grid_arrange_shared_legend(P1, 
                           P2+ylab(" "), 
                           P3+ylab(" "), 
                           P4+ylab(" "),
                           nrow = 1, position='right')


 


###################
## bootstrap CIs ##
###################
set.seed(1234)
system.time(
  csim.obj.sparse.boot  <- csim(y, Tr, X, type = "AIC",  rho.grid = rho.grid, nbasis.t=nbasis.t,
                                eff.aug = TRUE, trace = FALSE, n.max = 7, it.max = 100, boot.CI = TRUE, n.boot=500,  plots=FALSE)
)
csim.obj.sparse.boot$alpha.coef
csim.obj.sparse.boot$boot.results
for(s in 1:ncol(X)){
  cat("s=", s, "\n");
  print(boot.ci(csim.obj.sparse.boot$boot.results, type = "basic", index = s) )
}


set.seed(1234)
system.time(
  csim.obj.nonsparse.boot  <- csim(y, Tr, X, sparse =FALSE, rho.grid = rho.grid, nbasis.t=nbasis.t, eff.aug = TRUE,  n.max = 7,  it.max = 100,  
                                   boot.CI = TRUE, n.boot=500,  plots=FALSE)
)
for(s in 1:ncol(X)){
  cat("s=", s, "\n");
  print(boot.ci(csim.obj.nonsparse.boot$boot.results, type = "basic", index = s) )
}
csim.obj.nonsparse.boot$alpha.coef
csim.obj.nonsparse.boot$boot.results
#boot.ci(csim.obj.nonsparse.boot$boot.results, type = "basic", index = 11)


set.seed(1234)
mc.obj <- mc(y, Tr, X, eff.aug  = TRUE, use.lasso=TRUE, boot.CI = TRUE, n.boot=500)
mc.obj$alpha.coef
for(s in 1:ncol(X)){
  cat("s=", s, "\n");
  print(boot.ci(mc.obj$boot.results, type = "basic", index = s) )
}


set.seed(1234)
mc.obj2 <- mc(y, Tr, X, eff.aug  = TRUE, use.lasso=FALSE, boot.CI = TRUE, n.boot=500)
mc.obj2$alpha.coef
for(s in 1:ncol(X)){
  cat("s=", s, "\n");
  print(boot.ci(mc.obj2$boot.results, type = "basic", index = s) )
}


#save.file <- "csim-EMBARC.RData"
#save.image(save.file)


######################################################################
## END OF THE FILE
######################################################################
