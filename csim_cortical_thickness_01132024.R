# 1/13/2024
# Constrained SIM application to a depression (MDD) dataset (Section 7 of the paper)
# a higher dimensional setting with p=161 (n=103), incorporating pretreatment sMRI cortical thickness measures.

## required packages in "csim-main.R": 
# library(ggplot2)
# library(plyr)
# library(splines)
# library(boot)
# library(glmnet)
# library(MASS)
# library(magic)  # for adiag()
# library(mgcv) 
source('csim-main.R')
source('csim-other-functions.R') 

# load the (higher dimensinoal) data 
load("embarc_cortical_thickness_n103.rda")
colnames(data)
n <- nrow(data)

library(dplyr)
dat <- data %>% dplyr::select(-c(sex,Greg_FH,fatigue,hypersomnia,axis2,anger_attack,anxious)) %>% 
  mutate(dur_MDE = log(dur_MDE+1))  
colnames(dat)  
summary(dat)
Tr <- as.numeric(as.factor(dat[,"Stage1TX"]))
y <- dat[,"changeScore"]
X <- dat %>% dplyr::select(-c(ProjectSpecificId,Stage1TX,changeScore)) %>% as.matrix
dim(X)
colnames(X) 



## estimate the constrained single index model
rho.grid = c(0, 0.25, 0.5)
nbasis.t <- c(6,6,8) 
set.seed(1234)
# takes about 20 seconds, including tuning parameter selection 
csim.obj   <- csim(y, Tr, X, type = "AIC", rho.grid = rho.grid, nbasis.t=nbasis.t, eff.aug = FALSE, n.max = 7, trace =FALSE,  it.max = 100, plots=TRUE)
i.fx <- which(csim.obj$alpha.coef==1) 
ls(csim.obj)
#csim.obj$link.fn.plot

## Value (and computation time) comparison
results.aggregated2 = time.aggregated2 <- NULL   # will store results here
n.rep <- 500 # the number of replications for computing the "Value"

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
  csim.time <- system.time(
    csim.obj   <- csim(y.train, Tr.train, X.train, type = "AIC",
                       nbasis.t=nbasis.t, rho.grid=rho.grid, eff.aug = FALSE, i.fx=i.fx, n.max = 7, it.max = 100, plots=F)
  )
  csim.obj$alpha.coef
  ## performance assessment:
  pred.test <- predict(csim.obj, X.test)$pred.new
  csim.value  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  csim.value
  
  ## 2. fit the modificed covariate model with effiency augmentation
  mc.time <- system.time(
    mc.obj <- mc(y.train, Tr.train, X.train, eff.aug =FALSE, use.lasso=TRUE)
  )
  mc.obj$alpha.coef
  pred.test <- predict.mc(mc.obj, X.test)$pred.new
  mc.value  <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  mc.value
  
  ## 3. fit a system of K linear regression with lasso
  K.LR.time <- system.time(
    K.LR.obj <- K.LR(y.train, Tr.train, X.train, use.lasso=TRUE)
  )
  pred.test <- predict(K.LR.obj, X.test)$pred.new
  K.LR.value <- performance_measure(pred.test, y.test, Tr.test, X.test)$value
  K.LR.value
  
  ## 4. fit a system of K separate sparse addtive models, one for each treatment group
  K.SAM.time <- system.time(
    K.SAM.obj <- K.SAM(y.train, Tr.train, X.train)
  )
  pred.test <- predict(K.SAM.obj, X.test)$pred.new
  pred.test
  K.SAM.value <- performance_measure(pred.test, y.test, Tr.test, X.test)$value+1
  K.SAM.value 
  
  all.placebo <- mean(y.test[Tr.test == 1])
  all.placebo
  
  all.drug <- mean(y.test[Tr.test == 2])
  all.drug
  
  results.ii <- c(csim.value, mc.value,  K.LR.value, K.SAM.value,  all.placebo, all.drug)
  time.ii <- c(csim.time[3], mc.time[3], K.LR.time[3], K.SAM.time[3])
  
  print(results.ii)
  print(time.ii)
  results.aggregated2 <- rbind(results.aggregated2, results.ii)
  time.aggregated2 <- rbind(time.aggregated2, time.ii)
}


results.aggregated <- results.aggregated2[complete.cases(results.aggregated2), ]
colnames(results.aggregated) <-  c("CSIM(VS)","MC(VS)", "K-LR", "K-SAM", "all PBO", "all DRUG")
apply(results.aggregated, 2, mean)
apply(results.aggregated, 2, sd)
apply(results.aggregated, 2, median)

apply(time.aggregated2, 2, mean)
apply(time.aggregated2, 2, median)
apply(time.aggregated2, 2, sd)

#save.file <- "csim-EMBARC2_01132024.RData"
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
  results.aggregated[1:n.rep,2], # MC 
  results.aggregated[1:n.rep,3], # K-LR
  results.aggregated[1:n.rep,4], # K-SAM
  results.aggregated[1:n.rep,5],
  results.aggregated[1:n.rep,6])
Method <- factor(c(
  rep("CSIM", n.rep), rep("MC", n.rep), 
  rep("K-LR", n.rep), rep("K-SAM", n.rep),
  rep("All PBO", n.rep), rep("All DRUG", n.rep)),
  levels = c("CSIM", "MC",
             "K-LR", "K-SAM",
             "All PBO", "All DRUG")) 
value.data1 <- data.frame(Value, Method)


#mean(y[Tr==1])
#mean(y[Tr==2])
dodge <- position_dodge(width = 0.6)
P1 <- ggplot(value.data1, aes(x=Method, y=Value, fill=Method)) +
  geom_boxplot(width=.2, outlier.colour=NA, #position = dodge,
               aes(x=Method, y=Value, fill=Method, shape = Method, color = Method) ) +
  scale_fill_manual(values=c("lightcyan", "lightcyan", "lightcyan",
                             "lightcyan", #"lightcyan", "lightcyan",
                             "royalblue", "firebrick3", "lightcyan", "lightcyan")) +   #geom_violin(trim=TRUE, fill="white", position = dodge) +
  labs(x=" ", y = "Value") + ylim(1.5, 16) + #ylim(2.5, 13) +
  geom_hline(yintercept = 7.26,  colour="firebrick3", linetype="dashed", size = 0.7) +
  geom_hline(yintercept = 6.99,  colour="royalblue", linetype="dashed", size = 0.7) +
  stat_summary(fun.y=mean, geom="point", size=1, color="red") +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5))

P1 + theme_light(base_size=14)+ theme(legend.position="none")


######################################################################
## END OF THE FILE
######################################################################