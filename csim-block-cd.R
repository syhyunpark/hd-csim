
######################################################################
# To optimize the current nonzero (active) single index coefficients,
# block.cd implements the block-coordinate descent of Radchenko (2015), given the current index u = alpha'X;
# return: w.a, the updated (active) coefficients; u, the updated single index; it, the number of iterations
######################################################################
block.cd <- function(w.a, ind.m, ind.old,
                     u, y, Tr, X.a, d.link.fn.obj,
                     it.max =60, eps=10^-4, trace=FALSE,
                     nbasis.t = c(6,6,8),
                     rho.grid = c(0, 0.25, 0.5),
                     linear.link.fn= FALSE,
                     ortho.constr = TRUE)
{

  w.a.old <- w.a + 5
  it <- 0;
  chng.fl <- FALSE;
  w.a.start <- w.a;

  while((max(abs(w.a- w.a.old))>eps)&(max(abs(w.a- w.a.start))>eps/100)|(it==0))
  {
    w.a.old <- w.a;
    it <- it + 1
    for(j in 1:length(w.a))
    {
      gr.m <- crossprod(d.link.fn.obj$d.link.fn1, X.a[,ind.m]);
      gr1 <- crossprod(d.link.fn.obj$d.link.fn1, X.a[,j])
      tmp1 <- apply(as.matrix(X.a[,ind.old]), 2, function(v){crossprod(d.link.fn.obj$d.link.fn1, v)})
      tmp2 <- -sign(w.a[ind.old]);
      thr1 <- tmp1*tmp2
      thr2 <- -sign(w.a[ind.m])*crossprod(d.link.fn.obj$d.link.fn1, X.a[,ind.m]);
      rm(tmp1);
      rm(tmp2)
      thr.entr <- max(c(thr1,thr2))

      # it all begins when the gr1 exceeds the threshold.
      if( ( (abs(gr1)>=thr.entr) | (w.a[j]!=0) ) & (j!=ind.m) )
      {
        if(w.a[j]==0) {sn= -sign(gr1)} else{sn= sign(w.a[j])}
        gr.cmb <- crossprod(d.link.fn.obj$d.link.fn1, (X.a[,j] - sign(w.a[ind.m])*sn*X.a[,ind.m]))
        hes.cmb <- crossprod((X.a[,j] - sign(w.a[ind.m])*sn*X.a[,ind.m])^2, d.link.fn.obj$d.link.fn2)
        # compute "delta_j"
        dj <- -gr.cmb/hes.cmb
        old.w.aj <- w.a[j];
        old.w.am <- w.a[ind.m]
        w.a[j] <- w.a[j]+dj
        if(sign(w.a[j]*old.w.aj) < 0)
        {
          dj <- -old.w.aj;
          w.a[j] <- 0;
          w.a[ind.m] <- w.a[ind.m]*(abs(w.a[ind.m]) + abs(old.w.aj))/abs(w.a[ind.m]);
          dm <- w.a[ind.m]- old.w.am
        }
        else{
          dm <- -dj*sign(w.a[ind.m])*sign(w.a[j]);
          w.a[ind.m] <- w.a[ind.m] - dj*sign(w.a[ind.m])*sign(w.a[j])
        }
        if(sign(w.a[ind.m]*old.w.am)<=0)
        {
          dm <- -old.w.am;
          w.a[ind.m] <- 0;
          chng.fl <- TRUE;
          if(trace) cat("Ch")
          if(old.w.aj!=0)
          {
            w.a[j] <- old.w.aj*(abs(old.w.aj) + abs(old.w.am))/abs(old.w.aj)
          }
          else{
            w.a[j] <- sign(w.a[j])*abs(old.w.am)
          }
          dj <- w.a[j] - old.w.aj
        }

        u <- u + dj*X.a[,j]+ dm*X.a[,ind.m]

        # update the gt and their 1st derivatives
        link.fn.obj <- fit.link.fn.gcv(y, Tr, u= u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link.fn = linear.link.fn, ortho.constr = ortho.constr)
        d.link.fn.obj <- deriv.link.fn(link.fn.obj)

        if(trace) cat("wa1=",w.a,"\n");
        if(trace) cat("grads",gr1, gr.m,"\n")
        if(chng.fl)
        {
          ind.m <- which.max(abs(w.a));
          chng.fl <- FALSE
        }
      }
    }
    if(it > it.max)
    {
      if(trace) cat("\n","DID NOT CONVERGE","\n"); break
    }
  }

  results <- list(w.a=w.a, u=u, it=it)
  return(results)
}




######################################################################
# a wrapper function to fit the (B-spline approximated) link functions gt, given the current index u = alpha'X;
# an optimal smoothing parameter, rho.opt, is chosen by minimizing GCV.
# return: beta.t.coef, beta.0.coef, smoother, resid, and MSE.
######################################################################
fit.link.fn.gcv <- function(y, Tr, u, nbasis.t = NULL, rho.grid = c(0, 0.25, 0.5), linear.link.fn = FALSE, ortho.constr=TRUE, ini=FALSE)
{

  if(is.null(nbasis.t))
  {
    n <- length(y);
    nt <- summary(as.factor(Tr));
    K <- length(nt);
    for(t in 1:K)
    {
      nbasis.t[t] <- floor(nt[t]^{1/5.5})  + 4
    }
    nbasis.t[K+1] <- floor(n^{1/5.5}) + 6
  }

  if(length(rho.grid) >1)
  {
    smoother <-  smoother.fn(Tr, u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link.fn = linear.link.fn)
    rho.opt <- fit.link.fn(y,  smoother, ortho.constr=ortho.constr, ini=ini)$rho.opt
  }else{
    rho.opt <- rho.grid
  }
  smoother <-  smoother.fn(Tr, u,  nbasis.t=nbasis.t, rho.grid = rho.opt, linear.link.fn =linear.link.fn)
  link.fn.obj <- fit.link.fn(y, smoother, ortho.constr=ortho.constr, ini=ini)
  return(link.fn.obj)
}



######################################################################
# a subfunction to construct (B-spline) smoother matrices, given the current index u = alpha'X.
# return: the QR decomposed design matrices (and the knot sequences used in constructing the B-spline design matrices);
######################################################################
smoother.fn <- function(Tr, u, nbasis.t=c(6,6,8), rho.grid = c(0, 0.25, 0.5), linear.link.fn = FALSE)
{

  # create a list, dat.list, grouped by the treatment indicator
  dat <- data.frame(Tr=Tr, u=u)
  dat_list <- dlply(dat, .(dat$Tr), function(dat) as.matrix(dat[,-1]) )

  K <- length(unique(Tr))
  u.t = design.t <- vector("list", K+1)
  for(t in 1:K)
  {
    u.t[[t]] <-  dat_list[[t]][,1]    # data points from the tth treatment group
  }
  u.t[[K+1]] <- u
  u.min = min(u)
  u.max = max(u)

  # construct treatment group-specific design matrices
  design.t = knots.t  <- vector("list", K+1)
  for(t in 1:(K+1))
  {
    if(linear.link.fn)   # if linear.link.fn==TRUE, construct the linear model design matrix
    {
      nbasis.t[t] <- 2
      design.t[[t]] <- cbind(1, u.t[[t]])
    }else{
      knots.t[[t]] <- c(rep(u.min, 3),  quantile(u.t[[t]], probs = seq(0, 1, length = nbasis.t[t] -2)),  rep(u.max,3))
      design.t[[t]] <-  splineDesign(knots.t[[t]], x= u.t[[t]], outer.ok = TRUE)
    }
  }

  # construct the block-diagonal matrix consist of the treatment-specific design matrices, to approximate E(Y| u, T)
  design.t.block <- NULL
  for(t in 1:K)
  {
    design.t.block <- c( design.t.block, list(design.t[[t]]))
  }
  Bt <- Reduce(adiag,  design.t.block)  # Bt is the block-diagonal design matrix

  # QR decomposition of the design matrix Bt, given each value of the smoothness tuning parameter, rho
  Bt.qr <- vector("list", length(rho.grid))
  ncol.Bt <- ncol(Bt)
  for(r in seq_along(rho.grid)) # a ridge-type smoothing (equivalent to a regular least squares with added observations)
  {
    Bt.qr[[r]] <- qr( rbind(Bt, diag(sqrt(rho.grid[r]), ncol.Bt) ) )
  }

  # compute effective degrees of freedom of smoothers, so that later we use GCV to select an optimal smoothing parameter
  edf <- vector("list", length=length(rho.grid))
  #edf <- NULL
  #if(length(rho.grid) > 1)
  #{
    svd.Bt <- svd(Bt)
    for(r in seq_along(rho.grid))
    {
      edf[[r]] <-  sum(svd.Bt$d^2/(svd.Bt$d^2 +rho.grid[[r]] )) #/K
    }
  #}

  # QR decomposition of the design matrix B0
  B0= design.t[[K+1]]
  B0.qr  <- qr(B0)

  results <- list(Bt.qr= Bt.qr, B0.qr= B0.qr, Bt=Bt, B0= B0,
                  u.t = u.t, u.min = u.min, u.max = u.max,
                  edf= edf, rho.grid = rho.grid, K=K,
                  knots.t = knots.t, nbasis.t = nbasis.t, ncol.Bt = ncol.Bt,
                  linear.link.fn=linear.link.fn)
  return(results)
}




######################################################################
# a subfunction to fit the (B-spline approximated) link functions gt, given the current index u = alpha'X.
# an optimal smoothing parameter, rho.opt, is chosen by minimizing GCV.
# return: beta.t.coef, beta.0.coef, smoother, resid, and MSE.
######################################################################
fit.link.fn <- function(y, smoother, ortho.constr = TRUE, ini=FALSE)
{

  options(warn=-1)
  # a ridge-type regularization (equivalent to an OLS with added 0s)
  y.aug <- c(y, rep(0, smoother$ncol.Bt))
  n <- length(y)

  # pick an optimal regularization (smoothing) paramter by GCV
  if(length(smoother$rho.grid) >1)
  {
    GCV <- numeric()
    for(s in seq_along(smoother$rho.grid))
    {
      GCV[s] <- sum((y - qr.fitted(smoother$Bt.qr[[s]], y.aug)[1:n] )^2) /(1 - smoother$edf[[s]] / n )^2
    }
    rho.index.opt <- which.min(GCV)
  }else{
    rho.index.opt <- 1
  }

  proj.Vt <- qr.fitted(smoother$Bt.qr[[rho.index.opt]], y.aug)
  if(ortho.constr)
  {
    beta.0.coef <- qr.coef(smoother$B0.qr, proj.Vt[1:n])
    proj.V0 <- drop(smoother$B0 %*% beta.0.coef)
    y.hat <- proj.Vt  - proj.V0
    beta.t.coef <- qr.coef(smoother$Bt.qr[[1]], y.hat)
  }else{
    beta.0.coef <- rep(0, ncol(smoother$B0))
    proj.V0 <- rep(0, n)
    y.hat <- proj.Vt
    beta.t.coef <- qr.coef(smoother$Bt.qr[[1]], y.hat)
  }

  working.resid <-  y - proj.Vt[1:n]
  resid <- y - y.hat[1:n]
  if(ini)
  {
    MSE  <- mean(resid^2)
  }else{
    MSE  <- mean(working.resid^2)
  }

  results <- list(MSE= MSE, resid = resid, working.resid = working.resid,
                  y.hat = y.hat[1:n],
                  working.y.hat = proj.Vt[1:n],
                  proj.V0 = proj.V0,
                  smoother = smoother,
                  beta.t.coef = beta.t.coef, beta.0.coef = beta.0.coef,
                  rho.opt = smoother$rho.grid[rho.index.opt])

  class(results) <- c("gt", "list")
  return(results)
}



######################################################################
# a subfunction to compute the 1st derivatives of the estimated link functions,
# evaluated at the current index u = alpha'X.
# return: d.link.fn, a n x 1 vector of the 1st derivatives of the estimated link functions
# evaluated at the current index u = alpha'X; also, some other related quantities, d.link.fn1 and d.link.fn2.
######################################################################
deriv.link.fn <- function(link.fn.obj)
{

  smoother <- link.fn.obj$smoother
  K <- smoother$K
  u.t <- smoother$u.t
  knots.t <- smoother$knots.t

  t.ind  <- unlist(lapply(1:K, function(x) rep(x, smoother$nbasis.t[-(K+1)][x])))
  beta.t.coef <- split(link.fn.obj$beta.t.coef, t.ind)

  d.design.t <- vector("list", K)
  d.link.fn <- NULL
  for(t in 1:K)
  {
    if(smoother$linear.link.fn)
    {
      d.link.fn <- c(d.link.fn, rep(beta.t.coef[[t]][2], length(u.t[[t]])) )
    }else{
      d.design.t[[t]] <- splineDesign(knots.t[[t]], x=u.t[[t]], derivs=rep(1, length(u.t[[t]])), outer.ok = TRUE)  # compute the 1st derivative of the design functions
      d.link.fn <- c(d.link.fn, d.design.t[[t]] %*% beta.t.coef[[t]])
    }
  }
  rm(d.design.t)
  d.link.fn1 <- - link.fn.obj$working.resid * d.link.fn
  d.link.fn2 <- d.link.fn1^2

  return(list(d.link.fn=d.link.fn, d.link.fn1=d.link.fn1, d.link.fn2=d.link.fn2))
}


######################################################################
# a function to implement an ordinary (not a block-) coordinate descent given the current index u = alpha'X;
# this function can be used to fit an un-regularized constrained single index model.
# return: w.a, the updated (active) coefficients; u, the updated single index; it, the number of iterations
######################################################################
forw.cd <- function(w.a, u, y, Tr, X.a, nbasis.t = c(6,6,8), rho.grid = c(0,0.25,0.5),
                    it.max=60, eps=10^-4, trace=FALSE, linear.link.fn = FALSE, ortho.constr= TRUE, i.fx=NULL)
{

  i.m <- i.fx;
  w.a.old <- w.a + 5;
  it <- 0

  while(max(abs(w.a-w.a.old)) > eps)
  {
    w.a.old <- w.a;
    it <- it + 1;
    if(is.null(i.fx))  i.m <- which.max(w.a);

    for(j in 1:length(w.a))
    {
      if(j!=i.m)
      {
        link.fn.obj <- fit.link.fn.gcv(y, Tr, u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link.fn = linear.link.fn, ortho.constr=ortho.constr)
        d.link.fn.obj <- deriv.link.fn(link.fn.obj)
        gr1 <- crossprod(d.link.fn.obj$d.link.fn1, X.a[,j]);
        hes1 <- crossprod((X.a[,j])^2, d.link.fn.obj$d.link.fn2)
        dj <- -gr1/hes1;
        w.a[j] <- w.a[j] + dj;
        u <- u + dj*X.a[,j]
        if(trace) cat("w.a=",w.a,"\n")
        ####
        if(trace) cat("predictor=", j, "grad= ", gr1, "\n")
      }
    }
    if(it > it.max)
    {
      if(trace) cat("\n","DID NOT CONVERGE","\n");
      break
    }
  }

  list(w.a=w.a, u=u, it=it)
}
#######################################################################################################



######################################################################
## END OF THE FILE
######################################################################
