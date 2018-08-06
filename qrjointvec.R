###################################################################
library(quantreg)
library(splines)
library(kernlab)
library(Matrix)
library(MCMCpack)


qrjointvec <- function(x, y, nsamp = 1e3, thin = 10, cens = NULL,
                    wt = NULL, incr = 0.01, par = "prior", nknots = 6,
                    hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(0.1,0.1,1)),
                    shrink = FALSE, prox.range = c(.2,.95), acpt.target = 0.15,
                    ref.size = 3, blocking = "std5", temp = 1, expo = 2,
                    blocks.mu, blocks.S, fix.nu = FALSE, fbase = c("t", "logistic", "unif"), verbose = TRUE){

    # Set up base functions
    fbase.choice <- match(fbase[1], c("t", "logistic", "unif"))
    if(is.na(fbase.choice)) stop("Only 't', 'logistic' or 'unif' is allowed for the choice of fbase")
    if(fbase.choice == 1){
        q0 <- function(u, nu = Inf) return(1 / (dt(qt(unitFn(u), df = nu), df = nu) * qt(.9, df = nu)))
        Q0 <- function(u, nu = Inf) return(qt(unitFn(u), df = nu) / qt(.9, df = nu))
        F0 <- function(x, nu = Inf) return(pt(x*qt(.9, df = nu), df = nu))
    } else if(fbase.choice == 2){
        fix.nu <- 1
        q0 <- function(u, nu = Inf) return(1 / dlogis(qlogis(unitFn(u))))
        Q0 <- function(u, nu = Inf) return(qlogis(unitFn(u)))
        F0 <- function(x, nu = Inf) return(plogis(x))
    } else {
        fix.nu <- 1
        q0 <- function(u, nu = Inf) return(1 / (dunif(qunif(u, -1,1), -1,1)))
        Q0 <- function(u, nu = Inf) return(qunif(u, -1,1))
        F0 <- function(x, nu = Inf) return(punif(x, -1,1))
    }
    base.bundle <- list(q0 = q0, Q0 = Q0, F0 = F0)
    

    n <- nrow(x)
    p <- ncol(x)
    nresp <- ncol(y)
    x <- scale(x, chull.center(x))
    
    if(is.null(cens)) cens <- rep(0, nresp*n)
    if(is.null(wt)) wt <- rep(1,n)
    
    # tau.g:  sequence of tau probabilities to be spead 'incr' apart, 
    #         supplemented at upper and lower regions by additionally fine grid of taus
    # L:      number of tau locations
    # mid:    indexes prob .5, (e.g. median) location in tau grid
    # reg.ix: indices of tau belonging to non-supplemented grid	

    Ltail <- ceiling(2*log(n,2) + log(incr,2))
    if(Ltail > 0){
        tau.tail <- incr / 2^(Ltail:1)
        tau.g <- c(0, tau.tail, seq(incr, 1 - incr, incr), 1 - tau.tail[Ltail:1], 1)
        L <- length(tau.g); mid <- which.min(abs(tau.g - 0.5)); reg.ix <- (1:L)[-c(1 + 1:Ltail, L - 1:Ltail)]
    } else {
        tau.g <- seq(0, 1, incr)
        L <- length(tau.g); mid <- which.min(abs(tau.g - 0.5)); reg.ix <- (1:L)
    }
    
    tau.kb <- seq(0,1,len = nknots)
    tau.k <- tau.kb
	
    # set priors to defaults if not specified
    a.sig <- hyper$sig; if(is.null(a.sig)) a.sig <- c(.1, .1)
    a.lam <- hyper$lam; if(is.null(a.lam)) a.lam <- c(6, 4)
    a.kap <- hyper$kap; if(is.null(a.kap)) a.kap <- c(0.1, 0.1, 1); a.kap <- matrix(a.kap, nrow = 3); nkap <- ncol(a.kap); a.kap[3,] <- log(a.kap[3,])
    hyper.reduced <- c(a.sig, c(a.kap))
	
    # Create grid-discretized lambdas; retain sufficient overlap to get good mixing.
    # User specifies 1) range for reasonable correlations and 2) parameters for 
    # beta distribution dictating the correlation probabilities.
    # Code translates to lambda scale & gives discretized probs
    prox.grid <- proxFn(max(prox.range), min(prox.range), 0.5)
    ngrid <- length(prox.grid)
    lamsq.grid <- lamFn(prox.grid)^2
    prior.grid <- -diff(pbeta(c(1, (prox.grid[-1] + prox.grid[-ngrid])/2, 0), a.lam[1], a.lam[2]))
    lp.grid <- log(prior.grid)  # log probabilities of the lambda prior
    
    d.kg <- abs(outer(tau.k, tau.g, "-"))^expo	# dist between orig tau grid and reduced tau grid
    d.kk <- abs(outer(tau.k, tau.k, "-"))^expo  # dist between points on reduced tau grid
    gridmats <- matrix(NA, nknots*(L + nknots)+2, ngrid) # setup to hold stuff; (long) x dim of lambda grid
    K0 <- 0
    t1 <- Sys.time()
    for(i in 1:ngrid){
        K.grid <- exp(-lamsq.grid[i] * d.kg); K.knot <- exp(-lamsq.grid[i] * d.kk);	diag(K.knot) <- 1 + 1e-10	
        R.knot <- chol(K.knot); A.knot <- solve(K.knot, K.grid)   # Chol decomp of covariance (not yet inverted)
        # gridmats: Concatenation of first A_0g (L x m), then R_0g (m x m), 
        # then log det of R, then log lambda prior.  Done for each lambda separately	
        gridmats[,i] <- c(c(A.knot), c(R.knot), sum(log(diag(R.knot))), lp.grid[i])	
        # K0: Prior correlation, weighted by prior on different lambda values
        K0 <- K0 + prior.grid[i] * K.knot
    }
    t2 <- Sys.time()
    #cat("Matrix calculation time per 1e3 iterations =", round(1e3 * as.numeric(t2 - t1), 2), "\n")
	
    ## create list of parameters to be passed to c
    # n: number of observations
    # p: number of covariates in design matrix, not inluding intercept
    # L: number of tau grid points
    # mid: index of tau in grid acting as median.  (indexed as one fewer in C than R)
    # nknots: number of knots in interpolation approximation of each w_j
    # ngrid: number of points in lambda grid/discretized prior (lambda dictates
    #        length-scale in GP, ie how related different taus are to eachother)
    
    niter <- nsamp * thin
    dimpars <- c(n, p, L, mid - 1, nknots, ngrid, ncol(a.kap), niter, thin, nresp, nsamp)
    nparmg <- (nknots+1) * (p+1) + 2
    parmg <- matrix(0, nrow = nparmg, ncol = nresp)
    # First supported option: using prior to initialize MC iteration
    if(par[1] == "prior") {
        for(i in 1:nresp){
        if(fix.nu) parmg[nparmg, i] <- nuFn.inv(fix.nu) # nu in last slot of par
        
        # Get rq beta coefficients for each tau in grid. Then use 5th degree b-spline
        # to estimate curve for each beta coefficient
        # "dither" slightly perturbs responses in order to avoid potential degeneracy
        # of dual simplex algorithm (present with large number of ties in y) and "hanging"
        # of the rq Fortran code
        beta.rq <- sapply(tau.g, function(a) return(coef(suppressWarnings(rq(dither(y[,i]) ~ x, tau = a, weights = wt)))))   # THIS COULD BE DONE WITH tau=tau.g
        v <- bs(tau.g, df = 5)
		
        # over tau and per coefficient get smoothed fits through data (ie independently estimated quantiles);
        # ultimate goal: reasonable estimates at median, at delta and at 1-delta
        rq.lm <- apply(beta.rq, 1, function(z) return(coef(lm(z ~ v))))
		
        # Combine b-spline estimates to get estimate of coefficient when tau
        # is delta (lower quantile) and when tau is 1 - delta (upper quantile)
        delta <- tau.g[2]
        tau.0 <- tau.g[mid]
        rq.tau0 <- c(c(1, predict(v, tau.0)) %*% rq.lm)
        rq.delta <- c(c(1, predict(v, delta)) %*% rq.lm)
        rq.deltac <- c(c(1, predict(v, 1 - delta)) %*% rq.lm)
        
        parmg[nknots*(p+1) + 1:(p+1),i] <- as.numeric(rq.tau0) # median coefficient quantiles  placed in p + 1 positions just before last 2
        sigma <- 1
        parmg[(nknots+1)*(p+1) + 1,i] <- sigFn.inv(sigma, a.sig) # prior for sigma.  
        
        for(j in 1:(p+1)){
            kapsq <- sum(exp(a.kap[3,]) * (a.kap[2,] / a.kap[1,]))
            lam.ix <- sample(length(lamsq.grid), 1, prob = prior.grid)
            R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),lam.ix], nknots, nknots)
            z <- sqrt(kapsq) * c(crossprod(R, rnorm(nknots)))
            parmg[(j - 1) * nknots + 1:nknots,i] <- z - mean(z) # centered and placed in sets of length nknots one after another for p+1 covariates
        }
        
        # get betas associated with these prior params & corresponding quantile estimate (centered around median)
        beta.hat <- estFn(parmg, x, y[,i], gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, FALSE, base.bundle = base.bundle)
        qhat <- tcrossprod(cbind(1, x), beta.hat)
        
        infl <- max(max((y[,i] - qhat[,mid])/(qhat[,ncol(qhat) - 1] - qhat[,mid])), max((qhat[,mid] - y[,i])/(qhat[,mid] - qhat[,2])))
        oo <- .C("INIT", par = as.double(parmg[,i]), x = as.double(x), y = as.double(y[,i]), cens = as.integer(cens[(i-1)*n+1:n]), wt = as.double(wt),
        shrink = as.integer(shrink), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridpars = as.double(gridmats), 
        tau.g = as.double(tau.g), siglim = as.double(sigFn.inv(c(1.0 * infl * sigma, 10 * infl * sigma), a.sig)),
        fbase.choice = as.integer(fbase.choice))
        
        parmg[,i] <- oo$par
        }
    } else if (par[1] == "RQ"){
		
        # Initialize at a model space approximation of the estimates from rq
        par <- rep(0, (nknots+1) * (p+1) + 2)
        
        beta.rq <- sapply(tau.g, function(a) return(coef(suppressWarnings(rq(dither(y)~ x, tau = a, weights = wt)))))
        v <- bs(tau.g, df = 5)
        rq.lm <- apply(beta.rq, 1, function(z) return(coef(lm(z ~ v)))) # smooth though non-monotonic estimates
        
        delta <- tau.g[2]     # smallest but one
        tau.0 <- tau.g[mid]   # median or mid-most value of taus evaluated
        rq.tau0 <- c(c(1, predict(v, tau.0)) %*% rq.lm)       # estimate coef median MLEs 
        rq.delta <- c(c(1, predict(v, delta)) %*% rq.lm)      # estimate coef upper quantile MLEs
        rq.deltac <- c(c(1, predict(v, 1 - delta)) %*% rq.lm) # estimate coef lower quantile MLEs
        
        par[nknots*(p+1) + 1:(p+1)] <- as.numeric(rq.tau0)
        nu <- ifelse(fix.nu, fix.nu, nuFn(0))
        # Solving for sigma in zeta(tau) - F_0( (beta0(delta) - gam0)/sigma)
        sigma <- min((rq.delta[1] - rq.tau0[1]) / Q0(delta, nu), (rq.deltac[1] - rq.tau0[1]) / Q0(1 - delta, nu))
        par[(nknots+1)*(p+1) + 1]  <- sigFn.inv(sigma, a.sig)
        
        epsilon <- 0.1 * min(diff(sort(tau.k)))
        tau.knot.plus <- pmin(tau.k + epsilon, 1)
        tau.knot.minus <- pmax(tau.k - epsilon, 0)
        beta.rq.plus <- cbind(1, predict(v, tau.knot.plus)) %*% rq.lm
        beta.rq.minus <- cbind(1, predict(v, tau.knot.minus)) %*% rq.lm
        zeta0.plus <- F0((beta.rq.plus[,1] - rq.tau0[1]) / sigma, nu)
        zeta0.minus <- F0((beta.rq.minus[,1] - rq.tau0[1]) / sigma, nu)
        zeta0.dot.knot <- (zeta0.plus - zeta0.minus) / (tau.knot.plus - tau.knot.minus)
        w0.knot <- log(pmax(epsilon, zeta0.dot.knot)) / shrinkFn(p)
        w0.knot <- (w0.knot - mean(w0.knot))
        
        w0PP <- ppFn0(w0.knot, gridmats, L, nknots, ngrid)
        w0 <- w0PP$w
        
        zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))
        zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2)
        zeta0.tot <- zeta0[L-2]
        zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
        zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
        zeta0.dot[c(1,L)] <- 0
        
        par[1:nknots] <- w0.knot
        beta0.dot <- sigma * q0(zeta0, nu) * zeta0.dot
        
        tilt.knot <- tau.g[tilt.ix <- sapply(tau.k, function(a) which.min(abs(a - zeta0)))]
        
        tau.knot.plus <- pmin(tilt.knot + epsilon, 1)
        tau.knot.minus <- pmax(tilt.knot - epsilon, 0)
        beta.rq.plus <- cbind(1, predict(v, tau.knot.plus)) %*% rq.lm
        beta.rq.minus <- cbind(1, predict(v, tau.knot.minus)) %*% rq.lm
        
        beta.dot.knot <- (beta.rq.plus[,-1,drop=FALSE] - beta.rq.minus[,-1,drop=FALSE]) /  (tau.knot.plus - tau.knot.minus)
        
        par[nknots + 1:(nknots*p)] <- c(beta.dot.knot)
        beta.hat <- estFn(par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, FALSE, base.bundle = base.bundle)
        qhat <- tcrossprod(cbind(1, x), beta.hat)
        infl <- max(max((y - qhat[,mid])/(qhat[,ncol(qhat) - 1] - qhat[,mid])), max((qhat[,mid] - y)/(qhat[,mid] - qhat[,2])))
        
        oo <- .C("INIT", par = as.double(par), x = as.double(x), y = as.double(y), cens = as.integer(cens),
        wt = as.double(wt), shrink = as.integer(shrink), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridpars = as.double(gridmats),
        tau.g = as.double(tau.g), siglim = as.double(sigFn.inv(c(1.0 * infl * sigma, 10.0 * infl * sigma), a.sig)), fbase.choice = as.integer(fbase.choice))
        
        par <- oo$par
    } else if (par[1] == "noX") { # Initialization without regard to X, will definitely work
        fit0 <- qde(y, nsamp = 1e2, thin = 10, cens = cens, wt = wt, nknots = nknots, hyper = hyper, prox.range = prox.range, fbase = fbase, fix.nu = fix.nu, verbose = FALSE)
        par <- rep(0, (nknots + 1) * (p + 1) + 2)
        par[c(1:nknots, nknots * (p + 1) + 1, (nknots + 1) * (p + 1) + 1:2)] <- fit0$par
    }
    
	# Choose a blocking scheme for updates.
    if(blocking == "single"){                                                # ONE BLOCK
        blocks <- list(rep(TRUE, nparmg))                                        # Update all simultaneously
    } else if(blocking == "single2"){                                        # TWO BLOCKS
        blocks <- list(rep(TRUE, nparmg), rep(FALSE, nparmg))                      # First, update all simultaneously
        blocks[[2]][nknots*(p+1) + 1:(p+3)] <- TRUE                            # Second, update lambda_nots, sigma, and nu simultaneously
    } else if(blocking == "single3"){                                        # THREE BLOCKS
        blocks <- list(rep(TRUE, nparmg), rep(FALSE, nparmg), rep(FALSE, nparmg))    # First, update all simultaneously
        blocks[[2]][nknots*(p+1) + 1:(p+1)] <- TRUE                            # Second, update all lambda_nots
        blocks[[3]][(nknots+1)*(p+1) + 1:2] <- TRUE                            # Last, update sigma and nu
    } else if(blocking == "std0"){                                           # P + 1 BLOCKS
        blocks <- replicate(p + 1, rep(FALSE, nparmg), simplify = FALSE)         # Update everything related to covariate (GP, lamda_not) simulatneous with sigma, nu
        for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
    } else if(blocking == "std1"){                                           # P + 2 BLOCKS
        blocks <- replicate(p + 2, rep(FALSE, nparmg), simplify = FALSE)         # First, Lowrank GPS + sigma, nu for each covariate
        for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
        blocks[[p + 2]][nknots*(p+1) + 1:(p+3)] <- TRUE                        # Then all covariate medians, sigma, nu simultaneously 
    } else if(blocking == "std2"){                                           # P + 3 BLOCKS
        blocks <- replicate(p + 3, rep(FALSE, nparmg), simplify = FALSE)         # Just like P+2, with additional update that is only simga and nu
        for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
        blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
        blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
    } else if(blocking == "std3"){                                           # ALSO P + 3 BLOCKS
        blocks <- replicate(p + 3, rep(FALSE, nparmg), simplify = FALSE)
        for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots)] <- TRUE        # Update GPs related to covariates
        blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE                        # Update covariate medians
        blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE                        # Lastly, update sigma and nu
    } else if(blocking == "std4"){
        blocks <- replicate(p + 3, rep(FALSE, nparmg), simplify = FALSE)         # ALSO P + 3 BLOCKS
        for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1)] <- TRUE  # Update GPs and medians
        blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE                        # Update medians
        blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE                        # Update sigma and nu
    } else if(blocking == "std5"){
        blocks <- replicate(p + 4, rep(FALSE, nparmg), simplify = FALSE)        # Same as previous with additional update to ALL parameters
        for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1)] <- TRUE
        blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
        blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
        blocks[[p + 4]][1:nparmg] <- TRUE
    } else {                                                                # Univariate updates
        blocks <- replicate(nparmg, rep(FALSE, nparmg), simplify = FALSE)
        for(i in 1:nparmg) blocks[[i]][i] <- TRUE
    }
    
    nblocks <- length(blocks)
    if(fix.nu) for(j in 1:nblocks) blocks[[j]][(nknots+1) * (p+1) + 2] <- FALSE
    
    blocks.ix <- rep(c(unlist(lapply(blocks, which))) - 1, nresp)           # Index of locations of updates
    blocks.size <- sapply(blocks, sum)                                      # Size of MVN in each block update
    if(missing(blocks.mu)) blocks.mu <- rep(0, nresp*sum(blocks.size))
    if(missing(blocks.S)){
        sig.nu <- c(TRUE, !as.logical(fix.nu))
        blocks.S <- rep(list(lapply(blocks.size, function(q) diag(1, q))), nresp)
        if(substr(blocking, 1, 3) == "std"){
          for(i in 1:nresp) {for(j in 1:(p+1)) {blocks.S[[i]][[j]][1:nknots, 1:nknots] <- K0}}
            if(as.numeric(substr(blocking, 4,5)) > 1){
              for(i in 1:nresp){
                blocks.S[[i]][[p + 2]] <- summary(suppressWarnings(rq(dither(y[,i]) ~ x, tau = 0.5, weights = wt)), se = "boot", cov = TRUE)$cov
                blocks.S[[i]][[p + 3]] <- matrix(c(1, 0, 0, .1), 2, 2)[sig.nu,sig.nu]
              }
            }
            if(as.numeric(substr(blocking, 4,5)) == 5){
                slist <- list(); length(slist) <- p + 3
                for(i in 1:(p+1)) slist[[i]] <- K0
                slist[[p+3]] <- matrix(c(1, 0, 0, .1), 2, 2)[sig.nu,sig.nu]
                for(i in 1:nresp){
                slist[[p+2]] <- summary(suppressWarnings(rq(dither(y[,i]) ~ x, tau = 0.5, weights = wt)), se = "boot", cov = TRUE)$cov
                blocks.S[[i]][[p + 4]] <- as.matrix(bdiag(slist))
                }
            }
        }
        blocks.S <- unlist(blocks.S)
    }
    
    imcmc.par <- c(nblocks, ref.size, verbose, max(10, niter/1e4), rep(0, nresp*nblocks))
    dmcmc.par <- c(temp, 0.999, rep(acpt.target, nresp*nblocks), rep(2.38 / sqrt(blocks.size), nresp))
    
    # Updated in C code
    # parsamp:  holds all parameters from all iterations in one long vector
    # acptsamp: holds all acceptance rates for each block at each iteration in one long vector
    # lpsamp:   holds all evaluations of log posterior at each iteration 
    
    C <- diag(nresp)    

    tm.c <- system.time(oo <- .C("BJQRvec", par = as.double(parmg), C = as.double(C), x = as.double(t(x)), y = as.double(y), cens = as.integer(cens), wt = as.double(wt),
    shrink = as.integer(shrink), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridmats = as.double(gridmats),
    tau.g = as.double(tau.g), muV = as.double(blocks.mu), SV = as.double(blocks.S), blocks = as.integer(blocks.ix), 
    blocks.size = as.integer(blocks.size), dmcmcpar = as.double(dmcmc.par), 
    imcmcpar = as.integer(imcmc.par), parsamp = double(nsamp * nresp*nparmg), 
    acptsamp = double(nsamp *nresp* nblocks), lpsamp = double(nsamp), Csamp = double(nsamp*nresp*(nresp-1)/2), acptCsamp = double(nsamp), lcopsumsamp = double(nsamp), fbase.choice = as.integer(fbase.choice)))
    if(verbose) cat("elapsed time:", round(tm.c[3]), "seconds\n")
    
   oo$x <- x; oo$y <- y; # oo$xnames <- colnames(x); oo$terms <- mt;
   oo$gridmats <- gridmats; oo$prox <- prox.grid; oo$reg.ix <- reg.ix;
   oo$runtime <- tm.c[3]
   class(oo) <- "qrjointvec"
   return(oo)
}

#########################################################################

coef.qrjointvec <- function(object, burn.perc = 0.5, nmc = 200, plot = FALSE, show.intercept = TRUE, reduce = TRUE, ...){
  nsamp <- object$dim[11]                      # number of iterations retains in sample
  nresp <- object$dim[10]
  n <- object$dim[1]; p <- object$dim[2]; L <- object$dim[3]; mid <- object$dim[4] + 1; nknots <- object$dim[5]; ngrid <- object$dim[6]
  nparmg <- (p + 1)*(nknots+1) + 2
  pars <- vector('list', nresp)
  pars.mat <- matrix(object$parsamp, nrow = nresp*nsamp, byrow = T) # reorder parameters into npar x nsamp matrix
  for(i in 1:nresp){
    pars[[i]] <- t(pars.mat[seq(i, nresp*nsamp, by = nresp),])
  }
  
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1])) # indices over which to summarizing (exclude burn; keep nmc samples)
  
  a.sig <- object$hyper[1:2]; a.kap <- matrix(object$hyper[-c(1:2)], nrow = 3)
  tau.g <- object$tau.g; reg.ix <- object$reg.ix
  x.ce <- outer(rep(1, L), attr(object$x, "scaled:center")); x.sc <- outer(rep(1,L), attr(object$x, "scaled:scale"))
  
  base.bundle <- list()
  if(object$fbase.choice == 1){
    base.bundle$q0 <- function(u, nu = Inf) return(1 / (dt(qt(unitFn(u), df = nu), df = nu) * qt(.9, df = nu)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qt(unitFn(u), df = nu) / qt(.9, df = nu))
    base.bundle$F0 <- function(x, nu = Inf) return(pt(x*qt(.9, df = nu), df = nu))
  } else if(object$fbase.choice == 2){
    base.bundle$q0 <- function(u, nu = Inf) return(1 / dlogis(qlogis(unitFn(u))))
    base.bundle$Q0 <- function(u, nu = Inf) return(qlogis(unitFn(u)))
    base.bundle$F0 <- function(x, nu = Inf) return(plogis(x))
  } else {
    base.bundle$q0 <- function(u, nu = Inf) return(1 / (dunif(qunif(u, -1,1), -1,1)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qunif(u, -1,1))
    base.bundle$F0 <- function(x, nu = Inf) return(punif(x, -1,1))
  }
  
  # use estFn to turn posterior on parameters into posterior betas;
  # return 3-D array { L x (p+1) x nsim }
  beta.samp <- vector('list', nresp)
  for(i in 1:nresp)
  {beta.samp[[i]] <- sapply(ss, function(p1) estFn(pars[[i]][,p1], object$x, object$y[,i], object$gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce, x.ce, x.sc, base.bundle), simplify="array")}
  
  if(reduce) tau.g <- tau.g[reg.ix]
  L <- length(tau.g)
  for(i in 1:nresp){
  if(plot){
    nr <- ceiling(sqrt(p+show.intercept)); nc <- ceiling((p+show.intercept)/nr)
    cur.par <- par(no.readonly = TRUE)
    par(mfrow = c(nr, nc))
  }
  
  # store summaries in 3-D array {L x (p+1) x 3}

  beta.hat <- array(0,c(L,p+1,3))
  plot.titles <- c("Intercept", paste('X', 1:p, sep=''))
  beta.hat[,1,] <- getBands(beta.samp[[i]][,1,], plot = (plot & show.intercept), add = FALSE, x = tau.g, xlab = "tau", ylab = "Coefficient", bty = "n")
  if(plot & show.intercept) {title(main = plot.titles[1]); points(seq(0,1,0.001), tau.intercept(seq(0,1,0.001)), type = 'l')}
  for(j in 2:(p+1)){
    beta.hat[,j,] <- getBands(beta.samp[[i]][,j,], plot = plot, add = FALSE, x = tau.g, xlab = "tau", ylab = "Coefficient", bty = "n")
    if(plot) {
      title(main = plot.titles[j])
      abline(h = 0, lty = 2, col = 4)
      points(seq(0,1,0.001), tau.slope(seq(0,1,0.001)), type = 'l')
    }
  }	
  if(plot) suppressWarnings(par(cur.par,no.readonly = TRUE))  # return R parameters to pre-plot settings
  dimnames(beta.hat) <- list(tau=tau.g, beta=plot.titles, summary=c("b.lo", "b.med", "b.hi"))
  dimnames(beta.samp[[i]]) <- list(tau=tau.g, beta=plot.titles, iter=1:length(ss))
  invisible(list(beta.samp = beta.samp[[i]], beta.est = beta.hat))
  
  mid.red <- which(tau.g == object$tau.g[mid])
  parametric.list <- rbind(beta.samp[[i]][mid.red, , ,drop=TRUE],
                           sigma = sigFn(pars[[i]][nknots * (p+1) + (p+1) + 1,ss], a.sig),	
                           nu = nuFn(pars[[i]][(nknots + 1) * (p+1) + 2,ss]))	
  dimnames(parametric.list)[[1]][1 + 0:p] <- c("Intercept", paste('X', 1:p, sep=''))	
  gamsignu <- t(apply(parametric.list, 1, quantile, pr = c(0.5, 0.025, 0.975)))	
  dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")	
  invisible(list(beta.samp = beta.samp[[i]], beta.est = beta.hat, parametric = gamsignu))
}
}
#########################################################################



coef.qde <- function(object, burn.perc = 0.5, nmc = 200, reduce = TRUE, ...){
    niter <- object$dim[7]
    nsamp <- object$dim[9]
    pars <- matrix(object$parsamp, ncol = nsamp)
    ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
    
    n <- object$dim[1]; p <- 0; L <- object$dim[2]; mid <- object$dim[3] + 1; nknots <- object$dim[4]; ngrid <- object$dim[5]
    a.sig <- object$hyper[1:2]; a.kap <- matrix(object$hyper[-c(1:2)], nrow = 3)
    tau.g <- object$tau.g; reg.ix <- object$reg.ix
    
    base.bundle <- list()
    if(object$fbase.choice == 1){
        base.bundle$q0 <- function(u, nu = Inf) return(1 / (dt(qt(unitFn(u), df = nu), df = nu) * qt(.9, df = nu)))
        base.bundle$Q0 <- function(u, nu = Inf) return(qt(unitFn(u), df = nu) / qt(.9, df = nu))
        base.bundle$F0 <- function(x, nu = Inf) return(pt(x*qt(.9, df = nu), df = nu))
    } else if(object$fbase.choice == 2){
        base.bundle$q0 <- function(u, nu = Inf) return(1 / dlogis(qlogis(unitFn(u))))
        base.bundle$Q0 <- function(u, nu = Inf) return(qlogis(unitFn(u)))
        base.bundle$F0 <- function(x, nu = Inf) return(plogis(x))
    } else {
        base.bundle$q0 <- function(u, nu = Inf) return(1 / (dunif(qunif(u, -1,1), -1,1)))
        base.bundle$Q0 <- function(u, nu = Inf) return(qunif(u, -1,1))
        base.bundle$F0 <- function(x, nu = Inf) return(punif(x, -1,1))
    }
    
    beta.samp <- apply(pars[,ss], 2, function(p1) c(estFn.noX(p1, object$y, object$gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce, base.bundle)))
    
    if(reduce) tau.g <- tau.g[reg.ix]
    L <- length(tau.g)
    b <- beta.samp[1:L,]
    beta.hat <- getBands(b, plot = FALSE, add = FALSE, x = tau.g, ...)

    parametric.list <- list(gam0 = pars[nknots + 1,ss], sigma = sigFn(pars[nknots + 2,ss], a.sig), nu = nuFn(pars[nknots + 3,ss]))
    gamsignu <- t(sapply(parametric.list, quantile, pr = c(0.5, 0.025, 0.975)))
    dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")
    invisible(list(beta.samp = beta.samp, beta.est = beta.hat, parametric = gamsignu))
}


#########################################################################

# Function creates a set of plots to assess convergence of markov chains

summary.qrjoint <- function(object, ntrace = 1000, burn.perc = 0.5, plot.dev = TRUE, more.details = FALSE, ...){
	thin <- object$dim[9]	
	nsamp <- object$dim[10]
	pars <- matrix(object$parsamp, ncol = nsamp)
	ss <- unique(pmax(1, round(nsamp * (1:ntrace/ntrace))))
    post.burn <- (ss > nsamp * burn.perc)
	dimpars <- object$dim
	dimpars[8] <- length(ss)
	
	n <- object$dim[1]; p <- object$dim[2]; ngrid <- object$dim[6]
	  # Calcuate deviance
    sm <- .C("DEV", pars = as.double(pars[,ss]), x = as.double(object$x), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt),
			 shrink = as.integer(object$shrink), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g),
			 devsamp = double(length(ss)), llsamp = double(length(ss)*n), pgsamp = double(length(ss)*ngrid*(p+1)), rpsamp = double(length(ss)*n),fbase.choice = as.integer(object$fbase.choice))
	deviance <- sm$devsamp
	ll <- matrix(sm$llsamp, ncol = length(ss))
	rp <- matrix(sm$rpsamp, ncol = length(ss))
	fit.waic <- waic(ll[,post.burn])
	pg <- matrix(sm$pgsamp, ncol = length(ss))
	prox.samp <- matrix(NA, p+1, length(ss))
	for(i in 1:(p+1)){
		prox.samp[i,] <- object$prox[apply(pg[(i-1)*ngrid + 1:ngrid,], 2, function(pr) sample(length(pr), 1, prob = pr))]
	}	
	
	# Trace plots for beta coefficients
	if(more.details) {
	  cur.par <- par(no.readonly = TRUE)
	  par(mfrow = c(2,2), mar = c(5,4,3,2)+.1)
	  }
	if(plot.dev){
		plot(thin * ss, deviance, ty = "l", xlab = "Markov chain iteration", ylab = "Deviance", bty = "n", main = "Fit trace plot", ...)
		grid(col = "gray")
	}
	
	if(more.details){
		ngrid <- length(object$prox)
		prior.grid <- exp(object$gridmats[nrow(object$gridmats),])
		lam.priorsamp <- lamFn(sample(object$prox, ntrace, replace = TRUE, prob = prior.grid))
		lam.prior.q <- quantile(lam.priorsamp, pr = c(.025, .5, .975))
		lam.samp <- lamFn(prox.samp)
		a <- min(lamFn(object$prox))
		b <- diff(range(lamFn(object$prox))) * 1.2
		plot(thin * ss, lam.samp[1,], ty = "n", ylim = a + c(0, b * (p + 1)), bty = "n", ann = FALSE, axes = FALSE)
		axis(1)
		for(i in 1:(p+1)){
			abline(h = b * (i-1) + lamFn(object$prox), col = "gray")
			abline(h = b * (i - 1) + lam.prior.q, col = "red", lty = c(2,1,2))
			lines(thin * ss, b * (i-1) + lam.samp[i,], lwd = 1, col = 4)
			if(i %% 2) axis(2, at = b * (i-1) + lamFn(object$prox[c(1,ngrid)]), labels = round(object$prox[c(1,ngrid)],2), las = 1, cex.axis = 0.6) 
			mtext(substitute(beta[index], list(index = i - 1)), side = 4, line = 0.5, at = a + b * (i - 1) + 0.4*b, las = 1)			
		}
		title(xlab = "Markov chain iteration", ylab = "Proxmity posterior", main = "Mixing over GP scaling")					
		
		# Plot of geweke tests for convergence on all parameters estimated; includes
		# Benjamini Hochberg line for false discovery rates
		theta <- as.mcmc(t(matrix(object$parsamp, ncol = nsamp)[,ss[post.burn]]))
		gg <- geweke.diag(theta, .1, .5)
		zvals <- gg$z

		pp <- 2 * (1 - pnorm(abs(zvals)))
		plot(sort(pp), ylab = "Geweke p-values", xlab = "Parameter index (reordered)", main = "Convergence diagnosis", ty = "h", col = 4, ylim = c(0, 0.3), lwd = 2)
		abline(h = 0.05, col = 2, lty = 2)
		abline(a = 0, b = 0.1 / length(pp), col = 2, lty = 2)
		mtext(c("BH-10%", "5%"), side = 4, at = c(0.1, 0.05), line = 0.1, las = 0, cex = 0.6)
		
		npar <- length(object$par)
		image(1:npar, 1:npar, cor(theta), xlab = "Parameter index", ylab = "Parameter index", main = "Parameter correlation")
		suppressWarnings(par(cur.par,no.readonly = TRUE))
	}
	invisible(list(deviance = deviance, pg = pg, prox = prox.samp, ll = ll, rp = rp, waic = fit.waic))
}

summary.qde <- function(object, ntrace = 1000, burn.perc = 0.5, plot.dev = TRUE, more.details = FALSE, ...){
    thin <- object$dim[8]
    nsamp <- object$dim[9]
    pars <- matrix(object$parsamp, ncol = nsamp)
    ss <- unique(pmax(1, round(nsamp * (1:ntrace/ntrace))))
    post.burn <- (ss > nsamp * burn.perc)
    dimpars <- object$dim
    dimpars[7] <- length(ss)
    
    n <- object$dim[1]; p <- 0; ngrid <- object$dim[5]
    sm <- .C("DEV_noX", pars = as.double(pars[,ss]), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt),
			 hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g),
             devsamp = double(length(ss)), llsamp = double(length(ss)*n), pgsamp = double(length(ss)*ngrid), rpsamp=double(length(ss)*n),fbase.choice = as.integer(object$fbase.choice))
    deviance <- sm$devsamp
    ll <- matrix(sm$llsamp, ncol = length(ss))
    rp <- matrix(sm$rpsamp, ncol = length(ss))
    fit.waic <- waic(ll[,post.burn])
    pg <- matrix(sm$pgsamp, ncol = length(ss))
    prox.samp <- object$prox[apply(pg[1:ngrid,], 2, function(pr) sample(length(pr), 1, prob = pr))]
    
    if(more.details) {
      cur.par <- par(no.readonly = TRUE)
      par(mfrow = c(2,2), mar = c(5,4,3,2)+.1)
    }
    if(plot.dev){
        plot(thin * ss, deviance, ty = "l", xlab = "Markov chain iteration", ylab = "Deviance", bty = "n", main = "Fit trace plot", ...)
        grid(col = "gray")
    }
    
    if(more.details){
        ngrid <- length(object$prox)
        prior.grid <- exp(object$gridmats[nrow(object$gridmats),])
        lam.priorsamp <- lamFn(sample(object$prox, ntrace, replace = TRUE, prob = prior.grid))
        lam.prior.q <- quantile(lam.priorsamp, pr = c(.025, .5, .975))
        lam.samp <- lamFn(prox.samp)
        a <- min(lamFn(object$prox))
        b <- diff(range(lamFn(object$prox))) * 1.2
        plot(thin * ss, lam.samp, ty = "n", ylim = a + c(0, b), bty = "n", ann = FALSE, axes = FALSE)
        axis(1)
        for(i in 1:1){
            abline(h = b * (i-1) + lamFn(object$prox), col = "gray")
            abline(h = b * (i - 1) + lam.prior.q, col = "red", lty = c(2,1,2))
            lines(thin * ss, b * (i-1) + lam.samp, lwd = 1, col = 4)
            if(i %% 2) axis(2, at = b * (i-1) + lamFn(object$prox[c(1,ngrid)]), labels = round(object$prox[c(1,ngrid)],2), las = 1, cex.axis = 0.6)
            mtext(substitute(beta[index], list(index = i - 1)), side = 4, line = 0.5, at = a + b * (i - 1) + 0.4*b, las = 1)
        }
        title(xlab = "Markov chain iteration", ylab = "Proxmity posterior", main = "Mixing over GP scaling")
        
        theta <- as.mcmc(t(matrix(object$parsamp, ncol = nsamp)[,ss[post.burn]]))
        gg <- geweke.diag(theta, .1, .5)
        zvals <- gg$z
        
        pp <- 2 * (1 - pnorm(abs(zvals)))
        plot(sort(pp), ylab = "Geweke p-values", xlab = "Parameter index (reordered)", main = "Convergence diagnosis", ty = "h", col = 4, ylim = c(0, 0.3), lwd = 2)
        abline(h = 0.05, col = 2, lty = 2)
        abline(a = 0, b = 0.1 / length(pp), col = 2, lty = 2)
        mtext(c("BH-10%", "5%"), side = 4, at = c(0.1, 0.05), line = 0.1, las = 0, cex = 0.6)
        
        npar <- length(object$par)
        image(1:npar, 1:npar, cor(theta), xlab = "Parameter index", ylab = "Parameter index", main = "Parameter correlation")
        suppressWarnings(par(cur.par,no.readonly = TRUE)) 
    }
    invisible(list(deviance = deviance, pg = pg, prox = prox.samp, ll = ll, rp=rp, waic = fit.waic))
}

#########################################################################
# A function to calculate sampled predictions on original or new data set
# object:  Object of class inheriting from "qrjoint"
# newdata: An optional data frame containing variables on which to predict. If omitted,
#          the fitted values are used.
# at:      (optional) tau.grid values (e.g. .5, .9) to keep and return
# summarize:  Logical - medians of MC posterior


predict.qrjoint <- function(object, newdata=NULL, summarize=TRUE, burn.perc = 0.5, nmc = 200, ...){
  p <- object$dim[2];
  betas <- coef(object, burn.perc=burn.perc, nmc=nmc, plot=FALSE)
  nsamp <- dim(betas$beta.samp)[3]
  L <- dim(betas$beta.samp)[1]
  
  if(is.null(newdata)) {# if predicting on original data, restore design matrix to original center and scale
    Xpred <- cbind(1, sapply(1:p, function(r) object$x[,r]*attr(object$x,'scaled:scale')[r] + attr(object$x, 'scaled:center')[r]))
  } else{
    
    Xpred <- model.matrix(object$terms, data=newdata)
    # Would be good to add error catching for NA's in newdata dataframe
  }
  npred <- dim(Xpred)[1]
  
  pred <- array(NA, c(npred, L, nsamp))
  for (i in 1:nsamp){ pred[,,i] <- tcrossprod(Xpred, betas$beta.samp[,,i]) }
  dimnames(pred) <- list(obs=rownames(Xpred), tau=round(object$tau.g[object$reg.ix],4), samp=1:nsamp)
  
  # Posterior median as estimate for each observation at each tau
  if(summarize){ 
    pred <- apply(pred,c(1,2), quantile, probs=.5)
  }
  return(pred)
} 

predict.qde <- function(object, burn.perc = 0.5, nmc = 200, yRange = range(object$y), yLength = 401, ...){
    thin <- object$dim[8]
    nsamp <- object$dim[9]
    pars <- matrix(object$parsamp, ncol = nsamp)
    ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
    dimpars <- object$dim
    dimpars[7] <- length(ss)

    yGrid <- seq(yRange[1], yRange[2], len = yLength)
    n <- object$dim[1]; p <- 0; ngrid <- object$dim[5]
    dimpars[1] <- yLength
    
    pred <- .C("PRED_noX", pars = as.double(pars[,ss]), yGrid = as.double(yGrid), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), ldenssamp = double(length(ss)*yLength),fbase.choice = as.integer(object$fbase.choice))
    dens <- matrix(exp(pred$ldenssamp), ncol = length(ss))
    return(list(y = yGrid, fsamp = dens, fest = t(apply(dens, 1, quantile, pr = c(.025, .5, .975)))))
}


# Function to construct beta covariate coefficient functionals from posterior outputs

estFn <- function(par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce = TRUE, x.ce = 0, x.sc = 1, base.bundle){
	
	n <- length(y); p <- ncol(x)
	wKnot <- matrix(par[1:(nknots*(p+1))], nrow = nknots)
	w0PP <- ppFn0(wKnot[,1], gridmats, L, nknots, ngrid)
	w0 <- w0PP$w
	wPP <- apply(wKnot[,-1,drop=FALSE], 2, ppFn, gridmats = gridmats, L = L, nknots = nknots, ngrid = ngrid, a.kap = a.kap)
	wMat <- matrix(sapply(wPP, extract, vn = "w"), ncol = p)
	
	zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))          # subract max for numerical stability
	zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2) # take integral of derivative to get original function.
	zeta0.tot <- zeta0[L-2] 
	zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
	zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
	zeta0.dot[c(1,L)] <- 0
    zeta0.ticks <- pmin(L-1, pmax(1, sapply(zeta0, function(u) sum(tau.g <= u))))
    zeta0.dists <- (zeta0 - tau.g[zeta0.ticks]) / (tau.g[zeta0.ticks+1] - tau.g[zeta0.ticks])
    vMat <- apply(wMat, 2, transform.grid, ticks = zeta0.ticks, dists = zeta0.dists)
	
	reach <- nknots*(p+1)
	gam0 <- par[reach + 1]; reach <- reach + 1
	gam <- par[reach + 1:p]; reach <- reach + p
	sigma <- sigFn(par[reach + 1], a.sig); reach <- reach + 1
	nu <- nuFn(par[reach + 1]);
	
	# Intercept quantiles
	b0dot <- sigma * base.bundle$q0(zeta0, nu) * zeta0.dot
	beta0.hat <- rep(NA, L)
	beta0.hat[mid:L] <- gam0 + trape(b0dot[mid:L], tau.g[mid:L], L - mid + 1)
	beta0.hat[mid:1] <- gam0 + trape(b0dot[mid:1], tau.g[mid:1], mid)
	
	# Working towards other covariate coefficient quantiles
	# First four lines get 'aX' shrinkage factor
	vNorm <- sqrt(rowSums(vMat^2))
	a <- tcrossprod(vMat, x)
	aX <- apply(-a, 1, max)/vNorm
	aX[is.nan(aX)] <- Inf
	aTilde <- vMat / (aX * sqrt(1 + vNorm^2))
	ab0 <- b0dot * aTilde
	
	beta.hat <- kronecker(rep(1,L), t(gam))
	beta.hat[mid:L,] <- beta.hat[mid:L,] + apply(ab0[mid:L,,drop=FALSE], 2, trape, h = tau.g[mid:L], len = L - mid + 1)
	beta.hat[mid:1,] <- beta.hat[mid:1,] + apply(ab0[mid:1,,drop=FALSE], 2, trape, h = tau.g[mid:1], len = mid)
	beta.hat <- beta.hat / x.sc
	beta0.hat <- beta0.hat - rowSums(beta.hat * x.ce)
	betas <- cbind(beta0.hat, beta.hat)
	if(reduce) betas <- betas[reg.ix,,drop = FALSE]
	return(betas)
}

estFn.noX <- function(par, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce = TRUE, base.bundle){
    
    n <- length(y); p <- 0
    wKnot <- par[1:nknots]
    w0PP <- ppFn0(wKnot, gridmats, L, nknots, ngrid)
    w0 <- w0PP$w
    
    zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))
    zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2)
    zeta0.tot <- zeta0[L-2]
    zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
    zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
    zeta0.dot[c(1,L)] <- 0
    zeta0.ticks <- pmin(L-1, pmax(1, sapply(zeta0, function(u) sum(tau.g <= u))))
    zeta0.dists <- (zeta0 - tau.g[zeta0.ticks]) / (tau.g[zeta0.ticks+1] - tau.g[zeta0.ticks])
    
    reach <- nknots
    gam0 <- par[reach + 1]; reach <- reach + 1
    sigma <- sigFn(par[reach + 1], a.sig); reach <- reach + 1
    nu <- nuFn(par[reach + 1]);
    
    b0dot <- sigma * base.bundle$q0(zeta0, nu) * zeta0.dot
    beta0.hat <- rep(NA, L)
    beta0.hat[mid:L] <- gam0 + trape(b0dot[mid:L], tau.g[mid:L], L - mid + 1)
    beta0.hat[mid:1] <- gam0 + trape(b0dot[mid:1], tau.g[mid:1], mid)
    
    betas <- beta0.hat
    if(reduce) betas <- betas[reg.ix]
    return(betas)
}


chull.center <- function (x, maxEPts = ncol(x) + 1, plot = FALSE){
    sx <- as.matrix(apply(x, 2, function(s) punif(s, min(s), max(s))))
    dd <- rowSums(scale(sx)^2)
    ix.dd <- order(dd, decreasing = TRUE)
    sx <- sx[ix.dd, , drop = FALSE]
    x.chol <- inchol(sx, maxiter = maxEPts)   # uses default rbfdot (radial basis kernal function, "Gaussian")
    ix.epts <- ix.dd[pivots(x.chol)]          # picks off extreme points
    x.choose <- x[ix.epts, , drop = FALSE]
    xCent <- as.numeric(colMeans(x.choose))   # takes means of extreme points
    attr(xCent, "EPts") <- ix.epts            # returns list of extrememum points
    if (plot) {
        n <- nrow(x)
        p <- ncol(x)
        xjit <- x + matrix(rnorm(n * p), n, p) %*% diag(0.05 * apply(x, 2, sd), p)
        xmean <- colMeans(x)
        x.ept <- x[ix.epts, ]
        M <- choose(p, 2)
        xnames <- dimnames(x)[[2]]
        if (is.null(xnames)) xnames <- paste("x", 1:p, sep = "")
        par(mfrow = c(ceiling(M/ceiling(sqrt(M))), ceiling(sqrt(M))), mar = c(2, 3, 0, 0) + 0.1)
        xmax <- apply(x, 2, max)
        for (i in 1:(p - 1)) {
            for (j in (i + 1):p) {
                plot(x[, i], x[, j], col = "gray", cex = 1, ann = FALSE, ty = "n", axes = FALSE, bty = "l")
                title(xlab = xnames[i], line = 0.3)
                title(ylab = xnames[j], line = 0.3)
                ept <- chull(x[, i], x[, j])
                polygon(x[ept, i], x[ept, j], col = gray(0.9), border = "white")
                points(xjit[, i], xjit[, j], pch = ".", col = gray(0.6))
                points(xmean[i], xmean[j], col = gray(0), pch = 17, cex = 1)
                points(xCent[i], xCent[j], col = gray(0), pch = 1, cex = 2)
                points(x.ept[, i], x.ept[, j], col = gray(.3), pch = 10, cex = 1.5)
            }
        }
    }
    return(xCent)
}


waic <- function(logliks, print = TRUE){
	lppd <- sum(apply(logliks, 1, logmean))
	p.waic.1 <- 2 * lppd - 2 * sum(apply(logliks, 1, mean))
	p.waic.2 <- sum(apply(logliks, 1, var))
	waic.1 <- -2 * lppd + 2 * p.waic.1
	waic.2 <- -2 * lppd + 2 * p.waic.2
	if(print) cat("WAIC.1 =", round(waic.1, 2), ", WAIC.2 =", round(waic.2, 2), "\n")
	invisible(c(WAIC1 = waic.1, WAIC2 = waic.2))
}

# Returns W-tilde_0 (not explictly mentioned in the paper)
# Using t-distribution with three degrees of freedom as prior distribution
# get L-dimensional vector that is needed for likelihood evaluation.

ppFn0 <- function(w.knot, gridmats, L, nknots, ngrid){
	w.grid <- matrix(NA, L, ngrid)
	lpost.grid <- rep(NA, ngrid)
	for(i in 1:ngrid){
		A <- matrix(gridmats[1:(L*nknots),i], nrow = nknots)
		R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),i], nrow = nknots)
		r <- sum.sq(backsolve(R, w.knot, transpose = TRUE))
		w.grid[,i] <- colSums(A * w.knot)
		lpost.grid[i] <- -(0.5*nknots+0.1)*log1p(0.5*r/0.1) - gridmats[nknots*(L+nknots)+1,i] + gridmats[nknots*(L+nknots)+2,i]		
	}
	lpost.sum <- logsum(lpost.grid)
	post.grid <- exp(lpost.grid - lpost.sum)
	w <- c(w.grid %*% post.grid)
	return(list(w = w, lpost.sum = lpost.sum))
}

# Returns W-tilde_j used in the paper, L-dimensional vector needed for likelihood eval.
# Also returns the log posterior sum.

ppFn <- function(w.knot, gridmats, L, nknots, ngrid, a.kap){
	w.grid <- matrix(NA, L, ngrid)
	lpost.grid <- rep(NA, ngrid)
	for(i in 1:ngrid){
		A <- matrix(gridmats[1:(L*nknots),i], nrow = nknots)
		R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),i], nrow = nknots)
		r <- sum.sq(backsolve(R, w.knot, transpose = TRUE))
		w.grid[,i] <- colSums(A * w.knot)
		lpost.grid[i] <- (logsum(-(nknots/2+a.kap[1,])*log1p(0.5*r/ a.kap[2,]) + a.kap[3,] + lgamma(a.kap[1,]+nknots/2)-lgamma(a.kap[1,])-.5*nknots*log(a.kap[2,]))
						  - gridmats[nknots*(L+nknots)+1,i] + gridmats[nknots*(L+nknots)+2,i])		
	}
	lpost.sum <- logsum(lpost.grid)
	post.grid <- exp(lpost.grid - lpost.sum)
	w <- c(w.grid %*% post.grid)
	return(list(w = w, lpost.sum = lpost.sum))
}


lamFn <- function(prox) return(sqrt(-100*log(prox)))
nuFn <- function(z) return(0.5 + 5.5*exp(z/2)) 
nuFn.inv <- function(nu) return(2*log((nu - 0.5)/5.5))
sigFn <- function(z, a.sig) return(exp(z/2)) 
sigFn.inv <- function(s, a.sig) return(2 * log(s))
unitFn <- function(u) return(pmin(1 - 1e-10, pmax(1e-10, u)))

sum.sq <- function(x) return(sum(x^2))
extract <- function(lo, vn) return(lo[[vn]])       # pull out the "vn"th list item
logmean <- function(lx) return(max(lx) + log(mean(exp(lx - max(lx)))))
logsum <- function(lx) return(logmean(lx) + log(length(lx)))
shrinkFn <- function(x) return(1) ##(1/(1 + log(x)))
trape <- function(x, h, len = length(x)) return(c(0, cumsum(.5 * (x[-1] + x[-len]) * (h[-1] - h[-len]))))

getBands <- function(b, col = 2, lwd = 1, plot = TRUE, add = FALSE, x = seq(0,1,len=nrow(b)), remove.edges = TRUE, ...){
	colRGB <- col2rgb(col)/255
	colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.2)
	b.med <- apply(b, 1, quantile, pr = .5)
	b.lo <- apply(b, 1, quantile, pr = .025)
	b.hi <- apply(b, 1, quantile, pr = 1 - .025)
	L <- nrow(b)
	ss <- 1:L; ss.rev <- L:1
	if(remove.edges){
		ss <- 2:(L-1); ss.rev <- (L-1):2
	}
	if(plot){
		if(!add) plot(x[ss], b.med[ss], ty = "n", ylim = range(c(b.lo[ss], b.hi[ss])), ...)  
		polygon(x[c(ss, ss.rev)], c(b.lo[ss], b.hi[ss.rev]), col = colTrans, border = colTrans)
		lines(x[ss], b.med[ss], col = col, lwd = lwd)  
	}
	invisible(cbind(b.lo, b.med, b.hi))
}

# Function to calculate kulback-liebler divergence between two GPS with alternate
# lambdas, detailed in appendix B.2
# Uses formula for kullback liebler divergence of two multiva gauassian dist'ns,
# each of dimension k-knots
# solve(K2,K1) does K2^{-1}%*%K1... and sum/diag around it gives determinant

klGP <- function(lam1, lam2, nknots = 11){
	tau <- seq(0, 1, len = nknots)
	dd <- outer(tau, tau, "-")^2
	K1 <- exp(-lam1^2 * dd); diag(K1) <- 1 + 1e-10; R1 <- chol(K1); log.detR1 <- sum(log(diag(R1)))
	K2 <- exp(-lam2^2 * dd); diag(K2) <- 1 + 1e-10; R2 <- chol(K2); log.detR2 <- sum(log(diag(R2)))
	return(log.detR2-log.detR1 - 0.5 * (nknots - sum(diag(solve(K2, K1)))))
}

# Function to create non-evenly-spaced grid for lambda discretized grid points
# Involves using kulback-liebler divergence.  Makes sure that prior remains sufficiently
# overlapped for neighboring lambda values, preventing (hopefully) poor mixing of
# Markov chain sampler.

proxFn <- function(prox.Max, prox.Min, kl.step = 1){
	prox.grid <- prox.Max
	j <- 1
	while(prox.grid[j] > prox.Min){
		prox1 <- prox.grid[j]
		prox2 <- prox.Min
		kk <- klGP(lamFn(prox1), lamFn(prox2))
		while(kk > kl.step){
			prox2 <- (prox1 + prox2)/2
			kk <- klGP(lamFn(prox1), lamFn(prox2))
		}
		j <- j + 1
		prox.grid <- c(prox.grid, prox2)
	}
	return(prox.grid)
}

transform.grid <- function(w, ticks, dists){
    return((1-dists) * w[ticks] + dists * w[ticks+1])
}
