
##' Compute a confidence interval for the grand mean at a user-specified level.
##' 
##' @param yi a vector containing the primary study measurements
##' @param vi a vector of the same length as yi containing the variances of the of the primary study measurements contained in yi
##' @param c0 a vector containing the mixing parameters for the test statistics; defaults to .2 if length(y)>=6 and .6 otherwise
##' @param level the level of the confidence interval; defaults to .05
##' @param plot indicator whether to plot the contour of the confidence region; defaults to TRUE
##' @param tau2.bounds upper and lower bounds for the range of population variance values for constructing the confidence region; if NULL, value will be calculated from tau2.alpha
##' @param resolution resolution of the population variance values for constructing the confidence region; defaults to 1e2
##' @param Z a matrix of length(yi) rows with each row consisting of standard normal samples to be used in the monte carlo estimation of the null distribution of the test statistic; if NULL, B values
##' will be sampled per row
##' @param B the number of monte carlo replicates per primary study observation to be used; defaults to 300
##' @param tau2.alpha the level of the exact CI with which to bounds on population variance when constructing the confidence region
##' @return a matrix with length(c0) rows and each row containing the lower and upper endpoints of the confidence interval for the given mixing parameter
##' @seealso \code{\link{rma.exact}} for computing entire confidence regions
##' @examples
##' 
##' set.seed(1)
##' K <- 4
##' c0 <- c(.5,1)
##' yi=rnorm(K)*sqrt(vi+tau2)+mu0
##' vi <- (seq(1, 5, length=K))^2
##' Z <- matrix(rnorm(K*5e3),nrow=K)
##' rma.exact.fast(yi,vi,resolution=5e2)
##' @export       
rma.exact.fast <- function(yi,vi,c0=.2*(length(yi)>=6)+.6*(length(yi)<6),level=.05,plot=TRUE,tau2.bounds=NULL,resolution=1e2,Z=NULL,B=3e3,tau2.alpha=.995) {
    
    resolution.tau2 <- resolution
    K <- length(yi)
    if(is.null(Z)) {
        Z <- matrix(rnorm(B*K),nrow=K)
    } else {
        B <- ncol(Z)
    }
    if(is.null(tau2.bounds))  tau2.bounds <- tau2.ci(yi,vi,level=tau2.alpha)
    tau2 <- sort(runif(resolution.tau2)*(diff(tau2.bounds))+tau2.bounds[1])

    inv.vi <- 1/vi
    sum.inv.vi <- sum(inv.vi)
    mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
    tau2.DL.obs <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
    mu.DL.obs <- sum(yi/(vi+tau2.DL.obs))/sum(1/(vi+tau2.DL.obs))

    weights <- outer(vi,tau2,`+`)
    Y <- lapply(1:K, function(k) outer(sqrt(weights[k,]),Z[k,]))
    mu.fixed <- lapply(1:K,function(k)Y[[k]]/vi[k])
    mu.fixed <- Reduce(`+`,mu.fixed)/sum.inv.vi
    tau2.DL <- lapply(1:K,function(k)(Y[[k]]-mu.fixed)^2/vi[k])
    tau2.DL <- (Reduce(`+`,tau2.DL) - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
    tau2.DL[tau2.DL<0] <- 0
    weights.DL <- lapply(1:K,function(k) tau2.DL+vi[k])
    weights.DL.inv <- lapply(weights.DL,function(x)1/x)
    mu.DL <- lapply(1:K,function(k) Y[[k]]*weights.DL.inv[[k]])
    weights.DL.inv.summed <- Reduce(`+`,weights.DL.inv)
    mu.DL <- Reduce(`+`,mu.DL) / weights.DL.inv.summed

    T1 <- mu.DL^2 * weights.DL.inv.summed
    T2 <- lapply(1:K,function(k)( Y[[k]] - mu.DL)^2 * weights.DL.inv[[k]])
    T2 <- Reduce(`+`,T2) + Reduce(`+`,lapply(weights.DL,log))
    T2 <- -T2 + outer(colSums(log(weights)),colSums(Z^2),`+`)
    T2 <- T2/2

    ss <- sum(1/(vi+tau2.DL.obs))
    
    CI.by.c0 <- lapply(c0,function(c0.i) {
        
        T0 <- T1 + T2*c0.i
        q <- apply(T0,1,function(cdf)quantile(cdf,probs=1-level))

        u <- (c0.i/2)*colSums(1/weights) + ss
        v <- -2*mu.DL.obs*ss - c0.i*colSums(yi/weights)
        w <- (c0.i/2)*(colSums(yi^2/weights) + colSums(log(weights/(tau2.DL.obs+vi))) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi))) + mu.DL.obs^2*ss - q
        ## print(    cbind(lower=(-v-sqrt(v^2-4*u*w))/(2*u),upper=(-v+sqrt(v^2-4*u*w))/(2*u)))
        det <- v^2-4*u*w
        det[det<0] <- NA
        cbind(lower=(-v-sqrt(det))/(2*u),upper=(-v+sqrt(det))/(2*u))
    })
    ## print(CI.by.c0)

    CIs <- lapply(CI.by.c0,function(CI.by.tau2)c(min(CI.by.tau2[,'lower'],na.rm=TRUE),max(CI.by.tau2[,'upper'],na.rm=TRUE)))
    CIs <- matrix(unlist(CIs),nrow=2)
    rownames(CIs) <- c('lower','upper')
    colnames(CIs) <- c0
    CIs <- t(CIs)
    
    CI.by.c0 <- abind::abind(CI.by.c0,along=3)   
    dimnames(CI.by.c0)[[1]] <- round(tau2,2)
    dimnames(CI.by.c0)[[3]] <- c0

    colors <- cm.colors(length(c0))
    if(plot){
        plot(x=NULL,y=NULL,xlim=range(CI.by.c0,na.rm=TRUE),ylim=range(tau2),xlab=expression(mu),ylab=expression(tau^2))
        for(i in 1:length(c0)) {
            lines(CI.by.c0[,'lower',i],tau2,col=colors[i])
            lines(CI.by.c0[,'upper',i],tau2,col=colors[i])
        }
        if(length(c0)>1)
            legend('topright',col=colors,lty=1,legend=c0,title=expression(c[0]))

    }

    ## return(structure(CIs,CI.by.c0=CI.by.c0))
    return(CIs)
}


rma.exact.c0 <- function(yi,vi,c0=1,mu.bounds=NULL,tau2.bounds=NULL,resolution=1e2,Z=NULL,B=3e3,resolution.mu=resolution,resolution.tau2=resolution,mu.alpha=.995,tau2.alpha=.995) {

    K <- length(yi)
    if(is.null(Z)) {
        Z <- matrix(rnorm(B*K),nrow=K)
    } else {
        B <- ncol(Z)
    }
    if(is.null(mu.bounds)) {
        fit=metafor::rma.uni(yi=yi, vi=vi, method="DL")
        mu.bounds <- fit$b+c(-1,1)*fit$se*qt(mu.alpha,df=K-1)
    }
    if(is.null(tau2.bounds))  tau2.bounds <- tau2.ci(yi,vi,level=tau2.alpha)
    mu <- sort(runif(resolution.mu)*(diff(mu.bounds))+mu.bounds[1])
    tau2 <- sort(runif(resolution.tau2)*(diff(tau2.bounds))+tau2.bounds[1])

    ## print(tau2.bounds)
    inv.vi <- 1/vi
    sum.inv.vi <- sum(inv.vi)
    mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
    tau2.DL.obs <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
    mu.DL.obs <- sum(yi/(vi+tau2.DL.obs))/sum(1/(vi+tau2.DL.obs))

    ## T1.obs <- (mu.DL.obs-grid.points$mu)^2 / sum(1/(tau2.DL.obs+vi))
    T1.obs <- (mu.DL.obs-mu)^2 * sum(1/(tau2.DL.obs+vi))
    T2.obs <- lapply(1:K, function(i) outer((yi[i]-mu)^2,tau2+vi[i],`/`))
    T2.obs <- Reduce(`+`,T2.obs)
    T2.obs <- t(T2.obs) + colSums(log(outer(vi,tau2,`+`)))
    T2.obs <- t(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))
    ## T2.obs <- t(T2.obs)

    weights <- outer(vi,tau2,`+`)
    Y <- lapply(1:K, function(k) outer(sqrt(weights[k,]),Z[k,]))
    mu.fixed <- lapply(1:K,function(k)Y[[k]]/vi[k])
    mu.fixed <- Reduce(`+`,mu.fixed)/sum.inv.vi
    tau2.DL <- lapply(1:K,function(k)(Y[[k]]-mu.fixed)^2/vi[k])
    tau2.DL <- (Reduce(`+`,tau2.DL) - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
    tau2.DL[tau2.DL<0] <- 0
    ## weights.DL <- 1/outer(tau2.DL,vi,`+`)
    ## mu.DL <- lapply(1:K,function(k)Y[[k]]*weights.DL[,,k])
    ## mu.DL <- Reduce(`+`,mu.DL) / apply(weights.DL,1:2,sum)
    weights.DL <- lapply(1:K,function(k) tau2.DL+vi[k])
    weights.DL.inv <- lapply(weights.DL,function(x)1/x)
    mu.DL <- lapply(1:K,function(k) Y[[k]]*weights.DL.inv[[k]])
    weights.DL.inv.summed <- Reduce(`+`,weights.DL.inv)
    mu.DL <- Reduce(`+`,mu.DL) / weights.DL.inv.summed

    T1 <- mu.DL^2 * weights.DL.inv.summed
    T2 <- lapply(1:K,function(k)( Y[[k]] - mu.DL)^2 * weights.DL.inv[[k]])
    T2 <- Reduce(`+`,T2) + Reduce(`+`,lapply(weights.DL,log))
    T2 <- -T2 + outer(colSums(log(weights)),colSums(Z^2),`+`)

    p.vals <- lapply(c0,function(c0.i) {
        T0.obs <- t(T1.obs + T2.obs*c0.i)

        T0 <- T1 + T2*c0.i

        p.val <- cbind(T0.obs,T0)
        p.val <- t(apply(p.val,1,rank))
        p.val <- p.val[,1:length(mu)]
        p.val <- p.val - t(apply(T0.obs,1,rank))
        p.val <- 1 - p.val / B
        ## p.val <- structure(t(p.val),dimnames=c(
        dimnames(p.val) <- list(tau2=round(tau2,2),mu=round(mu,2))
        t(p.val)
    })


    p.vals <- abind::abind(p.vals,along=3)
    dimnames(p.vals)[[3]] <- c0

    return(structure(p.vals,class=c('RMA.Exact','matrix'),mu=mu,tau2=tau2))
}

##' Compute a confidence region for grand mean.
##' 
##' @param yi a vector containing the primary study measurements
##' @param vi a vector of the same length as yi containing the variances of the of the primary study measurements contained in yi
##' @param c0 a vector containing the mixing parameters for the test statistics; defaults to .2 if length(y)>=6 and .6 otherwise
##' @param level the level of the confidence interval; defaults to .05
##' @param mu.bounds upper and lower bounds for the range of population effect values for constructing the confidence region; if NULL, value will be calculated from mu.alpha
##' @param tau2.bounds upper and lower bounds for the range of population variance values for constructing the confidence region; if NULL, value will be calculated from tau2.alpha
##' @param resolution resolution of the population mean and variance values within the bounding box; defaults to 1e2 for each of the two dimensions
##' @param resolution.mu resolution of the population mean values within the bounding box; defaults to resolution
##' @param resolution resolution of the population variance values within the bounding box; defaults to resolution
##' @param Z a matrix of length(yi) rows with each row consisting of standard normal samples to be used in the monte carlo estimation of the null distribution of the test statistic; if NULL, B values
##' will be sampled per row
##' @param B the number of monte carlo replicates per primary study observation to be used; defaults to 300
##' @param mu.alpha the level of the exact CI for constructing the bounds on the population mean dimension of the bounding box
##' @param tau2.alpha  the level of the exact CI for constructing the bounds on the population variance dimension of the bounding box
##' @param test.stat (currently for internal use)
##' @param ... (currently for internal use)
##' @return an object of class  RMA.Exact
##' @seealso \code{\link{rma.exact.fast}} for computing confidence intervals at specified levels, \code{\link{plot.RMA.Exact}}, \code{\link{confint.RMA.Exact}}
##' @examples
##' set.seed(1)
##' 
##' K <- 5
##' c0 <- 1
##' mu0 <- 0
##' tau2 <- 12.5
##' vi <- (seq(1, 5, length=K))^2
##' yi=rnorm(K)*sqrt(vi+tau2)+mu0
##' rma0 <- rma.exact(yi=yi,vi=vi)
##' plot(rma0)
##' confint(rma0)
##' 
##' ## multiple c0 values
##' c0 <- c(0,.25,1)
##' tau2 <- 12.5
##' vi <- (seq(1, 5, length=K))^2
##' yi=rnorm(K)*sqrt(vi+tau2)+mu0
##' rma0 <- rma.exact(yi=yi,vi=vi,c0=c0)
##' plot(rma0)
##' confint(rma0)
##' 
##' ## setting tau2.bounds and other parameters to non-default values
##' Z <- matrix(rnorm(K*5e3),nrow=K)
##' B <- ncol(Z)
##' resolution <- 3e2
##' rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=resolution,c0=c0,tau2.bounds=c(1,120),resolution.tau2=1e3,resolution.mu=1e2)
##' plot(rma0)
##' 
##' c0 <- 1:4
##' rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=resolution,c0=c0,tau2.bounds=c(1,450),resolution.tau2=1e3,resolution.mu=1e2)
##' plot(rma0)
##' confint(rma0,levels=c(.05))
##' @export
rma.exact <- function(yi,vi,c0=1,mu.bounds=NULL,tau2.bounds=NULL,resolution= 1e2,Z=NULL,B=3e3,resolution.mu=resolution,resolution.tau2=resolution,mu.alpha=.995,tau2.alpha=.995,test.stat=NULL,...) {
    if(is.null(test.stat)) return(rma.exact.c0(yi=yi,vi=vi,c0=c0,mu.bounds=mu.bounds,tau2.bounds=tau2.bounds,resolution=resolution,Z=Z,B=B,resolution.mu=resolution.mu,resolution.tau2=resolution.tau2,mu.alpha=mu.alpha,tau2.alpha=tau2.alpha))

    K <- length(yi)
    if(is.null(Z)) {
        Z <- matrix(rnorm(B*K),nrow=K)
    } else {
        B <- ncol(Z)
    }
    if(is.null(mu.bounds)) {
        fit=metafor::rma.uni(yi=yi, vi=vi, method="DL")
        mu.bounds <- fit$b+c(-1,1)*fit$se*qt(mu.alpha,df=K-1)
    }
    if(is.null(tau2.bounds))  tau2.bounds <- tau2.ci(yi,vi,level=tau2.alpha)
    mus <- sort(runif(resolution.mu)*(diff(mu.bounds))+mu.bounds[1])
    tau2s <- sort(runif(resolution.tau2)*(diff(tau2.bounds))+tau2.bounds[1])

    ## print(tau2.bounds)
    T.distr <- sapply(tau2s,function(tau2) {
        Y <- Z*sqrt(tau2+vi)
        apply(Y,2,function(yi)test.stat(yi,vi,mu=0,tau2=tau2,...))
    })
    ## T.distr <- t(T.distr)

    yi.translates <- outer(yi,mus,`-`)
    T.obs <- sapply(tau2s,function(tau2) {
        apply(yi.translates,2,function(yi) {
            test.stat(yi=yi,vi=vi,tau2=tau2,...)
            })
    })

    p.val <- rbind(T.obs,T.distr)
    p.val <- t(apply(p.val,2,rank))
    p.val <- p.val[,1:length(mus)]
    p.val <- p.val - t(apply(T.obs,2,rank))
    p.val <- 1 - p.val / B
    dimnames(p.val) <- list(tau2=round(tau2s,2),mu=round(mus,2))
    p.val <- t(p.val)


    return(structure(p.val,class=c('RMA.Exact','matrix'),mu=mus,tau2=tau2s,T.distr=T.distr,T.obs=T.obs))
}


##' Combine confidence regions from two RMA.Exact objects
##'
##' @param rma1 the first RMA.Exact object
##' @param rma2 the second RMA.Exact object
##' @return an object of class RMA.Exact
##' @export
`+.RMA.Exact` <- function(rma1,rma2) {
    structure(rbind(rma1,rma2))
}


##' Plot a confidence region given by an RMA.Exact object
##' 
##' @param rma0 an object of class RMA.Exact
##' @param levels the significance levels to plot
##' @param mfrow the dimensions of the array of plotting windows for use when rma0 contains regions for multiple weight parameters; defaults to a single row with as many as columns as regions
##' contained in rma0
##' @param ... additional parameters passed down to graphics::contour (passed from there down to plot.window)
##' @return undefined
##' @examples
##' see ?RMA.Exact
##' @export
plot.RMA.Exact <- function(rma0,levels=c(.01,.02,.05,.1,.15,.2),mfrow=c(1,dim(rma0)[3]),...) {
    if(length(dim(rma0))==2) {
        dim(rma0) <- c(dim(rma0),1)
        op <- par(no.readonly = TRUE)
    } else {
        op <- par(mfrow=mfrow)
    }
    mu <- attr(rma0,'mu')
    tau2 <- attr(rma0,'tau2')
    apply(rma0,3,function(rma0.i) {
        image(mu,tau2,structure(rma0.i,class='matrix'),col=cm.colors(12),xlab=expression(mu),ylab=expression(tau^2))
        contour(mu,tau2,rma0.i,add=TRUE,levels=levels,...)
    })
    par(op)
}



##' Compute a confidence interval for the population mean from an RMA.Exact object
##'
##' A warning will be issued if the reported endpoints of the confidence interval are near the bounding box.
##' 
##' @param rma0 an object of class RMA.Exact
##' @param levels the significance levels at which the compute the confidence intervals
##' @param ... (currently for internal use)
##' @return a matrix with a row corresponding to each weight parameter c0 stored in rma0, and columns containing the upper and lower interval endpoints, the population variance values at whcih those
##' endpoints were obtained, etc.
##' @examples
##' see ?RMA.Exact
##' @export
confint.RMA.Exact <- function(rma0,levels=.05,...) {
    mu <- attr(rma0,'mu')
    tau2 <- attr(rma0,'tau2')
    CI.by.c0 <- lapply(1:dim(rma0)[3],function(i) {
        contour.xy <- do.call(rbind, lapply(contourLines(mu,tau2,rma0[,,i],levels=levels),as.data.frame))
        colnames(contour.xy) <- c('level','mu','tau2')
        CIs <- sapply(split(contour.xy,contour.xy$level),function(df)c(mu.lower=min(df$mu),mu.upper=max(df$mu),tau2.lower=df$tau2[which.min(df$mu)],tau2.upper=df$tau2[which.max(df$mu)]))
        CIs <- t(CIs)
        CIs <- cbind(CIs,length=CIs[,'mu.upper']-CIs[,'mu.lower'],level=as.numeric(rownames(CIs)),c0=as.numeric(dimnames(rma0)[[3]][i]))
        as.data.frame(CIs)
    })

    CIs <- do.call(rbind,CI.by.c0)
    rownames(CIs) <- NULL

    if (sum(CIs[,'tau2.upper']>quantile(tau2,.95)))
        warning('CI near edges of sampled region; increase tau^2 bounds')
    if (sum(CIs[,'mu.upper']>quantile(mu,.99)) | sum(CIs[,'mu.lower']<quantile(mu,.01)))
        warning('CI near edges of sampled region; increase mu bounds')

    return(CIs)

}


pivot <- function(tau2,yi,vi) { ## notation as in ms
    K <- length(yi)
    weights <- (1/vi)/sum(1/vi)
    W <- t(matrix(weights,nrow=K,ncol=K))
    W <- W[-K,,drop=FALSE]
    diag(W) <- diag(W)-1
    Sigma <- diag(vi+tau2)
    WY <- W%*%yi
    t(WY) %*% solve(W%*%Sigma%*%t(W)) %*% WY
}

tau2.ci <- function(yi,vi,level=.995) {
    K <- length(yi)
    q.lower <- qchisq(1-level, df=K-1)
    if (q.lower > pivot(0,yi,vi)) {
        warning('Empty CI')
        stopifnot(q.lower < pivot(0,yi,vi))
        ## return(c(CI.lower=0,CI.upper=0))##TODO
    }
    tau2.bound <- 9*var(yi)
    while(pivot(tau2.bound,yi,vi)>=q.lower) tau2.bound <- 2*tau2.bound
       CI.upper <- uniroot(function(tau2)pivot(tau2,yi,vi) - q.lower, interval=c(0,tau2.bound))$root
    return(c(CI.lower=0,CI.upper=CI.upper))
}







## Lu's version


rma.null=function(tau0, thetahat, varhat, error.mat, c0=1)
   {K=length(thetahat)
    B=length(error.mat[1,])
    L=length(c0)

    theta.star=error.mat*sqrt(varhat+tau0)
    weight.ini=1/varhat
    mu.ini=apply(theta.star*weight.ini,2,sum)/sum(weight.ini)
    tau.sim=pmax((apply(t(t(theta.star)-mu.ini)^2*weight.ini,2,sum)-(K-1))/(sum(weight.ini)-sum(1/varhat^2)/sum(weight.ini)), 0)
    tau.sim.mat=t(matrix(tau.sim, B, K))
    weight.dl=1/(tau.sim.mat+varhat)
    mu.sim=apply(theta.star*weight.dl,2,sum)/apply(weight.dl,2,sum)

    test.null=matrix(0, B, L)
        
    for(ll in 1:L)
       {if(min(c0)<Inf) test2.null=apply(theta.star^2/(varhat+tau0)+log(varhat+tau0), 2, sum)-apply( (t(t(theta.star)-mu.sim))^2*weight.dl-log(weight.dl), 2, sum)
        test1.null=mu.sim^2*apply(weight.dl, 2, sum)    
        if(c0[ll]<Inf)
          {test.null[,ll]=test2.null+test1.null*c0[ll]}
        if(c0[ll]==Inf)
          {test.null[,ll]=test1.null}    
        }

    return(test.null)
   }


######################################################

tau.ci=function(thetahat, varhat, level=0.999)
               {test.tau=function(tau2, thetahat, varhat, cov)
                                 {K=length(varhat)
                                  weight=1/(varhat+tau2)/sum(1/(varhat+tau2))
                                  weight.mat=matrix(0, (K-1), K)
                                  for(i in 1:(K-1))
                                     {weight.mat[i,]=weight
                                      weight.mat[i,i]=weight.mat[i,i]-1
                                      }

                                  varmat=diag(varhat+tau2, K, K)
                                  sigma=weight.mat%*%varmat%*%t(weight.mat)
                                  test.stat=t(weight.mat%*%thetahat)%*%solve(sigma)%*%(weight.mat%*%thetahat)
                                  c0=qchisq(cov, K-1)
                                  return(as.numeric(test.stat)-c0)
                                  }  
                                  
              if(test.tau(0, thetahat, varhat, 1-level)<=0)
                ci.upptau=0                 
                
              if(test.tau(0, thetahat, varhat, 1-level)>0)  
                {lower=0
                 upper=(max(thetahat+3*sqrt(varhat))-min(thetahat-3*sqrt(varhat)))^2/2
                 while(test.tau(upper, thetahat, varhat, 1-level)>0)                           {lower=upper
                       upper=1.5*upper}

                 ci.upptau=uniroot(test.tau, c(lower, upper), thetaha=thetahat, varhat=varhat, cov=1-level)$root
                 }
               
                return(c(0, ci.upptau))
               }


#############################################################
rma.exact.lu=function(thetahat, varhat,  B=1000, num.tau=300, c0=1, level=0.998){

c0=sort(unique(c0))
L=length(c0)
K=length(varhat)
set.seed(100)
error.mat=matrix(rnorm(B*K), K, B)


fit.obs=rma.uni(yi=thetahat, vi=varhat, method="DL") 
muhat=as.numeric(fit.obs$b)
tauhat=fit.obs$tau2

sigmainv.hat=sum(1/(tauhat+varhat))
llk.obs=sum(-(thetahat-muhat)^2/(varhat+tauhat)-log(varhat+tauhat))

citau=tau.ci(thetahat, varhat, level=level)
tau.sup=citau[2]
tau.inf=citau[1]


step.tau=qnorm((1:num.tau)/(2*num.tau+1)+0.5)
step.tau=(tau.sup-tau.inf)/sum(step.tau)*step.tau
step.tau=c(step.tau, step.tau[num.tau])

##################################################################################

tau.up=tauhat+cumsum(step.tau)
tau.down=tauhat-cumsum(step.tau)

tau.grd=c(rev(tau.down), tauhat, tau.up)
tau.grd=tau.grd[tau.grd>=0 & tau.grd<=tau.sup]

ntau=length(tau.grd)


test.cut=matrix(NA, ntau, L)
for(b in 1:ntau)
   {tau.new=tau.grd[b]
    test.null=rma.null(tau0=tau.new, thetahat=thetahat, varhat=varhat, error.mat=error.mat, c0=c0)
    for(ll in 1:L)
       test.cut[b, ll]=quantile(test.null[,ll], 0.95)
    }

aa=bb=cc=matrix(0, ntau, L)

thetahat.mat=t(matrix(rep(thetahat, ntau), K, ntau))
varhat.mat=t(matrix(rep(varhat, ntau), K, ntau))

cc=t(matrix(rep(c0*muhat^2*sigmainv.hat, ntau), L, ntau))+llk.obs+apply(thetahat.mat^2/(varhat.mat+tau.grd)+log(varhat.mat+tau.grd),1,sum)-test.cut
bb=-2*(t(matrix(rep(c0*muhat*sigmainv.hat, ntau), L, ntau))+apply(thetahat.mat/(varhat.mat+tau.grd),1,sum))
aa=t(matrix(rep(c0*sigmainv.hat, ntau), L, ntau))+apply(1/(varhat.mat+tau.grd),1,sum)

if(max(c0)==Inf)
  {cc[,L]=muhat^2*sigmainv.hat-test.cut[,L]
   bb[,L]=-2*muhat*sigmainv.hat
   aa[,L]=sigmainv.hat
   }

delta=bb^2-4*aa*cc
delta[delta<0]=0
delta=sqrt(delta)
mu.min=(-bb-delta)/2/aa
mu.max=(-bb+delta)/2/aa

   
ci.low=apply(mu.min,2,min)
ci.upp=apply(mu.max,2,max)

return(list(ci=cbind(c0, ci.low, ci.upp), tau.grd=tau.grd, mu.min=mu.min, mu.max=mu.max))
}

