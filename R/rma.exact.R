## require(metafor)
## require(abind)


##' Compute a confidence interval for the grand mean at a user-specified level.
##' 
##' @param yi a vector containing the primary study measurements
##' @param vi a vector of the same length as yi containing the variances of the of the primary study measurements contained in yi
##' @param c0 a vector containing the mixing parameters for the test statistics; defaults to 1 
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
rma.exact.fast <- function(yi,vi,c0=1,level=.05,plot=TRUE,tau2.bounds=NULL,resolution=1e2,Z=NULL,B=3e3,tau2.alpha=.995) {
    
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
        cbind(lower=(-v-sqrt(v^2-4*u*w))/(2*u),upper=(-v+sqrt(v^2-4*u*w))/(2*u))
    })
    ## print(CI.by.c0)

    CIs <- lapply(CI.by.c0,function(CI.by.tau2)c(min(CI.by.tau2[,'lower']),max(CI.by.tau2[,'upper'])))
    CIs <- matrix(unlist(CIs),nrow=2)
    rownames(CIs) <- c('lower','upper')
    colnames(CIs) <- c0
    CIs <- t(CIs)
    
    CI.by.c0 <- abind::abind(CI.by.c0,along=3)   
    dimnames(CI.by.c0)[[1]] <- round(tau2,2)
    dimnames(CI.by.c0)[[3]] <- c0

    colors <- cm.colors(length(c0))
    if(plot){
        plot(x=NULL,y=NULL,xlim=range(CI.by.c0),ylim=range(tau2),xlab=expression(mu),ylab=expression(tau^2))
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


rma.exact.c0 <- function(yi,vi,c0=1,mu.bounds=NULL,tau2.bounds=NULL,resolution=1e2,Z=NULL,B=ncol(Z),resolution.mu=resolution,resolution.tau2=resolution,mu.alpha=.995,tau2.alpha=.995) {

    K <- length(yi)
    if(is.null(Z)) Z <- matrix(rnorm(B*K),nrow=K)
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
##' @param c0 a vector containing the mixing parameters for the test statistics; defaults to 1 
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
rma.exact <- function(yi,vi,c0=1,mu.bounds=NULL,tau2.bounds=NULL,resolution=1e2,Z=NULL,B=3e3,resolution.mu=resolution,resolution.tau2=resolution,mu.alpha=.995,tau2.alpha=.995,test.stat=NULL,...) {
    if(is.null(test.stat)) return(rma.exact.c0(yi=yi,vi=vi,c0=c0,mu.bounds=mu.bounds,tau2.bounds=tau2.bounds,resolution=resolution,Z=Z,B=B,resolution.mu=resolution.mu,resolution.tau2=resolution.tau2,mu.alpha=mu.alpha,tau2.alpha=tau2.alpha))

    K <- length(yi)
    if(is.null(Z)) Z <- matrix(rnorm(B*K),nrow=K)
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


## plot.RMA.Exact <- function(rma0,levels=c(.01,.02,.05,.1,.15,.2),...) {
##     mu <- as.numeric(rownames(rma0))
##     tau2 <- as.numeric(colnames(rma0))
##     image(mu,tau2,structure(rma0,'class'='matrix'),col=cm.colors(12),xlab=expression(mu),ylab=expression(tau^2))
##     contour(mu,tau2,rma0,add=T,levels=levels,...)
## }

## ## support for non-square grid points
## plot.RMA.Exact <- function(rma0,levels=c(.01,.02,.05,.1,.15,.2),...) {
##     stopifnot(length(unique(rma0$c0))==1)
##     df.cnt <- getContourLines(as.data.frame(select(rename(rma0,x=mu,y=tau2,z=p.val),x,y,z)), levels=levels)#,...)
##     df.cnt$z <- as.factor(df.cnt$z)
##     df.cnt <- rename(df.cnt,mu=x,tau2=y,p.val=z)
##     plt <- ggplot(df.cnt,aes(x=mu,y=tau2,group=Group,color=p.val))+geom_path()+theme_bw()+scale_color_manual(values=brewer.pal(length(levels),'PuBu'))
##     plt+labs(x=expression(mu),y=expression(tau^2))
## }

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








## confint.RMA.Exact <- function(rma0,alpha=.05,...) {

##     CI.list <- lapply(split(rma0,rma0$c0),function(rma00) {
##         df.cnt <- getContourLines(as.data.frame(select(rename(rma00,x=mu,y=tau2,z=p.val),x,y,z)),levels=alpha)#,...)
##         ## CI <- aggregate(x ~ z, data=df.cnt, function(x0)c(mu.lower=min(x0),mu.upper=max(x0),length=max(x0)-min(x0)))
##         CI <- group_by(df.cnt,z) %>% summarize(mu.lower=min(x),mu.upper=max(x),tau2.lower=y[which.min(x)],tau2.upper=y[which.max(x)],length=max(x)-min(x)) %>% rename(alpha=z) %>% as.data.frame
##         CI$c0 <- rma00$c0[1]
##         return(CI)


##         ## df.cnt$y <- round(df.cnt$y,2) ## TODO better smoothing
##         ## CI <- aggregate(x ~ y + z, data=df.cnt,function(x)c(min.mu=min(x),max.mu=max(x),length=max(x)-min(x),n=length(x)))
##         ## CI <- data.frame(tau2=CI$y,p.val=CI$z,as.data.frame(CI$x),c0=unique(rma0$c0))

##         ## ## plot(c(CI$min.mu,CI$max.mu),c(CI$tau2,CI$tau2))
##         ## ## plot(CI$length,CI$tau2)
##         ## group_by(CI,p.val) %>% slice(which.max(length)) %>% as.matrix

##     })

##     do.call(rbind,CI.list)
## }



## ## using tensorA--more readable code and some speedup
## rma.exact.tensor <- function(mu,tau2,yi,vi,Z,c0=1) {
##     K <- length(yi)
##     B <- ncol(Z)

##     ## vi <- as.tensor(vi,dims=c(K=K))
##     ## yi <- as.tensor(yi,dims=c(K=K))
##     inv.vi <- as.tensor(1/vi,dims=c(K=K))
##     sum.inv.vi <- sum(inv.vi)
##     mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
##     tau2.DL.obs <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
##     mu.DL.obs <- sum(yi/(vi+tau2.DL.obs))/sum(1/(vi+tau2.DL.obs))

##     grid.points <- expand.grid(mu,tau2)
##     N <- length(mu)*length(tau2)
##     colnames(grid.points) <- c('mu','tau2')
##     ## T1.obs <- (mu.DL.obs-grid.points$mu)^2 / sum(1/(tau2.DL.obs+vi))
##     ## TODO move next few lines to end, where these vars are used
##     T1.obs <- (mu.DL.obs-grid.points$mu)^2 * sum(1/(tau2.DL.obs+vi))
##     T2.obs <- outer(yi,grid.points$mu,`-`)^2 / outer(vi,grid.points$tau2,`+`) + log(outer(vi,grid.points$tau2,`+`))
##     ## T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) / 2
##     T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) #/ 2


##     ZZ <- to.tensor(as.numeric(Z),c(K=K,B=B))
##     w.unnormalized <- as.tensor(outer(as.numeric(vi),grid.points$tau2,`+`),dims=c(K=K,N=N))
##     w <- w.unnormalized/margin.tensor(w.unnormalized,1)
##     Y <- ZZ*sqrt(w.unnormalized) + as.tensor(grid.points$mu,dims=c(N=N))
##     mu.fixed <- margin.tensor(Y*(inv.vi/sum.inv.vi),1)
##     tau2.DL <- (margin.tensor((Y-mu.fixed)^2*inv.vi,1) - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
##     tau2.DL[tau2.DL<0] <- 0
##     ## tau2.DL <- (apply((Y-spread(mu.fixed,K))^2*inv.vi,2:3,sum) - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
##     ## tau2.DL <- apply(tau2.DL,1:2,function(x)max(x,0))
##     w.hat.inv.unnormalized <- 1/(tau2.DL+as.tensor(vi,dims=c(K=K)))
##     w.hat.inv.sums <- margin.tensor(w.hat.inv.unnormalized,3)
##     w.hat.inv <- w.hat.inv.unnormalized/w.hat.inv.sums
##     mu.DL <- margin.tensor(Y*w.hat.inv,1)


##     ## T1 <- (mu.DL-grid.points$mu)^2 * w.hat.inv.sums
##     T1 <- (mu.DL-as.tensor(grid.points$mu,dims=c(N=N)))^2 * w.hat.inv.sums
##     T2 <- -(Y - mu.DL)^2 * w.hat.inv.unnormalized + log(w.hat.inv.unnormalized) + ZZ^2 + log(w.unnormalized)
##     ## T2 <- -(Y - spread(mu.DL,K))^2 * w.hat.inv.unnormalized + log(w.hat.inv.unnormalized) + ZZ^2 + replicate(B,log(w.unnormalized))
##     ## T2 <- apply(T2,2:3,sum) / 2
##     ## T2 <- apply(T2,2:3,sum)
##     T2 <- margin.tensor(T2,1)

##     T0.obs <- rowSums(T1.obs<=t(T1)) + rowSums(T2.obs<=t(T2))/c0 ##[[mean or sum affects c0 scale, match to ms]]
##     T0 <- B+1-t(apply(T1,2,rank)) + (B+1-t(apply(T2,2,rank)))/c0

##     ## for(n in ls())assign(paste0(n,'.old'),get(n))
##     structure(matrix(rowMeans(T0.obs >=T0),nrow=length(mu)), dimnames=list(round(mu,8),round(tau2,8)), mu=mu, tau2=tau2, class=c('RMA.Exact','matrix'))

## }

## ## optimized--remaining bottlenecks are in calls to rank() and to log()
## rma.exact.fast <- function(mu,tau2,yi,vi,Z,c0=1) {
##     K <- length(yi)
##     B <- ncol(Z)

##     ## vi <- as.tensor(vi,dims=c(K=K))
##     ## yi <- as.tensor(yi,dims=c(K=K))
##     inv.vi <- as.tensor(1/vi,dims=c(K=K))
##     sum.inv.vi <- sum(inv.vi)
##     mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
##     tau2.DL.obs <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
##     mu.DL.obs <- sum(yi/(vi+tau2.DL.obs))/sum(1/(vi+tau2.DL.obs))

##     grid.points <- expand.grid(mu,tau2)
##     N <- length(mu)*length(tau2)
##     colnames(grid.points) <- c('mu','tau2')
##     ## T1.obs <- (mu.DL.obs-grid.points$mu)^2 / sum(1/(tau2.DL.obs+vi))
##     T1.obs <- (mu.DL.obs-grid.points$mu)^2 * sum(1/(tau2.DL.obs+vi))
##     T2.obs <- outer(yi,grid.points$mu,`-`)^2 / outer(vi,grid.points$tau2,`+`) + log(outer(vi,grid.points$tau2,`+`))
##     ## T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) / 2
##     T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) #/ 2

##     inv.vi <- as.numeric(inv.vi)
##     ZZ <- to.tensor(as.numeric(Z),c(K=K,B=B))
##     w.unnormalized <- as.tensor(outer(as.numeric(vi),grid.points$tau2,`+`),dims=c(K=K,N=N))
##     w <- w.unnormalized/margin.tensor(w.unnormalized,1)

##     for (i in 1:K)
##         assign(paste0('Y',i), (outer(sqrt(w.unnormalized[i,]),ZZ[i,])+grid.points$mu))
##     mu.fixed <- Reduce(`+`,lapply(1:K,function(k)get(paste0('Y',k))*(inv.vi[k]/sum.inv.vi)))
##     tau2.DL <- Reduce(`+`,lapply(1:K,function(k)
##                                  (get(paste0('Y',k))-mu.fixed)^2*inv.vi[k]
##                                  ))
##     tau2.DL <- (tau2.DL - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
##     tau2.DL[tau2.DL<0] <- 0
##     for(k in 1:K)assign(paste0('w.hat.inv.unnormalized',k),1/(tau2.DL+vi[k]))
##     w.hat.inv.sums <- Reduce(`+`,lapply(1:K,function(k)get(paste0('w.hat.inv.unnormalized',k))))
##     for(k in 1:K)assign(paste0('w.hat.inv',k),get(paste0('w.hat.inv.unnormalized',k))/w.hat.inv.sums)
##     mu.DL <- Reduce(`+`,lapply(1:K,function(k)get(paste0('Y',k))*get(paste0('w.hat.inv',k))))

##     T1 <- (mu.DL-grid.points$mu)^2 * w.hat.inv.sums
##     class(ZZ) <- 'matrix'
##     class(w.unnormalized) <- 'matrix'
##     for(k in 1:K)
##         assign(paste0('T2',k),
##                t(-(get(paste0('Y',k)) - mu.DL)^2 * get(paste0('w.hat.inv.unnormalized',k)) + log(get(paste0('w.hat.inv.unnormalized',k))) + log(w.unnormalized[k,])) + ZZ[k,]^2
##                )##bottleneck on log call
##     T2 <- t(Reduce(`+`,lapply(1:K,function(k)get(paste0('T2',k)))))
##     T0.obs <- rowSums(T1.obs<=T1) + rowSums(T2.obs<=T2)/c0
##     T0 <- B+1-t(apply(T1,1,rank)) + (B+1-t(apply(T2,1,rank)))/c0 ##bottleneck on rank

##     structure(matrix(rowMeans(T0.obs >=T0),nrow=length(mu)), dimnames=list(round(mu,8),round(tau2,8)), mu=mu, tau2=tau2, class=c('RMA.Exact','matrix'))

## }

## ## support for multiple c0 values
## rma.exact.fast <- function(mu=NULL,tau2=NULL,yi=NULL,vi,Z,c0=1,resolution=1e2) {
##     K <- length(yi)
##     B <- ncol(Z)

##     if(is.null(mu)) {
##         fit=rma.uni(yi=yi, vi=vi, method="DL")
##         mu.bounds <- fit$b+c(-1,1)*fit$se*qnorm(.999)
##         mu <- sort(runif(resolution)*(diff(mu.bounds))+mu.bounds[1])
##     }
##     if(is.null(tau2)) {
##         tau2.bounds <- tau.ci(thetahat=yi,varhat=vi)
##         tau2 <- sort(runif(resolution)*(diff(tau2.bounds))+tau2.bounds[1])
##     }

##     ## vi <- as.tensor(vi,dims=c(K=K))
##     ## yi <- as.tensor(yi,dims=c(K=K))
##     inv.vi <- as.tensor(1/vi,dims=c(K=K))
##     sum.inv.vi <- sum(inv.vi)
##     mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
##     tau2.DL.obs <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
##     mu.DL.obs <- sum(yi/(vi+tau2.DL.obs))/sum(1/(vi+tau2.DL.obs))

##     grid.points <- expand.grid(mu,tau2)
##     N <- length(mu)*length(tau2)
##     colnames(grid.points) <- c('mu','tau2')
##     ## T1.obs <- (mu.DL.obs-grid.points$mu)^2 / sum(1/(tau2.DL.obs+vi))
##     T1.obs <- (mu.DL.obs-grid.points$mu)^2 * sum(1/(tau2.DL.obs+vi))
##     T2.obs <- outer(yi,grid.points$mu,`-`)^2 / outer(vi,grid.points$tau2,`+`) + log(outer(vi,grid.points$tau2,`+`))
##     ## T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) / 2
##     T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) #/ 2

##     inv.vi <- as.numeric(inv.vi)
##     ZZ <- to.tensor(as.numeric(Z),c(K=K,B=B))
##     w.unnormalized <- as.tensor(outer(as.numeric(vi),grid.points$tau2,`+`),dims=c(K=K,N=N))
##     w <- w.unnormalized/margin.tensor(w.unnormalized,1)

##     for (i in 1:K)
##         assign(paste0('Y',i), (outer(sqrt(w.unnormalized[i,]),ZZ[i,])+grid.points$mu))
##     mu.fixed <- Reduce(`+`,lapply(1:K,function(k)get(paste0('Y',k))*(inv.vi[k]/sum.inv.vi)))
##     tau2.DL <- Reduce(`+`,lapply(1:K,function(k)
##                                  (get(paste0('Y',k))-mu.fixed)^2*inv.vi[k]
##                                  ))
##     tau2.DL <- (tau2.DL - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
##     tau2.DL[tau2.DL<0] <- 0
##     for(k in 1:K)assign(paste0('w.hat.inv.unnormalized',k),1/(tau2.DL+vi[k]))
##     w.hat.inv.sums <- Reduce(`+`,lapply(1:K,function(k)get(paste0('w.hat.inv.unnormalized',k))))
##     for(k in 1:K)assign(paste0('w.hat.inv',k),get(paste0('w.hat.inv.unnormalized',k))/w.hat.inv.sums)
##     mu.DL <- Reduce(`+`,lapply(1:K,function(k)get(paste0('Y',k))*get(paste0('w.hat.inv',k))))

##     T1 <- (mu.DL-grid.points$mu)^2 * w.hat.inv.sums
##     class(ZZ) <- 'matrix'
##     class(w.unnormalized) <- 'matrix'
##     for(k in 1:K)
##         assign(paste0('T2',k),
##                t(-(get(paste0('Y',k)) - mu.DL)^2 * get(paste0('w.hat.inv.unnormalized',k)) + log(get(paste0('w.hat.inv.unnormalized',k))) + log(w.unnormalized[k,])) + ZZ[k,]^2
##                )##bottleneck on log call
##     T2 <- t(Reduce(`+`,lapply(1:K,function(k)get(paste0('T2',k)))))
##     ## T0.obs <- rowSums(T1.obs<=T1) + rowSums(T2.obs<=T2)/c0
##     T0.obs.a <- rowSums(T1.obs<=T1)
##     T0.obs.b <- rowSums(T2.obs<=T2)
##     for(i in 1:length(c0))
##         assign(paste0('T0.obs',i),T0.obs.a + T0.obs.b/c0[i])
##     ## T0 <- B+1-t(apply(T1,1,rank)) + (B+1-t(apply(T2,1,rank)))/c0 ##bottleneck on rank
##     T0.a <- B+1-t(apply(T1,1,rank))
##     T0.b <- (B+1-t(apply(T2,1,rank)))
##     for(i in 1:length(c0))
##         assign(paste0('T0',i),T0.a + T0.b/c0[i])

##     ans.list <- lapply(1:length(c0),function(i)
##                        structure(matrix(rowMeans(get(paste0('T0.obs',i)) >=get(paste0('T0',i))),nrow=length(mu),dimnames=list(round(mu,8),round(tau2,8))), mu=mu, tau2=tau2,
##                                  class=c('RMA.Exact','matrix'))
##                        )
##     if(length(ans.list)==1) ans.list <- ans.list[[1]]
##     ans.list
## }

## ## combining test statistics directly, allowing for non-square grid points
## ## attributes for yi,vi (needed for zoom)
## rma.exact <- function(grid.points=NULL,yi,vi,Z=NULL,c0=1,resolution=1e2,B=ncol(Z)) {
##     K <- length(yi)

##     if(is.null(Z)) Z <- matrix(rnorm(B*K),nrow=K)
##     ## B <- ncol(Z)

##     if(is.null(grid.points)) {
##         ## if(is.null(mu)) {
##         fit=rma.uni(yi=yi, vi=vi, method="DL")
##         mu.bounds <- fit$b+c(-1,1)*fit$se*qt(.99,df=K-1)
##         mu <- sort(runif(resolution)*(diff(mu.bounds))+mu.bounds[1])
##         ## }
##         ## if(is.null(tau2)) {
##         ## tau2.bounds <- tau.ci(thetahat=yi,varhat=vi,level=.99)
##         tau2.bounds <- tau2.ci(yi,vi,level=.9)
##         ## print(tau2.bounds)
##         tau2 <- sort(runif(resolution)*(diff(tau2.bounds))+tau2.bounds[1])
##         ## }

##         grid.points <- expand.grid(mu,tau2)
##         colnames(grid.points) <- c('mu','tau2')
##     }

##     N <- nrow(grid.points)
##     mu <- grid.points$mu
##     tau2 <- grid.points$tau2

##     ## vi <- as.tensor(vi,dims=c(K=K))
##     ## yi <- as.tensor(yi,dims=c(K=K))
##     inv.vi <- as.tensor(1/vi,dims=c(K=K))
##     sum.inv.vi <- sum(inv.vi)
##     mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
##     tau2.DL.obs <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
##     mu.DL.obs <- sum(yi/(vi+tau2.DL.obs))/sum(1/(vi+tau2.DL.obs))

##     ## T1.obs <- (mu.DL.obs-grid.points$mu)^2 / sum(1/(tau2.DL.obs+vi))
##     T1.obs <- (mu.DL.obs-grid.points$mu)^2 * sum(1/(tau2.DL.obs+vi))
##     T2.obs <- outer(yi,grid.points$mu,`-`)^2 / outer(vi,grid.points$tau2,`+`) + log(outer(vi,grid.points$tau2,`+`))
##     T2.obs <- (colSums(T2.obs) - sum((yi-mu.DL.obs)^2/(tau2.DL.obs+vi) + log(tau2.DL.obs+vi))) #/ 2

##     inv.vi <- as.numeric(inv.vi)
##     ZZ <- to.tensor(as.numeric(Z),c(K=K,B=B))
##     w.unnormalized <- as.tensor(outer(as.numeric(vi),grid.points$tau2,`+`),dims=c(K=K,N=N))
##     w <- w.unnormalized/margin.tensor(w.unnormalized,1)

##     for (i in 1:K)
##         assign(paste0('Y',i), (outer(sqrt(w.unnormalized[i,]),ZZ[i,])+grid.points$mu))
##     mu.fixed <- Reduce(`+`,lapply(1:K,function(k)get(paste0('Y',k))*(inv.vi[k]/sum.inv.vi)))
##     tau2.DL <- Reduce(`+`,lapply(1:K,function(k)
##                                  (get(paste0('Y',k))-mu.fixed)^2*inv.vi[k]
##                                  ))
##     tau2.DL <- (tau2.DL - (K-1)) / (sum.inv.vi - sum(inv.vi^2)/sum.inv.vi)
##     tau2.DL[tau2.DL<0] <- 0
##     for(k in 1:K)assign(paste0('w.hat.inv.unnormalized',k),1/(tau2.DL+vi[k]))
##     w.hat.inv.sums <- Reduce(`+`,lapply(1:K,function(k)get(paste0('w.hat.inv.unnormalized',k))))
##     for(k in 1:K)assign(paste0('w.hat.inv',k),get(paste0('w.hat.inv.unnormalized',k))/w.hat.inv.sums)
##     mu.DL <- Reduce(`+`,lapply(1:K,function(k)get(paste0('Y',k))*get(paste0('w.hat.inv',k))))

##     T1 <- (mu.DL-grid.points$mu)^2 * w.hat.inv.sums
##     class(ZZ) <- 'matrix'
##     class(w.unnormalized) <- 'matrix'
##     for(k in 1:K)
##         assign(paste0('T2',k),
##                t(-(get(paste0('Y',k)) - mu.DL)^2 * get(paste0('w.hat.inv.unnormalized',k)) + log(get(paste0('w.hat.inv.unnormalized',k))) + log(w.unnormalized[k,])) + ZZ[k,]^2
##                )##bottleneck on log call
##     T2 <- t(Reduce(`+`,lapply(1:K,function(k)get(paste0('T2',k)))))
##     for(i in 1:length(c0))
##         assign(paste0('T0.obs',i),T1.obs + T2.obs/c0[i])
##     for(i in 1:length(c0))
##         assign(paste0('T0',i),T1 + T2/c0[i])

##     ## ans.list <- lapply(1:length(c0),function(i)
##     ##                    structure(matrix(rowMeans(get(paste0('T0.obs',i)) <=get(paste0('T0',i))),nrow=length(mu),dimnames=list(round(mu,8),round(tau2,8))), mu=mu, tau2=tau2,
##     ##                              class=c('RMA.Exact','matrix'))
##     ##                    )
##     ## if(length(ans.list)==1) ans.list <- ans.list[[1]]
##     ans.list <- lapply(1:length(c0),function(i)
##                        structure(data.frame(mu=grid.points$mu,tau2=grid.points$tau2,p.val=rowMeans(get(paste0('T0.obs',i)) <=get(paste0('T0',i))),c0=c0[i]),
##                                  class=c('RMA.Exact','data.frame'),yi=yi,vi=vi)
##                        )
##     ans <- do.call(rbind,ans.list)
## }


## ## df <- data.frame(x=grid.points$mu,y=grid.points$tau2,z=rowMeans(get(paste0('T0.obs',i)) <=get(paste0('T0',i))))
## ## df.cnt <- getContourLines(df,nlevels=3)
## ## ggplot(df.cnt,aes(x=x,y=y,group=Group,color=z))+geom_path()+theme_bw()


## ## helper to run base contourLines if grid of p-values is available, otherwise contoureR
## contourLines.RMA.Exact <- function(rma0,levels) {
##     mu <- unique(rma0$mu) %>% sort
##     tau2 <- unique(rma0$tau2) %>% sort
##     if(length(mu)*length(tau2)==nrow(rma0)) {
##         rma0.mtx <- matrix(rma0$p.val, nrow=length(mu))
##         contour.lst <- contourLines(mu,tau2,rma0.mtx,levels=levels)
##         return(do.call(rbind, lapply(contour.lst,function(df)data.frame(x=df$x,y=df$y,z=df$level))))
##     } else {
##         return(getContourLines(as.data.frame(select(rename(rma.lower,x=mu,y=tau2,z=p.val),x,y,z)),levels=levels) %>% select(-LID,-GID,-PID,-Group))
##     }
## }

## zoom <- function(rma0,alpha=.05,scale=1.5,resolution=5e1,plot=TRUE) {
##     yi <- attr(rma0,'yi')
##     vi <- attr(rma0,'vi')
##     sd.mu <- sd(rma0$mu)
##     sd.tau2 <- sd(rma0$tau2)
##     CI <- confint(rma0,alpha=alpha)
##     stopifnot(nrow(CI)==1)
##     tau2.lower <- pmax(0,runif(5e1)*scale*sd.tau2 + (CI[,'tau2.lower']-scale*sd.tau2/2)) %>% unique %>% sort
##     mu.lower <- (runif(5e1)*scale*sd.mu + (CI[,'mu.lower']-scale*sd.mu/2)) %>% sort
##     rma.lower <- rma.exact(grid.points=expand.grid(mu=mu.lower,tau2=tau2.lower),yi=yi,vi=vi,B=1e3)
##     tau2.upper <- pmax(0,runif(5e1)*scale*sd.tau2 + (CI[,'tau2.upper']-scale*sd.tau2/2)) %>% unique %>% sort
##     mu.upper <- (runif(5e1)*scale*sd.mu + (CI[,'mu.upper']-scale*sd.mu/2)) %>% sort
##     rma.upper <- rma.exact(grid.points=expand.grid(mu=mu.upper,tau2=tau2.upper),yi=yi,vi=vi,B=1e3)


##     if(plot) {
##         ## upper.lines <- contourLines.RMA.Exact(as.data.frame(select(rename(rma.upper,x=mu,y=tau2,z=p.val),x,y,z)),levels=alpha) %>% mutate(group='upper')
##         ## ## plot(y ~ x, data=cc)
##         ## lower.lines <- contourLines.RMA.Exact(as.data.frame(select(rename(rma.lower,x=mu,y=tau2,z=p.val),x,y,z)),levels=alpha) %>% mutate(group='lower')
##         upper.lines <- contourLines.RMA.Exact(rma.upper,levels=alpha) %>% mutate(group='upper')
##         ## plot(y ~ x, data=cc)
##         lower.lines <- contourLines.RMA.Exact(rma.lower,levels=alpha) %>% mutate(group='lower')
##         new.lines <- rbind(upper.lines,lower.lines) %>% mutate(p.val=paste0('new ',alpha))
##         plt <- plot(rma0)
##         nlevels <- length(unique(plt$data$p.val))

##         suppressMessages(plt <- plt+geom_path(data=new.lines,aes(x=x,y=y,group=group,color=p.val))+scale_color_manual(values=c(brewer.pal(nlevels,'PuBu'),'red')))
##         print(plt)
##     }

##     ## plot(y~x,data=lower.lines)
##     ## new.lines <- getContourLines(rbind(as.data.frame(select(rename(rma.upper,x=mu,y=tau2,z=p.val),x,y,z)),as.data.frame(select(rename(rma.lower,x=mu,y=tau2,z=p.val),x,y,z))),levels=.05)
##     ## plot(y~x,data=new.lines)

##     ## lower.lines <- contourLines(mu,tau2,matrix(rma.lower$p.val,nrow=length(mu)),levels=alpha)
##     ## lower.lines <- mutate(as.data.frame(lower.lines),p.val=paste0('new ',alpha),group='upper')
##     ## ## plt <- plt+geom_path(data=lower.lines,aes(x=x,y=y,group=p.val,color=p.val))+scale_color_manual(values=c(cm.colors(6),'red'))
##     ## upper.lines <- contourLines(mu.upper,tau2,matrix(rma.upper$p.val,nrow=length(unique(rma.upper$mu))),levels=alpha)
##     ## upper.lines <- mutate(as.data.frame(upper.lines),p.val=paste0('new ',alpha),group='lower')
##     ## new.lines <- rbind(upper.lines,lower.lines)
##     ## ## invisible(plt <- plt+geom_path(data=new.lines,aes(x=x,y=y,group=group,color=p.val))+scale_color_manual(values=c(cm.colors(6),'red')))
##     ## ## return(plt)
##     ## plt <- plot(rma0)
##     ## nlevels <- length(unique(plt$data$p.val))

##     ## suppressMessages(plt <- plt+geom_path(data=new.lines,aes(x=x,y=y,group=group,color=p.val))+scale_color_manual(values=c(cm.colors(nlevels),'red')))
##     ## plt

##     return(rma.lower+rma.upper)
## }
