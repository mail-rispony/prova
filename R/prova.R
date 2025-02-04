a <- function(){
  
  
  bspline<-function(x, ndx, xlr=NULL, knots=NULL, deg=3, deriv=0){
    #x: vettore di dati
    #xlr: il vettore di c(xl,xr)
    #ndx: n.intervalli in cui dividere il range
    #deg: il grado della spline
    require(splines)
    if(is.null(knots)){
      if(is.null(xlr)){
        xl<-min(x)-.01*diff(range(x))
        xr<-max(x)+.01*diff(range(x))
      } else {
        if(length(xlr)!=2) stop("quando fornito, xlr deve avere due componenti")
        xl<-xlr[1]
        xr<-xlr[2]
      }
      dx<-(xr-xl)/ndx
      knots<-seq(xl-deg*dx,xr+deg*dx,by=dx)
    }
    B<-splineDesign(knots,x,ord=deg+1,derivs=rep(deriv,length(x)))
    B<-list(B=B, knots=knots)
    B#the B-spline base matrix
  }#end_fn
  
  blockdiag <- function(...) {
    args <- list(...)
    nc <- sapply(args,ncol)
    cumnc <- cumsum(nc)
    ##  nr <- sapply(args,nrow)
    ## NR <- sum(nr)
    NC <- sum(nc)
    rowfun <- function(m,zbefore,zafter) {
      cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
            matrix(0,ncol=zafter,nrow=nrow(m)))
    }
    ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
    for (i in 2:length(args)) {
      ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
    }
    ret
  }
  ##########################################
  my.lm.fit<-function(y,x,n.int,deg, d, spar=0, xlr=NULL, 
                      quant=FALSE, plot.it=0, conf.level=0, ...){
    #n.int: no. of intervals
    #quant: if TRUE, the knots are computed via quantiles
    #   otherwise equally-spaced values
    #plot.it: 0 =no plot, 1= plots data and fitted values
    #         2: adds fitted values only
    #... arguments to pass to points()
    n<-length(y)
    if(quant){
      pr<-seq(0,1,l=n.int+1)
      k<- quantile(x, probs=pr)
    } else {
      #equally spaced knots
      k<- seq(min(x), max(x), l=n.int+1)
    }
    #browser()
    #B<-model.matrix(~0+factor(cut(x,
    #                             k,include.lowest=T, labels=F)))
    B <- bspline(x, ndx=n.int, deg=deg, xlr=xlr)
    knots<-B$knots
    B<-B$B
    #sbrowser()
    # B=NULL
    # for(i in 1:(length(k)-1)) {
    #   B[[length(B)+1]]<-1*(x>=k[i] & x<k[i+1])
    # }
    # B<-matrix(unlist(B), n, length(B))
    #o=lm.fit(y=y,x=B)
    #b<-o$coefficients
    
    #browser()
    
    D<-diff(diag(ncol(B)), dif=d)
    #Xs<- rbind(B, sqrt(spar)*D)
    #ys<-c(y, rep(0, nrows(D)))
    #lm.fit(y=ys, x=Xs)
    P<-crossprod(D)
    BtB=crossprod(B)
    invBtB.spar=solve(BtB+spar*P)
    b=drop(invBtB.spar%*%t(B)%*%y)
    fitted.values=drop(B%*%b)
    edf<- sum(diag(BtB %*% invBtB.spar))
    s2<-sum((y-fitted.values)^2)/(n-edf)
    var.b<- s2*invBtB.spar %*% BtB %*% invBtB.spar
    
    if(plot.it>0){
      xnew<- seq(min(x,xlr+.001),max(x,xlr-.001),l=300)
      #Bnew<-model.matrix(~0+factor(cut(xnew,
      #                              k,include.lowest=T, labels=F)))
      #matplot(xnew,  Bnew%*%diag(o$coefficients), type="l", lwd=3)
      Bnew<- bspline(xnew, ndx=n.int, deg=deg, xlr=xlr)#, knots=knots )
      Bnew<-Bnew$B
      hat.mu.new<- drop(Bnew%*%b)
      if(plot.it==1) plot(x,y, xlim=range(knots))
      lines(xnew, hat.mu.new, ...)
      #browser()
      if(conf.level>0){
        #var.hat.mu.new<-Bnew %*% var.b %*% t(Bnew)
        #se.hat.mu.new <- sqrt(diag(var.hat.mu.new))
        se.hat.mu.new <-sqrt(rowSums((Bnew%*%var.b)*Bnew))
        z<- -qnorm((1-conf.level)/2)
        inf<- hat.mu.new-z*se.hat.mu.new
        sup<- hat.mu.new+z*se.hat.mu.new
        matlines(xnew, cbind(inf, sup), lty=3,...)
        
      }
      #abline(v=knots, lty=3)
      #browser()
    }
    hatsigma<- sqrt(sum((y-fitted.values)^2)/n)
    ll<- sum(dnorm(y, fitted.values, sd=hatsigma,log=TRUE))
    aic= -2*ll + 2*edf
    bic= -2*ll + edf*log(n)
    #browser()
    r<-list(edf=edf, aic=aic, bic=bic, ress=y-fitted.values,
            pen=drop(t(b)%*%P%*%b), s2=s2)
    r
  }
  
  ##########################################
  
  
  my.glm.fit<-function(y,x,n.int,deg, d, spar=0, xlr=NULL, family= ExpoFamily("gaussian"),
                       quant=FALSE, plot.it=0, conf.level=0, P=0, beta_start=NULL, ...){
    #n.int: no. of intervals
    #quant: if TRUE, the knots are computed via quantiles
    #   otherwise equally-spaced values
    #plot.it: 0 =no plot, 1= plots data and fitted values
    #         2: adds fitted values only
    #... arguments to pass to points()
    n<-length(y)
    if(quant){
      pr<-seq(0,1,l=n.int+1)
      k<- quantile(x, probs=pr)
    } else {
      #equally spaced knots
      k<- seq(min(x), max(x), l=n.int+1)
    }
    #browser()
    #B<-model.matrix(~0+factor(cut(x,
    #                             k,include.lowest=T, labels=F)))
    B <- bspline(x, ndx=n.int, deg=deg, xlr=xlr)
    knots<-B$knots
    B<-B$B
    #browser()
    
    D<-diff(diag(ncol(B)), dif=d)
    #Xs<- rbind(B, sqrt(spar)*D)
    #ys<-c(y, rep(0, nrows(D)))
    #lm.fit(y=ys, x=Xs)
    P<-crossprod(D)
    
    o=IWLSp(B, y, P, beta_start, family, 
            nsteps = 100L, eps = 1e-6, verbose = FALSE)
    o
    
  }
  
  ######################################################
  
  
  ExpoFamily <- function(distribution = c("gaussian", "poisson", "binomial", "inverse", "gamma"),
                         link) {
    distribution <- match.arg(distribution)
    if (missing(link)) {
      link <- switch(distribution,
                     gaussian = "identity",
                     poisson = "log",
                     binomial = "logit",
                     inverse = "inverse",
                     gamma = "inverse")
    } else {
      link_validi <- c("identity", "log", "logit", "inverse", "sqrt")
      link <- match.arg(link, link_validi)
    }
    
    mu_fun <- switch(link,
                     identity = function(eta) eta,
                     log = function(eta) exp(eta),
                     sqrt = function(eta) eta ^ 2,
                     logit = function(eta) 1 / (1 + exp(-eta)),
                     inverse = function(eta) 1 / eta)
    
    dmu_deta <- switch(link,
                       identity = function(eta) rep(1, length(eta)),
                       log = function(eta) exp(eta),
                       sqrt = function(eta) 2 * eta,
                       logit = function(eta) exp(-eta) / (1 + exp(-eta))^2,
                       inverse = function(eta) -1 / (eta ^ 2))
    
    var_fun <- switch(distribution,
                      gaussian = function(mu) rep(1, length(mu)),
                      poisson = function(mu) mu,
                      binomial = function(mu) mu * (1 - mu),
                      inverse = function(mu) mu ^ 3,
                      gamma = function(mu) mu ^ 2)
    
    deviance_residual <- function(y, mu, family) {
      mu <- pmax(mu, .Machine$double.eps)
      
      switch(family$distribution,
             gaussian = sum((y - mu) ^ 2),
             binomial = -2 * sum(y * log(ifelse(mu == 0, 1, mu)) + (1 - y) * log(ifelse(1 - mu == 0, 1, 1 - mu))),
             poisson = 2 * sum(ifelse(y == 0, 0, y * log(y / mu)) - (y - mu)),
             gamma = -2 * sum(log(ifelse(y == 0, 1, y / mu)) - (y - mu) / mu),
             inverse = sum((y - mu)^2 / (y * (mu^2)))
      )
    }
    
    list(mu_fun = mu_fun, var_fun = var_fun, dmu_deta = dmu_deta,
         deviance_residual = deviance_residual, distribution = distribution)
  }
  ######################################################
  
  IWLSp <- function(X, y, P, beta_start, family, 
                    nsteps, eps, verbose = FALSE) {
    mu_fun <- family$mu_fun
    var_fun <- family$var_fun
    dmu_deta <- family$dmu_deta
    distribution <- family$distribution
    eta <- drop(X %*% beta_start)
    mu <- mu_fun(eta)
    residuals <- y - mu  
    for (i in seq_len(nsteps)) {
      dmu <- dmu_deta(eta)
      V <- var_fun(mu)
      y_tilde <- eta + ((y - mu) / dmu)
      W <- dmu^2 / V
      H <- crossprod(X, W * X) +P
      alpha <- crossprod(X, W * y_tilde)
      beta_new <- solve(H, alpha)
      Norm2_db <- sqrt(sum((beta_new - beta_start)^2))
      if (verbose) cat("Step", i, " ||db|| = ", Norm2_db, "\n")
      if (Norm2_db <= eps) {
        conv <- TRUE
        break
      } else {
        conv <- FALSE
        beta_start <- beta_new
      }
      
      eta <- drop(X %*% beta_start)
      mu <- mu_fun(eta)
      residuals <- y - mu  
    }
    
    dmu_final <- dmu_deta(eta)
    V_final <- var_fun(mu)
    W_final <- dmu_final^2 / V_final
    H_final <- crossprod(X, W_final * X)+P
    hatM <- crossprod(X, W * X) %*% solve(H_final)
    list(beta = beta_start, H = H_final, conv = conv, fitted.values = mu, 
         hatM=hatM, residuals = residuals, V_final = V_final, distribution = distribution)
  }
  
  ######################################################
  
  
  
  
  
  
  #==============================================
  
  my.fit.lm.VC<-function(y, Xsmooth, Zint=NULL, 
                         Xlin=NULL, spar, intc=TRUE, ndx, deg, d) { 
    n<-length(y)
    Xsmooth<- as.matrix(Xsmooth)
    if(missing(spar)){
      spar<-rep(0, ncol(Xsmooth))
    } else {
      if(length(spar)!=ncol(Xsmooth)) stop("smdsjdm")
    }
    #X<-cbind(Xsmooth, Xlin)
    p<-B<-P<-Bplot<-NULL
    #browser()
    if(is.null(Zint)) {
      id.vc=FALSE
      Zint<-matrix(1, n, ncol(Xsmooth))
    } else {
      id.vc=TRUE
      Zint<-as.matrix(Zint)
      if(ncol(Zint)!=ncol(Xsmooth)) stop("if provided, Zint... ")
    }
    #sbrowser()
    C<-contr.treatment(ndx+deg) #ndx+deg
    for(j in 1:ncol(Xsmooth)){
      Bplot[[length(Bplot)+1]]<- BB<-bspline(Xsmooth[,j], ndx=ndx, deg=deg)$B[,-1]
      B[[length(B)+1]]<- diag(Zint[,j])%*%BB
      p[[length(p)+1]] <- ncol(B[[j]])
      P[[length(P)+1]]<- spar[j]*crossprod(diff(diag(p[[j]]+1), d=d)%*%C)
    }
    
    if(!is.null(Xlin)) {
      Xlin<-cbind(1, as.matrix(Xlin))
    } else {
      Xlin<- matrix(1, nrow=n, ncol=1)
    }
    if(id.vc) {
      id.ok=which(colSums(Zint)!=n)
      Zint<-Zint[,id.ok]
      Xlin<-  cbind(Xlin, Zint)
    }
    
    X<-cbind(Xlin, do.call(cbind, B))
    P<-c(list(matrix(0, ncol(Xlin), ncol(Xlin))), P)
    P=do.call(blockdiag, P)
    if(!intc) {
      X<-X[,-1]
      P<-P[-1,-1]
    }
    
    b=drop(solve(crossprod(X)+P,crossprod(X,y)))
    
    var_b <- var(y)*solve(crossprod(X)+P,crossprod(X,X))%*%solve(crossprod(X)+P)
    se_b <- sqrt(diag(var_b))
    
    stat_test <- vector("numeric", length=ncol(X))
    pvalue <- vector("numeric", length=ncol(X))
    for (j in 1:ncol(X)) {
      stat_test[j] <- b[j]/se_b[j]
      pvalue[j] <- 2*(1-pt(stat_test[j], df=nrow(X)-ncol(X)))
    }
    hat.mu<- drop(X%*%b)
    H=crossprod(X) %*% solve(crossprod(X)+P)
    edf.val<- diag(H)
    pen <- t(b)%*%P%*%b/(length(y)-sum(edf.val))
    #browser()
    start<-0
    edfL<-fitj<- vector("list", ncol(Xsmooth))
    e=y-hat.mu
    var_e <- (t(e)%*%e)/(length(y)-sum(edf.val))
    #browser()
    for(j in 1:ncol(Xsmooth)){
      id.ok <- ncol(Xlin)+(1:p[[j]])+start -(!intc)
      start<-start+ p[[j]]
      bj<- b[id.ok]
      muj<-drop(Bplot[[j]]%*%bj)
      muj<-muj[order(Xsmooth[,j])]
      fitj[[j]]<-cbind(sort(Xsmooth[,j]), muj, muj+e)
      edfL[[j]]<- sum(edf.val[id.ok])
    }
    #  browser()
    r<-list(pen= pen,var_e=var_e, b=cbind(b,pvalue), mu=hat.mu, part.eff=fitj, e=e, edfL=edfL, edf=sum(edf.val))
    class(r)<-"nonpar24"
    r
  }
  
  
  
  
  
  
  str(dd)
  
  y <- dd$n.dec
  z <- dd$INFLU
  z1 <- ifelse(dd$INFLU==1,1,0) #influenza si 
  z2 <- ifelse(dd$INFLU==0, 1, 0)#influenza no
  x1 <- dd$PM10
  x2 <- dd$MTSPL03
  x3 <- dd$giorno
  
  par(mfrow=c(2,1))
  plot(x1, y, ylab="N.dec", xlab="PM10")
  plot(x2,y, ylab="N.dec", xlab="MTSPL10")
  boxplot(dd$n.dec~dd$INFLU)
  # dal boxplot si evince ch wein termini di mediana si hanno piu decessi in presenza di influenza 
  
  hfs <- function(spar0, max_it){
    for (i in 1:max_it) {
      #browser()
      cat(spar0,"\n")
      m1 <- my.fit.lm.VC(y, Xsmooth = x3#, Zint = cbind(z1, z2)
                         , Xlin = cbind(z1,x1)
                         , spar=spar0, ndx=10, deg=3,d=2, intc = T)
      spar1 <- drop(m1$var_e)/drop(m1$pen)
      if(abs(spar1-spar0)<0.001) break
      spar0<- drop(spar1)
    }
    r <- my.fit.lm.VC(y, Xsmooth = x3#, Zint = cbind(z1, z2)
                      , Xlin = cbind(z1,x1)
                      , spar=spar1, ndx=10, deg=3,d=2, intc = T)
    r
    
  }
  
  
  m1 <- my.fit.lm.VC(y, Xsmooth = x3#, Zint = cbind(z1, z2)
                     , Xlin = cbind(z1,x1)
                     , spar=0.1, ndx=10, deg=3,d=2, intc = T)
  m1$b
  
  eff_part1 <- m1$part.eff[[1]]
  mu_part1 <- eff_part1[,2]
  pseudo_part1 <- eff_part1[,3]
  
  
  # plot effetto parziale 1
  plot(sort(x3),pseudo_part1[order(sort(x3))], xlab="PM10", ylab="Pseudo dati")
  lines(sort(x3), mu_part1[order(sort(x3))], col=2, lwd=2, xlab="PM10", ylab="Pseudo dat")
  
  a <-hfs(spar0 = 0.1, max_it=20)
  a
  
  
  eff_part1 <- a$part.eff[[1]]
  mu_part1 <- eff_part1[,2]
  pseudo_part1 <- eff_part1[,3]
  
  
  
  plot(sort(x3),pseudo_part1[order(sort(x3))], xlab="PM10", ylab="Pseudo dati")
  lines(sort(x3), mu_part1[order(sort(x3))], col=2, lwd=2, xlab="PM10", ylab="Pseudo dat")
  
  
}

b <- function(){
  
  bspline<-function(x, ndx, xlr=NULL, deg=3, deriv=0){
    #x: vettore di dati
    #xlr: il vettore di c(xl,xr)
    #ndx: n.intervalli in cui dividere il range
    #deg: il grado della spline
    require(splines)
    if(is.null(xlr)){
      xl<-min(x)-.01*diff(range(x))
      xr<-max(x)+.01*diff(range(x))
    } else {
      if(length(xlr)!=2) stop("quando fornito, xlr deve avere due componenti")
      xl<-xlr[1]
      xr<-xlr[2]
    }
    dx<-(xr-xl)/ndx
    knots<-seq(xl-deg*dx,xr+deg*dx,by=dx)
    B<-splineDesign(knots,x,ord=deg+1,derivs=rep(deriv,length(x)))
    B<-list(B=B, knots=knots)
    B#the B-spline base matrix
  }#end_f
  
  bspline.fit <- function(x1,x2,x3,x,y,ndx=30,p=3,d=2,lambda=0.1,xlr=NULL,plotdata=T){
    B <- bspline(x=x,ndx=ndx,deg = p,xlr=xlr)$B[,-1]
    L <- outer(x1,0:1,"^")
    L1 <-outer(x2,1,"^")
    L2 <- outer(x3,1,"^")
    Linear<- cbind(L,L1,L2)
    nc<- ncol(Linear)
    X <- cbind(Linear,B)
    
    D <- diff(diag(ncol(B)), d=d)
    P <- lambda*crossprod(D)
    
    zeros<- matrix(0,nrow=nc,ncol=nc)
    require(Matrix)
    
    P1 <- as.matrix(bdiag(zeros,P))
    
    id1<- 1:ncol(L)
    id2 <- (ncol(L)+1):(ncol(L1)+2)
    id3 <- (ncol(L)+2):(ncol(L1)+3)
    id4 <- (ncol(L2)+4):ncol(X)
    
    invXtX.P <- solve(crossprod(X)+P1)
    b <- drop(invXtX.P %*% crossprod(X,y))
    
    b.B <- b[id4]
    
    fit.mu <- X %*% b
    e <- y - fit.mu
    
    fit.mu.L <- X[,id1] %*% b[id1] 
    fit.mu.L1 <- X[,c(1,id2)] %*% b[c(1,id2)]
    fit.mu.L2 <- X[,c(1,id3)] %*% b[c(1,id3)]
    fit.mu.B <- X[,c(1,id4)] %*% b[c(1,id4)] 
    
    # edf
    XtX <- crossprod(X)
    H <- invXtX.P %*% XtX
    edf <- sum(diag(H))
    edf.L <- sum(diag(H)[id1])
    edf.L1 <- sum(diag(H)[id2])
    edf.L2 <- sum(diag(H)[id3])
    edf.B <- sum(diag(H)[id4])
    
    s2 <- sum(e^2)/(length(y)-edf)
    li <- dnorm(y,fit.mu,sqrt(s2 * (length(y)-edf)/length(y)), log = T)
    AIC <- -2*sum(li)+2*edf
    
    if(plotdata){
      par(mfrow = c(1, 3))
      m.resp.L <- e + fit.mu.L
      plot(sort(x1),m.resp.L[order(x1)],lty=2,lwd=3,col="brown",main = "Reletionship between hours\n and income")
      lines(sort(x1),fit.mu.L[order(x1)],lwd=2,col='darkgreen')
      
      m.resp.L1 <-e + fit.mu.L1
      plot(sort(x2),m.resp.L1[order(x2)],lty=2,lwd=3,col="brown",main = "Reletionship between hours\n and age")
      lines(sort(x2),fit.mu.L1[order(x2)],lwd=2,col='darkgreen')
      
      #m.resp.L2 <-e+fit.mu.L2
      #plot(sort(x3),m.resp.L2[order(x3)])
      #lines(sort(x3),fit.mu.L2,lty=2,col='darkgreen')
      
      res.p.B <-e + fit.mu.B
      plot(sort(x), res.p.B[order(x)],lty=2,lwd=3,col="brown",main = "Reletionship between hours\n and education")
      lines(sort(x), fit.mu.B[order(x)],lwd=2,col='darkgreen')
      
    }
    
    r<-list(coeff=b, b.B=b.B, mu=fit.mu, mu.B=fit.mu.B,e = e, s2=s2, edf=edf, edf.B = edf.B, aic=AIC,lambda=lambda)
    r
    
    
    
  }
  
  hfs.alg <- function(x1,x2,x3,x,y,lambda0=0.1,itmax=20,d=2,p=3,ndx=30,xlr=NULL){
    for(i in 1:itmax){
      o<- bspline.fit(x1=x1,x2=x2,x3=x3,x=x,y=y,d=d,p=p,lambda = lambda0,ndx =ndx,xlr=xlr,plotdata =F )
      sigma.s2.e <- o$s2
      sigma.b <- sum(diff(o$b.B, d=d)^2)/o$edf.B
      lambda1 <- sigma.s2.e/sigma.b
      if((abs(lambda1-lambda0)< .001||lambda1 > 1e6 )) break
      lambda0 <- lambda1
      cat("niter", i, "lambda", lambda0,  "\n")
    }
    
    r <- bspline.fit(x1=x1,x2=x2,x3=x3,x=x, y=y, ndx=ndx, lambda=lambda1, d=d, p=p,plotdata = T)
    
  }
  library(tidyverse)
  dati<-read_delim(file.choose(),delim = ",")
  dati1<-dati[,-c(1,6,7,8,10,11,12,13)]
  
  #education as smoother
  x<-education
  x1<-income
  x2<-age
  x3<-nonwhite
  y<-hours
  
  summary(dati1)
  str(dati1)
  
  bspline.fit(x1,x2,x3,x,y)
  hfs.alg(x1,x2,x3,x,y)  
  
  #income as smoother
  x1<-education
  x<-income
  x2<-age
  x3<-nonwhite
  y<-hours
  
  bspline.fit(x1,x2,x3,x,y)
  hfs.alg(x1,x2,x3,x,y) 
  
  #age as smoother
  x2<-education
  x1<-income
  x<-age
  x3<-nonwhite
  y<-hours
  
  bspline.fit(x1,x2,x3,x,y)
  hfs.alg(x1,x2,x3,x,y) 
  #chiaramente occorre aggiustare i titoli dei plot
  
}
