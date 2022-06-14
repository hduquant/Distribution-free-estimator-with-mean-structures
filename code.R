library(mvtnorm)
library(MASS)
library(matrixcalc)
library(doParallel)
library(expm)
library(semTools)
library(lavaan)

methodlist<-c("gammaN.model")

  s=k=1
  a=0
  
  ms=m*(m+1)/2
  ps=p*(p+1)/2
  q=2*m+2*p+ms+ps+1+p+1
  p2=p*2
  pks=(p+p*k+s)*(p+p*k+s+1)/2
  df=(2*p+1)*(2*p+2)/2+2*p+1-q

syntax <- function(Timepoint){
  
  i.com      <- paste("1*t",1:Timepoint,
                      sep = "", 
                      collapse = " + ")
  
  i.int      <- "i ~ iint*1+beta1*x"
  i.equation <- paste("i =~ ", i.com, sep = "")
  
  s.com      <- paste(0:(Timepoint-1), "*t", 1:Timepoint, 
                      sep = "", 
                      collapse = " + ")
  
  s.int      <- "s ~ sint*1+beta2*x"
  
  s.equation <- paste("s =~ ", s.com, sep = "")
  
  vcov.equation  <- paste("i ~~ s2i*i",
                          "s ~~ s2s*s",
                          "i~~is*s", sep = "\n")
  
  r.equation <- paste("t",1:Timepoint, "~~t", 1:Timepoint, 
                      sep = "", 
                      collapse = '\n')
  
  xvarying.equation <- paste("t",1:Timepoint, "~xv", 1:Timepoint, 
                             sep = "", 
                             collapse = '\n')
  
  xvarying.equation2=NULL
  for (i in 1:Timepoint){
    
    xvarying.equation2[[i]] <- paste("xv",1:(Timepoint-i+1), "~~xv", i:Timepoint, 
                                     sep = "", 
                                     collapse = '\n')
  }
  
  xvarying.equation2 <- paste( c(xvarying.equation2),  sep = "", 
                               collapse = '\n')
  
  x.equation  <- paste("x ~~ x", sep = "\n")
  
  x.equation2 <- paste("x", "~~0*xv", 1:Timepoint, 
                       sep = "", 
                       collapse = '\n')
  
  x.inter  <- paste("x ~ 1", sep = "\n")
  xv.inter <- paste("xv",1:Timepoint,  "~1", 
                    sep = "", 
                    collapse = '\n')
  
  model <- paste(i.equation, i.int,
                 s.equation, s.int,
                 xvarying.equation,xvarying.equation2, 
                 x.equation, x.equation2,
                 xv.inter, x.inter,
                 vcov.equation, r.equation,sep = "\n")
  
  return(model)
}

la.model <- syntax(p)

# DATA is the simulated data
 DATA <- cbind(y,xall)
 z<-DATA
 
 si1=t(z)
 si2=apply(z, 1, function(data) vech((data-colMeans(z))%*%t(data-colMeans(z))))
 si<-rbind(si1,si2)
 
 #get sample covariance and sample mean 
 ps <- p*(p+1)/2
 Scov <-cov(cbind(y,xall))
 vscov <- vech(Scov[1:p, 1:p])
 zbar <- colMeans(DATA)
 
 # this is the s in the pdf tile
 betav <- rbind(as.matrix(zbar),as.matrix(vscov))
 
 # define raw data in the same order as used by Scov
 DATA <- cbind(y,xall)
 
 mean_t <- apply(DATA,2,mean)
 z_c <- DATA-matrix(1,n,1)%*%mean_t
 k <- 1
 sigmaele <- matrix(NA,nrow=(2*p+1)*(2*p+2)/2,ncol=n)
 for(i in 1:(2*p+1)){
   for(j in i:(2*p+1)){
     sigmaele[k,] <- z_c[,i]*z_c[,j]
     k=k+1
   }
 }
 
 ## First, let's calculate V22, the covariance of two covariances.
 sigmaijkl <- c()
 gammaadf1 <- matrix(NA,nrow=ps+p+p^2+(p+1)*(p+2)/2,ncol=ps+p+p^2+(p+1)*(p+2)/2)
 k=1
 for(i in 1:(ps+p+p^2+(p+1)*(p+2)/2)){
   for(j in i:(ps+p+p^2+(p+1)*(p+2)/2)){
     sigmaijkl[k] <- sum(sigmaele[i,]*sigmaele[j,])/n
     gammaadf1[i,j] <- sigmaijkl[k]-sum(sigmaele[i,])*sum(sigmaele[j,])/n^2
     gammaadf1[j,i] <- gammaadf1[i,j]
     k=k+1
   }
 }
 
 
 kk=1
 index <- matrix(NA,nrow=p+1+p,ncol=p+1+p)
 for(i in 1:(2*p+1)){
   for(j in i:(2*p+1)){
     index[i,j]=index[j,i]=kk
     kk=kk+1
   }}
 
 # gammaadf1u is Browne's unbiased DF estimator
 gammaadf1u=matrix(NA,nrow=ps+p+p^2+(p+1)*(p+2)/2,ncol=ps+p+p^2+(p+1)*(p+2)/2)
 for(i in 1:(2*p+1)){
   for(j in i:(2*p+1)){
     for (k in  1:(2*p+1)){
       for (l in 1:(2*p+1)){
         if ( index[k,l]>=index[i,j] & l>=k){
           #    print(c(i,j,k,l))
           gammaadf1u[index[k,l],index[i,j]] <- n*(n-1)/(n-2)/(n-3)*(sum(sigmaele[index[k,l],]*sigmaele[index[i,j],])/n-
                                                                       sum(sigmaele[index[k,l],])*sum(sigmaele[index[i,j],])/n^2)-
             n/(n-2)/(n-3)*(sum(sigmaele[index[k,i],])*sum(sigmaele[index[l,j],])/n^2+
                              sum(sigmaele[index[k,j],])*sum(sigmaele[index[i,l],])/n^2-
                              2/(n-1)*sum(sigmaele[index[k,l],])*sum(sigmaele[index[i,j],])/n^2 )
           gammaadf1u[index[i,j],index[k,l]] <- gammaadf1u[index[k,l],index[i,j]]
         }
       }
     }
   }
 }
 
 ## Second, let's calculate V12, the covariance of a mean and a covariance.
 # gammaadf2 is Peter's ADF estimator
 xy_c<- z_c
 gammaadf2=matrix(NA,nrow=2*p+1,ncol=ps+p+p^2+(p+1)*(p+2)/2)
 for(i in 1:(2*p+1)){
   for(j in 1:(ps+p+p^2+(p+1)*(p+2)/2)){
     gammaadf2[i,j] <- sum(xy_c[,i]*sigmaele[j,])/n
   }
 }
 
 # gammaadf2u is my unbiased DF estimator
 indexbi=t(combn(n, 2))
 xy <- cbind(y,xall)
 gammaadf2u <- matrix(NA,nrow=2*p+1,ncol=ps+p+p^2+(p+1)*(p+2)/2)
 for(i in 1:(2*p+1)){
   for(k in 1:(2*p+1)){
     for (l in k:(2*p+1)){
       k.rst <- mean(xy[,i]*xy[,k]*xy[,l])
       k.r.st <- (sum(xy[indexbi[,1],i]*xy[indexbi[,2],k]*xy[indexbi[,2],l])+
                    sum(xy[indexbi[,2],i]*xy[indexbi[,1],k]*xy[indexbi[,1],l]))/n/(n-1)
       k.s.rt <- (sum(xy[indexbi[,1],k]*xy[indexbi[,2],i]*xy[indexbi[,2],l])+
                    sum(xy[indexbi[,2],k]*xy[indexbi[,1],i]*xy[indexbi[,1],l]))/n/(n-1)
       k.t.rs <- (sum(xy[indexbi[,1],l]*xy[indexbi[,2],k]*xy[indexbi[,2],i])+
                    sum(xy[indexbi[,2],l]*xy[indexbi[,1],k]*xy[indexbi[,1],i]))/n/(n-1)
       k.t.s.r <- (sum(xy[,i])*sum(xy[,k])*sum(xy[,l])- k.r.st*n*(n-1)- k.s.rt*n*(n-1)-
                     k.t.rs*n*(n-1)-k.rst*n)/n/(n-1)/(n-2)
       gammaadf2u[i,index[k,l]] <- ( k.rst-k.r.st- k.s.rt-  k.t.rs + 2*k.t.s.r)*(n-1)/n
     }
   }
 }
 
 # Third, let's calculate V11, the covariance of two means. The ADF and unbiased DF estimators are the same.
 gammaadf3 <- Scov*(n-1)/n
 
 # Finally, we calculate gamma_ADF and gamma_DF^unbiased
 
 # gammaadf is gamma_ADF
 gammaadf <- bdiag(  gammaadf3,  gammaadf1)
 gammaadf[1:(2*p+1),(2*p+2):ncol(gammaadf)] <- gammaadf2
 gammaadf[(2*p+2):ncol(gammaadf),1:(2*p+1)] <- t(gammaadf2)
 gammaadf <- as.matrix(gammaadf)
 
 # gammaadfu is gamma_DF^unbiased
 gammaadfu <- bdiag(  gammaadf3,  gammaadf1u)
 gammaadfu[1:(2*p+1),(2*p+2):ncol(gammaadfu)] <- gammaadf2u
 gammaadfu[(2*p+2):ncol(gammaadfu),1:(2*p+1)] <- t(gammaadf2u)
 gammaadfu <- as.matrix(gammaadfu)
 
#### Above is the SAME for each method  
###############
############################################# 
  
 datalavaan<-cbind(y,x,xv)
 colnames(  datalavaan)<-c(paste0("t",c(1:p)),"x",paste0("xv",c(1:p)))
 
 fites=try(lavaan(la.model, 
                  data = datalavaan, estimator = "DLS",
                  verbose = FALSE,
                  estimator.args = list(dls.a = 1, dls.GammaNT = "model"),
                  fixed.x=FALSE,conditional.x = FALSE,se = "standard"))

  ######
     theta0=coef(fites)
     vdsig<-lavTech(fites,"Delta")[[1]]
     Sigma <- lavTech(fites, "implied")[[1]]$cov
     GammaNT <- 2 * lavaan:::lav_matrix_duplication_ginv_pre_post(Sigma %x% Sigma)
     gammaNm<-bdiag(Sigma,GammaNT)
     gammaNm<-as.matrix(gammaNm)
     weight<- try(solve( gammaNm))
     dtd<-try(solve(t(vdsig) %*% weight %*%vdsig))
     if( class( dtd)=="try-error"){
       break
     }
     r.Var<-dtd%*%t(vdsig)%*%weight%*%gammaadf%*%weight%*%vdsig%*%dtd
     r.SE <- sqrt(diag(r.Var)/(n-1))
     r.Varu<-dtd%*%t(vdsig)%*%weight%*%gammaadfu%*%weight%*%vdsig%*%dtd
     r.SEu <- sqrt(diag(r.Varu)/(n-1))
     Var<-dtd
     SE <- sqrt(diag(Var)/(n-1))
     
     U = lavTech(fites,"UfromUGamma")[[1]]
    ugamma <-  as.matrix(U %*% gammaadf)
    ugammau <-  as.matrix(U %*% gammaadfu)
    rank<-qr(ugamma)$rank
    ranku<-qr(ugammau)$rank
    
    ###### calculate a2    
    U5=sqrtm(U)
    wi=apply(si, 2, function(si)  U5%*%si)
    wi=Re(wi)
    wbar=rowMeans(wi)
    yi=wi-wbar
    H= yi%*%t(yi)
    
    Dmatrix<-function(i){
      t(yi[,i])%*%yi[,i]
    } 
    D<-lapply(1:n, Dmatrix)
    D<-unlist(D)
    a2c<-1/(n*(n-1)*(n-2)*(n-3))*((n-2)*(n-1)*sum(diag(H%*%H))-
                                    n*(n-1)*sum(D^2)+sum(diag(H))^2 )
    
    ######
    Tstats<-fites@ test[[1]]$stat
    Fstats= Tstats/(n-1)
    
    #noncorrect
    pvalue1<-pchisq(Tstats,df,lower.tail =FALSE)
    
    #### unbiased gamma   
    #recaled correction
    c1 <- sum(diag(ugammau ))/df
    rTstats1 <- Tstats/c1
    pvalue2<-pchisq(rTstats1,df,lower.tail =FALSE)
    
    #rank deficient correction
    c2<-sum(diag(ugammau ))/ranku
    rTstats2 <- Tstats/c2
    pvalue3<-pchisq(rTstats2,ranku,lower.tail =FALSE) # should be ranku
    
    #mean-variance correction
    ugamma2u <- ugammau %*% ugammau
    c3 <-sum(diag( ugamma2u))/sum(diag(ugammau )) 
    df2<-(sum(diag(ugammau )))^2/ sum(diag( ugamma2u))
    rTstats3 <- Tstats/c3
    pvalue4<-pchisq(rTstats3,df2,lower.tail =FALSE)
    
    ###### mean-variance correction a2    
    a2c.u<-1/(n*(n-1)*(n-2)*(n-3))*((n-2)*(n-1)^3*sum(diag(U%*%gammaadfu%*%U%*%gammaadfu))-
                                      n*(n-1)*sum(D^2)+sum(diag((n-1)*U%*%gammaadfu))^2 )
    c13 <-a2c.u/sum(diag(ugammau )) 
    df2.2.u<-(sum(diag(ugammau )))^2/ a2c.u
    rTstats13 <- Tstats/c13
    pvalue5<-pchisq(rTstats13,df2.2.u,lower.tail =FALSE)
    
    ms<-Mardia.tran1(xy)$skew
    mk<-Mardia.tran1(xy)$kurt
    
    ######## biased gamma
    #recaled correction
    c1 <- sum(diag(ugamma ))/df
    rTstats1 <- Tstats/c1
    pvalue6<-pchisq(rTstats1,df,lower.tail =FALSE)
    
    #rank deficient correction; rank
    c2<-sum(diag(ugamma ))/rank
    rTstats2 <- Tstats/c2
    pvalue7<-pchisq(rTstats2,rank,lower.tail =FALSE)
    
    #mean-variance correction
    ugamma2 <- ugamma %*% ugamma
    c3 <-sum(diag( ugamma2))/sum(diag(ugamma )) 
    df2<-(sum(diag(ugamma )))^2/ sum(diag( ugamma2))
    rTstats3 <- Tstats/c3
    pvalue8<-pchisq(rTstats3,df2,lower.tail =FALSE)
    
    #mean-variance correction + correct a2
    c4 <-a2c/sum(diag(ugamma )) 
    df2.2<-(sum(diag(ugamma )))^2/ a2c
    rTstats4 <- Tstats/c4
    pvalue9<-pchisq(rTstats4,df2.2,lower.tail =FALSE)
    
    #cor2; df
    c5<-0.5*(c1+c2)
    rTstats5 <- Tstats/c5
    pvalue10<-pchisq(rTstats5,df,lower.tail =FALSE)
    
    #cor3; df
    c6<-0.5*(1/c1+1/c2)
    rTstats6 <- Tstats*c6
    pvalue11<-pchisq(rTstats6,df,lower.tail =FALSE)
    
    #cor4
    pvalue12<-(pvalue6+pvalue8)/2
    
    #cor5
    pvalue13<-(pvalue6+pvalue9)/2
    
 