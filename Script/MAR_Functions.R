## Function to get p value 

getpval<-function(m.ci){
  est<-c(m.ci$par$B,m.ci$par$Q,m.ci$par$U,m.ci$par$R)
  lo<-c(m.ci$par.lowCI$B,m.ci$par.lowCI$Q,m.ci$par.lowCI$U,m.ci$par.lowCI$R)
  up<-c(m.ci$par.upCI$B,m.ci$par.upCI$Q,m.ci$par.upCI$U,m.ci$par.upCI$R)
  se<-c(m.ci$par.se$B,m.ci$par.se$Q,m.ci$par.se$U,m.ci$par.se$R)
  restable<-data.frame(est=est,se=se,lo=lo,up=up)
  row.names(restable)<-c(attr(m.ci$start$B,which="dimnames")[[1]],attr(m.ci$start$Q,which="dimnames")[[1]],attr(m.ci$start$U,which="dimnames")[[1]],attr(m.ci$start$R,which="dimnames")[[1]])
  
  pv.fct<-function(vec){
    z<-abs(vec[1]/vec[2])
    pv<-exp(-0.717*z - 0.416*z*z)
    pv[pv>1]<-1
    return(pv)}
  
  restable$p.val<-round(apply(restable,1,pv.fct),3)
  return(restable)
}

# Utilitary function for extracting parameter values from a MAR model fit
getparam.fct<-function(m.ci,Bzero,C=FALSE,Czero,ObsErr=F){
  nsp<-ncol(Bzero)
  nenv<-ncol(Czero)
  
  Bfit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
  Bfit[,,1]<-Bzero
  Bfit[,,1][is.na(m.ci$call$model$B==0)]<-m.ci$par$B
  Bfit[,,2]<-Bzero
  Bfit[,,2][is.na(m.ci$call$model$B==0)]<-m.ci$par.lowCI$B
  Bfit[,,3]<-Bzero
  Bfit[,,3][is.na(m.ci$call$model$B==0)]<-m.ci$par.upCI$B
  Bfit[,,4]<-Bzero
  Bfit[,,4][is.na(m.ci$call$model$B==0)]<-m.ci$par.se$B
  
  if (C==TRUE){
    Cfit<-array(NA,dim=c(nsp,nenv,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
    Cfit[,,1]<-Czero
    Cfit[,,1][is.na(m.ci$call$model$C==0)]<-m.ci$par$U
    Cfit[,,2]<-Czero
    Cfit[,,2][is.na(m.ci$call$model$C==0)]<-m.ci$par.lowCI$U
    Cfit[,,3]<-Czero
    Cfit[,,3][is.na(m.ci$call$model$C==0)]<-m.ci$par.upCI$U
    Cfit[,,4]<-Czero
    Cfit[,,4][is.na(m.ci$call$model$C==0)]<-m.ci$par.se$U
    
  }
  
  Sigmafit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
  Sigmafit[,,1]<-diag(nsp)*c(m.ci$par$Q)
  Sigmafit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$Q)
  Sigmafit[,,3]<-diag(nsp)*c(m.ci$par.upCI$Q)
  Sigmafit[,,4]<-diag(nsp)*c(m.ci$par.se$Q)


##  fittedparam
  if (C==FALSE){ Cfit<-0}
  fittedparam<-list(Bfit,Cfit,Sigmafit)
  
  if (ObsErr==T){
    ObsErrfit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
    ObsErrfit[,,1]<-diag(nsp)*c(m.ci$par$R)
    ObsErrfit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$R)
    ObsErrfit[,,3]<-diag(nsp)*c(m.ci$par.upCI$R)
    ObsErrfit[,,4]<-diag(nsp)*c(m.ci$par.se$R)
    fittedparam<-list(Bfit,Cfit,Sigmafit,ObsErrfit)}
  
  return(fittedparam)
}

###
### model selection
###

## Backward selection 
# C matrix
Cbackward.fct<-function(newX,Cstart,B,listparam){
  
  Caic<-matrix(NA,nrow=nrow(Cstart),ncol=ncol(Cstart))
  model.gen=list(Z=listparam$Z,
                 A=listparam$A,
                 R=listparam$R,
                 B=B,
                 C=Cstart,
                 c=listparam$c,
                 U=listparam$U,
                 Q=listparam$Q,
                 x0=listparam$x0,
                 V0=listparam$V0)
  m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
  AICprev<-m.try$AICc
  test<-T
  
  while(test==T){
    for (ci in 1:nrow(Cstart)){
      for (cj in 1:ncol(Cstart)){
        Ctest<-Cstart
        Ctest[ci,cj]<-0
        model.gen=list(Z=listparam$Z,
                       A=listparam$A,
                       R=listparam$R,
                       B=B,
                       C=Ctest,
                       c=listparam$c,
                       U=listparam$U,
                       Q=listparam$Q,
                       x0=listparam$x0,
                       V0=listparam$V0)
        m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
        Caic[ci,cj]<-m.try$AICc
      }
    }
    
    AICnew<-min(Caic)
    test<-AICnew<AICprev
    if (test==T){Cstart[Caic==min(Caic)]<-0
    AICprev<-AICnew}
  }
  return(Cstart)
}
# B matrix 
Bbackward.fct<-function(newX,C,Bstart,listparam,fixed.param=NULL){
  # Model without matrix C
  if (C==0){
    Baic<-matrix(NA,nrow=nrow(Bstart),ncol=ncol(Bstart))
    model.gen=list(Z=listparam$Z,
                   A=listparam$A,
                   R=listparam$R,
                   B=Bstart,
                   #C=C,
                   #c=listparam$c,
                   U=listparam$U,
                   Q=listparam$Q,
                   x0=listparam$x0,
                   V0=listparam$V0)
    m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
    AICprev<-m.try$AICc
    test<-T
    
    while(test==T){
      for (bi in 1:nrow(Bstart)){
        for (bj in 1:ncol(Bstart)){
          Btest<-Bstart
          
          if (is.null(fixed.param)==F){
            testfix<-NULL
          for (fp in 1:nrow(fixed.param)){
            testfix<-c(testfix,FALSE %in% (c(bi,bj)==fixed.param[fp,]))  
          }
           }
          
          if ((FALSE %in% testfix)==F) {
            Btest[bi,bj]<-0
            print(paste(c(bi,bj),"removed"))
          }
          model.gen=list(Z=listparam$Z,
                         A=listparam$A,
                         R=listparam$R,
                         B=Btest,
                         #C=C,
                         #c=listparam$c,
                         U=listparam$U,
                         Q=listparam$Q,
                         x0=listparam$x0,
                         V0=listparam$V0)
          m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05),silent=T)
          Baic[bi,bj]<-m.try$AICc
        }
      }
      
      AICnew<-min(Baic)
      test<-AICnew<AICprev
      if (test==T) {Bstart[Baic==min(Baic)]<-0
      AICprev<-AICnew}
    }
    return(Bstart)
  }
  reformat.fct<-function(mat){
    nr<-nrow(mat)
    nc<-ncol(mat)
    mat[mat=="0"]<-list(0)
    mat<-matrix(mat,nr,nc)
    return(mat)
    
  }
  # model with both matrix C and matrix B
  Baic<-matrix(NA,nrow=nrow(Bstart),ncol=ncol(Bstart))
  model.gen=list(Z=listparam$Z,
                 A=listparam$A,
                 R=listparam$R,
                 B=Bstart,
                 C=C,
                 c=listparam$c,
                 U=listparam$U,
                 Q=listparam$Q,
                 x0=listparam$x0,
                 V0=listparam$V0)
  m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
  AICprev<-m.try$AICc
  test<-T
  
  while(test==T){
    for (bi in 1:nrow(Bstart)){
      for (bj in 1:ncol(Bstart)){
        Btest<-Bstart
        Btest[bi,bj]<-0
        model.gen=list(Z=listparam$Z,
                       A=listparam$A,
                       R=listparam$R,
                       B=Btest,
                       C=C,
                       c=listparam$c,
                       U=listparam$U,
                       Q=listparam$Q,
                       x0=listparam$x0,
                       V0=listparam$V0)
        m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
        Baic[bi,bj]<-m.try$AICc
      }
    }
    
    AICnew<-min(Baic)
    test<-AICnew<AICprev
    if (test==T){Bstart[Baic==min(Baic)]<-0
    AICprev<-AICnew}
  }
  return(Bstart)
}
reformat.fct<-function(mat){
  nr<-nrow(mat)
  nc<-ncol(mat)
  mat[mat=="0"]<-list(0)
  mat<-matrix(mat,nr,nc)
  return(mat)
}

### Forward selection 
# C matrix 
Cforward.fct<-function(newX,Cstart,Cfull,B,listparam,deltaAIC=2){
  
  Caic<-matrix(NA,nrow=nrow(Cstart),ncol=ncol(Cstart))
  model.gen=list(Z=listparam$Z,
                 A=listparam$A,
                 R=listparam$R,
                 B=B,
                 C=Cstart,
                 c=listparam$c,
                 U=listparam$U,
                 Q=listparam$Q,
                 x0=listparam$x0,
                 V0=listparam$V0)
  m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
  AICprev<-m.try$AICc
  test<-T
  
  while(test==T){
    
    for (ci in 1:nrow(Cstart)){
      for (cj in 1:ncol(Cstart)){
        Ctest<-Cstart
        Ctest[ci,cj]<-Cfull[ci,cj][[1]]
        Ctest<-reformat.fct(Ctest)
        model.gen=list(Z=listparam$Z,
                       A=listparam$A,
                       R=listparam$R,
                       B=B,
                       C=Ctest,
                       c=listparam$c,
                       U=listparam$U,
                       Q=listparam$Q,
                       x0=listparam$x0,
                       V0=listparam$V0)
        m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
        Caic[ci,cj]<-m.try$AICc}
    }
    AICnew<-min(Caic)
    test<-AICnew+deltaAIC<AICprev
    if (test==T){Cstart[Caic==min(Caic)]<-Cfull[Caic==min(Caic)][[1]]
    AICprev<-AICnew}
  }
  Cstart<-reformat.fct(Cstart)
  return(Cstart)
}

# B matrix 
Bforward.fct<-function(newX,Bstart,Bfull,C,listparam,deltaAIC=2){
  # model without matrix C
  if (C==0) {
    Baic<-matrix(NA,nrow=nrow(Bstart),ncol=ncol(Bstart))
    model.gen=list(Z=listparam$Z,
                   A=listparam$A,
                   R=listparam$R,
                   B=Bstart,
                   #C=C,
                   #c=listparam$c,
                   U=listparam$U,
                   Q=listparam$Q,
                   x0=listparam$x0,
                   V0=listparam$V0)
    m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
    AICprev<-m.try$AICc
    test<-T
    
    while(test==T){
      
      for (bi in 1:nrow(Bstart)){
        for (bj in 1:ncol(Bstart)){
          Btest<-Bstart
          Btest[bi,bj]<-Bfull[bi,bj][[1]]
          Btest<-reformat.fct(Btest)
          model.gen=list(Z=listparam$Z,
                         A=listparam$A,
                         R=listparam$R,
                         B=Btest,
                         #C=C,
                         #c=listparam$c,
                         U=listparam$U,
                         Q=listparam$Q,
                         x0=listparam$x0,
                         V0=listparam$V0)
          m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
          Baic[bi,bj]<-m.try$AICc}
      }
      AICnew<-min(Baic)
      test<-AICnew+deltaAIC<AICprev
      if (test==T){Bstart[Baic==min(Baic)]<-Bfull[Baic==min(Baic)][[1]]
      AICprev<-AICnew}
    }
    Bstart<-reformat.fct(Bstart)
    return(Bstart)
  }
  # model with both matrix C and matrix B
  Baic<-matrix(NA,nrow=nrow(Bstart),ncol=ncol(Bstart))
  model.gen=list(Z=listparam$Z,
                 A=listparam$A,
                 R=listparam$R,
                 B=Bstart,
                 C=C,
                 c=listparam$c,
                 U=listparam$U,
                 Q=listparam$Q,
                 x0=listparam$x0,
                 V0=listparam$V0)
  m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
  AICprev<-m.try$AICc
  test<-T
  
  while(test==T){
    
    for (bi in 1:nrow(Bstart)){
      for (bj in 1:ncol(Bstart)){
        Btest<-Bstart
        Btest[bi,bj]<-Bfull[bi,bj][[1]]
        Btest<-reformat.fct(Btest)
        model.gen=list(Z=listparam$Z,
                       A=listparam$A,
                       R=listparam$R,
                       B=Btest,
                       C=C,
                       c=listparam$c,
                       U=listparam$U,
                       Q=listparam$Q,
                       x0=listparam$x0,
                       V0=listparam$V0)
        m.try<-MARSS(newX,model=model.gen,control=list(conv.test.slope.tol=0.05))
        Baic[bi,bj]<-m.try$AICc}
    }
    AICnew<-min(Baic)
    test<-AICnew+deltaAIC<AICprev
    if (test==T){Bstart[Baic==min(Baic)]<-Bfull[Baic==min(Baic)][[1]]
    AICprev<-AICnew}
  }
  Bstart<-reformat.fct(Bstart)
  return(Bstart)
}

###
#### Model simulation 
###
get.ci.sim.fct<- function(data,var.num){
  ci.inf<-apply(data[var.num,,],c(1),quantile,probs=0.025)
  ci.sup<-apply(data[var.num,,],c(1),quantile,probs=0.975)
  med<-apply(data[var.num,,],c(1),quantile,probs=0.5)
  return(cbind(ci.inf,ci.sup,med))
  }

plot.simu.fct<-function(data,var.num,main,y.lab){
  ci.inf<-apply(data[var.num,,],c(1),quantile,probs=0.025)
  ci.sup<-apply(data[var.num,,],c(1),quantile,probs=0.975)
  med<-apply(data[var.num,,],c(1),quantile,probs=0.5)
  #med<-apply(data[var.num,,],c(1),mean)
  plot(dimnames(data)[[2]], med,type='b',pch=18, xlab='Year',ylab=y.lab,main= main,ylim=range(c(ci.inf-0.2,ci.sup+0.3)),lwd=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  points(dimnames(data)[[2]],ci.inf,type="l",lty=2,lwd=1.8,col='darkgrey')
  points(dimnames(data)[[2]],ci.sup,type="l",lty=2,lwd=1.8,col='darkgrey')
}
