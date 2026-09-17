# Analytic score for the B-spline likelihood on the SAME fixed quadrature/grid.
# Avoids a full likelihood evaluation per parameter in each nlminb gradient.
# It differentiates the positive hazard interpolation and its exact integral.
derivative.curve0914 <- function(b,curve) {
  dh<-b$B*curve$h;dt<-diff(b$grid)
  ds<-(dh[-1,,drop=FALSE]-dh[-nrow(dh),,drop=FALSE])/dt
  area<-(dh[-1,,drop=FALSE]+dh[-nrow(dh),,drop=FALSE])*dt/2
  dH<-rbind(rep(0,ncol(dh)),apply(area,2,cumsum))
  list(dh=dh,dH=dH,ds=ds)
}
eval.derivative0914 <- function(dc,idx) {
  j<-as.vector(idx$j);dx<-as.vector(idx$dx)
  list(dh=dc$dh[j,,drop=FALSE]+dc$ds[j,,drop=FALSE]*dx,
       dH=dc$dH[j,,drop=FALSE]+dc$dh[j,,drop=FALSE]*dx+dc$ds[j,,drop=FALSE]*(dx^2/2))
}
moment.theta0914 <- function(K,theta,m) {
  if(theta==0)return(K*0)
  x<-theta*K;f<-x/2-2*x^2/3+3*x^3/4
  large<-x>1e-4
  f[large]<-log1p(x[large])/x[large]-1/(1+x[large])
  K*f-m*x/(1+x)+(m==2)*theta/(1+theta)
}
gradient.bs0914 <- function(d,basis,set=settings0914(),nq=set$quadrature) {
  X<-as.matrix(d[COV0914]);ill<-which(d$d1==1);healthy<-which(d$d1==0)
  di<-d$d2[ill];dh<-d$d2[healthy];nd<-nodes0914(d,nq)
  nB<-vapply(basis,function(b)ncol(b$B),integer(1));end<-cumsum(nB);begin<-c(0,head(end,-1))
  nb<-sum(nB);n0<-nb+as.integer(set$frailty);p<-ncol(X)
  ixu<-lapply(basis[1:2],function(b)index0914(b$grid,nd$u));ixs<-index0914(basis[[3]]$grid,nd$s)
  ixy<-lapply(basis[1:2],function(b)index0914(b$grid,d$Y2[healthy]))
  function(par) {
    bad<-function()rep(NaN,length(par))
    theta<-if(set$frailty)exp(par[nb+1]) else 0
    eta<-lapply(0:2,function(k)as.vector(X%*%par[n0+k*p+seq_len(p)]));E<-lapply(eta,exp)
    cu<-lapply(1:3,function(k)curves0914(basis[[k]],par[begin[k]+seq_len(nB[k])]))
    if(any(vapply(cu,is.null,logical(1)))||!is.finite(theta))return(bad())
    dc<-lapply(1:3,function(k)derivative.curve0914(basis[[k]],cu[[k]]))
    uu<-lapply(1:2,function(k)evalcurve0914(cu[[k]],ixu[[k]]));uu[[3]]<-evalcurve0914(cu[[3]],ixs)
    du<-lapply(1:2,function(k)eval.derivative0914(dc[[k]],ixu[[k]]));du[[3]]<-eval.derivative0914(dc[[3]],ixs)
    yy<-lapply(1:2,function(k)evalcurve0914(cu[[k]],ixy[[k]]));dy<-lapply(1:2,function(k)eval.derivative0914(dc[[k]],ixy[[k]]))
    Ki<-lapply(1:3,function(k)uu[[k]]$H*E[[k]][ill]);K<-Reduce(`+`,Ki)
    Ky<-lapply(1:2,function(k)yy[[k]]$H*E[[k]][healthy]);Khealthy<-Reduce(`+`,Ky)
    if(any(!is.finite(K))||any(!is.finite(Khealthy)))return(bad())
    li<-uu[[1]]$logh+eta[[1]][ill]+moment0914(K,theta,1+di)+nd$logw
    death<-which(di==1);if(length(death))li[death,]<-li[death,,drop=FALSE]+uu[[3]]$logh[death,,drop=FALSE]+eta[[3]][ill][death]
    norm<-rowlogsum0914(li);W<-exp(li-norm)
    if(any(!is.finite(W)))return(bad())
    mk<--(1+(1+di)*theta)/(1+theta*K)
    my<--(1+dh*theta)/(1+theta*Khealthy)
    scores<-numeric(length(par))
    for(k in 1:3){
      g<-as.vector(crossprod(du[[k]]$dH,as.vector(W*mk*E[[k]][ill])))
      hweight<-if(k==1)W/exp(uu[[1]]$logh) else if(k==3)W*di/exp(uu[[3]]$logh) else NULL
      if(!is.null(hweight))g<-g+as.vector(crossprod(du[[k]]$dh,as.vector(hweight)))
      if(k<3){g<-g+as.vector(crossprod(dy[[k]]$dH,my*E[[k]][healthy]))
        if(k==2)g<-g+as.vector(crossprod(dy[[2]]$dh,dh/exp(yy[[2]]$logh)))}
      scores[begin[k]+seq_len(nB[k])]<-g
      subj<-numeric(nrow(d));subj[ill]<-rowSums(W*mk*Ki[[k]])+if(k==1)1 else if(k==3)di else 0
      if(k<3)subj[healthy]<-my*Ky[[k]]+if(k==2)dh else 0
      scores[n0+(k-1)*p+seq_len(p)]<-as.vector(crossprod(X,subj))
    }
    if(set$frailty)scores[nb+1]<-sum(moment.theta0914(Khealthy,theta,dh))+sum(W*moment.theta0914(K,theta,1+di))
    -scores+2*set$lambda*par
  }
}

# Exact differentiation of the Weibull quadrature, including its parameter-
# dependent endpoint transformations and their Jacobians.
gradient.wb0914 <- function(d,set=settings0914(),nq=set$quadrature) {
  X<-as.matrix(d[COV0914]);ill<-which(d$d1==1);healthy<-which(d$d1==0)
  di<-d$d2[ill];dh<-d$d2[healthy];n0<-6+as.integer(set$frailty);p<-ncol(X)
  lv<-log(gauss0914(nq)$x)
  function(par){
    lk<-par[c(1,3,5)];a<-exp(par[c(2,4,6)]);theta<-if(set$frailty)exp(par[7]) else 0
    nd<-nodes0914(d,nq,a,TRUE)
    eta<-lapply(0:2,function(k)as.vector(X%*%par[n0+k*p+seq_len(p)]));E<-lapply(eta,exp)
    Ki<-list(exp(lk[1]+a[1]*nd$logu)*E[[1]][ill],exp(lk[2]+a[2]*nd$logu)*E[[2]][ill],exp(lk[3]+a[3]*nd$logs)*E[[3]][ill])
    K<-Reduce(`+`,Ki);ly<-log(d$Y2[healthy]);Ky<-list(exp(lk[1]+a[1]*ly)*E[[1]][healthy],exp(lk[2]+a[2]*ly)*E[[2]][healthy]);Khealthy<-Reduce(`+`,Ky)
    lh1<-log(a[1])+lk[1]+(a[1]-1)*nd$logu
    lh3<-log(a[3])+lk[3]+(a[3]-1)*nd$logs
    li<-lh1+eta[[1]][ill]+moment0914(K,theta,1+di)+nd$logw
    death<-which(di==1);if(length(death))li[death,]<-li[death,,drop=FALSE]+lh3[death,,drop=FALSE]+eta[[3]][ill][death]
    W<-exp(li-rowlogsum0914(li));mk<--(1+(1+di)*theta)/(1+theta*K);my<--(1+dh*theta)/(1+theta*Khealthy)
    scores<-numeric(length(par))
    for(k in 1:3){
      event<-if(k==1)rep(1,length(ill)) else if(k==3)di else rep(0,length(ill))
      lt<-if(k==3)nd$logs else nd$logu
      sk<-event+mk*Ki[[k]]
      sa<-event*(1+a[k]*lt)+mk*Ki[[k]]*a[k]*lt
      if(k%in%c(1,3)){
        du<-ds<-jw<-matrix(0,length(ill),2*nq)
        if(k==1){rows<-which(d$L1[ill]==0);cols<-seq_len(nq)
          if(length(rows)){
            du[rows,cols]<-matrix(-lv/a[1],length(rows),nq,byrow=TRUE)
            ds[rows,cols]<- -nd$u[rows,cols,drop=FALSE]*du[rows,cols,drop=FALSE]/nd$s[rows,cols,drop=FALSE]
            jw[rows,cols]<-matrix(-1-lv/a[1],length(rows),nq,byrow=TRUE)
          }
        }else{rows<-which(d$R1[ill]==d$Y2[ill]&di==1);cols<-nq+seq_len(nq)
          if(length(rows)){
            ds[rows,cols]<-matrix(-lv/a[3],length(rows),nq,byrow=TRUE)
            du[rows,cols]<- -nd$s[rows,cols,drop=FALSE]*ds[rows,cols,drop=FALSE]/nd$u[rows,cols,drop=FALSE]
            jw[rows,cols]<-matrix(-1-lv/a[3],length(rows),nq,byrow=TRUE)
          }
        }
        sa<-sa+(a[1]-1)*du+di*(a[3]-1)*ds+mk*(Ki[[1]]*a[1]*du+Ki[[2]]*a[2]*du+Ki[[3]]*a[3]*ds)+jw
      }
      scores[2*k-1]<-sum(W*sk);scores[2*k]<-sum(W*sa)
      subj<-numeric(nrow(d));subj[ill]<-rowSums(W*sk)
      if(k<3){ey<-if(k==2)dh else rep(0,length(healthy));yk<-ey+my*Ky[[k]]
        scores[2*k-1]<-scores[2*k-1]+sum(yk)
        scores[2*k]<-scores[2*k]+sum(ey*(1+a[k]*ly)+my*Ky[[k]]*a[k]*ly)
        subj[healthy]<-yk}
      scores[n0+(k-1)*p+seq_len(p)]<-as.vector(crossprod(X,subj))
    }
    if(set$frailty)scores[7]<-sum(moment.theta0914(Khealthy,theta,dh))+sum(W*moment.theta0914(K,theta,1+di))
    -scores+2*set$lambda*par
  }
}
