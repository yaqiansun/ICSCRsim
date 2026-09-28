# Analytic gradients (scores) of the penalized negative log-likelihood, used by nlminb and for
# the Hessian in inference(). gradient.bs(): B-spline model; gradient.wb(): Weibull model.
# Both differentiate exactly the likelihood that objective() computes, on the same quadrature
# nodes and grid (including the Weibull endpoint transformations and their Jacobians).
# derivative.curve() / eval.derivative(): derivatives of the interpolated h and H with respect
# to the spline coefficients. moment.theta(): derivative of moment() with respect to log(theta).
derivative.curve <- function(b,curve) {
  dh<-b$B*curve$h;dt<-diff(b$grid)
  ds<-(dh[-1,,drop=FALSE]-dh[-nrow(dh),,drop=FALSE])/dt
  area<-(dh[-1,,drop=FALSE]+dh[-nrow(dh),,drop=FALSE])*dt/2
  dH<-rbind(rep(0,ncol(dh)),apply(area,2,cumsum))
  list(dh=dh,dH=dH,ds=ds)
}
eval.derivative <- function(dc,idx) {
  j<-as.vector(idx$j);dx<-as.vector(idx$dx)
  list(dh=dc$dh[j,,drop=FALSE]+dc$ds[j,,drop=FALSE]*dx,
       dH=dc$dH[j,,drop=FALSE]+dc$dh[j,,drop=FALSE]*dx+dc$ds[j,,drop=FALSE]*(dx^2/2))
}
moment.theta <- function(K,theta,m) {
  if(theta==0)return(K*0)
  x<-theta*K;f<-x/2-2*x^2/3+3*x^3/4
  large<-x>1e-4
  f[large]<-log1p(x[large])/x[large]-1/(1+x[large])
  K*f-m*x/(1+x)+(m==2)*theta/(1+theta)
}
gradient.bs <- function(d,basis,set=settings(),nq=set$quadrature) {
  X<-as.matrix(d[COVARIATES]);case12<-which(d$d1==1);case34<-which(d$d1==0)
  delta2_12<-d$d2[case12];delta2_34<-d$d2[case34];quad<-nodes(d,nq)
  nB<-vapply(basis,function(b)ncol(b$B),integer(1));end<-cumsum(nB);begin<-c(0,head(end,-1))
  nb<-sum(nB);n0<-nb+as.integer(set$frailty);p<-ncol(X)
  ixu<-lapply(basis[1:2],function(b)grid.index(b$grid,quad$u));ixs<-grid.index(basis[[3]]$grid,quad$s)
  ixy<-lapply(basis[1:2],function(b)grid.index(b$grid,d$Y2[case34]))
  function(par) {
    bad<-function()rep(NaN,length(par))
    theta<-if(set$frailty)exp(par[nb+1]) else 0
    eta<-lapply(0:2,function(k)as.vector(X%*%par[n0+k*p+seq_len(p)]));exp_eta<-lapply(eta,exp)
    cu<-lapply(1:3,function(k)make.curve(basis[[k]],par[begin[k]+seq_len(nB[k])]))
    if(any(vapply(cu,is.null,logical(1)))||!is.finite(theta))return(bad())
    dc<-lapply(1:3,function(k)derivative.curve(basis[[k]],cu[[k]]))
    uu<-lapply(1:2,function(k)eval.curve(cu[[k]],ixu[[k]]));uu[[3]]<-eval.curve(cu[[3]],ixs)
    du<-lapply(1:2,function(k)eval.derivative(dc[[k]],ixu[[k]]));du[[3]]<-eval.derivative(dc[[3]],ixs)
    yy<-lapply(1:2,function(k)eval.curve(cu[[k]],ixy[[k]]));dy<-lapply(1:2,function(k)eval.derivative(dc[[k]],ixy[[k]]))
    Ki<-lapply(1:3,function(k)uu[[k]]$H*exp_eta[[k]][case12]);K_u<-Reduce(`+`,Ki)
    Ky<-lapply(1:2,function(k)yy[[k]]$H*exp_eta[[k]][case34]);K2_y2<-Reduce(`+`,Ky)
    if(any(!is.finite(K_u))||any(!is.finite(K2_y2)))return(bad())
    log_integrand<-uu[[1]]$logh+eta[[1]][case12]+moment(K_u,theta,1+delta2_12)+quad$logw
    death<-which(delta2_12==1);if(length(death))log_integrand[death,]<-log_integrand[death,,drop=FALSE]+uu[[3]]$logh[death,,drop=FALSE]+eta[[3]][case12][death]
    norm<-rowlogsum(log_integrand);W<-exp(log_integrand-norm)
    if(any(!is.finite(W)))return(bad())
    mk<--(1+(1+delta2_12)*theta)/(1+theta*K_u)
    my<--(1+delta2_34*theta)/(1+theta*K2_y2)
    scores<-numeric(length(par))
    for(k in 1:3){
      g<-as.vector(crossprod(du[[k]]$dH,as.vector(W*mk*exp_eta[[k]][case12])))
      hweight<-if(k==1)W/exp(uu[[1]]$logh) else if(k==3)W*delta2_12/exp(uu[[3]]$logh) else NULL
      if(!is.null(hweight))g<-g+as.vector(crossprod(du[[k]]$dh,as.vector(hweight)))
      if(k<3){g<-g+as.vector(crossprod(dy[[k]]$dH,my*exp_eta[[k]][case34]))
        if(k==2)g<-g+as.vector(crossprod(dy[[2]]$dh,delta2_34/exp(yy[[2]]$logh)))}
      scores[begin[k]+seq_len(nB[k])]<-g
      subj<-numeric(nrow(d));subj[case12]<-rowSums(W*mk*Ki[[k]])+if(k==1)1 else if(k==3)delta2_12 else 0
      if(k<3)subj[case34]<-my*Ky[[k]]+if(k==2)delta2_34 else 0
      scores[n0+(k-1)*p+seq_len(p)]<-as.vector(crossprod(X,subj))
    }
    if(set$frailty)scores[nb+1]<-sum(moment.theta(K2_y2,theta,delta2_34))+sum(W*moment.theta(K_u,theta,1+delta2_12))
    -scores+2*set$lambda*par
  }
}

# Exact differentiation of the Weibull quadrature, including its parameter-
# dependent endpoint transformations and their Jacobians.
gradient.wb <- function(d,set=settings(),nq=set$quadrature) {
  X<-as.matrix(d[COVARIATES]);case12<-which(d$d1==1);case34<-which(d$d1==0)
  delta2_12<-d$d2[case12];delta2_34<-d$d2[case34];n0<-6+as.integer(set$frailty);p<-ncol(X)
  lv<-log(gauss.legendre(nq)$x)
  function(par){
    log_kappa<-par[c(1,3,5)];alpha<-exp(par[c(2,4,6)]);theta<-if(set$frailty)exp(par[7]) else 0
    quad<-nodes(d,nq,alpha,TRUE)
    eta<-lapply(0:2,function(k)as.vector(X%*%par[n0+k*p+seq_len(p)]));exp_eta<-lapply(eta,exp)
    Ki<-list(exp(log_kappa[1]+alpha[1]*quad$logu)*exp_eta[[1]][case12],exp(log_kappa[2]+alpha[2]*quad$logu)*exp_eta[[2]][case12],exp(log_kappa[3]+alpha[3]*quad$logs)*exp_eta[[3]][case12])
    K_u<-Reduce(`+`,Ki);ly<-log(d$Y2[case34]);Ky<-list(exp(log_kappa[1]+alpha[1]*ly)*exp_eta[[1]][case34],exp(log_kappa[2]+alpha[2]*ly)*exp_eta[[2]][case34]);K2_y2<-Reduce(`+`,Ky)
    log_lambda01_u<-log(alpha[1])+log_kappa[1]+(alpha[1]-1)*quad$logu
    log_lambda03_s<-log(alpha[3])+log_kappa[3]+(alpha[3]-1)*quad$logs
    log_integrand<-log_lambda01_u+eta[[1]][case12]+moment(K_u,theta,1+delta2_12)+quad$logw
    death<-which(delta2_12==1);if(length(death))log_integrand[death,]<-log_integrand[death,,drop=FALSE]+log_lambda03_s[death,,drop=FALSE]+eta[[3]][case12][death]
    W<-exp(log_integrand-rowlogsum(log_integrand));mk<--(1+(1+delta2_12)*theta)/(1+theta*K_u);my<--(1+delta2_34*theta)/(1+theta*K2_y2)
    scores<-numeric(length(par))
    for(k in 1:3){
      event<-if(k==1)rep(1,length(case12)) else if(k==3)delta2_12 else rep(0,length(case12))
      lt<-if(k==3)quad$logs else quad$logu
      sk<-event+mk*Ki[[k]]
      sa<-event*(1+alpha[k]*lt)+mk*Ki[[k]]*alpha[k]*lt
      if(k%in%c(1,3)){
        du<-ds<-jw<-matrix(0,length(case12),2*nq)
        if(k==1){rows<-which(d$L1[case12]==0);cols<-seq_len(nq)
          if(length(rows)){
            du[rows,cols]<-matrix(-lv/alpha[1],length(rows),nq,byrow=TRUE)
            ds[rows,cols]<- -quad$u[rows,cols,drop=FALSE]*du[rows,cols,drop=FALSE]/quad$s[rows,cols,drop=FALSE]
            jw[rows,cols]<-matrix(-1-lv/alpha[1],length(rows),nq,byrow=TRUE)
          }
        }else{rows<-which(d$R1[case12]==d$Y2[case12]&delta2_12==1);cols<-nq+seq_len(nq)
          if(length(rows)){
            ds[rows,cols]<-matrix(-lv/alpha[3],length(rows),nq,byrow=TRUE)
            du[rows,cols]<- -quad$s[rows,cols,drop=FALSE]*ds[rows,cols,drop=FALSE]/quad$u[rows,cols,drop=FALSE]
            jw[rows,cols]<-matrix(-1-lv/alpha[3],length(rows),nq,byrow=TRUE)
          }
        }
        sa<-sa+(alpha[1]-1)*du+delta2_12*(alpha[3]-1)*ds+mk*(Ki[[1]]*alpha[1]*du+Ki[[2]]*alpha[2]*du+Ki[[3]]*alpha[3]*ds)+jw
      }
      scores[2*k-1]<-sum(W*sk);scores[2*k]<-sum(W*sa)
      subj<-numeric(nrow(d));subj[case12]<-rowSums(W*sk)
      if(k<3){ey<-if(k==2)delta2_34 else rep(0,length(case34));yk<-ey+my*Ky[[k]]
        scores[2*k-1]<-scores[2*k-1]+sum(yk)
        scores[2*k]<-scores[2*k]+sum(ey*(1+alpha[k]*ly)+my*Ky[[k]]*alpha[k]*ly)
        subj[case34]<-yk}
      scores[n0+(k-1)*p+seq_len(p)]<-as.vector(crossprod(X,subj))
    }
    if(set$frailty)scores[7]<-sum(moment.theta(K2_y2,theta,delta2_34))+sum(W*moment.theta(K_u,theta,1+delta2_12))
    -scores+2*set$lambda*par
  }
}
