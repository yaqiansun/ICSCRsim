# Shared fitting implementation, 0914. All times are in days.
# No dependency beyond the recommended R packages survival and splines.
VERSION0914 <- "0914-scheme3-v1"
CORE_DIR0914 <- dirname(normalizePath(sys.frame(1)$ofile))
source(file.path(CORE_DIR0914, "gradient.R"), local=TRUE)
source(file.path(CORE_DIR0914, "visits.R"), local=TRUE)
COV0914 <- c("age_at_scan", "cs_sex", "outpatient_inpatient")

settings0914 <- function() list(
  version = VERSION0914, lambda = 1e-6, n.grid = 1000L, quadrature = 64L,
  quad.check.tol = 0.001, grid.check.tol = 0.01,
  nl.rel = 1e-10, nl.max = 3000L, nm.rel = 1e-9, nm.max = 15000L,
  restart.tol = 0.001, restart.max = 8L, gradient.tol = 0.1, score.tol = 0.05, curvature.tol = 1e-8,
  hessian.step = 0.001, se.check.tol = 0.1, theta.starts = c(0.01, 0.25),
  seed.offset = 1000L, visit.seed.offset = 914000L, n = 5000L, n.rep = 1000L,
  cens = c(17000, 20000), visit.mean = 7, frailty = TRUE)

ladder0914 <- list(
  M336 = list(degree = c(1L,1L,3L), internal = c(1L,1L,2L)),
  M436 = list(degree = c(2L,1L,3L), internal = c(1L,1L,2L)),
  M446 = list(degree = c(2L,2L,3L), internal = c(1L,1L,2L)),
  M456 = list(degree = c(2L,2L,3L), internal = c(1L,2L,2L)),
  M466 = list(degree = c(2L,3L,3L), internal = c(1L,2L,2L)),
  M467 = list(degree = c(2L,3L,3L), internal = c(1L,2L,3L)))
scenarios0914 <- list(
  knot02 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,0L,2L), beta1age=.05),
  knot12 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,1L,2L), beta1age=.05),
  knot22 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.05),
  knot23 = list(theta=.25, degree=c(2L,3L,3L), internal=c(1L,2L,2L), beta1age=.05),
  theta0001 = list(theta=.0001, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.05),
  theta05 = list(theta=.5, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.05),
  beta1age008 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.08),
  beta1age015 = list(theta=.25, degree=c(2L,3L,3L), internal=c(1L,2L,2L), beta1age=.15),
  beta1age03 = list(theta=.25, degree=c(2L,3L,3L), internal=c(1L,2L,2L), beta1age=.30))

hash0914 <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, version=2)
  unname(tools::md5sum(f))
}
atomic0914 <- function(x, file) {
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(file, ".tmp-", Sys.getpid()); saveRDS(x, tmp, version=2)
  if (file.exists(file)) stop("Refusing to overwrite completed result: ", file)
  if (!file.rename(tmp, file)) stop("Atomic save failed; retained ", tmp)
  invisible(file)
}
manifest0914 <- function(files, settings, context) {
  files<-unique(c(files,file.path(CORE_DIR0914,"gradient.R"),file.path(CORE_DIR0914,"visits.R")))
  list(version=VERSION0914, code=unname(tools::md5sum(files)),
    code.names=basename(files), settings=settings, context=context,
    R=paste(R.version$major,R.version$minor,sep="."),
    survival=as.character(utils::packageVersion("survival")))
}
check.manifest0914 <- function(old, expected) {
  if (!identical(old, expected)) stop("Code, input, R/package version, or settings differ. Use a NEW output directory; never resume mixed versions.")
  invisible(TRUE)
}
validate.data0914 <- function(d, real=FALSE) {
  req <- c("L1","R1","d1","Y2","d2",COV0914)
  if (!all(req %in% names(d))) stop("Missing input columns")
  if (!all(vapply(d[req],is.numeric,logical(1)))) stop("Model columns must be numeric")
  if (anyNA(d[req]) || any(!is.finite(as.matrix(d[setdiff(req,"R1")])))) stop("Non-finite model data")
  if (any(!d$d1 %in% 0:1 | !d$d2 %in% 0:1) || any(d$Y2<=0 | d$L1<0)) stop("Invalid indicators or times")
  a <- d$d1==1
  if (!any(a)) stop("No observed non-terminal events")
  if (any(!is.finite(d$R1[a]) | d$L1[a]>=d$R1[a] | d$R1[a]>d$Y2[a])) stop("Invalid observed event interval")
  if (any(d$R1[!a]!=Inf) || any(d$L1[!a]>d$Y2[!a])) stop("Invalid unobserved-event interval")
  if (real && any(d$L1[!a]!=d$Y2[!a])) stop("Real-data L1=Y2 convention violated")
  if (any(!d$cs_sex %in% 0:1 | !d$outpatient_inpatient %in% 0:1)) stop("Binary covariates must be coded 0/1")
  invisible(TRUE)
}

# Event generation is unchanged; independent full visit plans replace the old shortcut.
# Undetected illness is deliberately retained, measured, and reported separately.
generate0914 <- function(rep, scenario, set=settings0914()) {
  sc <- scenarios0914[[scenario]]; if(is.null(sc)) stop("Unknown scenario")
  set.seed(set$seed.offset+rep)
  n <- set$n
  X <- cbind(age_at_scan=round(rnorm(n,55,12),1), cs_sex=rbinom(n,1,.70),
             outpatient_inpatient=rbinom(n,1,.40))
  beta <- list(c(sc$beta1age,.15,-.74), c(.05,-.22,-.45), c(.03,.19,-.21))
  k <- exp(c(-9.23,-11.99,-7.81)); a <- exp(c(-.56,-.12,-.42))
  theta <- if(set$frailty) sc$theta else 0
  gamma <- if(theta>1e-6) rgamma(n,1/theta,1/theta) else rep(1,n)
  lp <- lapply(beta,function(b) as.vector(X%*%b))
  T1 <- rweibull(n,a[1],(k[1]*exp(lp[[1]])*gamma)^(-1/a[1]))
  T2 <- rweibull(n,a[2],(k[2]*exp(lp[[2]])*gamma)^(-1/a[2]))
  ii <- which(T1<T2)
  if(length(ii)) T2[ii] <- T1[ii]+rweibull(length(ii),a[3],(k[3]*exp(lp[[3]][ii])*gamma[ii])^(-1/a[3]))
  C <- runif(n,set$cens[1],set$cens[2]); Y2 <- pmin(T2,C); d2 <- as.integer(T2<=C)
  # Visit schedules receive only administrative endpoints, never event times.
  # Separate fixed seed preserves event draws and pairs knot scenarios.
  set.seed(set$visit.seed.offset+rep)
  plans <- lapply(C, visit.plan0914, mean.gap=set$visit.mean)
  obs <- t(vapply(seq_len(n), function(i)
    observe.visits0914(plans[[i]], T1[i], Y2[i]), numeric(5)))
  d1 <- as.integer(obs[,"d1"])
  d <- data.frame(id=seq_len(n), L1=obs[,"L1"], R1=obs[,"R1"],
    T1=T1, d1=d1, T2=T2, Y2=Y2, C=C, d2=d2,
    last.negative.visit=obs[,"last.negative.visit"],
    completed.visits=obs[,"completed.visits"],
    truly_ill=as.integer(T1<Y2), missed_ill=as.integer(T1<Y2 & d1==0), X)
  # Preserve latent T1/T2 for diagnostics. The model does not read these fields.
  # For Case 3, Inf expresses the agreed working interpretation, not a new truth.
  d$T1.working <- ifelse(d1==0 & d2==1, Inf, NA_real_)
  stopifnot(all(d$L1[d1==0]==d$Y2[d1==0]))
  d$cs_sex <- as.integer(d$cs_sex)
  d$outpatient_inpatient <- as.integer(d$outpatient_inpatient)
  validate.data0914(d); d
}
truth0914 <- function(scenario, frailty=TRUE) {
  sc <- scenarios0914[[scenario]]
  b <- c(sc$beta1age,.15,-.74,.05,-.22,-.45,.03,.19,-.21)
  names(b) <- unlist(lapply(1:3,function(k)paste0("beta",k,"_",COV0914)))
  if(frailty) c(b,log_theta=log(sc$theta)) else b
}

# Gauss-Legendre nodes on (0,1); no evaluation at singular endpoints.
.gl0914 <- new.env(parent=emptyenv())
gauss0914 <- function(n) {
  key<-as.character(n); if(exists(key,.gl0914,inherits=FALSE))return(get(key,.gl0914))
  v<-(1:(n-1))/sqrt(4*(1:(n-1))^2-1); J<-matrix(0,n,n)
  J[cbind(1:(n-1),2:n)]<-v; J[cbind(2:n,1:(n-1))]<-v
  z<-eigen(J,symmetric=TRUE); o<-order(z$values)
  ans<-list(x=(z$values[o]+1)/2,w=z$vectors[1,o]^2)
  assign(key,ans,.gl0914); ans
}
rowlogsum0914 <- function(x) {
  m<-apply(x,1,max); ans<-rep(-Inf,nrow(x)); ok<-is.finite(m)
  ans[ok]<-m[ok]+log(rowSums(exp(x[ok,,drop=FALSE]-m[ok]))); ans
}
moment0914 <- function(K,theta,m) {
  if(theta==0)return(-K)
  x<-theta*K; g<-log1p(x); ratio<-rep(1,length(x)); j<-x!=0
  ratio[j]<-g[j]/x[j]
  ans<--K*ratio-m*g+(m==2)*log1p(theta)
  dim(ans)<-dim(K); ans
}

basis0914 <- function(d,config,set=settings0914()) {
  ill<-d$d1==1; proxy<-ifelse(ill,(d$L1+ifelse(is.finite(d$R1),d$R1,d$Y2))/2,d$Y2)
  maxima<-c(max(d$Y2),max(d$Y2),max((d$Y2-d$L1)[ill]))
  samples<-list(proxy[ill],d$Y2[d$d2==1],(d$Y2-proxy)[ill])
  out<-lapply(1:3,function(k){
    ne<-config$internal[k]; deg<-config$degree[k]
    if(ne>0 && !length(samples[[k]]))stop("No events available for knot placement")
    knots<-if(ne)as.numeric(quantile(samples[[k]],ppoints(ne))) else numeric()
    if(anyDuplicated(knots)||any(knots<=0|knots>=maxima[k]))stop("Non-distinct or boundary interior knots")
    grid<-sort(unique(c(0,seq(0,maxima[k],length.out=set$n.grid),
      expm1(seq(0,log1p(maxima[k]),length.out=201)),knots)))
    full<-c(rep(0,deg+1),knots,rep(maxima[k],deg+1))
    B<-splines::splineDesign(full,grid,ord=deg+1,outer.ok=TRUE)
    # Projection avoids zero, while likelihood integration explicitly starts at zero.
    init.grid<-seq(min(1,maxima[k]/2),maxima[k],length.out=500)
    Bi<-splines::splineDesign(full,init.grid,ord=deg+1,outer.ok=TRUE)
    list(grid=grid,B=B,init.grid=init.grid,B.init=Bi,knots=knots,boundary=c(0,maxima[k]),degree=deg)
  }); names(out)<-paste0("t",1:3); out
}
# Positive piecewise-linear hazard and its EXACT piecewise-quadratic integral.
# This fixes both the missing first day and the mismatch between interpolants.
index0914 <- function(grid,x) {
  if(any(x<0|x>tail(grid,1)+1e-8))stop("Evaluation outside baseline support")
  j<-pmin(findInterval(x,grid),length(grid)-1L); j<-pmax(j,1L)
  list(j=j,dx=x-grid[j],dim=dim(x))
}
curves0914 <- function(b,phi) {
  lh<-as.vector(b$B%*%phi)
  if(any(!is.finite(lh))||any(lh>700|lh< -700))return(NULL)
  h<-exp(lh); dt<-diff(b$grid); slope<-diff(h)/dt
  H<-c(0,cumsum(dt*(head(h,-1)+tail(h,-1))/2))
  if(any(!is.finite(H)))return(NULL)
  list(h=h,H=H,slope=slope)
}
evalcurve0914 <- function(curve,idx) {
  j<-idx$j; dx<-idx$dx
  h<-curve$h[j]+curve$slope[j]*dx
  H<-curve$H[j]+curve$h[j]*dx+curve$slope[j]*dx^2/2
  dim(h)<-dim(H)<-idx$dim; list(logh=log(h),H=H)
}

nodes0914 <- function(d,nq,alpha=c(1,1,1),weibull=FALSE) {
  ii<-which(d$d1==1); lo<-d$L1[ii]; hi<-pmin(d$R1[ii],d$Y2[ii]); mid<-(lo+hi)/2
  gl<-gauss0914(nq); lv<-log(gl$x); w<-log(gl$w)
  pL<-ifelse(weibull & lo==0,1/alpha[1],1)
  pR<-ifelse(weibull & hi==d$Y2[ii] & d$d2[ii]==1,1/alpha[3],1)
  left.dx<-exp(outer(pL,lv)+log(mid-lo)); right.dx<-exp(outer(pR,lv)+log(hi-mid))
  uL<-left.dx+lo; sL<-d$Y2[ii]-uL
  sR<-right.dx+(d$Y2[ii]-hi); uR<-hi-right.dx
  loguL<-log(uL); loguL[lo==0,]<-(outer(pL,lv)+log(mid-lo))[lo==0,,drop=FALSE]
  logsR<-log(sR); edge<-hi==d$Y2[ii]
  logsR[edge,]<-(outer(pR,lv)+log(hi-mid))[edge,,drop=FALSE]
  weights<-function(p,len)outer(p-1,lv)+log(p)+log(len)+matrix(w,length(ii),nq,byrow=TRUE)
  list(u=cbind(uL,uR),s=cbind(sL,sR),logu=cbind(loguL,log(uR)),
       logs=cbind(log(sL),logsR),logw=cbind(weights(pL,mid-lo),weights(pR,hi-mid)),ii=ii)
}

objective0914 <- function(d,model="Weibull",basis=NULL,set=settings0914(),nq=set$quadrature) {
  validate.data0914(d); X<-as.matrix(d[COV0914]); p<-ncol(X); ill<-which(d$d1==1); healthy<-which(d$d1==0)
  iswb<-identical(model,"Weibull"); nB<-if(iswb)NULL else vapply(basis,function(b)ncol(b$B),integer(1))
  nb<-if(iswb)6L else sum(nB); n0<-nb+as.integer(set$frailty)
  if(!iswb){
    nd<-nodes0914(d,nq); ixu<-lapply(basis[1:2],function(b)index0914(b$grid,nd$u))
    ixs<-index0914(basis[[3]]$grid,nd$s)
    ixy<-lapply(basis[1:2],function(b)index0914(b$grid,d$Y2[healthy]))
  }
  function(par,details=FALSE) {
    fail<-function(reason)if(details)list(ok=FALSE,reason=reason,loglik=NA_real_,nll=Inf) else Inf
    if(length(par)!=n0+3*p||any(!is.finite(par)))return(fail("invalid parameters"))
    theta<-if(set$frailty)exp(par[nb+1L]) else 0
    if(!is.finite(theta))return(fail("frailty overflow"))
    eta<-lapply(0:2,function(k)as.vector(X%*%par[n0+k*p+seq_len(p)]))
    if(any(!is.finite(unlist(eta)))||any(abs(unlist(eta))>700))return(fail("linear predictor overflow"))
    E<-lapply(eta,exp); ll<-numeric(nrow(d))
    if(iswb){
      lk<-par[c(1,3,5)]; la<-par[c(2,4,6)]; a<-exp(la)
      if(any(!is.finite(a))||any(a==0))return(fail("invalid Weibull shape"))
      nd<-nodes0914(d,nq,a,TRUE)
      H<-function(k,lt)exp(lk[k]+a[k]*lt)
      logh<-function(k,lt)la[k]+lk[k]+(a[k]-1)*lt
      Khealthy<-H(1,log(d$Y2[healthy]))*E[[1]][healthy]+H(2,log(d$Y2[healthy]))*E[[2]][healthy]
      logh2<-logh(2,log(d$Y2[healthy]))
      K<-H(1,nd$logu)*E[[1]][ill]+H(2,nd$logu)*E[[2]][ill]+H(3,nd$logs)*E[[3]][ill]
      lh1<-logh(1,nd$logu); lh3<-logh(3,nd$logs)
    }else{
      ends<-cumsum(nB); begins<-c(0,head(ends,-1))
      curves<-lapply(1:3,function(k)curves0914(basis[[k]],par[begins[k]+seq_len(nB[k])]))
      if(any(vapply(curves,is.null,logical(1))))return(fail("baseline overflow"))
      yy<-lapply(1:2,function(k)evalcurve0914(curves[[k]],ixy[[k]]))
      uu<-lapply(1:2,function(k)evalcurve0914(curves[[k]],ixu[[k]])); ss<-evalcurve0914(curves[[3]],ixs)
      Khealthy<-yy[[1]]$H*E[[1]][healthy]+yy[[2]]$H*E[[2]][healthy];logh2<-yy[[2]]$logh
      K<-uu[[1]]$H*E[[1]][ill]+uu[[2]]$H*E[[2]][ill]+ss$H*E[[3]][ill]
      lh1<-uu[[1]]$logh;lh3<-ss$logh
    }
    if(any(!is.finite(K))||any(!is.finite(Khealthy)))return(fail("cumulative hazard overflow"))
    dh<-d$d2[healthy]; di<-d$d2[ill]
    ll[healthy]<-dh*(logh2+eta[[2]][healthy])+moment0914(Khealthy,theta,dh)
    li<-lh1+eta[[1]][ill]+moment0914(K,theta,1+di)
    # Avoid 0 * Inf when a non-terminal-only integrand has an endpoint limit.
    death<-which(di==1); if(length(death))li[death,]<-li[death,,drop=FALSE]+lh3[death,,drop=FALSE]+eta[[3]][ill][death]
    ll[ill]<-rowlogsum0914(li+nd$logw)
    if(any(!is.finite(ll)))return(fail("non-finite log contribution"))
    penalty<-set$lambda*sum(par^2); nll<--sum(ll)+penalty
    if(details)list(ok=TRUE,loglik=sum(ll),penalty=penalty,penalized.loglik=-nll,nll=nll,
      min.log.contribution=min(ll),floor.count=0L,clamp.count=0L,contributions=ll) else nll
  }
}

weibull.start0914 <- function(d,set=settings0914()) {
  mid<-ifelse(d$d1==1,(d$L1+ifelse(is.finite(d$R1),d$R1,d$Y2))/2,d$Y2)
  X<-as.matrix(d[COV0914]); ill<-d$d1==1
  fits<-list(survival::survreg(survival::Surv(mid,d$d1)~X,dist="weibull"),
    survival::survreg(survival::Surv(d$Y2,d$d2)~X,dist="weibull"),
    survival::survreg(survival::Surv(pmax(d$Y2[ill]-mid[ill],.Machine$double.eps),d$d2[ill])~X[ill,,drop=FALSE],dist="weibull"))
  aa<-vapply(fits,function(f)1/f$scale,numeric(1))
  base<-unlist(lapply(1:3,function(k)c(-aa[k]*coef(fits[[k]])[1],log(aa[k]))),use.names=FALSE)
  names(base)<-as.vector(rbind(paste0("log_kappa",1:3),paste0("log_alpha",1:3)))
  beta<-unlist(lapply(1:3,function(k)unname(-aa[k]*coef(fits[[k]])[-1])))
  names(beta)<-unlist(lapply(1:3,function(k)paste0("beta",k,"_",COV0914)))
  c(base,if(set$frailty)c(log_theta=log(set$theta.starts[1])),beta)
}
bs.start0914 <- function(w,basis,set=settings0914()) {
  phi<-unlist(lapply(1:3,function(k){
    a<-exp(w[paste0("log_alpha",k)]); b<-basis[[k]]
    lh<-w[paste0("log_alpha",k)]+w[paste0("log_kappa",k)]+(a-1)*log(b$init.grid)
    v<-as.vector(qr.solve(b$B.init,lh));setNames(v,paste0("phi",k,"_",seq_along(v)))
  }))
  c(phi,if(set$frailty)w["log_theta"],w[grep("^beta",names(w))])
}
gradient0914 <- function(q,fn,h=1e-4) vapply(seq_along(q),function(j){
  a<-b<-q;a[j]<-a[j]+h;b[j]<-b[j]-h;(fn(a)-fn(b))/(2*h)
},numeric(1))
inference0914 <- function(par,obj,set=settings0914(),gradient=NULL) {
  scale<-rep(1,length(par));scale[grepl("_age_at_scan$",names(par))]<-.02
  q<-par/scale;fn<-function(x)obj(x*scale)
  gr<-if(is.null(gradient))NULL else function(x)gradient(x*scale)*scale
  getH<-function(step)tryCatch(optimHess(q,fn,gr=gr,control=list(ndeps=rep(step,length(q)))),error=function(e)NULL)
  H1<-getH(set$hessian.step);H2<-getH(set$hessian.step/2)
  inspect<-function(H){
    if(is.null(H)||any(!is.finite(H)))return(NULL)
    H<-(H+t(H))/2;ch<-tryCatch(chol(H),error=function(e)NULL)
    if(is.null(ch))return(NULL)
    V<-chol2inv(ch);if(any(!is.finite(V))||any(diag(V)<=0))return(NULL)
    list(H=H,V=V,se=sqrt(diag(V))*scale,eigen=eigen(H,symmetric=TRUE,only.values=TRUE)$values)
  }
  a<-inspect(H1);b<-inspect(H2);pd<-!is.null(a)&&!is.null(b)
  curvature<-function(H){
    if(is.null(H)||any(!is.finite(H)))return(FALSE)
    ev<-eigen((H+t(H))/2,symmetric=TRUE,only.values=TRUE)$values
    min(ev)>= -set$curvature.tol*max(1,max(abs(ev)))
  }
  curvature.ok<-curvature(H1)&&curvature(H2)
  grad<-if(is.null(gr))gradient0914(q,fn) else gr(q)
  rel<-if(pd)max(abs(a$se-b$se)/pmax(a$se,b$se)) else Inf
  gstat<-if(pd&&all(is.finite(grad)))max(abs(grad)*sqrt(diag(b$V))) else Inf
  ok<-pd && rel<=set$se.check.tol && gstat<=set$gradient.tol
  se<-if(pd)setNames(b$se,names(par)) else setNames(rep(NA_real_,length(par)),names(par))
  V<-if(pd)b$V*outer(scale,scale) else NULL
  H<-if(!is.null(H2))H2/outer(scale,scale) else NULL
  list(ok=ok,pd=pd,curvature.ok=curvature.ok,H=H,H.scaled=H2,H.scaled.other=H1,V=V,se=se,
    gradient=grad,standardized.gradient=gstat,se.relative.change=rel,
    scaled.condition=if(pd)max(b$eigen)/min(b$eigen) else Inf,
    theta.wald.weak=if("log_theta"%in%names(se))!is.finite(se["log_theta"])||se["log_theta"]>2 else FALSE)
}

fit0914 <- function(d,model="Weibull",config=NULL,start=NULL,set=settings0914()) {
  t0<-proc.time()[3];iswb<-identical(model,"Weibull")
  basis<-if(iswb)NULL else basis0914(d,config,set)
  if(is.null(start))start<-weibull.start0914(d,set)
  if(!iswb && !any(grepl("^phi",names(start))))start<-bs.start0914(start,basis,set)
  obj<-objective0914(d,model,basis,set)
  gradient<-if(iswb)gradient.wb0914(d,set) else gradient.bs0914(d,basis,set)
  if(!is.finite(obj(start)))stop("Invalid starting values; no sentinel optimization is attempted")
  scale<-rep(1,length(start));scale[grepl("_age_at_scan$",names(start))]<-.02
  fn<-function(q)obj(setNames(q*scale,names(start)))
  gr<-function(q)gradient(setNames(q*scale,names(start)))*scale
  run<-function(s){
    q<-s/scale;trace<-numeric(); code<-1L;msg<-"";capped<-FALSE
    cat("  start log(theta)=",if(set$frailty)s["log_theta"] else NA_real_," objective=",fn(q),"\n")
    # Analytic score obtains an efficient initial solution before simplex checks.
    initial<-nlminb(q,fn,gradient=gr,control=list(eval.max=set$nl.max,iter.max=set$nl.max,rel.tol=set$nl.rel,x.tol=set$nl.rel))
    if(is.finite(initial$objective)&&initial$objective<=fn(q)){
      q<-initial$par;code<-initial$convergence;msg<-initial$message
    }
    if(!iswb){
      previous<-fn(q)
      for(pass in seq_len(set$restart.max)){
        z<-optim(q,fn,method="Nelder-Mead",control=list(maxit=set$nm.max,reltol=set$nm.rel))
        if(!is.finite(z$value)||z$value>previous+1e-8)break
        gain<-previous-z$value;q<-z$par;previous<-z$value;trace<-c(trace,-z$value)
        cat("    simplex pass",pass,"gain",gain,"\n")
        if(pass>=2 && gain<set$restart.tol)break
        if(pass==set$restart.max)capped<-TRUE
      }
    }
    z<-nlminb(q,fn,gradient=gr,control=list(eval.max=set$nl.max,iter.max=set$nl.max,rel.tol=set$nl.rel,x.tol=set$nl.rel))
    if(is.finite(z$objective)&&z$objective<=fn(q)){q<-z$par;code<-z$convergence;msg<-z$message}
    cat("  end objective=",fn(q)," code=",code,"\n")
    list(par=q*scale,nll=fn(q),code=code,message=msg,nm.trace=trace,nm.capped=capped)
  }
  starts<-list(start)
  if(set$frailty)for(th in set$theta.starts){s<-start;s["log_theta"]<-log(th);if(!isTRUE(all.equal(s,start)))starts[[length(starts)+1L]]<-s}
  trials<-lapply(starts,function(s)tryCatch(run(s),error=function(e)list(nll=Inf,code=99L,message=conditionMessage(e))))
  best<-trials[[which.min(vapply(trials,function(z)z$nll,numeric(1)))]]
  if(!is.finite(best$nll))stop("All optimization starts failed")
  est<-setNames(best$par,names(start));details<-obj(est,TRUE)
  check.obj<-objective0914(d,model,basis,set,nq=2L*set$quadrature)
  check.gradient<-if(iswb)gradient.wb0914(d,set,nq=2L*set$quadrature) else gradient.bs0914(d,basis,set,nq=2L*set$quadrature)
  quad.diff<-abs(check.obj(est)-obj(est))
  grid.diff<-0
  if(!iswb){s2<-set;s2$n.grid<-2L*set$n.grid;fine<-basis0914(d,config,s2)
    grid.diff<-abs(objective0914(d,model,fine,set,nq=2L*set$quadrature)(est)-check.obj(est))}
  inf<-inference0914(est,check.obj,set,check.gradient)
  point.ok<-best$code==0L && details$ok && is.finite(quad.diff)&&quad.diff<=set$quad.check.tol &&
    is.finite(grid.diff)&&grid.diff<=set$grid.check.tol && inf$curvature.ok && all(is.finite(inf$gradient))&&max(abs(inf$gradient))<=set$score.tol
  eligible<-point.ok && inf$ok
  np<-length(est)
  list(version=VERSION0914,model=model,config=config,est=est,start=start,
    loglik=details$loglik,penalized.loglik=details$penalized.loglik,penalty=details$penalty,
    npar=np,aic=2*np-2*details$loglik,code=best$code,message=best$message,
    point.ok=point.ok,inference.ok=eligible,inference=inf,
    quadrature.difference=quad.diff,grid.difference=grid.diff,
    nm.trace=best$nm.trace,nm.capped=best$nm.capped,
    starts=lapply(trials,function(z)z[c("nll","code","message","nm.trace","nm.capped")]),
    basis=basis,settings=set,minutes=(proc.time()[3]-t0)/60,
    final.diagnostics=details[setdiff(names(details),"contributions")])
}

estimate.table0914 <- function(f) {
  p<-f$est;s<-f$inference$se
  beta<-grepl("^beta",names(p));frailty<-names(p)=="log_theta"
  data.frame(Parameter=names(p),Est=unname(p),SE=unname(s),z=ifelse(frailty,NA_real_,unname(p/s)),p=ifelse(frailty,NA_real_,2*pnorm(-abs(unname(p/s)))),
    CI.lo=ifelse(frailty,NA_real_,unname(p-1.96*s)),CI.hi=ifelse(frailty,NA_real_,unname(p+1.96*s)),HR=ifelse(beta,exp(unname(p)),NA_real_),
    HR.lo=ifelse(beta,exp(unname(p-1.96*s)),NA_real_),HR.hi=ifelse(beta,exp(unname(p+1.96*s)),NA_real_),
    inference_usable=f$inference.ok & !frailty,
    role=ifelse(beta,"regression",ifelse(frailty,"frailty diagnostic; profile required","baseline")),row.names=NULL)
}
fit.row0914 <- function(f) data.frame(Model=f$model,p=f$npar,loglik=f$loglik,
  penalty=f$penalty,penalized_loglik=f$penalized.loglik,AIC=f$aic,code=f$code,
  point_ok=f$point.ok,inference_ok=f$inference.ok,hessian_pd=f$inference$pd,curvature_ok=f$inference$curvature.ok,
  standardized_gradient=f$inference$standardized.gradient,se_relative_change=f$inference$se.relative.change,
  scaled_condition=f$inference$scaled.condition,quadrature_difference=f$quadrature.difference,
  grid_difference=f$grid.difference,nm_capped=f$nm.capped,minutes=f$minutes,row.names=NULL)

summarize0914 <- function(results,truth,model) {
  fits<-lapply(results,function(z)z[[model]]); n<-length(fits)
  usable<-vapply(fits,function(f)!is.null(f)&&isTRUE(f$point.ok),logical(1))
  validci<-vapply(fits,function(f)!is.null(f)&&isTRUE(f$inference.ok),logical(1))
  out<-lapply(names(truth),function(p){
    est<-vapply(fits,function(f)if(is.null(f)||is.null(f$est[p]))NA_real_ else unname(f$est[p]),numeric(1))
    se<-vapply(fits,function(f)if(is.null(f)||is.null(f$inference$se[p]))NA_real_ else unname(f$inference$se[p]),numeric(1))
    use<-usable&is.finite(est); ci<-validci&is.finite(est)&is.finite(se)&se>0
    coverage<-abs(est-truth[p])<=1.96*se
    data.frame(Parameter=p,Truth=unname(truth[p]),n_attempted=n,n_point=sum(use),n_ci=sum(ci),
      Est=if(any(use))mean(est[use]) else NA_real_,Bias=if(any(use))mean(est[use])-truth[p] else NA_real_,
      SEe=if(sum(use)>1)sd(est[use]) else NA_real_,SEa=if(any(ci))mean(se[ci]) else NA_real_,
      SEe_CI_subset=if(sum(ci)>1)sd(est[ci]) else NA_real_,
      Coverage_conditional=if(any(ci))mean(coverage[ci]) else NA_real_,
      CI_success_fraction=if(n)sum(ci)/n else NA_real_,
      Covered_and_CI_available_fraction=if(n)sum(coverage[ci])/n else NA_real_,row.names=NULL)
  });do.call(rbind,out)
}
