# Model settings, data generation, likelihood, fitting and summaries. All times are in days.
# No dependency beyond the recommended R packages survival and splines.
VERSION <- "1.0"
CORE_DIR <- local({ # folder of this file, also when source() is called from another sourced script
  src <- Filter(function(e) !is.null(e$ofile), sys.frames())
  dirname(normalizePath(src[[length(src)]]$ofile))
})
source(file.path(CORE_DIR, "gradient.R"), local=TRUE)
source(file.path(CORE_DIR, "visits.R"), local=TRUE)
COVARIATES <- c("age", "sex", "setting")

# settings(): every numerical setting in one list; it is stored with each result (manifest).
#   Data: n, n.rep, cens (C ~ Uniform), visit.mean, seeds. Likelihood: lambda (ridge), quadrature
#   (nodes per half-interval), n.grid. Optimisation: nl.*, nm.*, restart.*, theta.starts.
#   Eligibility checks: quad.check.tol, grid.check.tol, score.tol, gradient.tol, se.check.tol,
#   hessian.step, curvature.tol (see fit() and inference()).
settings <- function() list(
  version = VERSION, lambda = 1e-6, n.grid = 1000L, quadrature = 64L,
  quad.check.tol = 0.001, grid.check.tol = 0.01,
  nl.rel = 1e-10, nl.max = 3000L, nm.rel = 1e-9, nm.max = 15000L,
  restart.tol = 0.001, restart.max = 8L, gradient.tol = 0.1, score.tol = 0.05, curvature.tol = 1e-8,
  hessian.step = 0.001, se.check.tol = 0.1, theta.starts = c(0.01, 0.25),
  seed.offset = 1000L, visit.seed.offset = 914000L, n = 5000L, n.rep = 1000L,
  cens = c(17000, 20000), visit.mean = 7, frailty = TRUE)

# ladder: the six B-spline candidates for the real data; the name gives the basis counts
# (M1, M2, M3), with M_k = degree_k + internal knots_k + 1.
ladder <- list(
  M336 = list(degree = c(1L,1L,3L), internal = c(1L,1L,2L)),
  M436 = list(degree = c(2L,1L,3L), internal = c(1L,1L,2L)),
  M446 = list(degree = c(2L,2L,3L), internal = c(1L,1L,2L)),
  M456 = list(degree = c(2L,2L,3L), internal = c(1L,2L,2L)),
  M466 = list(degree = c(2L,3L,3L), internal = c(1L,2L,2L)),
  M467 = list(degree = c(2L,3L,3L), internal = c(1L,2L,3L)))
# scenarios: simulation settings. knotXY = transition 2 with X interior knots and degree Y
# (transitions 1 and 3 fixed); theta = frailty variance; beta1age = transition-1 age effect.
# The seven scenarios in the study are knot02, knot12, knot22, knot23, theta0001, theta05 and
# beta1age008; beta1age015 and beta1age03 are defined but not part of the study.
scenarios <- list(
  knot02 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,0L,2L), beta1age=.05),
  knot12 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,1L,2L), beta1age=.05),
  knot22 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.05),
  knot23 = list(theta=.25, degree=c(2L,3L,3L), internal=c(1L,2L,2L), beta1age=.05),
  theta0001 = list(theta=.0001, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.05),
  theta05 = list(theta=.5, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.05),
  beta1age008 = list(theta=.25, degree=c(2L,2L,3L), internal=c(1L,2L,2L), beta1age=.08),
  beta1age015 = list(theta=.25, degree=c(2L,3L,3L), internal=c(1L,2L,2L), beta1age=.15),
  beta1age03 = list(theta=.25, degree=c(2L,3L,3L), internal=c(1L,2L,2L), beta1age=.30))

# hash(): MD5 fingerprint of an R object (used to record each generated dataset).
# atomic.save(): write to a temporary file, then rename; never overwrite a completed result.
# make.manifest() / check.manifest(): record code checksums, settings and R/package versions,
# and refuse to mix results produced by different code or settings.
# validate.data(): input checks before fitting (columns, 0/1 indicators, L1 < R1 <= Y2, ...).
hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, version=2)
  unname(tools::md5sum(f))
}
atomic.save <- function(x, file) {
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(file, ".tmp-", Sys.getpid()); saveRDS(x, tmp, version=2)
  if (file.exists(file)) stop("Refusing to overwrite completed result: ", file)
  if (!file.rename(tmp, file)) stop("Atomic save failed; retained ", tmp)
  invisible(file)
}
make.manifest <- function(files, settings, context) {
  files<-unique(c(files,file.path(CORE_DIR,"gradient.R"),file.path(CORE_DIR,"visits.R")))
  list(version=VERSION, code=unname(tools::md5sum(files)),
    code.names=basename(files), settings=settings, context=context,
    R=paste(R.version$major,R.version$minor,sep="."),
    survival=as.character(utils::packageVersion("survival")))
}
check.manifest <- function(old, expected) {
  if (!identical(old, expected)) stop("Code, input, R/package version, or settings differ. Use a NEW output directory; never resume mixed versions.")
  invisible(TRUE)
}
validate.data <- function(d, real=FALSE) {
  req <- c("L1","R1","d1","Y2","d2",COVARIATES)
  if (!all(req %in% names(d))) stop("Missing input columns")
  if (!all(vapply(d[req],is.numeric,logical(1)))) stop("Model columns must be numeric")
  if (anyNA(d[req]) || any(!is.finite(as.matrix(d[setdiff(req,"R1")])))) stop("Non-finite model data")
  if (any(!d$d1 %in% 0:1 | !d$d2 %in% 0:1) || any(d$Y2<=0 | d$L1<0)) stop("Invalid indicators or times")
  a <- d$d1==1
  if (!any(a)) stop("No observed non-terminal events")
  if (any(!is.finite(d$R1[a]) | d$L1[a]>=d$R1[a] | d$R1[a]>d$Y2[a])) stop("Invalid observed event interval")
  if (any(d$R1[!a]!=Inf) || any(d$L1[!a]>d$Y2[!a])) stop("Invalid unobserved-event interval")
  if (real && any(d$L1[!a]!=d$Y2[!a])) stop("Real-data L1=Y2 convention violated")
  if (any(!d$sex %in% 0:1 | !d$setting %in% 0:1)) stop("Binary covariates must be coded 0/1")
  invisible(TRUE)
}

# generate(): one simulated dataset for replication `rep` of a scenario: covariates,
# shared Gamma frailty, illness-death event times, administrative censoring, then
# visit schedules and the observed intervals (see visits.R).
# Returns one row per subject; the model uses L1, R1, d1, Y2, d2 and the covariates.
generate <- function(rep, scenario, set=settings()) {
  sc <- scenarios[[scenario]]; if(is.null(sc)) stop("Unknown scenario")
  # Event seed 1000 + rep: each dataset is reproducible, and the four knot scenarios
  # share the same data (paired comparison).
  set.seed(set$seed.offset+rep)
  n <- set$n
  # Covariates: age ~ N(55, 12), sex ~ Bernoulli(0.70), setting ~ Bernoulli(0.40).
  X <- cbind(age=round(rnorm(n,55,12),1), sex=rbinom(n,1,.70),
             setting=rbinom(n,1,.40))
  # True values (Weibull fit to the heart failure data). beta_k = (age, sex, setting) for
  # transitions 0->1, 0->2, 1->2; baseline cumulative hazard kappa_k * t^alpha_k, t in days.
  # Age is not centred, so kappa_k refers to age 0.
  beta <- list(c(sc$beta1age,.15,-.74), c(.05,-.22,-.45), c(.03,.19,-.21))
  k <- exp(c(-9.23,-11.99,-7.81)); a <- exp(c(-.56,-.12,-.42))
  theta <- if(set$frailty) sc$theta else 0
  # Shared frailty gamma ~ Gamma(1/theta, 1/theta): mean 1, variance theta; it multiplies
  # all three transition hazards of the same subject.
  gamma <- if(theta>1e-6) rgamma(n,1/theta,1/theta) else rep(1,n)
  # Linear predictors X %*% beta_k for the three transitions.
  lp <- lapply(beta,function(b) as.vector(X%*%b))
  # Competing Weibull times for illness (T1) and direct death (T2). This scale gives
  # S_k(t) = exp(-gamma * exp(lp_k) * kappa_k * t^alpha_k).
  T1 <- rweibull(n,a[1],(k[1]*exp(lp[[1]])*gamma)^(-1/a[1]))
  T2 <- rweibull(n,a[2],(k[2]*exp(lp[[2]])*gamma)^(-1/a[2]))
  # Illness first (T1 < T2): death = T1 + sojourn time from transition 3 (semi-Markov).
  ii <- which(T1<T2)
  if(length(ii)) T2[ii] <- T1[ii]+rweibull(length(ii),a[3],(k[3]*exp(lp[[3]][ii])*gamma[ii])^(-1/a[3]))
  # Administrative censoring C ~ Uniform(17000, 20000) days; Y2 = min(T2, C); d2 = I(T2 <= C).
  C <- runif(n,set$cens[1],set$cens[2]); Y2 <- pmin(T2,C); d2 <- as.integer(T2<=C)
  # Visit schedules use only C. The separate visit seed 914000 + rep means the visit draws
  # do not depend on how many random numbers the event step used.
  set.seed(set$visit.seed.offset+rep)
  plans <- lapply(C, visit.plan, mean.gap=set$visit.mean)
  # obs: one row per subject with (L1, R1, d1, last.negative.visit, completed.visits).
  obs <- t(vapply(seq_len(n), function(i)
    observe.visits(plans[[i]], T1[i], Y2[i]), numeric(5)))
  d1 <- as.integer(obs[,"d1"])
  d <- data.frame(id=seq_len(n), L1=obs[,"L1"], R1=obs[,"R1"],
    T1=T1, d1=d1, T2=T2, Y2=Y2, C=C, d2=d2,
    last.negative.visit=obs[,"last.negative.visit"],
    completed.visits=obs[,"completed.visits"],
    truly_ill=as.integer(T1<Y2), missed_ill=as.integer(T1<Y2 & d1==0), X)
  # Diagnostic only (not read by the model): latent T1, T2, C; truly_ill = onset before Y2;
  # missed_ill = onset before Y2 but d1 = 0; T1.working = Inf for Case 3 (death observed,
  # no illness recorded) and NA otherwise.
  d$T1.working <- ifelse(d1==0 & d2==1, Inf, NA_real_)
  # Check: every subject without detected illness has L1 = Y2.
  stopifnot(all(d$L1[d1==0]==d$Y2[d1==0]))
  d$sex <- as.integer(d$sex)
  d$setting <- as.integer(d$setting)
  # Final input checks, then return the dataset.
  validate.data(d); d
}
# true.values(): true regression coefficients and log(theta) of a scenario (for bias and coverage).
true.values <- function(scenario, frailty=TRUE) {
  sc <- scenarios[[scenario]]
  b <- c(sc$beta1age,.15,-.74,.05,-.22,-.45,.03,.19,-.21)
  names(b) <- unlist(lapply(1:3,function(k)paste0("beta",k,"_",COVARIATES)))
  if(frailty) c(b,log_theta=log(sc$theta)) else b
}

# gauss.legendre(n): Gauss-Legendre nodes and weights on (0,1), for the integral over the
# unknown illness time in Cases 1-2. Golub-Welsch method: the nodes (roots of the Legendre
# polynomial P_n) are the eigenvalues of a symmetric tridiagonal matrix J, and the weights are
# the squared first components of its eigenvectors. Nodes lie strictly inside (0,1).
.gl.cache <- new.env(parent=emptyenv())   # results are cached: computed once per n
gauss.legendre <- function(n) {
  key<-as.character(n); if(exists(key,.gl.cache,inherits=FALSE))return(get(key,.gl.cache))
  # J: zero diagonal, off-diagonals v_k = k / sqrt(4k^2 - 1) (Legendre three-term recurrence).
  v<-(1:(n-1))/sqrt(4*(1:(n-1))^2-1); J<-matrix(0,n,n)
  J[cbind(1:(n-1),2:n)]<-v; J[cbind(2:n,1:(n-1))]<-v
  z<-eigen(J,symmetric=TRUE); o<-order(z$values)
  # Map nodes from [-1,1] to (0,1): x = (t+1)/2; the weights then sum to 1.
  ans<-list(x=(z$values[o]+1)/2,w=z$vectors[1,o]^2)
  assign(key,ans,.gl.cache); ans
}
# rowlogsum(x): log(sum(exp(x))) for each row, computed as max + log(sum(exp(x - max)))
# (log-sum-exp) so that very small likelihood terms do not underflow to zero.
rowlogsum <- function(x) {
  m<-apply(x,1,max); ans<-rep(-Inf,nrow(x)); ok<-is.finite(m)
  ans[ok]<-m[ok]+log(rowSums(exp(x[ok,,drop=FALSE]-m[ok]))); ans
}
# moment(K, theta, m) = log E[gamma^m exp(-gamma K)] for gamma ~ Gamma(1/theta, 1/theta)
#   = log(c_m) - (1/theta + m) log(1 + theta K), with c_m = 1, 1, 1 + theta for m = 0, 1, 2
# (Appendix A). m = number of observed hazards: 2 (Case 1), 1 (Cases 2, 3), 0 (Case 4).
# (1/theta) log(1 + theta K) is written as K * log1p(theta K)/(theta K), which tends to K as
# theta -> 0, so the no-frailty limit (-K) is reached smoothly.
moment <- function(K,theta,m) {
  if(theta==0)return(-K)
  x<-theta*K; g<-log1p(x); ratio<-rep(1,length(x)); j<-x!=0
  ratio[j]<-g[j]/x[j]
  ans<--K*ratio-m*g+(m==2)*log1p(theta)
  dim(ans)<-dim(K); ans
}

# make.basis(): B-spline basis for each transition (log lambda0k(t) = sum_b eta_kb B_b(t), eq. 3.5).
# Interior knots at ppoints quantiles of approximate event times: transition 1 = interval
# midpoints (delta1 = 1), transition 2 = observed death times, transition 3 = Y2 - midpoint.
# One knot = median; two knots = about the 28th and 72nd percentiles. M_k = d_k + e_k + 1.
# Grid: 1000 equally spaced points, 201 log-spaced points (dense near 0) and the knots.
# (The log-spaced end point can exceed the boundary by ~1e-11, where the basis is 0; the
#  likelihood never evaluates beyond the boundary, so this point is not used.)
make.basis <- function(d,config,set=settings()) {
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
# The cumulative hazard is the exact integral of the interpolated hazard, starting at 0.
# grid.index(): grid cell j containing each time and the offset dx = t - grid[j].
grid.index <- function(grid,x) {
  if(any(x<0|x>tail(grid,1)+1e-8))stop("Evaluation outside baseline support")
  j<-pmin(findInterval(x,grid),length(grid)-1L); j<-pmax(j,1L)
  list(j=j,dx=x-grid[j],dim=dim(x))
}
# make.curve(): log hazard B %*% phi on the grid; h = exp(log hazard); H = cumulative trapezoid
# sums (exact for the piecewise-linear h); NULL if the log hazard exceeds +-700 (overflow).
make.curve <- function(b,phi) {
  lh<-as.vector(b$B%*%phi)
  if(any(!is.finite(lh))||any(lh>700|lh< -700))return(NULL)
  h<-exp(lh); dt<-diff(b$grid); slope<-diff(h)/dt
  H<-c(0,cumsum(dt*(head(h,-1)+tail(h,-1))/2))
  if(any(!is.finite(H)))return(NULL)
  list(h=h,H=H,slope=slope)
}
# eval.curve(): at any t, h = h_j + slope_j dx and H = H_j + h_j dx + slope_j dx^2/2.
eval.curve <- function(curve,idx) {
  j<-idx$j; dx<-idx$dx
  h<-curve$h[j]+curve$slope[j]*dx
  H<-curve$H[j]+curve$h[j]*dx+curve$slope[j]*dx^2/2
  dim(h)<-dim(H)<-idx$dim; list(logh=log(h),H=H)
}

# nodes(): quadrature points for each subject with delta1 = 1 (rows) and 2 x nq nodes (columns).
# (L, R*] with R* = min(R, Y2) is split at its midpoint; each half is mapped onto (0,1):
#   left  u = L  + (mid - L) x^pL,   right  u = R* - (R* - mid) x^pR,
# with Jacobians (mid - L) pL x^(pL-1) and (R* - mid) pR x^(pR-1) (Chapter 3).
# Weibull only: pL = 1/alpha1 if L = 0 and pR = 1/alpha3 if R* is an observed death time; the
# Jacobian then cancels the singular factor u^(alpha1-1) or (y2-u)^(alpha3-1). Otherwise p = 1.
# Returns u, s = y2 - u, their logs, and logw = log(Gauss-Legendre weight x Jacobian).
nodes <- function(d,nq,alpha=c(1,1,1),weibull=FALSE) {
  ii<-which(d$d1==1); lo<-d$L1[ii]; hi<-pmin(d$R1[ii],d$Y2[ii]); mid<-(lo+hi)/2
  gl<-gauss.legendre(nq); lv<-log(gl$x); w<-log(gl$w)
  pL<-ifelse(weibull & lo==0,1/alpha[1],1)
  pR<-ifelse(weibull & hi==d$Y2[ii] & d$d2[ii]==1,1/alpha[3],1)
  # outer(p, log x) gives a subjects x nodes table of p_i log x_j, i.e. log(x_j^p_i).
  left.dx<-exp(outer(pL,lv)+log(mid-lo)); right.dx<-exp(outer(pR,lv)+log(hi-mid))
  uL<-left.dx+lo; sL<-d$Y2[ii]-uL
  sR<-right.dx+(d$Y2[ii]-hi); uR<-hi-right.dx
  # Near a singular endpoint, log u (or log s) is computed directly from its parts for accuracy.
  loguL<-log(uL); loguL[lo==0,]<-(outer(pL,lv)+log(mid-lo))[lo==0,,drop=FALSE]
  logsR<-log(sR); edge<-hi==d$Y2[ii]
  logsR[edge,]<-(outer(pR,lv)+log(hi-mid))[edge,,drop=FALSE]
  # log(weight x Jacobian) = log w_j + log(half length) + log p + (p-1) log x_j.
  weights<-function(p,len)outer(p-1,lv)+log(p)+log(len)+matrix(w,length(ii),nq,byrow=TRUE)
  list(u=cbind(uL,uR),s=cbind(sL,sR),logu=cbind(loguL,log(uR)),
       logs=cbind(log(sL),logsR),logw=cbind(weights(pL,mid-lo),weights(pR,hi-mid)),ii=ii)
}

# objective(): builds the penalized negative log-likelihood; returns a function of the parameter
# vector. Names follow Chapter 3 and Appendix A: K2_y2 = K2(y2), K12_u = K12(u), K1_uy2 = K1(u, y2),
# log_lambda0k = log baseline hazard; case12 = subjects with delta1 = 1, case34 = delta1 = 0.
objective <- function(d,model="Weibull",basis=NULL,set=settings(),nq=set$quadrature) {
  validate.data(d); X<-as.matrix(d[COVARIATES]); p<-ncol(X); case12<-which(d$d1==1); case34<-which(d$d1==0)
  iswb<-identical(model,"Weibull"); nB<-if(iswb)NULL else vapply(basis,function(b)ncol(b$B),integer(1))
  nb<-if(iswb)6L else sum(nB); n0<-nb+as.integer(set$frailty)
  if(!iswb){
    quad<-nodes(d,nq); ixu<-lapply(basis[1:2],function(b)grid.index(b$grid,quad$u))
    ixs<-grid.index(basis[[3]]$grid,quad$s)
    ixy<-lapply(basis[1:2],function(b)grid.index(b$grid,d$Y2[case34]))
  }
  function(par,details=FALSE) {
    fail<-function(reason)if(details)list(ok=FALSE,reason=reason,loglik=NA_real_,nll=Inf) else Inf
    if(length(par)!=n0+3*p||any(!is.finite(par)))return(fail("invalid parameters"))
    theta<-if(set$frailty)exp(par[nb+1L]) else 0
    if(!is.finite(theta))return(fail("frailty overflow"))
    eta<-lapply(0:2,function(k)as.vector(X%*%par[n0+k*p+seq_len(p)]))
    if(any(!is.finite(unlist(eta)))||any(abs(unlist(eta))>700))return(fail("linear predictor overflow"))
    exp_eta<-lapply(eta,exp); loglik_i<-numeric(nrow(d))
    if(iswb){
      log_kappa<-par[c(1,3,5)]; log_alpha<-par[c(2,4,6)]; alpha<-exp(log_alpha)
      if(any(!is.finite(alpha))||any(alpha==0))return(fail("invalid Weibull shape"))
      quad<-nodes(d,nq,alpha,TRUE)
      Lambda0<-function(k,lt)exp(log_kappa[k]+alpha[k]*lt)                  # Lambda0k(t) = kappa_k t^alpha_k (lt = log t)
      log_lambda0<-function(k,lt)log_alpha[k]+log_kappa[k]+(alpha[k]-1)*lt  # log lambda0k(t) = log(alpha_k kappa_k t^(alpha_k - 1))
      K2_y2<-Lambda0(1,log(d$Y2[case34]))*exp_eta[[1]][case34]+Lambda0(2,log(d$Y2[case34]))*exp_eta[[2]][case34]
      log_lambda02_y2<-log_lambda0(2,log(d$Y2[case34]))
      K12_u<-Lambda0(1,quad$logu)*exp_eta[[1]][case12]+Lambda0(2,quad$logu)*exp_eta[[2]][case12]
      K1_uy2<-Lambda0(3,quad$logs)*exp_eta[[3]][case12]   # s = y2 - u (semi-Markov sojourn time)
      K_u<-K12_u+K1_uy2
      log_lambda01_u<-log_lambda0(1,quad$logu); log_lambda03_s<-log_lambda0(3,quad$logs)
    }else{
      ends<-cumsum(nB); begins<-c(0,head(ends,-1))
      curves<-lapply(1:3,function(k)make.curve(basis[[k]],par[begins[k]+seq_len(nB[k])]))
      if(any(vapply(curves,is.null,logical(1))))return(fail("baseline overflow"))
      yy<-lapply(1:2,function(k)eval.curve(curves[[k]],ixy[[k]]))
      uu<-lapply(1:2,function(k)eval.curve(curves[[k]],ixu[[k]])); ss<-eval.curve(curves[[3]],ixs)
      K2_y2<-yy[[1]]$H*exp_eta[[1]][case34]+yy[[2]]$H*exp_eta[[2]][case34];log_lambda02_y2<-yy[[2]]$logh
      K12_u<-uu[[1]]$H*exp_eta[[1]][case12]+uu[[2]]$H*exp_eta[[2]][case12]
      K1_uy2<-ss$H*exp_eta[[3]][case12]
      K_u<-K12_u+K1_uy2
      log_lambda01_u<-uu[[1]]$logh;log_lambda03_s<-ss$logh
    }
    if(any(!is.finite(K_u))||any(!is.finite(K2_y2)))return(fail("cumulative hazard overflow"))
    delta2_34<-d$d2[case34]; delta2_12<-d$d2[case12]
    # Cases 3-4: log f3 = log lambda2*(y2) + log[1+theta*K2]^(-1/theta-1); log f4 = log[1+theta*K2]^(-1/theta).
    loglik_i[case34]<-delta2_34*(log_lambda02_y2+eta[[2]][case34])+moment(K2_y2,theta,delta2_34)
    # Cases 1-2, at each node u: log lambda1*(u) + frailty term (m = 1 for Case 2, m = 2 for Case 1).
    log_integrand<-log_lambda01_u+eta[[1]][case12]+moment(K_u,theta,1+delta2_12)
    # Case 1 only: + log lambda3*(y2 | u). Added only for deaths, avoiding 0 * Inf at an endpoint.
    death<-which(delta2_12==1); if(length(death))log_integrand[death,]<-log_integrand[death,,drop=FALSE]+log_lambda03_s[death,,drop=FALSE]+eta[[3]][case12][death]
    # Gauss-Legendre: log sum_j w_j J_j g(u_j), approximately log of the integral of g(u) over (L, R].
    loglik_i[case12]<-rowlogsum(log_integrand+quad$logw)
    if(any(!is.finite(loglik_i)))return(fail("non-finite log contribution"))
    # Objective = -sum of log f_k + ridge penalty lambda * sum(par^2).
    penalty<-set$lambda*sum(par^2); nll<--sum(loglik_i)+penalty
    if(details)list(ok=TRUE,loglik=sum(loglik_i),penalty=penalty,penalized.loglik=-nll,nll=nll,
      min.log.contribution=min(loglik_i),floor.count=0L,clamp.count=0L,contributions=loglik_i) else nll
  }
}

# weibull.start(): starting values from three separate Weibull regressions (survreg), using the
# interval midpoint for the illness time: transition 1 (midpoint, d1), transition 2 (Y2, d2),
# transition 3 (Y2 - midpoint, d2; subjects with d1 = 1). AFT estimates are converted to the
# proportional hazards form: alpha = 1/sigma, log kappa = -alpha b0, beta = -alpha b.
# log(theta) starts at log(0.01). These are starting values only.
weibull.start <- function(d,set=settings()) {
  mid<-ifelse(d$d1==1,(d$L1+ifelse(is.finite(d$R1),d$R1,d$Y2))/2,d$Y2)
  X<-as.matrix(d[COVARIATES]); ill<-d$d1==1
  fits<-list(survival::survreg(survival::Surv(mid,d$d1)~X,dist="weibull"),
    survival::survreg(survival::Surv(d$Y2,d$d2)~X,dist="weibull"),
    survival::survreg(survival::Surv(pmax(d$Y2[ill]-mid[ill],.Machine$double.eps),d$d2[ill])~X[ill,,drop=FALSE],dist="weibull"))
  aa<-vapply(fits,function(f)1/f$scale,numeric(1))
  base<-unlist(lapply(1:3,function(k)c(-aa[k]*coef(fits[[k]])[1],log(aa[k]))),use.names=FALSE)
  names(base)<-as.vector(rbind(paste0("log_kappa",1:3),paste0("log_alpha",1:3)))
  beta<-unlist(lapply(1:3,function(k)unname(-aa[k]*coef(fits[[k]])[-1])))
  names(beta)<-unlist(lapply(1:3,function(k)paste0("beta",k,"_",COVARIATES)))
  c(base,if(set$frailty)c(log_theta=log(set$theta.starts[1])),beta)
}
# bs.start(): B-spline starting control points = least-squares projection (QR) of the fitted
# Weibull log baseline hazard onto the spline basis; beta and log(theta) are copied from w.
bs.start <- function(w,basis,set=settings()) {
  phi<-unlist(lapply(1:3,function(k){
    a<-exp(w[paste0("log_alpha",k)]); b<-basis[[k]]
    lh<-w[paste0("log_alpha",k)]+w[paste0("log_kappa",k)]+(a-1)*log(b$init.grid)
    v<-as.vector(qr.solve(b$B.init,lh));setNames(v,paste0("phi",k,"_",seq_along(v)))
  }))
  c(phi,if(set$frailty)w["log_theta"],w[grep("^beta",names(w))])
}
numeric.gradient <- function(q,fn,h=1e-4) vapply(seq_along(q),function(j){
  a<-b<-q;a[j]<-a[j]+h;b[j]<-b[j]-h;(fn(a)-fn(b))/(2*h)
},numeric(1))
# inference(): standard errors and numerical checks at the estimate. The Hessian (curvature of
# the objective) is computed twice by finite differences of the analytic gradient, with steps
# 0.001 and 0.0005. SE = sqrt(diag(H^-1)).
inference <- function(par,obj,set=settings(),gradient=NULL) {
  scale<-rep(1,length(par));scale[grepl("_age$",names(par))]<-.02
  q<-par/scale;fn<-function(x)obj(x*scale)
  gr<-if(is.null(gradient))NULL else function(x)gradient(x*scale)*scale
  getH<-function(step)tryCatch(optimHess(q,fn,gr=gr,control=list(ndeps=rep(step,length(q)))),error=function(e)NULL)
  H1<-getH(set$hessian.step);H2<-getH(set$hessian.step/2)
  # inspect(): symmetrise H; Cholesky succeeds only if H is positive definite (a true minimum).
  inspect<-function(H){
    if(is.null(H)||any(!is.finite(H)))return(NULL)
    H<-(H+t(H))/2;ch<-tryCatch(chol(H),error=function(e)NULL)
    if(is.null(ch))return(NULL)
    V<-chol2inv(ch);if(any(!is.finite(V))||any(diag(V)<=0))return(NULL)
    list(H=H,V=V,se=sqrt(diag(V))*scale,eigen=eigen(H,symmetric=TRUE,only.values=TRUE)$values)
  }
  a<-inspect(H1);b<-inspect(H2);pd<-!is.null(a)&&!is.null(b)
  # curvature(): no eigenvalue of H is materially negative.
  curvature<-function(H){
    if(is.null(H)||any(!is.finite(H)))return(FALSE)
    ev<-eigen((H+t(H))/2,symmetric=TRUE,only.values=TRUE)$values
    min(ev)>= -set$curvature.tol*max(1,max(abs(ev)))
  }
  curvature.ok<-curvature(H1)&&curvature(H2)
  grad<-if(is.null(gr))numeric.gradient(q,fn) else gr(q)
  # rel: largest relative change in SE between the two step sizes (must be <= 0.1).
  # gstat: standardized gradient max|gradient x SE| (must be <= 0.1).
  rel<-if(pd)max(abs(a$se-b$se)/pmax(a$se,b$se)) else Inf
  gstat<-if(pd&&all(is.finite(grad)))max(abs(grad)*sqrt(diag(b$V))) else Inf
  ok<-pd && rel<=set$se.check.tol && gstat<=set$gradient.tol
  se<-if(pd)setNames(b$se,names(par)) else setNames(rep(NA_real_,length(par)),names(par))
  V<-if(pd)b$V*outer(scale,scale) else NULL
  H<-if(!is.null(H2))H2/outer(scale,scale) else NULL
  list(ok=ok,pd=pd,curvature.ok=curvature.ok,H=H,H.scaled=H2,H.scaled.other=H1,V=V,se=se,
    gradient=grad,standardized.gradient=gstat,se.relative.change=rel,
    scaled.condition=if(pd)max(b$eigen)/min(b$eigen) else Inf,
    # theta.wald.weak flags SE(log theta) > 2, where a Wald interval for theta is unreliable.
    theta.wald.weak=if("log_theta"%in%names(se))!is.finite(se["log_theta"])||se["log_theta"]>2 else FALSE)
}

# fit(): minimise the penalized negative log-likelihood, then run the numerical checks.
fit <- function(d,model="Weibull",config=NULL,start=NULL,set=settings()) {
  t0<-proc.time()[3];iswb<-identical(model,"Weibull")
  basis<-if(iswb)NULL else make.basis(d,config,set)
  if(is.null(start))start<-weibull.start(d,set)
  if(!iswb && !any(grepl("^phi",names(start))))start<-bs.start(start,basis,set)
  obj<-objective(d,model,basis,set)
  gradient<-if(iswb)gradient.wb(d,set) else gradient.bs(d,basis,set)
  if(!is.finite(obj(start)))stop("Invalid starting values; no sentinel optimization is attempted")
  # The three age coefficients are optimised on the scale beta/0.02 to balance parameter scales.
  scale<-rep(1,length(start));scale[grepl("_age$",names(start))]<-.02
  fn<-function(q)obj(setNames(q*scale,names(start)))
  gr<-function(q)gradient(setNames(q*scale,names(start)))*scale
  run<-function(s){
    q<-s/scale;trace<-numeric(); code<-1L;msg<-"";capped<-FALSE
    cat("  start log(theta)=",if(set$frailty)s["log_theta"] else NA_real_," objective=",fn(q),"\n")
    # Stage 1: nlminb with the analytic gradient.
    initial<-nlminb(q,fn,gradient=gr,control=list(eval.max=set$nl.max,iter.max=set$nl.max,rel.tol=set$nl.rel,x.tol=set$nl.rel))
    if(is.finite(initial$objective)&&initial$objective<=fn(q)){
      q<-initial$par;code<-initial$convergence;msg<-initial$message
    }
    # Stage 2 (B-spline only): up to 8 Nelder-Mead runs, stopping after >= 2 runs once the gain
    # is < 0.001; a run that increases the objective is discarded.
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
    # Stage 3: final nlminb polish. Each stage is kept only if it does not increase the objective.
    z<-nlminb(q,fn,gradient=gr,control=list(eval.max=set$nl.max,iter.max=set$nl.max,rel.tol=set$nl.rel,x.tol=set$nl.rel))
    if(is.finite(z$objective)&&z$objective<=fn(q)){q<-z$par;code<-z$convergence;msg<-z$message}
    cat("  end objective=",fn(q)," code=",code,"\n")
    list(par=q*scale,nll=fn(q),code=code,message=msg,nm.trace=trace,nm.capped=capped)
  }
  # Starting values for theta: 0.01 and 0.25 (Weibull: 2 starts); for the B-spline model also
  # the Weibull estimate (3 starts). The result with the smallest objective is retained.
  starts<-list(start)
  if(set$frailty)for(th in set$theta.starts){s<-start;s["log_theta"]<-log(th);if(!isTRUE(all.equal(s,start)))starts[[length(starts)+1L]]<-s}
  trials<-lapply(starts,function(s)tryCatch(run(s),error=function(e)list(nll=Inf,code=99L,message=conditionMessage(e))))
  best<-trials[[which.min(vapply(trials,function(z)z$nll,numeric(1)))]]
  if(!is.finite(best$nll))stop("All optimization starts failed")
  est<-setNames(best$par,names(start));details<-obj(est,TRUE)
  # Numerical checks. point.ok: optimizer code 0, objective change <= 0.001 when the nodes are
  # doubled (64 -> 128) and <= 0.01 when the B-spline grid is doubled, non-negative curvature and
  # max|score| <= 0.05. inference.ok additionally requires positive definite Hessians, SE change
  # <= 10% and standardized gradient <= 0.1. Ineligible fits are kept and reported, not replaced.
  check.obj<-objective(d,model,basis,set,nq=2L*set$quadrature)
  check.gradient<-if(iswb)gradient.wb(d,set,nq=2L*set$quadrature) else gradient.bs(d,basis,set,nq=2L*set$quadrature)
  quad.diff<-abs(check.obj(est)-obj(est))
  grid.diff<-0
  if(!iswb){s2<-set;s2$n.grid<-2L*set$n.grid;fine<-make.basis(d,config,s2)
    grid.diff<-abs(objective(d,model,fine,set,nq=2L*set$quadrature)(est)-check.obj(est))}
  inf<-inference(est,check.obj,set,check.gradient)
  point.ok<-best$code==0L && details$ok && is.finite(quad.diff)&&quad.diff<=set$quad.check.tol &&
    is.finite(grid.diff)&&grid.diff<=set$grid.check.tol && inf$curvature.ok && all(is.finite(inf$gradient))&&max(abs(inf$gradient))<=set$score.tol
  eligible<-point.ok && inf$ok
  np<-length(est)
  list(version=VERSION,model=model,config=config,est=est,start=start,
    loglik=details$loglik,penalized.loglik=details$penalized.loglik,penalty=details$penalty,
    npar=np,aic=2*np-2*details$loglik,code=best$code,message=best$message,
    point.ok=point.ok,inference.ok=eligible,inference=inf,
    quadrature.difference=quad.diff,grid.difference=grid.diff,
    nm.trace=best$nm.trace,nm.capped=best$nm.capped,
    starts=lapply(trials,function(z)z[c("nll","code","message","nm.trace","nm.capped")]),
    basis=basis,settings=set,minutes=(proc.time()[3]-t0)/60,
    final.diagnostics=details[setdiff(names(details),"contributions")])
}

estimate.table <- function(f) {
  p<-f$est;s<-f$inference$se
  beta<-grepl("^beta",names(p));frailty<-names(p)=="log_theta"
  data.frame(Parameter=names(p),Est=unname(p),SE=unname(s),z=ifelse(frailty,NA_real_,unname(p/s)),p=ifelse(frailty,NA_real_,2*pnorm(-abs(unname(p/s)))),
    CI.lo=ifelse(frailty,NA_real_,unname(p-1.96*s)),CI.hi=ifelse(frailty,NA_real_,unname(p+1.96*s)),HR=ifelse(beta,exp(unname(p)),NA_real_),
    HR.lo=ifelse(beta,exp(unname(p-1.96*s)),NA_real_),HR.hi=ifelse(beta,exp(unname(p+1.96*s)),NA_real_),
    inference_usable=f$inference.ok & !frailty,
    role=ifelse(beta,"regression",ifelse(frailty,"frailty diagnostic; profile required","baseline")),row.names=NULL)
}
fit.row <- function(f) data.frame(Model=f$model,p=f$npar,loglik=f$loglik,
  penalty=f$penalty,penalized_loglik=f$penalized.loglik,AIC=f$aic,code=f$code,
  point_ok=f$point.ok,inference_ok=f$inference.ok,hessian_pd=f$inference$pd,curvature_ok=f$inference$curvature.ok,
  standardized_gradient=f$inference$standardized.gradient,se_relative_change=f$inference$se.relative.change,
  scaled_condition=f$inference$scaled.condition,quadrature_difference=f$quadrature.difference,
  grid_difference=f$grid.difference,nm_capped=f$nm.capped,minutes=f$minutes,row.names=NULL)

# summarize(): one row per parameter over all replications of a scenario.
# Point-estimate summaries (Est, Bias, SEe = SD of estimates) use fits with point.ok;
# SEa (mean model SE) and coverage use fits with inference.ok. Coverage = proportion with
# |est - truth| <= 1.96 SE (95% Wald interval). CI_success_fraction and
# Covered_and_CI_available_fraction use all attempted replications as the denominator.
summarize <- function(results,truth,model) {
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
