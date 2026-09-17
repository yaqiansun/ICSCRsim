argv <- commandArgs(FALSE)
script <- sub('^--file=','',argv[grepl('^--file=',argv)])
root <- dirname(dirname(dirname(normalizePath(script))))
source(file.path(root,'analysis0914/R/core.R'))
source(file.path(root,'analysis0914/R/io.R'))
checks <- 0L
check <- function(ok,label) {
  if(!isTRUE(ok)) stop('FAIL: ',label)
  checks <<- checks+1L
  cat('PASS:',label,'\n')
}
plan <- c(5,12,19,26)
a <- observe.visits0914(plan,3,20)
b <- observe.visits0914(plan,10,20)
c <- observe.visits0914(plan,18,15)
e <- observe.visits0914(plan,13,15)
f <- observe.visits0914(plan,Inf,3)
check(identical(unname(a[1:3]),c(0,5,1)), 'left-censored detected onset')
check(identical(unname(b[1:3]),c(5,12,1)), 'interval-censored detected onset')
check(identical(unname(c[1:3]),c(15,Inf,0)), 'no onset before end: censor at Y2')
check(identical(unname(e[1:3]),c(15,Inf,0)), 'undetected onset: same working observation, no terminal exam')
check(identical(unname(f[1:3]),c(3,Inf,0)), 'no completed visit still uses Y2 convention')
check(observe.visits0914(plan,12,12)['R1']==12, 'positive visit exactly at follow-up endpoint retained')
check(observe.visits0914(plan,12,12)['L1']==5, 'onset at visit has previous negative endpoint')
check(identical(names(formals(visit.plan0914)),c('administrative.end','mean.gap')),
      'schedule API has no event-time input')
set.seed(914001); v <- visit.plan0914(20000,7)
set.seed(914001); vv <- visit.plan0914(20000,7)
check(identical(v,vv) && all(diff(c(0,v))>0) && tail(v,1)>20000 &&
      sum(v>20000)==1, 'deterministic strictly increasing plan extends beyond C')
s <- settings0914(); s$n <- 200L
d <- generate0914(1,'knot22',s)
check(identical(d,generate0914(1,'knot22',s)), 'replication is deterministic')
check(identical(d,generate0914(1,'knot23',s)), 'knot scenarios share identical data')
check(all(d$L1[d$d1==0]==d$Y2[d$d1==0]) && all(d$R1[d$d1==0]==Inf),
      'all undetected cases use Y2/Inf')
check(all(is.infinite(d$T1.working[d$d1==0 & d$d2==1])), 'Case 3 working interpretation recorded')
if(file.exists(file.path(root,'analysis0910/R/core.R'))) {
  old <- new.env()
  source(file.path(root,'analysis0910/R/core.R'),local=old)
  os <- old$settings0910(); os$n <- s$n
  original <- old$generate0910(1,'knot22',os)
  fields <- c('id','T1','T2','Y2','C','d2',COV0914)
  check(identical(d[fields],original[fields]), 'latent event and covariate draws unchanged from 0910')
}
out <- file.path(root,'analysis0914/validation')
dir.create(out,showWarnings=FALSE)
scenarios <- c('knot02','knot12','knot22','knot23','theta0001','theta05','beta1age008')
rows <- lapply(scenarios,function(sc) {
  x <- generate0914(2,sc,s)
  check(isTRUE(validate.data0914(x)),paste('model input validation:',sc))
  a <- x$d1==1
  check(all(x$L1[a]<x$T1[a] & x$T1[a]<=x$R1[a] & x$R1[a]<=x$Y2[a]),
        paste('detected intervals:',sc))
  data.frame(scenario=sc,n=nrow(x),detected=sum(a),deaths=sum(x$d2),
             left=sum(a & x$L1==0),interval=sum(a & x$L1>0),right=sum(!a),
             type1=sum(a & x$d2==1),type2=sum(a & x$d2==0),
             type3=sum(!a & x$d2==1),type4=sum(!a & x$d2==0),
             latent_onset=sum(x$truly_ill),missed_onset=sum(x$missed_ill))
})
write.csv(do.call(rbind,rows),file.path(out,'generation_smoke_counts.csv'),row.names=FALSE)
# Evaluate the retained likelihood without optimization on newly generated data.
truth <- truth0914('knot22')
par <- c(as.vector(rbind(c(-9.23,-11.99,-7.81),c(-.56,-.12,-.42))),
         log(.25), unname(truth[1:9]))
check(is.finite(objective0914(d,set=s)(par)), 'retained Weibull likelihood accepts generated observations')
bs <- basis0914(d,scenarios0914$knot22,s)
bp <- c(rep(log(.0001),sum(vapply(bs,function(x)ncol(x$B),integer(1)))),
        log(.25),unname(truth[1:9]))
check(is.finite(objective0914(d,'Bspline',bs,s)(bp)), 'retained B-spline likelihood accepts generated observations')
cat('PASS:',checks,'generation checks; no model fitting performed.\n')
