argv<-commandArgs(FALSE);script<-sub("^--file=","",argv[grepl("^--file=",argv)])
root<-dirname(dirname(dirname(normalizePath(script))))
source(file.path(root,"analysis0914","R","core.R"));source(file.path(root,"analysis0914","R","io.R"))
checks<-0L
check<-function(ok,label){if(!isTRUE(ok))stop("FAIL: ",label);checks<<-checks+1L;cat("PASS:",label,"\n")}
near<-function(x,y,tol,label)check(max(abs(x-y))<tol,label)
set<-settings0914()
g<-gauss0914(64);near(sum(g$w),1,1e-14,"quadrature weights sum to one")
near(sum(g$w*g$x^6),1/7,1e-14,"quadrature integrates polynomial")
near(moment0914(c(.1,1,100),1e-20,2),-c(.1,1,100),1e-12,"frailty limit is stable near zero")

# Analytic constant-hazard independent illness-death probabilities/densities.
d<-data.frame(L1=c(0,3,0,3),R1=c(2,Inf,2,Inf),d1=c(1,0,1,0),Y2=c(3,3,3,3),d2=c(1,1,0,0),
              age_at_scan=rep(0,4),cs_sex=rep(0,4),outpatient_inpatient=rep(0,4))
k<-c(.2,.3,.4);p<-c(as.vector(rbind(log(k),rep(0,3))),rep(0,9));s<-set;s$frailty<-FALSE;s$lambda<-0
obj<-objective0914(d,set=s)
integ<-.2*exp(-.4*3)*(1-exp(-(.2+.3-.4)*2))/(.2+.3-.4)
truth<-c(.4*integ,.3*exp(-.5*3),integ,exp(-.5*3))
near(obj(p,TRUE)$contributions,log(truth),1e-12,"four patterns match analytic independent model")

# Independent adaptive integration, including origin and sojourn endpoint singularities.
d2<-d;d2$R1[c(1,3)]<-3
for(th in c(0,.0001,.25,2)){
  ss<-set;ss$frailty<-th>0;ss$lambda<-0;ss$quadrature<-128L
  shape<-c(.57,.89,.66);kap<-c(.2,.3,.4)
  par<-c(as.vector(rbind(log(kap),log(shape))),if(th>0)log(th),rep(0,9))
  ob<-objective0914(d2,set=ss)
  for(i in c(1L,3L)){
    fun<-function(u){K<-kap[1]*u^shape[1]+kap[2]*u^shape[2]+kap[3]*(3-u)^shape[3]
      h1<-shape[1]*kap[1]*u^(shape[1]-1);h3<-shape[3]*kap[3]*(3-u)^(shape[3]-1)
      h1*if(i==1)h3*exp(moment0914(K,th,2)) else exp(moment0914(K,th,1))}
    ref<-integrate(fun,0,3,rel.tol=1e-10,abs.tol=1e-12,subdivisions=1000)$value
    near(ob(par,TRUE)$contributions[i],log(ref),2e-5,paste("Weibull endpoint integration theta",th,"pattern",i))
  }
}

# A small, strictly positive individual likelihood is retained, never floored.
toy<-d[3,,drop=FALSE];toy$R1<-1;toy$Y2<-3000
ss<-set;ss$lambda<-0
pa<-c(rep(c(log(.01),0),3),log(.0001),rep(0,9))
ob<-objective0914(toy,set=ss)
check(is.finite(ob(pa))&&ob(pa)>30,"tiny likelihood remains finite below old floor")

# Exact derivative relationship H'=h, including the omitted first day.
b<-list(grid=c(0,1,3),B=diag(3));cu<-curves0914(b,log(c(2,4,8)))
near(evalcurve0914(cu,index0914(b$grid,c(0,1,3)))$H,c(0,3,15),1e-12,"cumulative hazard includes 0 to 1")
t<-c(.2,.8,1.2,2.8);h<-1e-6
deriv<-(evalcurve0914(cu,index0914(b$grid,t+h))$H-evalcurve0914(cu,index0914(b$grid,t-h))$H)/(2*h)
near(deriv,exp(evalcurve0914(cu,index0914(b$grid,t))$logh),1e-8,"interpolated cumulative derivative equals hazard")

# Constant B-spline hazards reproduce the same analytic model.
dd<-d;dd$L1[c(2,4)]<-dd$Y2[c(2,4)]
bs<-basis0914(dd,list(degree=c(1,1,1),internal=c(0,0,0)),s)
bp<-c(rep(log(k[1]),2),rep(log(k[2]),2),rep(log(k[3]),2),rep(0,9))
near(objective0914(dd,"Bspline",bs,s)(bp),-sum(log(truth)),1e-12,"B-spline and Weibull constant baselines agree")

# Inverse-diagonal positivity is not mistaken for positive definiteness.
H<-solve(matrix(c(1,2,2,1),2));fn<-function(x)sum(x*(H%*%x))/2
inf<-inference0914(setNames(c(0,0),c("a","b")),fn,set)
check(!inf$pd&&!inf$ok,"indefinite Hessian is rejected despite positive inverse diagonal")

# Point-estimate and CI denominators are distinct and explicit.
f<-list(point.ok=TRUE,inference.ok=TRUE,est=c(a=1),inference=list(se=c(a=.2)))
f2<-f;f2$inference.ok<-FALSE;f2$est<-c(a=3);f2$inference$se<-c(a=NA_real_)
su<-summarize0914(list(list(wb=f),list(wb=f2),list(wb=NULL)),c(a=1),"wb")
check(su$n_attempted==3&&su$n_point==2&&su$n_ci==1&&su$Est==2,"summary distinguishes point/CI denominators")

# Generation-specific comparisons are in test_generation.R.
for (f in list.files(file.path(root, 'analysis0914/R'), pattern='[.]R$', full.names=TRUE)) parse(f)
parse(file.path(root, 'analysis0914/run_simulation.R'))
cat('ALL', checks, 'CHECKS PASSED', '\n')
