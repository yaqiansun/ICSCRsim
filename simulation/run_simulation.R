#!/usr/bin/env Rscript
# A replicate contains both models fitted to the exact same generated dataset.
# Usage: run SCENARIO START END -> generate and fit replications START..END (each saved separately);
#        combine SCENARIO       -> summarise all saved replications into tables and figures.
# Each replication: generate() the data, fit the Weibull model, then the B-spline model starting
# from the Weibull estimate; the counts of Cases 1-4 and of missed illnesses are saved with the fits.
argv<-commandArgs(FALSE);script<-sub("^--file=","",argv[grepl("^--file=",argv)])
if(length(script)!=1L)stop("Run with Rscript, not source()")
here<-dirname(normalizePath(script));root<-dirname(here)
core<-file.path(root,"simulation","R","core.R");io<-file.path(root,"simulation","R","io.R")
source(core);source(io)
a<-parse.args(commandArgs(TRUE));mode<-if(length(a$pos))a$pos[1] else "help";set<-configure(a$opt)
if(mode=="help"){
  cat("Rscript simulation/run_simulation.R run SCENARIO START END [--out=ROOT] [--n=5000 --n.rep=1000]\n",
      "Rscript simulation/run_simulation.R combine SCENARIO [same settings and output root]\n",
      "Use a separate --out for a pilot. Each replicate is immutable, verified, and independently resumable.\n")
  quit(save="no")
}
if(length(a$pos)<2)stop("Specify scenario")
scenario<-a$pos[2];if(!scenario%in%names(scenarios))stop("Unknown scenario")
sc<-scenarios[[scenario]]
manifest<-make.manifest(c(core,io,normalizePath(script)),set,list(kind="simulation",scenario=scenario,configuration=sc))
baseout<-if(!is.null(a$opt$out))a$opt$out else file.path(here,"results")
out<-file.path(baseout,scenario);setup.output(out,manifest)
with.log(out,paste(mode,Sys.getenv("SLURM_ARRAY_TASK_ID"),sep="_"),function(){
  if(mode=="combine"){simulation.outputs(out,manifest,scenario,set);return(invisible(NULL))}
  if(mode!="run"||length(a$pos)<4)stop("Use run SCENARIO START END")
  first<-as.integer(a$pos[3]);last<-as.integer(a$pos[4])
  if(anyNA(c(first,last))||first<1||last<first||last>set$n.rep)stop("Invalid replicate range")
  for(rep in seq.int(first,last)){
    path<-file.path(out,"replicates",sprintf("rep_%04d.rds",rep))
    if(file.exists(path)){z<-readRDS(path);check.manifest(z$manifest,manifest)
      if(z$rep!=rep)stop("Replicate ID mismatch");cat("Rep",rep,"verified cached result\n");next}
    lock<-paste0(path,".lock");if(!dir.create(lock,showWarnings=FALSE))stop("Replicate already running, or stale lock: ",lock)
    tryCatch({
      cat("Rep",rep,"seed",set$seed.offset+rep,"\n");d<-generate(rep,scenario,set)
      errors<-list(wb="",bs="")
      w<-fit.safe(d,set=set)
      if(!successful.object(w)){errors$wb<-w$error;w<-NULL}
      start<-if(!is.null(w))w$est else NULL
      # A failed Weibull fit is explicitly logged; same midpoint survreg start is used as fallback.
      b<-fit.safe(d,"Bspline",sc,start,set)
      if(!successful.object(b)){errors$bs<-b$error;b<-NULL}
      obs<-data.frame(n=nrow(d),d11=sum(d$d1==1&d$d2==1),d01=sum(d$d1==0&d$d2==1),
        d10=sum(d$d1==1&d$d2==0),d00=sum(d$d1==0&d$d2==0),truly_ill=sum(d$truly_ill),
        missed_ill=sum(d$missed_ill),missed_ill_death=sum(d$missed_ill==1&d$d2==1),
        missed_ill_censored=sum(d$missed_ill==1&d$d2==0),left_censored=sum(d$d1==1&d$L1==0))
      result<-list(manifest=manifest,rep=as.integer(rep),seed=set$seed.offset+rep,visit.seed=set$visit.seed.offset+rep,data.hash=hash(d),
        wb=w,bs=b,errors=errors,observations=obs,maxima=c(max(d$Y2),max(d$Y2),max((d$Y2-d$L1)[d$d1==1])))
      atomic.save(result,path)
      print(obs);for(f in list(w,b))if(!is.null(f))print(fit.row(f))
      if(any(nzchar(unlist(errors))))print(errors)
    },finally=unlink(lock,recursive=TRUE))
  }
})
