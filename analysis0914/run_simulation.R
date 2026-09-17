#!/usr/bin/env Rscript
# A replicate contains both models fitted to the exact same generated dataset.
argv<-commandArgs(FALSE);script<-sub("^--file=","",argv[grepl("^--file=",argv)])
if(length(script)!=1L)stop("Run with Rscript, not source()")
here<-dirname(normalizePath(script));root<-dirname(here)
core<-file.path(root,"analysis0914","R","core.R");io<-file.path(root,"analysis0914","R","io.R")
source(core);source(io)
a<-args0914(commandArgs(TRUE));mode<-if(length(a$pos))a$pos[1] else "help";set<-configure0914(a$opt)
if(mode=="help"){
  cat("Rscript analysis0914/run_simulation.R run SCENARIO START END [--out=ROOT] [--n=5000 --n.rep=1000]\n",
      "Rscript analysis0914/run_simulation.R combine SCENARIO [same settings and output root]\n",
      "Use a separate --out for a pilot. Each replicate is immutable, verified, and independently resumable.\n")
  quit(save="no")
}
if(length(a$pos)<2)stop("Specify scenario")
scenario<-a$pos[2];if(!scenario%in%names(scenarios0914))stop("Unknown scenario")
sc<-scenarios0914[[scenario]]
manifest<-manifest0914(c(core,io,normalizePath(script)),set,list(kind="simulation",scenario=scenario,configuration=sc))
baseout<-if(!is.null(a$opt$out))a$opt$out else file.path(here,"results")
out<-file.path(baseout,scenario);setup.output0914(out,manifest)
with.log0914(out,paste(mode,Sys.getenv("SLURM_ARRAY_TASK_ID"),sep="_"),function(){
  if(mode=="combine"){simulation.outputs0914(out,manifest,scenario,set);return(invisible(NULL))}
  if(mode!="run"||length(a$pos)<4)stop("Use run SCENARIO START END")
  first<-as.integer(a$pos[3]);last<-as.integer(a$pos[4])
  if(anyNA(c(first,last))||first<1||last<first||last>set$n.rep)stop("Invalid replicate range")
  for(rep in seq.int(first,last)){
    path<-file.path(out,"replicates",sprintf("rep_%04d.rds",rep))
    if(file.exists(path)){z<-readRDS(path);check.manifest0914(z$manifest,manifest)
      if(z$rep!=rep)stop("Replicate ID mismatch");cat("Rep",rep,"verified cached result\n");next}
    lock<-paste0(path,".lock");if(!dir.create(lock,showWarnings=FALSE))stop("Replicate already running, or stale lock: ",lock)
    tryCatch({
      cat("Rep",rep,"seed",set$seed.offset+rep,"\n");d<-generate0914(rep,scenario,set)
      errors<-list(wb="",bs="")
      w<-fit.safe0914(d,set=set)
      if(!successful.object0914(w)){errors$wb<-w$error;w<-NULL}
      start<-if(!is.null(w))w$est else NULL
      # A failed Weibull fit is explicitly logged; same midpoint survreg start is used as fallback.
      b<-fit.safe0914(d,"Bspline",sc,start,set)
      if(!successful.object0914(b)){errors$bs<-b$error;b<-NULL}
      obs<-data.frame(n=nrow(d),d11=sum(d$d1==1&d$d2==1),d01=sum(d$d1==0&d$d2==1),
        d10=sum(d$d1==1&d$d2==0),d00=sum(d$d1==0&d$d2==0),truly_ill=sum(d$truly_ill),
        missed_ill=sum(d$missed_ill),missed_ill_death=sum(d$missed_ill==1&d$d2==1),
        missed_ill_censored=sum(d$missed_ill==1&d$d2==0),left_censored=sum(d$d1==1&d$L1==0))
      result<-list(manifest=manifest,rep=as.integer(rep),seed=set$seed.offset+rep,visit.seed=set$visit.seed.offset+rep,data.hash=hash0914(d),
        wb=w,bs=b,errors=errors,observations=obs,maxima=c(max(d$Y2),max(d$Y2),max((d$Y2-d$L1)[d$d1==1])))
      atomic0914(result,path)
      print(obs);for(f in list(w,b))if(!is.null(f))print(fit.row0914(f))
      if(any(nzchar(unlist(errors))))print(errors)
    },finally=unlink(lock,recursive=TRUE))
  }
})
