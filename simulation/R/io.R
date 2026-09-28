# Input/output utilities (no real patient rows are exported):
#   parse.args(), configure(): read command-line arguments and options into settings();
#   setup.output(), with.log(): create output folders, check the manifest, write log files;
#   savefit(), loadfit(), fit.safe(), successful.object(): save and load fits; a failed fit is
#     recorded as an error instead of stopping the simulation;
#   hazard(): fitted baseline hazard of transition k at times t (for figures);
#   real.outputs(): tables, AIC comparison and figures for the real-data fits;
#   simulation.outputs(): combine the replications of a scenario into summary tables and figures.
parse.args <- function(args) {
  opt<-list(); pos<-character()
  for(a in args){if(startsWith(a,"--")){
    p<-strsplit(substring(a,3),"=",fixed=TRUE)[[1]]
    if(length(p)!=2)stop("Options must be --name=value: ",a)
    opt[[p[1]]]<-p[2]
  }else pos<-c(pos,a)}
  list(pos=pos,opt=opt)
}
configure <- function(opt) {
  allowed<-c("n","n.rep","n.grid","quadrature","nm.max","restart.max","nl.max","lambda","theta.starts","out","data","report")
  unknown<-setdiff(names(opt),allowed)
  if(length(unknown))stop("Unknown option(s): ",paste(unknown,collapse=", "))
  set<-settings()
  for(k in intersect(names(opt),c("n","n.rep","n.grid","quadrature","nm.max","restart.max","nl.max"))) {
    set[[k]]<-as.integer(opt[[k]]);if(is.na(set[[k]])||set[[k]]<1)stop("Invalid --",k)
  }
  if(!is.null(opt$lambda)){set$lambda<-as.numeric(opt$lambda);if(!is.finite(set$lambda)||set$lambda<0)stop("Invalid penalty")}
  if(!is.null(opt$theta.starts)){set$theta.starts<-as.numeric(strsplit(opt$theta.starts,",",fixed=TRUE)[[1]])
    if(any(!is.finite(set$theta.starts)|set$theta.starts<=0))stop("Invalid theta starts")}
  if(set$quadrature<2 || set$n.grid<3)stop("quadrature must be >=2 and n.grid >=3")
  set
}
setup.output <- function(out,manifest) {
  dir.create(out,showWarnings=FALSE,recursive=TRUE)
  path<-file.path(out,"manifest.rds")
  # Directory lock for concurrent array tasks creating the same manifest.
  lock<-paste0(path,".lock")
  if(!file.exists(path)){
    owner<-dir.create(lock,showWarnings=FALSE)
    if(owner){on.exit(unlink(lock,recursive=TRUE),add=TRUE)
      if(!file.exists(path))atomic.save(manifest,path)
    }else{
      for(i in 1:100){if(file.exists(path))break;Sys.sleep(.1)}
      if(!file.exists(path))stop("Manifest lock has no completed manifest; inspect ",lock)
    }
  }
  check.manifest(readRDS(path),manifest)
  for(d in c("fits","tables","figures","replicates","logs"))dir.create(file.path(out,d),showWarnings=FALSE)
  invisible(out)
}
with.log <- function(out,label,fn) {
  path<-file.path(out,"logs",paste0(label,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid(),".log"))
  con<-file(path,"wt"); sink(con,split=TRUE);sink(con,type="message")
  on.exit({sink(type="message");sink();close(con)},add=TRUE)
  cat("Version:",VERSION,"\nStarted:",format(Sys.time()),"\n")
  tryCatch(fn(),error=function(e){cat("ERROR:",conditionMessage(e),"\n");stop(e)})
}
savefit <- function(f,path,manifest) {
  atomic.save(list(manifest=manifest,fit=f),path)
}
loadfit <- function(path,manifest) {
  z<-readRDS(path);check.manifest(z$manifest,manifest);z$fit
}
fit.safe <- function(...)tryCatch(fit(...),error=function(e)structure(list(error=conditionMessage(e)),class="fit_error"))
successful.object <- function(x)!is.null(x)&&!inherits(x,"fit_error")&&!is.null(x$est)

hazard <- function(f,k,t) {
  if(f$model=="Weibull"){
    a<-exp(f$est[paste0("log_alpha",k)]);kap<-exp(f$est[paste0("log_kappa",k)])
    return(a*kap*t^(a-1))
  }
  b<-f$basis[[k]];phi<-f$est[grep(paste0("^phi",k,"_"),names(f$est))]
  valid<-t>=0&t<=tail(b$grid,1);out<-rep(NA_real_,length(t))
  out[valid]<-exp(eval.curve(make.curve(b,phi),grid.index(b$grid,t[valid]))$logh);out
}
real.outputs <- function(out,manifest,report="auto") {
  names.all<-c("Weibull",names(ladder))
  files<-file.path(out,"fits",paste0(names.all,".rds"));exists<-file.exists(files)
  fits<-setNames(lapply(files[exists],loadfit,manifest=manifest),names.all[exists])
  if(!length(fits))stop("No completed fits")
  tab<-do.call(rbind,lapply(fits,fit.row));tab<-tab[order(tab$AIC),]
  eligible<-tab$point_ok
  tab$dAIC_all<-tab$AIC-min(tab$AIC)
  tab$dAIC_eligible<-if(any(eligible))tab$AIC-min(tab$AIC[eligible]) else NA_real_
  write.csv(tab,file.path(out,"tables","aic_comparison.csv"),row.names=FALSE)
  complete<-all(exists)
  candidate<-tab$Model!="Weibull" & eligible
  best.bs<-if(any(candidate))tab$Model[which(candidate)[1]] else NA_character_
  best.all<-if(any(eligible))tab$Model[which(eligible)[1]] else NA_character_
  selected<-best.all
  if(report!="auto"){
    if(!report%in%names(fits))stop("Requested reported model is missing")
    if(!isTRUE(fits[[report]]$point.ok))stop("Requested reported model did not pass numerical optimization checks")
    selected<-report
  }
  status<-data.frame(complete_ladder=complete,aic_best_eligible_model=best.all,aic_best_eligible_bspline=best.bs,
    reported_model=selected,reported_inference_usable=if(!is.na(selected))fits[[selected]]$inference.ok else FALSE,
    selection_rule=report,selection_provisional=!complete,
    note="Ordinary AIC at ridge-regularized estimates; frailty boundary and optimization diagnostics must be considered.")
  write.csv(status,file.path(out,"tables","selection.csv"),row.names=FALSE)
  for(nm in names(fits))write.csv(estimate.table(fits[[nm]]),file.path(out,"tables",paste0("estimates_",nm,".csv")),row.names=FALSE)
  # Explicit sensitivity table, coefficient ranges, and significance decisions.
  bnames<-grep("^beta",names(fits[[1]]$est),value=TRUE)
  sens<-data.frame(Parameter=bnames)
  for(nm in names(fits))sens[[nm]]<-fits[[nm]]$est[bnames]
  write.csv(sens,file.path(out,"tables","coefficient_sensitivity.csv"),row.names=FALSE)
  if(!is.na(selected)){
    writeLines(selected,file.path(out,"tables","selected_model.txt"))
  }else if(file.exists(file.path(out,"tables","selected_model.txt"))){
    stop("No model currently eligible but a selection file already exists; use a clean output directory")
  }
  comparison.bs<-if(!is.na(selected)&&selected!="Weibull")selected else best.bs
  if(!is.na(comparison.bs)){
    S<-fits[[comparison.bs]];W<-fits[["Weibull"]]
    pdf(file.path(out,"figures",paste0("baseline_hazards_Weibull_vs_",comparison.bs,".pdf")),width=12,height=4)
    on.exit(dev.off(),add=TRUE);par(mfrow=c(1,3),mar=c(4,4,1,1))
    for(k in 1:3){t<-seq(1,tail(S$basis[[k]]$grid,1),length.out=500);h<-hazard(S,k,t)
      hw<-if(!is.null(W))hazard(W,k,t) else rep(NA_real_,length(t))
      plot(t,h,type="l",lty=2,col="blue",xlab=if(k==3)"Sojourn time (days)" else "Time (days)",ylab="Baseline hazard",ylim=range(c(h,hw),finite=TRUE))
      if(!is.null(W))lines(t,hw);legend("topright",legend=c("Weibull",comparison.bs),col=c("black","blue"),lty=c(1,2),bty="n",cex=.8)}
  }
  print(tab,row.names=FALSE);print(status,row.names=FALSE)
  invisible(list(fits=fits,table=tab,selection=status))
}

simulation.outputs <- function(out,manifest,scenario,set) {
  files<-list.files(file.path(out,"replicates"),"^rep_[0-9]+\\.rds$",full.names=TRUE)
  if(!length(files))stop("No completed replicate files")
  results<-lapply(files,readRDS)
  for(z in results)check.manifest(z$manifest,manifest)
  ids<-vapply(results,`[[`,integer(1),"rep")
  if(anyDuplicated(ids)||any(ids<1|ids>set$n.rep))stop("Duplicate/out-of-range replicate IDs")
  results<-results[order(ids)];ids<-sort(ids)
  truth<-true.values(scenario,set$frailty)
  diag<-do.call(rbind,lapply(results,function(z){
    do.call(rbind,lapply(c("wb","bs"),function(m){f<-z[[m]]
      if(!successful.object(f))return(data.frame(rep=z$rep,model=m,error=z$errors[[m]],code=NA_integer_,point_ok=FALSE,inference_ok=FALSE,
        hessian_pd=FALSE,quad_difference=NA_real_,grid_difference=NA_real_,gradient=NA_real_,se_change=NA_real_))
      data.frame(rep=z$rep,model=m,error="",code=f$code,point_ok=f$point.ok,inference_ok=f$inference.ok,hessian_pd=f$inference$pd,
        quad_difference=f$quadrature.difference,grid_difference=f$grid.difference,gradient=f$inference$standardized.gradient,se_change=f$inference$se.relative.change)
    }))
  }))
  write.csv(diag,file.path(out,"tables","replicate_diagnostics.csv"),row.names=FALSE)
  summaries<-lapply(c("wb","bs"),function(m)summarize(results,truth,m));names(summaries)<-c("wb","bs")
  for(m in names(summaries))write.csv(summaries[[m]],file.path(out,"tables",paste0(m,"_summary.csv")),row.names=FALSE)
  event<-do.call(rbind,lapply(results,function(z)cbind(rep=z$rep,z$observations)))
  write.csv(event,file.path(out,"tables","observation_diagnostics.csv"),row.names=FALSE)
  write.csv(data.frame(n_expected=set$n.rep,n_completed=length(ids),complete=identical(ids,seq_len(set$n.rep))),file.path(out,"tables","completion.csv"),row.names=FALSE)
  saveRDS(list(manifest=manifest,settings=set,scenario=scenario,true.par=truth,rep_ids=ids,summary=summaries,diagnostics=diag),file.path(out,"summary.rds"))
  # All plotted curves and means use the same point-estimate eligibility mask.
  for(m in c("wb","bs")){
    ff<-lapply(results,`[[`,m);keep<-vapply(ff,function(f)successful.object(f)&&isTRUE(f$point.ok),logical(1));ff<-ff[keep]
    if(!length(ff))next
    upper<-vapply(1:3,function(k)min(vapply(results[keep],function(z)z$maxima[k],numeric(1))),numeric(1))
    pdf(file.path(out,"figures",paste0(m,"_hazards.pdf")),width=12,height=4)
    par(mfrow=c(1,3),mar=c(4,4,1,1))
    for(k in 1:3){t<-seq(1,upper[k],length.out=500);mat<-vapply(ff,hazard,numeric(length(t)),k=k,t=t)
      a<-exp(c(-.56,-.12,-.42))[k];kap<-exp(c(-9.23,-11.99,-7.81))[k];truth.h<-a*kap*t^(a-1)
      matplot(t,mat,type="l",lty=1,col="grey75",xlab=if(k==3)"Sojourn time (days)" else "Time (days)",ylab="Baseline hazard",ylim=range(c(mat,truth.h),finite=TRUE))
      lines(t,truth.h,col="red",lty=2,lwd=2);lines(t,rowMeans(mat),col="blue",lwd=2)
      legend("topright",legend=c(paste(length(ff),"eligible fits"),"Truth","Mean"),col=c("grey75","red","blue"),lty=c(1,2,1),bty="n",cex=.8)
    };dev.off()
  }
  print(summaries);invisible(summaries)
}
