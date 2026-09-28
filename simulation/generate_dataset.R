#!/usr/bin/env Rscript
# Generate and save one synthetic dataset only; no model is fitted.
argv <- commandArgs(FALSE)
script <- sub('^--file=','',argv[grepl('^--file=',argv)])
if(length(script)!=1L) stop('Run this entry point with Rscript')
here <- dirname(normalizePath(script))
source(file.path(here,'R/core.R')); source(file.path(here,'R/io.R'))
a <- parse.args(commandArgs(TRUE))
if(length(a$pos)<1) stop('Usage: generate_dataset.R SCENARIO [REP] [--n=5000 --out=DIRECTORY]')
scenario <- a$pos[1]
rep <- if(length(a$pos)>1) as.integer(a$pos[2]) else 1L
set <- configure(a$opt)
stopifnot(length(rep)==1L,!is.na(rep),rep>=1,rep<=set$n.rep)
out <- if(is.null(a$opt$out)) file.path(here,'generation_diagnostics') else a$opt$out
file <- file.path(out,paste0(scenario,'_rep_',sprintf('%04d',rep),'.rds'))
if(file.exists(file)) stop('Existing diagnostic: use a new output directory')
d <- generate(rep,scenario,set)
observations <- data.frame(n=nrow(d),left=sum(d$d1==1 & d$L1==0),
  interval=sum(d$d1==1 & d$L1>0),right=sum(d$d1==0),
  type1=sum(d$d1==1 & d$d2==1),type2=sum(d$d1==1 & d$d2==0),
  type3=sum(d$d1==0 & d$d2==1),type4=sum(d$d1==0 & d$d2==0),
  deaths=sum(d$d2),latent_onset=sum(d$truly_ill),missed_onset=sum(d$missed_ill))
manifest <- make.manifest(c(file.path(here,'R/core.R'),file.path(here,'R/io.R'),script),set,
                        list(kind='generation-only',scenario=scenario,rep=rep))
atomic.save(list(manifest=manifest,data=d,observations=observations,
               note='Synthetic diagnostic only; not fitted results.'),file)
print(observations)
cat('Saved:',normalizePath(file,winslash='/'),'\nNo models were fitted.\n')
