#!/usr/bin/env Rscript
# Explicit launch modes only. Scientific code is not modified by this wrapper.
arg <- commandArgs(trailingOnly=TRUE)
if(length(arg)!=1L || !arg %in% c("--smoke","--full")) {
  cat("Usage: Rscript --vanilla reproduce.R --smoke | --full\n")
  quit(save="no",status=if(length(arg)==0L) 0L else 2L)
}
cmd <- commandArgs(FALSE)
file <- sub("^--file=","",cmd[grepl("^--file=",cmd)])
stopifnot(length(file)==1L)
root <- dirname(normalizePath(file)); setwd(root)
stopifnot(requireNamespace("survival",quietly=TRUE),requireNamespace("splines",quietly=TRUE))
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
full <- arg=="--full"
out <- file.path(root,"results",if(full) "production" else "smoke")
if(dir.exists(out)) stop("Output exists. Preserve it; inspect the run before resuming with the underlying entry point.")
scenarios <- if(full) c("knot02","knot12","knot22","knot23","theta0001","theta05","beta1age008") else "knot22"
nrep <- if(full) "1000" else "1"
run <- function(...) {
  status <- system2(file.path(R.home("bin"),"Rscript"),
    shQuote(c("--vanilla","analysis0914/run_simulation.R",...)))
  if(status!=0L) stop("Execution failed; see preceding output. Remaining work stopped.")
}
dir.create(out,recursive=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
for(s in scenarios) {
  run("run",s,"1",nrep,"--n=5000",paste0("--n.rep=",nrep),paste0("--out=",out))
  run("combine",s,"--n=5000",paste0("--n.rep=",nrep),paste0("--out=",out))
}
