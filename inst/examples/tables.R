## Code that generated tables in the vignette.
## User may want to change how replications are
## generated in parallel.




library(doParallel)
library(reshape2)
library(foreach)
library(plyr)
library(sparseMVN)

run.par <- TRUE
if(run.par) registerDoParallel(cores=10) else registerDoParallel(cores=1)

N.vec <- c(200)
m.vec <- c(25, 250, 500)
p.vec <- c(2,5)
k.vec <- c(15)
reps <- 30

## uses functions in R/demo_funcs.R

R <- ldply(1:reps,function(i) run.compare(i,N.vec,m.vec,p.vec,k.vec),
           .parallel=run.par)

tmp <- melt(R,id.vars=c("rep","matrix","m","N","p","k"))
tab <- dcast(tmp,matrix+p+k+m+N~variable,
             fun.aggregate=mean)
