#!/usr/bin/env Rscript

library("simulation.utils")

# start a list of parameters
# ifyou 
params <- list()
params$sigma0 <- seq(0.01,1,0.1)
params$M1 <- 20
params$M2 <- c(5,20)
params$mu_m <- c(0.01)
params$mu_bh <- c(0,0.01)
params$mu_mb <- c(0,0.01)
params$baseline_survival <- c(0)
params$survival_strength_1 <- .1
params$survival_strength_2 <- 1
params$max_time <- 200000
params$fix_clutch <- c(0)
params$kin_comp <- 0

all.params <- expand.grid(params)

all.params <- all.params[all.params$mu_m > 0 | all.params$mu_bh > 0 | all.params$mu_mb > 0,] 

all.params$m_init <- NA

# make a lookup list with optima
c.coeff <- sort(unique(all.params$survival_strength_1))

for (c.coeff.i in c.coeff)
{
    x.opt <- optimise(f = W, M=10, c=c.coeff.i,mmin=1,interval=c(1,10)) 

    all.params[all.params$survival_strength_1 == c.coeff.i,"m_init"] <- x.opt$minimum
}

all.params$sigma1 <- all.params$sigma0

col.order.orig <- 1:ncol(all.params)

col.sigma1 <- which("sigma1" == names(all.params))
col.sigma0 <- which("sigma0" == names(all.params))

new.col.order <- c(col.order.orig[1:col.sigma0],col.sigma1,col.order.orig[(col.sigma0+1):(length(col.order.orig)-1)])

all.params <- all.params[,new.col.order]



make.batch.file(
                parameter_list=all.params
                ,executable_path="./build/Bethedging_clutch_size"
                ,n_replicates = 4)

