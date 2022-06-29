#!/usr/bin/env Rscript

library("simulation.utils")

# start a list of parameters
params <- list()
params$max_time <- 100000
params$d <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0) #Original (0.1, 0.2, 0.5)

all.params <- expand.grid(params)


make.batch.file(
                parameter_list=all.params
                ,executable_path="./build/simulation_main"
                ,n_replicates = 10)

