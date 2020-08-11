#script to run the PC algorithm with temoral covariates
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Causal_inference")) #change if necessary

args <- commandArgs(trailingOnly=TRUE) #intended to be run from terminal with seed passed as argument
seed <- args[1]
print(paste0("Seed: ", seed))

library(RCITcpp) #available at https://github.com/rarambepola/RCITcpp
library(foreach)
library(doParallel)


source("pcalg_both_steps.R")
load("prewhitened.RData")

#remove unwanted variables and subsample
res_mat <- res_mat[, -(1:3)]
res_mat <- res_mat[sample.int(632, 400), ]

#reformat for algorithm
res_mat_list <- list()
for(i in 1:dim(res_mat)[2]){
  res_mat_list[[i]] <- res_mat[, i]
}

names(res_mat_list) <- colnames(res_mat)
n_vars <- dim(res_mat)[2]

#make G_0
G_0 <- matrix(1, n_vars, n_vars)
timelags <- c(rep(0:3, 7), rep(4, 8), -1)
diag(G_0) <- 0
for(i in 1:n_vars){
  G_0[i, timelags[i] < timelags] <- 0
}


alpha_use <- 0.4 #test level

pr_index <- 37 #which index represents prevalence
temp_index <- c(1:28, 37) #only use temporal variables
#remove 3 month timelag and temperature suitability of plasmodium vivax
temp_index2 <- setdiff(c(setdiff(1:28, (1:7)*4), 37), grep("TSI_pv", names(res_mat_list)))

#run PC algorithm
ptm <- proc.time()
pc_both <- pcalg_both_steps(res_mat_list[temp_index2], alpha=alpha_use, 
                                   last.index=length(temp_index2),
                                   G_0=G_0[temp_index2, temp_index2],
                                   alpha_step_2 = runif(1, 0.7, 0.9),
                                   n_cl=4)
print(proc.time() - ptm)

save(list=c("pc_both"),
     file=paste0("Output/temporal_covariates_", seed, ".RData"))

