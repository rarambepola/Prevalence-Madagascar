#script to process the results of PC algorithm to produce potential feature sets
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Causal_inference")) #change if necessary

##temporal variables
temp_files <- list.files(path="Outputs", pattern="temp*", full.names = TRUE)

pc_list <- list()
for(temp_file in temp_files){
  load(temp_file)
  pc_list[[length(pc_list) + 1]] <- pc_both[[1]]
}

sum_mat <- Reduce('+', pc_list)
prop_mat <- sum_mat / length(temp_files)


temp_cov_list_full <- list()

n_vars <- dim(pc_list[[1]])[1]
par(mfrow=c(1, 2))
for(i in (15:50)*0.01){
  temp_cov_list_full[[length(temp_cov_list_full) + 1]] <- setdiff(colnames(graph.minimal.n(prop_mat >= i, n_vars)), "logit_pr")
}
temp_cov_list_full <- lapply(temp_cov_list_full, sort)
temp_cov_list <- unique(temp_cov_list_full)


##static variables
temp_files <- list.files(path="Outputs", pattern="static*", full.names = TRUE)

pc_list <- list()
for(temp_file in temp_files){
  load(temp_file)
  pc_list[[length(pc_list) + 1]] <- pc_both[[1]]
}

sum_mat <- Reduce('+', pc_list)
prop_mat <- sum_mat / length(temp_files)


static_cov_list_full <- list()
n_vars <- dim(pc_list[[1]])[1]
par(mfrow=c(1, 2))
for(i in (5:50)*0.01){
  static_cov_list_full[[length(static_cov_list_full) + 1]] <- setdiff(colnames(graph.minimal.n(prop_mat >= i, n_vars)), "logit_pr")
}
static_cov_list_full <- lapply(static_cov_list_full, sort)
static_cov_list <- unique(static_cov_list_full)


##combine sets
causal_sets_all <- list()
for(static_cov in static_cov_list){
  for(temp_cov in temp_cov_list){
    causal_sets_all[[length(causal_sets_all) + 1]] <- c(static_cov, temp_cov)
  }
}

save(causal_sets_all, file="causal_sets_all.RData")

