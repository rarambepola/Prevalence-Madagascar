

pcalg_both_steps <- function(obsDat, alpha, last.index, G_0 = NULL, supprMessages = FALSE,
                             n_rff=5, n_rffz=5,
                             alpha_step_2=NULL,
                             n_rff_step_2=NULL,
                             n_rffz_step_2=NULL,
                             filepaths=NULL,
                             n_cl=2){
  source(paste0(filepaths, "pcalg.R"))
  source(paste0(filepaths, "pcalg_laststep.R"))
  
  if(is.null(alpha_step_2)) alpha_step_2 <- alpha
  if(is.null(n_rff_step_2)) n_rff_step_2 <- n_rff
  if(is.null(n_rffz_step_2)) n_rffz_step_2 <- n_rffz
  
  pc_step_1 <- pcalg(obsDat = obsDat[-last.index], 
                     alpha = alpha,
                     G_0 = G_0[-last.index, -last.index],
                     supprMessages=supprMessages,
                     n_rff=n_rff,
                     n_rffz=n_rffz,
                     n_cl=n_cl)
  
  G_1 <- pc_step_1[[1]]
  G_0[-last.index, -last.index] <- G_1

  pc_step_2 <- pcalg.last(obsDat=obsDat,
                          alpha=alpha_step_2,
                          last.index=last.index,
                          G_0=G_0,
                          supprMessages=supprMessages,
                          n_rff=n_rff_step_2,
                          n_rffz=n_rffz_step_2,
                          n_cl=n_cl)
  
  if(!is.null(names(obsDat))){
    colnames(pc_step_2[[1]]) <- names(obsDat)
  }
  return(pc_step_2)
}