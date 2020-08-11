##################################################################################################
########################       Script for 02_raster_prep                 #########################
##################################################################################################

library(raster)
library(rgdal)
library(malariaAtlas)
library(doParallel)

hf.rast.files <-cbind(c(list.files("MDG_time_travel_rasters", pattern ="tif$", full.names = TRUE)))

new_folder <- ("Madagascar_time_travel_rasters_cropped/")
dir.create(new_folder)

MDG <- getShp(ISO = "MDG")  ## shapefile to crop to

######################
##### Step 1: crop and re-write files to remove water from the rasters.
######################

cores <- detectCores() 
cl <- makeCluster(10)   ##if running on personal computer this number will change
registerDoParallel(cl)
start_time <- Sys.time()
b_crop <- foreach(i = 1:nrow(hf.rast.files)) %dopar% {
  
  b <- raster::raster(hf.rast.files[i])
  b_crop <- raster::crop(b, MDG) 
  
  raster::writeRaster(b_crop, filename = file.path(paste0(new_folder, names(b), ".tif"))) 
  
}

end_time <- Sys.time()
end_time - start_time  ## 10mins
stopCluster(cl)

##alternative to above
test <- lapply(hf.rast.files, function(x) raster(x))
test_crop <- lapply(test, function(x) crop(x, MDG))
lapply(test_crop, function(x) writeRaster(x, filename=file.path(paste0(new_folder, names(x), ".tif"))))



## make sure that all of our masked rasters have the same number of cells
rasters_check <- sapply(seq_along(cropped), function(x) ncell(cropped[[x]]))
if(length(unique(rasters_check)) > 1){
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  warning("Check the following files: ", 
          paste(filenames[which(is.na(match(rasters_check, getmode(rasters_check))))], 
                collapse = "\n"), 
          sep = "")
  
}
rm(rasters_check)

#####################################
## Step 2: 1/travel time squared.
#####################################


hf.rast.files <-cbind(c(list.files("Madagascar_time_travel_rasters_cropped", pattern ="tif$", full.names = TRUE)))


nrow <- ncell(raster(hf.rast.files[1]))
ncol <- length(hf.rast.files)
mat <- matrix(nrow = nrow, ncol = ncol)  ## pre defining should speed things up a little.

### save travel time distances as a matrix prior to inverse
start_time <- Sys.time()
for(i in 1:nrow(hf.rast.files)){
  a <- raster(hf.rast.files[[i]])
  x <- as.matrix(a)
  
  ind <- i*1
  mat[, ind] <- unlist(x, use.names = TRUE)
  
}
end_time <- Sys.time()
end_time - start_time 


### inverse travel time distances
start_time <- Sys.time()
for(i in 1:nrow(hf.rast.files)){
  
  a <- raster(hf.rast.files[[i]])
  x <- as.matrix(a)
  
  x <- 1/(x^2)  
  x <- list(x)
  
  ind <- i*1
  mat_inv[, ind] <- unlist(x, use.names = TRUE)
  
  
}  

end_time <- Sys.time()
end_time - start_time ## 16 mins


### parralelise, check if any quicker.
cores <- detectCores() 
cl <- makeCluster(10)   ##if running on personal computer this number will change
registerDoParallel(cl)
start_time <- Sys.time()
mat_inv <- foreach(i = 1:nrow(hf.rast.files)) %dopar% {
  
  a <- raster(hf.rast.files[[i]])
  x <- as.matrix(a)
  
  x <- 1/(x^2)  
  x <- list(x)
  
  ind <- i*1
  mat[, ind] <- unlist(x, use.names = TRUE)
  
  
  print(paste("written", names(a)))
  
}  

end_time <- Sys.time()
end_time - start_time ## 
stopCluster(cl)



###############
## Step 3: Set cut offs
################


mat_inv_old <- mat_inv
catchments <- mat_inv^2


#######
# Step 4: Derive the proportion of the pixel population going to each health facility
#######

rm(proportion)
proportion <- matrix(nrow = 1434125, ncol = 2801)

## separate column names
column_names <- colnames(mat_inv)
colnames(mat) <- NULL
rownames(mat) <- NULL ## slowed down the next bit having these

start_time <- Sys.time()
proportion <- t(apply(mat_inv, 1, function(x) x/sum(x)))  ## na.rm = TRUE, should break if there are NA's
end_time <- Sys.time()
end_time - start_time   ###9 mins used about 70% of mem in bld1

### do any pixels have no people going to a hf? i.e. NA in all cols.
missing_any <- proportion[rowSums(is.na(proportion)) !=ncol(proportion), ]  # 55%


###############
# Step 5: Prepare population file
###############

pop <- raster("FB_HRSL_Pop_AfricaMostly_MGMatched.2015.Annual.Data.1km.sum.tif")
mdg_hf <- raster("mainland_raster.tif")
pop <- crop(pop, mdg_hf)
#pop_matrix <- matrix(pop, nrow = 1434125, ncol = 1)  
#pop_rast <- values(pop)


#############
## Step 6: Combine TS work
#############

ts <- raster("ts_logistic.tif")
ts <- crop(ts, mdg_hf)

pop_ts <- pop*ts
pop_ts_matrix <- matrix(pop_ts, nrow= 1434125, ncol = 1)
pop_ts_matrix <- as.numeric(pop_ts_matrix)

proportion_2 <- matrix(nrow = 1434125, ncol = 2801)
proportion_2 <- proportion*pop_ts_matrix  

#can't just use colsums
population <- as.data.frame(t(apply(proportion_2, 2, function(x) sum(x, na.rm = TRUE))))
pop_vec <- t(population)

