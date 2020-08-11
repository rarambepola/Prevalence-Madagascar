####file containing  pc algorithm
# library(FSIC)
library(foreach)
library(doParallel)
library(RCITcpp)

# current_wd <- getwd()
# setwd(paste0(Sys.getenv("HOME"), "/Map_madagascar"))
# source("CausalSelection/RCIT_cpp.R")
# setwd(current_wd)

# #use partial correlations for now
# fstest <- function(x, y, alpha=0.05, N_bootstrap=100, test=TRUE){
#   cor_ts <- cor(x, y, use="complete")
#   N_obs <- length(x)
#   cor_vec <- c()
#   for(i in 1:N_bootstrap){
#     cor_vec[i] <- cor(x[sample(N_obs)], y, use="complete")
#   }
#
#   p_val <- 1 - mean(cor_ts > cor_vec)
#   # print(p_val)
#   if(test){
#     #TRUE means dependent
#     return(p_val < alpha)
#   }else
#     return(p_val)
# }


# condtest <- function(x, y, z, alpha=0.05, N_bootstrap=100, test=TRUE){
#   fitxz <- lm(x ~ z)
#   fityz <- lm(y ~ z)
#   x <- fitxz$residuals
#   y <- fityz$residuals
#   cor_ts <- cor(x, y, use="complete")
#   N_obs <- length(x)
#   cor_vec <- c()
#   for(i in 1:N_bootstrap){
#     cor_vec[i] <- cor(x[sample(N_obs)], y, use="complete")
#   }
#
#   p_val <- 1 - mean(cor_ts > cor_vec)
#   # print(p_val)
#
#   if(test){
#     #TRUE means dependent
#     return(p_val < alpha)
#   }else
#     return(p_val)
# }

# ##try RCIT
# #use partial correlations for now
# fstest <- function(x, y, alpha=0.05, test=TRUE){
#   return(RCIT(x, y)$p < alpha)
# }
#
#
# condtest <- function(x, y, z, alpha=0.05, test=TRUE){
#   return(RCIT(x, y, z)$p < alpha)
# }

##try RCIT cpp version
#use partial correlations for now
fstest <- function(x, y, alpha=0.05, test=TRUE,
                   n_rff=5){
  return(RCITcpp::RIT_wrapper(x, y, n_rff = n_rff) < alpha)
}


condtest <- function(x, y, z, alpha=0.05, test=TRUE,
                     n_rff=5, n_rffz=5){
  return(RCITcpp::RCIT_wrapper(x, y, z, n_rff = n_rff, n_rffz = n_rffz) < alpha)
}



###first version is following the description in the pcalg package documentation
##G_0 is the initial adjacency matrix - usually a complete graph but made an input
##for added flexibility
pcalg <- function(obsDat, alpha, G_0 = NULL, supprMessages = FALSE,
                  n_rff=5,
                  n_rffz=5,
                  n_cl=2){
  nvar <- length(obsDat)

  #if G_0 is not supplied, make it a complete graph
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
  }

  #run step 1 of the algorithm to get the skeleton
  if(!supprMessages) print("Starting skeleton algorithm")
  skeleton <- pcskel(G_0, obsDat, alpha, supprMessages,
                     n_rff=n_rff, n_rffz=n_rffz,
                     n_cl=n_cl)
  if(!supprMessages) print("Skeleton found")
  #skeleton matrix
  G_skel <- skeleton[[1]]
  G_out <- G_skel
  G_temp <- G_out
  #separating sets
  S <- skeleton[[2]]

  #get all unshielded triples i-j-k, for each orient them if
  #j is not in the sepset of i,k
  #actually keep track of any that should be oriented and orient them
  #afterwards if there are no contradictions
  triples <- getUnTrpls(G_skel)

  #tripleMat is a matrix for keeping track of which edges should be oriented
  tripleMat <- matrix(0, ncol=nvar, nrow=nvar)
  for(triple in triples){
    i <- triple[1]
    j <- triple[2]
    k <- triple[3]
    #if j is not in sepset orient
    if(!(j %in% S[[deparse(c(i, k))]])){
      tripleMat[i, j] <- 1
      tripleMat[k, j] <- 1
    }
  }

  #go through edges in tripleMat, if i-j is > 0 but j-i is 0 there is no disagreement
  #and we can orient
  for(i in 1:nvar){
    for(j in 1:nvar){
      if(tripleMat[i, j] > 0 & tripleMat[j, i] == 0){
        G_temp[j, i] <- 0
      }
    }
  }

  ##apply rules 1-3 repeatedly until nothing changes
  changed <- TRUE
  while(changed){
    rule1.out <- rule1(G_temp)
    rule2.out <- rule2(rule1.out[[1]])
    rule3.out <- rule3(rule2.out[[1]])
    G_temp <- rule3.out[[1]]
    changed <- any(c(rule1.out[[2]], rule2.out[[2]], rule3.out[[2]]))
  }

  ##if at any time we deleted both directions of an edge, this is wrong
  ##so put them back (unless they were not there in the original
  ##graph - these are the extra if statements inside)
  for(i in 1:nvar){
    for(j in 1:i){
      if(anyEdge(G_out, i, j) && !anyEdge(G_temp, i, j)){
        if(G_0[i, j] == 1) G_temp[i, j] <- 1
        if(G_0[j, i] == 1) G_temp[j, i] <- 1
      }
    }
  }

  G_out <- G_temp

  return(list(G_out, G_skel, S))
}


##function that orients j - k to j -> k whenever exists i -> j st i and k
##are not adjacent (otherwise this would be a v-structure)
rule1 <- function(adjmat){
  changed <- FALSE
  nvar <- dim(adjmat)[1]
  for(i in 1:nvar){
    for(j in 1:nvar){
      if(dirEdge(adjmat, i, j)){
        for(k in (1:nvar)[-i]){
          if(undirEdge(adjmat, j, k) && !anyEdge(adjmat, i, k)){
            adjmat[k,j] <- 0
            changed <- TRUE
          }
        }
      }
    }
  }
  return(list(adjmat, changed))
}

##function that orients i-j to i->j whenever there is i->k->j (otherwise there
##would be a cycle)
rule2 <- function(adjmat){
  changed <- FALSE
  nvar <- dim(adjmat)[1]
  for(i in 1:nvar){
    for(k in 1:nvar){
      if(dirEdge(adjmat, i, k)){
        for(j in (1:nvar)[-i]){
          if(dirEdge(adjmat, k, j) && undirEdge(adjmat, i, j)){
            adjmat[j,i] <- 0
            changed <- TRUE
          }
        }
      }
    }
  }
  return(list(adjmat, changed))
}

##function that orients i-j to i->j whenever there are chains i-k->j and i-l->j
##such that k, l are not adjacent (otherwise a new v-structure or cycle would be formed)
rule3 <- function(adjmat){
  changed <- FALSE
  nvar <- dim(adjmat)[1]
  for(i in 1:nvar){
    for(j in 1:nvar){
      if(undirEdge(adjmat, i, j)){
        for(k in (1:nvar)[-c(i, j)]){
          if(undirEdge(adjmat, i, k) && dirEdge(adjmat, k, j)){
            for(l in (1:nvar)[-c(i, j, k)]){
              if(undirEdge(adjmat, i, l) && dirEdge(adjmat, l,l)){
                adjmat[j, i] <- 0
                CHANGED <- TRUE
              }
            }
          }
        }
      }
    }
  }
  return(list(adjmat, changed))
}


##function that produces the skeleton (i.e. undirected acyclic graph) using order-independent
##modifications to the algorithm.
pcskel <- function(G_0, obsDat, alpha, supprMessages = FALSE,
                   n_rff=5, n_rffz=5, n_cl=2){
  #G_old is the graph we will use for the adj sets during each loop
  #G_new is the graph with edges deleted as we go
  #need both since G_1 is only updated to G_2 each time the size of the
  #conditioning set grows
  G_old <- G_0
  G_new <- G_0

  #S is the list of seperation sets
  S <- list()

  ##Step 1:Test for unconditional independence
  #get list of edges from adjacency matrix
  edges <- getEdges(G_old)


  cl<-makeCluster(n_cl)
  registerDoParallel(cl)

  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }

  z <- foreach(i = 1:length(edges), .combine = comb,
               .init = list(list(), list()),
               .packages=c('RCITcpp'),
               .export=c('fstest')) %dopar%{
                 edge <- edges[[i]]
                 if(fstest(obsDat[[edge[1]]], obsDat[[edge[2]]], alpha=alpha, test=TRUE,
                           n_rff=n_rff)){
                   edgeremove <- NULL
                   sepset <- NULL
                 }else{
                   edgeremove <- edge
                   sepset <- NA
                 }
                 list(edgeremove, sepset)
               }

  edgeremove <- z[[1]]
  sepsets <- z[[2]]
  keep <- !unlist(lapply(edgeremove, is.null))
  edgeremove <- edgeremove[keep]
  sepsets <- sepsets[keep]
  if(length(edgeremove) > 0){
    for(k in 1:length(edgeremove)){
      edge <- edgeremove[[k]]
      sepset <- sepsets[[k]]
      i <- edge[1]
      j <- edge[2]
      G_new[i, j] <- 0
      G_new[j, i] <- 0
      S[[deparse(c(i, j))]] <- NA
      S[[deparse(c(j, i))]] <- NA
    }
  }

  #update adjacency matrix (actually could have just updated G_old for this
  #whole step, done to make it conceptually similar to later)
  G_old <- G_new

  if(length(getEdges(G_old)) == 0){
    stopCluster(cl)
    return(G_old)
  }

  ##get size of the maximum conditioning set...how long does it take??
  ptm <- proc.time()
  edges <- getEdges(G_old)
  maxsize <- 0
  belowMaxCondSize <- TRUE
  for(edge in edges){
    i <- edge[1]
    j <- edge[2]
    adjSet <- getAdj(G_old, i, j)
    maxsize <- max(maxsize, length(adjSet))
  }

  ##Step 2: Test for conditional independence
  condSetSize <- 1
  edges <- getEdges(G_old)
  #for each edge, check if the adjacency set is big enough (both directions)
  #to condition on
  #variable maxCondSize checks if the biggest possible size has been reached for
  #the conditioning sets
  belowMaxCondSize <- TRUE
  while(belowMaxCondSize){
    belowMaxCondSize <- FALSE
    print(paste0("condSetSize ", condSetSize))
    #go through each edge


    z <- foreach(i = 1:length(edges), .combine = comb,
                 .init = list(list(), list()),
                 .packages=c('RCITcpp'),
                 .export=c('getAdj', 'getSubsets', 'condtest')) %dopar%{
                   edge <- edges[[i]]
                   if(!supprMessages) print(paste0("testing edge ", edge[1], ",", edge[2], " with size ", condSetSize))
                   i <- edge[1]
                   j <- edge[2]
                   brokenEdge <- FALSE
                   #do the i, j direction
                   #get set adj(G,i)\j
                   adjSet <- getAdj(G_old, i, j)

                   #if the adjacency set is bigger than conditioning set size, test
                   #conditional independence
                   if(length(adjSet) >= condSetSize){
                     belowMaxCondSize <- TRUE
                     #test independence conditioning on  all subsets of right size
                     condSets <- getSubsets(adjSet, condSetSize)

                     for(condSet in condSets){
                       #make a matrix of conditioning data
                       condData <- do.call(cbind, obsDat[condSet])
                       #if they are conditionally independent, remove edge, add condSet to
                       #separation set list and stop trying this edge
                       if(!condtest(obsDat[[i]], obsDat[[j]], condData, alpha=alpha,
                                    n_rff=n_rff, n_rffz=n_rffz)){
                         edgeremove <- edge
                         sepset <- condSet
                         brokenEdge <- TRUE
                         break
                       }
                     }
                   }

                   if(!brokenEdge){
                     #do j, i direction
                     i <- edge[1]
                     j <- edge[2]

                     #do the i, j direction
                     #get set adj(G,i)\j
                     adjSet <- getAdj(G_old, j, i)

                     #if the adjacency set is bigger than conditioning set size, test
                     #conditional independence
                     if(length(adjSet) >= condSetSize){
                       #test independence conditioning on  all subsets of right size
                       condSets <- getSubsets(adjSet, condSetSize)

                       for(condSet in condSets){
                         #make a matrix of conditioning data
                         condData <- do.call(cbind, obsDat[condSet])
                         #if they are conditionally independent, remove edge, add condSet to
                         #separation set list and stop trying this edge
                         if(!condtest(obsDat[[i]], obsDat[[j]], condData, alpha=alpha)){
                           edgeremove <- edge
                           sepset <- condSet
                           brokenEdge <- TRUE
                           break
                         }
                       }
                     }
                   }
                   if(!brokenEdge){
                     edgeremove <- NULL
                     sepset <- NULL
                   }
                   list(edgeremove, sepset)
                 }
    #update the adjacency matrix with edges remove

    edgeremove <- z[[1]]
    # print(length(z[[1]]))
    sepsets <- z[[2]]
    # print(length(z[[2]]))
    # print(edgeremove)s
    keep <- !unlist(lapply(edgeremove, is.null))
    edgeremove <- edgeremove[keep]
    sepsets <- sepsets[keep]
    if(length(edgeremove) > 0){
      for(k in 1:length(edgeremove)){
        edge <- edgeremove[[k]]
        sepset <- sepsets[[k]]
        i <- edge[1]
        j <- edge[2]
        G_new[i, j] <- 0
        G_new[j, i] <- 0
        S[[deparse(c(i, j))]] <- sepset
        S[[deparse(c(j, i))]] <- sepset
      }
    }
    G_old <- G_new
    condSetSize <- condSetSize + 1
    edges <- getEdges(G_old)
    maxsize <- 0
    for(edge in edges){
      i <- edge[1]
      j <- edge[2]
      adjSet <- getAdj(G_old, i, j)
      maxsize <- max(maxsize, length(adjSet))
    }
    print(paste0("maxsize ", maxsize))
    if(maxsize >= condSetSize){
      belowMaxCondSize <- TRUE
    }
  }
  stopCluster(cl)
  return(list(G_old, S))
}



##function to get edges from (undirected) adjacency matrix
getEdges <- function(adjMat){
  n <- dim(adjMat)[1]
  edges <- list()
  for(i in 1:n){
    for(j in 1:n){
      if(adjMat[i, j] == 1){
        edges[[length(edges) + 1]] <- sort(c(i, j))
      }
    }
  }
  return(unique(edges))
}

##function to get set Adj(G, i)\j
getAdj <- function(adjMat, i, j){
  n <- dim(adjMat)[1]
  adjs <- c()
  for(k in (1:n)[-j]){
    if(adjMat[k, i] == 1){
      adjs <- c(adjs, k)
    }
  }
  return(adjs)
}

##function to make complete adjacency matrix
compltMat <- function(n){
  return(matrix(1, nrow=n, ncol=n) - diag(n))
}

##function that takes a set and a number and returns all
##subsets of that size
getSubsets <- function(set, nSubset){
  sslist <- list()
  if(nSubset == 1){
    for(i in 1:length(set)){
      sslist[[i]] <- c(set[i])
    }
  }else{
    for(i in 1:(length(set) - nSubset + 1)){
      sslist2 <- getSubsets(set[-(1:i)], nSubset - 1)
      for(j in 1:length(sslist2)){
        sslist[[length(sslist) + 1]] <- c(set[i], sslist2[[j]])
      }
    }
  }
  return(sslist)
}


##function that takes an adjacency matrix and returns
##unshielded triples
getUnTrpls <- function(adjmat){
  triples <- list()
  nvars <- dim(adjmat)[1]
  for(i in 1:nvars){
    for(j in 1:nvars){
      if(adjmat[i, j] == 1){
        for(k in (1:nvars)[-i]){
          if(adjmat[j, k] == 1){
            if(adjmat[i, k] == 0){
              if(i < k){
                triples[[length(triples) + 1]] <- c(i, j, k)
              }else{
                triples[[length(triples) + 1]] <- c(k, j, i)
              }
            }
          }
        }
      }
    }
  }
  triples <- unique(triples)
  return(triples)
}

##functions that returns boolean if:
##there is an edge i -> j
dirEdge <- function(adjmat, i, j){
  return(adjmat[i, j] == 1 && adjmat[j, i] == 0)
}
##there is an edge i <->j
undirEdge <- function(adjmat, i, j){
  return(adjmat[i, j] == 1 && adjmat[j, i] == 1)
}
##there is an edge i - j
anyEdge <- function(adjmat, i, j){
  return(adjmat[i, j] == 1 || adjmat[j, i] == 1)
}


if(FALSE){
  M <- 100
  x <- runif(M)
  y <- x + rnorm(M, sd = 0.1)
  z <- y + rnorm(M, sd = 0.1)

  print(pcalg(list(x,y,z), alpha = 0.01, compltMat(3)))

  print(pcalg(list(x,y,z), alpha = 0.01))

}
