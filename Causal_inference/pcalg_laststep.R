library(foreach)
library(doParallel)
library(RCITcpp)


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



pcalg.last <- function(obsDat, alpha, last.index, G_0 = NULL, supprMessages = FALSE,
                       n_rff=5,
                       n_rffz=5,
                       n_cl=2){
  nvar <- length(obsDat)
  #if G_0 is not supplied, make it a complete graph
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
    warnings("G_0 not supplied")
  }

  #run step 1 of the algorithm to get the skeleton
  if(!supprMessages) print("Starting skeleton algorithm")
  skeleton <- pcskel.last(G_0, obsDat, alpha, last.index, supprMessages,
                          n_rff,
                          n_rffz,
                          n_cl)
  if(!supprMessages) print("Skeleton found")
  #skeleton matrix
  #print(skeleton)

  return(list(skeleton[[1]]))
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
pcskel.last <- function(G_0, obsDat, alpha, last.index, supprMessages = FALSE,
                        n_rff=5,
                        n_rffz=5,
                        n_cl=2){
  #G_old is the graph we will use for the adj sets during each loop
  #G_new is the graph with edges deleted as we go
  #need both since G_1 is only updated to G_2 each time the size of the
  #conditioning set grows
  G_old <- G_0
  G_new <- G_0

  #S is the list of seperation sets
  S <- list()

  ##Step 1:Test for unconditional independence
  #get list of edges only connected to last element

  getEdgesLast <- function(G_0, last.index){
    adj.nodes <- unique(c(which(G_0[last.index, ] > 0), which(G_0[, last.index] > 0)))
    #print(G_0)
    #print(adj.nodes)
    edges <- list()
    if(length(adj.nodes) == 0) return(NULL)
    for(i in 1:length(adj.nodes)){
      edges[[i]] <- c(last.index, adj.nodes[i])
    }
    return(edges)
  }

  edges <- getEdgesLast(G_0, last.index)

  #print(edges)

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
  #print(edgeremove)
  sepsets <- z[[2]]
  keep <- !unlist(lapply(edgeremove, is.null))
  edgeremove <- edgeremove[keep]
  sepsets <- sepsets[keep]
  #print(edgeremove)
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

  ##get size of the maximum conditioning set...how long does it take??
  ptm <- proc.time()
  edges <- getEdgesLast(G_old, last.index)
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
  edges <- getEdgesLast(G_old, last.index)
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
                                    n_rff=n_rff,
                                    n_rffz=n_rffz)){
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
    #print("z[[1]]")
    #print(z[[1]])
    sepsets <- z[[2]]
    keep <- !unlist(lapply(edgeremove, is.null))
    edgeremove <- edgeremove[keep]
    sepsets <- sepsets[keep]
    #print(edgeremove)
    if(length(edgeremove) > 0){
      #print(edgeremove)
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
    edges <- getEdgesLast(G_old, last.index)
    if(!is.null(edges)){
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
