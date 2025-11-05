# Copyright (c) 2018 - 2019  Wenyu Chen [wenyuc@uw.edu]
# All rights reserved.  See the file COPYING for license terms.



###############
### Main method with top-down approach
###############
#' Estimate topological ordering and DAG using high dimensional top-down approach
#'
#' @param X An n-by-p data matrix.
#' @param J maximum number of parents to condition on.
#' @return Estimated Adjacency matrix and topological ordering.
#' @examples
#' X1<-rnorm(100)
#' X2<-X1+rnorm(100)
#' EqVarDAG_HD_TD(cbind(X1,X2),2)
#'
#' #$adj
#' #[,1] [,2]
#' #[1,]    0    1
#' #[2,]    0    0
#' #
#' #$TO
#' #[1] 1 2


###################
TD_simulation_fun = function(seed, n_real, p_real, d = 1, path = F ){
  #This algorithm is for simulations using Gaussian SEM Learning_algorithm.

  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################

  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$x[1:n_real,1:p_real]
  B = graph; B[B!=0] = 1
  result_GSEM = TD(data, q = d, sparsity_level = 3, graph = graph)
  return(result_GSEM)
}


TD<-function(data, q, sparsity_level = 3, graph = NULL){
  library(leaps)
  J = q
  X = as.matrix(data)
  TD_res = EqVarDAG_HD_TD(X,J,sparsity_level)
  Estimated_G = t(TD_res$adj)
  Ordering = TD_res$TO
  Runtime = TD_res$Runtime
  print(paste("TD takes:", Runtime))
  if( !is.null(graph) ){
    est_MEC = dag2cpdagAdj(Estimated_G)
    B = graph
    B[B!=0] =1
    MEC = B + t(B)
    evaluation_result_GSEM = evaluation_fun( B, Estimated_G )
    evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)

    S = NULL
    for( i in Ordering){
      S <- union(S, i)
      MEC[i,S] <- 0
    }
    evaluation_result_GSEM2 = evaluation_fun( B, t(MEC) )
  }else{
    evaluation_result_GSEM = NULL
    evaluation_result_GSEM_MEC = NULL
  }

  return(
    list( DAG_Evaluation = evaluation_result_GSEM,
          MEC_Evaluation = evaluation_result_GSEM_MEC,
          Oracle_Evaluation = evaluation_result_GSEM2,
          DAG = Estimated_G,
          Ordering = Ordering,
          Time = Runtime)
  )
}

EqVarDAG_HD_TD<-function(X,J, sparsity_level =3){
  # Input
  # X : n by p matrix of data
  # J : maximum number of parents to condition on
  # Output
  # adj: estimated adjacency matrix
  # TO : estimated topological ordering
  Runtime = proc.time()[3]
  n<-dim(X)[1]
  p<-dim(X)[2]
  rr<-rev(getOrdering(X,J)) # use J=5
  result<-matrix(0,p,p)

  # for (ii in 1:(p-1)){
  #   now<-rr[ii]
  #   this<-sort(rr[(ii+1):p])
  #   if (length(this)>1){
  #     # variable selection
  #     if (n>100){
  #       lassom<-glmnet::cv.glmnet(X[,this],X[,now]  )
  #       bfit<-coefficients(lassom)[-1]
  #     } else {
  #       lassom<-glmnet::glmnet(X[,this],X[,now] )
  #       bic<-n*log(colSums((predict(lassom,
  #                                   X[,this])-X[,now])^2)/n)+lassom$df*log(n)+
  #         2*lassom$df*log(p-ii)
  #       bfit<-coefficients(lassom)[,which(bic==min(bic))[1]][-1]
  #     }
  #     for (jj in 1:length(this)){
  #       if(bfit[jj]!=0)
  #         result[this[jj],now]<-1
  #     }
  #   } else {
  #     # deal with the last two nodes
  #     lmod<-summary(RcppEigen::fastLm(X[,now]~X[,this]))
  #     if (lmod$coef[2,4]<0.05) {
  #       result[this,now]<-1
  #     }
  #   }
  # }
  # Runtime = proc.time()[3] - Runtime

  #### I fixed (G.Park)####
  library(glmnet)
  Ordering = rev(rr)
  Estimated_G = Estimated_O = matrix(0, p ,p)
  RemNodes = 1:p
  for(i in p:2){
    Cand_Pa = setdiff(1:p,Ordering[i:p])
    if( length(Cand_Pa) > 1){
      glm_beta = glmnet(X[ ,Cand_Pa], X[,Ordering[i]],
                          family = "gaussian",alpha = 1,
                          lambda = sparsity_level * sqrt( log(p)/n ) )$beta[,1]
      Estimated_G[Ordering[i], Cand_Pa[ which(glm_beta !=0 ) ] ] = 1
    }else{
      glm_beta = glmnet(X[ ,rep(Cand_Pa,2)], X[,Ordering[i]], family = "gaussian",alpha = 1,
                          lambda = sparsity_level * sqrt( log(p)/n ) )$beta
      if( any( glm_beta !=0 ) ){
        Estimated_G[Ordering[i], Cand_Pa ] = 1
      }
    }
  }
  result = t(Estimated_G)
  Runtime = proc.time()[3] - Runtime
  return(list(adj=result,TO=rev(rr), Runtime = Runtime))
}

###############
### helper functions
###############
# compute best subset search
helper.func <- function(z, Y, Theta, J,mtd="exhaustive"){
  leaps::regsubsets(x = Y[, Theta, drop = F],
                    y = Y[, z, drop = F],
                    method=mtd,
                    nbest = 1,
                    nvmax = min(J, sum(Theta > 0 )), really.big = T)
}

# compute topological ordering
getOrdering <- function(Y, J){
  p <- dim(Y)[2]
  variances <- apply(Y, MAR = 2, sd)
  Theta <- rep(0, p)
  Theta[1] <- which.min(variances)
  out <- sapply(setdiff(1:p, Theta[1]),
                  function(z){
                    sum(resid(RcppEigen
                              ::fastLm(Y[, z] ~ Y[, Theta[1], drop = F]) )^2)})
  Theta[2] <- setdiff(1:p, Theta[1])[which.min(out)]
  for(i in 3:p){
    out <- lapply(setdiff(1:p, Theta),
                  function(jj)helper.func(jj, Y, Theta[seq(i-1)], J))
    nextRoot <- which.min(sapply(out,function(x){min(x$rss)}))
    Theta[i] <- setdiff(1:p, Theta)[nextRoot]
  }
  return(Theta)
}
