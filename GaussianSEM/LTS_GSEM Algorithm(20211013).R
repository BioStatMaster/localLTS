
###################
LTS_GSEM_simulation_fun = function(seed, n_real, p_real, d = 1, h_ratio = 0.5, thresh = qnorm(0.975), path = NULL ){
  #This algorithm is for simulations using Gaussian SEM Learning_algorithm. 
  
  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################
  
  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$x[1:n_real,1:p_real]
  result_GSEM = LTS_GSEM_Algorithm(data, lambda1 = 0.1^4, lambda2 = 0.5*sqrt(log(p_real)/(n_real*h_ratio)), alpha = h_ratio, thresh = thresh, graph = graph)
  return(result_GSEM)
}

#### High Dimensional Linaer SEM Learning using Lasso #################################
LTS_GSEM_Algorithm = function(data, lambda1 = 0.1^4, lambda2 = 0.1^3, alpha = 0.75, thresh = qnorm(0.975), graph = NULL){
  
  # lambda : penalty parameter
  # alpha : fraction of trimmed observations
  # thresh : truncation parameter for conditional variance
  
  library(robustHD)
  library(Rcpp)
robustHD::sparseLTS
methods("sparseLTS")
robustHD:::sparseLTS.formula
  ###################
  X = as.matrix(data)
  p = ncol(X)
  n = nrow(X)
  Ordering = NULL
  Estimated_G = matrix(0, p ,p)
  
  ####
  Runtime = proc.time()[3]
  
  ### Step (1): Ordering Estimation 
  Ordering = LTS_Ordering_Fun(X, lambda = lambda1, alpha, thresh = thresh)
  
  ### Step (2): Parents Estimation 
  for(m in 2:p){
    S = Ordering[1:(m-1)]
    LTS_Est = sparseLTS(x = X[,S], y = X[,Ordering[m]], lambda = lambda2, mode = "lambda", alpha, intercept = F)
    Estimated_G[Ordering[m], S[which(LTS_Est$coefficients != 0)]] = 1
  }
  Estimated_G
  
  ### Step (1): Inverse Covariance Estimation ###
  Runtime = proc.time()[3]- Runtime
  ####
  print(paste("LTS_GSEM takes: ", Runtime))
  ####
  
  evaluation_result_GSEM2 = NULL
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
    evaluation_result_GSEM2 = evaluation_fun( B, t( MEC) )
  }else{
    evaluation_result_GSEM = NULL
    evaluation_result_GSEM_MEC = NULL
    evaluation_result_GSEM2 = NULL
  }
  #evaluation_result_GSEM
  
  return(
    list( DAG_Evaluation = evaluation_result_GSEM, 
          MEC_Evaluation = evaluation_result_GSEM_MEC,  
          Oracle_Evaluation = evaluation_result_GSEM2, 
          DAG = Estimated_G,
          Ordering = Ordering, 
          Time = Runtime)
  )
}


LTS_Ordering_Fun = function(X, lambda = 0.1^4, alpha = 0.5, thresh = qnorm(0.975) ){
  p = ncol(X)
  Ordering = NULL;
  while(length(Ordering) < p-1 ){
    length(Ordering)
    RemNodes = setdiff(1:p, Ordering)
    Scores = sapply(RemNodes, function(j) Score_Fun(X = X, RemNodes = RemNodes, j, lambda = lambda, alpha = alpha, thresh = thresh ) )
    Scores
    Ordering = c(RemNodes[which.max(Scores)], Ordering)
  }
  Ordering = c(RemNodes[which.min(Scores)], Ordering)
  return(Ordering)
}

Score_Fun = function(X, RemNodes, j, lambda = 0.1^4, alpha = 0.75, thresh = 2.5){
  id1 = 1:nrow(X)
  S = setdiff(RemNodes, j)
  LTS_Est = sparseLTS(x = X[,S], y = X[,j], lambda = lambda, mode = "lambda", alpha = alpha, intercept = F)
  if( length(S) > 1){
    lm1 = sparseLTS(x = X[,S][,which(LTS_Est$coefficients != 0)], y = X[,j], lambda = 0, mode = "lambda", alpha = alpha)
  }else{
    if( LTS_Est$coefficients != 0 ){
      lm1 = sparseLTS(x = as.matrix(X[,S]), y = X[,j], lambda = 0, mode = "lambda", alpha = alpha)
    }else{
      lm1 = sparseLTS(x = rep(0,nrow(X) ), y = X[,j], lambda = 0, mode = "lambda", alpha = alpha, intercept = T)
    }
  }
  ResidSet = abs(lm1$residuals)
  
  Id_TrVar = which(ResidSet < thresh)
  if( length(Id_TrVar) > 0 ){
    var1 = mean(ResidSet[Id_TrVar]^2)
  }else{
    var1 = mean(ResidSet^2)
  }
  #var1 = mean(ResidSet[id1]^2)
  #sd1 = getScale(lm1) 
  #resid_id = which(ResidSet/sd1 > thresh)
  # while( length( resid_id ) > 0  ){
  #    id1 = setdiff(id1, resid_id)
  #    var1 = mean(ResidSet[id1]^2)
  #    resid_id = id1[which(ResidSet[id1]/sqrt(var1) > thresh)]
  # }
  return( var1 )
}





