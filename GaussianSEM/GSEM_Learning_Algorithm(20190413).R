

##### GSEM_simulation_fun: Simulations for recovering GSEM ####
GSEM_simulation_fun = function(seed, n_real, p_real, d = 1, max_degree = 1, direction = "forward", alpha = 0.001, path = F ){
  #This algorithm is for simulations using Gaussian SEM Learning_algorithm. 
  
  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################
  
  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$x[1:n_real,1:p_real]
  
  result_GSEM = GSEM_Algorithm(data, alpha = alpha, direction = direction, graph = graph, max_degree = max_degree)
  
  return(result_GSEM)
}

#####GSEM Learning Algorithm via Conditional Independence Test#########
GSEM_Algorithm = function(data, alpha = 0.05, direction ="forward", graph = NULL, max_degree = 1){
  
  library(bnlearn)
  ###################
  X = as.matrix( data )
  p = ncol(X)
  n = nrow(X)
  RemNode = 1:p 
  pi_GSEM = NULL
  Estimated_G = Estimated_O = matrix(0, p ,p)
  evaluation_result_GSEM = evaluation_result_GSEM_MEC = NULL
  ####
  Runtime = proc.time()[3]
  
  #### Step 1): Finding the Ordering ####
  if(direction == "forward"){
    Ordering = Forward_Learning_fun(X, max_degree = max_degree)
  }else if(direction =="backward"){
    Ordering = Backward_Learning_fun(X, max_degree = max_degree)
  }
  #### Step 2): Finding the Parents ####
  
  if(max_degree == 1){
    for(m in 2:p){
      Rho = cov2cor( solve( cov( X[,Ordering[1:m]] ) ) )[1:(m-1),m]
      Z_stat = abs( 1/2 * log( (1+Rho)/(1-Rho) ) * sqrt( n - 1 - m ) )
      parent_pvalue = pnorm(Z_stat,lower.tail = F) * 2
      Estimated_G[Ordering[m], Ordering[1:(m-1)][parent_pvalue < alpha]] = 1
    }
  }else if(max_degree == 2){
    for(m in 2:p){
      j = Ordering[m]
      S = Ordering[ 1:(m-1)]
      Z = NULL;
      for(l in S ){
        Z = cbind(Z, poly( X[,l], degree = max_degree, raw = T) )
      }
      Rho = cov2cor( solve( cov( cbind( Z, X[,j] ) ) ) )[-c(2*m-1),2*m-1]
      Z_stat = abs( 1/2 * log( (1+Rho)/(1-Rho) ) * sqrt( n - 2 - 2*m ) )
      parent_pvalue = pnorm(Z_stat,lower.tail = F) *2
      Estimated_G[Ordering[m], Ordering[1:(m-1)][unique( (which(parent_pvalue < alpha) + 1)  %/% 2 )]] = 1
    }
  }else{
    for(m in 2:p){
      j = Ordering[m]
      S = Ordering[ 1:(m-1)]
      Z = NULL;
      for(l in S ){
        Z = cbind(Z, poly( X[,l], degree = max_degree, raw = T) )
      }
      lm0 = summary( lm( X[,j] ~ Z ) )$coef
      parent_id = unique( (which( lm0[-1,4] < alpha ) + 1) %/% max_degree )
      S[parent_id]
      Estimated_G[j, S[parent_id]] = 1
    }
  }
  
  #library(corrplot)
  #corrplot(Estimated_G[Ordering, Ordering])
  
  ####
  Runtime = proc.time()[3] - Runtime
  
  print(paste("It takes: ", Runtime))
  
  ####
  est_MEC = dag2cpdagAdj(Estimated_G)
  
  
  ####
  if( !is.null(graph) ){
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
    evaluation_result_GSEM2 = NULL
  }
  
  return(
    list( DAG_Evaluation = evaluation_result_GSEM, 
          MEC_Evaluation = evaluation_result_GSEM_MEC,  
          Oracle = evaluation_result_GSEM2, 
          DAG = Estimated_G, 
          Ordering = Ordering, 
          Time = Runtime)
  )
}

#################################################

Forward_Learning_fun = function(X, max_degree = 1){
  #### Step 1): Finding the Ordering ####
  p = ncol(X)
  n = nrow(X)
  Ordering = rep(0, p)
  RemNodes = 1:p
  k = 1
  Sig = var(X)
  Ordering[k] = which.min( diag(Sig) )
  RemNodes = setdiff(RemNodes, Ordering)
  while( length(RemNodes) > 1 ){
    
    if( max_degree == 1){
      ###
      k = k + 1
      Ome = solve( Sig[Ordering[1:(k-1)],Ordering[1:(k-1)]] )
      scores = sapply(RemNodes, function(j) ( Sig[j,j] - Sig[j,Ordering[1:(k-1)]] %*% Ome %*% Sig[Ordering[1:(k-1)],j] ) )
      #sapply(RemNodes, function(j) 1/( var(X[,j]) - cov(X[,j], X[,Ordering[1:(k-1)]]) %*% Ome %*% cov(X[,Ordering[1:(k-1)]],X[,j])  ) )
      #scores = sapply(RemNodes, function(j) sum(lm(X[,j]~X[,Ordering[1:(k-1)]])$resid^2)/n )
      
    }else{
      k = k + 1
      Z = NULL;
      for(l in Ordering[1:(k-1)] ){
        Z = cbind(Z, poly( X[,l], degree = max_degree, raw = T) )
      }
      ###
      #scores = sapply(RemNodes, function(j) mean(lm(X[,j]~Z)$resid^2))
      scores = sapply(RemNodes, function(j) 1/solve(cov( cbind(X[,j], Z)))[1,1] )
    }
    
    #print(scores)
    Ordering_j = which( scores == min(scores) )
    if( length(Ordering_j) > 1){
      Ordering[k] = RemNodes[sample(Ordering_j, 1)]
    }else{
      Ordering[k] = RemNodes[Ordering_j]
    }
    RemNodes = setdiff(RemNodes, Ordering)
  }
  Ordering[p] = RemNodes
  return( Ordering )
}

Backward_Learning_fun = function(X, max_degree = 1){
  #### Step 1): Finding the Ordering ####
  p = ncol(X)
  n = nrow(X)
  Ordering = rep(0, p)
  RemNodes = 1:p
  k = p
  while( length(RemNodes) > 1 ){
    if(max_degree ==1){
      scores = 1/diag(solve(cov( X[,RemNodes])))
    }else{
      scores = NULL
      for(j in RemNodes){
        Z = NULL;
        for(l in setdiff(RemNodes,j) ){
          Z = cbind(Z, poly( X[,l], degree = max_degree, raw = T) )
        }
        scores = c( scores, 1/solve(cov( cbind(X[,j], Z)))[1,1] )
        #scores = c( scores, mean( lm(X[,j]~ Z)$resid^2) )
      }
    }
    #
    #scores = sapply(RemNodes, function(j) sum(lm(X[,j]~X[,setdiff(RemNodes,j)])$resid^2)/n )
    #print(scores)
    Ordering_j = which( scores == max(scores) )
    if( length(Ordering_j) > 1){
      Ordering[k] = RemNodes[sample(Ordering_j, 1)]
    }else{
      Ordering[k] = RemNodes[Ordering_j]
    }
    k = k - 1
    RemNodes = setdiff(RemNodes, Ordering)
  }
  Ordering[1] = RemNodes
  return( Ordering )
}

##################################
#####GSEM Learning Algorithm via Conditional Independence Test#########
GSEM_Algorithm_True = function(Sigma, direction ="forward", graph = NULL){
  
  library(bnlearn)
  p = ncol(Sigma)
  RemNode = 1:p 
  pi_GSEM = NULL
  Estimated_G = Estimated_O = matrix(0, p ,p)
  evaluation_result_GSEM = evaluation_result_GSEM_MEC = NULL
  ####
  Runtime = proc.time()[3]
  
  #### Step 1): Finding the Ordering ####
  if(direction == "forward"){
    Ordering = rep(0, p)
    RemNodes = 1:p
    k = 1
    Sig = Sigma
    Ordering[k] = which.min( diag(Sig) )
    RemNodes = setdiff(RemNodes, Ordering)
    while( length(RemNodes) > 1 ){
      ###
      k = k + 1
      Ome = solve( Sig[Ordering[1:(k-1)],Ordering[1:(k-1)]] )
      scores = sapply(RemNodes, function(j) ( Sig[j,j] - Sig[j,Ordering[1:(k-1)]] %*% Ome %*% Sig[Ordering[1:(k-1)],j] ) )
      #print(scores)
      Ordering_j = which( scores == min(scores) )
      if( length(Ordering_j) > 1){
        Ordering[k] = RemNodes[sample(Ordering_j, 1)]
      }else{
        Ordering[k] = RemNodes[Ordering_j]
      }
      RemNodes = setdiff(RemNodes, Ordering)
    }
    Ordering[p] = RemNodes
  }else if(direction =="backward"){
    Sig = Sigma
    Ordering = rep(0, p)
    RemNodes = 1:p
    k = p
    while( length(RemNodes) > 1 ){
      scores = 1/diag(solve( Sig[RemNodes,RemNodes]))
      Ordering_j = which( scores == max(scores) )
      if( length(Ordering_j) > 1){
        Ordering[k] = RemNodes[sample(Ordering_j, 1)]
      }else{
        Ordering[k] = RemNodes[Ordering_j]
      }
      k = k - 1
      RemNodes = setdiff(RemNodes, Ordering)
    }
    Ordering[1] = RemNodes
  }
  #print( proc.time()[3]  - Runtime )
  
  #### Step 2): Finding the Parents ####
  B = graph
  B[B!=0] =1
  MEC = B + t(B)
  S = NULL
  for( i in Ordering){
    S <- union(S, i)
    MEC[i,S] <- 0 
  }
  Estimated_G = t(MEC)
  
  #library(corrplot)
  #corrplot(Estimated_G[Ordering, Ordering])
  
  ####
  Runtime = proc.time()[3] - Runtime
  print(paste("Oracle takes: ", Runtime))
  ####
  est_MEC = dag2cpdagAdj(Estimated_G)
  ####
  if( !is.null(graph) ){
    B = graph
    B[B!=0] =1
    MEC = B + t(B)
    evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
    evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)
    evaluation_result_GSEM2 = NULL
  }
  
  return(
    list( DAG_Evaluation = evaluation_result_GSEM, 
          MEC_Evaluation = evaluation_result_GSEM_MEC,  
          Oracle = evaluation_result_GSEM2, 
          DAG = Estimated_G, 
          Ordering = Ordering, 
          Time = Runtime)
  )
}




##################################





#####GSEM Learning Algorithm via L1 regularized Regression #######################
GSEM_Algorithm2 = function(data, graph = NULL, lambda.lasso = 0.4, sparsity_level = 0.18){
  
  library(glmnet)
  tuning.parameter <- lambda.lasso
  
  ###################
  X = as.matrix( data )
  RemNode = 1:p 
  pi_GSEM = NULL
  ScoreMatrix = matrix(0, p, p)
  
  #### Finding Ordering ####
  Time_GSEM = proc.time()[3]
  
  #### ordering estimation ####
  Ordering = rep(0, p)
  Estimated_G = Estimated_O = matrix(0, p ,p)
  RemNodes = 1:p
  k = 1
  Ordering[k] = which.min( sapply(RemNodes, function(j) var(X[,j]) ) )
  RemNodes = setdiff(RemNodes, Ordering)
  while( length(RemNodes) >1 ){
    k = k + 1
    Ordering[k] = RemNodes[which.min( sapply(RemNodes, function(j) sum(lm(X[,j]~X[,Ordering[1:(k-1)]])$resid^2)/n ) )]
    RemNodes = setdiff(RemNodes, Ordering)
  }
  Ordering[p] = RemNodes

  #### paremnts estimation
  glm0 = glmnet(X[,rep(Ordering[1],2)], X[,Ordering[2]], family = "gaussian", alpha = 1, lambda = tuning.parameter )$beta
  if( sum( glm0[,1] != 0)   ){
    Estimated_G[Ordering[2], Ordering[1]] = 1
  }
  for(k in 3:p){
    glm0 = glmnet(X[, Ordering[1:(k-1)]], X[,Ordering[k]], family = "gaussian", alpha = 1, lambda = tuning.parameter )$beta
    if( sum( glm0[,1] != 0)   ){
      parents = which( glm0[,1] != 0 )
      Estimated_G[Ordering[k], Ordering[parents]] = 1
    }
  }
  
  if( !is.null(graph) ){
    B = graph
    B[B!=0] =1
    MEC = B + t(B)
    Oracle_DAG = estimated_graph_fun( MEC, Ordering )
    evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
    evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), dag2cpdagAdj(Estimated_G) )
    evaluation_result_GSEM_Oracle = evaluation_fun( B, Oracle_DAG ) 
  }
  
  return( list( DAG_Evaluation = evaluation_result_GSEM, 
                MEC_Evaluation = evaluation_result_GSEM_MEC,  
                Oracle_Evaluation = evaluation_result_GSEM_Oracle, 
                DAG = Estimated_G, 
                Ordering = Ordering , 
                Time = Time_GSEM) )
}  

####################



#RawData = GSEM_generator( n, p, d, noiseVar, beta_min, beta_max, graph_type, seed = 1)

#GSEM_Algorithm(data, graph, lambda.lasso = 0.4, sparsity_level = 1)$DAG_Evaluation



