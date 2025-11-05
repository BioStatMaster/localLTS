LISTEN_simulation_fun = function(seed, n_real, p_real, d = 1, reg_param = 0.001, hard_thresholding = 0.20, path = F ){
  setwd("D:/Dropbox/GroupMeetings/Algorithm(Rcode)")
  source("ComparisonAlg/LISTEN(ghoash)/LISTEN_Algorithm(20190515).R")
  
  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################
  
  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$x[1:n_real,1:p_real]
  
  result_GSEM = LISTEN_Algorithm(X = data, reg_param= reg_param, hard_thresholding = hard_thresholding, graph = graph)
  
  return(result_GSEM)
}

LISTEN_Algorithm = function(X, reg_param= 0.001, hard_thresholding = 0.20, graph = NULL){
  
  #### Start ###
  Runtime = proc.time()[3]
  result = LISTEN(X, reg_param = reg_param, hard_thresholding = hard_thresholding)
  Runtime = proc.time()[3] - Runtime
  #print(paste("It takes: ", Runtime))
  
  Estimated_G = result$DAG
  est_MEC = dag2cpdagAdj(Estimated_G)
  Ordering = result$Ordering
  evaluation_result_GSEM = evaluation_result_GSEM_MEC = NULL
  evaluation_result_GSEM2 = NULL
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
  }
  return( 
    list( DAG_Evaluation = evaluation_result_GSEM, 
          MEC_Evaluation = evaluation_result_GSEM_MEC,  
          Oracle_Evaluation = evaluation_result_GSEM2, 
          DAG = Estimated_G, 
          Ordering = Ordering, 
          Time = Runtime
    ) 
  )
}



LISTEN = function(X, reg_param = 0.001, hard_thresholding = 0.20){
  #   """
  # Learn the causal ordering, given the empirical
  # covariance marix and the estimated precision matrix.
  # :param Sigma: Empirical covariance matrix.
  # :param O_hat: Estimated precision matrix.
  # :param hard_thresholding: value below which entries in
  # the inverse covariance matrix are assumed to be 0.
  # :return:
  Sigma = cov(X)
  n = nrow(Sigma)
  z = rep(0,n)       # DAG order.
  vars = rep(0,n)    # variances.
  A = matrix(0, n, n) # DAG
  r = rep(0,n)                   # ratios.
  Rem_Node = 1:n
  count = n
  
  ## clime in high otherwise, covariance matrix ##
  O_hat = inv_cov_estimate_clime(X, threshold = reg_param)
  O = O_hat # light copy
  O_update = O # light copy2
  
  while( count >= 2){ #count check.
      # print "============================="
      # print r
      # print "============================="
    r = diag(O_update)
    tv = which.min(r) # This is the node label.
    z[count] = tv #
    r[tv] = Inf
    Rem_Node = setdiff(Rem_Node, tv)
    vars[tv] = Sigma[tv, tv]
    if(O[tv, tv] > 0){ 
      vars[tv] = 1/O[tv, tv] ### dimension 안맞음.
    }
    
    # Delete terminal vertex from precision matrix.
    v = O[tv,Rem_Node]
    if( O[tv, tv] > 0 && length(Rem_Node)>1 ){
      O[Rem_Node,Rem_Node] = inv_cov_estimate_clime(X[,Rem_Node], threshold = reg_param)
      #print(O[Rem_Node,Rem_Node] )
      #print(solve(cov(X[,Rem_Node])))
      #O[Rem_Node,Rem_Node] = O[Rem_Node,Rem_Node] - outer(v,v)/O[tv, tv]
      O[tv,] = O[,tv] = 0
    }else if(length(Rem_Node) == 1){
      O[Rem_Node,Rem_Node] = O[Rem_Node,Rem_Node] - outer(v,v)/O[tv, tv]
      O[tv,] = O[,tv] = 0
    }else{
      print("error")
      stop()
    }
    O_update[Rem_Node, Rem_Node] = O[Rem_Node, Rem_Node]
    O_update[tv,tv] = Inf
    A[tv,Rem_Node[ abs(v) * vars[tv] > hard_thresholding ] ] = 1
    count = count - 1
  }
  z[count] = Rem_Node
  return(list(Ordering = z, DAG = A))
}

### Inverse Covariance Matrix Estimation :: CLIME that requires threshold ###
inv_cov_estimate_clime = function(X, threshold = 0.1){
  
  p = ncol(X)
  n = nrow(X)
  
  #if( n <= p ){
  if( n/1000 <= p ){
    #install.packages("clime")
     library(clime)
    # re.clime <- clime( X, lambda.min=0.1^10, standardize=F,linsolver="simplex", pdmaxiter=100)
    # re.cv <- cv.clime(re.clime, fold = 10)
    # re.cv$lambdaopt
    # re.clime <- clime(X, lambda = re.cv$lambdaopt, standardize=F)
    # Omega_hat = re.clime$Omegalist[[1]]
    #re.clime <- clime(X, standardize=T, 0.001* sqrt( log(ncol(X))/nrow(X) ) )
    # re.clime <- clime(X, standardize=FALSE,  perturb=F, linsolver="simplex", 10)
     re.clime <- clime(X, threshold, standardize=FALSE)
     Omega_hat = re.clime$Omegalist[[1]]
    #install.packages("fastclime")
    #library(fastclime)
    #out1 = quiet(fastclime(cov(X), lambda.min = threshold , nlambda = 1))
    #out2 = fastclime.selector(out1$lambdamtx, out1$icovlist, 3* sqrt( log( p )/n ) )
    #Omega_hat = out2[[1]]
    #Omega_hat = out1$icovlist[[1]]
  }else{
    Omega_hat = solve(cov(X))
    #sum( (Omega_hat - solve(cov(X)) )^2 )
  }
  #Omega_hat[abs(Omega_hat) < threshold ] = 0
  return(Omega_hat)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 



