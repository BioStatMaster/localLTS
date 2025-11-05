#source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/inv_cov_estimate.R"))

# If the absolute value of the precision matrix
# is below this value, then assume that the entry is 0.
# Used for calculating support.
#INV_COVARIANCE_THRESHOLD = 0.0001

#dirname(rstudioapi::getActiveDocumentContext()$path)

ols_estimate = function(Sigma, S, i){
  #   """
  # :param Sigma: Empirical covariance matrix.
  # :param S: Support set.
  # :param i: variable for which the ols estimate has to be computed.
  # :return: ols estimate.
  # """
  theta = solve(Sigma[S,S])%*%(Sigma[S,i])
  theta = as.vector(t(theta))
  return(theta)
}

get_markov_blanket = function(O, i, 
                              inv_cov_threshold=INV_COVARIANCE_THRESHOLD){
  x = abs(O[i,]) > inv_cov_threshold
  x[i] = F
  mb = which(x) 
  return(mb)
}

get_dag_order = function(Sigma, O_hat,
                  inv_cov_threshold=INV_COVARIANCE_THRESHOLD){
  #   """
  # Learn the causal ordering, given the empirical
  # covariance marix and the estimated precision matrix.
  # :param Sigma: Empirical covariance matrix.
  # :param O_hat: Estimated precision matrix.
  # :param inv_cov_threshold: value below which entries in
  # the inverse covariance matrix are assumed to be 0.
  # :return:
  
  n = nrow(Sigma)
  z = rep(0,n)       # DAG order.
  vars = rep(0,n)    # variances.
  r = rep(0,n)                   # ratios.
  Rem_Node = 1:n
  
  count = n
  O = O_hat # light copy
  
  # Compute markov blanket for all nodes.
  # This is in terms of node labels.
  mbs = list()
  
  for(i in 1:n){
    mbs[[i]] = get_markov_blanket(O, i, inv_cov_threshold=inv_cov_threshold)
  } # mbd python list -> R list 로 변경
  
  for(i in 1:n){
    S = mbs[[i]]
    if( length(S) > 0){
      theta_i = ols_estimate(Sigma, S, i)
      Sp = which(abs(theta_i) > 0.001) #?
      if(prod(dim(Sp)) > 0){
        o = as.vector(O_hat[i,S])
        o = o[Sp]
        r[i] = max(abs(o/theta_i[Sp]))
        }
      }
  }
  
  # print r
  
  while( count >= 2){ #count check.
      # print "============================="
      # print r
      # print "============================="
    tv = which.min(r)  # This is the node label.
    z[count] = tv #
    r[tv] = Inf
    Rem_Node = setdiff(Rem_Node, tv)
    vars[tv] = Sigma[tv, tv] 
    if(O[tv, tv] > 0){ 
      vars[tv] = 1/O[tv, tv] ### dimension 안맞음.
    }
    
    # Delete terminal vertex from precision matrix.
    v = O[tv,Rem_Node]

    if( O[tv, tv] > 0){
      O[Rem_Node,Rem_Node] = O[Rem_Node,Rem_Node] - outer(v,v)/O[tv, tv]
      O[tv,] = O[,tv] = 0
    }else{
      print("error")
    }
    count = count - 1
    
    # Recompute ratios for variables in markov blanket of tv.
    for( j in mbs[[tv]]){
      # print "Recomputing ratio for:", j,
      mbs[[j]] = get_markov_blanket(O, j, inv_cov_threshold=inv_cov_threshold)
      # print "Updated markov blanket of ", j, "is:", mbs[j]
      
      r[j] = 0
      S = mbs[[j]]
      if(length(S) > 0){
        theta_j = ols_estimate(Sigma, S, j)
        Sp = abs(theta_j) > 0
        if( sum(Sp != 0) == 0 )next
        
        # print "theta = ", theta_j.round(5)
        O_j = O[j,][S]
        
        # print "Omega_j = ", O_j
        r[j] = max(abs(O_j[Sp]/theta_j[Sp]))
        # print " = ", r[j]
      }
    }
  }
  z[1] = which.min(r)
  vars[z[1]] = Sigma[z[1], z[1]]
  if( O[z[1], z[1]] > 1){
    vars[z[1]] =  1/O[z[1], z[1]]
  }
  
  return(list(z = z, vars = vars))
}


get_bn_from_order = function(Sigma, O_hat, dag_order, vars,
                      inv_cov_threshold=INV_COVARIANCE_THRESHOLD){
  # """
  # Learn from GBN given DAG order.
  # :param Sigma: n x n covariance matrix.
  # :param dag_order:
  # :param inv_cov_threshold: value below which the entries
  # in the inverse covariance matrix is assumed to be 0.
  # :return: (adjacency matrix, weights, noise variances)
  # """
  n = nrow(Sigma)
  A = matrix(0, ncol = n, nrow = n)
  W = matrix(0, ncol = n, nrow = n)
  mbs = list()
  
  for( i in  1:n) mbs[[i]] = get_markov_blanket(O_hat, i, inv_cov_threshold=inv_cov_threshold)

  # position of a node in the dag order.
  pos_in_dag_order = sort(dag_order, index.return=TRUE)$ix #
  for( i in 2:n){
    v = dag_order[i] # Variable for which we
    # want to parents and parameters.
    # print "learning parameters for node:", v
    S = intersect( mbs[[v]], dag_order[1:i] )
    if(length(S) > 0 ){
      theta_v = ols_estimate(Sigma, S, v)
      for( j in 1:length(theta_v) ) {
        u = S[j]
        if( abs(theta_v[j]) > inv_cov_threshold){
          A[v, u] = 1
          W[v, u] = theta_v[j]
        }
      }
    }
  }
  return(list(A = A, W = W, vars = vars))
}

  

learn_gbn = function(X, Sigma, reg_param = 0.05,
              inv_cov_threshold,
              use_marginal_variance_for_dag_order=FALSE){
  #   """
  # Learn a GBN from data.
  # :param Sigma: Sample covariance matrix.
  # :param reg_param:
  # :param inv_cov_threshold:
  # :return: An @GaussianBayesianNetwork instance.
  # """
  # Sigma = quantize(Sigma, threshold=inv_cov_threshold)
  tryCatch(
    {
      O_hat = inv_cov_estimate_clime(X, threshold = reg_param)
    },
    
    error = function(e){
      # print( e.message)
      message( "Failed to learn inverse covariance matrix.")
      return( NULL)
    }
  )
  
  tryCatch(
    {
      if(!use_marginal_variance_for_dag_order){
    result= get_dag_order(Sigma, O_hat, inv_cov_threshold= inv_cov_threshold)
    dag_order = result[[1]]
    v = result[[2]]
    # ordering estimation in population.
    # result = get_dag_order(Sigma = solve(diag(p)-B) %*% diag(noiseVar) %*% solve(diag(p) -t(B)), O_hat = (diag(p) -t(B)) %*% diag(1/noiseVar) %*% (diag(p)-B), inv_cov_threshold= INV_COVARIANCE_THRESHOLD)

    # 출력형식 확인.
    }else{
    # Use marginal variance to figure out DAG order.
    dag_order = sort(diag(Sigma), index.return = T)$ix
    v = rep(1,nrow(Sigma))
    }
      },
    error = function(e){
      message("Failed to learn DAG order.")
      # print(e.message)
      # tb.print_exc()
      return( NULL)
    }
  )
  
  # print "DAG order:", dag_order
  bn = get_bn_from_order(Sigma, O_hat, dag_order, v,
                         inv_cov_threshold=inv_cov_threshold)
  return(
    list( DAG = bn$'A', 
          Ordering = dag_order
        )
  )
}

GBN_Algorithm = function(X, Sigma, reg_param= 0.05,
                         inv_cov_threshold,
                         use_marginal_variance_for_dag_order=FALSE, graph = NULL){
  #### Start ###
  Runtime = proc.time()[3]
  result = learn_gbn(X, Sigma, reg_param,
                     inv_cov_threshold=inv_cov_threshold,
                     use_marginal_variance_for_dag_order=FALSE)
  Runtime = proc.time()[3] - Runtime
  print(paste("It takes: ", Runtime))
  
  Estimated_G = result$DAG
  est_MEC = dag2cpdagAdj(Estimated_G)
  Ordering = result$Ordering
  if( !is.null(graph) ){
    B = graph
    B[B!=0] =1
    MEC = B + t(B)
    evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
    evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)
  }
  
  
  return(
    list( DAG_Evaluation = evaluation_result_GSEM, 
          MEC_Evaluation = evaluation_result_GSEM_MEC,  
          Oracle_Evaluation = NULL, 
          DAG = Estimated_G, 
          Ordering = Ordering, 
          Time = Runtime)
  )
}


standardize = function(X){
  # standardize data.
  n = nrow(X)
  m = apply(X ,2, mean)
  s = apply(X , 2, function(X){sqrt(var(X)*(n-1)/n)})
  #python 은 분모를 n으로 사용.
  Z = (X - m)/s
  return( Z)
}

select_model = function(X, reg_params, thresholds){
  # """
  # Split data into train and test. Learn multiple GBNs
  # for all combinations regularization parameter and thresholds.
  # Return the model that has the highest test loglikelihood.
  # 
  # :param X: data (m x n), where m is #samples, n is #nodes.
  # :param reg_params: 
  # :param thresholds:
  # :return: The best model.
  # """
  X = as.matrix(X)
  Z = standardize(X)

  # split data into train and test.
  train_idx = sample( 1:nrow(Z ), round(0.5*nrow(Z)) )
  Z_train = Z[train_idx,] ; Z_test = Z[-train_idx,]

  print(paste("training data size:", paste(dim(Z_train),collapse =  'x')))
  print(paste("test data size:", paste(dim(Z_test) , collapse  = 'x')))
  
  # Compute sample covariance matrix.
  C_train = (1/nrow(Z_train)) * t(Z_train)%*%(Z_train)
  
  best_r = 0; best_t = 0 ;  min_neg_log_lik = Inf
  for( r in reg_params){
    for(t in thresholds){
      l = Inf
      bnh = learn_gbn(Z_train, C_train, r, inv_cov_threshold=t)
      if(!is.null(bnh))l = -bnh.log_likelihood(Z_test) ########## 이부분만 수정
      print( paste0("reg_param:",r,", threshold: ",t,", test negative log likelihood: ",round(l,3)))
      if( l < min_neg_log_lik){
        best_r = r ; best_t = t }
    }
  }
  C = (1./nrow(Z)) * t(Z)%*%Z
  bn = learn_gbn(Z, C, best_r * sqrt(1./2), t)
  
  return(list(bn = bn, best_r = best_r, best_t = best_t)) #  형식고려 
}

### Inverse Covariance Matrix Estimation :: CLIME that requires threshold ###
inv_cov_estimate_clime = function(X, threshold = 0.0001){
  if(1==2){
    #install.packages("fastclime")
    library(fastclime)
    omega1 = fastclime(as.matrix(X), lambda.min = 0.01, nlambda = 50)
    omega2 = fastclime.selector(out1$lambdamtx, out1$icovlist, 0.01)
    Omega_hat = omega2$icov
    ### (diag(p) -t(B)) %*% diag(1/noiseVar) %*% (diag(p)-B)
  }else{
    Omega_hat = solve(cov(X))
  }
  Omega_hat[abs(Omega_hat) < threshold ] = 0
  return(Omega_hat)
}

