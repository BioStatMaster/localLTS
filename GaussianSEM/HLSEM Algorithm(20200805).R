

###################
HLSEM_simulation_fun = function(seed, n_real, p_real, d = 1, path = F ){
  #This algorithm is for simulations using Gaussian SEM Learning_algorithm. 
  
  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################
  
  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$x[1:n_real,1:p_real]
  result_GSEM = HLSEM_Algorithm(data, sparsity_level = 3, CV = F, graph = graph)
  return(result_GSEM)
}

#### High Dimensional Linaer SEM Learning using Lasso #################################
HLSEM_Algorithm = function(data, sparsity_level = 1, CV = F, graph = NULL){
  #install.packages("glasso")
  library(glasso)
  library(bnlearn)
  library(glmnet)
  
  ###################
  X = as.matrix(data)
  X = apply( X, 2, function(x) x- mean(x) )
  Sigma = cov(X)
  p = ncol(X)
  n = nrow(X)
  RemNodes = 1:p 
  pi_GSEM = NULL
  Ordering = score = rep(0, p)
  Estimated_G = matrix(0, p ,p)
  
  ####
  Runtime = proc.time()[3]
  
  ### Step (1): Inverse Covariance Estimation ###
  while( length(RemNodes) > 1 ){
    i = length(RemNodes)
    library(huge)
    Y = ( X[,RemNodes] )
    if(CV){
      Y = scale(Y)
      MB = huge(Y, method = "mb")
      MB = huge.select(MB)
      MB = huge(Y, method = "mb", lambda = 1.0*MB$opt.lambda) # half for small.   
    }else{
      # sparsity_level
      MB = huge(Y, method = "mb", lambda = 3*sqrt( log(p)/n ) ) # half for small.
    }
    Omega = as.matrix( MB$path[[1]] )
    diag(Omega) <- 0
    
    for(j in RemNodes){
      S = RemNodes[ which( Omega[,which(j == RemNodes)] !=0 ) ]
      if( length(S) > 0){
      score[j] = Sigma[j,j] - Sigma[j,S] %*% solve(Sigma[S,S]) %*% Sigma[S,j]
      }else{
      score[j] = Sigma[j,j]
      }
    }
    ### Step (2): Ordering Estimation ###
    Ordering_id = which( score == max(score) )
    if( length(Ordering_id) > 1 ){
      Ordering[i] = sample(Ordering_id, 1)
    }else{
      Ordering[i] = Ordering_id
    }
    score[Ordering[i]] = -score[Ordering[i]]
    #### Step 3): Finding the Parents ####
    Cand_Pa = setdiff(RemNodes,Ordering[i])
    if( length( Cand_Pa ) > 1 ){
      if( CV ){
      cvglmfit = cv.glmnet( X[ ,Cand_Pa], X[,Ordering[i]], family = "gaussian", alpha = 1, nfolds = 10)
      min.lambda = cvglmfit$lambda.min;  se.lambda  = cvglmfit$lambda.1se
      glm_beta = glmnet(X[ ,Cand_Pa], X[,Ordering[i]], family = "gaussian",alpha = 1, 
                        lambda = sparsity_level * (se.lambda - min.lambda) + min.lambda )$beta[,1]
      }else{
        glm_beta = glmnet(X[ ,Cand_Pa], X[,Ordering[i]], family = "gaussian",alpha = 1, 
                        lambda = sparsity_level * sqrt( log(p)/n ) )$beta[,1]
      }
      Estimated_G[Ordering[i], Cand_Pa[ which(glm_beta !=0 ) ] ] = 1
    }else{
      
      if( CV ){
        cvglmfit = cv.glmnet( X[ ,rep(Cand_Pa,2)], X[,Ordering[i]], family = "gaussian", alpha = 1, nfolds = 5)
        min.lambda = cvglmfit$lambda.min;  se.lambda  = cvglmfit$lambda.1se
        glm_beta = glmnet(X[ ,rep(Cand_Pa,2)], X[,Ordering[i]], family = "gaussian",alpha = 1,
                          lambda = sparsity_level * (se.lambda - min.lambda) + min.lambda )$beta
      }else{
        glm_beta = glmnet(X[ ,rep(Cand_Pa,2)], X[,Ordering[i]], family = "gaussian",alpha = 1,
                          lambda = sparsity_level * sqrt( log(p)/n ) )$beta
      }
      if( any( glm_beta !=0 ) ){ 
        Estimated_G[Ordering[i], Cand_Pa ] = 1
      }
    }
    RemNodes = setdiff(RemNodes, Ordering[i])
  }
  Ordering[1] = RemNodes
  
  Runtime = proc.time()[3]- Runtime
  ####
  print(paste("HLSEM takes: ", Runtime))
  ####

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

# #### HGSEM Learning Algorithm #########
# HGSEM_Algorithm = function(data,  pc_alpha =0.01, 
#                            alpha = 0.001, sparsity_level = 0.85, graph = NULL){
#   #install.packages("glasso")
#   library(glasso)
#   library(glmnet)
#   library(bnlearn)
#   library(corrplot)
#   ###################
#   X = as.matrix( data )
#   X = apply(X, 2, function(x) x - mean(x) )
#   p = ncol(X)
#   n = nrow(X)
#   RemNode = 1:p 
#   pi_GSEM = NULL
#   
#   ####
#   Runtime = proc.time()[3]
#   
#   ### Step 1): Inverse Covariance Estimation ###
#   Omega = matrix(0,p,p)
#   ### graphical lasso based
#   if(1 == 1){
#     X = scale(X)
#     Omega = glasso( cov(X), rho = sparsity_level )$wi
#     #Omega = glasso( cov(X), rho = 0.50 )
#     #Omega = glasso( cov(X), rho = 0.10, w.init = Omega$w, wi.init = Omega$wi )$wi
#     
#     diag(Omega) <- 0 
#     Omega[Omega!=0] = 1
#     #corrplot(Omega)
#   }  
#   
#   ### graphical lasso based using new pacakge
#   if(1 == 2){
#     #install.packages("CovTools")
#     library(CovTools)
#     Omega = PreEst.glasso(X, method = list(type = "fixed", param = sparsity_level),
#                   parallel = FALSE)$C
#     diag(Omega) <- 0 
#     Omega[Omega!=0] = 1
#     corrplot(Omega)
#   }
#   
#   ### clime based using new pacakge
#   if(1 == 2){
#     #install.packages("flare.tiger")
#     library(flare.tiger)
#     Omega = flare.tiger(X) 
#     
#     diag(Omega) <- 0 
#     Omega[Omega!=0] = 1
#     corrplot(Omega)
#   }
#  
#   ### PC based
#   if(1 == 2){
#     MEC = pc.stable(data, alpha = pc_alpha)
#     Omega = t(amat(cpdag(MEC)))
#     Omega = t(Omega) + Omega 
#     Omega[Omega!=0] = 1
#   }
#   
#   ### lasso based
#   if(1 == 2){
#     X = scale(X)
#     for( i in 1:p){
#       #lambda_lasso = cv.glmnet( X[,setdiff(1:p,i)], X[,i] , family = "gaussian", alpha=1, nfold = 10, thresh = 1e-06)
#       #lambda_optimal = max( min(lambda_lasso$glmnet.fit$lambda) ,
#       #                      lambda_lasso$lambda.min + sparsity_level *( lambda_lasso$lambda.1se - lambda_lasso$lambda.min )
#       #  )
#       #lambda_optimal = lambda_lasso$lambda.min + sparsity_level *( lambda_lasso$lambda.1se - lambda_lasso$lambda.min )
#       #glmnet0 = lambda_lasso$glmnet.fit$beta[,min(which( lambda_lasso$glmnet.fit$lambda <= lambda_optimal) )]
#       
#       lambda_optimal = max(0, sparsity_level * sqrt( log(p)/n  ) )
#       glmnet0 = glmnet( X[,setdiff(1:p,i)], X[,i] , family = "gaussian", alpha=1, lambda = lambda_optimal, thresh = 1e-06)$beta
#       
#       #nlambda = length(lambda_lasso$glmnet.fit$lambda)
#       #glmnet0 = lambda_lasso$glmnet.fit$beta[, round(nlambda*sparsity_level,0) ]
#       
#       #glmnet0 = glmnet( X[,setdiff(1:p,i)], X[,i], family = "gaussian", alpha=1, lambda = sparsity_level*sqrt( log(p)/n )/1000 )$beta
#       pos001 = which( glmnet0 != 0 )
#       if( length(pos001) > 0 ){
#         Omega[i, setdiff(1:p,i)[pos001] ] = Omega[setdiff(1:p,i)[pos001], i] = 1 ;
#       }
#     }
#   } 
#   
#   print(paste("Step (1) takes: ", proc.time()[3]- Runtime))
#   
#   #### Step 2): Finding the Ordering ####
#   X = as.matrix( data )
#   X = apply(X, 2, function(x) x - mean(x) )
#   
#   Ordering = rep(0, p)
#   Estimated_G = matrix(0, p ,p)
#   RemNodes = 1:p
#   Ordering =rep(0, p)
#   ConditionSet = list(NULL)
#   for(i in 1:p){ ConditionSet[[i]] = 0 }
#   i = 1
#   score = apply(X, 2, function(x) var(x))
#   score_id = which(score == min(score))
#   if( length(score_id) == 1){
#     Ordering[i] = score_id
#   }else{
#     Ordering[i] = sample(score_id, 1)
#   }
# 
#   score[Ordering] = Inf
#   RemNodes = setdiff(RemNodes, Ordering)
#   while( length(RemNodes)>1 ){
#     
#     i = i + 1
#     # score update #
#     for(j in RemNodes){
#       k = intersect(Ordering[1:(i-1)], which( Omega[,j] !=0 ) )
#       
#       if( !identical( sort(ConditionSet[[j]]) , sort(k) ) && length(k) > 0 ){
#         ConditionSet[[j]] = k
#         ScoreData = data.frame( Y = X[,j], Z = X[,ConditionSet[[j]]] )
#         resid = lm(Y~ -1 + . , ScoreData)$residual
#         score[j] = sum(resid^2)/n
#       }
#     }
#     #Odering Selection
#     Ordering_id = which( score == min(score) )
#     if( length(Ordering_id) > 1 ){
#       Ordering[i] = sample(Ordering_id, 1)
#     }else{
#       Ordering[i] = Ordering_id
#     }
#     score[Ordering] = Inf
#     RemNodes = setdiff(RemNodes, Ordering[i])
#   }
# 
#   Ordering[p] = RemNodes  
#   k = intersect(Ordering[1:(p-1)], which( Omega[,p] !=0 ) )
#   if( !identical( sort(ConditionSet[[p]]) , sort(k) ) && length(k) > 0 ){
#     ConditionSet[[p]] = k
#   }
#   
#   print(paste("Step (2) takes: ", proc.time()[3]- Runtime))
#   
#   #### Step 3): Finding the Parents ####
#   
#   
#   
#   #### Algorithm Starts for finding Edges ####
#   Estimated_G = matrix(0, p, p)
# #  lambda.lasso2 = log(p)/40
# 
#   if(1 == 2){
#     for(m in 2:p){
#       j = Ordering[m]
#       if(length(ConditionSet[[j]]) >0 ){
#         ScoreData = data.frame( Y = X[,j], Z = X[,ConditionSet[[j]]] )
#         beta_j = lm(Y~ -1 + . , ScoreData)$coef
#         pa_j = ConditionSet[[j]][which( abs(beta_j) > 0.35 )]
#         if( length(pa_j) >0){
#           Estimated_G[j, pa_j] = 1
#         }
#       }
#     }
#   }
#   #### Finding Edges using CI ####
#   
#   if( 1 == 1){
#   #### Step 2): Finding the Parents ####
#   for(m in 2:p){
#     j = Ordering[m]
#     if( !all( ConditionSet[[j]] ==0) ){
#       for(k in ConditionSet[[j]]){
#         S = setdiff( ConditionSet[[j]], k ) 
#         if( length(S) > 0 ){
#           parent_pvalue = ci.test(X[,j], X[,k], X[,S], test = "cor")$p.value
#         }else{
#           parent_pvalue = ci.test(X[,j], X[,k], test = "cor")$p.value
#         }
#         if(parent_pvalue < alpha){
#           Estimated_G[j, k] = 1
#         }
#       }
#     }
#   }
#   }
#   
#   #### Finding Edges using Lasso ####
#   if( 1 == 2){
#     #sparsity_level = 0
#     lambda_lasso = cv.glmnet( X[,Ordering[rep(1,2)]], X[,Ordering[2]] , family = "gaussian", alpha=1, nfold = 10, thresh = 1e-06)
#     lambda_optimal = lambda_lasso$lambda.min + sparsity_level *( lambda_lasso$lambda.1se - lambda_lasso$lambda.min )
#     glmnet0 = lambda_lasso$glmnet.fit$beta[,which.min( lambda_lasso$glmnet.fit$lambda > lambda_optimal )]
#     pos001 = which( glmnet0 != 0 )
#     if( !is.null(pos001) ){
#       Estimated_G[Ordering[2],Ordering[1]] = 1 ;
#     }
#     for( i in 3:p){
#       lambda_lasso = cv.glmnet( X[,Ordering[1:(i-1)]], X[,Ordering[i]] , family = "gaussian", alpha=1, nfold = 5, thresh = 1e-06)
#       lambda_optimal = lambda_lasso$lambda.min + sparsity_level *( lambda_lasso$lambda.1se - lambda_lasso$lambda.min )
#       glmnet0 = lambda_lasso$glmnet.fit$beta[,which.min( lambda_lasso$glmnet.fit$lambda > lambda_optimal )]
#       #glmnet0 = glmnet( as.matrix(X[,Ordering[1:(i-1)]]), X[,Ordering[i]] , family = "gaussian", alpha=1, lambda = sparsity_level*sqrt( log(p)/n ) )$beta
#       pos001 = which( glmnet0 != 0 )
#       if( length(pos001) > 0 ){
#         Estimated_G[Ordering[i], Ordering[ pos001 ]] =1 ;
#       }
#     }
#   }  
#   
#   print(paste("Step (3) takes: ", proc.time()[3]- Runtime))
#   
#   print( sapply(1:p, function(j) length(ConditionSet[[j]]  ) ) )
#   
#   
#   ####
#   Runtime = proc.time()[3] - Runtime
#   
#   ####
#   est_MEC = dag2cpdagAdj(Estimated_G)
#   
#   ####
#   
#   if( !is.null(graph) ){
#     B = graph
#     B[B!=0] =1
#     MEC = B + t(B)
#     evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
#     evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)
#   }
#   
#   return(
#     list( DAG_Evaluation = evaluation_result_GSEM, 
#           MEC_Evaluation = evaluation_result_GSEM_MEC,  
#           Oracle_Evaluation = NULL, 
#           DAG = Estimated_G, 
#           Ordering = Ordering, 
#           Time = Runtime)
#   )
# }
# 
# ############## new HGSEM Algorithm using shur complement: real fast #######################################################
# 
# HGSEM_Algorithm2 = function(X, alpha, graph = NULL){
#   
#   Runtime = proc.time()[3]  
#   
#   sig0 = cov(X)
#   
#   pi= rep(0, p)
#   RemNode = 1:p
#   U_score = apply(X,2, var )
#   pi[1] = which.min(U_score)
#   
#   
#   RemNode = setdiff(RemNode, pi)
#   A = solve(sig0[pi[1],pi[1]]) 
#   Estimated_G = matrix(0, p, p)
#   
#   for(j in 2:p ){
#   #  j = j+1
#     S = pi[1:(j-1)]
#     U_score = sapply(RemNode, function(k) sig0[k,k] - sig0[k,S] %*% A %*% sig0[S,k] )
#     pi[j] = RemNode[which.min(U_score)]
#     RemNode = setdiff(RemNode, pi)
#     
#     # update inverse covariance matrix
#     K = 1/min(U_score)
#     #K = 1/min(U_score[1])
#     sig12 = -K %*% sig0[pi[j] ,S] %*% A
#     sig11 = A + A %*%  sig0[S,pi[j]] %*% K %*% sig0[pi[j],S] %*% A
#     updated_A = rbind( cbind(sig11, t(sig12)), cbind(sig12, K) ) 
#     
#     ###
#     r = -( sig12 ) / sqrt(K * diag(sig11) ) 
#     parents = S[-abs(log( (1+r)/(1-r) )/2*sqrt(n-j-3) ) < qnorm(alpha)]
#     Estimated_G[pi[j], parents] = 1
#     
#     ##
#     A = updated_A
# 
#   }
#   
#   
#   Ordering = pi
#   
#   print(paste("Step (1) takes: ", proc.time()[3]- Runtime))
#   
#   #### Algorithm Starts for finding Edges ####
#   #  lambda.lasso2 = log(p)/40
#   
#   Runtime = proc.time()[3]- Runtime
#   
#   ####
#   #est_MEC = dag2cpdagAdj(Estimated_G)
#   
#   ####
#   
#   
#   if( !is.null(graph) ){
#     B = graph
#     B[B!=0] =1
#     MEC = B + t(B)
#     evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
#     #evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)
#   }
#   
#   result = list( DAG_Evaluation = evaluation_result_GSEM, 
#                  MEC_Evaluation = NULL,  
#                  Oracle_Evaluation = NULL, 
#                  DAG = Estimated_G, 
#                  Ordering = Ordering, 
#                  Time = Runtime)
#   
#   #################
# }
# 
# 
# ######################
# 
# 
# ###############
# 
# 
# 
# ####HGSEM Learning Algorithm 3  #########
# HGSEM_Algorithm3 = function(data,  pc_alpha =0.01, 
#                            alpha = 0.001, sparsity_level = 0.85, graph = NULL){
#   #install.packages("glasso")
#   library(glasso)
#   library(glmnet)
#   library(bnlearn)
#   library(corrplot)
#   ###################
#   X = as.matrix( data )
#   X = apply(X, 2, function(x) x - mean(x) )
#   p = ncol(X)
#   n = nrow(X)
#   RemNode = 1:p 
#   pi_GSEM = NULL
#   
#   ####
#   Runtime = proc.time()[3]
#   
#   ### Step 1): Inverse Covariance Estimation ###
#   Omega = matrix(0,p,p)
#   ### graphical lasso based
#   if(1 == 1){
#     X = scale(X)
#     Omega = glasso( cov(X), rho = sparsity_level )$wi
#     #Omega = glasso( cov(X), rho = 0.50 )
#     #Omega = glasso( cov(X), rho = 0.10, w.init = Omega$w, wi.init = Omega$wi )$wi
#     
#     diag(Omega) <- 0 
#     Omega[Omega!=0] = 1
#     #corrplot(Omega)
#   }  
#   
#   
#   print(paste("Step (1) takes: ", proc.time()[3]- Runtime))
#   
#   #### Step 2): Finding the Ordering ####
#   X = as.matrix( data )
#   X = apply(X, 2, function(x) x - mean(x) )
#   
#   Ordering = rep(0, p)
#   Estimated_G = matrix(0, p ,p)
#   RemNodes = 1:p
#   Ordering =rep(0, p)
#   ConditionSet = list(NULL)
#   for(i in 1:p){ ConditionSet[[i]] = 0 }
#   i = 1
#   score = apply(X, 2, function(x) var(x))
#   score_id = which(score == min(score))
#   if( length(score_id) == 1){
#     Ordering[i] = score_id
#   }else{
#     Ordering[i] = sample(score_id, 1)
#   }
#   
#   score[Ordering] = Inf
#   RemNodes = setdiff(RemNodes, Ordering)
#   while( length(RemNodes)>1 ){
#     
#     i = i + 1
#     # score update #
#     for(j in RemNodes){
#       k = intersect(Ordering[1:(i-1)], which( Omega[,j] !=0 ) )
#       
#       if( !identical( sort(ConditionSet[[j]]) , sort(k) ) && length(k) > 0 ){
#         ConditionSet[[j]] = k
#         ScoreData = data.frame( Y = X[,j], Z = X[,ConditionSet[[j]]] )
#         
#         score[j] = var(X[,j]) - as.matrix( cov(X[,j], X[,k]) ) %*% as.matrix( solve( var(X[,k]) ) ) %*% t(as.matrix(cov(X[,j], X[,k])) )
#         #resid = lm(Y~ -1 + . , ScoreData)$residual
#         #score[j] = sum(resid^2)/n
#       }
#     }
#     #Odering Selection
#     Ordering_id = which( score == min(score) )
#     if( length(Ordering_id) > 1 ){
#       Ordering[i] = sample(Ordering_id, 1)
#     }else{
#       Ordering[i] = Ordering_id
#     }
#     score[Ordering] = Inf
#     RemNodes = setdiff(RemNodes, Ordering[i])
#   }
#   
#   Ordering[p] = RemNodes  
#   k = intersect(Ordering[1:(p-1)], which( Omega[,p] !=0 ) )
#   if( !identical( sort(ConditionSet[[p]]) , sort(k) ) && length(k) > 0 ){
#     ConditionSet[[p]] = k
#   }
#   
#   print(paste("Step (2) takes: ", proc.time()[3]- Runtime))
#   
#   #### Step 3): Finding the Parents ####
#   
#   
#   #### Finding Edges using CI ####
#   
#   if( 1 == 1){
#     #### Step 2): Finding the Parents ####
#     for(m in 2:p){
#       j = Ordering[m]
#       if( !all( ConditionSet[[j]] ==0) ){
#         for(k in ConditionSet[[j]]){
#           S = setdiff( ConditionSet[[j]], k ) 
#           if( length(S) > 0 ){
#             parent_pvalue = ci.test(X[,j], X[,k], X[,S], test = "cor")$p.value
#           }else{
#             parent_pvalue = ci.test(X[,j], X[,k], test = "cor")$p.value
#           }
#           if(parent_pvalue < alpha){
#             Estimated_G[j, k] = 1
#           }
#         }
#       }
#     }
#   }
#   
#   
#   print(paste("Step (3) takes: ", proc.time()[3]- Runtime))
#   
#   print( sapply(1:p, function(j) length(ConditionSet[[j]]  ) ) )
#   
#   
#   ####
#   Runtime = proc.time()[3] - Runtime
#   
#   ####
#   est_MEC = dag2cpdagAdj(Estimated_G)
#   
#   ####
#   
#   if( !is.null(graph) ){
#     B = graph
#     B[B!=0] =1
#     MEC = B + t(B)
#     evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
#     evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)
#   }
#   
#   return(
#     list( DAG_Evaluation = evaluation_result_GSEM, 
#           MEC_Evaluation = evaluation_result_GSEM_MEC,  
#           Oracle_Evaluation = NULL, 
#           DAG = Estimated_G, 
#           Ordering = Ordering, 
#           Time = Runtime)
#   )
# }
# 
# 
# ##############

