

##### PGSEM_simulation_fun: Simulations for recovering PGSEM ####
PGSEM_simulation_fun = function(seed, n_real, p_real, d = 1, direction = "forward", sparsity_level = 5, path = F ){
   #This algorithm is for simulations using Gaussian SEM Learning_algorithm. 
   
   ### Gaussian SEM Load #############
   setwd(path)
   GaussianDAG_filename = paste0("PGSEM_d", d, "_seed", seed, ".Rdata")
   load(GaussianDAG_filename)
   ###################################
   
   synthetic.graph =DAGsamples
   graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
   data = synthetic.graph$x[1:n_real,1:p_real]
   
   result_GSEM = PGSEM_Algorithm(data, direction ="forward",  sparsity_level = sparsity_level, CV = F, graph = graph)
   
   return(result_GSEM)
}


########## Ordering Estimation Function (Ordering Only) #########       
PGSEM_Forward_Learning_fun = function(data){
  
   X = as.matrix(data)
  
   #### Step 1): Finding the Ordering ####
   p = ncol(X)
   n = nrow(X)
   
   Ordering = rep(0, p)
   RemNodes = 1:p
   k = 1
   
   ## first element of the ordering ##
   scores = apply(X, 2, function(x) var(x)/abs( mean(x) ) ) 
   Ordering[k] = which.min( scores )
   
   ## m-th element of the ordering ##
   RemNodes = setdiff(RemNodes, Ordering)
   while( length(RemNodes) > 1 ){
      ###
      k = k + 1
      scores = NULL
      for(j in RemNodes){
         lm_result_j = lm(X[,j] ~ X[,Ordering[1:(k-1)]] )
         score_j = var( (lm_result_j$residuals/sqrt(abs(lm_result_j$fitted.values)) )[lm_result_j$fitted.values!=0])
         ## I should check, in the paper, how we address this issue. When fitted.values are 0.
         ## I think this is a main weakness of our paper.
         scores <- c(scores, score_j)
         print(scores)
      }
      Ordering_j = which( scores == min(scores) )
      if( length(Ordering_j) > 1){
         Ordering[k] = RemNodes[sample(Ordering_j, 1)]
         ## Is sampling the best solution? 
      }else{
         Ordering[k] = RemNodes[Ordering_j]
      }
      RemNodes = setdiff(RemNodes, Ordering)
   }
   Ordering[p] = RemNodes
   return( Ordering )
}




#### High Dimensional Linear SEM Learning using Lasso #################################
PGSEM_Algorithm = function(data, direction = "forward", sparsity_level = 1, CV = F, graph = NULL){
   #install.packages("glasso")
   library(glasso)
   library(bnlearn)
   library(glmnet)
   
   ###################
  
   X = as.matrix(data)
   p = ncol(X)
   n = nrow(X)
   RemNodes = 1:p 
   pi_GSEM = NULL
   Ordering = score = rep(0, p)
   Estimated_G = matrix(0, p ,p)
   
   #### Setting initial time of operating the code
   Runtime = proc.time()[3]
   
   #### Step 1): Finding the Ordering ####
   if(direction == "forward"){
      Ordering = PGSEM_Forward_Learning_fun(data = X)
   }else if(direction =="backward"){
      Ordering = PGSEM_Forward_Learning_fun(data = X)
   }
   print( proc.time()[3]  - Runtime )
   
   #### Step 2): Finding the Parents ####
   X = scale(X)
   ## scaling?
   ## Selecting a lambda should be adjusted after theoretical guaruntee.
   for(i in 2:p){
      Cand_Pa = Ordering[1:(i-1)]
      if( length( Cand_Pa ) > 1 ){
        ## Selecting lambda through CV
         if( CV ){
            cvglmfit = cv.glmnet( X[ ,Cand_Pa], X[,Ordering[i]], family = "gaussian", alpha = 1, nfolds = 10)
            min.lambda = cvglmfit$lambda.min;  se.lambda  = cvglmfit$lambda.1se
            glm_beta = glmnet(X[ ,Cand_Pa], X[,Ordering[i]], family = "gaussian", alpha = 1, 
                              lambda = sparsity_level * (se.lambda - min.lambda) + min.lambda )$beta[,1]
         }else{
            glm_beta = glmnet(X[ ,Cand_Pa], X[,Ordering[i]], family = "gaussian", alpha = 1, 
                              lambda = sparsity_level * sqrt( log(p)/n ) )$beta[,1]
         }
         Estimated_G[Ordering[i], Cand_Pa[ which(glm_beta !=0 ) ] ] = 1
      }else{
         if( CV ){
            cvglmfit = cv.glmnet( X[ ,rep(Cand_Pa,2)], 
                                  X[,Ordering[i]], 
                                  family = "gaussian", 
                                  alpha = 1, 
                                  nfolds = 5)
            min.lambda = cvglmfit$lambda.min;  se.lambda  = cvglmfit$lambda.1se
            glm_beta = glmnet(X[ ,rep(Cand_Pa,2)], 
                              X[,Ordering[i]], 
                              family = "gaussian",
                              alpha = 1,
                              lambda = sparsity_level * (se.lambda - min.lambda) + min.lambda )$beta
         }else{
            glm_beta = glmnet(X[ ,rep(Cand_Pa,2)], 
                              X[,Ordering[i]], 
                              family = "gaussian",
                              alpha = 1,
                              lambda = sparsity_level * sqrt( log(p)/n ) )$beta
         }
         if( any( glm_beta !=0 ) ){ 
            Estimated_G[Ordering[i], Cand_Pa ] = 1
         }
      }
   }
   Runtime = proc.time()[3]- Runtime
   ####
   print(paste("PGSEM takes: ", Runtime))
   ####
   
   
   
   ##### Evaulation #####
   if( !is.null(graph) ){
      est_MEC = dag2cpdagAdj(Estimated_G)
      B = graph
      B[B!=0] = 1
      MEC = B + t(B)
      evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
      evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC )
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
   
   return(
      list( DAG_Evaluation = evaluation_result_GSEM, 
            MEC_Evaluation = evaluation_result_GSEM_MEC,  
            Oracle_Evaluation = evaluation_result_GSEM2, 
            DAG = Estimated_G,
            Ordering = Ordering, 
            Time = Runtime)
   )
}






