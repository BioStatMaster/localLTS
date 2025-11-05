
### Graph Structure Generator ###
Graph_Generator = function( p, d, graph_type, beta_min,  beta_max, q= 0.01, h = 1, seed = 1, structure = NULL ){
  # p : number of nodes
  # d : maximum indegree
  # beta_min: minimum edge weight
  # beta_max: maximum edge weight
  # graph_type: 1: tree, 5: random, 6: quadratic random graph, 7 : dense graph
  # h : number of hubs
  set.seed(seed)
  B = B2 = matrix( 0, p ,p )
  
  # given structure
  if( graph_type == 0 ) { 
    B = structure
  }
  
  # tree
  if( graph_type == 1 ) { 
    for(i in 2:p){
      parents_i = (i-1)
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
    }
  }
  # hub graph 
  if( graph_type == 2 ) { 
    for(k in 2:3){
      for(i in k:p){
        parents_i = (i-k+1)
        B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
      }
    }
  }
  
  #Randomg graph with random number of parents
  if( graph_type == 3 ){
      new.seed = seed + 100
      set.seed(new.seed)
      B = B2 = matrix( 0, p ,p )
      B2[,1] = rep(1, p)
      iteration = 0
      while( max( rowSums( B2+t(B2)) ) >= d ) { 
        iteration = iteration + 1
        if( q > 0.0001 ){ q = q - 0.0001}
        B = B2 = B3 = matrix( 0, p ,p )
        for(i in 2:p){
          parents_i = sample( 1:(i-1), rbinom(1, (i-1), q ),replace = FALSE )
          B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
        }
        B2[B!=0] = 1
        B3 = B2
        #library(CombMSC)
        for(j in 3:p){
          S  = which( B2[j,] == 1 )
          if(length(S)>=3){
            SI = subsets(length(S), 2, v = S)
            for(k in 1:nrow(SI) ){
              B2[SI[k,1], SI[k,2]] = 1
              B2[SI[k,2], SI[k,1]] = 0
            } 
          }else if( length(S)==2){
            B2[S[1], S[2]] = 1
            B2[S[2], S[1]] = 0
          }
        }
      }
      #print( max( rowSums( B2+t(B2)) ) )
      #B2[B!=0]<-1 
      #print( max( rowSums( B3 ) ) )
  }
  #Randomg graph with random number of parents
  if( graph_type == 4 ) {
    for(i in 2:p){
      parents_i = sample( 1:(i-1), sample(0:min(i-1, d), 1),replace = FALSE )
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( x= c(-1,1), size = length(parents_i), replace = TRUE)
    }
  }
  #Randomg graph with fixed number of parents
  if( graph_type == 5 ) {
    for(i in 2:p){
      parents_i = sample( 1:(i-1), min(i-1, d),replace = FALSE )
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
    }
  }
  #Randomg quadratic graph with fixed number of parents
  if( graph_type == 6 ) {
    for(i in 2:p){
      parents_i = sample( 1:(i-1), min(i-1, d),replace = FALSE )
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
      if( sample( c(0), 1) == 0 ){
        B2[i,parents_i ] = runif( length(parents_i) , 0.20, 0.20 )* sample( c(-1,1), length(parents_i), replace = T) 
      }
      #B3[i,parents_i ] = runif( length(parents_i) , beta_min/3,  beta_max/3)* sample( c(-1,1), length(parents_i), replace = T)
    }
  }
  
  # Hub Graphs with the number of Hubs h #
  if( graph_type == 7 ) {
    Hubs = seq(from = 1, to = p, by = ceiling((p-1)/h))
    for(i in 1:(length(Hubs)-1)){
      B[(Hubs[i]+1):Hubs[i+1],Hubs[i]] = 
        runif(Hubs[i+1]-Hubs[i],beta_min,  beta_max)*
        sample(c(-1,1),Hubs[i+1]-Hubs[i], replace = T)
    }
    for(i in length(Hubs)){
      if(Hubs[i] < p ){
        B[(Hubs[i]+1):p, Hubs[i]] = 
          runif(p-Hubs[i],beta_min,  beta_max)*
          sample(c(-1,1),p-Hubs[i], replace = T)
      }
    }
  }
  
  # 5-topological-Layer Graph #
  if( graph_type == 8 ){
    q = 25
    A1 = 1:(q/5); 
    A2 = (q/5+1):(2*q/5); 
    A3 = (2*q/5+1):(3*q/5);
    A4 = (3*q/5+1):(4*q/5);
    A5 = (4*q/5+1):(p);
    
    B[A2,A1] <- runif( (q/5)^2 , beta_min,  beta_max)* sample( c(-1,1), length((q/5)^2), replace = T)
    B[A3,A2] <- runif( (q/5)^2 , beta_min,  beta_max)* sample( c(-1,1), length((q/5)^2), replace = T)
    B[A4,A3] <- runif( (q/5)^2 , beta_min,  beta_max)* sample( c(-1,1), length((q/5)^2), replace = T)
    B[A5,A4] <- runif( (q/5)*(p-4*q/5) , beta_min,  beta_max)* sample( c(-1,1), length((q/5)^2), replace = T)
  }
  
  # Star Graph #
  if( graph_type == 9 ){
    B[2:p,1] <- runif( p-1 , beta_min,  beta_max)* sample( c(-1,1), p-1, replace = T)
  }
  
  return( list(B = B, B2 = B2) )
}

### PGSEM generator function ###
LTS_GSEM_generator = function( n, p, d = 2, b = 0, outlier_node = NULL, var_Min = 1, var_Max = 1, dist = "Gaussian", beta_min, beta_max, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1, path = FALSE){
  
  ## b is # of outliers
  ## outlier_node: corrupted variables. 
  
  ### Generating data with p, n ###
  B = matrix( 0, p ,p )
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed = seed, structure)$B
  
  ### Sample with respect to the given distributions ###
  noiseVar = runif(p, min = var_Min, max = var_Max)
  x= matrix( 0, n, p)
  for( i in 1:p){
    non_zero_id = which( B[i,1:(i-1) ] != 0)
    theta_i = colSums( B[i,non_zero_id ]*t( x[,non_zero_id] ) )
    x[,i] = rnorm(n, theta_i, noiseVar[i])
    if( length(outlier_node) == 1 ){
      if( any(i == outlier_node) ){
        #x[sample(1:n, b),i] = rnorm(b, 10^3, noiseVar[i])
        x[1:b,i] = rnorm(b, 10^3, noiseVar[i])
      }
    }else if(length(outlier_node) ==p){
      x[seq(i,n, by = p),i] = rnorm(b, 10^3, noiseVar[i])
    }else if(length(outlier_node) > 1){
      temp_para = TRUE
      if( any(i == outlier_node) ){
        if( temp_para ){
          x[1:b,i] = rnorm(b, 10^3, noiseVar[i])
          temp_para = FALSE
        }else{
          x[sample(1:200, b),i] = rnorm(b, 10^3, noiseVar[i])
        }
      }
    }
  }
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B)
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  return( DAGsamples)
}



### PGSEM generator function ###
GSEM_generator_PGSEM = function( n, p, d = 2, beta_min, beta_max, graph_type, q = 0.025, structure = NULL, seed = 1, path = NULL, tau = 1){
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q = q, seed = seed, structure )$B
  
  ### all positives 
  B = abs(B)
  ### Sample with respect to the given distributions ###
  beta0 = runif(p, beta_min, beta_max)
  #beta0 = runif(p, 2, 2)
  x= matrix( 0, n, p)
  
  ### Nosie variance for Gaussian ###
  noiseVar = rep(0, p)
  
  for(j in 1:p){
    eta_j = colSums( B[j,1:(j-1) ]*t( x[,1:(j-1)] ) )
    theta_j = beta0[j] + eta_j 
    var_j = (tau^2)*abs(theta_j)
    x[,j] = theta_j + rnorm(n, 0, sqrt(var_j) )
  }
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("PGSEM_d",d,"_seed", seed, ".Rdata")
    save(DAGsamples, file = filename)
  }
  return(DAGsamples)
}



### GSEM generator function ###
GSEM_generator_LSEM = function( n, p, d =2, beta_min, beta_max, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1,path = NULL){
  
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed, structure )$B
  
  ### Sample with respect to the given distributions ###
  theta_i = 0 
  x= matrix( 0, n, p)
  
  ### Nosie variance for Gaussian ###
  noiseVar = rep(0, p)
  
  for( i in 1:p){
    eta_i = colSums( B[i,1:(i-1) ]*t( x[,1:(i-1)] ) )
    theta_i = eta_i
    if( i %% 4 == 1){
      var_max = 0.50
      x[,i] = theta_i + rnorm(n, 0, sqrt(var_max) )
      noiseVar[i] = var_max
    }else if( i %% 4 ==2){
      U_range = 1.25
      x[,i] = theta_i + runif(n, -U_range, U_range)
      noiseVar[i] = U_range^2/3
    }else if( i %% 4 ==3){
      T_df = 6
      x[,i] = theta_i + rt(n, 6)/sqrt(3)
      noiseVar[i] = T_df/ (T_df-2)/3
    }else if( i %% 4 ==0){
      x[,i] = theta_i + 2 * rbeta(n, 0.5, 0.5, ncp = 0) - 1
      noiseVar[i] = 1/2
    }
    #x[,i] = theta_i + runif(n, -U_range2, -U_range1)
    #x[,i] = theta_i +  rlaplace(n, m=0, s=1)
    #x[,i] = theta_i + rt(n, T_df)/2
    # x[,i] = theta_i + runif(n, -U_range2, -U_range1) + runif(n, U_range1, U_range2)
    # x[,i] = theta_i + runif(n, -U_range2, -U_range1) + runif(n, U_range1, U_range2)
    #x[,i] = theta_i + rnorm(n, 0, noiseVar[i]) + rnorm(n, -mu, noiseVar[i]) + rnorm(n, mu, noiseVar[i])
    #
    #x[,i] = theta_i +  rlaplace(n, m=-2, s= 0.5) +  rlaplace(n, m= 2, s=0.5)
    #
    #x[,i] = theta_i + rt(n, T_df)/2 + rlaplace(n, m=0, s=0.5)
  }
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  return( DAGsamples)
}

#### GSEM generator function ####
TSEM_generator = function( n, p, d =2, beta_min, beta_max, graph_type, q = 0.025, h = 1, seed = 1, structure = NULL, path = NULL){
  
  ### Generating data with p, n ###
  B = B2 = matrix( 0, p ,p )
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed, structure )$B
  noiseVar = rep(0, p)
  ### Sample with respect to the given distributions ###
  Time = proc.time()[3]
  x= matrix( 0, n, p)
  i = 1
  x[,i] = rt(n, df = 10)
  for( i in 1:p){
    non_zero_id = which( B[i,1:(i-1) ] != 0)
    theta_i = colSums( B[i,non_zero_id ]*t( x[,non_zero_id] ) )
    if( i %%3 ==0){
      x[,i] = theta_i + rt(n, df = 20) * sqrt(2)
      noiseVar[i] = 10/(10-2)
    }else if( i %%3 ==1){
      x[,i] = theta_i + rt(n, df = 10)* sqrt(2)
      noiseVar[i] = 15/(15-2)
    }else if( i %%3 ==2){
      x[,i] = theta_i + rt(n, df = 15)* sqrt(2)
      noiseVar[i] = 20/(20-2)
    }
     # variance is T_df/ (T_df-2)/4
  }
  # faster version #
  #x = rmvnorm(n, mean = rep(0, p), sigma = ( solve( diag(p) - t(B) ) %*% diag(noiseVar) %*% solve( diag(p) - B) ) )
  Time = proc.time()[3] - Time
  Time
  
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  } 
  return( DAGsamples)
}



### Sub-GSEM generator function ###
SubGSEM_generator = function(n, p, d =2, var_Min = 1, var_Max = 1, beta_min, beta_max, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1, path = NULL){
  
  #install.packages('truncnorm')
  library(truncnorm)
  
  ### Generating data with p, n ###
  B = B2 = matrix( 0, p ,p )
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed, structure )$B
  
  ### Sample with respect to the given distributions ###
  noiseVar = runif(p, min = var_Min, max = var_Max)
  
  ## link function: exponential; it can be any valid function
  x= matrix( 0, n, p)
  
  U_range = 2.5
  i = 1
  x[,i] = runif(n, -U_range, U_range)
  noiseVar[i] = U_range^2/3
  for( i in 2:p){
    eta_i = colSums( B[i,1:(i-1) ]*t( x[,1:(i-1)] ) )
    eta_i2 = colSums( B2[i,1:(i-1) ]*t( x[,1:(i-1)]^2 ) )
    beta_0 = 0
    theta_i = eta_i + eta_i2 + beta_0
    if( i %%3 ==0){
      x[,i] = rnorm(n, theta_i, noiseVar[i])
      noiseVar[i] = noiseVar[i]^2
    }else if( i %%3 == 1){
      x[,i] = theta_i + rtruncnorm(n, a=-2.5, b=2.5, mean = 0, sd = sqrt(10))
      noiseVar[i] = 1.915
    }else{
      x[,i] = theta_i + runif(n, -U_range,  U_range )
      noiseVar[i] = U_range^2/3
    }
  }
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )

  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  
  return( DAGsamples)
}

#### GSEM generator function ####
GSEM_generator = function( n, p, d =2, var_Min = 1, var_Max = 1, beta_min, beta_max, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1, path = NULL){
  
  ### Generating data with p, n ###
  B = B2 = matrix( 0, p ,p )
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed, structure )$B
  
  ### Sample with respect to the given distributions ###
  noiseVar = runif(p, min = var_Min, max = var_Max)
  
  Time = proc.time()[3]
  x= matrix( 0, n, p)
  i = 1
  x[,i] = rnorm(n, 0, noiseVar[i])
  for( i in 2:p){
    non_zero_id = which( B[i,1:(i-1) ] != 0)
    eta_i = colSums( B[i,non_zero_id ]*t( x[,non_zero_id] ) )
    #eta_i = colSums( B[i,1:(i-1) ]*t( x[,1:(i-1)] ) )
    #eta_i2 = colSums( B2[i,1:(i-1) ]*t( x[,1:(i-1)]^2 ) )
    #eta_i3 = colSums( B3[i,1:(i-1) ]*t( x[,1:(i-1)]^3 ) )
    beta_0 = 0
    theta_i = eta_i
    x[,i] = rnorm(n, theta_i, noiseVar[i])
  }
  # faster version #
  #x = rmvnorm(n, mean = rep(0, p), sigma = ( solve( diag(p) - t(B) ) %*% diag(noiseVar) %*% solve( diag(p) - B) ) )

  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  
  return( DAGsamples)
}

#### GSEM generator function with controlled eigenvalues ####
GSEM_generator2 = function( n, p, d =2, var_Min = 1, var_Max = 1, beta_min, beta_max, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1, path = NULL){
  
  ### Generating data with p, n ###
  B = B2 = matrix( 0, p ,p )
  
  for(ii in 1:(p%/%25)){
    B[1:25+25*(ii-1), 1:25+25*(ii-1)] = Graph_Generator( p = 25, d, graph_type, beta_min, beta_max, q, h, seed+100*(ii-1), structure )$B
  }
  
  ### Sample with respect to the given distributions ###
  noiseVar = runif(p, min = var_Min, max = var_Max)
  
  Time = proc.time()[3]
  x= matrix( 0, n, p)
  i = 1
  x[,i] = rnorm(n, 0, noiseVar[i])
  for( i in 2:p){
    non_zero_id = which( B[i,1:(i-1) ] != 0)
    eta_i = colSums( B[i,non_zero_id ]*t( x[,non_zero_id] ) )
    #eta_i = colSums( B[i,1:(i-1) ]*t( x[,1:(i-1)] ) )
    #eta_i2 = colSums( B2[i,1:(i-1) ]*t( x[,1:(i-1)]^2 ) )
    #eta_i3 = colSums( B3[i,1:(i-1) ]*t( x[,1:(i-1)]^3 ) )
    beta_0 = 0
    theta_i = eta_i
    x[,i] = rnorm(n, theta_i, noiseVar[i])
  }
  # faster version #
  #x = rmvnorm(n, mean = rep(0, p), sigma = ( solve( diag(p) - t(B) ) %*% diag(noiseVar) %*% solve( diag(p) - B) ) )
  
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  
  return( DAGsamples)
}


#### GSEM generator function ####
USEM_generator = function( n, p, d =2, beta_min, beta_max, U_lower = -1, U_upper = 1, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1, path = FALSE){
  
  ### Generating data with p, n ###
  B = B2 = matrix( 0, p ,p )
  B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed = seed, structure )$B
  
  ### Sample with respect to the given distributions ###
  Time = proc.time()[3]
  x= matrix( 0, n, p)
  i = 1
  x[,i] = runif(n, U_lower, U_upper)
  for( i in 2:p){
    non_zero_id = which( B[i,1:(i-1) ] != 0)
    eta_i = colSums( B[i,non_zero_id ]*t( x[,non_zero_id] ) )
    beta_0 = 0
    theta_i = eta_i
    x[,i] = theta_i + runif(n, U_lower, U_upper)
  }
  Time = proc.time()[3] - Time
  Time
  
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = rep(1/12 *(U_upper - U_lower)^2, p) )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  
  return( DAGsamples)
}





#### Sensitivity Analysis for Layer-wise learning  ####
# graph_type = 8; h = 5
# var_Min = sqrt(0.2);var_Max = sqrt(0.8)

GSEM_generator_Sensitivity = function( n, p, d =2, var_Min = 1, var_Max = 1, dist = "Gaussian", beta_min, beta_max, graph_type, q = 0.025, h = 1, structure = NULL, seed = 1, path = NULL){
  
  ### Generating data with p, n ###
  # p = 100;seed = 1

  trueB = matrix(1, p, p); est_B = matrix(0, p, p); seed2 = seed
  while( any(trueB!= est_B) ){
    set.seed(seed2)
    B = B2 = matrix( 0, p ,p )
    B = Graph_Generator( p, d, graph_type, beta_min, beta_max, q, h, seed2 )$B
    trueB = B; trueB[trueB!=0] = 1
    
    ### Sample with respect to the given distributions ###
    noiseVar = runif(p, min = var_Min, max = var_Max)
    
    ### Identifiability Check ###
  
    Ordering = rep(0, p)
    Rem = 1:p
    for(index in p:2){
      Omega =  t(diag(length(Rem))-B[Rem,Rem]) %*% diag(noiseVar[Rem]) %*% (diag(length(Rem))-B[Rem,Rem])
      Ordering[index] <- Rem[which.min( diag(Omega) )]
      Rem = setdiff(Rem, Ordering[index])
    }
    Ordering[1] <- Rem
    
    est_B = matrix(0, p, p)
    Rem = Ordering
    for(index in p:2){
      Omega =  t(diag(length(Rem))-B[Rem,Rem]) %*% diag(noiseVar[Rem]) %*% (diag(length(Rem))-B[Rem,Rem])
      est_B[Ordering[index], Ordering[1:(index-1)]] <- (Omega[length(Rem),-length(Rem)] != 0)
      Rem = setdiff(Rem, Ordering[index])
    }
    seed2 = seed2 + 100
  }
  ###
  
  Time = proc.time()[3]
  x= matrix( 0, n, p)
  i = 1
  x[,i] = rnorm(n, 0, noiseVar[i])
  for( i in 2:p){
    non_zero_id = which( B[i,1:(i-1) ] != 0)
    eta_i = colSums( B[i,non_zero_id ]*t( x[,non_zero_id] ) )
    #eta_i = colSums( B[i,1:(i-1) ]*t( x[,1:(i-1)] ) )
    #eta_i2 = colSums( B2[i,1:(i-1) ]*t( x[,1:(i-1)]^2 ) )
    #eta_i3 = colSums( B3[i,1:(i-1) ]*t( x[,1:(i-1)]^3 ) )
    beta_0 = 0
    theta_i = eta_i
    x[,i] = rnorm(n, theta_i, noiseVar[i])
  }
  # faster version #
  #x = rmvnorm(n, mean = rep(0, p), sigma = ( solve( diag(p) - t(B) ) %*% diag(noiseVar) %*% solve( diag(p) - B) ) )
  
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B, noise = noiseVar )
  
  ### path ###
  if( !is.null(path) ){
    setwd(path)
    filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
    save(DAGsamples, file = filename)
  }
  
  ### path ###
  # if( path ){
  #   #setwd("C:/2019GaussianSEM_SimulationResult/Data/Homo")
  #   setwd("D:/2021DGSEM_Simulation/Data")
  #   GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  #   #setwd("D:/2021DGSEM_Simulation/Data_DenseGraph")
  #   #GaussianDAG_filename = paste0("GaussianSEM_d",h,"_seed",seed,".Rdata")
  #   save(DAGsamples, file = GaussianDAG_filename)
  # }
  return( DAGsamples)
}
