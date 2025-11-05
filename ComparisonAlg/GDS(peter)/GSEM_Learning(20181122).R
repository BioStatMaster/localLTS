data = read.csv("C:/Users/gwpar/Dropbox/Future Work/microdust/data.csv")

head(data)

X = data[-1,-1]
SigmaHat <- cov(X)
pars <- list(SigmaHat = SigmaHat)
resGDS <- GDS(X, scoreName = "SEMSEV", pars, check = "checkUntilFirstMinK", output = FALSE)


a<-0
p <- 15
n<-1000
d<- 2/(p-1)


n = 200
p = 20
d = 5
beta_min = 0.25
beta_max = 0.75
graph_type = 5
a = 0.0
noiseVar <- runif(p, min = 1-a, max = 1+a)

RawData = GSEM_generator( n, p, d, noiseVar, dist, beta_min, beta_max, graph_type, structure = NULL, seed = 1)

data = RawData$x
graph = RawData$true_Matrix
graph[graph!=0] = 1
GSEM_Algorithm(data, graph, lambda.lasso = 0.4, sparsity_level = 1)$DAG_Evaluation

sapply(data,var)

### GSEM generator function ###
GSEM_generator = function( n, p, d, noiseVar, dist = "Gaussian", beta_min, beta_max, graph_type, structure = NULL, seed = 1){

  set.seed(seed)
  ## Generating data with p, n ##
  prob = d/ (p-1)
  B = matrix( 0, p ,p )

  ## various type of graphs ##
  ## 1: tree, 2: bipartite, 3: cycle, 4: random , 5: random2, 6 random, 6 manual ## 
  if( graph_type == 1 ) { rGRAPH = randomTree(p, n_level = p , lB = beta_min, uB =  beta_max);	B = -abs( (wgtMatrix(rGRAPH)) ) }
  if( graph_type == 2 ) { rGRAPH = randomBipartite(p, lB = beta_min, uB =  beta_max);	B = wgtMatrix(rGRAPH) }
  if( graph_type == 3 ) { rGRAPH = randomCycle(p, lB = beta_min, uB =  beta_max);	B = wgtMatrix(rGRAPH) }
  #expected number of parents
  if( graph_type == 4 ) { rGRAPH = randomDAG(p, prob , lB = beta_min, uB =  beta_max);	B = (wgtMatrix(rGRAPH)) }
  #fixed number of parents
  if( graph_type == 5 ) {
    for(i in 2:p){
      parents_i = sample( 1:(i-1), min(i-1, d),replace = FALSE )
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
    }
  }
  #fixed number of parents and ordering
  if( graph_type == 6 ) {
    for(i in 2:p){
      if( i <= d ){
        B[i, 1:(i-1)] = runif( (i-1) , beta_min,  beta_max)* sample( c(-1,1), (i-1), replace = T )
      }
      if( i > d ){
        B[i, i-1] = runif( 1 , beta_min,  beta_max)* sample( c(-1,1), 1)
        parents_i = sample( 1:(i-2), min(i-1, d-1), replace = FALSE )
        B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
      }
    }
  }
  ##B is given
  if( graph_type == 7 ) {
    for(i in 2:p){
      parents_i = which( structure[i,]!=0 ) #B is given
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max)* sample( c(-1,1), length(parents_i), replace = T)
    }
  }
  #random number of parents bounded by d
  if( graph_type == 8 ) {
    for(i in 2:p){
      parents_i = sample( 1:(i-1), sample(1:min(i-1, d), 1),replace = FALSE )
      B[i,parents_i ] = runif( length(parents_i) , beta_min,  beta_max) * sample( c(-1,1), length(parents_i), replace = T)
    }
  }

  ### Sample with respect to the given distributions ###
  
  ## link function: exponential; it can be any valid function
  x= matrix( 0, n, p) 
  
  i = 1
  x[,i] = rnorm(n, 0, noiseVar[i])
  for( i in 2:p){
    eta_i = colSums( B[i,1:(i-1) ]*t( x[,1:(i-1)] ) )
    beta_0 = runif(1, 0.5, 1)
    theta_i = eta_i + beta_0
    x[,i] = rnorm(n, theta_i, noiseVar[i])
  }
  x = as.data.frame(x)
  DAGsamples = list(x= x, true_Matrix= B )
  
  # Path
  #GaussianDAG_filename = paste0("GaussianSEM_p",p,"_d",d,"_seed",seed,".Rdata")
  #setwd("C:/Users/Uos/Dropbox/bobo/2018GSEM/")
  
  #return( save(DAGsamples, file = GaussianDAG_filename) )
  return( DAGsamples)
}

#####GSEM Learning Algorithm #######################
GSEM_Algorithm = function(data, graph = NULL, lambda.lasso = 0.4, sparsity_level = 0.18){
  
  library(glmnet)
  tuning.parameter <- lambda.lasso
  
  ###################
  X = as.matrix( data )
  RemNode = 1:p 
  pi_GSEM = NULL
  ScoreMatrix = matrix(0, p, p)
  
  #### Finding Ordering ####
  Time_GSEM = proc.time()[3]
  
  #### ordering estimation 
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



