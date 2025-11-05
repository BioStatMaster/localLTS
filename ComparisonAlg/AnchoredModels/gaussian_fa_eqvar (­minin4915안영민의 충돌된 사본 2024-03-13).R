
#################################################################
#### zhang 2017 gaussian equal variance with Factor Analysis ####
#################################################################

FA_EQVAR_simulation_fun = function(seed, n_real, p_real, nleaf, alpha, path = F ){
  #This algorithm is for simulations using anchored linear SEM Learning_algorithm. 
  
  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################
  
  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$z[1:n_real,1:p_real] ## It should be z ##
  result_FA_EQVAR = FA_EQVAR(data, alpha = alpha, nleaf = nleaf, graph = graph)
  return(result_FA_EQVAR)
}


###############################################
#### identifiable 하기 위해 필요한 리프 수 ####
###############################################

ident_leaf = function(n){
  c_n = ((8*n+1)^(1/2)-1)/ (2*n)
  l_n = ceiling(c_n*n)
  return(list(per = c_n, min = l_n))
}

######################
## X : 열이 variance 행이 sample인 데이터, l : 리프노드 수
## output은 아직 설정하지 않았습니다.

FA_EQVAR =  function(data, alpha=0.001, nleaf, graph=NULL){  
  library(pcalg)
  library(Matrix)
  library(corpcor)
  library(ppcor)

    ml = ident_leaf(ncol(data))$min
  #print(ml)

  if(nleaf==0){ ## 0이면 true leaf 또는 true leaf가 identifiable하기 위한 minimum leaf 넣어줌
    B = graph
    B[B!=0] =1
    nl = length(which(colSums(B == 0) == nrow(B)))# true leaf 넣어주기
    if(nl<ml){
      nleaf = ml
    }else{
      nleaf = nl
    }
  }else if(nleaf<ml){nleaf = ml}
  
  X = as.matrix(na.omit(data))
  #### FA 수행 (요인의 수 K 설정, 전체 variable 수 n에서 리프노드 수 l 빼기)
  K <- ncol(X) - nleaf
  
  Runtime = proc.time()[3]
  
  result <- factanal(X, factors = K, scores = "regression", rotation = "varimax")
  print('factor finish')
  # 결과 시각화
  # FA 결과 correlation에 관련. covariance 보기 위해 보정
  s = diag(var(X)) # observed variable들의 분산
  Estar_var =(result$uniquenesses)*s
  #plot(Estar_var, type = "b", main = "var(E*)") # E* 분산 plot
  
  L = diag(sqrt(s))%*%result$loadings
  #cal = L %*% t(L) # FA 결과 loadings로 나머지 분산 계산
  Xstar_var = cov(X)-diag(Estar_var) # observed variable 분산에서 E* 분산 빼서 계산
  
  leaf <- order(Estar_var, decreasing = TRUE)[1:nleaf] # E*분산이 큰 nleaf개를 leaf node로
  sub_Xstar_var = Xstar_var[-leaf,-leaf] # non-leaf node들의 분산. PC알고리즘 적용해야함
  print(sub_Xstar_var)
  #### non-leaf node의 CPDAG 찾기 (By PC algorithm)
  
  sub_Xstar_cor = cov2cor(sub_Xstar_var)
  suffStat = list(C = sub_Xstar_cor, n= nrow(X))
  pc.result = pc(suffStat, indepTest = gaussCItest, p = ncol(sub_Xstar_cor),alpha = alpha)
  #plot(pc.result,main="")
  
  adj_mat <- as.matrix(as(pc.result,"amat"))
  non_leaf = as.numeric(colnames(adj_mat))
  #graph1 <- graph_from_adjacency_matrix(t(adj_mat), mode = "directed")
  #plot(graph1, main = "Graph from PC algorithm")
  print('end without leaf')
  
  #### leaf node의의 parent 찾기####
  for (i in 1:length(leaf)){
    dep = rep(0,K)
    # partial correlation을을 통해 조건부 독립 확인
    leaf_data=X[,leaf[i]]
    for (k in 1:K){
      nonleaf_data = X[,non_leaf[k]]
      result <- pcor.test(leaf_data, nonleaf_data, X[,non_leaf[-k]], method = "spearman")

      if (!is.nan(result$p.value) & (result$p.value<=alpha)){
        dep[k]=1
      }
    }
    adj_mat = rbind(adj_mat,dep)
    rownames(adj_mat)[nrow(adj_mat)] = leaf[i]

  }
  # 비어있는 부분 adjacency matrix 만들기(leaf->nonleaf 가는 엣지 없으니까)
  zeros = matrix(0,nrow(adj_mat),nrow(adj_mat)-ncol(adj_mat))
  adj_mat = cbind(adj_mat,zeros)
  colnames(adj_mat) = rownames(adj_mat)

  #graph2 <- graph_from_adjacency_matrix(t(adj_mat), mode = "directed")
  #plot(graph2, main = "Graph from Adjacency Matrix")
  
  Runtime = proc.time()[3]- Runtime
  ####
  print(paste("Latent_PC_GSEM takes: ", Runtime))
  print(!is.null(graph))
  if( !is.null(graph) ){
    B = graph
    B[B!=0] =1
    evaluation_result_FA = NULL
    evaluation_result_FA_MEC = evaluation_fun( t(dag2cpdag(t(B))), adj_mat)
    evaluation_result_FA2 = NULL
  }else{
    evaluation_result_FA = NULL
    evaluation_result_FA_MEC = NULL
    evaluation_result_FA2 = NULL
  }
  return(
    list( DAG_Evaluation = evaluation_result_FA, 
          MEC_Evaluation = evaluation_result_FA_MEC,  
          Oracle_Evaluation = evaluation_result_FA2, 
          DAG = adj_mat,
          Ordering = NULL, 
          Time = Runtime)
  )
}



########################################
#### simulation data 생성(Gaussian) ####
# ########################################
# 
# N <- 20 # 노드 수
# sample_n <- 5000 # 샘플 수
# l <- 8 # 리프노드 수
# err_var = 0.25 # measurement error variance
# 
# # edge weight matrix
# B1 <- 1.2 * (matrix(runif(N * N), nrow = N, ncol = N) - 0.5)
# for (i in 1:N) {
#   B1[i, i:N] <- 0
# }
# B1[(N - l + 1):N, (N - l + 1):N] <- 0
# 
# # A_m = (I-B)^(-1)
# A_m <- solve(diag(N) - B1)
# 
# Err <- matrix(rnorm(N * sample_n), nrow = N, ncol = sample_n)
# 
# # Latent variable (True data) 생성
# X_tilde <- A_m %*% Err
# 
# # measurement error 생성 및 observed variable 계산
# noise_std <- rep(sqrt(err_var), N) # measurement error std(표준편차) 설정
# Noise <- diag(noise_std) %*% matrix(rnorm(N * sample_n), nrow = N, ncol = sample_n)
# X <- X_tilde + Noise # 행: variable, 열: sample (20 * 5000 matrix)
# 
# # identifiable 위해 필요한 리프 수
# leaf_need = ident_leaf(N)
# print(leaf_need$per)
# print(leaf_need$min)
# print(is.matrix(X))
# # 그래프 추출
# FA_EQVAR(t(X),alpha=0.001,nleaf=l,graph=NULL)
