if(!require('Rglpk')){
  install.packages('Rglpk')
  library(Rglpk)
}

#inv_cov_estimate_clime = function(Sm, reg_param , threshold = 0.01, average = TRUE){
#  Omega_hat = cvx_clime(Sm, reg_param, average = average) 
#  Omega_hat[abs(Omega_hat) < threshold ] = 0
#  return(Omega_hat)
#}

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

cvx_clime_update = function(Sigma, O, tv, reg_param, epsilon=0.001){
  #:param O: Precision matrix.
  #:param tv: Terminal vertex.
  #:param reg_param:
  #:return:
  
  tv_pars = abs(O[tv,]) > epsilon   # Parents of tv.
  tv_pars[tv] = FALSE
  if(sum(tv_pars != 0) == 0) return(O)
  
  for(v in which(tv_pars)){ # ?
    Sv = abs(O[v,]) > epsilon
    Sv[tv] = FALSE
    Sv  = Sv|tv_pars 
    if( sum(Sv != 0) == 0) next
    
    Sm = Sigma[Sv, ][, Sv]
    n = nrow(Sm)  # number of variables.
    sI = diag(rep(1,n))    # I.
    smI = -diag(rep(1,n))   # - I.
    S = as.matrix(Sm);mS = as.matrix(-Sm) # Sigma, - Sigma.
    sZ = matrix(0, ncol = n, nrow = n) # zeros n x n.
    z = rep(0,n)
    A = cbind(rbind(smI, sI, mS, S),
                   rbind(smI, smI, sZ, sZ))
    ln = rep(1,n) * reg_param # reg_param 이 np.array 의 형식으로 쓰여진거면 문제없지만 형식 체ㅋ
    b = as.matrix(cbind(z, z, ln, ln))
    c = as.matrix(cbind(z, rep(1,n)))
    i = sum(Sv[1:v] != 0)  # Index of vertex within support set Sv.
    beta = cvx_clime_get_column(n, c, A, b, i)
    beta = as.vector(t(beta))
    O[v, Sv] = beta
    O[Sv, v] = beta
  }
  O[tv, ] = 0
  O[, tv] = 0
  
    return(O)
  }

cvx_clime_get_column = function(n, c, A, b, i){
  b[2*n + i] = b[2*n + i] - 1  # 
  b[3*n + i] = b[3*n + i] + 1
  res = Rglpk_solve_LP(obj = c,
           mat = t(A), # Rgplk구_solve_LP 함수에서 ox.solver.lp 와 다른 인풋요구 
           rhs =  b,
           dir = rep('<=',nrow(t(A)))
           ) # 결과저장형식 알아보기
  b[2*n + i] = b[2*n + i] + 1
  b[3*n + i] = b[3*n + i] - 1
  
  if(res$status != 0){ # 0 is optimal.
    message("Failed to optimize.")
    print(res)
    beta_i = rep(0,n)
  }else{
    beta_i = res$solution
    beta_i = beta_i[1:n]
  }
  return(beta_i)
}

cvx_clime = function(Sigma, reg_param, average=TRUE){
  # Disable messages.
  # cx.solvers.options['show_progress'] = False # 진행상황 띄우기 끄기
  # cx.solvers.options['glpk'] = {'msg_lev': 'GLP_MSG_OFF'} # turns off the screen output subsequent calls lp with the 'glpk' option.
  
  n = nrow(Sigma)  # number of variables.
  sI = diag(rep(1,n))        # I.
  smI = -diag(rep(1,n))   # - I.
  S = as.matrix(Sigma); mS = as.matrix(-Sigma)                              # Sigma, - Sigma.
  sZ = matrix(0,ncol = n , nrow = n)                                   # zeros n x n.
  z = rep(0 , n)
  A = cbind(rbind(smI, sI, mS, S),
                 rbind(smI, smI, sZ, sZ))
  ln = rep(0,n) * reg_param #reg_param 의 형태 
  b = cbind(z, z, ln, ln)
  c = cbind(z, rep(1,n))
  

  O_hat = list()
  for(i in 1:n){
    O_hat[[i]] = cvx_clime_get_column(n, c, A, b, i)
  }
  O_hat = do.call(cbind, O_hat)
  
  if(average){
    return(0.5 * (O_hat + t(O_hat)))
  }else{
    # Take minimum.
    A = O_hat ; B = t(O_hat)}
  
  R = matrix(0, ncol = n, nrow = n)
  for(i in (1:n)){
    x = np.abs(A[i,]) ; y = np.abs(B[i,])
    R[i, x <= y] = A[i, x <= y]
    R[i, y < x] = B[i, y < x]
  }
  return(R)
}
  
  
  