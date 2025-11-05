inv_cov_estimate_clime = function(Sm, reg_param, threshold=None, average=True, solver='glpk'){
  Omega_hat = cvx_clime(Sm, reg_param, average=average, solver=solver)
  if(is.null(threshold)){
    M = norm(Omega_hat,type =  '1')
    threshold = 4 * M * reg_param
  }
  Omega_hat[abs(Omega_hat) < threshold , ] = 0
  return( Omega_hat)
}

cvx_clime_update = function(Sigma, O, tv, reg_param, col_sum=None, solver='glpk'){
  # """
  # :param O: Precision matrix.
  # :param tv: Terminal vertex.
  # :param reg_param:
  # :param col_sum: absolute column sum of O.
  # :param solver: LP solver.
  # :return:
  # """
    M = max(col_sum)
    col_sum[tv] = 0
    epsilon = 4 * M * reg_param
    tv_pars = abs(O[tv]) > epsilon   # Parents of tv.
    tv_pars[tv] = F

    if(sum(tv_pars!=0)== 0){
        O[tv, ] = 0
        O[, tv] = 0
        return(O)
    }

    for( v in which(tv_pars)){
        Sv = abs(O[v]) > epsilon
        Sv[tv] = F
        Sv = Sv|tv_pars
        if(sum(Sv!=0) == 0)next
            
        Sm = Sigma[Sv, Sv]
        n = nrow(Sm)  # number of variables.
        sI = diag(rep(1 ,n))         # I.
        smI = -diag(rep(1 ,n))   # - I.
        S = as.matrix(Sm); mS = as.matrix(-Sm)                                    # Sigma, - Sigma.
        sZ = matrix(0, ncol = n, nrow = n) # zeros n x n.
        z = rep(0,n)
        A = cbind(rbind(smI, sI, mS, S),
                  rbind(smI, smI, sZ, sZ))
        ln = rep(1,n) * reg_param # reg_param 이 np.array 의 형식으로 쓰여진거면 문제없음.
        b = as.matrix(cbind(z, z, ln, ln))
        c = as.matrix(cbind(z, rep(1,n)))
        i = sum(Sv[1:v] != 0)  # Index of vertex within support set Sv.
        beta = cvx_clime_get_column(n, c, A, b, i)
        beta = as.vector(t(beta))
        ##########
        col_sum[v] = col_sum[v] - sum(abs(O[v, Sv]))
        col_sum[v] = col_sum[v] + sum(abs(beta))
        O[v, Sv] = beta
        O[Sv, v] = beta
    }

    O[tv, ] = 0
    O[, tv] = 0
    return(O)
}


cvx_clime_get_column = function(n, c, A, b, i, solver='glpk'){
  b[2*n + i] = b[2*n + i] - 1
  b[3*n + i] = b[3*n + i] + 1
  res = Rglpk_solve_LP(obj = c,
                       mat = t(A), # Rgplk_solve_LP 함수에서 ox.solver.lp 와 다른 인풋요구 
                       rhs =  b,
                       dir = rep('<=',nrow(t(A)))
  ) # 결과저장형식 알아보기
  b[2*n + i] = b[2*n + i] + 1
  b[3*n + i] = b[3*n + i] - 1

  if(res$status != 0){ # 0 is optimal.
    print("Failed to optimize.")
    print(res)
    beta_i = matrix(0, nrow = n ,ncol = 1)
  }else{
    beta_i = res$solution
    beta_i = beta_i[1:n]
  }
  return(beta_i)
}

cvx_clime = function(Sigma, reg_param, average=True, solver='glpk'){
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
  
  if( average){
    return( 0.5 * (O_hat + O_hat.T))
  }else{
    # Take minimum.
    A = O_hat ; B = t(O_hat)
    R = matrix(0 , ncol = n ,nrow = n)
    for( i in 1:n ){
      x = abs(A[i,]) ; y = abs(B[i,])
      R[i, x <= y] = A[i, x <= y]
      R[i, y < x] = B[i, y < x]
    }
  return( R)
  }
}



