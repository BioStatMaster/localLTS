## Package
library(tensorA) # Einsum
library(purrr)   # Partial ftn

## Import Part
source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/inv_cov_estimate2.R'))

## l2_loss ftn

l2_loss = function(X,B)
{
  n = dim(X)[1]
  Xh = X %*% t(B)
  R = X- Xh
  r = R %e% R
  l = sum(r)
  l = l/n
  return(l)
}

## log_lik ftn
log_lik = function(X,B, noise_vars = NULL)
{
  # Log likelihood of the directed graph. 
  # :param X : data matrix.
  # :param B : autoregression matrix (weight matrix).
  # :param noise_vars : noise variances.
  
  ll = 0
  n = dim(X)[1]
  p = dim(X)[2]
  for (i in 1:p)
  {
    m = X %*% B[i,]
    r = X[,i] - X %*% B[i,]
    if(is.null(noise_vars)) s = apply(r, 2, function(x) sd(x) * sqrt((length(x) - 1) / length(x)) )
      else s = noise_vars[i,] # check
  }
  ll = ll -(0.5 * n * log(2*pi) + n * log(s) + 0.5 * sum((r/s)^2))
  return(ll)
}

## log_lik_u ftn
log_lik_u = function(X, O)
{
  # Log likelihood of the undirected graph.
  # :param X : data matrix.
  # :param O : Precision matrix.
  Z = X - apply(X,2,mean)
  ll = 0
  n = dim(X)[1]
  d = dim(X)[2]
  for (z in Z)
  {
    ll = ll + z %*% O %*% z
    ll = ll * -0.5
    ld = determinant(O)$modulus[1]
    ll = ll + n* 0.5 * (ld - d * log(2*pi))
    return(ll)
  }
}

## bic_score ftn
bic_score = function(l,B, samples)
{
  A = abs(B) > 0.001
  edges = sum(A != 0)
  # bic = log(samples) * edges + samples *log(l)
  bic = l - edges * log(samples) * 0.5
  return(bic)
}

## update_inv_cov_clime ftn
update_inv_cov_clime = function(reg_param, col_sum, Sm, Oh, tv, solver='glpk')
{
  cvx_clime_update(Sm, Oh, tv, reg_param, col_sum=col_sum, solver=solver)
}
 
## update_inv_cov_schur ftn 
update_inv_cov_schur = function(Sm, Oh, tv)
{
  Oh = Oh - 1 / Oh[tv,tv] * outer(Oh[tv,], Oh[tv,]) # check
  Oh[tv,] = 0
  Oh[,tv] = 0
}

## learn_weight_matrix ftn

learn_weight_matrix = function(Sm, Om, reg_param,
                               update_method='clime',
                               solver='glpk')
{
  # Learn the causal ordering, given the empirical
  # covariance marix and the estimated precision matrix.
  # :param Sm: Empirical covariance matrix.
  # :param Om: Estimated precision matrix.
  # :param reg_param : regularization parameter (only used for @update_method = 'clime'.
  #                                              :param update_method : Method to use to update the precision matrix when a
  #                                              terminal vertex is removed. Choices are 'clime' and 'schur'. If method is clime,
  #                                              then update a submatrix of the precision matrix using clime. If schur, then perform
  #                                              a schur complement update of the precision matrix.
  #                                              :param solver: LP solver.
  #                                              :return: B, vars : weight matrix and noise variance estimates.
  n = dim(Sm)[1]
  B = matrix(0,n,n)  # This is the SEM weight matrix.
  O = Om            # Keep a copy of the inv covariance matrix to update.
  vars = rep(0,n)
  if (update_method == 'clime')
  {
    col_sum = apply(abs(O),2,sum) # Absolute column sum. Used to calculate running L1 induced norm.
    update_inv_cov = partial(update_inv_cov_clime, reg_param, col_sum, solver=solver) # Check
  }else if(update_method == 'schur')
    {update_inv_cov = update_inv_cov_schur
  }else {message("Exception(Invalid update method. Acceptable values are: clime, schur.)")}
  
  for (i in 1:(n-1))
  {
    tv = which.min(diag(O))  # Terminal vertex.
    vars[tv] = 1/O[tv,tv] # Check
    b_tv = -(O[tv,]/O[tv,tv]) # Check
    b_tv[tv] = b_tv[tv] + 1
    B[tv] = b_tv
    update_inv_cov(Sm, O, tv)
    O[tv, tv] = Inf # check
    
    # Update the noise variance of the last vertex.
    tv = which.min(diag(O))
    vars[tv] = 1/O[tv,tv] # check
    
    return(list(B = B, vars = vars))
  }
}

## learn_sem ftn

learn_sem = function(Sigma, reg_param,
                     update_method='clime',
                     rho=0, 
                     solver='glpk')
{
  # Learn a SEM from data.
  # :param Sigma: Sample covariance matrix.
  # :param reg_param:
  #   :param update_method : Method to use to update the precision matrix when a
  # terminal vertex is removed. Choices are 'clime' and 'schur'. If method is clime,
  # then update a submatrix of the precision matrix using clime. If schur, then perform
  # a schur complement update of the precision matrix.
  # :param rho: If rho > 0 then add rho * I to the covariance matrix to make it positive definite.
  # :param solver: LP solver to use, defaults to GLPK.
  # :return: An weight matrix (autoregression matrix) and noise variances.
  
  tryCatch(
    {
      Sm = Sigma
      if (rho > 0){
        diag(Sm) = diag(Sm) + rho
        }
        O_hat = inv_cov_estimate_clime(Sm, reg_param, solver=solver)
    }
    error = function(e){
      message('Failed to learn inverse covariance matrix.')
      print(e)
      return(NULL)}
  )
  tryCatch(
    {
      B = learn_weight_matrix(Sm, O_hat, reg_param,
                              update_method=update_method, 
                              solver=solver)[[1]]
      vars = learn_weight_matrix(Sm, O_hat, reg_param,
                                 update_method=update_method, 
                                 solver=solver)[[2]]
      # t = 4 * apply(O_hat,1,norm)  * reg_param
      # B[abs(B) < t] = 0
      return (list(B = B,  vars = vars))
    }
    error = function(e){
      message('Failed to learn DAG order.')
      print(e)
      return(NULL)}
  )
}