import numpy as np
import numpy.linalg as la
from inv_cov_estimate import inv_cov_estimate_clime, cvx_clime_update
from functools import partial


def l2_loss(X, B):
    n = X.shape[0]  # Samples.
    Xh = X.dot(B.T)
    R = X - Xh
    r = np.einsum('i...,i...', R, R)
    l = np.sum(r)
    l /= float(n)
    return l


def log_lik(X, B, noise_vars=None):
    """
    Log likelihood of the directed graph. 
    :param X : data matrix.
    :param B : autoregression matrix (weight matrix).
    :param noise_vars : noise variances.
    """
    ll = 0
    n, p = X.shape
    for i in xrange(p):
        m = X.dot(B[i])
        r = X[:,i] - X.dot(B[i])
        if noise_vars is None:
            s = r.std()
        else:
            s = noise_vars[i]
        ll += -(0.5 * n * np.log(2*np.pi) + n * np.log(s) + 0.5 * np.sum((r/s)**2))
    return ll


def log_lik_u(X, O):
    """
    Log likelihood of the undirected graph.
    :param X : data matrix.
    :param O : Precision matrix.
    """
    Z = X - np.mean(X, 0)
    ll = 0
    n, d = X.shape
    for z in Z:
        ll += z.dot(O).dot(z)
    ll *= -0.5
    _, ld = la.slogdet(O)
    ll += n * 0.5 * (ld - d * np.log(2 * np.pi))
    return ll


def bic_score(l, B, samples):
    A = np.abs(B) > 0.001
    edges = np.count_nonzero(A)
    # bic = np.log(samples) * edges + samples * np.log(l)
    bic = l - edges * np.log(samples) * 0.5
    return bic


def update_inv_cov_clime(reg_param, col_sum, Sm, Oh, tv, solver='glpk'):
    cvx_clime_update(Sm, Oh, tv, reg_param, col_sum=col_sum, solver=solver)


def update_inv_cov_schur(Sm, Oh, tv):
    Oh -= 1./Oh[tv, tv] * np.outer(Oh[tv], Oh[tv])
    Oh[tv, :] = 0
    Oh[:, tv] = 0


def learn_weight_matrix(Sm, Om, reg_param,
                        update_method='clime',
                        solver='glpk'):
    """
        Learn the causal ordering, given the empirical
        covariance marix and the estimated precision matrix.
    :param Sm: Empirical covariance matrix.
    :param Om: Estimated precision matrix.
    :param reg_param : regularization parameter (only used for @update_method = 'clime'.
    :param update_method : Method to use to update the precision matrix when a
        terminal vertex is removed. Choices are 'clime' and 'schur'. If method is clime,
        then update a submatrix of the precision matrix using clime. If schur, then perform
        a schur complement update of the precision matrix.
    :param solver: LP solver.
    :return: B, vars : weight matrix and noise variance estimates.
    """
    n = Sm.shape[0]
    B = np.zeros((n, n))         # This is the SEM weight matrix.
    O = Om.copy()                # Keep a copy of the inv covariance matrix to update.
    vars = np.zeros(n)
    if update_method == 'clime':
        col_sum = np.sum(np.abs(O), 0)  # Absolute column sum. Used to calculate running L1 induced norm.
        update_inv_cov = partial(update_inv_cov_clime, reg_param, col_sum, solver=solver)
    elif update_method == 'schur':
        update_inv_cov = update_inv_cov_schur
    else:
        raise Exception("Invalid update method. Acceptable values are: clime, schur.")

    for i in xrange(n-1):
        tv = np.argmin(np.diag(O))  # Terminal vertex.
        vars[tv] = 1./O[tv,tv]
        b_tv = - (O[tv]/O[tv,tv])
        b_tv[tv] += 1.
        B[tv] = b_tv.copy()
        update_inv_cov(Sm, O, tv)
        O[tv, tv] = np.inf

    # Update the noise variance of the last vertex.
    tv = np.argmin(np.diag(O))
    vars[tv] = 1./O[tv, tv]

    return B, vars


def learn_sem(Sigma, reg_param,
              update_method='clime',
              rho=0, 
              solver='glpk'):
    """
        Learn a SEM from data.
    :param Sigma: Sample covariance matrix.
    :param reg_param:
    :param update_method : Method to use to update the precision matrix when a
        terminal vertex is removed. Choices are 'clime' and 'schur'. If method is clime,
        then update a submatrix of the precision matrix using clime. If schur, then perform
        a schur complement update of the precision matrix.
    :param rho: If rho > 0 then add rho * I to the covariance matrix to make it positive definite.
    :param solver: LP solver to use, defaults to GLPK.
    :return: An weight matrix (autoregression matrix) and noise variances.
    """
    
    try:
        Sm = Sigma
        if rho > 0:
            Sm[np.diag_indices_from(Sm)] += rho

        O_hat = inv_cov_estimate_clime(Sm, reg_param, solver=solver)

    except Exception as e:
        print e.message
        print "Failed to learn inverse covariance matrix."
        return None

    try:
        B, vars = learn_weight_matrix(Sm, O_hat, reg_param,
                                      update_method=update_method, 
                                      solver=solver)
        #t = 4 * la.norm(O_hat, 1) * reg_param
        #B[np.abs(B) < t] = 0
        return B, vars

    except Exception as e:
        print e.message
        print "Failed to learn DAG order."
        return None
