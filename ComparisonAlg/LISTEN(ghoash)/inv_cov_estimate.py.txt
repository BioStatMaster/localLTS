import numpy as np
import numpy.linalg as la
import cvxopt as cx


def inv_cov_estimate_clime(Sm, reg_param, threshold=None, average=True, solver='glpk'):
    Omega_hat = cvx_clime(Sm, reg_param, average=average, solver=solver)
    if threshold is None:
        M = la.norm(Omega_hat, 1)
        threshold = 4 * M * reg_param
    Omega_hat[np.abs(Omega_hat) < threshold] = 0

    return Omega_hat


def cvx_clime_update(Sigma, O, tv, reg_param, col_sum=None, solver='glpk'):
    """
    :param O: Precision matrix.
    :param tv: Terminal vertex.
    :param reg_param:
    :param col_sum: absolute column sum of O.
    :param solver: LP solver.
    :return:
    """
    M = col_sum.max()
    col_sum[tv] = 0
    epsilon = 4 * M * reg_param
    tv_pars = np.abs(O[tv]) > epsilon   # Parents of tv.
    tv_pars[tv] = False

    if np.count_nonzero(tv_pars) == 0:
        O[tv, :] = 0
        O[:, tv] = 0
        return O

    for v in np.where(tv_pars)[0]:
        Sv = np.abs(O[v]) > epsilon
        Sv[tv] = False
        Sv |= tv_pars
        if np.count_nonzero(Sv) == 0:
            continue

        Sm = Sigma[Sv, :][:, Sv]
        n = Sm.shape[0]  # number of variables.
        sI = cx.spmatrix(np.ones(n), np.arange(n), np.arange(n), (n, n))         # I.
        smI = cx.spmatrix(np.ones(n) * -1, np.arange(n), np.arange(n), (n, n))   # - I.
        S, mS = cx.matrix(Sm), cx.matrix(-Sm)                                    # Sigma, - Sigma.
        sZ = cx.spmatrix([0], [0], [0], (n,n))                                   # zeros n x n.
        z = np.zeros(n)
        A = cx.sparse([[smI, sI, mS, S],
                       [smI, smI, sZ, sZ]])
        ln = np.ones(n) * reg_param
        b = cx.matrix(np.hstack((z, z, ln, ln)))
        c = cx.matrix(np.hstack((z, np.ones(n))))
        i = np.count_nonzero(Sv[:v])  # Index of vertex within support set Sv.
        beta = cvx_clime_get_column(n, c, A, b, i, solver=solver)
        beta = beta.flatten()
        beta[np.abs(beta) < epsilon] = 0.
        col_sum[v] -= np.sum(np.abs(O[v, Sv]))
        col_sum[v] += np.sum(np.abs(beta))
        O[v, Sv] = beta
        O[Sv, v] = beta

    O[tv, :] = 0
    O[:, tv] = 0
    return O


def cvx_clime_get_column(n, c, A, b, i, solver='glpk'):
    b[2*n + i] -= 1
    b[3*n + i] += 1
    res = cx.solvers.lp(c, A, b, solver=solver)
    b[2*n + i] += 1
    b[3*n + i] -= 1

    if res['status'] != 'optimal':
        print "Failed to optimize." + str(res)
        beta_i = np.zeros((n,1))
    else:
        beta_i = np.array(res['x'])
        beta_i = beta_i[:n]
    return beta_i


def cvx_clime(Sigma, reg_param, average=True, solver='glpk'):
    n = Sigma.shape[0]  # number of variables.
    sI = cx.spmatrix(np.ones(n), np.arange(n), np.arange(n), (n, n))         # I.
    smI = cx.spmatrix(np.ones(n) * -1, np.arange(n), np.arange(n), (n, n))   # - I.
    S, mS = cx.matrix(Sigma), cx.matrix(-Sigma)                              # Sigma, - Sigma.
    sZ = cx.spmatrix([0], [0], [0], (n,n))                                   # zeros n x n.
    z = np.zeros(n)
    A = cx.sparse([[smI, sI, mS, S],
                   [smI, smI, sZ, sZ]])
    ln = np.ones(n) * reg_param
    b = cx.matrix(np.hstack((z, z, ln, ln)))
    c = cx.matrix(np.hstack((z, np.ones(n))))

    O_hat = [cvx_clime_get_column(n, c, A, b, i, solver=solver)
             for i in xrange(n)]
    O_hat = np.hstack(O_hat)

    if average:
        return 0.5 * (O_hat + O_hat.T)
    else:
        # Take minimum.
        A, B = O_hat, O_hat.T
        R = np.zeros((n, n))
        for i in xrange(n):
            x, y = np.abs(A[i]), np.abs(B[i])
            R[i, x <= y] = A[i, x <= y]
            R[i, y < x] = B[i, y < x]
        return R


