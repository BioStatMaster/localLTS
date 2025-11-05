import numpy as np
import numpy.linalg as la
from inv_cov_estimate import inv_cov_estimate_clime
import traceback as tb
from sklearn.model_selection import train_test_split
import itertools as it


# If the absolute value of the precision matrix
# is below this value, then assume that the entry is 0.
# Used for calculating support.
INV_COVARIANCE_THRESHOLD = 0.02


def ols_estimate(Sigma, S, i):
    """
    :param Sigma: Empirical covariance matrix.
    :param S: Support set.
    :param i: variable for which the ols estimate has to be computed.
    :return: ols estimate.
    """
    theta = la.inv(Sigma[S,:][:,S]).dot(Sigma[S,i])
    theta = np.array(theta).flatten()
    return theta


def get_markov_blanket(O, i, inv_cov_threshold=INV_COVARIANCE_THRESHOLD):
    x = np.abs(O[i]) > inv_cov_threshold
    x[i] = False
    mb = np.where(x)[0]
    return list(mb)


def get_dag_order(Sigma, O_hat,
                  inv_cov_threshold=INV_COVARIANCE_THRESHOLD):
    """
        Learn the causal ordering, given the empirical
        covariance marix and the estimated precision matrix.
    :param Sigma: Empirical covariance matrix.
    :param O_hat: Estimated precision matrix.
    :param inv_cov_threshold: value below which entries in
        the inverse covariance matrix are assumed to be 0.
    :return:
    """
    n = Sigma.shape[0]
    z = np.zeros(n, dtype=int)         # DAG order.
    vars = np.zeros(n, dtype=float)    # variances.
    r = np.zeros(n)                    # ratios.

    count = n
    O = O_hat.copy()

    # Compute markov blanket for all nodes.
    # This is in terms of node labels.
    mbs = [get_markov_blanket(O, i, inv_cov_threshold=inv_cov_threshold)
           for i in xrange(n)]

    for i in xrange(n):
        S = mbs[i]
        if len(S) > 0:
            theta_i = ols_estimate(Sigma, S, i)
            Sp = np.where(np.abs(theta_i) > 0.001)[0]
            if Sp.size > 0:
                o = O_hat[i,S]
                o = o[Sp]
                r[i] = np.max(np.abs(o/theta_i[Sp]))

    # print r

    while count >= 2:
        # print "============================="
        # print r
        # print "============================="

        tv = np.argmin(r)  # This is the node label.
        z[count - 1] = tv
        r[tv] = np.inf
        vars[tv] = Sigma[tv, tv]
        if O[tv, tv] > 0:
            vars[tv] = 1./O[tv, tv]

        # Delete terminal vertex from precision matrix.
        v = O[tv]
        if O[tv, tv] > 0:
            O -= np.outer(v,v)/O[tv, tv]
        else:
            O -= np.outer(v,v)
            O[tv, tv] = 0

        count -= 1

        # Recompute ratios for variables in markov blanket of tv.
        for j in mbs[tv]:
            # print "Recomputing ratio for:", j,
            mbs[j] = get_markov_blanket(O, j, inv_cov_threshold=inv_cov_threshold)
            # print "Updated markov blanket of ", j, "is:", mbs[j]

            r[j] = 0
            S = mbs[j]
            if len(S) > 0:
                theta_j = ols_estimate(Sigma, S, j)
                Sp = np.abs(theta_j) > 0.001
                if np.count_nonzero(Sp) == 0:
                    continue

                # print "theta = ", theta_j.round(5)
                O_j = O[j][S]

                # print "Omega_j = ", O_j
                r[j] = np.max(np.abs(O_j[Sp]/theta_j[Sp]))
            # print " = ", r[j]

    z[0] = np.argmin(r)
    vars[z[0]] = Sigma[z[0], z[0]]
    if O[z[0], z[0]] > 0:
        vars[z[0]] =  1./O[z[0], z[0]]

    return z, vars


def get_bn_from_order(Sigma, O_hat, dag_order, vars,
                      inv_cov_threshold=INV_COVARIANCE_THRESHOLD):
    """
        Learn from GBN given DAG order.
    :param Sigma: n x n covariance matrix.
    :param dag_order:
    :param inv_cov_threshold: value below which the entries
        in the inverse covariance matrix is assumed to be 0.
    :return: (adjacency matrix, weights, noise variances)
    """
    n = Sigma.shape[0]
    A = np.zeros((n, n))
    W = np.zeros((n, n))
    mbs = [get_markov_blanket(O_hat, i,
                              inv_cov_threshold=inv_cov_threshold)
           for i in xrange(n)]
    # position of a node in the dag order.
    pos_in_dag_order = np.argsort(dag_order)
    for i in xrange(1, n):
        v = dag_order[i] # Variable for which we
                         # want to parents and parameters.
        # print "learning parameters for node:", v
        S = filter(lambda l: pos_in_dag_order[l] < pos_in_dag_order[v],
                   mbs[v])
        theta_v = ols_estimate(Sigma, S, v)
        for j in xrange(theta_v.size):
            u = S[j]
            if abs(theta_v[j]) > inv_cov_threshold:
                A[v, u] = 1
                W[v, u] = theta_v[j]

    return (A, W, vars)


def learn_gbn(Sigma, reg_param,
              inv_cov_threshold=INV_COVARIANCE_THRESHOLD,
              use_marginal_variance_for_dag_order=False):
    """
        Learn a GBN from data.
    :param Sigma: Sample covariance matrix.
    :param reg_param:
    :param inv_cov_threshold:
    :return: An @GaussianBayesianNetwork instance.
    """
    # Sigma = quantize(Sigma, threshold=inv_cov_threshold)
    try:
        O_hat = inv_cov_estimate_clime(
                Sigma, reg_param, threshold=inv_cov_threshold
        )

    except Exception as e:
        print e.message
        print "Failed to learn inverse covariance matrix."
        return None

    try:
        if not use_marginal_variance_for_dag_order:
            dag_order, v = get_dag_order(Sigma, O_hat,
                                         inv_cov_threshold=inv_cov_threshold)
        else:
            # Use marginal variance to figure out DAG order.
            dag_order = list(np.argsort(np.diag(Sigma)))
            v = np.ones(Sigma.shape[0])
    except Exception as e:
        print e.message
        print "Failed to learn DAG order."
        tb.print_exc()
        return None

    # print "DAG order:", dag_order
    bn = get_bn_from_order(Sigma, O_hat, dag_order, v,
                           inv_cov_threshold=inv_cov_threshold)
    return bn


def standardize(X):
    # standardize data.
    m = X.mean(axis=0)
    s = X.std(axis=0)
    Z = (X - m)/s
    return Z


def select_model(X, reg_params, thresholds):
    """
        Split data into train and test. Learn multiple GBNs
        for all combinations regularization parameter and thresholds.
        Return the model that has the highest test loglikelihood.

        :param X: data (m x n), where m is #samples, n is #nodes.
        :param reg_params: 
        :param thresholds:
        :return: The best model.
    """

    Z = standardize(X)
    
    # split data into train and test.
    Z_train, Z_test = train_test_split(Z, test_size=0.5, random_state=17)

    print "training data size:", Z_train.shape
    print "test data size:", Z_test.shape

    # Compute sample covariance matrix.
    C_train = (1./Z_train.shape[0]) * Z_train.T.dot(Z_train)

    best_r, best_t, min_neg_log_lik = 0, 0, np.inf
    for r, t in it.product(reg_params, thresholds):
        l = np.inf
        bnh = learn_gbn(C_train, r, inv_cov_threshold=t)
        if bnh is not None:
            l = - bnh.log_likelihood(Z_test)
        print "reg_param: %f, threshold: %f, test negative log likelihood: %.3f" % (r, t, l)
        if l < min_neg_log_lik:
            best_r, best_t = r, t 
    
    C = (1./Z.shape[0]) * Z.T.dot(Z)
    bn = learn_gbn(C, best_r * np.sqrt(1./2), t)
    
    return bn, best_r, best_t

