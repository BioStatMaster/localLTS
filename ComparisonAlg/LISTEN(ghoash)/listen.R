source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/learn_sem.R'))

covariance  = function(X){ ### 계산 방법 오류가 있는 것같음.
  p = 0 # Number of variables.
  if( is.matrix(X) | is.data.frame(X) ){
    # X is an array. Compute covariance directly.
    n = nrow(X)
    mu = apply(X, 2 ,mean)
    C = t(X)%*%X - outer(mu, mu) #######
    C = C / n # 수정 C = t(X)%*%X/n - outer(mu, mu)  
    return( C)
  }
  stop('Invalid form')
}


compute_noise_vars = function(X, B){
  p = ncol(p)
  vars = rep(0,p)
  for( i in 1:p){
    m = X%*%B[i,]
    r = X[,i] - X%*%B[i,]
    vars[i] = sqrt(var(r)*(n-1)/n) # 파이썬과 계산 방식 차이로 인해 n-1/n 추가.
    }
  return( vars)
}

def read_data(fname):
  with open(fname) as fp:
    header = fp.readline()
    header = header.strip().split(",")
    try:
      v = float(header[0])
      fp.seek(0)
      data = np.loadtxt(fp, delimiter=',')
      header = NULL
    except ValueError:
      data = np.loadtxt(fp, delimiter=',')
  return header, data

def read_matrix(fname):
  with open(fname) as fp:
  try:
  header = fp.readline()
header = header.strip().split(",")
v = float(header[0])
fp.seek(0)
M = np.loadtxt(fp, delimiter=',')
except ValueError:
  hlength = len(header)
cols = range(1, hlength)
M = np.loadtxt(fp, delimiter=',', usecols=cols)

if M.shape[0] != M.shape[1]:
  raise Exception('Invalid data in file %s', fname)
return M


def write_matrix(fname, mat, var_names=None):
  if var_names is None:
  np.savetxt(fname, mat, fmt='%.6f', delimiter=',')
return
header = "," + ",".join(var_names)
with open(fname, 'w') as fp:
  fp.write("%s\n" % header)
for index in range(mat.shape[0]):
  line = ",".join(map(lambda x: '%.6f' % x, mat[index]))
fp.write("%s,%s\n" % (var_names[index], line))


def write_variances(fname, varss, var_names=None):
  with open(fname, 'w') as fp:
  if var_names is not None:
  header = ",".join(var_names)
fp.write("%s\n" % header)
line = ",".join(map(lambda x: '%.6f' % x, varss))
fp.write("%s\n" % line)


def read_variances(fname):
  with open(fname, 'rU') as fp:
  header = fp.readline()
header = fp.strip().split(",")
try:
  v = float(header[0])
return np.array(map(float, header))
except ValueError:
  return np.loadtxt(fp, delimiter=',')


def append_reg_param_to_fname(fname, reg_param):
  # Append lambda_<reg_param> in the file name and return new file name.
  splits = fname.split(".")
if len(splits) == 1:
  new_fname = "%s_lambda_%.6f.csv" % (fname, reg_param)
else:
  new_fname = "%s_lambda_%.6f.%s" % (".".join(splits[:-1]), reg_param, splits[-1])

return new_fname


def draw(M, fname, var_names=None, directed=False):
  if directed:
    G = nx.DiGraph()
  else:
    G = nx.Graph()
  
  if var_names is not None:
    nodes = var_names
  else:
    nodes = range(M.shape[0])
  
  G.add_nodes_from(nodes)
  # Get edges from adjacency matrix.
  temp = np.where(M)
  edges = temp[1], temp[0]
  edges = zip(*edges)
  if var_names is not None:
    edges = map(lambda (u,v): (var_names[u], var_names[v]), edges)
  G.add_edges_from(edges)
  
  nx.draw_spring(G, node_color='y',
                 with_labels=True)
  
  plt.savefig(fname, bbox_inches='tight')


##############
def main(): # 어떤 것인지 제대로 모르겠음.
  cx.solvers.options['show_progress'] = False

parser = argparse.ArgumentParser('listen', description='LISTEN: LInear STructural Equation model learNing.')
parser.add_argument('input', metavar='input', help='Data file.')
parser.add_argument('output', metavar='output', help='Output file. One output file is generated for each regularization parameter'
                    ' and the name of the file is appended by the regularization parameter.')
parser.add_argument('reg_params', metavar='lambda', help='Regularization parameter.', type=float, nargs='+')
parser.add_argument('--otype', default='sem', choices=['pre', 'sem', 'lik'], 
                    help='Output type. (i) pre: just compute precision matrix,'
                    ' (ii) sem: compute sem, i.e. the weight matrix and noise variances (optional, see --vars option),'
                    ' (iii) lik: just compute log-likelihood and bic score of the data (in input file) given either the precision matrix (--pre option)'
                    ' or the weight matrix (--weight option). If the precision matrix is specified then compute the undirected log-likelihood'
                    ' and bic score, while if the weight matrix is specified then compute the directed log-likelihood. Either one or both options'
                    ' must be specified. The reg_param argument is ignored if otype is lik.'
                    ' Defaults to sem.')
parser.add_argument('--pre', metavar='file', help='File containing precision matrix. If specified then use this precision matrix.', 
                    default=None)
parser.add_argument('--weight', metavar='file', help='File containing the weight matrix of the SEM.', 
                    default=None)
parser.add_argument('--vars', metavar='file', 
                    help='If otype is sem then the estimated noise variances are stored in this file after computing SEM.'
                    ' If otype is lik then use the noise variances in the file to compute the log-likelihhod. For otype lik'
                    ' if the --vars option is not provided then the likelihood is computed by estimating the noise variances from data file.'
                    ' Noise variances are only needed to compute the likelihood for the directed case.',
                    default=None)
parser.add_argument('--rho', metavar='N', type=float, help='Constant to add to diagonal of covariance matrix to make it invertible.'
                    'Default: 0', default=0)
parser.add_argument('--solver', choices=['glpk', 'mosek', None], default=None, 
                    help='LP solver to use. Defaults to (None) i.e. conelp solver in CVXOPT package.')
parser.add_argument('--draw', metavar='file', default=None, 
                    help='Draw and save the directed/undirected graph corresponding to the weight/precision matrix '
                    'in the specified file.')

args = parser.parse_args()

if args.solver == 'glpk':
  cx.solvers.options['glpk'] = {'msg_lev': 'GLP_MSG_OFF'}

inputf = args.input
reg_params = args.reg_params

C = None # Covariance matrix.

print "Reading data."
var_names, X = read_data(inputf)
print "Done."

if args.otype == 'lik':
  with open(args.output, 'w') as fp:
  if args.pre is None and args.weight is None:
  raise Exception('Must specify precision matrix or weight matrix for computing log likelihood.')
fp.write("Type,log_lik,bic\n")

if args.pre is not None:
  O = read_matrix(args.pre, var_names=var_names) # Inverse covariance matrix.
ll = log_lik_u(X, O)
bic = bic_score(ll, O, X.shape[0])
print "[Undirected] log_lik: %.6f, bic: %.6f" % (ll, bic)
fp.write("Undirected,%.6f,%.6f\n" % (ll, bic))

if args.weight is not None:
  B = read_matrix(args.weight, var_names=var_names) # Inverse covariance matrix.
noise_variances = None
if args.vars is not None:
  noise_variances = read_variances(args.vars)
ll = log_lik(X, B, noise_vars=noise_variances)
bic = bic_score(ll, B, X.shape[0])
print "[Directed] log_lik: %.6f, bic: %.6f" % (ll, bic)
fp.write("Directed,%.6f,%.6f\n" % (ll, bic))

return

for reg_param in reg_params:
  O = None # Inverse covariance matrix.
if args.pre is not None:
  O = read_matrix(args.pre)

if reg_param < 0:
  raise Exception('Regularization parameter should be greater than 0.')

if args.otype == 'pre':
  if O is None:
  if C is None:
  print "Computing covariance matrix."
t = time.time()
C = covariance(X)
t = time.time() - t
print "Done (%.3f sec)." % t

print "Computing precision matrix."
t = time.time()
O = inv_cov_estimate_clime(C, reg_param, solver=args.solver)
t = time.time() - t
print "Done (%.3f sec)." % t

print "Computing metrics."
ll = log_lik_u(X, O)
bic = bic_score(ll, O, X.shape[0])
print "reg_param: %.6f, log_lik: %.6f, bic: %.6f" % (reg_param, ll, bic)
write_matrix(append_reg_param_to_fname(args.output, reg_param), O, var_names=var_names)
if args.draw is not None:
  draw(B, args.draw, var_names=var_names, directed=True)

if args.otype == 'sem':
  if C is None:
  print "Computing covariance matrix."
t = time.time()
C = covariance(X)
t = time.time() - t
print "Done (%.3f sec)." % t

print "Learning SEM."
t = time.time()
if O is None:
  B, vars = learn_sem(C, reg_param, rho=args.rho, solver=args.solver)
else:
  B, vars = learn_weight_matrix(C, O, reg_param, solver=args.solver)
t = time.time() - t
# Re-estimate noise variances.
vars = compute_noise_vars(X, B)
print "Done (%.3f sec)." % t

print "Computing metrics."
ll = log_lik(X, B, noise_vars=vars)
bic = bic_score(ll, B, X.shape[0])
print "reg_param: %.6f, log_lik: %.6f, bic: %.6f" % (reg_param, ll, bic)
write_matrix(append_reg_param_to_fname(args.output, reg_param), B, var_names=var_names)
if args.vars:
  write_variances(append_reg_param_to_fname(args.vars, reg_param), vars, var_names=var_names)

if args.draw is not None:
  draw(B, args.draw, var_names=var_names, directed=True)

if __name__ == '__main__':
  main()







