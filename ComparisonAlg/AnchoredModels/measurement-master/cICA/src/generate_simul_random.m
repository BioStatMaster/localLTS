% generate random graph and data

% generate structure
parsG.N = 5; 
parsG.nPA = 2; 
parsG.min_PA = 1; 
parsG.degree = 3; 
parsG.graphFile = '../../datasets/simul_random/graph_structures0.mat';
generate_strcucture(parsG);

% generate data
parsD.N = parsG.N;
parsD.T = 500;
parsD.graphFile = parsG.graphFile;
parsD.dataFile = sprintf('../../datasets/simul_random/graph_structures0_%d.mat', parsD.T);
generate_data(parsD);
