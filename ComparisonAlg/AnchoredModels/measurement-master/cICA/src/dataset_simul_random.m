function [X, B, G] = dataset_simul_random(pars)
% simulation dataset, random graphs
id = pars.id;
T = pars.T;
filePath = pars.filePath;
load(fullfile(filePath,sprintf('graph_structures%d_%d.mat', id, T)));
X = Data_save;
B = B_save;
G = G_save;
