% Learning the causal structure with GES search using BIC score
% example 1
clear all,clc,close all
addpath(genpath(pwd))

D = 3;
T = 500;
X1 = randn(T,1);
X2 = 0.8*X1 + 0.5*randn(T,1);
X3 = 0.8*X2 + 0.5*randn(T,1);
X = [X1,X2,X3];

X=X-repmat(mean(X),size(X,1),1);
X=X*diag(1./std(X));

maxP = 2; % maximum number of parents when searching the graph
Record = GES(X,maxP);

%function [Record] = GES(X,maxP)
%INPUT:
% X: Data with T*D dimensions
% maxP: allowed maximum number of parents when searching the graph

% Output:
% Record.G: learned causal graph
% Record.update1: each update (Insert operator) in the forward step
% Record.update2: each update (Delete operator) in the backward step
% Record.G_step1: learned graph at each step in the forward step
% Record.G_step2: learned graph at each step in the backward step
% Record.score: the score of the learned graph
