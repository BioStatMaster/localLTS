function PVAL = indtest_corr(X,Y,Z,pars)

if isempty(Z) %unconditional HSIC
    [RHO,PVAL] = corr(X,Y);
else
    [RHO,PVAL] = partialcorr(X,Y,Z);
end