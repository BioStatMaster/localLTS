samplesize=1000;
alpha=0.05;
X=rand(samplesize,1)-0.5;
Y=rand(samplesize,1)-0.5;
Z=0.5*X+rand(samplesize,1)-0.5;
W=X+1.5*Z-1.9*Y+rand(samplesize,1)-0.5;
dag=pc([X Y Z W],'indtest_corr',[],2,alpha);
dag