
setwd("C:/Users/gwpar/Dropbox/GroupMeeting/Algorithm(Rcode, 20181220)")
source("GaussianSEM/LTS_GSEM Algorithm(20211013).R")
source("GaussianSEM/PGSEM_Algorithm.R")
source("GaussianSEM/HLSEM Algorithm(20200805).R")
source("GaussianSEM/GSEM_generator(20190930).R")
source("GaussianSEM/HGSEM Algorithm(20190220).R")
source("GaussianSEM/GSEM_Learning_Algorithm(20190413).R")
source("ComparisonAlg/GLS(ghoash)/learn_gbn.R")
source("ComparisonAlg/LISTEN(ghoash)/LISTEN_Algorithm(20190515).R")
source("ComparisonAlg/GES_Algorithm(20181219).R")
source("ComparisonAlg/GDS(peter)/GDS(20190515).R")
source("ComparisonAlg/OnCausal/EqVarDAG-master/R/EqVarDAG_HD_TD.R")
source("Evaluation_Algorithm(20181005).R")

############
#### Packages ####
if(!require(Rcpp)){
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(pcalg)){
  install.packages("pcalg")
  library(pcalg)
}
if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}
if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}
if(!require(graph)){
  install.packages("graph")
  library(graph)
}
if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}
if(!require(doParallel)){
  install.packages("doParallel")
  library(doParallel)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(plyr)){
  install.packages("plyr")
  library(plyr)
}
if(!require(hypergeo)){
  install.packages("hypergeo")
  library(hypergeo)
}
if(!require(glmnet)){
  install.packages("glmnet")
  library(glmnet)
}
if(!require(gamlss)){
  install.packages("gamlss")
  library(gamlss)
}
if(!require(bnlearn)){
  install.packages("bnlearn")
  library(bnlearn)
}
library(compiler)
if(!require(rmutil)){
  install.packages("rmutil")
  library(rmutil)
}
if(!require(CombMSC)){
  install.packages("CombMSC")
  library(CombMSC)
}

source("https://bioconductor.org/biocLite.R")
if(!require(Rgraphviz)){
  biocLite("Rgraphviz")
  n
  library(Rgraphviz)
}
if(!require(RBGL)){
  biocLite("RBGL")
  n
  library(RBGL)
}
if(!require(graph)){
  biocLite("graph")
  n
  library(graph)
}

#########################
#### Data Generation ####
########################

#### Simulation Settings for Toy Example #####
p = 20; n = 50; d = 1; seed = sample(1:100, 1);
beta_min = 1.00; beta_max = 1.00;
graph_type = 5
alpha = 1* max( c( 1 - pnorm( n^(1/3)/2 ), 0.1^100 * 1 ) )

outlier_node = 1
b = 1
synthetic.graph = LTS_GSEM_generator( n, p, d = 1, b = b, outlier_node = outlier_node, var_Min = 0.75, var_Max = 0.75, dist = "Gaussian", beta_min, beta_max, graph_type = 5, q = 0.025, h = 1, structure = NULL, seed = 1, path = NULL)
graph = synthetic.graph$true_Matrix
data = synthetic.graph$x
par(mfrow = c(3,1))

colSums(graph!=0)

boxplot(data)
barplot( sapply(data, mean))
barplot( sapply(data, var))

####### Simulation Start! ##############
res1 = res2 = res3 = res4 = res5 = NULL
res1 = LTS_GSEM_Algorithm(data, lambda1 = 0.1^4, lambda2 = 0.25*sqrt(log(p)/(n*0.5)), alpha = 0.5, thresh = 10, graph = graph)
res1


res2 = GSEM_Algorithm(data, alpha = alpha, direction ="backward", graph = graph, max_degree = 1)
res3 = HLSEM_Algorithm(data, sparsity_level = 2, CV = F, graph = graph)
res4 = HGSEM_Algorithm(data, alpha = alpha, sparsity_level = 2 , graph = graph)
res5 = TD(data, q = d, graph = graph)

res = cbind( RGSM = res1[[1]], US = res2[[1]], HLSM = res5[[1]], HGSM = res4[[1]], TD = res5[[1]] )
cbind( round( res, 3)[c(1,2,13:17),])





############# LTS_GSEM  #################
numCores <- 34
myCluster <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(myCluster)

P_real = 50
N_real = 200
Num_out = 1
Thresh = 2
H_ratio = 0.5

P_real = seq(25, 30, by = 5)
Num_out = c(1, seq(30, 30, by = 30))
P_real = seq(5, 25, by = 200)
Num_out = seq(30,30, by = 30)


N_real = seq(50, 500, by = 50)
P_real = c(10,20)
Num_out = c(1)
Model = "Outlier1";
Thresh = c(1, 1.5, 3)
H_ratio = 0.5

P_real = seq(20, 25, by = 5)

P_real = c(50)
Model = "Outlier1";
Num_out = c(1,30)
Thresh = 2
H_ratio = seq(0.5, 0.9, by =1.1)
N_real = seq(50, 500, by = 50)
d = 1

for(p_real in P_real){
  if(Model != 'Outlier1'){
    Model = paste0("OutlierAll",p_real)
  }
  for( model in Model ){
    for( num_out in Num_out){
      path  = paste0("D:/2022LTS_GSEM/Data/",model,"_",num_out)
      for(n_real in N_real){
        for(h_ratio in H_ratio){
          for(thresh in Thresh){
            print(n_real)
            evaluation_result = NULL
            evaluation_result = foreach::foreach(i = 1:100
                                                 , .combine = list, .multicombine = TRUE)  %dopar% {
                                                   tryCatch({
                                                     LTS_GSEM_simulation_fun(seed = i, n_real = n_real, p_real = p_real, d = d, h_ratio = h_ratio, thresh = thresh, path = path )
                                                   }, error = function(e) {
                                                     return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                                   })
                                                 }
            GaussianDAGresult_filename = paste0("LTS_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_h",h_ratio * 100,"_thresh",thresh,"_result.Rdata")
            if(model =="Outlier1"){
              setwd( paste0("D:/2022LTS_GSEM/Result/",model) )
            }else{
              setwd( paste0("D:/2022LTS_GSEM/Result/",'OutlierAll') )
            }
            save(evaluation_result, file = GaussianDAGresult_filename)
          }
        }
      }
    }
  }
}

##### Unifying US algorithm ####
P_real = c(50)
Num_out = c(1,30)
#Num_out = c(1, seq(30, 120, by = 30))
#Model = "OutlierAll";
N_real = seq(50, 500, by = 50)
d = 1

for(p_real in P_real){
  if(Model != 'Outlier1'){
    Model = paste0("OutlierAll",p_real)
  }
  for( model in Model ){
    for( num_out in Num_out){
      path  = paste0("D:/2022LTS_GSEM/Data/",model,"_",num_out)
      for(n_real in N_real){
        print(n_real)
        evaluation_result = NULL
        evaluation_result = foreach::foreach(i = 1:100
                                             , .combine = list, .multicombine = TRUE)  %dopar% {
                                               tryCatch({
                                                 GSEM_simulation_fun(seed = i, n_real, p_real, d = d, max_degree = 1, direction = 'backward', alpha = 1 - pnorm( n_real^(1/3)/2 ) , path = path)
                                               }, error = function(e) {
                                                 return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                               })
                                             }
        GaussianDAGresult_filename = paste0("USB_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(model =="Outlier1"){
          setwd( paste0("D:/2022LTS_GSEM/Result/",model) )
        }else{
          setwd( paste0("D:/2022LTS_GSEM/Result/","OutlierAll") )
        }
        save(evaluation_result, file = GaussianDAGresult_filename)
      }
    }
  }
}

##### HLSEM algorithm ####
for(p_real in P_real){
  if(Model != 'Outlier1'){
    Model = paste0("OutlierAll",p_real)
  }
  for( model in Model ){
    for( num_out in Num_out){
      path  = paste0("D:/2022LTS_GSEM/Data/",model,"_",num_out)
      for(n_index in N_real){
        print(n_index)
        n_real = n_index;  evaluation_result = NULL
        evaluation_result = foreach::foreach(i = 1:100
                                             , .combine = list, .multicombine = TRUE)  %dopar% {
                                               tryCatch({
                                                 HLSEM_simulation_fun(seed = i, n_real, p_real, d = d, path = path)
                                               }, error = function(e) {
                                                 return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                               })
                                             }
        GaussianDAGresult_filename = paste0("HLSM_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(model =="Outlier1"){
          setwd( paste0("D:/2022LTS_GSEM/Result/",model) )
        }else{
          setwd( paste0("D:/2022LTS_GSEM/Result/","OutlierAll") )
        }
        save(evaluation_result, file = GaussianDAGresult_filename)
      }
    }
  }
}

########## HGSEM #################################
for(p_real in P_real){
  if(Model != 'Outlier1'){
    Model = paste0("OutlierAll",p_real)
  }
  for( model in Model ){
    for( num_out in Num_out){
      path  = paste0("D:/2022LTS_GSEM/Data/",model,"_",num_out)
      for(n_index in N_real){
        print(n_index)
        n_real = n_index;  evaluation_result = NULL
        evaluation_result = foreach::foreach(i = 1:100
                                             , .combine = list, .multicombine = TRUE)  %dopar% {
                                               tryCatch({
                                                 HGSEM_simulation_fun(seed = i, n_real, p_real, d = d, alpha = 1 - pnorm( n_real^(1/3)/2 ), path = path)
                                               }, error = function(e) {
                                                 return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                               })
                                             }
        GaussianDAGresult_filename = paste0("HGSM_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(model =="Outlier1"){
          setwd( paste0("D:/2022LTS_GSEM/Result/",model) )
        }else{
          setwd( paste0("D:/2022LTS_GSEM/Result/","OutlierAll") )
        }
        save(evaluation_result, file = GaussianDAGresult_filename)
      }
    }
  }
}

##### TD ##########
for(p_real in P_real){
  if(Model != 'Outlier1'){
    Model = paste0("OutlierAll",p_real)
  }
  for( model in Model ){
    for( num_out in Num_out){
      path  = paste0("D:/2022LTS_GSEM/Data/",model,"_",num_out)
      for(n_index in N_real){
        print(n_index)
        n_real = n_index;
        evaluation_result = foreach::foreach(i = 1:100
                                             , .combine = list, .multicombine = TRUE)  %dopar% {
                                               tryCatch({
                                                 TD_simulation_fun(seed = i, n_real = n_real, p_real = p_real, d = 1, path = path )
                                               }, error = function(e) {
                                                 return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                               })
                                             }
        
        GaussianDAGresult_filename = paste0("TD_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(model =="Outlier1"){
          setwd( paste0("D:/2022LTS_GSEM/Result/",model) )
        }else{
          setwd( paste0("D:/2022LTS_GSEM/Result/","OutlierAll") )
        }
        save(evaluation_result, file = GaussianDAGresult_filename)
      }
    }
  }
}

############# GDS  #################
for(p_real in P_real){
  if(Model != 'Outlier1'){
    Model = paste0("OutlierAll",p_real)
  }
  for( model in Model ){
    for( num_out in Num_out){
      path  = paste0("D:/2022LTS_GSEM/Data/",model,"_",num_out)
      for(n_real in N_real){
        print(n_real)
        evaluation_result = foreach::foreach(i = 1:100
                                             , .combine = list, .multicombine = TRUE)  %dopar% {
                                               tryCatch({
                                                 GDS_simulation_fun(seed = i, n_real = n_real, p_real = p_real, d = d, path = path)
                                               }, error = function(e) {
                                                 return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                               })
                                             }
        GaussianDAGresult_filename = paste0("GDS_p",p_real,"_n",n_real,"_d",d,"_b",num_out,"_result.Rdata")
        if(model =="Outlier1"){
          setwd( paste0("D:/2022LTS_GSEM/Result/",model) )
        }else{
          setwd( paste0("D:/2022LTS_GSEM/Result/","OutlierAll") )
        }
        save(evaluation_result, file = GaussianDAGresult_filename)
      }
    }
  }
}


parallel::stopCluster(myCluster)


################################

#### Data Generation ####################
## 1: tree, 2: bipartite, 3: cycle, 4: random number of parents , 5: fixed number of parents random2, 6 fixed ordering random, 7 manual ##
numCores <- 34
myCluster <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(myCluster)

n = 1000; p = 25; beta_min = 0.75; beta_max = 1.0
model = "Random"; graph_type = 5;  
D_real = c(1)
Num_out = c(1)
outlier_node = 1:p

for(p in c(50)){
  outlier_node = 1:p
  for(j in 1:length(Num_out)){
    path  = paste0("D:/2022LTS_GSEM/Data/OutlierAll",p,"_",Num_out[j])
    data = foreach::foreach(i = 1:100
                            , .combine = list, .multicombine = TRUE) %dopar% {
                              tryCatch({
                                synthetic.graph = LTS_GSEM_generator( n, p, d = 1, b = Num_out[j], outlier_node = outlier_node, var_Min = 1, var_Max = 1, dist = "Gaussian", beta_min, beta_max, graph_type = 5, q = 0.025, h = 1, structure = NULL, seed = i, path = path)
                              }, error=function(e){})
                            }  
  }
}
parallel::stopCluster(myCluster)






### Random Graphs ###
n = 2000; p = 50; beta_min = 0.75; beta_max = 1.0
model = "Random"; graph_type = 5; D_real = c(1,2,3)
model =  "Tree"; graph_type = 1; D_real = 1

path  = paste0("D:/2021PGSEM/Data/",model)
for(d in D_real){
  data = foreach::foreach(i = 1:100
                          , .combine = list, .multicombine = TRUE) %dopar% {
                            tryCatch({
                              synthetic.graph = GSEM_generator_PGSEM( n, p, d, beta_min, beta_max, graph_type, q = 0.025, structure = NULL, seed = i, path = path)
                            }, error=function(e){})
                          }
}



########################



### TD ##########
for(p_real in P_real){
  for( model in Model ){
    path  = paste0("D:/2020HighLSEM/Data/",model,"P",p_real)
    for( d in c(5,8) ){
      for(n_index in N_real){
        print(n_index)
        n_real = n_index;
        evaluation_result = foreach::foreach(i = 1:100
                                             , .combine = list, .multicombine = TRUE)  %dopar% {
                                               tryCatch({
                                                 TD_simulation_fun(seed = i, n_real = n_real, p_real = p_real, d = d, path = path )
                                               }, error = function(e) {
                                                 return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                               })
                                             }

        GaussianDAGresult_filename = paste0("TD_GaussianDAG_p",p_real,"_n",n_real,"_d",d,"_result.Rdata")

        if(model =="GSEM"){
          setwd("D:/2020HighLSEM/Result/GSEM")
        }else if(model =="LSEM"){
          setwd("D:/2020HighLSEM/Result/LSEM")
        }
        
        save(evaluation_result, file = GaussianDAGresult_filename)
      }
    }
  }
}



#### PC, GES, Lingam ##########
for(p_real in P_real){
  for( model in Model ){
    path  = paste0("D:/2021PGSEM/Data/",model,"P",p_real)
    for( method in c("Lingam") ){
      for( d in D_real ){
        for(n_index in N_real){
          n_real = n_index;
          evaluation_result = foreach::foreach(i = 1:100
                                               , .combine = list, .multicombine = TRUE) %dopar% {
                                                 tryCatch({
                                                   GES_simulation_fun(seed = i, n_real = n_real, p_real = p_real, d = d, method = method, alpha = 1 - pnorm( n_real^(1/3)/2 ), model = model )
                                                 }, error = function(e) {
                                                   return(paste0("The variable '", i, "'", " caused the error: '", e, "'"))
                                                 })
                                               }
          
          GaussianDAGresult_filename = paste0(method,"_GaussianDAG_p",p_real,"_n",n_real,"_d",d,"_result.Rdata")
          if(model =="Random"){
            setwd("D:/2021PGSEM/Result/Random")
          }else if(model =="Tree"){
            setwd("D:/2021PGSEM/Result/Tree")
          }
          save(evaluation_result, file = GaussianDAGresult_filename)
        }
      }
    }
  }
}
#


## 1: tree, 2: bipartite, 3: cycle, 4: random , 5: fixed number of parents random2, 6 fixed ordering random, 7 manual ##

####
numCores <- 34
myCluster <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(myCluster)



##### US algorithm #########


