

##### Simulations for GES and MMHC ####
GES_simulation_fun = function(seed, n_real, p_real, d = 1, method = "GES", test = "cor", alpha = 0.05, model = "homoGSEM" , path){

  # setwd("C:/Users/gwpar/Dropbox/GroupMeeting/Algorithm(Rcode, 20181220)")
  # source("ComparisonAlg/GES_Algorithm(20181219).R")
  # 
  ###Poisson SEM #############
  # high dimensional PSEM
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonDAG")
  #p = 1000; n = 2000; d = 5; # result is for PoissonDAGresult
  #p = 100 ; n = 1000; d = 2;
  
  # high dimensional PSEM D1- D5
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonD1D5Exp")
  #p = 100 ; n = 1000; d = 5;
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonD1Exp")
  #p = 1000 ; n = 2000; d = 1;
  #p = 200 ; n = 2000; d = 1;
  
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonD1D10Exp")
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonD10Exp")
  #p = 1000 ; n = 2000; d = 10;
  
  # high dimensional PDAG: identity link
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonDAG2")
  #p = 500; n = 5000; d = 2;
  #setwd("C:/2018AISTATS_GSEM/R-code/PoissonD1D2Identity")
  #p = 100; d = 2;
  
  # Hyper Poisson DAG: all Poisson with log link
  #setwd("C:/PoissonDAG/PoissonD2Exp")
  #p = 1000; n = 2000; d = 2; 
  #PoissonDAG_filename = paste0("PoissonDAG_p",p,"_d",d,"_seed",seed,".Rdata")
  #load(PoissonDAG_filename)
  # 
  # if(model =="homoGSEM"){
  #   setwd("C:/2019GaussianSEM_SimulationResult/Data/Homo"); p = 200;
  # }
  # if(model =="heteroGSEM"){
  #   #setwd("C:/2019GaussianSEM_SimulationResult/Data/Hetero"); p = 200;
  #   setwd("C:/2020HGSEM_Simulation/Data/Hetero"); # p =3000;
  # }
  # if(model =="polyhomoGSEM"){
  #   setwd("C:/2019GaussianSEM_SimulationResult/Data/PolyHomo")
  #   p = 20;
  # }
  # if(model =="polyheteroGSEM"){
  #   setwd("C:/2019GaussianSEM_SimulationResult/Data/PolyHetero")
  #   p = 20;
  # }
  # if(model =="NonGSEM"){
  #   setwd("C:/2019GaussianSEM_SimulationResult/Data/NonGauss")
  #   p = 20;
  # }
  # if(model =="NonGSEM2"){
  #   setwd("C:/2019GaussianSEM_SimulationResult/Data/NonGauss2")
  #   p = 20;
  # }
  # if(model =="NonGPolySEM"){
  #   setwd("C:/2019GaussianSEM_SimulationResult/Data/NonGPolySEM")
  #   p = 20;
  # }
  # if(model =="Drton"){
  #   setwd("C:/2020HGSEM_Simulation/Data/Drton"); # p =3000;
  # }
  # if(model =="Park"){
  #   setwd("C:/2020HGSEM_Simulation/Data/Park");
  # }
  #setwd("C:/2019GaussianSEM_SimulationResult/Data/NonFaithful");p = 3
  #GaussianDAG_filename = paste0("GaussianSEM_p",p,"_d",d,"_seed",seed,".Rdata")
  
  ###Gaussian SEM Load #############
  setwd(path)
  GaussianDAG_filename = paste0("GaussianSEM_d",d,"_seed",seed,".Rdata")
  load(GaussianDAG_filename)
  ###################################
  
  synthetic.graph =DAGsamples
  graph = synthetic.graph$true_Matrix[1:p_real, 1:p_real]
  data = synthetic.graph$x[1:n_real,1:p_real]
  
  if( 1 == 1  ){
    GES_output = MEC_Alg(data, method, test = test, alpha = alpha)
    result_DAG = GES_output$DAG
    result_MEC = GES_output$MEC
    time_GES = GES_output$Runtime
  }else{
    result_DAG = NULL
    result_MEC = NULL
    time_GES = NULL
  }
  
  

  #### Evaluation ####
  B = graph
  B[B!=0] =1
  
  #A <- (B != result_DAG)
  #rbind( rowSums(B), colSums(B),  rowSums(A), colSums(A),  rowSums(result_DAG), colSums(result_DAG) )
  
  result_GES_DAG = evaluation_fun( B, result_DAG )
  result_GES_MEC = evaluation_fun( dag2cpdagAdj(B), result_MEC )

  # Return: 1.DAG, 2.MEC, 3.Oracle(known parents), 4.Estimated DAG, 5.Estimated Ordering, 6.Time
  
  return(list( DAG_Evaluation = result_GES_DAG, 
        MEC_Evaluation = result_GES_MEC,  
        Oracle_Evaluation = NULL, 
        MEC = result_MEC, 
        Ordering = NULL , 
        Time = time_GES)
  )
  
}

##### Comparison Algorithm: PC, GS, GES and MMHC ####
MEC_Alg = function(data, method = "GES", test = "cor", alpha = 0.05){
  #this algorithm has GES, MMHC, GS, PC methods. 
  
  library(bnlearn)
  library(pcalg)
  data = data.frame(data)
  colnames(data) = 1:ncol(data)
  p = ncol(data)
  est_MEC = est_DAG = matrix(0,p,p)
  MEC = NULL
  
  Runtime = proc.time()[3]
  if(method =="GES")  MEC = rsmax2(data)
  if(method =="MMHC") MEC = mmhc(data, restrict.args = list(test = test, alpha = alpha) )
  if(method =="GS") MEC = hc(data)
  if(method =="PC") MEC = pc.stable(data, alpha = alpha)
  if( !is.null(MEC) ){
    est_MEC = t(amat(cpdag(MEC)))
    vstructs(MEC)
    
    est_DAG = est_MEC
    est_DAG[est_MEC == t(est_MEC)] = 0
  }
  if(method =="LINGAM"){
    est_DAG = lingam(data)$Bpruned
    est_DAG[est_DAG !=0] = 1
    est_MEC = t( dag2cpdag( t(est_DAG) ) )
  } 
  Runtime = proc.time()[3] - Runtime
  
  return( list(DAG = est_DAG, MEC = est_MEC, Runtime = Runtime) )
}



