

############## Evaluation Methods: Ref, Eunho Yang  ###############
evaluation_fun = function(true_graph, estimated_graph){
  #Precision: the fraction of all predicted (directed) edges that are actually present in the true DAG
  #recall: the fraction of directed edges in the true DAG that the method was able to recover.
  #true positives: the numbers of correctly identified edges.
  #true negatives: the numbers of correctly identified absence of edges.
  #false positives: the numbers of edges were falsely added.
  #false negatives: the numbers of edges were falsely added falsely missing.
  diag(true_graph) = 0; diag(estimated_graph) = 0
  true_graph_edge = true_graph+t(true_graph); true_graph_edge[true_graph_edge!=0] = 1
  estimated_graph_edge = estimated_graph+t(estimated_graph); estimated_graph_edge[estimated_graph_edge!=0] = 1
  true_graph_total_edges = sum(true_graph)
  estimated_graph_total_edges = sum(estimated_graph)
  
  ### DAG ###
  true_vec = c(true_graph[lower.tri(true_graph)],true_graph[upper.tri(true_graph)])
  esti_vec = c(estimated_graph[lower.tri(estimated_graph)],estimated_graph[upper.tri(estimated_graph)])
  
  ### Skeleton ###
  true_edge_vec = true_graph_edge[lower.tri(true_graph_edge)]
  esti_edge_vec = estimated_graph_edge[lower.tri(estimated_graph_edge)]
  
  ### Evaluations ###
  TP = sum( (true_vec + esti_vec) == 2 ) # num of rocovered 1
  TN = sum( (true_vec + esti_vec) == 0 ) # num of rocovered 0
  FP = sum( true_vec < esti_vec) # num of fasely recovered 1
  FN = sum( true_vec > esti_vec) # num of fasely recovered 0
  
  precisition =  TP /(TP+FP) #positive predictive value, 
  recall = TP/(TP+FN) #sensitivity
  specificity = TN/(TN+FP) #specificity
  
  ### Evaluations for Edges ###
  TP_E = sum( (true_edge_vec + esti_edge_vec) == 2 ) # num of rocovered 1
  TN_E = sum( (true_edge_vec + esti_edge_vec) == 0 ) # num of rocovered 0
  FP_E = sum( true_edge_vec < esti_edge_vec) # num of fasely recovered 1
  FN_E = sum( true_edge_vec > esti_edge_vec) # num of fasely recovered 0
  
  precisition_E =  TP_E /(TP_E+FP_E) #positive predictive value, 
  recall_E = TP_E/(TP_E+FN_E) #sensitivity
  specificity_E = TN_E/(TN_E+FP_E) #specificity
  
  hamming_dist = FP + FN
  hamming_dist_E = FP_E + FN_E
  hamming_dist_ordering = hamming_dist - hamming_dist_E
  
  return( list(hamming_dist = hamming_dist,
            hamming_dist_E = hamming_dist_E, 
            precisition = precisition,
            recall = recall,
            specificity = specificity,
            precisition_E = precisition_E, 
            recall_E = recall_E,
            specificity_E = specificity_E,
            TP = TP, 
            TN = TN, 
            FP = FP, 
            FN = FN, 
            TP_E = TP_E, 
            TN_E = TN_E, 
            FP_E = FP_E, 
            FN_E = FN_E,
            True_NUM_Edges = true_graph_total_edges, 
            Esti_NUM_Edges = estimated_graph_total_edges) )
}


############## Estimated directed graph: estimated_graph_fun ####################
estimated_graph_fun = function(directed_graph_edges, ordering){
  for(i in 1:length(ordering)) directed_graph_edges[ordering[1:i],ordering[i]]<-0
  return( directed_graph_edges )
}


############## DAG2CPDAGAdj ##############
dag2cpdagAdj <- function(Adj){
  library(graph)
  #library(pcalg)
  library(bnlearn)
  
  d <- as(Adj, "graphNEL")
  cpd <- cpdag(as.bn(d) )
  result<- amat(cpd)
  
  # if pcalg is allowed, the following code works.
  #cpd <- dag2cpdag(d)
  #result <- as(cpd, "matrix")
  return(result)
}

