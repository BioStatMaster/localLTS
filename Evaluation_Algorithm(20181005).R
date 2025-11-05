

############## Evaluation Methods: Ref, Eunho Yang  ###############
evaluation_fun = function(true_graph, estimated_graph){
  #Precision: the fraction of all predicted (directed) edges that are actually present in the true DAG
  #recall: the fraction of directed edges in the true DAG that the method was able to recover.
  #true positives: the numbers of correctly identified edges.
  #true negatives: the numbers of correctly identified absence of edges.
  #false positives: the numbers of edges were falsely added.
  #false negatives: the numbers of edges were falsely added falsely missing.
  true_graph_edge = true_graph+t(true_graph)
  estimated_graph_edge = estimated_graph+t(estimated_graph)
  true_graph_total_edges = sum(true_graph)
  estimated_graph_total_edges = sum(estimated_graph)
  
  precisition =  sum( (true_graph + estimated_graph) == 2 )/sum(estimated_graph)
  recall = sum( (true_graph + estimated_graph) == 2 )/sum(true_graph)
  precisition_edge =  sum( (true_graph_edge + estimated_graph_edge) == 2 )/sum(estimated_graph_edge)
  recall_edge = sum( (true_graph_edge + estimated_graph_edge) == 2 )/sum(true_graph_edge)
  
  true_positives = sum( (true_graph + estimated_graph) == 2 )
  true_negatives = sum( (true_graph + estimated_graph) == 0 )
  false_positives = sum( true_graph < estimated_graph)
  false_negatives = sum( true_graph > estimated_graph)
  
  true_positives_edge = sum( (true_graph_edge + estimated_graph_edge) == 2 )/2
  true_negatives_edge = sum( (true_graph_edge + estimated_graph_edge) == 0 )/2
  false_positives_edge = sum( true_graph_edge < estimated_graph_edge)/2
  false_negatives_edge = sum( true_graph_edge > estimated_graph_edge)/2
  
  hamming_dist = sum(true_graph!=estimated_graph)
  hamming_dist_edge = sum( true_graph_edge != estimated_graph_edge ) /2
  hamming_dist_ordering = hamming_dist - hamming_dist_edge
  
  return( c(precisition = precisition, 
            recall = recall, 
            precisition_edge = precisition_edge, 
            recall_edge = recall_edge,
            true_positives = true_positives, 
            true_negatives = true_negatives, 
            false_positives = false_positives, 
            false_negatives = false_negatives, 
            true_positives_edge = true_positives_edge, 
            true_negatives_edge = true_negatives_edge, 
            false_positives_edge = false_positives_edge, 
            false_negatives_edge = false_negatives_edge, 
            hamming_dist = hamming_dist, 
            hamming_dist_edge = hamming_dist_edge, 
            hamming_dist_ordering  = hamming_dist_ordering, 
            true_graph_total_edges = true_graph_total_edges, 
            estimated_graph_total_edges = estimated_graph_total_edges) )
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

