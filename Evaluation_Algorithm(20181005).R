

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

############### Ploting functions for DAG ##############
# B[j, k] != 0 이면 edge k -> j 로 해석 (본문 코드와 동일한 규약)
as_igraph_from_B <- function(B, vertex_labels = NULL) {
  if (!is.matrix(B)) stop("B must be a matrix.")
  if (nrow(B) != ncol(B)) stop("B must be square.")
  p <- nrow(B)

  # edge list: from k to j whenever B[j,k] != 0
  ed <- which(B != 0, arr.ind = TRUE)
  if (nrow(ed) == 0) {
    g <- make_empty_graph(n = p, directed = TRUE)
  } else {
    # ed has rows (j, k); convert to vector (k, j, k, j, ...)
    edges_vec <- as.vector(t(cbind(ed[,2], ed[,1])))
    g <- graph(edges = edges_vec, n = p, directed = TRUE)
    E(g)$weight <- B[cbind(ed[,1], ed[,2])]
  }

  # labels
  if (is.null(vertex_labels)) vertex_labels <- as.character(seq_len(p))
  V(g)$name   <- vertex_labels
  V(g)$label  <- vertex_labels

  # 간단한 스타일 속성
  V(g)$size   <- 18
  V(g)$color  <- "grey95"
  V(g)$frame.color <- "grey60"

  if (!is.null(E(g)$weight)) {
    w <- E(g)$weight
    E(g)$color <- ifelse(w >= 0, "#1f78b4", "#e31a1c")  # blue(+) / red(-)
    # 가중치 절댓값으로 두께 스케일
    sc <- 1 + 3 * (abs(w) / (max(abs(w)) + 1e-12))
    E(g)$width <- sc
  } else {
    E(g)$color <- "#1f78b4"
    E(g)$width <- 1.5
  }

  E(g)$arrow.size <- 0.6
  E(g)$curved <- 0

  g
}

# DAG이면 위상정렬 기반 레이아웃(sugiyama), 아니면 KK 레이아웃으로 안전하게
plot_sem_graph <- function(B,
                           vertex_labels = NULL,
                           layout_method = c("auto","sugiyama","kk","tree"),
                           main = "SEM Graph (k → j if B[j,k] ≠ 0)") {
  layout_method <- match.arg(layout_method)
  g <- as_igraph_from_B(B, vertex_labels)

  # 사이클 유무 확인
  is_dag <- is_dag(g)

  # 레이아웃 선택
  lay <-
    if (layout_method == "sugiyama" || (layout_method == "auto" && is_dag)) {
      # igraph의 sugiyama는 레이어/방향을 잘 보이게 해줌
      layout_with_sugiyama(g)$layout
    } else if (layout_method == "tree" && is_dag) {
      layout_as_tree(g, mode = "out")
    } else if (layout_method == "kk") {
      layout_with_kk(g)
    } else if (is_dag) {
      layout_with_sugiyama(g)$layout
    } else {
      layout_with_kk(g)
    }

  plot(g,
       layout = lay,
       main   = main,
       vertex.label.family = "sans",
       vertex.label.cex    = 0.9,
       edge.arrow.size     = E(g)$arrow.size,
       edge.curved         = E(g)$curved)
  invisible(g)
}
