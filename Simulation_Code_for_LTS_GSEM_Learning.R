
# setwd("C:/Users/gwpar/Dropbox/GroupMeeting/Algorithm(Rcode, 20181220)")
source("GaussianSEM/LTS_GSEM Algorithm(20211013).R")
source("GaussianSEM/PGSEM_Algorithm.R")
source("GaussianSEM/HLSEM Algorithm(20200805).R")
source("GaussianSEM/GSEM_generator(20190930).R")
source("NEW.R")
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
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(plyr)){
  install.packages("plyr")
  library(plyr)
}
if(!require(future)){
  install.packages("future")
  library(future)
}
if(!require(future.apply)){
  install.packages("future.apply")
  library(future.apply)
}
if(!require(progressr)){
  install.packages("progressr")
  library(progressr)
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
synthetic.graph = CCLSM_generator(
  n = n,
  p = p,
  d = d,
  b = b,
  outlier_nodes = outlier_node,
  var_Min = 0.75,
  var_Max = 0.75,
  dist = "Gaussian",
  beta_min = beta_min,
  beta_max = beta_max,
  graph_type = graph_type,
  q = 0.025,
  h = 1,
  seed = seed
)
graph = synthetic.graph$true_Matrix
data = synthetic.graph$x
Y_mask = synthetic.graph$Y
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

res = cbind(
  LTS = res1[[1]],
  US = res2[[1]],
  HLSEM = res3[[1]],
  HGSEM = res4[[1]],
  TD = res5[[1]]
)
cbind( round( res, 3)[c(1,2,13:17),])

####### 단일 함수 기반 파이프라인 테스트 ########
toy_simulation <- CCLSM_simulation_fun(
  seed = seed,
  n_real = n,
  p_real = p,
  d = d,
  beta_min = beta_min,
  beta_max = beta_max,
  graph_type = graph_type,
  b = b,
  outlier_nodes = outlier_node,
  h_ratio = 0.5,
  thresh = 10,
  generator_output = synthetic.graph,
  save_data = TRUE
)
head(toy_simulation$evaluations)






############# Future 기반 병렬 시뮬레이션 #################

#' 시뮬레이션에서 사용할 평가 지표 이름 목록
evaluation_metric_names <- c(
  "precisition", "recall", "precisition_edge", "recall_edge",
  "true_positives", "true_negatives", "false_positives", "false_negatives",
  "true_positives_edge", "true_negatives_edge", "false_positives_edge", "false_negatives_edge",
  "hamming_dist", "hamming_dist_edge", "hamming_dist_ordering",
  "true_graph_total_edges", "estimated_graph_total_edges"
)

#' 초 단위 시간을 `HH:MM:SS.ss` 문자열로 변환한다.
format_duration <- function(seconds) {
  if (is.na(seconds) || is.infinite(seconds)) {
    return("NA")
  }
  seconds <- max(seconds, 0)
  hrs <- floor(seconds / 3600)
  mins <- floor((seconds %% 3600) / 60)
  secs <- seconds %% 60
  sprintf("%02d:%02d:%05.2f", hrs, mins, secs)
}

#' 알고리즘 수행 결과에서 평가 지표 한 행을 생성한다.
build_metric_row <- function(name, result) {
  metrics <- setNames(rep(NA_real_, length(evaluation_metric_names)), evaluation_metric_names)
  if (!is.null(result$DAG_Evaluation)) {
    metrics[names(result$DAG_Evaluation)] <- result$DAG_Evaluation
  }
  row <- data.frame(
    algorithm = name,
    Time = if (!is.null(result$Time)) result$Time else NA_real_,
    error = if (!is.null(result$error)) result$error else NA_character_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (metric_name in names(metrics)) {
    row[[metric_name]] <- metrics[[metric_name]]
  }
  row
}

#' 오류가 발생하더라도 기본 구조를 유지하도록 알고리즘 호출을 감싼다.
safe_algorithm_call <- function(expr) {
  tryCatch(expr, error = function(e) {
    list(
      DAG_Evaluation = NULL,
      MEC_Evaluation = NULL,
      Oracle_Evaluation = NULL,
      DAG = NULL,
      Ordering = NULL,
      Time = NA_real_,
      error = conditionMessage(e)
    )
  })
}

#' CCLSM 기반 데이터 생성과 알고리즘 평가를 한 번에 수행한다.
#'
#' @param generator_output `CCLSM_generator()`가 반환한 결과. 제공되면
#'   새로운 데이터를 생성하지 않고 이 객체를 그대로 사용한다.
CCLSM_simulation_fun <- function(
    seed,
    n_real,
    p_real,
    d = 1,
    beta_min = 0.8,
    beta_max = 1.2,
    graph_type = 5,
    var_Min = 1,
    var_Max = 1,
    outlier_nodes = NULL,
    b = 0,
    overlap_mode = "independent",
    overlap_fraction = 0,
    union_rows = NULL,
    Y_custom = NULL,
    ensure_disjoint_rest = TRUE,
    dist = "Gaussian",
    out_mean = 1e3,
    out_scale = 1,
    h_ratio = 0.5,
    thresh = qnorm(0.975),
    hlsem_sparsity = max(1, d + 1),
    hgsem_sparsity = max(1, d + 1),
    td_sparsity = 3,
    gsem_direction = "backward",
    gsem_alpha = NULL,
    gsem_max_degree = d,
    gds_startAt = "emptyGraph",
    save_data = FALSE,
    generator_output = NULL,
    ...
) {
  set.seed(seed)
  gen <- generator_output
  if (is.null(gen)) {
    gen <- CCLSM_generator(
      n = n_real,
      p = p_real,
      d = d,
      b = b,
      outlier_nodes = outlier_nodes,
      var_Min = var_Min,
      var_Max = var_Max,
      dist = dist,
      beta_min = beta_min,
      beta_max = beta_max,
      graph_type = graph_type,
      seed = seed,
      overlap_mode = overlap_mode,
      overlap_fraction = overlap_fraction,
      union_rows = union_rows,
      Y_custom = Y_custom,
      ensure_disjoint_rest = ensure_disjoint_rest,
      out_mean = out_mean,
      out_scale = out_scale,
      ...
    )
  }

  if (!all(c("x", "true_Matrix") %in% names(gen))) {
    stop("generator_output must contain 'x' and 'true_Matrix'")
  }

  data <- if (is.data.frame(gen$x)) gen$x else as.data.frame(gen$x)
  graph <- gen$true_Matrix
  n_real <- nrow(data)
  p_real <- ncol(data)

  if (is.null(gsem_alpha)) {
    gsem_alpha <- max(c(1 - pnorm(n_real^(1/3) / 2), 0))
  }
  lambda2 <- 0.5 * sqrt(log(p_real) / (n_real * h_ratio))

  lts_res <- safe_algorithm_call(
    LTS_GSEM_Algorithm(
      data,
      lambda1 = 0.1^4,
      lambda2 = lambda2,
      alpha = h_ratio,
      thresh = thresh,
      graph = graph
    )
  )
  us_res <- safe_algorithm_call(
    GSEM_Algorithm(
      data,
      alpha = gsem_alpha,
      direction = gsem_direction,
      graph = graph,
      max_degree = gsem_max_degree
    )
  )
  hlsem_res <- safe_algorithm_call(
    HLSEM_Algorithm(
      data,
      sparsity_level = hlsem_sparsity,
      CV = FALSE,
      graph = graph
    )
  )
  hgsem_res <- safe_algorithm_call(
    HGSEM_Algorithm(
      data,
      alpha = gsem_alpha,
      sparsity_level = hgsem_sparsity,
      graph = graph
    )
  )
  td_res <- safe_algorithm_call(
    TD(
      data,
      q = d,
      sparsity_level = td_sparsity,
      graph = graph
    )
  )
  gds_res <- safe_algorithm_call(
    GDS_Algorithm(
      data,
      scoreName = "SEMSEV",
      pars = list(regr.pars = list()),
      check = "checkUntilFirstMinK",
      output = FALSE,
      startAt = gds_startAt,
      graph = graph
    )
  )

  results <- list(
    LTS = lts_res,
    US = us_res,
    HLSEM = hlsem_res,
    HGSEM = hgsem_res,
    TD = td_res,
    GDS = gds_res
  )

  metrics <- do.call(
    rbind,
    lapply(names(results), function(name) {
      build_metric_row(name, results[[name]])
    })
  )
  metrics$seed <- seed

  list(
    seed = seed,
    data = if (save_data) data else NULL,
    graph = graph,
    contamination = list(Y = gen$Y, summary = gen$summary),
    evaluations = metrics,
    raw_results = results
  )
}

#' Future 기반 병렬 실행으로 여러 시드를 반복 평가한다.
run_parallel_CCLSM <- function(
    seeds,
    workers = future::availableCores(),
    plan = c("multisession", "multicore", "sequential"),
    progress_format = "[:bar] :percent | elapsed: :elapsed | eta: :eta | :message",
    ...
) {
  plan <- match.arg(plan)
  total_steps <- length(seeds)
  if (total_steps == 0) {
    return(list())
  }

  old_plan <- future::plan()
  on.exit({
    future::plan(old_plan)
  }, add = TRUE)

  if (plan == "multisession") {
    future::plan(future::multisession, workers = workers)
  } else if (plan == "multicore") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::sequential)
  }

  progressr::handlers(global = TRUE, progressr::handler_progress(format = progress_format))
  start_time <- Sys.time()

  progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)
    future.apply::future_lapply(seq_along(seeds), function(idx) {
      seed <- seeds[idx]
      res <- CCLSM_simulation_fun(seed = seed, ...)
      elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      remaining_est <- max(total_steps - idx, 0)
      eta_sec <- if (idx > 0) (elapsed_sec / idx) * remaining_est else NA_real_
      msg <- sprintf(
        "seed %d | elapsed %s | eta %s",
        seed,
        format_duration(elapsed_sec),
        format_duration(eta_sec)
      )
      p(message = msg)
      res
    }, future.seed = TRUE)
  })
}

#' 시뮬레이션 결과 리스트를 단일 데이터 프레임으로 변환한다.
bind_simulation_metrics <- function(sim_results) {
  if (length(sim_results) == 0) {
    return(data.frame())
  }
  metrics <- lapply(sim_results, function(res) {
    df <- res$evaluations
    df$seed <- res$seed
    df
  })
  do.call(rbind, metrics)
}

if (interactive()) {
  example_seeds <- 1:4
  parallel_results <- run_parallel_CCLSM(
    seeds = example_seeds,
    plan = "sequential",
    workers = 1,
    n_real = n,
    p_real = p,
    d = d,
    beta_min = beta_min,
    beta_max = beta_max,
    graph_type = graph_type,
    b = b,
    outlier_nodes = outlier_node,
    h_ratio = 0.5,
    thresh = 10
  )
  simulation_summary <- bind_simulation_metrics(parallel_results)
  print(head(simulation_summary))
}
