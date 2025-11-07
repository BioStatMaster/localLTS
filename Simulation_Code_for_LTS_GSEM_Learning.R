
# Set the working directory based on the active RStudio document when possible.
default_project_root <- "D:/2022LTS_GSEM"
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
  if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
    setwd(dirname(ctx$path))
  } else if (dir.exists(default_project_root)) {
    setwd(default_project_root)
  }
} else if (dir.exists(default_project_root)) {
  setwd(default_project_root)
}

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
#' Attempt to load an optional dependency without triggering installations.
#'
#' @param package Package name as a string.
#' @param attach Whether to call `library()` when the package is available.
load_optional_package <- function(package, attach = TRUE) {
  if (requireNamespace(package, quietly = TRUE)) {
    if (attach) {
      library(package, character.only = TRUE)
    }
    TRUE
  } else {
    warning(sprintf("Package '%s' is not available. Related functionality may be skipped.", package))
    FALSE
  }
}

optional_packages <- c(
  "Rcpp", "pcalg", "MASS", "gtools", "graph", "dplyr", "plyr", "future",
  "future.apply", "progressr", "hypergeo", "glmnet", "gamlss", "bnlearn",
  "rmutil", "CombMSC"
)
loaded_optional <- vapply(optional_packages, load_optional_package, logical(1))
rm(loaded_optional, optional_packages)

library(compiler)

# Bioconductor packages are also optional; skip attempts to download in offline environments.
load_optional_package("Rgraphviz")
load_optional_package("RBGL")
load_optional_package("graph")

############# 공통 유틸리티 함수 #############

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

#########################
#### Data Generation ####
########################

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

#' CCLSM 기반 데이터 생성과 알고리즘 평가를 한 번에 수행한다.
#' 
#' @param generator_output `CCLSM_generator()`가 반환한 결과. 제공되면
#'   새로운 데이터를 생성하지 않고 이 객체를 그대로 사용한다.
#' @param fixed_graph 사전 생성한 그래프 행렬. 제공 시 동일한 구조를
#'   유지한 채 데이터를 시뮬레이션한다.
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
    fixed_graph = NULL
) {
  set.seed(seed)
  gen <- generator_output
  if (is.null(gen)) {
    if (!is.null(fixed_graph)) {
      gen <- CCLSM_generate(
        n = n_real,
        B = fixed_graph,
        var_min = var_Min,
        var_max = var_Max,
        outlier_nodes = outlier_nodes,
        b = b,
        overlap_mode = overlap_mode,
        overlap_fraction = overlap_fraction,
        union_rows = union_rows,
        Y_custom = Y_custom,
        ensure_disjoint_rest = ensure_disjoint_rest,
        dist = dist,
        out_mean = out_mean,
        out_scale = out_scale,
        seed = seed
      )
    } else {
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
        out_scale = out_scale
      )
    }
  }

  if (!all(c("x", "true_Matrix") %in% names(gen))) {
    stop("generator_output must contain 'x' and 'true_Matrix'")
  }

  data <- if (is.data.frame(gen$x)) gen$x else as.data.frame(gen$x)
  graph <- gen$true_Matrix
  n_real <- nrow(data)
  p_real <- ncol(data)
  contamination_mask <- if ("y" %in% names(gen)) gen$y else gen$Y
  if (is.null(contamination_mask)) {
    contamination_mask <- matrix(0, n_real, p_real)
  }

  contamination_nodes <- if (!is.null(outlier_nodes) && length(outlier_nodes) > 0) {
    sort(unique(outlier_nodes))
  } else {
    sort(unique(which(colSums(contamination_mask) > 0)))
  }
  contamination_details <- summarise_contamination_nodes(graph, contamination_nodes)

  format_metadata_strings <- function(details, field) {
    if (length(details) == 0) {
      return("")
    }
    paste(vapply(details, function(entry) {
      values <- entry[[field]]
      if (length(values) == 0) {
        sprintf("%d: none", entry$node)
      } else {
        sprintf("%d: %s", entry$node, paste(values, collapse = ","))
      }
    }, character(1)), collapse = " | ")
  }

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
  metrics$contaminated_nodes <- if (length(contamination_nodes) == 0) "" else paste(contamination_nodes, collapse = ",")
  metrics$Pa_j <- format_metadata_strings(contamination_details, "Pa_j")
  metrics$ancestor_union <- format_metadata_strings(contamination_details, "ancestor_union")

  list(
    seed = seed,
    data = if (save_data) data else NULL,
    graph = graph,
    contamination = list(
      Y = contamination_mask,
      summary = gen$summary,
      node_details = contamination_details
    ),
    evaluations = metrics,
    raw_results = results
  )
}

#' Future 기반 병렬 실행으로 여러 시드를 반복 평가한다.
#'
#' @param sim_args `CCLSM_simulation_fun()`에 전달할 인자 목록.
run_parallel_CCLSM <- function(
    seeds,
    workers = future::availableCores(),
    plan = c("multisession", "multicore", "sequential"),
    progress_format = "[:bar] :percent | elapsed: :elapsed | eta: :eta | :message",
    sim_args = list()
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

  handler <- progressr::handler_progress(format = progress_format)
  progressr::handlers(handler)
  start_time <- Sys.time()

  progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)
    future.apply::future_lapply(seq_along(seeds), function(idx) {
      seed <- seeds[idx]
      args <- modifyList(sim_args, list(seed = seed))
      res <- do.call(CCLSM_simulation_fun, args)
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

# 데이터 생성 및 단일 시뮬레이션 예제 -----------------------------------------

#### Simulation Settings for Toy Example #####
p <- 20
n <- 50
d <- 1
seed <- sample(1:100, 1)
beta_min <- 1.00
beta_max <- 1.00
graph_type <- 5
alpha <- max(c(1 - pnorm(n^(1/3) / 2), 0))

outlier_node <- 1
b <- 1
synthetic_data <- CCLSM_generator(
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
graph <- synthetic_data$true_Matrix
toy_data <- synthetic_data$x
Y_mask <- synthetic_data$y
str(synthetic_data$summary)
par(mfrow = c(3, 1))

colSums(graph != 0)

boxplot(toy_data)
barplot(sapply(toy_data, mean))
barplot(sapply(toy_data, var))

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
  gsem_alpha = alpha,
  generator_output = synthetic_data,
  save_data = TRUE
)
toy_metrics <- toy_simulation$evaluations
print(head(toy_metrics))

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

#' Ensure evaluation metric vectors share the same ordering and length.
#'
#' @param metrics Numeric vector returned by `evaluation_fun()`.
#' @param metric_names Character vector of metric names to enforce.
#' @return Numeric vector aligned with `metric_names`.
pad_metrics <- function(metrics, metric_names) {
  template <- setNames(rep(NA_real_, length(metric_names)), metric_names)
  if (is.null(metrics)) {
    return(template)
  }
  metrics <- metrics[intersect(names(metrics), metric_names)]
  template[names(metrics)] <- as.numeric(metrics)
  template
}

#' 주어진 인접행렬에서 특정 노드의 부모 집합을 반환한다.
#'
#' @param B 가중치 인접행렬.
#' @param node 관심 노드 번호.
#' @return 부모 노드 인덱스 벡터.
get_parents <- function(B, node) {
  which(B[node, ] != 0)
}

#' 여러 노드에 도달할 수 있는 모든 조상 노드를 계산한다.
#'
#' @param B 가중치 인접행렬.
#' @param nodes 시작 노드 인덱스 벡터.
#' @return 조상 노드 인덱스 벡터(중복 제거).
get_ancestors <- function(B, nodes) {
  if (length(nodes) == 0) {
    return(integer(0))
  }
  p <- nrow(B)
  visited <- rep(FALSE, p)
  stack <- nodes
  while (length(stack) > 0) {
    current <- stack[[1]]
    stack <- stack[-1]
    parents <- which(B[current, ] != 0)
    new_nodes <- parents[!visited[parents]]
    if (length(new_nodes) > 0) {
      visited[new_nodes] <- TRUE
      stack <- c(stack, new_nodes)
    }
  }
  which(visited)
}

#' 오염 노드, 그 부모, 그리고 해당 조상들의 정보를 정리한다.
#'
#' @param B 가중치 인접행렬.
#' @param nodes 오염된 노드 인덱스 벡터.
#' @return 각 노드별 메타데이터 리스트.
summarise_contamination_nodes <- function(B, nodes) {
  if (length(nodes) == 0) {
    return(list())
  }
  lapply(as.integer(nodes), function(node) {
    parents <- get_parents(B, node)
    base_set <- unique(c(node, parents))
    ancestor_union <- sort(unique(c(base_set, get_ancestors(B, base_set))))
    list(
      node = node,
      Pa_j = parents,
      ancestor_union = ancestor_union
    )
  })
}

#' Build a single evaluation record for saving to disk.
#' 
#' @param algorithm_result Output list from an algorithm run (e.g. LTS).
#' @param metric_names Character vector describing metric order.
#' @return List with four numeric entries so legacy loaders can index
#'   `[[1]]`, `[[2]]`, `[[3]]`, or `[[4]]` safely.
build_evaluation_record <- function(algorithm_result, metric_names) {
  dag <- pad_metrics(algorithm_result$DAG_Evaluation, metric_names)
  mec <- pad_metrics(algorithm_result$MEC_Evaluation, metric_names)
  oracle <- pad_metrics(algorithm_result$Oracle_Evaluation, metric_names)
  list(
    DAG = dag,
    MEC = mec,
    Oracle = oracle,
    Estimated_DAG = if (!is.null(algorithm_result$DAG)) algorithm_result$DAG else NULL,
    Ordering = if (!is.null(algorithm_result$Ordering)) algorithm_result$Ordering else NULL
  )
}

#' Split simulation outputs by algorithm for downstream persistence.
#'
#' @param sim_results List of outputs produced by `CCLSM_simulation_fun()`.
#' @param metric_names Character vector describing metric order.
#' @return Named list where each entry contains the per-seed records for an
#'   algorithm.
split_results_by_algorithm <- function(sim_results, metric_names) {
  if (length(sim_results) == 0) {
    return(list())
  }
  algorithms <- names(sim_results[[1]]$raw_results)
  setNames(lapply(algorithms, function(algo) {
    lapply(sim_results, function(res) {
      build_evaluation_record(res$raw_results[[algo]], metric_names)
    })
  }), algorithms)
}

# Mapping from algorithm names to file-name prefixes expected by Plotting_LTS.
method_file_prefixes <- c(
  LTS = "LTS",
  US = "USB",
  HLSEM = "HLSM",
  HGSEM = "HGSM",
  TD = "TD",
  GDS = "GDS"
)

# Null-coalescing helper for optional list elements.
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

#' Construct the output filename for a method/parameter configuration.
construct_output_filename <- function(prefix, combo, num_out) {
  if (identical(prefix, "LTS")) {
    sprintf(
      "%s_p%d_n%d_d%d_b%d_h%d_thresh%s_result.Rdata",
      prefix,
      combo$p_real,
      combo$n_real,
      combo$d,
      num_out,
      round(combo$h_ratio * 100),
      formatC(combo$thresh, format = "fg", digits = 6)
    )
  } else {
    sprintf(
      "%s_p%d_n%d_d%d_b%d_result.Rdata",
      prefix,
      combo$p_real,
      combo$n_real,
      combo$d,
      num_out
    )
  }

  progressr::handlers(global = TRUE, progressr::handler_progress(format = progress_format))
  start_time <- Sys.time()

  progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)
    future.apply::future_lapply(seq_along(seeds), function(idx) {
      seed <- seeds[idx]
      args <- modifyList(sim_args, list(seed = seed))
      res <- do.call(CCLSM_simulation_fun, args)
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

#' Persist algorithm-specific results to disk.
#'
#' @param directory Output folder.
#' @param algorithm Name of the algorithm (e.g. "LTS").
#' @param records List of per-seed evaluation records.
#' @param combo List with fields `n_real`, `p_real`, `d`, `h_ratio`, `thresh`.
#' @param num_out Number of outlier rows used in the scenario.
#' @param seeds Vector of RNG seeds.
#' @param extra_metadata Additional parameters saved alongside the results.
#' @return Full path to the saved `.Rdata` file.
save_algorithm_results <- function(
    directory,
    algorithm,
    records,
    combo,
    num_out,
    seeds,
    extra_metadata) {
  prefix <- method_file_prefixes[[algorithm]]
  if (is.null(prefix)) {
    warning(sprintf("Unknown algorithm prefix for '%s'", algorithm))
    return(invisible(NULL))
  }
  filename <- construct_output_filename(prefix, combo, num_out)
  full_path <- file.path(directory, filename)
  evaluation_result <- records
  params <- modifyList(extra_metadata, list(algorithm = algorithm))
  save(evaluation_result, seeds, params, file = full_path)
  full_path
}

#' Save all algorithm results for a single parameter combination.
#'
#' @param sim_results List of outputs from `CCLSM_simulation_fun()`.
#' @param output_dir Directory to store `.Rdata` files.
#' @param combo List describing the simulation parameters.
#' @param num_out Number of contaminated rows used in the generator.
#' @param seeds Vector of RNG seeds.
#' @param metric_names Ordering for evaluation metrics.
#' @return Named character vector of saved file paths.
persist_simulation_results <- function(
    sim_results,
    output_dir,
    combo,
    num_out,
    seeds,
    metric_names,
    extra_metadata = list()) {
  per_algorithm <- split_results_by_algorithm(sim_results, metric_names)
  saved_paths <- vapply(names(per_algorithm), function(algo) {
    records <- per_algorithm[[algo]]
    if (length(records) == 0) {
      return(NA_character_)
    }
    metadata <- modifyList(extra_metadata, list(
      n_real = combo$n_real,
      p_real = combo$p_real,
      d = combo$d,
      h_ratio = combo$h_ratio,
      thresh = combo$thresh,
      num_out = num_out
    ))
    save_algorithm_results(
      directory = output_dir,
      algorithm = algo,
      records = records,
      combo = combo,
      num_out = num_out,
      seeds = seeds,
      extra_metadata = metadata
    )
  }, character(1), USE.NAMES = TRUE)
  saved_paths
}

#' Resolve scenario-specific values that may depend on the number of nodes.
resolve_spec <- function(spec, combo) {
  if (is.function(spec)) {
    spec(combo$p_real)
  } else {
    spec
  }
  metrics <- lapply(sim_results, function(res) {
    df <- res$evaluations
    df$seed <- res$seed
    df
  })
  do.call(rbind, metrics)
}

#' Collapse node-wise contamination counts to a single scalar for naming.
collapse_num_out <- function(values) {
  if (length(values) == 0) {
    return(0)
  }
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(0)
  }
  max(values)
}

#' Summarise and persist Plotting_LTS inputs for a scenario.
#'
#' @param scenario_label Character label for the contamination scenario.
#' @param scenario_dir Directory where results were written.
#' @param combos Data frame describing the parameter grid.
#' @param d Graph indegree setting.
#' @param num_out Number of outlier rows.
#' @param seeds Vector of RNG seeds used in the simulation.
#' @param h_ratios Distinct trimming ratios.
#' @param thresh_values Distinct threshold values for LTS.
emit_plotting_summary <- function(
    scenario_label,
    scenario_dir,
    combos,
    d,
    num_out_values,
    seeds,
    h_ratios,
    thresh_values) {
  plotting_inputs <- list(
    directory = scenario_dir,
    model = scenario_label,
    N = sort(unique(combos$n_real)),
    P = sort(unique(combos$p_real)),
    d = d,
    num_out = sort(unique(num_out_values)),
    h_ratio = sort(unique(h_ratios)),
    thresh = sort(unique(thresh_values)),
    seeds = seeds
  )
  message("[Plotting Inputs] model=", scenario_label,
          " | dir=", scenario_dir,
          " | N=", paste(plotting_inputs$N, collapse = ","),
          " | P=", paste(plotting_inputs$P, collapse = ","),
          " | h_ratio=", paste(plotting_inputs$h_ratio, collapse = ","),
          " | thresh=", paste(plotting_inputs$thresh, collapse = ","),
          " | num_out=", paste(plotting_inputs$num_out, collapse = ","))
  invisible(plotting_inputs)
}

#' Run simulations across a grid of parameters and save results for plotting.
#'
#' @param seeds Vector of RNG seeds used for replication.
#' @param n_values Candidate sample sizes.
#' @param p_values Candidate variable counts.
#' @param h_ratios Vector of trimming ratios for the LTS algorithm.
#' @param thresh_values Vector of threshold values for LTS.
#' @param scenarios List describing contamination scenarios. Each scenario should
#'   provide `label`, `b`, and `outlier_nodes` (scalars or functions of `p`).
#'   Optional `sim_args` entries override generator parameters.
#' @param common_args Additional arguments passed to `CCLSM_simulation_fun`.
#' @param output_root Root directory where scenario sub-folders will be created.
#' @param workers Number of future workers.
#' @param plan Future plan ("multisession", "multicore", or "sequential").
#' @return Nested list containing the saved file paths for each scenario.
run_simulation_grid <- function(
    seeds,
    n_values,
    p_values,
    h_ratios,
    thresh_values,
    scenarios,
    common_args,
    output_root,
    workers = future::availableCores(),
    plan = c("multisession", "multicore", "sequential")) {
  plan <- match.arg(plan)
  if (is.null(common_args$d)) {
    stop("common_args$d must be supplied for batch simulations")
  }
  combos <- expand.grid(
    n_real = n_values,
    p_real = p_values,
    h_ratio = h_ratios,
    thresh = thresh_values,
    stringsAsFactors = FALSE
  )
  combos$d <- common_args$d

  saved_paths <- list()
  for (scenario in scenarios) {
    scenario_label <- scenario$label
    scenario_dir <- file.path(output_root, scenario_label)
    dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
    num_out_values <- vapply(seq_len(nrow(combos)), function(i) {
      collapse_num_out(resolve_spec(scenario$b, combos[i, ]))
    }, numeric(1))

    scenario_paths <- list()
    for (idx in seq_len(nrow(combos))) {
      combo <- combos[idx, ]
      scenario_outliers <- resolve_spec(scenario$outlier_nodes, combo)
      scenario_b <- resolve_spec(scenario$b, combo)
      scenario_args <- modifyList(common_args, scenario$sim_args %||% list())
      sim_args <- modifyList(scenario_args, list(
        n_real = combo$n_real,
        p_real = combo$p_real,
        h_ratio = combo$h_ratio,
        thresh = combo$thresh,
        outlier_nodes = scenario_outliers,
        b = scenario_b
      ))

      sim_results <- run_parallel_CCLSM(
        seeds = seeds,
        workers = workers,
        plan = plan,
        sim_args = sim_args
      )

      combo_list <- as.list(combo)
      saved <- persist_simulation_results(
        sim_results = sim_results,
        output_dir = scenario_dir,
        combo = combo_list,
        num_out = collapse_num_out(scenario_b),
        seeds = seeds,
        metric_names = evaluation_metric_names
      )
      scenario_paths[[paste0("n", combo$n_real, "_p", combo$p_real, "_h", combo$h_ratio, "_t", combo$thresh)]] <- saved
    }

    emit_plotting_summary(
      scenario_label = scenario_label,
      scenario_dir = scenario_dir,
      combos = combos,
      d = common_args$d,
      num_out_values = num_out_values,
      seeds = seeds,
      h_ratios = h_ratios,
      thresh_values = thresh_values
    )
    saved_paths[[scenario_label]] <- scenario_paths
  }
  saved_paths
}

#' 고정된 그래프에서 노드별 오염을 순차적으로 평가하는 배치 실행기.
#'
#' @param sem_seed 그래프 생성을 위한 고정 시드.
#' @param seeds 데이터/알고리즘 반복을 위한 시드 벡터.
#' @param n_real 표본 크기.
#' @param p_real 변수 개수.
#' @param h_ratios LTS 절삭비 목록.
#' @param thresh_values LTS 임계값 목록.
#' @param b 각 노드에서 오염시킬 행 수.
#' @param contaminated_nodes 실험할 노드 인덱스(미지정 시 1:p_real).
#' @param common_args `CCLSM_simulation_fun()`에 전달할 공통 인자 리스트.
#' @param output_root 결과 저장 루트 경로.
#' @param workers Future 워커 수.
#' @param plan Future 실행 방식.
#' @return 노드별로 저장된 파일 경로 목록.
run_nodewise_contamination_batch <- function(
    sem_seed,
    seeds,
    n_real,
    p_real,
    h_ratios,
    thresh_values,
    b,
    contaminated_nodes = NULL,
    common_args,
    output_root,
    workers = future::availableCores(),
    plan = c("multisession", "multicore", "sequential")) {
  plan <- match.arg(plan)
  if (is.null(contaminated_nodes)) {
    contaminated_nodes <- seq_len(p_real)
  }

  d_val <- common_args$d
  if (is.null(d_val)) {
    stop("common_args$d must be supplied")
  }
  var_min <- common_args$var_Min %||% 1
  var_max <- common_args$var_Max %||% 1
  dist_choice <- common_args$dist %||% "Gaussian"
  beta_min <- common_args$beta_min %||% 0.8
  beta_max <- common_args$beta_max %||% 1.2
  graph_type <- common_args$graph_type %||% 5
  overlap_mode <- common_args$overlap_mode %||% "independent"
  overlap_fraction <- common_args$overlap_fraction %||% 0
  ensure_disjoint <- common_args$ensure_disjoint_rest %||% TRUE
  out_mean_val <- common_args$out_mean %||% 1e3
  out_scale_val <- common_args$out_scale %||% 1

  base_generator <- CCLSM_generator(
    n = n_real,
    p = p_real,
    d = d_val,
    b = 0,
    outlier_nodes = NULL,
    var_Min = var_min,
    var_Max = var_max,
    dist = dist_choice,
    beta_min = beta_min,
    beta_max = beta_max,
    graph_type = graph_type,
    seed = sem_seed,
    overlap_mode = overlap_mode,
    overlap_fraction = overlap_fraction,
    ensure_disjoint_rest = ensure_disjoint,
    out_mean = out_mean_val,
    out_scale = out_scale_val
  )
  fixed_graph <- base_generator$true_Matrix

  common_args_expanded <- modifyList(common_args, list(
    d = d_val,
    var_Min = var_min,
    var_Max = var_max,
    dist = dist_choice,
    beta_min = beta_min,
    beta_max = beta_max,
    graph_type = graph_type,
    overlap_mode = overlap_mode,
    overlap_fraction = overlap_fraction,
    ensure_disjoint_rest = ensure_disjoint,
    out_mean = out_mean_val,
    out_scale = out_scale_val
  ))

  combos <- expand.grid(
    h_ratio = h_ratios,
    thresh = thresh_values,
    stringsAsFactors = FALSE
  )
  combos$n_real <- n_real
  combos$p_real <- p_real
  combos$d <- common_args$d

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  node_results <- list()

  for (node in contaminated_nodes) {
    node_label <- sprintf("Node%d", node)
    node_dir <- file.path(output_root, node_label)
    dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)

    node_meta <- summarise_contamination_nodes(fixed_graph, node)
    node_detail <- if (length(node_meta) > 0) {
      node_meta[[1]]
    } else {
      list(node = node, Pa_j = integer(0), ancestor_union = node)
    }
    num_out_value <- if (length(b) == 1) {
      b
    } else {
      idx <- which(contaminated_nodes == node)
      if (length(idx) == 0) b[[1]] else b[[idx[1]]]
    }

    scenario_paths <- list()
    for (idx in seq_len(nrow(combos))) {
      combo_row <- combos[idx, ]
      sim_args <- modifyList(common_args_expanded, list(
        n_real = n_real,
        p_real = p_real,
        h_ratio = combo_row$h_ratio,
        thresh = combo_row$thresh,
        outlier_nodes = node,
        b = num_out_value,
        fixed_graph = fixed_graph
      ))

      sim_results <- run_parallel_CCLSM(
        seeds = seeds,
        workers = workers,
        plan = plan,
        sim_args = sim_args
      )

      combo_list <- as.list(combo_row)
      saved <- persist_simulation_results(
        sim_results = sim_results,
        output_dir = node_dir,
        combo = combo_list,
        num_out = num_out_value,
        seeds = seeds,
        metric_names = evaluation_metric_names,
        extra_metadata = list(
          sem_seed = sem_seed,
          contaminated_node = node,
          Pa_j = node_detail$Pa_j,
          ancestor_union = node_detail$ancestor_union
        )
      )
      scenario_paths[[paste0("h", combo_row$h_ratio, "_t", combo_row$thresh)]] <- saved
    }

    emit_plotting_summary(
      scenario_label = node_label,
      scenario_dir = node_dir,
      combos = combos,
      d = common_args$d,
      num_out_values = rep(num_out_value, nrow(combos)),
      seeds = seeds,
      h_ratios = h_ratios,
      thresh_values = thresh_values
    )

    node_results[[node_label]] <- scenario_paths
  }

  node_results
}

#########################
#### Batch Simulation ####
#########################

# Default output root aligns with Plotting_LTS expectations.
batch_output_root <- file.path(default_project_root, "Result")

# Common parameters shared across all batch jobs. Adjust as needed.
batch_common_args <- list(
  d = d,
  beta_min = beta_min,
  beta_max = beta_max,
  graph_type = graph_type,
  var_Min = 0.75,
  var_Max = 0.75,
  dist = "Gaussian",
  gsem_alpha = alpha,
  overlap_mode = "independent",
  ensure_disjoint_rest = TRUE
)

# Parameter grids for batch execution.
batch_seeds <- 1:10
batch_n_values <- c(50, 100, 150, 200)
batch_p_values <- c(20)
batch_h_ratios <- c(0.5, 0.7, 0.9)
batch_thresh_values <- c(2, 10)

# Contamination scenarios mirror the plotting folders (Outlier1 / OutlierAll).
batch_scenarios <- list(
  list(
    label = "Outlier1",
    b = 1,
    outlier_nodes = function(p) 1
  ),
  list(
    label = "OutlierAll",
    b = 1,
    outlier_nodes = function(p) seq_len(p)
  )
)

# Fraction-overlap scenarios contaminate every node with a shared `b` while
# varying the overlap fraction. Scenario labels encode the percentage so the
# resulting folders remain Plotting_LTS friendly (e.g. `OutlierAll_fraction100`).
fraction_overlap_values <- c(1, 0.75, 0.5, 0.25, 0)
fraction_b_value <- 1
fraction_scenarios <- lapply(fraction_overlap_values, function(frac) {
  list(
    label = sprintf("OutlierAll_fraction%03d", as.integer(round(frac * 100))),
    b = fraction_b_value,
    outlier_nodes = function(p) seq_len(p),
    sim_args = list(
      overlap_mode = "fraction",
      overlap_fraction = frac
    )
  )
})

# Toggle to launch the batch run when sourcing the script.
run_batch_now <- FALSE
if (isTRUE(run_batch_now)) {
  dir.create(batch_output_root, recursive = TRUE, showWarnings = FALSE)
  batch_results <- run_simulation_grid(
    seeds = batch_seeds,
    n_values = batch_n_values,
    p_values = batch_p_values,
    h_ratios = batch_h_ratios,
    thresh_values = batch_thresh_values,
    scenarios = batch_scenarios,
    common_args = batch_common_args,
    output_root = batch_output_root,
    workers = future::availableCores(),
    plan = "multisession"
  )
}

# Toggle for the fraction-overlap scenarios that contaminate every node using
# the same `b` value while sweeping overlap fractions.
run_fraction_batch_now <- FALSE
if (isTRUE(run_fraction_batch_now)) {
  dir.create(batch_output_root, recursive = TRUE, showWarnings = FALSE)
  fraction_batch_results <- run_simulation_grid(
    seeds = batch_seeds,
    n_values = batch_n_values,
    p_values = batch_p_values,
    h_ratios = batch_h_ratios,
    thresh_values = batch_thresh_values,
    scenarios = fraction_scenarios,
    common_args = batch_common_args,
    output_root = batch_output_root,
    workers = future::availableCores(),
    plan = "multisession"
  )
}

# 노드별 오염 시나리오 설정.
nodewise_sem_seed <- 777
nodewise_n <- 100
nodewise_p <- 20
nodewise_b <- 1
nodewise_nodes <- NULL  # NULL이면 1:nodewise_p 전체를 탐색한다.
nodewise_h_ratios <- batch_h_ratios
nodewise_thresh_values <- batch_thresh_values
nodewise_output_root <- file.path(batch_output_root, "Nodewise")
run_nodewise_now <- FALSE
if (isTRUE(run_nodewise_now)) {
  dir.create(nodewise_output_root, recursive = TRUE, showWarnings = FALSE)
  nodewise_results <- run_nodewise_contamination_batch(
    sem_seed = nodewise_sem_seed,
    seeds = batch_seeds,
    n_real = nodewise_n,
    p_real = nodewise_p,
    h_ratios = nodewise_h_ratios,
    thresh_values = nodewise_thresh_values,
    b = nodewise_b,
    contaminated_nodes = nodewise_nodes,
    common_args = batch_common_args,
    output_root = nodewise_output_root,
    workers = future::availableCores(),
    plan = "multisession"
  )
}

if (interactive()) {
  example_seeds <- 1:4
  parallel_results <- run_parallel_CCLSM(
    seeds = example_seeds,
    plan = "sequential",
    workers = 1,
    sim_args = list(
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
  )
  simulation_summary <- bind_simulation_metrics(parallel_results)
  print(head(simulation_summary))
}
