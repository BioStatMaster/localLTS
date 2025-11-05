library(rstudioapi)
install.packages("languageserver")
path <-setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## ============================================================
## CCLSM Generator with Overlap Controls (B 입력 필수)
## ============================================================
#TODO: 모듈화
CCLSM_generate <- function(
    n,
    B,                              # 반드시 입력 (Graph_Generator()의 B)
    var_min = 1, var_max = 1,       # 잡음 표준편차 범위
    outlier_nodes = NULL,           # 오염시킬 노드들 (예: c(2,5))
    b = 0,                          # 노드별 오염셀 수(스칼라 or 노드수 길이)
    overlap_mode = c("independent","shared","fraction","custom_union","custom_mask"),
    overlap_fraction = 0,           # fraction 모드에서 공유 비율 rho∈[0,1]
    union_rows = NULL,              # custom_union 모드에서 합집합 행
    Y_custom = NULL,                # custom_mask 모드에서 사용자 마스크 (n×p, 0/1)
    ensure_disjoint_rest = TRUE,    # fraction: 공유 이외는 겹치지 않게 배정 시도
    dist = c("Gaussian","t","Laplace","Uniform"),
    out_mean = 1e3, out_scale = 1,  # 오염분포 파라미터
    seed = 1
){
  set.seed(seed)
  overlap_mode <- match.arg(overlap_mode)
  dist <- match.arg(dist)
  
  if (!is.matrix(B)) stop("B must be a matrix")
  p <- ncol(B)
  
  ## ---------- 깨끗한 X ----------
  noise_sd <- runif(p, min = var_min, max = var_max)
  X <- matrix(0, n, p)
  for (i in 1:p) {
    pa <- if (i > 1) which(B[i, 1:(i-1), drop=FALSE] != 0) else integer(0)
    mu <- if (length(pa) == 0) 0 else X[, pa, drop=FALSE] %*% B[i, pa]
    X[, i] <- rnorm(n, mean = as.vector(mu), sd = noise_sd[i])
  }
  
  ## ---------- 오염 마스크 Y ----------
  Y <- matrix(0L, n, p)
  
  ## (A) custom_mask는 최우선 처리 (outlier_nodes 없어도 동작)
  if (overlap_mode == "custom_mask") {
    if (is.null(Y_custom)) stop("Provide Y_custom for custom_mask mode")
    if (!all(dim(Y_custom) == c(n, p))) stop("Y_custom must be n x p")
    if (!all(Y_custom %in% c(0,1))) stop("Y_custom must be 0/1")
    Y <- Y_custom
    
    ## (B) 나머지 모드들은 outlier_nodes가 있어야 의미가 있음
  } else if (!is.null(outlier_nodes) && length(outlier_nodes) > 0) {
    
    if (is.numeric(b)) {
      if (length(b) == 1L) b <- rep(b, length(outlier_nodes))
      if (length(b) != length(outlier_nodes)) stop("length(b) must equal length(outlier_nodes) or be scalar")
      if (any(b > n)) stop("b cannot exceed n")
    } else stop("b must be numeric")
    
    rlaplace <- function(m, loc = 0, scale = 1){
      u <- runif(m, -0.5, 0.5)
      loc - scale * sign(u) * log(1 - 2*abs(u))
    }
    gen_outliers <- function(m){
      switch(dist,
             "Gaussian" = rnorm(m, mean = out_mean, sd = out_scale),
             "t"        = out_mean + out_scale * rt(m, df = 3),
             "Laplace"  = rlaplace(m, loc = out_mean, scale = out_scale),
             "Uniform"  = runif(m, min = out_mean - out_scale, max = out_mean + out_scale)
      )
    }
    
    if (overlap_mode == "independent") {
      for (idx in seq_along(outlier_nodes)) {
        k <- outlier_nodes[idx]; bk <- b[idx]
        if (bk > 0) Y[sample.int(n, bk), k] <- 1L
      }
      
    } else if (overlap_mode == "shared") {
      b_shared <- min(b)
      S <- if (b_shared > 0) sample.int(n, b_shared) else integer(0)
      for (k in outlier_nodes) if (length(S)) Y[S, k] <- 1L
      for (idx in seq_along(outlier_nodes)) {
        k <- outlier_nodes[idx]; extra <- b[idx] - b_shared
        if (extra > 0) {
          rest <- setdiff(seq_len(n), S)
          Y[sample(rest, extra), k] <- 1L
        }
      }
      
    } else if (overlap_mode == "fraction") {
      if (overlap_fraction < 0 || overlap_fraction > 1) stop("overlap_fraction must be in [0,1]")
      base_b <- min(b)
      b_shared <- round(overlap_fraction * base_b)
      S <- if (b_shared > 0) sample.int(n, b_shared) else integer(0)
      for (k in outlier_nodes) if (length(S)) Y[S, k] <- 1L
      for (idx in seq_along(outlier_nodes)) {
        k <- outlier_nodes[idx]
        need <- b[idx] - b_shared
        if (need > 0) {
          if (ensure_disjoint_rest) {
            already <- which(rowSums(Y) > 0)
            candidates <- setdiff(seq_len(n), union(S, already))
            take <- min(length(candidates), need)
            sel  <- if (take > 0) sample(candidates, take) else integer(0)
            if (length(sel)) Y[sel, k] <- 1L
            rem <- need - take
            if (rem > 0) {
              rest <- setdiff(seq_len(n), c(S, sel))
              Y[sample(rest, rem), k] <- 1L
            }
          } else {
            rest <- setdiff(seq_len(n), S)
            Y[sample(rest, need), k] <- 1L
          }
        }
      }
      
    } else if (overlap_mode == "custom_union") {
      if (is.null(union_rows)) stop("Provide union_rows for custom_union mode")
      union_rows <- sort(unique(union_rows))
      if (!all(union_rows %in% seq_len(n))) stop("union_rows must be within 1..n")
      for (idx in seq_along(outlier_nodes)) {
        k <- outlier_nodes[idx]; bk <- b[idx]
        if (bk > length(union_rows)) stop("b cannot exceed length(union_rows)")
        if (bk > 0) Y[sample(union_rows, bk), k] <- 1L
      }
    }
    
    # gen_outliers는 아래에서 사용하려고 외부로 빼고자 하면
    # <<- 로 빼거나 함수 밖 재정의가 필요합니다. 여기서는 아래에서 다시 정의합니다.
  }
  
  ## ---------- 오염 적용 ----------
  # (gen_outliers가 위 블록 내부였던 분기 때문에 여기서 다시 정의)
  rlaplace <- function(m, loc = 0, scale = 1){
    u <- runif(m, -0.5, 0.5)
    loc - scale * sign(u) * log(1 - 2*abs(u))
  }
  gen_outliers <- function(m){
    switch(dist,
           "Gaussian" = rnorm(m, mean = out_mean, sd = out_scale),
           "t"        = out_mean + out_scale * rt(m, df = 3),
           "Laplace"  = rlaplace(m, loc = out_mean, scale = out_scale),
           "Uniform"  = runif(m, min = out_mean - out_scale, max = out_mean + out_scale)
    )
  }
  
  Z <- X
  idx_out <- which(Y == 1L, arr.ind = TRUE)
  if (nrow(idx_out) > 0) Z[idx_out] <- gen_outliers(nrow(idx_out))
  
  ## 요약
  per_node <- colSums(Y)
  per_row  <- rowSums(Y)
  union_size <- sum(per_row > 0)
  
  list(
    X_clean = as.data.frame(X),
    Z = as.data.frame(Z),
    Y = Y,
    B = B,
    summary = list(
      per_node   = per_node,
      per_row    = per_row,
      union_size = union_size
    )
  )
}


## ---------------------------
## 위상정렬 (B[child, parent] ≠ 0)
## ---------------------------
compute_topological_order <- function(B) {
  p <- ncol(B)
  A <- matrix(0L, p, p)               # A[u,v]=1 if u→v
  for (u in 1:p) {
    children <- which(B[, u] != 0)
    if (length(children) > 0) A[u, children] <- 1L
  }
  indeg <- colSums(A)
  order <- integer(0)
  S <- sort(which(indeg == 0))
  while (length(S) > 0) {
    u <- S[1]; S <- S[-1]
    order <- c(order, u)
    outs <- which(A[u, ] != 0)
    for (v in outs) {
      indeg[v] <- indeg[v] - 1L
      if (indeg[v] == 0L) S <- sort(c(S, v))
    }
  }
  if (length(order) != p) {
    warning("Topological order not found (cycle suspected). Fallback to 1:p.")
    order <- 1:p
  }
  order
}

## ---------------------------
## 조상 집합
## ---------------------------
get_ancestors <- function(B, nodes) {
  p <- ncol(B)
  A <- (B != 0)      # child x parent
  visited <- rep(FALSE, p)
  stack <- as.integer(unique(nodes))
  while (length(stack) > 0) {
    v <- stack[1]; stack <- stack[-1]
    parents <- which(A[v, ])
    parents <- parents[!visited[parents]]
    if (length(parents) > 0) {
      visited[parents] <- TRUE
      stack <- c(stack, parents)
    }
  }
  which(visited)
}

## ---------------------------
## T_j(r) 계산: T_j(r)={pi_1,...,pi_{p+1-r}}\{j}
## ---------------------------
get_Tj_from_graph <- function(B, j, r, topo_order = NULL) {
  p <- ncol(B)
  if (is.null(topo_order)) topo_order <- compute_topological_order(B)
  stopifnot(r >= 1, r <= p)
  prefix_len   <- p + 1 - r
  prefix_nodes <- topo_order[seq_len(prefix_len)]
  Tj <- setdiff(prefix_nodes, j)
  list(Tj_index = Tj, topo_order = topo_order,
       prefix_len = prefix_len, prefix_nodes = prefix_nodes)
}

get_Pa_from_graph <- function(B, j) {
  # B: adjacency/weight matrix where B[child, parent] != 0
  # j: scalar or vector of node indices
  if (!is.matrix(B)) stop("B must be a matrix")
  p <- ncol(B)
  if (any(j < 1 | j > p)) stop("j must be within 1..ncol(B)")
  
  if (length(j) == 1L) {
    pa_idx <- which(B[j, ] != 0)
    list(Pa_index = pa_idx,
         Pa_weights = if (length(pa_idx)) as.numeric(B[j, pa_idx]) else numeric(0))
  } else {
    # return a named list per requested node
    res <- lapply(j, function(jj) {
      pa_idx <- which(B[jj, ] != 0)
      list(Pa_index = pa_idx,
           Pa_weights = if (length(pa_idx)) as.numeric(B[jj, pa_idx]) else numeric(0))
    })
    names(res) <- as.character(j)
    res
  }

suppressPackageStartupMessages({ library(robustHD) })  # sparseLTS
LTS_GSEM_Algorithm(..)
eval_single_node_LTS <- function(Z, B_true, Y, j,
                                 Tj = NULL, r = NULL, topo_order = NULL,
                                 alpha = 0.75, lambda = 1e-4) {
  Z <- as.matrix(Z)
  n <- nrow(Z); p <- ncol(Z)
  # Pa_j(r) 결정
  Paj <- get_Pa_from_graph(B_true, j)
  Paj_idx <- Paj$Pa_index
  Paj_weights <- Paj$Pa_weights

  # T_j(r) 결정
  Tj_info <- NULL; Tj_source <- NULL
  if (!is.null(Tj)) {
    Tj <- as.integer(Tj); Tj_source <- "explicit_Tj"
  } else if (!is.null(r)) {
    Tj_info <- get_Tj_from_graph(B_true, j, r, topo_order)
    Tj <- Tj_info$Tj_index
    Tj_source <- sprintf("graph+r (r=%d)", r)
  } else {
    Tj <- setdiff(seq_len(p), j)              # fallback
    Tj_source <- "all_minus_j"
  }
  # M_j = {j} ∪ Pa_j(r)
  # S_j = {j} ∪ T_j ∪ Anc({j} ∪ T_j)
  Mj_core <- unique(c(j, Paj_idx))
  Panc_set <- get_ancestors(B_true, Mj_core)
  Mj_cols <- sort(unique(c(Mj_core, Panc_set)))
  Sj_core <- unique(c(j, Tj))
  Anc_set <- get_ancestors(B_true, Sj_core)
  Sj_cols <- sort(unique(c(Sj_core, Anc_set)))
  
  # 오염(행) 집합: S_j 열 중 하나라도 1인 행
  if (!is.null(Y)) {
    affected_idx <- which(rowSums(Y[, Sj_cols, drop = FALSE] != 0) > 0)
    affected_size <- length(affected_idx)
    cell_budget_total   <- sum(Y)
    cell_budget_per_node <- colSums(Y)
    cell_budget_per_row  <- rowSums(Y)
  } else {
    affected_idx <- integer(0); affected_size <- 0
    cell_budget_total <- NA_integer_
    cell_budget_per_node <- rep(NA_integer_, p)
    cell_budget_per_row  <- rep(NA_integer_, n)
  }
  
  # 이론 임계: n - floor(alpha n) + 1
  h_j <- floor(alpha * n)
  nodewise_threshold <- n - h_j + 1
  theory_hit <- affected_size >= nodewise_threshold
  
  # 진실/추정 계수 및 파라미터 오차
  beta_true <- as.numeric(B_true[j, Tj])
  if (length(Tj) == 0L) {
    beta_hat <- numeric(0)
    param_mse <- 0; param_mae <- 0
    support_true <- integer(0); support_hat <- integer(0)
  } else {
    x_mat <- Z[, Tj, drop = FALSE]  # 반드시 오염 데이터 Z
    y_vec <- Z[, j]
    fit <- robustHD::sparseLTS(
      x = x_mat, y = y_vec,
      alpha = alpha, lambda = lambda, mode = "lambda",
      intercept = FALSE
    )
    beta_hat <- as.numeric(fit$coefficients)
    param_mse <- sum((beta_hat - beta_true)^2)
    param_mae <- sum(abs(beta_hat - beta_true))
    support_true <- which(beta_true != 0)
    support_hat  <- which(beta_hat  != 0)
  }
  
  support_FP_idx <- setdiff(support_hat,  support_true)
  support_FN_idx <- setdiff(support_true, support_hat)
  
  list(
    j = j, r = if (is.null(r)) NA_integer_ else r,
    Paj_index = Paj_idx,
    Paj_weights = Paj_weights,
    Mj_cols = Mj_cols,
    Tj_source = Tj_source,
    topo_order = if (is.null(Tj_info)) topo_order else Tj_info$topo_order,
    prefix_len = if (is.null(Tj_info)) NA_integer_ else Tj_info$prefix_len,
    prefix_nodes = if (is.null(Tj_info)) NULL else Tj_info$prefix_nodes,
    Tj_index = Tj,
    Sj_cols = Sj_cols,
    beta_true = beta_true,
    beta_hat  = beta_hat,
    ParamMSE  = param_mse,
    ParamMAE  = param_mae,
    support_FP_idx = support_FP_idx,
    support_FN_idx = support_FN_idx,
    support_FP = length(support_FP_idx),
    support_FN = length(support_FN_idx),
    affected_idx = affected_idx,
    affected_size = affected_size,
    nodewise_threshold = nodewise_threshold,
    theory_hit = theory_hit,
    cell_budget_total = cell_budget_total,
    cell_budget_per_node = cell_budget_per_node,
    cell_budget_per_row  = cell_budget_per_row
  )
}

## CCLSM_generate() 결과를 바로 넣는 어댑터
eval_j_from_gen <- function(gen, j, r = 1, alpha = 0.75, lambda = 1e-4,
                            Tj = NULL, topo_order = NULL) {
  stopifnot(all(c("Z","Y","B") %in% names(gen)))
  eval_single_node_LTS(
    Z = gen$Z, B_true = gen$B, Y = gen$Y,
    j = j, Tj = Tj, r = r, topo_order = topo_order,
    alpha = alpha, lambda = lambda
  )
}

summarize_Y <- function(Y){
  stopifnot(is.matrix(Y))
  per_node <- colSums(Y)
  per_row  <- rowSums(Y)
  union_size <- sum(per_row > 0)
  pair_overlap <- t(Y) %*% Y
  list(per_node = per_node, per_row = per_row, union_size = union_size,
       pair_overlap = pair_overlap)
}

plot_Y_mask <- function(Y, main = "Y mask (white=contaminated)"){
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = c(4,4,2,1))
  image(t(Y[nrow(Y):1, ]), axes = FALSE, col = c("black","white"),
        main = main)
  axis(1, at = seq(0,1,length.out = ncol(Y)), labels = 1:ncol(Y),
       tick = FALSE, line = -1)
  axis(2, at = seq(0,1,length.out = 5),
       labels = round(seq(nrow(Y),1,length.out=5)), las = 1)
}

## 예시: 그래프 생성 → CCLSM 데이터 생성(독립/공유/부분공유/합집합고정) → 시각화
set.seed(1)
p <- 8; d <- 2
G <- Graph_Generator(p, d, graph_type=5, beta_min=0.8, beta_max=1.2, seed=1)
B <- G$B

n <- 200
nodes <- c(2,5); b_vec <- c(20,20)

## (A) independent
exA <- CCLSM_generate(n=n, B=B,
                      outlier_nodes = nodes, b = b_vec,
                      overlap_mode = "independent", seed = 11)
exA
print(exA$summary); plot_Y_mask(exA$Y, "independent")

## (B) shared
exB <- CCLSM_generate(n=n, B=B,
                      outlier_nodes = nodes, b = b_vec,
                      overlap_mode = "shared", seed = 12)
print(exB$summary); plot_Y_mask(exB$Y, "shared")

## (C) fraction (rho=0.5)
exC <- CCLSM_generate(n=n, B=B,
                      outlier_nodes = nodes, b = b_vec,
                      overlap_mode = "fraction", overlap_fraction = 0.5,
                      ensure_disjoint_rest = TRUE, seed = 13)
print(exC$summary); plot_Y_mask(exC$Y, "fraction (rho=0.5)")

## (D) custom_union (합집합 크기 정확히 제어)
U <- sort(sample.int(n, 30))  # |union|=30
exD <- CCLSM_generate(n=n, B=B,
                      outlier_nodes = nodes, b = b_vec,
                      overlap_mode = "custom_union", union_rows = U, seed = 14)
print(exD$summary); plot_Y_mask(exD$Y, "custom_union (|union|=30)")

## j=6, r=4에 대해 그래프에서 T_j(r) 자동 계산 후 평가
set.seed(101)
gen <- CCLSM_generate(n=300, B=B,
                      outlier_nodes=c(4,8), b=40,
                      overlap_mode="shared", seed=101)
res <- eval_j_from_gen(gen, j=6, r=4, alpha=0.75, lambda=1e-4)

str(res[c("j","r","Tj_source","Tj_index","Sj_cols",
          "ParamMSE","ParamMAE","support_FP","support_FN",
          "affected_size","nodewise_threshold","theory_hit",
          "cell_budget_total")])

simulate_node_break_curve <- function(
  n=400, B, j=10, r=5,
  alpha=0.75, lambda=1e-4,
  m_grid=seq(0, 200, by=10),
  reps=20, seed=1
){
  set.seed(seed)
  p <- ncol(B)

  Tinfo  <- get_Tj_from_graph(B, j, r)
  out <- vector("list", length(m_grid)*reps); idx <- 1L

  for (m in m_grid) for (rep in 1:reps) {
    Y_custom <- matrix(0L, n, p)
    if (m > 0) Y_custom[sample.int(n, m), j] <- 1L

    gen <- CCLSM_generate(n = n, B = B,
                          overlap_mode = "custom_mask", Y_custom = Y_custom,
                          dist = "Gaussian", out_mean = 1e3, out_scale = 1,
                          seed = seed + rep)

    res <- eval_single_node_LTS(Z=gen$Z, B_true=B, Y=gen$Y,
                                j=j, Tj=Tinfo$Tj_index,
                                alpha=alpha, lambda=lambda)

    out[[idx]] <- data.frame(
      m = m, rep = rep,
      ParamMSE = res$ParamMSE, ParamMAE = res$ParamMAE,
      FP = res$support_FP, FN = res$support_FN,
      affected_size = res$affected_size,
      nodewise_threshold = res$nodewise_threshold,
      theory_hit = res$theory_hit
    ); idx <- idx + 1L
  }
  do.call(rbind, out)
}

## 실행 & 집계 & 간단 ## S_j(r) 열들 중 하나의 열에 m개 행을 표시하는 0/1 마스크
## S_j(r) 열 중 하나에 m개 행을 찍는 유틸(기본: j열)
make_Y_for_node <- function(n, p, Sj_cols, m, target_cols) {
  Y <- matrix(0L, n, p)
  if (m > 0) {
    S <- sample.int(n, m)
    for (cc in target_cols) Y[S, cc] <- 1L
  }
  Y
}

## contam ∈ {"y_only","x_leverage","xy_aligned"}
simulate_node_break_curve <- function(
    n, B, j, r,
    alpha = 0.75, lambda = 1e-4,
    m_grid = seq(0, floor(n/2), by = 10),
    reps = 20, seed = 1,
    contam = c("y_only","x_leverage","xy_aligned"),
    k_pred = 1,
    align_values = TRUE, Mx = 1e6, My = 1e6
){
  contam <- match.arg(contam)
  set.seed(seed)
  p <- ncol(B)
  
  Tinfo <- get_Tj_from_graph(B, j, r)
  Sj    <- sort(unique(c(j, Tinfo$Tj_index, get_ancestors(B, c(j, Tinfo$Tj_index)))))
  
  out <- vector("list", length(m_grid) * reps); k <- 1L
  
  for (m in m_grid) for (rep in 1:reps) {
    ## 어떤 열을 더럽힐지 결정
    if (contam == "y_only") {
      target_cols <- j
      
    } else if (contam == "x_leverage") {
      if (length(Tinfo$Tj_index) == 0L) {
        target_cols <- integer(0)
      } else {
        target_cols <- sample(Tinfo$Tj_index, min(k_pred, length(Tinfo$Tj_index)))
      }
      
    } else if (contam == "xy_aligned") {
      if (length(Tinfo$Tj_index) == 0L) {
        target_cols <- j
      } else {
        preds <- sample(Tinfo$Tj_index, min(k_pred, length(Tinfo$Tj_index)))
        target_cols <- c(j, preds)
      }
    }
    
    ## Y 마스크 생성
    Yc <- make_Y_for_node(n, p, Sj_cols = Sj, m = m, target_cols = target_cols)
    
    ## 데이터 생성(오염 주입)
    gen <- CCLSM_generate(
      n = n, B = B,
      overlap_mode = "custom_mask", Y_custom = Yc,
      dist = "Gaussian", out_mean = 1e3, out_scale = 1,
      seed = seed + rep
    )
    
    ## xy_aligned의 경우 값 정렬까지 강제(더 강력한 ‘힌지’ 유도)
    if (m > 0 && contam == "xy_aligned" && align_values && length(target_cols) > 0) {
      rows <- which(rowSums(Yc[, target_cols, drop = FALSE]) > 0)
      pred_cols <- setdiff(target_cols, j)
      if (length(pred_cols) > 0) gen$Z[rows, pred_cols] <- Mx  # 큰 레버리지 X
      gen$Z[rows, j] <- My                                    # 큰 반응 y
    }
    
    ## j-회귀 평가
    res <- eval_single_node_LTS(
      Z = gen$Z, B_true = B, Y = gen$Y,
      j = j, Tj = Tinfo$Tj_index, alpha = alpha, lambda = lambda
    )
    
    out[[k]] <- data.frame(
      m = m, rep = rep,
      ParamMSE = res$ParamMSE, ParamMAE = res$ParamMAE,
      FP = res$support_FP, FN = res$support_FN,
      affected_size = res$affected_size,
      nodewise_threshold = res$nodewise_threshold,
      theory_hit = res$theory_hit
    ); k <- k + 1L
  }
  
  do.call(rbind, out)
}


suppressPackageStartupMessages(library(robustHD))

## 그래프 만들기
set.seed(123)
G  <- Graph_Generator(p = 12, d = 2, graph_type = 5, beta_min = 0.8, beta_max = 1.2, seed = 123)
B  <- G$B

## (1) y-only: 평평할 수 있음
df_y <- simulate_node_break_curve(
  n = 400, B = B, j = 10, r = 5,
  alpha = 0.75, m_grid = seq(0, 200, by = 20),
  reps = 20, seed = 1, contam = "y_only"
)

## (2) x-leverage: 임계 부근 상승
df_x <- simulate_node_break_curve(
  n = 400, B = B, j = 10, r = 5,
  alpha = 0.75, m_grid = seq(0, 200, by = 20),
  reps = 20, seed = 1, contam = "x_leverage", k_pred = 1
)

## (3) xy-aligned: 가장 뚜렷한 ‘힌지’
df_xy <- simulate_node_break_curve(
  n = 400, B = B, j = 10, r = 5,
  alpha = 0.55, m_grid = seq(0, 200, by = 20),
  reps = 20, seed = 1, contam = "xy_aligned",
  k_pred = 2, align_values = TRUE, Mx = 1e6, My = 1e6
)

## 임계 확인
unique(df_xy$nodewise_threshold)

## 평균 ParamMSE vs m
agg <- aggregate(ParamMSE ~ m, df_xy, mean)
plt<-plot(agg$m, agg$ParamMSE, type = "b", pch = 19,
     main = "Nodewise breakdown (xy-aligned leverage)",
     xlab = "affected_size = m", ylab = "mean ParamMSE")
abline(v = unique(df_xy$nodewise_threshold)[1], lty = 2)


## S_j(r) 내부 특정 열(들) 또는 외부 열(들)에
## 동일한 m개 행을 표시하는 0/1 마스크 생성
make_Y_for_cols <- function(n, p, cols, m, share_rows = TRUE) {
  Y <- matrix(0L, n, p)
  if (m <= 0 || length(cols) == 0L) return(Y)
  if (share_rows) {
    S <- sample.int(n, m)
    Y[S, cols] <- 1L
  } else {
    ## 분산 배치가 필요하면(보통은 share_rows=TRUE가 해석상 더 명확)
    S <- sample.int(n, m)
    chosen <- sample(cols, length(S), replace = TRUE)
    for (k in seq_along(S)) Y[S[k], chosen[k]] <- 1L
  }
  Y
}

## ------------------------------------------------------------
## B) S_j 내부 vs 외부 오염 비교
##   - 같은 m에 대해 inside(S_j) vs outside(S_j^c)
##   - 반환: 두 케이스 각각의 지표(ParamMSE/MAE, FP/FN, affected_size 등)
## ------------------------------------------------------------
compare_inside_vs_outside_Sj <- function(
  n, B, j, r,
  m = 120,
  alpha = 0.75, lambda = 1e-4,
  reps = 50, seed = 1,
  contam = c("y_only","x_leverage","xy_aligned"),
  k_pred = 1,             # x_leverage/xy_aligned에서 오염할 predictor 수
  align_values = TRUE,    # xy_aligned에서 값까지 강제 정렬
  Mx = 1e6, My = 1e6,     # 강한 값
  out_mean = 1e3, out_scale = 1
){
  set.seed(seed)
  contam <- match.arg(contam)
  p <- ncol(B)

  ## S_j(r) 계산
  Tinfo <- get_Tj_from_graph(B, j, r)
  Sj    <- sort(unique(c(j, Tinfo$Tj_index, get_ancestors(B, c(j, Tinfo$Tj_index)))))
  Sjc   <- setdiff(1:p, Sj)

  if (length(Sjc) == 0L) {
    warning("S_j의 보완집합이 비어 있어 outside 케이스를 구성할 수 없습니다.")
  }

  ## 내부/외부에서 실제 오염에 사용할 열들 선택 규칙
  pick_target_cols <- function(where = c("inside","outside")) {
    where <- match.arg(where)
    if (where == "inside") {
      if (contam == "y_only") {
        return(j)
      } else if (contam == "x_leverage") {
        if (length(Tinfo$Tj_index) == 0L) return(integer(0))
        return(sample(Tinfo$Tj_index, min(k_pred, length(Tinfo$Tj_index))))
      } else { # xy_aligned
        if (length(Tinfo$Tj_index) == 0L) return(j)
        preds <- sample(Tinfo$Tj_index, min(k_pred, length(Tinfo$Tj_index)))
        return(c(j, preds))
      }
    } else { # outside
      ## S_j^c에서 하나(또는 몇 개) 선택
      if (length(Sjc) == 0L) return(integer(0))
      ## y_only/xy_aligned라도 j를 더럽히지 않게 주의!
      ## outside는 반드시 S_j^c에서만 선택
      if (contam == "x_leverage") {
        return(sample(Sjc, min(k_pred, length(Sjc))))
      } else {
        ## y_only/xy_aligned라도 y=j는 금지 → S_j^c 중에서 선택
        return(sample(Sjc, 1))
      }
    }
  }

  run_case <- function(where) {
    resL <- vector("list", reps)
    for (rep in 1:reps) {
      target_cols <- pick_target_cols(where)
      ## 표시 마스크 (같은 m개 행을 같은 열들에 공유로 주는 편이 비교가 명확)
      Yc <- make_Y_for_cols(n, p, cols = target_cols, m = m, share_rows = TRUE)

      ## 데이터 생성 (custom_mask로 정확히 적용)
      gen <- CCLSM_generate(
        n = n, B = B,
        overlap_mode = "custom_mask", Y_custom = Yc,
        dist = "Gaussian", out_mean = out_mean, out_scale = out_scale,
        seed = seed + rep
      )

      ## xy_aligned면 값 자체도 강제 정렬(더 강한 ‘힌지’ 유도)
      if (m > 0 && contam == "xy_aligned" && align_values && length(target_cols) > 0) {
        rows <- which(rowSums(Yc[, target_cols, drop = FALSE]) > 0)
        pred_cols <- setdiff(target_cols, j)
        if (length(pred_cols) > 0) gen$Z[rows, pred_cols] <- Mx
        ## where == "outside"에서도 j는 포함되지 않음(설계상!)
        if ("outside" == where) {
          ## outside에서는 y=j를 건드리지 않음 (국소성 확인)
        } else {
          gen$Z[rows, j] <- My
        }
      }

      ## LTS 평가 (항상 같은 T_j(r) 사용)
      rj <- eval_single_node_LTS(
        Z = gen$Z, B_true = B, Y = gen$Y,
        j = j, Tj = Tinfo$Tj_index, alpha = alpha, lambda = lambda
      )

      resL[[rep]] <- c(
        ParamMSE = rj$ParamMSE, ParamMAE = rj$ParamMAE,
        FP = rj$support_FP, FN = rj$support_FN,
        affected_size = rj$affected_size,
        nodewise_threshold = rj$nodewise_threshold,
        theory_hit = as.integer(rj$theory_hit)
      )
    }
    as.data.frame(do.call(rbind, resL))
  }

  inside  <- run_case("inside")
  outside <- if (length(Sjc) > 0L) run_case("outside") else NULL

  list(
    Sj = Sj, Sjc = Sjc,
    inside = inside,
    outside = outside,
    threshold = n - floor(alpha*n) + 1
  )
}

## 그래프 준비
set.seed(202)
G <- Graph_Generator(p = 30, d = 3, graph_type = 5, beta_min = 0.8, beta_max = 1.2, seed = 202)
B <- G$B

## j=10, r=5, m=120에서 inside vs outside 비교
cmp <- compare_inside_vs_outside_Sj(
  n = 400, B = B, j = 10, r = 5, m = 120,
  alpha = 0.75, reps = 50, contam = "xy_aligned", k_pred = 1,
  align_values = TRUE, Mx = 1e6, My = 1e6, out_mean = 1e3
)

## 결과 요약
sapply(cmp[c("inside","outside")], colMeans, na.rm = TRUE)
cmp$threshold
