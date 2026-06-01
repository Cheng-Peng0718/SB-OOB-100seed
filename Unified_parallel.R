## ============================================================
## Unified_parallel.R — Job Array Version (Seed 1-100)
## 并行架构完全保留：ForkCluster + foreach %dorng%
## Usage: Rscript --vanilla Unified_parallel.R <seed>
## ============================================================

suppressPackageStartupMessages({
  library(rpart)
  library(mlbench)
  library(MASS)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(foreach)
  library(doParallel)
  library(doRNG)
})

## ============================================================
## 0) 从命令行读取 seed（由 SLURM array task ID 传入）
## ============================================================
args        <- commandArgs(trailingOnly = TRUE)
MASTER_SEED <- if (length(args) >= 1L) as.integer(args[1L]) else 1L
cat(sprintf("[INFO] MASTER_SEED = %d\n", MASTER_SEED))

## ============================================================
## 1) 固定参数
## ============================================================
N_RUNS   <- 50
N_BOOT   <- 100
VERBOSE  <- TRUE
FUN_FILE <- file.path(getwd(), "Unified_funs.R")

## ============================================================
## 2) 加载函数库
## ============================================================
if (!file.exists(FUN_FILE))
  stop("Unified_funs.R not found in: ", getwd())
source(FUN_FILE)

needed <- c("fit_cart","predict_cart","node_ids",
            "bag_cart_OOB","bag_cart_SB","oob_predict_train","bag_predict",
            "E1E2_node_prob","E1E2_node_err",
            "test_error_class","test_mse","lower_bound_class","lower_bound_reg",
            "gen_friedman1","gen_friedman2","gen_friedman3",
            "load_breast","load_pima","load_vehicle","load_satellite","load_dna",
            "load_boston","load_ozone")
miss <- setdiff(needed, ls(.GlobalEnv))
if (length(miss))
  stop("Missing from Unified_funs.R:\n", paste(miss, collapse = ", "))

## ============================================================
## 3) 输出目录（每个 seed 独立子目录）
## ============================================================
OUT_DIR <- file.path(getwd(), "results", sprintf("seed_%03d", MASTER_SEED))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Output dir: %s\n", OUT_DIR))

## ============================================================
## 4) 并行设置 — ForkCluster（限定在 agx/agg 节点上运行）
##    worker 数取 min(SLURM_CPUS_PER_TASK, N_RUNS) 避免空闲进程
## ============================================================
Sys.setenv(OMP_NUM_THREADS        = "1",
           OPENBLAS_NUM_THREADS   = "1",
           MKL_NUM_THREADS        = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS    = "1")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "0"))
n_workers  <- if (slurm_cpus > 0L) {
  cat(sprintf("[INFO] SLURM_CPUS_PER_TASK = %d\n", slurm_cpus))
  min(slurm_cpus, N_RUNS)   # 不超过任务数，避免空闲worker
} else {
  np <- parallel::detectCores(logical = FALSE)
  min(max(1L, if (is.na(np)) parallel::detectCores() else np), N_RUNS)
}
cat(sprintf("[INFO] Launching %d parallel workers (ForkCluster)\n", n_workers))

cl <- parallel::makeForkCluster(n_workers)
doParallel::registerDoParallel(cl)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
set.seed(MASTER_SEED)

## ============================================================
## 5) 工具函数 (v2)
## ============================================================
.save_exp <- function(res, expname) {
  cat(sprintf("\n[RESULTS] %s:\n", expname))
  print(res)
  fpath <- file.path(OUT_DIR, sprintf("%s.RData", expname))
  save(res, file = fpath)
  cat(sprintf("[INFO] Saved: %s\n", fpath))
  invisible(res)
}

## v2 新增：把原始 runs (N_RUNS × n_metrics × n_datasets) 也保存下来
## 以后可做层次模型 / 重新统计 / 任意 ad-hoc 分析
.save_raw_runs <- function(raw, expname) {
  fpath <- file.path(OUT_DIR, sprintf("%s_raw.RData", expname))
  save(raw, file = fpath)
  cat(sprintf("[INFO] Saved raw runs: %s (%d rows)\n", fpath, nrow(raw)))
  invisible(raw)
}

## v2 新增：U_b 摘要（验证 SB 真的固定了 distinct count）
.save_ub_summary <- function(ub_tab, expname) {
  fpath <- file.path(OUT_DIR, sprintf("%s_ub.RData", expname))
  save(ub_tab, file = fpath)
  cat(sprintf("[INFO] Saved U_b summary: %s (%d rows)\n", fpath, nrow(ub_tab)))
  invisible(ub_tab)
}


## 汇总 N_RUNS 行 -> 配对统计
## ------------------------------------------------------------------
## BUGFIX: 旧版在 summarise() 里先把 OOB / SB_OOB 用 mean() 改写成标量，
##         之后再调 wilcox.test(OOB, SB_OOB, paired=TRUE) —— 此时
##         OOB/SB_OOB 已是单一数值，p 值恒为 1.0，毫无信息。
##         dplyr summarise 是顺序求值的，覆盖名后下游表达式拿到的是新值。
##
## 修法：  ① 所有依赖向量的统计量（wilcox/t/CI/d_z/sign）放在前面，
##           此时 OOB、SB_OOB 还是 N_RUNS 长度的向量；
##         ② 再把 OOB、SB_OOB 用 mean() 折叠成标量。
##         （在外面定义辅助函数 .paired_stats 进一步降低重复出错的概率。）
## ------------------------------------------------------------------

## 辅助：从两个向量得到所有配对统计量，返回 1 行 tibble
.paired_stats <- function(x_oob, x_sb) {
  d <- x_sb - x_oob
  n <- sum(!is.na(d))
  out <- list(
    n_runs    = n,
    mean_diff = mean(d, na.rm = TRUE),
    sd_diff   = sd(d,   na.rm = TRUE),
    wilcox_p  = tryCatch(
      wilcox.test(x_oob, x_sb, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_),
    t_p       = tryCatch(
      t.test(x_oob, x_sb, paired = TRUE)$p.value,
      error = function(e) NA_real_),
    ci_lo     = tryCatch(t.test(d)$conf.int[1], error = function(e) NA_real_),
    ci_hi     = tryCatch(t.test(d)$conf.int[2], error = function(e) NA_real_),
    cohens_dz = if (isTRUE(sd(d, na.rm = TRUE) > 0))
                  mean(d, na.rm = TRUE) / sd(d, na.rm = TRUE) else NA_real_,
    n_pos     = sum(d >  0, na.rm = TRUE),
    n_neg     = sum(d <  0, na.rm = TRUE),
    n_zero    = sum(d == 0, na.rm = TRUE)
  )
  np <- out$n_pos; nn <- out$n_neg
  out$sign_p <- if ((np + nn) > 0)
    tryCatch(binom.test(min(np, nn), np + nn, p = 0.5)$p.value,
             error = function(e) NA_real_) else NA_real_
  tibble::as_tibble(out)
}

## 全局容器：每次 .summarise_runs 调用时，把原始 50 个 inner run
## 的配对数据 (OOB, SB_OOB, diff per metric per dataset) 收集起来。
## 主程序在每个 EXP 开始前设置 .CURRENT_EXP，结束后从这里把原始数据
## 取出来另存为 *_raw.RData。
.RAW_RUNS <- new.env(parent = emptyenv())

.summarise_runs <- function(runs, dataset_name, dataset_type) {
  ## --- v2 新增：先把原始 runs 存到全局容器 ---
  if (exists(".CURRENT_EXP", envir = .RAW_RUNS)) {
    cur <- get(".CURRENT_EXP", envir = .RAW_RUNS)
    raw_chunk <- runs %>%
      mutate(dataset = dataset_name, type = dataset_type, .before = 1)
    if (exists(cur, envir = .RAW_RUNS)) {
      assign(cur, bind_rows(get(cur, envir = .RAW_RUNS), raw_chunk),
             envir = .RAW_RUNS)
    } else {
      assign(cur, raw_chunk, envir = .RAW_RUNS)
    }
  }

  ## --- v2 修复 + 增强的配对统计 (已说明) ---
  runs %>%
    group_by(metric) %>%
    reframe(
      stats  = list(.paired_stats(OOB, SB_OOB)),
      OOB    = mean(OOB,    na.rm = TRUE),
      SB_OOB = mean(SB_OOB, na.rm = TRUE)
    ) %>%
    tidyr::unnest(stats) %>%
    mutate(dataset = dataset_name, type = dataset_type, .before = 1)
}

## ------------------------------------------------------------------
## 健全性检查：所有 wilcox_p 都 == 1 几乎一定是回归到旧 bug。
## main pipeline 跑完后，调一次 .sanity_check_pvalues(res) 即可。
## ------------------------------------------------------------------
.sanity_check_pvalues <- function(res) {
  if (!"wilcox_p" %in% names(res)) return(invisible(res))
  ps <- res$wilcox_p[!is.na(res$wilcox_p)]
  if (length(ps) > 0 && all(abs(ps - 1) < 1e-12)) {
    warning("[SANITY] 所有 wilcox_p 都 == 1，summarise_runs 可能又被改坏了！")
  }
  invisible(res)
}

## ============================================================
## EXP1 — Node class probability: E1_B / E2_B
## 并行：foreach %dorng% 跑 N_RUNS 次，每次独立生成数据+建模
## ============================================================
exp1_once_compare <- function(train_df, test_df, yname) {
  ref      <- fit_cart(train_df, yname, method = "class")
  nodes_tr <- paste0("node_", node_ids(ref, train_df))
  nodes_te <- paste0("node_", node_ids(ref, test_df))
  levs     <- levels(train_df[[yname]])

  p_star <- tapply(seq_len(nrow(test_df)), nodes_te, function(ix)
    as.numeric(prop.table(table(factor(test_df[[yname]][ix], levels = levs)))))
  q_star <- setNames(as.numeric(prop.table(table(nodes_te))), names(table(nodes_te)))

  .bag_probs <- function(bag) {
    pr <- oob_predict_train(bag, train_df)$probs
    tapply(seq_len(nrow(train_df)), nodes_tr, function(ix)
      colMeans(pr[ix, , drop = FALSE], na.rm = TRUE))
  }
  .met <- function(p_B) {
    cmn <- Reduce(intersect, list(names(p_star), names(p_B), names(q_star)))
    E1E2_node_prob(p_star[cmn], p_B[cmn], q_star[cmn])
  }

  p_B_oob <- .bag_probs(bag_cart_OOB(train_df, yname, n_boot = N_BOOT))
  p_B_sb  <- .bag_probs(bag_cart_SB(train_df,  yname, n_boot = N_BOOT))
  m_oob   <- .met(p_B_oob)
  m_sb    <- .met(p_B_sb)

  tibble(metric = c("E1_B","E2_B"),
         OOB    = c(m_oob$E1, m_oob$E2),
         SB_OOB = c(m_sb$E1,  m_sb$E2),
         diff   = SB_OOB - OOB)
}

run_exp1_all <- function() {
  syn_specs <- list(
    list(name="waveform",  gen=function(n) mlbench.waveform(n),  ntr=300, nte=5000),
    list(name="twonorm",   gen=function(n) mlbench.twonorm(n),   ntr=200, nte=5000),
    list(name="threenorm", gen=function(n) mlbench.threenorm(n), ntr=200, nte=5000),
    list(name="ringnorm",  gen=function(n) mlbench.ringnorm(n),  ntr=200, nte=5000)
  )
  res_all <- NULL
  for (i in seq_along(syn_specs)) {
    sp <- syn_specs[[i]]
    if (VERBOSE) cat("\n[EXP1] Synthetic:", sp$name, "\n")
    set.seed(MASTER_SEED + 10000L * i)
    ## ← 多核并行：N_RUNS 次分布在 n_workers 个 fork worker 上同时跑
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble","mlbench","rpart")) %dorng% {
      tr    <- sp$gen(sp$ntr); te <- sp$gen(sp$nte)
      df_tr <- as_tibble(tr$x); df_tr$y <- factor(tr$classes)
      df_te <- as_tibble(te$x); df_te$y <- factor(te$classes, levels = levels(df_tr$y))
      exp1_once_compare(df_tr, df_te, "y")
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, sp$name, "synthetic"))
  }

  reals <- list(
    list(name="breast-cancer", loader=load_breast,    y="Class"),
    list(name="diabetes",      loader=load_pima,      y="diabetes"),
    list(name="vehicle",       loader=load_vehicle,   y="Class"),
    list(name="satellite",     loader=load_satellite, y="classes"),
    list(name="dna",           loader=load_dna,       y="class")
  )
  for (j in seq_along(reals)) {
    it <- reals[[j]]
    if (VERBOSE) cat("\n[EXP1] Real:", it$name, "\n")
    df <- it$loader()
    set.seed(MASTER_SEED + 100000L + 10000L * j)
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble","mlbench","rpart")) %dorng% {
      idx   <- unlist(tapply(seq_len(nrow(df)), df[[it$y]], function(ix)
        sample(ix, size = floor(0.7 * length(ix)))))
      tr    <- df[idx, ]; te <- df[-idx, ]
      tr[[it$y]] <- factor(tr[[it$y]])
      te[[it$y]] <- factor(te[[it$y]], levels = levels(tr[[it$y]]))
      exp1_once_compare(tr, te, it$y)
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, it$name, "real"))
  }
  res_all
}

## ============================================================
## EXP2 — Regression node error: EB1 / EB2
## ============================================================
exp2_once_compare <- function(train_df, test_df, yname = "y") {
  ref      <- fit_cart(train_df, yname, method = "anova")
  nodes_tr <- paste0("node_", node_ids(ref, train_df))
  nodes_te <- paste0("node_", node_ids(ref, test_df))

  node_mu <- tapply(seq_len(nrow(train_df)), nodes_tr,
                    function(ix) mean(train_df[[yname]][ix]))
  e_star  <- tapply(seq_len(nrow(test_df)), nodes_te, function(ix) {
    nd <- as.character(nodes_te[ix][1])
    sqrt(mean((test_df[[yname]][ix] - node_mu[[nd]])^2))
  })
  q_star <- setNames(as.numeric(prop.table(table(nodes_te))), names(table(nodes_te)))

  .node_rmse <- function(bag) {
    preds <- oob_predict_train(bag, train_df)$preds
    if (anyNA(preds))
      preds[is.na(preds)] <- bag_predict(
        bag, train_df[is.na(preds), , drop = FALSE], type = "vector")
    res2 <- (train_df[[yname]] - preds)^2
    tapply(seq_len(nrow(train_df)), nodes_tr,
           function(ix) sqrt(mean(res2[ix])))
  }
  .met <- function(e_B) {
    cmn <- intersect(names(e_star), names(e_B))
    E1E2_node_err(e_star[cmn], e_B[cmn], q_star[cmn])
  }

  m_oob <- .met(.node_rmse(bag_cart_OOB(train_df, yname, n_boot = N_BOOT)))
  m_sb  <- .met(.node_rmse(bag_cart_SB(train_df,  yname, n_boot = N_BOOT)))

  tibble(metric = c("EB1","EB2"),
         OOB    = c(m_oob$E1, m_oob$E2),
         SB_OOB = c(m_sb$E1,  m_sb$E2),
         diff   = SB_OOB - OOB)
}

run_exp2_all <- function() {
  reg_specs <- list(
    list(name="friedman1", gen=gen_friedman1, ntr=100, nte=2000),
    list(name="friedman2", gen=gen_friedman2, ntr=100, nte=2000),
    list(name="friedman3", gen=gen_friedman3, ntr=100, nte=2000)
  )
  res_all <- NULL
  for (i in seq_along(reg_specs)) {
    sp <- reg_specs[[i]]
    if (VERBOSE) cat("\n[EXP2] Synthetic:", sp$name, "\n")
    set.seed(MASTER_SEED + 200000L + 10000L * i)
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble")) %dorng% {
      tr <- sp$gen(sp$ntr); te <- sp$gen(sp$nte)
      ## NOTE: gen_friedman1/2/3 已返回 tibble(...,y), 不要再 unpack
      exp2_once_compare(tr, te, "y")
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, sp$name, "synthetic"))
  }

  for (info in list(
    list(name="Boston", loader=load_boston, seed_off=240000L),
    list(name="Ozone",  loader=load_ozone,  seed_off=250000L)
  )) {
    if (VERBOSE) cat("\n[EXP2] Real:", info$name, "\n")
    df <- info$loader()
    set.seed(MASTER_SEED + info$seed_off)
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble")) %dorng% {
      idx <- sample.int(nrow(df), size = floor(0.7 * nrow(df)))
      exp2_once_compare(df[idx, ], df[-idx, ], "y")
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, info$name, "real"))
  }
  res_all
}

## ============================================================
## EXP3 — Within-node variability: R1..R4
## ============================================================
exp3_once_compare <- function(gen_fun, n_train, n_test = 5000) {
  tr    <- gen_fun(n_train); te <- gen_fun(n_test)
  df_tr <- as_tibble(tr$x); df_tr$y <- factor(tr$classes)
  df_te <- as_tibble(te$x); df_te$y <- factor(te$classes, levels = levels(df_tr$y))
  yname <- "y"; levs <- levels(df_tr$y)

  ref      <- fit_cart(df_tr, yname, method = "class")
  nodes_tr <- paste0("node_", node_ids(ref, df_tr))
  nodes_te <- paste0("node_", node_ids(ref, df_te))

  p_star_t <- tapply(seq_len(nrow(df_te)), nodes_te, function(ix)
    as.numeric(prop.table(table(factor(df_te$y[ix], levels = levs)))))
  q_star_t <- setNames(as.numeric(prop.table(table(nodes_te))), names(table(nodes_te)))
  p_R      <- tapply(seq_len(nrow(df_tr)), nodes_tr, function(ix)
    as.numeric(prop.table(table(factor(df_tr$y[ix], levels = levs)))))

  Z            <- model.matrix(~ y - 1, data = df_te); colnames(Z) <- levs
  p_star_t_mat <- t(sapply(as.character(nodes_te), function(nd) p_star_t[[nd]]))

  ENR <- sum(sapply(intersect(names(p_star_t), names(p_R)), function(nd)
    sum((p_star_t[[nd]] - p_R[[nd]])^2) * q_star_t[[nd]]))
  EV  <- mean(rowSums((Z - p_star_t_mat)^2))

  .R_metrics <- function(bag) {
    pr   <- oob_predict_train(bag, df_tr)$probs
    p_B  <- tapply(seq_len(nrow(df_tr)), nodes_tr, function(ix)
      colMeans(pr[ix, , drop = FALSE], na.rm = TRUE))
    pb_te <- bag_predict(bag, df_te, type = "prob")
    ENB   <- sum(sapply(intersect(names(p_star_t), names(p_B)), function(nd)
      sum((p_star_t[[nd]] - p_B[[nd]])^2) * q_star_t[[nd]]))
    EB_pw <- mean(rowSums((Z - as.matrix(pb_te[, levs, drop = FALSE]))^2))
    c(R1 = 100 * ENR / (ENR + EV),
      R2 = 100 * (ENB - ENR) / (ENR + EV),
      R3 = 100 * ENB / (ENB + EV),
      R4 = 100 * EB_pw / EV)
  }

  R_oob <- .R_metrics(bag_cart_OOB(df_tr, yname, n_boot = N_BOOT))
  R_sb  <- .R_metrics(bag_cart_SB(df_tr,  yname, n_boot = N_BOOT))

  tibble(metric = names(R_oob),
         OOB    = as.numeric(R_oob),
         SB_OOB = as.numeric(R_sb),
         diff   = SB_OOB - OOB)
}

run_exp3_all <- function() {
  syn_specs <- list(
    list(name="waveform",  gen=function(n) mlbench.waveform(n),  ntr=300),
    list(name="twonorm",   gen=function(n) mlbench.twonorm(n),   ntr=200),
    list(name="threenorm", gen=function(n) mlbench.threenorm(n), ntr=200),
    list(name="ringnorm",  gen=function(n) mlbench.ringnorm(n),  ntr=200)
  )
  res_all <- NULL
  for (i in seq_along(syn_specs)) {
    sp <- syn_specs[[i]]
    if (VERBOSE) cat("\n[EXP3] Synthetic:", sp$name, "\n")
    set.seed(MASTER_SEED + 300000L + 10000L * i)
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble","mlbench","rpart")) %dorng% {
      exp3_once_compare(sp$gen, sp$ntr, 5000)
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, sp$name, "synthetic"))
  }
  res_all
}

## ============================================================
## EXP4 — Generalization error: eTS / eOB / absdiff / ratio
## ============================================================
exp4_once_compare_class <- function(train_df, test_df, yname) {
  .run <- function(bag) {
    yhat  <- bag_predict(bag, test_df, type = "class")
    eTS   <- test_error_class(test_df[[yname]], yhat)
    probs <- oob_predict_train(bag, train_df)$probs
    cls   <- colnames(probs)[max.col(probs, ties.method = "first")]
    eOB   <- mean(train_df[[yname]] != factor(cls, levels = levels(train_df[[yname]])))
    lb    <- lower_bound_class(p_hat   = eTS,
                               N_train = nrow(train_df),
                               N_test  = nrow(test_df))
    c(eTS = eTS, eOB = eOB,
      absdiff = abs(eTS - eOB),
      ratio   = abs(eTS - eOB) / max(lb, .Machine$double.eps))
  }
  r_oob <- .run(bag_cart_OOB(train_df, yname, n_boot = N_BOOT))
  r_sb  <- .run(bag_cart_SB(train_df,  yname, n_boot = N_BOOT))
  tibble(metric = names(r_oob),
         OOB    = as.numeric(r_oob),
         SB_OOB = as.numeric(r_sb),
         diff   = SB_OOB - OOB)
}

exp4_once_compare_reg <- function(train_df, test_df, yname = "y") {
  .run <- function(bag) {
    yhat   <- bag_predict(bag, test_df, type = "vector")
    eTS    <- test_mse(test_df[[yname]], yhat)
    preds  <- oob_predict_train(bag, train_df)$preds
    if (anyNA(preds))
      preds[is.na(preds)] <- bag_predict(
        bag, train_df[is.na(preds), , drop = FALSE], type = "vector")
    eOB    <- mean((train_df[[yname]] - preds)^2)
    resid2 <- (test_df[[yname]] - yhat)^2
    lb     <- lower_bound_reg(resid2,
                              N_train = nrow(train_df),
                              N_test  = nrow(test_df))
    c(eTS = eTS, eOB = eOB,
      absdiff = abs(eTS - eOB),
      ratio   = abs(eTS - eOB) / max(lb, .Machine$double.eps))
  }
  r_oob <- .run(bag_cart_OOB(train_df, yname, n_boot = N_BOOT))
  r_sb  <- .run(bag_cart_SB(train_df,  yname, n_boot = N_BOOT))
  tibble(metric = names(r_oob),
         OOB    = as.numeric(r_oob),
         SB_OOB = as.numeric(r_sb),
         diff   = SB_OOB - OOB)
}

run_exp4_all <- function() {
  res_all <- NULL

  if (VERBOSE) cat("\n[EXP4] Classification: twonorm\n")
  set.seed(MASTER_SEED + 400000L)
  runs <- foreach(r = seq_len(N_RUNS),
                  .combine  = dplyr::bind_rows,
                  .export   = EXPORTS,
                  .packages = c("dplyr","tibble","mlbench","rpart")) %dorng% {
    tr    <- mlbench.twonorm(200); te <- mlbench.twonorm(5000)
    df_tr <- as_tibble(tr$x); df_tr$y <- factor(tr$classes)
    df_te <- as_tibble(te$x); df_te$y <- factor(te$classes, levels = levels(df_tr$y))
    exp4_once_compare_class(df_tr, df_te, "y")
  }
  res_all <- bind_rows(res_all, .summarise_runs(runs, "twonorm", "class"))

  if (VERBOSE) cat("\n[EXP4] Regression: friedman1\n")
  set.seed(MASTER_SEED + 410000L)
  runs <- foreach(r = seq_len(N_RUNS),
                  .combine  = dplyr::bind_rows,
                  .export   = EXPORTS,
                  .packages = c("dplyr","tibble")) %dorng% {
    tr <- gen_friedman1(100); te <- gen_friedman1(2000)
    exp4_once_compare_reg(tr, te, "y")
  }
  res_all <- bind_rows(res_all, .summarise_runs(runs, "friedman1", "reg"))
  res_all
}

## ============================================================
## EXP5 — Trees using OOB outputs: mse_original / mse_oob_outputs
## ============================================================
exp5_once_compare <- function(train_df, test_df, yname = "y") {
  ref     <- fit_cart(train_df, yname, method = "anova")
  mse_ref <- test_mse(test_df[[yname]],
                      predict_cart(ref, test_df, type = "vector"))
  .oob_mse <- function(bag) {
    preds <- oob_predict_train(bag, train_df)$preds
    if (anyNA(preds))
      preds[is.na(preds)] <- bag_predict(
        bag, train_df[is.na(preds), , drop = FALSE], type = "vector")
    df_new          <- train_df
    df_new[[yname]] <- preds
    refB <- fit_cart(df_new, yname, method = "anova")
    test_mse(test_df[[yname]], predict_cart(refB, test_df, type = "vector"))
  }
  mse_oob <- .oob_mse(bag_cart_OOB(train_df, yname, n_boot = N_BOOT))
  mse_sb  <- .oob_mse(bag_cart_SB(train_df,  yname, n_boot = N_BOOT))
  tibble(metric = c("mse_original","mse_oob_outputs"),
         OOB    = c(mse_ref, mse_oob),
         SB_OOB = c(mse_ref, mse_sb),
         diff   = SB_OOB - OOB)
}

run_exp5_all <- function() {
  specs <- list(
    list(name="friedman1", gen=gen_friedman1, ntr=200, nte=2000),
    list(name="friedman2", gen=gen_friedman2, ntr=200, nte=2000),
    list(name="friedman3", gen=gen_friedman3, ntr=200, nte=2000)
  )
  res_all <- NULL
  for (i in seq_along(specs)) {
    sp <- specs[[i]]
    if (VERBOSE) cat("\n[EXP5] Synthetic:", sp$name, "\n")
    set.seed(MASTER_SEED + 500000L + 10000L * i)
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble")) %dorng% {
      tr <- sp$gen(sp$ntr); te <- sp$gen(sp$nte)
      exp5_once_compare(tr, te, "y")
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, sp$name, "reg"))
  }

  for (info in list(
    list(name="Boston", loader=load_boston, seed_off=540000L),
    list(name="Ozone",  loader=load_ozone,  seed_off=550000L)
  )) {
    if (VERBOSE) cat("\n[EXP5] Real:", info$name, "\n")
    df <- info$loader()
    set.seed(MASTER_SEED + info$seed_off)
    runs <- foreach(r = seq_len(N_RUNS),
                    .combine  = dplyr::bind_rows,
                    .export   = EXPORTS,
                    .packages = c("dplyr","tibble")) %dorng% {
      idx <- sample.int(nrow(df), size = floor(0.9 * nrow(df)))
      exp5_once_compare(df[idx, ], df[-idx, ], "y")
    }
    res_all <- bind_rows(res_all, .summarise_runs(runs, info$name, "reg"))
  }
  res_all
}

## ============================================================
## 6) EXPORTS — 必须在所有函数定义之后、foreach 调用之前
## ============================================================
EXPORTS <- ls(.GlobalEnv)

## ============================================================
## 7) 运行全部实验，每个实验完成后立即保存
## ============================================================
cat(sprintf("\n=== SEED %d | N_RUNS=%d | N_BOOT=%d | Workers=%d ===\n",
            MASTER_SEED, N_RUNS, N_BOOT, n_workers))

## v2 辅助：跑一个实验，存 summary + raw runs，并做 sanity check
.run_and_save <- function(expname, runner) {
  assign(".CURRENT_EXP", expname, envir = .RAW_RUNS)
  res <- runner()
  .save_exp(res, expname)
  .sanity_check_pvalues(res)
  if (exists(expname, envir = .RAW_RUNS)) {
    raw <- get(expname, envir = .RAW_RUNS)
    .save_raw_runs(raw, expname)
  }
  res
}

exp1_res <- .run_and_save("exp1", run_exp1_all)
exp2_res <- .run_and_save("exp2", run_exp2_all)
exp3_res <- .run_and_save("exp3", run_exp3_all)
exp4_res <- .run_and_save("exp4", run_exp4_all)
exp5_res <- .run_and_save("exp5", run_exp5_all)

## ============================================================
## 8) 统一保存这个 seed 的全部结果
## ============================================================
all_results <- list(seed = MASTER_SEED,
                    exp1 = exp1_res,
                    exp2 = exp2_res,
                    exp3 = exp3_res,
                    exp4 = exp4_res,
                    exp5 = exp5_res)
final_path <- file.path(OUT_DIR, "all_results.RData")
save(all_results, file = final_path)
cat(sprintf("\n[INFO] Seed %d complete. Full results: %s\n",
            MASTER_SEED, final_path))