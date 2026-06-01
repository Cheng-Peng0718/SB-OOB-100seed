## ============================================================
## Unified_funs.R — functions shared by Unified_parallel.R
## ============================================================

suppressPackageStartupMessages({
  library(rpart)
  library(mlbench)
  library(MASS)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
})

## ---------------- Utilities ----------------

is_classification <- function(y) is.factor(y) || is.character(y)
bootstrap_indices <- function(n) sample.int(n, size = n, replace = TRUE)

fit_cart <- function(df, yname, method = c("class","anova")) {
  y <- df[[yname]]
  method <- if (is_classification(y)) "class" else "anova"
  f <- as.formula(paste0(yname, " ~ ."))
  model <- rpart::rpart(
    f, data = df, method = method,
    control = rpart.control(cp = 0, xval = 0, minsplit = 2)
  )
  model$method <- method
  model
}

predict_cart <- function(model, newdata, type = "prob") {
  # robustly restore method
  if (is.null(model$method))
    model$method <- ifelse(is.factor(model$y), "class", "anova")
  
  if (model$method == "class") {
    if (type %in% c("prob","matrix")) {
      out <- tryCatch({
        predict(model, newdata, type = "prob")
      }, error = function(e) {
        tryCatch(predict(model, newdata, type = "matrix"),
                 error = function(e2) predict(model, newdata))
      })
      return(as.data.frame(out))
    }
    if (type == "class") return(predict(model, newdata, type = "class"))
  } else {
    if (type %in% c("vector","prob","matrix")) {
      out <- tryCatch({
        predict(model, newdata, type = "vector")
      }, error = function(e) predict(model, newdata))
      return(as.numeric(out))
    }
  }
  stop("Unsupported type for predict_cart().")
}

## ============================================================
## Correct node_ids(): compute 'where' for *given* data
## ============================================================
node_ids <- function(model, df) {
  # remove response if present
  yname <- all.vars(formula(model))[1]
  if (yname %in% names(df)) df <- df[setdiff(names(df), yname)]
  # align columns/order to training terms
  tr_terms <- delete.response(terms(model))
  tr_vars  <- all.vars(tr_terms)
  missing_vars <- setdiff(tr_vars, names(df))
  for (v in missing_vars) df[[v]] <- NA
  df <- df[, tr_vars, drop = FALSE]
  # predict where for df
  out <- tryCatch(
    predict(model, df, type = "where"),
    error = function(e) rep(NA_integer_, nrow(df))
  )
  as.character(out)
}

## ============================================================
## Bagging schemes: Classical OOB & SB-OOB
## ============================================================

bag_cart_OOB <- function(df, yname, n_boot = 100) {
  n <- nrow(df)
  y <- df[[yname]]
  task <- if (is_classification(y)) "class" else "anova"
  models <- vector("list", n_boot)
  inbag <- matrix(FALSE, nrow = n, ncol = n_boot)
  distinct_counts <- integer(n_boot)   # v2 追踪 U_b
  total_draws     <- integer(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- bootstrap_indices(n)
    uniq <- unique(idx)
    inbag[uniq, b] <- TRUE
    distinct_counts[b] <- length(uniq)
    total_draws[b]     <- length(idx)
    models[[b]] <- fit_cart(df[idx, , drop = FALSE], yname, method = task)
  }
  list(models = models, inbag = inbag, yname = yname, task = task,
       distinct_counts = distinct_counts, total_draws = total_draws,
       scheme = "classical")
}

bag_cart_SB <- function(df, yname, n_boot = 100, k_distinct = NULL) {
  n <- nrow(df)
  if (is.null(k_distinct)) k_distinct <- ceiling(0.632 * n)
  y <- df[[yname]]
  task <- if (is_classification(y)) "class" else "anova"
  models <- vector("list", n_boot)
  inbag <- matrix(FALSE, nrow = n, ncol = n_boot)
  distinct_counts <- integer(n_boot)
  total_draws     <- integer(n_boot)
  for (b in seq_len(n_boot)) {
    U <- integer(0); idx <- integer(0)
    while (length(U) < k_distinct) {
      j <- sample.int(n, 1, replace = TRUE)
      idx <- c(idx, j)
      if (!(j %in% U)) U <- c(U, j)
    }
    inbag[U, b] <- TRUE
    distinct_counts[b] <- length(U)
    total_draws[b]     <- length(idx)
    models[[b]] <- fit_cart(df[idx, , drop = FALSE], yname, method = task)
  }
  list(models = models, inbag = inbag, yname = yname, task = task,
       distinct_counts = distinct_counts, total_draws = total_draws,
       k_distinct = k_distinct, scheme = "sequential")
}

## ============================================================
## OOB predictions (train) & bag predictions (test)
## ============================================================

oob_predict_train <- function(bag, df) {
  n <- nrow(df)
  y <- df[[bag$yname]]
  if (is.null(bag$task)) bag$task <- if (is_classification(y)) "class" else "anova"
  
  if (bag$task == "class") {
    levs <- levels(as.factor(y))
    prob_mat <- matrix(NA_real_, nrow = n, ncol = length(levs), dimnames = list(NULL, levs))
    for (i in seq_len(n)) {
      keep <- !bag$inbag[i, ]
      if (!any(keep)) next
      acc <- rep(0, length(levs)); cnt <- 0L
      for (m in bag$models[keep]) {
        pr <- tryCatch(predict_cart(m, df[i, , drop = FALSE], type = "prob"), error = function(e) NULL)
        if (!is.null(pr)) {
          pr <- as.numeric(pr[1, levs, drop = TRUE])
          acc <- acc + pr; cnt <- cnt + 1L
        }
      }
      if (cnt > 0) prob_mat[i, ] <- acc / cnt
    }
    list(type = "class", probs = prob_mat, levels = levs)
  } else {
    preds <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      keep <- !bag$inbag[i, ]
      if (!any(keep)) next
      vals <- vapply(bag$models[keep],
                     function(m) tryCatch(predict_cart(m, df[i, , drop = FALSE], type = "vector"),
                                          error = function(e) NA_real_),
                     numeric(1))
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0) preds[i] <- mean(vals)
    }
    list(type = "reg", preds = preds)
  }
}

bag_predict <- function(bag, newdata, type = c("prob","class","vector")) {
  type <- match.arg(type)
  if (bag$task == "class") {
    levs <- NULL; sum_prob <- NULL; n_models <- 0L
    for (m in bag$models) {
      pr <- tryCatch(predict_cart(m, newdata, type = "prob"), error = function(e) NULL)
      if (is.null(pr)) next
      if (is.null(levs)) { levs <- colnames(pr); sum_prob <- as.matrix(pr) }
      else { sum_prob <- sum_prob + as.matrix(pr[, levs, drop = FALSE]) }
      n_models <- n_models + 1L
    }
    prob_avg <- sum_prob / n_models
    if (type == "prob")  return(as.data.frame(prob_avg))
    if (type == "class") {
      cls <- levs[max.col(prob_avg, ties.method = "first")]
      return(factor(cls, levels = levs))
    }
  } else {
    mats <- lapply(bag$models, function(m) {
      tryCatch(predict_cart(m, newdata, type = "vector"),
               error = function(e) rep(NA_real_, nrow(newdata)))
    })
    preds <- rowMeans(do.call(cbind, mats), na.rm = TRUE)
    if (type == "vector") return(as.numeric(preds))
  }
  stop("Unsupported type for bag_predict().")
}

## ============================================================
## Metrics
## ============================================================

E1E2_node_prob <- function(p_star_list, p_hat_list, q_star) {
  nodes <- intersect(names(p_star_list), names(p_hat_list))
  if (!length(nodes)) return(list(E1 = NA_real_, E2 = NA_real_))
  J <- length(p_star_list[[nodes[1]]])
  E1 <- 0; E2 <- 0
  for (nd in nodes) {
    ps <- p_star_list[[nd]]; ph <- p_hat_list[[nd]]
    E1 <- E1 + q_star[[nd]] * sum(abs(ps - ph)) / J
    E2 <- E2 + q_star[[nd]] * sqrt(sum((ps - ph)^2) / J)
  }
  list(E1 = E1, E2 = E2)
}

E1E2_node_err <- function(e_star, e_hat, q_star) {
  nodes <- intersect(names(e_star), names(e_hat))
  if (!length(nodes)) return(list(E1 = NA_real_, E2 = NA_real_))
  E1 <- sum(q_star[nodes] * abs(e_star[nodes] - e_hat[nodes]))
  E2 <- sqrt(sum(q_star[nodes] * (e_star[nodes] - e_hat[nodes])^2))
  list(E1 = E1, E2 = E2)
}

test_error_class <- function(y_true, y_pred) mean(y_true != y_pred)
test_mse        <- function(y_true, y_pred) mean((y_true - y_pred)^2)

lower_bound_class <- function(p_hat, N_train, N_test) {
  sqrt(2 * (p_hat * (1 - p_hat) * (1 / N_train + 1 / N_test)) / pi)
}
lower_bound_reg <- function(resid2_test, N_train, N_test) {
  cvar <- var(resid2_test)
  sqrt(2 * (cvar * (1 / N_train + 1 / N_test)) / pi)
}

## ============================================================
## Data generators / loaders
## ============================================================

# Synthetic (used indirectly; some wrappers are handy)
gen_waveform  <- function(n) mlbench.waveform(n)
gen_twonorm   <- function(n) mlbench.twonorm(n)
gen_threenorm <- function(n) mlbench.threenorm(n)
gen_ringnorm  <- function(n) mlbench.ringnorm(n)

gen_friedman1 <- function(n, sd = 1) {
  d <- mlbench.friedman1(n, sd = sd)
  as_tibble(cbind(as.data.frame(d$x), y = d$y))
}
gen_friedman2 <- function(n, sd = 125) {
  d <- mlbench.friedman2(n, sd = sd)
  as_tibble(cbind(as.data.frame(d$x), y = d$y))
}
gen_friedman3 <- function(n, sd = 0.1) {
  d <- mlbench.friedman3(n, sd = sd)
  as_tibble(cbind(as.data.frame(d$x), y = d$y))
}

# Real data loaders
load_breast <- function() {
  data("BreastCancer", package = "mlbench")
  as_tibble(BreastCancer) %>%
    mutate(across(everything(), ~replace(., . %in% "?", NA))) %>%
    drop_na() %>%
    mutate(Class = factor(Class)) %>%
    select(-Id) %>%
    mutate(across(-Class, as.numeric))
}

load_pima <- function() {
  data("PimaIndiansDiabetes", package = "mlbench")
  as_tibble(PimaIndiansDiabetes) %>% mutate(diabetes = factor(diabetes))
}

load_vehicle <- function() {
  data("Vehicle", package = "mlbench")
  as_tibble(Vehicle) %>% mutate(Class = factor(Class))
}

load_satellite <- function() {
  data("Satellite", package = "mlbench")
  as_tibble(Satellite) %>% mutate(classes = factor(classes))
}

load_dna <- function() {
  data("DNA", package="mlbench")
  # Two possible shapes depending on mlbench version
  if (exists("DNA") && is.list(DNA) && all(c("x","classes") %in% names(DNA))) {
    df <- as_tibble(DNA$x); df$class <- factor(DNA$classes); return(df)
  }
  if (exists("DNA") && is.data.frame(DNA)) {
    df <- as_tibble(DNA)
    if (!"class" %in% names(df)) names(df)[ncol(df)] <- "class"
    df$class <- factor(df$class)
    # 保持 predictors 原始类型（字符/因子），不要强制数值化
    return(df)
  }
  stop("Unrecognized DNA dataset structure.")
}

load_boston <- function() MASS::Boston %>% as_tibble() %>% rename(y = medv)

load_ozone <- function() {
  data("Ozone", package = "mlbench")
  df <- as_tibble(Ozone); names(df)[ncol(df)] <- "y"; df
}