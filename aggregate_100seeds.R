## ============================================================
## aggregate_100seeds.R
##
## 跑完 100 个 seed (results/seed_001 ... seed_100) 之后, 这个脚本:
##   1. 把所有 seed 的 exp1..exp5 summary 拼成一个长表
##   2. 跨 seed 做配对统计 (Wilcoxon, t, Cohen's d_z, Pitman-Morgan)
##   3. 输出 paper Section 5 用到的 7 张 CSV
##
## 用法:
##   cd ~/My_Project
##   Rscript --vanilla aggregate_100seeds.R
##
## 输入: results/seed_001/exp{1..5}.RData ... results/seed_100/exp{1..5}.RData
##       (raw runs *_raw.RData 不强制需要; 主统计基于 summary RData)
##
## 输出 (写到 tables/):
##   all_seeds_flat.csv          摊平的长表 (6200 行)
##   agg_100seeds.csv            跨 seed 配对统计主表 (57 cells × 18 列)
##   table_exp{1..5}_100seeds.csv  按 EXP 拆分的 paper-ready 表
##   summary_by_exp.csv          每个 EXP 的总览
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
})

RESULTS_DIR <- "results"
OUT_DIR     <- "tables"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

## ---- 1) 把所有 seed 的 5 个 EXP summary 拼成长表 ----
seed_dirs <- list.files(RESULTS_DIR, pattern = "^seed_\\d{3}$", full.names = TRUE)
cat(sprintf("[INFO] Found %d seed dirs in %s\n", length(seed_dirs), RESULTS_DIR))
stopifnot(length(seed_dirs) > 0)

load_one_exp <- function(seed_dir, exp_name) {
  fp <- file.path(seed_dir, paste0(exp_name, ".RData"))
  if (!file.exists(fp)) return(NULL)
  env <- new.env()
  load(fp, envir = env)
  if (!"res" %in% ls(env)) return(NULL)
  df <- env$res
  df$exp  <- exp_name
  df$seed <- as.integer(sub(".*seed_(\\d+).*", "\\1", seed_dir))
  df
}

all_long <- map_dfr(seed_dirs, function(sd) {
  map_dfr(c("exp1","exp2","exp3","exp4","exp5"), ~ load_one_exp(sd, .x))
})

cat(sprintf("[INFO] Loaded %d rows total across %d seeds\n",
            nrow(all_long), length(unique(all_long$seed))))

write_csv(all_long, file.path(OUT_DIR, "all_seeds_flat.csv"))
cat(sprintf("[INFO] Wrote %s\n", file.path(OUT_DIR, "all_seeds_flat.csv")))

## ---- 2) 跨 seed 配对统计 ----
.cross_seed_stats <- function(oob, sb) {
  d <- sb - oob
  n <- length(d)
  out <- list(
    n_seeds   = n,
    OOB_mean  = mean(oob),
    SB_mean   = mean(sb),
    mean_diff = mean(d),
    sd_diff   = sd(d)
  )
  out$cohens_dz <- if (out$sd_diff > 0) out$mean_diff / out$sd_diff else NA_real_

  # Paired Wilcoxon
  out$wilcox_p <- tryCatch(
    if (any(d != 0))
      wilcox.test(oob, sb, paired = TRUE, exact = FALSE)$p.value
    else 1.0,
    error = function(e) NA_real_)

  # Paired t
  out$t_p <- tryCatch(t.test(oob, sb, paired = TRUE)$p.value,
                      error = function(e) NA_real_)

  # CI on diff
  out$ci_lo <- tryCatch(t.test(d)$conf.int[1], error = function(e) NA_real_)
  out$ci_hi <- tryCatch(t.test(d)$conf.int[2], error = function(e) NA_real_)

  # Sign test
  np <- sum(d > 0); nn <- sum(d < 0)
  out$n_pos  <- np; out$n_neg <- nn; out$n_zero <- sum(d == 0)
  out$sign_p <- if ((np + nn) > 0)
    tryCatch(binom.test(min(np, nn), np + nn, p = 0.5)$p.value,
             error = function(e) NA_real_) else NA_real_

  # Pitman-Morgan (paired variance test)
  out$sd_OOB <- sd(oob); out$sd_SB <- sd(sb)
  out$sd_ratio <- if (out$sd_OOB > 0) out$sd_SB / out$sd_OOB else NA_real_

  out$PM_p_two <- tryCatch({
    if (out$sd_OOB > 1e-15 && out$sd_SB > 1e-15) {
      sp_ <- oob + sb; sm_ <- sb - oob
      r <- cor(sp_, sm_)
      if (abs(r) < 0.9999) {
        tt <- r * sqrt((n - 2) / (1 - r^2))
        2 * (1 - pt(abs(tt), n - 2))
      } else NA_real_
    } else NA_real_
  }, error = function(e) NA_real_)

  tibble::as_tibble(out)
}

agg <- all_long %>%
  group_by(exp, dataset, type, metric) %>%
  summarise(stats = list(.cross_seed_stats(OOB, SB_OOB)),
            .groups = "drop") %>%
  unnest(stats)

write_csv(agg, file.path(OUT_DIR, "agg_100seeds.csv"))
cat(sprintf("[INFO] Wrote %s (%d cells)\n",
            file.path(OUT_DIR, "agg_100seeds.csv"), nrow(agg)))

## ---- 3) Per-EXP paper-ready tables ----
add_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

agg$sig_mean <- add_stars(agg$wilcox_p)
agg$sig_var  <- add_stars(agg$PM_p_two)

cols_paper <- c("dataset","type","metric",
                "OOB_mean","SB_mean","mean_diff","sd_diff",
                "cohens_dz","n_pos","n_neg","wilcox_p","sig_mean",
                "sd_OOB","sd_SB","sd_ratio","PM_p_two","sig_var")

for (e in c("exp1","exp2","exp3","exp4","exp5")) {
  fp <- file.path(OUT_DIR, sprintf("table_%s_100seeds.csv", e))
  agg %>% filter(exp == e) %>% select(all_of(cols_paper)) %>% write_csv(fp)
  cat(sprintf("[INFO] Wrote %s (%d rows)\n", fp, sum(agg$exp == e)))
}

## ---- 4) Cross-experiment summary ----
agg_clean <- agg %>%
  filter(!grepl("mse_original$", metric))   # mse_original 恒为 0

summary_by_exp <- agg_clean %>%
  group_by(exp) %>%
  summarise(
    n_cells        = n(),
    p_mean_lt05    = sum(wilcox_p < 0.05, na.rm = TRUE),
    SB_smaller_mean= sum(wilcox_p < 0.05 & mean_diff < 0, na.rm = TRUE),
    SB_larger_mean = sum(wilcox_p < 0.05 & mean_diff > 0, na.rm = TRUE),
    p_var_lt05     = sum(PM_p_two < 0.05, na.rm = TRUE),
    sd_ratio_lt1   = sum(PM_p_two < 0.05 & sd_ratio < 1, na.rm = TRUE),
    mean_sd_ratio  = mean(sd_ratio, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_by_exp, file.path(OUT_DIR, "summary_by_exp.csv"))
cat("\n[INFO] Cross-experiment summary:\n")
print(summary_by_exp)

## ---- 5) Real vs Synthetic dichotomy (EXP1, E1_B) ----
ex1 <- agg %>% filter(exp == "exp1", metric == "E1_B")
real_sr <- ex1$sd_ratio[ex1$type == "real"]
syn_sr  <- ex1$sd_ratio[ex1$type == "synthetic"]

cat("\n[INFO] Real vs Synthetic dichotomy (EXP1, E1_B):\n")
cat(sprintf("  real      (n=%d): mean sd_ratio = %.3f\n",
            length(real_sr), mean(real_sr)))
cat(sprintf("  synthetic (n=%d): mean sd_ratio = %.3f\n",
            length(syn_sr),  mean(syn_sr)))

if (length(real_sr) > 0 && length(syn_sr) > 0) {
  mw <- wilcox.test(real_sr, syn_sr, alternative = "less")
  cat(sprintf("  Mann-Whitney one-sided (real<syn): p = %.4f\n", mw$p.value))

  set.seed(0)
  obs_diff <- mean(real_sr) - mean(syn_sr)
  combined <- c(real_sr, syn_sr); n_r <- length(real_sr)
  B <- 100000
  perm_diffs <- replicate(B, {
    p <- sample(combined); mean(p[1:n_r]) - mean(p[-(1:n_r)])
  })
  cat(sprintf("  Permutation (B=%d) one-sided p = %.4f\n",
              B, mean(perm_diffs <= obs_diff)))
}

cat("\n[DONE] Aggregation complete. See tables/ for output.\n")
