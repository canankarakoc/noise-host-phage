# 02_growth.R — growth parameters from the blanked OD curves
#
# The original pipeline fitted four model families (Gompertz, Baranyi, logistic,
# Weibull) per clone and kept the best by AIC. That makes the resulting
# parameters non-comparable between clones, because mumax means something
# slightly different in each family, and `lag` was defined one way for Baranyi
# (the h0 parameter) and another way for the rest (an operational threshold).
#
# Here we use one definition for every curve:
#   PRIMARY   model-free rolling log-linear regression ("easy linear" method):
#             mumax = max slope of log(OD) over a sliding window with good fit
#             lag   = where that tangent crosses the initial OD
#             K     = maximum blanked OD reached
#   SECONDARY logistic and Gompertz fits, to confirm the ranking is not an
#             artefact of the estimator.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr)})

OUT_DIR <- Sys.getenv("OUT_DIR", "results")
prep <- readRDS(file.path(OUT_DIR, "prep.rds"))
clones <- prep$clones; anc <- prep$anc

WIN    <- 5      # window length in points (15 min apart) = 60 min
MINR2  <- 0.95   # minimum R^2 for a window to count

# ---- model-free estimator -------------------------------------------------
easylinear <- function(time, y) {
  ok <- is.finite(time) & is.finite(y) & y > 0
  time <- time[ok]; y <- y[ok]
  if (length(y) < WIN + 2) return(c(mumax=NA, lag=NA, K=NA, r2=NA, y0=NA))
  ly <- log(y)
  n  <- length(ly)
  best <- list(sl = -Inf, int = NA, r2 = NA)
  for (i in 1:(n - WIN + 1)) {
    idx <- i:(i + WIN - 1)
    fit <- stats::lm.fit(cbind(1, time[idx]), ly[idx])
    sl  <- fit$coefficients[2]
    if (!is.finite(sl) || sl <= 0) next
    ss_res <- sum(fit$residuals^2); ss_tot <- sum((ly[idx] - mean(ly[idx]))^2)
    r2 <- if (ss_tot > 0) 1 - ss_res/ss_tot else 0
    if (r2 >= MINR2 && sl > best$sl) best <- list(sl = sl, int = fit$coefficients[1], r2 = r2)
  }
  if (!is.finite(best$sl) || best$sl <= 0) return(c(mumax=NA, lag=NA, K=NA, r2=NA, y0=NA))
  y0  <- stats::median(y[1:min(3, n)])
  lag <- (log(y0) - best$int) / best$sl     # tangent crosses the starting OD
  c(mumax = unname(best$sl), lag = unname(lag), K = max(y), r2 = unname(best$r2), y0 = unname(y0))
}

# ---- parametric fits ------------------------------------------------------
fit_param <- function(time, y, model = c("logistic","gompertz")) {
  model <- match.arg(model)
  ok <- is.finite(time) & is.finite(y) & y > 0
  time <- time[ok]; y <- y[ok]
  if (length(y) < 6) return(c(mumax=NA, K=NA, rmse=NA))
  pred <- function(p) {
    y0 <- p[1]; mumax <- p[2]; K <- p[3]
    if (model == "logistic") (K*y0)/(y0 + (K - y0)*exp(-mumax*time))
    else K*exp(log(y0/K)*exp(-mumax*time))
  }
  obj <- function(p) {
    if (any(!is.finite(p)) || p[1] <= 0 || p[2] <= 0 || p[3] <= p[1]) return(1e6)
    yh <- pred(p); if (any(!is.finite(yh)) || any(yh <= 0)) return(1e6)
    sum((log(y) - log(yh))^2)                       # fit on the log scale
  }
  start <- c(max(min(y), 5e-4), 0.01, min(max(y), 0.7))
  res <- try(stats::optim(start, obj, method = "L-BFGS-B",
                          lower = c(1e-4, 5e-4, 0.05), upper = c(0.05, 0.1, 0.9)),
             silent = TRUE)
  if (inherits(res, "try-error") || res$value >= 1e6) return(c(mumax=NA, K=NA, rmse=NA))
  c(mumax = res$par[2], K = res$par[3], rmse = sqrt(res$value/length(y)))
}

estimate_all <- function(df) {
  df <- df[order(df$time), ]
  el <- easylinear(df$time, df$OD_cor)
  lo <- fit_param(df$time, df$OD_cor, "logistic")
  go <- fit_param(df$time, df$OD_cor, "gompertz")
  tibble(mumax_el = el["mumax"], lag_el = el["lag"], K_obs = el["K"], r2_el = el["r2"],
         mumax_log = lo["mumax"], K_log = lo["K"], rmse_log = lo["rmse"],
         mumax_gom = go["mumax"], K_gom = go["K"], rmse_gom = go["rmse"])
}

# ---- evolved clones, no-phage curves --------------------------------------
cl <- clones %>% filter(phage_treatment == "without_phage") %>% drop_na(OD_cor)

growth <- cl %>%
  group_by(microcosmID, sampleID, patch, noise, treatment, clones, day) %>%
  group_modify(~ estimate_all(.x)) %>%
  ungroup() %>%
  mutate(generation_el  = log(2)/mumax_el,
         generation_log = log(2)/mumax_log,
         generation_gom = log(2)/mumax_gom)

cat("=== GROWTH PARAMETERS ===\n")
cat("curves attempted:", nrow(growth), "\n")
cat("easylinear succeeded:", sum(!is.na(growth$mumax_el)),
    "| logistic:", sum(!is.na(growth$mumax_log)),
    "| gompertz:", sum(!is.na(growth$mumax_gom)), "\n")

# keep curves that actually grew and were well described
growth_f <- growth %>% filter(!is.na(mumax_el), K_obs > 0.2, r2_el >= MINR2)
cat("after filter (grew, K>0.2, R2>=", MINR2, "):", nrow(growth_f), "\n", sep="")

cat("\nagreement between estimators (Spearman on generation time):\n")
cat("  easylinear vs logistic:",
    round(cor(growth_f$generation_el, growth_f$generation_log, use="complete.obs", method="spearman"),3), "\n")
cat("  easylinear vs gompertz:",
    round(cor(growth_f$generation_el, growth_f$generation_gom, use="complete.obs", method="spearman"),3), "\n")

cat("\ngeneration time (min) by cell, model-free estimator:\n")
print(growth_f %>% group_by(noise, treatment) %>%
  summarise(n=n(), median=round(median(generation_el),1),
            mean=round(mean(generation_el),1), sd=round(sd(generation_el),1), .groups="drop") %>%
  as.data.frame(), row.names=FALSE)

cat("\ncarrying capacity (max blanked OD) by cell:\n")
print(growth_f %>% group_by(noise, treatment) %>%
  summarise(n=n(), median=round(median(K_obs),3), mean=round(mean(K_obs),3),
            sd=round(sd(K_obs),3), .groups="drop") %>% as.data.frame(), row.names=FALSE)

# ---- ancestral baselines, both strains ------------------------------------
anc_g <- anc %>% filter(phage_treatment == "without_phage") %>% drop_na(OD_cor) %>%
  group_by(strain, set, clones) %>% group_modify(~ estimate_all(.x)) %>% ungroup() %>%
  mutate(generation_el = log(2)/mumax_el) %>%
  filter(!is.na(mumax_el), K_obs > 0.2, r2_el >= MINR2)

cat("\nancestral baseline by strain (REL607 is the evolved lineages' ancestor):\n")
print(anc_g %>% group_by(strain) %>%
  summarise(n=n(), gen_median=round(median(generation_el),1),
            gen_min=round(min(generation_el),1), gen_max=round(max(generation_el),1),
            K_median=round(median(K_obs),3), .groups="drop") %>% as.data.frame(), row.names=FALSE)

write_csv(growth_f, file.path(OUT_DIR, "trait_growth.csv"))
write_csv(anc_g,    file.path(OUT_DIR, "ancestral_growth.csv"))
cat("\nwrote growth traits\n")
