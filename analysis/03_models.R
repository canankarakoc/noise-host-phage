# 03_models.R — noise x dispersal models for every response variable
#
# Two decisions made explicit here:
#  (1) CLONE vs MICROCOSM. Growth, phage sensitivity and plate MIC are measured
#      on ~8 clones per microcosm. Those clones are subsamples, not independent
#      replicates, so every clone-level model carries (1 | microcosmID). The
#      experimental unit is the microcosm (n = 36).
#  (2) DAY. Four microcosms (m15-m18, all white noise + isolation) were
#      contaminated at day 30 and are represented by their day-24 clones.
#      PRIMARY  = day 30, with day-24 clones standing in for those four.
#      POOLED   = days 24 and 30 together, with day as a fixed effect.
#      Day 18 is excluded throughout.
#
# Besides the console report, writes to OUT_DIR:
#   model_anova.csv      type III F tests (noise, dispersal, interaction)
#   model_emmeans.csv    marginal means by noise within dispersal, model scale
#                        (log2 for MIC, natural log for generation time)
#   model_contrasts.csv  Tukey contrasts within dispersal + dispersal main effect
#   paired_red_white.csv patch-paired red and white values (D below)
#   microcosm_cfu_variability.csv, microcosm_host_phage_corr.csv
# 05_figures.R reads these, so the figures show the same estimates as the text.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest)
  library(emmeans); library(car)
})
emm_options(lmerTest.limit = 10000, pbkrtest.limit = 10000)

# Type III tests need sum-to-zero contrasts. Under R's default treatment
# contrasts, car::Anova(type = 3) tests each "main effect" at the reference
# level of the other factor (noise within isolation, dispersal within the
# constant regime), not averaged over it. Marginal means and contrasts from
# emmeans are the same under either coding.
options(contrasts = c("contr.sum", "contr.poly"))

OUT_DIR <- Sys.getenv("OUT_DIR", "results")
prep <- readRDS(file.path(OUT_DIR, "prep.rds"))
CONTAM <- prep$CONTAM
NOISE_LEV <- c("constant", "white noise", "red noise")

fx <- function(d) {
  d <- d %>% mutate(
    noise = factor(noise, levels = NOISE_LEV),
    treatment = factor(treatment, levels = c("isolation","dispersal")),
    patch = factor(patch), microcosmID = factor(microcosmID))
  if ("day" %in% names(d)) d$day <- factor(d$day)
  d
}

growth <- fx(read_csv(file.path(OUT_DIR,"trait_growth.csv"), show_col_types=FALSE))
auc    <- fx(read_csv(file.path(OUT_DIR,"trait_phage_sensitivity.csv"), show_col_types=FALSE))
micp   <- fx(read_csv(file.path(OUT_DIR,"trait_mic_plate.csv"), show_col_types=FALSE))
# MIC_use == 0 means the well was already below the negative control at zero
# kanamycin, i.e. the clone did not grow at all. That is a failed well, not a
# minimum inhibitory concentration of zero.
cat("MIC plate wells with MIC_use == 0 (no growth even without drug), dropped:",
    sum(micp$MIC_use == 0), "of", nrow(micp), "\n")
print(micp %>% filter(MIC_use == 0) %>% count(noise, treatment) %>% as.data.frame())
micp <- micp %>% filter(MIC_use > 0)
mics   <- fx(read_csv(file.path(OUT_DIR,"trait_mic_stripe.csv"), show_col_types=FALSE))
sel    <- fx(read_csv(file.path(OUT_DIR,"trait_selection.csv"), show_col_types=FALSE))
dyn    <-    read_csv(file.path(OUT_DIR,"dynamics_long.csv"), show_col_types=FALSE)

# drop day 18 and build the two day-handling variants
prep_days <- function(d) {
  d <- d %>% filter(day %in% c("24","30"))
  primary <- bind_rows(
    d %>% filter(day == "30", !microcosmID %in% CONTAM),
    d %>% filter(day == "24",  microcosmID %in% CONTAM))
  list(primary = primary, pooled = d)
}

cat("##############################################################\n")
cat("# A. IS DAY CONFOUNDED WITH THE WHITE-NOISE / ISOLATION CELL?\n")
cat("##############################################################\n")
for (nm in c("growth","auc","micp")) {
  d <- get(nm) %>% filter(day %in% c("24","30"))
  tb <- table(interaction(d$noise, d$treatment, sep=" / "), d$day)
  cat("\n--", nm, ": rows by cell x day --\n"); print(tb)
}
cat("\nIn PRIMARY, every cell is day 30 except white noise/isolation, which is\n",
    "2 microcosms at day 30 and 4 at day 24. Any effect in that one cell is\n",
    "therefore partly a day effect.\n", sep="")

# ---------------------------------------------------------------
# generic reporting helper
# ---------------------------------------------------------------
TABLES <- list(anova = list(), emmeans = list(), contrasts = list(), pairs = list())

tidy_anova <- function(aov, key) {
  a <- as.data.frame(aov)
  fcol <- intersect(c("F", "F value"), names(a))[1]
  df_res <- if ("Df.res" %in% names(a)) a$Df.res else rep(a["Residuals", "Df"], nrow(a))
  tibble(response = key, term = rownames(a), F = a[[fcol]], df = a$Df,
         df_res = df_res, p = a[["Pr(>F)"]]) %>%
    filter(term %in% c("noise", "treatment", "noise:treatment"))
}

report <- function(model, label, data, unit = c("clone","microcosm"), key = NULL) {
  unit <- match.arg(unit)
  cat("\n==============================================================\n")
  cat("##", label, "\n")
  cat("   n =", nrow(data), unit, "observations across",
      n_distinct(data$microcosmID), "microcosms\n")
  cat("--------------------------------------------------------------\n")
  aov <- try(car::Anova(model, type = 3, test.statistic = "F"), silent = TRUE)
  if (inherits(aov, "try-error")) aov <- car::Anova(model, type = 3)
  print(aov)
  # partial eta^2 style effect size from F and df where available
  emm <- emmeans(model, ~ noise | treatment)
  ct  <- summary(pairs(emm), infer = c(TRUE, TRUE))
  dsp <- summary(pairs(emmeans(model, ~ treatment)), infer = c(TRUE, TRUE))
  cat("\nestimated marginal means:\n"); print(summary(emm))
  cat("\npairwise contrasts within each dispersal treatment (Tukey):\n")
  print(ct)
  cat("\nmain effect of dispersal:\n")
  print(dsp)
  if (!is.null(key)) {
    TABLES$anova[[key]] <<- tidy_anova(aov, key)
    TABLES$emmeans[[key]] <<- as_tibble(as.data.frame(summary(emm))) %>%
      mutate(response = key, .before = 1)
    TABLES$contrasts[[key]] <<- bind_rows(
      as_tibble(as.data.frame(ct)) %>% mutate(treatment = as.character(treatment)),
      as_tibble(as.data.frame(dsp)) %>% mutate(treatment = "averaged over noise")) %>%
      mutate(response = key, .before = 1)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------
# B. CLONE-LEVEL TRAITS
# ---------------------------------------------------------------
cat("\n\n##############################################################\n")
cat("# B. CLONE-LEVEL TRAITS  (random intercept for microcosm)\n")
cat("##############################################################\n")

# --- B1 generation time (log scale; right-skewed) --------------
g <- prep_days(growth)
m_gen <- lmer(log(generation_el) ~ noise*treatment + (1|microcosmID), data = g$primary)
report(m_gen, "B1. Generation time (model-free estimator, log min) - PRIMARY", g$primary,
       key = "generation_time")

cat("\n>> ROBUSTNESS: same model with the two parametric estimators\n")
for (v in c("generation_log","generation_gom")) {
  mm <- lmer(as.formula(paste0("log(",v,") ~ noise*treatment + (1|microcosmID)")), data = g$primary)
  cat("\n--", v, "--\n"); print(car::Anova(mm, type=3, test.statistic="F"))
  print(summary(emmeans(mm, ~ noise | treatment), type="response"))
}

# --- B2 carrying capacity --------------------------------------
m_K <- lmer(K_obs ~ noise*treatment + (1|microcosmID), data = g$primary)
report(m_K, "B2. Carrying capacity (max blanked OD) - PRIMARY", g$primary,
       key = "carrying_capacity")

# --- B3 phage sensitivity --------------------------------------
a <- prep_days(auc)
m_auc <- lmer(phage_sensitivity ~ noise*treatment + (1|microcosmID), data = a$primary)
report(m_auc, "B3. Phage sensitivity (dAUC; HIGHER = MORE sensitive) - PRIMARY", a$primary,
       key = "phage_sensitivity")

# --- B4 MIC plate assay (log2) ---------------------------------
mp <- prep_days(micp)
m_mic <- lmer(log2(MIC_use) ~ noise*treatment + (1|microcosmID), data = mp$primary)
report(m_mic, "B4. Kanamycin MIC, plate assay (log2 ug/ml) - PRIMARY", mp$primary,
       key = "mic_plate")
cat("\nNOTE: MIC is right-censored at 256 ug/ml for",
    sum(mp$primary$censored), "of", nrow(mp$primary), "wells.\n")

# ---------------------------------------------------------------
# C. MICROCOSM-LEVEL TRAITS
# ---------------------------------------------------------------
cat("\n\n##############################################################\n")
cat("# C. MICROCOSM-LEVEL TRAITS  (n = 36, one value per microcosm)\n")
cat("##############################################################\n")

s <- prep_days(sel)
m_sel <- lm(selection_coef ~ noise*treatment, data = s$primary)
report(m_sel, "C1. Selection coefficient vs ancestor - PRIMARY", s$primary,
       key = "selection_coef")

ms <- prep_days(mics)
m_str <- lm(log2(meanMIC) ~ noise*treatment, data = ms$primary)
report(m_str, "C2. Kanamycin MIC, gradient strips (log2) - PRIMARY", ms$primary,
       key = "mic_strip")

# --- C3 variability of bacterial density ------------------------
cv <- dyn %>% filter(cfu_measured) %>%
  group_by(microcosmID, sampleID, patch, noise, treatment) %>%
  summarise(n_obs = n(), mean_logCFU = mean(log10(CFU_ml)),
            sd_logCFU = sd(log10(CFU_ml)),
            cv_logCFU = sd_logCFU/mean_logCFU, .groups="drop") %>% fx()
m_cv <- lm(cv_logCFU ~ noise*treatment, data = cv)
report(m_cv, "C3. Variability of bacterial density (CV of log10 CFU/ml)", cv,
       key = "cv_logCFU")
cat("\nNOTE: CFU observations per microcosm range",
    min(cv$n_obs), "-", max(cv$n_obs), "; CV of a log-scaled quantity is not\n",
    "scale-free, so sd of log10 CFU is reported alongside it.\n")

# the sd of log10 CFU is the variability measure used in the manuscript
m_sd <- lm(sd_logCFU ~ noise*treatment, data = cv)
report(m_sd, "C3b. Variability of bacterial density (sd of log10 CFU/ml)", cv,
       key = "sd_logCFU")

# --- C4 host-phage cross-correlation ----------------------------
cc <- dyn %>% filter(cfu_measured, pfu_measured) %>%
  group_by(microcosmID, sampleID, patch, noise, treatment) %>%
  filter(n() >= 5) %>%
  summarise(n_paired = n(),
            corr_at_0 = cor(log10(PFU_ml), log10(CFU_ml)), .groups="drop") %>% fx()
m_cc <- lm(corr_at_0 ~ noise*treatment, data = cc)
report(m_cc, "C4. Host-phage cross-correlation at lag 0 (paired days only)", cc,
       key = "host_phage_corr")
cat("\nNOTE: computed only on days where BOTH counts exist (>=5 per microcosm),\n",
    "so no interpolated values enter the correlation. Paired days per microcosm: ",
    min(cc$n_paired), "-", max(cc$n_paired), ".\n", sep="")

# ---------------------------------------------------------------
# D. RED vs WHITE PAIRED BY PATCH
# ---------------------------------------------------------------
cat("\n\n##############################################################\n")
cat("# D. RED vs WHITE, PAIRED BY PATCH\n")
cat("# Red patch k is a reordering of the same 30 values as white patch k,\n")
cat("# so patch pairs them on mean, variance and value set exactly.\n")
cat("##############################################################\n")

paired_test <- function(d, yvar, label, key) {
  dd <- d %>% filter(noise %in% c("white noise","red noise")) %>%
    group_by(patch, treatment, noise) %>%
    summarise(y = mean(.data[[yvar]], na.rm=TRUE), .groups="drop") %>%
    pivot_wider(names_from = noise, values_from = y) %>%
    mutate(diff = `red noise` - `white noise`) %>% filter(is.finite(diff))
  TABLES$pairs[[key]] <<- dd %>%
    transmute(response = key, patch, treatment, white = `white noise`,
              red = `red noise`, diff)
  cat("\n--", label, "--\n")
  if (nrow(dd) < 3) { cat("  too few pairs\n"); return(invisible()) }
  tt <- t.test(dd$diff)
  cat(sprintf("  %d patch x dispersal pairs; mean(red - white) = %.3f  95%% CI [%.3f, %.3f]  t = %.2f, p = %.4f\n",
              nrow(dd), mean(dd$diff), tt$conf.int[1], tt$conf.int[2], tt$statistic, tt$p.value))
  print(as.data.frame(dd), row.names = FALSE)
}

paired_test(g$primary,  "generation_el",     "generation time (min)",    "generation_time")
paired_test(g$primary,  "K_obs",             "carrying capacity",        "carrying_capacity")
paired_test(a$primary,  "phage_sensitivity", "phage sensitivity (dAUC)", "phage_sensitivity")
paired_test(mp$primary %>% mutate(l2 = log2(MIC_use)), "l2", "MIC plate (log2)",  "mic_plate")
paired_test(ms$primary %>% mutate(l2 = log2(meanMIC)), "l2", "MIC strips (log2)", "mic_strip")
paired_test(s$primary,  "selection_coef",    "selection coefficient",    "selection_coef")
paired_test(cv,         "cv_logCFU",         "CV of log10 CFU",          "cv_logCFU")
paired_test(cv,         "sd_logCFU",         "sd of log10 CFU",          "sd_logCFU")

# ---------------------------------------------------------------
# E. POOLED-DAY SENSITIVITY CHECK
# ---------------------------------------------------------------
cat("\n\n##############################################################\n")
cat("# E. DOES POOLING DAYS 24 + 30 CHANGE THE ANSWER?\n")
cat("##############################################################\n")
pooled_check <- function(d, form, label) {
  m <- lmer(as.formula(form), data = d)
  cat("\n--", label, "--\n"); print(car::Anova(m, type=3, test.statistic="F"))
}
pooled_check(g$pooled,  "log(generation_el) ~ noise*treatment + day + (1|microcosmID)", "generation time")
pooled_check(g$pooled,  "K_obs ~ noise*treatment + day + (1|microcosmID)",              "carrying capacity")
pooled_check(a$pooled,  "phage_sensitivity ~ noise*treatment + day + (1|microcosmID)",  "phage sensitivity")
pooled_check(mp$pooled, "log2(MIC_use) ~ noise*treatment + day + (1|microcosmID)",      "MIC plate")
cat("\n-- selection coefficient --\n")
print(car::Anova(lmer(selection_coef ~ noise*treatment + day + (1|microcosmID), data=s$pooled),
                 type=3, test.statistic="F"))

# ---------------------------------------------------------------
# tables for the manuscript and 05_figures.R
# ---------------------------------------------------------------
write_csv(bind_rows(TABLES$anova),     file.path(OUT_DIR, "model_anova.csv"))
write_csv(bind_rows(TABLES$emmeans),   file.path(OUT_DIR, "model_emmeans.csv"))
write_csv(bind_rows(TABLES$contrasts), file.path(OUT_DIR, "model_contrasts.csv"))
write_csv(bind_rows(TABLES$pairs),     file.path(OUT_DIR, "paired_red_white.csv"))
write_csv(cv, file.path(OUT_DIR, "microcosm_cfu_variability.csv"))
write_csv(cc, file.path(OUT_DIR, "microcosm_host_phage_corr.csv"))
cat("\nwrote model tables to", OUT_DIR, "\n")

cat("\n\nDONE\n")
