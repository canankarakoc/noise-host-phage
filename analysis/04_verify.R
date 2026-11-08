# 04_verify.R — stress-test the results that matter
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest)
  library(emmeans); library(car)
})
emm_options(lmerTest.limit = 10000, pbkrtest.limit = 10000)
# sum-to-zero contrasts so type III F tests are main effects (see 03_models.R)
options(contrasts = c("contr.sum", "contr.poly"))
OUT_DIR <- Sys.getenv("OUT_DIR", "results")
prep <- readRDS(file.path(OUT_DIR,"prep.rds")); CONTAM <- prep$CONTAM
NOISE_LEV <- c("constant","white noise","red noise")
fx <- function(d){ d <- d %>% mutate(noise=factor(noise,levels=NOISE_LEV),
  treatment=factor(treatment,levels=c("isolation","dispersal")),
  patch=factor(patch), microcosmID=factor(microcosmID))
  if("day" %in% names(d)) d$day <- factor(d$day); d }

micp <- fx(read_csv(file.path(OUT_DIR,"trait_mic_plate.csv"),show_col_types=FALSE)) %>% filter(MIC_use>0)
sel  <- fx(read_csv(file.path(OUT_DIR,"trait_selection.csv"),show_col_types=FALSE))
auc  <- fx(read_csv(file.path(OUT_DIR,"trait_phage_sensitivity.csv"),show_col_types=FALSE))
grw  <- fx(read_csv(file.path(OUT_DIR,"trait_growth.csv"),show_col_types=FALSE))
dyn  <- read_csv(file.path(OUT_DIR,"dynamics_long.csv"),show_col_types=FALSE)

primary <- function(d) bind_rows(d %>% filter(day=="30", !microcosmID %in% CONTAM),
                                 d %>% filter(day=="24",  microcosmID %in% CONTAM))
day30only <- function(d) d %>% filter(day=="30", !microcosmID %in% CONTAM)
pooled   <- function(d) d %>% filter(day %in% c("24","30"))

hr <- function(x) cat("\n", strrep("=",64), "\n", x, "\n", strrep("=",64), "\n", sep="")

# ---------------------------------------------------------------
hr("V1. RED vs WHITE as a mixed model, blocked by patch")
cat("Red patch k uses the same 30 concentrations as white patch k, reordered.\n",
    "Blocking on patch therefore removes mean, variance and value set, leaving\n",
    "only temporal structure.\n", sep="")

rw_test <- function(d, resp, lab, clone_level = TRUE) {
  dd <- d %>% filter(noise %in% c("white noise","red noise")) %>%
    mutate(noise = droplevels(factor(noise, levels = c("white noise","red noise"))))
  f <- if (clone_level) paste0(resp," ~ noise*treatment + (1|patch) + (1|microcosmID)")
       else             paste0(resp," ~ noise*treatment + (1|patch)")
  m <- suppressMessages(lmer(as.formula(f), data = dd))
  cat("\n--", lab, "| n =", nrow(dd), "obs,", n_distinct(dd$microcosmID), "microcosms --\n")
  print(car::Anova(m, type=3, test.statistic="F"))
  ct <- summary(pairs(emmeans(m, ~ noise)), infer=c(TRUE,TRUE))
  print(ct)
  invisible(m)
}

cat("\n### PRIMARY (day 30, day-24 stand-ins for the 4 contaminated) ###")
rw_test(primary(micp), "log2(MIC_use)",     "MIC plate (log2 ug/ml)")
rw_test(primary(auc),  "phage_sensitivity", "phage sensitivity (dAUC)")
rw_test(primary(grw),  "log(generation_el)","generation time (log min)")
rw_test(primary(sel),  "selection_coef",    "selection coefficient", clone_level=FALSE)

cat("\n\n### POOLED (days 24 + 30) ###")
rw_test(pooled(micp), "log2(MIC_use)",     "MIC plate (log2 ug/ml)")
rw_test(pooled(sel),  "selection_coef",    "selection coefficient", clone_level=FALSE)

cat("\n\n### DAY 30 ONLY (contaminated microcosms dropped, no stand-ins) ###")
rw_test(day30only(micp), "log2(MIC_use)",  "MIC plate (log2 ug/ml)")
rw_test(day30only(sel),  "selection_coef", "selection coefficient", clone_level=FALSE)

# ---------------------------------------------------------------
hr("V2. IS THE RESULT CARRIED BY THE DAY-24 STAND-INS?")
cat("The 4 contaminated microcosms are all white noise / isolation, so that one\n",
    "cell is mostly day-24 clones. Refit the full model on day 30 only.\n", sep="")

cat("\n-- MIC plate, day 30 only (n =", nrow(day30only(micp)), "clones,",
    n_distinct(day30only(micp)$microcosmID), "microcosms) --\n")
m <- lmer(log2(MIC_use) ~ noise*treatment + (1|microcosmID), data = day30only(micp))
print(car::Anova(m, type=3, test.statistic="F"))
print(summary(emmeans(m, ~ noise | treatment)))

cat("\n-- selection coefficient, day 30 only --\n")
m2 <- lm(selection_coef ~ noise*treatment, data = day30only(sel))
print(car::Anova(m2, type=3))
print(summary(emmeans(m2, ~ noise | treatment)))
print(summary(pairs(emmeans(m2, ~ noise | treatment)), infer=c(TRUE,TRUE)))

# ---------------------------------------------------------------
hr("V3. DISPERSAL-ONLY CHECK (no contamination in these 18 microcosms)")
cat("All 18 dispersal microcosms have clean day-30 data, so this subset is free\n",
    "of the day confound entirely.\n", sep="")
dd <- day30only(micp) %>% filter(treatment == "dispersal")
m3 <- lmer(log2(MIC_use) ~ noise + (1|microcosmID), data = dd)
cat("\n-- MIC plate, dispersal only, day 30 (n =", nrow(dd), ") --\n")
print(car::Anova(m3, type=3, test.statistic="F"))
print(summary(pairs(emmeans(m3, ~ noise)), infer=c(TRUE,TRUE)))

di <- day30only(micp) %>% filter(treatment == "isolation")
m4 <- lmer(log2(MIC_use) ~ noise + (1|microcosmID), data = di)
cat("\n-- MIC plate, isolation only, day 30 (n =", nrow(di), ",",
    n_distinct(di$microcosmID), "microcosms; white noise has only 2 left) --\n")
print(car::Anova(m4, type=3, test.statistic="F"))
print(summary(pairs(emmeans(m4, ~ noise)), infer=c(TRUE,TRUE)))

# ---------------------------------------------------------------
hr("V4. IS CFU VARIABILITY CONFOUNDED WITH SAMPLING EFFORT?")
cv <- dyn %>% filter(cfu_measured) %>%
  group_by(microcosmID, sampleID, patch, noise, treatment) %>%
  summarise(n_obs=n(), mean_logCFU=mean(log10(CFU_ml)), sd_logCFU=sd(log10(CFU_ml)),
            cv_logCFU=sd_logCFU/mean_logCFU, .groups="drop") %>% fx()
cat("CFU observations per microcosm by cell:\n")
print(cv %>% group_by(noise,treatment) %>%
        summarise(mean_n=round(mean(n_obs),1), min_n=min(n_obs), max_n=max(n_obs),
                  .groups="drop") %>% as.data.frame(), row.names=FALSE)
cat("\ncorrelation(n_obs, cv_logCFU) =", round(cor(cv$n_obs, cv$cv_logCFU),3), "\n")
cat("\nsd of log10 CFU (scale-free alternative to CV):\n")
m5 <- lm(sd_logCFU ~ noise*treatment, data = cv)
print(car::Anova(m5, type=3))
print(summary(emmeans(m5, ~ noise | treatment)))
print(summary(pairs(emmeans(m5, ~ noise | treatment)), infer=c(TRUE,TRUE))[,c("contrast","treatment","estimate","p.value")])

# ---------------------------------------------------------------
hr("V5. EFFECT SIZES (variance explained, PRIMARY models)")
r2_lm <- function(m, lab) {
  a <- car::Anova(m, type=2); ss <- a[["Sum Sq"]]; nm <- rownames(a)
  tot <- sum(ss); pe <- ss/tot
  cat("\n--", lab, "-- (proportion of total SS)\n")
  for (i in seq_along(nm)) cat(sprintf("   %-18s %.3f\n", nm[i], pe[i]))
}
r2_lm(lm(selection_coef ~ noise*treatment, data = primary(sel)), "selection coefficient")
r2_lm(lm(cv_logCFU ~ noise*treatment, data = cv),                "CV of log10 CFU")
r2_lm(lm(sd_logCFU ~ noise*treatment, data = cv),                "sd of log10 CFU")

cat("\nFor the clone-level MIC model, variance components:\n")
mm <- lmer(log2(MIC_use) ~ noise*treatment + (1|microcosmID), data = primary(micp))
print(as.data.frame(VarCorr(mm)))
cat("\nproportion of residual+microcosm variance that is between-microcosm:",
    round(as.data.frame(VarCorr(mm))$vcov[1] / sum(as.data.frame(VarCorr(mm))$vcov), 3), "\n")
cat("-> clones within a microcosm are strongly non-independent, which is why\n",
    "   (1|microcosmID) matters.\n", sep="")

cat("\n\nDONE\n")
