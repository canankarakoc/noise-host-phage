# 01_traits.R — derive endpoint traits that need no curve fitting
# Noise-colour experimental evolution (E. coli x T4 x kanamycin)
# Reads data from DATA_DIR, writes tidy trait tables to OUT_DIR

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(pracma)
})

DATA_DIR <- Sys.getenv("DATA_DIR", "data")   # see data/README.md
OUT_DIR  <- Sys.getenv("OUT_DIR",  "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- helper: read a csv and strip whitespace from column names AND character cells
rd <- function(f, ...) {
  x <- read.csv(file.path(DATA_DIR, f), check.names = FALSE,
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", ...)
  names(x) <- trimws(names(x))
  x[] <- lapply(x, function(col) if (is.character(col)) trimws(col) else col)
  x
}

NOISE_LEV <- c("constant", "white noise", "red noise")

# ============================================================
# Design
# ============================================================
design <- rd("design.csv") %>%
  mutate(noise = factor(noise, levels = NOISE_LEV),
         treatment = factor(treatment, levels = c("isolation", "dispersal")),
         patch = factor(patch))

cat("=== DESIGN ===\n")
print(table(design$noise, design$treatment))

# Which microcosms are white noise + isolation? (contamination set)
wn_iso <- design %>% filter(noise == "white noise", treatment == "isolation") %>%
  arrange(patch)
cat("\nwhite noise + isolation microcosms:\n"); print(wn_iso)

# The two contamination definitions used in the original script:
#   growth section: patches p3-p6 of white noise / isolation
#   competition+MIC section: microcosmID m15, m16, m17, m18
contam_by_patch <- wn_iso %>% filter(patch %in% c("p3","p4","p5","p6")) %>% pull(microcosmID)
contam_hardcoded <- c("m15","m16","m17","m18")
cat("\ncontaminated by patch rule :", paste(sort(contam_by_patch), collapse=", "), "\n")
cat("contaminated hardcoded     :", paste(sort(contam_hardcoded), collapse=", "), "\n")
cat("IDENTICAL:", setequal(contam_by_patch, contam_hardcoded), "\n")
CONTAM <- contam_by_patch

# ============================================================
# 1. Competitive fitness (selection coefficient)
# ============================================================
comp <- rd("competition_alldays.csv")

sel <- comp %>%
  select(time, microcosmID, day, dilution_fold, pink_colonies, red_colonies) %>%
  drop_na() %>%
  mutate(sum = pink_colonies + red_colonies) %>%
  filter(sum > 1) %>%                                        # drop singletons
  filter((dilution_fold < 1e5 & time == "t0") |
         (dilution_fold < 1e8 & time == "t1")) %>%            # drop highest dilutions
  mutate(pink_dl = ifelse(pink_colonies == 0, 1, pink_colonies),
         red_dl  = ifelse(red_colonies  == 0, 1, red_colonies),
         red_cfu  = red_dl  * dilution_fold / 10,
         pink_cfu = pink_dl * dilution_fold / 10) %>%
  group_by(time, microcosmID, day) %>%
  summarise(ancestral = mean(red_cfu), evolved = mean(pink_cfu), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = c(ancestral, evolved)) %>%
  mutate(selection_coef = log(evolved_t1/ancestral_t1) - log(evolved_t0/ancestral_t0)) %>%
  left_join(design, by = "microcosmID") %>%
  filter(day != 18, !is.na(noise)) %>%
  select(microcosmID, sampleID, day, patch, noise, treatment, selection_coef)

cat("\n=== SELECTION COEFFICIENT ===\n")
cat("rows:", nrow(sel), " days:", paste(sort(unique(sel$day)), collapse=","), "\n")
print(sel %>% group_by(noise, treatment) %>%
        summarise(n = n(), mean = mean(selection_coef), sd = sd(selection_coef), .groups="drop"))

# ============================================================
# 2. Kanamycin MIC — plate dilution assay
# ============================================================
mic_raw <- rd("mic_od_combined.csv")
mic_cols <- grep("^MIC_", names(mic_raw), value = TRUE)
mic_conc <- as.numeric(sub("^MIC_", "", mic_cols))
cat("\n=== MIC (plate) ===\nconcentrations tested (ug/ml):", paste(mic_conc, collapse=", "), "\n")
MAXCONC <- max(mic_conc)

mic_proc <- mic_raw %>%
  mutate(row_id = row_number()) %>%
  rowwise() %>%
  mutate(
    od_vec       = list(as.numeric(c_across(all_of(mic_cols)))),
    fail_stamp   = stamp_control  > 0.1,
    fail_neg     = `negative control` > 0.1,
    is_monotonic = sum(diff(unlist(od_vec)) > 0.1, na.rm = TRUE) <= 2,
    MIC_idx      = { v <- unlist(od_vec)
                     i <- which(v < (`negative control` + 0.01))[1]; i },
    MIC_est      = if (!is.na(MIC_idx)) mic_conc[MIC_idx] else NA_real_,
    censored     = is.na(MIC_est),                 # never inhibited, even at max conc
    MIC_use      = ifelse(censored, MAXCONC, MIC_est)
  ) %>%
  ungroup()

cat("QC: fail_stamp", sum(mic_proc$fail_stamp), "| fail_neg", sum(mic_proc$fail_neg),
    "| non-monotonic", sum(!mic_proc$is_monotonic), "\n")

mic <- mic_proc %>%
  filter(!fail_stamp, !fail_neg, is_monotonic) %>%
  left_join(design, by = "sampleID") %>%
  filter(day != 18, !is.na(noise)) %>%
  select(day, sampleID, microcosmID, patch, noise, treatment, MIC_use, censored)

cat("rows after QC:", nrow(mic), "| right-censored at", MAXCONC, "ug/ml:",
    sum(mic$censored), sprintf("(%.1f%%)", 100*mean(mic$censored)), "\n")
print(mic %>% group_by(noise, treatment) %>%
        summarise(n = n(), median_MIC = median(MIC_use),
                  pct_censored = round(100*mean(censored),1), .groups="drop"))

# ============================================================
# 3. Kanamycin MIC — gradient test strips (technical reps averaged)
# ============================================================
stripe <- rd("MIC_stripe_test.csv") %>%
  rename(day = Day) %>%
  left_join(design, by = c("microcosmID", "sampleID")) %>%
  filter(day != 18, !is.na(noise)) %>%
  group_by(microcosmID, sampleID, patch, noise, treatment, day) %>%
  summarise(meanMIC = mean(MIC), n_tech = n(), .groups = "drop")

cat("\n=== MIC (stripe) ===\nrows:", nrow(stripe), "\n")
print(stripe %>% group_by(noise, treatment) %>%
        summarise(n = n(), mean = mean(meanMIC), sd = sd(meanMIC), .groups="drop"))

# ============================================================
# 4. Growth curves: blank-correct, assemble, phage sensitivity (AUC)
# ============================================================
plate_design <- rd("plate_design.csv")

gfiles <- sprintf("files/%02d.csv", 1:22)
od_long <- map_dfr(seq_along(gfiles), function(i) {
  x <- rd(gfiles[i])
  x %>% pivot_longer(-time, names_to = "well", values_to = "OD") %>% mutate(set = i)
})
cat("\n=== GROWTH CURVES ===\nOD readings:", nrow(od_long), "\n")

od <- od_long %>% left_join(plate_design, by = c("well", "set"))

blank_means <- od %>% filter(microcosmID == "blank") %>%
  group_by(set, time) %>% summarise(mean_blank = mean(OD, na.rm = TRUE), .groups = "drop")

od_corr <- od %>%
  left_join(blank_means, by = c("set", "time")) %>%
  mutate(OD_cor = pmax(OD - mean_blank, 0.001)) %>%
  filter(time > 0)

# evolved clones
clones <- od_corr %>%
  filter(!microcosmID %in% c("blank", "606", "607")) %>%
  left_join(design, by = c("microcosmID", "sampleID")) %>%
  filter(!is.na(noise)) %>%
  mutate(day = as.integer(day))

# ancestral controls — REL607 is the strain that was evolved; REL606 is the
# competition reference. Keep both so the baseline can be checked.
anc <- od_corr %>% filter(microcosmID %in% c("606", "607")) %>%
  mutate(strain = microcosmID, day = as.integer(day))

cat("clone readings:", nrow(clones), "| ancestral readings:", nrow(anc), "\n")
cat("clone days:", paste(sort(unique(clones$day)), collapse=","), "\n")
cat("ancestral wells per strain x phage:\n")
print(anc %>% distinct(strain, phage_treatment, set, clones) %>% count(strain, phage_treatment))

# --- phage sensitivity = AUC(no phage) - AUC(with phage), per clone
# keep only clones whose no-phage curve actually grew (max OD > 0.2)
grew <- clones %>% filter(phage_treatment == "without_phage") %>%
  group_by(microcosmID, clones, day) %>%
  summarise(max_OD = max(OD_cor), .groups = "drop") %>%
  filter(max_OD > 0.2)

auc <- clones %>%
  semi_join(grew, by = c("microcosmID", "clones", "day")) %>%
  group_by(microcosmID, sampleID, patch, noise, treatment, clones, day, phage_treatment) %>%
  arrange(time, .by_group = TRUE) %>%
  summarise(auc = trapz(time, OD_cor), .groups = "drop") %>%
  pivot_wider(names_from = phage_treatment, values_from = auc) %>%
  filter(!is.na(with_phage), !is.na(without_phage)) %>%
  mutate(phage_sensitivity = without_phage - with_phage)

anc_auc <- anc %>%
  group_by(strain, set, clones, phage_treatment) %>%
  arrange(time, .by_group = TRUE) %>%
  summarise(auc = trapz(time, OD_cor), .groups = "drop") %>%
  pivot_wider(names_from = phage_treatment, values_from = auc) %>%
  filter(!is.na(with_phage), !is.na(without_phage)) %>%
  mutate(phage_sensitivity = without_phage - with_phage)

cat("\n=== PHAGE SENSITIVITY (dAUC) ===\nclone-level rows:", nrow(auc), "\n")
print(auc %>% group_by(noise, treatment) %>%
        summarise(n = n(), mean = mean(phage_sensitivity),
                  sd = sd(phage_sensitivity), .groups="drop"))
cat("\nancestral dAUC by strain:\n")
print(anc_auc %>% group_by(strain) %>%
        summarise(n = n(), mean = mean(phage_sensitivity), sd = sd(phage_sensitivity),
                  min = min(phage_sensitivity), max = max(phage_sensitivity), .groups="drop"))

# ============================================================
# 5. Population dynamics: CFU, PFU, CV, extinction, cross-correlation
# ============================================================
cfu <- rd("dynamicsCFU.csv")
pfu <- rd("dynamicsPFU.csv")

ALL_DAYS <- seq(2, 30, by = 2)
grid <- design %>% distinct(microcosmID, sampleID, patch, noise, treatment) %>%
  tidyr::expand_grid(day = ALL_DAYS)

# CFU per ml: colonies in a 10 ul spot -> x100 for per-ml, x dilution
cfu_d <- cfu %>%
  mutate(CFU_ml = CFU * dilution * 100) %>%
  group_by(sampleID, day) %>%
  summarise(CFU_ml = mean(CFU_ml), n_dil = n(), .groups = "drop")

# PFU per ml: plaques in `volume` ul
pfu_d <- pfu %>%
  mutate(PFU_ml = count * dilution * (1000/volume)) %>%
  group_by(sampleID, day) %>%
  summarise(PFU_ml = mean(PFU_ml), pfu_detected = any(count > 0),
            n_dil = n(), .groups = "drop")

dyn <- grid %>%
  left_join(cfu_d %>% select(sampleID, day, CFU_ml), by = c("sampleID","day")) %>%
  left_join(pfu_d %>% select(sampleID, day, PFU_ml, pfu_detected), by = c("sampleID","day")) %>%
  mutate(
    cfu_measured = !is.na(CFU_ml),
    pfu_measured = !is.na(PFU_ml),
    # phage absent = measured but no plaques found
    phage_absent = pfu_measured & !pfu_detected
  )

cat("\n=== DYNAMICS ===\n")
cat("microcosm-days expected:", nrow(dyn), "\n")
cat("CFU measured:", sum(dyn$cfu_measured), "| PFU measured:", sum(dyn$pfu_measured), "\n")
cat("PFU measured-but-zero (true absence):", sum(dyn$phage_absent), "\n")
cat("PFU not measured (unknown):", sum(!dyn$pfu_measured), "\n")

write_csv(sel,    file.path(OUT_DIR, "trait_selection.csv"))
write_csv(mic,    file.path(OUT_DIR, "trait_mic_plate.csv"))
write_csv(stripe, file.path(OUT_DIR, "trait_mic_stripe.csv"))
write_csv(auc,    file.path(OUT_DIR, "trait_phage_sensitivity.csv"))
write_csv(anc_auc,file.path(OUT_DIR, "ancestral_phage_sensitivity.csv"))
write_csv(dyn,    file.path(OUT_DIR, "dynamics_long.csv"))
saveRDS(list(design = design, CONTAM = CONTAM, clones = clones, anc = anc,
             MAXCONC = MAXCONC),
        file.path(OUT_DIR, "prep.rds"))
cat("\nwrote trait tables to", OUT_DIR, "\n")
