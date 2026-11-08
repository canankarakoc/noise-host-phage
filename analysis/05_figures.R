# 05_figures.R — manuscript and supplementary figures
#
# Reads the tables written by 01-03 (results/) and the kanamycin schedules
# (schedules/), and writes every figure to figures/ as PDF (vector) and PNG
# (600 dpi). Model estimates come from results/model_emmeans.csv, so the
# figures show the same marginal means and 95% CIs as the text.
#
# Trait figures use the primary day rule from 03_models.R: day 30, with day-24
# clones standing in for the four contaminated microcosms (m15-m18). Stand-ins
# are drawn as triangles.
#
#   fig1_schedules         kanamycin schedules and run lengths
#   fig2_dynamics          bacterial and phage densities; variability of density
#   fig3_mic               kanamycin MIC; red vs white paired by patch
#   fig4_selection         competitive fitness against REL606
#   fig5_phage             phage sensitivity against the ancestor REL607
#   figS1_growth           generation time and carrying capacity
#   figS2_mic_strip        MIC from gradient strips
#   figS3_host_phage_corr  host-phage correlation on days with both counts

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

OUT_DIR <- Sys.getenv("OUT_DIR", "results")
FIG_DIR <- Sys.getenv("FIG_DIR", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

prep   <- readRDS(file.path(OUT_DIR, "prep.rds"))
CONTAM <- prep$CONTAM

NOISE_LEV <- c("constant", "white noise", "red noise")
DISP_LEV  <- c("isolation", "dispersal")
NOISE_LAB <- c(constant = "Constant", `white noise` = "White noise", `red noise` = "Red noise")
NOISE_SHORT <- c(constant = "Constant", `white noise` = "White", `red noise` = "Red")
DISP_LAB  <- c(isolation = "Isolation", dispersal = "Dispersal")

# Hue families of the earlier figures (orange, blue, crimson), re-stepped so
# all three pass colour-vision-deficiency separation and 3:1 contrast on white.
COL <- c(constant = "#c9821f", `white noise` = "#2f4b9a", `red noise` = "#b80422")
INK <- "grey15"; MUTED <- "grey45"; RULE <- "grey80"

theme_ms <- function(base = 7) {
  theme_classic(base_size = base) %+replace%
    theme(
      axis.line        = element_line(linewidth = 0.3, colour = INK),
      axis.ticks       = element_line(linewidth = 0.3, colour = INK),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text        = element_text(size = base, colour = INK),
      axis.title       = element_text(size = base, colour = INK),
      strip.background = element_blank(),
      strip.text       = element_text(size = base, face = "bold", colour = INK,
                                      margin = margin(1, 1, 3, 1)),
      panel.spacing    = unit(3, "mm"),
      legend.position  = "none",
      legend.title     = element_blank(),
      legend.text      = element_text(size = base, colour = INK),
      legend.key.size  = unit(3, "mm"),
      plot.tag         = element_text(size = base + 2, face = "bold"),
      plot.margin      = margin(2, 3, 2, 2)
    )
}

# PDF uses the base pdf() device (cairo_pdf needs XQuartz on macOS); its fonts
# are Latin-1, so Greek letters go through plotmath expressions, not Unicode.
save_fig <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w_mm, height = h_mm,
         units = "mm", device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w_mm, height = h_mm,
         units = "mm", dpi = 600, bg = "white")
  cat("wrote", name, "\n")
}

fx <- function(d) d %>% mutate(noise = factor(noise, levels = NOISE_LEV),
                               treatment = factor(treatment, levels = DISP_LEV))

# same rule as prep_days()$primary in 03_models.R
primary <- function(d) bind_rows(
  d %>% filter(day == 30, !microcosmID %in% CONTAM),
  d %>% filter(day == 24,  microcosmID %in% CONTAM)) %>%
  mutate(standin = microcosmID %in% CONTAM)

rd <- function(f) read_csv(file.path(OUT_DIR, f), show_col_types = FALSE)

emm_tab <- rd("model_emmeans.csv") %>% fx()
get_emm <- function(key, back = identity) {
  emm_tab %>% filter(response == key) %>%
    mutate(across(c(emmean, lower.CL, upper.CL), back))
}

# ---------------------------------------------------------------------------
# Shared panel: one trait by noise regime, faceted by connectivity.
#   mc     one row per microcosm: noise, treatment, y, standin
#   emm    model marginal means with 95% CI (emmean, lower.CL, upper.CL)
#   clones optional, one row per clone: noise, treatment, y
# Clones are faint dots, microcosm means are outlined circles (triangles for
# day-24 stand-ins), and the model mean with its 95% CI sits to the right.
# ---------------------------------------------------------------------------
compare_panel <- function(mc, emm, ylab, clones = NULL, ref = NULL) {
  p <- ggplot(mc, aes(noise, y))
  if (!is.null(ref)) {
    p <- p +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = ref$lo, ymax = ref$hi,
               fill = "grey90") +
      geom_hline(yintercept = ref$mid, colour = MUTED, linewidth = 0.3)
  }
  if (!is.null(clones)) {
    p <- p + geom_point(data = clones, aes(colour = noise),
                        position = position_jitter(width = 0.14, height = 0, seed = 1),
                        size = 0.5, alpha = 0.35, stroke = 0)
  }
  p +
    geom_point(aes(fill = noise, shape = standin), colour = "white", stroke = 0.3,
               size = 1.7, position = position_jitter(width = 0.07, height = 0, seed = 2)) +
    geom_linerange(data = emm, aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                   colour = INK, linewidth = 0.45, position = position_nudge(x = 0.32)) +
    geom_point(data = emm, aes(y = emmean), colour = INK, size = 1.3,
               position = position_nudge(x = 0.32)) +
    facet_wrap(~ treatment, labeller = as_labeller(DISP_LAB)) +
    scale_colour_manual(values = COL, guide = "none") +
    scale_fill_manual(values = COL, guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24),
                       labels = c(`FALSE` = "Day-30 clones",
                                  `TRUE` = "Day-24 clones (day-30 sample contaminated)"),
                       drop = TRUE) +
    guides(shape = guide_legend(override.aes = list(fill = MUTED, colour = "white",
                                                    size = 1.8))) +
    scale_x_discrete(labels = NOISE_SHORT, expand = expansion(add = c(0.5, 0.7))) +
    labs(x = NULL, y = ylab) +
    theme_ms()
}

log2_axis <- function(lims = c(-1, 8.2)) {   # MIC tested 0.5-256 ug/ml
  scale_y_continuous(breaks = -1:8, labels = format(2^(-1:8), drop0trailing = TRUE, trim = TRUE),
                     limits = lims)
}

# draw the most recently added layer first (reference lines under the data)
send_back <- function(p) {
  n <- length(p$layers); p$layers <- c(p$layers[n], p$layers[-n]); p
}

# ===========================================================================
# Figure 1 — schedules
# ===========================================================================
read_sched <- function(f) {
  x <- read.table(file.path("schedules", f), sep = "\t", dec = ",", header = TRUE)
  names(x) <- paste0("p", seq_len(ncol(x)))
  as_tibble(x) %>% mutate(day = row_number()) %>%
    pivot_longer(-day, names_to = "patch", values_to = "kan")
}
white <- read_sched("white_noise_mean2_sd1.txt") %>% mutate(noise = "white noise")
red   <- read_sched("fixed_corr_red_mean2_sd1_corr0.8.txt") %>% mutate(noise = "red noise")

# red patch k is a permutation of white patch k: check before plotting
stopifnot(all(sapply(paste0("p", 1:6), function(k)
  isTRUE(all.equal(sort(white$kan[white$patch == k]), sort(red$kan[red$patch == k]))))))

# Constant patches: summary of the pooled noise values, as in scenarios.R. The
# patch order (C1 = min ... C6 = max) follows the daily dosing sheet
# (1.5, 8.8, 13.0, 12.9, 16.4, 29.3 ul of 1 mg/ml stock in 6 ml).
s <- summary(c(white$kan, red$kan))
const <- tibble(patch = paste0("p", 1:6),
                lab   = c("min", "Q1", "median", "mean", "Q3", "max"),
                kan   = unname(c(s["Min."], s["1st Qu."], s["Median"], s["Mean"],
                                 s["3rd Qu."], s["Max."])))
const_series <- const %>% select(patch, kan) %>%
  tidyr::expand_grid(day = 1:30) %>% mutate(noise = "constant")

sched <- bind_rows(const_series, white, red) %>%
  mutate(noise = factor(noise, levels = NOISE_LEV),
         patch_lab = factor(sub("p", "Patch ", patch), levels = paste("Patch", 1:6)))

lag1 <- function(x) cor(x[-1], x[-length(x)])
longest_run <- function(x) { r <- rle(x > median(x)); max(r$lengths[r$values]) }
sched_stats <- bind_rows(white, red) %>%
  group_by(noise, patch) %>% arrange(day, .by_group = TRUE) %>%
  summarise(r1 = lag1(kan), run = longest_run(kan), .groups = "drop") %>%
  mutate(noise = factor(noise, levels = NOISE_LEV),
         patch_lab = factor(sub("p", "Patch ", patch), levels = paste("Patch", 1:6)))

f1a <- ggplot(sched, aes(day, kan, colour = noise)) +
  geom_line(linewidth = 0.4) +
  geom_text(data = sched_stats,
            aes(x = 1, y = 5.55, label = sprintf("italic(r)[1]=='%.2f'", r1)),
            parse = TRUE, hjust = 0, vjust = 1, size = 2.1, colour = INK) +
  geom_text(data = const %>% mutate(noise = factor("constant", levels = NOISE_LEV),
                                    patch_lab = factor(sub("p", "Patch ", patch),
                                                       levels = paste("Patch", 1:6))),
            aes(x = 1, y = 5.55, label = sprintf("%s  %.2f", lab, kan)),
            hjust = 0, vjust = 1, size = 2.1, colour = INK) +
  facet_grid(noise ~ patch_lab, labeller = labeller(noise = NOISE_LAB)) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = c(1, 10, 20, 30)) +
  scale_y_continuous(breaks = 0:5, limits = c(0, 5.6)) +
  labs(x = "Day", y = "Kanamycin (µg/ml)") +
  theme_ms() + theme(panel.spacing.x = unit(2, "mm"))

runs <- sched_stats %>% mutate(noise = droplevels(noise))
f1b <- ggplot(runs, aes(noise, run)) +
  geom_line(aes(group = patch), colour = RULE, linewidth = 0.4) +
  geom_point(aes(fill = noise), shape = 21, colour = "white", stroke = 0.3, size = 1.9) +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = NOISE_SHORT, expand = expansion(add = 0.35)) +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10.5)) +
  labs(x = NULL, y = "Longest run of days above\nthe patch median") +
  theme_ms()

fig1 <- (f1a | f1b) + plot_layout(widths = c(6, 1.1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold"))
save_fig(fig1, "fig1_schedules", 180, 95)

# ===========================================================================
# Figure 2 — population dynamics and variability of bacterial density
# ===========================================================================
dyn <- rd("dynamics_long.csv") %>% fx()
ND  <- 1.3   # where the one measured-but-zero phage count is drawn

dyn_panel <- function(d, ylab, breaks, labels, lims) {
  ggplot(d, aes(day, y, colour = noise)) +
    geom_line(aes(group = microcosmID), linewidth = 0.3, alpha = 0.75) +
    geom_point(size = 0.55, stroke = 0) +
    facet_grid(noise ~ treatment,
               labeller = labeller(noise = NOISE_LAB, treatment = DISP_LAB)) +
    scale_colour_manual(values = COL, guide = "none") +
    scale_x_continuous(breaks = c(2, 10, 20, 30)) +
    scale_y_continuous(breaks = breaks, labels = labels, limits = lims) +
    labs(x = "Day", y = ylab) +
    theme_ms() + theme(panel.spacing = unit(2, "mm"))
}

cfu_d <- dyn %>% mutate(y = log10(CFU_ml))          # NA where not measured: line breaks
f2a <- dyn_panel(cfu_d, expression("Bacteria ("*log[10]*" CFU ml"^-1*")"),
                 breaks = c(2, 4, 6, 8, 10), labels = c(2, 4, 6, 8, 10), lims = c(1.2, 10.6))

pfu_d <- dyn %>% mutate(y = ifelse(phage_absent, NA, log10(PFU_ml)))
f2b <- dyn_panel(pfu_d, expression("Phage ("*log[10]*" PFU ml"^-1*")"),
                 breaks = c(ND, 4, 6, 8, 10), labels = c("n.d.", 4, 6, 8, 10),
                 lims = c(1, 11.2)) +
  geom_point(data = dyn %>% filter(phage_absent), aes(y = ND), shape = 4,
             size = 1.2, stroke = 0.5)

vari <- rd("microcosm_cfu_variability.csv") %>% fx() %>%
  mutate(y = sd_logCFU, standin = FALSE)
f2c <- compare_panel(vari, get_emm("sd_logCFU"),
                     expression("SD of "*log[10]*" CFU ml"^-1)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

fig2 <- ((f2a | f2b) / (f2c | plot_spacer())) +
  plot_layout(heights = c(2.6, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold"))
save_fig(fig2, "fig2_dynamics", 180, 175)

# ===========================================================================
# Figure 3 — kanamycin resistance
# ===========================================================================
micp <- rd("trait_mic_plate.csv") %>% filter(MIC_use > 0) %>% fx() %>% primary() %>%
  mutate(y = log2(MIC_use))
mic_mc <- micp %>% group_by(microcosmID, noise, treatment, standin) %>%
  summarise(y = mean(y), .groups = "drop")
f3a <- compare_panel(mic_mc, get_emm("mic_plate"), "Kanamycin MIC (µg/ml)",
                     clones = micp) + log2_axis(c(1, 8.3)) + theme(legend.position = "bottom")
stopifnot(min(micp$y) >= 1)   # nothing in the primary data below 2 ug/ml

pairs_mic <- rd("paired_red_white.csv") %>% filter(response == "mic_plate") %>%
  pivot_longer(c(white, red), names_to = "noise", values_to = "y") %>%
  mutate(noise = factor(recode(noise, white = "white noise", red = "red noise"),
                        levels = c("white noise", "red noise")),
         treatment = factor(treatment, levels = DISP_LEV))
f3b <- ggplot(pairs_mic, aes(noise, y)) +
  geom_line(aes(group = patch), colour = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = noise), shape = 21, colour = "white", stroke = 0.3, size = 1.9) +
  facet_wrap(~ treatment, labeller = as_labeller(DISP_LAB)) +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = NOISE_SHORT, expand = expansion(add = 0.4)) +
  log2_axis(c(1, 8.3)) +
  labs(x = NULL, y = "Kanamycin MIC, microcosm mean (µg/ml)") +
  theme_ms()

fig3 <- (f3a | f3b) + plot_layout(widths = c(1.45, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold"))
save_fig(fig3, "fig3_mic", 180, 85)

# ===========================================================================
# Figure 4 — competitive fitness
# ===========================================================================
sel <- rd("trait_selection.csv") %>% fx() %>% primary() %>% mutate(y = selection_coef)
fig4 <- compare_panel(sel, get_emm("selection_coef"),
                      expression("Selection coefficient, "*italic(s)*" (per 24 h)")) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.3) +
  theme(legend.position = "bottom")
fig4 <- send_back(fig4)
save_fig(fig4, "fig4_selection", 100, 80)

# ===========================================================================
# Figure 5 — phage sensitivity against the ancestor
# ===========================================================================
anc_auc <- rd("ancestral_phage_sensitivity.csv") %>% filter(strain == 607)
auc <- rd("trait_phage_sensitivity.csv") %>% fx() %>% primary() %>%
  mutate(y = phage_sensitivity)
auc_mc <- auc %>% group_by(microcosmID, noise, treatment, standin) %>%
  summarise(y = mean(y), .groups = "drop")
ref_auc <- list(lo = min(anc_auc$phage_sensitivity), hi = max(anc_auc$phage_sensitivity),
                mid = median(anc_auc$phage_sensitivity))
fig5 <- compare_panel(auc_mc, get_emm("phage_sensitivity"),
                      expression(atop(displaystyle("Phage sensitivity, "*Delta*"AUC (OD600 × min)"),
                                      displaystyle("higher = more sensitive"))),
                      clones = auc, ref = ref_auc) +
  geom_hline(yintercept = 0, colour = RULE, linewidth = 0.3)
fig5 <- send_back(fig5) +
  annotate("text", x = 0.55, y = ref_auc$hi, label = "Ancestor REL607",
           hjust = 0, vjust = 1.4, size = 2.1, colour = MUTED) +
  theme(legend.position = "bottom")
save_fig(fig5, "fig5_phage", 100, 85)

# ===========================================================================
# Figure S1 — growth parameters
# ===========================================================================
anc_g <- rd("ancestral_growth.csv") %>% filter(strain == 607)
grw <- rd("trait_growth.csv") %>% fx() %>% primary()

gen <- grw %>% mutate(y = generation_el)
gen_mc <- gen %>% group_by(microcosmID, noise, treatment, standin) %>%
  summarise(y = exp(mean(log(y))), .groups = "drop")      # geometric mean, as modelled
fs1a <- compare_panel(gen_mc, get_emm("generation_time", back = exp), "Generation time (min)",
                      clones = gen,
                      ref = list(lo = min(anc_g$generation_el), hi = max(anc_g$generation_el),
                                 mid = median(anc_g$generation_el))) +
  scale_y_log10(breaks = c(10, 20, 30, 50, 100, 200)) +
  theme(legend.position = "bottom")

kk <- grw %>% mutate(y = K_obs)
kk_mc <- kk %>% group_by(microcosmID, noise, treatment, standin) %>%
  summarise(y = mean(y), .groups = "drop")
fs1b <- compare_panel(kk_mc, get_emm("carrying_capacity"), "Carrying capacity (max. OD600)",
                      clones = kk,
                      ref = list(lo = min(anc_g$K_obs), hi = max(anc_g$K_obs),
                                 mid = median(anc_g$K_obs))) +
  theme(legend.position = "bottom")

figS1 <- (fs1a | fs1b) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"), legend.position = "bottom")
save_fig(figS1, "figS1_growth", 180, 80)

# ===========================================================================
# Figure S2 — MIC from gradient strips
# ===========================================================================
strip <- rd("trait_mic_stripe.csv") %>% fx() %>% primary() %>% mutate(y = log2(meanMIC))
figS2 <- compare_panel(strip, get_emm("mic_strip"), "Kanamycin MIC, gradient strip (µg/ml)") +
  log2_axis(c(0, 7)) + theme(legend.position = "bottom")
save_fig(figS2, "figS2_mic_strip", 100, 80)

# ===========================================================================
# Figure S3 — host-phage correlation on days with both counts
# ===========================================================================
cc <- rd("microcosm_host_phage_corr.csv") %>% fx() %>% mutate(y = corr_at_0, standin = FALSE)
figS3 <- compare_panel(cc, get_emm("host_phage_corr"),
                       "Correlation of log10 CFU and log10 PFU") +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.3) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5))   # CIs from the linear model can pass -1
figS3 <- send_back(figS3)
save_fig(figS3, "figS3_host_phage_corr", 100, 70)

# ---------------------------------------------------------------------------
# numbers quoted in the figure legends
# ---------------------------------------------------------------------------
cat("\n--- for the legends ---\n")
cat("schedule lag-1 autocorrelation, white:", round(range(sched_stats$r1[sched_stats$noise == "white noise"]), 2),
    " red:", round(range(sched_stats$r1[sched_stats$noise == "red noise"]), 2), "\n")
cat("longest run above median, white:", range(runs$run[runs$noise == "white noise"]),
    " red:", range(runs$run[runs$noise == "red noise"]), "\n")
cat("MIC: clones", nrow(micp), "| microcosms", n_distinct(micp$microcosmID),
    "| censored at 256:", sum(micp$censored), "| red-white pairs:",
    n_distinct(paste(pairs_mic$patch, pairs_mic$treatment)),
    "| pairs with red > white:", sum(rd("paired_red_white.csv") %>%
                                       filter(response == "mic_plate") %>% pull(diff) > 0), "\n")
cat("phage sensitivity: clones", nrow(auc), "| REL607 wells", nrow(anc_auc),
    "| REL607 median", round(ref_auc$mid), "range", round(ref_auc$lo), "-", round(ref_auc$hi), "\n")
cat("  microcosm means above REL607 minimum:", sum(auc_mc$y > ref_auc$lo), "of", nrow(auc_mc),
    "| max microcosm mean", round(max(auc_mc$y)), "\n")
cat("  clones above REL607 minimum:", sum(auc$y > ref_auc$lo), "of", nrow(auc), "\n")
cat("growth: clones", nrow(grw), "| REL607 generation time median", round(median(anc_g$generation_el), 1),
    "range", round(range(anc_g$generation_el), 1), "| K median", round(median(anc_g$K_obs), 3),
    "range", round(range(anc_g$K_obs), 3), "\n")
cat("host-phage correlation: microcosms", nrow(cc), "| paired days", range(cc$n_paired), "\n")
cat("dynamics: PFU measured-but-zero:", sum(dyn$phage_absent), "\n")
