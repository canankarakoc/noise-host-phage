# Bacterial growth models / Virus infectivity 
# 6th May 2021 Karakoc 
# Klara Bachelor Thesis 

# Updating 2 May 2025

# install.packages("growthrates")
# see https://tpetzoldt.github.io/growthrates/doc/Introduction.html

library(growthrates)
library(tidyverse)
library(readr)
library(ggpubr)

set.seed(5621)

mycolors = c("#19798b", "#f3d567","#e74b47", "#ee9b43",  "#b80422", "#172767")

# Theme 
# panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
mytheme <- theme_bw() +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 16)) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    size = 1
  )) +
  theme(strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14), 
        strip.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    linewidth = 1
  )) +
  theme(
    axis.text.x.top = element_blank(), 
    axis.title.x.top = element_blank(),
    axis.text.y.right = element_blank(), 
    axis.title.y.right = element_blank()
  ) +
  theme(
    axis.title.x = element_text(margin = margin(10, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.text.y = element_text(margin = margin(0, 10, 0, 0))
  )

se <- function(data) {
  # Remove NA values before calculation
  data <- na.omit(data)
  sd(data) / sqrt(length(data))
}

# run from the repository folder; the data files go in data/ (see data/README.md)

# read data
# experimental design
exp_design <- read.table("data/design.csv", header = T, sep = ",")
# plate design
design     <- read.table("data/plate_design.csv", header = T, sep = ",")

b <- design[duplicated(design[, c("microcosmID", "sampleID", "day", "phage_treatment", "clones")]), ]

# OD data from plate reader 
files  <- list.files(path = "data/files", pattern = "*.csv", full.names = T)
# combine all the plate data with the plate identifier "set"
ldf <- list() 
# read.table function doesn't work here, because both separators and decimals are commas
for (i in 1:length(files)){
  ldf[[i]]    <- read.csv2(files[i], header = T, sep = ",", dec = ".")
  ldf[[i]]$set <- i
}
str(ldf[[1]]) 

# rbind list of sets as data frame 
df <- do.call(rbind.data.frame, ldf)

# transform data to longitudinal format 
longdata <- gather(df, well, OD, -time, -set)
# combine data with plate design
final_data <- left_join(longdata, design, by = c("well", "set"))

# subset the blank readings 
blank <- subset(final_data, final_data$microcosmID == "blank")

# mean the blanks by the plate id(set) and time 
blank_means <- blank %>%
  group_by(set, time) %>%
  mutate(mean_blank = mean(OD)) %>%
  select(time, set, mean_blank)

# extract blanks from readings, clean data, (median of clones) 
data_sum <- final_data %>%
  left_join(blank_means, by = c("set", "time")) %>%
  mutate(blanked_OD = OD - mean_blank) %>%
  mutate(OD_cor = if_else(blanked_OD < 0, 0.001, blanked_OD)) %>% #set negative readings to a small number
  filter(time > 0) %>%
  distinct_all() %>%
  full_join(exp_design, by = c("microcosmID", "sampleID")) %>% #join table with experimental design
  filter(microcosmID != "blank") %>%
  drop_na()  #it drops the ancestral control, i will plot it separately 

# HERE P4-P6 of isolation-white noise plate were contaminated on day 25-26. 
# So clones from day 24 can be used as backup frozen record. 

# Identify contaminated Day 30 patches replace with day 24 frozen record
missing_patches <- c("p3", "p4", "p5", "p6")

replacements <- data_sum %>%
  filter(day == 24,
         noise == "white noise",
         treatment == "isolation",
         patch %in% missing_patches)

data_sum_30  <- data_sum %>%
  filter(day == 30)

data_sum_24  <- data_sum %>%
  filter(day == 24)

final_clones      <- bind_rows(data_sum_30, replacements)
final_clones_ALL  <- bind_rows(data_sum_30, data_sum_24) # this is another option. Discuss. 
# in order to reduce the bias and get quality measurements we can pool last two frozen records 
# (we have data for last 3) 
# at the end we want to see the effect of patch environments.

# plot the raw growth data
ggplot(final_clones_ALL%>%drop_na(), aes(x=time, y = OD_cor, linetype = phage_treatment, color = clones))+
  geom_line()+
  facet_grid(patch~noise+treatment+day)+
  xlab("Time (minutes)")+
  ylab("Mean OD of clones")+
  scale_linetype_manual("Phage treatment", values = c("without_phage" = "solid", "with_phage" = "dashed"), 
                        labels= c("with_phage" = "With phage", "without_phage" = "Without phage"))+
  mytheme+
  theme(axis.text.x = element_text(angle = 90))

# Note the ones are not working. 

# plot ancestral control to visualize the control and the phage infection
control <- final_data %>%
  filter(microcosmID == "606") %>%
  left_join(blank_means, by = c("set", "time")) %>%
  mutate(blanked_OD = OD - mean_blank) %>%
  mutate(OD_cor = if_else(blanked_OD < 0, 0.001, blanked_OD)) %>%  # Set negative OD to small value
  #filter(time < 1000) %>%
  filter(phage_treatment == "without_phage")  %>%
  filter(set == 1) # set 2 has not reached saturation 

# below this is unrealistic that there is some physical limitation. they might cause bias. 
# since there are enough replicate wells, I remove them. 

# ancestral is phage sensitive. 

# plot the raw growth data
ggplot(control, aes(x=time, y = OD_cor, linetype = phage_treatment, color = clones))+
  geom_line(aes())+
  xlab("Time (minutes)")+
  ylab("Mean OD of clones")+
  scale_linetype_manual("Phage treatment", values = c("without_phage" = "solid", "with_phage" = "dashed"), 
                        labels= c("with_phage" = "With phage", "without_phage" = "Without phage"))+
  mytheme

# reshape data for growth models 
# to calculate intrinsic growth rate and carrying capacity, we need the growth curves without phages

# better parameter estimation?
data_mod_nophage <- final_clones_ALL %>% 
  filter(trimws(phage_treatment) != "with_phage") %>% 
  #filter(time < 1000) %>%  #OD reduces and messes up with the models
  drop_na()  %>%
  mutate(across(where(is.character), str_trim))

# Split dataset by clones and microcosmID
grouped_data <- data_mod_nophage %>%
  group_by(clones, microcosmID, day) %>%
  group_split()

# Define a helper function that forces p into the interior of (lower, upper)
adjust_p <- function(p, lower, upper) {
  p_adjusted <- p  # copy
  for(i in seq_along(p)) {
    if(p_adjusted[i] <= lower[i] || p_adjusted[i] >= upper[i]) {
      p_adjusted[i] <- (lower[i] + upper[i]) / 2
    }
  }
  return(p_adjusted)
}

fit_all_models <- function(df) {
  clone <- unique(df$clones)
  micro <- unique(df$microcosmID)
  day   <- unique(df$day)  # NEW
  
  if (nrow(df) < 5 || sd(df$OD_cor) < 0.01) return(NULL)
  
  # Guess parameters
  y0 <- max(min(df$OD_cor), 0.0005)
  K <- min(max(df$OD_cor), 0.65)
  rate <- max(diff(df$OD_cor) / diff(df$time), na.rm = TRUE)
  rate <- min(max(rate, 0.001), 0.05)
  
  try_models <- list()
  
  try_models$gompertz <- try(fit_growthmodel(FUN = grow_gompertz,
                                             p = c(y0 = y0, mumax = rate, K = K),
                                             time = df$time, y = df$OD_cor,
                                             lower = c(y0 = 0.0001, mumax = 0.0005, K = 0.3),
                                             upper = c(y0 = 0.05, mumax = 0.1, K = 0.7),
                                             transform = "log"
  ), silent = TRUE)
  
  try_models$baranyi <- try(fit_growthmodel(FUN = grow_baranyi,
                                            p = c(y0 = y0, mumax = rate, K = K, h0 = 1),
                                            time = df$time, y = df$OD_cor,
                                            lower = c(y0 = 0.0001, mumax = 0.0005, K = 0.3, h0 = 0),
                                            upper = c(y0 = 0.05, mumax = 0.1, K = 0.7, h0 = 10),
                                            transform = "log"
  ), silent = TRUE)
  
  try_models$logistic <- try(fit_growthmodel(FUN = grow_logistic,
                                             p = c(y0 = y0, mumax = rate, K = K),
                                             time = df$time, y = df$OD_cor,
                                             lower = c(y0 = 0.0001, mumax = 0.0005, K = 0.3),
                                             upper = c(y0 = 0.05, mumax = 0.1, K = 0.65),
                                             transform = "log"
  ), silent = TRUE)
  
  try_models$weibull <- try(fit_growthmodel(
    FUN = grow_weibull,
    p = c(y0 = y0, mumax = rate, K = K, m = 1),  # m = shape param
    time = df$time, y = df$OD_cor,
    lower = c(y0 = 0.0001, mumax = 0.0005, K = 0.3, m = 0.1),
    upper = c(y0 = 0.05,   mumax = 0.1,    K = 0.7, m = 5),
    transform = "log"
  ), silent = TRUE)
  
  # Filter out failures and select best by AIC
  fits <- keep(try_models, ~ !inherits(.x, "try-error"))
  if (length(fits) == 0) return(NULL)
  
  # best_name <- names(which.min(map_dbl(fits, AIC)))
  
  # New working block
  extract_aic <- function(fit_obj) {
    resid <- residuals(fit_obj)
    n <- length(resid)
    rss <- sum(resid^2)
    k <- length(coef(fit_obj))
    aic <- n * log(rss / n) + 2 * k
    return(aic)
  }
  
  # extract_aic defined above
  aic_values <- map_dbl(fits, extract_aic)
  best_name <- names(which.min(aic_values))
  best_fit <- fits[[best_name]]
  
  return(list(
    fit = best_fit,
    model = best_name,
    clone = clone,
    microcosmID = micro,
    day = day,  # NEW
    input_df = df
  ))
}

# Re-run the loop
fits <- purrr::map(grouped_data, fit_all_models)
fits <- purrr::compact(fits)

extract_fit_plot_data <- function(fit_obj) {
  df <- fit_obj$input_df
  model <- fit_obj$model
  label <- paste(fit_obj$clone, fit_obj$microcosmID, fit_obj$day, model, sep = "_")
  
  pred_time <- seq(min(df$time), max(df$time), length.out = 100)
  pred_OD <- predict(fit_obj$fit, newdata = data.frame(time = pred_time))
  
  raw_df <- df %>% mutate(label = label)
  pred_df <- data.frame(time = pred_time, OD_pred = pred_OD, label = label)
  
  list(raw = raw_df, pred = pred_df)
}

plot_data <- map(fits, extract_fit_plot_data)

raw_df  <- bind_rows(map(plot_data, "raw"))
pred_df <- bind_rows(map(plot_data, "pred"))

labels <- unique(raw_df$label)
chunks <- split(labels, ceiling(seq_along(labels)))

if (!dir.exists("figures")) dir.create("figures")

save_fit_plot <- function(fit_obj, outdir = "figures") {
  df <- fit_obj$input_df
  pred_time <- seq(min(df$time), max(df$time), length.out = 100)
  pred_OD <- predict(fit_obj$fit, newdata = data.frame(time = pred_time))
  
  pred_df <- data.frame(time = pred_time, OD_pred = pred_OD)
  
  p <- ggplot() +
    geom_point(data = df, aes(x = time, y = OD_cor), alpha = 0.6) +
    geom_line(data = pred_df, aes(x = time, y = OD_pred.y), color = "red") +
    labs(
      title = paste("Clone:", fit_obj$clone,
                    "| Microcosm:", fit_obj$microcosmID,
                    "| Model:", fit_obj$model, 
                    "| Day:", fit_obj$day),
      x = "Time",
      y = "OD600"
    ) +
    theme_minimal()
  
  filename <- paste0(outdir, "/", fit_obj$day, "_", fit_obj$microcosmID, "_", fit_obj$clone,".pdf")
  ggsave(filename, plot = p, width = 6, height = 4)
}

walk(fits, save_fit_plot)

extract_growth_parameters <- function(fit_obj) {
  coefs <- coef(fit_obj$fit)
  model <- fit_obj$model
  clone <- fit_obj$clone
  micro <- fit_obj$microcosmID
  day   <- fit_obj$day  # NEW
  
  # Extract common ones
  y0     <- coefs["y0"]
  mumax  <- coefs["mumax"]
  K      <- coefs["K"]
  
  # Helper function for operational lag (defined inside outer function)
  estimate_lag_operational <- function() {
    threshold <- y0 + 0.1 * (K - y0)
    time_seq <- seq(0, max(fit_obj$input_df$time), length.out = 500)
    preds <- predict(fit_obj$fit, newdata = data.frame(time = time_seq))
    
    # Check if threshold is ever crossed
    above_thresh <- which(preds >= threshold)
    
    if (length(above_thresh) == 0) {
      return(NA_real_)  # Explicit NA, or max(time_seq) if you prefer
    } else {
      return(time_seq[above_thresh[1]])
    }
  }
  
  # Lag time: model-specific or inferred
  if (model == "baranyi") {
    lag <- coefs["h0"]
  } else {
    lag <- estimate_lag_operational()
  }
  
  data.frame(
    clone = clone,
    microcosmID = micro,
    day = day,  # NEW
    model = model,
    y0 = y0,
    mumax = mumax,
    K = K,
    lag = lag
  )
}

param_df <- map_dfr(fits, extract_growth_parameters)

extract_fit_quality <- function(fit_obj) {
  resid <- residuals(fit_obj$fit)
  rmse <- sqrt(mean(resid^2))
  rss <- sum(resid^2)
  data.frame(
    clone = fit_obj$clone,
    microcosmID = fit_obj$microcosmID,
    day = fit_obj$day,  # NEW
    model = fit_obj$model,
    rmse = rmse,
    rss = rss
  )
}

quality_df <- map_dfr(fits, extract_fit_quality)
param_df <- map_dfr(fits, extract_growth_parameters)
param_df <- left_join(param_df, quality_df, by = c("day", "microcosmID", "clone", "model"))
param_df_filtered <- filter(param_df, rmse <= 0.1) #filter out off fits + visually checked the fits

param_df_extended <- param_df_filtered %>%
  left_join(exp_design, by = c("microcosmID")) %>%
  mutate(generation = 0.69 / mumax)

# Group data by clone
grouped_control <- control %>%
  group_by(set, clones) %>%
  group_split()

# Fit models
fit_cont <- map(grouped_control, fit_all_models)
fit_cont <- compact(fit_cont)

param_cont <- map_dfr(fit_cont, extract_growth_parameters)
param_cont$treatment <- "control"

quality_df_cont <- map_dfr(fit_cont, extract_fit_quality)
param_cont <- cbind.data.frame(param_cont, quality_df_cont[,5:6])
param_df_filtered_cont <- filter(param_cont, rmse <= 0.1) #filter out off fits + visually checked the fits

control_range <- param_df_filtered_cont %>%
  summarise(
    gen_min = min(0.69 / mumax, na.rm = TRUE),
    gen_max = max(0.69 / mumax, na.rm = TRUE),
    K_min = min(K, na.rm = TRUE),
    K_max = max(K, na.rm = TRUE),
    lag_min = min(lag, na.rm = TRUE),
    lag_max = max(lag, na.rm = TRUE), 
    gen_mean = mean(0.69 / mumax, na.rm = TRUE), 
    gen_se = se(0.69 / mumax), 
    K_mean = mean(K, na.rm = TRUE), 
    K_se = se(K), 
    lag_mean = mean(lag, na.rm = TRUE), 
    lag_se = se(lag), 
    gen_med = median(0.69 / mumax, na.rm = TRUE), 
    K_med = median(K), 
    lag_med = median(lag)
  )

gen_mean <- control_range$gen_mean
gen_se   <- control_range$gen_se
K_mean   <- control_range$K_mean
K_se     <- control_range$K_se
lag_mean <- control_range$lag_mean
lag_se   <- control_range$lag_se

param_summary <- param_df_extended %>%
  group_by(noise, treatment, day) %>%
  summarise(
    meangeneration = mean(generation, na.rm = TRUE),
    segeneration = se(generation),
    meanK = mean(K, na.rm = TRUE),
    seK = se(K),
    meanlag = mean(lag, na.rm = TRUE),
    selag = se(lag),
    .groups = "drop"
  )

param_df_extended$noise <- factor(param_df_extended$noise , levels = rev(c("constant", "white noise", "red noise")))

generation_time <- ggplot(param_df_extended, aes(y = noise, x = generation)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = control_range$gen_min, xmax = control_range$gen_max),
            fill = "gray95", inherit.aes = FALSE) +
  geom_vline(xintercept = control_range$gen_med, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = patch, shape = factor(day)), size = 2, alpha = .7, width = 0.2) +
  facet_wrap(~treatment, ncol = 3) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(17, 16)) +
  labs(x = "Generation time (min)", y = NULL) +
  mytheme 

carrying_capacity <- ggplot(param_df_extended, aes(y = noise, x = K)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = control_range$K_min, xmax = control_range$K_max),
            fill = "gray95", inherit.aes = FALSE) +
  geom_vline(xintercept = control_range$K_med, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = patch, shape = factor(day)), size = 2, alpha = .7, width = 0.2) +
  facet_wrap(~treatment, ncol = 3) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(17, 16)) +
  labs(x = "Carrying capacity (K)", y = NULL) +
  mytheme 

##############VIRUS INFECTIVITY###################
#filtering the curves not normal to prevent biased estimations 
control_max <- final_clones_ALL %>%
  filter(phage_treatment == "without_phage") %>%
  group_by(clones, patch, noise, treatment, day) %>%
  summarise(max_OD = max(OD_cor), .groups = "drop") %>%
  filter(max_OD > 0.2) 

filtered_data_inf <- final_clones %>%
  inner_join(control_max, by = c("clones", "patch", "noise", "treatment", "day")) 

library(pracma)
# Compute AUC per group
auc_data <- filtered_data_inf %>%
  group_by(clones, patch, noise, treatment, phage_treatment, day) %>%
  arrange(time) %>%
  summarise(
    auc = trapz(time, OD_cor),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = phage_treatment, values_from = auc) %>%
  mutate(phage_sensitivity = without_phage - with_phage)

control_all <- final_data %>%
  filter(microcosmID == "606") %>%
  left_join(blank_means, by = c("set", "time")) %>%
  mutate(blanked_OD = OD - mean_blank) %>%
  mutate(OD_cor = if_else(blanked_OD < 0, 0.001, blanked_OD)) %>%
  filter(set == 1)

# plot to see the curves - all sensitive!
ggplot(control_all, aes(x = time, y= OD_cor, linetype = phage_treatment)) +
  geom_line() +
  facet_wrap(~set+clones, ncol = 3) +
  scale_color_manual(values = mycolors) +
  labs(y = "OD", x = NULL) +
  mytheme

auc_data_control <- control_all %>%
  group_by(set, clones, phage_treatment) %>%
  arrange(time) %>%
  summarise(
    auc = trapz(time, OD_cor),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = phage_treatment, values_from = auc) %>%
  mutate(phage_sensitivity = without_phage - with_phage)

auc_data_control_sum <- auc_data_control %>%
  summarise(
    mean_sen = mean(phage_sensitivity), 
    se_sen = se(phage_sensitivity), 
    median_sen = median(phage_sensitivity), 
    min_sen = min(phage_sensitivity), 
    max_sen = max(phage_sensitivity)
  )

auc_data$noise <- factor(auc_data$noise , levels = c("constant", "white noise", "red noise"))

infectivity <- ggplot(auc_data, aes(y = noise, x = phage_sensitivity)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = auc_data_control_sum$min_sen, xmax = auc_data_control_sum$max_sen),
            fill = "gray95", inherit.aes = FALSE) +
  geom_vline(xintercept = auc_data_control_sum$median_sen, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = patch, shape = factor(day)), width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~treatment, ncol = 3) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(17, 16)) +
  labs(x = "Phage sensitivity (\u0394AUC)", y = NULL) +
  mytheme

figure3 <- ggarrange(generation_time,carrying_capacity,infectivity, ncol =1, labels = c("A","B","C"))

#####END#####

# Klara Thesis Competition assays
# Day30 & 24 & 18
# Author: Canan Karakoc
# 16 March 2021 
# Updated 7 May 2024

comp_data  <- read.table("data/competition_alldays.csv", header=T, sep=",", dec=",")
str(comp_data)
factor(comp_data$time)
factor(comp_data$day)
factor(comp_data$microcosmID)

data_clean <-  comp_data %>%
  dplyr::select(time, microcosmID, day, dilution_fold, pink_colonies, red_colonies) %>%
  drop_na() %>% #remove non-available data
  mutate(sum = pink_colonies+red_colonies) %>%
  filter(sum > 1) %>% #remove singletons
  filter((dilution_fold < 100000 & time == "t0") | (dilution_fold < 100000000 & time == "t1")) %>%#remove highest dilutions 
  mutate(pink_colonies_dl = ifelse(pink_colonies == 0, 1, pink_colonies)) %>%
  mutate(red_colonies_dl  = ifelse(red_colonies  == 0, 1, red_colonies)) %>%
  mutate(red_col_dil = (red_colonies_dl * dilution_fold)/10) %>% # calculate CFU per 1µL (10µl plated) 
  mutate(pink_col_dil = (pink_colonies_dl * dilution_fold)/10) %>%  
  group_by(time, microcosmID, day) %>%
  summarise(ancestral = mean(red_col_dil, na.rm = T), evolved = mean(pink_col_dil, na.rm = T) ) %>%
  ungroup() %>%
  pivot_wider(names_from = "time", 
              values_from = c("ancestral", "evolved")) %>%
  mutate(selection_coef = log(evolved_t1/ancestral_t1) - # calculate selection coefficient
           log(evolved_t0/ancestral_t0)) %>% 
  left_join(exp_design, by = c("microcosmID")) %>%
  dplyr::select(microcosmID, day, sampleID, patch, noise, treatment, selection_coef)
            
# Identify contaminated Day 30 patches
missing_microcosms <- c("m15", "m16", "m17", "m18")

replacements_selection <- data_clean %>%
  filter(day == 24,
         microcosmID %in% missing_microcosms)

data_clean_30  <- data_clean %>%
  filter(day == 30,
         !microcosmID %in% missing_microcosms)

data_clean_24  <- data_clean %>%
  filter(day == 24)

final_selection     <- bind_rows(data_clean_30, replacements_selection)
final_selection_ALL <- bind_rows(data_clean_30, data_clean_24)

# Plot 
final_selection_ALL$noise <- factor(final_selection_ALL$noise , levels = rev(c("constant", "white noise", "red noise")))

selection <- ggplot(final_selection_ALL, aes(y = noise, x = selection_coef)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = patch, shape = factor(day)), width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~treatment, ncol = 3) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(17, 16)) +
  labs(x = "Selection coefficient (s)", y = NULL) +
  mytheme

#########################MIC######################################
# Calculatig minimum inhibitory concentrations 
# 96-well plate tests
#libraryd(drc) # for EC50

micdata_raw <- read.table("data/mic_od_combined.csv", sep = ",", dec = ".", header = T)
glimpse(micdata_raw)

mic_cols <- names(micdata_raw)[grepl("^MIC_", names(micdata_raw))]

# Extract global concentration values (once)
mic_conc_vector <- as.numeric(gsub("MIC_", "", mic_cols))

# Process each row with robust logic
micdata_processed <- micdata_raw %>%
  group_by(plate) %>%
  mutate(blank_mean = mean(negative.control, na.rm = TRUE)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    mic_raw = list(c_across(all_of(mic_cols))),
    
    # Quality checks
    fail_stamp = stamp_control > 0.1,
    fail_negative = negative.control > 0.1,
    
    # Allow 2 upward jumps to tolerate noise
    is_monotonic = sum(diff(unlist(mic_raw)) > 0.1) <= 2,
    
    # MIC logic: compare to own negative control with tolerance
    MIC_estimate = {
      idx <- which(unlist(mic_raw) < (negative.control + 0.01))[1]
      if (!is.na(idx)) mic_conc_vector[idx] else NA_real_
    },
    
    is_resistant = is.na(MIC_estimate),
    MIC_estimate = ifelse(is_resistant, 256, MIC_estimate)  # assign a placeholder
  ) %>%
  ungroup() %>%
  filter(!fail_stamp, !fail_negative, is_monotonic)

micdata_final <- micdata_processed %>%
  left_join(exp_design, by = "sampleID") %>%
  dplyr::select(day, sampleID, MIC_estimate, microcosmID, patch, noise, treatment) %>%
  filter(!day == 18) 

# Identify contaminated Day 30 patches and replace them with fossil reccords
missing_microcosms <- c("m15", "m16", "m17", "m18")

replacements_selection_mic <- micdata_final %>%
  filter(day == 24,
         microcosmID %in% missing_microcosms)

micdata_final_30  <- micdata_final %>%
  filter(day == 30) %>%
  filter(!microcosmID %in% missing_microcosms)

micdata_final_24  <- micdata_final %>%
  filter(day == 24) 

micdata_final     <- bind_rows(micdata_final_30, replacements_selection_mic)
micdata_final_ALL <- bind_rows(micdata_final_30, micdata_final_24)

#plot
micdata_final_ALL$noise <- factor(micdata_final_ALL$noise , levels = c("constant", "white noise", "red noise"))

micPlot <- ggplot(micdata_final_ALL, aes(y = noise, x = MIC_estimate)) +
  geom_vline(xintercept = 6, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = patch, shape = factor(day)), width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~treatment, ncol = 3) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(17, 16)) +
  labs(x = "MIC (plate method)", y = NULL) +
  mytheme

# Check the filtering step 
micdata_raw %>%
  mutate(row_id = row_number()) %>%
  rowwise() %>%
  mutate(
    mic_raw = list(c_across(all_of(mic_cols))),
    max_OD = max(unlist(mic_raw), na.rm = TRUE),
    min_OD = min(unlist(mic_raw), na.rm = TRUE),
    is_flat = max_OD - min_OD < 0.05,
    fail_stamp = stamp_control > 0.1,
    fail_negative = negative.control > 0.1,
    is_monotonic = sum(diff(unlist(mic_raw)) > 0.1) <= 2
  ) %>%
  ungroup() %>%
  count(is_flat, fail_stamp, fail_negative, is_monotonic)


# MIC stripe data 
micdata_stripe <- read.table("data/MIC_stripe_test.csv", sep = ",", dec = ".", header = T)

micdata_stripe <- micdata_stripe %>%
 left_join(exp_design, by = c("microcosmID", "sampleID"))

  replacements_selection_stripe <- micdata_stripe %>%
    filter(Day == 24,
           microcosmID %in% missing_microcosms)
  
  stripe_data_final_30  <- micdata_stripe %>%
    filter(Day == 30) %>%
    filter(!microcosmID %in% missing_microcosms)
  
  stripe_data_final_24  <- micdata_stripe %>%
    filter(Day == 24)
  
stripe_data_final <- bind_rows(replacements_selection_stripe, stripe_data_final_30)  
stripe_data_final_all <- bind_rows(stripe_data_final_24, stripe_data_final_30)    

# Average the technical replicates
  stripe_data_final_all_sum <-stripe_data_final_all %>%
    group_by(microcosmID, patch, noise, treatment, Day) %>%
    summarize(meanMIC = mean(MIC))
  
#Plot 
stripe_data_final_all_sum$noise <- factor(stripe_data_final_all_sum$noise , levels = rev(c("constant", "white noise", "red noise")))    

micPlot_stripe <- ggplot(stripe_data_final_all_sum, aes(y = noise, x = meanMIC)) +
  geom_vline(xintercept = 6, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = patch, shape = factor(Day)), width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~treatment, ncol = 3) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(17, 16)) +
  labs(x = "MIC (plate method)", y = NULL) +
  mytheme

ggarrange(micPlot, micPlot_stripe, ncol = 1, labels = c("A","B","C"))

# Dynamics 
dynamics_cfu <- read.table("data/dynamicsCFU.csv", header = T, sep = ",")
dynamics_pfu <- read.table("data/dynamicsPFU.csv", header = T, sep = ",")

dynamics_cfu_cl <- dynamics_cfu %>%
  mutate(dil_factor = dilution*100) %>%
  mutate(CFU_ml = CFU * dil_factor) %>%
  group_by(day, sampleID) %>%
  summarize(CFU_sum = mean(CFU_ml)) %>%
  mutate(CFU_log = log10(CFU_sum))
  
# Define the sequence of days
all_days <- seq(2, 30, by = 2)

# Expand design × days to get all expected measurements
expected_grid <- exp_design %>%
  distinct(sampleID, patch, noise, treatment, microcosmID) %>%
  tidyr::expand_grid(day = all_days)

cfu_full <- expected_grid %>%
  left_join(dynamics_cfu_cl, by = c("sampleID", "day"))

cfu_filled <- cfu_full %>%
  mutate(CFU_sum_filled = if_else(is.na(CFU_sum), 100, CFU_sum),
         CFU_log_filled = log10(CFU_sum_filled),
         interpolated = is.na(CFU_sum))  # flag filled points

CFU_plot <- ggplot(cfu_filled, aes(x = day, y = CFU_log_filled)) +
  geom_line(size = 1, alpha = 0.6, aes( color = patch))+
  facet_grid(noise~treatment) +
  scale_color_manual(values = mycolors) +
  labs(x = "Days", y = "Colony forming units (CFU/mL)" ) +
  mytheme+
  scale_y_continuous(breaks = c(3,6,9))
  
# PFU 
dynamics_pfu_cl <- dynamics_pfu %>%
  mutate(dil_factor = dilution*100) %>%
  mutate(PFU_ml = count * dil_factor) %>%
  group_by(day, sampleID) %>%
  summarize(PFU_sum = mean(PFU_ml)) %>%
  mutate(PFU_log = log10(PFU_sum+1)) 

pfu_full <- expected_grid %>%
  left_join(dynamics_pfu_cl, by = c("sampleID", "day"))

# PFU interpolated
# a zero plaque count (R1D, day 22) is treated like a missing day (below detection); log10(0) breaks the cross-correlation
pfu_interp <- pfu_full %>%
  mutate(PFU_sum_filled = if_else(is.na(PFU_sum) | PFU_sum == 0, 100, PFU_sum),
         PFU_log_filled = log10(PFU_sum_filled),
         interpolated = is.na(PFU_sum))  # flag filled points

PFU_plot <- ggplot(pfu_interp, aes(x = day, y = PFU_log_filled)) +
  geom_line(size = 1, alpha = 0.6, aes( color = patch))+
  facet_grid(noise~treatment) +
  scale_color_manual(values = mycolors) +
  labs(x = "Days", y = "Plaque making units (PFU/mL)" ) +
  mytheme+
  scale_y_continuous(breaks = c(3,6,9))

ggarrange(CFU_plot, PFU_plot, ncol = 1, labels = c("A","B"))

# time extinction 
# Identify last real PFU observation
pfu_padded_l <- pfu_full %>%
  group_by(sampleID) %>%
  arrange(day) %>%
  mutate(last_obs_day = max(day[!is.na(PFU_sum)], na.rm = TRUE),
         PFU_sum_filled = if_else(day > last_obs_day, 100, PFU_sum),
         PFU_log_filled = log10(PFU_sum_filled)) %>%
  ungroup()

# Calculate TTE based on last day PFU_log was above threshold and after that remains below
tte_clean_df <- pfu_full %>%
  group_by(sampleID, patch, noise, treatment) %>%
  arrange(day) %>%
  summarize(
    time_to_extinction = max(day[!is.na(PFU_log)], na.rm = TRUE),
    .groups = "drop"
  )

tte_clean_df$noise <- factor(tte_clean_df$noise , levels = rev(c("constant", "white noise", "red noise")))
tte_plot <- ggplot(tte_clean_df, aes(y = noise, x = time_to_extinction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(size = 3, alpha = 0.7, aes(color = patch))+
  facet_grid(~treatment) +
  scale_color_manual(values = mycolors) +
  labs(y = NULL, x = "Time to extinction (days)" ) +
  mytheme

# stability 
cfu_stability <- cfu_filled %>%
  group_by(sampleID, patch, noise, treatment) %>%
  summarize(
    mean_CFU = mean(CFU_log_filled, na.rm = TRUE),
    sd_CFU = sd(CFU_log_filled, na.rm = TRUE),
    cv_CFU = sd_CFU / mean_CFU,
    .groups = "drop"
  )

cfu_stability$noise <- factor(cfu_stability$noise , levels = rev(c("constant", "white noise", "red noise")))
stab_plot <- ggplot(cfu_stability, aes(y = noise, x = cv_CFU)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(size = 3, alpha = 0.7, aes(color = patch))+
  facet_grid(~treatment) +
  scale_color_manual(values = mycolors) +
  labs(y = NULL, x = "Coefficient of variation" ) +
  mytheme

library(zoo)

pfu_padded_l_selected <- pfu_interp %>%
  dplyr::select(sampleID, day, PFU_sum_filled, PFU_log_filled)

pfu_cfu_df <- cfu_filled  %>%
  left_join(pfu_padded_l_selected, by = c("sampleID", "day"))


crosscorr_results <- pfu_cfu_df %>%
  group_by(sampleID, patch, noise, treatment) %>%
  summarize(
    corr_at_0 = cor(PFU_log_filled, CFU_log_filled, use = "complete.obs"),
    max_corr = max(ccf(PFU_log_filled, CFU_log_filled, lag.max = 10, plot = FALSE)$acf),
    lag_max_corr = ccf(PFU_log_filled, CFU_log_filled, lag.max = 10, plot = FALSE)$lag[which.max(ccf(PFU_log_filled, CFU_log_filled, lag.max = 10, plot = FALSE)$acf)],
    .groups = "drop"
  )

crosscorr_results$noise <- factor(crosscorr_results$noise , levels = rev(c("constant", "white noise", "red noise")))
cross <- ggplot(crosscorr_results, aes(y = noise, x = corr_at_0)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.7, aes(, color = patch)) +
  facet_grid(~treatment) +
  scale_color_manual(values = mycolors) +
  labs(x = "Cross-correlation (lag 0)", y = NULL) +
  mytheme


figure_Stab <- ggarrange(tte_plot,stab_plot,cross, ncol =1, labels = c("A","B","C"))


lm_corr <- lm(corr_at_0 ~ noise * treatment, data = crosscorr_results)
summary(lm_corr)
anova(lm_corr)

library(lme4)
lmm_corr <- lmer(corr_at_0 ~ noise * treatment + (1|patch), data = crosscorr_results)
summary(lmm_corr)
anova(lmm_corr)

library(lmPerm)
lm_corr_perm <- lmp(corr_at_0 ~ noise * treatment, data = crosscorr_results)
summary(lm_corr_perm)

library(lmerTest)

lmm_corr <- lmer(corr_at_0 ~ noise * treatment + (1|patch), data = crosscorr_results)
summary(lmm_corr)  # Now includes p-values
anova(lmm_corr, type = 3)  # Type 3 for interaction tests


