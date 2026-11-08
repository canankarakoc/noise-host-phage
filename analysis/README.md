# Analysis

Statistical analysis and figures for the noise-colour experiment. Run in order
from the repository folder. Each script writes to `results/` (tables) or
`figures/` (figures).

```r
# the data files go in data/ (see data/README.md); to read them from
# somewhere else, set DATA_DIR first, e.g. Sys.setenv(DATA_DIR = "~/Desktop/noise")
source("analysis/01_traits.R")   # endpoint traits that need no curve fitting
source("analysis/02_growth.R")   # growth parameters from the OD curves
source("analysis/03_models.R")   # noise x dispersal models for every response
source("analysis/04_verify.R")   # robustness checks
source("analysis/05_figures.R")  # manuscript and supplementary figures
```

To keep the model output, run the last three from a shell:

```sh
Rscript analysis/03_models.R > results/models_output.txt 2>&1
Rscript analysis/04_verify.R > results/verify_output.txt 2>&1
Rscript analysis/05_figures.R
```

Needs `dplyr`, `tidyr`, `readr`, `purrr`, `pracma`, `lme4`, `lmerTest`,
`pbkrtest`, `emmeans`, `car`, `ggplot2`, `patchwork`.

## Outputs

| Script | Writes |
|---|---|
| `01_traits.R` | `trait_selection.csv`, `trait_mic_plate.csv`, `trait_mic_stripe.csv`, `trait_phage_sensitivity.csv`, `ancestral_phage_sensitivity.csv`, `dynamics_long.csv`, `prep.rds` |
| `02_growth.R` | `trait_growth.csv`, `ancestral_growth.csv` |
| `03_models.R` | `model_anova.csv`, `model_emmeans.csv`, `model_contrasts.csv`, `paired_red_white.csv`, `microcosm_cfu_variability.csv`, `microcosm_host_phage_corr.csv` |
| `04_verify.R` | console only |
| `05_figures.R` | `figures/fig1_schedules` to `fig5_phage` and `figS1_growth` to `figS3_host_phage_corr`, each as PDF and 600 dpi PNG |

`model_emmeans.csv` holds marginal means on the model scale: log2 for MIC,
natural log for generation time. The figures read their means and 95%
confidence intervals from it, so they match the numbers in the text.

## What these do differently from the original `noise_main.R`

1. **Clones get a random intercept for microcosm.** Roughly eight clones come
   from each microcosm and are subsamples, not replicates; for MIC, 56% of the
   variance is between microcosms. The experimental unit is the microcosm,
   n = 36.
2. **Every response variable is modelled**, not just the cross-correlation.
3. **Type III tests use sum-to-zero contrasts**, so the noise and dispersal
   terms are main effects averaged over the other factor. With R's default
   treatment contrasts they would be simple effects at the reference level
   (noise within isolation, dispersal within the constant regime).
4. **No interpolated values enter the correlations.** The cross-correlation is
   computed only on days where both CFU and PFU were actually measured.
5. **One growth estimator for every curve** instead of selecting among four
   model families per clone by AIC. `02_growth.R` uses a rolling log-linear
   regression and reports logistic and Gompertz fits as a robustness check.
6. **Time to phage extinction is not computed.** Phage was detected on the last
   recorded day in all 36 microcosms, and recording effort is confounded with
   treatment (correlation 0.68 between last recorded day and number of
   assays), so the quantity is not estimable from these files. `01_traits.R`
   writes `dynamics_long.csv` with `pfu_measured` and `phage_absent` kept
   separate so this stays visible.
7. **Patch is used as a blocking factor only where it means the same thing**:
   the red-versus-white contrast, where red patch k is a reordering of the same
   30 concentrations as white patch k.

`03_models.R` reports a primary analysis (day 30, with day-24 clones standing
in for the four contaminated microcosms m15-m18) and a pooled-day sensitivity
check. `04_verify.R` adds the patch-blocked red-versus-white tests, a
day-30-only refit, a dispersal-only subset with no contamination, and effect
sizes.
