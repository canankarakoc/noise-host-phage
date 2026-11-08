# noise-host-phage

Experimental evolution of *Escherichia coli* under phage predation and fluctuating antibiotic stress, where the fluctuations differ in their temporal structure ("noise colour").

## Question

Natural environments fluctuate, and the fluctuations have structure. In white noise each day is independent of the last; in red noise consecutive days are similar, so conditions are more predictable. Does the colour of antibiotic fluctuations change how bacteria adapt, and does it change their interaction with a lytic phage? Does dispersal between populations alter the outcome?

## Design

*E. coli* B REL607 was evolved together with phage T4 for 30 days, with daily 10% transfers into 6 ml LB at 30 °C. Kanamycin was added every day under one of three regimes:

- **Constant**: six fixed concentrations, one per patch, spanning the range of the noise series (their minimum, quartiles, median, mean and maximum).
- **White noise**: independent daily concentrations from a truncated normal distribution (mean 2 µg/ml, sd 1 µg/ml).
- **Red noise**: each patch gets the same daily values as its white-noise counterpart, reordered so that consecutive days are correlated (lag-1 autocorrelation 0.8).

Each regime had six patches, either isolated or connected by dispersal (0.1% migration among patches), giving 36 microcosms. Bacteria and phage were counted every second day. Clones isolated on days 24 and 30 were assayed for growth with and without phage, kanamycin resistance (MIC) and competitive fitness against the marked ancestor REL606.

## My role

I designed and planned the experiments, generated the antibiotic schedules, and coded and ran all the analyses. Students helped with the lab work (2020-2021).

## Status

Manuscript in preparation; not yet published. The citation will be added here when it is.

## What the scripts do

`scenarios.R` generates the kanamycin schedules. It is seeded, so it reproduces the schedules used in the experiment, and it writes them to `schedules/`:

- `white_noise_mean2_sd1.txt` and `fixed_corr_red_mean2_sd1_corr0.8.txt`: the white- and red-noise series (30 days × 6 patches, µg/ml).
- `white.kan.mean2.sd1.txt` and `red.kan.mean2.sd1.c0.8.txt`: the same series as the daily volume (µl) of 1 mg/ml kanamycin stock added to each 6 ml culture.
- `redified_mean2_sd1_weight2.txt` and `red.kan.mean2.sd1.w2.txt`: an alternative red-noise series (fixed weighting instead of fixed autocorrelation) that was generated but not used.

It also prints the six constant concentrations. The constant patches C1 to C6 received the minimum, first quartile, median, mean, third quartile and maximum, in that order.

`analysis/` holds the analysis, five scripts run in order (`analysis/README.md` has the details):

1. **`01_traits.R`**: competitive fitness (selection coefficient against REL606), kanamycin MICs from the plate assay (with quality filters; clones not inhibited at 256 µg/ml are kept as right-censored) and from gradient strips, phage sensitivity (difference in area under the growth curve without and with phage), and bacterial and phage densities on the days they were measured.
2. **`02_growth.R`**: generation time and carrying capacity from the blanked OD curves, with one model-free estimator for every curve and logistic and Gompertz fits as a check.
3. **`03_models.R`**: a noise × dispersal model for every response, with a random intercept for microcosm where clones are the unit, type III F tests (Kenward-Roger), marginal means and Tukey contrasts within each dispersal treatment, and red versus white paired by patch.
4. **`04_verify.R`**: robustness checks: red versus white blocked by patch, day-30-only refits, the uncontaminated dispersal subset, and effect sizes.
5. **`05_figures.R`**: the manuscript and supplementary figures, written to `figures/`.

Four day-30 populations (white noise, isolation, patches 3-6; m15-m18) were contaminated. The primary analysis uses day 30 with their day-24 clones standing in, and `04_verify.R` repeats the key tests with day-30 data only.

## How to run

Last run in November 2026 with R 4.6.0, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, purrr 1.2.2, pracma 2.4.6, lme4 2.0.1, lmerTest 3.2.1, pbkrtest 0.5.5, emmeans 2.0.3, car 3.1.5, ggplot2 4.0.3 and patchwork 1.3.2. `scenarios.R` also needs tidyverse and truncnorm.

```r
install.packages(c("dplyr", "tidyr", "readr", "purrr", "pracma", "lme4", "lmerTest",
                   "pbkrtest", "emmeans", "car", "ggplot2", "patchwork",
                   "tidyverse", "truncnorm"))
```

1. `scenarios.R` needs no data. Set the working directory to the repository folder and run it.
2. For the analysis, put the data files in `data/` (see Data availability). `data/README.md` describes them.
3. From the repository folder, run `analysis/01_traits.R` to `analysis/05_figures.R` in order. Together they take under a minute. Tables go to `results/` and figures to `figures/`.

The original single-script analysis (`noise_main.R`) was replaced by `analysis/` in November 2026. It is still in the git history.

## Data availability

The data are not included while the work is unpublished. They are available from me on request.

## License

The code and the schedules are released under the MIT License (see `LICENSE`).
