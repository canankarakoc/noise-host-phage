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

It also prints the six constant concentrations.

`noise_main.R`, from top to bottom:

1. **Growth**: blanks the plate-reader OD, fits Gompertz, Baranyi, logistic and Weibull growth models to each clone (`growthrates`), keeps the best model by AIC, drops poor fits (RMSE > 0.1), and compares generation time, carrying capacity and lag with the ancestor. One plot per fitted curve (about 400) is saved to `figures/`.
2. **Phage sensitivity**: the difference in area under the growth curve without and with phage.
3. **Competitive fitness**: the selection coefficient of evolved clones against REL606.
4. **Kanamycin resistance**: MICs from the plate dilution assay (with quality filters) and from gradient test strips.
5. **Population dynamics**: bacterial and phage densities over time, time to phage extinction, variability of bacterial density (CV), and the correlation between phage and bacterial densities, tested with linear, mixed-effects and permutation models.

Four day-30 populations (white noise, isolation, patches 3-6) were contaminated, so their day-24 clones are used instead. Objects ending in `_ALL` pool days 24 and 30.

## How to run

Last run in October 2026 with R 4.6.0, growthrates 0.8.5, tidyverse 2.0.0, ggpubr 1.0.0, pracma 2.4.6, zoo 1.8.15, lme4 2.0.1, lmerTest 3.2.1, lmPerm 2.1.6 and truncnorm 1.0.9.

```r
install.packages(c("growthrates", "tidyverse", "ggpubr", "pracma", "zoo", "lme4", "lmerTest", "lmPerm", "truncnorm"))
```

1. `scenarios.R` needs no data. Set the working directory to the repository folder and run it.
2. For the analysis, put the data files in `data/` (see Data availability). `data/README.md` describes them.
3. Run `noise_main.R` from the repository folder, top to bottom. It takes a few minutes. The composite figures are `figure3` (growth and phage sensitivity) and `figure_Stab` (dynamics); the fitness and MIC panels are `selection`, `micPlot` and `micPlot_stripe`.

## Data availability

The data are not included while the work is unpublished. They are available from me on request.

## License

The code and the schedules are released under the MIT License (see `LICENSE`).
