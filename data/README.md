# Data files

The data are not included in this repository while the work is unpublished (see the main README). `noise_main.R` expects these files here. All are comma-separated.

**`design.csv`**: the 36 microcosms. Columns: `microcosmID` (m01-m36), `sampleID` (for example C1I = constant, patch 1, isolation), `patch` (p1-p6), `noise` (constant, white noise, red noise), `treatment` (isolation, dispersal).

**`files/01.csv` to `files/22.csv`**: plate-reader growth curves of the evolved clones, one file per 96-well plate. Columns: `time` (minutes), then one column per well (`A01`-`H12`, OD600). The file number is the plate number (`set` in `plate_design.csv`), so keep the names.

**`plate_design.csv`**: what is in each well of the 22 plates. Columns: `well`, `microcosmID` (or `blank`; 606 is the ancestral control), `sampleID`, `day` (isolation day of the clone), `phage_treatment` (with_phage, without_phage), `clones`, `set` (plate number).

**`competition_alldays.csv`**: competitions of evolved clones against the ancestor REL606. Columns: `time` (t0 or t1, 24 h), `microcosmID`, `sampleID`, `day`, `dilution_fold`, `volume_plated` (µl), `red_colonies` (REL606), `pink_colonies` (evolved).

**`mic_od_combined.csv`**: kanamycin MICs by dilution in 96-well plates. Columns: `well`, `MIC_0` to `MIC_256` (OD600 after 24 h at each kanamycin concentration, µg/ml), `negative control`, `stamp_control`, `sampleID`, `plate`, `day`.

**`MIC_stripe_test.csv`**: kanamycin MICs from gradient test strips, technical replicates in separate rows. Columns: `microcosmID`, `sampleID`, `MIC` (µg/ml), `Day`.

**`dynamicsCFU.csv`**: bacterial counts during the experiment, every second day. Columns: `sample_full`, `sampleID`, `day`, `dilution`, `CFU` (colonies in a 10 µl spot).

**`dynamicsPFU.csv`**: phage counts during the experiment, every second day. Columns: `day`, `microcosmID`, `sampleID`, `dilution`, `volume` (µl), `count` (plaques).
