# Why Does Wearing High-Heeled Footwear Stiffen the Achilles Tendon?

Code and data for the WCB 2026 poster (P2-212) by Andrew M. Thornton,
Jordyn N. Schroeder, Gregory S. Sawicki and Owen N. Beck
(Human Locomotion Lab, The University of Texas at Austin; Georgia Institute
of Technology).

> **Status:** work in progress. Items marked **[PLACEHOLDER]** will be filled in
> when the project is finished.

## Study overview

**Research question:** does Achilles tendon (AT) strain explain the increase
in AT stiffness caused by wearing high heels?

- **Participants:** 8 adults who were not habitual high-heel wearers. Each
  received custom high heels that set the ankle at 14° plantar flexion, for
  about 14 weeks.
  - **Users** (n = 4) wore the heels for more than 1,500 steps/day.
  - **Nonusers** (n = 4) wore them for fewer than 1,000 steps/day and act as
    the control group.
- **Visits:** Pre and Post intervention.
  - **Tendon stiffness** was measured with ultrasonography and dynamometry.
  - **Walking biomechanics** were recorded during treadmill walking at
    1.3 m/s in flat shoes and in the high heels. AT force and strain were
    estimated from the ankle moment, a subject-scaled AT moment arm, and the
    tendon's stiffness and length.

The statistics address three questions:

1. Did AT stiffness change from Pre to Post, and did the change differ
   between Users and Nonusers?
2. Did walking in heels change peak ground reaction force (GRF), ankle
   effective mechanical advantage (EMA), AT force and AT strain compared with
   flats?
3. Did the change in stiffness scale with daily steps in heels
   (dose-response)?

Methods, exclusions and limitations are documented in the code comments,
especially the header of [`R/WCB2026_Stats.R`](R/WCB2026_Stats.R).

## Repository layout

```
GA_Tech_Heels/
├── README.md
├── DATA.md                      data dictionary (files, columns, units, exclusions)
├── LICENSE                      [PLACEHOLDER]
├── CITATION.cff
├── GA_Tech_Heels.Rproj          open in RStudio so R paths resolve from the repo root
├── matlab/
│   ├── analyze_participant.m        summary values for one participant / visit / condition
│   ├── run_batch_analysis.m         runs analyze_participant.m for everyone -> batch_results_*.csv
│   ├── export_stance_curves.m       mean stance-phase time series for one participant / visit / condition
│   └── run_batch_stance_curves.m    runs export_stance_curves.m for everyone -> stance_curves.csv
├── R/
│   ├── build_participant_data_sheet.R   batch results -> analysis data sheet
│   └── WCB2026_Stats.R                  statistics, figures and tables
├── data/
│   ├── raw/          raw .mat files (not tracked; see "Getting the raw data")
│   ├── params/       per-visit tendon parameters used by the MATLAB pipeline
│   ├── processed/    MATLAB output
│   └── analysis/     R input
├── output/           figures and tables written by WCB2026_Stats.R (not tracked)
└── docs/
    └── WCB2026_poster.pdf
```

## Pipeline

Each step reads the previous step's output. All paths are relative to the
repository root, so the scripts run from a fresh clone without editing.

| Step | Run | Reads | Writes |
|---|---|---|---|
| 1 | `matlab/run_batch_analysis.m` | `data/raw/*.mat`, `data/params/participant_params_Day{1,3}.csv` | `data/processed/batch_results_{heels,flats}_{Day1,Day3}.csv` |
| 2 | `matlab/run_batch_stance_curves.m` | same as step 1 | `data/processed/stance_curves.csv` |
| 3 | `R/build_participant_data_sheet.R` | `data/processed/batch_results_*.csv`, existing `data/analysis/participant_data_sheet_R.csv` | `data/analysis/participant_data_sheet_R.csv`, `data/analysis/heels_flats_comp.csv` |
| 4 | `R/WCB2026_Stats.R` | `data/analysis/*.csv`, `data/processed/stance_curves.csv` | `output/figures/*.pdf`, `output/tables/baseline_characteristics.csv` |

- **Running only the statistics:** the processed and analysis files are
  included in the repository, so step 4 runs on its own without the raw data
  or MATLAB.
- **What step 3 keeps from the sheet:** it rebuilds only the heels / flats
  columns from the batch results. Group, step counts, anthropometrics and
  tendon properties are kept from the existing sheet.

### Running the MATLAB steps

1. Download the raw `.mat` files into `data/raw/` (see below).
2. In MATLAB, from any folder:
   ```matlab
   run("path/to/GA_Tech_Heels/matlab/run_batch_analysis.m")
   run("path/to/GA_Tech_Heels/matlab/run_batch_stance_curves.m")
   ```
   One run of each script covers both visits and both footwear conditions.
   Participants or visits with no data are skipped with a warning and left as
   `NaN` (expected: BO11 Pre, DF07 Pre flats, DF07 Post).

### Running the R steps

Open `GA_Tech_Heels.Rproj` in RStudio and source the scripts, or from the
repository root:

```bash
Rscript R/build_participant_data_sheet.R
```

```bash
Rscript R/WCB2026_Stats.R
```

`build_participant_data_sheet.R` prints every value it changes and saves a
dated backup of the files it overwrites (backups are ignored by git).

## Getting the raw data

The raw motion-capture and force data are too large for GitHub
(`MECH_Data.mat` is about 740 MB). They are archived at:

**[PLACEHOLDER: data archive link / DOI]**

Put these three files in `data/raw/`:

- `MECH_Data.mat`: ground reaction forces
- `JointMomentData.mat`: ankle moments
- `JointAngleData.mat`: ankle angles

Their structure is described in [`DATA.md`](DATA.md).

## Requirements

Versions used for the poster analysis:

- **MATLAB** R2021b (9.11) with the Signal Processing Toolbox 8.7 (`butter`,
  `filtfilt`, `findpeaks`)
- **R** 4.5.2 with:

  | Package | Version |
  |---|---|
  | tidyverse | 2.0.0 (dplyr 1.1.4, ggplot2 4.0.2, tidyr 1.3.1, purrr 1.2.0, readr 2.1.6) |
  | mosaic | 1.9.2 |
  | cowplot | 1.2.0 |
  | scales | 1.4.0 |

## Notes on the data

These are described in more detail in [`DATA.md`](DATA.md) and the script
comments.

- **AT length and cross-sectional area are not analyzed:** those
  measurements were not reliable enough. AT length is still used as the
  tendon slack length `Lo` in the strain calculation. The heels-vs-flats
  comparison is unaffected because both conditions share the same `Lo`, but
  absolute strain values and comparisons between visits inherit its error.
- **GJ03's Pre treadmill visit is excluded** because stance detection failed
  on that trial.
- **SE15's Post force-plate data were recorded at 2,000 Hz** instead of
  1,000 Hz. The MATLAB code detects this and downsamples before aligning the
  forces with the 1,000 Hz ankle moment data.
- **Stance only:** all per-trial summary values (peaks, means and strain
  impulse) are calculated within each stance, from heel strike to toe-off,
  then averaged over the 12–15 stances in the 15 s trial.

## Citation

If you use this code or data, please cite the poster (see
[`CITATION.cff`](CITATION.cff)):

> Thornton AM, Schroeder JN, Sawicki GS, Beck ON. Why does wearing high heeled
> footwear stiffen the Achilles tendon? World Congress of Biomechanics, 2026.
> **[PLACEHOLDER: full citation / DOI]**

## License

**[PLACEHOLDER]**: see [`LICENSE`](LICENSE).

## Contact

Andrew M. Thornton, Department of Kinesiology and Health Education, The
University of Texas at Austin. Email: amt4782@my.utexas.edu

## Acknowledgements

We thank Kinsey Herrin for helping fabricate the heel lifts.
