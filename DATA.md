# Data dictionary

This file describes every data file in the repository, plus the raw `.mat`
files archived separately. It covers their structure, columns, units, and the
known exclusions and quirks.

**Labels used throughout**

| Label | Meaning |
|---|---|
| `DayA_01` / `Day1` / `Pre` | Pre-intervention visit |
| `DayA_03` / `Day3` / `Post` | Post-intervention visit (about 14 weeks later) |
| `HEEL` / `heels` / `Heels` | Walking in the custom high heels (ankle at 14° plantar flexion) |
| `FLAT` / `flats` / `Flats` | Walking in flat shoes |
| `S13` | Treadmill trial at 1.3 m/s (the only speed analyzed) |

**Participants** in this study, by group:

| Group | Participants |
|---|---|
| User (> 1,500 steps/day in heels) | BO11, DF07, OJ08, SS11 |
| Nonuser (< 1,000 steps/day in heels) | GJ03, GP04, RL10, SE15 |

The raw `.mat` files also contain participants from a separate toe-walking
study (AA17, CP12, JS16, LQ17, SB06, SG01, TL05) and group averages
(`AllHeel`, `AllToe`). These are not used here.

## Missing data and exclusions

| Participant | Visit | Condition | Status |
|---|---|---|---|
| BO11 | Pre | both | No treadmill data recorded |
| DF07 | Pre | flats | No data recorded (heels only) |
| DF07 | Post | both | No treadmill data recorded |
| GJ03 | Pre | both | **Excluded**: stance detection failed (about 40 "stances" detected in 15 s). Set to `NA` in the analysis files by `build_participant_data_sheet.R` |
| OJ08, SE15 | n/a | n/a | No total step count |

**Other known issues**

- **SE15 Post (`DayA_03`) force-plate data were recorded at 2,000 Hz** (30,000
  samples for 15 s), while every other trial and all moment and angle data
  are at 1,000 Hz. `analyze_participant.m` and `export_stance_curves.m` detect
  this, low-pass filter at the native rate, and downsample to 1,000 Hz before
  aligning with the moment data.
- **AT length (`ATLength`, used as `Lo`) and cross-sectional area (`CSA`) are
  not reliable enough to analyze.** For example, length changes of up to 20%
  in about 100 days, recorded in 2.5–5 mm steps. They remain in the data
  files because `Lo` is part of the strain calculation and `CSA` of `E` (see
  below).

---

## `data/raw/`: raw MATLAB data (not tracked; see README)

Each file holds one struct with every participant as a top-level field.

| File | Path to the data used | Contents |
|---|---|---|
| `MECH_Data.mat` | `MECH_Data.(pid).(day).(cond).S13.allData.Force_Fy1` / `Force_Fz1` | Anterior–posterior (`Fy1`) and vertical (`Fz1`) GRF, N. The sign is negated in the code. 1,000 Hz, at least 15 s (SE15 Post: 2,000 Hz) |
| `JointMomentData.mat` | `JointMomentData.(pid).(day).(cond).AnkMom` | Ankle moment, Nm/kg, 1,000 Hz, 15,000 samples, aligned with the GRF |
| `JointAngleData.mat` | `JointAngleData.(pid).(day).(cond).AnkAng` | Ankle angle, degrees, 1,000 Hz, 15,000 samples |

`MECH_Data` trials also contain marker-based events (`heel_strikes`,
`gaitCycles`, ...). These are not used; stance is detected from the vertical
GRF with a 10 N threshold.

## `data/params/`: tendon parameters for the MATLAB pipeline

`participant_params_Day1.csv` (Pre) and `participant_params_Day3.csv` (Post)
have one row per participant. `participant_params_template.csv` shows the
expected format.

| Column | Units | Description |
|---|---|---|
| `ParticipantID` | n/a | Must match the struct field name in the `.mat` files |
| `lin_k` | N/mm | Linear AT stiffness at that visit |
| `rest_AtMA` | m | Resting AT moment arm (used to scale a generic moment arm–angle curve) |
| `CSA` | mm² | AT cross-sectional area (unreliable; used only for `E`) |
| `Lo` | mm | AT slack length (unreliable; used in the strain calculation) |

## `data/processed/`: MATLAB output

### `batch_results_{heels,flats}_{Day1,Day3}.csv` (from `run_batch_analysis.m`)

One row per participant. Most values are averaged over every stance detected
in the 15 s trial. `NaN` means no data.

| Column | Units | Description |
|---|---|---|
| `ParticipantID` | n/a | |
| `E` | N/mm² per % strain | Tendon stress ÷ strain at peak AT force, averaged over stances. Strain is in %, so multiply by 100 for MPa. Uses `CSA` (unreliable) |
| `mean_peak_Fmtu_EMA` | unitless | Effective mechanical advantage (AT moment arm ÷ GRF moment arm) at the instant of peak AT force, averaged over stances |
| `mean_pushoff_GRF_EMA` | unitless | EMA at the push-off (second) resultant GRF peak, the highest GRF peak in the second half of each stance (~78% of stance), averaged over stances. The largest GRF peak (`mean_peak_Fr`) is the loading peak (~25%) in most trials but the push-off peak in others, so it is not used as the EMA event |
| `mean_lin_strain_impulse` | % strain · s | Integral of AT strain over time from heel strike to toe-off of each stance, averaged over stances |
| `mean_peak_lin_strain` | % | Peak AT strain per stance, averaged over stances |
| `mean_lin_strain` | % | Mean AT strain within each stance (heel strike to toe-off), averaged over stances |
| `mean_peak_Fr` | N | Peak resultant GRF per stance, averaged over stances |
| `mean_peak_Fmtu` | N | Peak AT (muscle–tendon unit) force per stance, averaged over stances |
| `mean_Fmtu` | N | Mean AT force within each stance, averaged over stances |
| `mean_peak_ank_mom` | Nm/kg | Peak ankle moment per stance, averaged over stances |
| `mean_ank_mom` | Nm/kg | Mean (unfiltered) ankle moment within each stance, averaged over stances |

**How the quantities are calculated** (`analyze_participant.m`):

- **AT force:** `Fmtu = ankle moment / AT moment arm × body mass`. The AT
  moment arm depends on ankle angle; it uses the Maganaris resting values
  scaled to `rest_AtMA`.
- **Strain:** `strain = Fmtu / lin_k / Lo × 100`.
- **Peaks:** for each stance, the largest peak found by `findpeaks` with a
  minimum height of max/1.5.
- **Body mass:** estimated from the mean vertical GRF over all strides.
- **Filtering:** GRF and moments are low-pass filtered at 10 Hz (4th-order
  Butterworth, zero lag).
- **Stance only:** every value is calculated within each stance (heel strike
  to toe-off), or at an event within it (peak AT force, push-off GRF peak), and
  then averaged over stances; swing-phase samples are never included. Before 2026-10-08 the three `mean_*` values
  averaged every non-zero sample in the whole trial, including swing. The
  stance-only means are about 1.5 times larger.

### `stance_curves.csv` (from `run_batch_stance_curves.m`)

Mean time series over stance in long format: 101 rows per participant ×
visit × condition. Each stance is time-normalized to 0–100% before averaging.

| Column | Units | Description |
|---|---|---|
| `ParticipantID` | n/a | |
| `Visit` | n/a | `Pre` or `Post` |
| `Condition` | n/a | `HEEL` or `FLAT` |
| `stance_pct` | % | 0 = heel strike, 100 = toe-off |
| `strain_mean` | % | AT strain, mean over stances |
| `strain_sd` | % | SD of AT strain across stances |
| `grf_bw` | body weights | Resultant GRF ÷ body weight, mean over stances |
| `fmtu_bw` | body weights | AT force ÷ body weight, mean over stances. Slightly negative values early in stance reflect a net dorsiflexor ankle moment |
| `r_internal` | m | AT (internal) moment arm, mean over stances |
| `R_external` | m | GRF (external) moment arm = ankle moment ÷ resultant GRF, mean over stances |
| `ema_median` | unitless | Median over stances of `r_internal / R_external` |
| `n_stances` | n/a | Number of stances averaged |

EMA is undefined near heel strike and toe-off, where the ankle moment and
`R_external` pass through zero. `WCB2026_Stats.R` therefore computes group EMA
as mean `r_internal` ÷ mean `R_external` and shows it only where every
participant's `R_external` is above 2 cm (22–96% of stance).

## `data/analysis/`: input to `WCB2026_Stats.R`

### `participant_data_sheet_R.csv`

One row per participant × visit (16 rows). The heels / flats columns are
rebuilt from the batch results by `build_participant_data_sheet.R`; the other
columns are entered by hand.

| Column | Units | Description |
|---|---|---|
| `ParticipantID` | n/a | |
| `Time` | n/a | `Pre` or `Post` |
| `Compliance` | n/a | Group: `User` or `Nonuser` |
| `InterventionDays` | days | Length of the intervention period |
| `HTDSteps` | steps/day | Average daily steps in heels over the intervention (same value on both rows) |
| `TotalSteps` | steps/day | Average total daily steps over the intervention (same value on both rows) |
| `Height` | cm | |
| `Age` | years | |
| `Mass` | kg | Measured at that visit |
| `Sex` | n/a | `M` or `F` |
| `ATMomentArm` | m | Resting AT moment arm |
| `ATLength` | mm | AT slack length (**unreliable, not analyzed**) |
| `k_lin` | N/mm | Linear AT stiffness |
| `CSA` | mm² | AT cross-sectional area (**unreliable, not analyzed**) |
| `{heels,flats}_EMA` | unitless | = `mean_peak_Fmtu_EMA` |
| `{heels,flats}_EMA_pushoff` | unitless | = `mean_pushoff_GRF_EMA` |
| `{heels,flats}_strain_impulse` | % strain · s | = `mean_lin_strain_impulse` |
| `{heels,flats}_peak_strain` | % | = `mean_peak_lin_strain` |
| `{heels,flats}_mean_strain` | % | = `mean_lin_strain` (stance only) |
| `{heels,flats}_peak_Fmtu` | N | = `mean_peak_Fmtu` |
| `{heels,flats}_mean_Fmtu` | N | = `mean_Fmtu` |
| `{heels,flats}_peak_ank_mom` | Nm/kg | = `mean_peak_ank_mom` |
| `{heels,flats}_mean_ank_mom` | Nm/kg | = `mean_ank_mom` |
| `{heels,flats}_peak_Fr` | N | = `mean_peak_Fr` |

### `heels_flats_comp.csv`

A long-format copy of the main footwear variables, used for the heels vs flats
figures. It has one row per participant × condition × visit and is generated
by `build_participant_data_sheet.R`.

| Column | Units | Description |
|---|---|---|
| `ParticipantID`, `Condition` (`Flats` / `Heels`), `Time` | n/a | |
| `EMA`, `EMA_pushoff`, `strain_impulse`, `peak_strain`, `mean_strain` | as above | |
| `peak_Fr` | body weights | Peak resultant GRF ÷ (`Mass` × 9.81) |
| `peak_Fmtu` | body weights | Peak AT force ÷ (`Mass` × 9.81) |

## `output/`: generated by `WCB2026_Stats.R` (not tracked)

- `output/tables/baseline_characteristics.csv`: baseline characteristics for
  each participant, with mean (SD) by group.
- `output/figures/*.pdf`: poster figures. The script comments describe the
  contents of each figure and the caveats to put in its legend.
  `flats_vs_heels_stance_summary_fig.pdf` combines, for GRF, EMA, AT force and
  AT strain, the mean stance time series (left) with the heels vs flats box
  plot of the matching per-trial value (right).
