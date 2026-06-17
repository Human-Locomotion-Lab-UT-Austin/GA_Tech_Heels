# GA_Tech_Heels

WCB2026 Folder:
  Code:
  
| File | Purpose |
|---|---|
| `analyze_participant.m` | Function. Runs the full pipeline (GRF filtering, stance detection, body weight, ankle moment, internal/external moment arm, Fmtu, strain, EMA, Young's modulus) for **one** participant and returns the 6 summary values. |
| `run_batch_analysis.m` | Script. Loads the data, reads your parameter spreadsheet, calls `analyze_participant.m` for every participant, and saves a results table. |
| `participant_params_X_X.csv` | Spreadsheet containing necessary participant information | (Found in Data folder)

## Requirements

- MATLAB with the Signal Processing Toolbox (`butter`, `filtfilt`, `findpeaks`).
- Three `.mat` files, each containing **every participant** as a top-level
  field (e.g. `MECH_Data.SS11`, `MECH_Data.AA17`, ...):
  - `MECH_Data.mat` — GRF data, with path `(pid).`(day)`.`(cond)`.`(trial)`.allData.Force_Fy1` / `Force_Fz1`
  - `JointMomentData.mat` — ankle moment, with path `(pid).`(day)`.`(cond)`.AnkMom`
  - `JointAngleData.mat` — ankle angle, with path `(pid).DayA_01.`(cond)`.AnkAng`
- A parameter spreadsheet (`.csv` or `.xlsx`) with one row per participant.

> Edit the `day`/`cond`/`trial` variables near the top of
> the function based on which trials are being analyzed.
> `day` - DayA_01 (Pre-intervention), DayA_03 (Post-intervention)
> `cond` - Heels, Flats
> `trial` - 0.5, 0.9, 1.3, 1.7 (treadmill speed, m/s)

## Setting up the parameter spreadsheet

Copy `participant_params_template.csv`, rename it, and add one row per
participant with these exact column headers:

| Column | Description |
|---|---|
| `ParticipantID` | Must exactly match the struct field name in the `.mat` files (e.g. `SS11`) |
| `lin_k` | Linear tendon stiffness (N/mm) |
| `rest_AtMA` | Resting Achilles tendon moment arm (m) |
| `CSA` | Tendon cross-sectional area (mm²) |
| `Lo` | Tendon slack length (mm) |

You can keep it as `.csv` or save as `.xlsx` — `readtable` handles both, just
update the `param_file` path in `run_batch_analysis.m`.

## Running it

1. Open `run_batch_analysis.m`.
2. Edit the two paths at the top:
   ```matlab
   data_dir   = "/path/to/your/Data";
   param_file = fullfile(data_dir, "participant_params.csv");
   ```
3. Make sure `analyze_participant.m` is on your MATLAB path (same folder as
   `run_batch_analysis.m` is easiest).
4. Run `run_batch_analysis.m`.

## Output

A MATLAB `table` (`results_table`) with one row per participant:

`ParticipantID | E | mean_peak_Fmtu_EMA | mean_lin_strain_impulse | mean_peak_lin_strain | mean_lin_strain | mean_peak_Fr`

It's displayed in the Command Window and saved to `batch_results.csv` in
`data_dir`.

## Error handling

If a participant's data is missing, malformed, or has no detectable stance
phases, `run_batch_analysis.m` catches the error, prints a warning naming the
participant, and fills that row with `NaN` instead of stopping the batch.
Check the Command Window warnings after a run to see which participants (if
any) failed and why.
