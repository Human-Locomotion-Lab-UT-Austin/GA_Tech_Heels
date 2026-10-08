# =============================================================================
# build_participant_data_sheet.R
#
# Rebuilds the heels / flats columns of participant_data_sheet_R.csv from the
# four MATLAB batch result files (run_batch_analysis.m + analyze_participant.m),
# so the data sheet used by WCB2026_Stats.R is reproducible from one pipeline
# run. Also regenerates heels_flats_comp.csv (the long-format copy used for the
# heels vs flats figures) from the rebuilt sheet.
#
#   batch_results_{heels,flats}_Day1.csv -> Time == "Pre"
#   batch_results_{heels,flats}_Day3.csv -> Time == "Post"
#
# All other columns (group, step counts, anthropometrics, k_lin, ATLength,
# CSA, ...) are kept from the existing data sheet unchanged.
#
# Exclusions: visits whose GRF data are unreliable are set to NA for every
# heels / flats variable (see `exclusions` below).
#
# Existing files are copied to dated backups before they are overwritten, and
# a report of every changed value is printed. After checking the report, push
# both updated files to the GitHub repository that WCB2026_Stats.R reads from.
# =============================================================================

library(tidyverse)

data_dir   <- "/Users/andrewthornton/Documents/WCB2026/Data"
sheet_file <- file.path(data_dir, "participant_data_sheet_R.csv")
comp_file  <- file.path(data_dir, "heels_flats_comp.csv")
g <- 9.81 # m/s^2, converts mass to body weight for peak GRF in heels_flats_comp

# Data sheet variable (after the heels_ / flats_ prefix) <- batch file column
batch_columns <- c(
  EMA            = "mean_peak_Fmtu_EMA",      # EMA at peak AT force
  strain_impulse = "mean_lin_strain_impulse", # % strain * s
  peak_strain    = "mean_peak_lin_strain",    # %
  mean_strain    = "mean_lin_strain",         # %
  peak_Fmtu      = "mean_peak_Fmtu",          # N
  mean_Fmtu      = "mean_Fmtu",               # N
  peak_ank_mom   = "mean_peak_ank_mom",       # Nm/kg
  mean_ank_mom   = "mean_ank_mom",            # Nm/kg
  peak_Fr        = "mean_peak_Fr"             # N
)

# Visits excluded from every heels / flats variable, with the reason
exclusions <- tribble(
  ~ParticipantID, ~Time, ~reason,
  "GJ03",         "Pre", "GRF stance detection failed (about 40 'stances' detected in 15 s, peak strain ~1%)"
)


# 1. Read the batch files -------------------------------------------------------

batch_files <- expand_grid(condition = c("heels", "flats"), Time = c("Pre", "Post")) |>
  mutate(file = file.path(data_dir, sprintf("batch_results_%s_%s.csv", condition,
                                            if_else(Time == "Pre", "Day1", "Day3"))))

batch_long <- batch_files |>
  pmap_dfr(function(condition, Time, file) {
    read_csv(file, show_col_types = FALSE) |>
      select(ParticipantID, all_of(unname(batch_columns))) |>
      mutate(condition = condition, Time = Time)
  }) |>
  pivot_longer(all_of(unname(batch_columns)), names_to = "batch_column") |>
  mutate(
    variable = names(batch_columns)[match(batch_column, batch_columns)],
    value    = if_else(is.nan(value), NA_real_, value)
  )

# Apply exclusions
batch_long <- batch_long |>
  left_join(exclusions, by = c("ParticipantID", "Time")) |>
  mutate(value = if_else(is.na(reason), value, NA_real_)) |>
  select(-reason)


# 2. Rebuild the sheet ----------------------------------------------------------

old_sheet <- read_csv(sheet_file, show_col_types = FALSE)

footwear_wide <- batch_long |>
  mutate(column = paste(condition, variable, sep = "_")) |>
  select(ParticipantID, Time, column, value) |>
  pivot_wider(names_from = column, values_from = value)

footwear_cols <- names(old_sheet)[str_detect(names(old_sheet), "^(heels|flats)_")]
stopifnot(setequal(footwear_cols, setdiff(names(footwear_wide), c("ParticipantID", "Time"))))

new_sheet <- old_sheet |>
  select(-all_of(footwear_cols)) |>
  left_join(footwear_wide, by = c("ParticipantID", "Time")) |>
  select(all_of(names(old_sheet))) # keep the original column order

stopifnot(nrow(new_sheet) == nrow(old_sheet),
          identical(new_sheet$ParticipantID, old_sheet$ParticipantID),
          identical(new_sheet$Time, old_sheet$Time))


# 3. Report changes -------------------------------------------------------------

changes <- list(old = old_sheet, new = new_sheet) |>
  map(~ .x |>
        select(ParticipantID, Time, all_of(footwear_cols)) |>
        pivot_longer(all_of(footwear_cols), names_to = "column")) |>
  reduce(full_join, by = c("ParticipantID", "Time", "column"), suffix = c("_old", "_new")) |>
  filter(xor(is.na(value_old), is.na(value_new)) |
           abs(value_new - value_old) > 1e-9 * pmax(1, abs(value_old))) |>
  mutate(pct_change = (value_new - value_old) / abs(value_old) * 100)

cat("Changed values:", nrow(changes), "\n")
print(changes, n = Inf, width = Inf)


# 4. Back up the old sheet and write the new one --------------------------------

backup_file <- file.path(data_dir, sprintf("participant_data_sheet_R_backup_%s.csv", Sys.Date()))
if (!file.exists(backup_file)) invisible(file.copy(sheet_file, backup_file))
write_csv(new_sheet, sheet_file, na = "NA")
cat("Wrote", sheet_file, "\nBackup of the previous sheet:", backup_file, "\n")


# 5. Regenerate heels_flats_comp.csv ---------------------------------------------

# One row per participant x condition x visit, in data sheet row order with
# flats before heels; peak GRF expressed in body weights. This reproduces the
# previous heels_flats_comp.csv exactly when run on the previous sheet.
comp_vars <- c("EMA", "strain_impulse", "peak_strain", "mean_strain", "peak_Fr")

heels_flats_comp <- new_sheet |>
  select(ParticipantID, Time, Mass, matches(sprintf("^(heels|flats)_(%s)$", paste(comp_vars, collapse = "|")))) |>
  pivot_longer(matches("^(heels|flats)_"), names_to = c("Condition", ".value"),
               names_pattern = "(heels|flats)_(.*)") |>
  mutate(Condition = if_else(Condition == "heels", "Heels", "Flats"),
         peak_Fr   = peak_Fr / (Mass * g)) |>
  arrange(match(paste(ParticipantID, Time), paste(new_sheet$ParticipantID, new_sheet$Time)),
          Condition) |>
  select(ParticipantID, Condition, Time, all_of(comp_vars))

if (file.exists(comp_file)) {
  comp_backup <- file.path(data_dir, sprintf("heels_flats_comp_backup_%s.csv", Sys.Date()))
  if (!file.exists(comp_backup)) invisible(file.copy(comp_file, comp_backup))
}
write_csv(heels_flats_comp, comp_file, na = "")
cat("Wrote", comp_file, "\n")
