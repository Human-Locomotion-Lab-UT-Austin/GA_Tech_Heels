# =============================================================================
# WCB2026_Stats.R
#
# Statistical analysis and figures for the high-heel intervention study
# (WCB 2026). Participants were classified as high-heel Users or Nonusers and
# tested before (Pre) and after (Post) the intervention period. At each visit
# Achilles tendon (AT) mechanics were measured while walking in flat shoes
# (Flats) and in high heels (Heels).
#
# Research questions
#   Q1. Did AT stiffness change from Pre to Post within each group, and did the
#       change differ between Users and Nonusers?
#   Q2. Did walking in heels change peak ground reaction force (GRF), effective
#       mechanical advantage (EMA), peak AT force, peak AT strain, mean AT
#       strain, and AT strain impulse relative to flats?
#   Q3. Did the change in AT stiffness scale with daily steps in heels
#       (dose-response)?
#
# AT length and cross-sectional area were also recorded but are not analyzed:
# their measurements were not reliable enough (e.g. length changes of up to
# 20% within ~100 days, recorded in 2.5-5 mm steps).
#
# Statistical approach
#   * Baseline characteristics are reported descriptively for every
#     participant and by group, without significance tests.
#   * Q1 primary inference: two-sided permutation test on the difference in
#     mean Pre to Post change between Users and Nonusers (group labels
#     reassigned without replacement).
#   * Q1 within-group changes are reported descriptively (mean, 95% t-based
#     CI, number of participants who increased). With 4 participants per group
#     a sign-flip test has only 2^4 = 16 sign combinations, so the smallest
#     attainable two-sided P is 2/16 = 0.125 and the test cannot reach
#     significance. Its P value is still computed and shown alongside min_p.
#   * Q2 (paired heels vs flats): two-sided sign-flip permutation test on the
#     mean of the per-participant differences (H0: mean difference = 0).
#   * Permutation null distributions are enumerated exactly (every sign
#     combination or group assignment) when that is feasible, which is the case
#     for every test here; otherwise 10,000 random resamples are used.
#   * Q3: simple linear regression of the % change in stiffness on daily steps
#     in heels (primary). Supplementary: the heels - flats
#     difference in peak GRF averaged across visits, reported as r with a 95% CI.
#
# Inputs
#   data/analysis/participant_data_sheet_R.csv  one row per participant x timepoint
#   data/analysis/heels_flats_comp.csv          long format footwear comparison
#   data/processed/stance_curves.csv            mean AT strain, resultant GRF and EMA moment
#                                               arms over stance for each participant,
#                                               visit and condition (run_batch_stance_curves.m)
#
# Outputs
#   output/tables/baseline_characteristics.csv
#   output/figures/ (PDF):
#   users_vs_nonusers_stiffness_fig.pdf, flats_vs_heels_kinetics_fig.pdf,
#   users_vs_nonusers_steps_fig.pdf, flats_vs_heels_grf_ema_combined_fig.pdf,
#   flats_vs_heels_strain_stance_fig.pdf, flats_vs_heels_stance_summary_fig.pdf,
#   htd_steps_reg_fig.pdf,
#   supp_peak_grf_reg_fig.pdf
# =============================================================================


# 1. Setup --------------------------------------------------------------------

library(tidyverse) # data wrangling (dplyr, purrr, readr) and ggplot2
library(mosaic)    # histogram() for permutation null distributions
library(cowplot)   # plot_grid() and ggsave2() for multi-panel figures

set.seed(2026) # makes any Monte Carlo (non-exact) permutation P value reproducible

n_perm    <- 10000 # random resamples per test when exact enumeration is too large
max_exact <- 1e5   # largest null distribution that is enumerated exactly
tol       <- sqrt(.Machine$double.eps) # tolerance when comparing permuted and observed means

# Paths are relative to the repository root. Open GA_Tech_Heels.Rproj in
# RStudio, or run `Rscript R/WCB2026_Stats.R` from the repository root.
repo_dir <- "."
if (!dir.exists(file.path(repo_dir, "data", "analysis"))) {
  stop("Run this script from the repository root (e.g. open GA_Tech_Heels.Rproj).")
}
analysis_dir  <- file.path(repo_dir, "data", "analysis")
processed_dir <- file.path(repo_dir, "data", "processed")
fig_dir       <- file.path(repo_dir, "output", "figures")
table_dir     <- file.path(repo_dir, "output", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

g <- 9.81 # gravitational acceleration (m/s^2), converts mass to body weight


# 2. Load data ----------------------------------------------------------------

# One row per participant x timepoint (Pre/Post), with flats and heels outcomes
# in separate columns
d <- read_csv(file.path(analysis_dir, "participant_data_sheet_R.csv"), show_col_types = FALSE)

# Long format footwear comparison (participant x timepoint x condition), used
# for the heels vs flats figures
heels_flats_comp <- read_csv(file.path(analysis_dir, "heels_flats_comp.csv"), show_col_types = FALSE)

# Mean stance-phase time series (0 = heel strike, 100 = toe-off, 101 points)
# for each participant x visit x condition, averaged over every stance in the
# trial: AT strain (%), resultant GRF (BW), and internal (r) and external (R)
# moment arms for EMA
stance_curves <- read_csv(file.path(processed_dir, "stance_curves.csv"), show_col_types = FALSE)


# 3. Statistical helper functions ---------------------------------------------

# Two-sided P value from a permutation null distribution. min_p is the smallest
# P value the design can produce (the share of permutations at least as extreme
# as the most extreme one), which shows when a test cannot reach significance.
perm_p <- function(perm_diff, obs_diff) {
  list(p_val = mean(abs(perm_diff) >= abs(obs_diff) - tol),
       min_p = mean(abs(perm_diff) >= max(abs(perm_diff)) - tol))
}

# Two-sided sign-flip permutation test for paired data, plus descriptive
# statistics. x holds one difference per participant (e.g. Post - Pre, or
# heels - flats). Under H0 each difference is equally likely to be positive or
# negative, so the null distribution comes from flipping the sign of each
# difference: all 2^n combinations when feasible, otherwise nperm random draws.
sign_flip_test <- function(x, mu = 0, nperm = n_perm) {
  x <- x[!is.na(x)] # drops participants missing either measurement
  n <- length(x)
  obs_diff <- mean(x) - mu
  exact <- 2^n <= max_exact
  if (exact) {
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n))) # one row per sign combination
    perm_diff <- as.vector(signs %*% abs(x - mu)) / n
  } else {
    perm_diff <- replicate(nperm, mean(sample(c(-1, 1), n, replace = TRUE) * abs(x - mu)))
  }
  print(histogram(perm_diff, v = obs_diff))
  ci <- stats::t.test(x, mu = mu)$conf.int - mu # 95% CI for obs_diff (mean - mu)
  c(list(n = n, obs_diff = obs_diff, sd = sd(x), ci_low = ci[1], ci_high = ci[2],
         n_positive = sum(x > mu)),
    perm_p(perm_diff, obs_diff),
    list(method = if (exact) "exact" else "Monte Carlo"))
}

# Two-sided permutation test for a difference in group means (User - Nonuser),
# plus descriptive statistics. The null distribution comes from reassigning the
# group labels without replacement (keeping group sizes): every possible
# assignment when feasible, otherwise nperm random shuffles.
group_perm_test <- function(x, group, nperm = n_perm) {
  keep <- !is.na(x)
  x <- x[keep]
  group <- group[keep]
  n <- length(x)
  n_user <- sum(group == "User")
  mean_diff <- function(user_idx) mean(x[user_idx]) - mean(x[-user_idx])
  obs_diff <- mean_diff(which(group == "User"))
  exact <- choose(n, n_user) <= max_exact
  if (exact) {
    perm_diff <- apply(combn(n, n_user), 2, mean_diff) # one column per possible User set
  } else {
    perm_diff <- replicate(nperm, mean_diff(sample(n, n_user)))
  }
  print(histogram(perm_diff, nint = 25, v = obs_diff))
  ci <- stats::t.test(x[group == "User"], x[group == "Nonuser"])$conf.int # Welch 95% CI
  pooled_sd <- sqrt((var(x[group == "User"]) + var(x[group == "Nonuser"])) / 2)
  c(list(n = n, obs_diff = obs_diff, ci_low = ci[1], ci_high = ci[2],
         cohens_d = obs_diff / pooled_sd),
    perm_p(perm_diff, obs_diff),
    list(method = if (exact) "exact" else "Monte Carlo"))
}

# P value label for figures (journal style, 3 decimal places)
format_p <- function(p) {
  if (p < 0.001) "P < 0.001" else paste("P =", formatC(p, digits = 3, format = "f"))
}

# Descriptive label for a within-group change: mean and 95% CI
format_change <- function(test, unit = "%") {
  f <- function(v) formatC(v, digits = 1, format = "f", flag = "+")
  paste0("Mean change ", f(test$obs_diff), unit,
         " (95% CI ", f(test$ci_low), " to ", f(test$ci_high), ")")
}


# 4. Participant baseline characteristics --------------------------------------

# One row per participant at the Pre visit. Intervention length and step counts
# describe exposure during the intervention (one value per participant).
baseline <- d |>
  filter(Time == "Pre") |>
  transmute(
    Participant              = ParticipantID,
    Group                    = Compliance,
    Sex,
    `Age (y)`                = Age,
    `Height (cm)`            = Height,
    `Mass (kg)`              = Mass,
    `BMI (kg/m^2)`           = Mass / (Height / 100)^2,
    `AT moment arm (mm)`     = ATMomentArm * 1000,
    `AT stiffness (N/mm)`    = k_lin,
    `Intervention (days)`    = InterventionDays,
    `Steps in heels per day` = HTDSteps,
    `Total steps per day`    = TotalSteps
  ) |>
  arrange(Group, Participant)

# Decimal places shown for each numeric column
baseline_digits <- c(`Age (y)` = 0, `Height (cm)` = 1, `Mass (kg)` = 1, `BMI (kg/m^2)` = 1,
                     `AT moment arm (mm)` = 0, `AT stiffness (N/mm)` = 0,
                     `Intervention (days)` = 0, `Steps in heels per day` = 0,
                     `Total steps per day` = 0)

fmt_num <- function(x, digits) {
  if_else(is.na(x), "n/a", formatC(x, digits = digits, format = "f", big.mark = ","))
}

# Mean (SD) for a group, noting n when some participants have missing data
fmt_mean_sd <- function(x, digits) {
  n_total <- length(x)
  x <- x[!is.na(x)]
  out <- paste0(fmt_num(mean(x), digits), " (", fmt_num(sd(x), digits), ")")
  if (length(x) < n_total) paste0(out, " [n = ", length(x), "]") else out
}

summarize_baseline <- function(df, label) {
  tibble(Participant = label, Group = "",
         Sex = paste0(sum(df$Sex == "F"), " F / ", sum(df$Sex == "M"), " M"),
         !!!imap(baseline_digits, ~ fmt_mean_sd(df[[.y]], .x)))
}

# Table: individual participants, then mean (SD) by group and overall.
# No significance tests: with 4 per group they have almost no power, and a
# non-significant P would wrongly suggest the groups were comparable.
baseline_table <- bind_rows(
  baseline |> mutate(across(all_of(names(baseline_digits)), ~ fmt_num(.x, baseline_digits[[cur_column()]]))) |>
    filter(Group == "Nonuser"),
  summarize_baseline(filter(baseline, Group == "Nonuser"), "Nonusers, mean (SD)"),
  baseline |> mutate(across(all_of(names(baseline_digits)), ~ fmt_num(.x, baseline_digits[[cur_column()]]))) |>
    filter(Group == "User"),
  summarize_baseline(filter(baseline, Group == "User"), "Users, mean (SD)"),
  summarize_baseline(baseline, "All, mean (SD)")
)
print(baseline_table, n = Inf, width = Inf)
write_csv(baseline_table, file.path(table_dir, "baseline_characteristics.csv"))


# 5. Q1: Pre to Post change in AT stiffness ------------------------------------

# One row per participant: group, plus % change from Pre to Post in linear AT
# stiffness (k_lin)
diff_data <- d |>
  group_by(ParticipantID) |>
  summarize(
    compliance = Compliance[Time == "Pre"],
    k_lin      = (k_lin[Time == "Post"] - k_lin[Time == "Pre"]) / k_lin[Time == "Pre"] * 100
  )

users    <- filter(diff_data, compliance == "User")
nonusers <- filter(diff_data, compliance == "Nonuser")

# Baseline imbalance: Users started with lower stiffness than Nonusers (mean
# 239 vs 286 N/mm in the current data; see baseline_table). Because the change
# is expressed as % of the Pre value, a lower starting stiffness gives a larger
# % change for the same absolute change, and participants who start low tend
# to measure higher on retest (regression to the mean). Either could
# contribute to the larger increase in Users. The overlap is wide (SE15, a
# Nonuser, had the lowest baseline at 160 N/mm), and with 4 per group the
# design cannot adjust for baseline (e.g. ANCOVA), so this is a limitation to
# report rather than something the analysis corrects.

# Primary test for Q1: did the Pre to Post change differ between Users and
# Nonusers? Nonusers act as the control group.
stiffness_group_test <- group_perm_test(diff_data$k_lin, diff_data$compliance)

# Within-group changes, reported descriptively (mean, 95% CI, n increased).
# With n = 4 per group the sign-flip P cannot fall below 0.125 (see min_p), so
# these P values are not used for inference.
stiffness_user_test    <- sign_flip_test(users$k_lin)
stiffness_nonuser_test <- sign_flip_test(nonusers$k_lin)


# 6. Q2: Heels vs flats loading and strain -------------------------------------

# Heels - flats difference for each participant at each timepoint
footwear_diff_by_time <- d |>
  transmute(
    ParticipantID, Time,
    peak_Fr_diff        = (heels_peak_Fr - flats_peak_Fr) / (Mass * g),  # body weights (BW)
    EMA_pushoff_diff    = heels_EMA_pushoff - flats_EMA_pushoff,         # EMA at push-off GRF peak
    peak_Fmtu_diff      = (heels_peak_Fmtu - flats_peak_Fmtu) / (Mass * g), # BW
    peak_strain_diff    = heels_peak_strain - flats_peak_strain,         # % strain
    mean_strain_diff    = heels_mean_strain - flats_mean_strain,         # % strain
    strain_impulse_diff = heels_strain_impulse - flats_strain_impulse    # % strain * s
  )

# Primary analysis: one heels - flats difference per participant (one value per
# person, as in Q1), averaged across the timepoints that have both conditions
footwear_diff <- footwear_diff_by_time |>
  group_by(ParticipantID) |>
  summarize(across(ends_with("_diff"), ~ mean(.x, na.rm = TRUE)))

peak_GRF_test       <- sign_flip_test(footwear_diff$peak_Fr_diff)
EMA_test            <- sign_flip_test(footwear_diff$EMA_pushoff_diff)
peak_Fmtu_test      <- sign_flip_test(footwear_diff$peak_Fmtu_diff)
strain_impulse_test <- sign_flip_test(footwear_diff$strain_impulse_diff)
peak_strain_test    <- sign_flip_test(footwear_diff$peak_strain_diff)
mean_strain_test    <- sign_flip_test(footwear_diff$mean_strain_diff)

# Robustness: repeat each test within the Pre and Post visits separately, to
# check that the footwear effect appears at both visits and is not driven by
# one of them. Fewer participants have Pre data (n = 5), so min_p is larger.
footwear_by_time <- footwear_diff_by_time |>
  pivot_longer(ends_with("_diff"), names_to = "measure", values_to = "diff") |>
  group_by(measure, Time) |>
  summarize(as_tibble(sign_flip_test(diff)), .groups = "drop") |>
  arrange(measure, desc(Time))
print(footwear_by_time, width = Inf)


# 7. Summary of permutation tests ----------------------------------------------

# obs_diff units: % change for stiffness (Pre to Post); BW (GRF, AT force),
# unitless (EMA), % strain, and % strain * s for the heels - flats comparisons. ci_low / ci_high give the
# 95% CI for obs_diff. role marks which P values are used for inference.
perm_results <- list(
  "Stiffness change, Users - Nonusers" = stiffness_group_test,
  "Stiffness change, Users"            = stiffness_user_test,
  "Stiffness change, Nonusers"         = stiffness_nonuser_test,
  "Peak GRF, heels - flats"            = peak_GRF_test,
  "EMA at push-off GRF peak, heels - flats" = EMA_test,
  "Peak AT force, heels - flats"       = peak_Fmtu_test,
  "Strain impulse, heels - flats"      = strain_impulse_test,
  "Peak strain, heels - flats"         = peak_strain_test,
  "Mean strain, heels - flats"         = mean_strain_test
) |>
  map(as_tibble) |>
  bind_rows(.id = "test") |>
  mutate(role = if_else(str_detect(test, ", (Users|Nonusers)$"), "descriptive", "inferential")) |>
  select(test, role, n, obs_diff, sd, ci_low, ci_high, n_positive, cohens_d, p_val, min_p, method)
print(perm_results, width = Inf)


# 8. Q3: Predictors of the change in AT stiffness -----------------------------

# Primary predictor: daily steps in high heels, the dose of the intervention.
# Combined with Q2 (strain was lower in heels than flats in every participant),
# this tests whether heel use was associated with stiffening even though it
# lowered AT strain per step.
# Supplementary predictor: heels - flats difference in peak GRF (BW), the
# loading variable that was higher in heels (Q2), averaged across the Pre and
# Post visits (the same per-participant values as the Q2 test, footwear_diff).
# Reported as r with a 95% CI only. See the limitations listed with
# peak_grf_reg_fig in section 10 before interpreting it.
# Not analyzed: heels - flats differences in AT strain. Nonusers did not wear
# heels consistently, so a regression of stiffening on strain across all
# participants mixes exposed and unexposed tendons and is not interpretable;
# the strain finding is the Q2 heels vs flats comparison instead.
predictor_data <- d |>
  filter(Time == "Post") |> # step counts are identical on Pre and Post rows
  select(ParticipantID, htd_steps = HTDSteps) |>
  left_join(select(footwear_diff, ParticipantID, peak_Fr_diff), by = "ParticipantID")

reg_data <- left_join(diff_data, predictor_data, by = "ParticipantID")

# One simple linear regression per predictor: % change in stiffness ~ predictor
primary_predictor        <- "htd_steps"
supplementary_predictors <- "peak_Fr_diff"

reg_models <- c(primary_predictor, supplementary_predictors) |>
  set_names() |>
  map(~ lm(reformulate(.x, response = "k_lin"), data = reg_data))

# r with a 95% CI (Fisher z) for every predictor
reg_results <- reg_models |>
  map(function(m) {
    s <- summary(m)
    ci <- stats::cor.test(m$model[[2]], m$model[[1]])$conf.int
    tibble(n = nobs(m), slope = coef(m)[[2]],
           r = sign(coef(m)[[2]]) * sqrt(s$r.squared),
           r_ci_low = ci[1], r_ci_high = ci[2],
           r_squared = s$r.squared, p_val = coef(s)[2, 4])
  }) |>
  bind_rows(.id = "predictor") |>
  mutate(role = if_else(predictor == primary_predictor, "primary", "supplementary"), .after = predictor)
print(reg_results, width = Inf)


# 9. Figure styling -----------------------------------------------------------

col_light <- "#a9a9a9" # Pre / Flats / Nonusers
col_dark  <- "#CC6600" # Post / Heels / Users

# Okabe-Ito colorblind-safe palette, one fixed color per participant so each
# participant has the same color in every figure
all_ids <- sort(unique(c(d$ParticipantID, heels_flats_comp$ParticipantID)))
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
               "#F0E442", "#D55E00", "#CC79A7", "#000000")
stopifnot(length(all_ids) <= length(okabe_ito))
participant_colors <- setNames(okabe_ito[seq_along(all_ids)], all_ids)
participant_color_scale <- scale_color_manual(values = participant_colors, name = "Participant")

# Shared base theme: white background, no grid, black axes, inward ticks
theme_wcb <- function() {
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(-4, "pt"),
    plot.title       = element_text(hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5),
    legend.position  = "none"
  )
}

# Box plot of two paired conditions with each participant's values joined by a
# line. data needs box_x (1.5 / 2.5, box centers) and x_num (1.7 / 2.3, points)
paired_box_fig <- function(data, y, line_group, x_labels, title, subtitle, y_label, y_limits) {
  # ylim() silently drops values outside the limits (and changes the box
  # statistics), so stop if any data point would be hidden
  y_values <- pull(data, {{ y }})
  if (any(y_values < y_limits[1] | y_values > y_limits[2], na.rm = TRUE)) {
    stop(sprintf("%s: data range %.2f to %.2f is outside y_limits %.2f to %.2f",
                 title, min(y_values, na.rm = TRUE), max(y_values, na.rm = TRUE),
                 y_limits[1], y_limits[2]))
  }
  ggplot(data, aes(x = box_x, y = {{ y }})) +
    geom_boxplot(aes(group = box_x, fill = factor(box_x)), width = 0.18, color = "black",
                 linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
    geom_point(aes(x = x_num, color = ParticipantID), size = 3) +
    geom_line(aes(x = x_num, group = {{ line_group }}, color = ParticipantID)) +
    scale_x_continuous(breaks = c(1.5, 2.5), labels = x_labels, limits = c(1.3, 2.7)) +
    scale_fill_manual(values = c("1.5" = col_light, "2.5" = col_dark)) +
    participant_color_scale +
    ylim(y_limits) +
    labs(title = title, subtitle = subtitle, y = y_label) +
    theme_wcb()
}

# Scatter plot of % change in stiffness against one predictor, with the least
# squares line. The subtitle gives r and P for the primary predictor, and r with
# its 95% CI for supplementary predictors (see section 8)
reg_fig <- function(predictor, x_label, x_scale, title = NULL, axis_title_size = 26) {
  res <- filter(reg_results, .data$predictor == !!predictor)
  f <- function(v) formatC(v, digits = 2, format = "f")
  subtitle <- if (res$role == "primary") {
    paste0("r = ", f(res$r), ", ", format_p(res$p_val))
  } else {
    paste0("r = ", f(res$r), " (95% CI ", f(res$r_ci_low), " to ", f(res$r_ci_high), ")")
  }
  ggplot(reg_data, aes(x = .data[[predictor]], y = k_lin, color = ParticipantID)) +
    geom_point(size = 6, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = col_dark,
                linewidth = 1, fullrange = TRUE) +
    x_scale +
    scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
    coord_cartesian(ylim = c(-25, 35)) +
    participant_color_scale +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = expression(paste(Delta, "Tendon Stiffness (", Delta, "%)"))
    ) +
    theme_wcb() +
    theme(
      plot.title    = element_text(size = 30),
      plot.subtitle = element_text(size = 20),
      axis.title.y  = element_text(color = "black", size = axis_title_size),
      axis.title.x  = element_text(color = "black", size = axis_title_size, margin = margin(t = 10)),
      axis.ticks    = element_line(color = "black", linewidth = 0.4),
      axis.text.x   = element_text(color = "black", size = 18, margin = margin(t = 10)),
      axis.text.y   = element_text(color = "black", size = 18, margin = margin(r = 10))
    )
}

# Save a figure to fig_dir at a fixed size so output does not depend on the
# size of the RStudio plot pane
save_fig <- function(plot, filename, width, height) {
  ggsave2(filename, plot = plot, path = fig_dir, width = width, height = height, units = "in")
}


# 10. Figure: Q3 regressions --------------------------------------------------
# Built first because the HTD steps panel is reused in the stiffness figure

# Primary predictor: daily steps in heels (dose-response)
htd_steps_reg_fig <- reg_fig(
  "htd_steps", "Steps in High Heels per Day",
  scale_x_continuous(breaks = seq(0, 4000, by = 1000), limits = c(0, 4000)),
  title = "Steps in Heels"
)
htd_steps_reg_fig
save_fig(htd_steps_reg_fig, "htd_steps_reg_fig.pdf", width = 11.24, height = 7.86)

# Supplementary: heels - flats difference in peak GRF, averaged across the Pre
# and Post visits.
#
# Limitations for interpretation (state these in the figure legend):
#   1. Nonusers did not wear heels consistently, so their heels - flats GRF
#      difference was measured in the lab but was not a load their tendons
#      received during the intervention. Nonusers are 4 of the 7 points, and
#      two of them (GP04, RL10) combine some of the largest GRF differences
#      with the largest stiffness losses (about -15%), which pulls the slope
#      negative. The correlation may therefore reflect the group difference
#      rather than an effect of GRF.
#   2. Only 3 Users have GRF data (DF07 has none), too few to assess the
#      relationship within the participants who actually wore heels.
#   3. Values are not equally precise: 5 participants are averaged over two
#      visits, while BO11 and GJ03 have only Post data.
#   4. GRF is the external load, not the tendon load. AT strain was lower in
#      heels in every participant (Q2) even though peak GRF was higher, so a
#      larger GRF difference does not mean a larger tendon stimulus.
#   5. Lab walking trials do not capture daily loading (number of steps in
#      heels), and with n = 7 the 95% CI for r is wide and includes zero. The
#      figure is descriptive: it is not a test, and it does not show causation.
peak_grf_reg_fig <- reg_fig(
  "peak_Fr_diff",
  expression(paste(Delta, "Peak GRF (Heels - Flats) (BW)")),
  scale_x_continuous(breaks = seq(-0.05, 0.15, by = 0.05), limits = c(-0.05, 0.15)),
  title = "Peak GRF"
)
peak_grf_reg_fig
save_fig(peak_grf_reg_fig, "supp_peak_grf_reg_fig.pdf", width = 11.24, height = 7.86)


# 11. Figure: Q1 AT stiffness, Users vs Nonusers -------------------------------

# Point (x_num) and box (box_x) x positions for the Pre / Post figures
pre_post_positions <- function(data) {
  mutate(data,
         x_num = if_else(Time == "Pre", 1.7, 2.3),
         box_x = if_else(Time == "Pre", 1.5, 2.5))
}

theme_stiffness <- theme(
  plot.title    = element_text(size = 30),
  plot.subtitle = element_text(size = 18),
  axis.title.x  = element_text(size = 20, margin = margin(t = 14)),
  axis.title.y  = element_text(size = 20, margin = margin(r = 14)),
  axis.text     = element_text(size = 16, color = "black")
)

stiffness_users_fig <- d |>
  filter(Compliance == "User") |>
  pre_post_positions() |>
  paired_box_fig(k_lin, ParticipantID, c("Pre", "Post"), "Users",
                 format_change(stiffness_user_test), "Tendon Stiffness (N/mm)", c(100, 400)) +
  xlab("Timepoint") +
  theme_stiffness

stiffness_nonusers_fig <- d |>
  filter(Compliance == "Nonuser") |>
  pre_post_positions() |>
  paired_box_fig(k_lin, ParticipantID, c("Pre", "Post"), "Nonusers",
                 format_change(stiffness_nonuser_test), "Tendon Stiffness (N/mm)", c(100, 400)) +
  xlab("Timepoint") +
  theme_stiffness

# % change in stiffness by group, dashed line = no change
stiffness_diff_fig <- diff_data |>
  mutate(x_num = if_else(compliance == "Nonuser", 1.15, 1.35),
         box_x = if_else(compliance == "Nonuser", 1, 1.5)) |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_boxplot(aes(group = box_x, fill = factor(box_x)), width = 0.12, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, color = ParticipantID), size = 3) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_fill_manual(values = c("1" = col_light, "1.5" = col_dark)) +
  participant_color_scale +
  ylim(-20, 32) +
  labs(title = "Between Groups", subtitle = format_p(stiffness_group_test$p_val), x = "Group",
       y = expression(paste(Delta, "Tendon Stiffness (", Delta, "%)"))) +
  theme_wcb() +
  theme_stiffness

combined_stiffness_fig <- plot_grid(
  plot_grid(stiffness_users_fig, stiffness_nonusers_fig,
            labels = c("A", "B"), ncol = 2, align = "hv"),
  plot_grid(stiffness_diff_fig, htd_steps_reg_fig,
            labels = c("C", "D"), ncol = 2),
  nrow = 2
)
combined_stiffness_fig
save_fig(combined_stiffness_fig, "users_vs_nonusers_stiffness_fig.pdf", width = 12.65, height = 9.88)


# 12. Figure: Q2 heels vs flats loading and strain -----------------------------

heels_flats_comp <- heels_flats_comp |>
  mutate(x_num = if_else(Condition == "Flats", 1.7, 2.3),
         box_x = if_else(Condition == "Flats", 1.5, 2.5))

theme_footwear <- theme(
  plot.title   = element_text(size = 28),
  plot.subtitle = element_text(size = 18),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 22),
  axis.text.x  = element_text(size = 20, color = "black", margin = margin(t = 10)),
  axis.text.y  = element_text(size = 20, color = "black", margin = margin(r = 10))
)

# Figure vs test: the P values come from the Q2 tests, which use one heels -
# flats difference per participant averaged across the Pre and Post visits
# (footwear_diff). The figure instead shows both visits separately: each box
# pools every participant-visit (up to 2 per participant), and lines join each
# participant's flats and heels values within a visit. State this in the
# legend, since the boxes contain more points than the n used for the tests.
footwear_fig <- function(y, title, test, y_label, y_limits) {
  paired_box_fig(heels_flats_comp, {{ y }}, interaction(ParticipantID, Time),
                 c("Flats", "Heels"), title, format_p(test$p_val), y_label, y_limits) +
    theme_footwear
}

peak_forces_fig    <- footwear_fig(peak_Fr, "Peak GRFs", peak_GRF_test,
                                   "Resultant GRFs (BW)", c(1.1, 1.45))
peak_strain_fig    <- footwear_fig(peak_strain, "Peak Strain", peak_strain_test,
                                   "Peak Tendon Strain (%)", c(2.5, 10))
mean_strain_fig    <- footwear_fig(mean_strain, "Mean Strain", mean_strain_test,
                                   "Mean Tendon Strain (%)", c(1, 4.5))
strain_impulse_fig <- footwear_fig(strain_impulse, "Strain Impulse", strain_impulse_test,
                                   "Tendon Strain Impulse (% · s)", c(0.8, 3.1))

combined_strain_fig <- plot_grid(
  peak_strain_fig, peak_forces_fig, strain_impulse_fig, mean_strain_fig,
  labels = c("A", "B", "C", "D"), ncol = 2, align = "v"
)
combined_strain_fig
save_fig(combined_strain_fig, "flats_vs_heels_kinetics_fig.pdf", width = 10.61, height = 9.88)


# 13. Figure: daily step counts, Users vs Nonusers -----------------------------

# Step counts are one average per participant over the intervention period,
# repeated on the Pre and Post rows of d. Keep one row per participant so each
# person counts once in the group mean and SE.
step_counts <- d |>
  group_by(ParticipantID, Compliance) |>
  summarize(
    n_values   = n_distinct(HTDSteps) + n_distinct(TotalSteps),
    HTDSteps   = first(HTDSteps),
    TotalSteps = first(TotalSteps),
    .groups = "drop"
  )
stopifnot(all(step_counts$n_values == 2)) # Pre and Post rows must hold the same step counts

# Bar = group mean, error bar = SE, points = individual participants.
# Participants with no step data (NA) are left out of that panel.
steps_fig <- function(y, title, y_max, y_step) {
  step_counts |>
    mutate(x_num = if_else(Compliance == "Nonuser", 1.15, 1.35),
           box_x = if_else(Compliance == "Nonuser", 1, 1.5)) |>
    ggplot(aes(x = box_x, y = {{ y }})) +
    stat_summary(aes(group = box_x, fill = factor(box_x)), fun = mean, geom = "bar",
                 width = 0.12, color = "black", linewidth = 0.5) +
    stat_summary(aes(group = box_x), fun.data = mean_se, geom = "errorbar",
                 width = 0.06, color = "black", linewidth = 0.5) +
    geom_point(aes(x = x_num, color = ParticipantID), size = 5) +
    scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.05)),
                       breaks = scales::breaks_width(y_step)) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = col_light, "1.5" = col_dark)) +
    participant_color_scale +
    labs(title = title, x = "Group", y = title) +
    theme_wcb() +
    theme(
      plot.title   = element_text(size = 30),
      axis.text    = element_text(size = 18, color = "black"),
      axis.title.y = element_text(color = "black", size = 26, margin = margin(r = 14)),
      axis.title.x = element_text(size = 26, margin = margin(t = 14))
    )
}

# Dashed line: 1000 steps per day in heels
htd_steps_fig <- steps_fig(HTDSteps, "Steps in High Heels per Day", 3800, 500) +
  geom_hline(yintercept = 1000, color = "black", linewidth = 0.5, linetype = "dashed")
total_steps_fig <- steps_fig(TotalSteps, "Total Steps per Day", 16000, 2500)

combined_steps_fig <- plot_grid(
  htd_steps_fig, total_steps_fig,
  labels = c("A", "B"), ncol = 2, align = "h", axis = "tb"
)
combined_steps_fig
save_fig(combined_steps_fig, "users_vs_nonusers_steps_fig.pdf", width = 14.28, height = 9.88)


# 14. Stance-phase time series: participant averages ------------------------

# Curves come from run_batch_stance_curves.m, which uses the same processing as
# the summary values in the data sheet and averages every stance in each trial
# after time-normalizing it to 0-100% of stance. Every curve therefore starts at
# heel strike and ends at toe-off on the same 101-point grid.
#
# Curves are kept only for participant-visits that have both heels and flats
# strain in the data sheet, so the figures use the same data as the Q2 tests.
# This drops GJ03 Pre (excluded from the data sheet as unreliable) and DF07
# (heels only). Each participant's curves are averaged across the visits in
# stance_visits, then averaged across participants (sections 15-16), so every
# participant counts once. Pre only (n = 5), Post only (n = 7) and both visits
# (n = 7) give nearly identical strain curves (peak strain within 0.25%, peak
# at 79-80% of stance).
stance_visits <- c("Pre", "Post")

valid_visits <- d |>
  filter(!is.na(heels_peak_strain), !is.na(flats_peak_strain)) |>
  select(ParticipantID, Visit = Time)

# Shared x-axis for the stance figures. Stance runs from heel strike (0%) to
# toe-off (100%); both events are marked with dashed lines and labeled below
# the axis.
stance_events <- geom_vline(xintercept = c(0, 100), linetype = "dashed",
                            color = "grey50", linewidth = 0.5)
stance_x_scale <- scale_x_continuous(
  name = "Stance (%)",
  breaks = seq(0, 100, by = 25),
  labels = c("0%\nHeel strike", "25%", "50%", "75%", "100%\nToe-off"),
  expand = expansion(mult = c(0.02, 0.02))
)

participant_curves <- stance_curves |>
  semi_join(valid_visits, by = c("ParticipantID", "Visit")) |>
  filter(Visit %in% stance_visits) |>
  group_by(ParticipantID, Condition, stance_pct) |>
  summarize(across(c(strain_mean, grf_bw, fmtu_bw, r_internal, R_external), mean), .groups = "drop")


# 15. Figure: GRF and EMA over stance -------------------------------------------

# GRF is the mean across participants of resultant GRF / body weight.
#
# EMA = r / R, where r is the AT (internal) moment arm and R the GRF (external)
# moment arm, R = ankle moment / resultant GRF. Near heel strike and toe-off the
# ankle moment is close to zero or dorsiflexor, so R approaches or crosses zero
# and EMA becomes undefined (values jump to +/- infinity). To keep the average
# stable, EMA is computed as mean r / mean R across participants, and it is only
# shown where every participant's R exceeds min_R_external (the ankle moment is
# clearly plantarflexor). In that window this matches the mean of per-stance
# EMA medians within about 0.05. Outside it EMA is left blank.
min_R_external <- 0.02 # m

grf_ema_curves <- participant_curves |>
  group_by(Condition, stance_pct) |>
  summarize(
    n          = n(),
    grf_bw     = mean(grf_bw),
    ema        = mean(r_internal) / mean(R_external),
    ema_valid  = all(R_external > min_R_external),
    .groups = "drop"
  ) |>
  group_by(stance_pct) |>
  mutate(ema_valid = all(ema_valid)) |> # same window for both conditions
  ungroup() |>
  mutate(ema = if_else(ema_valid, ema, NA_real_),
         Source     = if_else(Condition == "HEEL", "Heels", "Flats"),
         GRF_series = paste(Source, "GRF"),
         EMA_series = paste(Source, "EMA"))
stopifnot(n_distinct(grf_ema_curves$n) == 1) # same participants at every point and in both conditions

ema_window <- range(grf_ema_curves$stance_pct[!is.na(grf_ema_curves$ema)])
cat("EMA shown from", ema_window[1], "to", ema_window[2], "% of stance\n")

# GRF (BW) is on the left axis and EMA on the right. EMA is multiplied by
# ema_scale so its peak sits just below the top of the shared y range
shared_ylim   <- c(0, 1.75)
shared_breaks <- seq(0, 1.75, by = 0.25)
ema_scale     <- shared_ylim[2] / (1.05 * max(grf_ema_curves$ema, na.rm = TRUE))

combined_grf_ema_fig <- grf_ema_curves |>
  ggplot(aes(x = stance_pct)) +
  stance_events +
  geom_line(aes(y = grf_bw, color = GRF_series), linewidth = 2) +
  geom_line(aes(y = ema * ema_scale, color = EMA_series), linewidth = 2, na.rm = TRUE) +
  scale_y_continuous(
    name = "Ground Reaction Force (BW)",
    breaks = shared_breaks,
    sec.axis = sec_axis(~ . / ema_scale, name = "Effective Mechanical Advantage (r/R)")
  ) +
  stance_x_scale +
  coord_cartesian(ylim = shared_ylim) +
  scale_color_manual(name = NULL, values = c(
    "Heels GRF" = "black", "Flats GRF" = "grey60",
    "Flats EMA" = col_light, "Heels EMA" = col_dark
  )) +
  theme_wcb() +
  theme(
    text                = element_text(size = 16),
    axis.title.y.left   = element_text(color = "black", size = 26, margin = margin(r = 14)),
    axis.title.y.right  = element_text(color = "black", size = 26, margin = margin(l = 14)),
    axis.title.x        = element_text(size = 26, margin = margin(t = 14)),
    axis.text           = element_text(color = "black", size = 16),
    axis.text.y.right   = element_text(color = "black", size = 14),
    axis.line.y.right   = element_line(color = col_dark, linewidth = 0.5),
    plot.margin         = margin(r = 20, l = 20),
    legend.position     = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background   = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.margin       = margin(t = 4, r = 6, b = 4, l = 6),
    legend.key          = element_rect(fill = "white", color = NA),
    legend.text         = element_text(size = 18)
  )
combined_grf_ema_fig
save_fig(combined_grf_ema_fig, "flats_vs_heels_grf_ema_combined_fig.pdf", width = 9.21, height = 7.03)


# 16. Figure: AT strain over stance ---------------------------------------------

strain_curves <- participant_curves |>
  group_by(Condition, stance_pct) |>
  summarize(n = n(), avg_lin_strain = mean(strain_mean), sd_lin_strain = sd(strain_mean), .groups = "drop")
stopifnot(n_distinct(strain_curves$n) == 1) # same participants at every point and in both conditions

strain_stance_fig <- ggplot(strain_curves, aes(x = stance_pct, y = avg_lin_strain, color = Condition)) +
  stance_events +
  geom_line(linewidth = 2) +
  labs(y = "Strain (%)", title = "Achilles Tendon Strain Over Stance") +
  stance_x_scale +
  scale_color_manual(name = NULL, values = c("FLAT" = col_light, "HEEL" = col_dark),
                     labels = c("FLAT" = "Flats", "HEEL" = "Heels")) +
  theme_wcb() +
  theme(
    plot.title         = element_text(size = 30),
    axis.title.y       = element_text(size = 26, color = "black", margin = margin(r = 10)),
    axis.title.x       = element_text(size = 26, margin = margin(t = 10)),
    axis.text.x        = element_text(size = 18, color = "black", margin = margin(t = 10)),
    axis.text.y        = element_text(size = 18, color = "black", margin = margin(r = 10)),
    legend.position    = "inside",
    legend.position.inside = c(0.15, 0.95),
    legend.justification = c("right", "top"),
    legend.background  = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.margin      = margin(t = 4, r = 6, b = 4, l = 6),
    legend.key         = element_rect(fill = "white", color = NA),
    legend.text        = element_text(size = 20),
    plot.margin        = margin(t = 5, r = 25, b = 5, l = 5)
  )
strain_stance_fig
save_fig(strain_stance_fig, "flats_vs_heels_strain_stance_fig.pdf", width = 12.32, height = 7.03)


# 17. Figure: stance time series with heels vs flats summaries -------------------

# Four rows: GRF, EMA, AT force and AT strain. Left column: mean time series over
# stance across participants (participant_curves, section 14), heels vs flats.
# Right column: heels vs flats box plot of the matching per-trial summary value
# with its Q2 test P value. The box plots have the same figure-vs-test caveat as
# section 12: boxes pool every participant-visit, while the tests use one value
# per participant. Forces are in body weights (BW).
#
# EMA: the time series is mean r / mean R, shown only where it is defined
# (section 15). The box plot is EMA at the push-off (second) GRF peak, the
# highest GRF peak in the second half of each stance (~78% of stance). The
# largest GRF peak is not used because it is the loading peak (~25%) in most
# trials but the push-off peak in others, so it is not a consistent event.
stance_means <- participant_curves |>
  group_by(Condition, stance_pct) |>
  summarize(grf = mean(grf_bw), fmtu = mean(fmtu_bw), strain = mean(strain_mean),
            .groups = "drop") |>
  left_join(select(grf_ema_curves, Condition, stance_pct, ema),
            by = c("Condition", "stance_pct"))

# Smaller text so eight panels fit on one figure
theme_summary_grid <- theme(
  plot.title    = element_text(size = 18),
  plot.subtitle = element_text(size = 13),
  axis.title.x  = element_text(size = 14, margin = margin(t = 6)),
  axis.title.y  = element_text(size = 14, margin = margin(r = 6)),
  axis.text     = element_text(size = 12, color = "black"),
  legend.text   = element_text(size = 13)
)

# Mean time series over stance, heels vs flats. The legend is drawn in one panel
# only, and the x-axis title on the bottom row only.
stance_panel <- function(y, title, y_label, show_legend = FALSE, show_x_title = FALSE) {
  ggplot(stance_means, aes(x = stance_pct, y = {{ y }}, color = Condition)) +
    stance_events +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    stance_x_scale +
    scale_color_manual(name = NULL, values = c("FLAT" = col_light, "HEEL" = col_dark),
                       labels = c("FLAT" = "Flats", "HEEL" = "Heels")) +
    labs(title = title, y = y_label) +
    theme_wcb() +
    theme_summary_grid +
    theme(
      axis.title.x           = if (show_x_title) element_text(size = 14) else element_blank(),
      legend.position        = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.02, 0.98),
      legend.justification   = c("left", "top"),
      legend.background      = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.key             = element_rect(fill = "white", color = NA),
      plot.margin            = margin(t = 5, r = 25, b = 5, l = 5) # room for "Toe-off"
    )
}

# y limits padded 8% beyond the data range, so no point is hidden
padded_limits <- function(x) {
  r <- range(x, na.rm = TRUE)
  r + c(-1, 1) * 0.08 * diff(r)
}

summary_box_panel <- function(y, title, test, y_label) {
  paired_box_fig(heels_flats_comp, {{ y }}, interaction(ParticipantID, Time),
                 c("Flats", "Heels"), title, format_p(test$p_val), y_label,
                 padded_limits(pull(heels_flats_comp, {{ y }}))) +
    theme_summary_grid +
    theme(axis.title.x = element_blank())
}

stance_summary_fig <- plot_grid(
  stance_panel(grf, "Ground Reaction Force", "GRF (BW)", show_legend = TRUE),
  summary_box_panel(peak_Fr, "Peak GRF", peak_GRF_test, "Peak GRF (BW)"),
  stance_panel(ema, "Effective Mechanical Advantage", "EMA (r/R)"),
  summary_box_panel(EMA_pushoff, "EMA at Push-off GRF Peak", EMA_test, "EMA (r/R)"),
  stance_panel(fmtu, "Achilles Tendon Force", "AT Force (BW)"),
  summary_box_panel(peak_Fmtu, "Peak AT Force", peak_Fmtu_test, "Peak AT Force (BW)"),
  stance_panel(strain, "Achilles Tendon Strain", "AT Strain (%)", show_x_title = TRUE),
  summary_box_panel(peak_strain, "Peak Strain", peak_strain_test, "Peak AT Strain (%)"),
  ncol = 2, rel_widths = c(1.6, 1), align = "hv", axis = "tblr",
  labels = LETTERS[1:8], label_size = 16
)
stance_summary_fig
save_fig(stance_summary_fig, "flats_vs_heels_stance_summary_fig.pdf", width = 12, height = 16)
