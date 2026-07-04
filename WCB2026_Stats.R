library(tidyverse)
library(infer)
library(ggplot2)
library(ggExtra)
library(mosaic)
library(cowplot)
library(lmtest)
library(lme4)
library(jsonlite)
#help me

f <- "https://raw.githubusercontent.com/Human-Locomotion-Lab-UT-Austin/GA_Tech_Heels/refs/heads/main/participant_data_sheet_R.csv"
d <- read_csv(f, col_names = TRUE)
head(d)

diff_data <- d |>
  group_by(ParticipantID) |>
  summarize(compliance = Compliance[Time == "Pre"], k_lin = (k_lin[Time == "Post"] - k_lin[Time == "Pre"])/k_lin[Time == "Pre"] * 100, at_length = (ATLength[Time == "Post"] - ATLength[Time == "Pre"])/ATLength[Time == "Pre"] * 100)


# did stiffness change (for each group & between groups)
# Users
users <- diff_data |>
  filter(compliance == "User")
x <- users$k_lin # change in k_lin
n <- length(x)
m <- mean(x) # average change in k_lin
mu <- 0 # expected change in k_lin (null hypothesis)
(obs_diff <- m - mu) # mean difference between our sample and expected

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(stiffness_user_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# Nonusers
nonusers <- diff_data |>
  filter(compliance == "Nonuser")
x <- nonusers$k_lin # change in k_lin
n <- length(x)
m <- mean(x) # average change in k_lin
mu <- 0 # expected change in k_lin (null hypothesis)
(obs_diff <- m - mu) # mean difference between our sample and expected

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(stiffness_nonuser_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# Between Groups
obs_diff <- diff_data |>
  summarize(diff = mean(k_lin[compliance == "User"]) - mean(k_lin[compliance == "Nonuser"])) |>
  pull(diff)

n_perm <- 10000
perm_diff <- vector()
for(i in 1:n_perm) {
  diff_data |>
    mutate(compliance = sample(compliance)) |> # shuffles compliance designations without replacement (maintains proportion)
    summarize(diff = mean(k_lin[compliance == "User"]) - mean(k_lin[compliance == "Nonuser"])) |>
    pull(diff) -> perm_diff[[i]]
}
histogram(perm_diff, nint = 25, v = obs_diff)
(stiffness_user_nonuser_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# did achilles tendon length change (for each group & between groups)
# Users
users <- diff_data |>
  filter(compliance == "User")
x <- users$at_length # change in at_length
n <- length(x)
m <- mean(x) # average change in at_length
mu <- 0 # expected change in at_length (null hypothesis)
(obs_diff <- m - mu) # mean difference between our sample and expected

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(length_user_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# Nonusers
nonusers <- diff_data |>
  filter(compliance == "Nonuser")
x <- nonusers$at_length # change in at_length
n <- length(x)
m <- mean(x) # average change in at_length
mu <- 0 # expected change in at_length (null hypothesis)
(obs_diff <- m - mu) # mean difference between our sample and expected

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(length_nonuser_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# Between Groups
obs_diff <- diff_data |>
  summarize(diff = mean(at_length[compliance == "User"]) - mean(at_length[compliance == "Nonuser"])) |>
  pull(diff)

n_perm <- 10000
perm_diff <- vector()
for(i in 1:n_perm) {
  diff_data |>
    mutate(compliance = sample(compliance)) |> # shuffles compliance designations without replacement (maintains proportion)
    summarize(diff = mean(at_length[compliance == "User"]) - mean(at_length[compliance == "Nonuser"])) |>
    pull(diff) -> perm_diff[[i]]
}
histogram(perm_diff, nint = 25, v = obs_diff)
(length_user_nonuser_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# were peak forces higher (by condition)
x <- d$heels_peak_Fr # peak ground reaction forces in heels
x <- x[!is.na(x)]
n <- length(x)
m <- mean(x, na.rm = TRUE) # average peak ground reaction forces in heels
mu <- mean(d$flats_peak_Fr, na.rm = TRUE) # expected ground reaction forces if equal to flats (null hypothesis)
(obs_diff <- m - mu) # mean difference between heels and flats

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(peak_GRF_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# did the difference in peak forces change from pre to post


# was strain impulse higher (by condition)
x <- d$heels_strain_impulse # peak ground reaction forces in heels
x <- x[!is.na(x)]
n <- length(x)
m <- mean(x, na.rm = TRUE) # average peak ground reaction forces in heels
mu <- mean(d$flats_strain_impulse, na.rm = TRUE) # expected ground reaction forces if equal to flats (null hypothesis)
(obs_diff <- m - mu) # mean difference between heels and flats

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(strain_impulse_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # strain impulse was significantly lower in heel condition than flat condition


# was peak strain higher (by condition)
x <- d$heels_peak_strain # peak ground reaction forces in heels
x <- x[!is.na(x)]
n <- length(x)
m <- mean(x, na.rm = TRUE) # average peak ground reaction forces in heels
mu <- mean(d$flats_peak_strain, na.rm = TRUE) # expected ground reaction forces if equal to flats (null hypothesis)
(obs_diff <- m - mu) # mean difference between heels and flats

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(peak_strain_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


# was mean strain higher (by condition)
x <- d$heels_mean_strain # peak ground reaction forces in heels
x <- x[!is.na(x)]
n <- length(x)
m <- mean(x, na.rm = TRUE) # average peak ground reaction forces in heels
mu <- mean(d$flats_mean_strain, na.rm = TRUE) # expected ground reaction forces if equal to flats (null hypothesis)
(obs_diff <- m - mu) # mean difference between heels and flats

nperm <- 10000 # number of permutation simulations
perm_diff <- vector() # set up a dummy vector to hold results for each permutation
for (i in 1:nperm) {
  # scramble the sign of individual observed - expected weights, and then take mean
  perm_diff[[i]] <- mean(sample(c(-1, 1), length(x), replace = TRUE) * abs(x - mu))
}
histogram(perm_diff, v = obs_diff)
(mean_strain_p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # mean strain was significantly lower in heel condition compared to flat condition


# what was the factor that best predicted change in stiffness (strain impulse, peak strain, mean strain, total steps, heel steps)
reg_data <- diff_data |>
  mutate(total_steps = d$TotalSteps[d$Time == "Post"], htd_steps = d$HTDSteps[d$Time == "Post"], 
         strain_impulse_diff = d$heels_strain_impulse[d$Time == "Post"] - d$flats_strain_impulse[d$Time == "Post"], 
         peak_strain_diff = d$heels_peak_strain[d$Time == "Post"] - d$flats_peak_strain[d$Time == "Post"], 
         mean_strain_diff = d$heels_mean_strain[d$Time == "Post"] - d$flats_mean_strain[d$Time == "Post"],
         peak_Fr_diff = d$heels_peak_Fr[d$Time == "Post"] - d$flats_peak_Fr[d$Time == "Post"])

ts_m <- lm(data = reg_data, k_lin ~ total_steps)
summary(ts_m) # not significant

htd_steps_model <- lm(data = reg_data, k_lin ~ htd_steps)
summary(htd_steps_model) # not significant

strain_impulse_model <- lm(data = reg_data, k_lin ~ strain_impulse_diff)
summary(strain_impulse_model) # not significant

peak_strain_model <- lm(data = reg_data, k_lin ~ peak_strain_diff)
summary(peak_strain_model)

# Graphs
users_full <- d |>
  filter(Compliance == "User") |>
  mutate(
    x_num = ifelse(Time == "Pre", 1.7, 2.3), 
    box_x = ifelse(Time == "Pre", 1.5, 2.5)
  )
nonusers_full <- d |>
  filter(Compliance == "Nonuser") |>
  mutate(
    x_num = ifelse(Time == "Pre", 1.7, 2.3), 
    box_x = ifelse(Time == "Pre", 1.5, 2.5)
  )
diff_data <- diff_data |>
  mutate(
    x_num = ifelse(compliance == "Nonuser", 1.15, 1.35),
    box_x = ifelse(compliance == "Nonuser", 1, 1.5)
  )

# Moments (optional)
stiffness_users_fig <- users_full |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_boxplot(aes(group = Time, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = k_lin), color = "black", size = 3) +
  geom_line(aes(x = x_num, group = ParticipantID), color = "black", alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Pre", "Post"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Users", subtitle = paste("P =", stiffness_user_p_val)) +
  xlab("Timepoint") +
  ylab(expression(paste("Achilles Tendon Stiffness ", italic("(N/mm)")))) +
  ylim(100, 400)
stiffness_users_fig

stiffness_nonusers_fig <- nonusers_full |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_boxplot(aes(group = Time, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = k_lin), color = "black", size = 3) +
  geom_line(aes(x = x_num, group = ParticipantID), color = "black", alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Pre", "Post"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Nonusers", subtitle = paste("P =", stiffness_nonuser_p_val)) +
  xlab("Timepoint") +
  ylab(expression(paste("Achilles Tendon Stiffness ", italic("(N/mm)")))) +
  ylim(100, 400)
stiffness_nonusers_fig

stiffness_diff_fig <- diff_data |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_boxplot(aes(group = compliance, fill = factor(box_x)), width = 0.12, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = k_lin), color = "black", size = 3) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_fill_manual(values = c("1" = "#FFCC99", "1.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Between Groups", subtitle = paste("P =", stiffness_user_nonuser_p_val)) +
  xlab("Group") +
  ylab(expression(paste("Change in Achilles Tendon Stiffness ", italic("(%)")))) +
  ylim(-20, 32)
stiffness_diff_fig

# Top row: two plots side by side
top_row <- plot_grid(
  stiffness_users_fig, stiffness_nonusers_fig,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

# Bottom row: one plot, sized to match (using NULL as an invisible spacer)
bottom_row <- plot_grid(
  NULL, stiffness_diff_fig, NULL,
  labels = c("", "C", ""),
  ncol = 3,
  rel_widths = c(0.5, 1, 0.5)  # spacer widths control centering
)

combined_stiffness_fig <- plot_grid(
  top_row, bottom_row,
  nrow = 2
)
combined_stiffness_fig


f <- "https://raw.githubusercontent.com/Human-Locomotion-Lab-UT-Austin/GA_Tech_Heels/refs/heads/main/HeelsFlatsComp_R.csv"
heels_flats_comp <- read_csv(f, col_names = TRUE)
heels_flats_comp <- heels_flats_comp |>
  mutate(
    x_num = ifelse(Condition == "Flats", 1.7, 2.3), 
    box_x = ifelse(Condition == "Flats", 1.5, 2.5)
  )
label_data <- heels_flats_comp |>
  group_by(box_x) |>
  summarise(
    mean_peak_Fr = mean(peak_Fr, na.rm = TRUE),
    mean_peak_strain = mean(peak_strain, na.rm = TRUE),
    mean_mean_strain = mean(mean_strain, na.rm = TRUE),
    mean_strain_impulse = mean(strain_impulse, na.rm = TRUE),
    peak_Fr_upper_whisker = boxplot.stats(peak_Fr)$stats[5],  # top whisker position
    peak_strain_upper_whisker = boxplot.stats(peak_strain)$stats[5],
    mean_strain_upper_whisker = boxplot.stats(mean_strain)$stats[5],
    strain_impulse_upper_whisker = boxplot.stats(strain_impulse)$stats[5]
  )

# Peak Forces
peak_forces_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=peak_Fr)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = peak_Fr, color = ParticipantID), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = ParticipantID), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Peak Ground Reaction Forces", subtitle = paste("P =", peak_GRF_p_val)) +
  geom_text(
    data = label_data,
    aes(x = box_x, y = peak_Fr_upper_whisker, label = paste("mean =", round(mean_peak_Fr, 1))),
    vjust = -1, size = 4, inherit.aes = FALSE) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Peak Ground Reaction Forces (N)"))) +
  ylim(500, 1250)
peak_forces_fig

# Peak Strain
peak_strain_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=peak_strain)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = peak_strain, color = ParticipantID), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = ParticipantID), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  geom_text(
    data = label_data,
    aes(x = box_x, y = peak_strain_upper_whisker, label = paste("mean =", round(mean_peak_strain, 1))),
    vjust = -1, size = 4, inherit.aes = FALSE) +
  ggtitle("Peak Strain", subtitle = paste("P =", peak_strain_p_val)) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Peak Achilles Tendon Strain ", italic("(%)")))) +
  ylim(0, 10)
peak_strain_fig

mean_strain_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=mean_strain)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = mean_strain, color = ParticipantID), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = ParticipantID), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  geom_text(
    data = label_data,
    aes(x = box_x, y = mean_strain_upper_whisker, label = paste("mean =", round(mean_mean_strain, 1))),
    vjust = -1, size = 4, inherit.aes = FALSE) +
  ggtitle("Mean Strain", subtitle = paste("P =", mean_strain_p_val)) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Mean Achilles Tendon Strain ", italic("(%)")))) +
  ylim(0, 3)
mean_strain_fig

strain_impulse_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=strain_impulse)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = strain_impulse, color = ParticipantID), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = ParticipantID), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  geom_text(
    data = label_data,
    aes(x = box_x, y = strain_impulse_upper_whisker, label = paste("mean =", round(mean_strain_impulse, 1))),
    vjust = -1, size = 4, inherit.aes = FALSE) +
  ggtitle("Strain Impulse", subtitle = paste("P =", strain_impulse_p_val)) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Achilles Tendon Strain Impulse", italic("(%)")))) +
  ylim(0, 3.1)
strain_impulse_fig

# Remove x-axis text and title from the top two plots
strain_impulse_fig_top <- strain_impulse_fig +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

peak_strain_fig_top <- peak_strain_fig +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

combined_strain_fig <- plot_grid(
  strain_impulse_fig_top, peak_strain_fig_top, mean_strain_fig, peak_forces_fig,
  labels = c("A", "B", "C"),
  ncol = 2,
  align = "v"
)
combined_strain_fig


# Users vs. Nonusers HTD Steps
htd_steps <- d |>
  mutate(
    x_num = ifelse(Compliance == "Nonuser", 1.15, 1.35), 
    box_x = ifelse(Compliance == "Nonuser", 1, 1.5)
  )

htd_steps_fig <- htd_steps |>
  ggplot(aes(x = box_x, y = HTDSteps)) +
  geom_boxplot(aes(group = box_x, fill = factor(box_x)), width = 0.12, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = HTDSteps), color = "black", size = 3) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_fill_manual(values = c("1" = "#FFCC99", "1.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Steps in Heels per Day") +
  xlab("Group") +
  ylab(expression(paste("Steps in Heels per Day "))) +
  ylim(-5, 3800)
htd_steps_fig  

total_steps_fig <- htd_steps |>
  ggplot(aes(x = box_x, y = TotalSteps)) +
  geom_boxplot(aes(group = box_x, fill = factor(box_x)), width = 0.12, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = TotalSteps), color = "black", size = 3) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_fill_manual(values = c("1" = "#FFCC99", "1.5" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Total Steps per Day") +
  xlab("Group") +
  ylab(expression(paste("Total Steps per Day "))) +
  ylim(-5, 20000)
total_steps_fig  

# GRF vs. EMA Time Series
f <- "/Users/andrewthornton/Documents/WCB2026/Data/SS11_FLAT_S13_stance1_GRF_EMA.csv"
flats_grf_ema <- read_csv(f, col_names = TRUE)

scale_factor <- max(flats_grf_ema$Fr_N, na.rm = TRUE) / max(flats_grf_ema$EMA, na.rm = TRUE)

flats_grf_ema_fig <- flats_grf_ema |>
  ggplot(aes(x = Time_s)) +
  geom_line(aes(y = Fr_N, color = "GRF"), linewidth = 0.8) +
  geom_line(aes(y = EMA * scale_factor, color = "EMA"), linewidth = 0.8) +
  scale_y_continuous(
    name = "Ground Reaction Force (N)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Effective Mechanical Advantage")
  ) +
  scale_color_manual(name = NULL, values = c("GRF" = "black", "EMA" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black"),
    legend.position = "top"
  ) +
  xlab(expression(paste("Time (s)"))) +
  ggtitle("Flats")
flats_grf_ema_fig

f <- "/Users/andrewthornton/Documents/WCB2026/Data/SS11_HEEL_S13_stance1_GRF_EMA.csv"
heels_grf_ema <- read_csv(f, col_names = TRUE)

scale_factor <- max(heels_grf_ema$Fr_N, na.rm = TRUE) / max(heels_grf_ema$EMA, na.rm = TRUE)

heels_grf_ema_fig <- heels_grf_ema |>
  ggplot(aes(x = Time_s)) +
  geom_line(aes(y = Fr_N, color = "GRF"), linewidth = 0.8) +
  geom_line(aes(y = EMA * scale_factor, color = "EMA"), linewidth = 0.8) +
  scale_y_continuous(
    name = "Ground Reaction Force (N)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Effective Mechanical Advantage")
  ) +
  scale_color_manual(name = NULL, values = c("GRF" = "black", "EMA" = "#CC6600")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black"),
    legend.position = "top"
  ) +
  xlab(expression(paste("Time (s)"))) +
  ggtitle("Heels")
heels_grf_ema_fig

combined_grf_ema_fig <- plot_grid(
  flats_grf_ema_fig, heels_grf_ema_fig,
  labels = c("A", "B"),
  ncol = 1,
  align = "h"
)
combined_grf_ema_fig


# Strain over a stride
strain_data <- fromJSON("/Users/andrewthornton/Documents/WCB2026/Data/strain_stride_data.json")

heel_data <- strain_data$HEEL

# common 0-100% stride grid to interpolate everyone onto
pct_grid <- seq(0, 100, length.out = 101)

interp_list <- map(heel_data, function(df) {
  approx(df$stride_percent, df$lin_strain_percent, xout = pct_grid)$y
})
interp_matrix <- do.call(cbind, interp_list)

avg_strain <- rowMeans(interp_matrix, na.rm = TRUE)

# average stride duration (s) across participants, used to rescale the x-axis
avg_stride_duration <- mean(map_dbl(heel_data, ~ max(.x$time_s)))
time_axis_s <- pct_grid / 100 * avg_stride_duration

plot_df <- data.frame(
  time_s = time_axis_s,
  avg_lin_strain = avg_strain
)

heels_strain_stride_fig <- plot_df |>
  ggplot(aes(x = time_s, y = avg_lin_strain)) +
  geom_line(color = "#CC6600", linewidth = 0.8) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  xlab("Time (s)") +
  ylab("Achilles Tendon Strain (%)") +
  ggtitle("Heels")
heels_strain_stride_fig



flat_data <- strain_data$FLAT

# common 0-100% stride grid to interpolate everyone onto
pct_grid <- seq(0, 100, length.out = 101)

interp_list <- map(heel_data, function(df) {
  approx(df$stride_percent, df$lin_strain_percent, xout = pct_grid)$y
})
interp_matrix <- do.call(cbind, interp_list)

avg_strain <- rowMeans(interp_matrix, na.rm = TRUE)

# average stride duration (s) across participants, used to rescale the x-axis
avg_stride_duration <- mean(map_dbl(flat_data, ~ max(.x$time_s)))
time_axis_s <- pct_grid / 100 * avg_stride_duration

plot_df <- data.frame(
  time_s = time_axis_s,
  avg_lin_strain = avg_strain
)

flats_strain_stride_fig <- plot_df |>
  ggplot(aes(x = time_s, y = avg_lin_strain)) +
  geom_line(color = "black", linewidth = 0.8) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  xlab("Time (s)") +
  ylab("Achilles Tendon Strain (%)") +
  ggtitle("Flats")
flats_strain_stride_fig


combined_strain_stride_fig <- plot_grid(
  flats_strain_stride_fig, heels_strain_stride_fig,
  labels = c("A", "B"),
  ncol = 1,
  align = "h"
)
combined_strain_stride_fig
