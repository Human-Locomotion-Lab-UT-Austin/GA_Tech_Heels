library(tidyverse)
library(infer)
library(ggplot2)
library(ggExtra)
library(mosaic)
library(cowplot)
library(lmtest)
library(lme4)
library(jsonlite)
install.packages("viridis")
library(viridis)
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
post_data <- d |>
  filter(Time == "Post") |>
  mutate(
    strain_impulse_diff = (heels_strain_impulse - flats_strain_impulse) / flats_strain_impulse * 100,
    peak_strain_diff = (heels_peak_strain - flats_peak_strain) / flats_peak_strain * 100,
    mean_strain_diff = (heels_mean_strain - flats_mean_strain) / flats_mean_strain * 100,
    peak_Fr_diff = (heels_peak_Fr / (Mass * 9.81)) - (flats_peak_Fr / (Mass * 9.81)),
    peak_Fmtu_diff = (heels_peak_Fmtu / (Mass * 9.81)) - (flats_peak_Fmtu / (Mass * 9.81)),
    mean_Fmtu_diff = (heels_mean_Fmtu / (Mass * 9.81))  - (flats_mean_Fmtu / (Mass * 9.81)),
    peak_ank_mom_diff = heels_peak_ank_mom - flats_peak_ank_mom,
    mean_ank_mom_diff = heels_mean_ank_mom - flats_mean_ank_mom
  ) |>
  select(ParticipantID, total_steps = TotalSteps, htd_steps = HTDSteps,
         strain_impulse_diff, peak_strain_diff, mean_strain_diff, peak_Fr_diff, peak_Fmtu_diff, mean_Fmtu_diff, peak_ank_mom_diff, mean_ank_mom_diff)

reg_data <- diff_data |>
  left_join(post_data, by = "ParticipantID")

ts_m <- lm(data = reg_data, k_lin ~ total_steps)
summary(ts_m) # not significant

htd_steps_model <- lm(data = reg_data, k_lin ~ htd_steps)
summary(htd_steps_model) # not significant

strain_impulse_model <- lm(data = reg_data, k_lin ~ strain_impulse_diff)
summary(strain_impulse_model) # not significant

peak_strain_model <- lm(data = reg_data, k_lin ~ peak_strain_diff)
summary(peak_strain_model) # not significant

peak_Fr_model <- lm(data = reg_data, k_lin ~ peak_Fr_diff)
summary(peak_Fr_model) # not significant

mean_strain_model <- lm(data = reg_data, k_lin ~ mean_strain_diff)
summary(mean_strain_model) # not significant

peak_Fmtu_model <- lm(data = reg_data, k_lin ~ peak_Fmtu_diff)
summary(peak_Fmtu_model) # not significant

mean_Fmtu_model <- lm(data = reg_data, k_lin ~ mean_Fmtu_diff)
summary(mean_Fmtu_model) # not significant

peak_ank_mom_model <- lm(data = reg_data, k_lin ~ peak_ank_mom_diff)
summary(peak_ank_mom_model) # not significant

mean_ank_mom_model <- lm(data = reg_data, k_lin ~ mean_ank_mom_diff)
summary(mean_ank_mom_model) # not significant

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

user_whisker_data <- users_full |>
  group_by(box_x) |>
  summarise(
    mean_k_lin = mean(k_lin, na.rm = TRUE),
    k_lin_whisker = boxplot.stats(k_lin)$stats[5],  # top whisker position
  )

nonuser_whisker_data <- nonusers_full |>
  group_by(box_x) |>
  summarise(
    mean_k_lin = mean(k_lin, na.rm = TRUE),
    k_lin_whisker = boxplot.stats(k_lin)$stats[5],  # top whisker position
  )
diff_whisker_data <- diff_data |>
  group_by(box_x) |>
  summarise(
    mean_k_lin = mean(k_lin, na.rm = TRUE),
    k_lin_whisker = boxplot.stats(k_lin)$stats[5],  # top whisker position
  )

f <- "https://raw.githubusercontent.com/Human-Locomotion-Lab-UT-Austin/GA_Tech_Heels/refs/heads/main/heels_flats_comp.csv"
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
htd_steps <- d |>
  mutate(
    x_num = ifelse(Compliance == "Nonuser", 1.15, 1.35), 
    box_x = ifelse(Compliance == "Nonuser", 1, 1.5)
  )
steps_label_data <- htd_steps |>
  group_by(box_x) |>
  summarise(
    mean_htd_steps = mean(HTDSteps, na.rm = TRUE),
    mean_total_steps = mean(TotalSteps, na.rm = TRUE),
    htd_steps_upper_whisker = boxplot.stats(HTDSteps)$stats[5],  # top whisker position
    total_steps_upper_whisker = boxplot.stats(TotalSteps)$stats[5]
  )
f <- "/Users/andrewthornton/Documents/WCB2026/Data/SS11_HEEL_S13_stance1_GRF_EMA.csv"
heels_grf_ema <- read_csv(f, col_names = TRUE)

f <- "/Users/andrewthornton/Documents/WCB2026/Data/SS11_FLAT_S13_stance1_GRF_EMA.csv"
flats_grf_ema <- read_csv(f, col_names = TRUE)

strain_data <- fromJSON("/Users/andrewthornton/Documents/WCB2026/Data/strain_stride_data.json")

heel_data <- strain_data$HEEL

# --- Master participant color palette (shared across ALL figures) ---
all_ids <- sort(unique(c(
  users_full$ParticipantID,
  nonusers_full$ParticipantID,
  diff_data$ParticipantID,
  heels_flats_comp$ParticipantID,
  htd_steps$ParticipantID
)))

okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#0072B2", # blue
  "#F0E442", # yellow
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#000000"  # black
)

participant_colors <- setNames(okabe_ito, all_ids)
participant_color_scale <- scale_color_manual(values = participant_colors, name = "Participant")

stiffness_users_fig <- users_full |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_boxplot(aes(group = Time, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = k_lin, color = factor(ParticipantID)), size = 3) +
  geom_line(aes(x = x_num, group = ParticipantID, color = factor(ParticipantID)), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Pre", "Post"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_text(size = 20, margin = margin(t = 14)),
    axis.title.y = element_text(size = 20, margin = margin(r = 14)),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
 # geom_text(
   # data = user_whisker_data,
   # aes(x = box_x, y = k_lin_whisker, label = paste("mean =", round(mean_k_lin, 0))),
   # vjust = -1, size = 3, inherit.aes = FALSE) +
  ggtitle("Users", subtitle = paste("P =", round(stiffness_user_p_val, 3))) +
  xlab("Timepoint") +
  ylab(expression(paste("Tendon Stiffness (N/mm)"))) +
  ylim(100, 400)
stiffness_users_fig

stiffness_nonusers_fig <- nonusers_full |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_boxplot(aes(group = Time, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = k_lin, color = factor(ParticipantID)), size = 3) +
  geom_line(aes(x = x_num, group = ParticipantID, color = factor(ParticipantID)), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Pre", "Post"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_text(size = 20, margin = margin(t = 14)),
    axis.title.y = element_text(size = 20, margin = margin(r = 14)),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
 # geom_text(
    #data = nonuser_whisker_data,
    #aes(x = box_x, y = k_lin_whisker, label = paste("mean =", round(mean_k_lin, 0))),
    #vjust = -1, size = 3, inherit.aes = FALSE) +
  ggtitle("Nonusers", subtitle = paste("P =", round(stiffness_nonuser_p_val, 3))) +
  xlab("Timepoint") +
  ylab(expression(paste("Tendon Stiffness (N/mm)"))) +
  ylim(100, 400)
stiffness_nonusers_fig

stiffness_diff_fig <- diff_data |>
  ggplot(aes(x = box_x, y = k_lin)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_boxplot(aes(group = compliance, fill = factor(box_x)), width = 0.12, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = k_lin, color = factor(ParticipantID)), size = 3) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_fill_manual(values = c("1" = "#FFCC99", "1.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_text(size = 20, margin = margin(t = 14)),
    axis.title.y = element_text(size = 20, margin = margin(r = 14)),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  #geom_text(
    #data = diff_whisker_data,
    #aes(x = box_x, y = k_lin_whisker, label = paste("mean =", round(mean_k_lin, 1))),
    #vjust = -1, size = 3, inherit.aes = FALSE) +
  ggtitle("Between Groups", subtitle = paste("P =", round(stiffness_user_nonuser_p_val, 3))) +
  xlab("Group") +
  ylab(expression(paste(Delta, "Tendon Stiffness (", Delta, "%)"))) +
  ylim(-20, 32)
stiffness_diff_fig

top_row <- plot_grid(
  stiffness_users_fig, stiffness_nonusers_fig,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

bottom_row <- plot_grid(
  NULL, stiffness_diff_fig, NULL,
  labels = c("", "C", ""),
  ncol = 3,
  rel_widths = c(0.5, 1, 0.5)
)

combined_stiffness_fig <- plot_grid(
  top_row, bottom_row,
  nrow = 2
)
combined_stiffness_fig

ggsave2("users_vs_nonusers_stiffness_fig.pdf", plot = combined_stiffness_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")




# Peak Forces
peak_forces_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=peak_Fr)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = peak_Fr, color = factor(ParticipantID)), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = factor(ParticipantID)), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 28, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, color = "black", margin = margin(t=10)),
    axis.text.y = element_text(size = 20, color = "black", margin = margin(r=10)),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Peak GRFs", subtitle = paste("P =", round(peak_GRF_p_val, 3))) +
 #geom_text(
    #data = label_data,
    #aes(x = box_x, y = peak_Fr_upper_whisker, label = paste("mean =", round(mean_peak_Fr, 1))),
    #vjust = -1, size = 4, inherit.aes = FALSE) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Resultant GRFs (BW)"))) +
  ylim(1.1, 1.45)

peak_forces_fig

# Peak Strain
peak_strain_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=peak_strain)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = peak_strain, color = factor(ParticipantID)), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = factor(ParticipantID)), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 28, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, color = "black", margin = margin(t=10)),
    axis.text.y = element_text(size = 20, color = "black", margin = margin(r=10)),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
    #data = label_data,
    #aes(x = box_x, y = peak_strain_upper_whisker, label = paste("mean =", round(mean_peak_strain, 1))),
    #vjust = -1, size = 4, inherit.aes = FALSE) +
  ggtitle("Peak Strain", subtitle = paste("P =", round(peak_strain_p_val, 3))) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Peak Tendon Strain (%)"))) +
  ylim(2.5, 10)
peak_strain_fig

mean_strain_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=mean_strain)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = mean_strain, color = factor(ParticipantID)), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = factor(ParticipantID)), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 28, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, color = "black", margin = margin(t=10)),
    axis.text.y = element_text(size = 20, color = "black", margin = margin(r=10)),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  #geom_text(
    #data = label_data,
    #aes(x = box_x, y = mean_strain_upper_whisker, label = paste("mean =", round(mean_mean_strain, 1))),
    #vjust = -1, size = 4, inherit.aes = FALSE) +
  ggtitle("Mean Strain", subtitle = paste("P =", round(mean_strain_p_val, 3))) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Mean Tendon Strain (%)"))) +
  ylim(0.8, 3)
mean_strain_fig

strain_impulse_fig <- heels_flats_comp |>
  ggplot(aes(x=box_x, y=strain_impulse)) +
  geom_boxplot(aes(group = Condition, fill = factor(box_x)), width = 0.18, alpha = 1, color = "black",
               linewidth = 0.5, whisker.linewidth = 0.5, staplewidth = 0.5) +
  geom_point(aes(x = x_num, y = strain_impulse, color = factor(ParticipantID)), size = 3) +
  geom_line(aes(x = x_num, group = interaction(ParticipantID, Time), color = factor(ParticipantID)), alpha = 1) +
  scale_x_continuous(breaks = c(1.5, 2.5), labels = c("Flats", "Heels"), limits = c(1.3, 2.7)) +
  scale_fill_manual(values = c("1.5" = "#FFCC99", "2.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 28, hjust = 0.5),
    plot.subtitle = element_text(hjust=0.5, size = 18),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, color = "black", margin = margin(t=10)),
    axis.text.y = element_text(size = 20, color = "black", margin = margin(r=10)),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  #geom_text(
    #data = label_data,
    #aes(x = box_x, y = strain_impulse_upper_whisker, label = paste("mean =", round(mean_strain_impulse, 1))),
    #vjust = -1, size = 4, inherit.aes = FALSE) +
  ggtitle("Strain Impulse", subtitle = paste("P =", round(strain_impulse_p_val, 3))) +
  xlab("Footwear Condition") +
  ylab(expression(paste("Tendon Strain Impulse (% ", "\u00B7", " s)"))) +
  ylim(0.8, 3.1)
strain_impulse_fig


combined_strain_fig <- plot_grid(
  peak_strain_fig, peak_forces_fig, strain_impulse_fig, mean_strain_fig,
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  align = "v"
)

combined_strain_fig

ggsave2("flats_vs_heels_kinetics_fig.pdf", combined_strain_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")


# Users vs. Nonusers HTD Steps
htd_steps_fig <- htd_steps |>
  ggplot(aes(x = box_x, y = HTDSteps)) +
  geom_hline(yintercept = 1000, color = "black", linewidth = 0.5, linetype = "dashed") +
  stat_summary(aes(group = box_x, fill = factor(box_x)), fun = mean, geom = "bar",
               width = 0.12, color = "black", linewidth = 0.5, alpha = 1) +
  stat_summary(aes(group = box_x), fun.data = mean_se, geom = "errorbar",
               width = 0.06, color = "black", linewidth = 0.5) +
  geom_point(aes(x = x_num, y = HTDSteps, color = factor(ParticipantID)), size = 5) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_y_continuous(limits = c(0, 3800), expand = expansion(mult = c(0, 0.05)), breaks = scales::breaks_width(500)) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("1" = "#FFCC99", "1.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 18, color = "black"),
    axis.title.y.left = element_text(color = "black", size = 26, margin = margin(r = 14)),
    axis.title.x = element_text(size = 26, margin = margin(t = 14)),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Steps in High Heels per Day") +
  xlab("Group") +
  ylab(expression(paste("Steps in High Heels per Day ")))
htd_steps_fig

total_steps_fig <- htd_steps |>
  ggplot(aes(x = box_x, y = TotalSteps)) +
  stat_summary(aes(group = box_x, fill = factor(box_x)), fun = mean, geom = "bar",
               width = 0.12, color = "black", linewidth = 0.5, alpha = 1) +
  stat_summary(aes(group = box_x), fun.data = mean_se, geom = "errorbar",
               width = 0.06, color = "black", linewidth = 0.5) +
  geom_point(aes(x = x_num, y = TotalSteps, color = factor(ParticipantID)), size = 5) +
  scale_x_continuous(breaks = c(1, 1.5), labels = c("Nonusers", "Users"), limits = c(0.9, 1.6)) +
  scale_y_continuous(limits = c(0, 16000), expand = expansion(mult = c(0, 0.05)), breaks = scales::breaks_width(2500)) +
  scale_fill_manual(values = c("1" = "#FFCC99", "1.5" = "#CC6600")) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 18, color = "black"),
    axis.title.y.left = element_text(color = "black", size = 26, margin = margin(r = 14)),
    axis.title.x = element_text(size = 26, margin = margin(t = 14)),
    axis.ticks.length = unit(-4, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  ggtitle("Total Steps per Day") +
  xlab("Group") +
  ylab(expression(paste("Total Steps per Day ")))
total_steps_fig

combined_steps_fig <- plot_grid(
  htd_steps_fig, total_steps_fig,
  labels = c("A", "B"),
  ncol = 2,
  align = "h",
  axis = "tb"
)
combined_steps_fig

ggsave2("users_vs_nonusers_steps_fig.pdf", combined_steps_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")


# --- ONE shared scale factor across both datasets ---
# --- Combine the two datasets into one, tagging each with its source ---
combined_data <- bind_rows(
  flats_grf_ema |> mutate(Source = "Flats"),
  heels_grf_ema |> mutate(Source = "Heels")
)

# --- Convert Time_s to % of stance, computed separately for each Source ---
combined_data <- combined_data |>
  group_by(Source) |>
  mutate(
    Pct_Stance = (Time_s - min(Time_s, na.rm = TRUE)) /
      (max(Time_s, na.rm = TRUE) - min(Time_s, na.rm = TRUE)) * 100
  ) |>
  ungroup()

# Build a single grouping variable so each of the 4 lines gets its own color
combined_data <- combined_data |>
  mutate(
    GRF_series = paste(Source, "GRF"),
    EMA_series = paste(Source, "EMA")
  )

combined_grf_ema_fig <- combined_data |>
  ggplot(aes(x = Pct_Stance)) +
  geom_line(aes(y = Fr_N / 510.12, color = GRF_series), linewidth = 2) +
  geom_line(aes(y = EMA * combined_scale_factor, color = EMA_series), linewidth = 2) +
  scale_y_continuous(
    name = "Ground Reaction Force (BW)",
    breaks = shared_breaks,
    sec.axis = sec_axis(~ . / combined_scale_factor, name = "Effective Mechanical Advantage (r/R)")
  ) +
  scale_x_continuous(
    name = "Stance (%)",
    breaks = seq(0, 100, by = 25),
    labels = scales::label_number(suffix = "%")
  ) +
  coord_cartesian(ylim = shared_ylim) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Heels GRF" = "black",
      "Flats GRF" = "grey60",
      "Flats EMA" = "#FFCC99",
      "Heels EMA" = "#CC6600"
    )
  ) +
  theme(
    text = element_text(size = 16),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 24, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(color = "black", linewidth = 0.5),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26, margin = margin(r = 14)),
    axis.title.y.right = element_text(color = "black", size = 26, margin = margin(l = 14)),
    axis.title.x = element_text(size = 26, margin = margin(t = 14)),
    axis.text.y.right = element_text(color = "black", size = 14),
    axis.line.y.right = element_line(color = "#CC6600", linewidth = 0.5),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 16),
    plot.margin = margin(r = 20, l = 20),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.margin = margin(t = 4, r = 6, b = 4, l = 6),
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 18)
  )

combined_grf_ema_fig

ggsave2("flats_vs_heels_grf_ema_combined_fig.pdf", combined_grf_ema_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")



# Strain over a stride
# common 0-100% stride grid to interpolate everyone onto
pct_grid <- seq(0, 100, length.out = 101)
interp_list <- map(heel_data, function(df) {
  approx(df$stride_percent, df$lin_strain_percent, xout = pct_grid)$y
})
interp_matrix <- do.call(cbind, interp_list)
heels_avg_strain <- rowMeans(interp_matrix, na.rm = TRUE)
heels_plot_df <- data.frame(
  stride_pct = pct_grid,
  avg_lin_strain = heels_avg_strain,
  condition = "HEEL"
)

flat_data <- strain_data$FLAT
interp_list <- map(flat_data, function(df) {
  approx(df$stride_percent, df$lin_strain_percent, xout = pct_grid)$y
})
interp_matrix <- do.call(cbind, interp_list)
flats_avg_strain <- rowMeans(interp_matrix, na.rm = TRUE)
flats_plot_df <- data.frame(
  stride_pct = pct_grid,
  avg_lin_strain = flats_avg_strain,
  condition = "FLAT"
)

combined_df <- rbind(heels_plot_df, flats_plot_df)

strain_stride_fig <- ggplot(combined_df, aes(x = stride_pct, y = avg_lin_strain, color = condition)) +
  geom_line(linewidth = 2) +
  labs(
    x = "Stride (%)",
    y = "Strain (%)",
    color = "Location",
    title = "Achilles Tendon Strain Over a Stride"
  ) +
  scale_x_continuous(breaks = seq(0, 100, by = 25), labels = scales::label_number(suffix = "%")) +
  scale_color_manual(name = NULL, values = c("FLAT" = "black", "HEEL" = "#CC6600"), labels = c("FLAT" = "Flats", "HEEL" = "Heels")) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(size= 26, color = "black", margin = margin(r = 10)),
    axis.title.x = element_text(size = 26, margin = margin(t = 10)),
    axis.ticks.length = unit(-4, "pt"),
    axis.text.x = element_text(size = 18, color = "black", margin = margin(t=10)),
    axis.text.y = element_text(size = 18, color = "black", margin = margin(r=10)),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.margin = margin(t = 4, r = 6, b = 4, l = 6),
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 20)
  )
strain_stride_fig

ggsave2("flats_vs_heels_strain_stride_fi g.pdf", strain_stride_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")



# LINEAR REGRESSIONS
# --- Plot 1: Total Steps ---
s1 <- summary(ts_m)
r2_1 <- formatC(s1$r.squared, digits = 3, format = "f")
p_1 <- s1$coefficients[2, 4]
p_label_1 <- if (p_1 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_1, digits = 3, format = "f"))

total_steps_reg_fig <- ggplot(reg_data, aes(x = total_steps, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(0, 20000, by = 1000),
    labels = function(x) ifelse(x %% 2000 == 0, x, ""),
    limits = c(0, 16000)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = "Total Steps per Day",
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(R^2 == .(r2_1) * ", " * .(p_label_1))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none"
  ) +
total_steps_reg_fig

# --- Plot 2: HTD Steps ---
s2 <- summary(htd_steps_model)
r_2 <- formatC(sign(coef(htd_steps_model)[2]) * sqrt(s2$r.squared), digits = 3, format = "f")
p_2 <- s2$coefficients[2, 4]
p_label_2 <- if (p_2 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_2, digits = 3, format = "f"))
htd_steps_reg_fig <- ggplot(reg_data, aes(x = htd_steps, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 6, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(0, 20000, by = 500),
    labels = function(x) ifelse(x %% 1000 == 0, x, ""),
    limits = c(0, 4000)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = "Steps in High Heels per Day",
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(r == .(r_2) * ", " * .(p_label_2))
  ) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text.x = element_text(color = "black", size = 18, margin = margin(t=10)),
    axis.text.y = element_text(color = "black", size = 18, margin = margin(r=10)),
    legend.position = "none"
  )
htd_steps_reg_fig

ggsave2("htd_steps_reg_fig.pdf", htd_steps_reg_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")


# --- Plot 3: Strain Impulse Difference ---
s3 <- summary(strain_impulse_model)
r_3 <- formatC(sign(coef(strain_impulse_model)[2]) * sqrt(s3$r.squared), digits = 3, format = "f")
p_3 <- s3$coefficients[2, 4]
p_label_3 <- if (p_3 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_3, digits = 3, format = "f"))

strain_impulse_reg_fig <- ggplot(reg_data, aes(x = strain_impulse_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 6, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-50, 50, by = 10),
    limits = c(-35, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = expression(paste(Delta, "Strain Impulse (Heels - Flats) (", Delta, "%)")),
    y = (expression(paste(Delta, "Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(r == .(r_3) * ", " * .(p_label_3))
  ) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 24),
    axis.title.x = element_text(color = "black", size = 24, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text.x = element_text(color = "black", size = 18, margin = margin(t=10)),
    axis.text.y = element_text(color = "black", size = 18, margin = margin(r=10)),
    legend.position = "none"
  ) +
  ggtitle("Strain Impulse")
strain_impulse_reg_fig

# --- Plot 4: Peak Strain Difference ---
s4 <- summary(peak_strain_model)
r_4 <- formatC(sign(coef(peak_strain_model)[2]) * sqrt(s4$r.squared), digits = 3, format = "f")
p_4 <- s4$coefficients[2, 4]
p_label_4 <- if (p_4 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_4, digits = 3, format = "f"))

peak_strain_reg_fig <- ggplot(reg_data, aes(x = peak_strain_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 6, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-30, 10, by = 5),
    limits = c(-15, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = expression(paste(Delta, "Peak Strain (Heels - Flats) (", Delta, "%)")),
    y = expression(paste(Delta, "Tendon Stiffness (", Delta, "%)")),
    subtitle = bquote(r == .(r_4) * ", " * .(p_label_4))
  ) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 24),
    axis.title.x = element_text(color = "black", size = 24, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text.x = element_text(color = "black", size = 18, margin = margin(t=10)),
    axis.text.y = element_text(color = "black", size = 18, margin = margin(r=10)),
    legend.position = "none"
  ) +
  ggtitle("Peak Strain")

peak_strain_reg_fig

# --- Plot 5: Peak Fr Difference ---
s5 <- summary(peak_Fr_model)
r2_5 <- formatC(s5$r.squared, digits = 3, format = "f")
p_5 <- s5$coefficients[2, 4]
p_label_5 <- if (p_5 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_5, digits = 3, format = "f"))

peak_Fr_reg_fig <- ggplot(reg_data, aes(x = peak_Fr_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 6, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    #breaks = seq(-1000, 1000, by = 0.01),
    #labels = function(x) ifelse(x %% 25 == 0, x, ""),
    limits = c(-.05, 0.13)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = expression(paste(Delta, "Peak GRFs (Heels - Flats) (%BW)")),
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(R^2 == .(r2_5) * ", " * .(p_label_5))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none"
  ) +
  ggtitle("Peak Ground Reaction Forces")
peak_Fr_reg_fig

# --- Plot 6: Mean Strain Difference ---
s6 <- summary(mean_strain_model)
r_6 <- formatC(sign(coef(mean_strain_model)[2]) * sqrt(s6$r.squared), digits = 3, format = "f")
p_6 <- s6$coefficients[2, 4]
p_label_6 <- if (p_6 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_6, digits = 3, format = "f"))

mean_strain_reg_fig <- ggplot(reg_data, aes(x = mean_strain_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 6, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-50, 10, by = 10),
    limits = c(-35, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = expression(paste(Delta, "Mean Strain (Heels - Flats) (", Delta, "%)")),
    y = (expression(paste(Delta, "Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(r == .(r_6) * ", " * .(p_label_6))
  ) +
  participant_color_scale +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 23),
    axis.title.x = element_text(color = "black", size = 23, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text.x = element_text(color = "black", size = 18, margin = margin(t=10)),
    axis.text.y = element_text(color = "black", size = 18, margin = margin(r=10)),
    legend.position = "none"
  ) +
  ggtitle("Mean Strain")
mean_strain_reg_fig

# --- Plot 7: Peak Fmtu Difference ---
s7 <- summary(peak_Fmtu_model)
r2_7 <- formatC(s7$r.squared, digits = 3, format = "f")
p_7 <- s7$coefficients[2, 4]
p_label_7 <- if (p_7 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_7, digits = 3, format = "f"))

peak_Fmtu_reg_fig <- ggplot(reg_data, aes(x = peak_Fmtu_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-1000, 1000, by = 50),
    labels = function(x) ifelse(x %% 100 == 0, x, ""),
    limits = c(-500, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = "Difference in Peak Achilles Tendon Forces (Heels - Flats)",
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(R^2 == .(r2_7) * ", " * .(p_label_7))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none"
  ) 
peak_Fmtu_reg_fig


# --- Plot 8: Mean Fmtu Difference ---
s8 <- summary(mean_Fmtu_model)
r2_8 <- formatC(s8$r.squared, digits = 3, format = "f")
p_8 <- s8$coefficients[2, 4]
p_label_8 <- if (p_8 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_8, digits = 3, format = "f"))

mean_Fmtu_reg_fig <- ggplot(reg_data, aes(x = mean_Fmtu_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-1000, 1000, by = 50),
    labels = function(x) ifelse(x %% 50 == 0, x, ""),
    limits = c(-250, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = "Difference in Mean Achilles Tendon Forces (Heels - Flats)",
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(R^2 == .(r2_8) * ", " * .(p_label_8))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none"
  )
mean_Fmtu_reg_fig

# --- Plot 9: Peak Ankle Moment Difference ---
s9 <- summary(peak_ank_mom_model)
r2_9 <- formatC(s9$r.squared, digits = 3, format = "f")
p_9 <- s9$coefficients[2, 4]
p_label_9 <- if (p_9 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_9, digits = 3, format = "f"))

peak_ank_mom_reg_fig <- ggplot(reg_data, aes(x = peak_ank_mom_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-1000, 0, by = .1),
    labels = function(x) ifelse(x %% .1 == 0, x, ""),
    limits = c(-0.5, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = "Difference in Peak Ankle Moments (Heels - Flats)",
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(R^2 == .(r2_9) * ", " * .(p_label_9))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none"
  ) 
peak_ank_mom_reg_fig

# --- Plot 10: Mean Ankle Moment Difference ---
s10 <- summary(mean_ank_mom_model)
r2_10 <- formatC(s10$r.squared, digits = 3, format = "f")
p_10 <- s10$coefficients[2, 4]
p_label_10 <- if (p_10 < 0.001) "p < 0.001" else paste0("p = ", formatC(p_10, digits = 3, format = "f"))

mean_ank_mom_reg_fig <- ggplot(reg_data, aes(x = mean_ank_mom_diff, y = k_lin, color = factor(ParticipantID))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#CC6600", fill = "#CC6600", alpha = 0.15, linewidth = 1, fullrange = TRUE) +
  scale_x_continuous(
    breaks = seq(-1000, 0, by = .05),
    labels = function(x) ifelse(x %% .1 == 0, x, ""),
    limits = c(-0.25, 0)
  ) +
  scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
  coord_cartesian(ylim = c(-25, 35)) +
  labs(
    x = "Difference in Mean Ankle Moments (Heels - Flats)",
    y = (expression(paste(Delta, "Achilles Tendon Stiffness (", Delta, "%)"))),
    subtitle = bquote(R^2 == .(r2_10) * ", " * .(p_label_10))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y.left = element_text(color = "black", size = 26),
    axis.title.x = element_text(color = "black", size = 26, margin = margin(t = 10)),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(-4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none"
  )
mean_ank_mom_reg_fig


failures_top_fig <- plot_grid(
  peak_strain_reg_fig, mean_strain_reg_fig,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

failures_bottom_fig <- plot_grid(
  NULL,strain_impulse_reg_fig,NULL, 
  labels = c("", "C", ""),
  ncol = 3,
  rel_widths = c(0.5,1,0.5),
  align = "hv"
)

failures_combined_fig <- plot_grid(
  failures_top_fig, failures_bottom_fig,
  ncol = 1,
  align = "v"
)

failures_combined_fig

ggsave2("failed_regressions_fig.pdf", failures_combined_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")


successes_combined_fig <- plot_grid(
  peak_Fr_reg_fig, htd_steps_reg_fig,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)
successes_combined_fig

ggsave2("successful_regressions_fig.pdf", successes_combined_fig, path = "/Users/andrewthornton/Documents/WCB2026/Figures")
