library(tidyverse)
library(infer)
library(ggplot2)
library(ggExtra)
library(mosaic)
library(cowplot)
library(lmtest)
library(lme4)
#help me

f <- "https://raw.githubusercontent.com/Human-Locomotion-Lab-UT-Austin/GA_Tech_Heels/refs/heads/main/participant_data_sheet_R.csv"
d <- read_csv(f, col_names = TRUE)
head(d)

diff_data <- d |>
  group_by(ParticipantID) |>
  summarize(compliance = Compliance[Time == "Pre"], k_lin = (k_lin[Time == "Post"] - k_lin[Time == "Pre"])/k_lin[Time == "Pre"], at_length = (ATLength[Time == "Post"] - ATLength[Time == "Pre"])/ATLength[Time == "Pre"])


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # strain impulse was significantly lower in heel condition than flat condition


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # not significant


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
(p_val <- mean(abs(perm_diff) >= abs(obs_diff))) # mean strain was significantly lower in heel condition compared to flat condition


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
  filter(Compliance == "User")
nonusers_full <- d |>
  filter(Compliance == "Nonuser")

p1 <- users_full |>
  ggplot(aes(x=Time, y=k_lin)) +
  geom_boxplot(alpha = 0.3, color = "black", fill = "white", linewidth = 0.5, whisker.linewidth = 1) +
  scale_x_discrete(limits = c("Pre", "Post")) +
  geom_point(color="orange", size=3) +
  geom_line(aes(group = ParticipantID), color = "orange", alpha = 1) +
  theme(
    legend.position="none",
    plot.title = element_text(size=20, hjust=0.5)
  ) +
  ggtitle("Users") +
  xlab("Timepoint") +
  ylab(expression(paste("Achilles Tendon Stiffness ", italic("(N/mm)")^2))) +
  ylim(100,400)
p1

p2 <- nonusers_full |>
  ggplot(aes(x=Time, y=k_lin)) +
  geom_boxplot(alpha = 0.3, color = "black", fill = "black", linewidth = 0.5, whisker.linewidth = 1) +
  scale_x_discrete(limits = c("Pre", "Post")) +
  geom_point(color="orange", size=3) +
  geom_line(aes(group = ParticipantID), color = "orange", alpha = 1) +
  theme(
    legend.position="none",
    plot.title = element_text(size=20, hjust=0.5)
  ) +
  ggtitle("Nonusers") +
  xlab("Timepoint") +
  ylab(expression(paste("Achilles Tendon Stiffness ", italic("(N/mm)")^2))) +
  ylim(100,400)
p2

f <- "https://raw.githubusercontent.com/Human-Locomotion-Lab-UT-Austin/GA_Tech_Heels/refs/heads/main/HeelsFlatsComp_R.csv"
heels_flats_comp <- read_csv(f, col_names = TRUE)

p3 <- heels_flats_comp |>
  ggplot(aes(x=Condition, y= strain_impulse)) +
  geom_boxplot(alpha = 0.3, color = "black", fill = "lightblue", linewidth = 0.5, whisker.linewidth = 1) +
  scale_x_discrete(limits = c("Flats", "Heels")) +
  geom_point(color="orange", size=3) +
  theme(
    legend.position="none",
    plot.title = element_text(size=20, hjust = 0.5)
  ) +
  ggtitle("Strain Difference Between Heels and Flats") +
  xlab("Condition") +
  ylab(expression(paste("Mean Achilles Tendon Strain", italic("(%)")))) +
  ylim(0,3)
p3

plot_grid(
  plot_grid(
    p1,
    p2,
    labels = c("A", "B"),
    label_size = 12,
    nrow = 1
  ),
  p3,
  labels = c("", "C"),
  label_size = 12,
  nrow = 2
)
