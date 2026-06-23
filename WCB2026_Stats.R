library(tidyverse)
library(infer)
library(ggplot2)
library(mosaic)

f <- "https://raw.githubusercontent.com/Human-Locomotion-Lab-UT-Austin/GA_Tech_Heels/refs/heads/main/participant_data_sheet_R.csv"
d <- read_csv(f, col_names = TRUE)
head(d)

diff_data <- d |>
  group_by(ParticipantID) |>
  summarize(compliance = Compliance[Time == "Pre"], k_lin = k_lin[Time == "Post"] - k_lin[Time == "Pre"], at_length = ATLength[Time == "Post"] - ATLength[Time == "Pre"])


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