# --- Starting Tree Bias Test (RF Distance Shuffle Test) ---
# This script evaluates whether the inferred phylogenetic trees
# are biased toward their corresponding guide trees used in PMSF inference.

# Required packages
library(ape)
library(ggplot2)

# ----------------------------------------------------------
# Example of loading guide trees and result trees:
# (Both must contain the same number of trees, in matching order)
# guide_trees <- read.tree("guide_trees.tre")   # e.g., 20 guide trees
# result_trees <- read.tree("result_trees.tre") # e.g., 20 inferred trees

# ----------------------------------------------------------
# Function: Perform Robinson-Foulds (RF) distance shuffle test
starting_tree_bias_test <- function(guide_trees, result_trees, reps = 999) {
  if (length(guide_trees) != length(result_trees)) {
    stop("guide_trees and result_trees must have the same length")
  }

  # Calculate the observed statistic: sum of RF distances for each guide-result pair
  rf_obs_each <- mapply(function(g, r) RF.dist(g, r), guide_trees, result_trees)
  obs <- sum(rf_obs_each)

  # Perform random shuffling of result_trees 999 times
  reps_dist <- numeric(reps)
  for (i in 1:reps) {
    shuffled <- sample(result_trees, length(result_trees), replace = FALSE)
    rf_rand_each <- mapply(function(g, r) RF.dist(g, r), guide_trees, shuffled)
    reps_dist[i] <- sum(rf_rand_each)
  }

  # One-tailed test: p-value (observed < random)
  position <- 1 + length(which(reps_dist <= obs))
  p.value <- position / (reps + 1)

  # Effect size: difference between observed and random mean
  effect.size <- obs - mean(reps_dist)

  # ----------------------------------------------------------
  # Plot the result distribution
  df <- data.frame(RF_sum = reps_dist)
  p <- ggplot(df, aes(x = RF_sum)) +
    geom_histogram(binwidth = diff(range(df$RF_sum)) / 30,
                   fill = "grey70", color = "black") +
    geom_vline(xintercept = obs, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Starting Tree Bias Test (RF Distance)",
      x = "Sum of RF Distances (Randomized Pairs)",
      y = "Frequency",
      subtitle = paste0(
        "Observed = ", round(obs, 2),
        ", Mean(Random) = ", round(mean(reps_dist), 2),
        ", p = ", signif(p.value, 3)
      )
    ) +
    theme_minimal(base_size = 14)

  # Return results
  return(list(
    obs = obs,
    random_dist = reps_dist,
    p.value = p.value,
    effect.size = effect.size,
    plot = p
  ))
}

# ----------------------------------------------------------
# Example of running:
# test_result <- starting_tree_bias_test(guide_trees, result_trees)
# print(test_result$p.value)
# test_result$plot
