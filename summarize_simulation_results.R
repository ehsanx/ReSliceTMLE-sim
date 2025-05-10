install.packages("devtools")
devtools::install_github("ehsanx/ReSliceTMLE")
library(ReSliceTMLE)
library(rsimsum)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(dplyr)
library(tidyr)
library(patchwork)
library(stringr)

ReSliceTMLE::load_dependencies()

sim_data <- read.csv("statin_sim_data_full.csv")

SL.glm.DCDR <- function(...) SL.glm(...)
SL.gam4.DCDR <- function(...) SL.gam(..., deg.gam=4)
SL.gam6.DCDR <- function(...) SL.gam(..., deg.gam=6)
SL.nnet.DCDR <- function(...) SL.nnet(..., size=4)
SL.mean.DCDR <- function(...) SL.mean(...)

library <- c("SL.glm.DCDR", "SL.gam4.DCDR", "SL.gam6.DCDR", "SL.nnet.DCDR", "SL.mean.DCDR")

results <- run_multiple_tmle_variants(
  data = sim_data,
  outcome_var = "Y",
  treatment_var = "statin",
  covariates = c("age", "ldl_log", "diabetes", "risk_score"),
  tmle_variants = c("vanilla", "cvq", "cvq_multiple", "fullcv", "fullcv_multiple", "singlecrossfit", "doublecrossfit"),
  num_repetitions = 100,
  Q.SL.library = library,
  g.SL.library = library,
  sim_id_col = "sim_id",  
  true_value = -0.108,
  parallel_strategy = "none",
  missing_handling = "complete_case",
  save_results = TRUE,
  results_dir = "./results",
  family = "binomial"
)

#####################################
# Reproduce the plots in manuscript #
#####################################

## For Figure 3: Comprehensive comparison of TMLE variants across multiple performance metrics.

# Map the method names
method_mapping <- c(
  "vanilla" = "VanillaTMLE",
  "cvq" = "CVQTMLE",
  "cvq_multiple" = "CVQTMLE_rep",
  "fullcv" = "CVTMLE",
  "fullcv_multiple" = "CVTMLE_rep", 
  "singlecrossfit" = "SCTMLE",
  "doublecrossfit" = "DCTMLE"
)
data <- results$raw_final %>%
  mutate(method = method_mapping[method])

# Set true value and calculate SE
true_value <- -0.1081508
data$true_value <- true_value
data$se <- case_when(
  data$method %in% c("CVQTMLE_rep", "CVTMLE_rep", "DCTMLE", "SCTMLE") ~ sqrt(data$mvd),
  TRUE ~ sqrt(data$vd)
)

# Create rsimsum object
results_sim <- simsum(
  data = data,
  estvarname = "rd",
  se = "se",
  true = true_value,
  methodvar = "method",
  ref = "VanillaTMLE"
)

# Extract summary statistics
summary_data <- summary(results_sim)$summ

# Define statistics to plot
plot_stats <- c("bias", "cover", "becover", "mse", "empse", "modelse", "relerror")

# Human-readable names for statistics for plot titles
human_readable_names <- c(
  "bias" = "Bias",
  "cover" = "Coverage",
  "becover" = "Bias-eliminated Coverage",
  "mse" = "Mean Squared Error",
  "empse" = "Empirical Standard Error",
  "modelse" = "Model-based Standard Error",
  "relerror" = "Relative % Error in SE"
)

# Create a function to generate a plot for each statistic
create_boxplot <- function(stat) {
  # Determine the reference line based on the statistic
  reference_line <- if (stat %in% c("cover", "becover")) 0.95 else 0
  
  # Filter data for the specific statistic
  plot_data <- summary_data[summary_data$stat == stat, ]
  
  # Create the plot
  p <- ggplot(
    plot_data,
    aes(x = method, y = est, ymin = lower, ymax = upper)
  ) +
    geom_hline(yintercept = reference_line, color = "red", linetype = "dashed", size = 0.7) +
    geom_point(size = 2) +
    geom_errorbar(width = 0.2) +
    theme_bw() +
    labs(
      title = paste("Boxplot of", human_readable_names[stat]),
      x = "Method",
      y = human_readable_names[stat]
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.title = element_text(size = 11),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    )
  
  # Adjust y-axis limits for better visualization if needed
  if (stat == "bias") {
    p <- p + ylim(min(plot_data$lower) * 1.1, max(plot_data$upper) * 1.1)
  } else if (stat == "mse") {
    p <- p + scale_y_continuous(labels = scales::scientific)
  }
  
  return(p)
}

# Create plots for all statistics
plot_list <- lapply(plot_stats, create_boxplot)

# Combine plots in a grid - 3x3 grid (with one empty spot)
combined_plot <- plot_grid(
  plotlist = plot_list,
  ncol = 3,
  nrow = 3,
  align = "hv",
  labels = "AUTO"
)

# Add an overall title
title <- ggdraw() +
  draw_label(
    "Comparison of Methods Across Multiple Statistics",
    fontface = "bold",
    size = 14,
    x = 0.5,
    hjust = 0.5
  )

# Combine the title and plot
final_plot <- plot_grid(
  title, combined_plot,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)

## For Figure 4: Performance comparison of different TMLE variants across increasing numbers of replicates. 

# Function to prepare data for simsum analysis
prepare_data_for_simsum <- function(data) {
  # Method mapping
  method_mapping <- list(
    "cvq_multiple" = "CVq",
    "singlecrossfit" = "SC",
    "doublecrossfit" = "DC",
    "fullcv_multiple" = "fullCV",
    "vanilla" = "VanillaTMLE"
  )
  
  # Extract and transform data
  processed_data <- data %>%
    # Apply method mapping
    mutate(method = unlist(method_mapping[method])) %>%
    select(sim_id, iteration, method, rd, vd)
  
  # Add true value and calculate standard error
  processed_data <- processed_data %>%
    mutate(
      true_value = true_value,
      se = sqrt(vd)
    )
  
  return(processed_data)
}

# Function to calculate metrics using simsum for a specific number of replicates
calculate_metrics_simsum <- function(processed_data, n_reps) {
  # For methods with repetitions, subset data to use only n_reps per simulation and calculate median
  rep_methods_data <- processed_data %>%
    filter(method %in% c("CVq", "SC", "DC", "fullCV")) %>%
    group_by(sim_id, method) %>%
    slice_head(n = n_reps) %>%
    # Calculate median estimate and variance for each sim_id
    summarise(
      rd = median(rd),
      vd = median(vd),
      se = sqrt(median(vd)),
      true_value = first(true_value),
      .groups = 'drop'
    )
  
  # For Vanilla TMLE, just use all data since there's only one observation per sim_id
  vanilla_method_data <- processed_data %>%
    filter(method == "VanillaTMLE") %>%
    group_by(sim_id) %>%
    slice_head(n = 1) %>%
    select(sim_id, method, rd, vd, se, true_value)
  
  # Combine the datasets
  data_subset <- bind_rows(rep_methods_data, vanilla_method_data)
  
  # Create simsum object with Vanilla TMLE as reference
  sim_results <- simsum(
    data = data_subset,
    estvarname = "rd",
    se = "se",
    true = "true_value",
    methodvar = "method",
    ref = "VanillaTMLE"  # Set Vanilla TMLE as reference
  )
  
  # Extract summary statistics
  summary_data <- summary(sim_results)$summ
  
  # Reshape to have one row per method with columns for each statistic
  results_wide <- summary_data %>%
    select(method, stat, est) %>%
    pivot_wider(
      names_from = stat,
      values_from = est
    ) %>%
    mutate(replicates = n_reps)
  
  return(results_wide)
}

# Function to process all methods and replicates
process_all_replicates <- function(data, max_reps = 100) {
  # Prepare data
  processed_data <- prepare_data_for_simsum(data)
  
  # Initialize results dataframe
  all_results <- data.frame()
  
  # Process each number of replicates
  for (reps in 1:max_reps) {
    if (reps <= max(processed_data$iteration, na.rm = TRUE)) {
      # Calculate metrics for the current number of replicates
      rep_results <- calculate_metrics_simsum(processed_data, reps)
      
      # Add to overall results
      all_results <- rbind(all_results, rep_results)
    }
  }
  
  return(all_results)
}

# Function to create individual plots for each method
create_method_plots <- function(results, method_name) {
  method_data <- results %>% filter(method == method_name)
  
  # Define colorblind-friendly colors for individual method plots
  color_map <- c(
    "CVq" = "#0072B2",     # Blue
    "DC" = "#E69F00",      # Orange
    "fullCV" = "#56B4E9",  # Light blue
    "SC" = "#9e70e4"       # Purple
  )
  
  # Get the color for this method
  method_color <- color_map[method_name]
  
  p1 <- ggplot(method_data, aes(x = replicates, y = bias)) +
    geom_line(color = method_color, linewidth = 1) +
    geom_point(color = method_color) +
    labs(title = paste(method_name, "- Bias"),
         x = "Number of Replicates",
         y = "Absolute Bias") +
    theme_minimal()
  
  p2 <- ggplot(method_data, aes(x = replicates, y = mse)) +
    geom_line(color = method_color, linewidth = 1) +
    geom_point(color = method_color) +
    labs(title = paste(method_name, "- MSE"),
         x = "Number of Replicates",
         y = "MSE") +
    theme_minimal()
  
  p3 <- ggplot(method_data, aes(x = replicates, y = relerror)) +
    geom_line(color = method_color, linewidth = 1) +
    geom_point(color = method_color) +
    labs(title = paste(method_name, "- Relative Error"),
         x = "Number of Replicates",
         y = "Relative Error (%)") +
    theme_minimal()
  
  p4 <- ggplot(method_data, aes(x = replicates, y = cover)) +
    geom_line(color = method_color, linewidth = 1) +
    geom_point(color = method_color) +
    labs(title = paste(method_name, "- Coverage"),
         x = "Number of Replicates",
         y = "Coverage Probability") +
    ylim(0.9, 1.0) +
    theme_minimal()
  
  p5 <- ggplot(method_data, aes(x = replicates, y = becover)) +
    geom_line(color = method_color, linewidth = 1) +
    geom_point(color = method_color) +
    labs(title = paste(method_name, "- BE-Coverage"),
         x = "Number of Replicates",
         y = "BE-Coverage") +
    ylim(0.9, 1.0) +
    theme_minimal()
  
  return((p1 + p2 + p3) / (p4 + p5 + plot_spacer()))
}

# Function to create combined plots for all methods
create_combined_plots <- function(results) {
  # Add a constant "VanillaTMLE" reference line
  vanilla_ref <- results %>%
    filter(method == "VanillaTMLE") %>%
    select(-replicates) %>%
    distinct()
  
  # Prepare data for plotting
  plot_data <- results %>%
    filter(method != "VanillaTMLE")  # Remove VanillaTMLE from line plots
  
  # Define a colorblind-friendly palette
  # Using a blue-orange-purple-teal palette that's distinguishable for most types of color blindness
  colorblind_palette <- c("CVq" = "#0072B2",    # Blue
                          "DC" = "#E69F00",      # Orange
                          "fullCV" = "#56B4E9",  # Light blue
                          "SC" = "#9e70e4")      # Purple
  
  vanilla_color <- "#D55E00"  # Brick red - distinct from the other colors
  
  # Create list of plots
  plot_list <- list()
  
  # Bias plot
  plot_list$bias <- ggplot(plot_data, aes(x = replicates, y = bias, color = method)) +
    geom_line(linewidth = 1) +
    geom_point() +
    geom_hline(yintercept = vanilla_ref$bias, linetype = "dashed", color = vanilla_color, linewidth = 1) +
    annotate("text", x = max(plot_data$replicates) * 0.8, y = vanilla_ref$bias * 1.1,
             label = "VanillaTMLE", color = vanilla_color) +
    scale_color_manual(values = colorblind_palette) +
    labs(title = "Bias Comparison",
         x = "Number of Replicates",
         y = "Absolute Bias") +
    theme_minimal()
  
  # MSE plot
  plot_list$mse <- ggplot(plot_data, aes(x = replicates, y = mse, color = method)) +
    geom_line(linewidth = 1) +
    geom_point() +
    geom_hline(yintercept = vanilla_ref$mse, linetype = "dashed", color = vanilla_color, linewidth = 1) +
    annotate("text", x = max(plot_data$replicates) * 0.8, y = vanilla_ref$mse * 1.1,
             label = "VanillaTMLE", color = vanilla_color) +
    scale_color_manual(values = colorblind_palette) +
    labs(title = "MSE Comparison",
         x = "Number of Replicates",
         y = "MSE") +
    theme_minimal()
  
  # Relative Error plot
  plot_list$rel_error <- ggplot(plot_data, aes(x = replicates, y = relerror, color = method)) +
    geom_line(linewidth = 1) +
    geom_point() +
    geom_hline(yintercept = vanilla_ref$relerror, linetype = "dashed", color = vanilla_color, linewidth = 1) +
    annotate("text", x = max(plot_data$replicates) * 0.8, y = vanilla_ref$relerror * 1.1,
             label = "VanillaTMLE", color = vanilla_color) +
    scale_color_manual(values = colorblind_palette) +
    labs(title = "Relative Error Comparison",
         x = "Number of Replicates",
         y = "Relative Error (%)") +
    theme_minimal()
  
  # Coverage plot
  plot_list$coverage <- ggplot(plot_data, aes(x = replicates, y = cover, color = method)) +
    geom_line(linewidth = 1) +
    geom_point() +
    geom_hline(yintercept = vanilla_ref$cover, linetype = "dashed", color = vanilla_color, linewidth = 1) +
    annotate("text", x = max(plot_data$replicates) * 0.8, y = min(vanilla_ref$cover * 0.99, 0.99),
             label = "VanillaTMLE", color = vanilla_color) +
    scale_color_manual(values = colorblind_palette) +
    labs(title = "Coverage Comparison",
         x = "Number of Replicates",
         y = "Coverage Probability") +
    ylim(0.9, 1.0) +
    theme_minimal()
  
  # BE-Coverage plot
  plot_list$be_coverage <- ggplot(plot_data, aes(x = replicates, y = becover, color = method)) +
    geom_line(linewidth = 1) +
    geom_point() +
    geom_hline(yintercept = vanilla_ref$becover, linetype = "dashed", color = vanilla_color, linewidth = 1) +
    annotate("text", x = max(plot_data$replicates) * 0.8, y = min(vanilla_ref$becover * 0.99, 0.99),
             label = "VanillaTMLE", color = vanilla_color) +
    scale_color_manual(values = colorblind_palette) +
    labs(title = "BE-Coverage Comparison",
         x = "Number of Replicates",
         y = "BE-Coverage") +
    ylim(0.9, 1.0) +
    theme_minimal()
  
  return((plot_list$bias + plot_list$mse + plot_list$rel_error) /
           (plot_list$coverage + plot_list$be_coverage + plot_spacer()))
}

# Function to create a comparative table of final results
create_comparison_table <- function(results, max_reps) {
  # Get results for the maximum number of replicates for comparison
  final_results <- results %>%
    filter(replicates == max_reps | method == "VanillaTMLE") %>%
    select(method, bias, mse, relerror, cover, becover) %>%
    arrange(mse)  # Sort by MSE for easy comparison
  
  return(final_results)
}

### Main execution
data <- results$raw_intermediate

# Extract vanilla data from raw_final
data_vanilla <- results$raw_final %>% 
  filter(method == "vanilla") %>%
  # Add iteration column to match structure
  mutate(iteration = 1) %>%
  # Make sure column names match
  select(r1, r0, rd, v1, v0, vd, method, sim_id, iteration)

# Combine datasets
combined_data <- bind_rows(
  data,
  data_vanilla
)

# Find the maximum number of iterations in the data
max_iterations <- max(data$iteration, na.rm = TRUE)

# Calculate metrics for all methods using the simsum approach
results_analysis <- process_all_replicates(combined_data, max_reps = max_iterations)

# Create individual method plots
cvq_plots <- create_method_plots(results_analysis, "CVq")
sc_plots <- create_method_plots(results_analysis, "SC")
dc_plots <- create_method_plots(results_analysis, "DC")
fullcv_plots <- create_method_plots(results_analysis, "fullCV")
vanilla_stats <- results_analysis %>% filter(method == "VanillaTMLE") %>% distinct()

# Create combined plots
combined_plots <- create_combined_plots(results_analysis)

