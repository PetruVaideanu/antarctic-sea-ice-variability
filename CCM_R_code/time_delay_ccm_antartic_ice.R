# This script is written for and tested with rEDM version 0.7.1
# Do not run with newer versions (e.g., ???1.13) unless adapted
library(rEDM)

# Define the cause and effect variables.
# Use Fig4 or Supp Info Fig S2 to select desired cause and effect variables
# In the figure captions one can find the Embedding Dimensions used

var1 <- 'amo_sic' # effect
var2 <- 'amo_index' # cause
l = 40 # maximum time delay (lag) to explore in both directions

# Set embedding dimensions previously determined for each direction
E1 = 5 # for var1 xmap var2
E2 = 5 # for var2 xmap var1

# Set the delay between time series values used in embedding
tau = 1

# Initialize placeholders (used for plotting)
tp1 <- NULL  # placeholder for optimal lag in first direction
tp2 <- NULL  # placeholder for optimal lag in second direction
lag_found1 = 0  # optionally use this to mark a vertical line at a meaningful lag

# Load the input data frame "variable_time_series" (see the xlsx file) containing both time series
cause_effect_var_df <- variable_time_series

# Create a vector of variable names
vars <- c(paste(var1), paste(var2))

# Generate all combinations of lib_column, target_column, and tp (lag)
params <- expand.grid(lib_column = vars, target_column = vars, tp = -l:l)

# Exclude combinations where cause and effect are the same variable
params <- params[params$lib_column != params$target_column, ]

# Assign the appropriate embedding dimension and delay to each row
params$E <- c(E2, E1)  # var2 xmap var1 uses E2; var1 xmap var2 uses E1
params$tau <- c(tau, tau)

# Run CCM for each parameter set and collect results in a single data frame
output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
  ccm(cause_effect_var_df, 
      E = params$E[i], 
      lib_sizes = NROW(cause_effect_var_df), 
      random_libs = TRUE, 
      num_samples = 100, 
      replace = FALSE, 
      lib_column = params$lib_column[i], 
      target_column = params$target_column[i],
      tp = params$tp[i], 
      tau = params$tau[i])
}))

# Add a direction label to each row for identification
output$direction <- paste(output$lib_column, "xmap", output$target_column)

# Load ggplot2 for visualization
library(ggplot2)

# Optional: separate and smooth if needed - here we select one direction only
output_smoothed <- subset(output, direction == paste(var1, "xmap", var2))
output_smoothed <- output_smoothed[order(output_smoothed$tp), ]

# Plot only the selected mapping direction
time_delay_ccm_fig <- ggplot(output_smoothed, aes(x = tp, y = pmax(0, rho), color = direction)) + 
  geom_line(size = 1.2) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
  scale_x_continuous(breaks = seq(min(output_smoothed$tp), max(output_smoothed$tp), 10)) +
  
  # Add vertical reference lines for lag 0 and a custom lag if desired
  geom_vline(xintercept = 0, size = 1) +
  geom_vline(xintercept = lag_found1, size = 1, linetype = "dashed") +
  
  # Set custom color (optional, single color used here)
  scale_color_manual(values = c("red")) +  
  
  # Axis and title formatting
  labs(
    x = "Lag",              # x-axis label
    y = "Cross Map Skill"   # y-axis label
  ) +
  
  # Apply styling for publication-ready figure
  theme(
    axis.text.x = element_text(size = 35, face = "bold"),
    axis.text.y = element_text(size = 35, face = "bold"),
    axis.title.x = element_text(size = 35, face = "bold"),
    axis.title.y = element_text(size = 35, face = "bold"),
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(-0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 2, fill = NA),
    legend.position = c(0.50, 0.89),
    legend.text = element_text(size = 25, face = "bold"),
    legend.background = element_rect(fill = "white"),
    legend.key = element_blank()
  ) + 
  
  guides(color = guide_legend(title = NULL))

# Display the plot
print(time_delay_ccm_fig)
