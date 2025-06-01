# Import the rEDM version 0.7.1 (important for compatibility)
library(rEDM)
library(ggplot2)
library(doParallel)
library(foreach)

# Define effect and cause variables
# Use Fig4 or Supp Info Fig S2 to select desired cause and effect variables
# In the figure captions one can find the Embedding Dimensions used

var1 <- 'amo_sic'     # This is the effect variable
var2 <- 'amo_index'   # This is the hypothesized cause

# Set up CCM parameters
num_surr = 100        # Number of surrogate series to generate
E = 5                 # Embedding dimension
tp = -8               # Prediction lag
tau = 1               # Embedding time delay

# Provide the input dataframe that contains both time series
cause_effect_var_df <- variable_time_series

# Extract the effect variable time series for CCM analysis
effect_surr <- cause_effect_var_df[[var1]]

# Generate surrogate time series for the cause variable
# Method 1: Ebisuzaki (Fourier-based)
cause_surr_ebs <- make_surrogate_data(cause_effect_var_df[[var2]], method = "ebisuzaki", num_surr = num_surr)

# Method 2: Block Bootstrap model (resampling with blocks)
cause_surr_btstr <- block_bootstrap_model(cause_effect_var_df[[var2]], nr_surr = num_surr, block_size = 10)

# Combine each set of surrogate cause series with the original effect series
surrogates_ebs <- data.frame(cause_surr_ebs, effect_surr)
surrogates_btstr <- data.frame(cause_surr_btstr, effect_surr)

# Perform CCM at full library size for each surrogate series (used to select top 5%)
ccms_surr_ebs <- do.call(rbind, lapply(seq_len(NCOL(surrogates_ebs)-1), function(i) {
  ccm(surrogates_ebs, E = E, lib_sizes = NROW(cause_effect_var_df), 
      random_libs = TRUE, num_samples = 100, replace = FALSE, 
      lib_column = num_surr + 1, target_column = i, tp = tp, 
      silent = TRUE, tau = tau)
}))

ccms_surr_btstr <- do.call(rbind, lapply(seq_len(NCOL(surrogates_btstr)-1), function(i) {
  ccm(surrogates_btstr, E = E, lib_sizes = NROW(cause_effect_var_df), 
      random_libs = TRUE, num_samples = 100, replace = FALSE, 
      lib_column = num_surr + 1, target_column = i, tp = tp, 
      silent = TRUE, tau = tau)
}))

# Select the top 5% most predictive surrogates based on CCM rho
five_percent <- 0.05 * num_surr
toptensurr_ebs <- order(-ccms_surr_ebs$rho)[five_percent]
toptensurr_btstr <- order(-ccms_surr_btstr$rho)[five_percent]

# Initialize parallel processing
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Compute CCM curves (rho vs library size) for selected top Ebisuzaki surrogates
ccm_for_surrogates_ebs <- foreach(lib_size = seq(2, NROW(cause_effect_var_df), 2), .combine = rbind, .packages = 'rEDM') %dopar% {
  ccm(surrogates_ebs, E = E, lib_sizes = lib_size, random_libs = TRUE, 
      num_samples = 100, replace = FALSE, 
      lib_column = num_surr + 1, target_column = toptensurr_ebs, 
      tp = tp, silent = TRUE, tau = tau)
}
ccm_surrogates_ebs <- ccm_means(ccm_for_surrogates_ebs)

# Compute CCM curves for selected top Bootstrap surrogates
ccm_for_surrogates_btstr <- foreach(lib_size = seq(2, NROW(cause_effect_var_df), 2), .combine = rbind, .packages = 'rEDM') %dopar% {
  ccm(surrogates_btstr, E = E, lib_sizes = lib_size, random_libs = TRUE, 
      num_samples = 100, replace = FALSE, 
      lib_column = num_surr + 1, target_column = toptensurr_btstr, 
      tp = tp, silent = TRUE, tau = tau)
}
ccm_surrogates_btstr <- ccm_means(ccm_for_surrogates_btstr)

# Compute the actual CCM curve between real variables (not surrogate)
sp_xmap_pdo <- foreach(lib_size = seq(2, NROW(cause_effect_var_df), 2), .combine = rbind, .packages = 'rEDM') %dopar% {
  ccm(cause_effect_var_df, E = E, lib_column = var1, target_column = var2, 
      tp = tp, lib_sizes = lib_size, random_libs = TRUE, 
      num_samples = 100, replace = FALSE, silent = FALSE, tau = tau)
}
s_xmap_p_means <- ccm_means(sp_xmap_pdo)

# Shut down the parallel backend
stopCluster(cl)

# Plotting the CCM skill curves for real and surrogate cases

# Create a combined data frame for plotting
plot_data <- data.frame(
  lib_size = rep(s_xmap_p_means$lib_size, 3),
  rho = c(
    pmax(0, s_xmap_p_means$rho),          # Actual data
    pmax(0, ccm_surrogates_ebs$rho),      # Ebisuzaki surrogates
    pmax(0, ccm_surrogates_btstr$rho)     # Bootstrap surrogates
  ),
  group = rep(c(paste(var1, "xmap", var2), "Ebisuzaki", "Bootstrap"), 
              each = length(s_xmap_p_means$lib_size))
)

# Factor the group column to control legend order
plot_data$group <- factor(plot_data$group, levels = c(paste(var1, "xmap", var2), "Ebisuzaki", "Bootstrap"))

# Define color and size for each line
color_mapping <- setNames(c("red", "gray", "darkgray"), 
                          c(paste(var1, "xmap", var2), "Ebisuzaki", "Bootstrap"))
size_mapping <- setNames(c(2.5, 1.5, 1.5), 
                         c(paste(var1, "xmap", var2), "Ebisuzaki", "Bootstrap"))

# Generate the plot
ggplot(plot_data, aes(x = lib_size, y = rho, color = group, size = group)) +
  
  # Shaded regions under surrogate lines
  geom_ribbon(data = subset(plot_data, group == "Ebisuzaki"), 
              aes(x = lib_size, ymin = 0, ymax = rho, fill = "Ebisuzaki"), 
              inherit.aes = FALSE, alpha = 1) +
  geom_ribbon(data = subset(plot_data, group == "Bootstrap"), 
              aes(x = lib_size, ymin = 0, ymax = rho, fill = "Bootstrap"), 
              inherit.aes = FALSE, alpha = 1) +
  
  geom_line() +  # Add line plots
  theme_bw() +   # Use white background
  
  # Apply manual color and fill settings
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = c("Ebisuzaki" = "gray", "Bootstrap" = "darkgray"), guide = "none") +
  scale_size_manual(values = size_mapping) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  
  # Customize plot theme
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", margin = margin(t = 0, b = 0)),
    axis.text.y = element_text(size = 35, face = "bold", margin = margin(l = 0, r = 0)),
    axis.title.x = element_text(size = 35, face = "bold", margin = margin(t = 0, b = 0)),
    axis.title.y = element_text(size = 35, face = "bold", margin = margin(l = 0, r = 0)),
    axis.ticks = element_line(size = 1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 2, fill = NA),
    legend.position = c(0.5, 0.89),
    legend.text = element_text(size = 25, face = "bold"),
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = "black", size = 1)
  ) +
  
  labs(
    x = "Library Size",
    y = "Cross Map Skill",
    color = "Legend",
    size = "Legend"
  ) +
  
  guides(
    color = guide_legend(override.aes = list(size = c(2.5, 1.5, 1.5)))
  )

