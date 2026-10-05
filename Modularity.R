# ==============================================================================
# Script: Bivariate Synchrony vs. Network Modularity (3-Level LMM Analysis)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries and Setup Environment
# ------------------------------------------------------------------------------
library(lme4)      # For Multilevel Linear Mixed-Effects Models (LMM)
library(MuMIn)     # For Nakagawa & Schielzeth's Marginal and Conditional R^2
library(tidyverse) # For data manipulation (dplyr, tidyr) and visualization (ggplot2)
library(ggplot2)   # Plotting tools
library(qgraph)    # Network plotting tools

# Set working directory

# ------------------------------------------------------------------------------
# 2. Data Ingestion and Merging
# ------------------------------------------------------------------------------
# Load macro-level network metrics (filtered for the 12-node full system)
net_data  <- read_csv("network_metrics_results.csv") %>% 
  filter(network_type == "full")

# Load classical dyadic/bivariate synchrony results
sync_data <- read_csv("classical_sync_results.csv")

# Merge datasets by session and window_id, removing missing cases
df <- inner_join(net_data, sync_data, by = c("session", "window_id")) %>% 
  drop_na()

# ==============================================================================
# Option A: Per-Channel Bivariate Scatter Plots (Panel View)
# ==============================================================================

# Reshape dataset to long format for faceted channel plotting
df_long <- df %>%
  pivot_longer(
    cols = c(movement_sync, head_sync, face_sync, eda_sync, hr_sync, rsa_sync),
    names_to = "channel",
    values_to = "pairwise_sync"
  ) %>%
  mutate(channel = case_when(
    channel == "movement_sync" ~ "Body Movement",
    channel == "head_sync"      ~ "Head Motion",
    channel == "face_sync"      ~ "Facial Expression",
    channel == "eda_sync"       ~ "EDA (Electrodermal)",
    channel == "hr_sync"        ~ "Heart Rate",
    channel == "rsa_sync"       ~ "RSA (Respiratory)",
    TRUE ~ channel
  ))

# Plot Modularity vs. Each Bivariate Synchrony Modality
p_channels <- ggplot(df_long, aes(x = pairwise_sync, y = modularity)) +
  geom_point(alpha = 0.2, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkred", se = FALSE) +
  facet_wrap(~channel, scales = "free_x", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Network Modularity vs. Individual Pairwise Synchrony Channels",
    x = "Pairwise Synchrony (|r|)",
    y = "Network Modularity"
  )

# Save faceted plot
ggsave("Modularity_vs_Individual_Channels.png", plot = p_channels, width = 10, height = 6)

# ==============================================================================
# Option B: Multivariate 3-Level Linear Mixed-Effects Model (LMM)
# Predicted vs. Actual Modularity & Residual Analysis
# ==============================================================================

# Step 1: Extract Dyad Identifier from Nested Session Column (e.g., "10_11" -> Dyad 10, Session 11)
df <- df %>%
  separate(session, into = c("dyad", "session_num"), sep = "_", remove = FALSE)

# Step 2: Fit 3-Level Linear Mixed-Effects Model (Time Windows nested in Sessions, nested in Dyads)
# Fixed Effects: All 6 pairwise synchrony channels
# Random Effects: Random intercepts for sessions nested within dyads (1 | dyad/session)
model_full_mixed <- lmer(
  modularity ~ movement_sync + head_sync + face_sync + eda_sync + hr_sync + rsa_sync + 
    (1 | dyad/session), 
  data = df
)

# Step 3: Extract Marginal R^2 (variance explained strictly by fixed effects / pairwise synchrony)
r2_mixed <- r.squaredGLMM(model_full_mixed)
marginal_r2 <- round(r2_mixed[1], 3) # Marginal R^2 isolating fixed effects

# Step 4: Compute Predicted Modularity derived purely from Fixed Effects (ignoring random baseline variation)
df$predicted_modularity_mixed <- predict(model_full_mixed, re.form = NA)

# Step 5: Plot Actual vs. Predicted Modularity relative to the Line of Identity (y = x)
p_pred_vs_act_mixed <- ggplot(df, aes(x = predicted_modularity_mixed, y = modularity)) +
  geom_point(alpha = 0.3, color = "midnightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Actual Modularity vs. Value Predicted by Bivariate Synchrony",
    subtitle = paste0("3-Level LMM Marginal R² = ", marginal_r2, 
                      " (Fixed Effects Only)"),
    x = "Predicted Modularity (from Fixed Effects of 6 Channels)",
    y = "Actual Network Modularity"
  )

# Save 3-Level LMM validation plot
ggsave("Predicted_vs_Actual_Modularity_3Level.png", plot = p_pred_vs_act_mixed, width = 8, height = 6)
