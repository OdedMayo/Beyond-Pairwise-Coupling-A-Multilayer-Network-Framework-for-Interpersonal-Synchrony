# ==============================================================================
# Dynamic Network Analysis: Network Reconfiguration (Regime Shifts) Pipeline
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(lme4)      # Added for Multilevel Models (LMM)
library(MuMIn)     # Added for Marginal R^2 calculation


# ------------------------------------------------------------------------------
# Step 1: Data Loading & Preprocessing
# ------------------------------------------------------------------------------
dynamic_data <- read_csv("dynamic_similarity_results.csv") %>% filter(network_type == "full")
sync_data    <- read_csv("classical_sync_results.csv")

# Define dyadic synchrony channels (must match column names in classical_sync_results.csv)
sync_cols  <- c("movement_sync", "eda_sync", "face_sync", "head_sync", "hr_sync", "rsa_sync")
delta_cols <- paste0("delta_", sync_cols)

# Compute absolute window-to-window changes (Delta |r|) for mean and individual channels
sync_dynamic <- sync_data %>%
  arrange(session, window_id) %>%
  group_by(session) %>%
  mutate(
    delta_mean_sync = abs(mean_sync - lag(mean_sync)),
    across(all_of(sync_cols), ~ abs(.x - lag(.x)), .names = "delta_{.col}")
  ) %>%
  ungroup()

# Merge network dynamics with delta synchrony metrics
df_task2 <- inner_join(dynamic_data, sync_dynamic, by = c("session", "window_id")) %>% 
  drop_na(network_reconfiguration, delta_mean_sync, all_of(delta_cols))


# ------------------------------------------------------------------------------
# Visualization A: Network Reconfiguration vs. Individual Channel Deltas (Facet Grid)
# ------------------------------------------------------------------------------
df_task2_long <- df_task2 %>%
  pivot_longer(cols = all_of(delta_cols), names_to = "delta_channel", values_to = "delta_value") %>%
  mutate(delta_channel = str_remove(delta_channel, "delta_"))

p_dyn_individual <- ggplot(df_task2_long, aes(x = delta_value, y = network_reconfiguration)) +
  geom_point(alpha = 0.2, color = "darkorange3") +
  geom_smooth(method = "lm", color = "darkblue", se = FALSE, linewidth = 1) +
  facet_wrap(~ delta_channel, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Structural Regime Shifts vs. Delta in Individual Synchrony Channels",
    x = "Absolute Change in Dyadic Synchrony (\u0394 |r|)",
    y = "Network Reconfiguration (Distance W_t-1 to W_t)"
  )

ggsave("RegimeShift_vs_Individual_Deltas.png", plot = p_dyn_individual, width = 10, height = 6)


# ------------------------------------------------------------------------------
# Visualization B: Network Reconfiguration vs. Multivariate 3-Level LMM
# ------------------------------------------------------------------------------

# Step 1: Extract Dyad Identifier from Nested Session Column
df_task2 <- df_task2 %>%
  separate(session, into = c("dyad", "session_num"), sep = "_", remove = FALSE)

# Step 2: Fit 3-Level Linear Mixed-Effects Model (Time Windows nested in Sessions, nested in Dyads)
model_dyn_multi <- lmer(
  network_reconfiguration ~ delta_movement_sync + delta_eda_sync + delta_face_sync + delta_head_sync + delta_hr_sync + delta_rsa_sync + 
    (1 | dyad/session), 
  data = df_task2
)

# Step 3: Extract Marginal R^2 (variance explained strictly by fixed effects)
r2_dyn_multi <- r.squaredGLMM(model_dyn_multi)
marginal_r2 <- round(r2_dyn_multi[1], 3)

# Step 4: Compute Predicted Reconfiguration purely from Fixed Effects
df_task2$predicted_reconfig <- predict(model_dyn_multi, re.form = NA)

# Step 5: Plot Actual vs. Predicted Reconfiguration
p_dyn_multi <- ggplot(df_task2, aes(x = predicted_reconfig, y = network_reconfiguration)) +
  geom_point(alpha = 0.3, color = "darkorange3", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Actual Regime Shifts vs. Value Predicted by Delta Channels",
    subtitle = paste0("3-Level LMM Marginal R² = ", marginal_r2, " (Fixed Effects Only)"),
    x = "Predicted Reconfiguration (from Fixed Effects of 6 Delta Channels)",
    y = "Actual Network Reconfiguration"
  )

ggsave("RegimeShift_vs_Multivariate_Delta.png", plot = p_dyn_multi, width = 8, height = 6)

message("All dynamic network reconfiguration (Regime Shifts) analyses and figures generated successfully!")