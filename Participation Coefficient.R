# ==============================================================================
# Dynamic Network Analysis: Node-Level Participation Coefficient (PC) Pipeline
# ==============================================================================

library(tidyverse)
library(igraph)
library(ggplot2)
library(lme4)      # Added for Multilevel Models (LMM)
library(MuMIn)     # Added for Marginal R^2 calculation



# ------------------------------------------------------------------------------
# Step 1: Data Loading and Channel Definition
# ------------------------------------------------------------------------------
full_dataset <- read_csv("full_dataset_merged.csv")
sync_df      <- read_csv("classical_sync_results.csv")

client_nodes    <- c("BodyC", "HeadC", "C_Valence", "C_Mean SC", "C_HR", "C_RSA")
therapist_nodes <- c("BodyT", "HeadT", "T_Valence", "T_Mean SC", "T_HR", "T_RSA")
nodes_full      <- c(client_nodes, therapist_nodes)

# Helper function to compute Participation Coefficient (standalone, package-independent)
calculate_pc <- function(g, memb) {
  vs <- V(g)
  pc <- numeric(length(vs))
  names(pc) <- names(vs)
  str <- strength(g, weights = E(g)$weight)
  
  for(i in seq_along(vs)) {
    v <- vs[i]
    if(str[i] == 0) {
      pc[i] <- 0
      next
    }
    neigh <- neighbors(g, v)
    if(length(neigh) == 0) { pc[i] <- 0; next }
    
    e_weights <- E(g)[v %--% neigh]$weight
    neigh_comms <- memb[names(neigh)]
    
    # Sum edge weights by community to quantify connectivity distribution
    df <- aggregate(e_weights ~ neigh_comms, FUN = sum)
    pc[i] <- 1 - sum((df$e_weights / str[i])^2)
  }
  return(pc)
}

# ------------------------------------------------------------------------------
# Step 2: Sliding-Window Processing Loop for Node-Level Metrics
# ------------------------------------------------------------------------------
window_length     <- 30    
step              <- 15    
min_valid_samples <- 15 
node_metrics_list <- list()

all_sessions <- unique(full_dataset$session)

for (sess in all_sessions) {
  session_data <- full_dataset %>% filter(session == sess)
  max_time     <- max(session_data$second, na.rm = TRUE)
  min_time     <- min(session_data$second, na.rm = TRUE)
  
  start_t   <- min_time
  window_id <- 1
  
  while (start_t + window_length <= max_time) {
    end_t <- start_t + window_length
    win_data <- session_data %>% filter(second >= start_t & second < end_t)
    
    actual_samples <- nrow(win_data)
    if (actual_samples >= min_valid_samples) {
      available_nodes <- intersect(nodes_full, colnames(win_data))
      
      # Quality filter: exclude flatlines and channels with >50% missing data
      valid_nodes <- c()
      for (n in available_nodes) {
        node_data <- win_data[[n]]
        if (sum(is.na(node_data))/length(node_data) < 0.5 && sd(node_data, na.rm=TRUE) > 0) {
          valid_nodes <- c(valid_nodes, n)
        }
      }
      
      if (length(valid_nodes) > 2) {
        win_filtered <- win_data %>% select(all_of(valid_nodes))
        cor_matrix <- abs(cor(win_filtered, use = "pairwise.complete.obs"))
        cor_matrix[is.na(cor_matrix)] <- 0
        diag(cor_matrix) <- 0
        
        g <- graph_from_adjacency_matrix(cor_matrix, mode = "undirected", weighted = TRUE)
        
        # Algorithmic community detection (Louvain method)
        comm <- cluster_louvain(g, weights = E(g)$weight)
        memb <- membership(comm)
        
        pc_vals <- calculate_pc(g, memb)
        
        for (v_name in names(pc_vals)) {
          node_metrics_list[[length(node_metrics_list) + 1]] <- tibble(
            session = sess,
            window_id = window_id,
            node = v_name,
            participation_coeff = pc_vals[v_name],
            community = memb[v_name]
          )
        }
      }
    }
    start_t <- start_t + step
    window_id <- window_id + 1
  }
}

# Export node-level results to CSV
node_metrics_df <- bind_rows(node_metrics_list)
write_csv(node_metrics_df, "node_metrics_results.csv")

# ------------------------------------------------------------------------------
# Step 3: Merging Node Metrics with Dyadic Synchrony & Generating Figures
# ------------------------------------------------------------------------------

# Define dyadic synchrony channel names (must match classical_sync_results.csv)
sync_cols <- c("movement_sync", "eda_sync", "face_sync", "head_sync", "hr_sync", "rsa_sync")

# Focus analysis on Client RSA node (C_RSA)
rsa_node_df <- node_metrics_df %>% filter(node == "C_RSA")

# Merge node metrics with dyadic synchrony data
df_merged <- inner_join(rsa_node_df, sync_df, by = c("session", "window_id")) %>% 
  drop_na(participation_coeff, all_of(sync_cols))


# ------------------------------------------------------------------------------
# Visualization A: RSA Participation Coefficient vs. Direct RSA Dyadic Synchrony
# ------------------------------------------------------------------------------
model_node_single <- lm(participation_coeff ~ rsa_sync, data = df_merged)
r2_node_single <- round(summary(model_node_single)$r.squared, 3)

p_node_single <- ggplot(df_merged, aes(x = rsa_sync, y = participation_coeff)) +
  geom_point(alpha = 0.4, color = "seagreen", size = 2) +
  geom_smooth(method = "lm", color = "darkred", se = FALSE, linetype = "dashed", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Role in Network vs. Dyadic Synchrony (Client RSA)",
    x = "Dyadic RSA Synchrony (|r|)",
    y = "RSA Participation Coefficient (Network Hubness)"
  )

ggsave("Node_PC_vs_Dyadic_Sync.png", plot = p_node_single, width = 7, height = 6)


# ------------------------------------------------------------------------------
# Visualization B: RSA Participation Coefficient vs. Individual Synchrony Channels
# ------------------------------------------------------------------------------
df_long <- df_merged %>%
  pivot_longer(cols = all_of(sync_cols), names_to = "channel", values_to = "sync_value")

p_node_individual <- ggplot(df_long, aes(x = sync_value, y = participation_coeff)) +
  geom_point(alpha = 0.2, color = "seagreen") +
  geom_smooth(method = "lm", color = "darkred", se = FALSE, linewidth = 1) +
  facet_wrap(~ channel, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Client RSA Hubness vs. Individual Synchrony Channels",
    x = "Dyadic Synchrony (|r|)",
    y = "RSA Participation Coefficient"
  )

ggsave("Node_PC_vs_Individual_Channels.png", plot = p_node_individual, width = 10, height = 6)


# ------------------------------------------------------------------------------
# Visualization C: RSA Participation Coefficient vs. Multivariate 3-Level LMM
# ------------------------------------------------------------------------------

# Step 1: Extract Dyad Identifier from Nested Session Column
df_merged <- df_merged %>%
  separate(session, into = c("dyad", "session_num"), sep = "_", remove = FALSE)

# Step 2: Fit 3-Level Linear Mixed-Effects Model (Time Windows nested in Sessions, nested in Dyads)
model_node_multi <- lmer(
  participation_coeff ~ movement_sync + eda_sync + face_sync + head_sync + hr_sync + rsa_sync + 
    (1 | dyad/session), 
  data = df_merged
)

# Step 3: Extract Marginal R^2 (variance explained strictly by fixed effects)
r2_node_multi <- r.squaredGLMM(model_node_multi)
marginal_r2 <- round(r2_node_multi[1], 3)

# Step 4: Compute Predicted Participation Coefficient purely from Fixed Effects
df_merged$predicted_pc <- predict(model_node_multi, re.form = NA)

# Step 5: Plot Actual vs. Predicted PC
p_node_multi <- ggplot(df_merged, aes(x = predicted_pc, y = participation_coeff)) +
  geom_point(alpha = 0.3, color = "seagreen", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed", linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Actual RSA Hubness vs. Value Predicted by Bivariate Synchrony",
    subtitle = paste0("3-Level LMM Marginal R² = ", marginal_r2, " (Fixed Effects Only)"),
    x = "Predicted Participation Coefficient (from Fixed Effects of 6 Channels)",
    y = "Actual RSA Participation Coefficient"
  )

ggsave("Node_PC_vs_Multivariate_Model.png", plot = p_node_multi, width = 8, height = 6)

message("All node-level (Participation Coefficient) analyses and figures generated successfully!")
