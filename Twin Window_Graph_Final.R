# ==============================================================================
# Extraction & Visual Comparison of Twin Network Windows Pipeline
# ==============================================================================

library(tidyverse)
library(igraph)
library(qgraph) # Specialized package for psychological/physiological network visualization

# Set working directory

# ------------------------------------------------------------------------------
# Step 1: Data Loading & Twin Windows Selection
# ------------------------------------------------------------------------------
# Load full time-series dataset
full_dataset <- read_csv("full_dataset_merged.csv")

# Load identified twin window pairs and select top-ranked pair
twin_pairs <- read_csv("twin_windows_multivariate.csv")
top_pair   <- twin_pairs %>% slice(1)

sess_A <- top_pair$session_A
win_A  <- top_pair$window_id_A
sess_B <- top_pair$session_B
win_B  <- top_pair$window_id_B

# ------------------------------------------------------------------------------
# Step 2: Define Network Nodes & Display Labels
# ------------------------------------------------------------------------------
client_nodes    <- c("BodyC", "HeadC", "C_Valence", "C_Mean SC", "C_HR", "C_RSA")
therapist_nodes <- c("BodyT", "HeadT", "T_Valence", "T_Mean SC", "T_HR", "T_RSA")
nodes_full      <- c(client_nodes, therapist_nodes)

# Human-readable labels for network plots
display_labels  <- c(
  "C_Body", "C_Head", "C_Face", "C_EDA", "C_HR", "C_RSA",  # Client
  "T_Body", "T_Head", "T_Face", "T_EDA", "T_HR", "T_RSA"   # Therapist
)

# ------------------------------------------------------------------------------
# Step 3: Function to Extract Adjacency Correlation Matrix
# ------------------------------------------------------------------------------
extract_window_matrix <- function(dataset, target_session, target_window_id, window_length = 30, step = 15) {
  sess_df <- dataset %>% filter(session == target_session)
  min_t   <- min(sess_df$second, na.rm = TRUE)
  start_t <- min_t + (target_window_id - 1) * step
  end_t   <- start_t + window_length
  
  win_data <- sess_df %>%
    filter(second >= start_t & second < end_t) %>%
    select(any_of(nodes_full))
  
  cor_m <- abs(cor(win_data, use = "pairwise.complete.obs"))
  cor_m[is.na(cor_m)] <- 0
  diag(cor_m) <- 0
  return(cor_m)
}

# Generate adjacency matrices for Window A and Window B
adj_A <- extract_window_matrix(full_dataset, sess_A, win_A)
adj_B <- extract_window_matrix(full_dataset, sess_B, win_B)

# ------------------------------------------------------------------------------
# Step 4: Extract Pairwise Synchrony Metrics for Figure Subtitle
# ------------------------------------------------------------------------------
get_pair_val <- function(m, n1, n2) {
  if (n1 %in% rownames(m) && n2 %in% colnames(m)) round(m[n1, n2], 2) else NA
}

body_v <- get_pair_val(adj_A, "BodyC", "BodyT")
head_v <- get_pair_val(adj_A, "HeadC", "HeadT")
face_v <- get_pair_val(adj_A, "C_Valence", "T_Valence")
eda_v  <- get_pair_val(adj_A, "C_Mean SC", "T_Mean SC")
hr_v   <- get_pair_val(adj_A, "C_HR", "T_HR")
rsa_v  <- get_pair_val(adj_A, "C_RSA", "T_RSA")

channel_str <- paste0("Body=", body_v, " | Head=", head_v, " | Face=", face_v, 
                      " | EDA=", eda_v, " | HR=", hr_v, " | RSA=", rsa_v)

# ------------------------------------------------------------------------------
# Step 5: High-Resolution Network Visualization
# ------------------------------------------------------------------------------
edge_threshold <- 0.25 # Edge filtering threshold for display clarity

png("Twin_Networks_Clean_Labels.png", width = 3600, height = 2000, res = 300)

# Outer margin specification (oma) reserves top area for main figure titles
par(mfrow = c(1, 2), oma = c(1, 1, 9, 1))

groups_list <- list(
  Client    = 1:6,
  Therapist = 7:12
)
group_colors <- c("#3498DB", "#E74C3C")

# --- Plot Network A ---
qgraph(
  adj_A,
  labels      = display_labels, 
  threshold   = edge_threshold,
  layout      = "spring",
  groups      = groups_list,
  color       = group_colors,
  vsize       = 9,
  label.cex   = 1.1,
  edge.width  = 3.5,
  maximum     = 1,
  mar         = c(3, 3, 10, 3), # Inner plot margins (Bottom, Left, Top, Right) to prevent label clipping
  title       = paste0("Window A (Sess ", sess_A, " | Win ", win_A, ")\nModularity Q = ", 
                       round(top_pair$modularity_A, 3), "\nFiltered Edges (|r| >= ", edge_threshold, ")"),
  title.cex   = 1.1
)

# --- Plot Network B ---
qgraph(
  adj_B,
  labels      = display_labels, 
  threshold   = edge_threshold,
  layout      = "spring",
  groups      = groups_list,
  color       = group_colors,
  vsize       = 9,
  label.cex   = 1.1,
  edge.width  = 3.5,
  maximum     = 1,
  mar         = c(3, 3, 10, 3), # Inner plot margins (Bottom, Left, Top, Right) to prevent label clipping
  title       = paste0("Window B (Sess ", sess_B, " | Win ", win_B, ")\nModularity Q = ", 
                       round(top_pair$modularity_B, 3), "\nFiltered Edges (|r| >= ", edge_threshold, ")"),
  title.cex   = 1.1
)

# Outer margin main headers
mtext("Network Structure Comparison Across Twin Windows", 
      side = 3, line = 5, outer = TRUE, cex = 1.5, font = 2)
mtext(paste0("Identical Pairwise Synchrony Channels (|r|): ", channel_str), 
      side = 3, line = 2, outer = TRUE, cex = 1.1, font = 1)

dev.off()

