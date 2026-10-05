# ==============================================================================
# Dynamic Network Analysis and Multimodal Synchrony Pipeline
# ==============================================================================

# Load required libraries
library(tidyverse)
library(igraph)
library(brainGraph) # Used for computing Global Efficiency

# ------------------------------------------------------------------------------
# Step 0: Data Loading and Preliminary Quality Checks
# ------------------------------------------------------------------------------

movement_data <- read_csv("MOV_anonymized.csv")
face_data     <- read_csv("FACE_anonymized.csv")
eda_data      <- read_csv("EDA_anonymized.csv")
hr_data       <- read_csv("HR_anonymized.csv")
rsa_data      <- read_csv("RSA_anonymized.csv")

# Validation function for each modality
validate_modality <- function(df, name) {
  required_cols <- c("filterindex", "c_code_anonymized", "StartTime")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if(length(missing_cols) > 0) {
    warning(paste("Warning: In file", name, "the following required columns are missing:", paste(missing_cols, collapse = ", ")))
  } else {
    message(paste("✓ File", name, "is valid and contains all required columns."))
  }
  
  # Print preliminary summary statistics
  df %>% 
    summarise(
      Modality = name,
      Total_Rows = n(),
      Unique_Sessions = n_distinct(filterindex),
      Unique_Patients = n_distinct(c_code_anonymized),
      Min_Time = min(StartTime, na.rm = TRUE),
      Max_Time = max(StartTime, na.rm = TRUE)
    ) %>% 
    print()
  
  # Row counts per session
  message(paste("Number of rows per session in", name, ":"))
  df %>% count(filterindex) %>% print()
  cat("\n-----------------------------------\n")
}

# Execute validation checks for all modalities
validate_modality(movement_data, "Movement")
validate_modality(face_data, "Face")
validate_modality(eda_data, "EDA")
validate_modality(hr_data, "HR")
validate_modality(rsa_data, "RSA")

# ------------------------------------------------------------------------------
# Step 1: Temporal Normalization (Zero-Alignment) and Downsampling (1 Hz)
# ------------------------------------------------------------------------------
downsample_and_zero_align <- function(df) {
  df %>%
    mutate(StartTime = StartTime - min(StartTime, na.rm = TRUE)) %>% 
    mutate(second = floor(StartTime)) %>%
    group_by(filterindex, c_code_anonymized, second) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
    select(-StartTime)
}

movement_1hz <- downsample_and_zero_align(movement_data)
face_1hz     <- downsample_and_zero_align(face_data)
eda_1hz      <- downsample_and_zero_align(eda_data)
hr_1hz       <- downsample_and_zero_align(hr_data)
rsa_1hz      <- downsample_and_zero_align(rsa_data)

# ------------------------------------------------------------------------------
# Step 2: Merging Modalities into a Unified Dataset
# ------------------------------------------------------------------------------
full_dataset <- movement_1hz %>%
  full_join(face_1hz, by = c("filterindex", "c_code_anonymized", "second")) %>%
  full_join(eda_1hz,  by = c("filterindex", "c_code_anonymized", "second")) %>%
  full_join(hr_1hz,   by = c("filterindex", "c_code_anonymized", "second")) %>%
  full_join(rsa_1hz,  by = c("filterindex", "c_code_anonymized", "second")) %>%
  rename(session = filterindex, patient = c_code_anonymized)

print("Column names in full_dataset:")
print(colnames(full_dataset))
write_csv(full_dataset, "full_dataset_merged.csv")

# ------------------------------------------------------------------------------
# Step 3: Define Network Nodes (Client and Therapist Channels)
# ------------------------------------------------------------------------------
client_nodes <- c("BodyC", "HeadC", "C_Valence", "C_Mean SC", "C_HR", "C_RSA")
therapist_nodes <- c("BodyT", "HeadT", "T_Valence", "T_Mean SC", "T_HR", "T_RSA")
nodes_full <- c(client_nodes, therapist_nodes)

# ------------------------------------------------------------------------------
# Step 4: Sliding Window Parameters and Data Containers
# ------------------------------------------------------------------------------
window_length     <- 30    
step              <- 15    
min_valid_samples <- 15 # Minimum valid seconds (samples) required for a valid window

network_metrics_results    <- tibble()
dynamic_similarity_results <- tibble()
classical_sync_results     <- tibble()
window_diagnostics_results <- tibble() # Data frame to track window quality diagnostics

all_sessions <- unique(full_dataset$session)

# ------------------------------------------------------------------------------
# Step 5: Main Sliding-Window Processing Loop
# ------------------------------------------------------------------------------
for (sess in all_sessions) {
  
  session_data <- full_dataset %>% filter(session == sess)
  max_time     <- max(session_data$second, na.rm = TRUE)
  min_time     <- min(session_data$second, na.rm = TRUE)
  
  prev_adj_matrices <- list(full = NULL, client_internal = NULL, therapist_internal = NULL, interpersonal = NULL)
  
  window_id <- 1
  start_t   <- min_time
  
  while (start_t + window_length <= max_time) {
    end_t <- start_t + window_length
    
    window_data <- session_data %>% 
      filter(second >= start_t & second < end_t) %>%
      select(second, any_of(nodes_full)) # Retain 'second' column for temporal continuity check
    
    # Initialize diagnostic tracking variables
    diag_status <- "Valid"
    diag_missing_nodes <- ""
    diag_dropped_flatlines <- ""
    actual_samples <- nrow(window_data)
    max_time_gap <- 0
    
    # Quality Check 1: Insufficient sample count check
    if (actual_samples < min_valid_samples) {
      diag_status <- "Dropped_Insufficient_Samples"
    } else {
      # Quality Check 2: Detection of large temporal gaps
      max_time_gap <- max(diff(sort(window_data$second)), na.rm = TRUE)
      if (max_time_gap > 10) { # Flag consecutive time gaps exceeding 10 seconds
        diag_status <- "Warning_Large_Time_Gap"
      }
      
      available_nodes <- intersect(nodes_full, colnames(window_data))
      missing_nodes <- setdiff(nodes_full, available_nodes)
      if(length(missing_nodes) > 0) diag_missing_nodes <- paste(missing_nodes, collapse = "|")
      
      if (length(available_nodes) > 2) {
        
        # Quality Check 3: Filter channels with zero variance (flatlines) or excessive missingness
        valid_nodes <- c()
        dropped_nodes <- c()
        
        for (node in available_nodes) {
          node_data <- window_data[[node]]
          na_ratio <- sum(is.na(node_data)) / length(node_data)
          node_sd <- sd(node_data, na.rm = TRUE)
          
          # Retain node if < 50% missing values and standard deviation > 0
          if (na_ratio < 0.5 && !is.na(node_sd) && node_sd > 0) {
            valid_nodes <- c(valid_nodes, node)
          } else {
            dropped_nodes <- c(dropped_nodes, node)
          }
        }
        
        if (length(dropped_nodes) > 0) diag_dropped_flatlines <- paste(dropped_nodes, collapse = "|")
        
        if (length(valid_nodes) > 2) {
          window_filtered <- window_data %>% select(all_of(valid_nodes))
          
          # Quality Check 4: Pairwise complete correlation computation
          raw_cor <- suppressWarnings(cor(window_filtered, use = "pairwise.complete.obs"))
          
          # Handle NAs in correlation matrix (occurs when pairwise sample overlap is insufficient)
          na_cor_count <- sum(is.na(raw_cor))
          if (na_cor_count > 0) {
            diag_status <- paste0("Warning_NAs_in_Correlation_Matrix(", na_cor_count, ")")
            raw_cor[is.na(raw_cor)] <- 0 # Zero-imputation to preserve matrix structure for graph metrics
          }
          
          abs_cor <- abs(raw_cor)
          diag(abs_cor) <- 0 
          
          c_curr <- intersect(client_nodes, valid_nodes)
          t_curr <- intersect(therapist_nodes, valid_nodes)
          
          networks <- list(full = abs_cor)
          if (length(c_curr) > 1) networks$client_internal <- abs_cor[c_curr, c_curr, drop = FALSE]
          if (length(t_curr) > 1) networks$therapist_internal <- abs_cor[t_curr, t_curr, drop = FALSE]
          if (length(c_curr) > 0 && length(t_curr) > 0) networks$interpersonal <- abs_cor[c_curr, t_curr, drop = FALSE]
          
          # ------------------------------------------------------------------------------
          # Step 6: Compute Graph Topological Metrics
          # ------------------------------------------------------------------------------
          for (net_name in names(networks)) {
            adj_m <- networks[[net_name]]
            
            if (nrow(adj_m) >= 3 && nrow(adj_m) == ncol(adj_m)) {
              g <- graph_from_adjacency_matrix(adj_m, mode = "undirected", weighted = TRUE, diag = FALSE)
              
              num_nodes <- vcount(g)
              num_edges <- ecount(g)
              
              density_val    <- edge_density(g)
              efficiency_val <- global_efficiency(g, weights = E(g)$weight, directed = FALSE)
              
              modularity_val <- suppressWarnings({
                comm <- cluster_louvain(g, weights = E(g)$weight)
                modularity(comm)
              })
              if(is.nan(modularity_val)) modularity_val <- 0
              
              node_strength      <- strength(g, weights = E(g)$weight)
              node_ev_centrality <- tryCatch(eigen_centrality(g, weights = E(g)$weight)$vector, error = function(e) rep(0, vcount(g)))
              
              network_metrics_results <- bind_rows(network_metrics_results, tibble(
                session         = sess,
                window_id       = window_id,
                start_time      = start_t,
                end_time        = end_t,
                network_type    = net_name,
                num_nodes       = num_nodes,
                num_edges       = num_edges,
                density         = density_val,
                efficiency      = efficiency_val,
                modularity      = modularity_val,
                mean_strength   = mean(node_strength, na.rm = TRUE),
                mean_centrality = mean(node_ev_centrality, na.rm = TRUE)
              ))
            }
            
            # --- Compute Dynamic Similarity / Network Reconfiguration ---
            if (!is.null(prev_adj_matrices[[net_name]])) {
              prev_m <- prev_adj_matrices[[net_name]]
              if (all(dim(adj_m) == dim(prev_m)) && all(rownames(adj_m) == rownames(prev_m))) {
                if (nrow(adj_m) == ncol(adj_m)) {
                  vec_current <- adj_m[lower.tri(adj_m)]
                  vec_prev    <- prev_m[lower.tri(prev_m)]
                } else {
                  vec_current <- as.vector(adj_m)
                  vec_prev    <- as.vector(prev_m)
                }
                
                if (length(vec_current) > 1 && sd(vec_current, na.rm = TRUE) > 0 && sd(vec_prev, na.rm = TRUE) > 0) {
                  net_sim  <- cor(vec_current, vec_prev, use = "pairwise.complete.obs")
                  net_dist <- 1 - net_sim
                } else {
                  net_sim  <- NA
                  net_dist <- NA
                }
                
                dynamic_similarity_results <- bind_rows(dynamic_similarity_results, tibble(
                  session                 = sess,
                  window_id               = window_id,
                  network_type            = net_name,
                  network_correlation     = net_sim,
                  network_reconfiguration = net_dist
                ))
              }
            }
            prev_adj_matrices[[net_name]] <- adj_m
          }
          
          # ------------------------------------------------------------------------------
          # Compute Classical Dyadic Synchrony (Absolute Bivariate Correlations)
          # ------------------------------------------------------------------------------
          get_cor <- function(n1, n2) {
            if(all(c(n1, n2) %in% rownames(raw_cor))) {
              return(raw_cor[n1, n2])
            } else { return(NA) }
          }
          
          mov_s  <- abs(get_cor("BodyC", "BodyT"))
          head_s <- abs(get_cor("HeadC", "HeadT"))
          face_s <- abs(get_cor("C_Valence", "T_Valence"))
          eda_s  <- abs(get_cor("C_Mean SC", "T_Mean SC"))
          hr_s   <- abs(get_cor("C_HR", "T_HR"))
          rsa_s  <- abs(get_cor("C_RSA", "T_RSA"))
          
          mean_s <- mean(c(mov_s, head_s, face_s, eda_s, hr_s, rsa_s), na.rm = TRUE)
          
          classical_sync_results <- bind_rows(classical_sync_results, tibble(
            session       = sess, 
            window_id     = window_id,
            movement_sync = mov_s, 
            head_sync     = head_s,
            face_sync     = face_s, 
            eda_sync      = eda_s,
            hr_sync       = hr_s, 
            rsa_sync      = rsa_s, 
            mean_sync     = mean_s
          ))
          
        } else { diag_status <- "Dropped_Not_Enough_Valid_Nodes" }
      } else { diag_status <- "Dropped_Not_Enough_Available_Nodes" }
    }
    
    # Store quality diagnostics log for current window
    window_diagnostics_results <- bind_rows(window_diagnostics_results, tibble(
      session = sess,
      window_id = window_id,
      start_time = start_t,
      actual_samples = actual_samples,
      max_time_gap = max_time_gap,
      status = diag_status,
      missing_nodes = diag_missing_nodes,
      dropped_nodes = diag_dropped_flatlines
    ))
    
    window_id <- window_id + 1
    start_t   <- start_t + step
  }
}

print("Preview of window diagnostics log - details of window validity:")
print(head(window_diagnostics_results))

# ------------------------------------------------------------------------------
# Step 7: Export Analysis Results to CSV Files
# ------------------------------------------------------------------------------

write_csv(network_metrics_results, "network_metrics_results.csv")
write_csv(dynamic_similarity_results, "dynamic_similarity_results.csv")
write_csv(classical_sync_results, "classical_sync_results.csv")
write_csv(window_diagnostics_results, "window_diagnostics_results.csv")

message("Files saved successfully to the working directory!")


# ==============================================================================
# Step 8: Descriptive Statistics and Sample Summary
# ==============================================================================

# 1. Session- and Patient-Level Descriptives
sample_summary <- full_dataset %>%
  group_by(session) %>%
  summarise(
    patient = first(patient),
    duration_sec = max(second, na.rm = TRUE) - min(second, na.rm = TRUE) + 1,
    duration_min = duration_sec / 60,
    .groups = "drop"
  )

cat("=== Sample Characteristics ===\n")
cat("Total Unique Patients:", n_distinct(sample_summary$patient), "\n")
cat("Total Unique Sessions:", n_distinct(sample_summary$session), "\n")
cat("Mean Session Duration (min):", round(mean(sample_summary$duration_min), 2), 
    "(SD =", round(sd(sample_summary$duration_min), 2), 
    ", Range:", round(min(sample_summary$duration_min), 2), "-", round(max(sample_summary$duration_min), 2), ")\n")
cat("Total Recorded Time (hours):", round(sum(sample_summary$duration_sec) / 3600, 2), "\n\n")

# 2. Sliding Window and Data Quality Diagnostics Summary
if (exists("window_diagnostics_results") && nrow(window_diagnostics_results) > 0) {
  window_summary <- window_diagnostics_results %>%
    group_by(session) %>%
    summarise(
      total_windows = n(),
      valid_windows = sum(status == "Valid" | str_detect(status, "Warning")),
      dropped_windows = sum(str_detect(status, "Dropped")),
      .groups = "drop"
    )
  
  cat("=== Windowing & Data Quality Summary ===\n")
  cat("Total Sliding Windows:", sum(window_summary$total_windows), "\n")
  cat("Mean Windows per Session:", round(mean(window_summary$total_windows), 2), 
      "(SD =", round(sd(window_summary$total_windows), 2), ")\n")
  
  cat("\nBreakdown of Window Statuses:\n")
  print(table(window_diagnostics_results$status))
}

# Export descriptive summary table
write_csv(sample_summary, "sample_summary_descriptives.csv")