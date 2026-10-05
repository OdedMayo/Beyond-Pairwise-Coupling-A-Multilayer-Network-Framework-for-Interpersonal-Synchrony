# Multimodal Dynamic Network Analysis of Dyadic Synchrony

This repository contains the full R analysis pipeline for processing, constructing, and analyzing dynamic psychophysiological and behavioral networks in patient-therapist dyads. 

The pipeline integrates multimodal time-series data (Body Movement, Head Motion, Facial Expression, Electrodermal Activity [EDA], Heart Rate [HR], and Respiratory Sinus Arrhythmia [RSA]) into a 12-node network model to examine dyadic interaction dynamics beyond bivariate synchrony.

---

## 📌 Repository Overview & Execution Workflow

To ensure reproducible data processing and analytical continuity, scripts should be executed in the following sequential order:

```text
┌─────────────────────────────────────────────────────────────────────────┐
│ 1. multi_modal_network.r                                                │
│    └─► Data validation, 1Hz downsampling, zero-alignment & merging      │
│    └─► Sliding-window network extraction & baseline graph metrics       │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
           ┌─────────────────────────┼─────────────────────────┐
           ▼                         ▼                         ▼
┌──────────────────────┐  ┌──────────────────────┐  ┌──────────────────────┐
│ 2. Modularity.R      │  │ 3. Participation     │  │ 4. Regime Shift.R    │
│                      │  │    Coefficient.R     │  │                      │
│ └─► Macro-level LMM  │  │ └─► Node-level       │  │ └─► Dynamic network  │
│     (3-Level LMM:    │  │     integration      │  │     reconfiguration  │
│     Network          │  │     (Client RSA      │  │     & structural     │
│     Modularity)      │  │     hubness)         │  │     rewiring         │
└──────────────────────┘  └──────────────────────┘  └──────────────────────┘
                                     │
                                     ▼
                          ┌──────────────────────┐
                          │ 5. Twin Window_      │
                          │    Graph_Final.R      │
                          │ └─► High-res visual  │
                          │     comparison of    │
                          │     twin networks    │
                          └──────────────────────┘
```

---

## 🛠 Requirements & Dependencies

The pipeline requires **R (v4.0.0 or higher)** and the following R packages:

```r
install.packages(c(
  "tidyverse",   # Data manipulation & ggplot2 visualization
  "igraph",      # Graph theory & network construction
  "qgraph",      # Network visualization
  "lme4",        # Multilevel Linear Mixed-Effects Models (LMM)
  "MuMIn"        # R-squared calculation for mixed models (GLMM)
))

# Note: brainGraph (for global efficiency calculation) may require Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("brainGraph")
```

---

## 📁 Data Architecture & File Descriptions

### Input Files (Raw Anonymized Data)
* `MOV_anonymized.csv`: Body movement time-series.
* `FACE_anonymized.csv`: Facial valence time-series.
* `EDA_anonymized.csv`: Electrodermal activity (Mean SC) time-series.
* `HR_anonymized.csv`: Heart rate time-series.
* `RSA_anonymized.csv`: Respiratory sinus arrhythmia time-series.
* `twin_windows_multivariate.csv`: Identified twin window candidate pairs for visual comparison.

### Core R Scripts

#### 1. `multi_modal_network.r` — Processing & Network Construction
* **Purpose**: Performs quality validation across all 5 modalities, temporal downsampling (1 Hz), zero-alignment, and time-series merging into `full_dataset_merged.csv`.
* **Sliding Window Parameters**: Window size = 30s, Step size = 15s, Minimum sample threshold = 15s.
* **Network Metrics Calculated**: Graph Density, Global Efficiency, Louvain Modularity ($Q$), Mean Node Strength, and Eigenvector Centrality across Client-internal, Therapist-internal, Interpersonal, and Full 12-node systems.
* **Outputs Generated**: `network_metrics_results.csv`, `dynamic_similarity_results.csv`, `classical_sync_results.csv`, `window_diagnostics_results.csv`, `sample_summary_descriptives.csv`.

#### 2. `Modularity.R` — Macro-Level Structure vs. Dyadic Synchrony
* **Purpose**: Evaluates whether global network modularity ($Q$) provides incremental structural information beyond classical bivariate synchrony.
* **Statistical Modeling**: Fits a 3-Level Linear Mixed-Effects Model (`lme4::lmer`) nesting sliding windows within sessions and dyads, isolating fixed effects via Marginal $R^2$ (`MuMIn`).
* **Outputs Generated**:
  * `Modularity_vs_Individual_Channels.png`: Panel view of modularity against individual synchrony channels.
  * `Predicted_vs_Actual_Modularity_3Level.png`: Predicted vs. actual modularity plot controlling for dyad/session random intercepts.

#### 3. `Participation Coefficient.R` — Micro/Meso Node-Level Analysis
* **Purpose**: Computes node-level Participation Coefficient (PC) to assess how individual modalities (e.g., Client RSA) integrate across network modules vs. pairwise synchrony.
* **Outputs Generated**:
  * `node_metrics_results.csv`: Node-level PC values across all sliding windows.
  * `Node_PC_vs_Dyadic_Sync.png`: Direct RSA dyadic synchrony vs. RSA Participation Coefficient.
  * `Node_PC_vs_Individual_Channels.png`: Faceted comparison of PC across modalities.
  * `Node_PC_vs_Multivariate_Model.png`: Actual vs. predicted PC from a multivariate synchrony model.

#### 4. `Regime Shift.R` — Dynamic Network Reconfiguration
* **Purpose**: Analyzes topological rewiring (network reconfiguration distance: $1 - r_{\text{adjacency}}$) across consecutive temporal windows relative to window-to-window changes in dyadic synchrony ($\Delta |r|$).
* **Outputs Generated**:
  * `RegimeShift_vs_Individual_Deltas.png`: Network rewiring vs. individual $\Delta |r|$ channels.
  * `RegimeShift_vs_Multivariate_Delta.png`: Actual regime shifts vs. multivariate model predictions.

#### 5. `Twin Window_Graph_24.9.R` — Topographical Comparison
* **Purpose**: Extracts and visualizes "twin network windows"—pairs of windows exhibiting identical pairwise synchrony profiles but distinct macro-topological modularity ($Q$).
* **Outputs Generated**:
  * `Twin_Networks_Clean_Labels.png`: High-resolution dual network plot rendered via `qgraph` highlighting topological divergence.

---

## 📐 Network Specification

The full network model comprises **12 nodes** representing 6 physiological and behavioral channels measured simultaneously from Client and Therapist:

| Node Label (Internal) | Display Label | Modality / Metric | Agent |
| :--- | :--- | :--- | :--- |
| `BodyC` | `C_Body` | Body Movement Dynamics | Client |
| `HeadC` | `C_Head` | Head Motion Dynamics | Client |
| `C_Valence` | `C_Face` | Facial Valence | Client |
| `C_Mean SC` | `C_EDA` | Skin Conductance Level | Client |
| `C_HR` | `C_HR` | Heart Rate | Client |
| `C_RSA` | `C_RSA` | Respiratory Sinus Arrhythmia | Client |
| `BodyT` | `T_Body` | Body Movement Dynamics | Therapist |
| `HeadT` | `T_Head` | Head Motion Dynamics | Therapist |
| `T_Valence` | `T_Face` | Facial Valence | Therapist |
| `T_Mean SC` | `T_EDA` | Skin Conductance Level | Therapist |
| `T_HR` | `T_HR` | Heart Rate | Therapist |
| `T_RSA` | `T_RSA` | Respiratory Sinus Arrhythmia | Therapist |

---

## 🚀 Quick Start Instructions

1. Clone or download this repository.
2. Ensure input CSV files are placed in the working directory (or update the `setwd()` path in the scripts).
3. Open R / RStudio and execute `multi_modal_network.r`.
4. Run subsequent analytical scripts (`Modularity.R`, `Participation Coefficient.R`, `Regime Shift.R`, `Twin Window_Graph_24.9.R`) as needed.