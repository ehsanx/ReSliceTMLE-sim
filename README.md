# ReSliceTMLE-sim

Code for replicating the simulation study comparing different variants of Targeted Maximum Likelihood Estimation (TMLE) in the manuscript. Includes scripts for simulation, data analysis, and reproducing the visualization of results.

---
### 📂 TMLE Variants
The simulation study evaluates the following TMLE variants:

| **Variant Name**       | **Description**                                           |
|-------------------------|------------------------------------------------------------|
| `vanilla`              | Standard (Vanilla) TMLE implementation                     |
| `cvq`                  | Cross-validated Q-learning TMLE                            |
| `cvq_multiple`         | CV-QTMLE with multiple repetitions                         |
| `fullcv`               | Full cross-validated TMLE                                  |
| `fullcv_multiple`      | Full CV-TMLE with multiple repetitions                     |
| `singlecrossfit`       | Single crossfitting TMLE                                   |
| `doublecrossfit`       | Double crossfitting TMLE                                   |

---
## ⚙️ How to Run Simulations

### 1. Setup
- **R required** (version ≥ 4.1 recommended)
- Install dependencies:
  ```r
  # Install devtools if not already installed
  install.packages("devtools")
  
  # Install ReSliceTMLE package from GitHub
  devtools::install_github("ehsanx/ReSliceTMLE")
  
  # Install other required packages
  install.packages(c(
    "rsimsum",
    "ggplot2",
    "gridExtra",
    "cowplot",
    "dplyr",
    "tidyr",
    "patchwork",
    "stringr"
  ))
  ```

### 2. Run the Simulation and Analysis
**Step 1: Execute the script**
```r
source("summarize_simulation_results.R")
```
This script performs the entire workflow:
- Loads the simulated data
- Runs different TMLE variants on the data
- Processes the results
- Generates visualizations for the manuscript figures

**Step 2: Adjust simulation parameters (optional)**
- To control the simulation, modify these parameters in the script:
  - `num_repetitions`: Number of repetitions for methods with multiple runs (default: 100)
  - `Q.SL.library` and `g.SL.library`: SuperLearner libraries for outcome and treatment models

---
## 📊 Visualization and Analysis

The script produces two main figures:

1. **Figure**: Comprehensive comparison of TMLE variants across multiple performance metrics
   - Includes bias, coverage, MSE, and other key metrics
   - Presents results with confidence intervals

2. **Figure**: Performance comparison of TMLE variants across increasing numbers of replicates
   - Shows how performance metrics converge with more replicates
   - Compares multiple TMLE variants against Vanilla TMLE reference

Both figures are displayed in the R graphics device and can be saved using standard R plotting commands.

---
## 📝 Notes
- The full simulation with 100 repetitions may take several hours to complete depending on your hardware
- For a quicker test run, reduce `num_repetitions` to a smaller value (e.g., 10)
- Results are saved to the `results/` directory by setting `save_results = TRUE` and `results_dir = "./results"` in the function call, allowing for analysis without rerunning simulations
- Simulation results are summarized using the `rsimsum` package, which calculates key performance metrics (bias, coverage, MSE, etc.) for each TMLE variant. The script includes functions to process both the final results (`results$raw_final`) and intermediate results with multiple iterations (`results$raw_intermediate`), allowing for comprehensive analysis of method performance across different numbers of replicates.
- The code is designed to work with the ReSliceTMLE package, which provides implementations of all TMLE variants
- The script handles both the simulation run and the visualization in a single workflow

---
## 📄 License
TBD

