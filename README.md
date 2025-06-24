# MultivariateEquivalency

**Multivariate Equivalency and Power Analysis Using Cramér’s V and Truncated Geometric Distributions**

---

## Overview

The `MultivariateEquivalency` package provides tools for determining **statistical equivalency between multivariate distributions**, particularly when features are binned into ordinal categories. It is designed to support:

- **Equivalency testing** using multivariate generalizations of Cramér’s V  
- **Power analysis** for determining optimal experimental design parameters (e.g., sample size, number of variables)  
- **Cost modeling** for estimating the monetary burden of proposed experiments  
- **Process stability test** of multivariate processes using Hotelling’s \( T^2 \)

This package is especially suited for high-dimensional data common in process monitoring and experimental validation contexts (e.g., additive manufacturing, quality control). For further information on the scope and underlying workings of this package, see Lynch et al. (paper). 

---

## Input Data Format

Both reference and candidate datasets should be data frames or CSV files with the following structure:

- Each row is an observation (e.g., a part, trial, or individual)  
- Each column is a numeric variable to be tested  
- No missing values (use `na.omit()` or imputation before analysis)  
- Variables are assumed to be continuous and will be binned into categories internally  

**Example structure:**

| Var1 | Var2 | Var3 |
|------|------|------|
| 1.23 | 0.45 | 0.91 |
| 1.17 | 0.49 | 1.10 |
| 1.28 | 0.51 | 0.85 |

---

## Installation

To install the latest development version from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install from GitHub
devtools::install_github("colinmichaellynch/MultivariateEquivalency", build_vignettes = TRUE)
```
## Running the Tutorial

This package includes a full tutorial vignette that walks through:

- Checking process stability  
- Computing bin edges  
- Estimating power or required sample size  
- Optimizing for cost  
- Performing multivariate equivalency tests

To open the tutorial:

```r
browseVignettes("MultivariateEquivalency")
# or
vignette("tutorial", package = "MultivariateEquivalency")
```

---

##  Example Usage

Example data used in the tutorial is located in:

```r
system.file("extdata", "referenceProcess.csv", package = "MultivariateEquivalency")
system.file("extdata", "candidateProcessEquivalent.csv", package = "MultivariateEquivalency")
system.file("extdata", "candidateProcessNotEquivalent.csv", package = "MultivariateEquivalency")

```

```r
# Load reference process data
ref_path <- system.file("extdata", "referenceProcess.csv", package = "MultivariateEquivalency")
reference <- read.csv(ref_path)

# Check for stability
multivariateProcessStability(reference)

# Compute bin edges
bins <- binLimits(reference, b = 8)

# Load candidate process
cand_path <- system.file("extdata", "candidateProcessEquivalent.csv", package = "MultivariateEquivalency")
candidate <- read.csv(cand_path)

# Run equivalency test
multivariateEquivalency(bins, candidate, b = 8)
```

---

## License
This package is released under the **MIT License**.
