Multi-Omics Integration: Ketogenic Diet & Melanoma
This repository houses the computational framework and analysis pipelines used to investigate the metabolic and transcriptomic effects of a ketogenic diet (KD) on melanoma xenografts.

🔬 Research Scope
We utilize two primary integration paradigms:

Correlation-based analysis: Identifying co-regulatory modules between datasets.

Multi-omics integration: Utilizing mixOmics/DIABLO to identify discriminant features and predictive signatures.

Data generated at SALK/PMU Salzburg; availability pending EU regulatory compliance.

🛠 Technical Stack
Package Development: KetoAnalysisTools (custom R package)

Integration: mixOmics

Visualization: pathview, ggplot2, and custom diagnostic plots.


## 📂 Repository Structure

```text
.
├── KetoAnalysisTools/        # R package source code
│   ├── DESCRIPTION           # Package metadata
│   ├── R/                    # Core integration functions
│   ├── man/                  # Documentation
│   └── vignettes/            # Analysis notebooks
├── analysis/                 # Research scripts & workflows
│   ├── 01_preprocessing.Rmd
│   ├── 02_integration.Rmd
│   └── outputs/              # Generated plots/tables
├── data/                     # Raw input data
├── .gitignore                # Project-wide ignore rules
└── README.md                 # This documentation
