# Generalized PCR Model for Stock Price Modelling

This repository contains a financial modelling project that applies a Generalized Principal Component Regression (GPCR) framework to forecast intrinsic stock prices and study overpricing/underpricing patterns in Hong Kong stock data.

The group project was completed as part of the course “RMSC4002 – Financial Data Analytics with Machine Learning” at The Chinese University of Hong Kong.

## 📅 Project Information
- Completion date: December 2020
- Languages: R

------------------------------------------------------------
PROJECT OVERVIEW
------------------------------------------------------------

We estimate intrinsic stock prices using:

- Daily financial ratios for ~180 HKEX stocks
- PCA compression
- Multivariate basis expansions
- A second-stage high-dimensional PCA
- Regression on GPCR components

Two main analyses are included:

------------------------------------------------------------
1. Intrinsic Price Backtest
------------------------------------------------------------

Script:
    src/gpcr_stock_pricing_backtest.R

This script performs:
- Rolling window GPCR estimation
- Intrinsic price generation
- Underpriced / overpriced signal detection
- Simple long-only backtest
- Calculation of passive, market, and active returns
- Diagnostic plots

------------------------------------------------------------
2. Overpricing Classification Tree
------------------------------------------------------------

Script:
    src/gpcr_stock_overpricing_tree.R

This script:
- Focuses on the latest time window
- Produces an overpricing indicator (1 = overpriced)
- Fits:
    * Logistic regression
    * Classification tree (rpart)
- Identifies which financial ratios drive overpricing signals

------------------------------------------------------------
REPOSITORY STRUCTURE
------------------------------------------------------------

generalized-pcr-stock-model/
├── src/
│   ├── gpcr_stock_pricing_backtest.R
│   └── gpcr_stock_overpricing_tree.R
├── report/
│   └── RMSC4002_Project.pdf
├── data/
│   └── README.md   (Bloomberg data instructions only)
├── .gitignore
└── README.md

------------------------------------------------------------
DATA USAGE (IMPORTANT)
------------------------------------------------------------

The data used in this project come from the Bloomberg Terminal and CANNOT be redistributed publicly.

Therefore:
- No CSV files are included in this repository.
- All data files are intentionally excluded via `.gitignore`.

To run the code yourself:
- Export stock ratio and price files from Bloomberg with the same data layout.
- Place them in the “data/” folder.
- Ensure the price CSV is alphabetically last.
- Follow the instructions in data/README.md.

------------------------------------------------------------
AUTHOR & CONTRIBUTION
------------------------------------------------------------

This was a group project. My contributions included:
- Implementing major parts of the GPCR modelling pipeline
- Writing the intrinsic-price backtesting script
- Developing the overpricing classification analysis
- Plotting, visualisation, and accuracy analysis
- Writing major sections of the final report

------------------------------------------------------------
FULL REPORT
------------------------------------------------------------

The complete methodology and analysis are documented in:

    report/RMSC4002_Project.pdf

------------------------------------------------------------
RUNNING THE CODE
------------------------------------------------------------

After placing your Bloomberg CSV files in the “data/” folder, run:

    source("src/gpcr_stock_pricing_backtest.R")
    source("src/gpcr_stock_overpricing_tree.R")

------------------------------------------------------------
PROJECT HIGHLIGHTS
------------------------------------------------------------

- Designed a two-stage PCA model for cross-sectional compression
- Estimated intrinsic stock prices via GPCR
- Built an interpretable long-only backtest
- Studied overpricing using logistic regression and classification trees
- Demonstrated financial feature importance diagnostics

------------------------------------------------------------
Feel free to explore the scripts or contact me for more details.
------------------------------------------------------------
