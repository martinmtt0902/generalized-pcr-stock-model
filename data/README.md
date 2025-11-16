# Data (Not Included)

This project uses proprietary financial datasets exported from the **Bloomberg Terminal**, including:

- Daily stock prices (`Z_PX_LAST_stock1.csv`)
- 9 financial ratios (e.g., Alpha, Bid–Ask Spread, Market Cap, P/B, P/S, Beta, Volume)

Because of Bloomberg licensing restrictions, these CSV files **cannot be shared** in a public repository and are intentionally excluded via `.gitignore`.

## To reproduce the analysis

1. Export ratio and price data from Bloomberg with the same shape:
   - Rows = dates  
   - First column = date  
   - Remaining columns = individual stocks  

2. Place all CSV files in the `data/` folder.

3. Ensure the **price CSV** is alphabetically last in the folder.

4. Then run:

```r
source("src/gpcr_stock_pricing_backtest.R")
source("src/gpcr_stock_overpricing_tree.R")
