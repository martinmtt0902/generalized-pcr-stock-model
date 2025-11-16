############################################################
## Generalized PCR Stock Pricing & Backtest Script
##
## Usage:
##   - Place all Bloomberg CSV files in the `data/` folder of
##     this repository, e.g.:
##         ALPHA_OVERRIDABLE_stock1.csv
##         AVERAGE_BID_ASK_SPREAD_stock1.csv
##         BETA_ADJ_OVERRIDABLE_stock1.csv
##         CUR_MKT_CAP_stock1.csv
##         EARN_YLD_HIST_stock1.csv
##         PX_CLOSE_1D_stock1.csv
##         PX_TO_BOOK_RATIO_stock1.csv
##         PX_TO_SALES_RATIO_stock1.csv
##         PX_VOLUME_stock1.csv
##         Z_PX_LAST_stock1.csv  (PRICE – must be last in alphabetical order)
##
##   - Run this script from the repository root with working
##     directory set to the repo root.
##
##   - This code builds a Generalized PCR (GPCR) model,
##     produces intrinsic prices, and runs a simple backtest.
############################################################

## ---------------------------------------------------------
## User settings
## ---------------------------------------------------------

# Relative path to the folder containing all CSV files
data_dir <- "data"

# Length of each moving window (number of time points)
l <- 7L

## ---------------------------------------------------------
## Load file list and basic dimensions
## ---------------------------------------------------------

files <- list.files(
  path       = data_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

if (length(files) < 2L) {
  stop("Expecting at least one ratio CSV and one price CSV in `data/`.")
}

# Sort file names to ensure deterministic order (alphabetical)
files <- sort(files)

# Assume all ratio CSVs share the same dimension and the last CSV is the price file
example_ratio <- read.csv(files[1], check.names = FALSE)
dimen         <- dim(example_ratio)

q <- dimen[2] - 1L       # number of stocks (all columns except date)
p <- length(files) - 1L  # number of ratios (all files except the last price file)
n <- dimen[1]            # number of time points

# Number of moving windows (periods)
per <- floor(n / l)

# Containers:
YN      <- array(NA_real_, c(per, q))  # observed prices y_n
YN_HAT  <- array(NA_real_, c(per, q))  # intrinsic prices y_n_hat
deci    <- array(NA_real_, c(per, q))  # 1 = underpriced (long), 0 = overpriced (short)

cat("q (stocks):", q, "| p (ratios):", p, "| n (time points):", n, "\n")
cat("Window length l:", l, "| Number of periods:", per, "\n\n")

## ---------------------------------------------------------
## Main loop over moving windows
## ---------------------------------------------------------

for (w in seq_len(per)) {
  
  cat("Window", w, "of", per, "...\n")
  
  ###########################################################
  ### Step (0): Build 3D array x[i, j, t] for this window ###
  ###########################################################
  # x: [stock i, ratio j, time within window t]
  x <- array(NA_real_, dim = c(q, p, l))
  
  for (k in seq_len(p)) {
    
    # Read each ratio sheet
    ratio_raw <- read.csv(files[k], check.names = FALSE)
    
    # Subset to this moving window and drop date column (assumed first)
    idx_rows   <- (w * l - l + 1):(w * l)
    ratio_block <- ratio_raw[idx_rows, -1L, drop = FALSE]
    
    # Fill x: (stock i, ratio k, time t)
    for (t in seq_len(l)) {
      for (i in seq_len(q)) {
        x[i, k, t] <- ratio_block[t, i]
      }
    }
  }
  
  ###########################################################
  ### Step (1): PCA and compressed matrices Z_t           ###
  ###########################################################
  
  # Z: [stock i, component (1..3), time t]
  Z <- array(NA_real_, dim = c(q, 3L, l))
  
  for (t in seq_len(l)) {
    
    pca <- princomp(x[, , t], cor = TRUE)
    
    pc1 <- pca$loadings[, 1]
    pc2 <- pca$loadings[, 2]
    pc3 <- pca$loadings[, 3]
    
    Z1 <- Z2 <- Z3 <- rep(0, q)
    
    for (i in seq_len(p)) {
      Z1 <- Z1 + pc1[i] * x[, i, t]
      Z2 <- Z2 + pc2[i] * x[, i, t]
      Z3 <- Z3 + pc3[i] * x[, i, t]
    }
    
    Z[, , t] <- cbind(Z1, Z2, Z3)  # compressed matrix Z_t
  }
  
  ##############################################################
  ### Step (2): Construction of multivariate basis functions ###
  ##############################################################
  
  N <- 7L    # max exponent / basis index
  K <- 100L  # number of basis functions used per PC
  
  # g_k: one multivariate basis function at index k, for matrix Z_t (q x 3)
  g_k <- function(k, Zt) {
    set.seed(k %% K)  # K possible Gamma choices; k%%K reused across blocks
    
    # Random exponents for each stock
    Gamma <- sample.int(N, size = nrow(Zt), replace = TRUE)
    
    # Decide which principal component (column of Zt) to use
    # out = (1,2,3,1,2,3) for 6K basis functions
    out  <- c(1L, 2L, 3L, 1L, 2L, 3L)
    ind  <- ceiling(k / K)
    pc_idx <- out[ind]
    
    # Polynomial-style basis: power of chosen PC
    Zt[, pc_idx]^Gamma
  }
  
  ##########################################################
  ### Step (3): Data matrix D_tilde & covariance matrix  ###
  ##########################################################
  
  # g(t, Z): returns a 6K x q matrix (transposed inside)
  g <- function(t, Z) {
    out <- array(NA_real_, dim = c(q, 6L * K))
    for (k in seq_len(6L * K)) {
      out[, k] <- g_k(k, Z[, , t])
    }
    t(out)  # final dimension: 6K x q
  }
  
  # High-dimensional data matrix D_tilde: [6K, q * l]
  D_tilde <- array(NA_real_, dim = c(6L * K, q * l))
  
  for (t in seq_len(l)) {
    D_tilde[, ((t - 1) * q + 1):(t * q)] <- g(t, Z)
  }
  
  # Sample covariance of D_tilde' (covariance of rows of g_k)
  S <- cov(t(D_tilde))  # dimension: 6K x 6K
  
  #########################################################
  ### Step (4): Spectral decomposition of S             ###
  #########################################################
  
  eig  <- eigen(S)
  H    <- eig$vectors  # eigenvectors (6K x 6K)
  # D <- diag(eig$values)  # eigenvalues (not used below)
  
  #######################################################
  ### Step (5): Principal components P_j(Z_t)          ###
  #######################################################
  
  # P: [stock i, PC index j, time t]
  P <- array(0, dim = c(q, 6L * K, l))
  
  for (t in seq_len(l)) {
    for (k in seq_len(6L * K)) {
      temp <- g_k(k, Z[, , t])  # length q
      for (j in seq_len(6L * K)) {
        P[, j, t] <- P[, j, t] + H[k, j] * temp
      }
    }
  }
  
  #########################################################
  ### Step (6): Regression coefficients (OLS)            ###
  #########################################################
  
  M <- K  # number of PCs to keep (M << 6K)
  
  # Read price file (last CSV in alphabetical order)
  price_raw <- read.csv(files[length(files)], check.names = FALSE)
  idx_rows  <- (w * l - l + 1):(w * l)
  price_block <- as.matrix(price_raw[idx_rows, -1L, drop = FALSE])
  
  # Build Y (vectorized prices) and PZ (design matrix)
  Y  <- rep(NA_real_, l * q)
  PZ <- array(NA_real_, dim = c(l * q, 6L * K))
  
  for (j in seq_len(l)) {
    row_idx <- ((j - 1) * q + 1):(j * q)
    Y[row_idx]    <- price_block[j, ]
    PZ[row_idx, ] <- P[, , j]
  }
  
  Y  <- as.numeric(Y)
  PZ <- PZ[, 1:M, drop = FALSE]
  
  XtX <- crossprod(PZ)
  XtY <- crossprod(PZ, Y)
  
  betahat <- solve(XtX, XtY)
  
  #########################################################
  ### Step (7): Estimation of y_n (intrinsic price)     ###
  #########################################################
  
  # PCs at last time in window
  Pn <- P[, 1:M, l, drop = FALSE]
  
  yn_hat <- as.numeric(Pn %*% betahat)    # intrinsic prices (length q)
  yn     <- as.numeric(price_block[l, ])  # actual prices at last time
  
  YN[w, ]     <- yn
  YN_HAT[w, ] <- yn_hat
  deci[w, ]   <- as.numeric(yn < yn_hat)  # 1 = underpriced (long), 0 = overpriced (short)
  
  cat("  Completed window", w, "\n\n")
}

#########################################################
#### Step (8): Trading strategy & backtest            ####
#########################################################

# Passive return: buy deci[1,]==1 portfolio at start, hold until last period
V0          <- sum(deci[1, ] * YN[1, ])
VT          <- sum(deci[1, ] * YN[per, ])
Pass.return <- (VT - V0) / V0

# Market return: buy all stocks at start, hold until last period
All.0      <- sum(YN[1, ])
All.T      <- sum(YN[per, ])
All.return <- (All.T - All.0) / All.0

# Active strategy P&L over time
PnL <- rep(NA_real_, per)
PnL[1] <- 0

for (i in 2:per) {
  PnL[i] <- PnL[i - 1] -
    sum(deci[i, ]     * YN[i, ]) +
    sum(deci[i - 1, ] * YN[i - 1, ])
}

Active.return <- PnL[per] / V0

# Returns: market, passive, active
c(All.return, Pass.return, Active.return)

# Back-test accuracy: how often direction matches decision
Diff <- YN[-1, ] - YN[-per, ]  # price difference between periods
acc  <- sum((Diff > 0) == deci[-per, ]) / ((per - 1) * q)
acc  # accuracy rate

#########################################################
#### Plots: intrinsic vs market price (last window)   ####
#########################################################

yn_last     <- yn
yn_hat_last <- yn_hat

par(mfrow = c(1, 3))

plot(
  yn_last, yn_hat_last,
  main = "Intrinsic Price vs Market Price",
  xlab = "Market Price",
  ylab = "Intrinsic Price"
)
abline(0, 1, col = 2, lwd = 2)

plot(
  yn_last, yn_hat_last,
  main = "Intrinsic vs Market (< $20)",
  xlab = "Market Price",
  ylab = "Intrinsic Price",
  xlim = c(0, 20),
  ylim = c(0, 20)
)
abline(0, 1, col = 2, lwd = 2)

plot(
  yn_last, yn_hat_last,
  main = "Intrinsic vs Market (< $5)",
  xlab = "Market Price",
  ylab = "Intrinsic Price",
  xlim = c(0, 5),
  ylim = c(0, 5)
)
abline(0, 1, col = 2, lwd = 2)
