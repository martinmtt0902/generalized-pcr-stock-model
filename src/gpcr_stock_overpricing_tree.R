############################################################
## GPCR Stock Overpricing Classification Tree Script
##
## Usage:
##   - Place your ratio CSV files and one price CSV file in
##     the `data/` folder of this repository.
##   - The price CSV (y_t) must be the LAST file in
##     alphabetical order.
##
##   - This script:
##       * uses the last l = 30 time points,
##       * computes GPCR intrinsic prices,
##       * defines an "overpriced" indicator,
##       * fits logistic regression and a classification tree.
############################################################

## ---------------------------------------------------------
## User settings
## ---------------------------------------------------------

data_dir <- "data"   # relative path for GitHub compatibility
l        <- 30L      # length of moving window (last l time points)

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

# Sort for deterministic order: last file will be treated as price file y_t
files <- sort(files)

# Assume all ratio CSVs share the same dimensions; use the first ratio file
example_ratio <- read.csv(files[1], check.names = FALSE)
dimen         <- dim(example_ratio)

q <- dimen[2] - 1L       # number of stocks (all columns except date)
p <- length(files) - 1L  # number of ratios (all files except the last price file)
n <- dimen[1]            # number of time points

cat("q (stocks):", q, "| p (ratios):", p, "| n (time points):", n, "\n")
cat("Using last", l, "time points.\n\n")

## ---------------------------------------------------------
## Rearrange data into x[i, j, t] for last l time points
## ---------------------------------------------------------

# x: [stock i, ratio j, time t] for last l dates
x <- array(NA_real_, dim = c(q, p, l))

for (k in seq_len(p)) {
  
  # Read each ratio sheet
  ratio_raw <- read.csv(files[k], check.names = FALSE)
  
  # Keep only the last l rows, drop date column (assumed first)
  idx_rows    <- (n - l + 1L):n
  ratio_block <- ratio_raw[idx_rows, -1L, drop = FALSE]
  
  for (t in seq_len(l)) {
    for (i in seq_len(q)) {
      x[i, k, t] <- ratio_block[t, i]
    }
  }
}

###########################################################
### Step (1): Apply PCA to obtain the first 3 loadings  ###
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
    # Evaluate 3 principal components for time t
    Z1 <- Z1 + pc1[i] * x[, i, t]
    Z2 <- Z2 + pc2[i] * x[, i, t]
    Z3 <- Z3 + pc3[i] * x[, i, t]
  }
  
  Z[, , t] <- cbind(Z1, Z2, Z3)  # compressed matrix Z_t
}

##############################################################
### Step (2): Construction of Multivariate Basis Function  ###
##############################################################

N <- 7L    # univariate basis index (Max exponent)
K <- 100L  # number of basis functions per PC (K << P^N_q)

# Define g_1, g_2, ..., g_{6K}
g_k <- function(k, Zt) {
  set.seed(k %% K)  # total of K possible Gamma choices
  
  # Random exponents for each stock
  Gamma <- sample.int(N, size = nrow(Zt), replace = TRUE)
  
  # Decide which principal component (column of Zt) to use
  # out = (1, 2, 3, 1, 2, 3) for 6K basis functions
  out    <- c(1L, 2L, 3L, 1L, 2L, 3L)
  ind    <- ceiling(k / K)
  pc_idx <- out[ind]
  
  # Polynomial-style basis: Z_t^{Gamma} on selected PC
  Zt[, pc_idx]^Gamma
}

##########################################################
### Step (3): Data Matrix D & Sample Covariance Matrix ###
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
  D_tilde[, ((t - 1L) * q + 1L):(t * q)] <- g(t, Z)
  if (t %% max(1L, round(l / 8L)) == 0L) {
    cat("  Building D_tilde, t =", t, "of", l, "\n")
  }
}

S <- cov(t(D_tilde))  # sample covariance of D_tilde

#########################################################
###       Step (4): Spectral Decomposition of D       ###
#########################################################

info <- eigen(S)
H    <- info$vectors  # eigenvectors
# D <- diag(info$values)  # eigenvalues (not used below)

#######################################################
###       Step (5): Principal Components of Z       ###
#######################################################

# P: [stock i, PC index j, time t]
P <- array(0, dim = c(q, 6L * K, l))

for (t in seq_len(l)) {
  for (k in seq_len(6L * K)) {
    temp <- g_k(k, Z[, , t])
    for (j in seq_len(6L * K)) {
      P[, j, t] <- P[, j, t] + H[k, j] * temp
    }
  }
  if (t %% 2L == 0L) {
    cat("  Computing PCs, t =", t, "of", l, "\n")
  }
}

#########################################################
###       Step (6): Regression coefficients           ###
#########################################################

M <- K  # choose M << 6K

# Price file is the last CSV in alphabetical order
y_raw <- read.csv(files[length(files)], check.names = FALSE)
y_mat <- as.matrix(y_raw[(n - l + 1L):n, -1L, drop = FALSE])

Y  <- rep(NA_real_, l * q)
PZ <- array(NA_real_, dim = c(l * q, 6L * K))

for (j in seq_len(l)) {
  idx <- ((j - 1L) * q + 1L):(j * q)
  Y[idx]    <- y_mat[j, ]
  PZ[idx, ] <- P[, , j]
}

Y  <- as.numeric(Y)
PZ <- PZ[, 1:M, drop = FALSE]

XtX <- crossprod(PZ)
XtY <- crossprod(PZ, Y)

betahat <- solve(XtX, XtY)

#########################################################
###            Step (7): Estimation of y_n            ###
#########################################################

# PCs at last time n (within window)
Pn <- P[, 1:M, l, drop = FALSE]

yn_hat <- as.numeric(Pn %*% betahat)    # intrinsic prices at time n
yn     <- as.numeric(y_mat[l, ])        # actual prices at time n

#########################################################
####         Step (8): Classification Tree           ####
#########################################################

# Overpricing indicator: 1 if overpriced, 0 otherwise
OP <- as.numeric(yn > yn_hat)

# Xn: ratio data matrix at time n (last slice)
Xn <- x[, , l]

# Combine predictors and response into a data frame
data <- cbind(Xn, OP)
colnames(data) <- c(
  "Overridable_alpha",
  "Average_Bid_Ask_Spread",
  "Overridable_Adjusted_Beta",
  "Log_Market_Value",
  "Earning_Yield",
  "Closing_Price_1D_Before",
  "Price_to_Book_Ratio",
  "Price_to_Sales_Ratio",
  "Log_Trading_Volume",
  "Overpriced"
)

data <- as.data.frame(data)
data$Overpriced <- factor(data$Overpriced, levels = c(0, 1))

## Logistic regression (for reference)
logit_fit <- glm(
  Overpriced ~ Overridable_alpha +
    Average_Bid_Ask_Spread +
    Overridable_Adjusted_Beta +
    Log_Market_Value +
    Earning_Yield +
    Closing_Price_1D_Before +
    Log_Trading_Volume +
    Price_to_Book_Ratio +
    Price_to_Sales_Ratio,
  family = binomial,
  data   = data
)

summary(logit_fit)

## Classification tree
library(rpart)

ctree <- rpart(
  Overpriced ~ Overridable_alpha +
    Average_Bid_Ask_Spread +
    Overridable_Adjusted_Beta +
    Log_Market_Value +
    Earning_Yield +
    Log_Trading_Volume +
    Price_to_Book_Ratio +
    Price_to_Sales_Ratio,
  data   = data,
  method = "class",
  control = rpart.control(maxdepth = 3)
)

plot(ctree, asp = 1)
text(ctree, cex = 0.8, use.n = TRUE)
print(ctree)
