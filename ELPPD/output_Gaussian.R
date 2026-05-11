rm(list = ls())



# ---------------------------------------------------------
# User inputs
# ---------------------------------------------------------

n <- 50

p <- 10

sigma2 <- 1

tau2 <- 1e2

reps <- 100

K.vec <- c(2, 5, 10, 20)



# ---------------------------------------------------------
# Exact ELPPD
# ---------------------------------------------------------

exact_elpd_gaussian <- function(X,
                                y,
                                sigma2,
                                tau2){
  
  n <- nrow(X)
  
  p <- ncol(X)
  
  elpd <- 0
  
  for(i in 1:n){
    
    Xi <- matrix(
      X[i, ],
      1,
      p
    )
    
    yi <- y[i]
    
    X_minus <- X[-i, , drop = FALSE]
    
    y_minus <- y[-i]
    
    V_post <- solve(
      crossprod(X_minus) / sigma2 +
        diag(1 / tau2, p)
    )
    
    m_post <- V_post %*%
      crossprod(X_minus, y_minus) /
      sigma2
    
    pred_mean <- Xi %*% m_post
    
    pred_var <- sigma2 +
      Xi %*% V_post %*% t(Xi)
    
    elpd <- elpd +
      dnorm(
        yi,
        pred_mean,
        sqrt(pred_var),
        log = TRUE
      )
  }
  
  as.numeric(elpd)
}



# ---------------------------------------------------------
# Simulate reference dataset
# ---------------------------------------------------------

set.seed(123)

X <- matrix(
  rnorm(n * p),
  n,
  p
)

beta_true <- rnorm(p)

y <- X %*% beta_true +
  rnorm(
    n,
    sd = sqrt(sigma2)
  )

y <- as.vector(y)



# ---------------------------------------------------------
# Exact ELPPD
# ---------------------------------------------------------

true_elpd <- exact_elpd_gaussian(
  X,
  y,
  sigma2,
  tau2
)



# ---------------------------------------------------------
# Storage
# ---------------------------------------------------------

g.list <- vector(
  "list",
  length(K.vec)
)

b.list <- vector(
  "list",
  length(K.vec)
)

z.list <- vector(
  "list",
  length(K.vec)
)



# ---------------------------------------------------------
# Read results
# ---------------------------------------------------------

for(k in seq_along(K.vec)){
  
  K <- K.vec[k]
  
  load(
    paste0(
      K,
      "_fold_elppd_results.RData"
    )
  )
  
  samp.var.g <- var(results$elppd_g)
  
  samp.var.b <- var(results$elppd_b)
  
  samp.var.z <- var(results$elppd_z)
  
  g.list[[k]] <- c(
    results$var_g,
    samp.var.g
  )
  
  b.list[[k]] <- c(
    results$var_b,
    samp.var.b
  )
  
  z.list[[k]] <- c(
    results$var_z,
    samp.var.z
  )
}



# ---------------------------------------------------------
# PDF output
# ---------------------------------------------------------

pdf(
  "elppd_variance_boxplots.pdf",
  width = 12,
  height = 8
)



# ---------------------------------------------------------
# Plot layout
# ---------------------------------------------------------

par(
  mfrow = c(1, 3),
  mar = c(4, 4, 2, 1),
  oma = c(0, 0, 1, 0)
)



# ---------------------------------------------------------
# Gelfand
# ---------------------------------------------------------

boxplot(
  g.list,
  names = K.vec,
  main = "Gelfand",
  xlab = "Number of folds",
  ylab = "Variance",
  outline = FALSE,
  frame.plot = TRUE
)

points(
  1:length(K.vec),
  sapply(g.list, tail, 1),
  pch = 19,
  cex = 1
)

legend(
  "topright",
  legend = c(
    "Asymptotic variance",
    "Sample variance"
  ),
  pch = c(NA, 19),
  lty = c(1, NA),
  bty = "n",
  cex = 0.8
)



# ---------------------------------------------------------
# Bronze
# ---------------------------------------------------------

boxplot(
  b.list,
  names = K.vec,
  main = "Bronze",
  xlab = "Number of folds",
  ylab = "Variance",
  outline = FALSE,
  frame.plot = TRUE
)

points(
  1:length(K.vec),
  sapply(b.list, tail, 1),
  pch = 19,
  cex = 1
)



# ---------------------------------------------------------
# Zanella
# ---------------------------------------------------------

boxplot(
  z.list,
  names = K.vec,
  main = "Zanella",
  xlab = "Number of folds",
  ylab = "Variance",
  outline = FALSE,
  frame.plot = TRUE
)

points(
  1:length(K.vec),
  sapply(z.list, tail, 1),
  pch = 19,
  cex = 1
)



# ---------------------------------------------------------
# Close device
# ---------------------------------------------------------

dev.off()