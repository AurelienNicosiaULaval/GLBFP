## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, message=FALSE-------------------------------------------------
library(GLBFP)

# Load the example dataset `ashua` from the package
data("ashua")  

df = data.frame(level = ashua$level, log_flow = log(ashua$flow))
head(df)

## -----------------------------------------------------------------------------
compute_bi_optim(df, m = c(1,1))

## -----------------------------------------------------------------------------
b_ashua <- c(0.19, 0.15)  # Bin width for each dimension 
x_ashua <- c(round(mean(df$level),2), round(mean(df$log_flow),2)) # single point to evaluate

out <- LBFP(x = x_ashua, data = df, b = b_ashua)
out
names(out)


## -----------------------------------------------------------------------------
ASH(x_ashua, df, b_ashua)
GLBFP(x_ashua, df, b_ashua)

## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
GLBFP(x_ashua, df, b_ashua, m = c(5, 8))

## -----------------------------------------------------------------------------

out <- LBFP_estimate(df, b = b_ashua, grid_size = 20) # default is grid_size = 20
out
names(out)
head(cbind(out$grid, out$densities))

## -----------------------------------------------------------------------------

## grid zoomed in mean +- 1*sd

spec_grid <- list(level = seq(mean(df$level)-sd(df$level),mean(df$level)+sd(df$level),length.out = 10),
                   log_flow =                   seq(mean(df$log_flow)-sd(df$log_flow),mean(df$log_flow)+sd(df$log_flow),length.out = 10))
spec_grid

spec_grid_ex <- as.matrix(expand.grid(spec_grid))

out_spec_grid <- LBFP_estimate(df, b = b_ashua, grid_points = spec_grid_ex) # default is grid_size = 20
out_spec_grid
names(out_spec_grid)
head(cbind(out_spec_grid$grid, out_spec_grid$densities))


## -----------------------------------------------------------------------------
out_GLBFP <- GLBFP::GLBFP_estimate(df, b = b_ashua, m = c(3,3), grid_size = 5) # default is grid_size = 20
out_GLBFP
names(out_GLBFP)
head(cbind(out_GLBFP$grid, out_GLBFP$densities))

## -----------------------------------------------------------------------------

plot(out, contour = TRUE)


## -----------------------------------------------------------------------------
plot(out, contour = FALSE)

## ----echo=TRUE, message=FALSE-------------------------------------------------
library(MASS)  # For generating multivariate normal data

# Define the parameters for the Markov chain
P <- matrix(c(0.7, 0.2, 0.05, 0.05,
              0.15, 0.7, 0.1, 0.05,
              0, 0.15, 0.8, 0.05,
              0.15, 0.05, 0.2, 0.6),
            nrow = 4, byrow = TRUE)

# Stationary distribution
pi_tilde <- c(2/9, 1/3, 1/3, 1/9)

# Means of the bivariate normals
mu_list <- list(c(-1.5, -1.5),
                c(-1.5, 1.5),
                c(1.5, -1.5),
                c(1.5, 1.5))

# Covariance matrices
V1 <- V2 <- matrix(c(0.15, 0.0415, 0.0415, 0.25), nrow = 2)
V3 <- V4 <- matrix(c(0.15, -0.0415, -0.0415, 0.25), nrow = 2)
V_list <- list(V1, V2, V3, V4)

# Number of realizations to generate
n <- 1000

# Initialize vectors to store the results
Y <- numeric(n)
X <- matrix(0, nrow = n, ncol = 2)

# Generate Y0 from the stationary distribution
set.seed(123)  # For reproducibility
Y[1] <- sample(1:4, size = 1, prob = pi_tilde)

# Generate X0 based on Y0
X[1, ] <- mvrnorm(1, mu = mu_list[[Y[1]]], Sigma = V_list[[Y[1]]])

# Generate the Markov chain and bivariate normals
for (t in 2:n) {
  # Generate Y_t based on Y_{t-1}
  Y[t] <- sample(1:4, size = 1, prob = P[Y[t - 1], ])
  
  # Generate X_t based on Y_t
  X[t, ] <- mvrnorm(1, mu = mu_list[[Y[t]]], Sigma = V_list[[Y[t]]])
}

# Estimate the density across a grid using estimate_LBFP


b <- c(0.075, 0.15)  # Bin width for each dimension
result <- LBFP_estimate(X, b)
result2 <- GLBFP_estimate(X, b , m=c(5,5))
plot(result)
plot(result2)



