library(tidyverse)
library(gurobi)
# library(slam)

N <- 9 # Total number of nodes
d_vars <- expand.grid(1:N, 1:N)[-(1+(N+1)*(0:(N-1))),] %>% mutate(name=paste0("w_", Var2, Var1)) %>% pull(name)
d_vars <- c(d_vars, paste0("S", 1:3))
n_borders <- 3 # < 8, Node 5 and 8 can't be a border point (in-land)
PO <- 0.1 # overdose proportional parameter
PS <- 5 # seizure proportionality parameter



set.seed(100)
{
seizure <- numeric(9)
seizure[c(2, 5, 8)] <- runif(3, 1, 100) %>% sort(decreasing=T)
seizure[1] <- seizure[2]* runif(1, 0.9, 1.1)
seizure[3] <- seizure[2]* runif(1, 0.9, 1.1)
seizure[4] <- seizure[5]* runif(1, 0.9, 1.1)
seizure[6] <- seizure[5]* runif(1, 0.9, 1.1)
seizure[7] <- seizure[8]* runif(1, 0.9, 1.1)
seizure[9] <- seizure[8]* runif(1, 0.9, 1.1)
O <- c(runif(9, 100, 1000)) %>% round(0)

price <- rep(0 , 9)
price[2] <- runif(1, 1, 100)
price[5] <- price[2]* (1 + runif(1, 0.1, 0.7))
price[8] <- price[5]* (1 + runif(1, 0.1, 0.7))
price[1] <- price[2]* runif(1, 0.9, 1.1)
price[3] <- price[2]* runif(1, 0.9, 1.1)
price[4] <- max(c(price[5]* runif(1, 0.9, 1.1), 1.1*price[1]))
price[6] <- max(c(price[5]* runif(1, 0.9, 1.1), 1.1*price[3]))
price[7] <- max(c(price[8]* runif(1, 0.9, 1.1), 1.1*price[5]))
price[9] <- max(c(price[8]* runif(1, 0.9, 1.1), 1.1*price[6]))

purity <- rep(0, 9)
purity[2] <- runif(1, 0.5, 0.8)
purity[5] <- min(c(purity[2]* (1 + runif(1, 0, 0.2)), 0.99))
purity[8] <- purity[5]* (1 + runif(1, 0.01, 0.2))
purity[1] <- purity[2]* runif(1, 0.9, 1.1)
purity[3] <- purity[2]* runif(1, 0.9, 1.1)
purity[4] <- max(c(purity[5]* runif(1, 0.9, 1.1), 1.1*purity[1]))
purity[6] <- max(c(purity[5]* runif(1, 0.9, 1.1), 1.1*purity[3]))
purity[7] <- max(c(purity[8]* runif(1, 0.9, 1.1), 1.1*purity[5]))
purity[9] <- max(c(purity[8]* runif(1, 0.9, 1.1), 1.1*purity[6]))
purity <- ifelse(purity > 0.99, 0.99, purity)

matrix(price, 3, 3, byrow = T)
matrix(purity, 3, 3, byrow = T)


# Constraints
  # Conservation of flow/consumption
A1 <- matrix(0, N, N*(N-1)+n_borders)
for (i in 1:N) {
  A_ref <- matrix(0, N, N)
  A_ref[i,] <- -1
  A_ref[, i] <- 1
  S_vec <- rep(0, n_borders)
  if (i < 4) S_vec[i] <- 1
  A1[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
}
b1 <- rep(0, N)
b1_alt <- PO*O

  # State seizure
A2 <- matrix(0, N, N*(N-1)+n_borders)
for (i in 1:N) {
  A_ref <- matrix(0, N, N)
  A_ref[, i] <- 1
  S_vec <- rep(0, n_borders)
  if (i < 4) S_vec[i] <- 1
  A2[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
}
b2 <- rep(0, N)
b2_alt <- seizure

  # Directionality of flow
less_than_price <- expand.grid(price, price)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 < Var1) %>% pull(less_than) %>% which
less_than_purity <- expand.grid(purity, purity)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 < Var1) %>% pull(less_than) %>% which
zero_d_vars_index <- unique(c(less_than_price, less_than_purity)) %>% sort
n_zero_d_vars <- length(zero_d_vars_index)
A_zero_d_vars <- matrix(0, n_zero_d_vars, N*(N-1)+n_borders)
for (i in 1:n_zero_d_vars) {
  A_zero_d_vars[i, zero_d_vars_index[i]] <- 1
}
b_zero_d_vars <- rep(0, n_zero_d_vars)
}

{
model <- list()
model$obj        <- rep(0, N*(N-1)+n_borders)
model$A          <- rbind(A1, A2, A_zero_d_vars)
model$rhs        <- c(b1_alt, b2_alt, b_zero_d_vars)

model$sense      <- rep("=", nrow(model$A))
model$sense[10:18] <- ">"
model$vtype      <- "C" # Continuous

result <- gurobi(model, list(Method=2))
summary(result)
}

{
model_alt <- list()
model_alt$obj        <- rep(0, N*(N-1)+n_borders)
model_alt$A          <- rbind(A1, A2)
model_alt$rhs        <- c(b1_alt, b2_alt)

model_alt$sense      <- rep("=", nrow(model_alt$A))
model_alt$sense[10:18] <- ">"
model_alt$vtype      <- "C" # Continuous

result <- gurobi(model_alt, list(Method=2))
result
}
data.frame(deicision_vars=d_vars, optimal_sols=result$x)

d_vars_relaxed <- c(d_vars, "epsilon_o", "epsilon_p")
A_epsilon <-  matrix(0, nrow(model$A), 2)
A_epsilon[1:9, 1] <- -PO
A_epsilon[10:18, 2] <- -PS

{
model_relaxed <- list()
model_relaxed$obj      <- rep(0, N*(N-1)+2)
model_relaxed$A        <- cbind(model$A, A_epsilon)
model_relaxed$rhs      <- c(b1, b2, b_zero_d_vars)
model_relaxed$sense    <- rep("=", nrow(model_alt$A))
model_relaxed$vtype    <- "C" # Continuous
model_relaxed$ub       <- c(rep(Inf, length(d_vars)), 1, 1)

model_relaxed  <- gurobi(model_relaxed)
summary(model_relaxed)
}

{
model_relaxed2 <- list()
model_relaxed2$obj          <- rep(0, N*(N-1)+2)
model_relaxed2$A            <- cbind(model$A, A_epsilon)
model_relaxed2$rhs          <- c(b1, b2_alt, b_zero_d_vars) # seizure constrains are now >= without PS
model_relaxed2$sense        <- rep("=", nrow(model_relaxed2$A))
model_relaxed2$sense[10:18] <- ">"
model_relaxed2$vtype        <- "C" # Continuous
model_relaxed2$ub           <- c(rep(Inf, length(d_vars)), 1, 1)

model_relaxed2  <- gurobi(model_relaxed2)
summary(model_relaxed2)
}
