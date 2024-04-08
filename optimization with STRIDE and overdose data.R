# setwd("/Users/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(readxl)
library(gurobi)
library(stringi)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

{
stride <- read.csv("STRIDE_Raw.csv") %>% as_tibble %>% filter(!is.na(Seize.Year) & !is.na(Seize.Month))
stride$Nt.Wt <- as.numeric(stride$Nt.Wt)
stride$Potency <- as.numeric(stride$Potency)
stride$Post.Price <- as.numeric(stride$Post.Price)
# for (i in c(1,2,6,7)) {
#   stride[,i] <- as.factor(stride[,i]) 
# }

stride %>% pull(Drug) %>% unique %>% sort
stride %>% filter(grepl("COCAINE", Drug)) %>% pull(Drug) %>% unique

# Cocaine Data
cocaine <- stride %>%
  filter(Drug %in% c("COCAINE", "COCAINE HYDROCHLORIDE") & Seize.Year != "NULL" & Seize.Month != "NULL" & State != "NULL") %>% 
  mutate(state=as.factor(State),
         MethAcq=as.factor(MethAcq),
         Drug=as.factor(Drug),
         Seize.Year=as.numeric(Seize.Year),
         Seize.Month=as.numeric(Seize.Month),
         adjusted_price=Post.Price/(Nt.Wt*Potency/100)) %>% 
  mutate(adjusted_price=ifelse(Potency == 0, NA, adjusted_price)) %>% 
  filter(Seize.Year > 1999) %>% 
  left_join(states %>% 
              rename(state=state_abbv) %>% 
              select(state, state_name) %>% 
              unique,
            by="state") %>% 
  select(-State, -state) %>% 
  rename(state=state_name) %>% 
  relocate(state)

overdose <- read.csv("Cocaine Network Optimization/Overdose Deaths 1999-2016.csv") %>% as_tibble
}

cocaine$MethAcq %>% unique
cocaine %>% filter(MethAcq=="S") %>% pull(Post.Price) %>% table
cocaine %>% filter(MethAcq=="P" & Nt.Wt >= 5 & Nt.Wt <= 1000) %>% pull(adjusted_price) %>% sort %>% tail

cocaine %>% 
  filter(adjusted_price & Nt.Wt >= 5 & Nt.Wt <= 1000) %>% 
  ggplot() +
  geom_boxplot(aes(x=MethAcq, y=adjusted_price))

cocaine %>% 
  filter(Nt.Wt >= 5 & Nt.Wt <= 1000) %>% 
  ggplot() +
  geom_point(aes(x=Nt.Wt*Potency/100, y=adjusted_price, color=MethAcq)) +
  labs(x="Cocaine Weights")

cocaine %>% 
  filter(MethAcq != "P" & Nt.Wt >= 5 & Nt.Wt <= 1000) %>% 
  ggplot() +
  geom_point(aes(x=Nt.Wt*Potency/100, y=adjusted_price, color=MethAcq)) +
  labs(x="Cocaine Weights")

cocaine %>% 
  filter(adjusted_price < 10000 & MethAcq == "X") %>% 
  ggplot() +
  geom_point(aes(x=Nt.Wt*Potency, y=adjusted_price)) +
  labs(x="Cocaine Weights")
# Prices only with MethAcq == "P" are used

cocaine %>%
  filter(Seize.Year > 2010 &  MethAcq == "S") %>% 
  group_by(state, Seize.Year) %>% 
  summarise(n_weight=sum(!is.na(Nt.Wt)),
            total_weight=sum(Nt.Wt, na.rm=T),
            avg_weight=mean(Nt.Wt, na.rm=T),
            med_weight=median(Nt.Wt, na.rm=T),
            total_cocaine_weight=sum(Nt.Wt*Potency/100, na.rm=T),
            avg_cocaine_weight=mean(Nt.Wt*Potency/100, na.rm=T),
            med_cocaine_weight=median(Nt.Wt*Potency/100, na.rm=T)) %>% view

cocaine %>%
  filter(Seize.Year > 2010 &  MethAcq == "P" & adjusted_price < 10000) %>% 
  group_by(state, Seize.Year) %>% 
  summarise(n_price=sum(!is.na(adjusted_price)),
            avg_price=mean(adjusted_price, na.rm=T),
            med_price=median(adjusted_price, na.rm=T)) %>% view

cocaine_annual_seizures <- cocaine %>%
  filter(MethAcq == "S") %>% 
  group_by(state, Seize.Year) %>% 
  summarise(n_weight=sum(!is.na(Nt.Wt)),
            total_weight=sum(Nt.Wt, na.rm=T),
            avg_weight=mean(Nt.Wt, na.rm=T),
            med_weight=median(Nt.Wt, na.rm=T),
            total_cocaine_weight=sum(Nt.Wt*Potency/100, na.rm=T),
            avg_cocaine_weight=mean(Nt.Wt*Potency/100, na.rm=T),
            med_cocaine_weight=median(Nt.Wt*Potency/100, na.rm=T),)

cocaine_annual_prices <- cocaine %>% 
  filter(MethAcq == "P" & Nt.Wt >= 5 & Nt.Wt <= 1000) %>% 
  group_by(state, Seize.Year) %>% 
  summarise(n_price=sum(!is.na(adjusted_price)),
            avg_price=mean(adjusted_price, na.rm=T),
            med_price=median(adjusted_price, na.rm=T))

cocaine_annual_purities <- cocaine %>% 
  filter(MethAcq == "P" & Nt.Wt >= 5 & Nt.Wt <= 1000) %>% 
  group_by(state, Seize.Year) %>% 
  summarise(n_purity=sum(!is.na(Potency)),
            avg_purity=mean(Potency, na.rm=T),
            med_purity=median(Potency, na.rm=T))
ex_year <- 2013
cocaine_annual_seizures %>% filter(Seize.Year == ex_year) %>% 
  select(state, Seize.Year, n_weight, total_weight, total_cocaine_weight) %>%  arrange(state) %>% as.data.frame
cocaine_annual_prices %>% filter(Seize.Year == ex_year) %>% arrange(state) %>% as.data.frame
# 1: California  2: Texas     3: Florida
# 4: Nevada      5: Missouri  6: West Virginia
# 7: Washington  8: Illinois  9: New York

nine_states <- c("California", "Texas", "Florida", "Nevada", "Missouri", "West Virginia", "Washington", "Illinois", "New York")
cocaine_annual_seizures %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% 
  select(state, Seize.Year, n_weight, total_weight, total_cocaine_weight) %>%  arrange(state) %>% as.data.frame
cocaine_annual_purities %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% arrange(state) %>% as.data.frame
overdose %>% filter(year == ex_year & state %in% nine_states) %>% arrange(state) %>% as.data.frame

cocaine %>%
  filter(MethAcq=="P" & Nt.Wt >= 5 & Nt.Wt <= 1000 & Seize.Year == 2013 & state %in% nine_states) %>% 
  ggplot() +
  geom_boxplot(aes(x=state, y=Potency))
  



# Optimization
ex_year <- 2013
nine_states_data <- cocaine_annual_prices %>%
  filter(Seize.Year == ex_year & state %in% nine_states) %>% 
  left_join(cocaine_annual_seizures %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% select(-Seize.Year), by="state") %>% 
  left_join(cocaine_annual_purities %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% select(-Seize.Year), by="state") %>% 
  left_join(overdose %>% filter(year == ex_year & state %in% nine_states) %>% select(-year), by="state")
nine_states_data <- nine_states_data[c(1,7,2, 5,4,9, 8,3,6),]
matrix(nine_states_data$med_price, 3, 3, byrow=T)
matrix(nine_states_data$med_purity, 3, 3, byrow=T)
# how to set source states? using price or purity?

# sources by price
N <- nrow(nine_states_data) # Total number of nodes
source_index <- c(1,2,8)
n_sources <- length(source_index) # < 8, Node 5 and 8 can't be a border point (in-land)
d_vars <- expand.grid(1:N, 1:N)[-(1+(N+1)*(0:(N-1))),] %>% mutate(name=paste0("w_", Var2, Var1)) %>% pull(name)
d_vars <- c(d_vars, paste0("S", source_index))
n_d_vars <- length(d_vars)
PO <- 0.1 # overdose proportional parameter
PS <- 5 # seizure proportionality parameter

{
  seizure <- nine_states_data$total_weight
  O <- nine_states_data$deaths
  price <- nine_states_data$med_price
  purity <- nine_states_data$med_purity
  
  # Constraints
  # Conservation of flow/consumption
  A1 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_ref <- matrix(0, N, N)
    A_ref[i,] <- -1
    A_ref[, i] <- 1
    S_vec <- rep(-O[i]/sum(O), n_sources)
    if (i %in% source_index) S_vec[which(source_index == i)] <- 1 + S_vec[which(source_index == i)]
      
    A1[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
  }
  b1 <- rep(0, N)
  
  A1_alt <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_ref <- matrix(0, N, N)
    A_ref[i,] <- -1
    A_ref[, i] <- 1
    S_vec <- rep(0, n_sources)
    if (i %in% source_index) S_vec[which(source_index == i)] <- 1
    
    A1_alt[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
  }
  b1_alt <- O
  
  # State seizure
  A2 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_ref <- matrix(0, N, N)
    A_ref[, i] <- 1
    S_vec <- rep(0, n_sources)
    if (i %in% source_index) S_vec[which(source_index == i)] <- 1
    A2[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
  }
  b2 <- PS*seizure
  b2_alt <- seizure
  
  # Directionality of flow
  less_than_price <- expand.grid(price, price)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 < Var1) %>% pull(less_than) %>% which
  less_than_purity <- expand.grid(purity, purity)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 < Var1) %>% pull(less_than) %>% which
  zero_d_vars_index <- unique(c(less_than_price, less_than_purity)) %>% sort
  n_zero_d_vars <- length(zero_d_vars_index)
  A_zero_d_vars <- matrix(0, n_zero_d_vars, n_d_vars)
  for (i in 1:n_zero_d_vars) {
    A_zero_d_vars[i, zero_d_vars_index[i]] <- 1
  }
  b_zero_d_vars <- rep(0, n_zero_d_vars)
}

{# Ax = PS*seizure with directionality
  model <- list()
  model$obj        <- rep(0, n_d_vars)
  model$A          <- rbind(A1, A2, A_zero_d_vars)
  model$rhs        <- c(b1, b2, b_zero_d_vars)
  
  model$sense      <- rep("=", nrow(model$A))
  model$vtype      <- "C" # Continuous
  
  result <- gurobi(model, list(Method=2))
  result
}# infeasible

{# Ax >= seizure with directionality
  model <- list()
  model$obj        <- rep(0, n_d_vars)
  model$A          <- rbind(A1, A2, A_zero_d_vars)
  model$rhs        <- c(b1, b2_alt, b_zero_d_vars)
  
  model$sense      <- rep(">", nrow(model$A))
  model$sense[10:18] <- ">"
  model$vtype      <- "C" # Continuous
  
  result <- gurobi(model, list(Method=2))
  result
}# infeasible

{# Ax = PS*seizure without directionality
  model_alt <- list()
  model_alt$obj        <- rep(0, n_d_vars)
  model_alt$A          <- rbind(A1, A2)
  model_alt$rhs        <- c(b1, b2)
  
  model_alt$sense      <- rep("=", nrow(model_alt$A))
  model_alt$vtype      <- "C" # Continuous
  
  result <- gurobi(model_alt, list(Method=2))
  result
}# infeasible

{# Ax >= seizure without directionality
  model_alt <- list()
  model_alt$obj        <- rep(0, n_d_vars)
  model_alt$A          <- rbind(A1, A2)
  model_alt$rhs        <- c(b1, b2_alt)
  
  model_alt$sense      <- rep("=", nrow(model_alt$A))
  model_alt$sense[10:18] <- ">"
  model_alt$vtype      <- "C" # Continuous
  
  params                <- list()
  params$Method         <- 2
  params$PoolSolutions  <- 1024
  params$PoolSearchMode <- 2
  
  result <- gurobi(model_alt, params)
  result
}# feasible
data.frame(deicision_vars=d_vars, optimal_sols=result$x)
coordinates <- data.frame(i=1:9, j=1:9, x=rep(1:3, 3), y=rep(1:3, each=3))
opt_sol <- data.frame(deicision_vars=d_vars, optimal_sols=result$x) %>% filter(optimal_sols > 0)
total_flow <- sum(opt_sol$optimal_sols[1:9])
opt_sol$i <- c(opt_sol$deicision_vars[1:9] %>% substr(3, 3) %>% as.numeric,
               opt_sol$deicision_vars[10:11] %>% substr(2, 2) %>% as.numeric)
opt_sol$j <- c(opt_sol$deicision_vars[1:9] %>% substr(4, 4) %>% as.numeric,
               opt_sol$deicision_vars[10:11] %>% substr(2, 2) %>% as.numeric)
opt_sol <- opt_sol %>%
  left_join(coordinates %>% select(-j), by="i") %>% 
  left_join(coordinates %>% select(-i), by="j") %>% 
  rename(source_x=x.x,
         source_y=y.x,
         destination_x=x.y,
         destination_y=y.y)

coordinates %>% ggplot(aes(x=x, y=y)) +
  geom_point() +
  geom_text(label=coordinates$i, nudge_x = 0.05, nudge_y = 0) +
  labs(x="", y="") +
  geom_segment(data=opt_sol[1:9,],
               aes(x=source_x, 
                   y=source_y, 
                   xend=destination_x,
                   yend=destination_y),
               linewidth = 5*opt_sol$optimal_sols[1:9]/total_flow,
               arrow=arrow(angle=10,
                           # length=unit(0.1, "cm"),
                           type="closed")
  ) +
  geom_text(data=opt_sol[10:11,],
            aes(x=source_x, y=source_y),
            label=opt_sol[10:11,]$optimal_sols %>% round(1),
            color="red",
            nudge_x = c(0.1, 0), nudge_y = c(-0.1, 0.1))

# sources by purity
N <- nrow(nine_states_data) # Total number of nodes
source_index <- c(1,4,7)
n_sources <- length(source_index) # < 8, Node 5 and 8 can't be a border point (in-land)
d_vars <- expand.grid(1:N, 1:N)[-(1+(N+1)*(0:(N-1))),] %>% mutate(name=paste0("w_", Var2, Var1)) %>% pull(name)
d_vars <- c(d_vars, paste0("S", source_index))
n_d_vars <- length(d_vars)
PO <- 0.1 # overdose proportional parameter
PS <- 5 # seizure proportionality parameter


{
  seizure <- nine_states_data$total_weight
  O <- nine_states_data$deaths
  price <- nine_states_data$med_price
  purity <- nine_states_data$med_purity
  
  # Constraints
  # Conservation of flow/consumption
  A1 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_ref <- matrix(0, N, N)
    A_ref[i,] <- -1
    A_ref[, i] <- 1
    S_vec <- rep(-O[i]/sum(O), n_sources)
    if (i %in% source_index) S_vec[which(source_index == i)] <- 1 + S_vec[which(source_index == i)]
    
    A1[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
  }
  b1 <- rep(0, N)
  
  A1_alt <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_ref <- matrix(0, N, N)
    A_ref[i,] <- -1
    A_ref[, i] <- 1
    S_vec <- rep(0, n_sources)
    if (i %in% source_index) S_vec[which(source_index == i)] <- 1
    
    A1_alt[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
  }
  b1_alt <- O
  
  # State seizure
  A2 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_ref <- matrix(0, N, N)
    A_ref[, i] <- 1
    S_vec <- rep(0, n_sources)
    if (i %in% source_index) S_vec[which(source_index == i)] <- 1
    A2[i,] <- c(as.vector(t(A_ref))[-(1+(N+1)*(0:(N-1)))], S_vec)
  }
  b2 <- PS*seizure
  b2_alt <- seizure
  
  # Directionality of flow
  less_than_price <- expand.grid(price, price)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 < Var1) %>% pull(less_than) %>% which
  less_than_purity <- expand.grid(purity, purity)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 < Var1) %>% pull(less_than) %>% which
  zero_d_vars_index <- unique(c(less_than_price, less_than_purity)) %>% sort
  n_zero_d_vars <- length(zero_d_vars_index)
  A_zero_d_vars <- matrix(0, n_zero_d_vars, n_d_vars)
  for (i in 1:n_zero_d_vars) {
    A_zero_d_vars[i, zero_d_vars_index[i]] <- 1
  }
  b_zero_d_vars <- rep(0, n_zero_d_vars)
}

{# Ax >= seizure without directionality
  model_alt <- list()
  model_alt$obj        <- rep(0, n_d_vars)
  model_alt$A          <- rbind(A1, A2)
  model_alt$rhs        <- c(b1, b2_alt)
  
  model_alt$sense      <- rep("=", nrow(model_alt$A))
  model_alt$sense[10:18] <- ">"
  model_alt$vtype      <- "C" # Continuous
  
  params                <- list()
  params$Method         <- 2
  params$PoolSolutions  <- 1024
  params$PoolSearchMode <- 2
  
  result <- gurobi(model_alt, params)
  result
}# feasible
data.frame(deicision_vars=d_vars, optimal_sols=result$x)
coordinates <- data.frame(i=1:9, j=1:9, x=rep(1:3, 3), y=rep(1:3, each=3))
opt_sol <- data.frame(deicision_vars=d_vars, optimal_sols=result$x) %>% filter(optimal_sols > 0)
total_flow <- sum(opt_sol$optimal_sols[1:9])
opt_sol$i <- c(opt_sol$deicision_vars[1:10] %>% substr(3, 3) %>% as.numeric,
               opt_sol$deicision_vars[11] %>% substr(2, 2) %>% as.numeric)
opt_sol$j <- c(opt_sol$deicision_vars[1:10] %>% substr(4, 4) %>% as.numeric,
               opt_sol$deicision_vars[11] %>% substr(2, 2) %>% as.numeric)
opt_sol <- opt_sol %>%
  left_join(coordinates %>% select(-j), by="i") %>% 
  left_join(coordinates %>% select(-i), by="j") %>% 
  rename(source_x=x.x,
         source_y=y.x,
         destination_x=x.y,
         destination_y=y.y)

coordinates %>% ggplot(aes(x=x, y=y)) +
  geom_point() +
  geom_text(label=coordinates$i, nudge_x = 0.05, nudge_y = 0) +
  labs(x="", y="") +
  geom_segment(data=opt_sol[1:10,],
               aes(x=source_x, 
                   y=source_y, 
                   xend=destination_x,
                   yend=destination_y),
               linewidth = 5*opt_sol$optimal_sols[1:10]/total_flow,
               arrow=arrow(angle=10,
                           # length=unit(0.1, "cm"),
                           type="closed")
  ) +
  geom_text(data=opt_sol[11,],
            aes(x=source_x, y=source_y),
            label=opt_sol[11,]$optimal_sols %>% round(1),
            color="red",
            nudge_x = 0.1, nudge_y = 0.1)




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
