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
VSRR <- read.csv("Cocaine Network Optimization/VSRR_Provisional_Drug_Overdose_Death_Counts (2015-2023).csv") %>%
  as_tibble %>% 
  mutate(state=State,
         Data.Value=as.numeric(Data.Value)) %>% 
  left_join(states %>% 
              rename(state=state_abbv) %>% 
              select(state, state_name) %>% 
              unique,
            by="state") %>% 
  select(-State, -state) %>% 
  rename(state=state_name) %>% 
  relocate(state)
VSRR_cocaine <- VSRR %>% 
  filter(Indicator == "Cocaine (T40.5)") %>% 
  group_by(state, Year) %>% 
  summarize(overdose_death=sum(Data.Value, na.rm=T))
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
  filter(Seize.Year == 2012 &  MethAcq == "P" & adjusted_price < 10000) %>% 
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
            max_weight=max(Nt.Wt, na.rm=T),
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
  filter(MethAcq == "P" & Nt.Wt >= 5 & Nt.Wt <= 100) %>% 
  group_by(state, Seize.Year) %>% 
  summarise(n_purity=sum(!is.na(Potency)),
            avg_purity=mean(Potency, na.rm=T),
            med_purity=median(Potency, na.rm=T))
ex_year <- 2012
cocaine_annual_seizures %>% filter(Seize.Year == ex_year) %>% 
  select(state, Seize.Year, n_weight, total_weight, total_cocaine_weight) %>%  arrange(state) %>% as.data.frame
cocaine_annual_prices %>% filter(Seize.Year == ex_year) %>% arrange(state) %>% as.data.frame
# 1: California  2: Texas     3: Florida
# 4: Nevada      5: Missouri  6: West Virginia
# 7: Washington  8: Illinois  9: New York

nine_states <- c("California", "Texas", "Florida", "Nevada", "Colorado", "Missouri", "Washington", "Illinois", "New York")
cocaine_annual_seizures %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% 
  select(state, Seize.Year, n_weight, total_weight, total_cocaine_weight) %>%  arrange(state) %>% as.data.frame
cocaine_annual_purities %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% arrange(state) %>% as.data.frame
overdose %>% filter(year == ex_year & state %in% nine_states) %>% arrange(state) %>% as.data.frame
VSRR_cocaine %>% filter(Year == 2015 & state %in% nine_states) %>% arrange(state) %>% as.data.frame
VSRR_cocaine %>% filter(Year == 2015) %>% arrange(state) %>% as.data.frame
VSRR_cocaine %>% filter(Year == 2015) %>%
  left_join(cocaine_annual_prices %>% filter(Seize.Year == 2012)) %>% arrange(desc(n_price), desc(overdose_death)) %>% 
  as.data.frame

cocaine %>%
  filter(MethAcq=="P" & Nt.Wt >= 5 & Nt.Wt <= 1000 & Seize.Year == ex_year & state %in% nine_states) %>% 
  ggplot() +
  geom_boxplot(aes(x=state, y=adjusted_price))

cocaine %>%
  filter(MethAcq=="P" & Nt.Wt >= 5 & Nt.Wt <= 10 & Seize.Year == ex_year & state %in% nine_states) %>% 
  ggplot() + ggtitle("5 ~ 10g") +
  geom_boxplot(aes(x=state, y=Potency))
  
cocaine %>%
    filter(MethAcq=="P" & Nt.Wt >= 5 & Nt.Wt <= 100 & Seize.Year == ex_year & state %in% nine_states) %>% 
    ggplot() + ggtitle("5 ~ 100g") +
  geom_boxplot(aes(x=state, y=Potency))


# Optimization
ex_year <- 2012
nine_states_data <- cocaine_annual_prices %>%
  filter(Seize.Year == ex_year & state %in% nine_states) %>% 
  left_join(cocaine_annual_seizures %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% select(-Seize.Year), by="state") %>% 
  left_join(cocaine_annual_purities %>% filter(Seize.Year == ex_year & state %in% nine_states) %>% select(-Seize.Year), by="state") %>% 
  left_join(overdose %>% filter(year == ex_year & state %in% nine_states) %>% select(-year), by="state")
nine_states_data <- nine_states_data[c(1,7,2, 5,4,9, 8,3,6),]
matrix(nine_states_data$med_price, 3, 3, byrow=T)
matrix(nine_states_data$med_purity, 3, 3, byrow=T)
matrix(nine_states_data$avg_price, 3, 3, byrow=T)
matrix(nine_states_data$avg_purity, 3, 3, byrow=T)
# how to set source states? using price or purity?

# all possible sources
N <- nrow(nine_states_data) # Total number of nodes
source_index <- 1:9
n_sources <- length(source_index) # < 8, Node 5 and 8 can't be a border point (in-land)
d_vars <- expand.grid(1:N, 1:N)[-(1+(N+1)*(0:(N-1))),] %>% mutate(name=paste0("w_", Var2, Var1)) %>% pull(name)
d_vars <- c(d_vars, paste0("S", source_index))
n_d_vars <- length(d_vars)
PO <- 0.1 # overdose proportional parameter
PS <- 5 # seizure proportionality parameter

{
  # seizure <- nine_states_data$total_weight
  seizure <- nine_states_data$max_weight
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
  less_than_price <- expand.grid(price, price)[-(1+(N+1)*(0:(N-1))),] %>% mutate(less_than=Var2 > Var1) %>% pull(less_than) %>% which
  greater_than_purity <- expand.grid(purity, purity)[-(1+(N+1)*(0:(N-1))),] %>% mutate(greater_than=Var2 < Var1) %>% pull(greater_than) %>% which
  zero_d_vars_index <- unique(c(less_than_price, greater_than_purity)) %>% sort
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
  model$obj        <- rep(0, n_d_vars - n_zero_d_vars)
  model$A          <- rbind(A1[,-zero_d_vars_index], A2[,-zero_d_vars_index])
  model$rhs        <- c(b1, b2_alt)
  
  model$sense      <- rep("=", nrow(model$A))
  model$sense[10:18] <- ">"
  model$vtype      <- "C" # Continuous
  
  result <- gurobi(model, list(Method=2))
  result
}# feasible

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
  
  result_alt <- gurobi(model_alt, params)
  result_alt
}# feasible

data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result$x)
coordinates <- data.frame(i=1:9, j=1:9, x=rep(1:3, 3), y=rep(1:3, each=3))
opt_sol <- data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result$x) %>% filter(optimal_sols > 0)
total_flow <- sum(opt_sol$optimal_sols[1:9])

flow_index <- grep("w", opt_sol$deicision_vars)
source_index <- grep("S", opt_sol$deicision_vars)

opt_sol$i <- c(opt_sol$deicision_vars[flow_index] %>% substr(3, 3) %>% as.numeric,
               opt_sol$deicision_vars[source_index] %>% substr(2, 2) %>% as.numeric)
opt_sol$j <- c(opt_sol$deicision_vars[flow_index] %>% substr(4, 4) %>% as.numeric,
               opt_sol$deicision_vars[source_index] %>% substr(2, 2) %>% as.numeric)
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
  geom_segment(data=opt_sol[flow_index,],
               aes(x=source_x, 
                   y=source_y, 
                   xend=destination_x,
                   yend=destination_y),
               linewidth = 5*opt_sol$optimal_sols[flow_index]/total_flow,
               arrow=arrow(angle=10,
                           # length=unit(0.1, "cm"),
                           type="closed")
  ) +
  geom_text(data=opt_sol[source_index,],
            aes(x=source_x, y=source_y),
            label=opt_sol[source_index,]$optimal_sols %>% round(1),
            color="red",
            nudge_x = c(0.1, 0), nudge_y = c(-0.1, 0.1))

{# Ax >= seizure with directionality and obj. func.
  model_obj_sum <- list()
  model_obj_sum$obj        <- c(rep(0, n_d_vars - n_zero_d_vars - n_sources), rep(1, n_sources))
  model_obj_sum$modelsense <- "min"
  model_obj_sum$A          <- rbind(A1[,-zero_d_vars_index], A2[,-zero_d_vars_index])
  model_obj_sum$rhs        <- c(b1, b2_alt)
  model_obj_sum$sense      <- rep("=", nrow(model_obj_sum$A))
  model_obj_sum$sense[10:18] <- ">"
  model_obj_sum$vtype      <- "C" # Continuous
  
  result_obj_sum <- gurobi(model_obj_sum, list(Method=2))
  result_obj_sum
}

data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_obj_sum$x)
coordinates <- data.frame(i=1:9, j=1:9, x=rep(1:3, 3), y=rep(1:3, each=3))
opt_sol <- data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_obj_sum$x) %>% filter(optimal_sols > 0)
total_flow <- sum(opt_sol$optimal_sols[1:9])

flow_index <- grep("w", opt_sol$deicision_vars)
source_index <- grep("S", opt_sol$deicision_vars)

opt_sol$i <- c(opt_sol$deicision_vars[flow_index] %>% substr(3, 3) %>% as.numeric,
               opt_sol$deicision_vars[source_index] %>% substr(2, 2) %>% as.numeric)
opt_sol$j <- c(opt_sol$deicision_vars[flow_index] %>% substr(4, 4) %>% as.numeric,
               opt_sol$deicision_vars[source_index] %>% substr(2, 2) %>% as.numeric)
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
  geom_segment(data=opt_sol[flow_index,],
               aes(x=source_x, 
                   y=source_y, 
                   xend=destination_x,
                   yend=destination_y),
               linewidth = 5*opt_sol$optimal_sols[flow_index]/total_flow,
               arrow=arrow(angle=10,
                           # length=unit(0.1, "cm"),
                           type="closed")
  ) +
  geom_text(data=opt_sol[source_index,],
            aes(x=source_x, y=source_y),
            label=opt_sol[source_index,]$optimal_sols %>% round(1),
            color="red",
            nudge_x = c(0.1, 0), nudge_y = c(-0.1, 0.1))
