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
  stride <- read.csv("STRIDE_Raw.csv") %>% as_tibble %>% filter(!is.na(Seize.Year) & !is.na(Seize.Month) & !(State %in% c("AK", "HI")))
  stride$Nt.Wt <- as.numeric(stride$Nt.Wt)
  stride$Potency <- as.numeric(stride$Potency)
  stride$Post.Price <- as.numeric(stride$Post.Price)
  # for (i in c(1,2,6,7)) {
  #   stride[,i] <- as.factor(stride[,i]) 
  # }
  
  stride %>% pull(Drug) %>% unique %>% sort
  stride %>% filter(grepl("COCAINE", Drug)) %>% pull(Drug) %>% unique
  
  neighbor_xlsx <- read_xlsx("Cocaine Network Optimization/states neighborhood.xlsx")
  neighbor <- list()
  for (i in 1:length(neighbor_xlsx$State)) {
    neighbor[[neighbor_xlsx$State[i]]] <- strsplit(neighbor_xlsx$Bordering_States[i], ", ")[[1]]
  }

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
            med_cocaine_weight=median(Nt.Wt*Potency/100, na.rm=T))

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

}

# Optimization
{
ex_year <- 2012
state_names <- names(neighbor)
states_data <- cocaine_annual_prices %>%
  filter(Seize.Year == ex_year & state %in% state_names) %>% 
  arrange(state) %>% 
  left_join(cocaine_annual_seizures %>% filter(Seize.Year == ex_year & state %in% state_names) %>% select(-Seize.Year), by="state") %>% 
  left_join(cocaine_annual_purities %>% filter(Seize.Year == ex_year & state %in% state_names) %>% select(-Seize.Year), by="state") %>% 
  left_join(overdose %>% filter(year == ex_year & state %in% state_names) %>% select(-year), by="state") %>% 
  filter(!is.na(sum(n_price, n_weight, deaths)))

# all possible sources
N <- nrow(states_data) # Total number of nodes
state_index <- tibble(state=states_data$state, state_index=1:N) %>% 
  arrange(state) %>% 
  left_join(states %>% 
              rename(state=state_name) %>% 
              group_by(state) %>% 
              summarise(long=mean(long), lat=mean(lat)), by="state")

d_vars <- c()
for (i in 1:N) {
  bordering_states_index <- state_index %>%
    filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(state_index)
  d_vars <- c(d_vars, paste0("w", i, ",", bordering_states_index, "."))
}

d_vars <- c(d_vars, paste0("S", 1:N))
d_vars <- c(d_vars, paste0("x", 1:N))
n_d_vars <- length(d_vars)
PO <- 100 # overdose proportional parameter
# PS <- 5 # seizure proportionality parameter
epsilon <- 0.9
}

{
  # seizure <- states_data$total_weight
  seizure <- states_data$max_weight
  log_seizure <- log(seizure)
  O <- states_data$deaths
  price <- states_data$med_price
  purity <- states_data$med_purity

    # Constraints
  # Conservation of flow/consumption
  A1 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A1[i, which(str_match(d_vars, paste0("w", i, ","))[,1] == paste0("w", i, ","))] <- -1
    A1[i, which(str_match(d_vars, paste0(",", i, "."))[,1] == paste0(",", i, "."))] <- 1
    A1[i, n_d_vars-2*N+i] <- 1
  }
  b1 <- PO*O/sum(O)
  b1_alt <- PO*O/sum(O)
  
  # State seizure
  A2 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A2[i, which(str_match(d_vars, paste0(",", i, "."))[,1] == paste0(",", i, "."))] <- 1
    A2[i, n_d_vars-2*N+i] <- 1
  }
  b2_ratio <- 100*seizure/sum(seizure)
  b2_log_ratio <- 100*log_seizure/sum(log_seizure)
  
  # Directionality of flow
  less_than_price <- c()
  greater_than_purity <- c()
  for (i in 1:(n_d_vars - 2*N)) {
    var <- d_vars[i]
    comma_index <- which(strsplit(var, "")[[1]] == ",")
    dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
    source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
    destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
    if (price[destination_index] < price[source_index]) less_than_price <- c(less_than_price, i)
    if (purity[destination_index] > purity[source_index]) greater_than_purity <- c(greater_than_purity, i)
  }
  zero_d_vars_index <- unique(c(less_than_price, greater_than_purity)) %>% sort
  n_zero_d_vars <- length(zero_d_vars_index)
  
  # Source indicator
  A3 <- matrix(0, 2*N, n_d_vars)
  for (i in 1:N) {
    A3[i, which(d_vars %in% c(paste0("S", i), paste0("x", i)))] <- c(1, -1)
    A3[i+N, which(d_vars %in% c(paste0("S", i), paste0("x", i)))] <- c(-0.01, 1)
  }
  b3 <- rep(0, 2*N)
  
  A_source_limit <- matrix(0, 1, n_d_vars)
  A_source_limit[1, grep("S", d_vars)] <- 1
  b_source_limit <- 100
}


{# Ax >= seizure with directionality and obj. func.
  
  model_obj <- list()
  model_obj$obj        <- c(rep(0, n_d_vars - n_zero_d_vars - N), rep(1, N))
  model_obj$modelsense <- "min"
  model_obj$A          <- rbind(A1[,-zero_d_vars_index], A2[,-zero_d_vars_index], A3[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_obj$sense      <- rep(">", nrow(model_obj$A))
  # model_obj$sense[1:N] <- "="
  model_obj$sense[nrow(model_obj$A)] <- "="
  model_obj$vtype      <- c(rep("C", length(model_obj$obj) - N), rep("B", N)) # Continuous and binary
  model_obj$rhs        <- c(b1, b2_ratio, b3, b_source_limit)
  result_obj_sum       <- gurobi(model_obj, list(Method=2))
  result_obj_sum
  
  # iterate multiple epsilons for the directionality of flow
  epsilton_results <- data.frame()
  for (epsilon in seq(0.5, 1, by=0.05)) {
    less_than_price <- c()
    greater_than_purity <- c()
    for (i in 1:(n_d_vars - 2*N)) {
      var <- d_vars[i]
      comma_index <- which(strsplit(var, "")[[1]] == ",")
      dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
      source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
      destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
      if (price[destination_index] < price[source_index]*epsilon) less_than_price <- c(less_than_price, i)
      if (purity[destination_index]*epsilon > purity[source_index]) greater_than_purity <- c(greater_than_purity, i)
    }
    zero_d_vars_index <- unique(c(less_than_price, greater_than_purity)) %>% sort
    n_zero_d_vars <- length(zero_d_vars_index)
    
    model_epsilon <- list()
    model_epsilon$obj        <- c(rep(0, n_d_vars - n_zero_d_vars - N), rep(1, N))
    model_epsilon$modelsense <- "min"
    model_epsilon$A          <- rbind(A1[,-zero_d_vars_index], A2[,-zero_d_vars_index], A3[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
    model_epsilon$sense      <- rep(">", nrow(model_epsilon$A))
    # model_epsilon$sense[1:N] <- "="
    model_epsilon$sense[nrow(model_epsilon$A)] <- "="
    model_epsilon$vtype      <- c(rep("C", length(model_epsilon$obj) - N), rep("B", N)) # Continuous and binary
    model_epsilon$rhs        <- c(b1, b2_ratio, b3, b_source_limit)
    result_epsilon_sum       <- gurobi(model_epsilon, list(Method=2))
    epsilton_results <- rbind(epsilton_results, c(epsilon, result_epsilon_sum$status))
  }
  names(epsilton_results) <- c("epsilon", "feasibility")
  epsilton_results
  
  # iterate multiple values for PS
  PS_results <- data.frame()
  for (PS in seq(1, 100, by=1)) {
    b2 <- PS*seizure
    b2_log <- PS*log_seizure
    model_obj$rhs <- c(b1, b2, b3, b_source_limit)
    result_obj_seizure <- gurobi(model_obj, list(Method=2))
    model_obj$rhs <- c(b1, b2_log, b3, b_source_limit)
    result_obj_log_seizure <- gurobi(model_obj, list(Method=2))
    PS_results <- rbind(PS_results, c(PS, result_obj_seizure$status, result_obj_log_seizure$status))
  }
  names(PS_results) <- c("PS", "seizure_solution", "log_seizure_solution")
  PS_results
}


## results with b2_ratio = 100*seizure/sum(seizure)
find_state_index <- function(d_var) {
  if (grepl("w", d_var)) {
    comma_index <- which(str_split(d_var, "")[[1]] == ",")
    dot_index <- which(str_split(d_var, "")[[1]] == ".")
    return(c(substr(d_var, 2, comma_index-1), substr(d_var, comma_index+1, dot_index-1)) %>% as.numeric)
  }
  
  if (grepl("S", d_var)) {
    return(rep(substr(d_var, 2, str_length(d_var)), 2) %>% as.numeric)
  }
}
find_state_index <- Vectorize(find_state_index)

result_obj_sum$pool %>% length
for (sol_i in 1:(result_obj_sum$pool %>% length)) {
  opt_sol <- data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_obj_sum$pool[[sol_i]]$xn) %>%
    filter(optimal_sols > 0 & !grepl("x", deicision_vars))
  opt_sol
  # write.csv(opt_sol,
  #           paste0("Cocaine Network Optimization/Results/optimal networks (epsilon=0.5, sum_x=", result_obj_sum$pool[[sol_i]]$objval, ").csv"),
  #           row.names=F)
  total_flow <- opt_sol %>%
    filter(grepl("S", deicision_vars)) %>%
    pull(optimal_sols) %>% sum
  total_flow
  
  opt_sol <- cbind(opt_sol, find_state_index(opt_sol$deicision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  opt_sol$long_i <- state_index$long[opt_sol$i]
  opt_sol$lat_i <- state_index$lat[opt_sol$i]
  opt_sol$long_j <- state_index$long[opt_sol$j]
  opt_sol$lat_j <- state_index$lat[opt_sol$j]
  
  result_flow_index <- grep("w", opt_sol$deicision_vars)
  result_source_index <- grep("S", opt_sol$deicision_vars)
  overdose_map_year <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, deaths),
                                 by="state")
  overdose_map_year %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=deaths),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
    coord_quickmap() +
    labs(x="", y="") +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(data=state_index,
               aes(x=long, y=lat)) +
    geom_segment(data=opt_sol[result_flow_index,],
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 5*opt_sol$optimal_sols[result_flow_index]/total_flow,
                 arrow=arrow(angle=10,
                             length=unit(0.3, "cm"),
                             type="closed")
    ) +
    geom_text(data=opt_sol[result_source_index,],
              aes(x=long_i, y=lat_i),
              label=opt_sol[result_source_index,]$optimal_sols %>% round(1),
              color="red",
              nudge_x = c(0.1, 0), nudge_y = c(-0.5, 0)) -> optimal_network_map
  # ggsave(paste0("Cocaine Network Optimization/Results/optimal networks (epsilon=0.5, sum_x=", result_obj_sum$pool[[sol_i]]$objval, ").png"),
  #        optimal_network_map, scale=1.5)
}
  
A1[, -zero_d_vars_index] %*% matrix(result_obj_sum$pool[[sol_i]]$xn, ncol=1)

cocaine %>%
  filter(Seize.Year == 2012 & !is.na(adjusted_price) & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  select(-MethAcq, -Drug) %>% 
  relocate(state, Seize.Year, Seize.Month) %>% 
  print(n=20)

cocaine %>%
  filter(Seize.Year == 2012 & !is.na(adjusted_price) & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  select(-MethAcq, -Drug) %>% 
  relocate(state, Seize.Year, Seize.Month) %>% 
  group_by(state) %>% 
  summarise(price_sd_max_wt1000=sd(adjusted_price)) %>% 
  print(n=43)

cocaine %>%
  filter(Seize.Year == 2012 & !is.na(adjusted_price) & Nt.Wt >= 5 & Nt.Wt <= 100) %>%
  select(-MethAcq, -Drug) %>% 
  relocate(state, Seize.Year, Seize.Month) %>% 
  group_by(state) %>% 
  summarise(price_sd_max_wt100=sd(adjusted_price)) %>% 
  print(n=43)

  # same results allowing S > 0 only for CA, FL, TX 
n_states <- nrow(state_index)
three_sources <- state_index %>% filter(state %in% c("California", "Florida", "Texas")) %>% pull(state_index)
sum(grepl("w", d_vars))
omitted_index <- c(zero_d_vars_index, (165:n_d_vars)[!((165:n_d_vars) %in% which(d_vars %in% paste0("S", three_sources)))])

model_obj <- list()
model_obj$obj        <- c(rep(0, n_d_vars - length(omitted_index)))
model_obj$modelsense <- "min"
model_obj$A          <- rbind(A1[,-omitted_index], A2[,-omitted_index], A3[,-omitted_index], A_source_limit[,-omitted_index])
model_obj$sense      <- rep(">", nrow(model_obj$A))
# model_obj$sense[1:N] <- "="
model_obj$sense[nrow(model_obj$A)] <- "="
model_obj$vtype      <- "C"
model_obj$rhs        <- c(b1, b2_ratio, b3, b_source_limit)
result_obj_sum       <- gurobi(model_obj, list(Method=2))
result_obj_sum


## results for log seizures
model_obj$obj <- rep(0, n_d_vars - n_zero_d_vars)

PS <- 1
b2_log <- PS*log_seizure
model_obj$rhs <- c(b1, b2_log, b3, b_source_limit)
result_obj_log_seizure <- gurobi(model_obj, list(Method=2))
result_obj_log_seizure$pool
opt_sol <- data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_obj_log_seizure$pool[[2]]$xn) %>%
  filter(optimal_sols > 0 & !grepl("x", deicision_vars))
opt_sol

A1[,-zero_d_vars_index] %*% matrix(result_obj_log_seizure$pool[[2]]$xn, ncol=1)