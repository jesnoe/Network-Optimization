# setwd("/Users/R")
# setwd("C:/Users/Users/Documents/R")
library(readxl)
library(gurobi)
library(stringi)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

{
stride <- read.csv("Cocaine Network Optimization/STRIDE_Raw.csv") %>% as_tibble %>% filter(!is.na(Seize.Year) & !is.na(Seize.Month) & !(State %in% c("AK", "HI")))
stride$Nt.Wt <- as.numeric(stride$Nt.Wt)
stride$Potency <- as.numeric(stride$Potency)
stride$Post.Price <- as.numeric(stride$Post.Price)
# for (i in c(1,2,6,7)) {
#   stride[,i] <- as.factor(stride[,i]) 
# }
stride$State <- ifelse(stride$State == "NB", "NE", stride$State)

neighbor_xlsx <- read_xlsx("Cocaine Network Optimization/states neighborhood.xlsx")
neighbor <- list()
for (i in 1:length(neighbor_xlsx$State)) {
  state_i <- neighbor_xlsx$State[i]
  neighbor[[state_i]] <- strsplit(neighbor_xlsx$Bordering_States[i], ", ")[[1]]
  if (state_i == "Maryland") neighbor[[state_i]] <- c(neighbor[[state_i]], "District of Columbia")
  if (state_i == "Virginia") neighbor[[state_i]] <- c(neighbor[[state_i]], "District of Columbia")
}
neighbor$`District of Columbia` <- c("Maryland", "Virginia")

# Cocaine Data
cocaine <- stride %>%
  filter(Drug %in% c("COCAINE", "COCAINE HYDROCHLORIDE") & Seize.Year != "NULL" & Seize.Month != "NULL" & State != "NULL") %>% 
  mutate(state=as.factor(State),
         MethAcq=as.factor(MethAcq),
         Drug=as.factor(Drug),
         Potency=ifelse(Potency > 100, Potency/10, Potency),
         Potency=ifelse(Potency > 0 & Potency < 1, Potency*100, Potency),
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

overdose <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.5 only 1999-2020.xlsx")
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

years1 <- 2000:2003
years2 <- 2004:2008
years3 <- 2009:2014
}

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

# no D.C. data
period <- years3
price_period <- cocaine %>%
  filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  group_by(state) %>% 
  summarise(avg_price=mean(adjusted_price, na.rm=T),
            med_price=median(adjusted_price, na.rm=T))

cocaine_prices <- cocaine %>% 
  filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000)

seizure_period <- cocaine %>%
  filter(!is.na(state) & Seize.Year %in% period) %>%
  group_by(state) %>% 
  summarise(max_weight=max(Nt.Wt, na.rm=T),
            max_cocaine_weight=max(Nt.Wt*Potency/100, na.rm=T))


names(overdose)[3] <- "deaths"
overdose_period <- overdose %>% 
  filter(!is.na(state) & year %in% period) %>%
  group_by(state) %>% 
  summarise(avg_death=mean(deaths, na.rm=T),
            med_death=median(deaths, na.rm=T))


states_data <- full_join(price_period, seizure_period, by="state") %>% 
  left_join(overdose_period, by="state") %>% 
  arrange(state) %>% 
  right_join(states %>% 
               rename(state=state_name) %>% 
               filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
               group_by(state) %>% 
               summarise(long=mean(long), lat=mean(lat)), by="state") %>% 
  mutate(states_index=1:length(state))

states_data$med_death <- ifelse(is.na(states_data$med_death), 0, states_data$med_death)  
states_data$max_weight <- ifelse(is.na(states_data$max_weight), 0, states_data$max_weight)


# Optimization
{
  states_data <- states_data %>% 
    filter(!is.na(med_death))
  N <- nrow(states_data) # Total number of nodes
  states_data$states_index <- 1:N
  
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_index)
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
  O <- states_data$med_death
  price <- states_data$med_price
  
  # Constraints
  # Conservation of flow/consumption
  A1 <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A1[i, which(str_match(d_vars, paste0("w", i, ","))[,1] == paste0("w", i, ","))] <- -1
    A1[i, which(str_match(d_vars, paste0(",", i, "."))[,1] == paste0(",", i, "."))] <- 1
    A1[i, n_d_vars-2*N+i] <- 1
  }
  b1 <- PO*O/sum(O)
  
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
  for (i in 1:(n_d_vars-2*N)) {
    var <- d_vars[i]
    comma_index <- which(strsplit(var, "")[[1]] == ",")
    dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
    source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
    destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
    if (is.na(price[destination_index]) | is.na(price[source_index])) next
    if (price[destination_index] < price[source_index]) less_than_price <- c(less_than_price, i)
  }
  zero_d_vars_index <- unique(less_than_price) %>% sort
  n_zero_d_vars <- length(zero_d_vars_index)
  
  # Source indicator
  A3 <- matrix(0, 2*N, n_d_vars)
  for (i in 1:N) {
    A3[i, which(d_vars %in% c(paste0("S", i), paste0("x", i)))] <- c(1/min(b1[b1>0]), -1)
    A3[i+N, which(d_vars %in% c(paste0("S", i), paste0("x", i)))] <- c(-0.01, 1)
  }
  b3 <- rep(0, 2*N)
  
  A_source_limit <- matrix(0, 1, n_d_vars)
  A_source_limit[1, grep("S", d_vars)] <- 1
  b_source_limit <- 100
}

## results with b2_ratio = 100*seizure/sum(seizure)

# results with varying alpha
# PoolSearchMode=2

alpha_results <- data.frame()
for (alpha in seq(1, 0, by=-0.1)) {
  model_alpha <- list()
  model_alpha$obj        <- c(rep(1-alpha, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha, N))
  model_alpha$modelsense <- "min"
  model_alpha$A          <- rbind(A1[,-zero_d_vars_index], A3[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_alpha$sense      <- rep(">", nrow(model_alpha$A))
  model_alpha$sense[1:N] <- "="
  model_alpha$sense[nrow(model_alpha$A)] <- "="
  model_alpha$vtype      <- c(rep("C", length(model_alpha$obj) - N), rep("B", N)) # Continuous and binary
  model_alpha$rhs        <- c(b1, b3, b_source_limit)
  result_alpha_sum       <- gurobi(model_alpha, list(Method=-1))
  objvals <- c()
  for (j in 1:length(result_alpha_sum$pool)) {
    objvals <- c(objvals, result_alpha_sum$pool[[j]]$objval)
  }
  alpha_results <- rbind(alpha_results, c(alpha, result_alpha_sum$status, sum(objvals == min(objvals))))
}
names(alpha_results) <- c("alpha", "feasibility", "n_opt_sols")
alpha_results


for (alpha in seq(0.99, 0.91, by=-0.01)) {
  model_alpha <- list()
  model_alpha$obj        <- c(rep(1-alpha, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha, N))
  model_alpha$modelsense <- "min"
  model_alpha$A          <- rbind(A1[,-zero_d_vars_index], A3[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_alpha$sense      <- rep(">", nrow(model_alpha$A))
  model_alpha$sense[1:N] <- "="
  model_alpha$sense[nrow(model_alpha$A)] <- "="
  model_alpha$vtype      <- c(rep("C", length(model_alpha$obj) - N), rep("B", N)) # Continuous and binary
  model_alpha$rhs        <- c(b1, b3, b_source_limit)
  result_alpha_sum       <- gurobi(model_alpha, list(Method=-1))
  
  opt_sol <- data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_alpha_sum$pool[[1]]$xn) %>%
    filter(optimal_sols > 0 & !grepl("x", deicision_vars))
  result_flow_index <- grep("w", opt_sol$deicision_vars)
  result_source_index <- grep("S", opt_sol$deicision_vars)
  max_flow <- max(opt_sol$optimal_sols[result_flow_index])
  
  opt_sol <- cbind(opt_sol, find_state_index(opt_sol$deicision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  opt_sol$long_i <- states_data$long[opt_sol$i]
  opt_sol$lat_i <- states_data$lat[opt_sol$i]
  opt_sol$long_j <- states_data$long[opt_sol$j]
  opt_sol$lat_j <- states_data$lat[opt_sol$j]
  
  overdose_map_year <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, med_death),
                                 by="state")
  overdose_map_year %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_death),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
    coord_quickmap() +
    labs(x="", y="", title=bquote(alpha == .(alpha))) +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(data=states_data,
               aes(x=long, y=lat)) +
    geom_segment(data=opt_sol[result_flow_index,],
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.2+2*opt_sol$optimal_sols[result_flow_index]/max_flow,
                 arrow=arrow(angle=10,
                             length=unit(0.3, "cm"),
                             type="closed")
    ) +
    geom_text(data=opt_sol[result_source_index,],
              aes(x=long_i, y=lat_i),
              label=opt_sol[result_source_index,]$optimal_sols %>% round(1),
              color="red",
              nudge_x = c(0.1, 0), nudge_y = c(-0.5, 0)) -> optimal_network_map
  ggsave(paste0("Cocaine Network Optimization/Results/optimal networks ",
                period[1], "-", period[length(period)],
                " (alpha=",
                alpha,
                ").png"),
         optimal_network_map, scale=1.5)
}


# epsilon_seizure test
alpha_i <- 0.9
epsilon_si <- 0.44
opt_val_i <- 18.96
epsilon_s_results <- data.frame()
for (epsilon_s in seq(0, 1, by=0.02)) {
  model_epsilon_s <- list()
  model_epsilon_s$obj        <- c(rep(1-alpha, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha, N))
  model_epsilon_s$modelsense <- "min"
  model_epsilon_s$A          <- rbind(A1[,-zero_d_vars_index], A2[,-zero_d_vars_index], A3[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_epsilon_s$sense      <- rep(">", nrow(model_epsilon_s$A))
  model_epsilon_s$sense[1:N] <- "="
  model_epsilon_s$sense[nrow(model_epsilon_s$A)] <- "="
  model_epsilon_s$vtype      <- c(rep("C", length(model_epsilon_s$obj) - N), rep("B", N)) # Continuous and binary
  model_epsilon_s$rhs        <- c(b1, epsilon_s*b2_ratio, b3, b_source_limit)
  result_epsilon_s_sum       <- gurobi(model_epsilon_s, list(Method=-1))
  objvals <- c()
  for (j in 1:length(result_epsilon_s_sum$pool)) {
    objvals <- c(objvals, result_epsilon_s_sum$pool[[j]]$objval)
  }
  epsilon_s_results <- rbind(epsilon_s_results, c(epsilon_s, result_epsilon_s_sum$status, sum(objvals == min(objvals))))
}
names(epsilon_s_results) <- c("epsilon_s", "feasibility", "n_opt_sols")
epsilon_s_results

for (alpha_i in seq(1, 0, by=-.1)) {
  for (epsilon_si in seq(.4, .1, by=-.1)) {
    model_epsilon_s <- list()
    model_epsilon_s$obj        <- c(rep(1-alpha_i, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha_i, N))
    model_epsilon_s$modelsense <- "min"
    model_epsilon_s$A          <- rbind(A1[,-zero_d_vars_index], A2[,-zero_d_vars_index], A3[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
    model_epsilon_s$sense      <- rep(">", nrow(model_epsilon_s$A))
    model_epsilon_s$sense[1:N] <- "="
    model_epsilon_s$sense[nrow(model_epsilon_s$A)] <- "="
    model_epsilon_s$vtype      <- c(rep("C", length(model_epsilon_s$obj) - N), rep("B", N)) # Continuous and binary
    model_epsilon_s$rhs        <- c(b1, epsilon_si*b2_ratio, b3, b_source_limit)
    result_epsilon_s_sum       <- gurobi(model_epsilon_s, list(Method=-1))
    
    opt_val_i <- round(result_epsilon_s_sum$objval, 2)
    
    opt_sol <- data.frame(decision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_epsilon_s_sum$pool[[1]]$xn) %>%
      filter(optimal_sols > 0 & !grepl("x", decision_vars))
    result_flow_index <- grep("w", opt_sol$decision_vars)
    result_source_index <- grep("S", opt_sol$decision_vars)
    max_flow <- max(opt_sol$optimal_sols[result_flow_index])
    
    opt_sol %>% write.csv(paste0("Cocaine Network Optimization/Results/optimal solution ",
                                 period[1], "-", period[length(period)],
                                 " (alpha=",
                                 alpha_i,
                                 ", epsilton_s=",
                                 epsilon_si,
                                 ", z=",
                                 opt_val_i,
                                 ").csv"),
                          row.names=F)
    
    # opt_sol <- read.csv("Cocaine Network Optimization/Results with Overdose deaths T40.5 1999-2016/optimal solution 2009-2014 (alpha=0.9, epsilton_s=0.44, z=18.96).csv") %>%
    #   rename(decision_vars = deicision_vars)
    # result_flow_index <- grep("w", opt_sol$decision_vars)
    # result_source_index <- grep("S", opt_sol$decision_vars)
    # max_flow <- max(opt_sol$optimal_sols[result_flow_index])
    opt_sol <- cbind(opt_sol, find_state_index(opt_sol$decision_vars) %>% t) %>% 
      rename(i="1", j="2")
    
    states <- get_urbn_map(map = "states", sf = TRUE)
    overdose_map_year <- left_join(states %>%
                                     rename(state=state_name) %>% 
                                     filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
                                     mutate(centroid = st_centroid(geometry)),
                                   states_data %>% select(state, med_death, states_index),
                                   by="state") %>% 
      arrange(states_index)
    
    overdose_map_year_coords <- st_coordinates(overdose_map_year$centroid)
    opt_sol$long_i <- overdose_map_year_coords[,1][opt_sol$i]
    opt_sol$lat_i <- overdose_map_year_coords[,2][opt_sol$i]
    opt_sol$long_j <- overdose_map_year_coords[,1][opt_sol$j]
    opt_sol$lat_j <- overdose_map_year_coords[,2][opt_sol$j]
    
    opt_sol_Si <- opt_sol %>% filter(grepl("S", decision_vars)) %>% rename(states_index = i)
    opt_sol_Si_sf <- left_join(overdose_map_year %>% filter(states_index %in% opt_sol_Si$states_index), opt_sol_Si %>% select(states_index, optimal_sols), by="states_index")
    
    overdose_map_year %>% ggplot() +
      geom_sf(aes(fill = med_death),
              color="black") +
      # geom_polygon(aes(x=long,
      #                  y=lat,
      #                  group=group,
      #                  fill=med_death),
      #              color="black") +
      scale_fill_viridis_c(na.value="white") +
      expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
      # coord_quickmap() +
      labs(x="", y="", fill="") +#title=bquote(paste(epsilon[s], "=", .(epsilon_si), ", ", alpha, "=", .(alpha_i), ", opt_val=", .(opt_val_i)) )) +
      theme_bw() + 
      theme(axis.ticks = element_blank(),
            axis.line =  element_blank(),
            axis.text = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_sf(data = overdose_map_year$centroid, color = "black", size = 2) +
      geom_segment(data=opt_sol[result_flow_index,],
                   aes(x=long_i, 
                       y=lat_i, 
                       xend=long_j,
                       yend=lat_j),
                   color="grey60",
                   linewidth = 0.2+2*opt_sol$optimal_sols[result_flow_index]/max_flow,
                   arrow=arrow(angle=10,
                               length=unit(0.3, "cm"),
                               type="closed")
      ) +
      geom_sf_text(data = opt_sol_Si_sf, aes(label=round(optimal_sols, 1)), color = "red", size =3) -> optimal_network_map
    ggsave(paste0("Cocaine Network Optimization/Results/optimal networks epsilton_s ",
                  period[1], "-", period[length(period)],
                  " (alpha=",
                  alpha_i,
                  ", epsilton_s=",
                  epsilon_si,
                  ", z=",
                  opt_val_i,
                  ").png"),
           optimal_network_map, scale=1)#width=20, height=10, unit="cm")
  }
}
