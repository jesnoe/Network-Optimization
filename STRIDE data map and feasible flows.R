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

years1 <- 2000:2003
years2 <- 2004:2008
years3 <- 2009:2014
}

cocaine %>%
  filter(Seize.Year %in% years1 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  select(-MethAcq, -Drug) %>% 
  relocate(state, Seize.Year, Seize.Month) %>% 
  group_by(state) %>% 
  summarise(n_price=sum(!is.na(adjusted_price)),
            price_sd_max_wt1000=sd(adjusted_price, na.rm=T)) %>% 
  print(n=49)

cocaine %>%
  filter(Seize.Year %in% years2 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  select(-MethAcq, -Drug) %>% 
  relocate(state, Seize.Year, Seize.Month) %>% 
  group_by(state) %>% 
  summarise(n_price=sum(!is.na(adjusted_price)),
            price_sd_max_wt1000=sd(adjusted_price, na.rm=T)) %>% 
  print(n=49)

cocaine %>%
  filter(Seize.Year %in% years3 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  select(-MethAcq, -Drug) %>% 
  relocate(state, Seize.Year, Seize.Month) %>% 
  group_by(state) %>% 
  summarise(n_price=sum(!is.na(adjusted_price)),
            price_sd_max_wt1000=sd(adjusted_price, na.rm=T)) %>% 
  print(n=49)

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
periods <- list(p1=years1, p2=years2, p3=years3)
for (period in periods) {
  price_period <- cocaine %>%
    filter(state != "District of Columbia") %>% 
    filter(!is.na(state) & Seize.Year %in% period & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
    group_by(state) %>% 
    summarise(avg_price=mean(adjusted_price, na.rm=T),
              med_price=median(adjusted_price, na.rm=T))
  
  seizure_period <- cocaine %>%
    filter(state != "District of Columbia") %>% 
    filter(!is.na(state) & Seize.Year %in% period) %>%
    group_by(state) %>% 
    summarise(max_weight=max(Nt.Wt, na.rm=T),
              max_cocaine_weight=max(Nt.Wt*Potency/100, na.rm=T))
  
  overdose_period <- overdose %>% 
    filter(!is.na(state) & year %in% period) %>%
    group_by(state) %>% 
    summarise(avg_death=mean(deaths, na.rm=T),
              med_death=median(deaths, na.rm=T))
  
  N <- nrow(price_period)
  states_data <- full_join(price_period, seizure_period, by="state") %>% 
    left_join(overdose_period, by="state") %>% 
    arrange(state) %>% 
    right_join(states %>% 
                 rename(state=state_name) %>% 
                 filter(!(state %in% c("Alaska", "District of Columbia", "Hawaii"))) %>% 
                 group_by(state) %>% 
                 summarise(long=mean(long), lat=mean(lat)), by="state") %>% 
    mutate(states_index=1:length(state))
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_index)
    d_vars <- c(d_vars, paste0("w", i, ",", bordering_states_index, "."))
  }
  n_d_vars <- length(d_vars)
  price <- states_data$med_price
  less_than_price <- c()
  for (i in 1:(n_d_vars)) {
    var <- d_vars[i]
    comma_index <- which(strsplit(var, "")[[1]] == ",")
    dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
    source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
    destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
    if (is.na(price[destination_index]) | is.na(price[source_index])) next
    if (price[destination_index] < price[source_index]) less_than_price <- c(less_than_price, i)
  }
  restricted_d_vars_index <- unique(less_than_price) %>% sort
  d_var_restricted <- d_vars[-restricted_d_vars_index]
  
  allowed_flows <- data.frame(deicision_vars=d_var_restricted)
  allowed_flows <- cbind(allowed_flows, find_state_index(allowed_flows$deicision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  allowed_flows$long_i <- states_data$long[allowed_flows$i]
  allowed_flows$lat_i <- states_data$lat[allowed_flows$i]
  allowed_flows$long_j <- states_data$long[allowed_flows$j]
  allowed_flows$lat_j <- states_data$lat[allowed_flows$j]
  
  allowed_flows_map <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, med_price),
                                 by="state")
  allowed_flows_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_price),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=allowed_flows_map$long, y=allowed_flows_map$lat) +
    coord_quickmap() +
    labs(x="", y="") +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(data=states_data,
               aes(x=long, y=lat)) +
    geom_segment(data=allowed_flows,
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="red",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> med_price_direction_map
  
  price <- states_data$avg_price
  less_than_price <- c()
  for (i in 1:(n_d_vars)) {
    var <- d_vars[i]
    comma_index <- which(strsplit(var, "")[[1]] == ",")
    dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
    source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
    destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
    if (is.na(price[destination_index]) | is.na(price[source_index])) next
    if (price[destination_index] < price[source_index]) less_than_price <- c(less_than_price, i)
  }
  restricted_d_vars_index <- unique(less_than_price) %>% sort
  d_var_restricted <- d_vars[-restricted_d_vars_index]
  
  allowed_flows <- data.frame(deicision_vars=d_var_restricted)
  allowed_flows <- cbind(allowed_flows, find_state_index(allowed_flows$deicision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  allowed_flows$long_i <- states_data$long[allowed_flows$i]
  allowed_flows$lat_i <- states_data$lat[allowed_flows$i]
  allowed_flows$long_j <- states_data$long[allowed_flows$j]
  allowed_flows$lat_j <- states_data$lat[allowed_flows$j]
  
  allowed_flows_map <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, avg_price),
                                 by="state")
  allowed_flows_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=avg_price),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=allowed_flows_map$long, y=allowed_flows_map$lat) +
    coord_quickmap() +
    labs(x="", y="") +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(data=states_data,
               aes(x=long, y=lat)) +
    geom_segment(data=allowed_flows,
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="red",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> avg_price_direction_map
  
  max_seizure_map <- left_join(states %>%
                                 rename(state=state_name) %>% 
                                 filter(!(state %in% c("Alaska", "Hawaii"))),
                               states_data %>% select(state, max_weight),
                               by="state")
  max_seizure_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=log(max_weight)),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=allowed_flows_map$long, y=allowed_flows_map$lat) +
    coord_quickmap() +
    labs(x="", y="", fill="Log max seizures") +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> max_seizure_map_period
  
  med_overdose_map <- left_join(states %>%
                                 rename(state=state_name) %>% 
                                 filter(!(state %in% c("Alaska", "Hawaii"))),
                               states_data %>% select(state, med_death),
                               by="state")
  med_overdose_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_death),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=allowed_flows_map$long, y=allowed_flows_map$lat) +
    coord_quickmap() +
    labs(x="", y="", fill="Median deaths") +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> med_overdose_map_period
  
  ggsave(paste0("Cocaine Network Optimization/Figs/median price restricted flows (", period[1], "-", period[length(period)], ").png"),
         med_price_direction_map, scale=1.5)
  ggsave(paste0("Cocaine Network Optimization/Figs/average price restricted flows (", period[1], "-", period[length(period)], ").png"),
         avg_price_direction_map, scale=1.5)
  ggsave(paste0("Cocaine Network Optimization/Figs/log max seizure map (", period[1], "-", period[length(period)], ").png"),
         max_seizure_map_period, scale=1.5)
  ggsave(paste0("Cocaine Network Optimization/Figs/median death map (", period[1], "-", period[length(period)], ").png"),
         med_overdose_map_period, scale=1.5)
}
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
  states_data <- tibble(state=states_data$state, states_data=1:N) %>% 
    arrange(state) %>% 
    left_join(states %>% 
                rename(state=state_name) %>% 
                group_by(state) %>% 
                summarise(long=mean(long), lat=mean(lat)), by="state")
  
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_data)
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

