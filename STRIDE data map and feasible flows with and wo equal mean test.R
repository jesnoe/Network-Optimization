# setwd("/Users/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(readxl)
library(gurobi)
library(stringi)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
# may need this for projected figures
# coordinate.HIDTA <- left_join(counties_sf, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
# names(coordinate.HIDTA)[1] <- "GEOID"
{
  stride <- read.csv("STRIDE_Raw.csv") %>% as_tibble %>% filter(!is.na(Seize.Year) & !is.na(Seize.Month) & !(State %in% c("AK", "HI")))
  stride$Nt.Wt <- as.numeric(stride$Nt.Wt)
  stride$Potency <- as.numeric(stride$Potency)
  stride$Post.Price <- as.numeric(stride$Post.Price)
  # for (i in c(1,2,6,7)) {
  #   stride[,i] <- as.factor(stride[,i]) 
  # }
  stride$State <- ifelse(stride$State == "NB", "NE", stride$State)
  
  price_t_test_summary <- read.csv("Cocaine Network Optimization/price t_test_summary (2009-2014).csv") %>% as_tibble
  purity_t_test_summary <- read.csv("Cocaine Network Optimization/purity t_test_summary (2009-2014).csv") %>% as_tibble
  
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

# price/purity maps without equal mean test
periods <- list(p1=years1, p2=years2, p3=years3)
for (period in periods) {
  price_period <- cocaine %>%
    filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
    group_by(state) %>% 
    summarise(avg_price=mean(adjusted_price, na.rm=T),
              med_price=median(adjusted_price, na.rm=T))
  
  purity_period <- cocaine %>%
    filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
    group_by(state) %>% 
    summarise(avg_purity=mean(Potency, na.rm=T),
              med_purity=median(Potency, na.rm=T))
  
  states_data <- full_join(price_period, purity_period, by="state") %>% 
    right_join(states %>% 
                 rename(state=state_name) %>% 
                 filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
                 group_by(state) %>% 
                 summarise(long=mean(long), lat=mean(lat)), by="state") %>% 
    arrange(state) %>% 
    mutate(states_index=1:length(state))
  N <- nrow(states_data)
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_index)
    d_vars <- c(d_vars, paste0("w", i, ",", bordering_states_index, "."))
  }
  n_d_vars <- length(d_vars)
  
  
  price <- states_data$med_price
  purity <- states_data$med_purity
  # Directionality of flow (no equal mean test)
  less_than_price_no_test <- c()
  greater_than_purity_no_test <- c()
  for (i in 1:(n_d_vars)) {
    var <- d_vars[i]
    comma_index <- which(strsplit(var, "")[[1]] == ",")
    dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
    source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
    destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
    if (is.na(price[destination_index]) | is.na(price[source_index])) next
    if (is.na(purity[destination_index]) | is.na(purity[source_index])) next
    if (price[destination_index] < price[source_index]) less_than_price_no_test <- c(less_than_price_no_test, i)
    if (purity[destination_index] > purity[source_index]) greater_than_purity_no_test <- c(greater_than_purity_no_test, i)
  }
  
  price_restricted_d_vars_index <- unique(less_than_price_no_test) %>% sort
  d_var_price_restricted <- d_vars[-price_restricted_d_vars_index]
  purity_restricted_d_vars_index <- unique(greater_than_purity_no_test) %>% sort
  d_var_purity_restricted <- d_vars[-purity_restricted_d_vars_index]
  d_var_intersect <- intersect(d_var_price_restricted, d_var_purity_restricted)
  
  allowed_flows <- data.frame(decision_vars=d_var_price_restricted)
  allowed_flows <- cbind(allowed_flows, find_state_index(allowed_flows$decision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  allowed_flows$long_i <- states_data$long[allowed_flows$i]
  allowed_flows$lat_i <- states_data$lat[allowed_flows$i]
  allowed_flows$long_j <- states_data$long[allowed_flows$j]
  allowed_flows$lat_j <- states_data$lat[allowed_flows$j]
  
  allowed_flows_map <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, med_price, med_purity),
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
                 color="grey60",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> med_price_direction_map
  
  ggsave(paste0("Cocaine Network Optimization/Figs/price map/median price restricted flows (", period[1], "-", period[length(period)], ").png"),
         med_price_direction_map, width=15, height=8, unit="cm")
  
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
    geom_segment(data=allowed_flows %>% filter(decision_vars %in% d_var_intersect),
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="grey60",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> intersect_direction_map
  ggsave(paste0("Cocaine Network Optimization/Figs/price map/median price, purity common flows (", period[1], "-", period[length(period)], ").png"),
         intersect_direction_map, width=15, height=8, unit="cm")
  
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
    geom_segment(data=allowed_flows %>% filter(!(decision_vars %in% d_var_intersect)),
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="grey60",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> price_only_direction_map
  ggsave(paste0("Cocaine Network Optimization/Figs/price map/median price only flows (", period[1], "-", period[length(period)], ").png"),
         price_only_direction_map, width=15, height=8, unit="cm")
    
  # purity maps
  allowed_flows <- data.frame(decision_vars=d_var_purity_restricted)
  allowed_flows <- cbind(allowed_flows, find_state_index(allowed_flows$decision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  allowed_flows$long_i <- states_data$long[allowed_flows$i]
  allowed_flows$lat_i <- states_data$lat[allowed_flows$i]
  allowed_flows$long_j <- states_data$long[allowed_flows$j]
  allowed_flows$lat_j <- states_data$lat[allowed_flows$j]
  
  allowed_flows_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_purity),
                 color="black") +
    scale_fill_viridis_c(na.value="white", limits=c(0,100)) +
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
    ) -> med_purity_direction_map
  
  ggsave(paste0("Cocaine Network Optimization/Figs/purity map/median purity restricted flows (", period[1], "-", period[length(period)], ").png"),
         med_purity_direction_map, width=15, height=8, unit="cm")
  
  allowed_flows_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_purity),
                 color="black") +
    scale_fill_viridis_c(na.value="white", limits=c(0,100)) +
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
    geom_segment(data=allowed_flows %>% filter(!(decision_vars %in% d_var_intersect)),
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="red",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> purity_only_direction_map
  ggsave(paste0("Cocaine Network Optimization/Figs/purity map/median purity only flows (", period[1], "-", period[length(period)], ").png"),
         purity_only_direction_map, width=15, height=8, unit="cm")
}
  # for period 2009-2014
d_var_price_restricted %>% length     # 115
d_var_purity_restricted %>% length # 115
intersect(d_var_price_restricted, d_var_purity_restricted) %>% length # 71

d_var_price_restricted[!(d_var_price_restricted %in% d_var_purity_restricted)] %>% length # 44
d_var_purity_restricted[!(d_var_purity_restricted %in% d_var_price_restricted)] %>% length # 44

price_opposite_index <- d_var_price_restricted[!(d_var_price_restricted %in% d_var_purity_restricted)]
price_allowed_flows <- data.frame(decision_vars=d_var_price_restricted)
price_allowed_flows <- cbind(price_allowed_flows, find_state_index(price_allowed_flows$decision_vars) %>% t) %>% 
  rename(i="1", j="2") %>% as_tibble
price_opposite_flows <- price_allowed_flows %>% filter(decision_vars %in% price_opposite_index)

purity_opposite_index <- d_var_purity_restricted[!(d_var_purity_restricted %in% d_var_price_restricted)]
purity_allowed_flows <- data.frame(decision_vars=d_var_purity_restricted)
purity_allowed_flows <- cbind(purity_allowed_flows, find_state_index(purity_allowed_flows$decision_vars) %>% t) %>% 
  rename(i="1", j="2") %>% as_tibble
purity_opposite_flows <- purity_allowed_flows %>% filter(decision_vars %in% purity_opposite_index)

price_purity_t_test_summary <- price_t_test_summary %>% 
  rename(price_p_value=p_value) %>% 
  mutate(purity_p_value=purity_t_test_summary$p_value,
         pair=paste(state_index, bordering_state_index))

opposite_dir_pairs <- price_purity_t_test_summary %>%
  filter((pair %in% paste(price_opposite_flows$i, price_opposite_flows$j)) | (pair %in% paste(price_opposite_flows$j, price_opposite_flows$i)))

insignificant_pairs <- tibble()
for(alpha_t_test in c(0.05, 0.1, 0.25, 0.5)) {
  price_insignificant <- opposite_dir_pairs %>% filter(price_p_value > alpha_t_test) %>% nrow
  purity_insignificant <- opposite_dir_pairs %>% filter(purity_p_value > alpha_t_test) %>% nrow
  both_insignificant <- opposite_dir_pairs %>% filter(price_p_value > alpha_t_test & purity_p_value > alpha_t_test) %>% nrow
  insignificant_pairs <- rbind(insignificant_pairs,
                               c(alpha_t_test, price_insignificant, purity_insignificant, both_insignificant))
}
names(insignificant_pairs) <- c("p_value", "n_price_relaxed", "n_purity_relaxed", "n_both_relaxed")
insignificant_pairs

same_dir_pairs <- price_purity_t_test_summary %>% filter(!(pair %in% opposite_dir_pairs$pair))
insignificant_pairs_same <- tibble()
for(alpha_t_test in c(0.05, 0.1, 0.25, 0.5)) {
  price_significant <- same_dir_pairs %>% filter(price_p_value > alpha_t_test) %>% nrow
  purity_significant <- same_dir_pairs %>% filter(purity_p_value > alpha_t_test) %>% nrow
  both_significant <- same_dir_pairs %>% filter(price_p_value > alpha_t_test & purity_p_value > alpha_t_test) %>% nrow
  insignificant_pairs_same <- rbind(insignificant_pairs_same,
                               c(alpha_t_test, price_significant, purity_significant, both_significant))
}
names(insignificant_pairs_same) <- c("p_value", "n_price_relaxed", "n_purity_relaxed", "n_both_relaxed")
insignificant_pairs_same

price_purity_t_test_summary %>% 
  filter((pair %in% paste(purity_opposite_flows$i, purity_opposite_flows$j)) | (pair %in% paste(purity_opposite_flows$j, purity_opposite_flows$i)))

price_purity_t_test_summary %>% filter(price_p_value > 0.05 & purity_p_value > 0.05) # 47
price_purity_t_test_summary %>% filter(price_p_value > 0.05 & purity_p_value <= 0.05) # 17
price_purity_t_test_summary %>% filter(price_p_value <= 0.05 & purity_p_value > 0.05) # 25


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
  geom_segment(data=allowed_flows %>% filter(decision_vars %in% d_var_price_restricted[!(d_var_price_restricted %in% d_var_purity_restricted)]),
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
# ggsave(paste0("Cocaine Network Optimization/Figs/price map/median price specific flows (", period[1], "-", period[length(period)], ").png"),
#        med_price_direction_map, width=15, height=8, unit="cm")

allowed_flows_map %>% ggplot() +
  geom_polygon(aes(x=long,
                   y=lat,
                   group=group,
                   fill=med_purity),
               color="black") +
  scale_fill_viridis_c(na.value="white", limits=c(0,100)) +
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
  geom_segment(data=allowed_flows %>% filter(decision_vars %in% d_var_purity_restricted[!(d_var_purity_restricted %in% d_var_price_restricted)]),
               aes(x=long_i, 
                   y=lat_i, 
                   xend=long_j,
                   yend=lat_j),
               linewidth = 0.3,
               color="red",
               arrow=arrow(angle=10,
                           length=unit(0.2, "cm"),
                           type="closed")
  ) -> med_purity_direction_map
# ggsave(paste0("Cocaine Network Optimization/Figs/purity map/median purity specific flows (", period[1], "-", period[length(period)], ").png"),
#        med_purity_direction_map, width=15, height=8, unit="cm")

price_t_test_summary %>% filter(mean_state < mean_bordering_state) # 68
price_retricted_pairs <- price_t_test_summary %>% 
  filter(mean_state < mean_bordering_state) %>% 
  select(state, bordering_state)
  
purity_t_test_summary %>% filter(mean_state > mean_bordering_state) # 70
purity_retricted_pairs <- purity_t_test_summary %>%
  filter(mean_state > mean_bordering_state) %>% 
  select(state, bordering_state)

intersect(price_retricted_pairs, purity_retricted_pairs) # 52


# price/purity maps with equal mean test
alpha_t_test <- 0.05
periods <- list(p1=years1, p2=years2, p3=years3)
for (period in periods) {
  price_period <- cocaine %>%
    filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
    group_by(state) %>% 
    summarise(avg_price=mean(adjusted_price, na.rm=T),
              med_price=median(adjusted_price, na.rm=T))
  
  purity_period <- cocaine %>%
    filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
    group_by(state) %>% 
    summarise(avg_purity=mean(Potency, na.rm=T),
              med_purity=median(Potency, na.rm=T))
  
  states_data <- full_join(price_period, purity_period, by="state") %>% 
    arrange(state) %>% 
    right_join(states %>% 
                 rename(state=state_name) %>% 
                 filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
                 group_by(state) %>% 
                 summarise(long=mean(long), lat=mean(lat)), by="state") %>% 
    mutate(states_index=1:length(state))
  N <- nrow(states_data)
  
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_index)
    d_vars <- c(d_vars, paste0("w", i, ",", bordering_states_index, "."))
  }
  n_d_vars <- length(d_vars)
  
  
  price <- states_data$med_price
  purity <- states_data$med_purity
  
  # Directionality of flow (equal mean test)
  relaxed_pairs_index <- tibble()
  price_relaxed_pairs_index <- tibble()
  purity_relaxed_pairs_index <- tibble()
  less_than_price <- c()
  greater_than_purity <- c()
  for (i in 1:(n_d_vars)) {
    var <- d_vars[i]
    comma_index <- which(strsplit(var, "")[[1]] == ",")
    dot_index <- which(strsplit(var, "")[[1]] == ".") - 1
    source_index <- substr(var, 2, comma_index - 1) %>% as.numeric
    destination_index <- substr(var, comma_index + 1, dot_index) %>% as.numeric
    if (is.na(price[destination_index]) | is.na(price[source_index])) {
      relaxed_pairs_index <- rbind(relaxed_pairs_index, c(source_index, destination_index))
      next
    }
    if (destination_index %in% (price_t_test_summary %>% filter(state_index == source_index) %>% pull(bordering_state_index))) {
      if (price_t_test_summary %>%
          filter(state_index==source_index & bordering_state_index == destination_index) %>%
          pull(p_value) > alpha_t_test) {
        price_relaxed_pairs_index <- rbind(price_relaxed_pairs_index, c(source_index, destination_index))
      } else {
        if (price[destination_index] < price[source_index]) less_than_price <- c(less_than_price, i)
      }
    } else {
      if (price[destination_index] < price[source_index]) less_than_price <- c(less_than_price, i)
    }
    
    if (destination_index %in% (purity_t_test_summary %>% filter(state_index == source_index) %>% pull(bordering_state_index))) {
      if (purity_t_test_summary %>%
          filter(state_index==source_index & bordering_state_index == destination_index) %>%
          pull(p_value) > alpha_t_test) {
        purity_relaxed_pairs_index <- rbind(purity_relaxed_pairs_index, c(source_index, destination_index))
      } else {
        if (purity[destination_index] > purity[source_index]) greater_than_purity <- c(greater_than_purity, i)
      }
    } else {
      if (purity[destination_index] > purity[source_index]) greater_than_purity <- c(greater_than_purity, i)
    }
  }
  names(relaxed_pairs_index) <- c("source_index", "destination_index")
  relaxed_pairs_index$var_name <- paste0("w", relaxed_pairs_index$source_index, ",", relaxed_pairs_index$destination_index, ".")
  relaxed_pairs_index$var_index <- apply(relaxed_pairs_index, 1, function(x) ifelse(x[3] %in% d_vars, which(x[3] == d_vars), NA))
  
  price_relaxed_pairs_index$var_name <- paste0("w", price_relaxed_pairs_index[,1], ",", price_relaxed_pairs_index[,2], ".")
  price_relaxed_pairs_index$var_index <- apply(price_relaxed_pairs_index, 1, function(x) ifelse(x[3] %in% d_vars, which(x[3] == d_vars), NA))
  price_relaxed_pairs_index_vec <- c(price_relaxed_pairs_index$var_index,
                                     which(d_vars %in% paste0("w", price_relaxed_pairs_index[,2], ",", price_relaxed_pairs_index[,1], ".")))
  relaxed_less_than_price <- less_than_price[which(!(less_than_price %in% price_relaxed_pairs_index_vec))]
  
  purity_relaxed_pairs_index$var_name <- paste0("w", purity_relaxed_pairs_index[,1], ",", purity_relaxed_pairs_index[,2], ".")
  purity_relaxed_pairs_index$var_index <- apply(purity_relaxed_pairs_index, 1, function(x) ifelse(x[3] %in% d_vars, which(x[3] == d_vars), NA))
  purity_relaxed_pairs_index_vec <- c(purity_relaxed_pairs_index$var_index,
                                      which(d_vars %in% paste0("w", purity_relaxed_pairs_index[,2], ",", purity_relaxed_pairs_index[,1], ".")))
  relaxed_greater_than_purity <- greater_than_purity[which(!(greater_than_purity %in% purity_relaxed_pairs_index_vec))]
  
  relaxed_pairs_index_vec <- intersect(price_relaxed_pairs_index_vec, purity_relaxed_pairs_index_vec) %>% unique
  relaxed_pairs_index_vec <- c(relaxed_pairs_index_vec, relaxed_pairs_index$var_index) %>% unique %>% sort
  
  relaxed_zero_d_vars_index <- unique(c(less_than_price, greater_than_purity)) %>% sort
  relaxed_zero_d_vars_index <- relaxed_zero_d_vars_index[which(!(relaxed_zero_d_vars_index %in% relaxed_pairs_index_vec))]
  
  d_var_price_restricted <- d_vars[-relaxed_less_than_price]
  d_var_purity_restricted <- d_vars[-relaxed_greater_than_purity]
  
  allowed_flows <- data.frame(decision_vars=d_var_price_restricted)
  allowed_flows <- cbind(allowed_flows, find_state_index(allowed_flows$decision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  allowed_flows$long_i <- states_data$long[allowed_flows$i]
  allowed_flows$lat_i <- states_data$lat[allowed_flows$i]
  allowed_flows$long_j <- states_data$long[allowed_flows$j]
  allowed_flows$lat_j <- states_data$lat[allowed_flows$j]
  
  allowed_flows_map <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, med_price, med_purity),
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
                 color="grey60",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             type="closed")
    ) -> med_price_direction_map
  
  ggsave(paste0("Cocaine Network Optimization/Figs/price map/relaxed median price restricted flows (", period[1], "-", period[length(period)], ").png"),
         med_price_direction_map, width=15, height=8, unit="cm")
  
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
    geom_segment(data=allowed_flows %>% filter(decision_vars %in% price_relaxed_pairs_index$var_name),
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="grey60",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             ends="both",
                             type="closed")
    ) -> med_price_relaxed_map
  
  ggsave(paste0("Cocaine Network Optimization/Figs/price map/relaxed price directions flows (", period[1], "-", period[length(period)], ").png"),
         med_price_relaxed_map, width=15, height=8, unit="cm")
  
  allowed_flows <- data.frame(decision_vars=d_var_purity_restricted)
  allowed_flows <- cbind(allowed_flows, find_state_index(allowed_flows$decision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  allowed_flows$long_i <- states_data$long[allowed_flows$i]
  allowed_flows$lat_i <- states_data$lat[allowed_flows$i]
  allowed_flows$long_j <- states_data$long[allowed_flows$j]
  allowed_flows$lat_j <- states_data$lat[allowed_flows$j]
  
  allowed_flows_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_purity),
                 color="black") +
    scale_fill_viridis_c(na.value="white", limits=c(0,100)) +
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
    ) -> med_purity_direction_map
  
  ggsave(paste0("Cocaine Network Optimization/Figs/purity map/relaxed median purity restricted flows (", period[1], "-", period[length(period)], ").png"),
         med_purity_direction_map, width=15, height=8, unit="cm")
  
  allowed_flows_map %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_purity),
                 color="black") +
    scale_fill_viridis_c(na.value="white", limits=c(0,100)) +
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
    geom_segment(data=allowed_flows %>% filter(decision_vars %in% purity_relaxed_pairs_index$var_name),
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.3,
                 color="red",
                 arrow=arrow(angle=10,
                             length=unit(0.2, "cm"),
                             ends="both",
                             type="closed")
    ) -> med_purity_relaxed_map
  
  ggsave(paste0("Cocaine Network Optimization/Figs/purity map/relaxed purity directions flows (", period[1], "-", period[length(period)], ").png"),
         med_purity_relaxed_map, width=15, height=8, unit="cm")
}
  # for 2009-2014
relaxed_pairs_index # 3
length(price_relaxed_pairs_index_vec) #  48
length(purity_relaxed_pairs_index_vec) # 70
intersect(price_relaxed_pairs_index_vec, purity_relaxed_pairs_index_vec) %>% length # 36
d_vars[relaxed_pairs_index_vec]
states_data

# for 2004-2014
relaxed_pairs_index # 3
nrow(price_relaxed_pairs_index) #  17
nrow(purity_relaxed_pairs_index) # 29
intersect(price_relaxed_pairs_index_vec, purity_relaxed_pairs_index_vec) %>% length # 22
d_vars[relaxed_pairs_index_vec]
states_data
