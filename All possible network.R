# setwd("C:/Users/User/Documents/R")
library(readxl)
library(slam)
library(gurobi) # need a license. Refer to https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation_guide.html
library(stringi)
library(urbnmapr) # install by running -> devtools::install_github("UrbanInstitute/urbnmapr")
library(tidyverse)
library(gridExtra)
library(lubridate)
library(sf)

# Put "Cocaine Network Optimization" folder in your working directory
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
  
  
  overdose <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.5 1999-2020.xlsx")
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
  t_test_summary <- read.csv("Cocaine Network Optimization/price t_test_summary (2009-2014).csv")
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

# no D.C. data (? data has D.C. data)
period <- years3
price_period <- cocaine %>%
  filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  group_by(state) %>% 
  summarise(avg_price=mean(adjusted_price, na.rm=T),
            med_price=median(adjusted_price, na.rm=T))

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
            med_death=median(deaths, na.rm=T),
            sd_death=sd(deaths, na.rm=T))


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

states_data <- states_data %>% 
  filter(!is.na(med_death))
N <- nrow(states_data) # Total number of nodes
states_data$states_index <- 1:N
{
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_index)
    d_vars <- c(d_vars, paste0("w", i, ",", bordering_states_index, "."))
  }
  
  d_vars <- c(d_vars, paste0("S", 1:N))
  d_vars <- c(d_vars, paste0("x", 1:N))
  n_d_vars <- length(d_vars)
}

d_vars <- d_vars[grepl("w", d_vars)]
d_vars <- find_state_index(d_vars) %>% t %>% as_tibble %>% filter(V1 <= V2)
names(d_vars) <- c("i", "j")

states_sf <- get_urbn_map(map = "states", sf = TRUE)
overdose_map_year <- left_join(states_sf %>%
                                 rename(state=state_name) %>% 
                                 filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
                                 mutate(centroid = st_centroid(geometry)),
                               states_data %>% select(state, states_index),
                               by="state") %>% 
  arrange(states_index)

overdose_map_year_coords <- st_coordinates(overdose_map_year$centroid)
d_vars$long_i <- overdose_map_year_coords[,1][d_vars$i]
d_vars$lat_i <- overdose_map_year_coords[,2][d_vars$i]
d_vars$long_j <- overdose_map_year_coords[,1][d_vars$j]
d_vars$lat_j <- overdose_map_year_coords[,2][d_vars$j]

overdose_map_year %>% ggplot() +
  geom_sf(color="black", fill="white") +
  expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
  labs(x="", y="", fill="") +
  theme_bw() + 
  theme(axis.ticks = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_sf(data = overdose_map_year$centroid, color = "black", size = 2) +
  geom_segment(data=d_vars,
               aes(x=long_i, 
                   y=lat_i, 
                   xend=long_j,
                   yend=lat_j),
               color="black",
               linewidth = 0.5,
               arrow=arrow(angle=12,
                           length=unit(0.3, "cm"),
                           ends="both",
                           type="closed")
  ) -> optimal_network_map
ggsave("Cocaine Network Optimization/Figs/all possible flows.png", optimal_network_map, scale=1)

set.seed(100)
overdose_map_year %>% ggplot() +
  geom_sf(color="black", fill="white") +
  expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
  labs(x="", y="", fill="") +
  theme_bw() + 
  theme(axis.ticks = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_sf(data = overdose_map_year$centroid, color = "black", size = 2) +
  geom_segment(data=d_vars[sample(1:nrow(d_vars), 50),],
               aes(x=long_i, 
                   y=lat_i, 
                   xend=long_j,
                   yend=lat_j),
               color="black",
               linewidth = 0.5,
               arrow=arrow(angle=12,
                           length=unit(0.3, "cm"),
                           type="closed")
  ) -> optimal_network_map
optimal_network_map
ggsave("Cocaine Network Optimization/Figs/random subflows seed=100.png", optimal_network_map, scale=1)

set.seed(98721)
overdose_map_year %>% ggplot() +
  geom_sf(color="black", fill="white") +
  expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
  labs(x="", y="", fill="") +
  theme_bw() + 
  theme(axis.ticks = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_sf(data = overdose_map_year$centroid, color = "black", size = 2) +
  geom_segment(data=d_vars[sample(1:nrow(d_vars), 50),],
               aes(x=long_i, 
                   y=lat_i, 
                   xend=long_j,
                   yend=lat_j),
               color="black",
               linewidth = 0.5,
               arrow=arrow(angle=12,
                           length=unit(0.3, "cm"),
                           type="closed")
  ) -> optimal_network_map
optimal_network_map
ggsave("Cocaine Network Optimization/Figs/random subflows seed=98721.png", optimal_network_map, scale=1)