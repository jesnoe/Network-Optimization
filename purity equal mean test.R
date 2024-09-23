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
  stride$State <- ifelse(stride$State == "NB", "NE", stride$State)
  
  price_t_test_summary <- read.csv("Cocaine Network Optimization/price t_test_summary (2009-2014).csv") %>% as_tibble
  
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

cocaine_purities <- cocaine %>% 
  filter(!is.na(state) & Seize.Year %in% c(years2, years3) & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000)
cocaine_purities %>% group_by(state) %>% summarize(n_pu=n()) %>% view
cocaine_purities %>% filter(state %in% c("Vermont", "South Dakota", "Oregon"))

purity_t_test_summary <- tibble()
for (i in 1:49) {
  state_i <- states_data$state[i]
  bordering_states <- neighbor[[state_i]]
  
  state_purities <- cocaine_purities %>% 
    filter(state == state_i) %>% 
    pull(Potency)
  
  if (length(state_purities) < 3) next
  
  for(bordering_state in bordering_states) {
    bordering_state_purities <- cocaine_purities %>% 
      filter(state == bordering_state) %>% 
      pull(Potency)
    
    if (length(bordering_state_purities) < 3) next
    price_t_test <- t.test(state_purities, bordering_state_purities)
    purity_t_test_summary <- rbind(purity_t_test_summary,
                            tibble(state_index=i,
                                   bordering_state_index=states_data$states_index[which(states_data$state == bordering_state)],
                                   state=state_i,
                                   bordering_state=bordering_state,
                                   mean_state=price_t_test$estimate[1],
                                   mean_bordering_state=price_t_test$estimate[2],
                                   p_value=price_t_test$p.value))
  }
}
purity_t_test_summary <- purity_t_test_summary %>% filter(state_index < bordering_state_index)
# write.csv(purity_t_test_summary, paste0("Cocaine Network Optimization/purity t_test_summary (", period[1], "-", period[length(period)], ").csv"), row.names=F)

purity_t_test_summary %>% select(state_index, bordering_state_index) %>% unique