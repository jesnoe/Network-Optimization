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
         Potency=ifelse(Potency > 100, Potency/10, Potency),
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
  filter(state != "District of Columbia") %>% 
  filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  group_by(state) %>% 
  summarise(avg_price=mean(adjusted_price, na.rm=T),
            med_price=median(adjusted_price, na.rm=T),
            sd_price=sd(adjusted_price, na.rm=T))

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


states_data <- full_join(price_period, seizure_period, by="state") %>% 
  left_join(overdose_period, by="state") %>% 
  arrange(state) %>% 
  right_join(states %>% 
               rename(state=state_name) %>% 
               filter(!(state %in% c("Alaska", "District of Columbia", "Hawaii"))) %>% 
               group_by(state) %>% 
               summarise(long=mean(long), lat=mean(lat)), by="state") %>% 
  mutate(states_index=1:length(state))

states_data$med_death <- ifelse(is.na(states_data$med_death), 0, states_data$med_death)  
states_data$max_weight <- ifelse(is.na(states_data$max_weight), 0, states_data$max_weight)


cocaine_prices <- cocaine %>% 
  filter(state != "District of Columbia" & !is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000)

delta_5 <- 5
delta_10 <- 10
delta_15 <- 15
delta_20 <- 20
equivalence_test_summary <- tibble()
for (i in 1:48) {
  state_i <- states_data$state[i]
  bordering_states <- neighbor[[state_i]]
  
  state_prices <- cocaine_prices %>% 
    filter(Seize.Year %in% period & state == state_i) %>% 
    pull(adjusted_price)
  
  if (length(state_prices) < 3) next
  
  for(bordering_state in bordering_states) {
    bordering_state_prices <- cocaine_prices %>% 
      filter(Seize.Year %in% period & state == bordering_state) %>% 
      pull(adjusted_price)
    
    if (length(bordering_state_prices) < 3) next
    price_t_test <- t.test(state_prices, bordering_state_prices, conf.level=0.9)
    reject_ij_5 <- ifelse(-delta_5 < price_t_test$conf.int[1] & delta_5 > price_t_test$conf.int[2], 1, 0)
    reject_ij_10 <- ifelse(-delta_10 < price_t_test$conf.int[1] & delta_10 > price_t_test$conf.int[2], 1, 0)
    reject_ij_15 <- ifelse(-delta_15 < price_t_test$conf.int[1] & delta_15 > price_t_test$conf.int[2], 1, 0)
    reject_ij_20 <- ifelse(-delta_20 < price_t_test$conf.int[1] & delta_20 > price_t_test$conf.int[2], 1, 0)
    equivalence_test_summary <- rbind(equivalence_test_summary,
                                      tibble(state_index=i,
                                             bordering_state_index=states_data$states_index[which(states_data$state == bordering_state)],
                                             state=state_i,
                                             bordering_state=bordering_state,
                                             mean_state=price_t_test$estimate[1],
                                             mean_bordering_state=price_t_test$estimate[2],
                                             LCB=price_t_test$conf.int[1],
                                             UCB=price_t_test$conf.int[2],
                                             reject_delta_5=reject_ij_5,
                                             reject_delta_10=reject_ij_10,
                                             reject_delta_15=reject_ij_15,
                                             reject_delta_20=reject_ij_20,
                                             p_value=price_t_test$p.value))
  }
}
equivalence_test_summary
equivalence_test_summary %>% filter(p_value > 0.05)
equivalence_test_summary %>% filter(reject_delta_20==1)
# write.csv(equivalence_test_summary, paste0("Cocaine Network Optimization/equivalence_test_summary (", period[1], "-", period[length(period)], ").csv"), row.names=F)

equivalence_test_summary %>% select(reject_delta_5:reject_delta_20) %>% apply(2, sum)
