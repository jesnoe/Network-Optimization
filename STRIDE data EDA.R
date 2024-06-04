# setwd("/Users/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
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

states_data$med_death <- ifelse(is.na(states_data$med_death), 0, states_data$med_death)  
states_data$max_weight <- ifelse(is.na(states_data$max_weight), 0, states_data$max_weight)

# monthly price
cocaine_monthly <- cocaine %>%
  filter(state != "District of Columbia") %>% 
  filter(!is.na(state) & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  group_by(state, Seize.Year, Seize.Month) %>% 
  summarise(n_price=sum(!is.na(adjusted_price)),
            med_price=median(adjusted_price, na.rm=T),
            avg_price=mean(adjusted_price, na.rm=T),
            sd_price=sd(adjusted_price, na.rm=T)
            ) %>% 
  relocate(Seize.Year, Seize.Month) %>% 
  arrange(Seize.Year, Seize.Month)

# cocaine_monthly %>% filter(state == "Arizona") %>% write.csv("Cocaine Network Optimization/monthly price (Arizona).csv", row.names=F)
# cocaine_monthly %>% filter(state == "California") %>% write.csv("Cocaine Network Optimization/monthly price (California).csv", row.names=F)

FL_months <- (cocaine_monthly %>% filter(state == "Florida") %>% select(Seize.Year, Seize.Month))[,-1]
GA_months <- (cocaine_monthly %>% filter(state == "Georgia") %>% select(Seize.Year, Seize.Month))[,-1]
FL_GA_price <- cocaine_monthly %>%
  filter(state %in% c("Florida", "Georgia") & Seize.Year %in% period) %>% 
  pivot_wider(names_from = "state", values_from = c("avg_price", "med_price"))

sum(FL_GA_price$med_price_Florida > FL_GA_price$med_price_Georgia, na.rm=T) # 21
sum(FL_GA_price$med_price_Florida < FL_GA_price$med_price_Georgia, na.rm=T) # 31
sum(FL_GA_price$avg_price_Florida > FL_GA_price$avg_price_Georgia, na.rm=T) # 29
sum(FL_GA_price$avg_price_Florida < FL_GA_price$avg_price_Georgia, na.rm=T) # 24


AZ_months <- (cocaine_monthly %>% filter(state == "Arizona") %>% select(Seize.Year, Seize.Month))[,-1]
CA_months <- (cocaine_monthly %>% filter(state == "California") %>% select(Seize.Year, Seize.Month))[,-1]
AZ_CA_price <- cocaine_monthly %>%
  select(Seize.Year, Seize.Month, state, med_price, avg_price) %>% 
  filter(state %in% c("Arizona", "California") & Seize.Year %in% period) %>% 
  pivot_wider(names_from = "state", values_from = c("avg_price", "med_price"))

cocaine %>% filter(state == "Arizona" & !is.na(adjusted_price)) %>% pull(adjusted_price) %>% sort
cocaine %>% filter(state == "California" & !is.na(adjusted_price)) %>% pull(adjusted_price) %>% sort

sum(AZ_CA_price$med_price_Arizona > AZ_CA_price$med_price_California, na.rm=T) # 1
sum(AZ_CA_price$med_price_Arizona < AZ_CA_price$med_price_California, na.rm=T) # 8
sum(AZ_CA_price$avg_price_Arizona > AZ_CA_price$avg_price_California, na.rm=T) # 1
sum(AZ_CA_price$avg_price_Arizona < AZ_CA_price$avg_price_California, na.rm=T) # 8

# price/purity vs. seizure weight
cocaine %>% 
  filter(Seize.Year %in% years3 & state == "Florida", !is.na(adjusted_price)) %>% 
  mutate(Seize.Year=as.factor(Seize.Year)) %>% 
  ggplot(aes(x=Nt.Wt,
             y=adjusted_price,
             color=Seize.Year)) +
  geom_point() + 
  lims(x=c(0, 300), y=c(0, 1000)) +
  labs(title="Florida") -> FL_price_plot
ggsave("Cocaine Network Optimization/Figs/price vs seizrue Florida.png", FL_price_plot, width=30, height=20, unit="cm")

cocaine %>% 
  filter(Seize.Year %in% years3 & state == "Georgia", !is.na(adjusted_price)) %>% 
  mutate(Seize.Year=as.factor(Seize.Year)) %>% 
  ggplot(aes(x=Nt.Wt,
             y=adjusted_price,
             color=Seize.Year)) +
  geom_point() + 
  lims(x=c(0, 300), y=c(0, 1000)) +
  labs(title="Georgia") -> GA_price_plot

ggsave("Cocaine Network Optimization/Figs/price vs seizrue Georgia.png", GA_price_plot, width=30, height=20, unit="cm")

 
