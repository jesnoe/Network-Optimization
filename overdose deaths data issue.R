# setwd("C:/Users/gkfrj/Documents/R")
library(readxl)
library(gurobi)
library(stringi)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
overdose_both <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.4 or T40.5 1999-2016.xlsx") %>% select(-population)
overdose_and <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.4 and T40.5 1999-2016.xlsx")
overdose_T40.5 <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.5 1999-2020.xlsx")
overdose_T40.4 <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.4 1999-2020.xlsx")
overdose_and_sup <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.4 and T40.5 suppressed 1999-2016.xlsx")

overdose_matched <- overdose_both %>% 
  full_join(overdose_T40.4, by=c("state", "year")) %>% 
  full_join(overdose_T40.5, by=c("state", "year")) %>%
  full_join(overdose_and, by=c("state", "year")) %>% 
  rename(deaths_or=deaths) %>% 
  mutate(deaths_T40.5_and=ifelse(is.na(deaths_T40.5), 0, deaths_T40.5)+ifelse(is.na(deaths_and), 0, deaths_and),
         deaths_union=ifelse(is.na(deaths_T40.5), 0, deaths_T40.5) +
           ifelse(is.na(deaths_T40.4), 0, deaths_T40.4) -
           ifelse(is.na(deaths_and), 0, deaths_and)) %>% 
  filter(!(state %in% c("Alaska", "Hawaii")))

overdose_matched
overdose_matched[which(!is.na(overdose_matched[,-(1:2)] %>% apply(1, sum))),] %>% 
  filter(deaths_or == deaths_union)
# write.csv(overdose_matched, "Cocaine Network Optimization/Overdose deaths matched 1999-2020.csv", row.names=F)
overdose_matched <- read_xlsx("Cocaine Network Optimization/Overdose deaths matched 1999-2020.xlsx")

# Suppressed = NA
sum(overdose_and_sup %>% full_join(overdose_and, by=c("state", "year")) %>% pull(deaths_and_suppressed) == "Suppressed")
overdose_and_sup %>% full_join(overdose_and, by=c("state", "year")) %>% pull(deaths_and) %>% is.na %>% sum
overdose_and_sup %>% full_join(overdose_and, by=c("state", "year")) %>% view

overdose_matched %>% arrange(deaths_T40.4, deaths_T40.5, deaths_T40.5_and) %>% view
overdose_matched %>% apply(1, function(x) sum(is.na(x[4:6]))) %>% table
overdose_matched[overdose_matched %>% apply(1, function(x) sum(is.na(x[4:6]))) > 1,] %>% view

overdose_matched$deaths_T40.5/overdose_matched$deaths_or

# Suppressed fraction
overdose_matched_ratio <- overdose_matched %>%
  mutate(deaths_or = as.numeric(deaths_or),
         deaths_T40.4 = as.numeric(deaths_T40.4),
         deaths_T40.5 = as.numeric(deaths_T40.5),
         deaths_and = as.numeric(deaths_and),
         deaths_T40.5_and = as.numeric(deaths_T40.5_and),
         T40.4_ratio = deaths_T40.4 / deaths_or,
         T40.5_ratio = deaths_T40.5 / deaths_or)
T40.4_NA_index <- which(is.na(overdose_matched_ratio$deaths_T40.4))
T40.5_NA_index <- which(is.na(overdose_matched_ratio$deaths_T40.5))
overdose_matched_ratio$deaths_T40.4[T40.4_NA_index] <- overdose_matched_ratio$deaths_or[T40.4_NA_index] + 
  overdose_matched_ratio$deaths_and[T40.4_NA_index] - overdose_matched_ratio$deaths_T40.5[T40.4_NA_index]
overdose_matched_ratio$deaths_T40.5[T40.5_NA_index] <- overdose_matched_ratio$deaths_or[T40.5_NA_index] + 
  overdose_matched_ratio$deaths_and[T40.5_NA_index] - overdose_matched_ratio$deaths_T40.4[T40.5_NA_index]
overdose_matched_ratio$T40.4_ratio[is.na(overdose_matched_ratio$T40.4_ratio)] <- 1 - overdose_matched_ratio$T40.5_ratio[is.na(overdose_matched_ratio$T40.4_ratio)]
overdose_matched_ratio$T40.5_ratio[is.na(overdose_matched_ratio$T40.5_ratio)] <- 1 - overdose_matched_ratio$T40.4_ratio[is.na(overdose_matched_ratio$T40.5_ratio)]
overdose_matched_ratio$n_NA <- overdose_matched_ratio %>% apply(1, function(x) sum(is.na(x[4:5])))
T40.4_NA_index <- which(is.na(overdose_matched_ratio$deaths_T40.4))
T40.5_NA_index <- which(is.na(overdose_matched_ratio$deaths_T40.5))

overdose_matched_ratio_summary <- overdose_matched_ratio %>%
  group_by(state) %>% 
  summarise(avg_T40.4_ratio = mean(T40.4_ratio, na.rm=T),
            sd_T40.4_ratio = sd(T40.4_ratio, na.rm=T),
            avg_T40.5_ratio = mean(T40.5_ratio, na.rm=T),
            sd_T40.5_ratio = sd(T40.5_ratio, na.rm=T)) %>% print(n=49)

overdose_matched_ratio <- overdose_matched_ratio %>% 
  left_join(overdose_matched_ratio_summary %>% select(-sd_T40.4_ratio), by="state")

overdose_matched_ratio$deaths_T40.4[T40.4_NA_index] <- ceiling(overdose_matched_ratio$deaths_or[T40.4_NA_index] * overdose_matched_ratio$avg_T40.4_ratio[T40.4_NA_index])
overdose_matched_ratio$deaths_T40.5[T40.5_NA_index] <- ceiling(overdose_matched_ratio$deaths_or[T40.5_NA_index] * overdose_matched_ratio$avg_T40.5_ratio[T40.5_NA_index])

overdose_matched_ratio$deaths_and <- overdose_matched_ratio$deaths_T40.4 +overdose_matched_ratio$deaths_T40.5 - overdose_matched_ratio$deaths_or
overdose_matched_ratio %>% filter(deaths_and < 0)
neg_death_and_index <- which(overdose_matched_ratio$deaths_and < 0)
overdose_matched_ratio$deaths_T40.5[neg_death_and_index] <- overdose_matched_ratio$deaths_T40.5[neg_death_and_index] -
  overdose_matched_ratio$deaths_and[neg_death_and_index]
overdose_matched_ratio$deaths_and[neg_death_and_index] <- 0

overdose_matched_ratio$deaths_T40.5_only <- overdose_matched_ratio$deaths_T40.5 - overdose_matched_ratio$deaths_and
overdose_matched_ratio %>%
  select(state, year, deaths_T40.5, deaths_T40.5_only) %>% 
  write.csv("Cocaine Network Optimization/Overdose deaths T40.5 only 1999-2020.csv", row.names=F)

overdose_matched_ratio %>% 
  select(state, year, deaths_T40.5, deaths_T40.5_only) %>% 
  mutate(T40.5_only_ratio=deaths_T40.5_only / deaths_T40.5) -> overdose_T40.5_only_summary
overdose_T40.5_only_summary$T40.5_only_ratio %>% summary
overdose_T40.5_only_summary$T40.5_only_ratio %>% sort
overdose_T40.5_only_summary %>% filter(T40.5_only_ratio == 0)
overdose_matched_ratio %>% filter(state == "Idaho")

overdose_T40.5_only_summary %>% filter(year %in% 2009:2014) %>% view
