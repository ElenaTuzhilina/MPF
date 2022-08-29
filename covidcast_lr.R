source("functions.R")
library(splines)
library(covidcast)
library(dplyr)
library(tidyr)
library(ggplot2)

#load raw data
time_stamp <- as.Date("2021-10-01")
aheads <- 0:27
q <- length(aheads)
lags <- 1:28
signals <- c("confirmed_7dav_incidence_prop", "smoothed_hh_cmnty_cli", "smoothed_cli")
scenario <- "non-missing" #can be "missing" as well

#train and test sets
df_train <- readRDS(paste0("Data/", time_stamp, "_prop_", scenario ,".rds"))
df_test <- readRDS(paste0("Data/", time_stamp, "_prop_test.rds"))

#filter 300 counties
top_geo <- find_top_geo(df_train, 300)

#preprocessing
train <- create_data(df_train, signals, aheads, lags, top_geo)
if(scenario == "non-missing") train <- remove_missing(train)
test <- create_data(df_test, signals, aheads, lags, top_geo)

#plot missing values
plot_missing(train$Y, train$info, test$Y, test$info)
ggsave(paste0("Plots/", scenario, "_missing.pdf"), height = 3, width = 6)

#fit MRR
model <- mrr(train$X, train$Y)
mae <- colMeans(abs(test$Y - test$X %*% model$B), na.rm = TRUE) %>% t() %>% data.frame()
colnames(mae) <- aheads
maes <- mae %>% mutate(method = "baseline", df = NA)
Ys <- list()
Ys[[q]] <- fit_merged(train, test, model, top_geo[1:5])

#fit the models with df = 1,...,6
dfs = 1:6
for(df in dfs){
  cat(df, " ")
  #create basis
  H <- compute_H(0, aheads, df)
  #fit MPF
  model <- mpf(train$X, train$Y, H)
  #compute test score
  mae <- colMeans(abs(test$Y - test$X %*% model$B), na.rm = TRUE) %>% t() %>% data.frame()
  colnames(mae) <- aheads
  maes <- rbind(maes, mae %>% mutate(method = "smooth", df = df))
  Ys[[df]] <- fit_merged(train, test, model, top_geo[1:5])
}

saveRDS(maes, paste0("Fits/lr_", scenario, "_maes.rds"))
saveRDS(Ys, paste0("Fits/lr_", scenario, "_fits.rds"))

#compute average score across aheads 
maes_avg <- maes %>% mutate(test_mse = rowMeans(select(maes, aheads))) %>%
  select(method, df, test_mse) 
maes_avg <- maes_avg %>% filter(method == "baseline") %>% select(method, test_mse) %>%
  data.frame(df = dfs) %>% 
  merge(maes_avg %>% filter(method == "smooth"), all = T)

saveRDS(maes_avg, paste0("Fits/lr_", scenario, "_maes_avg.rds"))
