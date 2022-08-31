source("functions.R")
library(splines)
library(covidcast)
library(dplyr)
library(tidyr)
library(ggplot2)
library(quantreg)


#load raw data
time_stamp <- as.Date("2021-10-01")
aheads <- 0:27
q <- length(aheads)
lags <- 1:28
signals <- c("confirmed_7dav_incidence_prop", "smoothed_hh_cmnty_cli", "smoothed_cli")
scenario <- "missing" #can be "missing" as well

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

#fit the models with df = 1,...,6
dfs = 1:6

################# LR ######################

#fit MRR
model <- mrr(train$X, train$Y)
mae <- colMeans(abs(test$Y - test$X %*% model$B), na.rm = TRUE)
maes_avg <- data.frame(mae = mean(mae), df = dfs, method = "baseline")
maes <- data.frame(mae, ahead = aheads, model = "baseline")
Ys <- list()
Ys[[q]] <- fit_merged(train, test, model, top_geo[1:5])

#fit MPF
for(df in dfs){
  cat(df, " ")
  #create basis
  H <- compute_H(0, aheads, df)
  #fit MPF
  model <- mpf(train$X, train$Y, H)
  #compute test score
  mae <- colMeans(abs(test$Y - test$X %*% model$B), na.rm = TRUE)
  maes_avg <- rbind(maes_avg, data.frame(mae = mean(mae), df = df, method = "smooth"))
  maes <- rbind(maes, data.frame(mae, ahead = aheads, model = paste("smooth", df)))
  Ys[[df]] <- fit_merged(train, test, model, top_geo[1:5])
}

saveRDS(Ys, paste0("Fits/lr_", scenario, "_fits.rds"))
saveRDS(maes, paste0("Fits/lr_", scenario, "_maes.rds"))
saveRDS(maes_avg, paste0("Fits/lr_", scenario, "_maes_avg.rds"))


################# QR+calibration ######################

train$info <- train$info %>% mutate(calibration = (time_stamp < max(train$stamps) - 28))
#fit QMRR
#no calibration
model <- qmrr(train$X, train$Y)
mcs_avg <- eval_qr(model, NULL, test, "baseline", dfs)
Bs <- list(qr = list(), qrc = list())
Bs[["qr"]][[q]] <- model$B
#calibration
model <- qmrr(train$X[!train$info$calibration,], train$Y[!train$info$calibration,])
calibrate <- correction(train$X[train$info$calibration,], train$Y[train$info$calibration,], model$B)
mcs_avg <- rbind(mcs_avg, eval_qr(model, calibrate, test, "baseline", dfs))
Bs[["qrc"]][[q]] <- model$B

#fit QMPF
for(df in dfs){
  cat(df, " ")
  #create basis
  H <- compute_H(0, aheads, df)
  #fit QMPF
  #no calibration
  model <- qmpf(train$X, train$Y, H)
  #compute test score
  mcs_avg <- rbind(mcs_avg, eval(model, NULL, test, "smooth", df))
  Bs[["qr"]][[q]] <- model$B
  #calibration
  model <- qmpf(train$X[!train$info$calibration,], train$Y[!train$info$calibration,], H)
  calibrate <- correction(train$X[train$info$calibration,], train$Y[train$info$calibration,], model$B)
  #compute test score
  mcs_avg <- rbind(mcs_avg, eval(model, calibrate, test, "smooth", df))
  Bs[["qrc"]][[q]] <- model$B
  saveRDS(Bs, paste0("Fits/qr_", scenario, "_coefs.rds"))
  saveRDS(mcs_avg, paste0("Fits/qr_", scenario, "_mcs_avg.rds"))
}

