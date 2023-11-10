source("functions.R")
library(splines)
library(covidcast)
library(dplyr)
library(tidyr)
library(ggplot2)
library(quantreg)

#Create "Fits" folder
#then run LR and QR parts for ""missing and then "non-missing" arguments
scenario <- "missing" 

#load raw data
time_stamp <- as.Date("2021-10-01")
aheads <- 0:27
q <- length(aheads)
lags <- 1:28
signals <- c("confirmed_7dav_incidence_prop", "smoothed_hh_cmnty_cli", "smoothed_cli")

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

#fit the models with df = 1,...,6
dfs = 1:6

################# LR ######################
score_loc <- function(scores){
  data.frame(scores, geo_value = test$info$geo_value) %>%
    group_by(geo_value) %>%
    summarise_at(vars(!starts_with("geo_value")), ~mean(., na.rm = T)) %>%
    select(-geo_value) %>% rowMeans(na.rm = T)
}

#fit MRR
model <- mrr(train$X, train$Y)
mae <- abs(test$Y - test$X %*% model$B)
maes_avg <- data.frame(mae = mean(colMeans(mae, na.rm = TRUE)), df = dfs, method = "baseline")
maes <- data.frame(mae = colMeans(mae, na.rm = TRUE), ahead = aheads, model = "baseline")
maes_loc <- data.frame(mae = score_loc(mae), df = 0, method = "baseline")

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
  mae <- abs(test$Y - test$X %*% model$B)
  maes_avg <- rbind(maes_avg, data.frame(mae = mean(colMeans(mae, na.rm = TRUE)), df = df, method = "smooth"))
  maes <- rbind(maes, data.frame(mae = colMeans(mae, na.rm = TRUE), ahead = aheads, model = paste("smooth", df)))
  maes_loc <- rbind(maes_loc, data.frame(mae = score_loc(mae), df = df, method = "smooth"))
  Ys[[df]] <- fit_merged(train, test, model, top_geo[1:5])
}

saveRDS(Ys, paste0("Fits/lr_", scenario, "_fits.rds"))
saveRDS(maes, paste0("Fits/lr_", scenario, "_maes.rds"))
saveRDS(maes_avg, paste0("Fits/lr_", scenario, "_maes_avg.rds"))
saveRDS(maes_loc, paste0("Fits/lr_", scenario, "_maes_loc.rds"))

################# QR ######################

train$info <- train$info %>% mutate(calibration = (time_stamp < max(train$stamps) - 28))
#fit QMRR
#no calibration
model <- qmrr(train$X, train$Y)
eval <- eval_qr(model, NULL, test)
mcs_avg <- merge(data.frame(eval$mcs_avg, method = "baseline"), data.frame(df = dfs))
mcs_loc <- data.frame(eval$mcs_loc, method = "baseline", df = 0)
Bs <- list(qr = list(), qrc = list())
Bs[["qr"]][[q]] <- model$B
#calibration
model <- qmrr(train$X[!train$info$calibration,], train$Y[!train$info$calibration,])
calibrate <- correction(train$X[train$info$calibration,], train$Y[train$info$calibration,], model$B)
eval <- eval_qr(model, calibrate, test)
mcs_avg <- rbind(mcs_avg, merge(data.frame(eval$mcs_avg, method = "baseline"), data.frame(df = dfs)))
mcs_loc <- rbind(mcs_loc, data.frame(eval$mcs_loc, method = "baseline", df = 0))
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
  eval <- eval_qr(model, NULL, test)
  mcs_avg <- rbind(mcs_avg, data.frame(eval$mcs_avg, method = "smooth", df = df))
  mcs_loc <- rbind(mcs_loc, data.frame(eval$mcs_loc, method = "smooth", df = df))
  Bs[["qr"]][[q]] <- model$B
  #calibration
  model <- qmpf(train$X[!train$info$calibration,], train$Y[!train$info$calibration,], H)
  calibrate <- correction(train$X[train$info$calibration,], train$Y[train$info$calibration,], model$B)
  #compute test score
  eval <- eval_qr(model, calibrate, test)
  mcs_avg <- rbind(mcs_avg, data.frame(eval$mcs_avg, method = "smooth", df = df))
  mcs_loc <- rbind(mcs_loc, data.frame(eval$mcs_loc, method = "smooth", df = df))
  Bs[["qrc"]][[q]] <- model$B
  saveRDS(Bs, paste0("Fits/qr_", scenario, "_coefs.rds"))
  saveRDS(mcs_avg, paste0("Fits/qr_", scenario, "_mcs_avg.rds"))
  saveRDS(mcs_loc, paste0("Fits/qr_", scenario, "_mcs_loc.rds"))
}

###################### Plots #########################
#run after completing computations for scenario = "missing" and then "non-missing"

time_stamp <- as.Date("2021-10-01")
aheads <- 0:27
q <- length(aheads)
lags <- 1:28
signals <- c("confirmed_7dav_incidence_prop", "smoothed_hh_cmnty_cli", "smoothed_cli")
dfs <- 1:6

################# LR ###################
#plot mae vs df 
maes_avg <- readRDS("Fits/lr_missing_maes_avg.rds")
maes_avg0 <- readRDS("Fits/lr_non-missing_maes_avg.rds")

rbind(data.frame(maes_avg, data = "all observed"),
      data.frame(maes_avg0, data = "fully-observed only")) %>%
  ggplot()+
  geom_line(mapping = aes(df, log(mae), color = data, linetype = method))+
  scale_linetype_manual(values=c("smooth" = "solid", "baseline" = "dashed"))+
  ylab("log(test mae)")+
  xlab("degrees-of-freedom")+
  theme_bw()+
  scale_x_continuous(breaks = dfs)

#find optimal df
df_opt <- maes_avg %>% filter(method == "smooth") %>% slice(which.min(mae)) %>% pull(df)

#test baseline vs optimal
maes_loc <- readRDS("Fits/lr_missing_maes_loc.rds")
maes_loc0 <- readRDS("Fits/lr_non-missing_maes_loc.rds")

b <- maes_loc %>% filter(df == 0) %>% pull(mae)
o <- maes_loc %>% filter(df == df_opt) %>% pull(mae)
wilcox.test(b[!is.na(b)], o[!is.na(b)], alternative = "greater", paired = TRUE)
prop.table(table(b[!is.na(b)] - o[!is.na(b)] > 0))

#compare fits for optimal df and MRR
Ys <- readRDS("Fits/lr_missing_fits.rds")

rbind(data.frame(Ys[[df_opt]], method = paste("smooth", df_opt)), data.frame(Ys[[q]], method = "baseline")) %>%
  ggplot()+
  geom_line(aes(time, cases, color = type, size = type, alpha = type, group = interaction(time_stamp, type)))+
  facet_grid(method~geo_value, scales = "free")+
  theme_bw()+
  theme(legend.position="top")+
  theme(legend.title = element_blank())+
  scale_size_manual(values=c("predicted train" = 0.4, "predicted test" = 0.4, "observed train" = 0.7, "unobserved test" = 0.7))+
  scale_color_manual(values=c("predicted train" = "dodgerblue1", "predicted test" = "tomato1", "observed train" = "dodgerblue4", "unobserved test" = "tomato4"))+
  scale_alpha_manual(values=c("predicted train" = 0.8, "predicted test" = 0.8, "observed train" = 1, "unobserved test" = 1))

#plot mae vs ahead
maes <- readRDS("Fits/lr_missing_maes.rds")

ggplot(maes)+
  geom_line(aes(ahead, log(mae), color = model), size = 0.3)+
  scale_x_continuous(breaks = seq(min(aheads), max(aheads), 5))+
  xlab("horizon")+
  ylab("log(test mae)")+
  theme_bw()

################# QR ###################
mcs_avg <- readRDS("Fits/qr_missing_mcs_avg.rds") %>% mutate(calibration = factor(calibration, levels = c("no calibration", "calibration")))

mcs_avg %>% filter(type == "mae") %>% 
  ggplot()+
  geom_line(mapping = aes(df, log(score), color = method, linetype = method))+
  scale_linetype_manual(values=c("smooth" = "solid", "baseline" = "dashed"))+
  scale_color_manual(values=c("smooth" = "black", "baseline" = "red"))+
  ylab("log(test mae)")+
  xlab("degrees-of-freedom")+
  theme_bw()+
  scale_x_continuous(breaks = dfs)+
  facet_wrap(~calibration)

mcs_loc <- readRDS("Fits/qr_missing_mcs_loc.rds") 

b <- mcs_loc %>% filter(df == 0, type == "mae", calibration == "no calibration") %>% pull(score)
o <- mcs_loc %>% filter(df == 3, type == "mae", calibration == "no calibration") %>% pull(score)
wilcox.test(b[!is.na(b)], o[!is.na(b)], alternative = "greater", paired = TRUE)
prop.table(table(b[!is.na(b)] - o[!is.na(b)] > 0))





