
########## method ########## 

compute_xy = function(X, Y, H){
  index = (rowSums(is.na(X))==0)
  X = X[index,]
  Y = Y[index,]
  y = c(Y)
  x = H %x% X
  x = x[!is.na(y),]
  y = y[!is.na(y)]
  w = rep(1, length(y))
  return(list(x = x, y = y, w = w))
}

mpf = function(X, Y, H){
  reshape = compute_xy(X, Y, H)
  theta = lm(reshape$y ~ reshape$x + 0, weights = reshape$w)$coef
  Theta = matrix(theta, nrow = ncol(H), byrow = TRUE)
  B = t(H %*% Theta)
  return(list(Theta = Theta, B = B))
} 

mrr = function(X, Y){
  B = matrix(0, ncol(X), ncol(Y))
  for(i in 1:ncol(Y)){
    y = Y[, i]
    x = X[!is.na(y),]
    y = y[!is.na(y)]
    B[, i] = lm(y ~ x + 0)$coef
  }
  return(list(B = B))
} 

########## basis ##########

compute_H = function(time_stamp, aheads, d){
  H = matrix(1/sqrt(length(aheads)), length(aheads), 1)
  if(d > 1 & d < length(aheads)) H = cbind(H, poly(time_stamp + aheads, df = d - 1, raw = F))
  if(d == length(aheads)) H = diag(d)
  return(H)
}

########## data preprocessing ##########

create_data = function(df, signals, aheads, lags, top_geo){
  list(Y = df %>% 
         filter(geo_value %in% top_geo) %>%
         select(paste0(signals[1], "+", aheads)) %>%
         as.matrix(),
       X = df %>% 
         filter(geo_value %in% top_geo) %>%
         select(as.vector(outer(signals, lags, paste, sep = "-"))) %>% 
         mutate("intercept" = 1) %>%
         as.matrix(),
       info = df %>% 
         filter(geo_value %in% top_geo) %>%
         select(geo_value, time_stamp),
       stamps = unique(df$time_stamp))
}

remove_missing = function(df){
  index <- (rowSums(is.na(df$Y))==0)
  list(Y = df$Y[index,],
       X = df$X[index,],
       info = df$info[index,],
       stamps = df$stamps[index])
}

reshape_y = function(Y, info){
  aheads <- find_aheads(Y)
  cbind(Y, info) %>% 
    gather(1:ncol(Y), key = "ahead", value = "cases") %>%
    mutate(ahead = factor(ahead, levels = paste0('confirmed_7dav_incidence_prop+', aheads), labels = aheads)) 
}

find_top_geo <- function(df, n){
  df <- na.omit(df)
  df %>% 
    group_by(geo_value) %>% 
    summarise_at(vars(contains("+")), ~sum(., na.rm = T))  %>%
    mutate(total_cases = rowSums(select(., contains("+")), na.rm = TRUE)) %>%
    arrange(desc(total_cases)) %>%
    top_n(n, total_cases) %>% 
    select(geo_value) %>%
    pull()
}

find_aheads = function(Y){
  as.numeric(substr(colnames(test$Y), nchar(signals[1])+2, nchar(signals[1])+3))
}

########## plots ##########

plot_missing = function(Y, info, Ytest, infotest){
  aheads <- find_aheads(Y)
  stamps <- unique(info$time_stamp)
  sub <- info$geo_value == info$geo_value[1]
  nas <- reshape_y(Y[sub,], info[sub,]) %>% mutate(missing = is.na(cases), set = "observed train")
  aheadstest <- find_aheads(Ytest)
  stampstest <- unique(infotest$time_stamp)
  subtest <- infotest$geo_value == info$geo_value[1]
  nastest <- reshape_y(Ytest[subtest,], infotest[subtest,]) %>% mutate(missing = is.na(cases), set = "unobserved test")
  ggplot(rbind(nas, nastest))+
    geom_tile(aes(as.numeric(ahead), time_stamp, fill = interaction(set, missing)))+ 
    #scale_fill_manual(values=c("dodgerblue4", "tomato4"))+
    scale_fill_manual(name = "legend",
                      values = c(
                        alpha("dodgerblue4", 1),
                        alpha("tomato4", 1),
                        alpha("dodgerblue4", 0.1)),
                      labels = c("observed train", 
                                 "unobserved test",
                                 "unobserved train"))+
    scale_x_continuous(breaks = seq(min(aheads), max(aheads), 5))+
    scale_y_continuous(breaks = max(stampstest) - seq(0, as.numeric(max(stampstest) - min(stamps)), 7))+
    xlab("ahead")+
    coord_flip()+
    theme_bw()+
    theme(legend.position="top")+
    theme(legend.title = element_blank())+
    theme(axis.text.x = element_text(size=8, angle=30, hjust = 1))
}

fit_merged = function(train, test, model, geo_values){
  Yhat_train <- train$X %*% model$B
  colnames(Yhat_train) <- colnames(train$Y)
  Yhat_test <- test$X %*% model$B
  colnames(Yhat_test) <- colnames(test$Y)
  Y <- rbind(data.frame(reshape_y(train$Y, train$info), "type" = "observed train"),
             data.frame(reshape_y(Yhat_train, train$info), "type" = "predicted train"),
             data.frame(reshape_y(test$Y, test$info), "type" = "unobserved test"),
             data.frame(reshape_y(Yhat_test, test$info), "type" = "predicted test")) %>% 
    mutate(type = factor(type, levels = c("predicted train", "predicted test", "observed train", "unobserved test"))) %>%
    mutate(ahead = as.numeric(ahead)) %>%
    mutate(time = time_stamp + ahead) %>%
    filter(geo_value %in% geo_values)
  return(Y)
}


