source("functions.R")
library(dplyr)
library(ggplot2)

set.seed(0)

n = 1000 #location
q = 30 #ahead
A = 1:q
p = 10 #number of predictors
d = 3 #spline basis size
m = p #assume all lags = 0

#generate the data
generate = function(n, q, p, d, A){
  X = matrix(rnorm(n * p), n, p)
  H = compute_H(0, A, d)
  Theta = matrix(rnorm(ncol(H) * m), ncol(H), m) #spline coefficients
  B = t(H %*% Theta) #smooth model coefficients
  Y = X %*% B #response
  W = (runif(length(Y)) <= 0.1) #10% of missing values
  Y[W] = NA
  return(list(Y = Y, X = X))
}

#run simulation for df = 1...6
result = c()
dfs = 1:6

for(rep in 1:100){
  cat(rep, " ")
  #generate the data
  gen = generate(1000, 30, 10, 3, 1:30)
  X = gen$X
  Y = gen$Y
  #train-test split
  index = sample(1:n, n*0.5)
  for(SNR in c(0.1, 0.5, 1, 2)){
    #compute sd for errors
    s = abs(mean(Y, na.rm = TRUE)/SNR)
    E = matrix(rnorm(length(Y), 0, s), dim(Y))
    for(d in dfs){
      #create basis
      H = compute_H(0, A, d)
      #find mof solution
      sol = mpf(X[index,], (Y+E)[index,], H)
      #calculate train-test scores
      train = mean(abs(Y[index,] - X[index,] %*% sol$B), na.rm = T)
      test = mean(abs(Y[-index,] - X[-index,] %*% sol$B), na.rm = T)
      res = data.frame(df = d, score = c(train, test), score_type = c("train", "test"))
      result = rbind(result, data.frame(res, s, SNR, model = "smooth", rep = rep))
    }
    #find the solution when applying separate regressions
    H = compute_H(0, A, q)
    sol = mpf(X[index,], (Y+E)[index,], H)
    train = mean(abs(Y[index,] - X[index,] %*% sol$B), na.rm = T)
    test = mean(abs(Y[-index,] - X[-index,] %*% sol$B), na.rm = T)
    res = data.frame(df = d, score = c(train, test), score_type = c("train", "test"))
    result = rbind(result, data.frame(df = rep(dfs, times = rep(2, length(dfs))), score = res$score, score_type = res$score_type, s, SNR, model = "baseline", rep = rep))
  }
}

#compute average score and confidence intervals
result_avg = result %>% 
  group_by(df, score_type, SNR, model) %>% 
  mutate(mean = mean(score), se = sd(score)/sqrt(100), lower = quantile(score, 0.05), upper = quantile(score, 0.95))
           
#plot train-test performance
ggplot(data = result_avg %>% filter(score_type == "test"), mapping = aes(x = df, y = log(mean, 10), color = model, fill = model, linetype = model))+
  geom_line(size = 0.5)+
  geom_ribbon(aes(x = df, y = log(mean, 10), ymax = log(mean +  se, 10), ymin = log(mean - se, 10)), alpha=0.3, color = NA)+
  facet_wrap(~SNR, ncol = 5, labeller = label_bquote(cols = SNR == .(SNR)))+ 
  xlab("degrees-of-freedom") + ylab("log(test mae)")+
  scale_x_continuous(breaks = dfs)+
  labs(color = "type", linetype = "type")+
  theme_bw()+
  theme(legend.position="none")+
  scale_color_manual(values=c("smooth" = "black", "baseline" = "red"))+
  scale_fill_manual(values=c("smooth" = "black", "baseline" = "red"))+
  scale_linetype_manual(values=c("smooth" = "solid", "baseline" = "dashed"))
ggsave("Plots/simulation.pdf", height = 2.5, width = 6)
