source("functions.R")
library(splines)
library(covidcast)
library(dplyr)
library(tidyr)
library(ggplot2)

#################LR plots###################
#plot mae vs df 
maes_avg <- readRDS("Fits/lr_missing_maes_avg.rds")
maes_avg0 <- readRDS("Fits/lr_non-missing_maes_avg.rds")

rbind(data.frame(maes_avg, data = "all observed"),
      data.frame(maes_avg0, data = "fully-observed only")) %>%
ggplot()+
  geom_line(mapping = aes(df, log(test_mse), color = data, linetype = method))+
  scale_linetype_manual(values=c("smooth" = "solid", "baseline" = "dashed"))+
  ylab("log(test mae)")+
  xlab("degrees-of-freedom")+
  theme_bw()+
  scale_x_continuous(breaks = dfs)
ggsave(paste0("Plots/lr_mse_vs_df.pdf"), height = 3, width = 6)

#find optimal df
df_opt <- maes_avg %>% filter(method == "smooth") %>% slice(which.min(test_mse)) %>% pull(df)

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
ggsave(paste0("Plots/lr_missing_compare_fits.pdf"), height = 3.5, width = 7)


#plot mae vs ahead
maes_avg <- readRDS("Fits/lr_missing_maes.rds")

maes %>% 
  gather(aheads+1, key = "ahead", value = "mse") %>% mutate(ahead = as.numeric(ahead)) %>%
  mutate(model = paste(method, ifelse(is.na(df), "", df))) %>%
ggplot()+
  geom_line(aes(ahead, log(mse), color = model), size = 0.3)+
  scale_x_continuous(breaks = seq(min(aheads), max(aheads), 5))+
  ylab("log(test mae)")+
  theme_bw()
ggsave(paste0("Plots/lr_mse_vs_ahead.pdf"), height = 3, width = 6)

