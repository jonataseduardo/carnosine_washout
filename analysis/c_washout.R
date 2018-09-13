library(data.table)
library(ggplot2)
library(brms)
library(rstan)
library(mice)

washout_raw <- fread('carnosine_data.csv')
washout <- melt(washout_raw, id = c(1,10), 
                variable.name = 'week_id')

week_info <- 
  data.table(week_id= c("PRE_S01", "POS_S02", "W01", 
                        "W02", "W03", "W04", "W05", "W06"),
             time_stamp = c(NA, 0, 1, 2, 4, 8, 12, 16)) 

washout <- washout[week_info, on = c("week_id")]


decay_sample <- 
  washout[!is.na(time_stamp) & group_id == 'T' & !is.na(value), 
          .(sample_id, value, time_stamp)]


prior1 <- prior(normal(2,3), nlpar = "b1") +
          prior(normal(-4,2), nlpar = "b2") +
          prior(normal(15,5), nlpar = "b3")

fit1 <- brm(bf(value ~ b1 * exp(b2 * time_stamp) + b3, 
               b1 ~ 1 + (1|sample_id),
               b2 ~ 1 + (1|sample_id),
               b3 ~ 1 + (1|sample_id),
               b1 + b2 + b3 ~ 1, 
               nl = TRUE),
            data = decay_sample, 
            prior = prior1,
            control = list(adapt_delta = 0.95))

summary(fit1)
plot(fit1)
plot(marginal_effects(fit1), points = TRUE)

id <- decay_sample[, .GRP, sample_id][, .(sample_id)]

meffects <- marginal_effects(fit1, conditions =  id, 
                             re_formula = NULL, method = "predict")

plot(meffects, points = TRUE)

ggplot(washout, aes(x = time_stamp, y = value, color = sample_id)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~ group_id)

warnings()
