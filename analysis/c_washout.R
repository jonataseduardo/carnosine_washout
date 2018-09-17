library(data.table)
library(ggplot2)
library(brms)

washout_raw <- fread('carnosine_data.csv')

week_info <- 
  data.table(week_id= c("PRE_S01", "POS_S02", "W01", 
                        "W02", "W03", "W04", "W05", "W06"),
             time_stamp = c(NA, 0, 1, 2, 4, 8, 12, 16)) 

washout <- washout[week_info, on = c("week_id")]

washout <- melt(washout_raw, id = c('group_id', 'sample_id'), 
                variable.name = 'week_id')


# First we will test if the samples in the target group are draw from the
# sample distribution from the control at time 0

ks.test(washout[is.na(time_stamp) & group_id == 'T', value], 
        washout[is.na(time_stamp) & group_id == 'C', value])

## Result
#D = 0.36364, p-value = 0.7473
#alternative hypothesis: two-sided
# There is no evidence that the two distribution are diferents from each
# other

# Withing individual variaility is lower then the group variability analisys 
washout[ group_id == 'C', var(value, na.rm = TRUE)]
washout[ group_id == 'C', var(value, na.rm = TRUE), by = sample_id]

#w_off <- 
#  washout_raw[,lapply(names(w_treated)[3:9], 
#                       function(x) get(x) - get('PRE_S01'))] 
#
#names(w_off) <- names(w_treated)[3:9]
#w_off[, sample_id := w_treated[, sample_id]]

decay_sample <- 
  washout[!is.na(value) & is.na(time_stamp), 
          .(sample_id, group_id, value, time_stamp)]

decay_sample[, sample_id := as.factor(sample_id)]
decay_sample[, group_id := as.factor(group_id)]






linear_fit <- 
  brm(value ~ time_stamp,
      data = decay_sample[group_id == 'T'])

summary(linear_fit)
plot(linear_fit)
plot(marginal_effects(linear_fit), points = TRUE)


mu_control <- 
  washout[ group_id == 'C' | week_id == 'PRE_S01', 
          mean(value, na.rm = TRUE)] 

sigma_control <- 
  washout[group_id == 'C' | week_id == 'PRE_S01', 
          sd(value, na.rm = TRUE)] 

mu_0 <- 
  washout[ group_id == 'T' & time_stamp == '0', 
          mean(value, na.rm = TRUE)] 
sigma_0 <- 
  washout[ group_id == 'T' & time_stamp == '0', 
          sd(value, na.rm = TRUE)] 


washout[,(max(value, na.rm = TRUE) - min(value, na.rm = TRUE))/16 ] 
(mu_0 - mu_control) / 16

3.276344 - 1.323023


sigma_0 + sigma_control

exp_prior <- 
  prior(normal(21.16837, 17.43511), nlpar = "b1") +
  #prior(normal(1.323023, 1.953321), lb = 0, nlpar = "b2" ) +
  prior(gamma(2,1), lb = 0, nlpar = "b2" ) +
  prior(normal(19.58739, 5.145628), nlpar = "b3") 

exponetial_fit <- 
  brm(bf(value  ~ b1 * exp(- b2 * time_stamp) + b3,  
         b1 ~ 1,
         b2 ~ 1,
         b3 ~ 1, 
         nl = TRUE,
         cmc = FALSE),
      data = decay_sample[group_id == 'T'& !is.na(time_stamp)], 
      prior = exp_prior,
      control = list(adapt_delta = 0.99))

summary(linear_fit)
summary(exponetial_fit)
plot(exponetial_fit)
18.53 * 0.47

gpl <- plot(marginal_effects(exponetial_fit), points = TRUE)

loo(exponetial_fit, linear_fit)

fit_classical <- 
  nls( value ~ SSasymp(time_stamp, Asym, R0, lrc), 
      data = decay_sample[group_id == 'T' & !is.na(time_stamp)])

summary(fit_classical)

gg <- decay_sample[group_id == 'T' & !is.na(time_stamp)]
ggplot(gg, aes(x = time_stamp, y = value)) +
  geom_point() + 
  stat_smooth(method = 'nls',
              formula = y ~ SSasymp(x, Asym, R0, lrc),
              se = FALSE
              )

ggplot(washout, aes(x = time_stamp, y = value, color = sample_id)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~ group_id)

