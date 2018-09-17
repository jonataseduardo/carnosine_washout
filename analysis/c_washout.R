library(data.table)
library(ggplot2)
library(brms)

washout_raw <- fread('carnosine_data.csv')

washout_0 <- 
  melt(washout_raw, id = c('group_id', 'sample_id'), 
       variable.name = 'week_id')

week_info <- 
  data.table(week_id= c("PRE_S01", "POS_S02", "W01", 
                        "W02", "W03", "W04", "W05", "W06"),
             time_stamp = c(NA, 0, 1, 2, 4, 8, 12, 16)) 

washout <- washout_0[week_info, on = c("week_id")][!is.na(value)]

# First we will test if the samples in the target group are draw from the
# sample distribution from the control at time 0

## Result
#D = 0.36364, p-value = 0.7473
#alternative hypothesis: two-sided
# There is no evidence that the two distribution are diferents from each
# other

# Withing individual variaility is lower then the group variability analisys 
washout[group_id == 'C', var(value)]
washout[group_id == 'C', var(value), by = sample_id]


t.test(washout[is.na(time_stamp) & group_id == 'T', value], 
       washout[is.na(time_stamp) & group_id == 'C', value])

wilcox.test(washout[is.na(time_stamp) & group_id == 'T', value], 
            washout[is.na(time_stamp) & group_id == 'C', value])

# Sencond we want to test if the distribution of the cardosine in the subjects
# are diferent for each time 

init_s <- washout[is.na(time_stamp) & group_id == 'T', value]

diff_time_tests <- 
  washout[!is.na(time_stamp) & time_stamp != 1 & group_id == 'T', 
          .(wc_pval = wilcox.test(init_s, value, alternative = 'less')$p.value,
            t_pval = t.test(init_s, value, alternative = 'less')$p.value), 
          by = time_stamp]

diff_time_tests

#Classical and Bayesian statist analysis for linear model 
classical_linear <- 
  lm(formula = value ~ time_stamp, 
     data = washout[group_id == 'T' & !is.na(time_stamp)])

summary(classical_linear)

bayes_linear <- 
  brm(value ~ time_stamp,
      data = washout[group_id == 'T' & !is.na(time_stamp)])

summary(bayes_linear)
plot(bayes_linear)
plot(marginal_effects(bayes_linear), points = TRUE)


classical_exp <- 
  nls( value ~ SSasymp(time_stamp, Asym, R0, lrc), 
      data = washout[group_id == 'T' & !is.na(time_stamp)])

summary(classical_exp)

ggplot(data = washout[group_id == 'T' & !is.na(time_stamp)], 
       aes(x = time_stamp, y = value)) +
  geom_point() + 
  stat_smooth(method = 'nls',
              formula = y ~ SSasymp(x, Asym, R0, lrc),
              se = FALSE
              )

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


exp_prior <- 
  prior(normal(21.16837, 17.43511), nlpar = "b1") +
  prior(normal(0.5, 1), lb = 0, nlpar = "b2" ) +
  prior(normal(19.58739, 5.145628), nlpar = "b3") 

bayes_exp <- 
  brm(bf(value  ~ b1 * exp(- b2 * time_stamp) + b3,  
         b1 ~ 1,
         b2 ~ 1,
         b3 ~ 1, 
         nl = TRUE,
         cmc = FALSE),
      data = washout[group_id == 'T'& !is.na(time_stamp)], 
      prior = exp_prior,
      control = list(adapt_delta = 0.99))

summary(bayes_exp)
plot(bayes_exp)
plot(marginal_effects(bayes_exp), points = TRUE)

loo(bayes_linear, bayes_exp)

##########################################
#Same analysis removing offset 
###########################################

w_off <- 
  washout_raw[,lapply(names(washout_raw)[3:9], 
                       function(x) get(x) - get('PRE_S01'))] 

names(w_off) <- names(washout_raw)[3:9]
w_off[, sample_id := washout_raw[, sample_id]]
w_off[, group_id := washout_raw[, group_id]]

lw_off <- 
  melt(w_off, id = c('group_id', 'sample_id'), 
       variable.name = 'week_id')

washout_off <- lw_off[week_info, on = c("week_id")][!is.na(value)]


#Classical and Bayesian statist analysis for linear model 
classical_linear_off <- 
  lm(formula = value ~ time_stamp, 
     data = washout_off[group_id == 'T'])

summary(classical_linear_off)
summary(classical_linear)

bayes_linear_off <- 
  brm(value ~ time_stamp,
      data = washout_off[group_id == 'T'])

prior_summary(bayes_linear)

summary(bayes_linear_off)
plot(bayes_linear_off)
plot(marginal_effects(bayes_linear_off), points = TRUE)

classical_exp_off <- 
  nls( value ~ SSasymp(time_stamp, Asym, R0, lrc), 
      data = washout_off[group_id == 'T'])

summary(classical_exp_off)

ggplot(data = washout_off[group_id == 'T'], 
       aes(x = time_stamp, y = value)) +
  geom_point() + 
  stat_smooth(method = 'nls',
              formula = y ~ SSasymp(x, Asym, R0, lrc),
              se = FALSE
              )

exp_prior_off <- 
  prior(normal(21.16837, 17.43511), nlpar = "b1") +
  prior(normal(0.5, 1), lb = 0, nlpar = "b2" )

bayes_exp_off <- 
  brm(bf(value  ~ b1 * exp(- b2 * time_stamp),  
         b1 ~ 1,
         b2 ~ 1,
         nl = TRUE,
         cmc = FALSE),
      data = washout_off[group_id == 'T'], 
      prior = exp_prior_off,
      control = list(adapt_delta = 0.99))

summary(bayes_exp_off)
plot(bayes_exp_off)
plot(marginal_effects(bayes_exp_off), points = TRUE)

loo(bayes_linear_off, bayes_exp_off)

ggplot(washout, aes(x = time_stamp, y = value, color = sample_id)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~ group_id)
