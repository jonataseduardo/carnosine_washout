library(data.table)
library(ggplot2)
library(brms)
library(propagate)

join_fit <- 
  function(classical_linear, bayes_linear, classical_exp, bayes_exp){
    ci_cl <- 
      as.data.table(
        predict(classical_linear, 
                newdata = data.table(time_stamp = 0:16), 
                interval = 'confidence'))
    ci_cl[, `:=`(fit_type = 'classical linear', time_stamp = 0:16)]

    ci_bl <- 
      as.data.table(
        fitted(bayes_linear,
               newdata = data.table(time_stamp = 0:16)
               ))[, c(1,3,4)]
    ci_bl[, `:=`(fit_type = 'bayesian linear', time_stamp = 0:16)]

    ci_ce <- 
      as.data.table(
        predictNLS(classical_exp, 
                   newdata = data.table(time_stamp = 0:16), 
                   interval = 'confidence')$summary
        )[, c('Prop.Mean.1', 'Sim.2.5%', 'Sim.97.5%')]
    ci_ce[, `:=`(fit_type = 'classical exponential', time_stamp = 0:16)]

    ci_be <- 
      as.data.table(
        fitted(bayes_exp,
               newdata = data.table(time_stamp = 0:16), 
               ))[, c(1,3,4)]
    ci_be[, `:=`(fit_type = 'bayesian exponential', time_stamp = 0:16)]


    estimates <- rbindlist(list(ci_cl, ci_bl, ci_ce, ci_be))
    return(estimates)
  }

  {
  washout_raw <- fread('carnosine_data.csv')

  washout_0 <- 
    melt(washout_raw, id = c('group_id', 'sample_id'), 
         variable.name = 'week_id')

  week_info <- 
    data.table(week_id= c("PRE_S01", "POS_S02", "W01", 
                          "W02", "W03", "W04", "W05", "W06"),
               time_stamp = c(NA, 0, 1, 2, 4, 8, 12, 16)) 

  washout <- washout_0[week_info, on = c("week_id")][!is.na(value)]

  ## Subtracting initial cardosine level from idividuals
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
  } 

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


#############################
# Fiting full data
###########################

## Classical and Bayesian linear fit
  {
  classical_linear <- 
    lm(formula = value ~ time_stamp, 
       data = washout[group_id == 'T' & !is.na(time_stamp)])

  bayes_linear <- 
    brm(value ~ time_stamp,
        data = washout[group_id == 'T' & !is.na(time_stamp)])

  }
## Classical and Bayesian Exponential Fit
  {
  classical_exp <- 
    nls( value ~ SSasymp(time_stamp, Asym , R0, lrc), 
        data = washout[group_id == 'T' & !is.na(time_stamp)])

  #estimating_priors
  ## values b3 prior
  washout[is.na(time_stamp) & group_id == 'T', .(mean(value), sd(value))]

  ## values for b1 prior
  washout[time_stamp == 0 & group_id == 'T', .(mean(value), sd(value))]
  40.75576 - 22.02241

  exp_prior <- 
    prior(normal(18.73335, 12.28948), nlpar = "b1") +
    #prior(normal(1, 5), lb = 0, nlpar = "b2" ) +
    prior(student_t(3, 0, 13), lb = 0, nlpar = "b2" ) +
    prior(normal(22.02241, 8.636423), nlpar = "b3") 

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
  }

##########################################
#Analysis removing offset 
###########################################

#Classical and Bayesian statist analysis for linear model 
  {

  classical_linear_off <- 
    lm(formula = value ~ time_stamp, 
       data = washout_off[group_id == 'T'])

  bayes_linear_off <- 
    brm(value ~ time_stamp,
        data = washout_off[group_id == 'T'])
  }

#Classical and Bayesian exponential model 
  {
  classical_exp_off <- 
    nls( value ~ SSasymp(time_stamp, Asym, R0, lrc), 
        data = washout_off[group_id == 'T'])

  ## values for b1 prior
  washout_off[time_stamp == 0 & group_id == 'T', .(mean(value), sd(value))]

  exp_prior_off <- 
    prior(normal(18.73336, 5.801748), nlpar = "b1") +
    #prior(normal(0.0, 5), lb = 0, nlpar = "b2" )
    prior(student_t(3, 0, 13), lb = 0, nlpar = "b2" )

  bayes_exp_off <- 
    brm(bf(value  ~ b1 * exp(- b2 * time_stamp),  
           b1 + b2 ~ 1,
           nl = TRUE,
           cmc = FALSE),
        data = washout_off[group_id == 'T'], 
        prior = exp_prior_off,
        control = list(adapt_delta = 0.99))

  }


estimates <- join_fit(classical_linear, bayes_linear, classical_exp, bayes_exp)

  ggplot(data = estimates, aes(x = time_stamp)) + 
  geom_line(aes(y = fit, color = fit_type), size = 2) + 
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = fit_type), alpha = 0.1) + 
  geom_point(data = washout[group_id == 'T'], aes(x = time_stamp, y = value)) + 
  theme_bw()



exp(coef(classical_exp_off)[['lrc']])
summary(bayes_exp)
plot(bayes_exp)
plot(marginal_effects(bayes_exp), points = TRUE)

loo(bayes_linear, bayes_exp)
loo(bayes_linear_off, bayes_exp_off)
summary(bayes_exp_off)

plot(bayes_exp_off)
plot(marginal_effects(bayes_exp_off), points = TRUE)
summary(bayes_linear_off)
plot(bayes_linear_off)
plot(marginal_effects(bayes_linear_off), points = TRUE)
