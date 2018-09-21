library(data.table)
library(ggplot2)
library(brms)
library(propagate)
library(knitr)
library(devtools)
#install_github('paul-buerkner/brms') 


make_labels <- 
  function(classical_linear, bayes_linear, classical_exp, bayes_exp){

    linear_fit_label <- 
      function(fl){
        paste0(format(fl[1], digits = 3), 
               ' - ', 
               format(-1 * fl[2], digits = 3), ' t')
      }

    c_exp_fit_label <- 
      function(fl){
        paste0(format(fl[2], digits = 3),
               ' Exp(-', 
               format(exp(fl[3]), digits = 2), 
               ' t) + ',
               format(fl[1], digits = 3))
      }

    b_exp_fit_label <- 
      function(fl){
        if(length(fl) == 3){
          paste0(format(fl[1], digits = 3),
                 ' Exp(-', 
                 format(fl[2], digits = 2), 
                 ' t) + ',
                 format(fl[3], digits = 3))
        }else{
          paste0(format(fl[1], digits = 3),
                 ' Exp(-', 
                 format(fl[2], digits = 2), 
                 ' t)')
        }
      }

    l_cl <- summary(classical_linear)$coefficients[,'Estimate']
    l_ce <- summary(classical_exp)$coefficients[,'Estimate']
    l_be <- summary(bayes_exp)$fixed[,'Estimate']
    l_bl <- summary(bayes_linear)$fixed[,'Estimate']

    label_dt <- 
      data.table(statistics = rep(c('classical', 'bayesian'), each = 2),
                 func = rep(c('linear', 'exponential'), 2),
                 label = c(linear_fit_label(l_cl),
                           c_exp_fit_label(l_ce),
                           linear_fit_label(l_bl),
                           b_exp_fit_label(l_be)))

    return(label_dt[])
  }

join_fit <- 
  function(classical_linear, bayes_linear, classical_exp, bayes_exp){

    ci_cl <- 
      as.data.table(
        predict(classical_linear, 
                newdata = data.table(time_stamp = 0:16), 
                interval = 'confidence'))
    ci_cl[, `:=`(statistics = 'classical', 
                 func = 'linear', 
                 time_stamp = 0:16)]

    ci_bl <- 
      as.data.table(
        fitted(bayes_linear,
               newdata = data.table(time_stamp = 0:16)
               ))[, c(1,3,4)]
    ci_bl[, `:=`(statistics = 'bayesian', 
                 func = 'linear', 
                 time_stamp = 0:16)]

    ci_ce <- 
      as.data.table(
        predictNLS(classical_exp, 
                   newdata = data.table(time_stamp = 0:16), 
                   interval = 'confidence')$summary
        )[, c('Prop.Mean.1', 'Sim.2.5%', 'Sim.97.5%')]
    ci_ce[, `:=`(statistics = 'classical', 
                 func = 'exponential', 
                 time_stamp = 0:16)]

    ci_be <- 
      as.data.table(
        fitted(bayes_exp,
               newdata = data.table(time_stamp = 0:16), 
               ))[, c(1,3,4)]
    ci_be[, `:=`(statistics = 'bayesian', 
                 func = 'exponential', 
                 time_stamp = 0:16)]

    label_dt <- 
      make_labels(classical_linear, bayes_linear, 
                  classical_exp, bayes_exp)

    estimates <- rbindlist(list(ci_cl, ci_bl, ci_ce, ci_be))
    return(estimates[label_dt, on = c('statistics', 'func')])
  }

  {
  washout_raw <- fread('carnosine_data.csv')

  washout_0 <- 
    melt(washout_raw, id = c('group_id', 'sample_id'), 
         variable.name = 'week_id')

  week_info <- 
    data.table(week_id= c("PRE_S01", "POS_S02", "W01", 
                          "W02", "W03", "W04", "W05", "W06"),
               time_stamp = c(NA, 0., 1., 2., 4., 8., 12., 16.)) 

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

##########################################
# Fiting full data
##########################################

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
  #40.75576 - 22.02241 = 18.73335

  exp_prior <- 
    prior(normal(18.73335, 12.28948), nlpar = "b1") +
    prior(student_t(3, 0, 13), lb = 0, nlpar = "b2" ) +
    prior(normal(22.02241, 8.636423), nlpar = "b3") 

  bayes_exp <- 
    brm(bf(value ~ b1 * exp(- b2 * time_stamp) + b3,  
           b2 ~ 1, 
           b1 ~ 1,
           b3 ~ 1, 
           nl = TRUE,
           sigma = TRUE, # sigma = TRUE and cmc = TRUE is needed to 
           cmc = TRUE),  # estimate b3 as offset (not in the documentation)
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
    brm(value ~ 1 + (time_stamp || sample_id),
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

  exp_prior_off2 <- 
    prior(normal(18.73335, 12.28948), nlpar = "b1") +
    prior(student_t(3, 0, 13), lb = 0, nlpar = "b2" ) +
    prior(normal(0, 8.636423), nlpar = "b3") 

  bayes_exp_off2 <- 
    brm(bf(value ~ b1 * exp(- b2 * time_stamp) + b3,  
           b2 ~ 1, 
           b1 ~ 1,
           b3 ~ 1, 
           nl = TRUE,
           sigma = TRUE,
           cmc = TRUE),
        data = washout_off[group_id == 'T'& !is.na(time_stamp)], 
        prior = exp_prior,
        control = list(adapt_delta = 0.99))

  summary(bayes_exp_off2)
  }


cl <- summary(classical_linear)$coefficients
ce <- summary(classical_exp)$coefficients
be <- summary(bayes_exp)$fixed
bl <- summary(bayes_linear)$fixed

clo <- summary(classical_linear_off)$coefficients
ceo <- summary(classical_exp_off)$coefficients
beo <- summary(bayes_exp_off)$fixed
blo <- summary(bayes_linear_off)$fixed

summary(bayes_exp)

kable(cl)
kable(ce)
kable(be)
kable(bl)
kable(clo)
kable(ceo)
kable(beo)
kable(blo)
  

plot(bayes_exp)
plot(marginal_effects(bayes_exp), points = TRUE)
loo(bayes_linear, bayes_exp)

plot(bayes_exp_off)
plot(marginal_effects(bayes_exp_off2), points = TRUE)

plot(bayes_linear_off)
plot(marginal_effects(bayes_linear_off), points = TRUE)
loo(bayes_linear_off, bayes_exp_off)


washout_plot <- 
  rbindlist(list(
       washout[, has_offset := TRUE], 
       washout_off[, has_offset := FALSE]))


estimates <- 
  rbindlist(list(
    join_fit(classical_linear, bayes_linear, 
             classical_exp, bayes_exp)[, has_offset := TRUE],
    join_fit(classical_linear_off, bayes_linear_off, 
             classical_exp_off, bayes_exp_off)[, has_offset := FALSE]))


{
ggplot(data = estimates[statistics == 'bayesian' & has_offset == FALSE], 
       aes(x = time_stamp)) + 
  theme_classic() + 
  labs(x = 'time in weeks', y = 'carnosine variation') + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.8, 0.9),
        text = element_text(size = 18),
        legend.key.width=unit(2,"line")
        ) + 
  scale_color_grey() + 
  scale_fill_grey() + 
  scale_linetype_discrete() + 
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = label), alpha = 0.3) + 
  geom_line(aes(y = fit, linetype = label), size = 2, color = 'black') + 
  geom_point(data = washout_plot[group_id == 'T' & has_offset == FALSE], 
             aes(x = time_stamp, y = value),
             size = 2) 
}
