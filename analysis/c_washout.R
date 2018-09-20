library(data.table)
library(ggplot2)
library(brms)
library(propagate)
library(knitr)

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
  #40.75576 - 22.02241 = 18.73335

  exp_prior <- 
    prior(normal(18.73335, 12.28948), nlpar = "b1") +
    #prior(normal(1, 5), lb = 0, nlpar = "b2" ) +
    prior(student_t(3, 0, 13), lb = 0, nlpar = "b2" ) +
    prior(normal(22.02241, 8.636423), nlpar = "b3") 

  bayes_exp <- 
    brm(bf(value  ~ b1 * exp(- b2 * time_stamp) + b3 + 0,  
           b1 + b2 + b3 ~ 1, 
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

estimates <- 
  join_fit(classical_linear, bayes_linear, 
           classical_exp, bayes_exp)

estimates_off <- 
  join_fit(classical_linear_off, bayes_linear_off, 
           classical_exp_off, bayes_exp_off)


ggplot(data = estimates_off[!grepl('classical', fit_type)], aes(x = time_stamp)) + 
geom_line(aes(y = fit, color = fit_type), size = 2) + 
geom_ribbon(aes(ymax = upr, ymin = lwr, fill = fit_type), alpha = 0.1) + 
geom_point(data = washout_off[group_id == 'T'], aes(x = time_stamp, y = value)) + 
theme_bw()


cl <- summary(classical_linear)
ce <- summary(classical_exp)
be <- summary(bayes_exp)
bl <- summary(bayes_linear)

clo <- summary(classical_linear_off)
ceo <- summary(classical_exp_off)
beo <- summary(bayes_exp_off)
blo <- summary(bayes_linear_off)

plot(bayes_exp)
plot(marginal_effects(bayes_exp), points = TRUE)
loo(bayes_linear, bayes_exp)

plot(bayes_exp_off)
plot(marginal_effects(bayes_exp_off), points = TRUE)

plot(bayes_linear_off)
plot(marginal_effects(bayes_linear_off), points = TRUE)
loo(bayes_linear_off, bayes_exp_off)

b1_e <- format(beo$fixed["b1_Intercept","Estimate"], digits = 3)
b2_e <- format(beo$fixed["b2_Intercept","Estimate"], digits = 2)
b1_l <- format(blo$fixed["Intercept","Estimate"], digits = 3)
b2_l <- format(abs(blo$fixed["time_stamp","Estimate"]), digits = 3)
le <- paste0(b1_e,' Exp(-', b2_e, ' t)')
ll <- paste0(b1_l, ' - ',  b2_l, ' t')

{
ggplot(data = estimates_off[!grepl('classical', fit_type)], aes(x = time_stamp)) + 
  theme_classic() + 
  labs(x = 'time in weeks', y = 'carnosine variation') + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.8, 0.9),
        text = element_text(size = 18),
        legend.key.width=unit(2,"line")
        ) + 
  scale_color_grey(labels = c(le,ll)) + 
  scale_fill_grey(labels = c(le,ll)) + 
  scale_linetype_discrete(labels = c(le,ll)) + 
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = fit_type), alpha = 0.3) + 
  geom_line(aes(y = fit, linetype = fit_type), size = 2, color = 'black') + 
  geom_point(data = washout_off[group_id == 'T'], 
             aes(x = time_stamp, y = value),
             size = 2) 
}
