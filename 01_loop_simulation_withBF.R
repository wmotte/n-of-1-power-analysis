#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
#
# LOOP version (frequentist only) ... 30 January 2025
#
# 1. Simulate single patient
# 2. Model freq.
#
################################################################################

library( 'brms' )
library( 'plyr' ) 
library( 'ggplot2' )
library( 'dplyr' )
library( 'tidyr' )
library('readxl')
library('writexl')

################################################################################
# BEGIN FUNCTIONS
################################################################################

# for plots
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

###
# Simulation
##
simulate_data <- function( patient, ncycles, n_seizures_per_week, lambda_delta_treatment_perc, lambda_delta_placebo_perc )
{
 
   # convert week freq to mean daily seizures
  lambda_baseline <- n_seizures_per_week / 7
  
  # four weeks per block
  days_per_group <- 7 * 4
  
  # get actual lambda reduction, based on given percentage reduction
  lambda_delta_placebo   <- ( lambda_delta_placebo_perc / 100 ) * lambda_baseline
  lambda_delta_treatment <- ( lambda_delta_treatment_perc / 100 ) * lambda_baseline
  
  # first sample baseline (only at beginning of n-of-1 trial
  #
  # baseline, with double number of days (for visualization)
  baseline_data <- rpois( days_per_group * 2, lambda_baseline )
  all <- data.frame( group = 'baseline', cycle = 0, y = baseline_data )
  
  # loop over cycles
  for( cycle in 1:ncycles )
  {
    
    # placebo effect
    placebo_data <- rpois( days_per_group, lambda_baseline + lambda_delta_placebo )
    dfP <- data.frame( group = 'placebo', cycle = cycle, y = placebo_data )
    
    # treatment effect
    treatment_data <- rpois( days_per_group, lambda_baseline + lambda_delta_treatment )
    dfT <- data.frame( group = 'treatment', cycle = cycle, y = treatment_data )
    
    # add to container
    all <- rbind( all, rbind( dfP, dfT ) )
  }
  
  # make factor
  all$group <- as.factor( all$group )    
  
  # add day
  all$day <- 1:nrow( all )
  
  # add for each block
  all$block_id <- paste0( all$group, '_', all$block )
  
  # add patient
  all$patient <- patient
  
  # add meta-info on used parameters
  all$used_lambda_baseline <- lambda_baseline
  all$used_lambda_delta_placebo <-lambda_delta_placebo
  all$used_lambda_delta_treatment <- lambda_delta_treatment
  all$n_seizures_per_week <- n_seizures_per_week
  all$lambda_delta_placebo_perc <- lambda_delta_placebo_perc
  all$lambda_delta_treatment_perc <- lambda_delta_treatment_perc
  
  return( all )
}


################################################################################
# END FUNCTIONS
################################################################################

# output dir
outdir <- 'MCIDsimulations_MW'
dir.create( outdir, showWarnings = FALSE )


####################################################################
########### read and prepare single-subject data ###################
####################################################################

# set seed for reproducibility
set.seed( 1234 )

# init container which stores all output
container <- NULL

# first loop variable
patients <- 1:10

# second loop variable
ncycles <- 1:3

clinically_relevant_reduction_50perc_RR <- 50
clinically_relevant_reduction_30perc_RR <- 30

# third loop variable
lambda_delta_treatment_percs <- c( 0, -20, -30, -40, -50, -60, -70, -90 )

patient <- patients[ 1 ]
ncycle <- ncycles[ 1 ]
lambda_delta_treatment_perc <- lambda_delta_treatment_percs[ 1 ]

# loop over patients
for( patient in patients )
{
  print( paste0( '*** patient ***: ', patient ) )
 
  # loop over cycles
  for( ncycle in ncycles )
  {
    print( ncycle )    
    
    # loop over treatment effects
    for( lambda_delta_treatment_perc in lambda_delta_treatment_percs )
    {
      print( lambda_delta_treatment_perc )
      
      # init
      n_seizures_per_week <- 14
      lambda_delta_placebo_perc <- -15

      # get simulated data, given patient id and number of cycles
      single <- simulate_data( patient = patient, ncycles = ncycle, n_seizures_per_week, lambda_delta_treatment_perc, lambda_delta_placebo_perc )
      
      # make a 'null-group' with 'placebo and treatment' combined, but separate from baseline
      single$group_null <- as.character( single$group )
      single[ single$group_null %in% c( 'treatment', 'placebo' ), 'group_null' ] <- 'post_baseline'

            # frequentist modellen 
      model <- glm( y ~ group, family = poisson( link = "log" ), data = single )
      
      # fit a null-model
      model_null <- glm( y ~ group_null, family = poisson( link = "log" ), data = single )
      
      # check with the likelihood ratio which model fits the data best [Analysis of Deviance]
      anova_res <- anova( model_null, model, test = "Chisq" )
      
      ################################################
      ###### MCID calculation for freq. model ########
      ##
      ## Percentage Reduction = ( 1 − IRR ) × 100 
      ## with IRR = Incidence Rate Ratio
      
      # get the IRRs and confidence intervals
      # exponentiate the coefficients to get relative risks [i.e., 'groupplacebo' and 'grouptreatment']
      params <- parameters::parameters( model, exponentiate = TRUE )

      # get relative risk for the "treatment" group
      params_treatment <- params %>% filter( Parameter == "grouptreatment" )
      RR_treatment <- params_treatment$Coefficient
      RR_treatment_low <- params_treatment$CI_low
      RR_treatment_high <- params_treatment$CI_high
      
      # calculate reduction percentage (1 - relative risk) * 100
      reduction_percentage <- round( ( 1 - RR_treatment ) * 100, 1 )
      
      # get CI values as well
      reduction_percentage_CI_high <- round( ( 1 - RR_treatment_low ) * 100, 1 )
      reduction_percentage_CI_low <- round( ( 1 - RR_treatment_high ) * 100, 1 )
      
      # check if reduction in treatment group is bigger than 50%/30%, based on median
      mcid_reached_50 <- reduction_percentage > 50
      mcid_reached_30 <- reduction_percentage > 30
      
      # fit Bayesian Poisson model (baseline, placebo, treatment)
      # refresh = 0 turns off iteration info print
      bmodel_H1 <- brms::brm( y ~ group, data = single,
                              iter = 25000, warmup = 5000,
                              family = poisson(),
                              thin = 5,
                              chains = 4,
                              save_pars = save_pars( all = TRUE ), refresh = 0 )
      
      # update Bayesian Poisson model [H0], now with null-group (i.e., baseline versus post-baseline)
      bmodel_H0 <- update( bmodel_H1, formula = y ~ group_null, newdata = single, recompile = TRUE, refresh = 0 )
      
      # Bayes Factor, through bridge sampling. H1 is clearly higher in evidence (log(BF) == 6.2)
      BF <- brms::bayes_factor( bmodel_H1, bmodel_H0, log = TRUE )
      BF_value <- as.numeric(BF$bf)
      
      #print( paste0( "*** log( BF_H1/H0) => ", round( as.numeric( BF$bf ), 1 ) ) )
      
      
      # Parameter      | Median |        95% CI |     pd |  Rhat |      ESS
      # -------------------------------------------------------------------
      # (Intercept)    |   0.76 | [0.55, 1.01] | 97.08% | 1.000 | 15378.00
      # groupplacebo   |   1.19 | [0.82, 1.74] | 81.41% | 1.000 | 15580.00
      # grouptreatment |   0.57 | [0.37, 0.89] | 99.33% | 1.000 | 15728.00
      parameters::parameters( bmodel_H1, exponentiate = TRUE )
      
      # Parameter               | Median |        95% CI |   pd |  Rhat |      ESS
      # --------------------------------------------------------------------------
      # (Intercept)             |   0.76 | [0.55, 1.01] | 97.06% | 1.000 | 16432.00
      # group_nullpost_baseline |   0.88 | [0.63, 1.26] | 75.59% | 1.000 | 15698.00
      parameters::parameters( bmodel_H0, exponentiate = TRUE )
      
      
      # posterior post-warm-up samples [16_000 samples]   
      draws <- as.data.frame( bmodel_H1 )
      
      # get risk-ratios for treatment group + 95% Credibility Interval (CrI)
      bayes_RR_treatment <- median( exp( draws$b_grouptreatment ) )
      bayes_RR_treatment_low <- quantile( exp( draws$b_grouptreatment ), 0.025 )
      bayes_RR_treatment_high <- quantile( exp( draws$b_grouptreatment ), 0.975 )
      
      bayes_reduction_percentage <- round( ( 1 - bayes_RR_treatment ) * 100, 1 )
      
      # get CI values for reduction as well
      bayes_reduction_percentage_CI_low <- round( ( 1 - bayes_RR_treatment_high ) * 100, 1 )
      bayes_reduction_percentage_CI_high <- round( ( 1 - bayes_RR_treatment_low ) * 100, 1 )
      
         # check if reduction in treatment group is bigger than 50/ 30%, based on median
      bayes_mcid_reached_50 <- bayes_reduction_percentage > 50
      bayes_mcid_reached_30 <- bayes_reduction_percentage > 30
      
      pdata <- data.frame( iter = 1:nrow( draws ), reduction_perc = ( 1 - exp( draws$b_grouptreatment ) ) * 100 )
      
      # booleans
      pdata$passed_50 <- pdata$reduction_perc > clinically_relevant_reduction_50perc_RR
      pdata$passed_30 <- pdata$reduction_perc > clinically_relevant_reduction_30perc_RR
      
      # probabilities
      prob_50 <- round( 100 * ( sum( pdata$passed_50 ) / nrow( pdata ) ), 1 )
      prob_30 <- round( 100 * ( sum( pdata$passed_30 ) / nrow( pdata ) ), 1 )   
      
    
      # make output (single row)
      output <- data.frame( patient = patient, 
                            ncycle = ncycle,
                            lambda_delta_treatment_perc = lambda_delta_treatment_perc,
                            anova_deviance = anova_res$Deviance[ 2 ],
                            anova_pvalue = anova_res$`Pr(>Chi)`[ 2 ],
                            RR_treatment = RR_treatment,
                            RR_treatment_low = RR_treatment_low,
                            RR_treatment_high = RR_treatment_high,
                            mcid_reached_30 = mcid_reached_30, 
                            mcid_reached_50 = mcid_reached_50,
                            BF = BF_value,
                            bayes_RR_treatment = bayes_RR_treatment,
                            bayes_RR_treatment_low = bayes_RR_treatment_low,
                            bayes_RR_treatment_high = bayes_RR_treatment_high,
                            bayes_mcid_reached_30 = bayes_mcid_reached_30, 
                            bayes_mcid_reached_50 = bayes_mcid_reached_50,
                            prob_30 = prob_30,
                            prob_50 = prob_50 )
      
      # add to container
      container <- rbind( container, output )
      
    } # closes lambda_delta_treatment loop
  } # closes ncycle loop
  
} # closes patient loop

# write to file
readr::write_csv( container, file = paste0( outdir, '/simulation10pt_SF14_P15.csv' ), quote = 'all' )
