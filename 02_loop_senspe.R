

######################################################################

###
# Get fitted data
##
get_fitted_data <- function()
{
  all <- as.data.frame( readr::read_csv( 'MCIDsimulations_MW/simulation10pt_SF14_P15.csv', show_col_types = FALSE ) )
  all$stratum <- ''
  
  # verschillende strata 
  all[ all$lambda_delta_treatment_perc ==   0, 'stratum' ] <- 'verum 0'
  all[ all$lambda_delta_treatment_perc == -20, 'stratum' ] <- 'verum 20'
  all[ all$lambda_delta_treatment_perc == -30, 'stratum' ] <- 'verum 30'
  all[ all$lambda_delta_treatment_perc == -40, 'stratum' ] <- 'verum 40'
  all[ all$lambda_delta_treatment_perc == -50, 'stratum' ] <- 'verum 50'
  all[ all$lambda_delta_treatment_perc == -60, 'stratum' ] <- 'verum 60'
  all[ all$lambda_delta_treatment_perc == -70, 'stratum' ] <- 'verum 70'
  all[ all$lambda_delta_treatment_perc == -90, 'stratum' ] <- 'verum 90'
  
  all$`...1` <- NULL
  
  return( all )
}

# get fitted data
data <- get_fitted_data()

subset_data <- data %>%
  filter(stratum %in% c("verum 0", "verum 20", "verum 30", "verum 40" , "verum 50", 
                        "verum 60", "verum 70", "verum 90"))

# Zet de 'stratum' om naar een binaire variabele voor de twee groepen
subset_data$stratum_binary <- ifelse(subset_data$stratum %in% c("verum 0"), 0, 
                                     ifelse(subset_data$stratum %in% c( "verum 20", "verum 30", "verum 40" , 
                                                                        "verum 50", "verum 60", "verum 70", "verum 90"), 1, NA))

#########################################################################


# Samenvattingen van percentages per uitkomstmaat van effectiviteit
summary <- subset_data %>%
  group_by(stratum, ncycle) %>%
  summarise(
    totaal = n(),  # Aantal patiënten per groep
    p_significant = (sum(anova_pvalue < 0.05) / totaal) * 100,  # Aantal met p < 0.05
    RR_freq_30 = (sum(RR_treatment < 0.70)/ totaal) * 100,  # Aantal met RR van minimaal 30%
    BF_above1 = (sum(BF > 1)/ totaal) * 100,
    BF_above2 = (sum(BF > 2)/ totaal) * 100,
    BF_above3 = (sum(BF > 3)/ totaal) * 100,
    BF_above10 = (sum(BF > 10)/ totaal) * 100,
    Cert_under20 = (sum(prob_30 < 20)/ totaal) * 100,
    Cert_20to80 = (sum(prob_30 < 80 & prob_30 > 20)/ totaal) * 100,   # aantal met een prob_30% van tussen 20-80% certainty 
    Cert_above80 = (sum(prob_30 > 80)/ totaal) * 100,
    Cert_above30 = (sum(prob_30 > 30)/ totaal) * 100, 
    Cert_above40 = (sum(prob_30 > 40)/ totaal) * 100,
    Cert_above50 = (sum(prob_30 > 50)/ totaal) * 100,
    Cert_above60 = (sum(prob_30 > 60)/ totaal) * 100,
    Cert_above70 = (sum(prob_30 > 70)/ totaal) * 100
  )

print(summary)

write_xlsx(summary, path = "summary_results_MCID_MW_SF14_p15.xlsx")

###########################################################################



######### PPV en NPV per cyclus, VOOR SPECIFIEKE THRESHOLDS ############
library(pROC)

# Prevalentie van de positieve klasse
prevalence <- 0.50  

# Specificeer de gewenste thresholds
desired_thresholds <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# Lijst om resultaten op te slaan per cyclus en per stratum
selectedresults_per_stratum_cycle <- list()

# Loop door elke cyclus
for (ncycle in unique(subset_data$ncycle)) {
  # Filter de data per cyclus
  subset_data_cycle <- subset(subset_data, ncycle == ncycle)
  
  # Loop door de verschillende strata
  for (current_stratum in unique(subset_data_cycle$stratum)) {
    
    # Maak "verum 0" de negatieve groep
    subset_data_stratum <- subset(subset_data_cycle, stratum == current_stratum | stratum == "verum 0")
    
    # Controleer of er voldoende data is (inclusief "verum 0")
    if (nrow(subset_data_stratum) > 1) {
      
      # Zet de stratum om naar een binaire variabele (met "verum 0" als negatieve groep)
      subset_data_stratum$stratum_binary <- ifelse(subset_data_stratum$stratum == "verum 0", 0, 1)
      
      # Verwijder rijen met NA in stratum_binary
      subset_data_stratum <- subset(subset_data_stratum, !is.na(stratum_binary))
      
      # Controleer of beide niveaus (0 en 1) aanwezig zijn
       if (length(unique(subset_data_stratum$stratum_binary)) == 2) {
        
        # Bereken de ROC voor deze stratum
        roc_result <- roc(subset_data_stratum$stratum_binary, subset_data_stratum$BF)
        
        # Interpoleer sensitiviteit en specificiteit voor de gewenste thresholds
        coords_df <- as.data.frame(coords(roc_result, "all", ret = c("threshold", "sensitivity", "specificity")))
        coords_df$threshold <- as.numeric(as.character(coords_df$threshold))
        
        # Interpoleren naar gewenste thresholds
        interpolated_results <- data.frame(
          threshold = desired_thresholds,
          sensitivity = approx(coords_df$threshold, coords_df$sensitivity, xout = desired_thresholds)$y,
          specificity = approx(coords_df$threshold, coords_df$specificity, xout = desired_thresholds)$y
        )
        
        # Bereken PPV en NPV
        interpolated_results$PPV <- (interpolated_results$sensitivity * prevalence) / 
          ((interpolated_results$sensitivity * prevalence) + 
             ((1 - interpolated_results$specificity) * (1 - prevalence)))
        
        interpolated_results$NPV <- (interpolated_results$specificity * (1 - prevalence)) / 
          (((1 - interpolated_results$sensitivity) * prevalence) + 
             (interpolated_results$specificity * (1 - prevalence)))
        
        # Voeg cyclus- en stratum-informatie toe
        interpolated_results$ncycle <- ncycle
        interpolated_results$stratum <- current_stratum
        
        # Voeg toe aan de lijst
        selectedresults_per_stratum_cycle[[paste(ncycle, current_stratum, sep = "_")]] <- interpolated_results
      
      }
    }
  }
} 

# Combineer alle resultaten in één dataframe
selectedthresholds_results <- do.call(rbind, selectedresults_per_stratum_cycle)

# Bekijk de resultaten
print(selectedthresholds_results)

# Exporteer naar Excel
# write_xlsx(selectedfinal_results, path = "thresholds_resultsmarit2_by_stratum.xlsx")

#################################

######### SENSITIVITEIT / SPECIFICITEIT VOOR PROB_30 ############

# Zet de 'stratum' om naar een binaire variabele voor de twee groepen
subset_data$stratum_binary <- ifelse(subset_data$stratum %in% c("verum 0"), 0, 
                                     ifelse(subset_data$stratum %in% c("verum 50"), 1, NA))
# Specificeer de gewenste thresholds
desired_thresholds <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

# Lijst om resultaten op te slaan per cyclus
selectedresults_per_cycle <- list()

# Prevalentie instellen (aanpassen indien nodig)
prevalence <- 0.50

# Loop door elke cyclus
for (ncycle in unique(subset_data$ncycle)) {
  # Filter de data per cyclus
  subset_data_cycle <- subset(subset_data, ncycle == ncycle)
  
  # Controleer of er voldoende data is
  if (nrow(subset_data_cycle) > 1) {
    # Bereken de ROC voor deze cyclus
    roc_result <- roc(subset_data_cycle$stratum_binary, subset_data_cycle$prob_30)
    
    # Interpoleer sensitiviteit en specificiteit voor de gewenste thresholds
    coords_df <- as.data.frame(coords(roc_result, "all", ret = c("threshold", "sensitivity", "specificity")))
    coords_df$threshold <- as.numeric(as.character(coords_df$threshold))
    
    # Interpoleren naar gewenste thresholds
    interpolated_results <- data.frame(
      threshold = desired_thresholds,
      sensitivity = approx(coords_df$threshold, coords_df$sensitivity, xout = desired_thresholds)$y,
      specificity = approx(coords_df$threshold, coords_df$specificity, xout = desired_thresholds)$y
    )
    
    # Bereken PPV en NPV
    interpolated_results$PPV <- (interpolated_results$sensitivity * prevalence) / 
      ((interpolated_results$sensitivity * prevalence) + 
         ((1 - interpolated_results$specificity) * (1 - prevalence)))
    
    interpolated_results$NPV <- (interpolated_results$specificity * (1 - prevalence)) / 
      (((1 - interpolated_results$sensitivity) * prevalence) + 
         (interpolated_results$specificity * (1 - prevalence)))
    
    # Voeg cyclusinformatie toe
    interpolated_results$cycle <- ncycle
    
    # Opslaan in de lijst
    selectedresults_per_cycle[[as.character(ncycle)]] <- interpolated_results
  }
}

# Combineer alle resultaten in één dataframe
selectedfinal_results2 <- do.call(rbind, selectedresults_per_cycle)

# Bekijk de resultaten
print(selectedfinal_results2)
