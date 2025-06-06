# Title: Revisions ----

# Author: Milan Wiedemann
# Status: Accepted

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                 #
# 1      Analysis - Sudden gains and interviewer assessed PTSD Symptoms (PSS-I)   #
# 1.1        Sample 1 - PTSD Symptoms (PSS-I)                                     #
# 1.2        Sample 2 - PTSD Symptoms (PSS-I)                                     #
#                                                                                 #
# 2      Analysis - Correlations at baseline                                      #
# 2.1        Sample 1                                                             #
# 2.1        Sample 2                                                             #
#                                                                                 #
# 3      Analysis - Stable reversals                                              #
# 3.1        Sample 1                                                             #
# 3.1        Sample 2                                                             #
#                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# 1. Analysis - Sudden gains and treatment outcome ----

# Load packages
library(nlme) # Version 3.1-142
library(multcomp) # Version 1.4-11
library(broom) # Version 0.5.3
library(MBESS) # Version 4.6.0

# Define contrast for group (sg vs nosg) differences at end of treatment
# b0 = Intercept
# b1 = Covariate
# b2 = Sudden gains group

contrast_outcome <- rbind("sg_vs_nosg_at_end" = c(0, 0, 1))

# 1.1. Sample 1 - PTSD Symptoms (PSS-I) ----

# Run linear mixed effect model
a1_model_01_sg_pssi_outcome <- lme(pssi_end ~ pssi_s0 + sg,
                                   random = ~ 1 | id,
                                   method = "ML",
                                   na.action = na.omit,
                                   data = data_a1_pssi)

# Show summary of the model
summary(a1_model_01_sg_pssi_outcome)

# Check asssumption of normality of the residuals
qqnorm(resid(a1_model_01_sg_pssi_outcome))

# Calculate contrasts specified in in the contrast matrix "contrast_outcome"
a1_model_01_sg_pssi_outcome_contrasts <- glht(a1_model_01_sg_pssi_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
summary(a1_model_01_sg_pssi_outcome_contrasts, test = adjusted("none"))

# Select estimates, standard error and p-value
table_a1_model_01_sg_pssi_outcome_full <- table_a1_model_01_sg_pssi_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Take absolute value of adjusted difference from the model
a1_pssi_end_diff <- abs(table_a1_model_01_sg_pssi_outcome_full$beta[[1]])

# Calculate pooled standard deviation at baseline
sd_a1_pssi_s0 <- sd(data_a1_pssi$pssi_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a1_pssi_end <- a1_pssi_end_diff / sd_a1_pssi_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a1_pssi_end_95ci <- MBESS::ci.smd(smd = cohens_d_a1_pssi_end, 
                                           n.1 = sum(!is.na(filter(data_a1_pssi, sg == 1)$pssi_end)), 
                                           n.2 = sum(!is.na(filter(data_a1_pssi, sg == 0)$pssi_end)), 
                                           conf.level = .95)

# 1.2. Sample 2 - PTSD Symptoms (PSS-I) ----

# Run linear mixed effect model
a2_model_01_sg_pssi_outcome <- lme(pssi_end ~ pssi_s0 + sg,
                                   random = ~ 1 | id,
                                   method = "ML",
                                   na.action = na.omit,
                                   data = data_a2_pssi)

# Show summary of the model
summary(a2_model_01_sg_pssi_outcome)

# Check asssumption of normality of the residuals
qqnorm(resid(a2_model_01_sg_pssi_outcome))

# Calculate contrasts specified in in the contrast matrix "contrast_outcome"
a2_model_01_sg_pssi_outcome_contrasts <- glht(a2_model_01_sg_pssi_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
summary(a2_model_01_sg_pssi_outcome_contrasts, test = adjusted("none"))

# Select estimates, standard error and p-value
table_a2_model_01_sg_pssi_outcome_full <- table_a2_model_01_sg_pssi_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Take absolute value of adjusted difference from the model
a2_pssi_end_diff <- abs(table_a2_model_01_sg_pssi_outcome_full$beta[[1]])

# Calculate pooled standard deviation at baseline
sd_a2_pssi_s0 <- sd(data_a2_pssi$pssi_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a2_pssi_end <- a2_pssi_end_diff / sd_a2_pssi_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a2_pssi_end_95ci <- MBESS::ci.smd(smd = cohens_d_a2_pssi_end, 
                                           n.1 = sum(!is.na(filter(data_a2_pssi, sg == 1)$pssi_end)), 
                                           n.2 = sum(!is.na(filter(data_a2_pssi, sg == 0)$pssi_end)), 
                                           conf.level = .95)

# 2. Analysis - Correlations at baseline ----

# Load packages
library(tidyverse) # Version 1.3.0
library(psych) # Version 1.8.12

# 2.1. Sample 1 ----
data_a1_byperson %>% 
  select(pds_s0, ptci22r_s0, mem4_s0, bdi_s0, bai_s0) %>% 
  corr.test()

# 2.2. Sample 2 ----
a2_corr_s0 <- data_a2_byperson %>% 
  select(pds_s0, ptci20_s0, mem5_s0, phq9_s0, gad7_s0) %>% 
  corr.test()

# 3. Analysis - Stable reversals ----

# Load packages
library(tidyverse) # Version 1.3.0
library(suddengains) # Version 0.4.0

# For more information see:
# Wucherpfennig, F., Rubel, J. A., Hofmann, S. G., & Lutz, W. (2017a). 
# Processes of change after a sudden gain and relation to treatment outcome - Evidence for an upward spiral. 
# Journal of Consulting and Clinical Psychology, 85(12), 1199–1210. 
# doi:10.1037/ccp0000263

# 3.1. Sample 1 ----

# Filter all patients who met the Tang and DeRubeis criteria for a reversal
a1_sg_reversals <- data_a1_byperson %>% 
  filter(sg_reversal_byperson == TRUE) %>% 
  select(id, sg_session_n, pds_s0:pds_s12)

# Check if reversals also met criteria for a sudden loss
create_bysg(data = a1_sg_reversals, 
            id_var_name = "id", 
            sg_var_list = c("pds_s1", "pds_s2", "pds_s3", "pds_s4", "pds_s5", "pds_s6",
                            "pds_s7", "pds_s8", "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
            tx_start_var_name = "pds_s0",
            tx_end_var_name = "pds_s12",
            sg_measure_name = "pds",
            sg_crit1_cutoff = -6.15, 
            identify = "sl")

# No sudden gainer in Sample 1 experienced a stable reversal.

# 3.2.  Sample 2 ----
# Filter all patients who met the Tang and DeRubeis criteria for a reversal
a2_sg_reversals <- data_a2_byperson %>% 
  filter(sg_reversal_byperson == TRUE) %>% 
  select(id, sg_session_n, pds_s0:pds_s12)

# Check if reversals also met criteria for a sudden loss
create_bysg(data = a2_sg_reversals, 
            id_var_name = "id", 
            sg_var_list = c("pds_s1", "pds_s2", "pds_s3", "pds_s4", "pds_s5", "pds_s6"
                            ,"pds_s7", "pds_s8", "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
            tx_start_var_name = "pds_s0",
            tx_end_var_name = "pds_s12",
            sg_measure_name = "pds",
            sg_crit1_cutoff = -6.15, 
            identify = "sl")

# Three sudden gainers in Sample 2 experienced a stable reversal.
# Following suggestions from a reviewer and the editor we repeated all analyses excluding patients who met the criteria for a stable reversal.
# The results did not differ.
