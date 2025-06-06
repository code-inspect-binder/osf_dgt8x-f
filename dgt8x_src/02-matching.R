# Title: Matching ----

# Author: Milan Wiedemann
# Status: Accepted

# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
# 1.    Matching                                  #
# 1.1.1     Sample 1: Calculate propensity score  #
# 1.2.1     Sample 1: Matching                    #
# 1.1.2     Sample 2: Calculate propensity score  #
# 1.2.2     Sample 2: Matching                    #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # 

# 1. Matching ----

# Load packages
library(tidyverse)
library(MatchIt)

# Load data ----
data_a1_byperson_matchit <- data_a1_byperson %>% 
  select(id, sg_crit123, age, months_since_trauma, 
         pds_s0, ptci22r_s0, bdi_s0, bai_s0, mem4_s0, 
         traumatype_index, sex, comorbid_depression)

data_a2_byperson_matchit <- data_a2_byperson %>% 
  select(id, sg_crit123, age, months_since_trauma, 
         pds_s0, ptci20_s0, phq9_s0, gad7_s0, mem5_s0, 
         traumatype_index, sex, comorbid_depression)

# 1.1.1 Sample 1: Calculate propensity score ----
a1_psmodel <- glm(sg_crit123 ~ age + sex + months_since_trauma + traumatype_index + comorbid_depression +
                               pds_s0 + bdi_s0 + bai_s0 + ptci22r_s0 + mem4_s0,
                family = binomial(),
                data = data_a1_byperson_matchit)

# Show summary of the model
summary(a1_psmodel)

# Extract fitted propensity scores and add to dataset
data_a1_byperson_matchit$pscore <- fitted(a1_psmodel)

# 1.2.1 Sample 1: Matching ----
a1_matchit_out <- matchit(sg_crit123 ~ age + sex + months_since_trauma + traumatype_index + comorbid_depression +
                                       pds_s0 + bdi_s0 + bai_s0 + ptci22r_s0 + mem4_s0 + pscore,
                          mahvars = c("age", "sex", "months_since_trauma", "traumatype_index", "comorbid_depression", 
                                      "pds_s0", "bdi_s0", "bai_s0", "ptci22r_s0", "mem4_s0", "pscore"),
                          caliper = .25,
                          distance = "mahalanobis",
                          method = "nearest",
                          replace = FALSE, 
                          ratio = 1,
                          data = data_a1_byperson_matchit)

# Show summary of match
summary(a1_matchit_out)

# 1.1.2 Sample 2: Calculate propensity score ----
a2_psmodel <- glm(sg_crit123 ~ age + sex + months_since_trauma +  traumatype_index + comorbid_depression +
                               pds_s0 + phq9_s0 + gad7_s0 + ptci20_s0 + mem5_s0,
                    family = binomial(),
                    data = data_a2_byperson_matchit)

# Show summary of the model
# summary(a2_psmodel)

# Extract fitted propensity scores and add to dataset
data_a2_byperson_matchit$pscore <- fitted(a2_psmodel)

# 1.2.2 Sample 2: Matching ----
a2_matchit_out <- matchit(sg_crit123 ~ age + sex + months_since_trauma + traumatype_index + comorbid_depression + 
                                       pds_s0 + phq9_s0 + gad7_s0 + ptci20_s0 + mem5_s0 + pscore,
                          mahvars = c("age", "sex", "months_since_trauma", "traumatype_index", "comorbid_depression", 
                                      "pds_s0", "phq9_s0", "gad7_s0", "ptci20_s0", "mem5_s0", "pscore"),
                          caliper = .25,
                          distance = "mahalanobis",
                          method = "nearest",
                          replace = FALSE,
                          ratio = 1,
                          data = data_a2_byperson_matchit)

# Show summary of match
summary(a2_matchit_out)