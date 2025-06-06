# Title: Analysis scripts for SG in CT-PTSD ----

# Author: Milan Wiedemann
# Status: Accepted

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                 #
# 1      Analysis - Sudden gains and treatment outcome                            #
# 1.1.1      Sample 1 - PTSD Symptoms                                             #
# 1.1.2      Sample 2 - PTSD Symptoms                                             #
# 1.2.1      Sample 1 - Depression Symptoms                                       #
# 1.2.2      Sample 2 - Depression Symptoms                                       #
# 1.3.1      Sample 1 - Anxiety Symptoms                                          #
# 1.3.2      Sample 2 - Anxiety Symptoms                                          #
#                                                                                 #
# 2      Analysis - Baseline predictors                                           #
# 2.1.1      Analysis (Sample 1) - Baseline predictors (uni)                      #
# 2.2.2      Analysis (Sample 1) - Baseline predictors (multi)                    #
# 2.3.1      Analysis (Sample 1) - Baseline predictors (multi using sig from uni) #
# 2.1.2      Analysis (Sample 2) - Baseline predictors (uni)                      #
# 2.2.1      Analysis (Sample 2) - Baseline predictors (multi)                    #
#                                                                                 #
# 3      Analysis - Changes around sudden gains                                   #
# 3.1.1      Analysis (Sample 1) - Negative appraisals                            #
# 3.1.2      Analysis (Sample 2) - Negative appraisals                            #
# 3.1.3      Analysis (Samples 1 & 2) - Negative appraisals, Pooled effects       #
#                                                                                 #
# 3.2.1      Analysis (Sample 1) - Memory characteristics                         #
# 3.2.2      Analysis (Sample 2) - Memory characteristics                         #
# 3.2.3      Analysis (Samples 1 & 2) - Memory characteristics, Pooled effects    #
#                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# 1. Analysis: Sudden gains and treatment outcome ----

# Define contrasts for group (sg vs nosg) differences at: (1) end of treatment and (2) follow-up
# b0 = Intercept
# b1 = Covariate
# b2 = Time (0 = end, 1 = fu)
# b3 = Sudden gains group
# b4 = time for sg group (0 = end, 1 = fu)

contrast_outcome <- rbind("sg_vs_nosg_at_end" = c(0, 0, 0, 1, 0),
                          "sg_vs_nosg_at_fu" = c(0, 0, 0, 1, 1))

# 1.1.1 Sample 1 - PTSD Symptoms ----

# Reshape data from wide to long format
data_a1_sg_pds_outcome <- data_a1_byperson %>% 
  select(id, 
         sg = sg_crit123, 
         pds_s0, pds_end, pds_fu) %>% 
  gather(c("pds_end", "pds_fu"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "pds_end", 0),
         time = replace(time, time == "pds_fu", 1),
         time_fct = factor(time, levels = c("0", "1")))

# Run linear mixed effect model
a1_model_01_sg_pds_outcome <- lme(value ~ pds_s0 + time_fct + sg + time_fct * sg,
                                  random = ~ 1 | id,
                                  method = "ML",
                                  na.action = na.omit,
                                  data = data_a1_sg_pds_outcome)

# Show summary of the model
# summary(a1_model_01_sg_pds_outcome)

# Check assumption of normality of the residuals
qqnorm(resid(a1_model_01_sg_pds_outcome))

# Calculate contrasts specified in in the contrast matrix "contrast_outcome"
a1_model_01_sg_pds_outcome_contrasts <- glht(a1_model_01_sg_pds_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a1_model_01_sg_pds_outcome_contrasts, test = adjusted("none"))

# Extract model parameters and assign to a new dataframe
table_a1_model_01_sg_pds_outcome_raw <- tidy(summary(a1_model_01_sg_pds_outcome_contrasts, test = adjusted("none")))

# Select estimates, standard error and p-value
table_a1_model_01_sg_pds_outcome_full <- table_a1_model_01_sg_pds_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Calculate Cohen's d
# Take absolute value of adjusted difference from the model
a1_pds_end_diff <- abs(table_a1_model_01_sg_pds_outcome_full$beta[[1]])
a1_pds_fu_diff <- abs(table_a1_model_01_sg_pds_outcome_full$beta[[2]])

# Calculate the (shared) standard deviation at baseline of the total sample
sd_a1_pds_s0 <- sd(data_a1_sg_pds_outcome$pds_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a1_pds_end <- a1_pds_end_diff / sd_a1_pds_s0
cohens_d_a1_pds_fu <- a1_pds_fu_diff / sd_a1_pds_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a1_pds_end_95ci <- MBESS::ci.smd(smd = cohens_d_a1_pds_end, 
                                          n.1 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 1)$pds_end)), 
                                          n.2 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 0)$pds_end)), 
                                          conf.level = .95)

cohens_d_a1_pds_fu_95ci <- MBESS::ci.smd(smd = cohens_d_a1_pds_fu, 
                                         n.1 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 1)$pds_fu)), 
                                         n.2 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 0)$pds_fu)), 
                                         conf.level = .95)

# 1.1.2 Sample 2 - PTSD Symptoms ----

# Reshape data from wide to long format
data_a2_sg_pds_outcome <- data_a2_byperson %>% 
  select(id, 
         sg = sg_crit123, 
         pds_s0, pds_end, pds_fu) %>% 
  gather(c("pds_end", "pds_fu"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "pds_end", 0),
         time = replace(time, time == "pds_fu", 1),
         time_fct = factor(time, levels = c("0", "1")))

# Run linear mixed effect model
a2_model_01_sg_pds_outcome <- lme(value ~ pds_s0 + time_fct + sg + time_fct * sg,
                                  random = ~ 1 | id,
                                  method = "ML",
                                  na.action = na.omit,
                                  data = data_a2_sg_pds_outcome)

# Show summary of the model
# summary(a2_model_01_sg_pds_outcome)

# Check assumption of normality of the residuals
qqnorm(resid(a2_model_01_sg_pds_outcome))

# Calculate contrasts specified in in the contrast matrix "contrast_outcome"
a2_model_01_sg_pds_outcome_contrasts <- glht(a2_model_01_sg_pds_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a2_model_01_sg_pds_outcome_contrasts, test = adjusted("none"))

# Extract model parameters and assign to a new dataframe
table_a2_model_01_sg_pds_outcome_raw <- tidy(summary(a2_model_01_sg_pds_outcome_contrasts, test = adjusted("none")))

# Select estimates, standard error and p-value
table_a2_model_01_sg_pds_outcome_full <- table_a2_model_01_sg_pds_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Calculate Cohen's d
# Take absolute value of adjusted difference from the model
a2_pds_end_diff <- abs(table_a2_model_01_sg_pds_outcome_full$beta[[1]])
a2_pds_fu_diff <- abs(table_a2_model_01_sg_pds_outcome_full$beta[[2]])

# Calculate the (shared) standard deviation at baseline of the total sample
sd_a2_pds_s0 <- sd(data_a2_sg_pds_outcome$pds_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a2_pds_end <- a2_pds_end_diff / sd_a2_pds_s0
cohens_d_a2_pds_fu <- a2_pds_fu_diff / sd_a2_pds_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a2_pds_end_95ci <- MBESS::ci.smd(smd = cohens_d_a2_pds_end, 
                                          n.1 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 1)$pds_end)), 
                                          n.2 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 0)$pds_end)), 
                                          conf.level = .95)

cohens_d_a2_pds_fu_95ci <- MBESS::ci.smd(smd = cohens_d_a2_pds_fu, 
                                         n.1 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 1)$pds_fu)), 
                                         n.2 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 0)$pds_fu)), 
                                         conf.level = .95)

# 1.2.1 Sample 1 - Depression symptoms ----

# Reshape data from wide to long format
data_a1_sg_bdi_outcome <- data_a1_byperson %>% 
  select(id, 
         sg = sg_crit123, 
         bdi_s0, bdi_end, bdi_fu) %>% 
  gather(c("bdi_end", "bdi_fu"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "bdi_end", 0),
         time = replace(time, time == "bdi_fu", 1),
         time_fct = factor(time, levels = c("0", "1")))

# Run linear mixed effect model
a1_model_01_sg_bdi_outcome <- lme(value ~ bdi_s0 + time_fct + sg + time_fct * sg,
                                  random = ~ 1 | id,
                                  method = "ML",
                                  na.action = na.omit,
                                  data = data_a1_sg_bdi_outcome)

# Show summary of the model
# summary(a1_model_01_sg_bdi_outcome)

# Check assumption of normality of the residuals
qqnorm(resid(a1_model_01_sg_bdi_outcome))

# Calculate contrasts specified in the matrix from the model
a1_model_01_sg_bdi_outcome_contrasts <- glht(a1_model_01_sg_bdi_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a1_model_01_sg_bdi_outcome_contrasts, test = adjusted("none"))

# Extract model parameters and assign to a new dataframe
table_a1_model_01_sg_bdi_outcome_raw <- tidy(summary(a1_model_01_sg_bdi_outcome_contrasts, test = adjusted("none")))

# Select only variables needed for the table
table_a1_model_01_sg_bdi_outcome_full <- table_a1_model_01_sg_bdi_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Calculate Cohen's d
# Take absolute value of adjusted difference from the model
a1_bdi_end_diff <- abs(table_a1_model_01_sg_bdi_outcome_full$beta[[1]])
a1_bdi_fu_diff <- abs(table_a1_model_01_sg_bdi_outcome_full$beta[[2]])

# Calculate shared standard deviation at baseline, see method described in @Freeman2017
sd_a1_bdi_s0 <- sd(data_a1_sg_bdi_outcome$bdi_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a1_bdi_end <- a1_bdi_end_diff / sd_a1_bdi_s0
cohens_d_a1_bdi_fu <- a1_bdi_fu_diff / sd_a1_bdi_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a1_bdi_end_95ci <- MBESS::ci.smd(smd = cohens_d_a1_bdi_end, 
                                          n.1 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 1)$bdi_end)), 
                                          n.2 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 0)$bdi_end)), 
                                          conf.level = .95)

cohens_d_a1_bdi_fu_95ci <- MBESS::ci.smd(smd = cohens_d_a1_bdi_fu, 
                                         n.1 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 1)$bdi_fu)), 
                                         n.2 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 0)$bdi_fu)), 
                                         conf.level = .95)

# 1.2.2 Sample 2 - Depression symptoms ----

# Reshape data from wide to long format
data_a2_sg_phq9_outcome <- data_a2_byperson %>% 
  select(id, 
         sg = sg_crit123, 
         phq9_s0, phq9_end, phq9_fu) %>% 
  gather(c("phq9_end", "phq9_fu"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "phq9_end", 0),
         time = replace(time, time == "phq9_fu", 1),
         time_fct = factor(time, levels = c("0", "1")))

# Run linear mixed effect model
a2_model_01_sg_phq9_outcome <- lme(value ~ phq9_s0 + time_fct + sg + time_fct * sg,
                                   random = ~ 1 | id,
                                   method = "ML",
                                   na.action = na.omit,
                                   data = data_a2_sg_phq9_outcome)

# Show summary of the model
# summary(a2_model_01_sg_phq9_outcome)

# Check assumption of normality of the residuals
qqnorm(resid(a2_model_01_sg_phq9_outcome))

# Calculate contrasts specified in the matrix from the model
a2_model_01_sg_phq9_outcome_contrasts <- glht(a2_model_01_sg_phq9_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a2_model_01_sg_phq9_outcome_contrasts, test = adjusted("none"))

# Extract model parameters and assign to a new dataframe
table_a2_model_01_sg_phq9_outcome_raw <- tidy(summary(a2_model_01_sg_phq9_outcome_contrasts, test = adjusted("none")))

table_a2_model_01_sg_phq9_outcome_full <- table_a2_model_01_sg_phq9_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Calculate Cohen's d
# Take absolute value of adjusted difference from the model
a2_phq9_end_diff <- abs(table_a2_model_01_sg_phq9_outcome_full$beta[[1]])
a2_phq9_fu_diff <- abs(table_a2_model_01_sg_phq9_outcome_full$beta[[2]])

# Calculate the (shared) standard deviation at baseline of the total sample
sd_a2_phq9_s0 <- sd(data_a2_sg_phq9_outcome$phq9_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a2_phq9_end <- a2_phq9_end_diff / sd_a2_phq9_s0
cohens_d_a2_phq9_fu <- a2_phq9_fu_diff / sd_a2_phq9_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a2_phq9_end_95ci <- MBESS::ci.smd(smd = cohens_d_a2_phq9_end, 
                                           n.1 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 1)$phq9_end)), 
                                           n.2 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 0)$phq9_end)), 
                                           conf.level = .95)

cohens_d_a2_phq9_fu_95ci <- MBESS::ci.smd(smd = cohens_d_a2_phq9_fu, 
                                          n.1 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 1)$phq9_fu)), 
                                          n.2 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 0)$phq9_fu)), 
                                          conf.level = .95)

# 1.3.1 Sample 1 - Anxiety symptoms ----

# Reshape data from wide to long format
data_a1_sg_bai_outcome <- data_a1_byperson %>% 
  select(id, 
         sg = sg_crit123, 
         bai_s0, bai_end, bai_fu) %>% 
  gather(c("bai_end", "bai_fu"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "bai_end", 0),
         time = replace(time, time == "bai_fu", 1),
         time_fct = factor(time, levels = c("0", "1")))

# Run linear mixed effect model
a1_model_01_sg_bai_outcome <- lme(value ~ bai_s0 + time_fct + sg + time_fct * sg,
                                  random = ~ 1 | id,
                                  method = "ML",
                                  na.action = na.omit,
                                  data = data_a1_sg_bai_outcome)

# Show summary of the model
# summary(a1_model_01_sg_bai_outcome)

# Check assumption of normality of the residuals
qqnorm(resid(a1_model_01_sg_bai_outcome))

# Calculate contrasts specified in the matrix from the model
a1_model_01_sg_bai_outcome_contrasts <- glht(a1_model_01_sg_bai_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a1_model_01_sg_bai_outcome_contrasts, test = adjusted("none"))

# Extract model parameters and assign to a new dataframe
table_a1_model_01_sg_bai_outcome_raw <- tidy(summary(a1_model_01_sg_bai_outcome_contrasts, test = adjusted("none")))

# Select only variables needed for the table
table_a1_model_01_sg_bai_outcome_full <- table_a1_model_01_sg_bai_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Calculate Cohen's d
# Take absolute value of adjusted difference from the model
a1_bai_end_diff <- abs(table_a1_model_01_sg_bai_outcome_full$beta[[1]])
a1_bai_fu_diff <- abs(table_a1_model_01_sg_bai_outcome_full$beta[[2]])

# Calculate the (shared) standard deviation at baseline of the total sample
sd_a1_bai_s0 <- sd(data_a1_sg_bai_outcome$bai_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a1_bai_end <- a1_bai_end_diff / sd_a1_bai_s0
cohens_d_a1_bai_fu <- a1_bai_fu_diff / sd_a1_bai_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a1_bai_end_95ci <- MBESS::ci.smd(smd = cohens_d_a1_bai_end, 
                                          n.1 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 1)$bai_end)), 
                                          n.2 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 0)$bai_end)), 
                                          conf.level = .95)

cohens_d_a1_bai_fu_95ci <- MBESS::ci.smd(smd = cohens_d_a1_bai_fu, 
                                         n.1 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 1)$bai_fu)), 
                                         n.2 = sum(!is.na(filter(data_a1_byperson, sg_crit123 == 0)$bai_fu)), 
                                         conf.level = .95)

# 1.3.2 Sample 2 - Anxiety symptoms ----

# Reshape data from wide to long format
data_a2_sg_gad7_outcome <- data_a2_byperson %>% 
  select(id, 
         sg = sg_crit123, 
         gad7_s0, gad7_end, gad7_fu) %>% 
  gather(c("gad7_end", "gad7_fu"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "gad7_end", 0),
         time = replace(time, time == "gad7_fu", 1),
         time_fct = factor(time, levels = c("0", "1")))

# Run linear mixed effect model
a2_model_01_sg_gad7_outcome <- lme(value ~ gad7_s0 + time_fct + sg + time_fct * sg,
                                   random = ~ 1 | id,
                                   method = "ML",
                                   na.action = na.omit,
                                   data = data_a2_sg_gad7_outcome)

# Show summary of the model
# summary(a2_model_01_sg_gad7_outcome)

# Check assumption of normality of the residuals
qqnorm(resid(a2_model_01_sg_gad7_outcome))

# Calculate contrasts specified in the matrix from the model
a2_model_01_sg_gad7_outcome_contrasts <- glht(a2_model_01_sg_gad7_outcome, contrast_outcome)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a2_model_01_sg_gad7_outcome_contrasts, test = adjusted("none"))

# Extract model parameters and assign to a new dataframe
table_a2_model_01_sg_gad7_outcome_raw <- tidy(summary(a2_model_01_sg_gad7_outcome_contrasts, test = adjusted("none")))

table_a2_model_01_sg_gad7_outcome_full <- table_a2_model_01_sg_gad7_outcome_raw %>% 
  dplyr::select(beta = estimate, se = std.error, p = p.value)

# Calculate Cohen's d
# Take absolute value of adjusted difference from the model
a2_gad7_end_diff <- abs(table_a2_model_01_sg_gad7_outcome_full$beta[[1]])
a2_gad7_fu_diff <- abs(table_a2_model_01_sg_gad7_outcome_full$beta[[2]])

# Calculate the (shared) standard deviation at baseline of the total sample
sd_a2_gad7_s0 <- sd(data_a2_sg_gad7_outcome$gad7_s0, na.rm = TRUE)

# Calculate Cohen's d using the absolute value of the baseline adjusted difference from the linear mixed effect model
cohens_d_a2_gad7_end <- a2_gad7_end_diff / sd_a2_gad7_s0
cohens_d_a2_gad7_fu <- a2_gad7_fu_diff / sd_a2_gad7_s0

# Calculate confidence intervals of Cohen's d
cohens_d_a2_gad7_end_95ci <- MBESS::ci.smd(smd = cohens_d_a2_gad7_end, 
                                           n.1 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 1)$gad7_end)), 
                                           n.2 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 0)$gad7_end)), 
                                           conf.level = .95)

cohens_d_a2_gad7_fu_95ci <- MBESS::ci.smd(smd = cohens_d_a2_gad7_fu, 
                                          n.1 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 1)$gad7_fu)), 
                                          n.2 = sum(!is.na(filter(data_a2_byperson, sg_crit123 == 0)$gad7_fu)), 
                                          conf.level = .95)

# 2. Analysis: Baseline predictors ----

# 2.1.1 Sample 1: Univariate logistic regressions ----

# Sample 1: Age
a1_sg_pred_model_age <- glm(sg_crit123 ~  age, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_age)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_age <- exp(coef(a1_sg_pred_model_age))[[2]]
a1_LL_age <- exp(cbind(OR = coef(a1_sg_pred_model_age), confint(a1_sg_pred_model_age)))[ , 2][[2]]
a1_UL_age <- exp(cbind(OR = coef(a1_sg_pred_model_age), confint(a1_sg_pred_model_age)))[ , 3][[2]]

# Extract p value from model
a1_p_age <- tidy(a1_sg_pred_model_age)[[2 , 5]]

# Sample 1: Gender
a1_sg_pred_model_sex <- glm(sg_crit123 ~  sex, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_sex)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_sex <- exp(coef(a1_sg_pred_model_sex))[[2]]
a1_LL_sex <- exp(cbind(OR = coef(a1_sg_pred_model_sex), confint(a1_sg_pred_model_sex)))[ , 2][[2]]
a1_UL_sex <- exp(cbind(OR = coef(a1_sg_pred_model_sex), confint(a1_sg_pred_model_sex)))[ , 3][[2]]

# Extract p value from model
a1_p_sex <- tidy(a1_sg_pred_model_sex)[[2 , 5]]

# Sample 1: Months since trauma
a1_sg_pred_model_months_since_trauma <- glm(sg_crit123 ~  months_since_trauma, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_months_since_trauma)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_months_since_trauma <- exp(coef(a1_sg_pred_model_months_since_trauma))[[2]]
a1_LL_months_since_trauma <- exp(cbind(OR = coef(a1_sg_pred_model_months_since_trauma), confint(a1_sg_pred_model_months_since_trauma)))[ , 2][[2]]
a1_UL_months_since_trauma <- exp(cbind(OR = coef(a1_sg_pred_model_months_since_trauma), confint(a1_sg_pred_model_months_since_trauma)))[ , 3][[2]]

# Extract p value from model
a1_p_months_since_trauma <- tidy(a1_sg_pred_model_months_since_trauma)[[2 , 5]]

# Sample 1: PTSD symptoms at baseline
a1_sg_pred_model_pds_s0 <- glm(sg_crit123 ~  pds_s0, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_pds_s0)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_pds_s0 <- exp(coef(a1_sg_pred_model_pds_s0))[[2]]
a1_LL_pds_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_pds_s0), confint(a1_sg_pred_model_pds_s0)))[ , 2][[2]]
a1_UL_pds_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_pds_s0), confint(a1_sg_pred_model_pds_s0)))[ , 3][[2]]

# Extract p value from model
a1_p_pds_s0 <- tidy(a1_sg_pred_model_pds_s0)[[2 , 5]]

# Sample 1: PTSD symptoms at baseline
a1_sg_pred_model_bdi_s0 <- glm(sg_crit123 ~  bdi_s0, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_bdi_s0)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_bdi_s0 <- exp(coef(a1_sg_pred_model_bdi_s0))[[2]]
a1_LL_bdi_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_bdi_s0), confint(a1_sg_pred_model_bdi_s0)))[ , 2][[2]]
a1_UL_bdi_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_bdi_s0), confint(a1_sg_pred_model_bdi_s0)))[ , 3][[2]]

# Extract p value from model
a1_p_bdi_s0 <- tidy(a1_sg_pred_model_bdi_s0)[[2 , 5]]

# Sample 1: Anxiety symptoms at baseline
a1_sg_pred_model_bai_s0 <- glm(sg_crit123 ~  bai_s0, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_bai_s0)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_bai_s0 <- exp(coef(a1_sg_pred_model_bai_s0))[[2]]
a1_LL_bai_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_bai_s0), confint(a1_sg_pred_model_bai_s0)))[ , 2][[2]]
a1_UL_bai_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_bai_s0), confint(a1_sg_pred_model_bai_s0)))[ , 3][[2]]

a1_p_bai_s0 <- tidy(a1_sg_pred_model_bai_s0)[[2 , 5]]

# Sample 1: Comorbid depression
a1_sg_pred_model_comorbid_depression <- glm(sg_crit123 ~  comorbid_depression, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_comorbid_depression)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_comorbid_depression <- exp(coef(a1_sg_pred_model_comorbid_depression))[[2]]
a1_LL_comorbid_depression <- exp(cbind(OR = coef(a1_sg_pred_model_comorbid_depression), confint(a1_sg_pred_model_comorbid_depression)))[ , 2][[2]]
a1_UL_comorbid_depression <- exp(cbind(OR = coef(a1_sg_pred_model_comorbid_depression), confint(a1_sg_pred_model_comorbid_depression)))[ , 3][[2]]

a1_p_comorbid_depression <- tidy(a1_sg_pred_model_comorbid_depression)[[2 , 5]]

# Sample 1: Negative appraisals at baseline
a1_sg_pred_model_ptci22r_s0 <- glm(sg_crit123 ~  ptci22r_s0, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_ptci22r_s0)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_ptci22r_s0 <- exp(coef(a1_sg_pred_model_ptci22r_s0))[[2]]
a1_LL_ptci22r_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_ptci22r_s0), confint(a1_sg_pred_model_ptci22r_s0)))[ , 2][[2]]
a1_UL_ptci22r_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_ptci22r_s0), confint(a1_sg_pred_model_ptci22r_s0)))[ , 3][[2]]

# Extract p value from model
a1_p_ptci22r_s0 <- tidy(a1_sg_pred_model_ptci22r_s0)[[2 , 5]]

# Sample 1: Memory characteristics at baseline
a1_sg_pred_model_mem4_s0 <- glm(sg_crit123 ~  mem4_s0, family = binomial(), data = data_a1_byperson)
# summary(a1_sg_pred_model_mem4_s0)

# Calculate odds ratio and 95% CI of the estimates
a1_OR_mem4_s0 <- exp(coef(a1_sg_pred_model_mem4_s0))[[2]]
a1_LL_mem4_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_mem4_s0), confint(a1_sg_pred_model_mem4_s0)))[ , 2][[2]]
a1_UL_mem4_s0 <- exp(cbind(OR = coef(a1_sg_pred_model_mem4_s0), confint(a1_sg_pred_model_mem4_s0)))[ , 3][[2]]

# Extract p value from model
a1_p_mem4_s0 <- tidy(a1_sg_pred_model_mem4_s0)[[2 , 5]]


# 2.2.1 Sample 1: Multivariate logistic regression ----
a1_sg_pred_model <- glm(sg_crit123 ~ 
                        # Demographic predictors 
                        age + sex + months_since_trauma +
                        # Pretreatment psychopathology
                        pds_s0 + bdi_s0 + bai_s0 + comorbid_depression +
                        # Pretreatment cognitive processes
                        ptci22r_s0 + mem4_s0, 
                        family = binomial(),
                        data = data_a1_byperson)

# Show summary of the model
# summary(a1_sg_pred_model)

# Show odds ratio of the estimates
# exp(coef(a1_sg_pred_model))

# Calculate odds ratio and 95% CI of the estimates
a1_OR_multi <- exp(cbind(OR = coef(a1_sg_pred_model), confint(a1_sg_pred_model)))[ , 1]
a1_LL_multi <- exp(cbind(OR = coef(a1_sg_pred_model), confint(a1_sg_pred_model)))[ , 2]
a1_UL_multi <- exp(cbind(OR = coef(a1_sg_pred_model), confint(a1_sg_pred_model)))[ , 3]

# Extract estimates from model
a1_b_terms_multi <- tidy(a1_sg_pred_model)[ , 1]

# Extract p value from model
a1_b_p_multi <- tidy(a1_sg_pred_model)[ , 5]

# 2.3.1 Sample 1: Multivariate logistic regression using only significant predictors from univariate logistig regressions ----
a1_sg_pred_model_sig <- glm(sg_crit123 ~ age  + bai_s0 + comorbid_depression,
                            family = binomial(),
                            data = data_a1_byperson)

# Show summary of the model
# summary(a1_sg_pred_model_sig)

# Show odds ratio of the estimates
# exp(coef(a1_sg_pred_model_sig))
# exp(cbind(OR = coef(a1_sg_pred_model_sig), confint(a1_sg_pred_model_sig)))

# 2.1.2 Sample 2: Univariate logistic regressions ----

# Sample 2: Age
a2_sg_pred_model_age <- glm(sg_crit123 ~  age, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_age)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_age <- exp(coef(a2_sg_pred_model_age))[[2]]
a2_LL_age <- exp(cbind(OR = coef(a2_sg_pred_model_age), confint(a2_sg_pred_model_age)))[ , 2][[2]]
a2_UL_age <- exp(cbind(OR = coef(a2_sg_pred_model_age), confint(a2_sg_pred_model_age)))[ , 3][[2]]

# Extract p value from model
a2_p_age <- tidy(a2_sg_pred_model_age)[[2 , 5]]

# Sample 2: Gender
a2_sg_pred_model_sex <- glm(sg_crit123 ~  sex, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_sex)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_sex <- exp(coef(a2_sg_pred_model_sex))[[2]]
a2_LL_sex <- exp(cbind(OR = coef(a2_sg_pred_model_sex), confint(a2_sg_pred_model_sex)))[ , 2][[2]]
a2_UL_sex <- exp(cbind(OR = coef(a2_sg_pred_model_sex), confint(a2_sg_pred_model_sex)))[ , 3][[2]]

# Extract p value from model
a2_p_sex <- tidy(a2_sg_pred_model_sex)[[2 , 5]]

# Sample 2: Months since trauma
a2_sg_pred_model_months_since_trauma <- glm(sg_crit123 ~  months_since_trauma, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_months_since_trauma)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_months_since_trauma <- exp(coef(a2_sg_pred_model_months_since_trauma))[[2]]
a2_LL_months_since_trauma <- exp(cbind(OR = coef(a2_sg_pred_model_months_since_trauma), confint(a2_sg_pred_model_months_since_trauma)))[ , 2][[2]]
a2_UL_months_since_trauma <- exp(cbind(OR = coef(a2_sg_pred_model_months_since_trauma), confint(a2_sg_pred_model_months_since_trauma)))[ , 3][[2]]

# Extract p value from model
a2_p_months_since_trauma <- tidy(a2_sg_pred_model_months_since_trauma)[[2 , 5]]

# Sample 2: PTSD symptoms at baseline
a2_sg_pred_model_pds_s0 <- glm(sg_crit123 ~  pds_s0, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_pds_s0)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_pds_s0 <- exp(coef(a2_sg_pred_model_pds_s0))[[2]]
a2_LL_pds_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_pds_s0), confint(a2_sg_pred_model_pds_s0)))[ , 2][[2]]
a2_UL_pds_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_pds_s0), confint(a2_sg_pred_model_pds_s0)))[ , 3][[2]]

# Extract p value from model
a2_p_pds_s0 <- tidy(a2_sg_pred_model_pds_s0)[[2 , 5]]

# Sample 2: Depression symptoms at baseline
a2_sg_pred_model_phq9_s0 <- glm(sg_crit123 ~  phq9_s0, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_phq9_s0)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_phq9_s0 <- exp(coef(a2_sg_pred_model_phq9_s0))[[2]]
a2_LL_phq9_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_phq9_s0), confint(a2_sg_pred_model_phq9_s0)))[ , 2][[2]]
a2_UL_phq9_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_phq9_s0), confint(a2_sg_pred_model_phq9_s0)))[ , 3][[2]]

# Extract p value from model
a2_p_phq9_s0 <- tidy(a2_sg_pred_model_phq9_s0)[[2 , 5]]

# Sample 2: Anxiety symptoms at baseline
a2_sg_pred_model_gad7_s0 <- glm(sg_crit123 ~  gad7_s0, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_gad7_s0)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_gad7_s0 <- exp(coef(a2_sg_pred_model_gad7_s0))[[2]]
a2_LL_gad7_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_gad7_s0), confint(a2_sg_pred_model_gad7_s0)))[ , 2][[2]]
a2_UL_gad7_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_gad7_s0), confint(a2_sg_pred_model_gad7_s0)))[ , 3][[2]]

# Extract p value from model
a2_p_gad7_s0 <- tidy(a2_sg_pred_model_gad7_s0)[[2 , 5]]

# Sample 2: Comorbid depression
a2_sg_pred_model_comorbid_depression <- glm(sg_crit123 ~  comorbid_depression, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_comorbid_depression)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_comorbid_depression <- exp(coef(a2_sg_pred_model_comorbid_depression))[[2]]
a2_LL_comorbid_depression <- exp(cbind(OR = coef(a2_sg_pred_model_comorbid_depression), confint(a2_sg_pred_model_comorbid_depression)))[ , 2][[2]]
a2_UL_comorbid_depression <- exp(cbind(OR = coef(a2_sg_pred_model_comorbid_depression), confint(a2_sg_pred_model_comorbid_depression)))[ , 3][[2]]

# Extract p value from model
a2_p_comorbid_depression <- tidy(a2_sg_pred_model_comorbid_depression)[[2 , 5]]

# Sample 2: Negative appraisals at baseline
a2_sg_pred_model_ptci20_s0 <- glm(sg_crit123 ~  ptci20_s0, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_ptci20_s0)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_ptci20_s0 <- exp(coef(a2_sg_pred_model_ptci20_s0))[[2]]
a2_LL_ptci20_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_ptci20_s0), confint(a2_sg_pred_model_ptci20_s0)))[ , 2][[2]]
a2_UL_ptci20_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_ptci20_s0), confint(a2_sg_pred_model_ptci20_s0)))[ , 3][[2]]

# Extract p value from model
a2_p_ptci20_s0 <- tidy(a2_sg_pred_model_ptci20_s0)[[2 , 5]]

# Sample 2: Memory characteristics at baseline
a2_sg_pred_model_mem5_s0 <- glm(sg_crit123 ~  mem5_s0, family = binomial(), data = data_a2_byperson)
# summary(a2_sg_pred_model_mem5_s0)

# Calculate odds ratio and 95% CI of the estimates
a2_OR_mem5_s0 <- exp(coef(a2_sg_pred_model_mem5_s0))[[2]]
a2_LL_mem5_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_mem5_s0), confint(a2_sg_pred_model_mem5_s0)))[ , 2][[2]]
a2_UL_mem5_s0 <- exp(cbind(OR = coef(a2_sg_pred_model_mem5_s0), confint(a2_sg_pred_model_mem5_s0)))[ , 3][[2]]

# Extract p value from model
a2_p_mem5_s0 <- tidy(a2_sg_pred_model_mem5_s0)[[2 , 5]]

# 2.2.2 Sample 2: Run multivariate logistic regression ----
a2_sg_pred_model <- glm(sg_crit123 ~ 
                        # Demographic predictors 
                        age + sex + months_since_trauma +
                        # Pretreatment psychopathology
                        pds_s0 + phq9_s0 + gad7_s0 + comorbid_depression +
                        # Pretreatment cognitive process
                        ptci20_s0 + mem5_s0, 
                        family = binomial(),
                        data = data_a2_byperson)

# Show summary of the model
# summary(a2_sg_pred_model)

# Show odds ratio of the estimates
# exp(coef(a2_sg_pred_model))

a2_OR_multi <- exp(cbind(OR = coef(a2_sg_pred_model), confint(a2_sg_pred_model)))[ , 1]
a2_LL_multi <- exp(cbind(OR = coef(a2_sg_pred_model), confint(a2_sg_pred_model)))[ , 2]
a2_UL_multi <- exp(cbind(OR = coef(a2_sg_pred_model), confint(a2_sg_pred_model)))[ , 3]

# Extract estimates from model
a2_b_terms_multi <- tidy(a2_sg_pred_model)[ , 1]

# Extract p value from model
a2_b_p_multi <- tidy(a2_sg_pred_model)[ , 5]

# 3 Analysis: Change processes around sudden gains ----

# Define contrasts for group (sg vs nosg) differences in change at four time points:
# (1) N-2 to N-1: "sg_pre2_to_pre1"
# (2) N-1 to N: "sg_pre1_to_0"   
# (3) N to N+1: "sg_0_to_post1"  
# (4) N+1 to N+2: "sg_post1_to_post2"

# beta0 - intercept at 0
# beta1 - group0, change from -2 to 0
# beta2 - group0, change from -1 to 0
# beta3 - group0, change from 0 to 1
# beta4 - group0, change from 0 to 2
# beta5 - group0, change from 0 to 3
# beta6 - group
# beta7 - group1, change from -2 to 0 
# beta8 - group1, change from -1 to 0
# beta9 - group1, change from 0 to 1
# beta10 - group1, change from 0 to 2
# beta11 - group1, change from 0 to 3

contrast_matrix <- rbind(
  # Contrasts for difference within sg group
  "sg_pre2_to_pre1" = c(0, -1, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0),
  "sg_pre1_to_0" = c(0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0),
  "sg_0_to_post1" = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
  "sg_post1_to_post2" = c(0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 0),
  # Contrasts for difference within matched group
  "matched_pre2_to_pre1" = c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "matched_pre1_to_0" = c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "matched_0_to_post1" = c(0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0),
  "matched_post1_to_post2" = c(0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
  # Contrasts for differences between sg and matched groups
  "sg_vs_matched_pre2_to_pre1" = c(0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0),
  "sg_vs_matched_pre1_to_0" = c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0),
  "sg_vs_matched_0_to_post1" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
  "sg_vs_matched_post1_to_post2" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0)
  ) 

# 3.1.1 Sample 1: Negative appraisals ----
# Reshape data from wide to long format
data_a1_ptci22r_long <- data_a1_byperson_matchit %>% 
  select(id, 
         sg = sg_crit123, 
         pre3 = sg_ptci22r_pre3, pre2 = sg_ptci22r_pre2, pre1 = sg_ptci22r_pre1,
         post1 = sg_ptci22r_post1, post2 = sg_ptci22r_post2, post3 = sg_ptci22r_post3) %>% 
  gather(c("pre3", "pre2", "pre1", "post1", "post2", "post3"), 
         key = "time", value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "pre3", -2),
         time = replace(time, time == "pre2", -1),
         time = replace(time, time == "pre1", 0),
         time = replace(time, time == "post1", 1),
         time = replace(time, time == "post2", 2),
         time = replace(time, time == "post3", 3),
         time_fct = factor(time, levels = c("0", "-2", "-1", "1", "2", "3"))) %>%
  mutate(time = as.numeric(time))

# Run linear mixed effect model
a1_model_01_ptci22r <- lme(value ~ time_fct + sg + time_fct * sg,
                           random = ~ 1 | id,
                           method = "ML",
                           na.action = na.omit,
                           data = data_a1_ptci22r_long)

# Show summary of the model
# summary(a1_model_01_ptci22r)

# Check assumption of normality of the residuals
qqnorm(resid(a1_model_01_ptci22r))

# Calculate contrasts specified in in the contrast matrix "contrast_matrix"
a1_model_01_contrasts_ptci22r <- glht(a1_model_01_ptci22r, contrast_matrix)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a1_model_01_contrasts_ptci22r, test = adjusted("none"))

# 3.1.2 Sample 2: Negative appraisals ----
# Reshape data from wide to long format
data_a2_ptci20_long <- data_a2_byperson_matchit %>% 
  select(id, 
         sg = sg_crit123, 
         pre3 = sg_ptci20_pre3, pre2 = sg_ptci20_pre2, pre1 = sg_ptci20_pre1,
         post1 = sg_ptci20_post1, post2 = sg_ptci20_post2, post3 = sg_ptci20_post3) %>% 
  gather(c("pre3", "pre2", "pre1",
           "post1", "post2", "post3"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "pre3", -2),
         time = replace(time, time == "pre2", -1),
         time = replace(time, time == "pre1", 0),
         time = replace(time, time == "post1", 1),
         time = replace(time, time == "post2", 2),
         time = replace(time, time == "post3", 3),
         time_fct = factor(time, levels = c("0", "-2", "-1", "1", "2", "3"))) %>% 
  mutate(time = as.numeric(time))

# Run linear mixed effect model
a2_model_01_ptci20 <- lme(value ~ time_fct + sg + time_fct * sg,
                          random = ~ 1 | id,
                          method = "ML",
                          na.action = na.omit,
                          data = data_a2_ptci20_long)

# Show summary of the model
# summary(a2_model_01_ptci20)

# Check assumption of normality of the residuals
qqnorm(resid(a2_model_01_ptci20))

# Calculate contrasts specified in in the contrast matrix "contrast_matrix"
a2_model_01_contrasts_ptci20 <- glht(a2_model_01_ptci20, contrast_matrix)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a2_model_01_contrasts_ptci20, test = adjusted("none"))

# 3.1.3 Samples 1 & 2: Negative appraisals, Pooled effects ----

# Calculate difference scores and SD of difference scores for Sample 1
results_meta$table_cog_a1 <- data_a1_byperson_matchit %>% 
  transmute(id = id, sg_crit123 = sg_crit123,
            ptci22r_diff_pre3_to_pre2 = sg_ptci22r_pre3 - sg_ptci22r_pre2,
            ptci22r_diff_pre2_to_pre1 = sg_ptci22r_pre2 - sg_ptci22r_pre1,
            ptci22r_diff_pre1_to_post1 = sg_ptci22r_pre1 - sg_ptci22r_post1,
            ptci22r_diff_post1_to_post2 = sg_ptci22r_post1 - sg_ptci22r_post2,
            ptci22r_diff_post2_to_post3 = sg_ptci22r_post2 - sg_ptci22r_post3) %>% 
  group_by(sg_crit123) %>% 
  summarise(n_total = n(),
            n_pre3_to_pre2 = sum(!is.na(ptci22r_diff_pre3_to_pre2)),
            n_pre2_to_pre1 = sum(!is.na(ptci22r_diff_pre2_to_pre1)),
            n_pre1_to_post1 = sum(!is.na(ptci22r_diff_pre1_to_post1)),
            n_post1_to_post2 = sum(!is.na(ptci22r_diff_post1_to_post2)),
            n_post2_to_post3 = sum(!is.na(ptci22r_diff_post2_to_post3)),
            sd_pre3_to_pre2 = sd(ptci22r_diff_pre3_to_pre2, na.rm = T),
            sd_pre2_to_pre1 = sd(ptci22r_diff_pre2_to_pre1, na.rm = T),
            sd_pre1_to_post1 = sd(ptci22r_diff_pre1_to_post1, na.rm = T),
            sd_post1_to_post2 = sd(ptci22r_diff_post1_to_post2, na.rm = T),
            sd_post2_to_post3 = sd(ptci22r_diff_post2_to_post3, na.rm = T)) %>% 
  mutate(sample = "sample1")

# Calculate difference scores and SD of difference scores for Sample 2
results_meta$table_cog_a2 <- data_a2_byperson_matchit %>% 
  transmute(id = id, sg_crit123 = sg_crit123,
            ptci20_diff_pre3_to_pre2 = sg_ptci20_pre3 - sg_ptci20_pre2,
            ptci20_diff_pre2_to_pre1 = sg_ptci20_pre2 - sg_ptci20_pre1,
            ptci20_diff_pre1_to_post1 = sg_ptci20_pre1 - sg_ptci20_post1,
            ptci20_diff_post1_to_post2 = sg_ptci20_post1 - sg_ptci20_post2,
            ptci20_diff_post2_to_post3 = sg_ptci20_post2 - sg_ptci20_post3) %>% 
  group_by(sg_crit123) %>% 
  summarise(n_total = n(),
            n_pre3_to_pre2 = sum(!is.na(ptci20_diff_pre3_to_pre2)),
            n_pre2_to_pre1 = sum(!is.na(ptci20_diff_pre2_to_pre1)),
            n_pre1_to_post1 = sum(!is.na(ptci20_diff_pre1_to_post1)),
            n_post1_to_post2 = sum(!is.na(ptci20_diff_post1_to_post2)),
            n_post2_to_post3 = sum(!is.na(ptci20_diff_post2_to_post3)),
            sd_pre3_to_pre2 = sd(ptci20_diff_pre3_to_pre2, na.rm = T),
            sd_pre2_to_pre1 = sd(ptci20_diff_pre2_to_pre1, na.rm = T),
            sd_pre1_to_post1 = sd(ptci20_diff_pre1_to_post1, na.rm = T),
            sd_post1_to_post2 = sd(ptci20_diff_post1_to_post2, na.rm = T),
            sd_post2_to_post3 = sd(ptci20_diff_post2_to_post3, na.rm = T)) %>% 
  mutate(sample = "sample2")

# Create dataframe with all information for meta-analysis
results_meta$table_cog_a1a2 <- tibble(
  study = as.character(c("Sample1.1", "Sample2.1", 
                         "Sample1.2", "Sample2.2", 
                         "Sample1.3", "Sample2.3", 
                         "Sample1.4", "Sample2.4")),
  time = as.character(c("pre3_to_pre2", "pre3_to_pre2", 
                        "pre2_to_pre1", "pre2_to_pre1", 
                        "pre1_to_post1", "pre1_to_post1", 
                        "post1_to_post2", "post1_to_post2")),
  # Sudden gains group
  n1i = c(results_meta$table_cog_a1$n_pre3_to_pre2[[2]], results_meta$table_cog_a2$n_pre3_to_pre2[[2]], 
          results_meta$table_cog_a1$n_pre2_to_pre1[[2]], results_meta$table_cog_a2$n_pre2_to_pre1[[2]], 
          results_meta$table_cog_a1$n_pre1_to_post1[[2]], results_meta$table_cog_a2$n_pre1_to_post1[[2]], 
          results_meta$table_cog_a1$n_post1_to_post2[[2]], results_meta$table_cog_a2$n_post1_to_post2[[2]]),
  m1i = c(table_a1_ptci22r_1$beta[[1]], table_a2_ptci20_1$beta[[1]],
          table_a1_ptci22r_1$beta[[2]], table_a2_ptci20_1$beta[[2]],
          table_a1_ptci22r_1$beta[[3]], table_a2_ptci20_1$beta[[3]],
          table_a1_ptci22r_1$beta[[4]], table_a2_ptci20_1$beta[[4]]),
  sd1i = c(results_meta$table_cog_a1$sd_pre3_to_pre2[[2]], results_meta$table_cog_a2$sd_pre3_to_pre2[[2]],   
           results_meta$table_cog_a1$sd_pre2_to_pre1[[2]], results_meta$table_cog_a2$sd_pre2_to_pre1[[2]],  
           results_meta$table_cog_a1$sd_pre1_to_post1[[2]], results_meta$table_cog_a2$sd_pre1_to_post1[[2]], 
           results_meta$table_cog_a1$sd_post1_to_post2[[2]], results_meta$table_cog_a2$sd_post1_to_post2[[2]]),
  # Matched control group
  n2i = c(results_meta$table_cog_a1$n_pre3_to_pre2[[1]], results_meta$table_cog_a2$n_pre3_to_pre2[[1]], 
          results_meta$table_cog_a1$n_pre2_to_pre1[[1]], results_meta$table_cog_a2$n_pre2_to_pre1[[1]], 
          results_meta$table_cog_a1$n_pre1_to_post1[[1]], results_meta$table_cog_a2$n_pre1_to_post1[[1]], 
          results_meta$table_cog_a1$n_post1_to_post2[[1]], results_meta$table_cog_a2$n_post1_to_post2[[1]]),
  m2i = c(table_a1_ptci22r_2$beta[[1]], table_a2_ptci20_2$beta[[1]],
          table_a1_ptci22r_2$beta[[2]], table_a2_ptci20_2$beta[[2]],
          table_a1_ptci22r_2$beta[[3]], table_a2_ptci20_2$beta[[3]],
          table_a1_ptci22r_2$beta[[4]], table_a2_ptci20_2$beta[[4]]),
  sd2i = c(results_meta$table_cog_a1$sd_pre3_to_pre2[[1]], results_meta$table_cog_a2$sd_pre3_to_pre2[[1]],   
           results_meta$table_cog_a1$sd_pre2_to_pre1[[1]], results_meta$table_cog_a2$sd_pre2_to_pre1[[1]],  
           results_meta$table_cog_a1$sd_pre1_to_post1[[1]], results_meta$table_cog_a2$sd_pre1_to_post1[[1]], 
           results_meta$table_cog_a1$sd_post1_to_post2[[1]], results_meta$table_cog_a2$sd_post1_to_post2[[1]])
  )

# Calculate SMD for negative appraisals 
results_meta$table_smd_cog_a1a2 <- escalc(measure = "SMD",
                                          m1i = m1i,
                                          sd1i = sd1i,
                                          n1i = n1i,
                                          m2i = m2i,
                                          sd2i = sd2i,
                                          n2i = n2i,
                                          data = results_meta$table_cog_a1a2)

# Add character variable indicating each sample, add spaces as these names have to be "unique"
results_meta$table_smd_cog_a1a2$sample <- c("Sample 1", "Sample 2",
                                            "Sample 1 ", "Sample 2 ",
                                            "Sample 1  ", "Sample 2  ",
                                            "Sample 1   ", "Sample 2   ")

# Run meta-analysis for negative appraisals for all time points 
results_meta$smd_cog_a1a2 <- rma(yi = yi, vi = vi,
                                 measure = "SMD",
                                 method = "FE",
                                 slab = sample,
                                 data = results_meta$table_smd_cog_a1a2)

# Run meta analyses for each time interval separately using the subset command
results_meta$smd_cog_a1a2_pre3_to_pre2 <- rma(yi = yi, vi = vi,
                                              measure = "GEN",
                                              method = "FE",
                                              level = 95,
                                              subset = (time == "pre3_to_pre2"),
                                              data = results_meta$table_smd_cog_a1a2)

results_meta$smd_cog_a1a2_pre2_to_pre1 <- rma(yi = yi, vi = vi,
                                              measure = "GEN",
                                              method = "FE",
                                              level = 95,
                                              subset = (time == "pre2_to_pre1"),
                                              data = results_meta$table_smd_cog_a1a2)

results_meta$smd_cog_a1a2_pre1_to_post1 <- rma(yi = yi, vi = vi,
                                               measure = "GEN",
                                               method = "FE",
                                               level = 95,
                                               subset = (time == "pre1_to_post1"),
                                               data = results_meta$table_smd_cog_a1a2)

results_meta$smd_cog_a1a2_post1_to_post2 <- rma(yi = yi, vi = vi,
                                                measure = "GEN",
                                                method = "FE",
                                                level = 95,
                                                subset = (time == "post1_to_post2"),
                                                data = results_meta$table_smd_cog_a1a2)

# 3.2.1 Sample 1: Memory characteristics ----
# Reshape data from wide to long format
data_a1_mem4_long <- data_a1_byperson_matchit %>% 
  select(id, 
         sg = sg_crit123, 
         pre3 = sg_mem4_pre3, pre2 = sg_mem4_pre2, pre1 = sg_mem4_pre1,
         post1 = sg_mem4_post1, post2 = sg_mem4_post2, post3 = sg_mem4_post3) %>% 
  gather(c("pre3", "pre2", "pre1",
           "post1", "post2", "post3"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "pre3", -2),
         time = replace(time, time == "pre2", -1),
         time = replace(time, time == "pre1", 0),
         time = replace(time, time == "post1", 1),
         time = replace(time, time == "post2", 2),
         time = replace(time, time == "post3", 3),
         time_fct = factor(time, levels = c("0", "-2", "-1", "1", "2", "3"))) %>%
  mutate(time = as.numeric(time))

# Run linear mixed effect model
a1_model_01_mem4 <- lme(value ~ time_fct + sg + time_fct * sg,
                        random = ~ 1 | id,
                        method = "ML",
                        na.action = na.omit,
                        data = data_a1_mem4_long)

# Show summary of the model
summary(a1_model_01_mem4)

# Check assumption of normality of the residuals
qqnorm(resid(a1_model_01_mem4))

# Calculate contrasts specified in in the contrast matrix "contrast_matrix"
a1_model_01_contrasts_mem4 <- glht(a1_model_01_mem4, contrast_matrix)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a1_model_01_contrasts_mem4, test = adjusted("none"))

# 3.2.2 Sample 2: Memory characteristics ----
# Reshape data from wide to long format
data_a2_mem5_long <- data_a2_byperson_matchit %>% 
  select(id, 
         sg = sg_crit123, 
         pre3 = sg_mem5_pre3, pre2 = sg_mem5_pre2, pre1 = sg_mem5_pre1,
         post1 = sg_mem5_post1, post2 = sg_mem5_post2, post3 = sg_mem5_post3) %>% 
  gather(c("pre3", "pre2", "pre1",
           "post1", "post2", "post3"), 
         key = "time", 
         value = "value") %>% 
  mutate(sg = as.factor(sg),
         time_chr = time) %>%
  mutate(time = replace(time, time == "pre3", -2),
         time = replace(time, time == "pre2", -1),
         time = replace(time, time == "pre1", 0),
         time = replace(time, time == "post1", 1),
         time = replace(time, time == "post2", 2),
         time = replace(time, time == "post3", 3),
         time_fct = factor(time, levels = c("0", "-2", "-1", "1", "2", "3"))) %>% 
  mutate(time = as.numeric(time))

# Run linear mixed effect model
a2_model_01_mem5 <- lme(value ~ time_fct + time_fct * sg,
                        random = ~ 1 | id,
                        method = "ML",
                        na.action = na.omit,
                        data = data_a2_mem5_long)

# Show summary of the model
# summary(a2_model_01_mem5)

# Check assumption of normality of the residuals
qqnorm(resid(a2_model_01_mem5))

# Calculate contrasts specified in in the contrast matrix "contrast_matrix"
a2_model_01_contrasts_mem5 <- glht(a2_model_01_mem5, contrast_matrix)

# Show contrasts without adjusting p-values for multiple comparisons
# summary(a2_model_01_contrasts_mem5, test = adjusted("none"))

# 3.2.3 Samples 1 & 2: Memory characteristics, Pooled effects ----

# Calculate difference scores and SD of difference scores for Sample 1
results_meta$table_mem_a1 <- data_a1_byperson_matchit %>% 
  transmute(id = id, sg_crit123 = sg_crit123,
            mem4_diff_pre3_to_pre2 = sg_mem4_pre3 - sg_mem4_pre2,
            mem4_diff_pre2_to_pre1 = sg_mem4_pre2 - sg_mem4_pre1,
            mem4_diff_pre1_to_post1 = sg_mem4_pre1 - sg_mem4_post1,
            mem4_diff_post1_to_post2 = sg_mem4_post1 - sg_mem4_post2,
            mem4_diff_post2_to_post3 = sg_mem4_post2 - sg_mem4_post3) %>% 
  group_by(sg_crit123) %>% 
  summarise(n_total = n(),
            n_pre3_to_pre2 = sum(!is.na(mem4_diff_pre3_to_pre2)),
            n_pre2_to_pre1 = sum(!is.na(mem4_diff_pre2_to_pre1)),
            n_pre1_to_post1 = sum(!is.na(mem4_diff_pre1_to_post1)),
            n_post1_to_post2 = sum(!is.na(mem4_diff_post1_to_post2)),
            n_post2_to_post3 = sum(!is.na(mem4_diff_post2_to_post3)),
            sd_pre3_to_pre2 = sd(mem4_diff_pre3_to_pre2, na.rm = T),
            sd_pre2_to_pre1 = sd(mem4_diff_pre2_to_pre1, na.rm = T),
            sd_pre1_to_post1 = sd(mem4_diff_pre1_to_post1, na.rm = T),
            sd_post1_to_post2 = sd(mem4_diff_post1_to_post2, na.rm = T),
            sd_post2_to_post3 = sd(mem4_diff_post2_to_post3, na.rm = T)
  ) %>% 
  mutate(sample = "sample1")

# Calculate difference scores and SD of difference scores for Sample 2
results_meta$table_mem_a2 <- data_a2_byperson_matchit %>% 
  transmute(id = id, sg_crit123 = sg_crit123,
            mem5_diff_pre3_to_pre2 = sg_mem5_pre3 - sg_mem5_pre2,
            mem5_diff_pre2_to_pre1 = sg_mem5_pre2 - sg_mem5_pre1,
            mem5_diff_pre1_to_post1 = sg_mem5_pre1 - sg_mem5_post1,
            mem5_diff_post1_to_post2 = sg_mem5_post1 - sg_mem5_post2,
            mem5_diff_post2_to_post3 = sg_mem5_post2 - sg_mem5_post3) %>% 
  group_by(sg_crit123) %>% 
  summarise(n_total = n(),
            n_pre3_to_pre2 = sum(!is.na(mem5_diff_pre3_to_pre2)),
            n_pre2_to_pre1 = sum(!is.na(mem5_diff_pre2_to_pre1)),
            n_pre1_to_post1 = sum(!is.na(mem5_diff_pre1_to_post1)),
            n_post1_to_post2 = sum(!is.na(mem5_diff_post1_to_post2)),
            n_post2_to_post3 = sum(!is.na(mem5_diff_post2_to_post3)),
            sd_pre3_to_pre2 = sd(mem5_diff_pre3_to_pre2, na.rm = T),
            sd_pre2_to_pre1 = sd(mem5_diff_pre2_to_pre1, na.rm = T),
            sd_pre1_to_post1 = sd(mem5_diff_pre1_to_post1, na.rm = T),
            sd_post1_to_post2 = sd(mem5_diff_post1_to_post2, na.rm = T),
            sd_post2_to_post3 = sd(mem5_diff_post2_to_post3, na.rm = T)
  ) %>% 
  mutate(sample = "sample2")

# Create dataframe with all information for meta-analysis
results_meta$table_mem_a1a2 <- tibble(
  study = as.character(c("Sample1.1", "Sample2.1", 
                         "Sample1.2", "Sample2.2", 
                         "Sample1.3", "Sample2.3", 
                         "Sample1.4", "Sample2.4")),
  time = as.character(c("pre3_to_pre2", "pre3_to_pre2", 
                        "pre2_to_pre1", "pre2_to_pre1", 
                        "pre1_to_post1", "pre1_to_post1", 
                        "post1_to_post2", "post1_to_post2")),
  # Sudden gains group
  n1i = c(results_meta$table_mem_a1$n_pre3_to_pre2[[2]], results_meta$table_mem_a2$n_pre3_to_pre2[[2]] , 
          results_meta$table_mem_a1$n_pre2_to_pre1[[2]], results_meta$table_mem_a2$n_pre2_to_pre1[[2]] , 
          results_meta$table_mem_a1$n_pre1_to_post1[[2]], results_meta$table_mem_a2$n_pre1_to_post1[[2]] , 
          results_meta$table_mem_a1$n_post1_to_post2[[2]], results_meta$table_mem_a2$n_post1_to_post2[[2]]),
  m1i = c(table_a1_mem4_1$beta[[1]], table_a2_mem5_1$beta[[1]],
          table_a1_mem4_1$beta[[2]], table_a2_mem5_1$beta[[2]],
          table_a1_mem4_1$beta[[3]], table_a2_mem5_1$beta[[3]],
          table_a1_mem4_1$beta[[4]], table_a2_mem5_1$beta[[4]]),
  sd1i = c(results_meta$table_mem_a1$sd_pre3_to_pre2[[2]], results_meta$table_mem_a2$sd_pre3_to_pre2[[2]],   
           results_meta$table_mem_a1$sd_pre2_to_pre1[[2]], results_meta$table_mem_a2$sd_pre2_to_pre1[[2]],  
           results_meta$table_mem_a1$sd_pre1_to_post1[[2]], results_meta$table_mem_a2$sd_pre1_to_post1[[2]], 
           results_meta$table_mem_a1$sd_post1_to_post2[[2]], results_meta$table_mem_a2$sd_post1_to_post2[[2]]),
  # Matched control group
  n2i = c(results_meta$table_mem_a1$n_pre3_to_pre2[[1]], results_meta$table_mem_a2$n_pre3_to_pre2[[1]], 
          results_meta$table_mem_a1$n_pre2_to_pre1[[1]], results_meta$table_mem_a2$n_pre2_to_pre1[[1]], 
          results_meta$table_mem_a1$n_pre1_to_post1[[1]], results_meta$table_mem_a2$n_pre1_to_post1[[1]], 
          results_meta$table_mem_a1$n_post1_to_post2[[1]], results_meta$table_mem_a2$n_post1_to_post2[[1]]),
  m2i = c(table_a1_mem4_2$beta[[1]], table_a2_mem5_2$beta[[1]],
          table_a1_mem4_2$beta[[2]], table_a2_mem5_2$beta[[2]],
          table_a1_mem4_2$beta[[3]], table_a2_mem5_2$beta[[3]],
          table_a1_mem4_2$beta[[4]], table_a2_mem5_2$beta[[4]]),
  sd2i = c(results_meta$table_mem_a1$sd_pre3_to_pre2[[1]], results_meta$table_mem_a2$sd_pre3_to_pre2[[1]],   
           results_meta$table_mem_a1$sd_pre2_to_pre1[[1]], results_meta$table_mem_a2$sd_pre2_to_pre1[[1]],  
           results_meta$table_mem_a1$sd_pre1_to_post1[[1]], results_meta$table_mem_a2$sd_pre1_to_post1[[1]], 
           results_meta$table_mem_a1$sd_post1_to_post2[[1]], results_meta$table_mem_a2$sd_post1_to_post2[[1]])
)

# Calculate SMD for memory characteristics
results_meta$table_smd_mem_a1a2 <- escalc(measure = "SMD",
                                          m1i = m1i,
                                          sd1i = sd1i,
                                          n1i = n1i,
                                          m2i = m2i,
                                          sd2i = sd2i,
                                          n2i = n2i,
                                          data = results_meta$table_mem_a1a2
)

# Add character variable indicating each sample, add spaces as these names have to be "unique"
results_meta$table_smd_mem_a1a2$sample <- c("Sample 1", "Sample 2",
                                            "Sample 1 ", "Sample 2 ",
                                            "Sample 1  ", "Sample 2  ",
                                            "Sample 1   ", "Sample 2   ")

# Run meta-analysis for memory characteristics for all time points 
results_meta$smd_mem_a1a2 <- rma(yi = yi, vi = vi,
                                 measure = "SMD",
                                 method = "FE",
                                 slab = sample,
                                 data = results_meta$table_smd_mem_a1a2)

# Run meta analyses for each time interval separately using the subset command
results_meta$smd_mem_a1a2_pre3_to_pre2 <- rma(yi = yi, vi = vi,
                                              measure = "GEN",
                                              method = "FE",
                                              level = 95,
                                              subset = (time == "pre3_to_pre2"),
                                              data = results_meta$table_smd_mem_a1a2)

results_meta$smd_mem_a1a2_pre2_to_pre1 <- rma(yi = yi, vi = vi,
                                              measure = "GEN",
                                              method = "FE",
                                              level = 95,
                                              subset = (time == "pre2_to_pre1"),
                                              data = results_meta$table_smd_mem_a1a2)

results_meta$smd_mem_a1a2_pre1_to_post1 <- rma(yi = yi, vi = vi,
                                               measure = "GEN",
                                               method = "FE",
                                               level = 95,
                                               subset = (time == "pre1_to_post1"),
                                               data = results_meta$table_smd_mem_a1a2)

results_meta$smd_mem_a1a2_post1_to_post2 <- rma(yi = yi, vi = vi,
                                                measure = "GEN",
                                                method = "FE",
                                                level = 95,
                                                subset = (time == "post1_to_post2"),
                                                data = results_meta$table_smd_mem_a1a2)
