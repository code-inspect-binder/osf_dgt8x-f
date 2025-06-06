# Title: Identify sudden gains (SG) ----

# Author: Milan Wiedemann
# Status: Accepted

# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
# 1.   Identify sudden gains                  # 
# 1.1.1      Sample 1: Define cut-off         #
# 1.2.1      Sample 1: Select cases           #
# 1.3.1      Sample 1: Identify sudden gain   #
#                                             #
# 1.1.2      Sample 2: Define cut-off         #
# 1.2.2      Sample 2: Select cases           #
# 1.3.2      Sample 2: Identify sudden gains  #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # # 

# 1. Identify sudden gains ----

# Load packages
library(tidyverse)
library(suddengains)

# 1.1.1 Sample 1: Define cut-of ----

# Load data: PDS session by session variables to identify SG and create bysg and byperson datasets 
data_a1 <- read_csv("data/data_audit1_pds.csv") %>%
  select(id, pds_s0:pds_s12, pds_end)

# Define cut-off for Sample 1: 6.15, see Method section in paper
define_crit1_cutoff(sd = 10.54, reliability = .83)$standard_error_difference

# 1.2.1 Sample 1: Select cases ----

# All cases with sufficient session-by-session data to apply @Tang1999 SG criteria will be selected
data_a1 <- select_cases(data = data_a1,
                        id_var_name = "id", 
                        sg_var_list = c("pds_s1", "pds_s2", "pds_s3", "pds_s4", 
                                        "pds_s5", "pds_s6", "pds_s7", "pds_s8", 
                                        "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
                        method = "pattern") %>%
  filter(sg_select == TRUE)

# 2.3.1 Sample 1: Identify sudden gains ----

# Identify sudden gains using all three SG criteria and create bysg dataset
data_a1_bysg <- create_bysg(data = data_a1, 
                            sg_crit1_cutoff = 6.15,
                            sg_crit2_pct = .25,
                            id_var_name = "id", 
                            tx_start_var_name = "pds_s0", 
                            tx_end_var_name = "pds_end",
                            sg_var_list = c("pds_s0", 
                                            "pds_s1", "pds_s2", "pds_s3", "pds_s4", 
                                            "pds_s5", "pds_s6", "pds_s7", "pds_s8", 
                                            "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
                            identify_sg_1to2 = TRUE,
                            sg_var_name = "pds") %>% 
                filter(sg_session_n > 1)

# Create byperson dataset
data_a1_byperson <- create_byperson(data = data_a1_bysg, 
                                    sg_crit1_cutoff = 6.15,
                                    sg_crit2_pct = .25,
                                    id_var_name = "id", 
                                    tx_start_var_name = "pds_s0", 
                                    tx_end_var_name = "pds_end",
                                    sg_var_list = c("pds_s0", 
                                                    "pds_s1", "pds_s2", "pds_s3", "pds_s4", 
                                                    "pds_s5", "pds_s6", "pds_s7", "pds_s8", 
                                                    "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
                                    identify_sg_1to2 = TRUE,
                                    sg_var_name = "pds",
                                    multiple_sg_select = "first",
                                    data_is_bysg = TRUE)

# 1.1.2 Sample 2: Define cut-off ----

# Load data: PDS session by session variables to identify SG and create bysg and byperson datasets 
data_a2 <- read_csv("data/data_audit2_pds.csv") %>%
  select(id, pds_s0:pds_s12, pds_end)

# Define cut-off for Sample 2: 6.15 see Method section in paper
define_crit1_cutoff(sd = 10.54, reliability = .83)$standard_error_difference

# 1.2.2 Sample 2: Select cases ----
# All cases with sufficient session-by-session data to apply @Tang1999 SG criteria will be selected
data_a2 <- select_cases(data = data_a2,
                        id_var_name = "id", 
                        sg_var_list = c("pds_s1", "pds_s2", "pds_s3", "pds_s4", 
                                        "pds_s5", "pds_s6", "pds_s7", "pds_s8", 
                                        "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
                        method = "pattern") %>%
  filter(sg_select == TRUE)

# 2.3.2 Sample 2: Identify sudden gains ----

# Identify sudden gains using all three SG criteria and create bysg dataset
data_a2_bysg <- create_bysg(data = data_a2, 
                            sg_crit1_cutoff = 6.15,
                            sg_crit2_pct = .25,
                            id_var_name = "id", 
                            tx_start_var_name = "pds_s0", 
                            tx_end_var_name = "pds_end",
                            sg_var_list = c("pds_s0", 
                                            "pds_s1", "pds_s2", "pds_s3", "pds_s4", 
                                            "pds_s5", "pds_s6", "pds_s7", "pds_s8", 
                                            "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
                            identify_sg_1to2 = TRUE,
                            sg_var_name = "pds") %>% 
  filter(sg_session_n > 1)

# Create byperson dataset
data_a2_byperson <- create_byperson(data = data_a2_bysg, 
                                    sg_crit1_cutoff = 6.15,
                                    sg_crit2_pct = .25,
                                    id_var_name = "id", 
                                    tx_start_var_name = "pds_s0", 
                                    tx_end_var_name = "pds_end",
                                    sg_var_list = c("pds_s0", 
                                                    "pds_s1", "pds_s2", "pds_s3", "pds_s4", 
                                                    "pds_s5", "pds_s6", "pds_s7", "pds_s8", 
                                                    "pds_s9", "pds_s10", "pds_s11", "pds_s12"),
                                    identify_sg_1to2 = TRUE,
                                    sg_var_name = "pds",
                                    multiple_sg_select = "first",
                                    data_is_bysg = TRUE)
