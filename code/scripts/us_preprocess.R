# @@@@@@@@@@@@@@@@@@@@@@@
# US: preprocessing-----
# @@@@@@@@@@@@@@@@@@@@@@@

# These are the variable that will be selected
to_select <- c(
  "date", "county", "state", "fips", "cases", "deaths", "total_pop", "perc_families",
  "perc_family_only_onep", "perc_edu_bachelor", "perc_withinternet", "perc_imm65","total_beds", 
  "ratio_beds","perc_alzheimer_dementia","perc_asthma", "perc_atrial_fibrillation",
  "perc_cancer_breast", "perc_cancer_colorectal", "perc_cancer_lung","perc_ch_obstructive_pulm", 
  "perc_chronic_kidney_disease", "perc_depression","perc_diabetes", "perc_heart_failure", 
  "perc_hypertension", "perc_ischemic_heart_disease", "perc_obesity", "perc_rheumatoid_arthritis",
  "perc_stroke", "perc_tobacco_use", "median_income", "pm2.5", "summer_temp", "summer_hum", 
  "winter_temp","winter_hum", "perc_age65_over", "median_age", "sex_ratio", "child_dependency",
  "perc_black", "perc_lat", "perc_white", "perc_asian", "perc_island","perc_other", 
  "perc_two_more_races", "total_tests", "days_f0", "perc_imm65"
)

# data freeze 
date_freeze <- "2020-06-10"

# get the data
df1_us_jhu <- getus_all() %>% 
 
  # Rode Islande has not deaths at the county level:
  # https://coronavirus.jhu.edu/us-map-faq
  # we are going to remove RI
  filter(state != "Rhode Island")

# we get RI from the NYT repository
df1_us_nyt <- getus_all(repo = "nyt") %>% 
  # fips 0 is for the state unassigned deaths
  filter(state == "Rhode Island")

# row bind and filter for till data freeze
df1_us <- bind_rows(df1_us_jhu, df1_us_nyt) %>%
      filter(date <= !!date_freeze)

# some cleaning
suppressWarnings(
df2 <- 
  df1_us %>%
  # calculate age65_over 
  mutate(perc_age65_over = `perc_65-69` + `perc_70-74` + `perc_75-79` + `perc_80-84` + perc_85_over) %>%
  mutate(urban = if_else(urban == "Urban", 1, 0)) %>%
  
  mutate(total_tests = positive + negative) %>%
  # total hospital beds normalized per population
  mutate(ratio_beds = total_beds/total_pop) %>% 
  
  # calculate day since first case
  # 0/0 generate warning
  mutate(f_date = case_when(cases >= 1 ~ date)) %>%
  group_by(fips) %>%
  mutate(f_date = min(f_date, na.rm = TRUE), days_f0 = as.numeric(date - f_date)) %>%
  ungroup() %>%
  mutate(days_f0 = if_else(is.finite(days_f0), days_f0, NA_real_)) %>%
  
  
  # perc races
  # mutate_at and mutate(across) crashes so I have to repeat code
  mutate(perc_black = total_black / total_pop * 100) %>%
  mutate(perc_white = total_white / total_pop * 100) %>%
  mutate(perc_lat = total_latino / total_pop * 100) %>%
  mutate(perc_asian = total_asian / total_pop * 100) %>%
  mutate(perc_island = total_pacific_islander / total_pop * 100) %>%
  mutate(perc_native = total_native / total_pop * 100) %>%
  mutate(perc_other = total_other_race / total_pop * 100) %>%
  mutate(perc_two_more_races = total_two_more_races / total_pop * 100) %>%
  
  # perc divided by 100
  mutate_at(vars(starts_with("perc")), function(x) x / 100) %>% 
  
  # family with one parent together
  mutate(perc_family_only_onep = perc_families_only_female + perc_families_only_male) 

)

# all_dat select and ready for model -------------
all_dat <- df2 %>%
  
  # select variables ----------
  dplyr::select(!!to_select) %>%
  
  # we remove NAs
  na.omit() %>% 
  # filter(cases >= 1) %>%
  # make variables numeric
  mutate_at(vars(
    -date, -county, -state, -fips), 
    as.numeric
  ) %>%
  # rename
  # add XX_ to all the variables excepts those
  rename_at(
    vars(
      -date, -county, -state, -fips, -cases, -deaths,
      -perc_imm65, -cases, -total_pop
    ),
    ~ paste0("XX_", .)
  ) %>%
  rename(ZZ_perc_imm65 = perc_imm65) %>%
  rename_at(vars(deaths), ~ paste0("YY_", .)) %>%
  rename_at(vars(cases), ~ paste0("NC_", .)) %>%
  rename_at(vars(total_pop), ~ paste0("NP_", .)) %>%
  # calculate logitZZ
  mutate(logitZZ_perc_imm65 = logit(ZZ_perc_imm65)) 



# These dataframea are used mainly for data exploration

#  Analyses are done with datataframe all_dat that contains all the dates
# (we fit model for each day)

# Direct adjustments and mixed models are done only on last day

today_us <- all_dat %>% 
  filter(date == max(date))

# all variables only last day
today_us_all <- df2 %>% 
  filter(date == max(date))

# selected variable and filtered only last day
today_us_filt <-  today_us %>% 
  filter(NC_cases >= 10)











