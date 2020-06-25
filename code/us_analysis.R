path_scripts <- ("code/scripts/")
us_scripts <- list(
  "libraries_functions.R",
  "us_preprocess.R"
)

lapply(us_scripts, function(x) source(paste0(path_scripts, x)))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=> direct adjustment-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

us_perc_imm65_deahts_pop_covariate <- glm(YY_deaths ~ ZZ_perc_imm65 
                                          + XX_perc_families 
                                          + XX_ratio_beds              
                                          + XX_perc_asthma             
                                          + XX_perc_cancer_breast         
                                          + XX_perc_cancer_lung           
                                          + XX_perc_ch_obstructive_pulm   
                                          + XX_perc_hypertension           
                                          + XX_median_income              
                                          + XX_summer_temp                
                                          + XX_summer_hum                 
                                          + XX_winter_temp                
                                          + XX_winter_hum                 
                                          + XX_perc_age65_over             
                                          + XX_sex_ratio
                                          + XX_days_f0,                    
                                          family  = quasipoisson(link="log"),
                                          offset= log(NP_total_pop),
                                          data = today_us_filt)


summary(us_perc_imm65_deahts_pop_covariate)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=>  Quintile-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Propensity score is divided in quintile and fed to the model
# This take some time: for each date of the last 60 days runs 10 regressions.

# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)

us_perc_imm65_deahts_pop_quint <- glm_pp(
  dat = all_dat,
  dates_tostudy = dates_tostudy,
  filt_cases = to_filt_cases ,
  form = form_ps,
  offset_f = "log(NP_total_pop)",
  var_dep = "YY_deaths",
  var_int = "ZZ_perc_imm65",
  include_model = TRUE,
  verbose = TRUE,
  nation = "us",
  quint = "TRUE",
  inter = "no"
)

today_us_perc_imm65_deahts_pop_quint <- extract_model(us_perc_imm65_deahts_pop_quint, filt_cases = 10)
summary(today_us_perc_imm65_deahts_pop_quint)

extract_mortality_ratio(us_perc_imm65_deahts_pop_quint, filt_cases = 10, incr = 0.1)

ggcovid(us_perc_imm65_deahts_pop_quint, "date", nation = "us")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=> Split by quint and do analysis in each stratum-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# This analysis is done only on the last day without using the function.
# Calculate the PP and quintile, than split the dataframe based on quintile
# and divide it in 5. Easier to do that that implement it in the function.


# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

PropScores_dat <- lm(form_ps, data = today_us_filt)

# dat quint.
dat_quint_split <- 
  today_us_filt %>% 
  mutate(PP = as.factor(ntile(fitted.values(PropScores_dat), 5)),
         PP_cont = fitted.values(PropScores_dat)
  )

# take the dataframe and we nest it by PP
# then we apply the glm 
us_perc_imm65_deahts_pop_pp_splitted <- 
  dat_quint_split %>% 
  nest(-PP) %>% 
  arrange(PP)  %>% 
  mutate(
    fit = map(data, ~ glm(
      YY_deaths ~ 
        ZZ_perc_imm65,
      family = quasipoisson(link = "log"),
      offset = log(NP_total_pop), data = .x)
    ),
    tidy = map(fit, tidy),
    glanced = map(fit, glance),
    augumented = map(fit, augment)
  )


us_perc_imm65_deahts_pop_pp_splitted %>% 
  unnest(tidy) %>% 
  select(PP, term, estimate, std.error, statistic, p.value)

summary(us_perc_imm65_deahts_pop_pp_splitted)

# check the pp_scores
summary(dat_quint_split$PP_cont)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=>  Quintile NO NewYork-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# quintiles without NY State

# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)

us_perc_imm65_deahts_pop_quint_nony <-
  
  all_dat %>% 
  filter(state != "New York") %>% 
  glm_pp(
    dat = .,
    dates_tostudy = max(dates_tostudy),
    filt_cases = 10,
    form = form_ps,
    offset_f = "log(NP_total_pop)",
    var_dep = "YY_deaths",
    var_int = "ZZ_perc_imm65",
    include_model = TRUE,
    verbose = TRUE,
    nation = "us",
    quint = "TRUE",
    inter = "no"
  )

today_us_perc_imm65_deahts_pop_quint_nony <- extract_model(us_perc_imm65_deahts_pop_quint_nony, filt_cases = 10)
summary(today_us_perc_imm65_deahts_pop_quint_nony)

#@@@@@@@@@@@@@@@@@@@@@@@@
# US: PM2.5 quintile-----
#@@@@@@@@@@@@@@@@@@@@@@@@

XX <- grep("XX", names(all_dat), value = TRUE)
# remove XX_pm2.5 added ZZ_perc_imm65
XX <- XX[XX != "XX_pm2.5"]
XX <- c("ZZ_perc_imm65", XX)

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "XX_pm2.5")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)

us_pm2.5_deaths_pop_quint <- 
  
  all_dat %>% 
  glm_pp(
    dates_tostudy = max(dates_tostudy),
    filt_cases = 10,
    form = form_ps, 
    offset_f = "log(NP_total_pop)",
    var_dep = "YY_deaths",
    var_int = "XX_pm2.5",
    include_model = TRUE,
    verbose = TRUE,
    nation = "us",
    quint = "TRUE",
    inter = "no"
  )

today_us_pm2.5_deaths_pop_quint <- extract_model(us_pm2.5_deaths_pop_quint, filt_cases = 10)
summary(today_us_pm2.5_deaths_pop_quint)
extract_mortality_ratio(us_pm2.5_deaths_pop_quint, filt_cases = 10, incr = 1)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=> PP Continuous-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)

us_perc_imm65_deahts_pop <- glm_pp(
  dat = all_dat,
  dates_tostudy = dates_tostudy,
  filt_cases = to_filt_cases,
  form = form_ps,
  offset_f = "log(NP_total_pop)",
  var_dep = "YY_deaths",
  var_int = "ZZ_perc_imm65",
  include_model = TRUE,
  verbose = TRUE,
  nation = "us"
  
)

today_us_perc_imm65_deahts_pop <- extract_model(us_perc_imm65_deahts_pop, filt_cases = 10)
summary(today_us_perc_imm65_deahts_pop)

# this extracts  mortality ratios 
extract_mortality_ratio(us_perc_imm65_deahts_pop, filt_cases = 10, incr = 0.1)

ggcovid(us_perc_imm65_deahts_pop, "date",  filt_cases = 10)
ggcovid(us_perc_imm65_deahts_pop, "cases")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac No New York PP Continuous-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)
#  No-New York
us_perc_imm65_deahts_pop_nonewyork  <-
  all_dat %>%
  filter(state != "New York") %>%
  glm_pp(
    dates_tostudy = dates_tostudy,
    filt_cases = to_filt_cases,
    form = form_ps,
    offset_f = "log(NP_total_pop)",
    var_dep = "YY_deaths",
    var_int = "ZZ_perc_imm65",
    include_model = TRUE,
    verbose = TRUE
  )

today_perc_imm65_deahts_pop_nonewyork <- extract_model(us_perc_imm65_deahts_pop_nonewyork)
summary(today_perc_imm65_deahts_pop_nonewyork)

extract_mortality_ratio(us_perc_imm65_deahts_pop_nonewyork, filt_cases = 10, incr = 0.1)

ggcovid(us_perc_imm65_deahts_pop_nonewyork, "date", nation = "us")
ggcovid(us_perc_imm65_deahts_pop_nonewyork, "cases", nation = "us")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=>  Quintile direct adjustment: state-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

PropScores_dat <- lm(form_ps, data = today_us_filt)

dat_quint <- 
  today_us_filt %>% 
  mutate(PP = as.factor(ntile(fitted.values(PropScores_dat), 5)),
         PP_cont = fitted.values(PropScores_dat)
  ) 

us_imm65_deaths_pop_quint_state <- glm(
                                      YY_deaths ~ ZZ_perc_imm65 + PP + state,
                                      family = quasipoisson(link="log"), 
                                      offset = log(NP_total_pop), 
                                      data = dat_quint)
summary(us_imm65_deaths_pop_quint_state)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=>  Quintile mixed model: state-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

us_imm65_deaths_pop_quint_mixed_state <- glmer(
  YY_deaths ~ ZZ_perc_imm65 + PP + (1|state),
  family = poisson(link="log"), 
  offset = log(NP_total_pop), 
  data = dat_quint)

summary(us_imm65_deaths_pop_quint_mixed_state)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=>  Quintile mixed model: state no NYC-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

dat_quint_no_nys <- dat_quint %>% 
  filter(state != "New York")

us_imm65_deaths_pop_quint_mixed_state_nonys <- glmer(
  YY_deaths ~ ZZ_perc_imm65 + PP + (1|state),
  family = poisson(link="log"), 
  offset = log(NP_total_pop), 
  data = dat_quint_no_nys)

summary(us_imm65_deaths_pop_quint_mixed_state_nonys)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac No Rhode Island PP Continuous-----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)

#  No-RI
us_perc_imm65_deahts_pop_noRI  <-
  all_dat %>%
  filter(state != "Rhode Isand") %>%
  glm_pp(
    dates_tostudy = dates_tostudy,
    filt_cases = to_filt_cases,
    form = form_ps,
    offset_f = "log(NP_total_pop)",
    var_dep = "YY_deaths",
    var_int = "ZZ_perc_imm65",
    include_model = TRUE,
    verbose = TRUE,
    quint = FALSE,
    inter = "no"
  )

today_perc_imm65_deahts_pop_noRI <- extract_model(us_perc_imm65_deahts_pop_noRI)
summary(today_perc_imm65_deahts_pop_noRI)

extract_mortality_ratio(us_perc_imm65_deahts_pop_noRI, filt_cases = 10, incr = 0.1)

ggcovid(us_perc_imm65_deahts_pop_noRI, "date", nation = "us")
ggcovid(us_perc_imm65_deahts_pop_noRI, "cases", nation = "us")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# US: infl_vac 65=> PP Continuous at least one death -----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# again we select the names that start we XX
XX <- grep("XX", names(all_dat), value = TRUE)
# we don't include total beds because we have the bed X population already in
XX <- XX[!XX %in% c("XX_total_beds")]

# formula for propensity score
form_ps <- reformulate(termlabels = XX, response = "logitZZ_perc_imm65")

dates_tostudy <- sort(unique(all_dat$date))[52:length(unique(all_dat$date))]
dates_tostudy <- dates_tostudy[dates_tostudy >= "2020-04-20"]

to_filt_cases <- seq(10, 100, 10)

us_perc_imm65_deahts_pop_d1 <- 
  
  all_dat %>% 
  filter(YY_deaths >= 1) %>% 
  glm_pp(
  dat = .,
  dates_tostudy = dates_tostudy,
  filt_cases = to_filt_cases,
  form = form_ps,
  offset_f = "log(NP_total_pop)",
  var_dep = "YY_deaths",
  var_int = "ZZ_perc_imm65",
  include_model = TRUE,
  verbose = TRUE,
  nation = "us"
  
)

today_us_perc_imm65_deahts_pop_d1 <- extract_model(us_perc_imm65_deahts_pop_d1, filt_cases = 10)
summary(today_us_perc_imm65_deahts_pop_d1)

extract_mortality_ratio(us_perc_imm65_deahts_pop_d1, filt_cases = 10, incr = 0.1)






