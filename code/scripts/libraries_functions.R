# Jun 4, 2020
# Baltimore, MD
# Claudio Zanettini
#
# Code of the analyses for paper xxxxx 


# All the data can be retrieved using the `covid19census`s package developed by the lab
# Information regarding the package can be found at https://github.com/c1au6i0/covid19census
# The script used to import and preprocess data from the different sources, and the raw-data can be found in the package
# github repository https://github.com/c1au6i0/covid19census/blob/dev/data-raw/import_raw.R
#
# Note that the script relays heavily on the `dplyr` functions `rename_at` and `summarize_at` that
# have been recently  replaced by `across` and `where` in dplyr v1

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Load libraries and defines functions -----
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Note: use package `renv` and `renv.lock` file in parent folder to install the 
# same versions of packages that we used.


# # Use this to install the covid19census developmental pacakge
# devtools::install_github("c1au6i0/covid19census")

library(oddsratio)
library(covid19census)
library(lme4)
library(tidyverse)
library(broom)
library(purrr)
library(knitr)
library(stringr)
library(xtable)
theme_set(theme_bw())


#' logit
#'
#' @param p numeric vector
#'
#' @return logit
logit <- function(p) {
  log(p) - log(1 - p)
}

# As described in the paper, the model was run using many different subset of data
# (e.s several dates X several filtering criteria). To avoid unnecessary repetition of code,
# a series of functions were developed to: calculate the propensity score and the glm model for
# a dataframe (`get_PPScore`),  apply that function (using different filtering criteria)
# to a date nested dataframe (`glm_pp`). Function `extract_model`, is a wrapper to extract the model 
# of last date, with filter criteria cases, and `extract_tidy_results` is used to extract tidy summary regarding var on interest  
# from object returned by `glm_pp`. Finally, `extract_mortality_ratio` calculates mortality ratios using the
# package `oddratio`. There are also some functions  to plot data.
# Of course, the dataframe used by the function needs to have certain structure
# (we did not create  a class), but at least there is some flexibility that provides the possibility of choosing 
# different parameters.

# NOTE: The functions could have be used in a mapply to run all the model in one line of code,
# but for the sake of clarity each main model was run inside a separate chunk of code
# (so yes there is repetition of code).


#' get_PPScore
#'
#' Calculate propensity score of variable of interest, `var_int` X confounders, and glm quasipoisson(link = "log")
#' model  with `var_dep` as response and `var_int` + propensity score as factors, and `offset_f`
#' as offset
#'
#' @param dat dataframe (need to have specific shape and names)
#' @param filt_cases threshold of number of cases x county to be included in the analysis
#' @param form  formula to use for propensity score, this is a `reformulate` expression
#' @param offset_f the offset  "log(NC_cases)"
#' @param var_dep dependent variable
#' @param var_int variable of interest (for example influenza)
#' @param include_model include the glm object in the dataframe
#' @param verbose boolean
#' @param quint if TRUE calculate the quintile of the PP score
#' @param nation "us" for us dataset. Return object changes depending on the nation (italy does not have counties)
#' @param inter c("yes", "no") calculate interaction of var of inter and PP score.
#' @param other_a if  "only_quint" it calculate only the interaction for a model in var_dep ~ var_int:as.factor(PP) where PP is 
#'     in quintile. The "only_quint" and "PP_main" are experimental features, the output cannot be extracted with extract_model.
#'     quint needs to be TRUE for "quint only"
#' @return dataframe of results. Each raw contains is a separate linear regression analysis:
#'  Most columns are self-explantory, `model` and `model_PP` are  glm objects

get_PPScore <- function(
  dat = .x,
  filt_cases = 1,
  form = form_ps,
  offset_f = "log(NC_cases)",
  var_dep = "YY_deaths",
  var_int = "ZZ_perc_imm65",
  quint = FALSE,
  include_model = FALSE,
  verbose = TRUE,
  nation = "us",
  inter = "no",
  PP_main = FALSE,
  other_a = FALSE
  ) {
  
  if(other_a == "quint_only" && quint != TRUE){
    stop("Interaction quint_only is only avaliable when quint = TRUE")
  }
  
  dat <- dat %>%
    filter(dat$NC_cases >= !!filt_cases)
  
  PropScoresLM <- lm(form, data = dat)
  
  
  
  if (quint == TRUE) {
    dat <- 
      dat %>% 
      mutate(PP = as.factor(ntile(fitted.values(PropScoresLM), 5)))
  } else {
    dat[, "PP"] <- fitted.values(PropScoresLM)
  }
  
  # this create the formula that will be used in the analysis
  if (inter == "yes") {
    form_glm <-   as.formula(paste(var_dep, "~", var_int, "*",  "PP"))
  } 
    
  if (inter == "no") {
    form_glm <- reformulate(termlabels = c(var_int, "PP"), response = var_dep)
  }
  
  if (other_a == "quint_only" && quint == TRUE) {
    form_glm <- as.formula(paste(var_dep, "~ ", var_int, ":",  "as.factor(PP)"))
  }
  
  if (other_a == "PP_main" && quint == TRUE) {
    form_glm <- as.formula(paste(var_dep, "~ ","as.factor(PP) + ", var_int, ":",  "as.factor(PP)"))
  }

  
  outcome_PP <-
    # we add PP on the fly to the dataframe 
    dat %>%
    glm(form_glm,
        family = quasipoisson(link = "log"),
        offset = eval(parse(text = offset_f)),
        data = .
    )
  
  # if we are looking only interactions with quintiles we can't filter them out
  if(other_a %in%  c("quint_only", "PP_main" )) {
    res <- tidy(outcome_PP, conf.int = TRUE)
  } else {
    res <- tidy(outcome_PP, conf.int = TRUE) %>%
      filter(term == !!var_int) 
  }
  
  res <- res %>% 
    mutate(
      deaths = sum(dat$YY_deaths),
      cases = sum(dat$NC_cases),
      formula = deparse(formula(outcome_PP)),
      offset_f = offset_f,
      quint = quint
    )
  
  
  if (nation == "us") {
    res <-
      res %>%
      mutate(
        n_counties = nrow(dat)
      )
  }
  
  if (include_model == TRUE) {
    res <- res %>%
      mutate(model = list(outcome_PP)) %>%
      mutate(pp_model = list(PropScoresLM))
  }
  
  if (verbose == TRUE) {
    message(
      cat(paste0(
        "\nDate ",
        first(dat$date2),
        ":\nApplying model ",
        first(res$formula),
        "\noffset: ",
        first(offset_f),
        "\nfilt_cases = ",
        filt_cases,
        "!"
      ))
    )
  }
  res
}


#' glm_pp
#'
#' Given a `covid-19` dataframe, apply the propensity glm analysis to subset of dates, varying inclusion criteria per county/region
#'
#' @param dat dataframe
#' @param filt_cases numeric vector of thresholds of cases X county  (default = 1)
#' @param form_ps  formula to use for propensity score, this is a `reformulate` expression
#' @param offset_f the offset, es "log(NC_cases)"
#' @param var_dep dependent variable, es. "YY_death"
#' @param var_int indipendent variable, es ZZ_perc_imm65"
#' @param include_model include a glm in the dataframe
#' @param verbose boolean
#' @param inter c("yes", "no") calculate interaction of var of inter and PP score.
#' @param other_a if  "only_quint" it calculate only the interaction for a model in var_dep ~ var_int:as.factor(PP) where PP is 
#'     in quintile. The "only_quint" and "PP_main" are experimental features, the output cannot be extracted with extract_model.
#'
#' @return dataframe
#'


glm_pp <- function(dat,
                   dates_tostudy = dates_tostudy,
                   filt_cases,
                   offset_f = "log(NC_cases)",
                   var_dep,
                   var_int,
                   include_model = FALSE,
                   verbose = TRUE,
                   nation = "us",
                   quint = FALSE,
                   inter = "no",
                   other_a  = "FALSE",
                   form_ps = form_ps,
                   ...) {
  
  
  suppressWarnings(
    zz_term <-
      dat %>%
      # filter(YY_deaths > 0) %>%
      filter(date %in% dates_tostudy) %>%
      mutate(date2 = date) %>%
      nest(-date) %>%
      # we repeat the date slice X number of case threeshold that we want to study
      slice(rep(1:n(), each = length(filt_cases))) %>%
      mutate(filt_cases_c = rep(filt_cases, nrow(.) / length(filt_cases))) %>%
      
      mutate(
        fit =
          map2(
            data,
            filt_cases_c,
            get_PPScore,
            form = !!form_ps,
            offset_f = !!offset_f,
            var_int = !!var_int,
            var_dep = !!var_dep,
            nation = !!nation,
            include_model = !!include_model,
            verbose = !!verbose,
            quint = !!quint,
            inter = !!inter,
            other_a = !!other_a
          )
      ) %>%
      unnest(fit)
  )
  
  zz_term_nd <- zz_term %>%
    select(-data)
  
  zz_term_nd
}


#' extract_model
#'
#' @param x onject return by glm_pp
#' @param fill_cases extract the model with that criterion
#' @param model "pp_model" to extract propensity score, "main" for the main model
#'
#' @return most recent model with filt cases = 1
extract_model <- function(x, filt_cases = 10, model = "main") {
  
  if(model == "main") {
    mod_extr <- "model"
  } 
  
  if(model == "pp_model"){
    mod_extr <- model
  } 
  
    mod <- x %>%
        filter(date == max(date), filt_cases_c == !!filt_cases) %>%
        select(!!mod_extr) %>%
        flatten()
    
  mod[[1]]
    
}


#' extract_tidy_results
#'
#' Extract tidy summary of last day available regarding var on interest  from object returned by
#' glm_pp 
#'
#' @param x  object return by glm_pp 
#' @param filt_cases case threshold
#' @return dataframe tidy summary
#'
#' @examples
extract_tidy_results <- function(x, filt_cases){
  
  mod <- extract_model(x, filt_cases = filt_cases)
  dat_f <- tidy(mod, conf.int = TRUE )[2,]
  dat_f[, "degrees of freedom"] <-  mod$df.residual
  
  dat_f
}


#' extract_mortality_ratio
#'
#' Extract mortality ratios from last day of aggregated models
#'
#' @param model_pp model as return by `glm pp`
#' @param incr increments
#' @param filt_cases used to extract model with that filtering threeshold.
#'
#' @return dataframe as returned by `oddratio::or_glm()`
#'
#' @examples
extract_mortality_ratio <- function(model_pp, 
                                    incr, 
                                    filt_cases){
  
  
  model_today <- extract_model(model_pp, filt_cases = filt_cases)
  
  incr_list <- as.list(incr)
  names(incr_list) <- unique(model_pp$term)
  
  or_glm(data = model_today$data, 
         model = model_today, 
         incr = incr_list) %>% 
    filter(predictor ==  unique(model_pp$term)) %>% 
    mutate(filt_cases =!!filt_cases) %>% 
    mutate(increment = as.numeric(increment))
  
}


#@@@@@@@@@@@@@@@@@@@@@
# graph functions-----
#@@@@@@@@@@@@@@@@@@@@@


#' ggcovid
#' 
#' Plot results of `glm_pp`
#'
#' @param dat a glm_pp object
#' @param x_var date or cases. Date plots the time-serie of the slope of the variable of interest in the model, 
#'     cases takes the last day and plot the slope of the variable of interest in the model changing inclusion criterion
#'     for country.
#' @param filt_cases  minimum number of cases to included in the model,
#' @param nation it or us
ggcovid <- function(dat, x_var, filt_cases = 10, nation = "us") {
  if (!x_var %in% c("date", "cases")) stop("x_var can be date or cases!")
  
  if (x_var == "date") {
    to_add <- paste0("cases >=", filt_cases)
    aes_x <- "date"
    dat <- dat %>%
      filter(filt_cases_c == !!filt_cases)
    x_title <- NULL
  }
  
  if (x_var == "cases") {
    to_add <- paste0("; date ", max(dat$date))
    aes_x <- "filt_cases_c"
    dat <- dat %>%
      filter(date == max(date))
    x_title <- "inclusion criteria (minimum cases)"
  }
  
  if (nation == "it") {
    var_facet <- c("estimate", "deaths", "cases")
  }
  
  
  if (nation == "us") {
    var_facet <- c("estimate", "number of counties")
    
    dat <- dat %>%
      rename(`number of counties` = n_counties)
  }
  
  
  # dependent variable used in title
  dep_var <- stringr::word(first(dat$formula), 1)
  

  dat %>%
    pivot_longer(
      cols = !!var_facet,
      names_to = "var_name"
    ) %>%
    mutate(var_name = factor(var_name, levels = !!var_facet)) %>%
    mutate(sign = if_else(p.value < 0.05, "p<0.05", "NS")) %>%
    mutate(sign = factor(sign, levels = c("p<0.05", "NS"))) %>%
    # mutate(std.error = if_else(var_name == "estimate", std.error, 0)) %>%
    
    mutate_at(vars(std.error, conf.low, conf.high), ~ if_else(var_name == "estimate", ., 0)) %>%
    ggplot(aes_string(aes_x, "value", color = "sign", group = "term")) +
    geom_line() +
    geom_point() +
    geom_linerange(aes(ymax = conf.high, ymin = conf.low), show.legend = FALSE) +
    facet_grid(vars(var_name), scales = "free") +
    scale_color_manual(values = c("NS" = "#00BFC4", "p<0.05" = "#F8766D"), drop = FALSE) +
    theme(
      # legend.position = "top",
      legend.title = element_blank(),
      legend.box.background = element_rect(colour = "black")
    ) +
    labs(
      title = paste0("Slope of ", first(dat$term), " over time"),
      subtitle = paste0(
        dep_var, "; offset = ", first(dat$offset_f),
        "\n", to_add
      ),
      x = x_title
    ) +
    expand_limits(y = 0)
}

#' lord of the ring function
#' 
#' Given a COVID-19 dataframe of data, graph a specific variable by state/region in a jitter plot.
#' It will be lapply to the dataframe to get with one function all the graphs needed. One function to graph 
#' them all (this is lame I know)
#'
#' @param dat COVID-19 dataframe
#' @param var_n the column to be plotted
#' @param model_dat boolean, is this for the data of the model or not? (model data have less country)
#' @param log_v boolean, do you want a 
#' @param path_f where to save data
#'
#' @return generate a ggplot
lord_rings <- function(var_n, 
                       dat, 
                       model_dat,
                       log_v = FALSE, 
                       nation,
                       path_f = getwd()
                       ){

  x <- var_n
  model_dat <- paste0("Model data: ", model_dat )
  tot_NA<- paste0("; Total NA: ", sum(is.na(dat[, var_n])))


  if(log_v ==  FALSE){
    log_name_file <- "abs_"
  } 
  
  if(log_v ==  TRUE){
    log_name_file <- "log_"
    var_n <- paste0("log10(", var_n, ")")
  }
  if(model_dat ==  FALSE){
    model_dat_name_file <- "_model_data.png"
  } else {
    model_dat_name_file <- "_all_data.png"
  }
  
  name_file <- paste0(path_f, "/", log_name_file, x, model_dat_name_file)
  message("Drawing graph, ", x, " just a moment!")
  
  
  # create colors to use and change title dynamically 
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  if (nation == "us"){

    n_rep <- ceiling(length(unique(dat$state)) / 8)
    colors_to_use <- rep(cbbPalette, n_rep)[1:length(unique(dat$state))]
    tot_observ <- paste0("; Counties: ", nrow(dat))
    
    plot_r <- dat %>%  
      # mutate(state = factor(state, levels = ordered_levels)) %>% 
      ggplot(aes_string(var_n, "state", colour = "state")) +
          geom_jitter(alpha = 0.5, height = 0.3) +
          theme(legend.title = element_blank(),
                legend.position = "none") +
          scale_color_manual(values = colors_to_use) +
          labs( 
            title = paste0(toupper(var_n)),
            subtitle = paste0(model_dat, tot_NA, tot_observ),
            y = NULL
          )
   
  }
  
  if (nation == "it"){
    
    tot_observ <- paste0("; Regions: ", nrow(dat))
    
    plot_r <- dat %>%  
      ggplot(aes_string(var_n, "region")) +
      geom_col(fill = "grey30") +
      theme(legend.title = element_blank(),
            legend.position = "none") +
      labs( 
        title = paste0(toupper(var_n)),
        subtitle = paste0(model_dat, tot_NA, tot_observ),
        y = NULL
      )

  }
  

  ggsave(name_file, units = "in", width = 6, height = 7)
}



#' the fifth element
#' 
#' Function to plot each variable based on quint of propensity score
#'
#' @param dat COVID-19 dataframe of one date, with var PP with quintiles
#' @param var_n the column to be plotted

#' @param log_v boolean, do you want a 
#' @param path_f where to save data
#'
#' @return generate a ggplot
fifth_element <- function(var_n, 
                          dat, 
                          log_v = FALSE, 
                          path_f = getwd()
){
  
  
  x <- var_n
  
  
  if(log_v ==  FALSE){
    log_name_file <- "abs_"
  } 
  
  if(log_v ==  TRUE){
    log_name_file <- "log_"
    var_n <- paste0("log10(", var_n, ")")
  }
  
  name_file <- paste0(path_f, "/", log_name_file, x, ".png")
  message("Drawing graph, ", x, " just a moment!")
  
  
  # plot_r <-
  dat %>%  
    # mutate(state = factor(state, levels = ordered_levels)) %>% 
    ggplot(aes_string(var_n, "PP")) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    
    theme(legend.title = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = colors_to_use) +
    labs( 
      y = "Quintile"
    )
  
  ggsave(name_file, units = "in", width = 6, height = 7)
}

































