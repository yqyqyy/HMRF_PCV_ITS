# -------------------------------------------------------------------------
# Author: QY Yu
# Date: 08 July 2021
# Purpose: HMRF PCV on Pneumonia ITS data analyses
# -------------------------------------------------------------------------


# Section I. Data reading and cleaning ------------------------------------
# 0. Source R packages and functions for data cleaning and analyses -------
source('R_packages_functions_for_PCV_ITS.R')


# 1.1 Read raw converted Pneumonia data during 2004-2017 (ACP) ------------
df <- mergefile('raw_data/acp_pneu/')


# 1.2 Read raw converted PP lab and culture data during 2004-2019 (PP) ----  
pp_lab <- mergefile('raw_data/PP_lab/') # 34,261
pp_cul <- mergefile('raw_data/PP_culture/') # 34,999


# 1.3 Read raw converted Fracture data during 2004-2019 (fracture) --------
ft_04_19 <- mergefile('raw_data/fracture/')


# 1.4 Read raw converted UTI data during 2004-2018 (uti) ------------------  
uti_04_18 <- mergefile('raw_data/uti/')


# 2.0 Function to clean raw combined data and remove NA data --------------
clean_data <- function(data, code = NULL){
  colnames(data) <- colrename(data)
  clean_dt <- data %>%
    mutate(discharge.date.yyyy.mm.dd = gsub('\n', '', discharge.date.yyyy.mm.dd),
           date.of.birth.yyyy.mm.dd = as.Date(date.of.birth.yyyy.mm.dd),
           admission.date.yyyy.mm.dd = as.Date(admission.date.yyyy.mm.dd),
           discharge.date.yyyy.mm.dd = as.Date(discharge.date.yyyy.mm.dd)) %>%
    mutate(code_col = ifelse('principal.diagnosis.code' %in% names(.), 1, 0),
           principal.diagnosis.code = ifelse(code_col == 1, gsub('\n', '', principal.diagnosis.code), code)) %>% 
    filter(!is.na(date.of.birth.yyyy.mm.dd),
           !is.na(discharge.date.yyyy.mm.dd),
           year(admission.date.yyyy.mm.dd) %in% c(2004:2017)) %>%
    select(reference.key, sex, date.of.birth.yyyy.mm.dd, admission.date.yyyy.mm.dd, discharge.date.yyyy.mm.dd, principal.diagnosis.code)
  clean_dt
}


# 2.1 Clean raw combined Pneumonia data and remove NA data ----------------
tot <- clean_data(data = df)


# 2.2 Clean raw combined PP lab & culture data and remove NA data ---------
names(pp_lab) <- colrename(pp_lab)
names(pp_cul) <- colrename(pp_cul)


### a) PP Lab
pp_lab_cl <- pp_lab %>% 
  mutate(date.of.birth.yyyy.mm.dd = as.Date(date.of.birth.yyyy.mm.dd),
         admission.date.yyyy.mm.dd = as.Date(admission.date.yyyy.mm.dd),
         discharge.date.yyyy.mm.dd = as.Date(discharge.date.yyyy.mm.dd),
         lis.reference.datetime = as.Date(lis.reference.datetime)) %>% 
  select(reference.key, sex, date.of.birth.yyyy.mm.dd, principal.diagnosis.code,
         admission.date.yyyy.mm.dd, discharge.date.yyyy.mm.dd, lis.reference.datetime, 
         lis.test.description, lis.result, lis.specimen) %>% 
  filter(!is.na(date.of.birth.yyyy.mm.dd), # 1
         !is.na(lis.reference.datetime)) # 11,769

pp_lab_cl_pos <- pp_lab_cl %>% 
  mutate(positive = ifelse(grepl('Positive|^Detected|^Streptococcus|^Untypeable|11A|12F|14|15[ABCF]|16F|19[AF]|20|22F|23[AF]|29|3|34|35[BC]|4|5|6[ABC]|7[CF]|8|9[NV]', lis.result), TRUE, FALSE)) %>%
  filter(positive == TRUE) %>% # 2886 positive
  select(reference.key, sex, date.of.birth.yyyy.mm.dd, principal.diagnosis.code,
         admission.date.yyyy.mm.dd, discharge.date.yyyy.mm.dd, lis.reference.datetime) # 2886 obs
pp_lab_cl_pos

# Test lab date in hospitalization period
# as.Date(substr(pp_lab_cl_pos$lis.reference.datetime,1,10)) %within% interval(pp_lab_cl_pos$admission.date.yyyy.mm.dd, pp_lab_cl_pos$discharge.date.yyyy.mm.dd) %>% table()
# -> all within admission period


### b) PP Culture
pp_cul_cl <- pp_cul %>% 
  mutate(date.of.birth.yyyy.mm.dd = as.Date(date.of.birth.yyyy.mm.dd),
         admission.date.yyyy.mm.dd = as.Date(admission.date.yyyy.mm.dd),
         discharge.date.yyyy.mm.dd = as.Date(discharge.date.yyyy.mm.dd),
         bacterial.culture.and.st.within.hospital.stay.ref.date.yyyy.mm.dd.hh.mm = as.Date(bacterial.culture.and.st.within.hospital.stay.ref.date.yyyy.mm.dd.hh.mm)) %>% 
  select(reference.key, sex, date.of.birth.yyyy.mm.dd, principal.diagnosis.code,
         admission.date.yyyy.mm.dd, discharge.date.yyyy.mm.dd, 
         bacterial.culture.and.st.within.hospital.stay.organism.description.corporate,
         bacterial.culture.and.st.within.hospital.stay.organism.description.local,
         bacterial.culture.and.st.within.hospital.stay.ref.date.yyyy.mm.dd.hh.mm) %>% 
  filter(!is.na(date.of.birth.yyyy.mm.dd), # 1
         !is.na(bacterial.culture.and.st.within.hospital.stay.ref.date.yyyy.mm.dd.hh.mm)) # 20,361

pp_cul_cl_pos <- pp_cul_cl %>% 
  # mutate(positive = TRUE) %>%
  select(reference.key, sex, date.of.birth.yyyy.mm.dd, principal.diagnosis.code,
         admission.date.yyyy.mm.dd, discharge.date.yyyy.mm.dd, 
         bacterial.culture.and.st.within.hospital.stay.ref.date.yyyy.mm.dd.hh.mm) %>% 
  rename(lis.reference.datetime = bacterial.culture.and.st.within.hospital.stay.ref.date.yyyy.mm.dd.hh.mm)
pp_cul_cl_pos

# Test lab date in hospitalization period
# as.Date(substr(pp_cul_cl_pos$lis.reference.datetime,1,10)) %within% interval(pp_cul_cl_pos$admission.date.yyyy.mm.dd, pp_cul_cl_pos$discharge.date.yyyy.mm.dd) %>% table()
# -> all within admission period


### c) Merge PP lab and culture data during 2004-2017
pp_lab_cul <- rbind(pp_lab_cl_pos, pp_cul_cl_pos) %>% filter(year(admission.date.yyyy.mm.dd) %in% c(2004:2017))


# 2.3 Clean raw combined Fracture data, and remove NA dob data ------------ 
ft_04_17 <- clean_data(data = ft_04_19, code = 810)


# 2.4 Clean raw combined UTI data, and remove NA dob data -----------------     
uti_04_17 <- clean_data(data = uti_04_18, code = 599)


# 3. Subset ACP and flu dataset for primary outcome analysis --------------  
### a) All-cause pneumonia (acp) data with code 480-486   
acp <- tot %>% filter(grepl('^48[0-6]', tot$principal.diagnosis.code))
table(acp$principal.diagnosis.code)


### b) Flu (flu) data with code 487 (episodes as confounder adjustment in model)  
flu <- tot %>% filter(grepl('^487', principal.diagnosis.code))
table(flu$principal.diagnosis.code)


# 4.0 Function to merge overlap time interval and change dob to age -------
merge_interval_episode <- function(data, time_interval){
  data %>% group_by(reference.key) %>% 
    arrange(reference.key, admission.date.yyyy.mm.dd, discharge.date.yyyy.mm.dd) %>% 
    mutate(indx = c(0, cumsum(as.numeric(lead(admission.date.yyyy.mm.dd)) >
                                (cummax(as.numeric(discharge.date.yyyy.mm.dd)+time_interval)))[-n()])) %>%
    group_by(reference.key, indx) %>% 
    arrange(reference.key, indx) %>% 
    summarise(sex = unique(sex), 
              dob = unique(date.of.birth.yyyy.mm.dd),
              code = paste(unique(principal.diagnosis.code), collapse = ','),
              year = year(min(admission.date.yyyy.mm.dd)), 
              month = month(min(admission.date.yyyy.mm.dd)),
              ad_date = min(admission.date.yyyy.mm.dd), 
              dc_date = max(discharge.date.yyyy.mm.dd)) %>%
    select(-indx) %>% 
    mutate(age = age(dob, ad_date))
}


# 4.1 Merge overlap interval for ACP data as episodes data ----------------
acpdt <- merge_interval_episode(data = acp, time_interval = 1)


# 4.2 Merge overlap interval for flu data as episodes data ----------------
fludt <- merge_interval_episode(data = flu, time_interval = 1)


# 4.3 Merge overlap interval for PP data as episodes data -----------------
ppdt <- merge_interval_episode(data = pp_lab_cul, time_interval = 1)


# 4.4 Merge overlap interval for Fracture data as episodes data -----------
ft_0417_dt <- merge_interval_episode(data = ft_04_17, time_interval = 1)


# 4.5 Merge overlap interval for UTI data as episodes data ----------------  
utidt <- merge_interval_episode(data = uti_04_17, time_interval = 1)


# 5. Remove raw data to clean up memory -----------------------------------
rm(df, pp_lab, pp_cul, ft_04_19, uti_04_18)


# -------------------------------------------------------------------------
# Section II. Age-standardization for overall data set --------------------
# 1. Age population per year data for standardization ---------------------
agepop <- hk_pop_agegp %>% filter(year %in% c(2004:2017)) %>% 
  mutate(agegp = cut(age, c(-1,4,19,49,64,74,84,85), 
                     labels = c('0-4', '5-19', '20-49', '50-64', '65-74', '75-84', '85+'))) %>% 
  select(year, sex, age, hk_pop, agegp)

### a) Combine age group to calculate the combined population / age incidence (proportion)  
agepop_gp <- agepop %>% group_by(year, agegp) %>% summarise(gppop = sum(hk_pop))
agepop_gp

### b) Using HK mid population in 2017 as reference (standard population)  
agepop_ref <- agepop_gp %>% filter(year == 2017)
agepop_ref


# 2.0 Function to get monthly standardized incidence ----------------------
monthly_std_inc <- function(data){
  dt <- data %>% mutate(agegp = cut(age, c(-1,4,19,49,64,74,84,Inf),
                                    labels = c('0-4', '5-19', '20-49', '50-64', '65-74', '75-84', '85+')))
  dt_agegp_inc <- 
    merge(dt %>% group_by(year, month, agegp, .drop = FALSE) %>% summarise(agegp_count = n()),
          agepop_gp, by = c('year', 'agegp')) %>%
    arrange(year, month, agegp) %>%
    mutate(inc = agegp_count/gppop)
  dt_agegp_adjcount <- 
    merge(dt_agegp_inc,
          agepop_ref %>% rename(year_ref = year, stdpop = gppop),
          by = 'agegp') %>% 
    mutate(adj_agegp_count = inc*stdpop) %>% 
    arrange(year, month, agegp)
  monthly_inc_rate <- dt_agegp_adjcount %>% 
    group_by(year, month) %>% 
    summarise(raw_case = sum(agegp_count), 
              cal_case = round(sum(adj_agegp_count)),
              pop = sum(stdpop)) %>% 
    ungroup() %>% 
    mutate(rate = cal_case/pop*10^5, agegp = 'all') %>% 
    arrange(year, month)
  monthly_inc_rate
}


# 2.1 Get ACP monthly age-standardized incidence (Primary outcome) --------
acp_mir <- monthly_std_inc(data = acpdt)


# 2.2 Get flu monthly age-standardized incidence (Reviewers' comments) ----
flu_mir <- monthly_std_inc(data = fludt)


# 2.3 Get PP monthly age-standardized incidence (Secondary outcome) -------
pp_mir <- monthly_std_inc(data = ppdt)


# 2.4 Get Fracture monthly standardized incidence (Negative outcome 1) ---- 
ft_mir <- monthly_std_inc(data = ft_0417_dt)


# 2.5 Get UTI monthly standardized incidence (Negative outcome 2) --------- 
uti_mir <- monthly_std_inc(data = utidt)


# -------------------------------------------------------------------------
# Section III. Subgroup data set incidence rate ---------------------------
# 1. Three subgroups for age population per year --------------------------
agepop_3g <- agepop %>% 
  mutate(age3g = ifelse(age %in% c(0:19), '0-19', ifelse(age %in% c(20:64), '20-64', '65+'))) %>% 
  group_by(year, age3g) %>% summarise(gppop = sum(hk_pop))
agepop_3g


# 2.0 Function to get monthly incidence for three subgroups ---------------
subgroup3_inc <- function(data){
  dt <- data %>% mutate(age3g = cut(age, c(-1,19,64,Inf), labels = c('0-19', '20-64', '65+')))
  monthly_inc_rate_age3g <- 
    merge(dt %>% group_by(year, month, age3g, .drop = FALSE) %>% summarise(cal_case = n()),
          agepop_3g, by = c('year', 'age3g'), all.x = T) %>%
    mutate(rate = cal_case/gppop*10^5, raw_case = cal_case) %>% 
    arrange(year, month) %>% 
    rename(agegp = age3g, pop = gppop)
  monthly_inc_rate_age3g
}


# 2.1 Three subgroups for All-cause pneumonia (0-19y, 20-64y, 65+y) -------
acp_mir_age3g <- subgroup3_inc(data = acpdt)


# 2.2 Three subgroups for flu (0-19y, 20-64y, 65+y) -----------------------
flu_mir_age3g <- subgroup3_inc(data = fludt)


# 2.3 Three subgroups for PP (0-19y, 20-64y, 65+y) ------------------------
pp_mir_age3g <- subgroup3_inc(data = ppdt)


# 2.4 Three subgroups for Fracture (0-19y, 20-64y, 65+y) ------------------
ft_mir_age3g <- subgroup3_inc(data = ft_0417_dt)


# 2.5 Three subgroups for UTI (0-19y, 20-64y, 65+y) -----------------------
uti_mir_age3g <- subgroup3_inc(data = utidt)


# 3. Seven subgroups for age population per year --------------------------
agepop_7g <- agepop %>% group_by(year, agegp) %>% summarise(gppop = sum(hk_pop))
agepop_7g


# 4.0 Function to get monthly incidence for seven subgroups ---------------
subgroup7_inc <- function(data){
  dt <- data %>% mutate(agegp = cut(age, c(-1,4,19,49,64,74,84,Inf), 
                                    labels = c('0-4', '5-19', '20-49', '50-64', '65-74', '75-84', '85+')))
  monthly_inc_rate_age7g <- 
    merge(dt %>% group_by(year, month, agegp, .drop = FALSE) %>% summarise(cal_case = n()),
          agepop_7g, by = c('year', 'agegp'), all.x = T) %>%
    mutate(rate = cal_case/gppop*10^5, raw_case = cal_case) %>% 
    arrange(year, month) %>% 
    rename(pop = gppop)
  monthly_inc_rate_age7g
}


# 4.1 Seven subgroups for All-cause pneumonia -----------------------------
acp_mir_age7g <- subgroup7_inc(data = acpdt)


# 4.2 Seven subgroups for flu ---------------------------------------------
flu_mir_age7g <- subgroup7_inc(data = fludt)


# 4.3 Seven subgroups for PP ----------------------------------------------
pp_mir_age7g <- subgroup7_inc(data = ppdt)


# 4.4 Seven subgroups for Fracture ----------------------------------------
ft_mir_age7g <- subgroup7_inc(data = ft_0417_dt)


# 4.5 Seven subgroups for UTI ---------------------------------------------
uti_mir_age7g <- subgroup7_inc(data = utidt)


# -------------------------------------------------------------------------
# Section IV. Data preparation --------------------------------------------
# 0. Function to set up the time series data set --------------------------
time_series_data <- function(data, subgroup, flu_data = flu_mir){
  tempdf <- data %>%
      filter(agegp %in% subgroup) %>%
      select(year, month, cal_case, pop, rate) %>%
      arrange(year, month) %>%
      mutate(month = factor(month),
             time = c(0: (length(year)-1)),
             timec = time - 96, # time centred at pcv introduction
             timeafter = ifelse(timec < 0, 0, timec)) %>%
      mutate(pcv = rep(c(0,NA,1),c(60,36,72))) %>%
      mutate(flu_ad = log(flu_data %>% filter(agegp %in% subgroup) %>% pull(raw_case)),
             flu_ad = ifelse(flu_ad == -Inf, 0, flu_ad)) %>%
      mutate(pcv3 = rep(c(0,NA,1),c(60,39,69)),
             timeafter3 = ifelse((time-99) < 0, 0, (time-99)),
             pcv6 = rep(c(0,NA,1),c(60,42,66)),
             timeafter6 = ifelse((time-102) < 0, 0, (time-102)),
             pcv12 = rep(c(0,NA,1),c(60,48,60)),
             timeafter12 = ifelse((time-108) < 0, 0, (time-108)))
  tempdf
}


# 1.1 Prepare overall data set for ACP data -------------------------------
pcv <- time_series_data(data = acp_mir, subgroup = 'all', flu_data = flu_mir)


# 1.2 Prepare 3 subgroup data set for ACP data ----------------------------   
age019 <- time_series_data(data = acp_mir_age3g, subgroup = '0-19', flu_data = flu_mir_age3g)
age2064 <- time_series_data(data = acp_mir_age3g, subgroup = '20-64', flu_data = flu_mir_age3g)
age65 <- time_series_data(data = acp_mir_age3g, subgroup = '65+', flu_data = flu_mir_age3g)


# 1.3 Prepare 7 subgroup data set for ACP data (for HMRF Reviewers) -------   
age04 <- time_series_data(data = acp_mir_age7g, subgroup = '0-4', flu_data = flu_mir_age7g)
age519 <- time_series_data(data = acp_mir_age7g, subgroup = '5-19', flu_data = flu_mir_age7g)
age2049 <- time_series_data(data = acp_mir_age7g, subgroup = '20-49', flu_data = flu_mir_age7g)
age5064 <- time_series_data(data = acp_mir_age7g, subgroup = '50-64', flu_data = flu_mir_age7g)
age6574 <- time_series_data(data = acp_mir_age7g, subgroup = '65-74', flu_data = flu_mir_age7g)
age7584 <- time_series_data(data = acp_mir_age7g, subgroup = '75-84', flu_data = flu_mir_age7g)
age85 <- time_series_data(data = acp_mir_age7g, subgroup = '85+', flu_data = flu_mir_age7g)


# 2.1 Prepare overall data set for PP data --------------------------------
pp_pcv <- time_series_data(data = pp_mir, subgroup = 'all', flu_data = flu_mir)


# 2.2 Prepare 3 subgroup data set for PP data -----------------------------
ppage019 <- time_series_data(data = pp_mir_age3g, subgroup = '0-19', flu_data = flu_mir_age3g)
ppage2064 <- time_series_data(data = pp_mir_age3g, subgroup = '20-64', flu_data = flu_mir_age3g)
ppage65 <- time_series_data(data = pp_mir_age3g, subgroup = '65+', flu_data = flu_mir_age3g)


# 2.3 Prepare 7 subgroup data set for PP data (for HMRF Reviewers) --------   
ppage04 <- time_series_data(data = pp_mir_age7g, subgroup = '0-4', flu_data = flu_mir_age7g)
ppage519 <- time_series_data(data = pp_mir_age7g, subgroup = '5-19', flu_data = flu_mir_age7g)
ppage2049 <- time_series_data(data = pp_mir_age7g, subgroup = '20-49', flu_data = flu_mir_age7g)
ppage5064 <- time_series_data(data = pp_mir_age7g, subgroup = '50-64', flu_data = flu_mir_age7g)
ppage6574 <- time_series_data(data = pp_mir_age7g, subgroup = '65-74', flu_data = flu_mir_age7g)
ppage7584 <- time_series_data(data = pp_mir_age7g, subgroup = '75-84', flu_data = flu_mir_age7g)
ppage85 <- time_series_data(data = pp_mir_age7g, subgroup = '85+', flu_data = flu_mir_age7g)


# 3.1 Prepare overall data set for flu data (common model trend test*) ----
influenza <- time_series_data(data = flu_mir, subgroup = 'all', flu_data = flu_mir)


# 3.2 Prepare 3 subgroup data set for flu (for reference*) ----------------
flu_age019 <- time_series_data(data = flu_mir_age3g, subgroup = '0-19', flu_data = flu_mir_age3g)
flu_age2064 <- time_series_data(data = flu_mir_age3g, subgroup = '20-64', flu_data = flu_mir_age3g)
flu_age65 <- time_series_data(data = flu_mir_age3g, subgroup = '65+', flu_data = flu_mir_age3g)


# 3.3 Prepare 7 subgroup data set for flu data (for reference*) -----------   
fluage04 <- time_series_data(data = flu_mir_age7g, subgroup = '0-4', flu_data = flu_mir_age7g)
fluage519 <- time_series_data(data = flu_mir_age7g, subgroup = '5-19', flu_data = flu_mir_age7g)
fluage2049 <- time_series_data(data = flu_mir_age7g, subgroup = '20-49', flu_data = flu_mir_age7g)
fluage5064 <- time_series_data(data = flu_mir_age7g, subgroup = '50-64', flu_data = flu_mir_age7g)
fluage6574 <- time_series_data(data = flu_mir_age7g, subgroup = '65-74', flu_data = flu_mir_age7g)
fluage7584 <- time_series_data(data = flu_mir_age7g, subgroup = '75-84', flu_data = flu_mir_age7g)
fluage85 <- time_series_data(data = flu_mir_age7g, subgroup = '85+', flu_data = flu_mir_age7g)


# 4.1 Prepare overall data set for Fracture data --------------------------
ftdt <- time_series_data(data = ft_mir, subgroup = 'all', flu_data = flu_mir)


# 4.2 Prepare 3 subgroup data set for Fracture data -----------------------
ftage019 <- time_series_data(data = ft_mir_age3g, subgroup = '0-19', flu_data = flu_mir_age3g)
ftage2064 <- time_series_data(data = ft_mir_age3g, subgroup = '20-64', flu_data = flu_mir_age3g)
ftage65 <- time_series_data(data = ft_mir_age3g, subgroup = '65+', flu_data = flu_mir_age3g)


# 4.3 Prepare 7 subgroup data set for Fracture data -----------------------   
ftage04 <- time_series_data(data = ft_mir_age7g, subgroup = '0-4', flu_data = flu_mir_age7g)
ftage519 <- time_series_data(data = ft_mir_age7g, subgroup = '5-19', flu_data = flu_mir_age7g)
ftage2049 <- time_series_data(data = ft_mir_age7g, subgroup = '20-49', flu_data = flu_mir_age7g)
ftage5064 <- time_series_data(data = ft_mir_age7g, subgroup = '50-64', flu_data = flu_mir_age7g)
ftage6574 <- time_series_data(data = ft_mir_age7g, subgroup = '65-74', flu_data = flu_mir_age7g)
ftage7584 <- time_series_data(data = ft_mir_age7g, subgroup = '75-84', flu_data = flu_mir_age7g)
ftage85 <- time_series_data(data = ft_mir_age7g, subgroup = '85+', flu_data = flu_mir_age7g)


# 5.1 Prepare overall data set for UTI data -------------------------------
uti <- time_series_data(data = uti_mir, subgroup = 'all', flu_data = flu_mir)


# 5.2 Prepare 3 subgroup data set for UTI data ----------------------------
utiage019 <- time_series_data(data = uti_mir_age3g, subgroup = '0-19', flu_data = flu_mir_age3g)
utiage2064 <- time_series_data(data = uti_mir_age3g, subgroup = '20-64', flu_data = flu_mir_age3g)
utiage65 <- time_series_data(data = uti_mir_age3g, subgroup = '65+', flu_data = flu_mir_age3g)


# 5.3 Prepare 7 subgroup data set for UTI data ----------------------------   
utiage04 <- time_series_data(data = uti_mir_age7g, subgroup = '0-4', flu_data = flu_mir_age7g)
utiage519 <- time_series_data(data = uti_mir_age7g, subgroup = '5-19', flu_data = flu_mir_age7g)
utiage2049 <- time_series_data(data = uti_mir_age7g, subgroup = '20-49', flu_data = flu_mir_age7g)
utiage5064 <- time_series_data(data = uti_mir_age7g, subgroup = '50-64', flu_data = flu_mir_age7g)
utiage6574 <- time_series_data(data = uti_mir_age7g, subgroup = '65-74', flu_data = flu_mir_age7g)
utiage7584 <- time_series_data(data = uti_mir_age7g, subgroup = '75-84', flu_data = flu_mir_age7g)
utiage85 <- time_series_data(data = uti_mir_age7g, subgroup = '85+', flu_data = flu_mir_age7g)


# -------------------------------------------------------------------------
# Section V. Data descriptive analyses ------------------------------------
# 1.0 Function to summarize data statistics -------------------------------  
summary_statistics <- function(inc_rate_dt, case_dt, clean_dt){
  statistics <- 
    list(min_rate = min(inc_rate_dt$rate),
         max_rate = max(inc_rate_dt$rate),
         avg_rate = mean(inc_rate_dt$rate),
         sd_rate = sd(inc_rate_dt$rate),
         min_case = min(case_dt$raw_case),
         max_case = max(case_dt$raw_case),
         avg_case = mean(case_dt$raw_case),
         sd_case = sd(case_dt$raw_case),
         sex_dist = clean_dt %>% select(reference.key, sex) %>% unique(.) %>% pull(sex) %>% table(.),
         sex_prop = clean_dt %>% select(reference.key, sex) %>% unique(.) %>% pull(sex) %>% table(.) %>% prop.table(.))
  statistics
}


# 1.1 Summary statistics of monthly ACP incidence/episodes (Primary) ------
acp_statistics <- summary_statistics(inc_rate_dt = pcv, case_dt = acp_mir, clean_dt = acpdt)


# 1.2 Summary statistics of monthly PP incidence/episodes (Secondary) -----  
pp_statistics <- summary_statistics(inc_rate_dt = pp_pcv, case_dt = pp_mir, clean_dt = ppdt)


# 1.3 Summary statistics of monthly flu incidence/episodes (reference) ----
flu_statistics <- summary_statistics(inc_rate_dt = influenza, case_dt = flu_mir, clean_dt = fludt)


# 2.0 Function to plot trend of monthly incidence (for reference*) --------
monthly_inc_plot <- function(data, min_ylim = 0, max_ylim = 100, y_lab = 'Std monthly incidence x 100,000', line_color = 3, title_no = "(a)"){
  plot(data$rate, type = "n", ylim = c(min_ylim, max_ylim), xlab = "Year", font.lab = 2,
       ylab = y_lab, bty = "n", xaxt = "n", yaxt = "n")
  rect(60,min_ylim,96,max_ylim, col = grey(0.9), border = F)
  points(data$rate, cex = 0.7)
  lines(data$rate, cex = 0.7, lwd = 2, col = line_color)
  axis(1, at = 0:14*12, labels = F, lwd = 3)
  axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, srt = 45, lwd = 3, font = 4, las = 2)
  axis(2, lwd = 3, las = 1, font = 4)
  par(adj = 0.02)
  title(title_no, font.lab = 1)
}


# 2.1 Scatter plot of monthly incidence of ACP (for reference*) -----------
monthly_inc_plot(data = pcv)


# 2.2 Scatter plot of monthly incidence of PP (for reference*) ------------
monthly_inc_plot(data = pp_pcv, max_ylim = 3, line_color = 6, title_no = "(b)")


# 2.3 Scatter plot of monthly incidence of Fracture (for reference*) ------
monthly_inc_plot(data = ftdt, min_ylim = 0, max_ylim = 50, line_color = 4, title_no = "(c)")


# 2.4 Scatter plot of monthly incidence of UTI (for reference*) -----------
monthly_inc_plot(data = uti, min_ylim = 0, max_ylim = 50, line_color = 7, title_no = "(d)")


# 3.0 Function to summarize annual incidence ------------------------------
annual_inc <- function(data){
  data %>% group_by(year) %>% summarise(yrcount = sum(cal_case), yrpop = unique(pop), yrrate = yrcount/yrpop*10^5)
}


# 3.1 Annual incidence of ACP, overall and by age (3 groups) --------------
acp_annual_inc <- 
  list(inc_all = annual_inc(pcv), inc_019 = annual_inc(age019), inc_2064 = annual_inc(age2064), inc_65 = annual_inc(age65),
       inc_04 = annual_inc(age04), inc_519 = annual_inc(age519), inc_2049 = annual_inc(age2049), 
       inc_5064 = annual_inc(age5064), inc_6574 = annual_inc(age6574), inc_7584 = annual_inc(age7584), inc_85 = annual_inc(age85))


# 3.2 Annual incidence of PP, overall and by age (3 groups) ---------------
pp_annual_inc <- 
  list(inc_all = annual_inc(pp_pcv), inc_019 = annual_inc(ppage019), inc_2064 = annual_inc(ppage2064), inc_65 = annual_inc(ppage65),
       inc_04 = annual_inc(ppage04), inc_519 = annual_inc(ppage519), inc_2049 = annual_inc(ppage2049), 
       inc_5064 = annual_inc(ppage5064), inc_6574 = annual_inc(ppage6574), inc_7584 = annual_inc(ppage7584), inc_85 = annual_inc(ppage85))


# 4.0 Function to plot trend of annual incidence --------------------------
annual_inc_plot <- function(data, min_ylim = 0, max_ylim = 1000, y_lab = 'Std annual incidence x 100,000', line_color = 3, title_no = "(a)"){
  plot(annual_inc(data)$yrrate, type = "n", ylim = c(min_ylim, max_ylim), xlab = "Year", font.lab = 2,
       ylab = y_lab, bty = "n", xaxt = "n", yaxt = "n")
  rect(6,min_ylim,8,max_ylim, col = grey(0.9), border = F)
  points(annual_inc(data)$yrrate, cex = 0.7)
  lines(annual_inc(data)$yrrate, cex = 0.7, lwd = 2, col = line_color)
  axis(1, at = 1:14, labels = F, lwd = 3)
  axis(1, at = 1:14, tick = F, labels = 2004:2017, srt = 45, lwd = 3, font = 4, las = 2)
  axis(2, lwd = 3, las = 1, font = 4)
  par(adj = 0.02)
  title(title_no, font.lab = 1)
}


# 4.1 Scatter plot for annual incidence of ACP ---------------------------
annual_inc_plot(data = pcv, min_ylim = 0, max_ylim = 1000, y_lab = 'Std annual incidence x 100,000')


# 4.2 Scatter plot for annual incidence of PP ----------------------------
annual_inc_plot(data = pp_pcv, min_ylim = 0, max_ylim = 30, y_lab = 'Std annual incidence x 100,000', line_color = 6, title_no = "(b)")


# 4.3 Scatter plot for annual incidence of Fracture ----------------------
annual_inc_plot(data = ftdt, min_ylim = 0, max_ylim = 500, y_lab = 'Std annual incidence x 100,000', line_color = 4, title_no = "(c)")


# 4.4 Scatter plot for annual incidence of UTI ----------------------------
annual_inc_plot(data = uti, min_ylim = 0, max_ylim = 500, y_lab = 'Std annual incidence x 100,000', line_color = 7, title_no = "(d)")


# 5.1 Three subgroups monthly patterns plot for ACP (Appendix ref) --------
plot(age65$rate, type = "n", ylim = c(0,500), xlab = "Year",  font.lab = 2,
     ylab = "Monthly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
rect(60,0,96,500, col = grey(0.9), border = F)
axis(1, at = 0:14*12, labels = F, lwd = 3)
axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, lwd = 3, font = 4, las = 2)
axis(2, lwd = 3, las = 1, font = 4)
par(adj = 0.02)
title("(a)", font.lab = 1)
points(age019$rate, cex = 0.7, col = 'royalblue', pch = 0)
lines(age019$rate, cex = 0.7, col = 'royalblue', lwd = 2)
points(age2064$rate, cex = 0.7, pch = 1, col = 'darkorange')
lines(age2064$rate, cex = 0.7, col = 'darkorange', lwd = 2)
points(age65$rate, cex = 0.7, pch = 2, col = 'deeppink')
lines(age65$rate, cex = 0.7, col = 'deeppink', lwd = 2)
legend(merge =  T, "topright", c("0-19y", "20-64y", "65+y"),
       lty = 1, lwd = 2, pch = (c(0:2)), col = c('royalblue','darkorange','deeppink'),
       inset = 0.0001, bty = "n", ncol = 3, text.font = 2, pt.cex = 1,
       x.intersp = 0.15, seg.len = 1.5, text.width = 15)


# 5.2 Three subgroups monthly patterns plot for PP (Appendix ref) ---------  
plot(ppage65$rate, type = "n", ylim = c(0,12), xlab = "Year",  font.lab = 2,
     ylab = "Monthly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
rect(60,0,96,12, col = grey(0.9), border = F)
axis(1, at = 0:14*12, labels = F, lwd = 3)
axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, lwd = 3, font = 4, las = 2)
axis(2, lwd = 3, las = 1, font = 4)
par(adj = 0.02)
title("(b)", font.lab = 1)
points(ppage019$rate, cex = 0.7, col = 'royalblue', pch = 0)
lines(ppage019$rate, cex = 0.7, col = 'royalblue', lwd = 2)
points(ppage2064$rate, cex = 0.7, pch = 1, col = 'darkorange')
lines(ppage2064$rate, cex = 0.7, col = 'darkorange', lwd = 2)
points(ppage65$rate, cex = 0.7, pch = 2, col = 'deeppink')
lines(ppage65$rate, cex = 0.7, col = 'deeppink', lwd = 2)
legend(merge =  T, "topright", c("0-19y", "20-64y", "65+y"),
       lty = 1, lwd = 2, pch = (c(0:2)), col = c('royalblue','darkorange','deeppink'),
       inset = 0.0001, bty = "n", ncol = 3, text.font = 2, pt.cex = 1,
       x.intersp = 0.15, seg.len = 1.5, text.width = 15)


# 5.3 Three subgroups monthly patterns plot for Fracture (Appendix ref) ----
plot(ftage65$rate, type = "n", ylim = c(0,150), xlab = "Year",  font.lab = 2,
     ylab = "Monthly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
rect(60,0,96,150, col = grey(0.9), border = F)
axis(1, at = 0:14*12, labels = F, lwd = 3)
axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, lwd = 3, font = 4, las = 2)
axis(2, lwd = 3, las = 1, font = 4)
par(adj = 0.02)
title("(c)", font.lab = 1)
points(ftage019$rate, cex = 0.7, col = 'royalblue', pch = 0)
lines(ftage019$rate, cex = 0.7, col = 'royalblue', lwd = 2)
points(ftage2064$rate, cex = 0.7, pch = 1, col = 'darkorange')
lines(ftage2064$rate, cex = 0.7, col = 'darkorange', lwd = 2)
points(ftage65$rate, cex = 0.7, pch = 2, col = 'deeppink')
lines(ftage65$rate, cex = 0.7, col = 'deeppink', lwd = 2)

legend(merge =  T, "topright", c("0-19y", "20-64y", "65+y"),
       lty = 1, lwd = 2, pch = (c(0:2)), col = c('royalblue','darkorange','deeppink'),
       inset = 0.0001, bty = "n", ncol = 3, text.font = 2, pt.cex = 1,
       x.intersp = 0.15, seg.len = 1.5, text.width = 15)


# -------------------------------------------------------------------------
# Section VI. Poisson model fitting ---------------------------------------
# 0.0 New data set for prediction use -------------------------------------
datanew <- data.frame(pop = mean(pcv$pop), pcv = rep(c(0,NA,1),c(60,36,72)), 
                      time = 0:167, timeafter = pcv$timeafter, month = rep(1:12,14), flu_ad = pcv$flu_ad,
                      pcv3 = rep(c(0,NA,1),c(60,39,69)), timeafter3 = pcv$timeafter3, 
                      pcv6 = rep(c(0,NA,1),c(60,42,66)), timeafter6 = pcv$timeafter6,
                      pcv12 = rep(c(0,NA,1),c(60,48,60)), timeafter12 = pcv$timeafter12)
datanew


# 0.1 Function to check Poisson model assumption (residual & ACF/PACF) ----
res_plot <- function(model){
  res <- residuals(model, type = "deviance")
  res_acf <- 
    list(res_plot = plot(res, ylim = c(-10,10), pch = 16, cex = 0.7, col = grey(0.4),
                         main = "Residuals over time", ylab = "Deviance residuals", xlab = "Time") %>%
           abline(h = 0, lty = 2, lwd = 2),
         acf_plot = acf(res),
         pacf_plo = pacf(res))
}


# 0.2 Function to fit the model with Newey-West method adjustment ---------
poisson_model_result <- function(tempdf, dataname = "", time_lag = "0m lag",
                                 intervention = tempdf$pcv, t_after = tempdf$timeafter, auto_lag = 3, adj_flu = 1){
  # get temp model, results
  if(adj_flu == 1){
    tempformula <- formula(cal_case ~ offset(log(pop)) + time + intervention + t_after:intervention + flu_ad + factor(month))
  } else {
    tempformula <- formula(cal_case ~ offset(log(pop)) + time + intervention + t_after:intervention + factor(month))
  }
  tempmodel <- glm(tempformula, family = poisson, data = tempdf)
  
  # res <- residuals(tempmodel, type = "deviance")
  # acf(res)
  # pacf(res)
  res_plot(model = tempmodel)
  
  # Newey-West standard errors used to account for autocorrelation up to lag 3 (package: sandwich)
  NeweyWest(tempmodel, lag = NULL, verbose = TRUE) # lag truncation parameter chosen: 3
  bw <- bwNeweyWest(tempmodel) # if set lag = NULL, floor(bwNeweyWest(x)) is used

  est <- exp(c(coef(tempmodel)["time"], coef(tempmodel)["intervention"], coef(tempmodel)["intervention:t_after"]))
  se <- c(sqrt(diag(NeweyWest(tempmodel, prewhite = F, lag = auto_lag)))["time"],
          sqrt(diag(NeweyWest(tempmodel, prewhite = F, lag = auto_lag)))["intervention"],
          sqrt(diag(NeweyWest(tempmodel, prewhite = F, lag = auto_lag)))["intervention:t_after"])
  
  lb <- est * exp(-1.96 * se)
  ub <- est * exp(1.96 * se)
  
  table <- cbind(round(est, digits = 6), round(lb, digits = 6), round(ub, digits = 6), round(se, digits = 6))
  p_val <- coeftest(tempmodel, vcov. = NeweyWest(tempmodel, prewhite = F, lag = auto_lag)) # Test SE get P-value
  
  list(Group = paste0("The group dataset is ", dataname, " (", time_lag, ")"), BW = bw, Table = table, Summary = p_val)
}


# 1.1 Main analysis for ACP adjusting flu and seasonality (month as factor) ----
poisson_model_result(tempdf = pcv, dataname = "pcv", time_lag = "0m lag", 
                     intervention = pcv$pcv, t_after = pcv$timeafter, auto_lag = 3, adj_flu = 1)


# 1.2 Predict data set for ACP model and draw an estimated plot ----
pred_poisson_acp <- predict(glm(cal_case ~ offset(log(pop)) + time + pcv + timeafter:pcv + flu_ad + factor(month),
                                family = poisson, data = pcv), 
                            type = "response",
                            transform(datanew, month = 5))/mean(pcv$pop)*10^5

plot(pcv$rate, type = "n", ylim = c(0,100), xlab = "Year", font.lab = 2,
     ylab = "Std monthly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
rect(60, 0, 96, 100, col = grey(0.9), border = F)
abline(v = 60, lwd = 3, pch = 2)
abline(v = 96, lwd = 3, pch = 2)
points(pcv$rate, cex = 0.7, col = 'orange')
axis(1, at = 0:14*12, labels = F, lwd = 3)
axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, lwd = 3, font = 4, las = 2)
axis(2, lwd = 3, las = 1, font = 4)
lines(1:168, pred_poisson_acp, col = 3, lwd = 3)
par(adj = 0.02)
title("(a)", font.lab = 1)


# 2.1 Main analysis for PP adjusting flu and seasonality (month as factor) ----
poisson_model_result(tempdf = pp_pcv, dataname = "pp_pcv", time_lag = "0m lag", 
                     intervention = pp_pcv$pcv, t_after = pp_pcv$timeafter, auto_lag = 1, adj_flu = 1)


# 2.2 Predict data set for PP model and draw an estimated plot ----
pred_poisson_pp <- predict(glm(cal_case ~ offset(log(pop)) + time + pcv + timeafter:pcv + flu_ad + factor(month),
                               family = poisson, data = pp_pcv),
                           type = "response",
                           transform(datanew, month = 4))/mean(pp_pcv$pop)*10^5

plot(pp_pcv$rate, type = "n", ylim = c(0,max(pp_pcv$rate*4/3)), xlab = "Year", font.lab = 2,
     ylab = "Std monthly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
rect(60, 0, 96, max(pp_pcv$rate*4/3), col = grey(0.9), border = F)
abline(v = 60, lwd = 3, pch = 2)
abline(v = 96, lwd = 3, pch = 2)
points(pp_pcv$rate, cex = 0.7, col = 'orange')
axis(1, at = 0:14*12, labels = F, lwd = 3)
axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, lwd = 3, font = 4, las = 2)
axis(2, lwd = 3, las = 1, font = 4)
lines(1:168, pred_poisson_pp, col = 6, lwd = 3)
par(adj = 0.02)
title("(b)", font.lab = 1)


# 3.1 Main analysis for Fracture using month adjustment (without flu) ---- 
poisson_model_result(tempdf = ftdt, dataname = "ftdt", time_lag = "0m lag", 
                     intervention = ftdt$pcv, t_after = ftdt$timeafter, auto_lag = 15, adj_flu = 0)


# 3.2 Predict data set for Fracture model and draw an estimated plot (with ACP together) ----
pred_poisson_ft <- predict(glm(cal_case ~ offset(log(pop)) + time + pcv + timeafter:pcv + factor(month),
                               family = poisson, data = ftdt),
                           type = "response",
                           transform(datanew, month = 3))/mean(ftdt$pop)*10^5

plot(ftdt$rate, type = "n", ylim = c(0,100), xlab = "Year", font.lab = 2,
     ylab = "Std monthly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
rect(60, 0, 96, 100, col = grey(0.9), border = F)
abline(v = 60, lwd = 3, pch = 2)
abline(v = 96, lwd = 3, pch = 2)
points(ftdt$rate, cex = 0.7, col = 'grey')
points(pcv$rate, cex = 0.7, col = 'orange')
axis(1, at = 0:14*12, labels = F, lwd = 3)
axis(1, at = 0:13*12+6, tick = F, labels = 2004:2017, lwd = 3, font = 4, las = 2)
axis(2, lwd = 3, las = 1, font = 4)
lines(1:168, pred_poisson_ft, col = 'blue', lwd = 3)
lines(1:168, pred_poisson_acp, col = 'green', lwd = 3)
par(adj = 0.02)
title("(a)", font.lab = 1)
legend("topright", c("All-cause pneumonia", "Limb fracture"),
       lty = 1, col = c('green','blue'), inset = -0.0005, bty = "n", cex = 0.7)


# 4.1 Main analysis for UTI using month adjustment (without flu) ---- 
poisson_model_result(tempdf = uti, dataname = "uti", time_lag = "0m lag", 
                     intervention = uti$pcv, t_after = uti$timeafter, auto_lag = 4, adj_flu = 0)


# -------------------------------------------------------------------------
# Section VII. Subgroup analyses ------------------------------------------
# 1.1 Subgroup analysis for ACP (Primary outcome - 3 groups) --------------
poisson_model_result(tempdf = age019, dataname = "age019", time_lag = "0m lag", intervention = age019$pcv, t_after = age019$timeafter, auto_lag = 1)
poisson_model_result(tempdf = age019, dataname = "age019", time_lag = "3m lag", intervention = age019$pcv3, t_after = age019$timeafter3, auto_lag = 1)
poisson_model_result(tempdf = age019, dataname = "age019", time_lag = "6m lag", intervention = age019$pcv6, t_after = age019$timeafter6, auto_lag = 1)
poisson_model_result(tempdf = age019, dataname = "age019", time_lag = "12m lag", intervention = age019$pcv12, t_after = age019$timeafter12, auto_lag = 1)

poisson_model_result(tempdf = age2064, dataname = "age2064", time_lag = "0m lag", intervention = age2064$pcv, t_after = age2064$timeafter, auto_lag = 3)
poisson_model_result(tempdf = age2064, dataname = "age2064", time_lag = "3m lag", intervention = age2064$pcv3, t_after = age2064$timeafter3, auto_lag = 3)
poisson_model_result(tempdf = age2064, dataname = "age2064", time_lag = "6m lag", intervention = age2064$pcv6, t_after = age2064$timeafter6, auto_lag = 3)
poisson_model_result(tempdf = age2064, dataname = "age2064", time_lag = "12m lag", intervention = age2064$pcv12, t_after = age2064$timeafter12, auto_lag = 3)

poisson_model_result(tempdf = age65, dataname = "age65", time_lag = "0m lag", intervention = age65$pcv, t_after = age65$timeafter, auto_lag = 1)
poisson_model_result(tempdf = age65, dataname = "age65", time_lag = "3m lag", intervention = age65$pcv3, t_after = age65$timeafter3, auto_lag = 0)
poisson_model_result(tempdf = age65, dataname = "age65", time_lag = "6m lag", intervention = age65$pcv6, t_after = age65$timeafter6, auto_lag = 1)
poisson_model_result(tempdf = age65, dataname = "age65", time_lag = "12m lag", intervention = age65$pcv12, t_after = age65$timeafter12, auto_lag = 1)


# 1.2 Subgroup analysis for ACP (Primary outcome - 7 groups) --------------
poisson_model_result(tempdf = age04, dataname = "age04", auto_lag = 4)
poisson_model_result(tempdf = age519, dataname = "age519", auto_lag = 2)
poisson_model_result(tempdf = age2049, dataname = "age2049", auto_lag = 3)
poisson_model_result(tempdf = age5064, dataname = "age5064", auto_lag = 1)
poisson_model_result(tempdf = age6574, dataname = "age6574", auto_lag = 6)
poisson_model_result(tempdf = age7584, dataname = "age7584", auto_lag = 3)
poisson_model_result(tempdf = age85, dataname = "age85", auto_lag = 4)


# 2.1 Subgroup analysis for PP (Secondary outcome - 3 groups) --------------
poisson_model_result(tempdf = ppage019, dataname = "ppage019", time_lag = "0m lag", intervention = ppage019$pcv, t_after = ppage019$timeafter, auto_lag = 8)
poisson_model_result(tempdf = ppage2064, dataname = "ppage2064", time_lag = "0m lag", intervention = ppage2064$pcv, t_after = ppage2064$timeafter, auto_lag = 1)
poisson_model_result(tempdf = ppage65, dataname = "ppage65", time_lag = "0m lag", intervention = ppage65$pcv, t_after = ppage65$timeafter, auto_lag = 3)


# 2.2 Subgroup analysis for PP (Primary outcome - 7 groups) ---------------
poisson_model_result(tempdf = ppage04, dataname = "ppage04", auto_lag = 12)
poisson_model_result(tempdf = ppage519, dataname = "ppage519", auto_lag = 9)
poisson_model_result(tempdf = ppage2049, dataname = "ppage2049", auto_lag = 0)
poisson_model_result(tempdf = ppage5064, dataname = "ppage5064", auto_lag = 2)
poisson_model_result(tempdf = ppage6574, dataname = "ppage6574", auto_lag = 4)
poisson_model_result(tempdf = ppage7584, dataname = "ppage7584", auto_lag = 5)
poisson_model_result(tempdf = ppage85, dataname = "ppage85", auto_lag = 2)


# 3. Subgroup analysis for Fracture (Negative outcome 1 - 3 groups) -------
poisson_model_result(tempdf = ftage019, dataname = "ftage019", time_lag = "0m lag", 
                     intervention = ftage019$pcv, t_after = ftage019$timeafter, auto_lag = 3, adj_flu = 0)
poisson_model_result(tempdf = ftage2064, dataname = "ftage2064", time_lag = "0m lag", 
                     intervention = ftage2064$pcv, t_after = ftage2064$timeafter, auto_lag = 8, adj_flu = 0)
poisson_model_result(tempdf = ftage65, dataname = "ftage65", time_lag = "0m lag", 
                     intervention = ftage65$pcv, t_after = ftage65$timeafter, auto_lag = 5, adj_flu = 0)


# 4. Subgroup analysis for UTI (Negative outcome 2 - 3 groups) ------------
poisson_model_result(tempdf = utiage019, dataname = "utiage019", time_lag = "0m lag", 
                     intervention = utiage019$pcv, t_after = utiage019$timeafter, auto_lag = 5, adj_flu = 0)
poisson_model_result(tempdf = utiage2064, dataname = "utiage2064", time_lag = "0m lag", 
                     intervention = utiage2064$pcv, t_after = utiage2064$timeafter, auto_lag = 2, adj_flu = 0)
poisson_model_result(tempdf = utiage65, dataname = "utiage65", time_lag = "0m lag", 
                     intervention = utiage65$pcv, t_after = utiage65$timeafter, auto_lag = 5, adj_flu = 0)


# -------------------------------------------------------------------------
# Section VIII. Sensitivity analyses --------------------------------------
# 1. Delay in post-PCV effect with seasonality adjustment -----------------
poisson_model_result(tempdf = pcv, dataname = "pcv", time_lag = "3m lag", intervention = pcv$pcv3, t_after = pcv$timeafter3, auto_lag = 3)
poisson_model_result(tempdf = pcv, dataname = "pcv", time_lag = "6m lag", intervention = pcv$pcv6, t_after = pcv$timeafter6, auto_lag = 3)
poisson_model_result(tempdf = pcv, dataname = "pcv", time_lag = "12m lag", intervention = pcv$pcv12, t_after = pcv$timeafter12, auto_lag = 2)


# 2. Change time interval for overall ACP (30d/60d) -----------------------
### a) Merge overlapping time interval define 30d/60d as interval between admission date and discharge date
# time interval of 30 days
acpdt_30d <- merge_interval_episode(data = acp, time_interval = 30)
fludt_30d <- merge_interval_episode(data = flu, time_interval = 30)

# time interval of 60 days
acpdt_60d <- merge_interval_episode(data = acp, time_interval = 60)
fludt_60d <- merge_interval_episode(data = flu, time_interval = 60)


### b) Get acp_30d/acp_60d monthly age-standardised incidence, and flu monthly admission cases
# time interval of 30 days
acp_30d_mir <- monthly_std_inc(data = acpdt_30d)
flu_30d_mir <- monthly_std_inc(data = fludt_30d)

# time interval of 60 days
acp_60d_mir <- monthly_std_inc(data = acpdt_60d)
flu_60d_mir <- monthly_std_inc(data = fludt_60d)


### c) Three subgroups for acp_30d/acp_60d (0-19y, 20-64y, 65+y), and for flu 
# time interval of 30 days
acp_30d_mir_age3g <- subgroup3_inc(data = acpdt_30d)
flu_30d_mir_age3g <- subgroup3_inc(data = fludt_30d)

# time interval of 60 days
acp_60d_mir_age3g <- subgroup3_inc(data = acpdt_60d)
flu_60d_mir_age3g <- subgroup3_inc(data = fludt_60d)


### d) Prepare overall data set for acp_30d/acp_60d data with pcv and delay pcv time indicator
# time interval of 30 days
pcv_30d <- time_series_data(data = acp_30d_mir, subgroup = 'all', flu_data = flu_30d_mir)
# time interval of 60 days
pcv_60d <- time_series_data(data = acp_60d_mir, subgroup = 'all', flu_data = flu_60d_mir)


### e) Prepare subgroup data set for pcv_30d/pcv_60d
# time interval of 30 days
age019_30d <- time_series_data(data = acp_30d_mir_age3g, subgroup = '0-19', flu_data = flu_30d_mir_age3g)
age2064_30d <- time_series_data(data = acp_30d_mir_age3g, subgroup = '20-64', flu_data = flu_30d_mir_age3g)
age65_30d <- time_series_data(data = acp_30d_mir_age3g, subgroup = '65+', flu_data = flu_30d_mir_age3g)

# time interval of 60 days
age019_60d <- time_series_data(data = acp_60d_mir_age3g, subgroup = '0-19', flu_data = flu_60d_mir_age3g)
age2064_60d <- time_series_data(data = acp_60d_mir_age3g, subgroup = '20-64', flu_data = flu_60d_mir_age3g)
age65_60d <- time_series_data(data = acp_60d_mir_age3g, subgroup = '65+', flu_data = flu_60d_mir_age3g)


### f) Main analysis adjusting flu episodes and seasonality (month as factor as indicator)
# time interval of 30 days
poisson_model_result(tempdf = pcv_30d, dataname = "pcv", time_lag = "0m lag", 
                     intervention = pcv$pcv, t_after = pcv$timeafter, auto_lag = 2)


# time interval of 60 days
poisson_model_result(tempdf = pcv_60d, dataname = "pcv", time_lag = "0m lag", 
                     intervention = pcv$pcv, t_after = pcv$timeafter, auto_lag = 2)


### g) Subgroup analysis for ACP_30d/60d with autocorrelationa adjustment (Newey-West)
poisson_model_result(age019_30d, auto_lag = 1)
poisson_model_result(age2064_30d, auto_lag = 2)
poisson_model_result(age65_30d, auto_lag = 1)

poisson_model_result(age019_60d, auto_lag = 0)
poisson_model_result(age2064_60d, auto_lag = 2)
poisson_model_result(age65_60d, auto_lag = 1)


# 3. Not include flu as a covariate in the model --------------------------
poisson_model_result(tempdf = pcv, dataname = "pcv", auto_lag = 3, adj_flu = 0)


# -------------------------------------------------------------------------
# Section IX. Invasive pneumococcal disease 2015-2019 data set ------------
# 1. Read IPD data 2015-2019 from CHP ----
ipd_chp <- read.csv('raw_data/chp_ipd_number.csv')
names(ipd_chp)[1] <- 'year'


# 2. Add 2018/2019 HK population ----
totalpop <- agepop_gp %>% group_by(year) %>% summarise(pop = sum(gppop)) %>% 
  mutate(year = as.character(year), pop = as.character(pop)) %>% 
  add_row(year = c('2018', '2019'), pop = c('7451000', '7507400'))


# 3. IPD yearly Crude incidence using raw HK pop ----
ipd_dt <- merge(ipd_chp %>% group_by(year) %>% summarise(ipd_ad = sum(ipd_ad)),
                totalpop %>% filter(year > 2014) %>% mutate(year = as.numeric(year), pop = as.numeric(pop)),
                by = 'year') %>% 
  mutate(rate = ipd_ad/pop*10^5)


# 4. Scatter plot of annual incidence of IPD ----
plot(ipd_dt$rate, type = "n", ylim = c(0,3), xlab = "Year", font.lab = 2,
     ylab = "Yearly incidence x 100,000", bty = "n", xaxt = "n", yaxt = "n")
points(ipd_dt$rate, cex = 0.7)
lines(ipd_dt$rate, cex = 0.7, lwd = 2, col = 8)
axis(1, at = 1:5, labels = F, lwd = 3)
axis(1, at = 1:5, tick = F, labels = 2015:2019, srt = 45, lwd = 3, font = 4, las = 1)
axis(2, lwd = 3, las = 1, font = 4)
title("Annual incidence of invasive pneumococcal disease", font.lab = 1)


# -------------------------------------------------------------------------



# -------------------------------------------------------------------------
# Update on 2 Aug 2021.
# -------------------------------------------------------------------------




