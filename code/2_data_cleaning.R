
options(scipen=20)

################################################################
# load and install binary packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("janitor", "tidyverse")
ipak(packages)


################################################################
# load and clean data

setwd("~/Documents/IQSS/water-quality")

dat <- read_csv("data_raw/WHA_7_19_2020.csv")

# sort out column names
dat <- dat %>% clean_names()

# sort out dates
dat <- dat %>%
    mutate(dateofvisit2 = mdy(dat$dateofvisit2),
                  date = dateofvisit2,
                  dayofyear = yday(date)) %>%
    arrange(dayofyear)

# convert strings to factors
facVars <- as.character(quote(c( 
watersourcetype2, # 8
watersourcetype, # 16
division, # 8
parish, # 51
villagetown, # 76
sampletaken, # 2
topography, # 3
landcover_un_developedor_developed, # 2
temp_c, # 60
functional_v1, # 2
knownfunctionalproblem, # 2
x1_1subcounty, # 4
x1_9_positionrespond, # 8
x2_1_responsible, # 10
x2_2_a_caretaker, # 2
x2_2_b_caretakerpaid, # 2
x2_4wuc_admin, # 1
x2_4wu_cfeecollect, # 1
x2_4wu_cfinmgmt, # 1
x2_4wu_ctechnical, # 1
x2_6_meetsched, # 3
x2_9_meetminutes, # 2
x2_10_minutesviewed, # 4
x2_11_wu_cpay80, # 2
x2_12_a_pubmeet, # 2
x2_13_a_meetaccount, # 2
x3_1_a_functional, # 2
x3_1_b_openallday, # 2
x3_1_c_opentimes, # 2
x3_1_d_primaryuse, # 5
x3_4_a_breaklastyear, # 2
x3_4_b_breaklast2weeks, # 2
x3_6_enoughwater, # 2
x4_0_yearbuilt, # 39
x4_3_a_authorizedlist, # 2
x4_3_c_updateduserlist, # 2
x4_3_d_listshown, # 4
x4_4_feecollectionsystem, # 2
x4_5_feeregular, # 2
x4_6_b_percentfeesother, # 0
x4_6_b_percentfeesdontknow, # 0
x4_11_a_feecollectrecords, # 2
x4_11_b_feerecordsshown, # 5
x4_12_reactive, # 5
x4_14_a_finmgmtrecords, # 2
x4_14_b_namenumber, # 8
x4_15_savingsstorage, # 6
x4_17_plan_oand_m, # 2
x4_18_enoughfundsrepair, # 2
x4_19_personrepair, # 2
x4_20_persontrained, # 2
x5_1_a_facilitator, # 2
x5_1_c_whofacilitated, # 4
x5_2_a_od_fdeclaration, # 3
x6_1_irritated, # 3
x6_2_trustworthy, # 3
x6_3_interviewquality, # 3
ttcge1, # 2
ttcge10, # 2
ttcge50, # 2
ttcge100, # 2
ttcge125, # 2
unmanaged # 2
)))[-1]

dat[, facVars] <- lapply(dat[, facVars], factor)
str(dat)

# view factor levels
lapply(dat[, facVars], table)
lapply(dat[, facVars], nlevels)

# missingness
missing <- sapply(dat[, facVars], function(x) (sum(is.na(x)) / length(x)) * 100)
sort(missing, dec=TRUE) %>% round()


####################################################################################
# merge in environmental variables by date

load("data_cleaned/dat_all_day_tempC_avg.Rdata")
load("data_cleaned/dat_all_day_relhumid_avg.Rdata")
load("data_cleaned/dat_all_day_precip2_mm.Rdata")

dat <- merge(dat, dat_all_day_tempC_avg, by = c("date", "dayofyear"), all.x = TRUE)
dat <- merge(dat, dat_all_day_relhumid_avg, by = c("date", "dayofyear"), all.x = TRUE)
dat <- merge(dat, dat_all_day_precip2_mm, by = c("date", "dayofyear"), all.x = TRUE)

str(dat)

####################################################################################
# predictors

# 37 predictor variables from main data, plus 3 from enviromental data
predictor_vars <- as.character(quote(c(   
avgprox3func,                 
avgprox3func_nwsc, 
avgprox3funcprotect,          
avgprox3funcprotect_nwsc, 
landcover_un_developedor_developed,
proxalt,
prox_nwsc, 
proxaltfunc,                  
proxaltfunc_nwsc, 
proxbore,                     
proxfuncbore, 
proxprotect,                  
proxfuncprotect, 
proxfuncprotect_nwsc,                           
prox_tc, 
proxmarsh,                 
proxfunct_ttc_125, 
proxfunctprotect_ttc_125,      
proxfunctprotect_ttc_125nwsc, 
riskofcontamination10,  
topography,
watersourcetype2,
x2_1_responsible,
x2_2_a_caretaker,
x2_3_inspectorvisits,
x3_1_b_openallday,
x3_1_d_primaryuse,
x4_2_regularuserhhs,
x4_3_a_authorizedlist,
x4_4_feecollectionsystem,   
x4_19_personrepair,                
x4_17_plan_oand_m,          
x4_18_enoughfundsrepair,            
x4_20_persontrained, 
x5_0_latrineaccess,
x5_1_a_facilitator,
x5_2_a_od_fdeclaration,
tempC_avg,
precip2_mm,
relhumid_avg
)))[-1]

# missingness
missing <- sapply(dat[, predictor_vars], function(x) (sum(is.na(x)) / length(x)) * 100)
sort(missing, dec=TRUE) %>% round()


####################################################################################
# variables for each model

# outcome model 1: two alternate responses
functionality_response <- c("functionalinferred",  "x3_1_a_functional")      
functionality_vars <- c(functionality_response, predictor_vars)
#
# outcome model 2: dichotomize at: 1) 11+, 2) 101+
ttc_response <- "ttc" 
ttc_vars <- c(ttc_response, predictor_vars)
#
# outcome model 3 
breakdown_response <- "x3_3_a_lengthbreakdown" 
breakdown_predictors <- predictor_vars[!(predictor_vars %in% c("x5_0_latrineaccess", "x5_1_a_facilitator"))]
breakdown_vars <- c(breakdown_response, breakdown_predictors)
# 
# outcome model 4
# actual_fees_response <- "x4_9_actualmonthlyfees" 
# actual_fees_vars <- c(actual_fees_response, predictor_vars)
#
# outcome model 5 
# reported_fees_response <- "reportedfeesdividedbyallusers" 
# reported_fees_vars <- c(reported_fees_response, predictor_vars)
#
# outcome model 6
enough_water_response <- "x3_6_enoughwater" 
enough_water_predictors <- "functionalinferred"
enough_water_vars <- c(enough_water_response, enough_water_predictors, predictor_vars)

# imputation model
all_vars <- Reduce(union, list(predictor_vars, "functionalinferred", "ttc" , "x3_6_enoughwater")) # 44

# save output
save(dat, 
        all_vars,
        predictor_vars,
        file = "cleaned_data.Rdata", 
        compress = "gzip")

