
options(scipen=20)

################################################################
# load and install binary packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("mice", "tidyverse")
ipak(packages)


################################################################
# load cleaned data

setwd("~/Documents/IQSS/water-quality")
load("data_cleaned/cleaned_data.Rdata")


################################################################
# multiple imputation (mice) with x2_1_responsible reference at level 2

# change ref level for x2_1_responsible
dat$x2_1_responsible <- relevel(dat$x2_1_responsible, ref = "2")

imp2 <- mice(dat[, all_impute_vars], m=5)

imp2_df <- complete(imp2, action = "long")
imp2_list <- split(imp2_df, f = imp2_df$.imp)
imp2_list <- lapply(imp2_list, function(x) x[, !(colnames(x) %in% c(".imp", ".id"))])

missing2 <- lapply(imp2_list[[1]][, predictor_vars], function(x) which(is.na(x)))
# missing values for rows 100 and 221 for "proxfuncprotect_nwsc" and "proxfunctprotect_ttc_125nwsc"

save(imp2_list, file = "data_cleaned/mice_imp2_data.Rdata", compress = "gzip")
        
        
################################################################
# multiple imputation (mice) with x2_1_responsible reference at level 8

# change ref level for x2_1_responsible
dat$x2_1_responsible <- relevel(dat$x2_1_responsible, ref = "8")

imp8 <- mice(dat[, all_impute_vars], m=5)

imp8_df <- complete(imp8, action = "long")
imp8_list <- split(imp8_df, f = imp8_df$.imp)
imp8_list <- lapply(imp8_list, function(x) x[, !(colnames(x) %in% c(".imp", ".id"))])

missing8 <- lapply(imp8_list[[1]][, predictor_vars], function(x) which(is.na(x)))
# missing values for rows 100 and 221 for "proxfuncprotect_nwsc" and "proxfunctprotect_ttc_125nwsc"

save(imp8_list, file = "data_cleaned/mice_imp8_data.Rdata", compress = "gzip")

