
options(scipen=20)

################################################################
# load and install binary packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("glmnet", "geepack", "mice", "tidyverse")
ipak(packages)


####################################################################################
# set up predictor design matrix
 
design_matrix <- function(dat, xvars = predictor_vars) {
	dat <- dat[-c(100, 221), xvars]
    xfactors <- formula( paste("~ ", paste(colnames(dat[, sapply(dat, is.factor)]), collapse = " + ")))
    xfactors_dm <- model.matrix(xfactors, data = dat)[, -1]
    xnumeric_vars <- colnames(dat[, sapply(dat, is.numeric)])
    xnumeric <- dat[, xnumeric_vars]
    xnumeric_std <- lapply(xnumeric, scale)
    X_mat <- as.matrix(data.frame(xnumeric_std, xfactors_dm))
	return(X_mat)
}

load("data_cleaned/mice_imp2_data.Rdata")
X_mat2_list <- lapply(imp2_list, design_matrix)


####################################################################################
####################################################################################
# 1) Model for functionality (ref 2)

cv_glmmod_func <- list()
best_lambda_func <- list()
best_coefs_df_func <- list()
non_zero_func <- list()

# remove "openallday"
X_mat2_list_func <- X_mat2_list
for (i in 1:length(X_mat2_list_func)) {
	X_mat2_list_func[[i]] <- X_mat2_list_func[[i]][, !(colnames(X_mat2_list_func[[i]]) %in% "x3_1_b_openallday1")]
}

# alpha=1 for lasso only
for (i in seq_along(X_mat2_list_func)) {
	print(i)
    # use CV to select best lambda
    cv_glmmod_func[[i]] <- cv.glmnet(X_mat2_list_func[[i]], y=factor(imp2_list[[i]][-c(100, 221), "functionalinferred"]),
                                     alpha=1, family="binomial", nfolds=10, parallel=FALSE)
    plot(cv_glmmod_func[[i]])
    best_lambda_func[[i]] <- cv_glmmod_func[[i]]$lambda.min
    best_coefs_df_func[[i]] <- coef(cv_glmmod_func[[i]], s = "lambda.min") %>% as.matrix() %>% as.data.frame()
    non_zero_func[[i]] <- best_coefs_df_func[[i]][, "1"] != 0
}

best_coefs_df_func

non_zero_df_func <- as.data.frame(non_zero_func)
colnames(non_zero_df_func) <- paste0("x", 1:5)
row_idx_func <- which(rowSums(non_zero_df_func) == 5)
keep_vars_func <- rownames(best_coefs_df_func[[1]])[row_idx_func] # intersection of non-zero coefs across all imputations

# -------------------------------------------------------------------------
# glm

# combine new design matrix
glm_mat2_func <- list()

for (i in seq_along(X_mat2_list_func)) {
    glm_mat2_func[[i]] <- data.frame(X_mat2_list_func[[i]], functionalinferred=imp2_list[[i]][-c(100, 221), "functionalinferred"])
}

# specification of intersection of non-zero coefs across all imputations
func_form <- formula( paste("functionalinferred ~ ", paste(keep_vars_func[-1], collapse = " + ")))

# fit complete-data model	
functional_models2 <- lapply(glm_mat2_func, function(x) {
	geeglm(func_form, family = poisson(link = "log"), data = x, id = 1:352, corstr = "exchangeable") 
	})

# pool results
functional_models2_pooled <- pool(functional_models2)

# summarize results
functional_results2 <- summary(functional_models2_pooled) 

# calculate RR point and interval estimates
functional_results2 <- functional_results2 %>%
    mutate(RR = exp(estimate),
           RR_lcl = exp(estimate - (1.96*std.error)),
           RR_ucl = exp(estimate + (1.96*std.error))
           )

write.csv(functional_results2, file = "results/functional_results_ref2.csv", row.names = FALSE)


####################################################################################
####################################################################################
# 2) Model for TTC above certain contamination thresholds (dichotomous)

# dichotomize at 11+ and 101+
for (i in seq_along(imp2_list)) {
    imp2_list[[i]][, "ttc_11"] <- factor(ifelse(imp2_list[[i]][, "ttc"] >= 11, "high", "low"), levels = c("low", "high"))
    imp2_list[[i]][, "ttc_101"] <- factor(ifelse(imp2_list[[i]][, "ttc"] >= 101, "high", "low"), levels = c("low", "high"))
}

####################################################################################
# ttc11 (ref 2)

cv_glmmod_ttc11 <- list()
best_lambda_ttc11 <- list()
best_coefs_df_ttc11 <- list()
non_zero_ttc11 <- list()

# alpha=1 for lasso only
for (i in seq_along(X_mat2_list)) {
	print(i)
    # use CV to select best lambda
    cv_glmmod_ttc11[[i]] <- cv.glmnet(X_mat2_list[[i]], y=factor(imp2_list[[i]][-c(100, 221), "ttc_11"]),
                                      alpha=1, family="binomial", nfolds=10, parallel=FALSE)
    plot(cv_glmmod_ttc11[[i]])
    best_lambda_ttc11[[i]] <- cv_glmmod_ttc11[[i]]$lambda.min
    best_coefs_df_ttc11[[i]] <- coef(cv_glmmod_ttc11[[i]], s = "lambda.min") %>% as.matrix() %>% as.data.frame()
    non_zero_ttc11[[i]] <- best_coefs_df_ttc11[[i]][, "1"] != 0
}

best_coefs_df_ttc11

non_zero_df_ttc11 <- as.data.frame(non_zero_ttc11)
colnames(non_zero_df_ttc11) <- paste0("x", 1:5)
row_idx_ttc11 <- which(rowSums(non_zero_df_ttc11) == 5)
keep_vars_ttc11 <- rownames(best_coefs_df_ttc11[[1]])[row_idx_ttc11] # intersection of non-zero coefs across all imputations

# -------------------------------------------------------------------------
# glm

# combine new design matrix
glm_mat2_ttc11 <- list()

for (i in seq_along(X_mat2_list)) {
    glm_mat2_ttc11[[i]] <- data.frame(X_mat2_list[[i]], ttc_11=as.numeric(imp2_list[[i]][-c(100, 221), "ttc_11"])-1)
}

# specification of intersection of non-zero coefs across all imputations
ttc11_form <- formula( paste("ttc_11 ~ ", paste(keep_vars_ttc11[-1], collapse = " + ")))

# fit complete-data model
ttc11_models2 <- lapply(glm_mat2_ttc11, function(x) {
	geeglm(ttc11_form, family = poisson(link = "log"), data = x, id = 1:352, corstr = "exchangeable")  
	})

# pool results
ttc11_models_pooled2 <- pool(ttc11_models2)

# summarize results
ttc11_results2 <- summary(ttc11_models_pooled2) 

# calculate RR point and interval estimates
ttc11_results2 <- ttc11_results2 %>%
    mutate(RR = exp(estimate),
           RR_lcl = exp(estimate - (1.96*std.error)),
           RR_ucl = exp(estimate + (1.96*std.error))
           )

write.csv(ttc11_results2, file = "results/ttc11_results_ref2.csv", row.names = FALSE)


####################################################################################
# ttc101 (ref 2)

cv_glmmod_ttc101 <- list()
best_lambda_ttc101 <- list()
best_coefs_df_ttc101 <- list()
non_zero_ttc101 <- list()

# alpha=1 for lasso only
for (i in seq_along(X_mat2_list)) {
	print(i)
    # use CV to select best lambda
    cv_glmmod_ttc101[[i]] <- cv.glmnet(X_mat2_list[[i]], y=factor(imp2_list[[i]][-c(100, 221), "ttc_101"]),
                                       alpha=1, family="binomial", nfolds=10, parallel=FALSE)
    plot(cv_glmmod_ttc101[[i]])
    best_lambda_ttc101[[i]] <- cv_glmmod_ttc101[[i]]$lambda.min
    best_coefs_df_ttc101[[i]] <- coef(cv_glmmod_ttc101[[i]], s = "lambda.min") %>% as.matrix() %>% as.data.frame()
    non_zero_ttc101[[i]] <- best_coefs_df_ttc101[[i]][, "1"] != 0
}

best_coefs_df_ttc101

non_zero_df_ttc101 <- as.data.frame(non_zero_ttc101)
colnames(non_zero_df_ttc101) <- paste0("x", 1:5)
row_idx_ttc101 <- which(rowSums(non_zero_df_ttc101) == 5)
keep_vars_ttc101 <- rownames(best_coefs_df_ttc101[[1]])[row_idx_ttc101] # intersection of non-zero coefs across all imputations

# -------------------------------------------------------------------------
# glm

# combine new design matrix
glm_mat2_ttc101 <- list()

for (i in seq_along(X_mat2_list)) {
    glm_mat2_ttc101[[i]] <- data.frame(X_mat2_list[[i]], ttc_101=as.numeric(imp2_list[[i]][-c(100, 221), "ttc_101"])-1)
}

# specification of intersection of non-zero coefs across all imputations
ttc101_form <- formula( paste("ttc_101 ~ ", paste(keep_vars_ttc101[-1], collapse = " + ")))

# fit complete-data model
ttc101_models2 <- lapply(glm_mat2_ttc101, function(x) {
	geeglm(ttc101_form, family = poisson(link = "log"), data = x, id = 1:352, corstr = "exchangeable")   
	})

# pool results
ttc101_models_pooled2 <- pool(ttc101_models2)

# summarize results
ttc101_results2 <- summary(ttc101_models_pooled2) 

# calculate RR point and interval estimates
ttc101_results2 <- ttc101_results2 %>%
    mutate(RR = exp(estimate),
                  RR_lcl = exp(estimate - (1.96*std.error)),
                  RR_ucl = exp(estimate + (1.96*std.error))
                  )

write.csv(ttc101_results2, file = "results/ttc101_results_ref2.csv", row.names = FALSE)


####################################################################################
####################################################################################
# 6) Model for "x3_6_enoughwater": whether water source produces enough water to meet community needs year round
# ref 2

cv_glmmod_enoughwater <- list()
best_lambda_enoughwater <- list()
best_coefs_df_enoughwater <- list()
non_zero_enoughwater <- list()

# alpha=1 for lasso only
for (i in seq_along(X_mat2_list)) {
	print(i)
    # use CV to select best lambda
    cv_glmmod_enoughwater[[i]] <- cv.glmnet(X_mat2_list[[i]], y=factor(imp2_list[[i]][-c(100, 221), "x3_6_enoughwater"]),
                                            alpha=1, family="binomial", nfolds=10, parallel=FALSE)
    plot(cv_glmmod_enoughwater[[i]])
    best_lambda_enoughwater[[i]] <- cv_glmmod_enoughwater[[i]]$lambda.min
    best_coefs_df_enoughwater[[i]] <- coef(cv_glmmod_enoughwater[[i]], s = "lambda.min") %>% as.matrix() %>% as.data.frame()
    non_zero_enoughwater[[i]] <- best_coefs_df_enoughwater[[i]][, "1"] != 0
}

best_coefs_df_enoughwater

non_zero_df_enoughwater <- as.data.frame(non_zero_enoughwater)
colnames(non_zero_df_enoughwater) <- paste0("x", 1:5)
row_idx_enoughwater <- which(rowSums(non_zero_df_enoughwater) == 5)
keep_vars_enoughwater <- rownames(best_coefs_df_enoughwater[[1]])[row_idx_enoughwater] # intersection of non-zero coefs across all imputations

# -------------------------------------------------------------------------
# glm

# combine new design matrix
glm_mat2_enoughwater <- list()

for (i in seq_along(X_mat2_list)) {
    glm_mat2_enoughwater[[i]] <- data.frame(X_mat2_list[[i]], x3_6_enoughwater=as.numeric(imp2_list[[i]][-c(100, 221), "x3_6_enoughwater"])-1)
}

# specification of intersection of non-zero coefs across all imputations
enoughwater_form <- formula( paste("x3_6_enoughwater ~ ", paste(keep_vars_enoughwater[-1], collapse = " + ")))

# fit complete-data model
enoughwater_models2 <- lapply(glm_mat2_enoughwater, function(x) {
	geeglm(enoughwater_form, family = poisson(link = "log"), data = x, id = 1:352, corstr = "exchangeable")  
	})

# pool results
enoughwater_models_pooled2 <- pool(enoughwater_models2)

# summarize results
enoughwater_results2 <- summary(enoughwater_models_pooled2) 

# calculate RR point and interval estimates
enoughwater_results2 <- enoughwater_results2 %>%
    mutate(RR = exp(estimate),
                  RR_lcl = exp(estimate - (1.96*std.error)),
                  RR_ucl = exp(estimate + (1.96*std.error))
                  )

write.csv(enoughwater_results2, file = "results/enoughwater_results_ref2.csv", row.names = FALSE)


####################################################################################
####################################################################################
# 3) Model for length of last breakdown - ref 2

# take out zeros (i.e., no breakdown)

x3_3_a_lengthbreakdown <- dat$x3_3_a_lengthbreakdown[-c(100, 221)]

keep_obs <- which(!is.na(x3_3_a_lengthbreakdown))
x3_3_a_lengthbreakdown_ln <- log(x3_3_a_lengthbreakdown[!is.na(x3_3_a_lengthbreakdown)]+1)

X_mat2_list_sub <- list()

for (i in seq_along(X_mat2_list)) {
	X_mat2_list_sub[[i]] <- X_mat2_list[[i]][keep_obs, !(colnames(X_mat2_list[[i]]) %in% c("x5_1_a_facilitator1", "x5_0_latrineaccess", "x3_1_b_openallday1"))]
}

####################################################################################
cv_glmmod_breakdown <- list()
best_lambda_breakdown <- list()
best_coefs_df_breakdown <- list()
non_zero_breakdown <- list()

# alpha=1 for lasso only
for (i in seq_along(X_mat2_list)) {
	print(i)
    # use CV to select best lambda
    cv_glmmod_breakdown[[i]] <- cv.glmnet(X_mat2_list_sub[[i]], y=x3_3_a_lengthbreakdown_ln,
                                          alpha=1, family="gaussian", nfolds=10, parallel=FALSE)
    plot(cv_glmmod_breakdown[[i]])
    best_lambda_breakdown[[i]] <- cv_glmmod_breakdown[[i]]$lambda.min
    best_coefs_df_breakdown[[i]] <- coef(cv_glmmod_breakdown[[i]], s = "lambda.min") %>% as.matrix() %>% as.data.frame()
    non_zero_breakdown[[i]] <- best_coefs_df_breakdown[[i]][, "1"] != 0
}

best_coefs_df_breakdown

non_zero_df_breakdown <- as.data.frame(non_zero_breakdown)
colnames(non_zero_df_breakdown) <- paste0("x", 1:5)
row_idx_breakdown <- which(rowSums(non_zero_df_breakdown) == 5)
keep_vars_breakdown <- rownames(best_coefs_df_breakdown[[1]])[row_idx_breakdown] # intersection of non-zero coefs across all imputations

# -------------------------------------------------------------------------
# glm

# combine new design matrix
glm_mat2_breakdown <- list()

for (i in seq_along(X_mat2_list)) {
    glm_mat2_breakdown[[i]] <- data.frame(X_mat2_list_sub[[i]], x3_3_a_lengthbreakdown_ln=x3_3_a_lengthbreakdown_ln)
}

# specification of intersection of non-zero coefs across all imputations
breakdown_form <- formula( paste("x3_3_a_lengthbreakdown_ln ~ ", paste(keep_vars_breakdown[-1], collapse = " + ")))

# fit complete-data model
breakdown_models2 <- lapply(glm_mat2_breakdown, function(x) {
	lm(breakdown_form, data = x) 
	})

# pool results
breakdown_models_pooled2 <- pool(breakdown_models2)

# summarize results
breakdown_results2 <- summary(breakdown_models_pooled2) 

# calculate interval estimates
breakdown_results2 <- breakdown_results2 %>%
    mutate(lcl = estimate - (1.96*std.error),
           ucl = estimate + (1.96*std.error)
           )

write.csv(breakdown_results2, file = "results/breakdown_results_ref2.csv", row.names = FALSE)

