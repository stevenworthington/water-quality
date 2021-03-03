
options(scipen=20)

################################################################
# load and install binary packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("stringi", "lubridate", "forecast", "data.table", "tibbletime", "tidyverse")
ipak(packages)


################################################################
setwd("~/Documents/IQSS/water-quality")

dat <- fread("data_raw/TAHMOGuluweather2018_15min.csv")
colnames(dat) <- c("time", "precip_mm",  "precip2_mm",  "relhumid_avg",
                               "tempC_avg", "tempC_max", "tempC_min")
# str(dat)

dat <- dat %>%
    mutate(time = parse_date_time(time, order = "mdy HM"),
                month = factor(months(time), levels = month.name),
                week = factor(stringi::stri_datetime_fields(time)$WeekOfMonth, levels = 1:5),
                day = yday(time)) %>%
    filter(time <= "2018-07-22 20:15:00")

setDF(dat)

# str(dat)


################################################################
################################################################
# plots

# -------------------------------------------------------------------------------------------------------
# plot week

dat_week <- dat %>%
    group_by(month, week) %>%
    summarize_all(mean, na.rm = TRUE)

dat_week_long <- dat_week %>% 
    pivot_longer(cols = c("precip_mm", "precip2_mm", "relhumid_avg", "tempC_avg", "tempC_max", "tempC_min"), names_to = "variable")

grid <- expand.grid(month.name, 1:5)
grid <- grid[order(grid$Var1), ]
dat_week_long <- dat_week_long %>%
    mutate(month_week = factor(paste(month, week), levels = paste(grid$Var1, grid$Var2)))

ggplot(dat_week_long, aes(x = month_week, y = value, group = 1)) +
    geom_line() +
    geom_smooth() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, angle = 90))

temp <- ts(dat_week[, "tempC_avg"], 
             freq = 34, 
             start = c(2018, 01, 01),
             end = c(2018, 07, 23)) 
plot(temp)

# -------------------------------------------------------------------------------------------------------
# plot day

dat_day <- dat %>%
    group_by(month, week, day) %>%
    summarize_all(mean, na.rm = TRUE)

dat_day_long <- dat_day %>% 
    pivot_longer(cols = c("precip_mm", "precip2_mm", "relhumid_avg", "tempC_avg", "tempC_max", "tempC_min"), names_to = "variable")
    
ggplot(dat_day_long, aes(x = day, y = value, group = 1)) +
    geom_line() +
    geom_smooth() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_classic() 

# -------------------------------------------------------------------------------------------------------
# plot time

dat_time <- dat %>% 
    pivot_longer(cols = c("precip_mm", "precip2_mm", "relhumid_avg", "tempC_avg", "tempC_max", "tempC_min"), names_to = "variable")
    
ggplot(dat_time, aes(x = time, y = value)) +
    geom_path() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_classic() 


################################################################
################################################################
# ARIMA models

# ----------------------------------------
# relhumid_avg

# 864 = 15mins*24hours*9days = 9 days of missing data
# 672 = 15mins*24hours*7days = 1 week
# 96 = 15mins*24hours = 1 day

# multiseasonal time series data object
mseries_relhumid_avg <- msts(dat[, "relhumid_avg"], seasonal.periods = c(96, 672))

# estimate ARIMA model with seasonality
z_relhumid_avg <- fourier(mseries_relhumid_avg, K = c(10, 10))
arima_fit_relhumid_avg <- auto.arima(mseries_relhumid_avg, xreg = z_relhumid_avg, ic = "aicc", lambda = "auto", biasadj = TRUE,
    seasonal = FALSE, stepwise = TRUE, stationary = TRUE, approximation = TRUE, parallel = FALSE)

# forecast ARIMA model    
zf_relhumid_avg <- fourier(mseries_relhumid_avg, K = c(10, 10), h = 864)    
arima_forecast_relhumid_avg <- forecast(arima_fit_relhumid_avg, xreg = zf_relhumid_avg, h = 864)
plot(arima_forecast_relhumid_avg)

# put forecast in a data frame with time variable 
newdat_relhumid_avg <- data.frame(  
    time = seq(from = parse_date_time("07-23-2018 00:30", order = "mdy HM"),
               to = parse_date_time("07-31-2018 23:45", order = "mdy HM"), 
               by = "15 min"),
    relhumid_avg = arima_forecast_relhumid_avg$mean[-c(1:2)]
)

# need to append forecast to original data to the get daily average 
arima_relhumid_avg <- rbind(dat[, c("time", "relhumid_avg")], newdat_relhumid_avg)

# create a function for the rolling cumulative sum
rollsum_3 <- tibbletime::rollify(sum, window = 3)

# create daily averages and cumulative sums over a lagged 3 day window   
arima_daily_relhumid_avg <- arima_relhumid_avg %>%
    mutate(day = yday(time),
                date = date(time),
                dayofyear = yday(date)) %>%
    group_by(day) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(relhumid_avg_lag = lag(relhumid_avg),
                relhumid_avg_cumsum = rollsum_3(relhumid_avg_lag))

# save
save(arima_daily_relhumid_avg, file = "data_cleaned/arima_daily_relhumid_avg.Rdata")

# plot daily average
plot_arima_relhumid_avg <- ggplot(arima_daily_relhumid_avg, aes(x = day, y = relhumid_avg_cumsum, group = 1)) +
    geom_line() +
    geom_vline(xintercept = 202, linetype = "dashed", color = "darkred")+
    geom_smooth() +
    theme_classic() 
ggsave(plot_arima_relhumid_avg, file = "results/arima_relhumid_avg_cumsum.pdf", height = 3, width = 7)    


# -------------------------------------------------------------------------------------------------------
# tempC_avg

# multiseasonal time series data object
mseries_tempC_avg <- msts(dat[, "tempC_avg"], seasonal.periods = c(96, 672))

# estimate ARIMA model with seasonality
z_tempC_avg <- fourier(mseries_tempC_avg, K = c(10, 10))
arima_fit_tempC_avg <- auto.arima(mseries_tempC_avg, xreg = z_tempC_avg, , ic = "aicc", lambda = "auto", biasadj = TRUE,
    seasonal = FALSE, stepwise = TRUE, stationary = TRUE, approximation = TRUE, parallel = FALSE)

# forecast ARIMA model
zf_tempC_avg <- fourier(mseries_tempC_avg, K = c(10, 10), h = 864)    
arima_forecast_tempC_avg <- forecast(arima_fit_tempC_avg, xreg = zf_tempC_avg, h = 864)
plot(arima_forecast_tempC_avg)
 
# put forecast in a data frame with time variable 
newdat_tempC_avg <- data.frame(  
    time = seq(from = parse_date_time("07-23-2018 00:30", order = "mdy HM"),
               to = parse_date_time("07-31-2018 23:45", order = "mdy HM"), 
               by = "15 min"),
    tempC_avg = arima_forecast_tempC_avg$mean[-c(1:2)]
)

# need to append forecast to original data to the get daily average
arima_tempC_avg <- rbind(dat[, c("time", "tempC_avg")], newdat_tempC_avg)

# create a function for the rolling cumulative sum
rollsum_3 <- tibbletime::rollify(sum, window = 3)

# create daily averages and cumulative sums over a 3 lagged day window    
arima_daily_tempC_avg <- arima_tempC_avg %>%
    mutate(day = yday(time),
                   date = date(time),
                   dayofyear = yday(date)) %>%
    group_by(day) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(tempC_avg_lag = lag(tempC_avg),
                tempC_avg_cumsum = rollsum_3(tempC_avg_lag))

# save
save(arima_daily_tempC_avg, file = "data_cleaned/arima_daily_tempC_avg.Rdata")

# plot daily average
plot_arima_tempC_avg <- ggplot(arima_daily_tempC_avg, aes(x = day, y = tempC_avg_cumsum, group = 1)) +
    geom_line() +
    geom_vline(xintercept = 202, linetype = "dashed", color = "darkred")+
    geom_smooth() +
    theme_classic() 
ggsave(plot_arima_tempC_avg, file = "results/arima_tempC_avg_cumsum.pdf", height = 3, width = 7)    


# -------------------------------------------------------------------------------------------------------
# precip2_mm

# multiseasonal time series data object
mseries_precip2_mm <- msts(dat[, "precip2_mm"], seasonal.periods = c(96, 672))

# estimate ARIMA model with seasonality
z_precip2_mm <- fourier(mseries_precip2_mm, K = c(10, 10))
arima_fit_precip2_mm <- auto.arima(mseries_precip2_mm, xreg = z_precip2_mm, seasonal = FALSE)

# forecast ARIMA model
zf_precip2_mm <- fourier(mseries_precip2_mm, K = c(10, 10), h = 864)
arima_forecast_precip2_mm <- forecast(arima_fit_precip2_mm, xreg = zf_precip2_mm, h = 864)
plot(arima_forecast_precip2_mm)

# put forecast in a data frame with time variable 
newdat_precip2_mm <- data.frame(  
    time = seq(from = parse_date_time("07-23-2018 00:30", order = "mdy HM"),
               to = parse_date_time("07-31-2018 23:45", order = "mdy HM"), 
               by = "15 min"),
    precip2_mm = arima_forecast_precip2_mm$mean[-c(1:2)]
)

# need to append forecast to original data to the get daily average 
arima_precip2_mm <- rbind(dat[, c("time", "precip2_mm")], newdat_precip2_mm)

# create a function for the rolling cumulative sum
rollsum_3 <- tibbletime::rollify(sum, window = 3)

# create daily averages and cumulative sums over a lagged 3 day window    
arima_daily_precip2_mm <- arima_precip2_mm %>%
    mutate(day = yday(time),
                   date = date(time),
                   dayofyear = yday(date)) %>%
    group_by(day) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(precip2_mm_lag = lag(precip2_mm),
                precip2_mm_cumsum = rollsum_3(precip2_mm_lag))

# save
save(arima_daily_precip2_mm, file = "data_cleaned/arima_daily_precip2_mm.Rdata")

# plot daily average
plot_arima_precip2_mm <- ggplot(arima_daily_precip2_mm, aes(x = day, y = precip2_mm_cumsum, group = 1)) +
    geom_line() +
    geom_vline(xintercept = 202, linetype = "dashed", color = "darkred")+
    geom_smooth() +
    theme_classic() 
ggsave(plot_arima_precip2_mm, file = "results/arima_precip2_mm_cumsum.pdf", height = 3, width = 7)    

