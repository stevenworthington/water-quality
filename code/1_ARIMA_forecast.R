
options(scipen=20)

################################################################
# load and install binary packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("parallel", "stringi", "lubridate", "forecast", "data.table", "tidyverse")
ipak(packages)


################################################################
# do computations with multiple processors

# number of cores: 
nc <- detectCores() 
# create clusters 
cl <- makeCluster(rep("localhost", nc))

# -------------------------------------------------------------------------------------------------------
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

# -------------------------------------------------------------------------------------------------------
# plot week
dat_week <- dat %>%
    group_by(month, week) %>%
    summarize_all(mean, na.rm = TRUE)

dat_week_long <- reshape2::melt(dat_week, id.vars = c("time", "month", "week", "day"))

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

dat_day_long <- reshape2::melt(dat_day, id.vars = c("time", "month", "week", "day"))

ggplot(dat_day_long, aes(x = day, y = value, group = 1)) +
    geom_line() +
    geom_smooth() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_classic() 

# -------------------------------------------------------------------------------------------------------
# plot time
dat_time <- reshape2::melt(dat, id.vars = "time")

ggplot(dat_time, aes(x = time, y = value)) +
    geom_path() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_classic() 

################################################################
# ARIMA

# ----------------------------------------
# relhumid_avg

# 864 = 15mins*24hours*9days = 9 days of missing data
# 672 = 15mins*24hours*7days = 1 week
# 96 = 15mins*24hours = 1 day

mseries_relhumid_avg <- msts(dat[, "relhumid_avg"], seasonal.periods = c(96, 672))

z_relhumid_avg <- fourier(mseries_relhumid_avg, K = c(10, 10))
arima_fit_relhumid_avg <- auto.arima(mseries_relhumid_avg, xreg = z_relhumid_avg, ic = "aicc", lambda = "auto", biasadj = TRUE,
    seasonal = FALSE, stepwise = TRUE, stationary = TRUE, approximation = TRUE, parallel = FALSE)
    
zf_relhumid_avg <- fourier(mseries_relhumid_avg, K = c(10, 10), h = 864)    
arima_forecast_relhumid_avg <- forecast(arima_fit_relhumid_avg, xreg = zf_relhumid_avg, h = 864)
plot(arima_forecast_relhumid_avg)

# need to append this to original the get daily average 
newdat_relhumid_avg <- data.frame(  
    time = seq(from = parse_date_time("07-23-2018 00:30", order = "mdy HM"),
               to = parse_date_time("07-31-2018 23:45", order = "mdy HM"), 
               by = "15 min"),
    relhumid_avg = arima_forecast_relhumid_avg$mean[-c(1:2)]
)
 
dat_all_relhumid_avg <- rbind(dat[, c("time", "relhumid_avg")], newdat_relhumid_avg)
    
dat_all_day_relhumid_avg <- dat_all_relhumid_avg %>%
    mutate(day = yday(time),
                   date = date(time),
                   dayofyear = yday(date)) %>%
    group_by(day) %>%
    summarize_all(mean, na.rm = TRUE)

save(dat_all_day_relhumid_avg, file = "data_cleaned/dat_all_day_relhumid_avg.Rdata")

arima_relhumid_avg <- ggplot(dat_all_day_relhumid_avg, aes(x = day, y = relhumid_avg, group = 1)) +
    geom_line() +
    geom_vline(xintercept = 202, linetype = "dashed", color = "darkred")+
    geom_smooth() +
    theme_classic() 
ggsave(arima_relhumid_avg, file = "results/arima_relhumid_avg.pdf", height = 3, width = 7)    


# -------------------------------------------------------------------------------------------------------
# tempC_avg

mseries_tempC_avg <- msts(dat[, "tempC_avg"], seasonal.periods = c(96, 672))

z_tempC_avg <- fourier(mseries_tempC_avg, K = c(10, 10))
arima_fit_tempC_avg <- auto.arima(mseries_tempC_avg, xreg = z_tempC_avg, , ic = "aicc", lambda = "auto", biasadj = TRUE,
    seasonal = FALSE, stepwise = TRUE, stationary = TRUE, approximation = TRUE, parallel = FALSE)

zf_tempC_avg <- fourier(mseries_tempC_avg, K = c(10, 10), h = 864)    
arima_forecast_tempC_avg <- forecast(arima_fit_tempC_avg, xreg = zf_tempC_avg, h = 864)
plot(arima_forecast_tempC_avg)
 
# need to append this to original the get daily average 
newdat_tempC_avg <- data.frame(  
    time = seq(from = parse_date_time("07-23-2018 00:30", order = "mdy HM"),
               to = parse_date_time("07-31-2018 23:45", order = "mdy HM"), 
               by = "15 min"),
    tempC_avg = arima_forecast_tempC_avg$mean[-c(1:2)]
)
 
dat_all_tempC_avg <- rbind(dat[, c("time", "tempC_avg")], newdat_tempC_avg)
    
dat_all_day_tempC_avg <- dat_all_tempC_avg %>%
    mutate(day = yday(time),
                   date = date(time),
                   dayofyear = yday(date)) %>%
    group_by(day) %>%
    summarize_all(mean, na.rm = TRUE)

save(dat_all_day_tempC_avg, file = "data_cleaned/dat_all_day_tempC_avg.Rdata")

arima_tempC_avg <- ggplot(dat_all_day_tempC_avg, aes(x = day, y = tempC_avg, group = 1)) +
    geom_line() +
    geom_vline(xintercept = 202, linetype = "dashed", color = "darkred")+
    geom_smooth() +
    theme_classic() 
ggsave(arima_tempC_avg, file = "results/arima_tempC_avg.pdf", height = 3, width = 7)    


# -------------------------------------------------------------------------------------------------------
# precip2_mm

mseries_precip2_mm <- msts(dat[, "precip2_mm"], seasonal.periods = c(96, 672))

z_precip2_mm <- fourier(mseries_precip2_mm, K = c(10, 10))
arima_fit_precip2_mm <- auto.arima(mseries_precip2_mm, xreg = z_precip2_mm, seasonal = FALSE)

zf_precip2_mm <- fourier(mseries_precip2_mm, K = c(10, 10), h = 864)
arima_forecast_precip2_mm <- forecast(arima_fit_precip2_mm, xreg = zf_precip2_mm, h = 864)
plot(arima_forecast_precip2_mm)

# need to append this to original the get daily average 
newdat_precip2_mm <- data.frame(  
    time = seq(from = parse_date_time("07-23-2018 00:30", order = "mdy HM"),
               to = parse_date_time("07-31-2018 23:45", order = "mdy HM"), 
               by = "15 min"),
    precip2_mm = arima_forecast_precip2_mm$mean[-c(1:2)]
)
 
dat_all_precip2_mm <- rbind(dat[, c("time", "precip2_mm")], newdat_precip2_mm)
    
dat_all_day_precip2_mm <- dat_all_precip2_mm %>%
    mutate(day = yday(time),
                   date = date(time),
                   dayofyear = yday(date)) %>%
    group_by(day) %>%
    summarize_all(mean, na.rm = TRUE)

save(dat_all_day_precip2_mm, file = "data_cleaned/dat_all_day_precip2_mm.Rdata")

arima_precip2_mm <- ggplot(dat_all_day_precip2_mm, aes(x = day, y = precip2_mm, group = 1)) +
    geom_line() +
    geom_vline(xintercept = 202, linetype = "dashed", color = "darkred")+
    geom_smooth() +
    theme_classic() 
ggsave(arima_precip2_mm, file = "results/arima_precip2_mm.pdf", height = 3, width = 7)    

