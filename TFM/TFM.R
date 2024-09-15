############################# MIBEL data set ###################################
########################### Preliminary details ################################
rm(list=ls())
options(digits=4)
setwd("C:/Users/migue/Documents/MIGUEL/MÁSTER/MASTER DATA SCIENCE/4º SEMICUATRIMESTRE/TFM")
suppressMessages(library(fda))
suppressMessages(library(fda.usc))
suppressMessages(library(ftsa))
suppressMessages(library(forecast))
suppressMessages(library(viridis)) # for colors
color_1 <- "#4686FBFF"

######################### 2. DATABASE CONSTRUCTION #################################
############## Load the data set and plot the time series of prices ############
X_ts <- read.csv(file='Data.csv',header=FALSE,dec=".")
X_ts <- as.matrix(X_ts) 
n_0 <- sum(X_ts==0) # There are 16 missing values
set.seed(1)
X_ts[X_ts==0] <- abs(rnorm(n_0)) # Input missing values as a normal distribution
n_X_ts <- length(X_ts)
n_X_ts # 53712 observations
n_hours <- 24
n_days <- n_X_ts/24
n_days # There are 2238 days
dates_all <- seq(as.Date("2017-01-01"),as.Date("2023-02-16"),by="days") # Year-month-day format
days_plot <- round(seq(from=1,to=length(dates_all),length.out=6))
X_ts_dates <- dates_all[days_plot]

# Missing values
zero_indices <- which(X_ts == 0)
days_with_zero_values <- ceiling(zero_indices / n_hours)
dates_with_zero_values <- dates_all[days_with_zero_values]
dates_with_zero_values

# Plot time series (Figure 2.1)
plot.ts(X_ts,axes=F,xlab="Day",ylab="Price (€/MWh)",col=color_1,lwd=1)
axis(2)
axis(1,labels=X_ts_dates,at=c(1,n_X_ts*(1:5)/5))
box()

################# Obtain the cumulative intraday returns (CIDRs) ###############
X_days <- matrix(X_ts,nrow=n_hours,ncol=n_days) # Rows 24h y columns 2238 days
X_days <- rbind(X_days,c(X_days[1,2:n_days],X_days[n_hours,n_days]+rnorm(1))) # we add another row to make it 0h to 24h, and repeat the next day
X_days_CIDRs <- 100 * (log(X_days) - as.matrix(rep(1,n_hours+1),nrow=n_hours+1,ncol=1) %*% log(X_days[1,])) # first row all 0

######## Plot the CIDRs using colors to distinguish the day of the week ########
# Define the dates and weekdays corresponding to the CIDRs and plot all of them
days_CIDRs <- dates_all
days_CIDRs[1] # first day 2017/01/01
days_CIDRs[n_days] # last day 2023/02/16
week_days_CIDRs <- c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday") # January 1, 2017 was a Sunday
week_days_CIDRs <- rep(week_days_CIDRs,n_days/7+1)
week_days_CIDRs <- week_days_CIDRs[1:n_days]
week_days_CIDRs[1] # first day Sunday
week_days_CIDRs[n_days] # last day Thursday
X_days_CIDRs_colors <- character(n_days)
X_days_CIDRs_colors[week_days_CIDRs=="Sunday"] <- turbo(7)[1]
X_days_CIDRs_colors[week_days_CIDRs=="Monday"] <- turbo(7)[2]
X_days_CIDRs_colors[week_days_CIDRs=="Tuesday"] <- turbo(7)[3]
X_days_CIDRs_colors[week_days_CIDRs=="Wednesday"] <- turbo(7)[4]
X_days_CIDRs_colors[week_days_CIDRs=="Thursday"] <- turbo(7)[5]
X_days_CIDRs_colors[week_days_CIDRs=="Friday"] <- turbo(7)[6]
X_days_CIDRs_colors[week_days_CIDRs=="Saturday"] <- turbo(7)[7]

# Plot (Figure 2.2)
matplot(1:(n_hours+1),axes=F,X_days_CIDRs,type="l",lty=1,lwd=1.8,col=X_days_CIDRs_colors,
        xlab="Hours",ylab="Cumulative Intra-Daily Returns (CIDRs)")
axis(2)
axis(1,labels=0:24,at=1:25)
box()
legend(x = "bottomright", 
       legend = c("Sunday", "Monday", "Tuesday", "Wednesday", 
                  "Thursday", "Friday", "Saturday"),
       col = turbo(7), lwd = 5, cex = 0.5, pt.cex = 0.8) 

# Plot the CIDRs using colors to distinguish the past time (Figure 2.3)
library(viridis)
dates_all <- as.Date(dates_all)  
date_range <- as.numeric(dates_all) 
color_scale <- viridis(n_days, option = "D")  
date_colors <- color_scale[as.numeric(cut(dates_all, breaks = n_days))]

matplot(1:(n_hours+1), axes=F, X_days_CIDRs, type="l", lty=1, lwd=1.8, col=date_colors,
        xlab="Hours", ylab="Cumulative Intra-Daily Returns (CIDRs)")
axis(2)
axis(1, labels=0:24, at=1:25)
box()
legend("bottomright", 
       legend = c("Older Dates", "Recent Dates"),
       fill = c(color_scale[1], color_scale[n_days]),  
       border = "black", cex = 0.8, pt.cex = 0.8) 

################### Smooth the data #########################################
# Define the number of observations per curves and the number of consecutive curves
m <- n_hours+1 # Number of observations per curves (25)
n <- n_days # Number of consecutive curves  (2238 days)

# Define some arguments to define the functional object
argvals_CIDRs <- 0:(m-1)
rangeval_CIDRs <- c(0,m-1)

# Select the number of basis B-SPLINES with GCV
nbasis_CIDRs <- 4:20
l_nbasis_CIDRs <- length(nbasis_CIDRs)
gcv_CIDRs <- 1:l_nbasis_CIDRs
for (i in nbasis_CIDRs){
  basisobj_CIDRs <- create.bspline.basis(rangeval=rangeval_CIDRs,nbasis=i,norder=4)
  CIDRs_smooth <- smooth.basis(argvals=argvals_CIDRs,y=X_days_CIDRs,fdParobj=basisobj_CIDRs)
  gcv_CIDRs[i-3] <- mean(CIDRs_smooth$gcv)
}
# Plot GCV (Figure 2.4)
plot(nbasis_CIDRs,gcv_CIDRs,pch=19,col=color_1,
     main="GCV criterion for OLS",xlab="k",ylab="Value of GCV")  
nbasis_CIDRs <- which.min(gcv_CIDRs)+3 # k=15 the best, we add the +3 because the bspline is cubic
nbasis_CIDRs
# K=15
basisobj_CIDRs <- create.bspline.basis(rangeval=rangeval_CIDRs,nbasis=nbasis_CIDRs,norder=4)
CIDRs_smooth <- smooth.basis(argvals=argvals_CIDRs,y=X_days_CIDRs,fdParobj=basisobj_CIDRs)

# Plot the CIDR corresponding to February, 16, 2023 (last day) (Figure 2.5)
# Graphic 1: With linear interpolation
par(fig=c(0, 1, 0.55, 1), mar=c(2.1, 4.1, 4.1, 2.1))  
matplot(1:(n_hours+1), axes=F, X_days_CIDRs[,n_days], type="l", lty=1, lwd=5,
        col=color_1, ylab="CIDRs", ylim=c(-22, 15))
axis(2) 
axis(1, labels=c(0,4,8,12,16,20,24), at=seq(1, n_hours+1, length.out=7))  
box()
mtext("With linear interpolation", side=3, line=1.5, cex=1.2)  
mtext("Hours", side=1, line=2.5, cex=1.2) 

# Graphic 2: Smoothed with OLS and 15 B-splines
par(fig=c(0, 1, 0, 0.45), new=TRUE, mar=c(4.1, 4.1, 2.1, 2.1))  
grid_smooth <- seq(0, m-1, length.out=138)
matplot(axes=F, eval.fd(grid_smooth, CIDRs_smooth$fd[n_days]), type="l", lty=1, lwd=5,
        col=color_1, ylab="CIDRs curve", ylim=c(-22, 15))
axis(2) 
axis(1, labels=c(0,4,8,12,16,20,24), at=seq(1, 138, length.out=7))  
box()
mtext("Smoothed with OLS and 15 B-splines", side=3, line=1.5, cex=1.2) 
mtext("Hours", side=1, line=2.5, cex=1.2)  

######################## 3. DESCRIPTIVE ANALYSIS OF CIDRs #########################
# Plot all the CIDRs (Figure 3.1)
plot(CIDRs_smooth,lty=1,lwd=2,xlab="Hours",ylab="CIDRs",col=color_1)
title("Smoothed CIDRs with OLS and 15 B-splines")
axis(1, at=seq(0, 24, by=1), labels=seq(0, 24, by=1))

# Stationary test by Horváth, Kokoszka and Rice (2014)
T_stationary(CIDRs_smooth$y)

######################## Sample functional characteristics ####################
# Sample functional mean (Figure 3.2)
color_2 <- "darkorchid2"
mean_CIDRs <- mean.fd(CIDRs_smooth$fd)
plot(CIDRs_smooth,lty=1,lwd=2,xlab="Hours",ylab="CIDRs",col=color_1)
title("Sample functional mean")
lines(mean_CIDRs,col=color_2,lwd=4)
axis(1, at=seq(0, 24, by=1), labels=seq(0, 24, by=1))

# Sample functional sd (Figure 3.3)
sd_CIDRs <- sd.fd(CIDRs_smooth$fd)
plot(sd_CIDRs,lty=1,lwd=2,xlab="Hours",ylab="CIDRs",col=color_1)
title("Sample functional standard deviation")
lines(sd_CIDRs,col=color_2,lwd=4)
axis(1, at=seq(0, 24, by=1), labels=seq(0, 24, by=1))

# Sample functional covariance (Figure 3.4)
cov_CIDRs <- var.fd(CIDRs_smooth$fd)
points_CIDRs <- seq(0,24,length.out=100)
cov_points_CIDRs <- eval.bifd(points_CIDRs,points_CIDRs,cov_CIDRs)
persp(points_CIDRs,points_CIDRs,cov_points_CIDRs,phi=30,theta=30,expand=.5,col=color_1,
      ltheta=120,shade=0.5,ticktype="detailed",xlab="t",ylab="s",zlab="",
      r=40,d=.1,border=color_2,main="Sample covariance function of the CIDRs")
contour(points_CIDRs,points_CIDRs,cov_points_CIDRs,lwd=2,col=color_1,
        main="Contour plot of the sample covariance function of the CIDRs")


############################ Outliers ###########################################
# TIME CONSUMING!!!
library(fda.usc)
tt <- 0:24
CIDRs_smooth_fdausc <- eval.fd(tt,CIDRs_smooth$fd)
fdataobj_CIDRs <- fdata(t(CIDRs_smooth_fdausc),tt)
out_trimming <- outliers.depth.trim(fdataobj_CIDRs)
out_trimming$outliers

######################## 5. DEVELOPMENT OF THE STUDY #############################
################################ 5.1. FPCs ################################
######### Compute functional principal components (FPCs) and FPC scores ########
fpcs_CIDRs_smooth <- pca.fd(CIDRs_smooth$fd,nharm=15,harmfdPar=fdPar(CIDRs_smooth$fd))
fpcs_CIDRs_psi <- fpcs_CIDRs_smooth$harmonics
CIDRs_scores <- fpcs_CIDRs_smooth$scores

####################### Determine the truncation level #########################
# Have a look at eigenvalues and the explained variability
table_fpcs_CIDRs_smooth <- cbind(fpcs_CIDRs_smooth$values[1:15],fpcs_CIDRs_smooth$varprop,cumsum(fpcs_CIDRs_smooth$varprop))
colnames(table_fpcs_CIDRs_smooth) <- c("Eigenvalues","Proportion of variability","Accumulated proportion of variability")
kable(table_fpcs_CIDRs_smooth)

################## Plot sample eigenvalues (Figure 5.1) ######################
par(mfrow=c(1,2), mar=c(4, 4, 2, 1))  
plot(1:15, table_fpcs_CIDRs_smooth[,1], pch=19, col=color_1, type="b", 
     xlab="Eigenvalue number", ylab="Eigenvalue", lwd=3, 
     main="Sample eigenvalues")
plot(1:15, table_fpcs_CIDRs_smooth[,3], pch=19, col=color_1, type="b", 
     xlab="Eigenvalue number", ylab="Eigenvalue", lwd=3, 
     main="Proportion of variability")

############### Plot the first six sample functional eigenfunctions (Figure 5.2) ###############
par(mfrow=c(3,2), mar=c(4.1, 4.1, 2.1, 1.1))
plot(axes=F,eval.fd(grid_smooth,fpcs_CIDRs_smooth$harmonics[1]),type="l",lwd=5,col=color_1,lty=1,
     xlab="Hour",ylab="Value", main="First FCP")
abline(h=0,lwd=3,lty=3)
axis(2)
axis(1,labels=c(0,4,8,12,16,20),at=floor(c(1,5,9,13,17,21)*138/24),cex.axis=1.3)
box()
plot(axes=F,eval.fd(grid_smooth,fpcs_CIDRs_smooth$harmonics[2]),type="l",lwd=5,col=color_1,lty=1,
     xlab="Hour",ylab="Value", main="Second FCP")
abline(h=0,lwd=3,lty=3)
axis(2)
axis(1,labels=c(0,4,8,12,16,20),at=floor(c(1,5,9,13,17,21)*138/24),cex.axis=1.3)
box()
plot(axes=F,eval.fd(grid_smooth,fpcs_CIDRs_smooth$harmonics[3]),type="l",lwd=5,col=color_1,lty=1,
     xlab="Hour",ylab="Value", main="Third FCP")
abline(h=0,lwd=3,lty=3)
axis(2)
axis(1,labels=c(0,4,8,12,16,20),at=floor(c(1,5,9,13,17,21)*138/24),cex.axis=1.3)
box()
plot(axes=F,eval.fd(grid_smooth,fpcs_CIDRs_smooth$harmonics[4]),type="l",lwd=5,col=color_1,lty=1,
     xlab="Hour",ylab="Value", main="Fourth FCP")
abline(h=0,lwd=3,lty=3)
axis(2)
axis(1,labels=c(0,4,8,12,16,20),at=floor(c(1,5,9,13,17,21)*138/24),cex.axis=1.3)
box()
plot(axes=F,eval.fd(grid_smooth,fpcs_CIDRs_smooth$harmonics[5]),type="l",lwd=5,col=color_1,lty=1,
     xlab="Hour",ylab="Value", main="Fifth FCP")
abline(h=0,lwd=3,lty=3)
axis(2)
axis(1,labels=c(0,4,8,12,16,20),at=floor(c(1,5,9,13,17,21)*138/24),cex.axis=1.3)
box()
plot(axes=F,eval.fd(grid_smooth,fpcs_CIDRs_smooth$harmonics[6]),type="l",lwd=5,col=color_1,lty=1,
     xlab="Hour",ylab="Value", main="Sixth FCP")
abline(h=0,lwd=3,lty=3)
axis(2)
axis(1,labels=c(0,4,8,12,16,20),at=floor(c(1,5,9,13,17,21)*138/24),cex.axis=1.3)
box()

############ Plot the first six FPC scores (Figure 6.2)(APPENDIX) #############
scores_plot <- round(seq(from=1,to=length(dates_all),length.out=5))
scores_dates <- dates_all[scores_plot]

par(mfrow=c(3,2), mar=c(4.1, 4.1, 2.1, 1.1))
plot.ts(axes=F,CIDRs_scores[,1],type="l",lwd=1,col=color_1,lty=1,
        xlab="",ylab="", main="FCP1 scores")
axis(2)
axis(1,labels=scores_dates,at=c(1,n*(1:4)/4),cex.axis=1.3)
box()
plot.ts(axes=F,CIDRs_scores[,2],type="l",lwd=1,col=color_1,lty=1,
        xlab="",ylab="", main="FCP2 scores")
axis(2)
axis(1,labels=scores_dates,at=c(1,n*(1:4)/4),cex.axis=1.3)
box()
plot.ts(axes=F,CIDRs_scores[,3],type="l",lwd=1,col=color_1,lty=1,
        xlab="",ylab="", main="FCP3 scores")
axis(2)
axis(1,labels=scores_dates,at=c(1,n*(1:4)/4),cex.axis=1.3)
box()
plot.ts(axes=F,CIDRs_scores[,4],type="l",lwd=1,col=color_1,lty=1,
        xlab="",ylab="", main="FCP4 scores")
axis(2)
axis(1,labels=scores_dates,at=c(1,n*(1:4)/4),cex.axis=1.3)
box()
plot.ts(axes=F,CIDRs_scores[,5],type="l",lwd=1,col=color_1,lty=1,
        xlab="",ylab="", main="FCP5 scores")
axis(2)
axis(1,labels=scores_dates,at=c(1,n*(1:4)/4),cex.axis=1.3)
box()
plot.ts(axes=F,CIDRs_scores[,6],type="l",lwd=1,col=color_1,lty=1,
        xlab="",ylab="", main="FCP6 scores")
axis(2)
axis(1,labels=scores_dates,at=c(1,n*(1:4)/4),cex.axis=1.3)
box()

# Plot ACFs of the first two FPC scores (Figure 5.3)
par(mfrow=c(1,2), mar=c(5, 4, 4, 2) + 0.1, oma=c(0, 0, 2, 0)) 
# Function to plot ACF excluding lag 0
plot_acf_exclude_zero <- function(series, title) {
  acf_result <- acf(series, plot=FALSE)
  # Remove lag 0
  lag_values <- acf_result$lag[-1]
  acf_values <- acf_result$acf[-1]
  # Plot
  plot(lag_values, acf_values, type="h", lwd=3, col=color_1,
       main=title, ylab="Autocorrelation", xlab="Lag", ylim=c(0, 0.5))
  abline(h=0, lwd=3, lty=3)
}

# Plot the ACF for FPC1 scores
plot_acf_exclude_zero(CIDRs_scores[,1], "ACF of the FPC1 scores")
# Plot the ACF for FPC2 scores
plot_acf_exclude_zero(CIDRs_scores[,2], "ACF of the FPC2 scores")



############################## 5.2 Predictive models #######################################
############################ Univariate models #################################
##################### Define hyper-parameters grid search #############################
arima_params_grid <- expand.grid(order_p = c(1, 2, 3), order_d = c(0,1), order_q = c(1, 2, 3), seasonal_order_p = c(0, 1), seasonal_order_d = c(0, 1), seasonal_order_q = c(0, 1), seasonal_period = c(7), include_mean = c(FALSE))
stl_params_grid <- expand.grid(s_window = c("periodic"), s_degree = c(1))
hw_params_grid <- expand.grid(seasonal = c("additive"))
ann_params_grid <- expand.grid(p = c(3, 5, 7, 10, 15), P = c(1), k = c(3, 5, 7, 9))
prophet_params_grid <- expand.grid(dummy = 1)  

# Function to obtain the forecast with ARIMA
forecast_arima <- function(ts_data, params, h=1) {
  model <- arima(ts_data, order=c(params$order_p, params$order_d, params$order_q), seasonal=list(order=c(params$seasonal_order_p, params$seasonal_order_d, params$seasonal_order_q), period=params$seasonal_period), include.mean=params$include_mean)
  forecast(model, h=h)
}

# Function to obtain the forecast with STL
forecast_stl <- function(ts_data, params, h=1) {
  model <- stl(ts_data, s.window=params$s_window, s.degree=params$s_degree)
  return(forecast(model, h=h))
}

# Function to obtain the forecast with HW-Additive
forecast_hw <- function(ts_data, params, h=1) {
  model <- HoltWinters(ts_data, seasonal=params$seasonal)
  forecast(model, h=h)
}

# Function to obtain the forecast with NNAR
forecast_ann <- function(ts_data, params, h=1) {
  model <- nnetar(ts_data, p=params$p, P=params$P, k=params$k)
  forecast(model, h=h)
}

# Function to obtain the forecast with Prophet
library(prophet)
# Exists also library(fable.prophet)
forecast_prophet <- function(ts_data, params, h=1) {
  # Create a dataframe compatible with Prophet
  df <- data.frame(ds = time(ts_data), y = as.numeric(ts_data))
  # Adjust the Prophet model
  model <- prophet(df, growth="linear")
  future <- make_future_dataframe(model, periods=h)
  forecast <- predict(model, future)
  return(forecast[(nrow(forecast)-h+1):nrow(forecast), "yhat"])
}

# Function to obtain the error with MISFE
calculate_error <- function(for_CIDRs, i, CIDRs_smooth, n_for, metric="MISFE") {
  error <- for_CIDRs[i] - CIDRs_smooth$fd[n-n_for+i]
  if (metric == "MISFE") {
    return(norm.fdata(fdata(error, argvals=error$argvals, rangeval=error$rangeval), metric=metric.lp, lp=2))
  } else {
    return(norm.fdata(fdata(error, argvals=error$argvals, rangeval=error$rangeval), metric=metric.lp, lp=1))
  }
}

# Function to perform hyper-parameters search grid and return the best ones
hyperparameter_search <- function(ts_scores, params_grid, forecast_function, CIDRs_smooth, n_for, metric="MISFE") {
  best_error <- Inf
  best_params <- NULL
  for (params in 1:nrow(params_grid)) {
    current_params <- params_grid[params, ]
    message("Evaluando parámetros: ", paste(current_params, collapse=", "))
    for_CIDRs <- CIDRs_smooth$fd[(n-n_for):n]
    for (i in 1:n_for) {
      message("Iteración: ", i)
      tryCatch({
        for_scores <- lapply(ts_scores, function(ts) forecast_function(ts, params=current_params, h=1))
        # Verify not NA
        if (any(sapply(for_scores, function(x) any(is.na(x$mean))))) {
          message("for_scores contiene NA en la iteración ", i, " con parámetros: ", paste(current_params, collapse=", "))
          next
        }
        for_CIDRs$coefs[,i] <- mean.fd(CIDRs_smooth$fd)$coefs + 
          sum(sapply(1:6, function(j) as.numeric(for_scores[[j]]$mean) * fpcs_CIDRs_smooth$harmonics[j]$coefs))
      }, error = function(e) {
        message("Error en la iteración ", i, " con parámetros: ", paste(current_params, collapse=", "))
        message("Mensaje de error: ", e$message)
        next
      })
    }
    current_error <- mean(sapply(1:n_for, function(i) calculate_error(for_CIDRs, i, CIDRs_smooth, n_for, metric)))
    if (!is.na(current_error) && current_error < best_error) {
      best_error <- current_error
      best_params <- current_params
    }
  }
  return(list(best_params=best_params, best_error=best_error))
}


# Define the objects that will contain the 47 predicted curves
n_for <- 47
ts_scores <- lapply(1:6, function(j) {
  ts(CIDRs_scores[1:(n-n_for+i-1), j], start=c(2017, as.numeric(format(dates_all[1], "%j"))), frequency=365)
})

# Perform the hyperparameter search for each model (TIME CONSUMING!!!)
best_arima <- hyperparameter_search(ts_scores, arima_params_grid, forecast_arima, CIDRs_smooth, n_for, metric="MISFE")
best_ann <- hyperparameter_search(ts_scores, ann_params_grid, forecast_ann, CIDRs_smooth, n_for, metric="MISFE")

best_stl <- hyperparameter_search(ts_scores, stl_params_grid, forecast_stl, CIDRs_smooth, n_for, metric="MISFE")
best_hw <- hyperparameter_search(ts_scores, hw_params_grid, forecast_hw, CIDRs_smooth, n_for, metric="MISFE")
best_prophet <- hyperparameter_search(ts_scores, prophet_params_grid, forecast_prophet, CIDRs_smooth, n_for, metric="MISFE")

# Print the best hyper-parameters
print(best_arima) # order_p=2, seasonal_order_p=1 
print(best_stl)
print(best_hw)
print(best_ann) # p=5, P=1, k=5
print(best_prophet)

############################## Apply univariate models #######################################
# Obtain the score time series and forecast the next day 
# There are 6 years and 47 days, we are going to forecast the last 47 days
# Define the initial data set for estimation
# Note that there are 47 days to forecast (31 days in January and 16 days in February) 
# Now, we have to repeat the procedure starting with the computation of the FPCs

library(prophet)
library(forecast)
library(fda)

# Define the objects that will contain the 47 predicted curves
n_for <- 47
for_CIDRs_arima <- for_CIDRs_stl <- for_CIDRs_hw_ad <- for_CIDRs_ann <- for_CIDRs_prophet <- CIDRs_smooth$fd[(n-n_for):n]

for (i in 1 : n_for){
  
  print(i)
  
  # Obtain the FPC scores for the corresponding day
  
  fpcs_CIDRs_smooth <- pca.fd(CIDRs_smooth$fd[1:(n-n_for+i-1)],nharm=6,harmfdPar=fdPar(CIDRs_smooth$fd[1:(n-n_for+i-1)]))
  CIDRs_scores <- fpcs_CIDRs_smooth$scores
  
  ts_score_1 <- ts(CIDRs_scores[1:(n-n_for+i-1),1],start=c(2017,as.numeric(format(dates_all[1],"%j"))),frequency=365)
  ts_score_2 <- ts(CIDRs_scores[1:(n-n_for+i-1),2],start=c(2017,as.numeric(format(dates_all[1],"%j"))),frequency=365)
  ts_score_3 <- ts(CIDRs_scores[1:(n-n_for+i-1),3],start=c(2017,as.numeric(format(dates_all[1],"%j"))),frequency=365)
  ts_score_4 <- ts(CIDRs_scores[1:(n-n_for+i-1),4],start=c(2017,as.numeric(format(dates_all[1],"%j"))),frequency=365)
  ts_score_5 <- ts(CIDRs_scores[1:(n-n_for+i-1),5],start=c(2017,as.numeric(format(dates_all[1],"%j"))),frequency=365)
  ts_score_6 <- ts(CIDRs_scores[1:(n-n_for+i-1),6],start=c(2017,as.numeric(format(dates_all[1],"%j"))),frequency=365)

  # Helper function to prepare data for Prophet
  prepare_prophet_df <- function(ts_data) {
    start_date <- as.Date("2017-01-01")
    time_indices <- seq(start_date, by = "days", length.out = length(ts_data))
    data.frame(ds = time_indices, y = as.numeric(ts_data))
  }
  
  # Obtain the forecasts using ARIMA models
  
  for_score_1_arima <- forecast(arima(ts_score_1,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=7),include.mean=F),h=1)
  for_score_2_arima <- forecast(arima(ts_score_2,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=7),include.mean=F),h=1)
  for_score_3_arima <- forecast(arima(ts_score_3,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=7),include.mean=F),h=1)
  for_score_4_arima <- forecast(arima(ts_score_4,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=7),include.mean=F),h=1)
  for_score_5_arima <- forecast(arima(ts_score_5,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=7),include.mean=F),h=1)
  for_score_6_arima <- forecast(arima(ts_score_6,order=c(2,0,0),seasonal=list(order=c(1,0,0),period=7),include.mean=F),h=1)

  # Obtain the forecasted CIDR with ARIMA
  
  for_CIDRs_arima$coefs[,i] <- mean.fd(CIDRs_smooth$fd)$coefs + 
   as.numeric(for_score_1_arima$mean) * fpcs_CIDRs_smooth$harmonics[1]$coefs +
   as.numeric(for_score_2_arima$mean) * fpcs_CIDRs_smooth$harmonics[2]$coefs +
   as.numeric(for_score_3_arima$mean) * fpcs_CIDRs_smooth$harmonics[3]$coefs +
   as.numeric(for_score_4_arima$mean) * fpcs_CIDRs_smooth$harmonics[4]$coefs + 
   as.numeric(for_score_5_arima$mean) * fpcs_CIDRs_smooth$harmonics[5]$coefs +
   as.numeric(for_score_6_arima$mean) * fpcs_CIDRs_smooth$harmonics[6]$coefs

  # Obtain the forecasts using STL decomposition
  
  for_score_1_stl <- forecast(stl(ts_score_1,s.window="periodic",s.degree=1),h=1)
  for_score_2_stl <- forecast(stl(ts_score_2,s.window="periodic",s.degree=1),h=1)
  for_score_3_stl <- forecast(stl(ts_score_3,s.window="periodic",s.degree=1),h=1)
  for_score_4_stl <- forecast(stl(ts_score_4,s.window="periodic",s.degree=1),h=1)
  for_score_5_stl <- forecast(stl(ts_score_5,s.window="periodic",s.degree=1),h=1)
  for_score_6_stl <- forecast(stl(ts_score_6,s.window="periodic",s.degree=1),h=1)

  # Obtain the forecasted CIDR with STL decomposition
  
  for_CIDRs_stl$coefs[,i] <- mean.fd(CIDRs_smooth$fd)$coefs + 
   as.numeric(for_score_1_stl$mean) * fpcs_CIDRs_smooth$harmonics[1]$coefs +
   as.numeric(for_score_2_stl$mean) * fpcs_CIDRs_smooth$harmonics[2]$coefs +
   as.numeric(for_score_3_stl$mean) * fpcs_CIDRs_smooth$harmonics[3]$coefs +
   as.numeric(for_score_4_stl$mean) * fpcs_CIDRs_smooth$harmonics[4]$coefs + 
   as.numeric(for_score_5_stl$mean) * fpcs_CIDRs_smooth$harmonics[5]$coefs + 
   as.numeric(for_score_6_stl$mean) * fpcs_CIDRs_smooth$harmonics[6]$coefs

  # Obtain the forecast using Holt-Winters with additive seasonality (exponential smoothing)
  
  for_score_1_hw_ad <- forecast(HoltWinters(ts_score_1,seasonal="additive"),h=1)
  for_score_2_hw_ad <- forecast(HoltWinters(ts_score_2,seasonal="additive"),h=1)
  for_score_3_hw_ad <- forecast(HoltWinters(ts_score_3,seasonal="additive"),h=1)
  for_score_4_hw_ad <- forecast(HoltWinters(ts_score_4,seasonal="additive"),h=1)
  for_score_5_hw_ad <- forecast(HoltWinters(ts_score_5,seasonal="additive"),h=1)
  for_score_6_hw_ad <- forecast(HoltWinters(ts_score_6,seasonal="additive"),h=1)

  # Obtain the forecasted CIDR with Holt-Winters with additive seasonality
  
  for_CIDRs_hw_ad$coefs[,i] <- mean.fd(CIDRs_smooth$fd)$coefs + 
   as.numeric(for_score_1_hw_ad$mean) * fpcs_CIDRs_smooth$harmonics[1]$coefs +
   as.numeric(for_score_2_hw_ad$mean) * fpcs_CIDRs_smooth$harmonics[2]$coefs +
   as.numeric(for_score_3_hw_ad$mean) * fpcs_CIDRs_smooth$harmonics[3]$coefs +
   as.numeric(for_score_4_hw_ad$mean) * fpcs_CIDRs_smooth$harmonics[4]$coefs +
   as.numeric(for_score_5_hw_ad$mean) * fpcs_CIDRs_smooth$harmonics[5]$coefs +
   as.numeric(for_score_6_hw_ad$mean) * fpcs_CIDRs_smooth$harmonics[6]$coefs

  # Obtain the forecast using autoregressive neural networks
  
  for_score_1_ann <- forecast(nnetar(ts_score_1,5,1,5),h=1)
  for_score_2_ann <- forecast(nnetar(ts_score_2,5,1,5),h=1)
  for_score_3_ann <- forecast(nnetar(ts_score_3,5,1,5),h=1)
  for_score_4_ann <- forecast(nnetar(ts_score_4,5,1,5),h=1)
  for_score_5_ann <- forecast(nnetar(ts_score_5,5,1,5),h=1)
  for_score_6_ann <- forecast(nnetar(ts_score_6,5,1,5),h=1)

  # Obtain the forecasted CIDR using autoregressive neural networks
  
  for_CIDRs_ann$coefs[,i] <- mean.fd(CIDRs_smooth$fd)$coefs + 
    as.numeric(for_score_1_ann$mean) * fpcs_CIDRs_smooth$harmonics[1]$coefs +
    as.numeric(for_score_2_ann$mean) * fpcs_CIDRs_smooth$harmonics[2]$coefs +
    as.numeric(for_score_3_ann$mean) * fpcs_CIDRs_smooth$harmonics[3]$coefs +
    as.numeric(for_score_4_ann$mean) * fpcs_CIDRs_smooth$harmonics[4]$coefs +
    as.numeric(for_score_5_ann$mean) * fpcs_CIDRs_smooth$harmonics[5]$coefs +
    as.numeric(for_score_6_ann$mean) * fpcs_CIDRs_smooth$harmonics[6]$coefs
  
  # Obtain the forecast using Facebook Prophet
  
  # Prepare data for Prophet with the function previously created
  ts_score_1_df <- prepare_prophet_df(ts_score_1)
  ts_score_2_df <- prepare_prophet_df(ts_score_2)
  ts_score_3_df <- prepare_prophet_df(ts_score_3)
  ts_score_4_df <- prepare_prophet_df(ts_score_4)
  ts_score_5_df <- prepare_prophet_df(ts_score_5)
  ts_score_6_df <- prepare_prophet_df(ts_score_6)
  
  # Adjust and predict with library(prophet)
  model_1 <- prophet::prophet(ts_score_1_df, weekly.seasonality = TRUE)
  model_2 <- prophet::prophet(ts_score_2_df, weekly.seasonality = TRUE)
  model_3 <- prophet::prophet(ts_score_3_df, weekly.seasonality = TRUE)
  model_4 <- prophet::prophet(ts_score_4_df, weekly.seasonality = TRUE)
  model_5 <- prophet::prophet(ts_score_5_df, weekly.seasonality = TRUE)
  model_6 <- prophet::prophet(ts_score_6_df, weekly.seasonality = TRUE)
  
  future_1 <- make_future_dataframe(model_1, periods=1)
  future_2 <- make_future_dataframe(model_2, periods=1)
  future_3 <- make_future_dataframe(model_3, periods=1)
  future_4 <- make_future_dataframe(model_4, periods=1)
  future_5 <- make_future_dataframe(model_5, periods=1)
  future_6 <- make_future_dataframe(model_6, periods=1)
  
  forecast_1 <- predict(model_1, future_1)
  forecast_2 <- predict(model_2, future_2)
  forecast_3 <- predict(model_3, future_3)
  forecast_4 <- predict(model_4, future_4)
  forecast_5 <- predict(model_5, future_5)
  forecast_6 <- predict(model_6, future_6)
  
  for_score_1_prophet <- forecast_1$yhat[nrow(forecast_1)]
  for_score_2_prophet <- forecast_2$yhat[nrow(forecast_2)]
  for_score_3_prophet <- forecast_3$yhat[nrow(forecast_3)]
  for_score_4_prophet <- forecast_4$yhat[nrow(forecast_4)]
  for_score_5_prophet <- forecast_5$yhat[nrow(forecast_5)]
  for_score_6_prophet <- forecast_6$yhat[nrow(forecast_6)]
  
  # Obtain the forecasted CIDR using Prophet
  for_CIDRs_prophet$coefs[, i] <- mean.fd(CIDRs_smooth$fd)$coefs + 
    for_score_1_prophet * fpcs_CIDRs_smooth$harmonics[1]$coefs +
    for_score_2_prophet * fpcs_CIDRs_smooth$harmonics[2]$coefs +
    for_score_3_prophet * fpcs_CIDRs_smooth$harmonics[3]$coefs +
    for_score_4_prophet * fpcs_CIDRs_smooth$harmonics[4]$coefs +
    for_score_5_prophet * fpcs_CIDRs_smooth$harmonics[5]$coefs +
    for_score_6_prophet * fpcs_CIDRs_smooth$harmonics[6]$coefs
}


################################ Multivariate models ###########################
# VAR (Vector Autorregresive)
library(ftsa)

# Define the objects that will contain the 47 predicted curves
n_for <- 47
for_CIDRs_var <- CIDRs_smooth$fd[(n - n_for):n]

for (i in 1:n_for) {
  
  print(paste("Predicción para el día:", i))
  
  fpcs_CIDRs_smooth <- pca.fd(CIDRs_smooth$fd[1:(n - n_for + i - 1)], nharm = 6, harmfdPar = fdPar(CIDRs_smooth$fd[1:(n - n_for + i - 1)]))
  CIDRs_scores <- fpcs_CIDRs_smooth$scores
  
  # Verify NA in FPCs scores
  if (any(is.na(CIDRs_scores))) {
    stop(paste("Valores NA encontrados en las puntuaciones principales en la iteración", i))
  }
  
  # Verify column names valid
  if (is.null(colnames(CIDRs_scores))) {
    colnames(CIDRs_scores) <- paste0("Score_", 1:ncol(CIDRs_scores))
  }
  
  # Convert scores to functional series (fts)
  fts_scores <- fts(x = 1:nrow(CIDRs_scores), y = CIDRs_scores)
  
  # Perform prediction using farforecast function
  forecast_var <- tryCatch({
    farforecast(fts_scores, h = 1, var_type = "const", Dmax_value = 6, Pmax_value = 1, level = 95, PI = FALSE)
  }, error = function(e) {
    warning(paste("Error en la predicción VAR en la iteración", i, ":", e$message))
    NULL
  })
  
  # Verify if the prediction was made correctly
  if (!is.null(forecast_var) && !is.null(forecast_var$mean)) {
    predicted_scores <- forecast_var$mean
    
    for_CIDRs_var$coefs[, i] <- mean.fd(CIDRs_smooth$fd)$coefs +
      as.numeric(predicted_scores[1]) * fpcs_CIDRs_smooth$harmonics[1]$coefs +
      as.numeric(predicted_scores[2]) * fpcs_CIDRs_smooth$harmonics[2]$coefs +
      as.numeric(predicted_scores[3]) * fpcs_CIDRs_smooth$harmonics[3]$coefs +
      as.numeric(predicted_scores[4]) * fpcs_CIDRs_smooth$harmonics[4]$coefs +
      as.numeric(predicted_scores[5]) * fpcs_CIDRs_smooth$harmonics[5]$coefs +
      as.numeric(predicted_scores[6]) * fpcs_CIDRs_smooth$harmonics[6]$coefs
    
  } else {
    warning("Predicciones no válidas en la iteración ", i)
  }
}



############ Compute the MISFE and the MIAFE using the 47 forecasted curves with all the methods considered ##################################

MISFE <- MIAFE <- matrix(0,nrow=n_for,ncol=6)
colnames(MISFE) <- colnames(MIAFE) <- c("ARIMA","STL","HW-Additive","NNAR", "Prophet","VAR")

for (i in 1 : n_for){
  
  err_CIDRs_arima <- for_CIDRs_arima[i] - CIDRs_smooth$fd[n-n_for+i]
  MISFE[i,1] <- norm.fdata(fdata(err_CIDRs_arima,argvals=err_CIDRs_arima$argvals,rangeval=err_CIDRs_arima$rangeval),
                           metric=metric.lp,lp=2)
  MIAFE[i,1] <- norm.fdata(fdata(err_CIDRs_arima,argvals=err_CIDRs_arima$argvals,rangeval=err_CIDRs_arima$rangeval),
                           metric=metric.lp,lp=1)
  
  err_CIDRs_stl <- for_CIDRs_stl[i] - CIDRs_smooth$fd[n-n_for+i]
  MISFE[i,2] <- norm.fdata(fdata(err_CIDRs_stl,argvals=err_CIDRs_stl$argvals,rangeval=err_CIDRs_stl$rangeval),
                           metric=metric.lp,lp=2)
  MIAFE[i,2] <- norm.fdata(fdata(err_CIDRs_stl,argvals=err_CIDRs_stl$argvals,rangeval=err_CIDRs_stl$rangeval),
                           metric=metric.lp,lp=1)
  
  err_CIDRs_hw_ad <- for_CIDRs_hw_ad[i] - CIDRs_smooth$fd[n-n_for+i]
  MISFE[i,3] <- norm.fdata(fdata(err_CIDRs_hw_ad,argvals=err_CIDRs_hw_ad$argvals,rangeval=err_CIDRs_hw_ad$rangeval),
                           metric=metric.lp,lp=2)
  MIAFE[i,3] <- norm.fdata(fdata(err_CIDRs_hw_ad,argvals=err_CIDRs_hw_ad$argvals,rangeval=err_CIDRs_hw_ad$rangeval),
                           metric=metric.lp,lp=1)
  
  err_CIDRs_ann <- for_CIDRs_ann[i] - CIDRs_smooth$fd[n-n_for+i]
  MISFE[i,4] <- norm.fdata(fdata(err_CIDRs_ann,argvals=err_CIDRs_ann$argvals,rangeval=err_CIDRs_ann$rangeval),
                           metric=metric.lp,lp=2)
  MIAFE[i,4] <- norm.fdata(fdata(err_CIDRs_ann,argvals=err_CIDRs_ann$argvals,rangeval=err_CIDRs_ann$rangeval),
                           metric=metric.lp,lp=1)
  
  err_CIDRs_prophet <- for_CIDRs_prophet[i] - CIDRs_smooth$fd[n-n_for+i]
  MISFE[i,5] <- norm.fdata(fdata(err_CIDRs_prophet,argvals=err_CIDRs_prophet$argvals,rangeval=err_CIDRs_prophet$rangeval),
                           metric=metric.lp,lp=2)
  MIAFE[i,5] <- norm.fdata(fdata(err_CIDRs_prophet,argvals=err_CIDRs_prophet$argvals,rangeval=err_CIDRs_prophet$rangeval),
                           metric=metric.lp,lp=1)
  err_CIDRs_var <- for_CIDRs_var[i] - CIDRs_smooth$fd[n-n_for+i]
  MISFE[i,6] <- norm.fdata(fdata(err_CIDRs_var,argvals=err_CIDRs_var$argvals,rangeval=err_CIDRs_var$rangeval),
                           metric=metric.lp,lp=2)
  MIAFE[i,6] <- norm.fdata(fdata(err_CIDRs_var,argvals=err_CIDRs_var$argvals,rangeval=err_CIDRs_var$rangeval),
                           metric=metric.lp,lp=1)
}  

# Plot comparative (Figure 5.4)
par(mfrow=c(1,2))
boxplot(MISFE,col=color_1,main="MISFE")
boxplot(MIAFE,col=color_1,main="MIAFE")

# Table 5.2
apply(MISFE, 2, median)

######## Compare the forecasts for January 9th, 2023 and 16th February 2023 #######################
# Figure 5.3
par(mfrow=c(1,2), mar=c(4.5, 5, 4, 2))  
# January 9th, 2023
plot(seq(0, 23, length.out=1000), eval.fd(seq(0,23,length.out=1000), CIDRs_smooth$fd[n-38]), 
     col="black", lwd=3, ylim=c(-100,400), type="l",
     main="Forecasts for January 9th, 2023", xlab="Hour of the day", ylab="CIDR", 
     cex.main=1.2, cex.lab=1.1, xaxt="n")  # Oculta el eje x predeterminado
axis(1, at=seq(0, 24, by=4), labels=c(0, 4, 8, 12, 16, 20, 24))

lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_arima[n_for-38]), col="blue", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_stl[n_for-38]), col="red", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_hw_ad[n_for-38]), col="green", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_ann[n_for-38]), col="violet", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_prophet[n_for-38]), col="orange", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_var[n_for-38]), col="yellow", lwd=3)

legend("topleft", legend=c("Original CIDR", "ARIMA", "STL", "HW-Additive", "ANN", "Prophet", "VAR"), 
       col=c("black", "blue", "red", "green", "violet", "orange", "yellow"), lwd=3, 
       cex=0.5, text.width = 2, pt.lwd = 0.5)  # Ajuste de tamaños de leyenda y símbolos

# 16th February 2023
plot(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), CIDRs_smooth$fd[n]), col="black", lwd=3, ylim=c(-100,400), type="l",
     main="Forecasts for February 16th, 2023", xlab="Hour of the day", ylab="CIDR", 
     cex.main=1.2, cex.lab=1.1, xaxt="n")  # Tamaño del título y etiquetas
axis(1, at=seq(0, 24, by=4), labels=c(0, 4, 8, 12, 16, 20, 24))

lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_arima[n_for]), col="blue", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_stl[n_for]), col="red", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_hw_ad[n_for]), col="green", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_ann[n_for]), col="violet", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_prophet[n_for]), col="orange", lwd=3)
lines(seq(0, 23, length.out=1000),eval.fd(seq(0,23,length.out=1000), for_CIDRs_var[n_for]), col="yellow", lwd=3)

legend("topleft", legend=c("Original CIDR", "ARIMA", "STL", "HW-Additive", "ANN", "Prophet", "VAR"), 
       col=c("black", "blue", "red", "green", "violet", "orange", "yellow"), lwd=3, 
       cex=0.7, text.width = 3.5, pt.lwd = 2)  # Ajuste de tamaños de leyenda y símbolos


########################### Prediction bands Prophet 16th February 2023 ##############################
# Figure 5.6
n_total <- 2238
alpha <- 0.05  # Significance level for prediction band
t_seq <- seq(0, 23, length.out=1000)

# Prophet prediction for the last day
theta_pred <- eval.fd(t_seq, for_CIDRs_prophet[n_for])  

# Obtain the functional residuals
theta_obs <- eval.fd(t_seq, CIDRs_smooth$fd[n_for])  # Real observations
residuals <- theta_obs - theta_pred  

# Functional standard deviation of residuals
residual_mean <- mean(residuals)
sigma_eps <- sqrt(1 / (n_total - 1) * sum((residuals - residual_mean)^2))

# q_L y q_U
q_L <- qnorm(alpha / 2)  
q_U <- qnorm(1 - alpha / 2) 

# Lower and upper bands
lower_band <- pred_next + q_L * sigma_eps  
upper_band <- pred_next + q_U * sigma_eps  


plot(t_seq, eval.fd(t_seq, CIDRs_smooth$fd[n]), type="l", col="black", lwd=3, ylim=c(min(lower_band), max(upper_band)),
     main="Forecasts for February 16th, 2023", xlab="Hour of the day", ylab="CIDR", 
     cex.main=1.2, cex.lab=1.1, xaxt="n") 
axis(1, at=seq(0, 24, by=4), labels=c(0, 4, 8, 12, 16, 20, 24))

lines(t_seq, pred_next, col="orange", lwd=3) # Prophet Forecast
lines(t_seq, lower_band, col="gray", lwd=2, lty=2)  # Lower band
lines(t_seq, upper_band, col="gray", lwd=2, lty=2)  # Upper band
polygon(c(t_seq, rev(t_seq)), c(upper_band, rev(lower_band)), 
        col=rgb(0.9, 0.9, 0.9, 0.5), border=NA)  # Band shading

legend("topleft", legend=c("Original CIDR", "Prophet Forecast", "95% Prediction Band"), 
       col=c("black", "orange", "gray"), lwd=3, 
       cex=0.7, text.width = 3.5, pt.lwd = 2)




