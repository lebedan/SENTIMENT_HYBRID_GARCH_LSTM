
garch_x<-function(file,garch_train_size){
  
suppressWarnings({

library(xts)

library(rugarch)

# Assuming 'data' is your pandas DataFrame, convert it to an R data frame

data <- read.csv(file)
data$Date<-as.Date(data$Date, format="%Y-%m-%d")

data[, -1] <- sapply(data[, -1], as.numeric)

data[is.na(data)] <- 0
endog_var=c("Log.Returns")
exog_vars=c("News.Positive.Sentiment", "News.Negative.Sentiment", "Twitter.Positive.Sentiment","Twitter.Negative.Sentiment")


# Define the range for number of lags to test
max_lags <- 21
bic_values <- numeric(max_lags)

# Function to create lagged variables
create_lagged_data <- function(data, exog_vars, num_lags) {
  for (i in 1:ncol(exog_vars)) {
    for (lag in 1:num_lags) {
      col_name <- paste0("lag", lag, "_", names(exog_vars)[i])  # Create lagged column name
      data[[col_name]] <- c(rep(NA, lag), exog_vars[[i]][1:(nrow(exog_vars)-lag)])  # Create lagged column values
    }
  }
  data[is.na(data)] <- 0
  return(data)
}

# Loop through different values of num_lags
for (num_lags in 1:max_lags) {
  # Create lagged variables
  lagged_data <- create_lagged_data(data, data[exog_vars], num_lags)
  lagged_exog_vars <- grep("lag", names(lagged_data), value = TRUE)
  
  # Scale the exogenous variables
  exog <- scale(lagged_data[lagged_exog_vars])
  
  # Define the endogenous variable
  endog <- xts(lagged_data[endog_var], order.by = lagged_data$Date)
  
  # Fit GARCH model
  garchx <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                             external.regressors = exog),
                       mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                       distribution.model = "norm")
  
  garch_fit <- tryCatch(ugarchfit(spec = garchx, data = endog), error = function(e) NULL)
  
  # Calculate BIC
  if (is.null(garch_fit)) {
    bic_values[num_lags] <- Inf
  } else {
    tryCatch({bic_values[num_lags] <- infocriteria(garch_fit)[2]
    }, error = function(e) {
      bic_values[num_lags] <- Inf
    })
  }
}

# Select the number of lags with the lowest BIC
optimal_lags <- which.min(bic_values)
if(optimal_lags==Inf) {
  optimal_lags <-1
} 

print(paste("Optimal number of lags:", optimal_lags))

# Create lagged variables with the optimal number of lags
final_lagged_data <- create_lagged_data(data, data[exog_vars], optimal_lags)
final_lagged_exog_vars <- grep("lag", names(final_lagged_data), value = TRUE)
final_exog <- scale(final_lagged_data[final_lagged_exog_vars])

# Define the endogenous variable
final_endog <- xts(final_lagged_data[endog_var], order.by = final_lagged_data$Date)

# Fit the final model
garchx_final <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                                 external.regressors = final_exog),
                           mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                           distribution.model = "norm")

#garchx.res <- ugarchfit(garchx, endog, xreg =exog, solver = 'hybrid')
expanding_forecast <- ugarchroll(garchx_final, data = endog, xreg = final_exog,
                                 n.ahead = 1, forecast.length = 1,
                                 refit.every = 1, n.start = garch_train_size,solver = 'hybrid')
sigma<-expanding_forecast@forecast$density[,"Sigma"]
write.csv(expanding_forecast@forecast$density, file = "~/Balázs/TDK-Data-Analysis-main/Only GARCH_X/Python_Input.csv", row.names = TRUE)
})
}
#garch_x("~/Balázs/TDK-Data-Analysis-main/Only GARCH_X/R_data.csv",256)
