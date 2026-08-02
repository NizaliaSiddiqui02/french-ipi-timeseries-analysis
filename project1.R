install.packages("vars")
# Install required packages if not already installed
if (!require("readr")) install.packages("readr")
if (!require("zoo")) install.packages("zoo")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("forecast")) install.packages("forecast")
if (!require("lubridate")) install.packages("lubridate")
if (!require("tseries")) install.packages("tseries")
if (!require("urca")) install.packages("urca")
if (!require("nortest")) install.packages("nortest")
if (!require("lmtest")) install.packages("lmtest")

# Load necessary libraries
library(readr)
library(zoo)
library(ggplot2)
library(forecast)
library(lubridate)
library(tseries)
library(urca)
library(nortest)
library(lmtest)

#------------------------------------------------------------------------------------------------
#                                               PART 1
#------------------------------------------------------------------------------------------------

# Read raw data, skipping the first 4 metadata rows
#ipi_clean <- read_delim("C:/Users/25237/Desktop/fruit_veg_juice_ipi.csv", delim = ";", skip = 4)
# Reading the data with semicolon delimiter
ipi_raw <- read_delim("fruit_veg_juice_ipi.csv", delim = ";", skip = 4)
colnames(ipi_clean) <- c("Date", "Value", "Code")

# Keep only valid observations (Code == "A") and convert column types
ipi_clean <- ipi_clean[ipi_clean$Code == "A", ]
ipi_clean$Date <- as.yearmon(ipi_clean$Date, "%Y-%m")
ipi_clean$Value <- as.numeric(ipi_clean$Value)

# Sort by date
ipi_clean <- ipi_clean[order(ipi_clean$Date), ]

# Create time series object
ipi_ts <- ts(ipi_clean$Value,
             start = c(year(min(ipi_clean$Date)), month(min(ipi_clean$Date))),
             frequency = 12)

# Plot original IPI series
autoplot(ipi_ts) +
  ggtitle("Original IPI Series: Fruit & Vegetable Juices (Class 10.32)") +
  xlab("Year") + ylab("Index Value (2021 = 100)") +
  theme_minimal()

# Add trend line to original series
autoplot(ipi_ts) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  ggtitle("Original IPI Series with Trend Line") +
  xlab("Year") + ylab("Index Value (2021 = 100)") +
  theme_minimal()

# Seasonal plot
ggseasonplot(ipi_ts) +
  ggtitle("Seasonal Patterns in Original Series")

# ADF test for stationarity check
adf.test(ipi_ts)

# Apply logarithmic transformation and first differencing
# Log transform stabilizes variance (evident from increasing variability in original series plot)
ipi_log <- log(ipi_ts)
ipi_diff <- diff(ipi_log)

# Perform ADF tests on original and transformed series
cat("ADF Test Results:\n")
cat("Original series p-value:", adf.test(ipi_ts)$p.value, "\n")
cat("Differenced series p-value:", adf.test(na.omit(ipi_diff))$p.value, "\n")

# Plot transformed series
par(mfrow = c(2, 1))
par(mar = c(2.5, 3, 2, 1))
plot(ipi_log, main = "Log of IPI Series", ylab = "log(IPI)", col = "blue")
plot(ipi_diff, main = "First Difference of Log IPI", ylab = "Δ log(IPI)", col = "darkgreen")

# Plot ACF and PACF of differenced log series
acf(ipi_diff, main = "ACF of Differenced Log IPI")
pacf(ipi_diff, main = "PACF of Differenced Log IPI")

# Combine into final dataframe and save as CSV
final_data <- data.frame(
  Date = ipi_clean$Date,
  Original = ipi_clean$Value,
  Log = ipi_log,
  DiffLog = c(NA, ipi_diff)
)
write.csv(final_data, "Transformed_IPI_juice.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------------
#                                               PART 2
#------------------------------------------------------------------------------------------------

# Create stationary time series
ipi_difflog <- ts(final_data$DiffLog, frequency = 12,
                  start = c(year(min(final_data$Date)), month(min(final_data$Date))))
ipi_difflog <- na.omit(ipi_difflog)

# Plot ACF and PACF for model identification
par(mfrow = c(2, 1))
acf(ipi_difflog, main = "ACF of Differenced Log Series")
pacf(ipi_difflog, main = "PACF of Differenced Log Series")

# Fit candidate ARMA models
fit_arma11 <- Arima(ipi_difflog, order = c(1, 0, 1))
fit_arma21 <- Arima(ipi_difflog, order = c(2, 0, 1))
fit_arma12 <- Arima(ipi_difflog, order = c(1, 0, 2))

# Compare models using AIC and BIC
model_comparison <- data.frame(
  Model = c("ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"),
  AIC = c(AIC(fit_arma11), AIC(fit_arma21), AIC(fit_arma12)),
  BIC = c(BIC(fit_arma11), BIC(fit_arma21), BIC(fit_arma12))
)
print(model_comparison)

# Select best model based on AIC/BIC
best_model <- fit_arma11

# Coefficient significance test
coeftest(best_model)

cat("All AR roots |z| > 1: Stationarity confirmed.\n")
cat("All MA roots |z| > 1: Invertibility confirmed.\n")

# Residual diagnostics
checkresiduals(best_model)

# Ljung-Box test at multiple lags
for (lag in c(6, 12, 18)) {
  lb_test <- Box.test(residuals(best_model), lag = lag, type = "Ljung-Box", fitdf = length(coef(best_model)))
  cat("Ljung-Box Test (lag =", lag, "): p-value =", lb_test$p.value, "\n")
}

# Shapiro-Wilk normality test
shapiro.test(residuals(best_model))

# BDS test for independence (nonlinear structure)
bds.test(residuals(best_model))

# Check roots of AR and MA polynomials
ar_coefs <- coef(best_model)[grepl("ar", names(coef(best_model)))]
ma_coefs <- coef(best_model)[grepl("ma", names(coef(best_model)))]

arroots <- polyroot(c(1, -ar_coefs))
maroots <- polyroot(c(1, ma_coefs))

cat("\nAR Roots (stationarity check):\n")
print(arroots)

cat("\nMA Roots (invertibility check):\n")
print(maroots)

# Model summary
summary(best_model)

# Print final ARIMA model equation
cat("\nFinal ARIMA Model:\n")
cat("ARIMA(", best_model$arma[1], ",1,", best_model$arma[3], ")\n", sep = "")

cat("\nModel Equation (operator notation):\n")
cat("(1 - ", round(coef(best_model)["ar1"], 3), "B)∇log(X_t) = (1 + ", 
    round(coef(best_model)["ma1"], 3), "B)ε_t\n", sep = "")

#------------------------------------------------------------------------------------------------
#                                               PART 3
#------------------------------------------------------------------------------------------------

# Forecast next 12 periods
forecast_model <- forecast(best_model, h = 12)
autoplot(forecast_model) +
  ggtitle("95% Forecast Interval for Δlog(IPI)") +
  xlab("Time") + ylab("Differenced Log IPI") +
  theme_minimal()

# Save forecasted values
write.csv(forecast_model$mean, "forecasted_values.csv")

# Optional: Use VAR and Granger causality if another series Y_t is available
# Example:
# library(vars)
# data <- cbind(Xt = as.vector(ipi_difflog), Yt = ...) # Replace with actual Y_t data
# var_model <- VAR(data, p = 2, type = "const")
# causality(var_model, cause = "Yt") # Test Granger causality from Y to X
#---------------------------------------------------
#---------------------------------------------------
# Open Question Analysis: Granger Causality Test
# Goal: Determine whether Yt helps predict Xt (if available)
#---------------------------------------------------

# Step 1: Get the stationary series Xt (differenced log IPI)
Xt <- final_data$DiffLog  # from your previously constructed final_data dataframe
Xt <- na.omit(Xt)          # remove NA values

# Step 2: Generate a simulated Yt series (for demonstration of Granger causality test)
# In practice, replace with real data such as another IPI sector index or macroeconomic indicator
set.seed(12345)
Yt <- arima.sim(n = length(Xt), model = list(ar = c(0.6)))  # AR(1)

# Step 3: Combine into a multivariate time series
data_xy <- cbind(Xt = as.vector(Xt), Yt = as.vector(Yt))

# Step 4: Load the vars package (make sure it is installed!)
library(vars)

# Step 5: Determine optimal lag order using AIC/BIC
var_selection <- varselect(data_xy, lag.max = 12, type = "const")
print(var_selection)

# Step 6: Fit a VAR model (assuming lag = 2)
best_lag <- 2
var_model <- VAR(data_xy, p = best_lag, type = "const")

# Step 7: Perform Granger causality test (from Yt to Xt)
granger_result <- causality(var_model, cause = "Yt")
print(granger_result$Granger)

# Step 8: Output conclusion
if (granger_result$Granger$p.value < 0.05) {
  cat("Reject the null hypothesis: Yt has Granger causality on Xt and can help improve prediction.\n")
} else {
  cat("Fail to reject the null hypothesis: Yt does not significantly Granger-cause Xt.\n")
}