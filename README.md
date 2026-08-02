# 📈 French Industrial Production Index (IPI) Time Series Analysis

An end-to-end time series analysis and forecasting project built using **R**, focusing on the French Industrial Production Index (IPI) for the **Preparation of Fruit and Vegetable Juices** subsector (INSEE / NAF Rev. 2, Class 10.32).

## 📄 Documentation & Code
- **Full Report**: [`docs/Report.pdf`](docs/Report.pdf)
- **Main R Script**: [`project1.R`](project1.R)
- **Raw Data**: [`data/fruit_veg_juice_ipi.csv`](data/fruit_veg_juice_ipi.csv)

---

##  Project Overview
The objective of this study is to model production dynamics in the French juice manufacturing industry using monthly index data ($2021 = 100$). The project covers:
1. **Data Preprocessing & Cleaning**: Filtering validated observations (`Code A`) and converting dates to `yearmon` format.
2. **Stationarity Transformations**: Applying logarithmic transformations and first-differencing ($\Delta \log(X_t)$) to address non-stationarity.
3. **Statistical Testing**: Conducting Augmented Dickey-Fuller (ADF) tests to verify stationarity.
4. **ARMA Model Identification & Selection**: Evaluating $ARIMA(p,0,q)$ candidate models using AIC and BIC.
5. **Diagnostics & Validation**: Coefficient significance tests (`coeftest`) and Ljung-Box test for residual autocorrelation.
6. **Forecasting**: Generating a 12-period forward forecast with 95% confidence intervals.
7. **Multivariate Extension**: Implementing a Vector Autoregression (VAR) framework to test Granger causality.



## 🛠️ Key Results

| Metric / Test | Result / Best Value | Conclusion |
| :--- | :--- | :--- |
| **Original ADF $p$-value** | $0.9763$ | Non-stationary |
| **Differenced Log ADF $p$-value** | $0.0132$ | **Stationary** ($p < 0.05$) |
| **Selected Model** | **$\text{ARIMA}(1,0,1)$** | Lowest AIC ($-1253.380$) & BIC ($-136.54$) |
| **$MA(1)$ Parameter** | $-0.5858$ ($p = 2.15 \times 10^{-12}$) | Highly statistically significant |



## 🚀 Getting Started

### Prerequisites
Ensure you have **R (>= 4.0)** installed along with the following packages:
```R
install.packages(c("readr", "zoo", "ggplot2", "forecast", "lubridate", "tseries", "urca", "nortest", "lmtest", "vars"))
