#load necessary libraries
library(KFAS)         # Kalman Filtering & State-Space Models
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tseries)      # ADF Test

# Load data (replace path with actual file)
data <- read.csv("C:/Users/sanga/Desktop/ito/GSPC.csv", stringsAsFactors = FALSE)
data$Date <- as.Date(data$Date, format = "%m-%d-%Y")
data <- data %>% arrange(Date)
head(data)

# Compute Log Returns (Initial values set to 0)
data <- data %>% mutate(Return = c(0, diff(log(data$GSPC.Adjusted))))
head(data)

# Function to calculate first-order autocorrelation using a moving window
moving_window_autocorr <- function(returns, window_size) {
  n <- length(returns)
  autocorrs <- numeric(n - window_size + 1)
  
  for (i in 1:(n - window_size + 1)) {
    window <- returns[i:(i + window_size - 1)]
    autocorrs[i] <- cor(window[-length(window)], window[-1], use = "complete.obs")
  }
  
  return(autocorrs)
}

# Apply moving window autocorrelation function
window_size <- 100  # Using a 100-day moving window
autocorrs <- moving_window_autocorr(data$Return, window_size)
head(data)

# Add the autocorrelations to the data frame (aligned with the end of each window)
data$Autocorrelation <- c(rep(NA, window_size - 1), autocorrs)
head(data)
tail(data)


# Plot the moving window autocorrelation
ggplot(data, aes(x = Date)) +
  geom_line(aes(y = Autocorrelation, color = "Moving Window Autocorrelation")) +
  labs(title = "Moving Window Autocorrelation of S&P 500 Returns",
       x = "Date",
       y = "First-Order Autocorrelation",
       color = "Metric") +
  theme_minimal()

# Load necessary package
library(urca)

# Initialize vector to store BIC values
bic_values <- numeric(6)

# Loop through lags 0 to 5
for (i in 0:5) {
  test <- ur.df(data$Return, type = "trend", lags = i)
  
  # Extract residuals from the ADF regression
  residuals <- test@res
  n <- length(residuals)
  sigma2 <- sum(residuals^2) / n
  
  # Calculate BIC manually
  k <- i + 2  # i lags + intercept + trend
  bic <- n * log(sigma2) + k * log(n)
  
  bic_values[i + 1] <- bic
  cat("Lag:", i, "- BIC:", round(bic, 3), "\n")
}

best_lag <- which.min(bic_values) - 1
cat("Best lag based on BIC:", best_lag, "\n")

# Perform ADF test with constant and trend
# Use the ur.df function from 'urca' package
# Set lags = 0 for lag order (as per SBIC selection)
adf_result <- ur.df(data$Return, type = "trend", lags = 0)

# Show test results
summary(adf_result)

##  or #####

# Augmented Dickey-Fuller Test (with trend & constant, SBIC criterion)
adf_result <- adf.test(data$Return, k = 0)  # Lag 0 chosen based on SBIC
print(adf_result)

### or ######

# Define State-Space Model for Time-Varying AR(1)
n <- nrow(data)

# Observation Equation: Return_t = AR_t * Return_{t-1} + Noise
Zt <- array(0, dim = c(1, 1, n))  
for (i in 2:n) {
  Zt[,,i] <- data$Return[i - 1]
}

Ht <- matrix(NA)  # Observation variance (estimated)

# Transition Equation: AR_t = AR_{t-1} + Noise (Random Walk)
Tt <- array(1, dim = c(1, 1, n))  # Unit root assumption
Qt <- matrix(NA)  # Process variance (estimated)

# Initial state: Market starts fully efficient (AR = 0)
a1 <- matrix(0)  # Initial state at zero
P1 <- matrix(0)  # Initial variance at zero (as per the paper)
P1inf <- matrix(0)  # Ensuring proper estimation

# Build State-Space Model
model <- SSModel(
  data$Return ~ SSMcustom(
    Z = Zt,  
    T = Tt,  
    Q = Qt,  
    a1 = a1,  
    P1 = P1,  
    P1inf = P1inf
  ), 
  H = Ht
)

# Estimate Model with Kalman Filter
fit <- fitSSM(model, inits = c(0.1, 0.1), method = "BFGS")
smoothed_states <- KFS(fit$model, smoothing = "state")

# Extract the time-varying AR(1) coefficients
ar_coefficients <- as.numeric(smoothed_states$alphahat)

# Ensure ar_coefficients aligns with data
if (length(ar_coefficients) > nrow(data)) {
  ar_coefficients <- tail(ar_coefficients, nrow(data))
} else if (length(ar_coefficients) < nrow(data)) {
  ar_coefficients <- c(rep(0, nrow(data) - length(ar_coefficients)), ar_coefficients)
}

# Assign AR coefficients to data
data$AR_Coefficient <- ar_coefficients
head(data)
tail(data)

# Plot the Time-Varying AR(1) Coefficients (Relative Market Inefficiency)
ggplot(data, aes(x = Date, y = AR_Coefficient)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "SP500", 
       x = "Time", 
       y = "AR Coefficient") 
theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Plot Time-Varying AR(1) Coefficients (Market Inefficiency Indicator)
ggplot(data, aes(x = Date, y = AR_Coefficient)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Time-Varying AR(1) Coefficient of S&P 500",
       x = "Time", y = "AR Coefficient") +
  theme_minimal()

a1 <- matrix(0)  # Ensure AR coefficient starts at zero
P1 <- matrix(1e-6)  # Small initial variance to enforce stability
Tt <- array(1, dim = c(1, 1, n))  # AR coefficient follows a random walk
Qt <- matrix(0.0001)  # Small process variance to avoid overfitting

cat("Data rows:", nrow(data), "\n")
cat("AR Coefficient rows:", length(ar_coefficients), "\n")
ar_coefficients <- tail(ar_coefficients, nrow(data))
ar_coefficients <- c(rep(0, nrow(data) - length(ar_coefficients)), ar_coefficients)
data$AR_Coefficient <- ar_coefficients
ggplot(data, aes(x = Date, y = AR_Coefficient)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "SP500 - Time-Varying AR(1) Coefficient", 
       x = "Time", 
       y = "AR Coefficient") +
  theme_minimal()






