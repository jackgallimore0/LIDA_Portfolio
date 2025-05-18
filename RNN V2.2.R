# Version 2.2 
# This version will allow variation of paremters
# takes data from the new splined data
# All output data frames will remain the same

test2 <- RNN2()

RNN2 <- function(range = c(3:5), x = 1, lay = 10, drop = 0.1, batch = 64, epoc = 50) {
  
  cat("#################\n",
      "Summary of RNN2 function parameters:\n",
      "Range: ", paste(range, collapse = ", "), "\n",
      "X: ", x, "\n",
      "Layers: ", lay, "\n",
      "Dropout rate: ", drop, "\n",
      "Batch size: ", batch, "\n",
      "Epochs: ", epoc, "\n",
      "#################\n")
  
  set.seed(123)  # Set seed for R
  tensorflow::set_random_seed(123)  # Set seed for TensorFlow
  
  # Load mortality data for age groups from Excel file
  year_5 <- read_excel("Splined_Mortality_Zoe.xlsx")
  ages <- colnames(year_5)
  
  # Define the age groups to use for training
  #range <- c(14:16)
  
  # Select the specific age group from the defined range for prediction
  #x <- 1
  
  # Filter the data to include only the specified age groups
  year_5 <- year_5[, range]
  
  # Create a logical vector for training (70% TRUE, 30% FALSE)
  istrain0 <- c(rep(TRUE, round(nrow(year_5) * 0.7)), rep(FALSE, round(nrow(year_5) * 0.3)))
  
  # Add the training vector to the data and convert to data frame
  year_5$train <- istrain0
  year_5 <- data.frame(year_5)
  test_for <- c("1", "2", "3")
  colnames(year_5) <- c(test_for, "train")
  
  # Create a matrix of the mortality data for the specified age groups
  xdata0 <- data.matrix(year_5[, test_for])
  
  # Scale the data
  xdata0 <- scale(xdata0)
  
  # Function to create lagged matrix
  lagm <- function(x, k = 1) {
    n <- nrow(x)
    pad <- matrix(NA, k, ncol(x))
    rbind(pad, x[1:(n - k), ])
  }
  
# Create a data frame with lagged variables
arframe0 <- data.frame(
  x = xdata0[, test_for[x]],
  L1 = lagm(xdata0, 1),
  L2 = lagm(xdata0, 2),
  L3 = lagm(xdata0, 3),
  L4 = lagm(xdata0, 4),
  L5 = lagm(xdata0, 5)
)

# Remove the first 5 rows with NA values due to lagging
arframe0 <- arframe0[-(1:5), ]
istrain0 <- istrain0[-(1:5)]

# Fit a linear model to the training data
arfit0 <- lm(x ~ ., data = arframe0)

# Predict using the linear model on the test data
arpred0 <- predict(arfit0, arframe0[!istrain0, ])

# Calculate variance of the test data
v00 <- var(arframe0[!istrain0, "x"])

# Number of rows in the data frame
n0 <- nrow(arframe0)

# Convert data frame to matrix excluding the response variable
xrnn0 <- data.matrix(arframe0[, -1])

# Format the data into a 3D array for RNN input
xrnn0 <- array(xrnn0, c(n0, 3, 5))

# Reverse the order of the 3rd dimension (lags)
xrnn0 <- xrnn0[, , 5:1]

# Reformat the array with shifting windows
xrnn0 <- aperm(xrnn0, c(1, 3, 2))

# Define the RNN model
model0 <- keras_model_sequential() %>%
  layer_simple_rnn(units = lay,
                   input_shape = list(5, 3),
                   dropout = drop, recurrent_dropout = 0.1) %>%
  layer_dense(units = 1)

# Compile the model
model0 %>% compile(
  optimizer = optimizer_rmsprop(),
  loss = "mse"
)

# Train the model
history <- model0 %>% fit(
  xrnn0[istrain0, , ], arframe0[istrain0, "x"],
  batch_size = batch, epochs = epoc,
  validation_data = list(xrnn0[!istrain0, , ], arframe0[!istrain0, "x"])
)

# Make predictions with the trained model
kpred0 <- predict(model0, xrnn0[!istrain0, , ])
actual_mort <- year_5[!year_5$train, test_for[x]]

Training_Loss = history$metrics$loss
Validation_Loss = history$metrics$val_loss
Epoch = 1:length(history$metrics$loss)

#Stores the data required for graph
Loss_data <- data.frame(
  Epoch = Epoch,
  Training_Loss = Training_Loss,
  Validation_Loss = Validation_Loss
)

#Data frame for the actual and the prediction
Comparison_table <- data.frame(
  Month = Ghana_Dates$Month[!istrain0],
  Year = Ghana_Dates$Year[!istrain0],
  Actual_mortality = actual_mort,
  Scaled_prediction = c(kpred0)
)

Comparison_table$Unscaled_prediction <- (kpred0 * sd(actual_mort)) + mean(actual_mort)

# Squared error MSE and RMSE
Comparison_table$Squared_Error <- (Comparison_table$Actual_mortality - Comparison_table$Unscaled_prediction)^2

mean_actual <- mean(Comparison_table$Actual_mortality)
MSE <- mean(Comparison_table$Squared_Error)
RMSE <- MSE / (mean_actual^2)

#Summary of loss and loss val trend and total mse
model1 <- lm(Training_Loss ~ Epoch)
model2 <- lm(Validation_Loss ~ Epoch)

slope_TL <- as.numeric(coef(model1)[2])
slope_VL <- as.numeric(coef(model2)[2])

#List of MSE and slope values
results <- list(
  MSE = MSE,
  RMSE = RMSE,
  slope_TL = slope_TL,
  slope_VL = slope_VL
)

#This result can be returned in a function
#can get data out like result$loss_data$Training_Loss
output <- list(
  Loss_data = Loss_data,
  Comparison_table = Comparison_table,
  results = results
)

return(output)
}
