#---------------------------------------------------------------------------
#Linear Regression using Stepwise search for minimum AIC models
#---------------------------------------------------------------------------
library(MASS)
library(dplyr)

df<-read.csv("C:/Users/rocka/OneDrive/Documents/output_file.csv")

na_count <- colSums(is.na(df))
print(na_count)

df <- rename(df, response=CCS)

SEED<-1234

# Search for minimum AIC model using all the data
modeltry <- lm(response ~ .^2, data = df)
length(coef(modeltry))

step.model <- stepAIC(modeltry, trace=0)
length(coef(step.model))
round(coef(step.model), 3)
AIC(modeltry, step.model)

# Validation-set approach to estimate MSE
# Create indices for training and test data
set.seed(SEED)
nobs <- nrow(df)
train <- sample(1:nobs, nobs/2)
test <- (-train)

# Computing test-set mean-sum-of-squares for use in computing
# test-set R-squared
mean <- mean(df$response[test])
MSS <- mean((df$response[test]-mean)^2)

#Searching for minimum AIC model using the training data
step.model <- stepAIC(lm(response~.^2, data=df[train,]),
                      trace=0)

#Using the test data to estimate MSE, RMSE, and test-set R-squared
truth <- df[test,]$response
pred <- predict(step.model, newdata=df[test,])
(MSE <- mean((truth-pred)^2))
(RMSE <- sqrt(MSE))
(R_sq <- 1 - MSE/MSS)

#Doing a final build of the model using all of the data
names(step.model)
model.final <- lm(step.model$model, data=df)

(coef <- round(coef(model.final), 5))


# number of effects (other than the intercept) in the final model
(length(coef) -1)
