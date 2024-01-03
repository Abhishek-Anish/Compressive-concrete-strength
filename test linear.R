library(glmnet) 
library(ggcorrplot) 
library(tidyverse)
library(dplyr)
library(GGally)
library(leaps)
library(gridExtra)
library(coefplot)
library(car)

df<-read.csv("C:/Users/rocka/OneDrive/Documents/output_file.csv")

na_count <- colSums(is.na(df))
print(na_count)

#Correlation Plots

varsnum <-select(df, where(is.numeric)) |>
  relocate(CCS)
#+ fig.width=10.25, fig.height=10
ggcorrplot(corr=cor(varsnum), type="lower",
           ggtheme=ggplot2::theme_gray, lab=TRUE,
           title="Correlations for continuous variables")

gg <- ggpairs(varsnum,
              switch="y", upper=list(continuous="points"),
              lower=list(continuous="cor"))
#+ fig.width=8.25, fig.height=8
gg + labs(title="Scatterplots (and correlations) for continuous variables")

source("C:/Users/rocka/OneDrive/Documents/rss_regress_funcs_v5.R")

df <- rename(df, response=CCS)

lmodel<-lm(response~ . ,df)
summary(lmodel)

vif_results <- vif(lmodel)

print(vif_results) 


par(mar=c(4, 4, 2, 1))  
par(mfrow=c(2,2))
plot(lmodel)


normalQQ(lmodel,title2=". . .Using normalQQ()")

#-----------------------------------------------------------------------------------
# Linear regression using Best subset
#-----------------------------------------------------------------------------------



SEED<-1234
set.seed(SEED)

npredictors<-8

regfit.full <- regsubsets(response~., data=df) #nvmax not essential(8)
summary(regfit.full)
reg.summary <- summary(regfit.full)
names(reg.summary)

reg.summary$rsq

best.plot <- function(varName, varLabel, minmax=" ") {
  gg <- ggplot(data.frame(varName), aes(x=seq_along(varName), y=varName)) +
    geom_line() +
    labs(x="Number of variables", y=varLabel, title="Best subsets")
  if (minmax=="min") {
    gg <- gg + geom_point(aes(x=which.min(varName), y=min(varName)),
                          color="red") +
      geom_vline(aes(xintercept=which.min(varName)), linetype=
                   "dotted")
  }
  if (minmax=="max") {
    gg <- gg + geom_point(aes(x=which.max(varName), y=max(varName)),
                          color="red") +
      geom_vline(aes(xintercept=which.max(varName)), linetype=
                   "dotted")
  }
  return(gg)
}

results <- with(reg.summary, data.frame(rss,adjr2,cp,bic))
names(results) 

grid.arrange(best.plot(results$rss, "RSS"),
             best.plot(results$adjr2, "Adjusted RSq"
                       , "max"),
             best.plot(results$cp, "Cp", minmax="min"),
             best.plot(results$bic, "BIC", minmax="min"),
             ncol=2)

#Use plot.regsubsets() to plot accuracy measures vs selected predictors
dev.off()
plot(regfit.full,scale="r2")

plot(regfit.full,scale="adjr2")

plot(regfit.full,scale="Cp")

plot(regfit.full,scale="bic")



# Display coefficients for best 6- and 8-variables subsets
round(coef(regfit.full,6), 3)
round(coef(regfit.full,8), 3)

# Compute training-data MSEs for best 6 & 8 variable models
# for comparison with later cross-validated results

(MSE6 <- results$rss[6]/nrow(df))

(MSE8 <- results$rss[8]/nrow(df))

# Results so far:
#                     # predictors  training-set MSE
#   Best BIC & Cp model:      6         107.6148
#   Best Adj Rsq model:       8         107.1972



# Use a validation set to choose among models

# Use sampling with replacement to select an approx 50% training sample
set.seed(SEED)
train <- sample(c(TRUE,FALSE), nrow(df), replace=TRUE)
str(train)
test <- (!train)

# Compute test-set mean sum-of-squares for use in computing
# test-set R-squared
mean <- mean(df$response[test])
MSS <- mean((df$response[test]-mean)^2)

# Use the training sample to determine best subsets (for RSS) 
regfit.best <- regsubsets(response~., data=df[train,], nvmax=npredictors)

# Use the validation set to compute test-set MSEs.  There is no predict()
# function for regsubsets(,) so need code to compute predictions. This is
# done by creating an X matrix and then using matrix multiplication.
test.mat=model.matrix(response~.,data=df[test,]) # Create X matrix
val.errors=rep(NA,npredictors)                    # Reserves space to hold MSE's
for(i in 1:npredictors){
  coefi <- coef(regfit.best,id=i) # Coefficients for selected variables
  pred <- test.mat[,names(coefi)]%*%coefi # multiply X matrix by coefficients
  val.errors[i] <- mean((df$response[test]-pred)^2)  # mean squared error
}

val.errors
which.min(val.errors)  

round(coef(regfit.best, 5),3)  # Display coefficients for this model

# Repeat without having to view the console output to see # of variables
(pmin <- which.min(val.errors))
round(coef(regfit.best,pmin), 3)

# Display MSE, root-mean-squared error (RMSE), and test-set R-squared
(MSE <-val.errors[pmin])
(RMSE <- sqrt(MSE))
(Rsq <- 1 - MSE/MSS)

# Finally, use the full data set to calculate the coefficients for the 
# best 5 variable model
reg.best <- regsubsets(response~., data=df, nvmax=npredictors)
round(coef(reg.best,5),3)

# Results so far:
#                     # predictors  training-set MSE   Test MSE
#   Best BIC & Cp model:      6         107.6148
#   Best Adj Rsq model:       8         107.1972
#
#   Validation set:
#     Best test-set model     5                         111.8174 

# 4. ---------------------------------------------------------------------------
#  Use cross validation to choose among models.

# First, create a user-defined predict function for regsubsets().
# This function uses syntax not taught in class.
predict.regsubsets <- function(object, newdata, id, ...){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form,newdata)
  coefi <- coef(object, id=id)
  xvars <- names(coefi)
  mat[,xvars] %*% coefi
}

# Create random assignments for  k=10 folds
k <- 10
n <- nrow(df)
set.seed(SEED)
folds=sample(rep(1:k, length=n))

# Reserve space for a matrix to store the cross-validated MSEs
# for k test folds and all of the best models.  Initialize to NA. 
cv.errors=matrix(NA,k,npredictors,
                 dimnames=list(NULL, paste(1:npredictors)))

# For each test fold, use the k-1 other folds to determine best models and
# use the user-defined predict() function to compute each model's test-fold MSE.
for(j in 1:k){
  best.fit=regsubsets(response~.,data=df[folds!=j,],nvmax=npredictors)
  for(i in 1:npredictors){
    pred=predict(best.fit,df[folds==j,], id=i)
    cv.errors[j,i] <- mean( (df$response[folds==j] - pred)^2)
  }
}

# Compute the cross-validated MSE for each model.
mean.cv.errors = rep(0, npredictors)
for (i in 1:npredictors) {
  mean.cv.errors[i] <- mean(cv.errors[,i])
}
round(mean.cv.errors,0)

# Use the user-defined best.plot function to display results
best.plot(mean.cv.errors, "Cross-validated MSE", "min")

# Display cross-validated MSE and RMSE for best 6-variable model
(MSE <- mean.cv.errors[6])
(RMSE <- sqrt(MSE))

mean <- mean(df$response)
MSS <- mean((df$response-mean)^2)
(Rsq <- 1 - MSE/MSS)

# # Finally, use the full data set to calculate the coefficients for the 
# best 6 variable model
reg.best <- regsubsets(response~., data=df, nvmax=npredictors)
round(coef(reg.best,6),3)

# Results so far:
#                     # predictors  training-set MSE   Test MSE
#   Best BIC & Cp model:      6         107.6148
#   Best Adj Rsq model:       8         107.1972
#
#   Validation set:
#   Best test-set model       5                         111.8174 
#   K-fold CV best model      6                         109.1922


#--------------------------------------------------------------------------------
#Lasso regression
#--------------------------------------------------------------------------------
y<-df$response
x<-model.matrix(response~.,df)[,-1]

# Create training set 
set.seed(SEED)
train <- sample(1:nrow(x), size=nrow(x)/2)
test <- (-train)

# Recompute test-set mean sum-of-squares for use in computing
# test-set R-squared
mean <- mean(df$response[test])
MSS <- mean((df$response[test]-mean)^2)

# Use training set to perform lasso regression for 100 lambda values
# from 10^10 to 1/10^189.  alpha=1 for lasso(default)
grid <- 10^seq(10,-2,length=100)
lasso.mod <- glmnet(x[train,],y[train], alpha=1, lambda=grid)

# Cross validation of the training set determines the lambda minimizing MSE.
# Also, determines lambda for the 1-standard-error rule.
set.seed(SEED)
cv.out <- cv.glmnet(x[train,],y[train],alpha=1, nfold=10)
(bestlam <- cv.out$lambda.min)
(bestlam.1se <- cv.out$lambda.1se)

plot(cv.out)

# Plot coefficient paths using plot() function in glmnet package
par(mfrow=c(1,1))
plot(cv.out$glmnet.fit, xvar="lambda", label=TRUE) # function in glmnet package
abline(v=log(c(bestlam, bestlam.1se)), lty=2)

# Determine test-set MSE for lambda minimizing cross-validated training-set MSE
y.test <- y[test]
lasso.pred.min=predict(lasso.mod, s=bestlam, newx=x[test,])
(MSE.lasso.min <- mean((lasso.pred.min-y.test)^2))
(RMSE <- sqrt(MSE.lasso.min))
(Rsq <- 1 - MSE.lasso.min/MSS)

#Plot residuals vs fitted for best-lambda lasso
#+ fig.height=4
ggplot(data.frame(y.test,lasso.pred.min),
       aes(x=lasso.pred.min,
           y=lasso.pred.min-y.test)) +
  geom_point() +
  labs(x="fitted", y="residual",
       title="Residual vs fitted for best lambda lasso")

# Determine test-set MSE for lambda from 1-standard-error rule
lasso.pred.1se=predict(lasso.mod, s=bestlam.1se, newx=x[test,])
(MSE.lasso.1se <-mean((lasso.pred.1se-y.test)^2))
(RMSE <- sqrt(MSE.lasso.1se))
(Rsq <- 1 - MSE.lasso.1se/MSS)


# Finally, use the full data set to perform lasso regression and
# display non-zero coefficients for best lambda determined from
# cross validation of the training set
out <- glmnet(x,y, alpha=1, lambda=grid)
lasso.coef.min <- predict(out, type="coefficients",
                          s=bestlam)[1:(ncol(x)+1),]
round(lasso.coef.min[lasso.coef.min != 0], 3)

# Number of predictors with non-zero coefficients for best lambda
(nvars <- length(lasso.coef.min[lasso.coef.min!=0])-1)

# Use coefplot() in coefplot package to plot coefficients for best lambda
coefplot(out, lambda=bestlam, sort='magnitude',
         title="Nonzero coefficients for best lambda")

# Also, display non-zero coefficients for lambda from 1se rule
lasso.coef.1se <- predict(out,type="coefficients",s=bestlam.1se)[1:(ncol(x)+1),]
round(lasso.coef.1se[lasso.coef.1se!=0], 3)

# Number of predictors with non-zero coefficients for 1se rule
(nvars <- length(lasso.coef.1se[lasso.coef.1se!=0])-1)

# Use coefplot() in coefplot package to plot coefficients for 1se rule
coefplot(out, lambda=bestlam.1se, sort='magnitude',
         title="Nonzero coefficients for 1se rule")

# Results so far:
#                     # predictors  training-set MSE   Test MSE
#   Best BIC & Cp model:      6         107.6148
#   Best Adj Rsq model:       8         107.1972
#
#   Validation set:
#   Best test-set model       5                         111.8174 
#   K-fold CV best model      6                         109.1922
#
#   Best-lambda Lasso         8                         113.5209
#   1se-lambda Lasso          7                         114.7777
