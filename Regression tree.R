#-------------------------------------------------------------------------------------
#Regression Tree
#-------------------------------------------------------------------------------------

library(rpart)
library(rpart.plot)
library(DMwR2)
library(dplyr)
library(tidyverse)

df<-read.csv("C:/Users/rocka/OneDrive/Documents/output_file.csv")


df <- rename(df, response=CCS)

SEED<-1234

set.seed(SEED)
rpart.modela<-rpart(response ~ .,data=df,method="anova", cp=0.001)

plotcp(rpart.modela)

rpart.modela.1se <- rt.prune(rpart.modela, se=1)

rpart.plot(rpart.modela.1se , extra=1, roundint=FALSE,
           digits=3, main="1se regression tree for concrete compressive strength (in MPa)",
           cex=0.5)

var(df$response)

printcp(rpart.modela.1se)

(MSE<-278.81*0.24033)
(rmse<-sqrt(MSE))

# Use a validation set to estimate test-set MSE
set.seed(7)
train = sample(1:nrow(df), nrow(df)/2)

set.seed(SEED)
rpart.concrete.train <- rpart(response~
                                ., data = df[train,],
                              method="anova", cp=0)
rpart.concrete.train <- rt.prune(rpart.concrete.train, se=1)

Yhat <- predict(rpart.concrete.train, newdata = df[-train,])
concrete.test <- df[-train,"response"]

(MSE <- mean((Yhat-concrete.test)^2))

(RMSE <- sqrt(MSE))

#repeating with different seed for creating training set
set.seed(5)
train <- sample(1:nrow(df), nrow(df)/2)

set.seed(SEED)
rpart.concrete.train <- rpart(response~
                                ., data = df[train,],
                              method="anova", cp=0)
rpart.concrete.train <- rt.prune(rpart.concrete.train, se=1)

Yhat <- predict(rpart.concrete.train, newdata = df[-train,])
concrete.test <- df[-train,"response"]

(MSE <- mean((Yhat-concrete.test)^2))

(RMSE <- sqrt(MSE))

ggplot(data.frame(Yhat, concrete.test), aes(x=Yhat ,y=concrete.test)) +
  geom_point() +
  geom_abline(slope=1,intercept=0) +
  labs(x="predicted CSS",
       y="test-set CSS",
       title="regression tree")

#  Results:                                 MSE           RMSE
#     Cross-validated MSE                  67.00641     8.185744
#     Validation set MSE (set.seed(7))     68.21776     8.259404
#     Validation set MSE (set.seed(5))     64.34217     8.021357
