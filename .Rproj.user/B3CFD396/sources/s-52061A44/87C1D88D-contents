# The purpose of this script is to train a bunch a models and stack them


# file paths
scripts <- 'C:/Users/tommy/Google Drive/Coursework/1machine_learning/scripts/'
data <- "C:/Users/tommy/Google Drive/Coursework/1machine_learning/data/"
raw <- "C:/Users/tommy/Google Drive/Coursework/1machine_learning/data/raw/"

# load custom functions
source(paste0(scripts,'helper_functions.R'))

# load packages
LoadPackages(c("ggplot2", "caret", 'ranger', 'glmnet', 
               'kernlab', 'stringr', 'xgboost', 'plyr', 'fastICA', 
               'earth', 'Cubist', "dplyr", 'corrplot'))

# load data
submit_data <- readRDS(paste0(data,"submit_data.rds"))
train <- readRDS(paste0(data,"train.rds")) %>% select(-id, -set)
test <- readRDS(paste0(data,"test.rds"))

#  create list with indexes for each cv fold
#  this makes the models more comaparable bc each one will use the same folds
set.seed(888)
cv_folds <- createFolds(train$target, k = 5, list = TRUE, returnTrain = TRUE)

# create object to control model tuning
# save predictions to get out of sample predictions for stacking
ctrl <- trainControl(method = "cv", index = cv_folds, savePredictions = 'final',
                     verboseIter = T)

# RandomForest
set.seed(888)
rf_grid <- expand.grid(mtry = c(2,5,9,12), 
                       min.node.size = c( 60, 80, 85, 90, 100),
                       splitrule = c("variance", 'extratrees'))

rf <- train(target ~ ., data = train, tuneGrid = rf_grid, num.tree = 1000, 
            importance = 'permutation', metric = "MAE", method = "ranger", 
            trControl = ctrl)

rf$results %>% arrange(MAE)

saveRDS(rf, paste0(data, "rf.rds"))

# Elastice net regression
set.seed(888)
enet_grid <- expand.grid(alpha = seq(0,1, length = 15), 
                         lambda = c(seq(0.0001, 1, length =  300),
                                    seq(1,5, length = 20))) %>%
  sample_frac(.2) 

enet <- train(target ~ . , data = train, tuneGrid = enet_grid,
              preProcess = c("center", "scale"), metric = "MAE",
              method = "glmnet", trControl = ctrl)

enet$results %>% arrange(MAE)

saveRDS(enet, paste0(data, "enet.rds"))

# Support Vector Regression
set.seed(888)
svm_grid <- expand.grid(sigma = 2^seq(-10, -7, length = 12), 
                        C = seq(.8, 1.5, by = .1))

svm <- train(target ~ ., data = train, tuneGrid = svm_grid,
             metric = "MAE", method = "svmRadial", 
             preProcess = c("center", "scale"), trControl = ctrl)

svm$results %>% arrange(MAE) %>% head()

saveRDS(svm, paste0(data, "svm.rds"))

# k nearest neighbor
set.seed(888)
knn_grid <- expand.grid(k = seq(10,120, by = 5))

knn <- train(target ~ ., data = train, tuneGrid = knn_grid,
             metric = "MAE", method = "knn", 
             preProcess = c("center", "scale"), trControl = ctrl)

knn$results %>% arrange(MAE)

saveRDS(knn, paste0(data, "knn.rds"))

# xgboost
set.seed(888)
xgb_grid <- expand.grid(nrounds = c(50, 100, 150), 
                        max_depth = c(2, 5, 10, 15), 
                        eta = seq(.01,.3, by = .04), 
                        gamma = c(0, 1, 3), 
                        colsample_bytree = c(.5, .7, .9), 
                        min_child_weight = seq(5, 55, by = 10), 
                        subsample = c(.5, .7, .9)) %>%
  sample_frac(.06)

xgb <- train(target ~ ., data = train, tuneGrid = xgb_grid,
             metric = "MAE", method = "xgbTree",trControl = ctrl)

xgb$results %>% arrange(MAE) %>% head()

saveRDS(xgb, paste0(data, "xgb.rds"))


# Independent Component Regression
set.seed(888)
icr <- train(target ~ . , data = train, tuneLength = 20, # picks grid for u (got lazy)
             metric = "MAE", method = "icr", 
             preProcess = c("center", "scale"), trControl = ctrl)

icr$results %>% arrange(MAE)

# Projection Pursuit Regression
set.seed(888)
ppr <- train(target ~ . , data = train, tuneLength = 20,
             metric = "MAE", method = "ppr", 
             preProcess = c("center", "scale"), trControl = ctrl)

ppr$results %>% arrange(MAE)

# Bagged MARS
set.seed(888)
bagEarth <- train(target ~ . , data = train, tuneLength = 20, 
                  metric = "MAE", method = "bagEarth", 
                  preProcess = c("center", "scale"), trControl = ctrl)

bagEarth$results %>% arrange(MAE)

# cubist
set.seed(888)
cubist <- train(target ~ . , data = train, tuneLength = 20,
                metric = "MAE", method = "cubist", 
                preProcess = c("center", "scale"), trControl = ctrl)

cubist$results %>% arrange(MAE)

# null model (predicts the mean)
null <- train(target ~ . , data = train, 
              metric = "MAE", method = "null", trControl = ctrl)

null$results %>% arrange(MAE)

# Create data frame for out of sample predictions from the models
rf_oos <- rf$pred %>% arrange(rowIndex) %>% pull(pred)
enet_oos <- enet$pred %>% arrange(rowIndex) %>% pull(pred)
svm_oos <- svm$pred %>% arrange(rowIndex) %>% pull(pred)
knn_oos <- knn$pred %>% arrange(rowIndex) %>% pull(pred)
xgb_oos <- xgb$pred %>% arrange(rowIndex) %>% pull(pred)
icr_oos <- icr$pred %>% arrange(rowIndex) %>% pull(pred)
ppr_oos <- ppr$pred %>% arrange(rowIndex) %>% pull(pred)
bagEarth_oos <- bagEarth$pred %>% arrange(rowIndex) %>% pull(pred)
cubist_oos <- cubist$pred %>% arrange(rowIndex) %>% pull(pred)

train_stack <- data.frame(target = train$target,
                          RF_oos = rf_oos,
                          Enet_oos = enet_oos,
                          SVR_oos = svm_oos,
                          KNN_oos = knn_oos,
                          XGB_oos = xgb_oos,
                          ICR_oos = icr_oos,
                          PPR_oos = ppr_oos,
                          BagMars_oos = bagEarth_oos,
                          Cubist_oos = cubist_oos)

# correlation plot
ts2 <- train_stack
names(ts2) <- names(ts2) %>% str_replace_all("_oos", "")
corrplot(cor(ts2 %>% select(-target)), method = 'circle',
          tl.col = "black", tl.srt = 30, tl.cex = 1, mar = c(2,0,1,0))


# remove randomforest predictions b/c they are correlated with others
train_stack <- train_stack %>% select(-RF_oos)

# xgboost with oos predictions from other models as meta features
set.seed(888)
xgb_stack <- train(target ~ ., data = train_stack, tuneGrid = xgb_grid,
                   metric = "MAE", method = "xgbTree",trControl = ctrl)

xgb_stack$results %>% arrange(MAE) %>% head()

saveRDS(xgb_stack, paste0(data, "xgb_stack.rds"))

################################################################################
# Results on test set for best models

# xgboost MAE for test set
MAE(predict(xgb, test), test$target)
# svm MAE for test set
MAE(predict(svm, test), test$target)
