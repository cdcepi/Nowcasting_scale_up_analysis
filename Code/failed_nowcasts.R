library(tidyverse)
library(here)
library(randomForest)
library(caret)

load(here("Data/all_nobbs_mean.Rdata"))
dat <-all_nobbs_mean %>%
  mutate(failed_log_score=ifelse(failed_log_score=="Log scores <= -9.2", 1, 0))%>% #, 
  ungroup()%>%
  dplyr::select(pe_full, pe_part,  
                failed_log_score,
                range_50_cases, range_95_cases,
                mean_count)%>%
  drop_na()

## Random Forest
set.seed(100)
train <-sample(nrow(dat), 0.8*nrow(dat), replace=FALSE) 

train.data <-dat[train,]
test.data <-dat[-train,]

model <-randomForest(formula=as.factor(failed_log_score) ~ ., 
                     data=train.data, importance=TRUE, 
                     ntree=1000)
saveRDS(model, file="./Data/RF_full_model.RDS")

print(model) 
varImpPlot(model)
importance(model) 
class_prediction <-predict(object=model,  
                           newdata=test.data,  
                           type="class") #Return classification labels


#Confusion matrix for the test set
cm <-confusionMatrix(data=class_prediction, #Predicted classes
                     reference=as.factor(test.data$failed_log_score), #Actual classes
                     positive="1")  

print(cm)
probs <-predict(model, data.frame(test.data), type="prob")[,2]
ROC <-pROC::roc(response = test.data$failed_log_score, predictor = probs)
ROC$auc

######################################################################################
######################################################################################

# Data management to address correlation with high counts
dat_scaled <-all_nobbs_mean %>%
  mutate(failed_log_score=ifelse(failed_log_score=="Log scores <= -9.2", 1, 0))%>% #, 
  ungroup()%>%
  dplyr::select(pe_full, pe_part,  
                failed_log_score,
                range_50, range_95,
                range_50_cases, range_95_cases,
                mean_count)%>%
  drop_na() %>%
  mutate(range_50=ifelse(range_50==0, 0.33, range_50),
         range_50=log(range_50), 
         range_95=log(range_95), 
         mean_count_log=log(mean_count),
         mean_count_log_std=(log(mean_count)-mean(log(mean_count)))/sd(log(mean_count)),
         range95_count_log_std=(log(range_95_cases)-mean(log(range_95_cases)))/sd(log(range_95_cases)),
         range_95_cases_log=log(range_95_cases)) 

train_scaled <-sample(nrow(dat_scaled), 0.8*nrow(dat_scaled), replace=FALSE) #80/20 split for testing and training

train.data_scaled <-dat_scaled[train_scaled,]
test.data_scaled <-dat_scaled[-train_scaled,]

## Regression model
mod <-glm(failed_log_score ~ pe_part + range_95_cases_log + mean_count_log,
           data=train.data_scaled, family = binomial)
summary(mod)
probs_reg <-predict(mod, data.frame(test.data_scaled), type="response")

ROC_reg <-pROC::roc(response = test.data_scaled$failed_log_score, predictor = probs_reg)
ROC_reg$auc

class_prediction_reg <-predict(object=mod,
                           newdata=test.data_scaled,
                           type="response")
class_prediction_reg <-as.factor(ifelse(class_prediction_reg <= 0.5, 0, 1))

#Confusion matrix for the test set
cm_reg <-confusionMatrix(data=class_prediction_reg, #Predicted classes
                     reference=as.factor(test.data_scaled$failed_log_score), #Actual classes
                     positive="1")
print(cm_reg)
saveRDS(mod, "./Data/glm_part.RDS")
