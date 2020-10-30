## ############################## ##
##
## Script Name: Random Forest Classification of digitised polygons
## 
## Author: JNCC
## 
## Date Created: 2019-11-15
## 
## Date Modified: 2020-10-30
## 
## Licence: MIT Licence
##
## 
## Abstract:
## functions to train RF model 'RFclass' and then to predict across variables.
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## caret_6.0-84        lattice_0.20-38    
## randomForest_4.6-14 dplyr_0.8.5         raster_3.0-7       
## sp_1.4-2  
##
## ############################## ##


#' #### Train the random forest classification ####
#'
#' @param train.df a dataframe containing the training data for the random forest classification
#' @param vars filepath for the rasterstack of predictor variables
#' @param varnames vector containing the variable names within the stack in order
#' @param prescol column noting the class in the training data
#' @param max_tries number of iterations to run through the model
#' @param resamp cv resample iterations.
#' @param prop.test proportion of data to set aside for testing
#' @param out.folder folder to save outputs
#'
#' @return
#' @export
#'
#' @examples
#' 
RFclass <- function(train.df,vars, varnames, prescol="barecat",max_tries=5,resamp=5,prop.test=0.25,out.folder){
  
  #load packages
  library(raster)
  library(dplyr)
  library(randomForest)
  library(caret)
  
  
  #set up
  if (!file.exists(paste(out.folder,"Outputs",sep=""))){
    dir.create(paste(out.folder,"Outputs",sep=""))
  }
  if (!file.exists(paste(out.folder,"Outputs/modpredict",sep=""))){
    dir.create(paste(out.folder,"Outputs/modpredict",sep=""))
  }
  
  ##### build varstack and check variable factors stay as factors ####
  #get var layers
  varstack <- stack(vars)
  names(varstack) <-varnames
  #check factors stay as factors
  list.factors <- names(varstack)[is.factor(varstack)]
  lvls <- NULL
  if(length(list.factors)>=1){
    for(k in 1:length(list.factors)){
      lvls[[list.factors[k]]] <- levels(varstack[[list.factors[k]]])[[1]]
    }
  }    
  varfactors <- names(varstack)
  
  #load in data
  trainvals <- read.csv(train.df)
  all_evals <- NULL
  all_imp <- NULL
  all_models <-NULL
  
  #### iterate over number of max_tries ####
  for (i in 1:max_tries){
  #### split out into training and test data ####
  training.data <- NULL
  test.data <- NULL
  # split out by category
  trainvals_split <- trainvals %>% dplyr::select(prescol,varfactors) %>% rename(response =prescol)
  ntrain <- trainvals_split %>% dplyr::group_by(response) %>% dplyr::summarise(count=n()) 
  samp <- min(ntrain$count)
  for(s in unique(ntrain$response)){
    traincatdf <- trainvals_split %>% dplyr::filter(response==s)
    traincatsamp <- traincatdf[sample(nrow(traincatdf),samp),]
    #split into training and test
    trainssamp <- sample(nrow(traincatsamp), ceiling(prop.test * (nrow(traincatsamp))))
    vals_train <- traincatsamp[-trainssamp,] 
    vals_test <- traincatsamp[trainssamp,] 
    training.data <- rbind(training.data, vals_train)
    test.data <- rbind(test.data, vals_test)
  }
  print(paste0("Training data: ", nrow(training.data), ". Test data: ", nrow(test.data)))
  #make sure classes are factors
  training.data$response <- as.factor(training.data$response)
  test.data$response <- as.factor(test.data$response)
 
  #### random forest classification ####
  print("Running random forest...")
  nameRF<- paste("RF",i,sep="")
  # tuning mtrys
  mtry <- tryCatch(randomForest::tuneRF(x = training.data[, 2:(ncol(training.data))], y = training.data$response, trace = FALSE,plot = FALSE), error = function(err) NA)
  best.m <- tryCatch(mtry[mtry[, 2] == min(mtry[, 2]), 1],error = function(err) NA)
  tunegrid <- expand.grid(.mtry=best.m)
  #train random forest model
  if (resamp==0){
    RF <- caret::train(response ~ ., training.data,method = "rf",tuneGrid=tunegrid,trControl = trainControl(method = "none",verboseIter = T))
  } else {
    RF <- caret::train(response ~ ., training.data,method = "rf",tuneGrid=tunegrid,trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = T))
  }
  model <- list(RF$finalModel)
  names(model)  <-  nameRF
  save(model, file = paste(out.folder, "Outputs/modpredict/RFmodel_",max_tries,"_",i, sep = ""))
  #model evaluation
  RF_pred <- predict(RF, test.data)
  RF_eval<-list(postResample(pred = RF_pred, obs = test.data$response)) #tset rmse
  names(RF_eval) <- nameRF
  #model importance values
  RFimp <-tibble::rownames_to_column(data.frame(randomForest::importance(RF$finalModel)),"var")
  names(RFimp) <- paste(nameRF,names(RFimp),sep="_")
  all_evals <- append(all_evals, RF_eval)
  all_imp <- append(all_imp, RFimp)
  all_models <- append(all_models, model)
  
  message("Run ", i, " completed.")
  }
  # save models and evals before prediction
  utils::write.csv(all_evals, file = paste(out.folder, "Outputs/evals_",max_tries,".csv", sep = ""))
  utils::write.csv(all_imp, file = paste(out.folder, "Outputs/imps_",max_tries,".csv", sep = ""))
  
  #get best model
  evals <- data.frame(mod=names(all_evals),t(data.frame(all_evals)))
  bestmod <- evals$mod[evals$Accuracy==max(evals$Accuracy)]
  RF <-all_models[[bestmod]]
  
   cat("Predicting random forest...                            \n ")
  RFpredimage <-predict(varstack, RF, filename=paste0(out.folder, "Outputs/Prediction_Class.tif"), progress='text', format='GTiff', datatype='FLT4S', type='response', overwrite=TRUE)
  RFpredimage <-predict(varstack, RF, filename=paste0(out.folder, "Outputs/Prediction_prob.tif"), progress='text', format='GTiff', datatype='FLT4S', type='prob', overwrite=TRUE)
  
  
  #plot output classification
  classout <- raster::raster(paste0(out.folder, "Outputs/Prediction_Class.tif"))
  plot(classout)
  #return class error matrix and importance
  print(RF)
  print(RF$importance)
  #write out info
  sink(paste0(out.folder, "Outputs/modinfo.txt"))
  print(RF)
  print(RF$importance)
  sink()
  
}

