## ############################## ##
##
## Script Name: Regression modelling 
## 
## Author: JNCC
## 
## Date Created: 2019-11-15
## 
## Date Modified:2020-10-30
## 
## Licence: MIT Licence
##
## 
## Abstract: Script to run through Bare peat regression modelling. the script lets you input a variety of models to test, the number of points to sample and iterations to resample data and run the models. The end result is a mean prediction across the model iterations as a raster tif.
## 
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
##  party_1.3-3          strucchange_1.5-2   
## sandwich_2.5-1       zoo_1.8-6           
## modeltools_0.2-22    mvtnorm_1.0-11      
## purrr_0.3.2          doParallel_1.0.14   
## iterators_1.0.10     foreach_1.4.4       
## quantregForest_1.3-7 RColorBrewer_1.1-2  
## brnn_0.7             Formula_1.2-3       
## caret_6.0-84         ggplot2_3.1.1       
## lattice_0.20-38      ranger_0.11.2       
## dismo_1.1-4          dplyr_0.8.1         
## rgdal_1.4-4          raster_3.0-7        
## randomForest_4.6-14  maptools_0.9-5      
## sp_1.3-1  
##
## ############################## ##

#' Regression modelling 
#'
#' @param training a dataframe of training data with extracted variable values
#' @param vars filepath for the rasterstack of predictor variables
#' @param varnames vector containing the variable names within the stack in order
#' @param out.folder filepath to save the outputs in
#' @param models vector containing models you wish to run out of c("RF","QRF","BRNN","SVM","CIF","BRT")
#' @param max_tries number of iterations to run through the model
#' @param prop.test proportion of data to set aside for testing
#' @param nsamp number of points to sample from the training dataframe
#' @param resamp cv resample iterations.
#'
#' @return
#' @export
#'
#' @examples
#' 
RFReg <- function(training, vars,varnames,out.folder,models=c("RF","QRF","BRNN","SVM","CIF","BRT"),max_tries=5,prop.test=0.25,nsamp=1000,resamp=5,stratified=T,fillsamp=F){
  
  #Load libraries
  require(maptools)
  require(sp)
  require(randomForest)
  require(raster)
  require(rgdal)
  library(dplyr)
  library(dismo)
  library(ranger)
  library(caret)
  library(brnn)
  library(quantregForest)
  library(doParallel)
  library(purrr)
  library(party)
  library(bst)
  
  #set up
  all_evals <- NULL
  all_imp <- NULL
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

  #get training data
  trainvals <- read.csv(training, header=TRUE) 
  
  #### iterate over number of max_tries ####
  for (i in 1:max_tries){
    ptm <- proc.time()
    
    #### split out into training and test data ####
    if (stratified==F){
      trainvals <- trainvals %>% dplyr::select("bare",varfactors)
      if(!is.na(nsamp)){ # if a number to sample is given sample that number otherwise all values
        trainno <- sample(nrow(trainvals),nsamp)
        trains <- trainvals[trainno,]
        } else{
          trains <- trainvals
          }
     #split into training and test
      trainssamp <- sample(nrow(trains), ceiling(prop.test * (nrow(trains))))
      training.data <- trains[-trainssamp,] %>% rename(response = bare)
      test.data <- trains[trainssamp,] %>% rename(response = bare)
      print(paste("Training data:", nrow(training.data),"Test.data",nrow(test.data)))
      
      } else {
      #sample by category 
        trainvals <- trainvals %>% dplyr::select("bare","barecat",varfactors)
        ntrain <- trainvals %>% dplyr::group_by(barecat) %>% dplyr::summarise(count=n())
        nsampdf <- data.frame(barecat=ntrain$barecat,totn=ntrain$count,nsamp=c(nsamp/length(unique(trainvals$barecat)))) %>% dplyr::mutate(sampled=ifelse(totn>nsamp,nsamp,totn))

        #update sampling numbers if fill in data selected
        if(fillsamp==T ){
          while(sum(nsampdf$sampled) < nsamp){
            remaining <- sum(nsampdf$nsamp-nsampdf$sampled)
            catnotfull<-nsampdf %>% dplyr::filter(nsamp==sampled) %>% dplyr::select(barecat)
            addval <- round(remaining/length(catnotfull$barecat),0) #get number to add on
            nsampdf <- nsampdf %>% mutate(nsamp=ifelse(nsampdf$barecat %in% catnotfull$barecat,nsamp+addval,nsamp))
            #see if sample size available
            nsampdf <- nsampdf %>% dplyr::mutate(sampled=ifelse(totn>nsamp,nsamp,totn))

          }
        }

        #sample data        
        training.data <- NULL
        test.data <- NULL
        for(s in unique(nsampdf$barecat)){
          traincatdf <- trainvals %>% dplyr::filter(trainvals$barecat==s)
          samp <- nsampdf %>% dplyr::filter(nsampdf$barecat==s) %>% dplyr::select(sampled) %>% as.numeric()
          traincatsamp <- traincatdf[sample(nrow(traincatdf),samp),]
          #split into training and test
          trainssamp <- sample(nrow(traincatsamp), ceiling(prop.test * (nrow(traincatsamp))))
          vals_train <- traincatsamp[-trainssamp,] %>% dplyr::rename(response = bare) %>% dplyr::select(-barecat)
          vals_test <- traincatsamp[trainssamp,] %>% dplyr::rename(response = bare) %>% dplyr::select(-barecat)
          training.data <- rbind(training.data, vals_train)
          test.data <- rbind(test.data, vals_test)
        }
        
        
    
     print(paste("Training data:", nrow(training.data),"Test.data",nrow(test.data)))
    all.data <- rbind(training.data,test.data)
    #record sampling information
    sink(paste0(out.folder, "Outputs/modpredict/trainingdata_",max_tries,"_",i,".txt"))
    cat("Sampled from binned data:\n")
    print(nsampdf)
    sink()
    #plot histogram of training data
    hist(all.data$response)
    ggplot2::ggplot(all.data, ggplot2::aes(response)) +
      ggplot2::geom_histogram(bins=10)
    ggplot2::ggsave(paste0(out.folder, "Outputs/modpredict/trainingdata_",max_tries,"_",i,".png"))
    }
    
    
    #### Run Random Forest ####
    if ("RF" %in% models){
      cat("tuning mtrys...                                      ")
      nameRF<- paste("RF",i,sep="")
      mtry <- tryCatch(randomForest::tuneRF(x = training.data[, 2:(ncol(training.data))], y = training.data$response, trace = FALSE,plot = FALSE), error = function(err) NA)
      best.m <- tryCatch(mtry[mtry[, 2] == min(mtry[, 2]), 1],error = function(err) NA)
      tunegrid <- expand.grid(.mtry=best.m)
      cat("Training random forest...                             \n ")
      if (resamp==0){
      RF <- caret::train(response ~ ., training.data,method = "rf",tuneGrid=tunegrid,trControl = trainControl(method = "none",verboseIter = T))
      } else {
        RF <- caret::train(response ~ ., training.data,method = "rf",tuneGrid=tunegrid,trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = T))
      }
      model <- list(RF$finalModel)
      names(model)  <- nameRF
      save(model, file = paste(out.folder, "Outputs/modpredict/RFmodel_",max_tries,"_",i, sep = ""))
      #model evaluation
      RF_pred <- predict(RF, test.data)
      RF_eval<-list(postResample(pred = RF_pred, obs = test.data$response)) #tset rmse
      names(RF_eval) <- nameRF
      utils::write.csv(RF_eval, file = paste(out.folder, "Outputs/modpredict/RFevals_",max_tries,"_",i,".csv", sep = ""))
      #model importance values
      RFimp <-tibble::rownames_to_column(data.frame(randomForest::importance(RF$finalModel)),"var")
      names(RFimp) <- paste(nameRF,names(RFimp),sep="_")
      utils::write.csv(RFimp, file = paste(out.folder, "Outputs/modpredict/RFimps_",max_tries,"_",i,".csv", sep = ""))
      #model prediction
      cat("Predicting random forest...                            \n ")
      RFpredimage <-predict(varstack, RF$finalModel, filename=paste(out.folder, "Outputs/modpredict/RFpredout",i,".tif",sep=""), progress='text', format='GTiff', datatype='FLT4S', type='response', overwrite=TRUE)
      all_evals <- append(all_evals, RF_eval)
      all_imp <- append(all_imp, RFimp)
      message("Random forest Regression model run ",i, " completed.") 
    }
    
    #### Run Quantile regression ####
    if ("QRF" %in% models){
      cat("tuning mtrys...                                        \n ")
      nameQRF<- paste("QRF",i,sep="")
      mtry <- tryCatch(randomForest::tuneRF(x = training.data[, 2:(ncol(training.data))], y = training.data$response, trace = FALSE,plot = FALSE), error = function(err) NA)
      best.m <- tryCatch(mtry[mtry[, 2] == min(mtry[, 2]), 1],error = function(err) NA)
      tunegrid <- expand.grid(.mtry=best.m)
      cat("Training Quantile regression forest...                  \n ")
      if (resamp==0){
      QRF <- caret::train(response ~ ., training.data,method = "qrf",tuneGrid=tunegrid,trControl = trainControl(method = "none",verboseIter = T))
      } else {
        QRF <- caret::train(response ~ ., training.data,method = "qrf",tuneGrid=tunegrid,trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = T))
      }
      modelQ <- list(QRF$finalModel)
      names(modelQ)  <- nameQRF
      save(modelQ, file = paste(out.folder, "Outputs/modpredict/QRFmodel_",max_tries,"_",i, sep = ""))
      #model evaluation
      QRF_pred <- predict(QRF, test.data)
      QRF_eval<-list(postResample(pred = QRF_pred, obs = test.data$response)) #tset rmse
      names(QRF_eval) <- nameQRF
      utils::write.csv(QRF_eval, file = paste(out.folder, "Outputs/modpredict/QRFevals_",max_tries,"_",i,".csv", sep = ""))
      #model importance values
      QRFimp <-tibble::rownames_to_column(data.frame(randomForest::importance(QRF$finalModel)),"var")
      names(QRFimp) <- paste(nameQRF,names(QRFimp),sep="_")
      utils::write.csv(QRFimp, file = paste(out.folder, "Outputs/modpredict/QRFimps_",max_tries,"_",i,".csv", sep = ""))
      # model prediction
      cat("Predicting Quantile regregression forest...               \n ")
      QRFpredimage <-predict(varstack,QRF,what = c(0.1, 0.5, 0.9),filename=paste(out.folder, "Outputs/modpredict/QRFpredout",i,".tif",sep=""), progress='text', format='GTiff',overwrite=TRUE)
      all_evals <- append(all_evals, QRF_eval)
      all_imp <- append(all_imp, QRFimp)
      message("Quantile forest Regression model run ",i, " completed.")
    }
    
    #### Bayesian Regularized Neural Networks ####
    if ("BRNN" %in% models){
      cat("Calculating Bayesian Regularized Neural Networks...         \n ")
      nameBRNN<- paste("BRNN",i,sep="")
      if (resamp==0){
        BRNN <- caret::train(response ~ ., training.data,method = "brnn",trControl = trainControl(method = "none",verboseIter = F),verbose=F)
      } else {
      BRNN <- caret::train(response ~ ., training.data,method = "brnn",trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = F),verbose=F)
      }
      modelNN <- list(BRNN$finalModel)
      names(modelNN)  <- nameBRNN
      save(modelNN, file = paste(out.folder, "Outputs/modpredict/BRNNmodel_",max_tries,"_",i, sep = ""))
      #model evaluation
      BRNN_pred <- predict(BRNN, test.data)
      BRNN_eval<-list(postResample(pred = BRNN_pred, obs = test.data$response)) #tset rmse
      names(BRNN_eval) <- nameBRNN
      utils::write.csv(BRNN_eval, file = paste(out.folder, "Outputs/modpredict/BRNNevals_",max_tries,"_",i,".csv", sep = ""))
      # model prediction
      cat("Predicting Bayesian Regularized Neural Networks...              ")
      BRNNpredimage <-predict(varstack,BRNN,filename=paste(out.folder, "Outputs/modpredict/BRNNpredout",i,".tif",sep=""), progress='text', format='GTiff',   overwrite=TRUE)
      all_evals <- append(all_evals, BRNN_eval)
      message("Bayesian Regularized Neural Networks model run ",i, " completed.")
    }
    
    #### Support Vector Machines ####
  if ("SVM" %in% models){
    cat("Calculating Support Vector Machines...              ")
    nameSVM<- paste("SVM",i,sep="")
    if(resamp==0){
      SVM <- caret::train(response ~ ., training.data,method = "svmLinear",trControl = trainControl(method = "none",verboseIter = F),verbose=F)
    } else{
      SVM <- caret::train(response ~ ., training.data,method = "svmLinear",trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = F),verbose=F)
    }
    modelSVM <- list(SVM$finalModel)
    names(modelSVM)  <- nameSVM
    save(modelSVM, file = paste(out.folder, "Outputs/modpredict/SVMmodel_",max_tries,"_",i, sep = ""))
    #model evaluation
    SVM_pred <- predict(SVM, test.data)
    SVM_eval<-list(caret::postResample(pred = SVM_pred, obs = test.data$response)) #test rmse
    names(SVM_eval) <- nameSVM
    utils::write.csv(SVM_eval, file = paste(out.folder, "Outputs/modpredict/SVM_eval_",max_tries,"_",i,".csv", sep = ""))
    all_evals <- append(all_evals, SVM_eval)
    # model prediction
    cat("Predicting Support Vector Machines ...              ")
    SVMpredimage <-predict(varstack,SVM,filename=paste(out.folder, "Outputs/modpredict/SVMpredout",i,".tif",sep=""), progress='text', format='GTiff',   overwrite=TRUE)
    message("Support Vector Machines model run ",i, " completed.")
  }
    
    #### Conditional Inference Forest ####
    if ("CIF" %in% models){
      cat("Calculating conditional inference forest...              ")
      nameCIF<- paste("CIF",i,sep="")
      if(resamp==0){
        CIF <- caret::train(response ~ ., training.data,method = "cforest",trControl = trainControl(method = "none",verboseIter = F))
      } else{
        CIF <- caret::train(response ~ ., training.data,method = "cforest",trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = F))
      }
      modelCIF <- list(CIF$finalModel)
      names(modelCIF)  <- nameCIF
      save(modelCIF, file = paste(out.folder, "Outputs/modpredict/CIFmodel_",max_tries,"_",i, sep = ""))
      #model evaluation
      CIF_pred <- predict(CIF, test.data)
      CIF_eval<-list(caret::postResample(pred = CIF_pred, obs = test.data$response)) #tset rmse
      names(CIF_eval) <- nameCIF
      utils::write.csv(CIF_eval, file = paste(out.folder, "Outputs/modpredict/CIF_eval_",max_tries,"_",i,".csv", sep = ""))
      all_evals <- append(all_evals, CIF_eval)
      # model prediction
      cat("Predicting conditional inference forest ...              ")
      CIFpredimage <-predict(varstack,CIF,filename=paste(out.folder, "Outputs/modpredict/CIFpredout",i,".tif",sep=""), progress='text', format='GTiff',   overwrite=TRUE)
      message("conditional inference forest model run ",i, " completed.")
    }
    
    #### Boosted Regression Tree ####
    if ("BRT" %in% models){
      cat("Calculating boosted regression tree...              ")
      nameBRT<- paste("BRT",i,sep="")
      if(resamp==0){
        BRT <- caret::train(response ~ ., training.data,method = "bstTree",trControl = trainControl(method = "none",verboseIter = F))
      } else{
        BRT <- caret::train(response ~ ., training.data,method = "bstTree",trControl = trainControl(method = "repeatedcv", number = resamp,repeats = resamp,verboseIter = F))
      }
      modelBRT <- list(BRT$finalModel)
      names(modelBRT)  <- nameBRT
      save(modelBRT, file = paste(out.folder, "Outputs/modpredict/BRTmodel_",max_tries,"_",i, sep = ""))
      #model evaluation
      BRT_pred <- predict(BRT, test.data)
      BRT_eval<-list(caret::postResample(pred = BRT_pred, obs = test.data$response)) #tset rmse
      names(BRT_eval) <- nameBRT
      utils::write.csv(BRT_eval, file = paste(out.folder, "Outputs/modpredict/BRT_eval_",max_tries,"_",i,".csv", sep = ""))
      all_evals <- append(all_evals, BRT_eval)
      # model prediction
      cat("Predicting boosted regression tree ...              ")
      BRTpredimage <-predict(varstack,BRT,filename=paste(out.folder, "Outputs/modpredict/BRTpredout",i,".tif",sep=""), progress='text', format='GTiff',   overwrite=TRUE)
      message("boosted regression tree model run ",i, " completed.")
    }
    
    message("Run ", i, " completed.")
    print(proc.time() - ptm)
    gc()
  }
  
  # save models and evals before prediction
  utils::write.csv(all_evals, file = paste(out.folder, "Outputs/evals_",max_tries,".csv", sep = ""))
  utils::write.csv(all_imp, file = paste(out.folder, "Outputs/imps_",max_tries,".csv", sep = ""))

  #### mean across the predictions ####
  if(max_tries>1){
   predimages <- list.files(paste(out.folder, "Outputs/modpredict/",sep=""))
   predimages <- predimages[which(grepl(".tif",predimages))] 
    for (j in 1:length(models)){
      mod <- models[j]
      modimages <- predimages[which(grepl(mod,predimages))] 
      if(mod =="RF"){
        modimages <- modimages[which(!grepl("QRF",modimages))] 
        }
      imagestack <- raster::raster(paste(out.folder, "Outputs/modpredict/",modimages[1],sep=""))
      for (k in 2:length(modimages)){
        image <- raster::raster(paste(out.folder, "Outputs/modpredict/",modimages[k],sep=""))
        imagestack <- raster::addLayer(imagestack,image)
        }
      meanpred<- raster::calc(imagestack,mean,progress='text')
      raster::writeRaster(meanpred,paste(out.folder, "Outputs/pred",mod,".tif",sep = ""),overwrite=T)
      } 
  }
  print("mean predictions done.")
}



