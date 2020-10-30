## ############################## ##
##
## Script Name: Results evaluation and visulisation
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
## Abstract: Functions to evaluate model performance of classification and regression models and generate plots. Functions include:
## -  spec.plot        - to generate a plot of spectral information vs data classes
## -  ind.comp.plot   - create plots and perform correlation tests between the indices results for APGB vs Sentinel imagery
## -  ind.site.plot    - create plots and perform correlation tests to compare the indices results between different sites
## -  mod.eval         - script to evaluate a models performance across multiple methods and model runs
## -  imp.eval         - script to evaluate variable importance across multiple methods and model runs
## -  extrap.eval      - to evaluate accuracy of extrapolated data in non-training regions
## -  pnt_analysis     -  to extracting point data from imagery tiles and indices
## -  extrap_rmse      -  to evaluate the accuracy of data in extrapolated regions using continuous values and rmse statistic
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## sf_0.8-0      raster_3.0-7 
## sp_1.3-1      knitr_1.23   
## dplyr_0.8.1   stringr_1.4.0
## ggplot2_3.1.1 plyr_1.8.4   
## tidyr_0.8.3   cowplot_1.0.0
##
## ############################## ##


 #### Plot summary of data - spec bands ####
#'
#' @param dat dataframe of spectral bands and response
#' @param resp string, name of the classes to plot band response against
#' @param vars string, names of variables
#'
#' @return
#' @export
#'
#' @examples
spec.plot <- function(dat, resp,vars){
  #load dependencies
  library(cowplot)
  library(tidyr)
  library(plyr)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  #read data and format df
  dat.df <- read.csv(dat) 
  trim.df <- dat.df  %>%   dplyr::mutate_if(is.character, stringr::str_trim) %>%
    dplyr::mutate_if(is.character,as.numeric) 
  trim.df[[resp]] <- as.factor(trim.df[[resp]])
  names(trim.df)[which(names(trim.df)==resp)] <- "class"
  #step is specific to colnames
   new_df <- trim.df %>% dplyr::group_by(class) %>% plyr::summarise(mean = mean(3))
   alldf <-data.frame()
   #get stats per class
   for (i in 1:length(vars)){
     varname <- vars[i]
      vardf<- tapply(trim.df[[vars[i]]], trim.df$class, function(x){c(MEAN=mean(x), MIN=min(x),MAX=max(x) )},simplify = F)
     for (j in 1:length(vardf)){
       statdf <- raster::as.data.frame(vardf[j]) %>% t() %>% data.frame() %>% dplyr::mutate(var=varname) %>% dplyr::mutate(class=names(vardf)[j]) 
       alldf <- rbind(alldf,statdf)
     }
   }
   alldf$class <- plyr::revalue(alldf$class,c("0"="vegetated", "1"="bare"))
   #generate plot
   ggplot2::ggplot(alldf, ggplot2::aes(x=var, y=MEAN,group=class, colour=class,fill=class)) +
     ggplot2::geom_ribbon(ggplot2::aes(ymin=MIN,ymax=MAX),alpha=.3,colour=NA,show.legend = T)+ 
     ggplot2::geom_line(ggplot2::aes(y=MEAN)) +
     ggplot2::ylim(0,max(alldf$MAX)) +
     ggplot2::scale_fill_manual(values=c(vegetated = "#228B22",bare = "#8B4513")) +
     ggplot2::scale_colour_manual(values=c(vegetated = "#228B22",bare = "#8B4513")) + 
     ggplot2::labs(group=NULL,colour="",x="Band",y="Reflectance") 
   return(alldf)
}
  
   

#### correlation test and t tests ####
#'
#' @param APGB  path to the APGB data results
#' @param Sentinel path to the Sentinel data results
#' @param indices indices to produce comparisons for
#' @param out.path output path to store results
#' @param suffix naming convention
#'
#' @return
#' @export
#'
#' @examples
#' 
ind.comp.plot <- function(APGB,Sentinel,indices,out.path,suffix="train"){
  #load dependencies
  library(dplyr)
  library(plyr)
  library(ggplot2)
  
  #load in the data and process for comparison
  APGB_train <- read.csv(APGB) %>% dplyr::mutate(Imagery="APGB") %>% dplyr::mutate(class=as.factor(class))
  Sent_train <- read.csv(Sentinel) %>% dplyr::mutate(Imagery="Sentinel") %>% dplyr::mutate(class=as.factor(class))
  alldat <- rbind(APGB_train, Sent_train)
  #rename for displaying bare and vegetated
  alldat$class <- plyr::revalue(alldat$class,c("0"="vegetated", "1"="bare"))
  #group for the plots
  alldat$Group <- paste(alldat$Imagery,alldat$class)
  
  #iterate through the indices 
  for (i in 1:length(indices)){
    
    #create plots
    ggplot2::ggplot(alldat, ggplot2::aes_string(x="Imagery",y=indices[i] ,group="Group",colour="Group")) +
      ggplot2::geom_boxplot(outlier.shape = NA,show.legend=FALSE) +
      ggplot2::scale_colour_manual(values=c("APGB bare" = "#8B4513","APGB vegetated" = "#228B22","Sentinel vegetated" = "#228B22","Sentinel bare" = "#8B4513")) 
    ggplot2::ggsave(paste0(out.path,indices[i],"_",suffix,"comparison.png"))
    
    #group for generating summary stats
    summstat <- alldat %>% dplyr::select(Group,index=indices[i]) %>% 
      dplyr::group_by(Group) %>% 
      dplyr::summarise(mean = mean(index), min=min(index), max=max(index), sd=stats::sd(index), IQR = stats::IQR(index), Q1 = stats::quantile(index, 0.25), Q3 = stats::quantile(index, 0.75), var=stats::var(index))
    
    ## APGB bare vs veg
    APGB.ind <- APGB_train %>% dplyr::select(class,index=indices[i]) 
    vegbare_APGB <- stats::t.test(index ~ class, data = APGB.ind, paired=T)
    
    ## Sent bare vs veg
    Sent.ind <- Sent_train %>% dplyr::select(class,index=indices[i]) 
    Sent_APGB <- stats::t.test(index ~ class, data = Sent.ind, paired=T) 
    
    # correlation sent vs apgb
    corrindex <- stats::cor.test(APGB.ind$index, Sent.ind$index, method = "pearson")
    
    #t-test APGB vs Sentinel
    all.img <- alldat %>% dplyr::select(Imagery,index=indices[i]) 
    all.img.test <- stats::t.test(index ~ Imagery, data = all.img, paired=T)
    #bare peat 
    bare.img <- alldat %>% dplyr::filter(class=="bare") %>% select(Imagery,index=indices[i]) 
    bare.img.test <- stats::t.test(index ~ Imagery, data = bare.img, paired=T) 
    #vegetated peat
    veg.img <- alldat %>% dplyr::filter(class=="vegetated") %>% dplyr::select(Imagery,index=indices[i]) 
    veg.img.test <- stats::t.test(index ~ Imagery, data = veg.img, paired=T) 
    
    #list all stats
    tests <- list(df_stats = summstat,APGB_bare_veg = vegbare_APGB,Sentinel_bare_veg = Sent_APGB,correlation_APGB_Sentinel=corrindex, all_APGB_Sentinel = all.img.test, bare_APGB_Sentinel = bare.img.test, veg_APGB_Sentinel = veg.img.test)
    
    # write outputs to a text file
    sink(paste0(out.path,indices[i],"_",suffix,"stat.txt"))
    print(tests)
    sink()
  }
}

#### comparisons of indices between different sites ####
#'
#' @param img.dat path to results data
#' @param indices indices to draw comparisons on
#' @param out.path output path
#'
#' @return
#' @export
#'
#' @examples
#' 
ind.site.plot <- function(img.dat,indices, out.path){
  
  #load dependencies
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  #get all the data
  alldat <- NULL
  for (i in 1: length(img.dat)){
    img.df <- read.csv(img.dat[[i]])
    img.df$Imagery <- names(img.dat[i])
    alldat <- rbind(alldat, img.df)
  }
  alldat$class <- as.factor(alldat$class)
  alldat$class <- plyr::revalue(alldat$class,c("0"="vegetated", "1"="bare"))
  
  #iterate through each indices
  for (j in 1:length(indices)){
    ##  plot per indices
    #overall value comparison
    ggplot2::ggplot(alldat, ggplot2::aes_string(x="Imagery",y=indices[j])) +
      ggplot2::geom_boxplot(outlier.shape = NA,show.legend=FALSE) 
    ggplot2::ggsave(paste0(out.path,indices[j],"_overall_comparison.png"))
    #bare vs veg
    g1 <- alldat %>% dplyr::mutate(Group = paste0(Imagery,"_",class))
    ggplot2::ggplot(g1, ggplot2::aes_string(x="Imagery",y=indices[j],group="Group",color="Group")) +
      ggplot2::geom_boxplot(outlier.shape = NA,show.legend=FALSE) +
      ggplot2::scale_colour_manual(values=c("train_bare" = "#8B4513","train_vegetated" = "#228B22","evalshp1_vegetated" = "#228B22","evalshp1_bare" = "#8B4513","evalshp2_vegetated" = "#228B22","evalshp2_bare" = "#8B4513")) 
    ggplot2::ggsave(paste0(out.path,indices[j],"_bareVveg_comparison.png"))
    
    ## summary stats
    summstat <- g1 %>% dplyr::select(Group,index=indices[j]) %>% dplyr::group_by(Group) %>% dplyr::summarise(mean = mean(index), min=min(index), max=max(index), sd=stats::sd(index), IQR = stats::IQR(index), Q1 = stats::quantile(index, 0.25), Q3 = stats::quantile(index, 0.75), var=stats::var(index)) %>% dplyr::mutate(index=indices[j])
    
    #t-test train vs evalshp2
    all.img <- alldat %>% dplyr::select(Imagery,index=indices[j]) %>% dplyr:: filter(Imagery!="evalshp2")
    all.eval1.test <- stats::t.test(index ~ Imagery, data = all.img, paired=T) 
    bare.img <- alldat %>%  dplyr::filter(Imagery!="evalshp2") %>%  dplyr::filter(class=="bare") %>%  dplyr::select(Imagery,index=indices[j]) 
    bare.eval1.test <- stats::t.test(index ~ Imagery, data = bare.img, paired=T) 
    veg.img <- alldat %>%  dplyr::filter(Imagery!="evalshp2") %>%  dplyr::filter(class=="vegetated") %>%  dplyr::select(Imagery,index=indices[j]) 
    veg.eval1.test <- stats::t.test(index ~ Imagery, data = veg.img, paired=T) 
    
    #t-test train vs evalshp1
    all.img <- alldat %>% select(Imagery,index=indices[j]) %>%  dplyr::filter(Imagery!="evalshp1")
    all.eval2.test <- stats::t.test(index ~ Imagery, data = all.img, paired=T) 
    bare.img <- alldat %>%  dplyr::filter(Imagery!="evalshp1") %>%  dplyr::filter(class=="bare") %>%  dplyr::select(Imagery,index=indices[j]) 
    bare.eval2.test <- stats::t.test(index ~ Imagery, data = bare.img, paired=T) 
    veg.img <- alldat %>%  dplyr::filter(Imagery!="evalshp1") %>%  dplyr::filter(class=="vegetated") %>%  dplyr::select(Imagery,index=indices[j]) 
    veg.eval2.test <- stats::t.test(index ~ Imagery, data = veg.img, paired=T) 
    
    #list all stats
    tests <- list(df_stats = summstat,train_vs_evalshp1_all = all.eval1.test,train_vs_evalshp1_bare = bare.eval1.test,train_vs_evalshp1_veg = veg.eval1.test,train_vs_evalshp2_all = all.eval2.test,train_vs_evalshp2_bare = bare.eval2.test,train_vs_evalshp2_veg = veg.eval2.test)
    
    # write outputs to a text file
    sink(paste0(out.path,indices[j],"_stat.txt"))
    print(tests)
    sink()
  }
  
}


 #### looking at model performance ####
#'
#' @param results.fold link to folder where results are stored
#' @param max_tries number of max tries the model ran for
#'
#' @return
#' @export
#'
#' @examples
mod.eval <- function(results.fold,max_tries=5){
  #load dependencies
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(knitr)
  #read in data
  evals <- read.csv(paste(results.fold,"evals_",max_tries,".csv",sep=""))
  ## RMSE
  #tidy values
  RMSE <- evals %>% dplyr::filter(X=="RMSE") %>% tidyr::gather(mod,RMSE,2:ncol(evals)) %>% dplyr::select(-X) %>% dplyr::mutate(group=stringr::str_extract(mod, "[aA-zZ]+"))
  #create plot
  ggplot2::ggplot(RMSE, ggplot2::aes(x=group, y=RMSE, fill=group)) +
    ggplot2::geom_boxplot(show.legend = FALSE)  + 
    ggplot2::scale_fill_brewer(palette="Blues") + theme_classic()
  ggplot2::ggsave(paste(results.fold,"RMSE.PNG",sep=""),plot = last_plot())
  #find best model
  best <- RMSE[which(RMSE$RMSE==min(RMSE$RMSE)),]
  bestmods <- RMSE %>% dplyr::group_by(group) %>% dplyr::summarise(meangrp = mean(RMSE))
  print(paste("Best model:",best$mod,round(best$RMSE,3),". Mean RSME values per model:"))
  print(knitr::kable(bestmods))
  ## Rsquared
  #tidy values
  Rsq <- evals %>% dplyr::filter(X=="Rsquared") %>% tidyr::gather(mod,Rsquared,2:ncol(evals)) %>% dplyr::select(-X) %>% dplyr::mutate(group=stringr::str_extract(mod, "[aA-zZ]+"))
  #create plot
  ggplot2::ggplot(Rsq, ggplot2::aes(x=group, y=Rsquared, fill=group)) +
    ggplot2::geom_boxplot(show.legend = FALSE)  + 
    ggplot2::scale_fill_brewer(palette="Blues") + theme_classic()
  ggplot2::ggsave(paste(results.fold,"Rsq.PNG",sep=""),plot = last_plot())
  #find best model
  bestRsq <- Rsq[which(Rsq$Rsquared==max(Rsq$Rsquared)),]
  bestmodsRsq <- Rsq %>% dplyr::group_by(group) %>% dplyr::summarise(meangrp = mean(Rsquared))
  print(paste("Best model:",bestRsq$mod,round(bestRsq$Rsquared,3),". Mean Rsq values per model:"))
  print(knitr::kable(bestmodsRsq))
  #return results to environment
  list2env(list(RMSE=RMSE),envir = .GlobalEnv)
}


 #### importance values ####
#'
#' @param results.fold link to folder where results are stored
#' @param max_tries number of max tries the model ran for
#'
#' @return
#' @export
#'
#' @examples
imp.eval <- function(results.fold,max_tries=5){
  #load dependencies
  library(ggplot2)
  library(tidyr)
  library(stringr)
  #read in importance sheet
  imp <- read.csv(paste(results.fold,"imps_",max_tries,".csv",sep=""))
  imp$X <-imp[,2]
  #format df
  imp.df <- imp %>% dplyr::select(which(!grepl("var",names(imp)))) %>% 
    tidyr::gather(Model,IncNodePurity,2:((ncol(imp)+1)/2) ) %>% 
    dplyr::mutate(Model = stringr::str_replace(Model, pattern="_IncNodePurity", "")) %>% 
    dplyr::mutate(group=stringr::str_extract(Model, "[aA-zZ]+")) %>% 
    dplyr::rename(Variable=X)
  #generate plot
  ggplot2::ggplot(imp.df, ggplot2::aes(x=Variable, y=IncNodePurity, fill=group)) +
    ggplot2::geom_boxplot()  + 
    ggplot2::scale_fill_brewer(palette="Greens") + theme_classic()
  ggplot2::ggsave(paste(results.fold,"VarImportance.PNG",sep=""),plot = last_plot())
  # get most important per model
  mostimp <- imp.df %>% dplyr::group_by(group,Variable) %>% 
    dplyr::summarise(maxvar = max(IncNodePurity)) 
  QRFvar <- mostimp %>% dplyr::filter(group=="QRF") %>% dplyr::filter(maxvar==max(maxvar))
  RFvar <- mostimp %>% dplyr::filter(group=="RF") %>% dplyr::filter(maxvar==max(maxvar))
  print(paste("Most Important Variable: QRF=",QRFvar$Variable,QRFvar$maxvar,", RF=", RFvar$Variable, RFvar$maxvar))
}

## evaluate accuracy of extrapolated data in non-training regions ####
#'
#' @param results folder of tests and predicted modelled images
#' @param training bare peat training data; polygons or points, will use raster::extract function
#'@param threshval = threshold value to class predicted raster as bare peat present
#'@param out_folder = output folder
#' @return
#' @export
#'
#' @examples
extrap.eval <- function(results,training,threshval=0.1,out_folder){
  #load packages
  library(dplyr)
  library(raster)
  library(sf)
  library(caret)
  
  #store all results
  evaluation_results <- NULL
  
  ## prepare shapefile into points
  #read in shapefile
  shape <- sf::st_read(training,quiet=T)
  train_shp <- shape %>% dplyr::mutate(id=1:nrow(shape))
  
  ##get predicted rasters
  #get files or folders
  results_fold <- list.dirs(results,recursive=F,full.names=F)
  #iterating through each results folder
  for(i in 1:(length(results_fold))){
    foldername <- results_fold[i]
    if (!file.exists(paste0(results,foldername,"/predRF.tif"))){ 
      #if directory of test folders list all tests in folder
      test_fold <-list.files(paste0(results,foldername),full.names = T,recursive=F)
    } else {
      #otherwise folder is test folder
      test_fold <- paste0(results,foldername)
    }
    #iterate through each test folder
    for (j in 1:length(test_fold)){
      #read in predicted raster
      testfiles <- list.files(test_fold[j],include.dirs = F) 
      #get the output predicted maps
      pred_elements <- testfiles[which(grepl(pattern="pred",testfiles) & grepl(pattern="tif",testfiles))]
      #iterate through every predicted map
      for (k in 1:length(pred_elements)){
        #read in predicted map
        predimg <- raster::raster(paste0(test_fold[j],"/", pred_elements[k]))
        #extract values at the training points/polygons
        ind.vals<- raster::extract(predimg, train_shp, df=TRUE,progress='text')
        #convert to dataframe
        train_shp_info <- train_shp %>% sf::st_as_sf() %>% raster::as.data.frame() %>% dplyr::select(-geometry)
        all_dat <- merge(ind.vals,train_shp_info,by.x="ID",by.y="id")
        names(all_dat) <- c("ID","pred_bare","bare")
        #revalue those which are above your threshold as 1 (present), below as 0 (absent)
        new_dat <- all_dat %>% dplyr::mutate(pred_bare = ifelse(pred_bare>threshval,1,0)) %>% dplyr::mutate_if(is.numeric,as.factor)
        #generate confusion matrix to see if present=present, absent=absent
        conf <- caret::confusionMatrix(new_dat$pred_bare,new_dat$bare)
        overall <- conf$overall
        byclass <- conf$byClass
        #compile the results in a data frame
        eval_out <- data.frame(Category = foldername,
                               Test = basename(test_fold)[j],
                               Method = gsub(gsub(pred_elements[k],pattern=".tif",replacement=""),pattern = "pred",replacement=""),
                               Accuracy=overall[1],
                               Kappa=overall[2],
                               ConfInt_95_lower = overall[3],
                               ConfInt_95_upper = overall[4],
                               No_Info_Rate = overall[5],
                               AccuracyPValue=overall[6],
                               McnemarPValue=overall[7],
                               Sensitivity = byclass[1],
                               Specificity = byclass[2],
                               Precision = byclass[5],
                               Prevalence = byclass[8],
                               Detection_Rate = byclass[9],
                               Detection_Prevalence  = byclass[10],
                               Balanced_Accuracy= byclass[11],
                               Positive_class=conf[["positive"]]
        )
        #combine all the results from the different tests
        evaluation_results <- rbind(evaluation_results,eval_out)
        #report on progress
        print(paste("results folder ", i," : ","test ", j, " of ", length(test_fold), " : ", k, "of", length(pred_elements), "methods."))
      }
    }
  }
  #write out results
  write.csv(evaluation_results,paste0(out_folder,'evaluation_results.csv'))
  
}


#' Point analysis for extracting point data from imagery tiles and indices
#'
#' @param img_lst list of imagery files to iterate through
#' @param points_shp shapefile
#' @param name name to write out to
#' @param tiles_fold folder path to indices layers
#'
#' @return
#' @export
#'
#' @examples
pnt_analysis <- function(img_lst, points_shp,name, tiles_fold){
  
  library(raster)
  library(sf)
  library(dplyr)
  
  ## loop and run point analysis
  images_lst <- NULL
  for (i in 1:length(img_lst)){
    combo <- raster::brick(paste0(tiles_fold,img_lst[i],".tif"))
    bright <- raster::raster(paste0(tiles_fold,"Indices/",img_lst[i],"/Brightness/Brightness_",img_lst[i],".tif"))
    NDVI <- raster::raster(paste0(tiles_fold,"Indices/",img_lst[i],"/NDVI/NDVI_",img_lst[i],".tif"))
    RG <- raster::raster(paste0(tiles_fold,"Indices/",img_lst[i],"/RG/RG_",img_lst[i],".tif"))
    RB <- raster::raster(paste0(tiles_fold,"Indices/",img_lst[i],"/RB/RB_",img_lst[i],".tif"))
    allstk <- raster::stack(combo,bright,NDVI,RG,RB)
    images_lst <- append(images_lst,allstk)
  } 
  if(length(images_lst)==1){
    allimg <- images_lst[[1]]
  } else{
    names(images_lst) <- NULL
    images_lst$fun <- mean
    allimg <- do.call(raster::mosaic,images_lst)
  }
  # read in points
  eval_pnt <- sf::st_read(points_shp,quiet=T)
  #extract point layer
  vals <- raster::extract(allimg,eval_pnt,df=T)
  names(vals)<-c("ID","Red_band","Blue_band","Green_band","NIR_band","Brightness","NDVI","RG","RB")
  # compile statistics
  df <- data.frame(vals,bare=eval_pnt$bare) %>% dplyr::group_by(bare) %>% 
    dplyr::summarise(nir_min = min(NIR_band),nir_max = max(NIR_band),nir_1st = raster::quantile(NIR_band,0.25,na.rm=TRUE),nir_3rd = raster::quantile(NIR_band,0.75,na.rm=TRUE),                             brightness_min = min(Brightness),brightness_max = max(Brightness),brightness_1st = raster::quantile(Brightness,0.25,na.rm=TRUE),brightness_3rd = raster::quantile(NIR_band,0.75,na.rm=TRUE),
                     NDVI_min = min(NDVI),NDVI_max = max(NDVI),NDVI_1st = raster::quantile(NDVI,0.25,na.rm=TRUE),NDVI_3rd = raster::quantile(NDVI,0.75,na.rm=TRUE),
                     RG_min = min(RG),RG_max = max(RG),RG_1st = raster::quantile(RG,0.25,na.rm=TRUE),RG_3rd = raster::quantile(RG,0.75,na.rm=TRUE),
                     RB_min = min(RB),RB_max = max(RB),RB_1st = raster::quantile(RB,0.25,na.rm=TRUE),RB_3rd = raster::quantile(RB,0.75,na.rm=TRUE))
  if(!dir.exists(paste0(tiles_fold,name))){
    dir.create(paste0(tiles_fold,name))
  }
  write.csv(df,paste0(tiles_fold,name,"/",name,"_stats.csv"))
  return(allimg)
}



#' extrap_rmse function
#'
#' @param results.fold folder path where results from RFReg.R are saved
#' @param eval_dat evaluation data filepath
#' @param name filename to save results to
#'
#' @return
#' @export
#'
#' @examples
extrap_rmse <- function(results.fold,eval_dat,name="extrap_eval"){
  #load packages
  library(dplyr)
  library(caret)
  library(randomForest)
  library(raster)
  library(brnn)
  library(bst)
  library(quanregforest)
  #read in eval data
  eval_df <- read.csv(eval_dat)
  
  #get list of models
  modlist <- data.frame(file=list.files(paste0(results.fold,"modpredict"))) %>% filter(grepl(file,pattern="model"))
  eval_results <- NULL
  for (i in 1:nrow(modlist)){
    #load model
    load(paste0(results.fold,"modpredict/",as.character(modlist[i,])))
    if(exists("modelQ")){ #QRF model
      model <- modelQ[[1]]
      group="QRF"
      rm(modelQ)
    } else if(exists("modelSVM")){ #SVM model
      model <- modelSVM[[1]]
      group="SVM"
      rm(modelSVM)
    } else if (exists("modelNN")){ #BRNN model
      model <- modelNN[[1]]
      group="BRNN"
      rm(modelNN)
    } else if (exists("modelCIF")){ #CIF model
      model <- modelCIF[[1]]
      group="CIF"
      rm(modelCIF)
    } else if (exists("modelBRT")){ #BRT model
      model <- modelBRT[[1]]
      group="BRT"
      rm(modelBRT)
    } else{ #RF model
      model <- model[[1]]
      group="RF"
    }
    
    if (group=="SVM"){
      varnames = names( data.frame(model@xmatrix))
    } else {
      varnames = names(eval_df)[names(eval_df) %in% model$xNames]
    }
    
    eval_cols <- eval_df %>% dplyr::rename(response = bare) %>% dplyr::select(c("response",varnames))
    mod_pred <- predict(model, eval_cols[,2:(length(varnames)+1)])
    mod_eval<-list(postResample(pred = mod_pred, obs = eval_cols$response)) # rmse
    mod_df <-data.frame(Group=group,Run=basename(as.character(modlist[i,])),RMSE=mod_eval[[1]][1],Rsquared=mod_eval[[1]][2],MAE=mod_eval[[1]][3],row.names = NULL)
    eval_results <- rbind(eval_results,mod_df)
  }
  
  #find the best model results and print
  best <- eval_results[which(eval_results$RMSE==min(eval_results$RMSE)),]
  bestmods <- eval_results %>% dplyr::group_by(Group) %>% dplyr::summarise(meangrp = mean(RMSE))
  print(paste("Best model:",best$Run,round(best$RMSE,3),". Mean RSME values per model:"))
  print(knitr::kable(bestmods))
  print(paste("Best model:",best$Run,round(best$Rsquared,3),". Mean Rsquared values per model:"))
  bestrsq <- eval_results %>% dplyr::group_by(Group) %>% dplyr::summarise(meangrp = mean(Rsquared))
  print(knitr::kable(bestrsq))
  #create plot and save
  ggplot2::ggplot(eval_results, ggplot2::aes(x=Group, y=RMSE, fill=group)) +
    ggplot2::geom_boxplot(show.legend = FALSE)  + 
    ggplot2::scale_fill_brewer(palette="Blues") + theme_classic()
  ggplot2::ggsave(paste(results.fold,name,"_RMSE.PNG",sep=""),plot = last_plot())
  
  ggplot2::ggplot(eval_results, ggplot2::aes(x=Group, y=Rsquared, fill=group)) +
    ggplot2::geom_boxplot(show.legend = FALSE)  + 
    ggplot2::scale_fill_brewer(palette="Blues") + theme_classic()
  ggplot2::ggsave(paste(results.fold,name,"_Rsquared.PNG",sep=""),plot = last_plot())
  #write out results
  write.csv(eval_results,paste0(results.fold,name,".csv"))
} 



