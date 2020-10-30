## ############################## ##
##
## Script Name:script to assess bare vs vegetated against multiple indices using digitised polygons as training data 
## 
## Author: JNCC
## 
## Date Created: 2019-11-14
## 
## Date Modified: 2020-10-30
## 
## Licence: MIT Licence
##
## 
## Abstract: Script to extract information about the indices values for the multispectral imagery and compare results for given classes providing in the training data. Here training data can either be supplied as polygons, which will then be samples to retrieve points, or training points to extract from. 
## 
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## ggplot2_3.1.1   fasterize_1.0.0
## tidyr_0.8.3     rgdal_1.4-4    
## rgeos_0.4-3     sf_0.7-4       
## dplyr_0.8.1     raster_3.0-7   
## sp_1.3-1  
##
## ############################## ##

#' Indices profiling
#'
#' @param train.path path to the training polygons
#' @param NDVI.path optional - path to the NDVI file, values to extract values from, this is mostly as this is in a different folder for ease processing, if in same folder as indices path then leave null
#' @param Ind.path path to the imagery indices folder
#' @param out.path output folder to store results
#' @param nsamps number of sample points to assess, default will take minimum number from classes present
#' #' @param nf if indices were created with a new folder per image
#'
#' @return
#' @export
#'
#' @examples
#' 
indi.profile <- function(train.path, Ind.path, out.path, nsamps=NULL,nf=F){
  
  #load in dependencies
  library(sp)
  library(sf)
  library(rgeos)
  library(rgdal)
  library(dplyr)
  library(tidyr)
  library(raster)
  library(fasterize)
  library(ggplot2)
  
  #read in indices
  if(nf==T){
    imagery.list <- list.files(Ind.path)
  } else {
    imagery.list <- basename(Ind.path)
  }
  
  for (i in 1:length(imagery.list)){ #for every image 
    # get the indices folder and no. indices
    if(nf==T){
      Ind.list <- list.files(paste0(Ind.path,imagery.list[i]))
      images <- list.files(paste0(Ind.path,imagery.list[i],"/",Ind.list[1]))
    } else {
      Ind.list <- list.files(Ind.path)
      images <- list.files(paste0(Ind.path,"/",Ind.list[1]))
    }      
    # get list of images in the indices folder (where nf==F this will produce multiple images within indices folders if supplied)
    #images <- list.files(paste0(Ind.path,imagery.list[i],"/",Ind.list[1]))
    #iterate through images
    for (j in 1:length(images)){
      #get image name
      name <- images[j] %>% gsub(pattern=Ind.list[1],replacement="")
      #grab all indices for that image
      indat.list <- list()
      for (n in 1:length(Ind.list)){
        if(nf==T){
          ind.1 <- raster::raster(paste0(Ind.path,imagery.list[i],"/",Ind.list[n],"/",Ind.list[n],name))
        }else{
          ind.1 <- raster::raster(paste0(Ind.path,"/",Ind.list[n],"/",Ind.list[n],name))
        }
        item <- list(ind.1)
        names(item) <- Ind.list[n]
        indat.list <- append(indat.list, item)
        }
      #read in training polygons
      training <- rgdal::readOGR(train.path)
      if(class(training)=="SpatialPolygonsDataFrame"){
        train.sf <- sf::read_sf(train.path)
        #get extraction points
        r <- raster::raster(ext=raster::extent(indat.list[[1]]), res=raster::res(indat.list[[1]]))
        train.r <- fasterize::fasterize(train.sf, r, "bare")
        peat.out <- raster::overlay(train.r,indat.list[["NDVI"]],fun=function(x,y){ 
          x[is.na(y[])] <- NA
          return(x)
        })
        peat.points <- raster::rasterToPoints(peat.out, spatial = T, na.rm = T)
        #get no. to sample per class
        classes <- data.frame(class = unique(peat.points[[1]]), num = 0)
        for (r in 1:length(classes$class)) {
          cl <- classes$class[r]
          count <- nrow(peat.points[peat.points[[1]] == cl, ])
          classes$num[r] <- count
        }
        if(is.null(nsamps)){
          nsamps <- min(classes$num)
        }
        #vegetated
        class_data<- peat.points[peat.points$layer==classes$class[1],]
        veg_class <-  class_data[sample(1:length(class_data), nsamps), ]
        #bare
        bare_data<- peat.points[peat.points$layer==classes$class[2],]
        bare_class <-  bare_data[sample(1:length(bare_data), nsamps), ]
        allpoints <- rbind(veg_class,bare_class)
      } else{
        allpoints <- training
      }
      
    # extract the indices values
    if(class(training)=="SpatialPolygonsDataFrame"){
      classVal <- raster::extract(peat.out, allpoints, df=TRUE)
      } else {
        allpoints$id <- 1:nrow(allpoints)
        classVal <- data.frame(allpoints$id,allpoints$layer)
      }
      
      for (k in 1:length(indat.list)){
        ind <- indat.list[[k]]
        ind.vals<- raster::extract(ind, allpoints, df=TRUE)
        colnames(ind.vals)[2] <- names(indat.list[k])
        classVal <- cbind(classVal,ind.vals[,2])
        print(paste(k,"of",length(indat.list)))
      }

  colnames(classVal) <- c("ID", "class",as.character(names(indat.list)))
  indi.data <- classVal[,2:ncol(classVal)]
  indi.data$class <- as.factor(indi.data$class)
  
  #write out
  if (!file.exists(paste(out.path,"Indices_analysis",sep=""))){
    dir.create(paste(out.path,"Indices_analysis",sep=""))
  }
  
  newname <- gsub(name,pattern=".tif",replacement="")
  #save training points
  if(class(training)=="SpatialPolygonsDataFrame"){
    raster::shapefile(allpoints, paste(out.path,"Indices_analysis/train",newname,".shp",sep=""),overwrite=T)
  }
  #write out results
  write.csv(indi.data, paste(out.path,"Indices_analysis/resp",newname,".csv",sep=""))
  
  list2env(list(indi.data=indi.data),.GlobalEnv)
      
      
    }
    
  }
}

