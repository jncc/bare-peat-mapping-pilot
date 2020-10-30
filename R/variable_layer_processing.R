## ############################## ##
##
## Script Name: Functions for extracting values from large variable layers
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
## Abstract: Various for extracting data from variable layers using different methods and inputs for deriving points to extract from.These include:
##  - var.prep - from a satellite image will prepare rasters of all the layers, generate indices specified, prepare additional variable layers supplied, all in the same resolution, crs, extent, and will mask out ares specified in a mask rasterlayer.  
##  - var.extract - take a training raster as the input, derives points and tterates through variable layers and extract values at points
##  - var.extract.pnts - simple extracts points from a list of variables
##  - var.extract.iterate - simple extract points from a folder to variable layers
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## fasterize_1.0.0 sf_0.7-4        dplyr_0.8.1     rgdal_1.4-4     raster_3.0-7   
## rgeos_0.4-3     sp_1.3-1 
##
## ############################## ##


#' #### variable layer prep ####
#'# function which takes a satellite image and list of variable layer locations. It first can create the indicies layer from the sat image, then iterates through all the layers to mask out the peatland boundaries layer and write all the prepared files out to the specified folder 
#' @param sat.img filepath to satellite image
#' @param indices character vector of indices to calculate, passed to the run_Indices.R function
#' @param layer.list list of addition layers to include in the analysis, expects filepaths
#' @param mask shapefile to mask the variables to 
#' @param msk.field character,field in the shapefile with which to mask with
#' @param nirband band number for the near infrared band
#' @param rband band number for the red band
#' @param bband band number for the blue band
#' @param gband band number for the green band
#' @param swirband band number for the short wave infrared band
#' @param out.folder folderpath where to save the outputs
#' @param gDAL.path path to where the gdal repository can be found
#'
#' @return
#' @export
#'
#' @examples
var.prep <- function(sat.img,indices=NA,layer.list=NA,mask,msk.field="DN",nirband=8,rband=3,bband=1,gband=2,swirband=NA,out.folder, gDAL.path){
  
  #load packages
  library(sp)
  library(rgeos)
  library(raster)
  library(rgdal)
  library(dplyr)
  library(sf)
  library(fasterize)
  #create folder for outputs
  if (!file.exists(paste(out.folder, "RFR_layers", sep = ""))) {
    dir.create(paste(out.folder, "RFR_layers", sep = ""))
  }
  ## run indices if required ##
  if(!is.na(indices)){
    source("run_Indices.R"
)
    runIndices(image.path=sat.img, out.path=paste(out.folder,"RFR_layers/",sep=""), nir=nirband,r=rband, b=bband, g=gband,swir=swirband, indices=indices, nf=F)
  }
  ## load in image
  sat.image <- raster::brick(sat.img)
  ## mask non peat from layers ## 
  masked<- sf::read_sf(mask)
  r <- raster::raster(ext=extent(sat.image[[1]]), res=res(sat.image[[1]]),crs=crs(sat.image[[1]]))
  msk <- fasterize::fasterize(masked, r, msk.field)
  #create folder for outputs
  if (!file.exists(paste0(out.folder, "RFR_layers/rfrlayers/"))){
    dir.create(paste0(out.folder, "RFR_layers/rfrlayers/"))
  }
  ## mask sentinel indiv layers ##
  for (i in 1:sat.image@data@nlayers){
    layer <- sat.image[[i]]
    sat.imagemsk <- raster::overlay(layer,msk,fun=function(x,y){ 
      x[is.na(y[])] <- NA
      return(x)
    })
    raster::writeRaster(sat.imagemsk,paste(out.folder, "RFR_layers/rfrlayers/satimgband",i,".tif", sep = ""),overwrite=T)
  }
  ## mask indices layers ##
  if(!is.na(indices)){
    list.s <- list()
    imgname<- basename(sat.img)
    for (j in 1:length(indices)){
      index <- indices[j]
      index.l <-list(index = paste(out.folder,"RFR_layers/Indices/",index,"/",index,"_",imgname,sep=""))
      names(index.l) <- index
      list.s <- append(list.s,index.l)
    }
    for (k in 1:length(list.s)){
      layer <- raster::raster(list.s[[k]])
      layer.msk <- raster::overlay(layer,msk,fun=function(x,y){ 
        x[is.na(y[])] <- NA
        return(x)
      })
      index.name <- names(list.s[k])
      raster::writeRaster(layer.msk,paste(out.folder, "RFR_layers/rfrlayers/",index.name,".tif", sep = ""))
    }
  }
  ## prep and mask other variable layer inputs ##
  if (!is.na(layer.list)){
    gdalUtils::gdal_setInstallation(search_path=gDAL.path)
    for (l in 1: length(layer.list)){
      gdalUtils::gdalwarp(srcfile = layer.list[[l]][1], dstfile = paste(out.folder,"RFR_layers/rfrlayers/",names(layer.list[l]), "_tmp.tif", sep=""),  r = layer.list[[l]][2], tr=res(sat.image), t_srs=crs(sat.image), overwrite = T, verbose = F)
      
      
      gdalUtils::gdalwarp(srcfile = layer.list[[l]][1], dstfile = paste(out.folder,"RFR_layers/rfrlayers/",names(layer.list[l]), "_tmp.tif", sep=""),  r = layer.list[[l]][2], tr=res(sat.image), t_srs=crs(sat.image), overwrite = T, verbose = T)
      #align with the raster - extent
      gdalUtils::align_rasters(unaligned = paste0(out.folder,"RFR_layers/rfrlayers/",names(layer.list[l]), "_tmp.tif"), 
                               reference = sat.img, 
                               dstfile = paste0(out.folder,"RFR_layers/rfrlayers/",names(layer.list[l]), ".tif"),
                               nThreads = 'ALL_CPUS', r = 'near', co = 'COMPRESS=LZW', overwrite = T, verbose = F)
      
      file.remove(paste0(out.folder,"RFR_layers/rfrlayers/",names(layer.list[l]), "_tmp.tif"))
      varlayer <- raster::raster(paste0(out.folder,"RFR_layers/rfrlayers/",names(layer.list[l]), ".tif"))
      origin(varlayer) <- raster::origin(sat.image)
      varcropmsk <- raster::overlay(varlayer,msk,fun=function(x,y){ 
        x[is.na(y[])] <- NA
        return(x)
      })
      var.name <- names(layer.list[l])
      raster::writeRaster(varcropmsk,paste(out.folder, "RFR_layers/rfrlayers/",var.name,".tif", sep=""),overwrite=T)

    }
  }
}



#' #### Extract values from variables using matrices ####
#'# Iterate through and extract points from folder of variable imagery
#' @param var.folder folder path to the variable layers
#' @param train.r rasterlayer of training values
#' @param out.folder output folder to write to
#'
#' @return
#' @export
#'
#' @examples
var.extract<- function(var.folder,train.r,out.folder,img.name="combotrain"){
  library(raster)
  library(dplyr)
  # load in training raster
  train <- raster::raster(train.r)
  RFR.layers <- list.files(var.folder)
  #get just the tifs
  tifs <- RFR.layers[which(grepl(pattern=".tif",RFR.layers))]
  #crop and stack
  varstack <- raster::stack(paste0(var.folder,tifs[1]))
  varstack <- raster::crop(varstack,extent(train))
  for (i in 2:length(tifs)){
    r <- raster::raster(paste0(var.folder,tifs[i]))
    r.cro <- raster::crop(r,varstack)  
    varstack <- raster::addLayer(varstack,r.cro)
  }
  origin(varstack) <- raster::origin(train)
  names(train) <- "bare"
  trainvar <- raster::addLayer(varstack,train)
  #df to write into
  out.df<- data.frame(raster::getValues(trainvar[[1]]))
  names(out.df)<- names(trainvar[[1]])
  for (j in 2:length(trainvar@layers)){
    layer.dat<- data.frame(raster::getValues(trainvar[[j]]))
    names(layer.dat)<- names(trainvar[[j]])
    out.df <- cbind(out.df,layer.dat)
  }
  train.complete<- out.df[which(complete.cases(out.df)),]
  #split out presence and absences to get an even number of the two, ensures more sampling of values>0 in the regression.
  presence <- train.complete %>% dplyr::filter(bare > 0)
  absence <- train.complete %>% dplyr::filter(bare == 0)
  if(nrow(presence)!=nrow(absence)){ #assumes more absence than presence
    absamp <- sample(nrow(absence), nrow(presence)) #randomly sample absences
    absence <- absence[absamp,]
    }
  combotrain <- rbind(presence, absence)
  if (!file.exists(paste(out.folder, "trainingvals", sep = ""))) {
    dir.create(paste(out.folder, "trainingvals", sep = ""))
  }
    write.csv(combotrain,paste0(out.folder, "trainingvals/",img.name,".csv"))
}


#' #### extract points from variables ####
#'
#' @param points  training points
#' @param varlist list of variables
#'
#' @return
#' @export
#'
#' @examples
var.extract.pnts <- function(points,varlist){
  library(raster)
  varstack <- raster::raster(varlist[[1]])
  for (i in 2:length(varlist)){
    layer <- raster::raster(varlist[[i]])
    varstack <- raster::addLayer(varstack,layer)
  }
  pointy <- rgdal::readOGR(points)
  ind.vals<- raster::extract(varstack, pointy, df=TRUE)
  point_sf <- pointy %>% sf::st_as_sf() %>% raster::as.data.frame() %>% dplyr::select(-geometry)
  all_dat <- merge(ind.vals,point_sf,by.x="ID",by.y="id")
  return(all_dat)
}



#' #### extract points from variables - iterate through folders ####
#'# Iterate through and extract points from folder of variable imagery
#creates training dataframe of presence and absences
#' @param points training points to extract values at
#' @param var.folder folder path of variable layers
#' @param out output folder to write to
#'
#' @return
#' @export
#'
#' @examples
var.extract.iterate <- function(points, var.folder,out){
  library(sp)
  library(rgeos)
  library(raster)
  library(rgdal)
  library(dplyr)
  var.list <- list.files(var.folder)
  points.df <- read.csv(points)
  sp::coordinates(points.df)<-c("x","y")
  vars <- NULL
  cpts <- data.frame(sp::coordinates(points.df))
  for (i in 1:length(var.list)){
    r.item <- raster::brick(paste(var.folder, var.list[i], sep=""))
    var.pnt <- raster::extract(r.item,points.df, na.rm=T, df=T,quick=F)
    var.points <- cbind(cpts,var.pnt[-1])
    var.comp <- var.points[complete.cases(var.points)==T,]
    names(var.comp)<-c("x","y",1:(length(names(var.comp))-2))
    vars <- rbind(vars,var.comp)
  }
  write.csv(vars,paste(out,basename(var.folder),"_vartrain.csv",sep=""))
}
