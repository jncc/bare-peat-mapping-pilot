## ############################## ##
##
## Script Name: Thresholding bare peat based on indices
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
## Abstract: function for thresholding based on indices. This will generate the specified indices and threshold based on a supplied function using these layers. The function must use the indices in alphabetical order in order to read in correctly.
## 
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
##  dplyr_0.8.1  rgeos_0.4-3  rgdal_1.4-4 
##  raster_3.0-7 sp_1.3-1  
##
## ############################## ##
#


#' thresholding bare peat based on indices
#'
#' @param Img.path path to the folder with indices values or list of file paths
#' @param out.path output folder to write output to
#' @param ind.name name of indices, will assume that the file naming convention is ind.name_... and will replace with BARE_ when writing out
#' @param c.fun function with which to calculate across the indices stack
#' @param nir near infrared band in the multispectral imagery
#' @param r red band in the multispectral imagery
#' @param g green band in the multispectral imagery
#' @param b blue band in the multispectral imagery
#' @param swir short wave infrared band in the multispectral imagery
#' @param start starting image to process, defaults to the first image
#' @param end  ending image to process, if NULL will default to the last image in the folder
#'
#' @return
#' @export
#'
#' @examples
barethresh <- function(Img.path,out.path,spec.bands=NA,ind.name = NA, c.fun, nir=4, r=1,g=2,b=3,swir=NA,start=1,end=NA){
  
  #load dependencies
  library(raster)
  library(rgdal)
  require(rgeos)
  require(dplyr)
  
  #create bare peat folder
  if (!file.exists(paste(out.path,"bare",sep=""))){
    dir.create(paste(out.path,"/bare",sep=""))
  }
  
  #creates a temporary file to write to
  if (!file.exists(paste0(out.path,"bare","/temp"))){
    dir.create(paste0(out.path,"bare","/temp"))
  }
  
  # list all images
  if(class(Img.path)=="character"){
    Img.list <- list.files(Img.path)
  } else {
    Img.list <-Img.path
  }
  
  # get the last image number to loop through if none supplied
  if(is.na(end)){ 
    end <- length(Img.list)
  }
  
  #load up run_Indices function if needed
  if(!is.na(ind.name[1])){
    source("./Functions/run_Indices.R")
  }
  
  #loop through images
  for (i in start:end){
    # get spec bands if included in thresholding rule
    if (!is.na(spec.bands[1])){  
      if(class(Img.path)=="character"){
        spec.1 <- raster::stack(paste0(Img.path,Img.list[i]))
      } else{
        spec.1 <- raster::stack(Img.list[[i]])
      }
      #loop through spectral bands required
      all.band <- list()
      for (j in spec.bands){
        spec.band <- list(spec.1[[j]])
        names(spec.band) <- paste0("sb", j)
        all.band <-append(all.band,spec.band)
      }
    } else{
      all.band <- list()
    }
    
    if(!is.na(ind.name[1])){
      #generates the indices with the runIndices function in run_Indices.R
      if(class(Img.path)=="character"){
        runIndices(image.path=paste0(Img.path,Img.list[i]),out.path=paste0(out.path,"bare","/temp/"), nir=4, r=1,g=2,b=3,swir=swir, indices=ind.name,nf=F)
        name <- basename(paste0(Img.path,Img.list[i]))
      } else{
        runIndices(image.path=Img.list[[i]],out.path=paste0(out.path,"bare","/temp/"), nir=4, r=1,g=2,b=3,swir=swir, indices=ind.name,nf=F)
        name <- paste0(names(Img.list[i]),".tif")
      }
      #takes first indices and lists files in folder
      Indices.list <- list.files(paste0(out.path,"bare/temp/Indices/"))
      indat.list <- list()
      for (n in 1:length(ind.name)){
        ind.1 <- raster::raster(paste0(out.path,"bare/temp/Indices/",Indices.list[n],"/",Indices.list[n],"_",name))
        item <- list(ind.1)
        names(item) <- Indices.list[n]
        indat.list <- append(indat.list, item)
      }
      
    } else {
      indat.list <- list()
    }
    
    #join the lists and stack
    all.list <- append(all.band,indat.list)
    all.stack <- raster::stack(all.list)
    #Apply function to raster stack
    r.class <- raster::overlay(all.stack, fun=c.fun)
    #save it
    raster::writeRaster(r.class,paste(out.path,"bare/BARE_",name,sep=""), overwrite=T)
    print(paste(i, "of", length(Img.list), "done"))
    #remove the temporary folder
    unlink(paste0(out.path,"bare/temp"), recursive =TRUE)
  }
  
}
