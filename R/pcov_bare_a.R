## ############################## ##
##
## Script Name: Percentage cover of high res pixels compared to lower resolution imagery
## 
## Author: JNCC
## 
## Date Created: 2019-11-15
## 
## Date Modified:  2020-10-30
## 
## Licence: MIT Licence
##
## 
## Abstract:
## # calculating the percentage cover of a low res image from resampling a higher resolution image. A matrix approach to speed up processing time as opposed to spatial calculations. data retrieval is based on cell sizes and alginment between the low resolution and high resolution images.
##outputs: will produce two csvs, one with percentage covers of all classes, another with percentage cover of just the class of interest (here bare peat) denoted by the class variable.
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## prioritizr_4.1.3 proto_1.0.0      sf_0.7-4        
## fasterize_1.0.0  tibble_2.1.3     plyr_1.8.4      
## tidyr_0.8.3      dplyr_0.8.1      raster_3.0-7    
## sp_1.3-1 
## 
##
## ############################## ##


#' Percentage cover of high res pixels compared to lower resolution imagery
#'
#' @param inClassImage high resolution imagery filepath
#' @param inPredImage low resolution imagery filepath
#' @param out.path output folder
#' @param mask optional, mask image outputs
#' @reclass whether you wish to reclassify the training data layer, i.e. if more than one class is present
#' @classcol if reclass = T, the class you wish to keep
#'
#' @return
#' @export
#'
#' @examples
pcov.resample <-  function(inClassImage, inPredImage, out.path, mask=NA,reclass=T,classcol=NA){
  #load packages
  library(raster)
  library(dplyr)
  library(tidyr)
  library(plyr)
  library(tibble)
  library(fasterize)
  library(sf)
  library(prioritizr)
  
  ## Load the classified image
  classImage <- raster::raster(inClassImage)
  predImage <- raster::raster(inPredImage, band = 1)
  raster::NAvalue(predImage) <- 0
  
  ## mask out peat soils from sat image with peat mask - if mask provided
  if (!is.na(mask)){
    mask.r <- sf::read_sf(mask)
    r <-raster::raster(ext = raster::extent(predImage),crs = raster::crs(predImage),res = raster::res(predImage))
    mask.rast <- fasterize::fasterize(mask.r, r)
    predImage.peat <-raster::overlay(predImage,mask.rast,fun = function(x, y) {
      x[is.na(y[])] <- NA
          return(x)
        }
      )
    } else {
      predImage.peat <- predImage
    }
  
  ## calculate cell size
  predImageRes <- raster::res(predImage.peat)[1]
  classImageRes <- raster::res(classImage)[1]
  halfCell <- predImageRes / 2 # distance from the center of cell to the edge=assuming a square pixel
  
  ## get common extent between images
  #line up rasters
  origin(predImage.peat) <- raster::origin(classImage)
  # get common extent
  commonExt <- raster::intersect(raster::extent(predImage.peat), raster::extent(classImage))
  # crop low res image to common extent
  sampleCellscrp <- raster::crop(predImage.peat, commonExt)
  #writeRaster(sampleCellscrp,paste(out.path,"sampleCellscrp.tif",sep=""))
  #reclassify class image
  if(reclass == T){
    minval <- classImage@data@min
    maxval <- classImage@data@max
    classmat <- data.frame(from = c(NA, seq(minval:maxval)), to = 0)
    classmat$to[which(classmat$from == classcol)] <- 1
    classed <- reclassify(classImage, rcl = classmat, progress = "window")
    raster::writeRaster(classed, paste0(out.path, "pcov_classed.tif"), overwrite =T)
    } else {
      classImage[is.na(classImage)] <- 0
      classed <- classImage
      raster::writeRaster(classed, paste0(out.path, "pcov_classed.tif"), overwrite =T)
    }
  
  ## create output matrix
  out <- NULL
  # testing - get cell number from Sentinel image
  # x <- cellFromXY(sampleCellscrp,xy=matrix(c(402755,402935), ncol=6))
  # ext.x <- extentFromCells(sampleCellscrp,x)
  # class.crp <- crop(classed,ext.x)
  # z <-rowColFromCell(sampleCellscrp,x)
  ##
  
  # get difference in resolution
  resdif <- predImageRes / classImageRes #40
  #number of rows in low res image
  nrowno <-(sampleCellscrp@extent@ymax - sampleCellscrp@extent@ymin) / predImageRes #1900
  #number of cols in low res image
  ncolno <-(sampleCellscrp@extent@xmax - sampleCellscrp@extent@xmin) / predImageRes #2600
  outputMatrix <- matrix(nrow = nrowno, ncol = ncolno)
  
  ## iterate through to get the values from high res image
  for (i in 1:nrowno) {
    for (j in 1:ncolno) {
      # retreive values in high res block
      vals <-raster::getValuesBlock(
          classed,
          row = 1 + (resdif * (i - 1)),
          col = (resdif * (j - 1)) + 1,
          ncol = resdif,
          nrow = resdif
        )
      #get a percentage cover
      pcnt <- sum(vals) / (resdif ^ 2)
      #write back to matrix
      outputMatrix[i, j] <- pcnt
      cat("row ", i, ", col ", j, " done. \r", sep = "")
      flush.console()
    }
  }
  
  ## write output as raster
  r <-raster::raster(
      outputMatrix,
      xmn = raster::xmin(sampleCellscrp),
      xmx = raster::xmax(sampleCellscrp),
      ymn = raster::ymin(sampleCellscrp),
      ymx = raster::ymax(sampleCellscrp)
    )
  raster::writeRaster(r, paste(out.path, "pcov_training.tif", sep = ""), overwrite =
                        T)
  write.csv(outputMatrix, paste(out.path, "pcov_training.csv", sep = ""))
}
