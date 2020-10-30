## ############################## ##
##
## Script Name: script to calculate vegetation indices from multispectral imagery
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
## Abstract: This script can calculate the several indices, taking as an input either the path to the multispectral image or a folder containing many images. The indices it accepts are:
##            - EVI - Enhanced Vegetation Index
##            - GLI - Green Leaf Index
##             - GNDVI - Green normalized difference index
##             - RDVI - Renormalized Difference Vegetation Index
##             - SAVI - Soil Adjusted Vegetation Index
##             - SBL - Soil Background line
##             - RB - Red/blue
##             - RG - Red/Green
##             - Brightness - Brightness (mean of RGB bands)
##             - NBR - Normalized Difference NIR/SWIR Normalized burn ratio
##
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## dplyr_0.8.1  raster_3.0-7
## sp_1.3-1  
##
## ############################## ##

#' Run Indices
#'
#' @param image.path path to folder with multispectral imagery
#' @param out.path path to save NDVI output files
#' @param nir near infra-red band
#' @param r red band
#' @param b blue band
#' @param g green band
#' @param swir optional - Short wave infrared band
#' @param indices a vector of indices to retrieve e.g. c("EVI", "GLI","GNDVI","RDVI","SAVI","SBL"). The function will accept EVI,GLI,GNDVI,RDVI,SAVI,SBL,RB,RG,Brightness,NBR
#' @param nf logical, if you wish to create new folder per image
#'
#' @return
#' @export
#'
#' @examples
#' 
runIndices <- function(image.path, out.path, nir,r, b, g, swir=NA,indices, nf=T){
  require(raster)
  require(dplyr)
  #new folder and checks
  if (!file.exists(paste(out.path,"Indices",sep=""))){
    dir.create(paste(out.path,"Indices",sep=""))
  }
  if(is.na(nir)|is.na(r)|is.na(b)|is.na(g)){
    stop("need to supply min band numbers.")
  }
  #test if file or folder
  if (file_test("-f",image.path)==T){
    # get file and break up into raster bands
    filename <- image.path
    granule <- raster::brick(filename)
    nir.band <- granule[[nir]]
    r.band <- granule[[r]]
    g.band <- granule[[g]]
    b.band <- granule[[b]]
    if(!is.na(swir)){
      swir.band <- granule[[swir]]
    }
    name <- basename(filename) %>% gsub(pattern=".tif", replacement="")
    if (nf == T){
      if (!file.exists(paste(out.path,"Indices/",name,sep=""))){
        dir.create(paste(out.path,"Indices/",name,sep=""))
      }
    }
    ## calculate indices
    #EVI - Enhanced vegetation index
    if("EVI" %in% indices){
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/EVI",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/EVI",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/EVI",sep=""))){
          dir.create(paste(out.path,"Indices/EVI",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){ # if file is small run calculation directly
        EVI <- 2.5*((nir.band- r.band)/(nir.band+(6* r.band)-(7.5*b.band))+1)
        if (nf == T){
          raster::writeRaster(EVI, paste0(out.path, "Indices/",name,"/EVI/EVI_",name,".tif"), overwrite=T)
        } else { 
          raster::writeRaster(EVI, paste0(out.path, "Indices/EVI/EVI_",name,".tif"), overwrite=T)
        }
      } else { # if large run calculation on blocks
        bs <- raster::blockSize(granule, n=1) #get block size
        #EVI
        if(nf==T){
          named <- paste0(out.path, "Indices/",name,"/EVI/EVI_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/EVI/EVI_",name,".tif")
        }
        EVI.out <- raster::raster(granule)
        EVI.out <- raster::writeStart(EVI.out, named) # set up writing to new raster
        for (n in 1:bs$n) { #iterate through the blocks, run calc and writing out to those blocks on the new raster
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
          EVI <- 2.5*((nir.val-r.val)/(nir.val+(6*r.val)-(7.5*b.val))+1)
          EVI.out <- raster::writeValues(EVI.out, EVI, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        EVI.out <- raster::writeStop(EVI.out)
      }
      print(paste("EVI done"))
    }
    
    if("GLI" %in% indices){
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/GLI",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/GLI",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/GLI",sep=""))){
          dir.create(paste(out.path,"Indices/GLI",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        GLI <- ((2*g.band)-r.band-b.band)/((2*g.band)+r.band+b.band)
        if (nf==T){
          raster::writeRaster(GLI, paste0(out.path, "Indices/",name,"/GLI/GLI_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(GLI, paste0(out.path, "Indices/GLI/GLI_",name,".tif"), overwrite=T)
        }
      } else{
        bs <- raster::blockSize(granule, n=1)
        if (nf==T){
          named <- paste0(out.path, "Indices/",name,"/GLI/GLI_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/GLI/GLI_",name,".tif")
        }
        GLI.out <- raster::raster(granule)
        GLI.out <- raster::writeStart(GLI.out, named)
        for (n in 1:bs$n) {
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
          b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
          GLI <- ((2*g.val)-r.val-b.val)/((2*g.val)+r.val+b.val)
          GLI.out <- raster::writeValues(GLI.out, GLI, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        GLI.out <- raster::writeStop(GLI.out)
      }      
      print(paste("GLI done"))
    }
    
    if("GNDVI" %in% indices){
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/GNDVI",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/GNDVI",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/GNDVI",sep=""))){
          dir.create(paste(out.path,"Indices/GNDVI",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        GNDVI <- (nir.band-g.band)/(nir.band+g.band)
        if (nf==T) {        
          raster::writeRaster(GNDVI, paste0(out.path, "Indices/",name,"/GNDVI/GNDVI_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(GNDVI, paste0(out.path, "Indices/GNDVI/GNDVI_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/GNDVI/GNDVI_",name,".tif")
        } else{
          named <- paste0(out.path, "Indices/GNDVI/GNDVI_",name,".tif")
        }
        
        GNDVI.out <- raster::raster(granule)
        GNDVI.out <- raster::writeStart(GNDVI.out, named)
        for (n in 1:bs$n) {
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
          GNDVI <- (nir.val-g.val)/(nir.val+g.val)
          GNDVI.out <- raster::writeValues(GNDVI.out, GNDVI, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        GNDVI.out <- raster::writeStop(GNDVI.out)
      }
      print(paste("GNDVI done"))
    }
    
    if("RDVI" %in% indices){
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/RDVI",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/RDVI",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/RDVI",sep=""))){
          dir.create(paste(out.path,"Indices/RDVI",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        RDVI <- (nir.band-r.band)/((nir.band+r.band)^0.5)
        if (nf == T) {
          raster::writeRaster(RDVI, paste0(out.path, "Indices/",name,"/RDVI/RDVI_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(RDVI, paste0(out.path, "Indices/RDVI/RDVI_",name,".tif"), overwrite=T)
        }       
        
      } else {
        bs <- raster::blockSize(granule, n=1)
        if  (nf == T){
          named <- paste0(out.path, "Indices/",name,"/RDVI/RDVI_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/RDVI/RDVI_",name,".tif")
        }
        
        RDVI.out <- raster::raster(granule)
        RDVI.out <- raster::writeStart(RDVI.out, named)
        for (n in 1:bs$n) {
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          RDVI <- (nir.val-r.val)/((nir.val+r.val)^0.5)
          RDVI.out <- raster::writeValues(RDVI.out, RDVI, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        RDVI.out <- raster::writeStop(RDVI.out)
      }
      print(paste("RDVI done"))
    }
    
    if("SAVI" %in% indices){
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/SAVI",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/SAVI",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/SAVI",sep=""))){
          dir.create(paste(out.path,"Indices/SAVI",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        SAVI <- ((nir.band-r.band)/(nir.band+r.band+0.5))*(1+0.5)
        if (nf == T){
          raster::writeRaster(SAVI, paste0(out.path, "Indices/",name,"/SAVI/SAVI_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(SAVI, paste0(out.path, "Indices/SAVI/SAVI_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/SAVI/SAVI_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/SAVI/SAVI_",name,".tif")
        }
        SAVI.out <- raster::raster(granule)
        SAVI.out <- raster::writeStart(SAVI.out, named)
        for (n in 1:bs$n) {
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          SAVI <- ((nir.val-r.val)/(nir.val+r.val+0.5))*(1+0.5)
          SAVI.out <- raster::writeValues(SAVI.out, SAVI, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        SAVI.out <- raster::writeStop(SAVI.out)
      }
      print(paste("SAVI done"))
    }
    
    if("SBL" %in% indices){
      if (nf ==T){
        if (!file.exists(paste(out.path,"Indices/",name,"/SBL",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/SBL",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/SBL",sep=""))){
          dir.create(paste(out.path,"Indices/SBL",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        SBL <- (nir.band-(2.4*r.band))
        if (nf == T){
          raster::writeRaster(SBL, paste0(out.path, "Indices/",name,"/SBL/SBL_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(SBL, paste0(out.path, "Indices/SBL/SBL_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/SBL/SBL_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/SBL/SBL_",name,".tif")
        }
        SBL.out <- raster::raster(granule)
        SBL.out <- raster::writeStart(SBL.out, named)
        for (n in 1:bs$n) {
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          SBL <- (nir.val-(2.4*r.val))
          SBL.out <- raster::writeValues(SBL.out, SBL, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        SBL.out <- raster::writeStop(SBL.out)
      }
      print(paste("SBL done"))
    }
    
    if("RB" %in% indices){
      if ( nf==T){
        if (!file.exists(paste(out.path,"Indices/",name,"/RB",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/RB",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/RB",sep=""))){
          dir.create(paste(out.path,"Indices/RB",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        RB <- (r.band/b.band)
        if (nf ==T){
          raster::writeRaster(RB, paste0(out.path, "Indices/",name,"/RB/RB_",name,".tif"), overwrite=T)
        } else{
          raster::writeRaster(RB, paste0(out.path, "Indices/RB/RB_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/RB/RB_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/RB/RB_",name,".tif")
        }
        RB.out <- raster::raster(granule)
        RB.out <- raster::writeStart(RB.out, named)
        for (n in 1:bs$n) {
          b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          RB <- (r.val/b.val)
          RB.out <- raster::writeValues(RB.out, RB, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        RB.out <- raster::writeStop(RB.out)
      }
      print(paste("RB done"))
    }
    
    if("RG" %in% indices){
      if (nf ==T){
        if (!file.exists(paste(out.path,"Indices/",name,"/RG",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/RG",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/RG",sep=""))){
          dir.create(paste(out.path,"Indices/RG",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        RG <- (r.band/g.band)
        if (nf==T){
          raster::writeRaster(RG, paste0(out.path, "Indices/",name,"/RG/RG_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(RG, paste0(out.path, "Indices/RG/RG_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf ==T){
          named <- paste0(out.path, "Indices/",name,"/RG/RG_",name,".tif") 
        } else {
          named <- paste0(out.path, "Indices/RG/RG_",name,".tif") 
        }
        RG.out <- raster::raster(granule)
        RG.out <- raster::writeStart(RG.out, named)
        for (n in 1:bs$n) {
          g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          RG <- (r.val/g.val)
          RG.out <- raster::writeValues(RG.out, RG, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        RG.out <- raster::writeStop(RG.out)
      }
      print(paste("RG done"))
    }
    
    if("Brightness" %in% indices) {
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/Brightness",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/Brightness",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/Brightness",sep=""))){
          dir.create(paste(out.path,"Indices/Brightness",sep=""))
        }
      }
      
      if (file.size(filename)<1000000000){
        bright <- ((r.band+g.band+b.band)/3)
        if (nf == T){
          raster::writeRaster(bright, paste0(out.path, "Indices/",name,"/Brightness/Brightness_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(bright, paste0(out.path, "Indices/Brightness/Brightness_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/Brightness/Brightness_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/Brightness/Brightness_",name,".tif")
        }
        bright.out <- raster::raster(granule)
        bright.out <- raster::writeStart(bright.out, named)
        for (n in 1:bs$n) {
          g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
          bright <- ((r.val+g.val+b.val)/3)
          bright.out <- raster::writeValues(bright.out, bright, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        bright.out <- raster::writeStop(bright.out)
      }
      print(paste("bright done"))
    }
    
    if("NDVI" %in% indices) {
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/NDVI",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/NDVI",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/NDVI",sep=""))){
          dir.create(paste(out.path,"Indices/NDVI",sep=""))
        }
      }
      if (file.size(filename)<1000000000){
        NDVI <- (nir.band-r.band)/(nir.band+r.band)
        if (nf == T){
          raster::writeRaster(NDVI, paste0(out.path, "Indices/",name,"/NDVI/NDVI_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(NDVI, paste0(out.path, "Indices/NDVI/NDVI_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/NDVI/NDVI_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/NDVI/NDVI_",name,".tif")
        }
        NDVI.out <- raster::raster(granule)
        NDVI.out <- raster::writeStart(NDVI.out, named)
        for (n in 1:bs$n) {
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
          NDVI <- (nir.val-r.val)/(nir.val+r.val)
          NDVI.out <- raster::writeValues(NDVI.out, NDVI, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        NDVI.out <- raster::writeStop(NDVI.out)
      }
      print(paste("NDVI done"))
    }
    # NBR
    if("NBR" %in% indices) {
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,"/NBR",sep=""))){
          dir.create(paste(out.path,"Indices/",name,"/NBR",sep=""))
        }
      } else {
        if (!file.exists(paste(out.path,"Indices/NBR",sep=""))){
          dir.create(paste(out.path,"Indices/NBR",sep=""))
        }
      }
      if (file.size(filename)<1000000000){
        NBR <- (nir.band-swir.band)/(nir.band+swir.band)
        if (nf == T){
          raster::writeRaster(NBR, paste0(out.path, "Indices/",name,"/NBR/NBR_",name,".tif"), overwrite=T)
        } else {
          raster::writeRaster(NBR, paste0(out.path, "Indices/NBR/NBR_",name,".tif"), overwrite=T)
        }
      } else {
        bs <- raster::blockSize(granule, n=1)
        if (nf == T){
          named <- paste0(out.path, "Indices/",name,"/NBR/NBR_",name,".tif")
        } else {
          named <- paste0(out.path, "Indices/NBR/NBR_",name,".tif")
        }
        NBR.out <- raster::raster(granule)
        NBR.out <- raster::writeStart(NBR.out, named)
        for (n in 1:bs$n) {
          swir.val <- raster::getValues(swir.band, row=bs$row[n], nrows=bs$nrows[n] )
          nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
          NBR <- (nir.val-swir.val)/(nir.val+swir.val)
          NBR.out <- raster::writeValues(NBR.out, NBR, bs$row[n])
          print(paste("Processed", n, "in", bs$n))
        }
        NBR.out <- raster::writeStop(NBR.out)
      }
      print(paste("NBR done"))
    }
  } else {
    imagery.list <- list.files(image.path) # if folder supplied get list of images to iterate through
    #iterate through files in folder
    for (i in 1:length(imagery.list)){  
      filename <- paste(image.path, imagery.list[i], sep="")
      granule <- raster::brick(filename)
      nir.band <- granule[[nir]]
      r.band <- granule[[r]]
      g.band <- granule[[g]]
      b.band <- granule[[b]]
      if(!is.na(swir)){
        swir.band <- granule[[swir]]
      }
      name <- basename(filename) %>% gsub(pattern=".tif", replacement="")
      if (nf == T){
        if (!file.exists(paste(out.path,"Indices/",name,sep=""))){
          dir.create(paste(out.path,"Indices/",name,sep=""))
        }
      }
      #EVI - Enhanced vegetation index
      if("EVI" %in% indices){
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/EVI",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/EVI",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/EVI",sep=""))){
            dir.create(paste(out.path,"Indices/EVI",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          EVI <- 2.5*((nir.band- r.band)/(nir.band+(6* r.band)-(7.5*b.band))+1)
          if (nf == T){
            raster::writeRaster(EVI, paste(out.path, "Indices/",name,"/EVI/EVI_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(EVI, paste(out.path, "Indices/EVI/EVI_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          #EVI
          if(nf==T){
            named <- paste(out.path, "Indices/",name,"/EVI/EVI_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/EVI/EVI_",imagery.list[i],sep="")
          }
          EVI.out <- raster::raster(granule)
          EVI.out <- raster::writeStart(EVI.out, named)
          for (n in 1:bs$n) {
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
            EVI <- 2.5*((nir.val-r.val)/(nir.val+(6*r.val)-(7.5*b.val))+1)
            EVI.out <- raster::writeValues(EVI.out, EVI, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          EVI.out <- raster::writeStop(EVI.out)
        }
        print(paste("EVI done"))
      }
      
      if("GLI" %in% indices){
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/GLI",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/GLI",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/GLI",sep=""))){
            dir.create(paste(out.path,"Indices/GLI",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          GLI <- ((2*g.band)-r.band-b.band)/((2*g.band)+r.band+b.band)
          if (nf==T){
            raster::writeRaster(GLI, paste(out.path, "Indices/",name,"/GLI/GLI_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(GLI, paste(out.path, "Indices/GLI/GLI_",imagery.list[i],sep=""), overwrite=T)
          }
        } else{
          bs <- raster::blockSize(granule, n=1)
          if (nf==T){
            named <- paste(out.path, "Indices/",name,"/GLI/GLI_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/GLI/GLI_",imagery.list[i],sep="")
          }
          GLI.out <- raster::raster(granule)
          GLI.out <- raster::writeStart(GLI.out, named)
          for (n in 1:bs$n) {
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
            b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
            GLI <- ((2*g.val)-r.val-b.val)/((2*g.val)+r.val+b.val)
            GLI.out <- raster::writeValues(GLI.out, GLI, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          GLI.out <- raster::writeStop(GLI.out)
        }      
        print(paste("GLI done"))
      }
      
      if("GNDVI" %in% indices){
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/GNDVI",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/GNDVI",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/GNDVI",sep=""))){
            dir.create(paste(out.path,"Indices/GNDVI",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          GNDVI <- (nir.band-g.band)/(nir.band+g.band)
          if (nf==T) {        
            raster::writeRaster(GNDVI, paste(out.path, "Indices/",name,"/GNDVI/GNDVI_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(GNDVI, paste(out.path, "Indices/GNDVI/GNDVI_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/GNDVI/GNDVI_",imagery.list[i],sep="")
          } else{
            named <- paste(out.path, "Indices/GNDVI/GNDVI_",imagery.list[i],sep="")
          }
          
          GNDVI.out <- raster::raster(granule)
          GNDVI.out <- raster::writeStart(GNDVI.out, named)
          for (n in 1:bs$n) {
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
            GNDVI <- (nir.val-g.val)/(nir.val+g.val)
            GNDVI.out <- raster::writeValues(GNDVI.out, GNDVI, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          GNDVI.out <- raster::writeStop(GNDVI.out)
        }
        print(paste("GNDVI done"))
      }
      
      if("RDVI" %in% indices){
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/RDVI",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/RDVI",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/RDVI",sep=""))){
            dir.create(paste(out.path,"Indices/RDVI",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          RDVI <- (nir.band-r.band)/((nir.band+r.band)^0.5)
          if (nf == T) {
            raster::writeRaster(RDVI, paste(out.path, "Indices/",name,"/RDVI/RDVI_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(RDVI, paste(out.path, "Indices/RDVI/RDVI_",imagery.list[i],sep=""), overwrite=T)
          }       
          
        } else {
          bs <- raster::blockSize(granule, n=1)
          if  (nf == T){
            named <- paste(out.path, "Indices/",name,"/RDVI/RDVI_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/RDVI/RDVI_",imagery.list[i],sep="")
          }
          
          RDVI.out <- raster::raster(granule)
          RDVI.out <- raster::writeStart(RDVI.out, named)
          for (n in 1:bs$n) {
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            RDVI <- (nir.val-r.val)/((nir.val+r.val)^0.5)
            RDVI.out <- raster::writeValues(RDVI.out, RDVI, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          RDVI.out <- raster::writeStop(RDVI.out)
        }
        print(paste("RDVI done"))
      }
      
      if("SAVI" %in% indices){
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/SAVI",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/SAVI",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/SAVI",sep=""))){
            dir.create(paste(out.path,"Indices/SAVI",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          SAVI <- ((nir.band-r.band)/(nir.band+r.band+0.5))*(1+0.5)
          if (nf == T){
            raster::writeRaster(SAVI, paste(out.path, "Indices/",name,"/SAVI/SAVI_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(SAVI, paste(out.path, "Indices/SAVI/SAVI_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/SAVI/SAVI_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/SAVI/SAVI_",imagery.list[i],sep="")
          }
          SAVI.out <- raster::raster(granule)
          SAVI.out <- raster::writeStart(SAVI.out, named)
          for (n in 1:bs$n) {
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            SAVI <- ((nir.val-r.val)/(nir.val+r.val+0.5))*(1+0.5)
            SAVI.out <- raster::writeValues(SAVI.out, SAVI, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          SAVI.out <- raster::writeStop(SAVI.out)
        }
        print(paste("SAVI done"))
      }
      
      if("SBL" %in% indices){
        if (nf ==T){
          if (!file.exists(paste(out.path,"Indices/",name,"/SBL",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/SBL",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/SBL",sep=""))){
            dir.create(paste(out.path,"Indices/SBL",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          SBL <- (nir.band-(2.4*r.band))
          if (nf == T){
            raster::writeRaster(SBL, paste(out.path, "Indices/",name,"/SBL/SBL_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(SBL, paste(out.path, "Indices/SBL/SBL_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/SBL/SBL_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/SBL/SBL_",imagery.list[i],sep="")
          }
          SBL.out <- raster::raster(granule)
          SBL.out <- raster::writeStart(SBL.out, named)
          for (n in 1:bs$n) {
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            SBL <- (nir.val-(2.4*r.val))
            SBL.out <- raster::writeValues(SBL.out, SBL, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          SBL.out <- raster::writeStop(SBL.out)
        }
        print(paste("SBL done"))
      }
      
      if("RB" %in% indices){
        if ( nf==T){
          if (!file.exists(paste(out.path,"Indices/",name,"/RB",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/RB",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/RB",sep=""))){
            dir.create(paste(out.path,"Indices/RB",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          RB <- (r.band/b.band)
          if (nf ==T){
            raster::writeRaster(RB, paste(out.path, "Indices/",name,"/RB/RB_",imagery.list[i],sep=""), overwrite=T)
          } else{
            raster::writeRaster(RB, paste(out.path, "Indices/RB/RB_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/RB/RB_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/RB/RB_",imagery.list[i],sep="")
          }
          RB.out <- raster::raster(granule)
          RB.out <- raster::writeStart(RB.out, named)
          for (n in 1:bs$n) {
            b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            RB <- (r.val/b.val)
            RB.out <- raster::writeValues(RB.out, RB, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          RB.out <- raster::writeStop(RB.out)
        }
        print(paste("RB done"))
      }
      
      if("RG" %in% indices){
        if (nf ==T){
          if (!file.exists(paste(out.path,"Indices/",name,"/RG",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/RG",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/RG",sep=""))){
            dir.create(paste(out.path,"Indices/RG",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          RG <- (r.band/g.band)
          if (nf==T){
            raster::writeRaster(RG, paste(out.path, "Indices/",name,"/RG/RG_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(RG, paste(out.path, "Indices/RG/RG_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf ==T){
            named <- paste(out.path, "Indices/",name,"/RG/RG_",imagery.list[i],sep="") 
          } else {
            named <- paste(out.path, "Indices/RG/RG_",imagery.list[i],sep="") 
          }
          RG.out <- raster::raster(granule)
          RG.out <- raster::writeStart(RG.out, named)
          for (n in 1:bs$n) {
            g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            RG <- (r.val/g.val)
            RG.out <- raster::writeValues(RG.out, RG, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          RG.out <- raster::writeStop(RG.out)
        }
        print(paste("RG done"))
      }
      
      if("Brightness" %in% indices) {
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/Brightness",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/Brightness",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/Brightness",sep=""))){
            dir.create(paste(out.path,"Indices/Brightness",sep=""))
          }
        }
        
        if (file.size(filename)<1000000000){
          bright <- ((r.band+g.band+b.band)/3)
          if (nf == T){
            raster::writeRaster(bright, paste(out.path, "Indices/",name,"/Brightness/Brightness_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(bright, paste(out.path, "Indices/Brightness/Brightness_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/Brightness/Brightness_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/Brightness/Brightness_",imagery.list[i],sep="")
          }
          bright.out <- raster::raster(granule)
          bright.out <- raster::writeStart(bright.out, named)
          for (n in 1:bs$n) {
            g.val <- raster::getValues(g.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            b.val <- raster::getValues(b.band, row=bs$row[n], nrows=bs$nrows[n] )
            bright <- ((r.val+g.val+b.val)/3)
            bright.out <- raster::writeValues(bright.out, bright, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          bright.out <- raster::writeStop(bright.out)
        }
        print(paste("bright done"))
      }
      
      if("NDVI" %in% indices) {
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/NDVI",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/NDVI",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/NDVI",sep=""))){
            dir.create(paste(out.path,"Indices/NDVI",sep=""))
          }
        }
        if (file.size(filename)<1000000000){
          NDVI <- (nir.band-r.band)/(nir.band+r.band)
          if (nf == T){
            raster::writeRaster(NDVI, paste(out.path, "Indices/",name,"/NDVI/NDVI_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(NDVI, paste(out.path, "Indices/NDVI/NDVI_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/NDVI/NDVI_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/NDVI/NDVI_",imagery.list[i],sep="")
          }
          NDVI.out <- raster::raster(granule)
          NDVI.out <- raster::writeStart(NDVI.out, named)
          for (n in 1:bs$n) {
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            r.val <- raster::getValues(r.band, row=bs$row[n], nrows=bs$nrows[n] )
            NDVI <- (nir.val-r.val)/(nir.val+r.val)
            NDVI.out <- raster::writeValues(NDVI.out, NDVI, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          NDVI.out <- raster::writeStop(NDVI.out)
        }
        print(paste("NDVI done"))
      }
      # NBR
      if("NBR" %in% indices) {
        if (nf == T){
          if (!file.exists(paste(out.path,"Indices/",name,"/NBR",sep=""))){
            dir.create(paste(out.path,"Indices/",name,"/NBR",sep=""))
          }
        } else {
          if (!file.exists(paste(out.path,"Indices/NBR",sep=""))){
            dir.create(paste(out.path,"Indices/NBR",sep=""))
          }
        }
        if (file.size(filename)<1000000000){
          NBR <- (nir.band-swir.band)/(nir.band+swir.band)
          if (nf == T){
            raster::writeRaster(NBR, paste(out.path, "Indices/",name,"/NBR/NBR_",imagery.list[i],sep=""), overwrite=T)
          } else {
            raster::writeRaster(NBR, paste(out.path, "Indices/NBR/NBR_",imagery.list[i],sep=""), overwrite=T)
          }
        } else {
          bs <- raster::blockSize(granule, n=1)
          if (nf == T){
            named <- paste(out.path, "Indices/",name,"/NBR/NBR_",imagery.list[i],sep="")
          } else {
            named <- paste(out.path, "Indices/NBR/NBR_",imagery.list[i],sep="")
          }
          NBR.out <- raster::raster(granule)
          NBR.out <- raster::writeStart(NBR.out, named)
          for (n in 1:bs$n) {
            swir.val <- raster::getValues(swir.band, row=bs$row[n], nrows=bs$nrows[n] )
            nir.val <- raster::getValues(nir.band, row=bs$row[n], nrows=bs$nrows[n] )
            NBR <- (nir.val-swir.val)/(nir.val+swir.val)
            NBR.out <- raster::writeValues(NBR.out, NBR, bs$row[n])
            print(paste("Processed", n, "in", bs$n))
          }
          NBR.out <- raster::writeStop(NBR.out)
        }
        print(paste("NBR done"))
      }
      #finished image indices
      print(paste(i, "of", length(imagery.list), "done"))
    }
  }
}
