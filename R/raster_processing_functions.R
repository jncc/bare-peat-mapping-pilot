## ############################## ##
##
## Script Name: RasterLayer processing functions
##
## Author: JNCC
##
## Contact (if different from above):
##
## Date Created: 2019-11-13
##
## Date Modified: 2020-10-30
##
## Licence: MIT Licence
##
##
## Abstract: Various functions for processing rasters:
##  - r.process - crops a raster to a shape, performs reprojecting the raster and reclassifys
##  - shp.buff - buffer a shapefile and output as a raster or shapefile
##  - combo.img - Function designed to combine RGB imagery with additional NIR band from infrared image
##  - gdb.tile - mosaics data from a geodatabase of tiles, this is based on providing a shapefile of the area of interest you wish to retreive data for and uses the geodatabase index to retreive just those tiles
##  - topo.gen - creates hillshade, aspect and slope layers using gdalUtils::gdaldem
##  - mosaic.img - function to mosaic either all the images in a given folder or a list of rasters.
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## prioritizr_4.1.3   proto_1.0.0
## gdalUtils_2.0.1.14 stringi_1.4.3
## dplyr_0.8.1        fasterize_1.0.0
## APfun_0.1.3        sf_0.7-4
## rgeos_0.4-3        raster_3.0-7
## rgdal_1.4-4        sp_1.3-1
##
## ############################## ##


#' #### raster processing####
# will crop to a shapefile if a shape is provided, assumes same projection.
# reclassify to a given matrix if one is provided#'
#' @param layer The Rasterlayer you wish to process
#' @param AOI Optional - a shapefile to crop the extent to
#' @param reproj Optional - crs you wish to reproject to
#' @param reclass Optional - matrix to reclassify the raster with, see ?reclassify for details.
#' @param poly logical, if you wish to write out as a polygon
#' @param out.path folder to save the output to
#' @param out.name file name to save the output as
#' @param OSGeoPath Path of OSGeo, required for polygonizing
#'
#' @return
#' @export
#'
#' @examples
#'
r.process <- function(layer,AOI=NA,reproj=NA ,reclass=NA,poly=F,out.path,out.name=NA, OSGeoPath){
  #get packages
  library(rgdal)
  library(raster)
  library(rgeos)
  library(sf)
  library(sp)
  library(APfun)
  library(fasterize)
  #load in raster
  layer.r <- raster::brick(layer)
  if(layer.r@data@nlayers==1){
    layer.r <- raster::raster(layer)
  }
  # if AOI provided clip the raster to extent
  if(!is.na(AOI)){
    AOIcrp <- raster::shapefile(AOI)
   layer_out <- raster::crop(layer.r,AOIcrp)
   cat("cropped..        \r")
  }
  # if reproj is supplied, then reproject
  if(!is.na(reproj)){
    if(class(reproj)!="CRS"){stop("reproj needs to be CRS.")}
    layer_out <- raster::projectRaster(layer_out, crs=reproj)
    cat("reprojected..    \r")
  }
  # if reclass matrix provided reclassify the raster
  if(class(reclass)!="logical"){
    layer_out <- raster::reclassify(layer_out, reclass)
    cat("reclassed..    \r")
  }
  #write out files as raster
  if(is.na(out.name)){
    out.name <- gsub(basename(layer),pattern=".tif",replacement="")
  }
  raster::writeRaster(layer_out, paste(out.path,out.name,"_clp.tif",sep=""))
  # polygonise if needed
  if(poly==T){
    APfun::APpolygonize(paste(out.path,out.name,"_clp.tif",sep=""), readToMemory = TRUE, outFile = paste(out.path,out.name,"_clp.shp",sep=""),OSGeoPath = OSGeoPath)
  }
}


#' #### buffer shapefile layer ####
#' function will buffer then shapefile and return as raster if required
#' @param shape shapefile
#' @param dist distance to buffer around shape
#' @param temp.r raster template to rasterize the output
#' @param field field in the shapefile to rasterize on
#' @param shp logical, if true write out as a shapefile
#' @param tif  logical, if true write out as a raster
#' @param out.path output folder to write to
#' @param out.name output file name
#'
#' @return
#' @export
#'
#' @examples
shp.buff <- function(shape,dist=0,temp.r,field,shp=T,tif=F,out.path,out.name=NA){
  #get packages
  library(sf)
  library(fasterize)
  library(dplyr)
  #load in shape
  shape.sf <- sf::read_sf(shape)
  #buffer
  shape.buff <- shape.sf  %>% sf::st_buffer(dist = 10) %>% as('Spatial')
  if(is.na(out.name)){
    out.name <- gsub(basename(shape),pattern=".shp",replacement="")
  }
  if(shp==T){ #write out as shapefile
    raster::shapefile(shape.buff,paste(out.path,out.name,".shp",sep=""))
  }
  if(tif==T){
    shape.buffsf <- sf::st_as_sf(shape.buff)
    peat.buff.r <- fasterize::fasterize(shape.buffsf, temp.r, field)
    raster::writeRaster(peat.buff.r,paste(out.path,out.name,".tif",sep=""))
  }
}



#' #### combo.img ####
#' Function designed to combine RGB imagery with additional NIR band from infrared image. It resamples image b to image a resolution then combines the two
# adds just the first layer from the b.path image but can be modified to add all or select layers
#' @param a.path path to RGB image folder
#' @param b.path path to infrared image folder
#' @param nir band in b image which contains the Near-infrared band to join
#' @param out.path output folder to save to
#' @param start image to start on, will use list.files of the .tif files in folder
#' @param end image to end on, will use list.files of the .tif files in folder
#'
#' @return
#' @export
#'
#' @examples
combo.img <- function(a.path, b.path,nir=1,out.path, start=1,end=NULL){
  library(raster)
  library(rgeos)
  library(rgdal)
  #get all files
  aerial.list <- list.files(a.path)
  aerial.listed <- aerial.list[which(grepl(pattern=".tif",aerial.list))]
  infra.list <- list.files(b.path)
  infra.listed <- infra.list[which(grepl(pattern=".tif",infra.list))]
  #makes ouput folder
  if (!file.exists(paste(out.path,"combined",sep=""))){
    dir.create(paste(out.path,"combined",sep=""))
  }
  #iterates through resample to 0.5 res and combine with infra-red band
  if(is.null(end)){
    end <- length(aerial.listed)
  }
  for (i in start:end){
    granule <- raster::stack(paste(a.path, aerial.listed[i], sep=""))
    inf.gran <- raster::stack(paste(b.path, infra.listed[i], sep=""))
    #resample granule to same resolution
    gran.samp <- raster::resample(inf.gran,granule,  method = "bilinear")
    gran.all <- raster::addLayer(granule, gran.samp[[nir]])
    name <- aerial.listed[i]
    raster::writeRaster(gran.all, paste(out.path,"combined/",name,sep=""), overwrite=T)
    print(paste(i,"of",length(aerial.list),"done."))
  }

}



#' ##### geodatabase mosaicing ####
#' Function to retrieve geodatabase tiles from a given AOI and then mosaic them together outputing a raster
#' @param AOI Area of Interest to retrieve tiles
#' @param gdbindex.path path to the geodatabase index as a shapefile of tile references. This will assume the tiles are in the same folder path
#' @param out.path output folder to write to
#'
#' @return
#' @export
#'
#' @examples
gdb.tile <- function(AOI, gdbindex.path, out.path){
  library(raster)
  library(sp)
  library(rgeos)
  library(rgdal)
  library(sf)
  library(dplyr)
  library(stringi)
  library(stringr)

  #load and crop grid ref to aoi
  aoi <- sf::st_read(AOI)
  grid <- sf::st_read(gdbindex.path)
  aoi.grid <- st_intersection(grid,aoi)
  files <- aoi.grid %>% sf::st_as_sf() %>% raster::as.data.frame() %>%
    dplyr::select(LETTERS,FULL_PATH) %>% dplyr::mutate_all(as.character)
  #retrieve last part of file path
  for (i in 1:nrow(files)){
    cut.off<- stringr::str_locate(pattern=files$LETTERS[i],files$FULL_PATH[i])
    files$name[i]<- substr(files$FULL_PATH[i], cut.off[1,1], nchar(files$FULL_PATH[i]))
  }
  #get full path names
  files.load <- files %>% dplyr::mutate(fullname = paste(dirname(gdbindex.path),"/",name,sep=""))
  files.load$fullname <- gsub("\\\\","/",files.load$fullname)
  #loop through and grab tiles
  tiles <- list()
  for (i in 1:nrow(files.load)){
    r <- raster::raster(files.load$fullname[i])
    tiles <- append(tiles, list(r))
  }
  #mosaic
  tiles$fun <- mean
  rast.mosaic <- do.call(mosaic,tiles)
  if (!file.exists(paste(out.path,"DTM",sep=""))){
    dir.create(paste(out.path,"/DTM",sep=""))
  }
  raster::writeRaster(rast.mosaic,paste(out.path,"DTM/","DTM_mosaic.tif",sep=""))
}


#' #### topology layers ####
#'
#' @param dem The full path to the dem file.
#' @param out.path folder to store outputs in
#' @param h logical, TRUE to run hillslope
#' @param a logical, TRUE to run aspect
#' @param s logical, TRUE to run slope
#' @param GDAL.path path to link to GDAL
#'
#' @return
#' @export
#'
#' @examples
topo.gen <- function(dem,out.path,h=T,a=T,s=T,GDAL.path){
  require(gdalUtils)
  gdalUtils::gdal_setInstallation(GDAL.path)
  if(h==T){
    gdalUtils::gdaldem(mode="hillshade",input_dem=dem,output = paste(out.path,"DTM_hillshade.tif",sep=""),verbose=T)
    }
  if(a==T){
    gdalUtils::gdaldem(mode="aspect",input_dem=dem,output = paste(out.path,"DTM_aspect.tif",sep=""),verbose=T)
    }
  if(s==T){
    gdalUtils::gdaldem(mode="slope",input_dem=dem,output = paste(out.path,"DTM_slope.tif",sep=""),verbose=T)
  }
}



#' #### mosaic image ####
#'#mosaicing function. takes either a file path or image list of rasters
#' @param combo.path path to folder of rasters to mosaic together
#' @param img.list list of rasterlayers or rasterstacks you wish to mosaic together
#' @param out.path folder to save the output
#' @param name name of the output raster
#'
#' @return
#' @export
#'
#' @examples
mosaic.img <- function(combo.path=NULL,img.list=NULL,out.path, name="mosiac.tif"){
  library(raster)
  #check either have combo.path or img.list
  if(is.null(combo.path)&is.null(img.list)){
    stop("need to supply either combo.path OR img.list")
  } else if(!is.null(combo.path)&!is.null(img.list)){
    stop("only supply either combo.path OR img.list")
  }
  if(!is.null(combo.path)){
    combo.list <- list.files(combo.path)
    for (i in 1:length(combo.list)){
      filename <- paste(combo.path, combo.list[i], sep="")
      granule <- raster::raster(filename)
      img.list <- append(img.list,list(granule))
    }
  }
  names(img.list) <- NULL
  img.list$fun <- mean
  rast.mosaic <- do.call(raster::mosaic,c(img.list,progress="window"))
  raster::writeRaster(rast.mosaic,paste(out.path,name,sep=""))
}
