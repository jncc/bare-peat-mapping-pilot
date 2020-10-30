## ############################## ##
##
## Script Name: Image overlap
## 
## Author: JNCC
## 
## Date Created: 2019-11-14
## 
## Date Modified:2020-10-30
## 
## Licence: MIT Licence
##
## 
## Abstract:
## Script to pull out overlapping tiles between imagery and return as a list of rasters
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## dplyr_0.8.1  rgeos_0.4-3  rgdal_1.4-4 
## raster_3.0-7 sp_1.3-1 doParallel_1.0.14
##
## ############################## ##
# 

#' Image overlap
#'
#' @param P.path the path to larger image to overlap to
#' @param A.path the path to imagery to find overlaps with
#' @param img.no img.no is the img number in the list from P.path you want to use in the overlap. defaults to 1 where only 1 image present in folder
#' @param polys optional polys to restrict tiles to a mask shapefile 
#' @param parallel whether to run with parallel processing
#'
#' @return
#' @export
#'
#' @examples
overlap.img <- function(P.path,A.path,img.no=1, polys=NULL,parallel=F){
  #load packages
  library(raster)
  library(rgeos)
  library(parallel)
  #get first image
  if (file_test("-f",P.path)==F){
    p.list <- list.files(P.path)
    filename <- paste(P.path, p.list[img.no], sep="")
    granule <- raster::brick(filename)
  } else {
    granule <- raster::brick(P.path)
  }
  ukgrid = raster::crs("+init=epsg:27700")
  raster::crs(granule) <- ukgrid
  # compare to aerial, list those which overlap
  a.list <- list.files(A.path)
  aerial.overlap <- list()
  #set up parallel if needed
  if (parallel == T) { 
    cl <- parallel::makeCluster(parallel::detectCores() - 1)  
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl,expr=paste("library",c("raster", "sf", "parallel","foreach","doParallel"), (";")))
    message("Parallel clustering established.")
  }
  
  if(parallel == F){
   for (i in 1:length(a.list)){
    a.name <- paste(A.path, a.list[i], sep="")
    aerial <- raster::brick(a.name)
    raster::crs(aerial) <- ukgrid
    # find if the tile intersects then add to list
    over.a <- raster::intersect(raster::extent(aerial),raster::extent(granule))
    name <-  gsub(a.list[i],pattern=".tif",replacement="")
    if(!is.null(over.a)){
      item <- list(name=aerial)
      names(item)<- name
      aerial.overlap <- append(aerial.overlap,item)
    }
    } 
  } else {
    foreach::foreach(i = 1:length(a.list)) %dopar% {
      a.name <- paste(A.path, a.list[i], sep="")
      aerial <- raster::brick(a.name) 
      raster::crs(aerial) <- ukgrid
      over.a <- raster::intersect(raster::extent(aerial),raster::extent(granule))
      name <-  gsub(a.list[i],pattern=".tif",replacement="")
      if(!is.null(over.a)){
        item <- list(name=aerial)
        names(item)<- name
        aerial.overlap <- append(aerial.overlap,item)
      }
    }
  }
  
  if (!is.null(polys)){ # finds if it also intersects with a shapefile
    poly <- raster::shapefile(polys)
    raster::crs(poly)<-NA
    poly.overlap <- list()
    for (j in 1:length(aerial.overlap)){
      r.item <- aerial.overlap[[j]]
      rast.over<- rgeos::gIntersects(as(extent(r.item),'SpatialPolygons'),poly)
      if(isTRUE(rast.over)){
        item <- list(r.item)
        names(item) <- names(aerial.overlap[j])
        poly.overlap <- append(poly.overlap,item)
      }
    }
    list2env(list(overlap.list=poly.overlap),.GlobalEnv)
  } else{
    list2env(list(overlap.list=aerial.overlap),.GlobalEnv)
  }
}

