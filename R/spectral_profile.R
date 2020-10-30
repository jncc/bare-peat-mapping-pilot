## ############################## ##
##
## Script Name: Spectral profiling
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
## Abstract:
## script to assess the spectral profile for bare and vegetated peat using digitised polygons as training data 
##
##
## R version 3.6.0 (2019-04-26)
## Dependencies:
## ggplot2_3.1.1   fasterize_1.0.0
## tidyr_0.8.3     dplyr_0.8.1    
## rgdal_1.4-4     sf_0.7-4       
## rgeos_0.4-3     raster_3.0-7   
## sp_1.3-1  
##
## ############################## ##

#' Spectral profiling
#'
#' @param train.path path to the training polygons
#' @param sat.path path to the imagery to extract band data from
#' @param out.path output folder to store results
#' @param img.no image number in the imagery list to use, default will take first image
#' @param nsamps number of sample points to assess, default will take minimum number from classes present
#'
#' @return
#' @export
#'
#' @examples
spec.profile <- function(train.path, sat.path, out.path, img.no=1,nsamps=NULL){
  
  # load package dependencies
  library(sp)
  library(sf)
  library(rgeos)
  library(plyr)
  library(rgdal)
  library(dplyr)
  library(tidyr)
  library(raster)
  library(fasterize)
  library(ggplot2)
  
  #read in data
  training <- rgdal::readOGR(train.path)
  granule <- raster::brick(sat.path)
  #if polygons supplied then randomly sample points
  if(class(training)=="SpatialPolygonsDataFrame"){
    train.sf <- sf::read_sf(train.path)
    #convert to raster the to points to get most possible points
    r <- raster::raster(ext=extent(granule), res=res(granule))
    train.r <- fasterize::fasterize(train.sf, r, "bare")
    peat.out <- raster::overlay(train.r,granule[[1]],fun=function(x,y){ 
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
    ## randomly sample points and extract values ##
    #vegetated
    class_data<- peat.points[peat.points$layer==classes$class[1],]
    veg_class <-  class_data[sample(1:length(class_data), nsamps), ]
    #bare
    bare_data<- peat.points[peat.points$layer==classes$class[2],]
    bare_class <-  bare_data[sample(1:length(bare_data), nsamps), ]
    allpoints <- rbind(veg_class,bare_class)
    } else { # otherwise if points use points
      allpoints <- training
    }
  
  #extracting the values from the points
  satValues <- raster::extract(granule, allpoints, df=TRUE)
  if(class(training)=="SpatialPolygonsDataFrame"){
    classVal <- raster::extract(peat.out, allpoints, df=TRUE)
  } else {
    allpoints$id <- 1:nrow(allpoints)
    classVal <- data.frame(allpoints$id,allpoints$layer)
  }
  Sat.data <- cbind(satValues[,2:5],classVal[,2])
  colnames(Sat.data) <- c("band_1", "band_2", "band_3", "band_4", "class")
  Sat.data$class <- as.factor(Sat.data$class)
 
  
  #save training points
  if (!file.exists(paste(out.path,"Spectral_analysis",sep=""))){
    dir.create(paste(out.path,"Spectral_analysis",sep=""))
  }
  name <- gsub(basename(sat.path),pattern="_osgb_vmsk_rad_sref",replacement="")
  name <- gsub(name,pattern=".tif",replacement="")
  
  if(class(training)=="SpatialPolygonsDataFrame"){
    raster::shapefile(allpoints, paste0(out.path,"Spectral_analysis/train_",name,".shp"),overwrite=T)
  }
  
  #create facet plot
  satlongdat <- Sat.data %>% dplyr::mutate(bare = ifelse(class==0, "Vegetated peat", "Bare peat")) %>% tidyr::gather(Bands,Value,1:4)
  
  p <- ggplot(data = satlongdat, aes(x=Bands, y=Value)) + 
    geom_boxplot(aes(fill=bare)) +
    ggplot2::scale_fill_manual(values=c("Bare peat" = "#8B4513","Vegetated peat" = "#228B22")) 
  p + facet_wrap( ~ Bands, scales="free")
  ggplot2::ggsave(paste0(out.path,"Spectral_analysis/facet_",name,".png"))
  
  ##plot spectral profile
  # summarise dataframe into min mean max per group
  satdat<- Sat.data %>% dplyr::group_by(class) %>% 
    dplyr::summarise(b1=mean(band_1), b2=mean(band_2),b3=mean(band_3), b4=mean(band_4),b1min=min(band_1),b2min=min(band_2),b3min=min(band_3),b4min=min(band_4),b1max=max(band_1),b2max=max(band_2),b3max=max(band_3),b4max=max(band_4))
  g1 <- satdat %>% tidyr::gather(key="band",value="val",2:5) %>% dplyr::select(1,10,11)
  g2 <- satdat %>% tidyr::gather(key="band",value="min",6:9) %>% dplyr::select(1,10,11)
  g3 <- satdat %>% tidyr::gather(key="band",value="max",10:13) %>% dplyr::select(1,10,11)
  data.gath <- cbind(g1,min=g2$min,max=g3$max)
  data.gath$class <- plyr::revalue(data.gath$class,c("0"="vegetated", "1"="bare"))
  
  #plot 1 - with ribbons
  ggplot2::ggplot(data.gath, aes(x=band, y=val,group=class, colour=class,fill=class)) +
    ggplot2::geom_ribbon(aes(ymin=min,ymax=max),alpha=.3,colour=NA,show.legend = T)+ 
    ggplot2::geom_line(aes(y=val)) +
    ggplot2::ylim(0,max(data.gath$max)) +
    ggplot2::scale_fill_manual(values=c(vegetated = "#228B22",bare = "#8B4513")) +
    ggplot2::scale_colour_manual(values=c(vegetated = "#228B22",bare = "#8B4513")) + 
    ggplot2::labs(group=NULL,colour="",x="Band",y="Reflectance") 
  #write out results
  write.csv(Sat.data, paste0(out.path,"Spectral_analysis/resp_",name,".csv"))
  ggplot2::ggsave(paste0(out.path,"Spectral_analysis/resp_",name,".png"))
  
 #plot 2
  ggplot2::ggplot(data.gath, aes(x=band, y=val,group=class, colour=class,fill=class)) +
    ggplot2::geom_line(aes(y=val)) +
    ggplot2::ylim(0,max(data.gath$max)) +
    ggplot2::scale_fill_manual(values=c(vegetated = "#228B22",bare = "#8B4513")) +
    ggplot2::scale_colour_manual(values=c(vegetated = "#228B22",bare = "#8B4513")) + 
    ggplot2::labs(group=NULL,colour="",x="Band",y="Reflectance") 
  ggplot2::ggsave(paste0(out.path,"Spectral_analysis/resp2_",name,".png"))
 #return data to environment
 list2env(list(sat.dat=Sat.data),.GlobalEnv)
  }

