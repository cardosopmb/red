#####RED - IUCN Redlisting Tools
#####Version 1.2.1 (2017-09-27)
#####By Pedro Cardoso
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P.(in prep.) An R package to facilitate species red list assessments according to the IUCN criteria.
#####Changed from v1.2.0:
#####Added vignette
#####Centroid at outlier graph
#####on.exit at red.setup() to avoid unexpected behavior if function fails

#Todo:
# rli*: add trees (for phylogenetic and functional diversity)
# map.habitat: automatically select habitat or range
# map.change: new function and automated download of future climate layers and forest change across years
# map.sdm: multicore p multiple runs
# kbas: implement

#####data origins:
#####climate -> Fick, S.E. & Hijmans, R.J. (2017) Worldclim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, in press.
#####altitude -> Farr, T. G., et al. (2007), The Shuttle Radar Topography Mission, Rev. Geophys., 45, RG2004
#####landcover -> Tuanmu, M.-N. & Jetz, W. (2014) A global 1-km consensus land-cover product for biodiversity and ecosystem modeling. Global Ecology and Biogeography, 23: 1031-1045.

#####RED Stats:
#####library("cranlogs")
#####day <- cran_downloads(package = "red", from = "2016-08-19", to = "2017-06-21")
#####group <- matrix(day$count, 100, byrow = TRUE)
#####plot(rowSums(group), type = "n")
#####lines(rowSums(group))

#####required packages
library("BAT")
library("dismo")
library("geosphere")
library("graphics")
library("grDevices")
library("jsonlite")
library("maptools")
library("raster")
library("rgdal")
library("rgeos")
library("sp")
library("stats")
library("utils")
#' @import graphics
#' @import jsonlite
#' @import maptools
#' @import rgdal
#' @import rgeos
#' @import sp
#' @import stats
#' @import utils
#' @importFrom geosphere areaPolygon
#' @importFrom grDevices chull dev.copy dev.off pdf
#' @importFrom raster area cellStats clump crop extent extract getValues layerStats mask raster rasterize rasterToPoints rasterToPolygons reclassify res sampleRandom scalebar terrain trim writeRaster xmax xmin

###############################################################################
##############################AUX FUNCTIONS####################################
###############################################################################

raster::rasterOptions(maxmemory = 2e+09)

longlat2utm <- function(longlat){
  longlat = as.matrix(longlat)
  minlong = min(longlat[,1])
  zone = floor((minlong + 180) / 6) + 1
  res = rgdal::project(longlat, paste("+proj=utm +zone=",zone," ellps=WGS84",sep=''))
  return(res)
}

utm2longlat <- function(utm, zone){
  if(class(utm) == "RasterLayer"){
    if(!is.null(zone))
      raster::crs(utm) <- paste("+proj=utm +zone=", zone, sep="")
    res <- raster::projectRaster(utm, crs = "+proj=longlat +datum=WGS84", method='ngb')
  } else {
    utm <- SpatialPoints(utm, CRS(paste("+proj=utm +zone=", zone,sep="")))
    res <- as.data.frame(spTransform(utm,CRS(paste("+proj=longlat"))))
  }
  return(res)
}

##warn if maxent.jar is not available
warnMaxent <- function(){
  warning("RED could not find maxent.jar.
  1. Download the latest version of maxent from:
  https://biodiversityinformatics.amnh.org/open_source/maxent/
  2. Move the file maxent.jar to the java directory inside dismo package
  (there should be a file named dismo.jar already there)
  3. Install the latest version of java runtime environment (JRE) with the same architecture (32 or 64 bits) as your version of R:
  http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html")
}

##detect which layers are categorical by checking if all values are integers and if the max is less than 50 (may fail, just an attempt)
find.categorical <- function(layers){
  categorical = c()
  for(l in 1:(dim(layers)[3])){
    lay <- raster::as.matrix(layers[[l]])
    lay[is.na(lay)] <- 0
    if(sum(floor(lay)) == sum(lay) && max(lay) < 50)
      categorical = c(categorical, l)
  }
  return(categorical)
}

##basic function to calculate the rli of any group of species
rli.calc <- function(spData, tree = NULL, boot = FALSE, runs = 1000){
  if(max(spData) > 5){                                ##if letters are given, convert to [0,1]
    spData <- replace(spData, which(spData == "EX" ), 0)
    spData <- replace(spData, which(spData == "EW" ), 0)
    spData <- replace(spData, which(spData == "RE" ), 0)
    spData <- replace(spData, which(spData == "CR" ), 0.2)
    spData <- replace(spData, which(spData == "EN" ), 0.4)
    spData <- replace(spData, which(spData == "VU" ), 0.6)
    spData <- replace(spData, which(spData == "NT" ), 0.8)
    spData <- replace(spData, which(spData == "LC" ), 1)
    spData <- replace(spData, which(spData == "DD" ), 2)
    spData <- as.numeric(spData)
    spData <- subset(spData, spData < 2)
  } else if (max(spData) > 1){                       ##if a scale [0,5] is given, convert to [0,1]
    spData <- 1 - spData / 5
  }
  if(is.null(tree)){                           ##if not weighted by PD or FD
    if(!boot){                                 ##if no bootstrap to be made
      return (mean(spData))
    } else {
      run <- rep(NA, runs)
      for(i in 1:runs){
        rnd <- sample(spData, length(spData), replace = TRUE) ##bootstrap
        run[i] <- mean(rnd)
      }
      res <- matrix(quantile(run, c(0.025, 0.5, 0.975)), nrow = 1)
      colnames(res) <- c("LowCL", "Median", "UpCL")
      return(res)
    }
  } else {                                     ##if weighted by PD or FD
    comm <- matrix(1, nrow = 2, ncol = length(spData))
    contrib <- BAT::contribution(comm, tree, relative = TRUE)[1,]
    if(!boot){                                 ##if no bootstrap to be made
      return (sum(contrib * spData))
    } else {
      for(i in 1:runs){
        rnd <- sample(spData, length(spData), replace = TRUE) ##bootstrap
        for(i in 1:length(rnd))
          run[i] <- run[i] + contrib[rnd[i]] * spData[rnd[i]]
      }
      return(quantile(run, c(0.025, 0.5, 0.975)))
    }
  }
}

##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################

#' Setup GIS directory.
#' @description Setup directory where GIS files are stored.
#' @details Writes a txt file in the red directory allowing the package to always access the world GIS files directory.
#' @export
red.setDir <- function(){
  redFile <- paste(find.package("red"), "/red.txt", sep = "")
  dir <- readline("Input directory for storing world gis layers:")
  dir <- paste(dir, "/", sep = "")
  dput(dir, redFile)
}

#' Read GIS directory.
#' @description Read directory where GIS files are stored.
#' @details Reads a txt file pointing to where the world GIS files are stored.
#' @export
red.getDir <- function(){
  redFile <- paste(find.package("red"), "/red.txt", sep = "")
  if (file.exists(redFile)){       #if there is already a file read from it
    dir <- dget(redFile)
  } else {
    warning(paste(redFile, "not found, please run red.setDir()"))
    return()
  }
  return(dir)
}

#' Download and setup GIS files.
#' @description Setup red to work with species distribution modelling and layers available online.
#' @details Please check that you have at least 50Gb free in your disk (and a fast internet connection) to download all files. In the end of the process "only" 17.4Gb will be left though. This function will:
#' 1. Check if maxent.jar is available in the dismo package directory.
#' 2. Ask user input for GIS directory.
#' 3. Download global bioclim and elevation files (20) from http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip.
#' 4. Download landcover files (12) from http://data.earthenv.org/consensus_landcover/without_DISCover/.
#' 5. Unzip all files and delete the originals.
#' 6. Create a new layer (1) with the dominant land cover at each cell.
#' 7. Resample all files (33) to approximately 10x10km (for use with widespread species) grid cells.
#' Sit back and enjoy, this should take a while.
#' @export
red.setup <- function(){

  ##test if maxent.jar is in the right directory
  if(!file.exists(paste(.libPaths()[[1]], "/dismo/java/maxent.jar", sep=""))){
    warnMaxent()
    return()
  }

  oldwd = getwd()
  on.exit(expr = setwd(oldwd))
  gisdir = red.setDir()
  setwd(gisdir)

  ##basic setup
  pb <- txtProgressBar(min = 0, max = 33, style = 3)

  ##download and process bioclim
  download.file("http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip", "bioclim2.zip")
  unzip(zipfile = "bioclim.zip")
  file.remove("bioclim.zip")
  for(i in 1:19){
    setTxtProgressBar(pb, i)
    if(i < 10)
      rast <- raster(paste("wc2.0_bio_30s_0", i, ".tif", sep=""))
    else
      rast <- raster(paste("wc2.0_bio_30s_", i, ".tif", sep=""))
    rast <- crop(rast, c(-180, 180, -56, 90))
    writeRaster(rast, paste("red_1km_", i, ".tif", sep=""))
    rast <- aggregate(rast, 10)
    writeRaster(rast, paste("red_10km_", i, ".tif", sep=""))
    if(i < 10)
      file.remove(paste("wc2.0_bio_30s_0", i, ".tif", sep=""))
    else
      file.remove(paste("wc2.0_bio_30s_", i, ".tif", sep=""))
    gc()
  }

  ##download and process altitude
  setTxtProgressBar(pb, 20)
  download.file("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/alt_30s_bil.zip", "alt_30s_bil.zip")
  unzip(zipfile = "alt_30s_bil.zip")
  file.remove("alt_30s_bil.zip")
  rast <- raster("alt.bil")
  rast <- crop(rast, c(-180, 180, -56, 90))
  writeRaster(rast, "red_1km_20.tif")
  rast <- aggregate(rast, 10)
  writeRaster(rast, "red_10km_20.tif")
  file.remove("alt.bil")
  file.remove("alt.hdr")
  gc()

  ##download and process land cover
  altmask1 = raster("red_1km_20.tif")
  altmask10 =  raster("red_10km_20.tif")
  for(i in 5:12){
    setTxtProgressBar(pb, (i+20))
    download.file(paste("http://data.earthenv.org/consensus_landcover/without_DISCover/Consensus_reduced_class_", i, ".tif", sep=""), destfile = paste("Consensus_reduced_class_", i, ".tif", sep=""), mode = "wb")
    rast <- raster(paste("Consensus_reduced_class_", i, ".tif", sep=""))
    rast <- mask(rast, altmask1)
    writeRaster(rast, paste("red_1km_", (i+20), ".tif", sep=""))
    rast <- aggregate(rast, 10)
    #maskLayer <- sum(altmask, rast)
    #maskLayer[!is.na(maskLayer)] <- 1
    rast <- mask(rast, altmask10)
    writeRaster(rast, paste("red_10km_", (i+20), ".tif", sep=""))
    file.remove(paste("Consensus_reduced_class_", i, ".tif", sep=""))
    gc()
  }
  remove(rast)

  ##create new rasters with most common landcover at each cell
  setTxtProgressBar(pb, 33)
  max1 <- raster()
  max10 <- raster()
  for(i in 21:32){
    rast <- raster(paste("red_1km_", i, ".tif", sep=""))
    max1 <- raster::stack(max1, rast)
    rast <- raster(paste("red_10km_", i, ".tif", sep=""))
    max10 <- raster::stack(max10, rast)
  }
  max1 <- which.max(max1)
  writeRaster(max1, "red_1km_33.tif")
  max10 <- which.max(max10)
  writeRaster(max10, "red_10km_33.tif")
  remove(max1, max10)
  gc()
  setwd(oldwd)


  ##Now the files should be named as:
  ##red_1km_1.tif
  ##...
  ##red_10km_33.tif
  ##Where 1 to 19 are the corresponding bioclim variables, 20 is altitude, 21 to 32 are landcover proportion and 33 is most common landcover per cell

  #download country borders (not working Feb. 2017)
  #download.file("http://biogeo.ucdavis.edu/data/gadm2.6/countries_gadm26.rds", destfile = paste("worldcountries.rds"), mode = "wb")
}

#' Download taxon records from GBIF.
#' @description Downloads species or higher taxon data from GBIF and outputs non-duplicate records with geographical coordinates.
#' @param taxon Taxon name.
#' @details As always when using data from multiple sources the user should be careful and check if records "make sense". This can be done by either ploting them in a map (e.g. using red::map.draw()) or using red::outliers().
#' @return A data.frame with longitude and latitude, plus species names if taxon is above species.
#' @examples records("Nephila senegalensis")
#' @export
records <- function(taxon){
  taxon = unlist(strsplit(taxon, split = " ")[[1]])
  dat <- dismo::gbif(taxon[1], paste(taxon[2], "*", sep = ""))
  dat <- dat[c("species","lon","lat")] #filter columns
  dat <- dat[!(is.na(dat$lon) | is.na(dat$lat)),] #filter rows
  dat <- unique(dat)       #delete duplicate rows
  colnames(dat) <- c("Species", "long", "lat")
  if (length(taxon) == 1){      #if genus
    dat[which(is.na(dat[,1])),1] <- paste(taxon, "sp.")
  } else {                #if species
    dat <- dat[,-1]
  }
  return(dat)
}

#' Move records to closest non-NA cell.
#' @description Identifies and moves presence records to cells with environmental values.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster.
#' @param buffer Maximum distance in map units that a record will move. If 0 all NA records will be changed.
#' @details Often records are in coastal or other areas for which no environmental data is available. This function moves such records to the closest cells with data so that no information is lost during modelling.
#' @return A matrix with new coordinate values.
#' @examples rast <- raster::raster(matrix(c(rep(NA,100), rep(1,100), rep(NA,100)), ncol = 15))
#' pts <- cbind(runif(100, 0, 0.55), runif(100, 0, 1))
#' raster::plot(rast)
#' points(pts)
#' pts <- move(pts, rast)
#' raster::plot(rast)
#' points(pts)
#' @export
move <- function(longlat, layers, buffer = 0){
  layers <- layers[[1]]
  values <- extract(layers, longlat)   #get values of each record
  suppressWarnings(
    for(i in which(is.na(values))){    #if a value is NA, move it
      distRaster = raster::distanceFromPoints(layers, longlat[i,])
      distRaster = mask(distRaster, layers)
      vmin = raster::minValue(distRaster)
      if(buffer <= 0 || buffer > vmin){
        vmin = rasterToPoints(distRaster, function(x) x == vmin)
        longlat[i,] = vmin[1,1:2]
      }
    }
  )
  return(longlat)
}

#' Visual detection of outliers.
#' @description Draws plots of sites in geographical (longlat) and environmental (2-axis PCA) space.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster. It can be any set of environmental layers thought to allow the identification of environmental outliers.
#' @details Erroneous data sources or errors in transcriptions may introduce outliers that can be easily detected by looking at simple graphs of geographical or environmental space.
#' @return A data.frame with coordinate values and distance to centroid in pca is returned. Two plots are drawn for visual inspection. The environmental plot includes row numbers for easy identification of possible outliers.
#' @examples data(red.records)
#' data(red.layers)
#' outliers(red.records, red.layers[[1:3]])
#' @export
outliers <- function(longlat, layers){
  if(dim(layers)[3] == 33)      #if layers come from raster.read
    pca <- raster.reduce(layers[[1:19]], n = 2)
  else
    pca <- raster.reduce(layers, n = 2)
  ##extract pca values from longlat
  pca <- as.data.frame(raster::extract(pca, longlat))
  goodRows <-  which(!is.na(pca[,1]))
  pca <- pca[goodRows,]
  longlat <- longlat[goodRows,]
  par(mfrow = c(1,2))
  map.draw(longlat, layers[[1]], spName = "Geographical")
  raster::plot(pca, main = "Environmental", type = "n")
  centroid = colMeans(pca)
  text(centroid[1], centroid[2], label = "X")
  for(i in 1:nrow(pca)){
    text(pca[i,1], pca[i,2], label = row.names(longlat)[i])
  }

  ##build new matrix ordered by distance to centroid
  dist2centroid = apply(pca, 1, function(x) dist(rbind(x, centroid)))
  out = as.data.frame(cbind(longlat, dist2centroid))
  out = out[order(-dist2centroid),]
  return(out)
}

#' Spatial thinning of occurrence records.
#' @description Thinning of records with minimum distances either absolute or relative to the species range.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param distance Distance either in relative terms (proportion of maximum distance between any two records) or in raster units.
#' @param relative If TRUE, represents the proportion of maximum distance between any two records. If FALSE, is in raster units.
#' @param runs Number of runs
#' @details Clumped distribution records due to ease of accessibility of sites, emphasis of sampling on certain areas in the past, etc. may bias species distribution models.
#' The algorithm used here eliminates records closer than a given distance to any other record. The choice of records to eliminate is random, so a number of runs are made and the one keeping more of the original records is chosen.
#' @return A matrix of species occurrence records separated by at least the given distance.
#' @examples records <- matrix(sample(100), ncol = 2)
#' par(mfrow=c(1,2))
#' graphics::plot(records)
#' records <- thin(records, 0.1)
#' graphics::plot(records)
#' @export
thin <- function(longlat, distance = 0.01, relative = TRUE, runs = 100){
  longlat = longlat[!duplicated(longlat),]                #first, remove duplicate rows
  nSites = nrow(longlat)
  if(nSites < 4)
    return(longlat)

  ##if relative, calculate maxDist between any two points
  if(relative){
    maxDist = 0
    for(x in 1:(nSites-1)){
      for(y in (x+1):nSites){
        maxDist = max(maxDist,((longlat[x,1]-longlat[y,1])^2+(longlat[x,2]-longlat[y,2])^2)^.5)
      }
    }
    distance = maxDist*distance
  }

  listSites = matrix(longlat[1,], ncol=2, byrow = TRUE)
  for (r in 1:runs){
    longlat = longlat[sample(nSites),]       ##shuffle rows (sites)
    rndSites = longlat[1,]                   ##start with first random site
    for(newSite in 2:nSites){
      for(oldSite in 1:(newSite-1)){
        addSite = TRUE
        dist = ((longlat[newSite,1]-longlat[oldSite,1])^2+(longlat[newSite,2]-longlat[oldSite,2])^2)^.5
        if(dist < distance){
          addSite = FALSE
          break
        }
      }
      if(addSite)
        rndSites = rbind(rndSites, longlat[newSite,])
    }
    if(nrow(rndSites) > nrow(listSites))
      listSites = rndSites
  }
  return(as.matrix(listSites))
}

#' Read and buffer raster layers.
#' @description Read raster layers of environmental or other variables and crop them to a given extent around the known occurrences.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster.
#' @param ext Either extent of map or buffer around the known records used to crop layers. If buffer, it is relative to the maximum distance between any two records.
#' @details If layers are not given, the function will read either 30 arc-second (approx. 1km) or 5 arc-minutes (approx. 10km) resolution rasters from worldclim (Fick & Hijmans 2017) and landcover (Tuanmu & Jetz 2014) if red.setup() is run previously.
#' @return A RasterStack object (If no layers are given: Variables 1-19 = bioclim, 20 = elevation, 21-32 = proportion landcover, 33 = most common landcover).
#' @references Fick, S.E. & Hijmans, R.J. (2017) Worldclim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, in press.
#' @references Tuanmu, M.-N. & Jetz, W. (2014) A global 1-km consensus land-cover product for biodiversity and ecosystem modeling. Global Ecology and Biogeography, 23: 1031-1045.
#' @examples data(red.layers)
#' data(red.records)
#' par(mfrow=c(1,2))
#' raster::plot(red.layers[[1]])
#' points(red.records)
#' croppedLayers <- raster.read(red.records, red.layers, 0.1)
#' raster::plot(croppedLayers[[1]])
#' points(red.records)
#' @export
raster.read <- function(longlat, layers = NULL, ext = 1){

  xmin = min(longlat[,1])
  xmax = max(longlat[,1])
  xlen = xmax - xmin
  ymin = min(longlat[,2])
  ymax = max(longlat[,2])
  ylen = ymax - ymin

  if(is.null(layers)){          ##if no layers are provided read the ones available
    gisdir = red.getDir()

    ##calculate species range and buffer around it
    if(eoo(longlat) < 200000){
      layers <- raster::stack(raster::raster(paste(gisdir, "red_1km_1.tif", sep = "")))
      for(i in 2:33)
        layers <- raster::stack(layers, raster::raster(paste(gisdir, "red_1km_", i, ".tif", sep = "")))
    } else {
      layers <- raster::stack(raster::raster(paste(gisdir, "red_10km_1.tif", sep = "")))
      for(i in 2:33)
        layers <- raster::stack(layers, raster::raster(paste(gisdir, "red_10km_", i, ".tif", sep = "")))
    }
    ##determine longitude limits of species to check if crop and paste are needed around longitude 180 for Pacific species
    if(xmin < -90 && xmax > 90 && sum(longlat[longlat[,1] < 90 && longlat[,1] > -90,]) != 0){
      ##crop and merge layers
      rightHalf = crop(layers, c(0,180,raster::extent(layers)@ymin,raster::extent(layers)@ymax))
      raster::extent(rightHalf) <- c(-180,0,raster::extent(layers)@ymin,raster::extent(layers)@ymax)
      leftHalf = crop(layers, c(-180,0,raster::extent(layers)@ymin,raster::extent(layers)@ymax))
      raster::extent(leftHalf) <- c(0,180,raster::extent(layers)@ymin,raster::extent(layers)@ymax)
      layers <- merge(rightHalf, leftHalf)
      ##modify longlat
      for(i in 1:nrow(longlat))
        if(longlat[i,1] > 0)
          longlat[i,1] = longlat[i,1] - 180
      else
        longlat[i,1] = longlat[i,1] + 180
    }
  }

  if(length(ext) == 4)                          ##if absolute extent is given crop and return, else calculate buffer
    return(crop(layers, ext))

  if(xlen == 0)      ##in case some dimensions are inexistent consider equal to extent
    xlen = ext
  if(ylen == 0)
    ylen = ext

  ##calculate new extent of layers and crop
  ext = max(1, ((xlen + ylen) * ext))
  xmin <- max(raster::extent(layers)@xmin, xmin-ext)
  xmax <- min(raster::extent(layers)@xmax, xmax+ext)
  ymin <- max(raster::extent(layers)@ymin, ymin-ext)
  ymax <- min(raster::extent(layers)@ymax, ymax+ext)
  layers <- crop(layers, c(xmin,xmax,ymin,ymax))
  return(layers)
}

#' Uniformize raster layers.
#' @description Crop raster layers to minimum size possible and uniformize NA values across layers.
#' @param layers Raster* object as defined by package raster.
#' @details Excludes all marginal rows and columns with only NA values and change values to NA if they are NA in any of the layers.
#' @return A Raster* object, same class as layers.
#' @examples data(red.layers)
#' raster::plot(raster.clean(red.layers))
#' @export
raster.clean <- function(layers){

  ##apply mask to have NAs everywhere where any layer has NAs
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  layers <- mask(layers, maskLayer)

  ##crop by excluding external rows and columns with NAs only
  layers <- trim(layers)

  return(layers)
}

#' Reduce dimensionality of raster layers.
#' @description Reduce the number of layers by either performing a PCA on them or by eliminating highly correlated ones.
#' @param layers Raster* object as defined by package raster.
#' @param method Either Principal Components Analysis ("pca", default) or Pearson's correlation ("cor").
#' @param n Number of layers to reduce to.
#' @param thres Value for pairwise Pearson's correlation above which one of the layers (randomly selected) is eliminated.
#' @details Using a large number of explanatory variables in models with few records may lead to overfitting. This function allows to avoid it as much as possible.
#' If both n and thres are given, n has priority. If method is not recognized and layers come from raster.read function, only landcover is reduced by using only the dominating landuse of each cell.
#' @return A RasterStack object.
#' @export
raster.reduce <- function(layers, method = "pca", n = NULL, thres = NULL){
  ##method = "pca, cor", if unrecognized method only reduce landcover but not climate
  out <- raster::stack()

  if(dim(layers)[3] == 33){          ##check if layers are obtained with raster.read
    out <- raster::stack(layers[[33]])
    layers = layers[[1:19]]
  }
  if(method == "cor"){                       ##if correlation
    if(is.null(n)){
      if(is.null(thres))
        thres = 0.7
      for(i in 1:dim(layers)[3]){                  ##delete layers until none are correlated above threshold
        cor = as.matrix(as.dist(layerStats(layers, 'pearson', na.rm = TRUE)[[1]]))
        if(max(cor) < thres)
          break
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    } else {
      while (dim(layers)[3] > n){                   ##delete layers until reaching n layers
        cor = abs(as.matrix(as.dist(layerStats(layers, 'pearson', na.rm = TRUE)[[1]])))
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    }
  } else if(method == "pca"){                                  ##if pca
    if(is.null(n))
      n = 3
    if(sum(!is.na(getValues(layers[[1]]))) > 2000)
      sr <- sampleRandom(layers, 1000)
    else
      sr <- sampleRandom(layers, as.integer(sum(!is.na(getValues(layers[[1]])))/2))
    pca <- prcomp(sr)
    layers <- raster::predict(layers, pca, index = 1:n)
    for(i in 1:n)
      names(layers[[i]]) <- paste("pca",i)
  }
  out <- raster::stack(layers, out)
  return(out)
}

#' Create longitude layer.
#' @description Create a layer depicting longitude based on any other.
#' @param layers Raster* object as defined by package raster.
#' @details Using longitude (and latitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data(red.layers)
#' raster::plot(raster.long(red.layers))
#' @export
raster.long <- function(layers){
  if(dim(layers)[3] > 1)
    layers <- layers[[3]]
  x <- rasterToPoints(layers)[,1:2]
  long <- rasterize(x, layers, x[,1])
  long <- mask(long, layers)
  names(long) <- "longitude"
  return(long)
}

#' Create latitude layer.
#' @description Create a layer depicting latitude based on any other.
#' @param layers Raster* object as defined by package raster.
#' @details Using latitude (and longitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data(red.layers)
#' raster::plot(raster.lat(red.layers[[1]]))
#' @export
raster.lat <- function(layers){
  if(dim(layers)[3] > 1)
    layers <- layers[[3]]
  x <- rasterToPoints(layers)[,1:2]
  lat <- rasterize(x, layers, x[,2])
  lat <- mask(lat, layers)
  names(lat) <- "latitude"
  return(lat)
}

#' Create eastness layer.
#' @description Create a layer depicting eastness based on an elevation layer.
#' @param dem RasterLayer object of elevation (a digital elevation model - DEM) as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples data(red.layers)
#' raster::plot(raster.east(red.layers[[3]]))
#' @export
raster.east <- function(dem){
  asp <- terrain(dem, opt = "aspect")
  return(sin(asp))
}

#' Create northness layer.
#' @description Create a layer depicting northness based on an elevation layer.
#' @param dem RasterLayer object of elevation (a digital elevation model - DEM) as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples data(red.layers)
#' raster::plot(raster.north(red.layers[[3]]))
#' @export
raster.north <- function(dem){
  asp <- terrain(dem, opt = "aspect")
  return(cos(asp))
}

#' Predict species distribution.
#' @description Prediction of potential species distributions using maximum entropy (maxent).
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layers Predictor variables, a Raster* object as defined by package raster.
#' @param error Vector of spatial error in longlat (one element per row of longlat) in the same unit as longlat. Used to move any point randomly within the error radius.
#' @param categorical Vector of layer indices of categorical (as opposed to quantitative) data. If NULL the package will try to find them automatically based on the data.
#' @param thres Threshold of logistic output used for conversion of probabilistic to binary (presence/absence) maps. If 0 this will be the value that maximizes the sum of sensitivity and specificity.
#' @param testpercentage Percentage of records used for testing only. If 0 all records will be used for both training and testing.
#' @param mcp Used for a precautionary approach. If TRUE, all areas predicted as present but outside the minimum convex hull polygon encompassing all occurrence records are converted to absence. Exceptions are cells connected to other areas inside the polygon.
#' @param eval If TRUE, build a matrix with AUC, Kappa, TSS, EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @param runs If <= 0 no ensemble modelling is performed. If > 0, ensemble modelling with n runs is made. For each run, a new random sample of occurrence records (if testpercentage > 0), background points and predictive variables (if subset > 0) are chosen. In the ensemble model, each run is weighted as max(0, (runAUC - 0.5)) ^ 2.
#' @param subset Number of predictive variables to be randomly selected from layers for each run if runs > 0. If <= 0 all layers are used on all runs. Using a small number of layers is usually better than using many variables for rare species, with few occurrence records (Lomba et al. 2010, Breiner et al. 2015).
#' @details Builds maxent (maximum entropy) species distribution models (Phillips et al. 2004, 2006; Elith et al. 2011) using function maxent from R package dismo (Hijmans et al. 2017). Dismo requires the MaxEnt species distribution model software, a java program that can be downloaded from http://biodiversityinformatics.amnh.org/open_source/maxent. Copy the file 'maxent.jar' into the 'java' folder of the dismo package. That is the folder returned by system.file("java", package="dismo"). You need MaxEnt version 3.3.3b or higher. Please note that this program (maxent.jar) cannot be redistributed or used for commercial or for-profit purposes.
#' @return List with either one or two raster objects (depending if ensemble modelling is performed, in which case the second is a probabilistic map from all the runs) and, if eval = TRUE, a matrix with AUC, Kappa, TSS, EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model). Aggregate values are taken from maps after transformation of probabilities to incidence, with presence predicted for cells with ensemble values > 0.5.
#' @references Breiner, F.T., Guisan, A., Bergamini, A., Nobis, M.P. (2015) Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218.
#' @references Hijmans, R.J., Phillips, S., Leathwick, J., Elith, J. (2017) dismo: Species Distribution Modeling. R package version 1.1-4. https://CRAN.R-project.org/package=dismo
#' @references Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J., Guisan, A. (2010) Overcoming the rare species modelling paradox: a novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143: 2647-2657.
#' @references Phillips, S.J., Dudik, M., Schapire, R.E. (2004) A maximum entropy approach to species distribution modeling. Proceedings of the Twenty-First International Conference on Machine Learning. p. 655-662.
#' @references Phillips, S.J., Anderson, R.P., Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling, 190: 231-259.
#' @references Elith, J., Phillips, S.J., Hastie, T., Dudik, M., Chee, Y.E., Yates, C.J. (2011) A statistical explanation of MaxEnt for ecologists. Diversity and Distributions, 17: 43-57.
#' @export
map.sdm <- function(longlat, layers, error = NULL, categorical = NULL, thres = 0, testpercentage = 0, mcp = TRUE, eval = TRUE, runs = 0, subset = 0){

  raster::rasterOptions(maxmemory = 2e+09)

  ##if ensemble is to be done
  if(runs > 0){
    if(eval)
      runEval = matrix(NA, nrow = 1, ncol = 7)
    runMap <- rasterize(longlat, layers[[1]], field = 0, background = 0)
    pb <- txtProgressBar(min = 0, max = runs, style = 3)
    totalAUC = 0
    for(i in 1:runs){
      if(subset > 0 && subset < dim(layers)[3]){
        runLayers <- layers[[sample.int(dim(layers)[3], subset)]]
        thisRun <- map.sdm(longlat, runLayers, error, categorical, thres, testpercentage, mcp, eval, runs = 0, subset = 0)
      } else {
        thisRun <- map.sdm(longlat, layers, error, categorical, thres, testpercentage, mcp, eval, runs = 0, subset = 0)
      }
      runAUC = 1
      if(eval){
        runAUC <- thisRun[[2]][1]
        runAUC <- max(0, (runAUC - 0.5)) ^ 2      #weight the map by its AUC above 0.5 to the square
        runEval <- rbind(runEval, thisRun[[2]])
        thisRun <- thisRun[[1]]
      }
      totalAUC = totalAUC + runAUC
      runMap <- runMap + (thisRun * runAUC)
      setTxtProgressBar(pb, i)
    }
    runMap <- raster::calc(runMap, function(x) {x/totalAUC})
    upMap <- reclassify(runMap, matrix(c(0,0.025,0,0.025,1,1), ncol = 3, byrow = TRUE))
    consensusMap <- reclassify(runMap, matrix(c(0,0.499,0,0.499,1,1), ncol = 3, byrow = TRUE))
    downMap <- reclassify(runMap, matrix(c(0,0.975,0,0.975,1,1), ncol = 3, byrow = TRUE))
    if(mcp && aoo(consensusMap) >= 4)
      consensusMap <- map.habitat(longlat, consensusMap, mcp = TRUE, eval = FALSE)

    if(eval){
      runEval <- runEval[-1,]
      clEval <- matrix(NA, nrow = 3, ncol = 7)
      colnames(clEval) <- c("AUC", "Kappa", "TSS", "EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
      rownames(clEval) <- c("UpCL", "Consensus", "LowCL")
      clEval[1,] <- apply(runEval, 2,  quantile, probs= 0.975, na.rm = TRUE)
      clEval[2,] <- apply(runEval, 2,  quantile, probs= 0.5, na.rm = TRUE)
      clEval[3,] <- apply(runEval, 2,  quantile, probs= 0.025, na.rm = TRUE)
      clEval[1:3,4] <- eoo(longlat)
      clEval[1:3,6] <- aoo(longlat)
      clEval[1,5] <- eoo(upMap)
      clEval[1,7] <- aoo(upMap)
      clEval[2,5] <- eoo(consensusMap)
      clEval[2,7] <- aoo(consensusMap)
      clEval[3,5] <- eoo(downMap)
      clEval[3,7] <- aoo(downMap)
      return(list(consensusMap, runMap, clEval))
    } else {
      return (consensusMap)
    }
  }

  #if there is error randomly move points within its radius
  if(!is.null(error)){
    for(i in 1:nrow(longlat)){
      #move up to given error (angular movement converted to x and y)
      rndAngle = sample(1:360, 1)
      rndDist = runif(1, 0, error[i])
      longlat[i,1] = longlat[i,1] + rndDist * cos(rndAngle)
      longlat[i,2] = longlat[i,2] + rndDist * sin(rndAngle)
    }
  }
  longlat <- move(longlat, layers)  #move all records falling on NAs

  nPoints = min(1000, sum(!is.na(as.vector(layers[[1]])), na.rm = TRUE)/4)
  bg <- dismo::randomPoints(layers, nPoints)                                ##extract background points

  ##if no categorical variables are given try to figure out which
  if(is.null(categorical))
    categorical <- find.categorical(layers)

  llTrain <- longlat
  llTest <- longlat
  if(testpercentage > 0){
    testRecords <- sample(1:nrow(longlat), ceiling(nrow(longlat)*testpercentage/100))
    llTrain <- longlat[-testRecords,]
    llTest <- longlat[testRecords,]
  }
  mod <- dismo::maxent(layers, llTrain, a = bg, factors = categorical) ##build model
  p <- raster::predict(mod, layers)                                     ##do prediction

  e <- dismo::evaluate(p = llTrain, a = bg, model = mod, x = layers)         ##do evaluation of model
  if(thres == 0)
    thres <- dismo::threshold(e)$spec_sens                                   ##extract threshold from evaluation
  p <- reclassify(p, matrix(c(0,thres,0,thres,1,1), nrow=2, byrow = TRUE))  ##convert to presence/absence

  if(mcp && aoo(p) >= 4)
    p <- map.habitat(longlat, p, mcp = TRUE, eval = FALSE)

  if(eval){
    e <- dismo::evaluate(p = llTest, a = bg, model = mod, x = layers, tr = thres)                  ##do evaluation of model with threshold
    auc <- e@auc
    kappa <- e@kappa
    sensitivity <- as.numeric(e@TPR/(e@TPR+e@FNR))
    specificity <- as.numeric(e@TNR/(e@TNR+e@FPR))
    tss <- sensitivity + specificity - 1
    eooRaw <- eoo(longlat)
    aooRaw <- aoo(longlat)
    aooModel <- aoo(p)
    if(aooModel > 8)
      eooModel <- eoo(p)
    else
      eooModel = aooModel
    txtEval <- matrix(c(auc, kappa, tss, eooRaw, eooModel, aooRaw, aooModel), nrow = 1)
    colnames(txtEval) <- c("AUC", "Kappa", "TSS", "EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
    return(list(p, txtEval))
  } else {
    return(p)
  }
}

#' Map species distribution of habitat specialist.
#' @description Mapping of all habitat patches where the species is known to occur.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence (1/0) of a single habitat type.
#' @param move If TRUE, identifies and moves presence records to closest cells with suitable habitat. Use when spatial error might put records outside the correct patch.
#' @param mcp If TRUE, all habitat patches inside the minimum convex hull polygon encompassing all occurrence records are converted to presence.
#' @param eval If TRUE, build a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @details In many cases a species has a very restricted habitat and we generally know where it occurs. In such cases using the distribution of the known habitat patches may be enough to map the species.
#' @return One raster object and, if eval = TRUE, a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @export
map.habitat <- function(longlat, layer, move = TRUE, mcp = FALSE, eval = TRUE){
  if(move){
    moveLayer <- layer
    moveLayer[moveLayer == 0] <- NA
    longlat <- move(longlat, moveLayer)
    remove(moveLayer)
  }

  if(mcp){
    vertices <- chull(longlat)
    vertices <- c(vertices, vertices[1])
    vertices <- longlat[vertices,]
    poly = Polygon(vertices)
    poly = Polygons(list(poly),1)
    poly = SpatialPolygons(list(poly))    ##minimum convex polygon
    patches <- raster::clump(layer, gaps=FALSE)       ##individual patches, numbered
    selPatches <- raster::unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside polygon
  } else {
    patches <- raster::clump(layer, gaps=FALSE)       ##individual patches, numbered
    selPatches <- raster::unique(extract(patches, longlat, df = TRUE, weights = TRUE)$clumps) ##which patches have the species
  }
  selPatches <- selPatches[!is.na(selPatches)]
  allPatches <- raster::unique(patches)
  allPatches <- as.data.frame(cbind(allPatches, rep(0, length(allPatches))))
  colnames(allPatches) <- c("patches", "selected")
  allPatches[selPatches, 2] <- 1
  patches <- raster::subs(patches, allPatches)
  layer <- mask(layer, patches, maskvalue = 0, updatevalue = 0)
  if(eval){
    eooRaw <- eoo(longlat)
    eooModel <- eoo(layer)
    aooRaw <- aoo(longlat)
    aooModel <- aoo(layer)
    txtEval <- matrix(c(eooRaw, eooModel, aooRaw, aooModel), nrow = 1)
    colnames(txtEval) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)")
    return(list(layer, txtEval))
  } else {
    return(layer)
  }
}

#' Map recorded distribution of species.
#' @description Mapping of all cells where the species is known to occur.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layers Raster* object as defined by package raster. Any raster with the relevant extent and cell size can be used.
#' @param eval If TRUE, build a matrix with EOO and AOO calculated from occurrence records only.
#' @details To be used if either information on the species is very scarce (and it is not possible to model the species distribution) or, on the contrary, complete (and there is no need to model the distribution).
#' @return One raster object and, if EVAL = TRUE, a matrix with EOO and AOO.
#' @examples
#' data(red.records)
#' data(red.layers)
#' raster::plot(map.points(red.records, red.layers, eval = FALSE))
#' points(red.records)
#' @export
map.points <- function(longlat, layers, eval = TRUE){
  p <- rasterize(longlat, layers[[1]], field = 1, background = 0)
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  p <- mask(p, maskLayer)
  if(eval){
    eooRaw <- eoo(longlat)
    aooRaw <- aoo(longlat)
    txtEval <- matrix(c(eooRaw, aooRaw), nrow = 1)
    colnames(txtEval) <- c("EOO", "AOO")
    return(list(p, txtEval))
  } else {
    return(p)
  }
}

## #' Predict species distribution change in time.
## #' @description Prediction and projection in time of potential species distribution using maximum entropy (maxent).
## #' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
## #' @param layers Raster* object as defined by package raster, used for modelling current conditions.
## #' @param projectLayers Raster* object as defined by package raster, used for projection into the future.
## #' @param bg Background data as a matrix of longitude and latitude (two columns). If not defined 1000 points will be randomly selected.
## #' @param categorical Vector of layer indices of categorical (as opposed to quantitative) data. If NULL the package will try to find them automatically based on the data itself.
## #' @param thres Threshold of logistic output used for conversion of probabilistic to binary (presence/absence) maps. If 1 this will be the value that maximizes the sum of sensitivity and specificity.
## #' @param polygon Used for a precautionary approach. If TRUE, all areas predicted as present but outside the minimum convex hull polygon encompassing all occurrence records are converted to absence. Only cells connected to other areas inside the polygon are kept for both present and future projections.
## #' @param runs If <= 0 no ensemble modelling is performed. If > 0, ensemble modelling with n runs is made. For each run a bootstrap (random sampling with replacement) of all occurrence records and a new set of 1000 background points are chosen.
## #' @details Builds maxent (maximum entropy) species distribution models (Phillips et al. 2004, 2006; Elith et al. 2011) using function maxent from R package dismo (Hijmans et al. 2017) for both the present and future.
## #' @return Either three or six RasterLayer (depending if ensemble modelling is performed, in which case the second trio are probabilistic maps from all the runs) and a matrix with Present AOO, Future AOO, Gain, Keep and Loss.
## #' @references Hijmans, R.J., Phillips, S., Leathwick, J., Elith, J. (2017) dismo: Species Distribution Modeling. R package version 1.1-4. https://CRAN.R-project.org/package=dismo
## #' @references Phillips, S.J., Dudik, M., Schapire, R.E. (2004) A maximum entropy approach to species distribution modeling. Proceedings of the Twenty-First International Conference on Machine Learning. p. 655-662.
## #' @references Phillips, S.J., Anderson, R.P., Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling 190:231-259.
## #' @references Elith, J., Phillips, S.J., Hastie, T., Dudik, M., Chee, Y.E., Yates, C.J. (2011) A statistical explanation of MaxEnt for ecologists. Diversity and Distributions 17:43-57.
## #' @export
## map.change <- function(longlat, layers, projectLayers, bg = NULL, categorical = NULL, thres = 1, polygon = FALSE, runs = 0){
##   options(warn=-1)
#   ##if ensemble modelling is to be done
#   if(runs > 0){
#     runEval = matrix(NA, nrow = 1, ncol = 5)
#     runMap <- rasterize(longlat, layers[[1]], field = 0, background = 0)
#     runMap <- raster::stack(runMap, runMap, runMap)
#     runMap01 <- runMap
#     pb <- txtProgressBar(min = 0, max = runs, style = 3)
#     for(i in 1:runs){
#       runData <- longlat[sample(nrow(longlat), replace = TRUE),]
#       thisRun <- map.change(runData, layers, projectLayers, bg, categorical, thres, polygon, runs = 0)
#       runEval <- rbind(runEval, thisRun[[4]])
#       thisRun <- list(thisRun[[1]], thisRun[[2]], thisRun[[3]])
#       for(j in 1:3)
#         runMap[[j]] <- runMap[[j]] + thisRun[[j]] / runs
#       setTxtProgressBar(pb, i)
#     }
#     for(i in 1:2)
#       runMap01[[i]] <- reclassify(runMap[[i]], matrix(c(0,0.5,0,0.5,1,1), ncol = 3, byrow = TRUE))
#     runMap01[[3]] <- runMap01[[2]] * 2 - runMap01[[1]] ##gain = 2, kept = 1, loss = -1, never exists = 0
#     runEval <- runEval[-1,]
#     clEval <- matrix(NA, nrow = 4, ncol = 5)
#     colnames(clEval) <- c("Present AOO", "Future AOO", "Gain", "Keep", "Loss")
#     rownames(clEval) <- c("Aggregate", "LowCL", "Median", "UpCL")
#     clEval[1,1] <- aoo(runMap01[[1]])
#     clEval[1,2] <- aoo(runMap01[[2]])
#     clEval[1,3] <- cellStats((raster::area(runMap01[[3]]) * subs(runMap01[[3]], as.data.frame(matrix(c(0,0,1,0,-1,0,2,1), ncol = 2, byrow = TRUE)))),sum)
#     clEval[1,4] <- cellStats((raster::area(runMap01[[3]]) * subs(runMap01[[3]], as.data.frame(matrix(c(0,0,1,1,-1,0,2,0), ncol = 2, byrow = TRUE)))),sum)
#     clEval[1,5] <- cellStats((raster::area(runMap01[[3]]) * subs(runMap01[[3]], as.data.frame(matrix(c(0,0,1,0,-1,1,2,0), ncol = 2, byrow = TRUE)))),sum)
#     clEval[2,] <- apply(runEval, 2,  quantile, probs= 0.025, na.rm = TRUE)
#     clEval[3,] <- apply(runEval, 2,  quantile, probs= 0.5, na.rm = TRUE)
#     clEval[4,] <- apply(runEval, 2,  quantile, probs= 0.975, na.rm = TRUE)
#     options(warn=0)
#     return(list(runMap01, runMap, clEval))
#   }
#
#   ##if no background points are given randomly sample them
#   if(is.null(bg))
#     bg <- randomPoints(layers, 1000)                                ##extract background points (to use as absence)
#   ##if no categorical variables are given try to figure out which
#   if(is.null(categorical))
#     categorical = find.categorical(layers)
#
#   model <- dismo::maxent(layers, longlat, a = bg, factors = categorical) ##build model
#   present <- raster::predict(model, layers)                             ##do prediction for the present
#   future <- raster::predict(model, projectLayers)
#   e <- dismo::evaluate(longlat, bg, model, layers)                       ##do evaluation of model
#   if(thres >= 1)
#     thres <- threshold(e)$spec_sens                                   ##extract threshold from evaluation
#   present <- reclassify(present, matrix(c(0,thres,0,thres,1,1), nrow=2, byrow = TRUE))  ##convert to presence/absence
#   future <- reclassify(future, matrix(c(0,thres,0,thres,1,1), nrow=2, byrow = TRUE))  ##convert to presence/absence
#
#   if(polygon){                          ##if species is limited in dispersal ability
#     vertices <- chull(longlat)
#     vertices <- c(vertices, vertices[1])
#     vertices <- longlat[vertices,]
#     poly = Polygon(vertices)
#     poly = Polygons(list(poly),1)
#     poly = SpatialPolygons(list(poly))    ##original EOO, before modelling
#     ##present
#     patches <- clump(present, gaps=FALSE)       ##individual patches, numbered, present
#     selPatches <- unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside original EOO
#     present <- replace(present, !(patches %in% selPatches), 0)
#     ##future
#     patches <- clump(future, gaps=FALSE)       ##individual patches, numbered, future
#     selPatches <- unique(extract(patches, poly, df = TRUE, weights = TRUE)$clumps) ##which patches are inside original EOO
#     future <- replace(future, !(patches %in% selPatches), 0)
#   }
#
#   spDiff <- future * 2 - present ##gain = 2, kept = 1, loss = -1, never exists = 0
#   txtChange <- rep(NA,5)
#   names(txtChange) <- c("Present AOO", "Future AOO", "Gain", "Keep", "Loss")
#   txtChange[1] <- cellStats((raster::area(present) * present), sum)
#   txtChange[2] <- cellStats((raster::area(future) * future), sum)
#   txtChange[3] <- cellStats((raster::area(spDiff) * subs(spDiff, as.data.frame(matrix(c(0,0,1,0,-1,0,2,1), ncol = 2, byrow = TRUE)))),sum)
#   txtChange[4] <- cellStats((raster::area(spDiff) * subs(spDiff, as.data.frame(matrix(c(0,0,1,1,-1,0,2,0), ncol = 2, byrow = TRUE)))),sum)
#   txtChange[5] <- cellStats((raster::area(spDiff) * subs(spDiff, as.data.frame(matrix(c(0,0,1,0,-1,1,2,0), ncol = 2, byrow = TRUE)))),sum)
#   return(list(present, future, spDiff, txtChange))
# }

#' Species distributions made easy (multiple species).
#' @description Single step for prediction of multiple species distributions. Output of maps (in pdf format), klms (for Google Earth) and relevant data (in csv format).
#' @param longlat data.frame of taxon names, longitude and latitude or eastness and northness (three columns in this order) of each occurrence record.
#' @param layers If NULL analyses are done with environmental layers read from data files of red.setup(). If a Raster* object as defined by package raster, analyses use these.
#' @param habitat Raster* object as defined by package raster. Habitat extent layer (0/1) used instead of layers if any species is an habitat specialist.
#' @param zone UTM zone if data is in metric units. Used only for correct placement of kmls and countries.
#' @param error Vector of spatial error in longlat (one element per row of longlat) in the same unit as longlat. Used to move any point randomly within the error radius.
#' @param move If TRUE, identifies and moves presence records to closest cells with environmental data. Use when spatial error might put records outside such data.
#' @param dem RasterLayer object. It should be a digital elevation model for calculation of elevation limits of the species. If NULL, dem from red.setup() is used if possible, otherwise it will be 0.
#' @param pca Number of pca axes for environmental data reduction. If 0 (default) no pca is made.
#' @param filename Name of output csv file with all results. If NULL it is named "Results_All.csv".
#' @param mapoption Vector of values within options: points, habitat and sdm; each value corresponding to the function to be used for each species (map.points, map.habitat, map.sdm). If a single value, all species will be modelled according to it. If NULL, the function will perform analyses using map.points. Species values must be in same order as latlong.
#' @param testpercentage Percentage of records used for testing only. If 0 all records will be used for both training and testing.
#' @param mintest Minimim number of total occurrence records of any species to set aside a test set. Only used if testpercentage > 0.
#' @param runs If <= 0 no ensemble modelling is performed. If > 0, ensemble modelling with n runs is made. For each run, a new random sample of occurrence records (if testpercentage > 0), background points and predictive variables (if subset > 0) are chosen. In the ensemble model, each run is weighted as max(0, (runAUC - 0.5)) ^ 2.
#' @param subset Number of predictive variables to be randomly selected from layers for each run if runs > 0. If <= 0 all layers are used on all runs. Using a small number of layers is usually better than using many variables for rare species, with few occurrence records (Lomba et al. 2010, Breiner et al. 2015).
#' @return Outputs maps in asc, pdf and kml format, plus a file with EOO, AOO and a list of countries where the species is predicted to be present if possible to extract.
#' @references Breiner, F.T., Guisan, A., Bergamini, A., Nobis, M.P. (2015) Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218.
#' @references Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J., Guisan, A. (2010) Overcoming the rare species modelling paradox: a novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143: 2647-2657.
#' @export
map.easy <- function(longlat, layers = NULL, habitat = NULL, zone = NULL, error = NULL, move = TRUE, dem = NULL, pca = 0, filename = NULL, mapoption = NULL, testpercentage = 0, mintest = 20, runs = 0, subset = 0){

  try(dev.off(), silent = TRUE)
  spNames <- unique(longlat[,1])
  nSp <- length(spNames)

  if(is.null(mapoption))
    mapoption = rep("points", nSp)
  else if(length(mapoption) == 1)
    mapoption = rep(mapoption, nSp)
  else if(length(mapoption) != nSp)
    return(warning("Number of species different from length of mapoption"))

  if("sdm" %in% mapoption){
    if(!file.exists(paste(.libPaths()[[1]], "/dismo/java/maxent.jar", sep=""))){
      warnMaxent()
      return()
    }
  }

  if (all(mapoption == rep("points", nSp))){
    res <- matrix(NA, nrow = nSp, ncol = 5)
    colnames(res) <- c("EOO", "AOO", "Min elevation", "Max elevation", "Countries")
  } else if (("sdm" %in% mapoption) && runs > 0) {
    res <- matrix(NA, nrow = nSp, ncol = 11)
    colnames(res) <- c("EOO (raw)", "EOO (LowCL)", "EOO (Consensus)", "EOO (UpCL)", "AOO (raw)", "AOO (LowCL)", "AOO (Consensus)", "AOO (UpCL)", "Min elevation", "Max elevation", "Countries")
  } else {
    res <- matrix(NA, nrow = nSp, ncol = 7)
    colnames(res) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)", "Min elevation", "Max elevation", "Countries")
  }
  rownames(res) <- spNames

  if(is.null(layers))
    newLayers <- TRUE
  else
    newLayers <- FALSE
  if(is.null(dem))
    newDem <- TRUE
  else
    newDem <- FALSE

  rad = 0.1
  for(s in 1:nSp){
    cat("\nSpecies", s, "of", nSp, "-", toString(spNames[s]),"\n")
    spData <- longlat[longlat[,1] == spNames[s], -1]
    if(!is.null(error)){
      spError <- error[longlat[,1] == spNames[s]]
      if(max(spError) > 1)
        rad <- spError/100000
      else
        rad <- spError
    } else {
      spError <- NULL
    }
    if(newLayers){
      layers <- raster.read(spData)
      if(newDem)
        dem <- layers[[20]]
      if(pca > 0)
        layers <- raster.reduce(layers, n = pca)
    }

    if(mapoption[s]  == "sdm" && aoo(move(spData, layers)) > 8){
      if(move)
        spData <- move(spData, layers)
      if(testpercentage > 0)
        p <- map.sdm(spData, layers, spError, testpercentage = testpercentage, mcp = TRUE, runs = runs, subset = subset)
      else
        p <- map.sdm(spData, layers, spError, testpercentage = 0, mcp = TRUE, runs = runs, subset = subset)
    } else if (mapoption[s] == "habitat"){
      p <- map.habitat(spData, habitat, move)
    } else {
      mapoption[s] = "points"
      p <- map.points(spData, layers)
    }
    writeRaster(p[[1]], paste(toString(spNames[s]), ".asc", sep=""), overwrite = TRUE)
    map.draw(spData, p[[1]], spNames[s], sites = FALSE, print = TRUE)
    if(mapoption[s] != "points"){
      kml(p[[1]], zone = zone, paste(toString(spNames[s]), ".kml", sep=""), mapoption = "aoo")
      countryList <- countries(p[[1]], zone = zone)
      if(is.null(dem))
        elev <- c(0, 0)
      else
        elev <- elevation(p[[1]], dem)
    } else {
      kml(spData, zone = zone, paste(toString(spNames[s]), ".kml", sep=""), mapoption = "points", rad = rad)
      countryList <- countries(spData, zone = zone)
      if(is.null(dem))
        elev <- c(0, 0)
      else
        elev <- elevation(spData, dem)
    }

    if(mapoption[s]  == "sdm" && aoo(spData) > 8 && runs > 0){
      writeRaster(p[[2]], paste(toString(spNames[s]), "_prob.asc", sep=""), overwrite = TRUE)
      map.draw(spData, p[[2]], paste(toString(spNames[s]), "_prob", sep = ""), legend = TRUE, print = TRUE)
    }

    ##write output values to csv
    spRes = p[[length(p)]]
    if(ncol(res) == 5){      #colnames(res) <- c("EOO", "AOO", "Min elevation", "Max elevation", "Countries")
      res[s,] <- c(spRes, elev, toString(countryList))
    }
    if(ncol(res) == 7){      #colnames(res) <- c("EOO (raw)", "EOO (model)", "AOO (raw)", "AOO (model)", "Min elevation", "Max elevation", "Countries")
      if(length(spRes) == 7)
        res[s,] <- c(spRes[4:7], elev, toString(countryList))
      else  #if length(spRes) < 7
        res[s,] <- c(spRes[c(1,1,2,2)], elev, toString(countryList))
    }
    if(ncol(res) == 11){     #colnames(res) <- c("EOO (raw)", "EOO (LowCL)", "EOO (Consensus)", "EOO (UpCL)", "AOO (raw)", "AOO (LowCL)", "AOO (Consensus)", "AOO (UpCL)", "Min elevation", "Max elevation", "Countries")
      if(length(spRes) == 2)
        res[s,] <- c(spRes[c(1,1,1,1,2,2,2,2)], elev, toString(countryList))
      else if(length(spRes) == 4)
        res[s,] <- c(spRes[c(1,2,2,2,3,4,4,4)], elev, toString(countryList))
      else if(is.null(dim(spRes)))
        res[s,] <- c(spRes[4:7], elev, toString(countryList))
      else   #if matrix
        res[s,] <- c(spRes[2,4], spRes[3:1,5], spRes[2,6], spRes[3:1,7], elev, toString(countryList))
    }
    write.csv(res[s,], paste(toString(spNames[s]), ".csv", sep = ""))
    if(mapoption[s]  == "sdm" && aoo(spData) > 8){
      if(runs > 0)
        write.csv(p[[3]], paste(toString(spNames[s]), "_detail.csv", sep = ""))
      else
        write.csv(p[[2]], paste(toString(spNames[s]), "_detail.csv", sep = ""))
    }


  }
  if(is.null(filename))
    write.csv(res, "Results_All.csv")
  else
    write.csv(res, toString(filename))
  return(as.data.frame(res))
}

#' Map creation.
#' @description Creates maps ready to print in pdf or other formats.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence map for the species.
#' @param spName String of species name.
#' @param borders If TRUE country borders are drawn.
#' @param scale If TRUE a distance scale in km is drawn.
#' @param legend If TRUE the legend for the map is drawn.
#' @param sites If TRUE the record locations are drawn.
#' @param mcp If TRUE the minimum convex polygon representing the Extent of Occurrence is drawn.
#' @param print If TRUE a pdf is saved instead of the output to the console.
#' @examples data(red.records)
#' data(red.range)
#' par(mfrow = c(1,2))
#' map.draw(red.records, layer = red.range, mcp = TRUE)
#' @export
map.draw <- function(longlat = NULL, layer, spName,  borders = FALSE, scale = TRUE, legend = FALSE, sites = TRUE, mcp = FALSE, print = FALSE){
  worldborders <- NULL
  data(worldborders, envir = environment())
  if (borders){
    layer[layer == 0] <- NA
    raster::plot(layer, main = spName, legend = legend, xlab = "longitude", ylab = "latitude", col = "forestgreen")
    lines(worldborders)
  } else {
    raster::plot(layer, main = spName, legend = legend, colNA = "lightblue", xlab = "longitude", ylab = "latitude")
  }
  if (scale){
    width = (xmax(layer) - xmin(layer))
    d = round(width/10^(nchar(width)-1))*10^(nchar(width)-2)
    scalebar(d = d, type="bar", divs = 2)
  }
  if (sites && !is.null(longlat))
    points(longlat)
  if (mcp){
    e <- rasterToPoints(layer, fun = function(dat){dat == 1})   ##convert raster to points
    vertices <- chull(e[,1], e[,2])
    vertices <- c(vertices, vertices[1])
    vertices <- e[vertices,c(1,2)]
    poly <- SpatialPolygons(list(Polygons(list(Polygon(vertices)),1)))
    raster::plot(poly, add = TRUE)
  }
  if(print){
    dev.copy(device = pdf, file = paste(toString(spName), ".pdf", sep=""))
    dev.off()
  }
}

#' Extent of Occurrence (EOO).
#' @description Calculates the Extent of Occurrence of a species based on either records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @details EOO is calculated as the minimum convex polygon covering all known or predicted sites for the species.
#' @return A single value in km2.
#' @examples data(red.records)
#' data(red.range)
#' eoo(red.records)
#' eoo(red.range)
#' @export
eoo <- function(spData){
  if(class(spData) == "RasterLayer"){
    if (raster::xmax(spData) <= 180) {  #if longlat data
      e <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
      vertices <- chull(e[,1], e[,2])
      if(length(vertices) < 3) return(0)
      vertices <- c(vertices, vertices[1])
      vertices <- e[vertices,c(1,2)]
      area = geosphere::areaPolygon(vertices)/1000000
    } else {
      spData[spData < 1] <- NA
      spData <- rasterToPoints(spData)
      vertices <- chull(spData)
      if(length(vertices) < 3) return(0)
      vertices <- c(vertices, vertices[1])
      vertices <- spData[vertices,]
      area = 0
      for(i in 1:(nrow(vertices)-1))
        area = area + (as.numeric(vertices[i,1])*as.numeric(vertices[(i+1),2]) - as.numeric(vertices[i,2])*as.numeric(vertices[(i+1),1]))
      area = abs(area/2000000)
    }
  } else if (ncol(spData) == 2){
    vertices <- chull(spData)
    if(length(vertices) < 3) return(0)
    vertices <- c(vertices, vertices[1])
    vertices <- spData[vertices,]
    if(max(spData) <= 180) {  #if longlat data
      area = geosphere::areaPolygon(vertices)/1000000
    } else { #if square data in meters
      area = 0
      for(i in 1:(nrow(vertices)-1))
        area = area + (as.numeric(vertices[i,1])*as.numeric(vertices[(i+1),2]) - as.numeric(vertices[i,2])*as.numeric(vertices[(i+1),1]))
      area = abs(area/2000000)
    }
  } else {
    return(warning("Data format not recognized"))
  }
  return(area)
}

#' Area of Occupancy (AOO).
#' @description Calculates the Area of Occupancy of a species based on either known records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @details AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as required by IUCN.
#' @return A single value in km2.
#' @examples data(red.range)
#' aoo(red.range)
#' @export
aoo <- function(spData){
  if (class(spData) == "RasterLayer"){ #if rasterlayer with 0/1
    if(raster::maxValue(spData) == 0){  #if no data (empty raster)
      area = 0
    } else if (raster::xmax(spData) <= 180) {  #if longlat data
      #area = cellStats((raster::area(spData) * spData), sum)          #old version using cell area instead of 2x2 grid
      if(res(spData)[1] > 0.05)   #if resolution is 10km convert to 1km
        spData = disaggregate(spData, fact = 10)
      spData[spData < 1] <- NA
      spData <- rasterToPoints(spData)
      if(nrow(unique(spData)) == 1){
        area = 4
      } else {
        spData <- longlat2utm(spData[,-3])
        spData = floor(spData/2000)
        ncells = nrow(unique(spData))
        area = ncells * 4
      }
    } else { #if square data in meters
      spData[spData < 1] <- NA
      spData <- rasterToPoints(spData)
      spData = floor(spData/2000)
      ncells = nrow(unique(spData))
      area = ncells * 4
    }
  } else if (ncol(spData) == 2){
    if (max(spData) <= 180) {  #if longlat data
      #layer = raster.read(spData, ext = 0.1)[[1]]                     #old version using cell area instead of 2x2 grid
      #gisdir = red.getDir()
      #layer = raster(paste(gisdir, "red_2km_1.tif", sep = ""))
      #spData <- spData[!is.na(extract(layer, spData)),]
      #layer = rasterize(spData, layer, field = 1)
      #area = cellStats((raster::area(layer) * layer), sum)
      spData <- longlat2utm(spData)
      spData = floor(spData/2000)
      ncells = nrow(unique(spData))
      area = ncells * 4
    } else { #if square data in meters
      spData = floor(spData/2000)
      ncells = nrow(unique(spData))
      area = ncells * 4
    }
  } else {
    return(warning("Data format not recognized!"))
  }
  return(area)
}

#' Elevation limits.
#' @description Calculates the elevation (or depth) limits (range) of a species based on either known records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @param dem RasterLayer object. Should be a digital elevation model (DEM) of the relevant area. If not given the function will try to read it from base data, only works with longlat data.
#' @details Maximum and minimum elevation are calculated based on the DEM.
#' @return A vector with two values (min and max) in meters above (or below) sea level.
#' @examples data(red.records)
#' data(red.range)
#' data(red.layers)
#' dem = red.layers[[3]]
#' elevation(red.records, dem)
#' elevation(red.range, dem)
#' @export
elevation <- function(spData, dem = NULL){
  if(class(spData) != "RasterLayer"){ #if no rasterlayer is given but just a matrix of longlat.
    if(is.null(dem) && max(spData) <= 180){
      gisdir = red.getDir()
      dem <- raster::raster(paste(gisdir, "red_1km_20.tif", sep =""))
      dem <- crop(dem, c(min(spData[,1])-0.1, max(spData[,1]+0.1), min(spData[,2])-0.1, max(spData[,2])+0.1))
    }
    spData = rasterize(spData, dem, field = 1, background = NA) #create a layer of presence based on the dem
  } else if (is.null(dem)){
    gisdir = red.getDir()
    dem <- raster::raster(paste(gisdir, "red_1km_20.tif", sep = ""))
    dem <- crop(dem, spData)
  }
  spData[spData == 0] <- NA
  spData <- raster::overlay(spData, dem, fun = function(x,y){(x*y)})
  out <- c(raster::minValue(spData), raster::maxValue(spData))
  names(out) <- c("Min", "Max")
  return(out)
}

#' Countries of occurrence.
#' @description Extracts the names or ISO codes of countries of occurrence of a species based on either records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @param zone UTM zone if data is in metric units.
#' @param ISO Outputs either country names (FALSE) or ISO codes (TRUE).
#' @details Country boundaries and designations are based on data(worldborders) from package maptools.
#' @return A vector with country names or codes.
#' @examples data(red.records)
#' data(red.range)
#' countries(red.records)
#' countries(red.range, ISO = TRUE)
#' @export
countries <- function(spData, zone = NULL, ISO = FALSE){
  if ((class(spData) == "RasterLayer" && raster::xmax(spData) > 180) || (class(spData) != "RasterLayer" && max(spData) > 180))   ##if need to project to longlat
    spData <- utm2longlat(spData, zone)

  worldborders <- NULL
  data(worldborders, envir = environment())
  if(class(spData) == "RasterLayer")
    spData <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
  countryList <- sp::over(sp::SpatialPoints(spData), sp::SpatialPolygons(worldborders@polygons))
  if(ISO)
    countryList <- unique(worldborders@data[countryList,])$ISO2
  else
    countryList <- unique(worldborders@data[countryList,])$NAME
  countryList <- sort(as.vector(countryList[!is.na(countryList)]))
  return(countryList)
}

#' Output kml files.
#' @description Creates kml files for Google Maps as required by IUCN guidelines.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (0/1 values).
#' @param zone UTM zone if data is in metric units.
#' @param filename The name of file to save, should end with .kml.
#' @param mapoption Type of representation, any of "points", "eoo" or "aoo".
#' @param rad radius of circles in degrees if mapoption is "points". It can be the same value for all points or a vector with length equal to number of records in spData representing associated error. The default is about 10km (0.1 degrees) as per IUCN guidelines.
#' @return A kml with polygon or circles around records.
#' @export
kml <- function(spData, zone = NULL, filename, mapoption = "aoo", rad = 0.1){
  if ((class(spData) == "RasterLayer" && raster::xmax(spData) > 180) || (class(spData) != "RasterLayer" && max(spData) > 180))   ##if need to project to longlat
    spData <- utm2longlat(spData, zone)

  if(mapoption == "aoo" && class(spData) == "RasterLayer"){
    spData[spData != 1] <- NA
    spData <- rasterToPolygons(spData, dissolve = TRUE)
    writeOGR(spData, filename, layer = filename, overwrite_layer = TRUE, driver = "KML")
  } else if(mapoption == "points" || (class(spData) == "RasterLayer" && aoo(spData) <= 8) || nrow(spData) < 3){
    poly = list()
    for(i in 1:nrow(spData)){
      pts = seq(0, 2 * pi, length.out = 100)
      if(length(rad) == 1)
        xy = cbind(spData[i, 1] + rad * sin(pts), spData[i, 2] + rad * cos(pts))
      else
        xy = cbind(spData[i, 1] + rad[i] * sin(pts), spData[i, 2] + rad[i] * cos(pts))
      poly[[i]] = Polygon(xy)
    }
    poly = Polygons(poly,1)
    kmlPolygon(poly, filename, border = "red")
  } else {
    if (class(spData) == "RasterLayer"){
      e <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
      vertices <- chull(e[,1], e[,2])
      vertices <- c(vertices, vertices[1])
      vertices <- e[vertices,c(1,2)]
    } else {
      vertices <- chull(spData)
      vertices <- c(vertices, vertices[1])
      vertices <- spData[vertices,]
    }
    poly = Polygon(vertices)
    poly = Polygons(list(poly),1)
    kmlPolygon(poly, filename, border = "red")
  }
}

#' Red List Index.
#' @description Calculates the Red List Index (RLI) for a group of species.
#' @param spData Either a vector with species assessment categories for a single point in time or a matrix with two points in time in different columns (species x date). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index and this approach is here implemented.
#' @return Either a vector (if no two dates are given) or a matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
#' @references Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.
#' @references Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
#' @examples rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX"), ncol = 2, byrow = TRUE)
#' colnames(rliData) <- c("2000", "2010")
#' rli(rliData[,1])
#' rli(rliData[,1], boot = TRUE)
#' rli(rliData)
#' rli(rliData, boot = TRUE)
#' @export
rli <- function (spData, boot = FALSE, runs = 1000){
  ##RLI with phylogenetic or functional data to be implemented soon
  ##rli <- function (spData, tree = NULL, boot = FALSE, runs = 1000){
  tree = NULL   ##to add soon

  ##if only one point in time is given
  if(is.null(dim(spData)))
    return(rli.calc(spData, tree, boot, runs))  ##return either 1 or 3 values

  ##if two points in time are given
  ts <- apply(spData, 2, function(x) rli.calc(x, tree, boot = FALSE))
  sl <- ts[2] - ts[1]
  if(!boot){
    res <- matrix(c(ts, sl), nrow = 1)
    colnames(res) <- c(colnames(spData), "Change")
    rownames(res) <- c("Raw")
    return(res)
  } else {
    tr <- apply(spData, 2, function(x) rli.calc(x, tree, boot, runs))
    p = 0
    rndSl = rep(NA, runs)
    for(r in 1:runs){
      rndSl[r] <- rli.calc(spData[,2], tree, boot, 1)[2] - rli.calc(spData[,1], tree, boot, 1)[2]
      if(sign(sl) < sign(rndSl[r]) || sign(sl) > sign(rndSl[r]))
        p = p + 1
    }
    p = p / runs
    rndSl = quantile(rndSl, c(0.025, 0.5, 0.975))
    res <- matrix(c(ts[1], tr[,1], ts[2], tr[,2], sl, rndSl), nrow = 4, ncol = 3)
    colnames(res) <- c(colnames(spData), "Change")
    rownames(res) <- c("Raw", "LowCL", "Median", "UpCL")
    return(list("Values" = res, "P_change" = p))
  }
}

#' Red List Index for multiple groups.
#' @description Calculates the Red List Index (RLI) for multiple groups of species.
#' @param spData A matrix with group names (first column) and species assessment categories for one or two points in time (remaining columns). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to 5 (Extinct/Extinct in the Wild).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index and this approach is here implemented.
#' @return A matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
#' @references Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.
#' @references Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
#' @examples rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX"), ncol = 2, byrow = TRUE)
#' colnames(rliData) <- c("2000", "2010")
#' rliData <- cbind(c("Arthropods","Arthropods","Birds","Birds","Birds"), rliData)
#' rli.multi(rliData[,1:2])
#' rli.multi(rliData[,1:2], boot = TRUE)
#' rli.multi(rliData)
#' rli.multi(rliData, boot = TRUE)
#' @export
rli.multi <- function (spData, boot = FALSE, runs = 1000){
  ##RLI with phylogenetic or functional data to be implemented soon
  ##rli.multi <- function (spData, tree = NULL, boot = FALSE, runs = 1000){
  tree = NULL ##to add soon xxx

  groups <- unique(spData[,1])
  nGroups <- length(groups)
  if(ncol(spData) == 2 && !boot){
    res <- matrix(NA, nrow = nGroups, ncol = 1)
  } else if((ncol(spData) == 2 && boot) || (ncol(spData) == 3 && !boot)){
    res <- matrix(NA, nrow = nGroups, ncol = 3)
  } else {
    res <- matrix(NA, nrow = nGroups, ncol = 13)
    colnames(res) <- c(paste(colnames(spData)[2], "(raw)"), paste(colnames(spData)[2], "(lowCL)"), paste(colnames(spData)[2], "(median)"), paste(colnames(spData)[2], "(upCL)"), paste(colnames(spData)[3], "(raw)"), paste(colnames(spData)[3], "(lowCL)"), paste(colnames(spData)[3], "(median)"), paste(colnames(spData)[3], "(upCL)"), "Change (raw)", "Change (lowCL)", "Change (median)", "Change (upCL)", "p (change)")
  }
  row.names(res) <- groups
  for(g in 1:nGroups){
    if(is.null(tree)){
      v <- rli(spData[spData[,1] == groups[g],-1], boot = boot, runs = runs)
      if(ncol(res) < 13){
        res[g,] <- v
        colnames(res) <- colnames(v)
      } else {
        res[g,1:4] <- v$Values[,1]
        res[g,5:8] <- v$Values[,2]
        res[g,9:12] <- v$Values[,3]
        res[g,13] <- v$P_change
      }
    } else {
    }
  }
  return(res)
}

#' Sampled Red List Index.
#' @description Calculates accumulation curve of confidence limits in sampled RLI.
#' @param spData A vector with species assessment categories for a single point in time. Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param p p-value of confidence limits (in a two-tailed test).
#' @param runs Number of runs for smoothing accumulation curves.
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Yet, in many groups, it is not possible to assess all species due to huge diversity and/or lack of resources. In such case, the RLI is estimated from a randomly selected sample of species - the Sampled Red List Index (SRLI; Stuart et al. 2010).
#' This function allows to calculate how many species are needed to reach a given maximum error of the SRLI around the true value of the RLI (with all species included) for future assessments of the group.
#' @return A vector with the accumulation of the error of the SRLI around the true value of the RLI (with all species included).
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PLoS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PLoS One, 2: e140.
#' @references Stuart, S.N., Wilson, E.O., McNeely, J.A., Mittermeier, R.A. & Rodriguez, J.P. (2010) The barometer of Life. Science 328, 117.
#' @examples rliData <- c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX")
#' rli.sampled(rliData)
#' @export
rli.sampled <- function (spData, p = 0.05, runs = 1000){
  ##RLI with phylogenetic or functional data to be implemented soon
  ##rli.sampled <- function (spData, tree = NULL, p = 0.05, runs = 1000){
  tree = NULL ##to add soon xxx

  nSpp <- length(spData)
  accum <- rep(NA, nSpp)
  for(n in 1:nSpp){               #test with n species from the entire set
    diff = rep(NA, runs)          #try runs times each species
    for(r in 1:runs){             #do r runs for each n species
      diff[r] = abs(rli.calc(spData, tree, FALSE, 1) - rli.calc(sample(spData, n), tree, FALSE, 1))  #calculate absolute difference between true and sampled rli for each run
    }
    accum[n] = quantile(diff, (1-p))
  }
  return(accum)   #returns the accumulation curve of confidence limit of sampled RLI
}

#' Occurrence records for Hogna maderiana (Walckenaer, 1837).
#'
#' Occurrence records for Hogna maderiana (Walckenaer, 1837).
#'
#' @docType data
#' @keywords datasets
#' @name red.records
#' @usage data(red.records)
#' @format Matrix of longitude and latitude (two columns) of occurrence records for Hogna maderiana (Walckenaer, 1837), a spider species from Madeira Island.
NULL

#' Geographic range for Hogna maderiana (Walckenaer, 1837).
#'
#' Geographic range for Hogna maderiana (Walckenaer, 1837).
#'
#' @docType data
#' @keywords datasets
#' @name red.range
#' @usage data(red.range)
#' @format RasterLayer object as defined by package raster of range for Hogna maderiana (Walckenaer, 1837), a spider species from Madeira Island.
NULL

#' Environmental layers for Madeira.
#'
#' Average annual temperature, total annual precipitation, altitude and landcover for Madeira Island (Fick & Hijmans 2017, Tuanmu & Jetz 2014).
#'
#' @docType data
#' @keywords datasets
#' @name red.layers
#' @usage data(red.layers)
#' @format RasterStack object as defined by package raster.
#' @references Fick, S.E. & Hijmans, R.J. (2017) Worldclim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, in press.
#' @references Tuanmu, M.-N. & Jetz, W. (2014) A global 1-km consensus land-cover product for biodiversity and ecosystem modeling. Global Ecology and Biogeography, 23: 1031-1045.
NULL

#'
#'
#' World country borders.
#'
#' World country borders.
#'
#' @docType data
#' @keywords datasets
#' @name worldborders
#' @usage data(worldborders)
#' @format SpatialPolygonsDataFrame.
NULL
