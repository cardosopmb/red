#####RED - IUCN Redlisting Tools
#####Version 1.5.0 (2020-05-04)
#####By Pedro Cardoso
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P.(2017) An R package to facilitate species red list assessments according to the IUCN criteria. Biodiversity Data Journal 5: e20530 doi: 10.3897/BDJ.5.e20530
#####Changed from v1.4.0:
#####added function rli.predict to interpolate and extrapolate linearly beyond the years assessed
#####added new options in functions rli and rli.multi on how to deal with DD species when bootstrapping

#####required packages
library("BAT")
library("dismo")
library("gdistance")
library("geosphere")
library("graphics")
library("grDevices")
library("jsonlite")
library("maptools")
library("methods")
library("raster")
library("rgdal")
library("rgeos")
library("sp")
library("stats")
library("utils")
#' @import gdistance
#' @import graphics
#' @import jsonlite
#' @import maptools
#' @import rgdal
#' @import rgeos
#' @import sp
#' @import stats
#' @import utils
#' @importFrom BAT contribution
#' @importFrom geosphere areaPolygon
#' @importFrom grDevices chull dev.copy dev.off pdf
#' @importFrom methods slot
#' @importFrom raster area cellStats clump crop extent extract getValues layerStats mask raster rasterize rasterToPoints rasterToPolygons reclassify res sampleRandom scalebar terrain trim writeRaster xmax xmin

raster::rasterOptions(maxmemory = 2e+09)
globalVariables(c("worldborders"))

###############################################################################
##############################AUX FUNCTIONS####################################
###############################################################################

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
    lay <- as.vector(lay)
    lay <- lay[!is.na(lay)]
    if(sum(floor(lay)) == sum(lay) && length(unique(lay)) < 50)
      categorical = c(categorical, l)
  }
  return(categorical)
}

##basic function to calculate the rli of any group of species
rli.calc <- function(spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000){
  if(all(is.na(spData)))
    return(NA)
  spData <- rli.convert(spData)                ##call function to convert spData to a 0-1 scale

  if(is.null(tree)){                           ##if not weighted by PD or FD
    if(!boot){                                 ##if no bootstrap to be made
      return (mean(spData, na.rm = TRUE))
    } else {
      run <- rep(NA, runs)
      if(!dd){
        for(i in 1:runs){
          rnd <- sample(spData, replace = TRUE) ##bootstrap with all species
          run[i] <- mean(rnd, na.rm = TRUE)
        }
      } else {                                       ##bootstrap with only DD species
        nDD = sum(is.na(spData))                     ##number of DD species
        rliBase = sum(spData, na.rm = TRUE)
        for(i in 1:runs){
          rnd <- sample(spData[!is.na(spData)], nDD, replace = TRUE)
          run[i] <- (rliBase + sum(rnd)) / length(spData)
        }
      }
      res <- matrix(quantile(run, c(0.025, 0.5, 0.975)), nrow = 1)
      colnames(res) <- c("LowCL", "Median", "UpCL")
      return(res)
    }
  } else {                                     ##if weighted by PD or FD, still to work, not available at the moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    comm <- matrix(1, nrow = 2, ncol = length(spData))
    contrib <- BAT::contribution(comm, tree, relative = TRUE)[1,]
    contrib <- contrib/sum(contrib[!is.na(spData)]) #needed to standardize the contribution by the total contribution of species living in the community
    if(!boot){                                 ##if no bootstrap to be made
      return(sum(spData * contrib, na.rm = TRUE))
    } else {
      run <- rep(NA, runs)
      for(i in 1:runs){
        rndSpp <- sample(length(spData), replace = TRUE)
        rndComm <- spData[rndSpp]
        rndContrib <- contrib[rndSpp]/sum(contrib[rndSpp])
        run[i] <- sum(rndComm * rndContrib, na.rm = TRUE)
      }
        res <- matrix(quantile(run, c(0.025, 0.5, 0.975)), nrow = 1)
        colnames(res) <- c("LowCL", "Median", "UpCL")
        return(res)
    }
  }
}

##function to convert strings to numbers in the RLI
rli.convert <- function(spData){
  if(!is.numeric(spData)){                                ##if letters are given, convert to [0,1]
    spData <- replace(spData, which(spData == "EX" ), 0)
    spData <- replace(spData, which(spData == "EW" ), 0)
    spData <- replace(spData, which(spData == "RE" ), 0)
    spData <- replace(spData, which(spData == "CR" ), 0.2)
    spData <- replace(spData, which(spData == "CR(PE)" ), 0.2)
    spData <- replace(spData, which(spData == "EN" ), 0.4)
    spData <- replace(spData, which(spData == "VU" ), 0.6)
    spData <- replace(spData, which(spData == "NT" ), 0.8)
    spData <- replace(spData, which(spData == "LC" ), 1)
    spData <- replace(spData, which(spData == "DD" ), NA)
    spData <- as.numeric(spData)
  } else if (all(spData == floor(spData))){  #if all integers, a scale [0,5] is given, convert to [0,1]
    spData <- 1 - spData/5
  }
  return(spData)
}

##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################

#' Setup GIS directory.
#' @description Setup directory where GIS files are stored.
#' @param gisPath Path to the directory where the gis files are stored.
#' @details Writes a txt file in the red directory allowing the package to always access the world GIS files directory.
#' @export
red.setDir <- function(gisPath = NULL){
  if(is.null(gisPath))
    gisPath <- readline("Input directory for storing world gis layers:")
  gisPath <- paste(gisPath, "/", sep = "")
  redFile <- paste(find.package("red"), "/red.txt", sep = "")
  dput(gisPath, redFile)
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
    if(nSites < 40){ #if limited number of sites use all data
      maxDist = 0
      for(x in 1:(nSites-1)){
        for(y in (x+1):nSites){
          maxDist = max(maxDist,((longlat[x,1]-longlat[y,1])^2+(longlat[x,2]-longlat[y,2])^2)^.5)
        }
      }
    } else { #if many sites use hypothenusa of square encompassing all of them
      horiDist = max(longlat[,1]) - min(longlat[,1])
      vertDist = max(longlat[,2]) - min(longlat[,2])
      maxDist = (horiDist^2 + vertDist^2)^0.5
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

#' Create distance layer.
#' @description Creates a layer depicting distances to records using the minimum, average, distance to the minimum convex polygon or distance taking into account a cost surface.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster to serve as model to create distance layer. Cost surface in case of param ="cost".
#' @param type text string indicating whether the output should be the "minimum", "average", "mcp" or "cost" distance to all records. "mcp" means the distance to the minimum convex polygon encompassing all records.
#' @details Using distance to records in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data(red.layers)
#' alt = red.layers[[3]]
#' data(red.records)
#' par(mfrow=c(3,2))
#' raster::plot(alt)
#' points(red.records)
#' raster::plot(raster.distance(red.records, alt))
#' raster::plot(raster.distance(red.records, alt, type = "average"))
#' raster::plot(raster.distance(red.records, alt, type = "mcp"))
#' raster::plot(raster.distance(red.records, alt, type = "cost"))
#' @export
raster.distance <- function(longlat, layers, type = "minimum"){
  if(dim(layers)[3] > 1)
    layers <- layers[[1]]
  layers[!is.na(layers)] <- 0
  if(type == "average"){
    for(d in 1:nrow(longlat)){
      layers <- layers + raster::distanceFromPoints(layers, longlat[d,])
    }
    layers <- layers/nrow(longlat)
    names(layers) <- "average distance"
  } else if (type == "mcp"){
    vertices <- chull(longlat)
    vertices <- c(vertices, vertices[1])
    vertices <- longlat[vertices,]
    poly = Polygon(vertices)
    poly = Polygons(list(poly),1)
    poly = SpatialPolygons(list(poly))    ##minimum convex polygon
    longlat = rasterToPoints(rasterize(poly, layers))[,1:2]
    layers <- mask(raster::distanceFromPoints(layers, longlat), layers)
    names(layers) <- "mcp distance"
  } else if (type == "cost"){
    layers <- transition(layers, function(x) 1/mean(x), 8)
    layers <- geoCorrection(layers)
    layers <- accCost(layers, as.matrix(longlat))
    names(layers) <- "cost distance"
  } else {
    layers <- mask(raster::distanceFromPoints(layers, longlat), layers)
    names(layers) <- "minimum distance"
  }
  return(layers)
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
#' @param year Vector of sampling years in longlat (one element per row of longlat). Used to exclude old records with a given probability proportional to time passed since sampling (never excluded only for current year).
#' @param idconf Vector of identification confidence in longlat (one element per row of longlat). Used to exclude uncertain records with a given probability. Can be on any scale where max values are certain (e.g. from 1 - very uncertain to 10 - holotype).
#' @param categorical Vector of layer indices of categorical (as opposed to quantitative) data. If NULL the package will try to find them automatically based on the data.
#' @param thres Threshold of logistic output used for conversion of probabilistic to binary (presence/absence) maps. If 0 this will be the value that maximizes the sum of sensitivity and specificity.
#' @param testpercentage Percentage of records used for testing only. If 0 all records will be used for both training and testing.
#' @param mcp Used for a precautionary approach. If TRUE, all areas predicted as present but outside the minimum convex hull polygon encompassing all occurrence records are converted to absence. Exceptions are cells connected to other areas inside the polygon.
#' @param points If TRUE, force map to include cells with presence records even if suitable habitat was not identified.
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
map.sdm <- function(longlat, layers, error = NULL, year = NULL, idconf = NULL, categorical = NULL, thres = 0, testpercentage = 0, mcp = TRUE, points = FALSE, eval = TRUE, runs = 0, subset = 0){

  raster::rasterOptions(maxmemory = 2e+09)
  origLonglat = longlat

  ##if ensemble is to be done
  if(runs > 0){
    longlat = origLonglat
    
    #if there is spatial error randomly move points within its radius
    if(!is.null(error)){
      for(i in 1:nrow(longlat)){
        #move up to given error (angular movement converted to x and y)
        rndAngle = sample(1:360, 1)
        rndDist = runif(1, 0, error[i])
        longlat[i,1] = longlat[i,1] + rndDist * cos(rndAngle)
        longlat[i,2] = longlat[i,2] + rndDist * sin(rndAngle)
      }
    }
    
    #if there is year
    if(!is.null(year)){
      for(i in 1:nrow(longlat)){
        if(year[i] < sample(min(year):as.integer(substr(Sys.Date(), 1, 4)), 1))
          longlat = longlat[-i,]
      }
    }

    #if there is idconf
    if(!is.null(idconf)){
      for(i in 1:nrow(longlat)){
        if(idconf[i] < sample(1:max(idconf), 1))
          longlat = longlat[-i,]
      }
    }

    if(eval)
      runEval = matrix(NA, nrow = 1, ncol = 7)
    runMap <- rasterize(longlat, layers[[1]], field = 0, background = 0)
    pb <- txtProgressBar(min = 0, max = runs, style = 3)
    totalAUC = 0
    for(i in 1:runs){
      if(subset > 0 && subset < dim(layers)[3]){
        runLayers <- layers[[sample.int(dim(layers)[3], subset)]]
        thisRun <- map.sdm(longlat, runLayers, error = NULL, year = NULL, idconf = NULL, categorical, thres, testpercentage, mcp, points, eval, runs = 0, subset = 0)
      } else {
        thisRun <- map.sdm(longlat, layers, error = NULL, year = NULL, idconf = NULL, categorical, thres, testpercentage, mcp, points, eval, runs = 0, subset = 0)
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

  longlat <- move(longlat, layers)  #move all records falling on NAs

  nPoints = min(1000, sum(!is.na(as.vector(layers[[1]])), na.rm = TRUE)/4)
  bg <- dismo::randomPoints(layers, nPoints)                                ##extract background points

  ##if no categorical variables are given try to figure out which are
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
  
  if(points)
    p <- max(p, map.points(longlat, p, eval = FALSE))

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
#' @param points If TRUE, force map to include cells with presence records even if suitable habitat was not identified.
#' @param eval If TRUE, build a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @details In many cases a species has a very restricted habitat and we generally know where it occurs. In such cases using the distribution of the known habitat patches may be enough to map the species.
#' @return One raster object and, if eval = TRUE, a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
#' @export
map.habitat <- function(longlat, layer, move = TRUE, mcp = FALSE, points = FALSE, eval = TRUE){

  if(points)
    layer <- max(layer, map.points(longlat, layer, eval = FALSE))
  
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

#' Species distributions made easy (multiple species).
#' @description Single step for prediction of multiple species distributions. Output of maps (in pdf format), klms (for Google Earth) and relevant data (in csv format).
#' @param longlat data.frame of taxon names, longitude and latitude or eastness and northness (three columns in this order) of each occurrence record.
#' @param layers If NULL analyses are done with environmental layers read from data files of red.setup(). If a Raster* object as defined by package raster, analyses use these.
#' @param habitat Raster* object as defined by package raster. Habitat extent layer (0/1) used instead of layers if any species is an habitat specialist.
#' @param zone UTM zone if data is in metric units. Used only for correct placement of kmls and countries.
#' @param thin boolean defining if species data should be thinned before modeling (only for SDMs).
#' @param error Vector of spatial error in longlat (one element per row of longlat) in the same unit as longlat. Used to move any point randomly within the error radius.
#' @param move If TRUE, identifies and moves presence records to closest cells with environmental data. Use when spatial error might put records outside such data.
#' @param dem RasterLayer object. It should be a digital elevation model for calculation of elevation limits of the species. If NULL, dem from red.setup() is used if possible, otherwise it will be 0.
#' @param pca Number of pca axes for environmental data reduction. If 0 (default) no pca is made.
#' @param filename Name of output csv file with all results. If NULL it is named "Results_All.csv".
#' @param mapoption Vector of values within options: points, habitat and sdm; each value corresponding to the function to be used for each species (map.points, map.habitat, map.sdm). If a single value, all species will be modelled according to it. If NULL, the function will perform analyses using map.points. Species values must be in same order as latlong.
#' @param testpercentage Percentage of records used for testing only. If 0 all records will be used for both training and testing.
#' @param mintest Minimim number of total occurrence records of any species to set aside a test set. Only used if testpercentage > 0.
#' @param points If TRUE, force map to include cells with presence records even if suitable habitat was not identified.
#' @param runs If <= 0 no ensemble modelling is performed. If > 0, ensemble modelling with n runs is made. For each run, a new random sample of occurrence records (if testpercentage > 0), background points and predictive variables (if subset > 0) are chosen. In the ensemble model, each run is weighted as max(0, (runAUC - 0.5)) ^ 2.
#' @param subset Number of predictive variables to be randomly selected from layers for each run if runs > 0. If <= 0 all layers are used on all runs. Using a small number of layers is usually better than using many variables for rare species, with few occurrence records (Lomba et al. 2010, Breiner et al. 2015).
#' @return Outputs maps in asc, pdf and kml format, plus a file with EOO, AOO and a list of countries where the species is predicted to be present if possible to extract.
#' @references Breiner, F.T., Guisan, A., Bergamini, A., Nobis, M.P. (2015) Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218.
#' @references Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J., Guisan, A. (2010) Overcoming the rare species modelling paradox: a novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143: 2647-2657.
#' @export
map.easy <- function(longlat, layers = NULL, habitat = NULL, zone = NULL, thin = TRUE, error = NULL, move = TRUE, dem = NULL, pca = 0, filename = NULL, mapoption = NULL, testpercentage = 0, mintest = 20, points = FALSE, runs = 0, subset = 0){

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
      if(thin)
        spData <- thin(spData)
      if(testpercentage > 0)
        p <- map.sdm(spData, layers, spError, testpercentage = testpercentage, mcp = TRUE, points = points, runs = runs, subset = subset)
      else
        p <- map.sdm(spData, layers, spError, testpercentage = 0, mcp = TRUE, points = points, runs = runs, subset = subset)
    } else if (mapoption[s] == "habitat"){
      p <- map.habitat(spData, habitat, move, points = points)
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
    points(longlat, pch = 19)
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
#' @param spData spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (either 0/1 or probabilistic values).
#' @details EOO is calculated as the minimum convex polygon covering all known or predicted sites for the species.
#' @return A single value in km2 or a vector with lower confidence limit, consensus and upper confidence limit (probabilities 0.975, 0.5 and 0.025 respectively).
#' @examples data(red.records)
#' data(red.range)
#' eoo(red.records)
#' eoo(red.range)
#' @export
eoo <- function(spData){
  if(class(spData) == "RasterLayer"){
    if(!all(raster::as.matrix(spData) == floor(raster::as.matrix(spData)), na.rm = TRUE)){ #if probabilistic map
      upMap <- reclassify(spData, matrix(c(0,0.025,0,0.025,1,1), ncol = 3, byrow = TRUE))
      consensusMap <- reclassify(spData, matrix(c(0,0.499,0,0.499,1,1), ncol = 3, byrow = TRUE))
      downMap <- reclassify(spData, matrix(c(0,0.975,0,0.975,1,1), ncol = 3, byrow = TRUE))
      area <- c(eoo(downMap), eoo(consensusMap), eoo(upMap))
    } else {
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
  return(round(area))
}

#' Area of Occupancy (AOO).
#' @description Calculates the Area of Occupancy of a species based on either known records or predicted distribution.
#' @param spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (either 0/1 or probabilistic values).
#' @details AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as required by IUCN.
#' @return A single value in km2 or a vector with lower confidence limit, consensus and upper confidence limit (probabilities 0.975, 0.5 and 0.025 respectively).
#' @examples data(red.range)
#' aoo(red.range)
#' @export
aoo <- function(spData){
  if (class(spData) == "RasterLayer"){ #if rasterlayer
    if(raster::maxValue(spData) == 0){  #if no data (empty raster)
      area = 0
    } else if(!all(raster::as.matrix(spData) == floor(raster::as.matrix(spData)), na.rm = TRUE)){ #if probabilistic map
      upMap <- reclassify(spData, matrix(c(0,0.025,0,0.025,1,1), ncol = 3, byrow = TRUE))
      consensusMap <- reclassify(spData, matrix(c(0,0.499,0,0.499,1,1), ncol = 3, byrow = TRUE))
      downMap <- reclassify(spData, matrix(c(0,0.975,0,0.975,1,1), ncol = 3, byrow = TRUE))
      area <- c(aoo(downMap), aoo(consensusMap), aoo(upMap))
    } else {
      if (raster::xmax(spData) <= 180) {  #if longlat data
        if(res(spData)[1] > 0.05){ #if resolution is > 1km use area of cells rounded to nearest 4km
          area = round(cellStats((raster::area(spData) * spData), sum)/4)*4
        } else {
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
        }
      } else { #if square data in meters
        spData[spData < 1] <- NA
        spData <- rasterToPoints(spData)
        spData = floor(spData/2000)
        ncells = nrow(unique(spData))
        area = ncells * 4
      }
    }
  } else if (ncol(spData) == 2){
    if (max(spData) <= 180) {  #if longlat data
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
  return(round(area))
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
  return(round(out))
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
#' @param smooth Smooths the kml lines as per IUCN guidelines. Higher values represent smoother polygons.
#' @param rad radius of circles in degrees if mapoption is "points". It can be the same value for all points or a vector with length equal to number of records in spData representing associated error. The default is about 10km (0.1 degrees) as per IUCN guidelines.
#' @return A kml with polygon or circles around records.
#' @export
kml <- function(spData, zone = NULL, filename, mapoption = "aoo", smooth = 0, rad = 0.1){
  if ((class(spData) == "RasterLayer" && raster::xmax(spData) > 180) || (class(spData) != "RasterLayer" && max(spData) > 180))   ##if need to project to longlat
    spData <- utm2longlat(spData, zone)
  
  if(mapoption == "aoo" && class(spData) == "RasterLayer"){
    spData[spData != 1] <- NA
    spData <- rasterToPolygons(spData, dissolve = TRUE)

    #simplify
    if(smooth > 0){
      trytol <- c(seq(0.001,0.01,0.001),seq(0.02,0.1,0.01),seq(0.2,1,0.1),2:10,seq(20,100,10),seq(200,1000,100),seq(2000,10000,1000),seq(20000,100000,10000),seq(200000,1000000,100000))
      for (i in trytol){
        if(class(try(gSimplify(spData, tol = (1 / i)), silent = TRUE)) != "try-error"){
          spData <- gSimplify(spData, tol = (smooth / (i*10)))
          break
        }
      }

      #cut to coast
      spData <- gIntersection(worldborders, spData)

      #round
      smooth = smooth * 100
      polys = methods::slot(spData@polygons[[1]], "Polygons")

      spline.poly <- function(xy, vertices, k=3, ...) {
        # Assert: xy is an n by 2 matrix with n >= k.

        # Wrap k vertices around each end.
        n <- dim(xy)[1]
        if (k >= 1) {
          data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
        } else {
          data <- xy
        }

        # Spline the x and y coordinates.
        data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
        x <- data.spline$x
        x1 <- data.spline$y
        x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y

        # Retain only the middle part.
        cbind(x1, x2)[k < x & x <= n+k, ]
      }

      spData <- SpatialPolygons(
        Srl = lapply(1:length(polys),
                     function(x){
                       p <- polys[[x]]

                       #applying spline.poly function for smoothing polygon edges
                       px <- methods::slot(polys[[x]], "coords")[,1]
                       py <- methods::slot(polys[[x]], "coords")[,2]
                       bz <- spline.poly(methods::slot(polys[[x]], "coords"),smooth, k=3)
                       bz <- rbind(bz, bz[1,])
                       methods::slot(p, "coords") <- bz

                       # create Polygons object
                       poly <- Polygons(list(p), ID = x)
                     }
        )
      )
      spData <- SpatialPolygonsDataFrame(spData, data=data.frame(ID = 1:length(spData)))
      kmlPolygons(spData, filename, name = filename, col = '#FFFFFFAA', border = "red", lwd = 2)
    } else {
      kmlPolygon(spData, filename, name = filename, col = '#FFFFFFAA', border = "red", lwd = 2)
    }
    
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
    kmlPolygon(poly, filename, name = filename, col = '#FFFFFFAA', border = "red", lwd = 2)
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
    kmlPolygon(poly, filename, name = filename, col = '#FFFFFFAA', border = "red", lwd = 2)
  }
}

#' Red List Index.
#' @description Calculates the Red List Index (RLI) for a group of species.
#' @param spData Either a vector with species assessment categories for a single point in time or a matrix with two points in time in different columns (species x date). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param tree An hclust or phylo object (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param dd bootstrap among all species (FALSE) or Data Deficient species only (TRUE).
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
#' To incorporate Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
#' The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
#' Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index and this approach is here implemented.
#' @return Either a vector (if no two dates are given) or a matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
#' @references Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.
#' @references Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
#' @examples rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","DD","DD"), ncol = 2, byrow = TRUE)
#' colnames(rliData) <- c("2000", "2010")
#' rli(rliData[,1])
#' rli(rliData[,1], boot = TRUE)
#' rli(rliData)
#' rli(rliData, boot = TRUE, dd = TRUE)
#' @export
rli <- function (spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000){

  ##if only one point in time is given
  if(is.null(dim(spData)))
    return(rli.calc(spData, tree, boot, dd, runs))  ##return either 1 or 3 values

  ##if two points in time are given
  ts <- apply(spData, 2, function(x) rli.calc(x, tree, boot = FALSE))
  sl <- (ts[2] - ts[1]) / (as.numeric(colnames(spData))[2] - as.numeric(colnames(spData))[1])
  if(!boot){
    res <- matrix(c(ts, sl), nrow = 1)
    colnames(res) <- c(colnames(spData), "Change/year")
    rownames(res) <- c("Raw")
    return(res)
  } else {
    tr <- apply(spData, 2, function(x) rli.calc(x, tree, boot, dd, runs))
    p = 0
    rndSl = rep(NA, runs)
    for(r in 1:runs){
      rndSl[r] <- rli.calc(spData[,2], tree, boot, dd, runs = 1)[2] - rli.calc(spData[,1], tree, boot, dd, runs = 1)[2]
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
#' @param tree A list of hclust or phylo objects, each corresponding to a tree per group (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @param boot If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.
#' @param dd bootstrap among all species (FALSE) or Data Deficient species only (TRUE).
#' @param runs Number of runs for bootstrapping
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to 5 (Extinct/Extinct in the Wild).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
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
rli.multi <- function (spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000){

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
    if(is.null(tree))
      v <- rli(spData[spData[,1] == groups[g],-1], tree = NULL, boot = boot, dd = dd, runs = runs)
    else
      v <- rli(spData[spData[,1] == groups[g],-1], tree[[g]], boot = boot, dd = dd, runs = runs)
    if(ncol(res) < 13){
      res[g,] <- v
      colnames(res) <- colnames(v)
    } else {
      res[g,1:4] <- v$Values[,1]
      res[g,5:8] <- v$Values[,2]
      res[g,9:12] <- v$Values[,3]
      res[g,13] <- v$P_change
    }
  }
  return(res)
}

#' Prediction of Red List Index.
#' @description Linearly interpolates and extrapolates RLI values to any years.
#' @param rliValue Should be a vector with RLI values and names as the corresponding year numbers.
#' @param from Starting year of the sequence to predict.
#' @param to Ending year of the sequence to predict.
#' @param rliPlot Plots the result
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' @return A matrix with the RLI values and confidence limits.
#' @examples rliValue <- c(4.5, 4.3, 4.4, 4.2, 4.0)
#' names(rliValue) <- c(2000, 2004, 2008, 2011, 2017)
#' rli.predict(rliValue, 1990, 2020)
#' @export
rli.predict <- function(rliValue, from = NA, to = NA, rliPlot = FALSE){
  year = as.numeric(c(names(rliValue)))
  rliTable = data.frame(rliValue, year)
  if(is.na(from))
    from = min(year)
  if(is.na(to))
    to = max(year)
  newYear = data.frame(year = seq(from = from, to = to, by = 1))
  lmOut = predict(lm(rliValue ~ year, data = rliTable), newYear, interval = "confidence", level = 0.95)
  res = lmOut[,c(2,1,3)]
  colnames(res) = c("LowCL", "Fitted RLI", "UpCL")
  rownames(res) = newYear$year
  
  if(rliPlot){
    plot(year, rliValue, xlab="Year", ylab="Fitted RLI", xlim = c(from, to), ylim = c(0,5))
    abline(lm(rliValue ~ year, data = rliTable), col = "red")
    matlines(newYear, lmOut[,2:3], col = "blue", lty = 2)
  }
  
  return(res)
}

#' Sampled Red List Index.
#' @description Calculates accumulation curve of confidence limits in sampled RLI.
#' @param spData A vector with species assessment categories for a single point in time. Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param tree An hclust or phylo object (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
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
rli.sampled <- function (spData, tree = NULL, p = 0.05, runs = 1000){

  nSpp <- length(spData)
  accum <- rep(NA, nSpp)
  for(n in 1:nSpp){               #test with n species from the entire set
    diff = rep(NA, runs)          #try runs times each species
    for(r in 1:runs){             #do r runs for each n species
      rndComm = rep(NA, nSpp)
      rndSpp = sample(nSpp, n)
      rndComm[rndSpp] = spData[rndSpp]
      diff[r] = abs(rli.calc(spData, tree, FALSE, FALSE, runs = 1) - rli.calc(rndComm, tree, FALSE, FALSE, runs = 1))  #calculate absolute difference between true and sampled rli for each run
    }
    accum[n] = quantile(diff, (1-p))
  }
  return(accum)   #returns the accumulation curve of confidence limit of sampled RLI
}

#' Mapping the Red List Index.
#' @description Creates a map for the red list index according to species distribution and threat status.
#' @param spData Either a vector with species assessment categories for a single point in time or a matrix with two points in time in different columns (species x date). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).
#' @param layers Species distributions (0/1), a Raster* object as defined by package raster.
#' @param layers2 Species distributions (0/1) on the second point in time, a Raster* object as defined by package raster. If there are two dates but no layers2, the distributions are assumed to be kept constant in time.
#' @param tree An hclust or phylo object (used when species are weighted by their unique contribution to phylogenetic or functional diversity).
#' @details The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
#' The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to Extinct/Extinct in the Wild (5).
#' Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
#' Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
#' @return A RasterLayer with point values  (if a single date is given) or change per cell (if two dates are given).
#' @references Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.
#' @references Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.
#' @examples sp1 <- raster::raster(matrix(c(1,1,1,0,0,0,0,0,NA), ncol = 3))
#' sp2 <- raster::raster(matrix(c(1,0,0,1,0,0,1,0,NA), ncol = 3))
#' sp3 <- raster::raster(matrix(c(1,0,0,0,0,0,0,0,NA), ncol = 3))
#' sp4 <- raster::raster(matrix(c(0,1,1,1,1,1,1,1,NA), ncol = 3))
#' layers <- raster::stack(sp1, sp2, sp3, sp4)
#' spData <- c("CR","EN","VU","LC")
#' raster::plot(rli.map(spData, layers))
#' @export
rli.map <- function (spData, layers, layers2 = NULL, tree = NULL){
    
    if(!is.null(dim(spData))){             #if to calculate change call this same function twice
        if(is.null(layers2)){
            layers2 <- layers
        }
        map1 <- rli.map(spData[,1], layers = layers, tree = tree)
        map2 <- rli.map(spData[,2], layers = layers2, tree = tree)
        return(map2 - map1)
    }

  #convert rasters to array
  layers = raster::as.array(layers)
  
  #get data for each cell (row by row)
  cells = matrix(NA, (nrow(layers) * ncol(layers)), dim(layers)[3])
  i = 0
  for (r in 1:nrow(layers)){
    for(c in 1:ncol(layers)){
        i = i+1
        cells[i,] = layers[r,c,]
    }
  }

  #RLI of each cell
  rliCells = rep(NA, nrow(cells))
  for (i in 1:nrow(cells)){
      rliNA <- ifelse(cells[i,] == 1, spData, NA) #only consider species present in each cell
      rliCells[i] = rli.calc(rliNA, tree = tree)
  }

  #create RLI map
  rliMap = raster::raster(matrix(rliCells, nrow = nrow(layers), byrow = T))
  return(rliMap)
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
