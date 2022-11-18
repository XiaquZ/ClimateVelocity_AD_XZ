##Load packages
library("sp", lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/") #error when loading raster
library("iterators",lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/") # Otherwise error when loading doParallel
library("codetools",lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/") # Otherwise error when loading raster
library("sf",lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")# Otherwise error when loading gdalUtilities
library("dplyr", lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('raster',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('igraph',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('Matrix',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('rgdal',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")#Otherwise rgeos and gdistance package cannot be loaded
library('gdistance',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('rgeos',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('foreach',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('doParallel',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('snow',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library("snowfall", lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
library('data.table',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")


##Set climate bin width.
bin.width <- 0.25   
## +- 0.25 degrees C 
## total bin diameter = 0.5 degrees C 

#cost.penalty <- 2 ## Two penalty units per degree C dissimilarity from temperature of interest

## Get current and future climate and resistance raster
pre1 <- raster("/lustre1/scratch/348/vsc34871/input/ForestMicroMAT2020_be_EPSG3035.tif"); names(pre1) <- 'pre1' ## Mean annual temperature (1995); units are degrees C 
fut30 <- raster("/lustre1/scratch/348/vsc34871/input/ForestMicroMAT_be_EPSG3035_ssp126.tif"); names(fut30) <- 'fut30' ## Mean annual temperature (2085); units are degrees C 

####Create a stack for only pre and fut data, so 2 layers.
the.stack_ssp126 <-  stack(pre1,fut30)
plot(the.stack_ssp126)

####Create initial climate data frame.####
clim_ssp126 <- the.stack_ssp126 %>%
  getValues() %>%
  data.frame(cid = 1:ncell(the.stack_ssp126)) %>%
  na.omit()

#Get the coordination for each pixel.
clim_ssp126[,c("x","y")] <- xyFromCell(the.stack_ssp126, clim_ssp126$cid) 

##Round the pre1 and fut30 to keep 1 decimal.
colnames(clim_ssp126)[1] <- "pre1" #Change the column name of data frame
colnames(clim_ssp126)[2] <- "fut30"
clim_ssp126$pre1 <- round(clim_ssp126$pre1,1) 
clim_ssp126$fut30 <- round(clim_ssp126$fut30, 1)
head(clim_ssp126)
clim_ssp126 <- as.data.table(clim_ssp126)

#### Make a foreach parallel processing around the code you want to run for each temperature value####

##Some input indices.
n=1 #The number of climate variable, in this case, only one climate variable (mean annual temp)
geoTol = 180000 #Searching radius for the analog pixel. (2 km yr-1)
tdiff = 90 #Time interval between present and future.

##Prepare data
dat <- na.omit(data.table(clim_ssp126))
# matrix with the future climatic values for all cells
fut <- dat[, seq(2, (2*n), by = 2), with=FALSE] #Get future predicted temperature.

## set things up for parallel processing.
cores = detectCores()
ncores = cores[1]-1
cuts <- cut(1:nrow(dat), ncores, labels = FALSE)
cl <- makeCluster(ncores)
registerDoParallel(cl)

result <- foreach(x = 1:ncores, .combine = rbind, 
                  .multicombine = TRUE) %dopar% { #Change the raster package to terra?
                    library("sp", lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/") #error when loading raster
                    library("iterators",lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/") # Otherwise error when loading doParallel
                    library("codetools",lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/") # Otherwise error when loading raster
                    library('raster',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('igraph',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('Matrix',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('rgdal',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('gdistance',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('rgeos',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('foreach',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('doParallel',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('snow',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library("snowfall", lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    library('data.table',lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
                    
                    a <- x
                    Dat <- dat[cuts == a,]
                    
                    resu <- data.table(focal = Dat$cid, target = as.integer(NA), climDis = as.double(NA), geoDis_m = as.double(NA), ang = as.double(NA), vel_km_yr = as.double(NA))
                    i <- 0
                    while(i <= nrow(Dat)){
                      i <- i+1
                      
                      ## for each focal cell subset neighbor cells # Added by AD
                      pres_XY <- as.numeric(Dat[i,(ncol(Dat)-1):ncol(Dat)])  # Added by AD
                      dif_XY <- data.table(sweep(dat[,(ncol(dat)-1):ncol(dat)], 2, pres_XY, "-"))  # Added by AD
                      # Keep only cells with X,Y differences < than k threshold # Added by AD
                      upperXY = colnames(dif_XY)  # Added by AD
                      l_XY <- lapply(upperXY, function(x) call("<", call("abs", as.name(x)), 10000))  # Added by AD 
                      #Comment from Xiaqu: If I understood correctly, the "0.1" here is the threshold of the search radius, and the 0.1 is around 11 km.
                      
                      ii_XY = Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1=c1, .c2=c2)), l_XY)  # Added by AD
                      #Xiaqu: The line above is a bit difficult to understand 
                      
                      anacid_XY <- dat$cid[dif_XY[eval(ii_XY), which=TRUE]]  # cids analogue cells  # Added by AD
                      
                      # for each focal cell subset target cell analogues (within ClimTol)
                      pres <- as.numeric(Dat[i, seq(1, (2*n), by = 2), with=FALSE])
                      fut <- dat[anacid_XY, seq(2, (2*n), by = 2), with=FALSE] # Added by AD
                      dif <- data.table(sweep(fut, 2, pres, "-"))
                      
                      # Identify future analogue cells
                      # Ohlemuller et al 2006 / Hamann et al 2015
                      upper = colnames(dif)
                      l <- lapply(upper, function(x) call("<", call("abs", as.name(x)), bin.width[grep(x, colnames(dif))]))
                      ii = Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1=c1, .c2=c2)), l)
                      anacid <- dat$cid[dif[eval(ii), which=TRUE]]  # cids analogue cells
                      
                      
                      # LOCATE CLOSEST ANALOGUE
                      if(length(anacid)>0){
                        # check which of those are within distance and get the analogue at minimum distance
                        x <- dat$x[dat$cid %in% anacid]
                        match <- as.data.frame(cbind(Dat$x[i],Dat$y[i], dat$x[dat$cid %in% anacid], dat$y[dat$cid %in% anacid]))
                        colnames(match)[1] <- "from.x" #Change the column name of data frame
                        colnames(match)[2] <- "from.y"
                        colnames(match)[3] <- "to.x" #Change the column name of data frame
                        colnames(match)[4] <- "to.y"
                        match$d <- sqrt((match$from.x-match$to.x)^2 + (match$from.y-match$to.y)^2) # in x/y units
                        
                        
                        an <- anacid[match$d < geoTol]       # cids analogue cells within search radius
                        dis <- match$d[match$d < geoTol]     # distance to candidate analogues
                        if (length(an) > 0){
                          resu[i, target := an[which.min(dis)]]   # cid of geographically closest climate analogue
                          resu[i, climDis := mean(as.numeric(dif[which(anacid == resu[i, target]),]))]  # mean clim difference for the closest analogue
                          resu[i, geoDis_m := round(min(dis), digits = 1)]
                          #resu[i, ang := geosphere::bearing(Dat[i, c("x","y")], dat[cid == resu[i, target], c("x","y")])]
                          resu[i, vel_km_yr := round((resu$geoDis_m[i]/1000)/tdiff, digits = 1)] #velocity in km yr-1.
                        }}
                    }
                    return(resu)
                  }
stopCluster(cl)

#Plot the results of velocity.
Eucli_VoCC <- raster(pre1)
Eucli_VoCC[result$focal] <-  result$vel_km_yr
library("terra", lib.loc = "/vsc-hard-mounts/leuven-data/348/vsc34871/Rlib/4.0.2-foss-2018a/")
Eucli_VoCC <- rast(Eucli_VoCC)
Eucli_VoCC <- round(Eucli_VoCC, digits = 1)
plot(Eucli_VoCC)

#write.csv(result, paste('/vsc-hard-mounts/leuven-data/348/vsc34871/VoCC/Output/', 'FVoCC_EuclideanDist_BE_Forest_SSP126', '.csv', sep=''), row.names=F)
writeRaster(Eucli_VoCC, filename = "/vsc-hard-mounts/leuven-data/348/vsc34871/VoCC/Output/FVoCC_EuclideanDist_BE_Forest_SSP126.tif", overwrite=TRUE)


