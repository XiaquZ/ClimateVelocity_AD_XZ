##Load packages
library(raster)
library(gdistance)
library(parallel)
library(doParallel)
library(dplyr)
library(rgdal)
library(data.table)

##Set climate bin width.
climTol <- 0.25   
## +- 0.25 degrees C 
## total bin diameter = 0.5 degrees C 

#cost.penalty <- 2 ## For LCP method: Two penalty units per degree C dissimilarity from temperature of interest

## Get current and future climate and resistance raster
## Get current and future climate and resistance raster
setwd('D:/VoCC_data')
pre1 <- raster('./micro2020Forest_BE_WGS84.tif'); names(pre1) <- 'pre1' ## Mean annual temperature (1995); units are degrees C 
fut30 <- raster('./InputSSP126/micro2100_ssp126_Forest_BE_WGS84.tif'); names(fut30) <- 'fut30' ## Mean annual temperature (2085); units are degrees C 
#resistance.mask <- raster('C:/Users/u0142858/OneDrive - KU Leuven/KUL/PhD/My Project/WP1_Mapping_CB/OutputsVoCC/OutputTIF/OutputSSP126/resistance.mask_EU_CHELSAbio1_ssp126.tif')

####Create a stack for only pre and fut data, so 2 layers.
the.stack_ssp126 <-  stack(pre1,fut30)
plot(the.stack_ssp126)

####Create initial climate data frame.####
clim <- the.stack_ssp126 %>%
  getValues() %>%
  data.frame(cid = 1:ncell(the.stack_ssp126)) %>%
  na.omit()

#Get the coordination for each pixel.
clim[,c("x","y")] <- xyFromCell(the.stack_ssp126, clim$cid) 

##Round the pre1 and fut30 to keep 1 decimal.
colnames(clim)[1] <- "pre1" #Change the column name of data frame
colnames(clim)[2] <- "fut30"
clim$pre1 <- round(clim$pre1,1) 
clim$fut30 <- round(clim$fut30, 1)
head(clim)
clim <- as.data.table(clim)

##Some input indices.
n=1 #The number of climate variable, in this case, only one climate variable (mean annual temp)
geoTol = 180 #Searching radius for the analog pixel. (2 km yr-1)
tdiff = 90 #Time interval between present and future.

##Prepare data
dat <- na.omit(data.table(clim))
# matrix with the future climatic values for all cells
#fut <- dat[, seq(2, (2*n), by = 2), with=FALSE] #Get future predicted temperature.


##################################################################################
#' Distance-based velocity based on geographically closest climate analogue
#' Initial dVoCC code modified by Aggeliki Doxa for Xiaqu Zhou analysis
##################################################################################

dVoCC_AD_GreatCircle_Single <- function(clim, n, tdiff, method = "Single", climTol, geoTol, distfun = "GreatCircle", trans = NA, lonlat = TRUE){
  
  
  dat <- na.omit(data.table(clim))
  
  # set things up for parallel processing
  cores = detectCores()
  ncores = cores[1]-1
  cuts <- cut(1:nrow(dat), ncores, labels = FALSE)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  result <- foreach(x = 1:ncores, .combine = rbind, fill = TRUE, .packages = c('raster', 'gdistance', 'sp', 'rgeos', 'geosphere', 'rgdal', 'data.table'), .multicombine = TRUE) %dopar% {
    
    a <- x
    Dat <- dat[cuts == a,]
    
    
    resu <- data.table(focal = Dat$cid, target = as.integer(NA), climDis = as.double(NA), geoDis = as.double(NA), vel = as.double(NA))
    i <- 0
    while(i <= nrow(Dat)){
      i <- i+1
      ## for each focal cell subset neighbor cells # Added by AD
      pres_XY <- as.numeric(Dat[i,(ncol(Dat)-1):ncol(Dat)])  # Added by AD
      dif_XY <- data.table(sweep(dat[,(ncol(dat)-1):ncol(dat)], 2, pres_XY, "-"))  # Added by AD
      # Keep only cells with X,Y differences < than k threshold # Added by AD
      upperXY = colnames(dif_XY)  # Added by AD
      l_XY <- lapply(upperXY, function(x) call("<", call("abs", as.name(x)), 0.02))  # Added by AD 
      #Comment from Xiaqu: If I understood correctly, the "0.1" here is the threshold of the search radius, and the 0.1 is around 11 km.
      
      ii_XY = Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1=c1, .c2=c2)), l_XY)  # Added by AD
      #Xiaqu: The line above is a bit difficult to understand, but it also appeared in the code of Garcia et al, 2019. 
      #Just I did not think about using this scripts to subset neighbor cells...
      
      anacid_XY <- dat$cid[dif_XY[eval(ii_XY), which=TRUE]]  # cids analogue cells  # Added by AD
      
      # for each focal cell subset target cell analogues (within ClimTol) WITHIN THE SPECIFIED DISTANCE # Added by AD
      pres <- as.numeric(Dat[i, seq(1, (2*n), by = 2), with=FALSE])
      fut <- dat[anacid_XY, seq(2, (2*n), by = 2), with=FALSE] # Added by AD
      dif <- data.table(sweep(fut, 2, pres, "-"))
      
      # Identify future analogue cells
      #      if(method == "Single"){     # Ohlemuller et al 2006 / Hamann et al 2015
      upper = colnames(dif)
      l <- lapply(upper, function(x) call("<", call("abs", as.name(x)), climTol[grep(x, colnames(dif))]))
      ii = Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1=c1, .c2=c2)), l)
      anacid <- dat$cid[dif[eval(ii), which=TRUE]]  # cids analogue cells
      #     }
      
      # LOCATE CLOSEST ANALOGUE
      if(length(anacid)>0){ 
        d <- (geosphere::distHaversine(cbind(Dat$x[i],Dat$y[i]), cbind(dat$x[dat$cid %in% anacid], dat$y[dat$cid %in% anacid])))/1000   # in km
        
        an <- anacid[d < geoTol]       # cids analogue cells within search radius
        dis <- d[d < geoTol]      # distance to candidate analogues
        if (length(an) > 0){
          resu[i, target := an[which.min(dis)]]   # cid of geographically closest climate analogue
          resu[i, climDis := mean(as.numeric(dif[which(anacid == resu[i, target]),]))]  # mean clim difference for the closest analogue
          resu[i, geoDis_m := round(min(dis), digits = 1)] # Xiaqu: Round the digital number to reduce file size.
          # resu[i, ang := geosphere::bearing(Dat[i, c("x","y")], dat[cid == resu[i, target], c("x","y")])]
          resu[i, vel := resu$geoDis_m[i]/tdiff]  # Xiaqu: Round the digital number to reduce file size.
        }
      }
      
      
    }
    # setTxtProgressBar(pb, i) #added by AD
    return(resu)
  }
  
  # close(pb) #added by AD
  stopCluster(cl)
  
  return(result)
  
}
GC_VoCC <- dVoCC_AD_GreatCircle_Single(clim, n=1, tdiff = 90, method = "Single", climTol = 0.25, geoTol = 180, distfun = "GreatCircle",  trans = NA, lonlat = TRUE)
head(GC_VoCC)
r1 <- raster(pre1)
r1[GC_VoCC$focal] <- GC_VoCC$vel
plot(r1)
