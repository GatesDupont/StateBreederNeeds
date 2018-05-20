  # Gates Dupont      #
  # GLD44@cornell.edu #
  # 2018              #
  
  library(raster)
  library(jsonlite)
  library(ggplot2)
  APIkey = "xxxxxxxxxxxx"
  
  
  #--------FUNCTIONS--------
  
  #----Pull frequencies from internal API----
  freq.API.pull = function(lat, lng, thr){
    url = paste("https://ebird.org/ws2.0/product/geo/freqlist?lat=", 
                lat, "&lng=", lng, "&key=", APIkey, "&m=6", sep = "")
    data = fromJSON(readLines(url, warn=FALSE))
    #data = data$frequencies
    return(data)
  }
  
  #----Pull taxonomy data from internal API-----
  tax.API.pull = function(species){
    species = paste(species, collapse=',')
    url = paste0("https://ebird.org/ws2.0/ref/taxonomy/ebird?fmt=json&species=", species)
    data = fromJSON(readLines(url, warn=FALSE))
    data = data.frame(data$taxonOrder, data$sciName, data$comName, data$speciesCode)
    colnames(data) = c("TaxonomicOrder", "ScientificName", "CommonName", "SpeciesCode")
    return(data)
  }
  
  
  
  #--------LOADING DATA--------
  
  #----Loading my needs list and breeders list----
  needs = read.csv("~/Desktop/MassNeeds_MayJune.csv")
  breeds = read.csv("~/Desktop/BreedersMA2.csv")
  needs.rm = setdiff(needs$ComName,breeds$ComName)
  
  #----Removing species that don't breed in MA----
  needs2 = needs
  for(i in 1:length(needs.rm)){
    needs2 = needs2[needs2$ComName != (needs.rm[i]),]
  }
  needs=needs2; rownames(needs) = NULL
  
  #----Load Massachusetts Polygon----
  us <- raster::getData('GADM', country = 'US', level = 1)
  ma <- us[us$NAME_1 == "Massachusetts",]
  
  
  
  #--------MAKING SAMPLING POINTS--------
  
  #----Making a grid for MA----
  proj4string(ma) <- "+proj=longlat +datum=WGS84 +no_defs"
  aea.proj = paste("+proj=aea +lat_1=", ma@bbox["y",][1], 
                   " +lat_2=", ma@bbox["y",][2], 
                   " +lat_0=", (ma@bbox["y",][2] + ma@bbox["y",][1])/2, 
                   " +lon_0=", (ma@bbox["x",][2] + ma@bbox["x",][1])/2, 
                   " +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m", sep = "")
  ma = spTransform(ma, CRS(aea.proj))
  
  # Expanding grid
  ma.ext = extent(-168115.1, 168978.2, -90110.46, 91864.72)
  ma.ext = as(ma.ext, "SpatialPolygons")
  sp::proj4string(ma.ext) = aea.proj
  
  grid = makegrid(ma.ext, CRS(aea.proj), cellsize = (20*1000)) #km to m -- 20km grid
  grid <- SpatialPoints(grid, proj4string = CRS(proj4string(ma)))
  #grid <- grid[ma, ]
  API.points = spTransform(grid, CRS("+proj=longlat +datum=WGS84 +no_defs"))
  #At this point, I have my needs, and my API sampling points.
  
  
  
  #--------PULLING ALL SPECIES FREQUENCIES FROM API--------
  
  #----Pulling frequency for all species----
  point.data = list()
  for(i in (1:length(API.points@coords[,"x1"]))){
    point.data[[i]] = freq.API.pull(lat = API.points@coords[,"x2"][i],
                                    lng = API.points@coords[,"x1"][i])
  }
  
  #----Looping through to get only max freq for May-June (5 week period) records----
  good.data = list()
  for(j in 1:length(point.data)){
    allData = as.data.frame(point.data[[j]])
    species = c()
    max.freq = c()
    for(i in unique(point.data[[j]]$speciesCode)){
      species = c(species, i)
      to.no.100 =  unlist(allData[allData$speciesCode == i,][2])
      to.no.100 = to.no.100[23:27]
      no.100 = to.no.100[to.no.100 != 100]
      max.freq = c(max.freq, max(no.100))
    }
    good.data[[j]] = data.frame(species, max.freq)
  }
  
  #----Working towards selecting only my needs----
  
  options(scipen=999) # Turning off scientific notation for matrix
  
  final.data.all.sp = list()
  for(j in 1:length(good.data)){
    mydata = as.data.frame(good.data[[j]])
    
    # Cleaning filter-level frequencies of 100 
    m = as.matrix(mydata)
    m[m == -Inf & is.na(m)] = 0
    rawOutput0 = as.data.frame(m)
    
    # Preparing to combine taxnomony
    rawOutput = rawOutput0[order(rawOutput0$species),] 
    erdTax0 = tax.API.pull(rawOutput$species)
    erdTax = erdTax0[order(erdTax0$SpeciesCode),]
    
    # Final output, ordered, with combined taxonomy
    output = cbind(erdTax, rawOutput)
    output$species = NULL # Dropping duplicate column originally used to check row alignment
    final.data.all.sp[[j]] = output[order(output$TaxonomicOrder),]
  }
  # I now have frequencies for all MA species in May-June, with Common Names
  
  
  
  #--------SELECTING MY STATE NEEEDS BREEDERS--------
  
  #----Selecting the birds that I do not want included (already have OR non-breeders)----
  species.rm = c()
  for(j in 1:length(final.data.all.sp)){
    species.rm = c(species.rm, as.character(final.data.all.sp[[j]]$CommonName))
  }
  species.rm = unique(species.rm)
  species.rm = setdiff(species.rm, needs$ComName)
  
  #----Removing all other birds at each point----
  point.needs = list()
  for(j in 1:length(final.data.all.sp)){
    df = final.data.all.sp[[j]]
    for(i in 1:length(species.rm)){
      df = df[df$CommonName != species.rm[i],]
    }
    point.needs[[j]] = df
  }
  
  
  
  #--------SUMMING FREQUENCIES--------
  
  #----Summing frequency of needed species at each point----
  counter = c()
  for(j in 1:length(point.needs)){
    counter = c(counter, mean(as.numeric(point.needs[[j]]$max.freq)))
  }
  counter = counter/100
  
  output = data.frame(x = API.points@coords[,2], y = API.points@coords[,2], z = counter)
  
  #----FINAL PLOT----
  plot(rasterFromXYZ(data.frame(grid@coords, output$z)), 
       col=colorRampPalette(c("yellow", "red"))(100),
       frame.plot=F, axes=F, box=F, legend=T, zlim=c(0,.08))
  lines(ma)
  
  # Save the plot and use a Layer Mask of Massachusetts in Photoshop
