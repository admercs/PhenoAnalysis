read_phenoanalysis <- function(modeldir,experiments,sitename,varname,xmin,xmax,ymin,ymax,dx,dy,
                               date.start=c(),date.end=c()) {

  library(ncdf)
 
  # derive number of output grid points
  nx <- round((xmax-xmin)/dx)
  ny <- round((ymax-ymin)/dy)

  # define number of experiments
  ne <- length(experiments)
  
  # set NetCDF varname
  if (varname == "SOS") {
    varname.read <- "LAI"
  } else if (varname == "EOS") {
    varname.read <- "LAI"
  } else {
    varname.read <- varname
  }

  # prepare output data
  outdata <- vector(mode="list",length=(ne+1))
 
  for (e in 1:ne) {
    
    # flags
    first <- TRUE

    indir <- paste(modeldir,"/",experiments[e],"/output/",sep="")

    if (!file.exists(indir)) {
      print(paste("Phenoanalysis Output directory :",indir,"does not exist."))
      browser()
    }

    # derive start and end year
    year.start <- as.numeric(format(date.start,"%Y",tz="GMT"))
    year.end <- as.numeric(format(date.end,"%Y",tz="GMT"))
   
    # create file suffix
    suffix <- ".nc"
    
    # define input granularity
    instep <- -1
    
    # create file prefix
    prefix <- paste(sitename,'.analysis.',sep='')
    
    # find files matching that prefix
    pattern <- glob2rx(paste(prefix,"*",suffix,sep=""))
    infiles <- list.files(path=indir,pattern=pattern)
    
    if (length(infiles) == 0) {
      print("No phenoanalysis output files found in input directory.")
      browser()
    }

    for (year in year.start:year.end) {

      syear <- formatC(year,width=4,flag="0",format="d")

      pattern <- glob2rx(paste(prefix,"*",syear,"*",suffix,sep=""))
      infiles <- list.files(path=indir,pattern=pattern)

      if (length(infiles) != 0) {
        
        for (f in 1:length(infiles)) {
          
          ncid <- open.ncdf(paste(indir,infiles[f],sep=''))
          
          # get time
          time <- get.var.ncdf(ncid,"time")
          orig <- att.get.ncdf(ncid,"time","units")[[2]]
          orig1 <- strsplit(orig, " ")[[1]]
          orig2 <- strsplit(orig1[[3]], "-")[[1]]
          orig3 <- strsplit(orig1[[4]], ":")[[1]]
          timestep <- tolower(orig1[1])
          year0 <- as.double(orig2[1])
          month0 <- as.double(orig2[2])
          day0 <- as.double(orig2[3])
          hour0 <- as.double(orig3[1])
          minute0 <- as.double(orig3[2])
          second0 <- as.double(orig3[3])
          
          if (timestep == "years") {
            temptime <- ISOdatetime(year0+round(time),month0,day0,hour0,minute0,second0,tz="GMT")
          } else if (timestep == "months") {
            # FIX: currently only for single year (time: 1..12)
            temptime <- ISOdatetime(year0,month0+round(time),day0,hour0,minute0,second0,tz="GMT")
          } else if (timestep == "days") {
            temptime <- ISOdatetime(year0,month0,day0,hour0,minute0,second0,tz="GMT") + round(time * 8640)*10
          } else if (timestep == "hours") {
            temptime <- ISOdatetime(year0,month0,day0,hour0,minute0,second0,tz="GMT") + round(time * 360)*10
          } else if (timestep == "minutes") {
            temptime <- ISOdatetime(year0,month0,day0,hour0,minute0,second0,tz="GMT") + round(time * 6)*10
          } else if (timestep == "seconds") {
            temptime <- ISOdatetime(year0,month0,day0,hour0,minute0,second0,tz="GMT") + round(time)
          }
          
          # evaluate timestep length
          if (length(temptime) > 1) {
            dt <- as.double(temptime[2]) - as.double(temptime[1])
          } else {
            if (timestep == "years") {
              dt <- NULL
            } else if (timestep == "months") {
              dt <- NULL
            } else if (timestep == "days") {
              dt <- 86400
            } else if (timestep == "hours") {
              dt <- 3600
            } else if (timestep == "minutes") {
              dt <- 60
            } else if (timestep == "seconds") {
              dt <- 1
            }
            
          }

          # derive time bounds to read
          tidx <- which((temptime>=date.start) & (temptime<=date.end))

          if (tidx[1] != -1) {
            t0 <- tidx[1]
            t1 <- tidx[length(tidx)]
            
            # derive / subset geographic extent
            if (first) {
              
              lon <- get.var.ncdf(ncid,"lon")
              lat <- get.var.ncdf(ncid,"lat")
              
              if (length(lon) == 1) {
                # single longitude
                if (nx > 1) stop("Input file has single longitude and output grid has many longitudes")
                if (mean(c(xmin,xmax)) != lon[1]) stop("Input file longitude does not match requested longitude")
                x0 <- 1
                x1 <- 1            
              } else {
                
              }
              
              if (length(lat) == 1) {
                # single latitude
                if (ny > 1) stop("Input file has single latitude and output grid has many latitudes")
                if (mean(c(ymin,ymax)) != lat[1]) stop("Input file latitude does not match requested latitude")
                y0 <- 1
                y1 <- 1
              } else {
                
              }

              
              # get units
              units <- att.get.ncdf(ncid,varname.read,"units")[[2]]
              
              # get long name
              longname <- att.get.ncdf(ncid,varname.read,"long_name")[[2]]
              
              # get version
              version <- att.get.ncdf(ncid,varname.read,"version")[[2]]
              
              # get production date
              proddate <- att.get.ncdf(ncid,varname.read,"prod_date")[[2]]
              
              
            } # read geographic domain

            tempdata <- get.var.ncdf(ncid,varname.read,start=c(x0,y0,t0),count=c(x1-x0+1,y1-y0+1,t1-t0+1))

            if (first) {
              alltime <- temptime[tidx]
              alldata <- tempdata
            } else {
              alltime <- c(alltime,temptime[tidx])
              alldata <- c(alldata,tempdata)
            }

            first <- FALSE
            
          } # time step ok

        } # files

      } # files available for year

    } # years

    outdata[[1]] <- alltime
    outdata[[e+1]] <- alldata
  
    names(outdata)[1] <- "Time"
    names(outdata)[e+1] <- varname.read
    attr(outdata[[e+1]],"units") <- units
    attr(outdata[[e+1]],"longname") <- longname
    attr(outdata[[e+1]],"experiment") <- experiments[e]
    attr(outdata[[e+1]],"version") <- version
    attr(outdata[[e+1]],"proddate") <- proddate
 
  } # experiments
 
  return(outdata)

} # function
