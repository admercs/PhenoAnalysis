PRO read_modis_composite, datadir, year, has_obs, gridmask, $
                          lonmin, lonmax, latmin, latmax, dlon, dlat, pollon, pollat, $
                          obsvar, obsvarerror, obsvariance, obsvarabsdev, obsdays,$
                          obsvarname, obsvarunits, obsvarlongname

  has_obs = 0

  IF (obsvarname EQ 'FPAR') OR (obsvarname EQ 'LAI') THEN BEGIN

     indir = datadir[0] + '/MOD15A2.005_composite/'

     ncvarname = obsvarname

     syear = string(year,format='(I4)')
     infile = indir +'MOD15A2.global_0.5.'+syear+'.nc'
     temp = file_search(infile)
     IF temp EQ '' THEN BEGIN
        print,'MODIS composite file ',infile,' not found. stopping.'
        stop
     ENDIF ELSE has_obs = 1

;; define output dimensions
     nlon = round((lonmax - lonmin)/dlon) ; number of longitudes in output
     nlat = round((latmax - latmin)/dlat) ; number of latitudes in output

;; inquire geographic domain of input file
     ncid = ncdf_open(infile)

     varid = ncdf_varid(ncid,'lon')
     IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlon')
     ncdf_varget,ncid,varid,lons_modis

     varid = ncdf_varid(ncid,'lat')
     IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlat')
     ncdf_varget,ncid,varid,lats_modis

     varid = ncdf_varid(ncid,'time')
     ncdf_varget,ncid,varid,time_modis

     ntime = n_elements(time_modis)

     dlon_modis = abs(lons_modis[1]-lons_modis[0])
     dlat_modis = abs(lats_modis[1]-lats_modis[0])

     latmin_modis = min(lats_modis) - 0.5d0*dlat_modis
     lonmin_modis = min(lons_modis) - 0.5d0*dlon_modis
     latmax_modis = max(lats_modis) + 0.5d0*dlat_modis
     lonmax_modis = max(lons_modis) + 0.5d0*dlon_modis

     nlon_modis = n_elements(lons_modis)
     nlat_modis = n_elements(lats_modis)

     xmin_modis = -1L
     xmax_modis = -1L
     ymin_modis = -1L
     ymax_modis = -1L

;; define precision for getting rid of rounding errors below
     pr = 1.d-6

     FOR i=nlon_modis-1L,0L,-1L DO IF (round((lons_modis[i]-0.5d0*dlon_modis)/pr,/L64) LE round(lonmin/pr,/L64)) $
        AND (xmin_modis EQ -1L) THEN xmin_modis = i
     FOR i=nlat_modis-1L,0L,-1L DO IF (round((lats_modis[i]-0.5d0*dlat_modis)/pr,/L64) LE round(latmin/pr,/L64)) $
        AND (ymin_modis EQ -1L) THEN ymin_modis = i
     FOR i=0L,nlon_modis-1L,1L DO IF (round((lons_modis[i]+0.5d0*dlon_modis)/pr,/L64) GE round(lonmax/pr,/L64)) $
        AND (xmax_modis EQ -1L) THEN xmax_modis = i
     FOR i=0L,nlat_modis-1L,1L DO IF (round((lats_modis[i]+0.5d0*dlat_modis)/pr,/L64) GE round(latmax/pr,/L64)) $
        AND (ymax_modis EQ -1L) THEN ymax_modis = i

     IF ((xmin_modis EQ (-1L)) OR (ymin_modis EQ (-1L)) OR (xmax_modis EQ (-1L)) OR (ymax_modis EQ (-1L))) THEN BEGIN
        print,'MODIS composite dataset does not cover requested area. stopping.'
        print,'We need: ',latmin,latmax,lonmin,lonmax
        print,'We get:  ',latmin_modis,latmax_modis,lonmin_modis,lonmax_modis
        stop
     ENDIF

     nx_modis = xmax_modis - xmin_modis + 1L
     ny_modis = ymax_modis - ymin_modis + 1L

     ;; define output variables
     obsvar = fltarr(nlon,nlat,ntime)
     obsvarerror = fltarr(2,nlon,nlat,ntime)
     obsvariance = fltarr(nlon,nlat,ntime)
     obsvarabsdev = fltarr(nlon,nlat,ntime)
     obsdays = time_modis
     
     ncdf_varget,ncid, ncvarname, data, offset = [xmin_modis,ymin_modis,0], count = [nx_modis,ny_modis,ntime]      
     ncdf_attget,ncid, ncvarname,'_FillValue',nodata
     ncdf_attget,ncid, ncvarname,'long_name',tmp
     obsvarlongname=strtrim(string(tmp),2)
     ncdf_attget,ncid, ncvarname,'units',tmp
     obsvarunits=strtrim(string(tmp),2)

     badvalues = where(data EQ nodata,badcount)
     IF badcount GT 0 THEN data[badvalues] = !values.f_nan

     ncdf_varget,ncid, strtrim(ncvarname,2)+'_MAX_MIN', error, $
                 offset = [0,xmin_modis,ymin_modis,0], count = [2,nx_modis,ny_modis,ntime]  
     ncdf_attget,ncid, ncvarname,'_FillValue',nodata
     
     badvalues = where(error EQ nodata,badcount)
     IF badcount GT 0 THEN error[badvalues] = !values.f_nan

;; scale pft's to output dimension if needed
     IF (nx_modis ne nlon) OR (ny_modis ne nlat) THEN BEGIN
        
        IF (nx_modis EQ (nx_modis/nlon)*nlon) AND (ny_modis EQ (ny_modis/nlat)*nlat) THEN BEGIN
           ;; simple rebinning
           obsvar[*] = rebin(data,nlon,nlat,ntime)
           obsvarerror[*] = rebin(error,2,nlon,nlat,ntime)
        ENDIF ELSE BEGIN
           IF (nlon*nx_modis GT 20000) OR (nlat*ny_modis GT 20000) THEN BEGIN
              ;; NN rebinning
              FOR t=0,ntime-1 DO BEGIN
                 obsvar[*,*,t] = congrid(reform(data[*,*,t]),nlon,nlat)
                 FOR n=0,1 DO BEGIN
                    obsvarerror[n,*,*,t] = congrid(reform(error[n,*,*,t]),nlon,nlat)
                 ENDFOR
              ENDFOR
           ENDIF ELSE BEGIN
              ;; find common denominator, then do magnifying
              FOR t=0,ntime-1 DO BEGIN
                 obsvar[*,*,t] = rebin(congrid(reform(data[*,*,t]),nlon*nx_modis,nlat*ny_modis),nlon,nlat)
                 FOR n=0,1 DO BEGIN
                    obsvarerror[n,*,*,t] = rebin(congrid(reform(error[n,*,*,t]),nlon*nx_modis,nlat*ny_modis),$
                                                 nlon,nlat)
                 ENDFOR
              ENDFOR
           ENDELSE
        ENDELSE
        
     ENDIF ELSE BEGIN
        
        obsvar[*] = data
        obsvarerror[*] = error
        
     ENDELSE
     
     ncdf_close,ncid

  ENDIF


END
