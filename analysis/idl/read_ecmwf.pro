;+
; :Description:
;
;   This code read daily ECMWF ERA-40, ERA-INTERIM or Operational
;   data on the chosen output grid. If the elevation of the output
;   grid is provided, then the ECMWF data is downscaled with a given
;   lapse rate (only makes sense for near-surface temperature fields!)
;   or by use of exponential downscaling with a scale height (for
;   instance for integrated water vapor etc.)
;
; :Categories:
;   File reading
;
; :Params:
;   ecmwfpath: in, required, type=string
;      Full directory to ECMWF data
;   varname: in, required, type=string
;      Variable name = NetCDF variable name = File name:
;      varname.YYYYMM.nc
;   lon: in, required, type="float or fltarr"
;      longitude needed output grid (can be irregular array, regular array or
;      vector) [degrees E]
;   lat: in, required, type="float or fltarr"
;      latitude needed output grid (can be irregular array, regular array or
;      vector) [degrees N]
;   year: in, required, type=integer
;      year [YYYY]
;   month: in, required, type=integer
;      month [MM, 1-12]
;   day: in, required, type=integer
;      day of month [DD, 1-31]
;
; :Keywords:
;   dlon: in, optional, type=float
;      longitude grid distance for regular grids [degrees]
;   dlat: in, optional, type=float
;      latitude grid distance for regular grids [degrees]
;   elevation: in, optional, type="float or integer"
;      elevation [m a.s.l] of output grid (needs to have the same dimension as
;      lon/lat). If keyword is set, then either lapse_rate or
;      scale_height needs to be set also
;   lapse_rate: in, optional, type=float
;      lapse rate for linear downscaling [m]
;   scale_height: in, optional, type=float
;      scale height for exponential downscaling [m]
;   missing: in, optional, type=float
;      missing data value on output
;
; :Returns: 
;   read (and downscaled) ECMWF field in desired output grid
; 
; :Bugs:
;   For efficiency reasons the regridding information from the ECMWF
;   to the output grid is only calculated once and stored in a common
;   data space. If different grids are used in a single executable
;   where this code is called from the common data space needs to be
;   extended for handling different output grids simultaneously.
;
;   The dateline is not handled correctly. For bilinear interpolation
;   across the dateline the respective region across the dateline
;   should be read. Right now, the dateline is simply curtailed at
;   -179.5 and 179.5, respectively.
;
; :Todo:
;   Adapt for other types of ECMWF data (monthly forecasts, ensemble
;   data etc.)
;
; :Requires:
;   IDL 7.1
;
; :Examples:
;
; :History:
;   2012/03/29 : included exponential downscaling in addition to
;   linear downscaling
;
;   2012/05/08 : corrected bug with undefined scale_fit and lapse_fit
;
;   2012/06/25 : always use default lapse rate and scale height
;
;   2013/12/27 : provide full directory to ECMWF data, this script
;   will find out whether the data is stored in daily, monthly or
;   yearly files
;   2014/01/06 : replaced congrid downscaling with bilinear /
;   trilinear interpolation for regular grids
;
; :Author:
;   Reto Stockli (Blue Marble Research)
;    
; :Copyright:
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;
;   You should have received a copy of the GNU General Public License
;   along with this program.  If not, see http://www.gnu.org/licenses/.
;
;-
FUNCTION read_ecmwf,ecmwfpath, varname, lon, lat, year, month, day, $
                    dlon=dlon, dlat=dlat, elevation=elevation, lapse_rate=lapse_rate, $
                    scale_height=scale_height, missing=missing

  COMMON ECMWF_DATA, xmin_ecmwf, xmax_ecmwf, ymin_ecmwf, ymax_ecmwf, $
     nx_ecmwf, ny_ecmwf, z_ecmwf, x_ecmwf, y_ecmwf, z_ecmwf_grid, $
     dlon_ecmwf, dlat_ecmwf

  ;; parameters
  ntday = 4L ;; number of ECMWF analysis time steps per day
  nt = ntday ;; number of ECMWF time steps to read

  ;; check whether we need to produce a regular grid or not
  IF keyword_set(dlon) AND keyword_set(dlat) THEN regular=1 ELSE regular=0

  ;; check for keywords
  IF n_elements(missing) EQ 0 THEN missing = -9999.0

  syear = string(year,format='(I4.4)')
  smonth = string(month,format='(I2.2)')
  sday = string(day,format='(I2.2)')

  ;; define output dimensions
  dims = size(lon,/dimensions)
  IF n_elements(dims) GT 1 THEN ll2d=1 ELSE ll2d=0

  IF regular THEN BEGIN
     IF ll2d THEN BEGIN
        nlon = (size(lon,/dimensions))[0]
        nlat = (size(lon,/dimensions))[1]
     ENDIF ELSE BEGIN
        nlon = (size([lon],/dimensions))[0]
        nlat = (size([lat],/dimensions))[0]
     ENDELSE

     lonmin = min(lon)-0.5d0*dlon
     lonmax = max(lon)+0.5d0*dlon
     latmin = min(lat)-0.5d0*dlat
     latmax = max(lat)+0.5d0*dlat
  ENDIF ELSE BEGIN
     IF ll2d THEN BEGIN
        nlon = (size(lon,/dimensions))[0]
        nlat = (size(lon,/dimensions))[1]
     ENDIF ELSE BEGIN
        nlon = (size([lon],/dimensions))[0]
        nlat = 1
     ENDELSE

     lonmin = min(lon)
     lonmax = max(lon)
     latmin = min(lat)
     latmax = max(lat)
  ENDELSE
  
  ;; build output array
  IF regular OR ll2d THEN BEGIN
     outdata = make_array(nlon,nlat,nt,/float,value=missing)
  ENDIF ELSE BEGIN
     outdata = make_array(nlon,nt,/float,value=missing)        
  ENDELSE

  ;; inquire geographic sub-domain of ecmwf file and read elevation at
  ;; first read (all successive reads in this process will use the same domain!!!)

  IF n_elements(xmin_ecmwf) EQ 0L THEN BEGIN

     infile = file_search(ecmwfpath+'Z.nc')
     IF infile EQ '' THEN BEGIN
        print,'ecmwf file '+ecmwfpath+'Z.nc'+' containing ECMWF orography (geopotential height) not found. stopping.'
        stop
     ENDIF
     ncid = ncdf_open(infile[0])

     varid = ncdf_varid(ncid,'lon')
     IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlon')
     ncdf_varget,ncid,varid,lons_ecmwf

     varid = ncdf_varid(ncid,'lat')
     IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlat')
     ncdf_varget,ncid,varid,lats_ecmwf

     dlon_ecmwf = lons_ecmwf[1]-lons_ecmwf[0]
     dlat_ecmwf = lats_ecmwf[1]-lats_ecmwf[0]

     latmin_ecmwf = min(lats_ecmwf) - 0.5d0*dlat_ecmwf
     lonmin_ecmwf = min(lons_ecmwf) - 0.5d0*dlon_ecmwf
     latmax_ecmwf = max(lats_ecmwf) + 0.5d0*dlat_ecmwf
     lonmax_ecmwf = max(lons_ecmwf) + 0.5d0*dlon_ecmwf

     nlon_ecmwf = n_elements(lons_ecmwf)
     nlat_ecmwf = n_elements(lats_ecmwf)

     xmin_ecmwf = -1L
     xmax_ecmwf = -1L
     ymin_ecmwf = -1L
     ymax_ecmwf = -1L

     ;; define precision for getting rid of rounding errors below
     pr = 1.d-6

     IF ~regular THEN BEGIN
        FOR i=nlon_ecmwf-1L,0L,-1L DO IF (round(lons_ecmwf[i]/pr,/L64) LE round(lonmin/pr,/L64)) $
           AND (xmin_ecmwf EQ -1L) THEN xmin_ecmwf = i
        FOR i=nlat_ecmwf-1L,0L,-1L DO IF (round(lats_ecmwf[i]/pr,/L64) LE round(latmin/pr,/L64)) $
           AND (ymin_ecmwf EQ -1L) THEN ymin_ecmwf = i
        FOR i=0L,nlon_ecmwf-1L,1L DO IF (round(lons_ecmwf[i]/pr,/L64) GE round(lonmax/pr,/L64)) $
           AND (xmax_ecmwf EQ -1L) THEN xmax_ecmwf = i
        FOR i=0L,nlat_ecmwf-1L,1L DO IF (round(lats_ecmwf[i]/pr,/L64) GE round(latmax/pr,/L64)) $
           AND (ymax_ecmwf EQ -1L) THEN ymax_ecmwf = i
     ENDIF ELSE BEGIN
        FOR i=nlon_ecmwf-1L,0L,-1L DO IF (round(lons_ecmwf[i]/pr,/L64) LE round((lonmin+0.5d0*dlon)/pr,/L64)) $
           AND (xmin_ecmwf EQ -1L) THEN xmin_ecmwf = i
        FOR i=nlat_ecmwf-1L,0L,-1L DO IF (round(lats_ecmwf[i]/pr,/L64) LE round((latmin+0.5d0*dlat)/pr,/L64)) $
           AND (ymin_ecmwf EQ -1L) THEN ymin_ecmwf = i
        FOR i=0L,nlon_ecmwf-1L,1L DO IF (round(lons_ecmwf[i]/pr,/L64) GE round((lonmax-0.5d0*dlon)/pr,/L64)) $
           AND (xmax_ecmwf EQ -1L) THEN xmax_ecmwf = i
        FOR i=0L,nlat_ecmwf-1L,1L DO IF (round(lats_ecmwf[i]/pr,/L64) GE round((latmax-0.5d0*dlat)/pr,/L64)) $
           AND (ymax_ecmwf EQ -1L) THEN ymax_ecmwf = i
     ENDELSE

     ;; try to recover from out of bounds in polar regions and dateline
     IF (xmax_ecmwf EQ (-1L)) AND (max(lon) GT max(lons_ecmwf)) AND (lonmax_ecmwf EQ 180.d0) THEN xmax_ecmwf = nlon_ecmwf-1L
     IF (xmin_ecmwf EQ (-1L)) AND (min(lon) LT min(lons_ecmwf)) AND (lonmin_ecmwf EQ -180.d0) THEN xmin_ecmwf = 0L
     IF (ymax_ecmwf EQ (-1L)) AND (max(lat) GT max(lats_ecmwf)) AND (latmax_ecmwf EQ 90.d0) THEN ymax_ecmwf = nlat_ecmwf-1L
     IF (ymin_ecmwf EQ (-1L)) AND (min(lat) LT min(lats_ecmwf)) AND (latmin_ecmwf EQ -90.d0) THEN ymin_ecmwf = 0L

     IF ((xmin_ecmwf EQ (-1L)) OR (ymin_ecmwf EQ (-1L)) OR (xmax_ecmwf EQ (-1L)) OR (ymax_ecmwf EQ (-1L))) THEN BEGIN
        print,'ecmwf dataset does not cover requested area. stopping.'
        print,'We need: ',lonmin,lonmax,latmin,latmax
        print,'We get:  ',lonmin_ecmwf,lonmax_ecmwf,latmin_ecmwf,latmax_ecmwf
        print,xmin_ecmwf,xmax_ecmwf,ymin_ecmwf,ymax_ecmwf
        print,lons_ecmwf[xmin_ecmwf],lons_ecmwf[xmax_ecmwf],lats_ecmwf[ymin_ecmwf],lats_ecmwf[ymax_ecmwf]
        stop
     ENDIF

     nx_ecmwf = xmax_ecmwf - xmin_ecmwf + 1L
     ny_ecmwf = ymax_ecmwf - ymin_ecmwf + 1L

     ;; adjust for boundary conditions in bilinear interpolation
     lontmp = (lon>(lonmin_ecmwf+0.5d0*dlon_ecmwf))<(lonmax_ecmwf-0.5d0*dlon_ecmwf)
     lattmp = (lat>(latmin_ecmwf+0.5d0*dlat_ecmwf))<(latmax_ecmwf-0.5d0*dlat_ecmwf) 

     ;; derive bilinear interpolation indices
     x_ecmwf = (lontmp - lons_ecmwf[xmin_ecmwf])/dlon_ecmwf
     y_ecmwf = (lattmp - lats_ecmwf[ymin_ecmwf])/dlat_ecmwf

;     print,x_ecmwf
;     print,y_ecmwf

     ;; read topography
     ncdf_varget,ncid, 'Z', data, offset = [xmin_ecmwf,ymin_ecmwf,0], count = [nx_ecmwf,ny_ecmwf,1]

     ;; close NetCDF file
     ncdf_close,ncid

     ;; convert topography from m2/s2 to m
     data /= 9.81

     ;; maintain original grid, even if 1x1
     data = reform(data,nx_ecmwf,ny_ecmwf)

     ;; store ECMWF topography on original grid
     z_ecmwf_grid = reform(data,nx_ecmwf,ny_ecmwf)

     ;; regrid topography to output grid
     IF regular THEN BEGIN

        IF (dlon GE dlon_ecmwf) AND (dlat GE dlat_ecmwf) THEN BEGIN
           ;; rebinning for upscaling regular grid
           z_ecmwf = rebin_new(data,nlon,nlat)
        ENDIF ELSE BEGIN
           ;; bilinear interpolation for downscaling regular grid
           z_ecmwf = interpolate(data,x_ecmwf,y_ecmwf)
        ENDELSE
        
     ENDIF ELSE BEGIN
        ;; use bilinear interpolation for irregular grid
        z_ecmwf = interpolate(data,x_ecmwf,y_ecmwf)
     ENDELSE

  ENDIF

  prefix = varname

  ;; check for daily file
  t_ecmwf = 0L
  nt_ecmwf = nt
  infile = file_search(ecmwfpath+prefix+'.'+syear+smonth+sday+'.nc')
  IF infile[0] EQ '' THEN BEGIN
     ;; check for monthly file
     t_ecmwf = (long(day) - 1L) * ntday
     nt_ecmwf = nt
     infile = file_search(ecmwfpath+prefix+'.'+syear+smonth+'.nc')
  ENDIF

  IF infile[0] NE '' THEN BEGIN

     ncid = ncdf_open(infile[0])
     ncdf_varget,ncid, varname, data, offset = [xmin_ecmwf,ymin_ecmwf,t_ecmwf], $
                 count = [nx_ecmwf,ny_ecmwf,nt_ecmwf]
     ncdf_close,ncid

     ;; try fitting new lapse rate or scale height
     IF n_elements(elevation) EQ 0L AND (stddev(z_ecmwf_grid) GT 200.) THEN BEGIN

        ;; fit lapse rate [m]
        IF keyword_set(lapse_rate) THEN BEGIN
           lapse_fit = make_array(nt,/float,value=missing)
           FOR t = 0,nt-1 DO BEGIN
              fit = regress(z_ecmwf_grid[*],(data[*,*,t])[*],correlation=correlation,ftest=ftest)
;;              print,correlation[0], ftest[0], fit[0]
              IF (abs(correlation) GT 0.5) AND (ftest GT 50.) THEN lapse_fit[t] = fit[0]
           ENDFOR
        ENDIF

        ;; fit scale height [m]
        IF keyword_set(scale_height) THEN BEGIN
           scale_fit = make_array(nt,/float,value=missing)
           FOR t = 0,nt-1 DO BEGIN
              fit = regress(z_ecmwf_grid[*],(alog(data[*,*,t]))[*],const=const,correlation=correlation,ftest=ftest)
;;              print,correlation[0], ftest[0], -1./fit[0],const
              IF (correlation LT -0.5) AND (ftest GT 50.) THEN scale_fit[t] = -1./fit[0]
           ENDFOR
        ENDIF

     ENDIF

     ;; regrid topography to output grid
     FOR t=0,nt-1 DO BEGIN

        IF regular THEN BEGIN

           IF (dlon GE dlon_ecmwf) AND (dlat GE dlat_ecmwf) THEN BEGIN
              ;; rebinning for upscaling regular grid
              IF regular OR ll2d THEN BEGIN
                 outdata[*,*,t] = rebin_new(data[*,*,t],nlon,nlat)
              ENDIF ELSE BEGIN
                 outdata[*,t] = rebin_new(data[*,*,t],nlon,nlat)
              ENDELSE
           ENDIF ELSE BEGIN
              ;; bilinear interpolation for downscaling regular grid
              IF regular OR ll2d THEN BEGIN
                 outdata[*,*,t] = interpolate(data[*,*,t],x_ecmwf,y_ecmwf)
              ENDIF ELSE BEGIN
                 outdata[*,t] = interpolate(data[*,*,t],x_ecmwf,y_ecmwf)
              ENDELSE
           ENDELSE
           
        ENDIF ELSE BEGIN
           ;; use bilinear interpolation for irregular grid
           IF regular OR ll2d THEN BEGIN
              outdata[*,*,t] = interpolate(data[*,*,t],x_ecmwf,y_ecmwf)
           ENDIF ELSE BEGIN
              outdata[*,t] = interpolate(data[*,*,t],x_ecmwf,y_ecmwf)
           ENDELSE
        ENDELSE

     ENDFOR

     ;; downscaling
     IF n_elements(elevation) NE 0L THEN BEGIN

        ;; linearly downscale output variable
        IF keyword_set(lapse_rate) THEN BEGIN

           FOR t=0,nt-1 DO BEGIN
;              IF keyword_set(lapse_fit) THEN BEGIN
;                 IF lapse_fit[t] NE missing THEN lapse = lapse_fit[t] ELSE lapse = lapse_rate
;              ENDIF ELSE BEGIN
              lapse = lapse_rate
;              ENDELSE

              IF regular OR ll2d THEN BEGIN
                 outdata[*,*,t] = outdata[*,*,t] + (float(elevation) - z_ecmwf) * lapse
              ENDIF ELSE BEGIN
                 outdata[*,t] = outdata[*,t] + (float(elevation) - z_ecmwf) * lapse
              ENDELSE
           ENDFOR

        ENDIF

        ;; exponentially downscale output variable
        IF keyword_set(scale_height) THEN BEGIN

           FOR t=0,nt-1 DO BEGIN
;              IF keyword_set(scale_fit) THEN BEGIN
;                 IF scale_fit[t] NE missing THEN scale = scale_fit[t] ELSE scale = scale_height
;              ENDIF ELSE BEGIN
              scale = scale_height
;              ENDELSE

              IF regular OR ll2d THEN BEGIN
                 outdata[*,*,t] *= exp((z_ecmwf - float(elevation))/scale)
              ENDIF ELSE BEGIN
                 outdata[*,t] *= exp((z_ecmwf - float(elevation))/scale)
              ENDELSE
           ENDFOR

        ENDIF

     ENDIF

  ENDIF ELSE BEGIN
     print,'ecmwf file containing ',varname,' for YYYYMMDD ',syear,smonth,sday,' not available..'
     stop
  ENDELSE

  RETURN, outdata

END
