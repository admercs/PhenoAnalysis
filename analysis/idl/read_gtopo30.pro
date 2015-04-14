;+
; :Description:
;   This code read the GTOPO30 Tiles and
;   composites them for a given output region and resolution.
;
; :Categories:
;   Geographic Processing
;
; :Params:
;   gtopodir : in, required, type=string 
;     absolute path (directory) to the GTOPO DEM binary files
;   lonmin: in, required, type="float or double"
;     western edge of requested domain (degrees E)
;   latmin: in, required, type="float or double"
;      southern edge of requested domain (degrees N)
;   lonmax: in, required, type="float or double"
;      eastern edge of requested domain (degrees E)
;   latmax: in, required, type="float or double"
;      northern edge of requested domain (degrees N)
;   dlon: in, required, type="float or double"
;      output resolution west-east (degrees)
;   dlat: in, required, type="float or double"
;      output resolution south-north (degrees)
;   topo: out, required, type="intarr"
;      topography on chosen grid (2D grid)
;;
; :Requires:
;   IDL 7.1
;
; :Examples:
;
; :Uses:
;   rebin_new
;
; :History:
;   2010/11/22 Reto Stockli (Blue Marble Research)
;
;   2011/03/18 : added standard idldoc-compatible rst-style header
;
;   2014/01/06 : fixed handling of Antarctic tiles and missing values
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
PRO read_gtopo30, gtopodir, lonmin, latmin, lonmax, latmax, dlon, dlat, topo

  ;; how do we fill the missing areas (ocean?)
  topo_fillvalue = 0

  ;; output dimensions
  nlon = round((lonmax - lonmin)/dlon)
  nlat = round((latmax - latmin)/dlat)
  
  ;; GTOPO30 files: geographical file structure
  gtopo_nx = 4800L
  gtopo_ny = 6000L
  gtopo_nx2 = 7200L
  gtopo_ny2 = 3600L

  gtopo_nlon = 40L
  gtopo_nlat = 50L
  gtopo_nlon2 = 60L
  gtopo_nlat2 = 30L

  gtopo_lon0 = -180.0
  gtopo_lat0 = 90.0
  gtopo_lon1 = 180.0
  gtopo_lat1 = -90.0

  gtopo_lonlen = round(gtopo_nx/gtopo_nlon)
  gtopo_latlen = round(gtopo_ny/gtopo_nlat)
  gtopo_dlon = double(gtopo_nlon)/double(gtopo_nx) 
  gtopo_dlat = double(gtopo_nlat)/double(gtopo_ny) 

  ;; north of -60S

  ;; western longitudes for each input file x index
  gtopo_lonind0 = lindgen(round((gtopo_lon1-gtopo_lon0)/gtopo_nlon))*gtopo_nlon*gtopo_lonlen + $
                  round(gtopo_lon0*double(gtopo_lonlen)) 
  ;; eastern longitudes for each input file x index
  gtopo_lonind1 = lindgen(round((gtopo_lon1-gtopo_lon0)/gtopo_nlon))*gtopo_nlon*gtopo_lonlen + $
                  round(gtopo_lon0*double(gtopo_lonlen)) + gtopo_nlon*gtopo_lonlen 
  ;; southern latitude for each input file y index  
  gtopo_latind0 = - lindgen(round((gtopo_lat0-gtopo_lat1)/gtopo_nlat))*gtopo_nlat*gtopo_latlen + $
                  round(gtopo_lat0*double(gtopo_latlen)) - gtopo_nlat*gtopo_latlen 
  ;; northern latitude for each input file y index
  gtopo_latind1 = - lindgen(round((gtopo_lat0-gtopo_lat1)/gtopo_nlat))*gtopo_nlat*gtopo_latlen + $
                  round(gtopo_lat0*double(gtopo_latlen)) 

  ;; south of -60S

  ;; western longitudes for each input file x index
  gtopo_lonind2 = lindgen(round((gtopo_lon1-gtopo_lon0)/gtopo_nlon2))*gtopo_nlon2*gtopo_lonlen +$
                  round(gtopo_lon0*double(gtopo_lonlen)) 
  ;; eastern longitudes for each input file x index
  gtopo_lonind3 = lindgen(round((gtopo_lon1-gtopo_lon0)/gtopo_nlon2))*gtopo_nlon2*gtopo_lonlen + $
                  round(gtopo_lon0*double(gtopo_lonlen)) + gtopo_nlon2*gtopo_lonlen 

  gtopo_latind0[3] = -90L*gtopo_latlen
  gtopo_latind1[3] = -60L*gtopo_latlen

;; GTOPO30 missing data value 
  gtopo_nodata = -9999

  lonind0 = round(lonmin*double(gtopo_lonlen))
  lonind1 = round(lonmax*double(gtopo_lonlen))
  latind0 = round(latmin*double(gtopo_latlen))
  latind1 = round(latmax*double(gtopo_latlen))

  nlon_gtopo = lonind1 - lonind0
  nlat_gtopo = latind1 - latind0
  
  print,'We expect : ',nlon_gtopo,'x',nlat_gtopo,' GTOPO30 grid points on input'
  
  ;; GTOPO30 input data
  gtopo30 = intarr(nlon_gtopo,nlat_gtopo)

  print,'searching for: ',float(lonind0)*gtopo_dlon,float(lonind1)*gtopo_dlon,$
        float(latind0)*gtopo_dlat,float(latind1)*gtopo_dlat

  IF latind1 LE -60 THEN BEGIN
     cx0 = where((lonind0 GE gtopo_lonind2) AND (lonind0 LT gtopo_lonind3),x0count)
     cx1 = where((lonind1 GT gtopo_lonind2) AND (lonind1 LE gtopo_lonind3),x1count)
  ENDIF ELSE BEGIN
     cx0 = where((lonind0 GE gtopo_lonind0) AND (lonind0 LT gtopo_lonind1),x0count)
     cx1 = where((lonind1 GT gtopo_lonind0) AND (lonind1 LE gtopo_lonind1),x1count)
  ENDELSE
  cy1 = where((latind0 GE gtopo_latind0) AND (latind0 LT gtopo_latind1),y0count)
  cy0 = where((latind1 GT gtopo_latind0) AND (latind1 LE gtopo_latind1),y1count)

  IF (x0count EQ 0L) OR (x1count EQ 0L) OR (y0count EQ 0L) OR (y1count EQ 0L) THEN BEGIN
     print,"GTOPO30 Bounds reached. Try SRTM."
     stop
  ENDIF 

  FOR cy=cy0[0],cy1[0] DO BEGIN
     FOR cx=cx0[0],cx1[0] DO BEGIN

        
        IF latind1 LE -60L*gtopo_latlen THEN BEGIN
           scx = string(abs(round(double(gtopo_lonind2[cx])*gtopo_dlon)),format='(I3.3)')
           IF gtopo_lonind2[cx] LE 0 THEN scx='W'+scx ELSE scx='E'+scx
        ENDIF ELSE BEGIN
           scx = string(abs(round(double(gtopo_lonind0[cx])*gtopo_dlon)),format='(I3.3)')
           IF gtopo_lonind0[cx] LE 0 THEN scx='W'+scx ELSE scx='E'+scx
        ENDELSE
        scy = string(abs(round(double(gtopo_latind1[cy])*gtopo_dlat)),format='(I2.2)')
        IF gtopo_latind1[cy] LE 0 THEN scy='S'+scy ELSE scy='N'+scy

        gtopofile = scx+scy+'.DEM'

        print,'Searching for: ',gtopofile
        infile = file_search(gtopodir+gtopofile,count=count)
        IF count GT 0L THEN BEGIN
           
           print,'Reading: ',gtopodir+gtopofile

           ;; cut desired sub-area
           IF latind1 LE -60L*gtopo_latlen THEN BEGIN
              x0 = lonind0 - gtopo_lonind2[cx]
              x1 = lonind1 - gtopo_lonind2[cx]-1
              y0 = latind0 - gtopo_latind0[cy]
              y1 = latind1 - gtopo_latind0[cy]-1

              dx0 = (x0>0) - x0
              dx1 = x1 - (x1<(gtopo_nx2-1))
              dy0 = (y0>0) - y0
              dy1 = y1 - (y1<(gtopo_ny2-1))
           ENDIF ELSE BEGIN
              x0 = lonind0 - gtopo_lonind0[cx]
              x1 = lonind1 - gtopo_lonind0[cx]-1
              y0 = latind0 - gtopo_latind0[cy]
              y1 = latind1 - gtopo_latind0[cy]-1
              
              dx0 = (x0>0) - x0
              dx1 = x1 - (x1<(gtopo_nx-1))
              dy0 = (y0>0) - y0
              dy1 = y1 - (y1<(gtopo_ny-1))
           ENDELSE
           
           x0 = x0 + dx0
           x1 = x1 - dx1
           y0 = y0 + dy0
           y1 = y1 - dy1

           nx = x1-x0+1
           ny = y1-y0+1

           ;; geotiff is saved with LL corner = 1/1 
           IF latind1 LE -60L*gtopo_latlen THEN BEGIN
              nxtot = ulong(gtopo_nx2)
              y0 = gtopo_ny2 - 1 - y0
              y1 = gtopo_ny2 - 1 - y1
           ENDIF ELSE BEGIN
              nxtot = ulong(gtopo_nx)
              y0 = gtopo_ny - 1 - y0
              y1 = gtopo_ny - 1 - y1
           ENDELSE

;           IF nx GT 0 THEN BEGIN

              gtopoline = intarr(nx)
              
              openr,lun,gtopodir+gtopofile,/get_lun,/swap_if_little_endian
              FOR y=0,ny-1 DO BEGIN

                 yr = y

                 IF latind1 LE -60L*gtopo_latlen THEN BEGIN
                    ;; southernmost line in Antarctica is missing. Why?
                    IF (y1+y) EQ 3599 THEN yr -= 1
                 ENDIF

                 point_lun,lun,2UL*(ulong(x0) + ulong(y1+yr)*nxtot)
                 readu,lun,gtopoline

                 IF (latind1 LE -60L*gtopo_latlen) AND (lonind1 GE 21599L) THEN BEGIN
                    ;; easternmost 2 columns are missing in
                    ;; Antarctica. Why?
                    gtopoline[(nx-2)>0:nx-1] = gtopoline[(nx-3)>0]
                 ENDIF

                 IF (latind1 LE -60L*gtopo_latlen) AND (lonind0 LE -21600L) THEN BEGIN
                    ;; westernmost 1 column is missing in
                    ;; Antarctica. Why?
                    gtopoline[0] = gtopoline[1<(nx-1)]
                 ENDIF
                 
                 gtopo30[dx0:dx0+nx-1,dy0+ny-y-1] = gtopoline

              ENDFOR
              close,lun
              free_lun,lun

;           ENDIF
                        
        ENDIF

     ENDFOR
  ENDFOR

  dlon_raw = gtopo_dlon
  dlat_raw = gtopo_dlat

  idx = where(gtopo30 EQ gtopo_nodata,count)
  IF count GT 0L THEN gtopo30[idx] = topo_fillvalue

  ;; rescale GTOPO30 if needed
  IF (nlon NE nlon_gtopo) OR (nlat NE nlat_gtopo) THEN BEGIN

     IF (nlon GT nlon_gtopo) OR (nlat GT nlat_gtopo) THEN BEGIN
        ;; magnifying: bilinear interpolation
        topo = congrid(gtopo30,nlon,nlat,/interp,/center)

     ENDIF ELSE BEGIN
        ;; minifying: two methods
        IF ((nlon_gtopo EQ (nlon_gtopo/nlon)*nlon) AND (nlat_gtopo EQ (nlat_gtopo/nlat)*nlat)) OR $
           ((float(nlon)/float(nlon_gtopo) EQ float(nlon/nlon_gtopo)) AND $
            (float(nlat)/float(nlat_gtopo) EQ float(nlat/nlat_gtopo))) THEN BEGIN
           ;; simple rebinning -- fast
           topo = rebin(gtopo30,nlon,nlat) 
        ENDIF ELSE BEGIN
           ;; new dimensions are not integer factor of old dimensions
           topo = fix(rebin_new(gtopo30,nlon,nlat))
        ENDELSE

     ENDELSE

  ENDIF ELSE BEGIN
     ;; no regridding
     topo = gtopo30
  ENDELSE

END
