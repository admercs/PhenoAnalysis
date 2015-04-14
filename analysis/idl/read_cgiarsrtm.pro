;+
; :Description:
;   This code read the CGIAR SRTM Version 4 GeoTiff Tiles and
;   composites them for a given output region and resolution.
;
;   For more information on the CGIAR srtm please see:
;
;   http://srtm.csi.cgiar.org/
;
;   For downloading the CGIAR srtm data and bypassing their web
;   interface please go to:
;
;   ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/
;
;   ftp://xftp.jrc.it/pub/srtmV4/tiff/
;
; :Categories:
;   Geographic Processing
;
; :Params:
;   srtmdir : in, required, type=string 
;     absolute path (directory) to the CGIAR SRTM GeoTiff files
;   lonmin: in, required, type="float or double"
;     western edge of requested domain (degrees E)
;   latmin: in, required, type="float or double"
;      southern edge of requested domain (degrees N)
;   lonmax: in, required, type="float or double"
;      eastern edge of requested domain (degrees E)
;   latmax: in, required, type="float or double"
;      northern edge of requested domain (degrees N)
;   dlon: in, required, type="float or double"
;      requested resolution west-east (degrees)
;   dlat: in, required, type="float or double"
;      reqested resolution south-north (degrees)
;   topo: out, required, type="intarr"
;      topography on chosen grid (2D grid)
;
; :Keywords:
;   out_lonmin: out, optional, type="float or double"
;     western edge of extracted domain (degrees E)
;   out_latmin: out, optional, type="float or double"
;      southern edge of extracted domain (degrees N)
;   out_lonmax: out, optional, type="float or double"
;      eastern edge of extracted domain (degrees E)
;   out_latmax: out, optional, type="float or double"
;      northern edge of extracted domain (degrees N)
;   out_dlon: out, optional, type="float or double"
;      output resolution west-east (degrees)
;   out_dlat: out, optional, type="float or double"
;      output resolution south-north (degrees)
;
; :Bugs:
;  The code likely fails when reading SRTM data from 180.0E since the
;  1/2 pixel shift would require to remap the easternmost pixel around
;  the dateline    
;
; :Uses:
;   rebin_new
;
; :Requires:
;   IDL 7.1
;
; :Examples:
;
; :History:
;   2010/11/22 Reto Stockli (Blue Marble Research)
;
;   2011/03/18 : added standard idldoc-compatible rst-style header
;
;   2012/07/09 : added cgiar_lon/lat/min/max keywords which give the
;   actual extracted domain (that may differ from the requrested
;   domain by +/- grid size of the CGIAR SRTM
;
;   2013/07/16 : corrected CGIAR-SRTM by 1/2 pixel shift in N and E
;   direction
; 
;   2014/01/08 : fixed integer bug and implemented new regridding
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
PRO read_cgiarsrtm, srtmdir, lonmin, latmin, lonmax, latmax, dlon, dlat, topo, $
                    out_lonmin=out_lonmin,out_lonmax=out_lonmax, $
                    out_latmin=out_latmin,out_latmax=out_latmax, $
                    out_dlon=out_dlon, out_dlat=out_dlat

  ;; how do we fill the missing areas (ocean?)
  topo_fillvalue = 0

  ;; SRTM CGIAR files: geographical file structure
  cgiar_nx = 6000L             ; number of columns in GeoTIFF
  cgiar_ny = 6000L             ; number of lines in GeoTIFF

  cgiar_nlon = 5L              ; longitude spread for GeoTIFF
  cgiar_nlat = 5L              ; latitude spread for GeoTIFF

  ;; derived CGIAR geographical parameters
  cgiar_lonlen = round(cgiar_nx/cgiar_nlon) ; number of SRTM columns per degree lon
  cgiar_latlen = round(cgiar_ny/cgiar_nlat) ; number of SRTM lines per degree lat
  cgiar_dlon = double(cgiar_nlon)/double(cgiar_nx) ; longitude spacing of SRTM dataset pixels
  cgiar_dlat = double(cgiar_nlat)/double(cgiar_ny) ; latitude spacing of SRTM dataset pixels

  ;; origin of CGIAR tiles
  cgiar_lon0 = -180.d0           ; starting longitude of input file x index
  cgiar_lon1 =  180.d0           ; ending longitude of input file x index
  cgiar_lat0 = -60.d0            ; starting longitude of input file x index
  cgiar_lat1 =  60.d0            ; ending latitude of input file y index

  ;; western longitudes for each input file x index
  cgiar_lonind0 = lindgen(round((cgiar_lon1-cgiar_lon0)/cgiar_nlon))*cgiar_nlon*cgiar_lonlen + $
                  round(cgiar_lon0*double(cgiar_lonlen))
  ;; eastern longitudes for each input file x index
  cgiar_lonind1 = lindgen(round((cgiar_lon1-cgiar_lon0)/cgiar_nlon))*cgiar_nlon*cgiar_lonlen + $
                  round(cgiar_lon0*double(cgiar_lonlen)) + cgiar_nlon*cgiar_lonlen
  ;; southern latitude for each input file y index  
  cgiar_latind0 = - lindgen(round((cgiar_lat1-cgiar_lat0)/cgiar_nlat))*cgiar_nlat*cgiar_latlen + $
                  round(cgiar_lat1*double(cgiar_latlen)) - cgiar_nlat*cgiar_latlen 
  ;; northern latitude for each input file y index
  cgiar_latind1 = - lindgen(round((cgiar_lat1-cgiar_lat0)/cgiar_nlat))*cgiar_nlat*cgiar_latlen + $
                  round(cgiar_lat1*double(cgiar_latlen)) 

  ;; SRTM missing data value (Version 4 of CGIAR SRTM)
  cgiar_nodata = -32768

  ;; introduce 1/2 pixel shift due to changed CGIAR SRTM pixel
  ;; origin. Please see: http://srtm.csi.cgiar.org/SRTM_FAQ.asp
  lonind0 = round(((lonmin + 0.5d0*cgiar_dlon)>cgiar_lon0)*double(cgiar_lonlen))
  lonind1 = round(((lonmax + 0.5d0*cgiar_dlon)<cgiar_lon1)*double(cgiar_lonlen))
  latind0 = round(((latmin + 0.5d0*cgiar_dlat)>cgiar_lat0)*double(cgiar_latlen))
  latind1 = round(((latmax + 0.5d0*cgiar_dlat)<cgiar_lat1)*double(cgiar_latlen))

  nlon_cgiar = lonind1 - lonind0
  nlat_cgiar = latind1 - latind0
  
  ;; generate output array
  srtm = make_array(nlon_cgiar,nlat_cgiar,/integer,value=topo_fillvalue)

  ;; adjust output bounds to integer SRTM pixel bounds and grid distance by keeping the
  ;; requested amount of output pixels
  nlon = round((lonmax - lonmin)/dlon)
  nlat = round((latmax - latmin)/dlat)
  out_lonmin=float(double(lonind0)/double(cgiar_lonlen)-0.5d0*cgiar_dlon)
  out_lonmax=float(double(lonind1)/double(cgiar_lonlen)-0.5d0*cgiar_dlon)
  out_latmin=float(double(latind0)/double(cgiar_latlen)-0.5d0*cgiar_dlat)
  out_latmax=float(double(latind1)/double(cgiar_latlen)-0.5d0*cgiar_dlat)
  out_dlon = (out_lonmax-out_lonmin)/float(nlon)
  out_dlat = (out_latmax-out_latmin)/float(nlat)

  print,'Reading CGIAR SRTM with bounds: ',out_lonmin,'-',out_lonmax,'/', $
        out_latmin,'-',out_latmax,format="(A,F12.6,A,F12.6,A,F12.6,A,F12.6)"
  print,'Number of SRTM grid points : ',nlon_cgiar,'x',nlat_cgiar,' on input',format="(A,I6,A,I6,A)"
  print,'Number of SRTM grid points : ',nlon,'x',nlat,' on output',format="(A,I6,A,I6,A)"
  
  print,'Searching for global indices: ',lonind0,'-',lonind1,' / ',latind0,'-',latind1,format="(A,I8,A,I8,A,I8,A,I8)"

  cx0 = where((lonind0 GE cgiar_lonind0) AND (lonind0 LT cgiar_lonind1),x0count) + 1L
  cx1 = where((lonind1 GT cgiar_lonind0) AND (lonind1 LE cgiar_lonind1),x1count) + 1L
  cy0 = where((latind1 GT cgiar_latind0) AND (latind1 LE cgiar_latind1),y0count) + 1L
  cy1 = where((latind0 GE cgiar_latind0) AND (latind0 LT cgiar_latind1),y1count) + 1L

  print,'Searching for global tiles:   ',cx0[0],'-',cx1[0],' / ',cy0[0],'-',cy1[0],format="(A,I8,A,I8,A,I8,A,I8)"

  IF (x0count EQ 0L) OR (x1count EQ 0L) OR (y0count EQ 0L) OR (y1count EQ 0L) THEN BEGIN
     print,"CGIAR SRTM Bounds reached. Stopping."
     stop
  ENDIF

  FOR cy=cy0[0],cy1[0] DO BEGIN
     FOR cx=cx0[0],cx1[0] DO BEGIN

        scx = string(cx,format='(I2.2)')
        scy = string(cy,format='(I2.2)')

        srtmfile = 'srtm_'+scx+'_'+scy        

        infile = file_search(srtmdir+srtmfile+'.tif',count=count)
        IF count EQ 0L THEN BEGIN
           infile = file_search(srtmdir+srtmfile+'.tif.gz',count=count)
           
           IF count GT 0L THEN BEGIN
              spawn,'gunzip -f '+srtmdir+srtmfile+'.tif.gz'
           ENDIF
        ENDIF
        
        IF count EQ 0L THEN BEGIN
           infile = file_search(srtmdir+srtmfile+'.zip',count=count)
           
           IF count GT 0L THEN BEGIN
              spawn,'cd '+srtmdir+'; unzip -o '+srtmfile+'.zip'
              spawn,'\rm -f '+srtmdir+srtmfile+'.zip'
           ENDIF
        ENDIF
        
        infile = file_search(srtmdir+srtmfile+'.tif',count=count)
        IF count GT 0L THEN BEGIN
           
           print,'Reading: ',srtmdir+srtmfile+'.tif'
           
           ;; cut desired sub-area
           x0 = lonind0 - cgiar_lonind0[cx-1]
           x1 = lonind1 - cgiar_lonind0[cx-1] - 1L
           y0 = latind0 - cgiar_latind0[cy-1]
           y1 = latind1 - cgiar_latind0[cy-1] - 1L
           
           dx0 = (x0>0L) - x0
           dx1 = x1 - (x1<(cgiar_nx-1L))
           dy0 = (y0>0L) - y0
           dy1 = y1 - (y1<(cgiar_ny-1L))
           
           x0 = x0 + dx0
           x1 = x1 - dx1
           y0 = y0 + dy0
           y1 = y1 - dy1
           
           nx = x1-x0+1L
           ny = y1-y0+1L
           
           ;; read full dataset instead of sub-rectangle (why: see just below)
           ;; reverse: geotiff is saved with LL corner = 1/1 
           data = reverse(read_tiff(srtmdir+srtmfile+'.tif'),2)

           ;; gaps in CGIAR SRTM are extensive water bodies, fill with standard fill value
           ;; sometimes we have a 255 fill values in deep ocean
           ;; areas. Andy tells me that it happens with GeoTiff if no
           ;; values within the file are above the byte range.
           cgiar_nodata2 = 255          
           
           IF max(data) EQ cgiar_nodata2 THEN BEGIN
              data[where(data EQ cgiar_nodata2)] = cgiar_nodata
           ENDIF
 
           srtm[dx0:dx0+nx-1L,dy0:dy0+ny-1L] = data[x0:x1,y0:y1]

       ENDIF ELSE BEGIN
           print,'Not Found: ',srtmdir+srtmfile+'.tif'
        ENDELSE

     ENDFOR
  ENDFOR

  idx = where(srtm EQ cgiar_nodata,count)
  IF count GT 0L THEN srtm[idx] = topo_fillvalue

  ;; fix for some integer issues due to grid-point edge definition of
  ;; CGIAR SRTM data
  IF (nlon GT 1000L) AND (nlon_cgiar EQ (nlon-1L)) THEN BEGIN
     srtm = [srtm[0L,*],srtm]
     nlon_cgiar = nlon
  ENDIF

  IF (nlon GT 1000L) AND (nlon_cgiar EQ (nlon+1L)) THEN BEGIN
     srtm = srtm[1L:nlon_cgiar-1L,*]
     nlon_cgiar = nlon
  ENDIF

  IF (nlat GT 1000L) AND (nlat_cgiar EQ (nlat-1L)) THEN BEGIN
     srtm = [[srtm[*,0L]],[srtm]]
     nlat_cgiar = nlat
  ENDIF

  IF (nlat GT 1000L) AND (nlat_cgiar EQ (nlat+1L)) THEN BEGIN
     srtm = srtm[*,1L:nlat_cgiar-1L]
     nlat_cgiar = nlat
  ENDIF

  ;; rescale SRTM if needed
  IF (nlon NE nlon_cgiar) OR (nlat NE nlat_cgiar) THEN BEGIN

     IF ((nlon_cgiar EQ (nlon_cgiar/nlon)*nlon) AND (nlat_cgiar EQ (nlat_cgiar/nlat)*nlat)) OR $
        ((float(nlon)/float(nlon_cgiar) EQ float(nlon/nlon_cgiar)) AND $
         (float(nlat)/float(nlat_cgiar) EQ float(nlat/nlat_cgiar))) THEN BEGIN
        ;; simple rebinning -- fast
        topo = rebin(srtm,nlon,nlat) 
     ENDIF ELSE BEGIN
        ;; use new rebinning for non-integer dimensional scaling
        topo = fix(rebin_new(srtm,nlon,nlat))
     ENDELSE

  ENDIF ELSE BEGIN
     ;; no regridding
     topo = srtm
  ENDELSE

END
