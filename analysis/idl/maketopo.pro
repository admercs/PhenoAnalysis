;+
; :Description:
;   This script creates a subgrid - topography distribution 
;   ((flag: histo=1) or a mean composite topography map (histo=0) from the 90 m
;   CGIAR SRTM v4.1 and 1 km GTOPO30 topography at a chosen output resolution.
;
;   Note: The RADARSAT-II Antarctica support was removed. We'll
;   wait for either the IceSat GLAS or the ASTER DEM for polar areas
;   and use GTOPO30 in the meantime.
;
;   Dependencies: 
;   - CGIAR SRTM Dataset in GeoTIFF format. Download from: ftp://srtm.csi.cgiar.org/SRTM_v41/SRTM_Data_GeoTIFF/
;   - GTOPO30 Dataset in Binary DEM format. Download from: ftp://edcftp.cr.usgs.gov/pub/data/gtopo30/global/
; 
; :Categories:
;   Geographical Data Processing
;
; :Params:
;
; :Keywords:
;   basedir: in, required, type="string"
;     directory where all the input and output data resides
;   topocase: in, required, type="string"
;     pft case (experiment) name, required for TOPO directory within basedir
;   lonmin0 : in, required, type="double"
;     western (edge, not center) boundary of output dataset
;   lonmax0 : in, required, type="double"
;     eastern (edge, not center) boundary of output dataset
;   latmin0 : in, required, type="double"
;     southern (edge, not center) boundary of output dataset
;   latmax0 : in, required, type="double"
;     northern (edge, not center) boundary of output dataset
;   dlon0 : in, required, type="double"
;     longitude spacing of output grid
;   lat0 :  in, required, type="double"
;     latitude spacing of output grid
;   region : in, required, type="string"
;     name of the area of interest (becomes part of output file name)
;   xproc : in, optional, type="integer"
;     longitudinal processing subtile number
;   yproc : in, optional, type="integer"
;     latitudinal processing subtile number
;   nxprocs : in, optional, type="integer" 
;     number of longitudinal subtiles in total
;   nyprocs : in, optional, type="integer"
;     number of latitudinal subtiles in total
;   histo: in, optional, type="boolean"
;      produce elevation histogram (1) or mean elevation map (0)
;
; :Returns: 
; 
; :Uses:
;   read_cgiarsrtm
;   read_gtopo30
;   rebin_new
;
; :Bugs:
;
; :Todo:
;
; :Requires:
;   IDL 7.1
;
; :Examples:
;
; :History:
;   2007-2008: Initial development at Colorado State University
;
;   2009: Continued development at Blue Marble Research
;
;   2012/09/28 : revised code, included IDL/RST documentation
;
;   2012/10/01 : removed lonlen, latlen, dlon0, dlat0 keywords and
;   replaced with nx0, ny0, dlon, dlat keywords
;
;   2013/12/30 : rewritten interface
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
PRO maketopo,basedir=basedir, topocase=topocase, lonmin0=lonmin0,lonmax0=lonmax0,$
             latmin0=latmin0,latmax0=latmax0,dlon0=dlon0,dlat0=dlat0,region=region, $
             xproc=xproc,yproc=yproc,nxprocs=nxprocs,nyprocs=nyprocs,histo=histo


  ;; read command line arguments (useful for IDL VM or RT mode)
  ;; they are the same as the program arguments above
  args = command_line_args(count=nargs)
  IF nargs GT 0 THEN BEGIN
     FOR a=0,nargs-1 DO BEGIN
        temp = strtrim(strsplit(args[a],'=',/extract),2)
        IF n_elements(temp) EQ 2 THEN BEGIN
           CASE temp[0] OF
              'basedir' : basedir = temp[1]
              'topocase' : topocase = temp[1]
              'lonmin0' :  lonmin0 = double(temp[1])
              'lonmax0' :  lonmax0 = double(temp[1])
              'latmin0' :  latmin0 = double(temp[1])
              'latmax0' :  latmax0 = double(temp[1])
              'dlon0' :  dlon0 = double(temp[1])
              'dlat0' :  dlat0 = double(temp[1])
              'region' : region = temp[1]
              'xproc' : xproc = fix(temp[1])
              'yproc' : yproc = fix(temp[1])
              'nxprocs' : nxprocs = fix(temp[1])
              'nyprocs' : nyprocs = fix(temp[1])
              'histo' : histo = fix(temp[1])
              ELSE : BEGIN
                 print,'Command Line Argument ',args[a],' not implemented.'
                 stop
              END
           ENDCASE
        ENDIF ELSE BEGIN
           print,'Command Line Argument ',args[a],' not understood.'
           print,'Please specify it like: argument=value .'
           stop
        ENDELSE
     ENDFOR
  ENDIF

;; check for cluster processing information
  IF not keyword_set(nxprocs) THEN nxprocs = 1
  IF not keyword_set(nyprocs) THEN nyprocs = 1
  IF not keyword_set(xproc) THEN xproc = 1
  IF not keyword_set(yproc) THEN yproc = 1

;; version
  version = "v1.90"
  
  ;; set directories
  
  IF topocase EQ "" THEN topodir = basedir + '/topo/' ELSE topodir = basedir + '/topo-' + topocase + '/'
  gtopodir = basedir + '/topo/GTOPO30/'
  srtmdir = basedir + '/topo/SRTM_CGIAR/GeoTiff/'

;; report fpe's
  !except=2
  
;; flags
  verbose = 1 ;; print stuff or not
  graphics = 0 

;; elevation ranges (for distributed - histogram elevation mode)
  minhgt = 0.
  maxhgt = 4000.
  nhgt = 40

;; how do we fill the missing areas (ocean?)
  topo_fillvalue = 0

;; read topography at its native resolution (approx 100 m at the
;; equator)
  dlon_raw = 1.0d/1200.0d
  dlat_raw = 1.0d/1200.0d

;; we'll work through the chosen output domain with tiles defined by nx0 and ny0
;; in order to save memory.

;; code profiling information
  time_total = 0.
  time_read = 0.
  time_calc = 0.
  time_write = 0.

  time_total -= systime(1)

;; initialize random numbers
  seed = long(double(xproc+100)*double(yproc+1)*systime(1,/seconds)/1.d7)

;; -------- START MAIN CODE ---------

  ;; prepare geographic information of full grid
  region=strtrim(region,2)

  nlon0 = round((lonmax0 - lonmin0)/dlon0)
  nlat0 = round((latmax0 - latmin0)/dlat0)
  dlon0 = (lonmax0 - lonmin0)/double(nlon0)
  dlat0 = (latmax0 - latmin0)/double(nlat0)

;; calculate grid spacing of output sub-grid
  tx0 = long(nlon0)*long(xproc-1L)/long(nxprocs)
  tx1 = long(nlon0)*long(xproc)/long(nxprocs)-1L
  ty0 = long(nlat0)*long(yproc-1L)/long(nyprocs)
  ty1 = long(nlat0)*long(yproc)/long(nyprocs)-1L
  lonmin = double(tx0) * (lonmax0 - lonmin0) / double(nlon0) + lonmin0
  lonmax = double(tx1+1L) * (lonmax0 - lonmin0) / double(nlon0) + lonmin0
  latmin = double(ty0) * (latmax0 - latmin0) / double(nlat0) + latmin0
  latmax = double(ty1+1L) * (latmax0 - latmin0) / double(nlat0) + latmin0

  dlon = dlon0
  dlat = dlat0
  nlon = round((lonmax - lonmin)/dlon)
  nlat = round((latmax - latmin)/dlat)

;; create lon/lat arrays (centered on pixel)
  lon = ((dindgen(nlon) + 0.5d0)*(lonmax-lonmin)/double(nlon) + lonmin) # replicate(1.d0,nlat) ; W->E
  lat = replicate(1.d0,nlon) # ((dindgen(nlat) + 0.5d0)*(latmax-latmin)/double(nlat) + latmin) ; S->N

  ;; North Pole (do not change)
  pollon = 180.0
  pollat = 90.0

  print,'Longitude Bounds: ',lonmin,'-',lonmax
  print,'Latitude Bounds:  ',latmin,'-',latmax
  print,'We expect : ',nlon0,'x',nlat0,' grid points in full dataset'
  print,'We expect : ',nlon,'x',nlat,' grid points in this dataset'

;; define graphics output if needed
  IF graphics THEN BEGIN
     device,decomposed=0,retain=2,set_font='HELVETICA',/tt_font
     loadct,0
  ENDIF

;; define elevation levels (histogram mode)
  hgtlevels = (findgen(nhgt)+0.5)/float(nhgt)*(maxhgt-minhgt) + minhgt
  
;; define netcdf names
  IF keyword_set(histo) THEN BEGIN
     nvar = 1
     ncvarid=lonarr(nvar)
     ncvarnames = ['hgt']
     ncvarlongnames = ['Subgrid-Scale elevation distribution']
     ncvarunits = '%'
     ncrange=[0B,100B]
     ncmissing = 129B
     topofile='hgt.'+region+'.nc'
  ENDIF ELSE BEGIN
     nvar = 1
     ncvarid=lonarr(nvar)
     ncvarnames = ['Z']
     ncvarlongnames = ['Elevation above sea level']
     ncvarunits = 'm'
     ncrange=[-10000s,10000s]
     ncmissing = -32767s
     topofile='topo.'+region+'.nc'
  ENDELSE

  ;; initialize NetCDF file
  
  ;; create lock
  file_lock,topodir,topofile,/lock,id=100L*xproc+yproc

  time_write -= systime(1)

  ;; ---------- BEGIN NETCDF DEFINE ----------
  exist = 1B
  temp = file_search(topodir+topofile)
  
  IF strtrim(temp) EQ '' THEN BEGIN
     print,topodir+topofile+' does not exist; creating new one.'
     exist = 0B
  ENDIF
  
  IF exist EQ 0B THEN BEGIN

     ;; calculate current time, format for NetCDF output
     caldat,systime(/julian,/utc),prod_month,prod_day,prod_year,prod_hour,prod_minute,prod_second
     
     prod_date = string(prod_year,format='(I4)')+'-'+string(prod_month,format='(I2.2)')+'-'+ $
                 string(prod_day,format='(I2.2)')+' '+string(prod_hour,format='(I2.2)')+':'+ $
                 string(prod_minute,format='(I2.2)')+':'+string(round(prod_second),format='(I2.2)')
     
     ;; create NetCDF file (header and dimensions etc.)
     ncid = NCDF_CREATE(topodir+topofile,/CLOBBER)
     NCDF_CONTROL,ncid,/NOFILL
     londim = NCDF_DIMDEF(ncid, 'lon',nlon0)
     latdim = NCDF_DIMDEF(ncid, 'lat',nlat0)
     IF keyword_set(histo) THEN hgtdim = NCDF_DIMDEF(ncid, 'Z',nhgt)
     
     ;; create lon/lat variables for a regular grid
     
     lonid = NCDF_VARDEF(ncid, 'lon', [londim], /DOUBLE)
     NCDF_ATTPUT,ncid,lonid,'long_name','Longitude'
     NCDF_ATTPUT,ncid,lonid,'units','degrees_east'
     
     latid = NCDF_VARDEF(ncid, 'lat', [latdim], /DOUBLE)
     NCDF_ATTPUT,ncid,latid,'long_name','Latitude'
     NCDF_ATTPUT,ncid,latid,'units','degrees_north'
     
     IF keyword_set(histo) THEN BEGIN
        hgtid = NCDF_VARDEF(ncid, 'Z', [hgtdim], /FLOAT)
        NCDF_ATTPUT,ncid,hgtid,'long_name','Elevation above sea level'
        NCDF_ATTPUT,ncid,hgtid,'positive','up'
        NCDF_ATTPUT,ncid,hgtid,'units','m'
     ENDIF

     ;; create data variables
     FOR v = 0,nvar-1 DO BEGIN
        
        ;; 2D data variable
        IF keyword_set(histo) THEN BEGIN
           ncvarid[v] = NCDF_VARDEF(ncid,ncvarnames[v], [londim,latdim,hgtdim], /BYTE)
        ENDIF ELSE BEGIN
           ncvarid[v] = NCDF_VARDEF(ncid,ncvarnames[v], [londim,latdim], /SHORT)
        ENDELSE
        ;; Attributes for data variables
        NCDF_ATTPUT, ncid, ncvarid[v], '_FillValue',ncmissing
        NCDF_ATTPUT, ncid, ncvarid[v], 'valid_range',ncrange           
        NCDF_ATTPUT, ncid, ncvarid[v], 'long_name',ncvarlongnames[v]
        NCDF_ATTPUT, ncid, ncvarid[v], 'units', ncvarunits
        NCDF_ATTPUT, ncid, ncvarid[v], 'version', version
        NCDF_ATTPUT, ncid, ncvarid[v], 'prod_date', prod_date
        
     ENDFOR

     ;; create global attributes
     IF keyword_set(histo) THEN BEGIN
        NCDF_ATTPUT, ncid, /GLOBAL, 'title','Fractional Elevation Distribution Map'
     ENDIF ELSE BEGIN
        NCDF_ATTPUT, ncid, /GLOBAL, 'title','Elevation Map'
     ENDELSE
     NCDF_ATTPUT, ncid, /GLOBAL, 'institution','Blue Marble Research'
     NCDF_ATTPUT, ncid, /GLOBAL, 'source','SRTM (CGIAR) + GTOPO30 (USGS)'
     NCDF_ATTPUT, ncid, /GLOBAL, 'Author','Reto Stockli'
     NCDF_ATTPUT, ncid, /GLOBAL, 'references','Bliss, N.B., and Olsen, L.M., 1996. Development of a 30-arc-second digital elevation model of South America. In: Pecora Thirteen, Human Interactions with the Environment - Perspectives from Space, Sioux Falls, South Dakota, August 20-22, 1996. AND Jarvis A., H.I. Reuter, A.  Nelson, E. Guevara, 2008, Hole-filled  seamless SRTM data V4, International  Centre for Tropical  Agriculture (CIAT), available  from http://srtm.csi.cgiar.org.'
     NCDF_ATTPUT, ncid, /GLOBAL, 'Conventions','CF-1.4'
     
     NCDF_CONTROL, ncid, /ENDEF ; Put file in data mode. 

     ;; ----------- END NETCDF DEFINE ---------------

     ;; write elevation levels
     IF keyword_set(histo) THEN NCDF_VARPUT, ncid, hgtid, hgtlevels
     
     NCDF_CLOSE, ncid

  ENDIF

  ;; remove lock
  file_lock,topodir,topofile,/unlock,id=100L*xproc+yproc

  time_write += systime(1)

  time_read -= systime(1)

  ;; read extended SRTM topography, re-adjust geographical boundaries
  ;; with actual SRTM grid spacing and grid bounds
  ;; GTOPO30 processing <60S or >60N
  IF (latmin LT (-60.d0)) AND (latmax GT 60.d0) THEN BEGIN
     ;; most difficult case: fill N and S boundary with GTOPO30
     read_gtopo30,gtopodir,lonmin,latmin,lonmax,-60.d0,dlon_raw,dlat_raw,gtopo_south
     read_gtopo30,gtopodir,lonmin,60.d0,lonmax,latmax,dlon_raw,dlat_raw,gtopo_north
     read_cgiarsrtm,srtmdir,lonmin,-60.d0,lonmax,60.d0,dlon_raw,dlat_raw,stopo_middle
     topo = [[gtopo_south],[stopo_middle],[gtopo_north]]
     gtopo_south = 0
     stopo_middle = 0
     gtopo_north = 0
  ENDIF
  IF (latmin LT (-60.d0)) AND (latmax LE 60.d0) THEN BEGIN
     ;; only fill S boundary with GTOPO30
     IF latmax LE (-60.d0) THEN BEGIN
        read_gtopo30,gtopodir,lonmin,latmin,lonmax,latmax,dlon_raw,dlat_raw,topo
     ENDIF ELSE BEGIN
        read_gtopo30,gtopodir,lonmin,latmin,lonmax,-60.0d,dlon_raw,dlat_raw,gtopo_south
        read_cgiarsrtm,srtmdir,lonmin,-60.d0,lonmax,latmax,dlon_raw,dlat_raw,stopo_middle
        topo = [[gtopo_south],[stopo_middle]]
        gtopo_south = 0
        stopo_middle = 0
     ENDELSE
  ENDIF
  IF (latmin GE (-60.d0)) AND (latmax GT 60.d0) THEN BEGIN
     ;; only fill N boundary with GTOPO30
     IF latmin GE 60.d0 THEN BEGIN
        read_gtopo30,gtopodir,lonmin,latmin,lonmax,latmax,dlon_raw,dlat_raw,topo
     ENDIF ELSE BEGIN
        read_gtopo30,gtopodir,lonmin,60.d0,lonmax,latmax,dlon_raw,dlat_raw,gtopo_north
        read_cgiarsrtm,srtmdir,lonmin,latmin,lonmax,60.d0,dlon_raw,dlat_raw,stopo_middle
        topo = [[stopo_middle],[gtopo_north]]
        stopo_middle = 0
        gtopo_north = 0
     ENDELSE
  ENDIF
  IF (latmin GE (-60.d0)) AND (latmax LE 60.d0) THEN BEGIN
     ;; simple case: SRTM only
     read_cgiarsrtm,srtmdir,lonmin,latmin,lonmax,latmax,dlon_raw,dlat_raw,topo
  ENDIF

  time_read += systime(1)
  
  time_calc -= systime(1)

  ;; build output arrays
  IF keyword_set(histo) THEN BEGIN
     ncvar = bytarr(nlon,nlat,nhgt)
  ENDIF ELSE BEGIN
     ncvar = intarr(nlon,nlat)
  ENDELSE
  
  IF keyword_set(histo) THEN BEGIN
     ;; make a histogram of elevation bins
     scx = round(dlon/dlon_raw)
     scy = round(dlat/dlat_raw)
     FOR y=0,nlat-1 DO BEGIN
        FOR x=0,nlon-1 DO BEGIN
           a = histogram((topo[x*scx:(x+1)*scx-1,y*scy:(y+1)*scy-1]>minhgt)<(maxhgt-1.0), $
                         min=minhgt,max=maxhgt-(hgtlevels[1]-hgtlevels[0]),nbins=nhgt)

           sum = 0
           pct = round(a/total(a)*100.)
           FOR h=0,nhgt-1 DO BEGIN
              sum = sum + pct[h]
              ncvar[x,y,h] = pct[h] - ((sum - 100)>0)
              sum = sum < 100
           ENDFOR

        ENDFOR
     ENDFOR
     
     print,'check for <100% total everywhere'
     less = where(total(ncvar,3) LT 100,lesscount)
     WHILE (lesscount GT 0) DO BEGIN
        v = fix(randomu(seed)*nhgt)
        temp = fix(ncvar[*,*,v])
        temp[less] = (temp[less] + 1*(temp[less] NE 0.))<100
        ncvar[*,*,v] = byte(temp)
        less = where(total(ncvar,3) LT 100,lesscount)
     ENDWHILE

     IF verbose THEN BEGIN
        FOR h = 0,nhgt-1 DO print,'level: ',hgtlevels[h],' %: ',mean(ncvar[*,*,h])
        print,'total: ',mean(total(ncvar,3))
        print,'min: ',min(ncvar), ' max: ',max(ncvar)
     ENDIF
  ENDIF ELSE BEGIN

     dims = size(topo,/dimensions)
     nlon_raw = dims[0]
     nlat_raw = dims[1]

     IF (nlon NE nlon_raw) OR (nlat NE nlat_raw) THEN BEGIN

        IF ((nlon_raw EQ (nlon_raw/nlon)*nlon) AND (nlat_raw EQ (nlat_raw/nlat)*nlat)) OR $
           ((float(nlon)/float(nlon_raw) EQ float(nlon/nlon_raw)) AND $
            (float(nlat)/float(nlat_raw) EQ float(nlat/nlat_raw))) THEN BEGIN
           ;; simple rebinning -- fast
           ncvar[*] = rebin(topo,nlon,nlat) 
        ENDIF ELSE BEGIN
           ;; use new rebinning for non-integer dimensional scaling
           ncvar[*] = fix(rebin_new(topo,nlon,nlat))
        ENDELSE

     ENDIF ELSE BEGIN
        ;; no regridding
        ncvar[*] = topo
     ENDELSE           

     IF verbose THEN BEGIN
        print,'Topography: min: ',min(ncvar), ' max: ',max(ncvar)
     ENDIF
  ENDELSE
  
  time_calc += systime(1)
  
  ;; --------- NETCDF WRITE DATA SECTION ---------------

  time_write -= systime(1)
  
  ;; Now write the data to the NetCDF file(s)

  ;; only open NetCDF file[s] if they are not in use by another process
  file_lock,topodir,topofile,/lock,id=100L*xproc+yproc

  ;; open file
  ncid = NCDF_OPEN(topodir+topofile,/write)

  latid = ncdf_varid(ncid,'lat')
  lonid = ncdf_varid(ncid,'lon')
  ncvarid = ncdf_varid(ncid,ncvarnames[0])

  ;; write longitude, latitude
  NCDF_VARPUT, ncid, lonid, reform(lon[*,0]), offset=[tx0], count=[nlon]
  NCDF_VARPUT, ncid, latid, reform(lat[0,*]), offset=[ty0], count=[nlat]
  
  ;; write data
  IF keyword_set(histo) THEN BEGIN
     NCDF_VARPUT, ncid, ncvarid,ncvar, offset=[tx0,ty0,0],count=[nlon,nlat,nhgt]
  ENDIF ELSE BEGIN
     NCDF_VARPUT, ncid, ncvarid,ncvar, offset=[tx0,ty0],count=[nlon,nlat]
  ENDELSE
  
  ;; close NetCDF file
  NCDF_CLOSE, ncid 

  ;; remove lock
  file_lock,topodir,topofile,/unlock,id=100L*xproc+yproc

  time_write += systime(1)

  time_total += systime(1)

  print
  print,"Total wall clock by process: "
  print,'Read         : ',time_read
  print,'Calc         : ',time_calc
  print,'Write        : ',time_write
  print,'TOTAL        : ',time_total

END
