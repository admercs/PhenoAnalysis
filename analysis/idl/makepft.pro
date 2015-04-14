;+
; :Description:
;   The code creates a regional subgridscale PFT distribution maps according to the
;   algorithms presented in P.J. Lawrence et al. 2007, extended here
;   for a high resolution 1 km subgrid-scale PFT map
;
;   Dependencies: 
;   In the directory basedir/ the following sub-directories with datasets need to exist:
;   1. modis/MOD15A2.005/       --> global MODIS FPAR/LAI dataset
;   2. modis/MCD12Q1.005/       --> global MODIS land cover dataset
;   3. modis/MOD44B.004/        --> Hansen et al. global MODIS-based vcf dataset
;   4. avhrr/vcf/               --> DeFries et al. global tree cover / leaf type dataset
;   5. crops/                   --> Ramankutty & Foley global crop dataset
;   6. ecmwf/era-interim/       --> ERA INTERIM temperature and precipitation climatology
;   7. topo-pftcase/            --> topography dataset, has to be created by IDL code maketopo.pro
;   8. pft-pftcase/             --> output directory
;
; :Categories:
;   Geographic / Biophysical Processing
;
; :Params:
;
; :Keywords:
;   basedir: in, required, type="string"
;     directory where all the input and output data resides
;   pftcase: in, required, type="string"
;     pft case (experiment) name, required for PFT directory within basedir
;   manyfiles: in, required, type="boolean"
;     create a single big NetCDF file or one file per pft (when single
;     NetCDF 3 File would exceed 2GB)
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
;
; :Returns: 
; 
; :Uses:
;
; :Bugs:
;
; :Todo:
;   Transfer code to Fortran 95
;   Re-visit climatological splitting of pft's by use of a
;   gradient of temperature / precip requirements
;
; :Requires:
;   IDL 7.1
;
; :Examples:
;
; :History:
;   2007 : first version as part of my Colorado State University
;   employment
;   2014 : rewrite for Collection 5 MODIS datasets and ECMWF
;   climatology
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
PRO makepft,basedir=basedir, pftcase=pftcase, manyfiles=manyfiles, lonmin0=lonmin0,lonmax0=lonmax0,$
            latmin0=latmin0,latmax0=latmax0,dlon0=dlon0,dlat0=dlat0,region=region, $
            xproc=xproc,yproc=yproc,nxprocs=nxprocs,nyprocs=nyprocs

  COMMON MODISTILE,vmodis,hmodis,xmodis,ymodis,xidmodis,yidmodis,modis_lon,modis_lat, $
     npmodis,hmin,hmax,vmin,vmax

  ;; read command line arguments (useful for IDL VM or RT mode)
  ;; they are the same as the program arguments above
  args = command_line_args(count=nargs)
  IF nargs GT 0 THEN BEGIN
     FOR a=0,nargs-1 DO BEGIN
        temp = strtrim(strsplit(args[a],'=',/extract),2)
        IF n_elements(temp) EQ 2 THEN BEGIN
           CASE temp[0] OF
              'basedir' : basedir = temp[1]
              'pftcase' : pftcase = temp[1]
              'manyfiles' : manyfiles = fix(temp[1])
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
  IF ~keyword_set(nxprocs) THEN nxprocs = 1
  IF ~keyword_set(nyprocs) THEN nyprocs = 1
  IF ~keyword_set(xproc) THEN xproc = 1
  IF ~keyword_set(yproc) THEN yproc = 1

;; version
  version = "v1.90"

;; graphics 
  graphics = 0

;; report fpe's
  !except=2


  IF graphics THEN BEGIN
     device,decomposed=0,retain=2,set_font='HELVETICA',/tt_font
     loadct,0
  ENDIF
  
;; flags
  wsc = 5                       ; scaling factor for display
  verbose = 1                   ; print stuff or not

;; we'll work through the whole globe in tiles here!
;; don't pick too high dlat0 and dlon0 tiles else we run out of memory
 
;; start and end year for surface climatology and LAI
  year_start = 2001
  year_end = 2010

;; number of months in a year
  nmonth = 12 

;; default nodata value
  nodata = -9999.

;; initialize random numbers
  seed = long(double(xproc+100)*double(yproc+1)*systime(1,/seconds)/1.d7)

;; code profiling information
  time_total = 0.
  time_netcdf = 0.
  time_climread = 0.
  time_climcalc = 0.
  time_lairead = 0.
  time_laicalc = 0.
  time_cropread = 0.
  time_avhrrread  = 0.
  time_lcread = 0.
  time_vcfread = 0.
  time_calc = 0.

  time_total -= systime(1)

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

;; netcdf names

  nvar = 35
  ncvarnames = ['bar_all',$
                'enf_tem',$
                'enf_bor',$
                'dnf_bor',$
                'ebf_tro',$
                'ebf_tem',$
                'dbf_tro',$
                'dbf_tem',$
                'dbf_bor',$
                'ebs_all',$
                'dbs_tem',$
                'dbs_bor',$
                'c3g_arc',$
                'c3g_nar',$
                'c4g_all',$
                'cro_brl',$
                'cro_cas',$
                'cro_cot',$
                'cro_grn',$
                'cro_mze',$
                'cro_mil',$
                'cro_oil',$
                'cro_oth',$
                'cro_pot',$
                'cro_pul',$
                'cro_rap',$
                'cro_ric',$
                'cro_rye',$
                'cro_sor',$
                'cro_soy',$
                'cro_sgb',$
                'cro_sgc',$
                'cro_sun',$
                'cro_wht',$
                'wat_all']

  ncvarlongnames = ['Bare Soil (-)',$
                    'Evergreen Needleleaf Trees (Temperate)', $
                    'Evergreen Needleleaf Trees (Boreal)', $
                    'Deciduous Needleleaf Trees (Boreal)', $
                    'Evergreen Broadleaf Trees (Tropical)', $
                    'Evergreen Broadleaf Trees (Temperate)', $
                    'Deciduous Broadleaf Trees (Tropical)', $
                    'Deciduous Broadleaf Trees (Temperate)', $
                    'Deciduous Broadleaf Trees (Boreal)',$
                    'Evergreen Broadleaf Shrubs (All)',$
                    'Deciduous Broadleaf Shrubs (Temperate)',$
                    'Deciduous Broadleaf Shrubs (Boreal)',$
                    'C3 Grass (Arctic)',$
                    'C3 Grass (Non-Arctic)',$
                    'C4 Grass (All)',$
                    'Crops (Barley)','Crops (Cassava)','Crops (Cotton)',$
                    'Crops (Groundnuts)','Crops (Maize)','Crops (Millet)', $
                    'Crops (Oilpalm)','Crops (Other)','Crops (Potato)','Crops (Pulses)', $
                    'Crops (Rape)','Crops (Rice)','Crops (Rye)','Crops (Sorghum)',$
                    'Crops (Soy)','Crops (Sugarbeets)','Crops (Sugarcane)',$
                    'Crops (Sunflower)','Crops (Wheat)',$
                    'Water (All)']
  
  ncvarunits = replicate('%',nvar)

;; ---------- BEGIN NETCDF DEFINE ----------

  time_netcdf -= systime(1)

  IF pftcase EQ "" THEN pftdir = basedir + '/pft/' ELSE pftdir = basedir + '/pft-' + pftcase + '/'

  IF manyfiles THEN BEGIN
     pftfile=ncvarnames+'.'+region+'.nc' 
  ENDIF ELSE BEGIN
     pftfile='pft.'+region+'.nc'
  ENDELSE

  IF manyfiles THEN BEGIN
     nfiles = nvar
  ENDIF ELSE BEGIN
     nfiles = 1
  ENDELSE

  ncid = lonarr(nfiles)
  lonid = lonarr(nfiles)
  latid = lonarr(nfiles)
  ncvarid = lonarr(nvar)

  FOR file = 0, nfiles-1 DO BEGIN
     ;; create NetCDF file (header and dimensions etc.)
     ;; if xproc or yproc are 0 independently if it exists.

     ;; only open NetCDF file[s] if they are not in use by another process
     file_lock,pftdir,pftfile[file],/lock,id=100L*xproc+yproc

     exist = 1B
     temp = file_search(pftdir+pftfile[file])
     
     IF strtrim(temp) EQ '' THEN BEGIN
        print,pftdir+pftfile[file]+' does not exist; creating new one.'
        exist = 0B
     ENDIF
     
     IF exist EQ 0B THEN BEGIN
        ;; we need to create a new file with new id's

        ;; calculate current time, format for NetCDF output
        caldat,systime(/julian,/utc),prod_month,prod_day,prod_year,prod_hour,prod_minute,prod_second
        
        prod_date = string(prod_year,format='(I4)')+'-'+string(prod_month,format='(I2.2)')+'-'+ $
                    string(prod_day,format='(I2.2)')+' '+string(prod_hour,format='(I2.2)')+':'+ $
                    string(prod_minute,format='(I2.2)')+':'+string(round(prod_second),format='(I2.2)')

        ncid[file] = NCDF_CREATE(pftdir+pftfile[file],/CLOBBER)
        NCDF_CONTROL,ncid[file],/NOFILL
        londim = NCDF_DIMDEF(ncid[file], 'lon',nlon0)
        latdim = NCDF_DIMDEF(ncid[file], 'lat',nlat0)

        ;; create lon/lat variables

        ;; regular grid        
        lonid[file] = NCDF_VARDEF(ncid[file], 'lon', [londim], /DOUBLE)
        NCDF_ATTPUT,ncid[file],lonid[file],'long_name','Longitude'
        NCDF_ATTPUT,ncid[file],lonid[file],'units','degrees_east'
        
        latid[file] = NCDF_VARDEF(ncid[file], 'lat', [latdim], /DOUBLE)
        NCDF_ATTPUT,ncid[file],latid[file],'long_name','Latitude'
        NCDF_ATTPUT,ncid[file],latid[file],'units','degrees_north'

        missing = 129B ;; signed byte value corresponds to 255B for unsigned byte in NetCDF
        scale = 1.0
        offset = 0.0

        ;; create data variables
        IF manyfiles THEN BEGIN
           v=file

           ;; define PFT variables
           
           ;; 2D data variable
           ncvarid[v] = NCDF_VARDEF(ncid[file],ncvarnames[v], [londim,latdim], /BYTE)
           
           ;; Attributes for data variables
           NCDF_ATTPUT, ncid[file], ncvarid[v], '_FillValue',missing
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'valid_range',[0B,100B]           
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'long_name',ncvarlongnames[v]
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'units', ncvarunits[v]
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'version', version
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'prod_date', prod_date

        ENDIF ELSE BEGIN
           FOR v = 0,nvar-1 DO BEGIN
              ;; define PFT variables
              
              ;; 2D data variable
              ncvarid[v] = NCDF_VARDEF(ncid[file],ncvarnames[v], [londim,latdim], /BYTE)
              
              ;; Attributes for data variables
              NCDF_ATTPUT, ncid[file], ncvarid[v], '_FillValue',missing
              NCDF_ATTPUT, ncid[file], ncvarid[v], 'valid_range',[0B,100B]
              NCDF_ATTPUT, ncid[file], ncvarid[v], 'long_name',ncvarlongnames[v]
              NCDF_ATTPUT, ncid[file], ncvarid[v], 'units', ncvarunits[v]
              NCDF_ATTPUT, ncid[file], ncvarid[v], 'version', version
              NCDF_ATTPUT, ncid[file], ncvarid[v], 'prod_date', prod_date
              
           ENDFOR

        ENDELSE

        ;; create global attributes
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Title','Fractional plant functional type distribution'
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Institution','Blue Marble Research'
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Source','MOD44B, MOD12Q, AVHRR VCF, Crop database'
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Author','Reto Stockli'
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'References','Reference Publication: P.J. Lawrence et al. (2007). Representing a new MODIS consistent land surface in the Community Land Model (CLM 3.0). JGR, Vol. 112, G01023'
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Conventions','CF-1.4'
        
        NCDF_CONTROL, ncid[file], /ENDEF ; Put file in data mode. 

        ;; close NetCDF file
        NCDF_CLOSE, ncid[file] 

        print,'Finished initializing NetCDF output file: ',pftdir+pftfile[file]

     ENDIF

     ;; remove lock
     file_lock,pftdir,pftfile[file],/unlock,id=100L*xproc+yproc

  ENDFOR

  time_netcdf += systime(1)

;; ----------- END NETCDF DEFINE ---------------


  ;; create high resolution climate fields
  time_climread -= systime(1)

  ;; read topography           
  IF pftcase EQ "" THEN topodir = 'topo' ELSE topodir = 'topo-' + pftcase
  read_topodata, basedir, topodir, region, lon, lat, topo, dlon=dlon, dlat=dlat

  ecmwfdir = basedir + '/ecmwf/era-interim/'

  ;; create high resolution climate fields
  thighres = make_array(nlon,nlat,nmonth,/float,value=0.)
  phighres = make_array(nlon,nlat,nmonth,/float,value=0.)
  nhighres = make_array(nlon,nlat,nmonth,/long,value=0L)

;;  goto,bla

  FOR year = year_start,year_end DO BEGIN

     ;; testing for leap year
     monthlen = [31,28,31,30,31,30,31,31,30,31,30,31]
     IF ((((year MOD 4) EQ 0) AND ((year MOD 100) NE 0)) OR $
         (((year MOD 100) EQ 0) AND ((year MOD 400) EQ 0))) THEN BEGIN
        leapyear = 1
        monthlen[1] = 29
        ndayofyear = total(monthlen)
     ENDIF ELSE BEGIN
        leapyear = 0
        ndayofyear = total(monthlen)
     ENDELSE
     
     FOR month = 1,12 DO BEGIN

        ;; number of days in current month
        ndayofmonth = monthlen[month-1]
        
        FOR day = 1,ndayofmonth DO BEGIN
           tmp = read_ecmwf(ecmwfdir, '2T', lon, lat, year, month, day, $
                            dlon=dlon, dlat=dlat, elevation=topo, lapse_rate=-0.0065, $
                            missing=!values.f_nan)
           
           thighres[*,*,month-1] += total(tmp,3)

           badidx = where(finite(tmp,/nan),badcount)
           IF badcount GT 0L THEN BEGIN
              print, "Missing ECMWF data. Stopping."
              stop
           ENDIF

           tmp = read_ecmwf(ecmwfdir, 'LSP', lon, lat, year, month, day, $
                            dlon=dlon, dlat=dlat, elevation=topo, $
                            missing=!values.f_nan)
           
           phighres[*,*,month-1] += total(tmp,3)*1000.

           badidx = where(finite(tmp,/nan),badcount)
           IF badcount GT 0L THEN BEGIN
              print, "Missing ECMWF data. Stopping."
              stop
           ENDIF

           tmp = read_ecmwf(ecmwfdir, 'CP', lon, lat, year, month, day, $
                            dlon=dlon, dlat=dlat, elevation=topo, $
                            missing=!values.f_nan)
           
           phighres[*,*,month-1] += total(tmp,3)*1000.

           badidx = where(finite(tmp,/nan),badcount)
           IF badcount GT 0L THEN BEGIN
              print, "Missing ECMWF data. Stopping."
              stop
           ENDIF

           nhighres[*,*,month-1] += 4L

        ENDFOR

     ENDFOR

  ENDFOR
  
  ;; create monthly mean temperature fields (K)
  thighres /= float(nhighres)

  ;; convert from K to degree Celsius
  thighres -= 273.16

  ;; create monthly total precipitation fields (mm per month)
  phighres /= float(year_end - year_start + 1)

  bla:

  time_climread += systime(1)
  
  time_climcalc -= systime(1)

  ;; minimum monthly temperature
  tmin = min(thighres,dimension=3)

  ;; maximum monthly temperature
  tmax = max(thighres,dimension=3)

  ;; GDD (daily sums of temperature above threshold!)
  gddmin = 5.0 ;; degC
  gdd = total((congrid(thighres,nlon,nlat,365,/interp,/center) - gddmin)>0.,3)

  ;; annual total precipitation
  pann = total(phighres,3)

  ;; winter total precipitation
  ;; NH:  NOV-APR, SH: MAY-OCT
  nhidx = [0,1,2,3,10,11]
  shidx = [4,5,6,7,8,9]
  pwin = total(phighres[*,*,nhidx],3)*(lat GT 0.) + total(phighres[*,*,shidx],3)*(lat LE 0.)

  ;; minimum monthly precipitation
  pmin = min(phighres,dimension=3)

  ;; maximum precipitation for month where T > 22degC
  pmaxtg22 = max(phighres * (thighres GT 22.0),dimension=3)

  ;; test if we can expect any pixels with C4 climate for
  ;; any month, if yes, read the MODIS FPAR distribution:

  c4potidx = where((thighres GT 22.0) AND (phighres GT 25.0),c4potcount)

  time_climcalc += systime(1)

  IF c4potcount GT 0 THEN BEGIN

     ;; read MODIS LAI for entire year to distribute C3/C4 grasses

     time_lairead -= systime(1)

     lai = make_array(nlon,nlat,nmonth,/float,value=!values.f_nan)
     laimin = make_array(nlon,nlat,/float,value=!values.f_nan)
     laimax = make_array(nlon,nlat,/float,value=!values.f_nan)

     lai_version = 5

     ;; initalize MODIS grid
     read_modistile, basedir, 0, 0, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'LAI', obsvarunits, obsvarlongnames, lai_version, temp_data, temp_error, temp_mask, $
                     pollon=pollon, pollat=pollat,/initialize


     lai_year_start = year_start
     lai_year_end   = year_end
     nplai = npmodis
     ntlai = 46
     dtlai = 8
     nylai = lai_year_end - lai_year_start + 1
     lai_data = make_array(nplai,ntlai,nylai,/float,value=!values.f_nan)
     

     FOR lai_year = lai_year_start,lai_year_end DO BEGIN
        FOR t=0,ntlai-1 DO BEGIN
           
           lai_day = 1+t*dtlai

           read_modistile, basedir, lai_year, lai_day, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                           'LAI', obsvarunits, obsvarlongnames, lai_version, temp_data, temp_error, temp_mask, $
                           pollon=pollon, pollat=pollat
           
           lai_data[*,t,lai_year - lai_year_start] = temp_data
           
        ENDFOR
     ENDFOR

     time_lairead += systime(1)

     time_laicalc -= systime(1)

     idx1 = round(findgen(nmonth)*ntlai/float(nmonth))
     idx2 = round((findgen(nmonth)+1.0)*ntlai/float(nmonth)-1.0)

     ;; scale LAI with min/max values
     scale_minmax = 0B


     IF total(finite(lai_data)) LT 1.0 THEN BEGIN
        print,"no valid LAI data found"
;        stop
     ENDIF

     ;; regridding
     xyindex = yidmodis*nlon + xidmodis
     hist = histogram(xyindex,omin=omin,reverse_indices = R)
     FOR i=0L,n_elements(hist)-1L DO BEGIN
        IF R[i] NE R[i+1] THEN BEGIN
           FOR m = 0,nmonth-1 DO BEGIN
              IF total(finite(lai_data[R[R[i]:R[i+1]-1],idx1[m]:idx2[m],*])) GT 2.0 THEN $
                 lai[(i+omin)+nlon*nlat*m] = $
                 mean(lai_data[R[R[i]:R[i+1]-1],idx1[m]:idx2[m],*],/nan)
           ENDFOR

           IF scale_minmax THEN BEGIN

              temparr = lai_data[R[R[i]:R[i+1]-1],*,*]
              
              goodlai = where(finite(temparr),goodlaicount)
              IF goodlaicount GT 5 THEN BEGIN
                 arr = temparr[goodlai]
                 ;; LAImin 5% Quantile, LAImax 98% quantile
                 laimin[i+omin]  = (arr) [ (sort(arr))[5  * (goodlaicount-1) / 100] ]
                 laimax[i+omin]  = (arr) [ (sort(arr))[98 * (goodlaicount-1) / 100] ]
              ENDIF

           ENDIF

        ENDIF
     ENDFOR

     temp_data = 0
     temp_error = 0
     temp_mask = 0
     xmodis = 0
     ymodis = 0
     hmodis = 0
     vmodis = 0
     xidmodis = 0
     yidmodis = 0
     xyindex = 0
     hist = 0
     R = 0

     ;; close holes in LAI
     idx=where(finite(lai,/nan),count,ncomplement=count2)
     IF (count2 EQ 0) THEN BEGIN
        ;; no MODIS data available
        lai[*] = 0.
     ENDIF ELSE BEGIN
        ;; MODIS data available, check for gaps
        IF (count GT 0) THEN BEGIN
           ;; gaps available, fill gaps
           lai[idx]=nodata

           FOR m = 0,nmonth-1 DO BEGIN
              lai_temp = reform(lai[*,*,m])

              gaps = where(lai_temp EQ nodata, gapcount)
              dia = 4
              WHILE (gapcount GT 0L) AND (dia LT 32) DO BEGIN
                 
                 print,'Month: ',m+1,'. Trying to fill ',gapcount, $
                       ' holes in monthly LAI with radius: ',dia
                 
                 lai_temp[gaps] = (fill_dist(lai_temp,dia,nodata))[gaps]
                 
                 dia=dia*2
                 gaps = where(lai_temp EQ nodata, gapcount)
                 
              ENDWHILE
              
              ;; remaining gaps: probably water, fill with 0.
              IF gapcount GT 0L THEN lai_temp[gaps] = 0.
              
              lai[*,*,m] = lai_temp

           ENDFOR
           
        ENDIF
     ENDELSE

     IF scale_minmax THEN BEGIN

        ;; close holes in LAImin/LAImax
        idx=where(finite(laimin,/nan),count,ncomplement=count2)
        IF (count2 EQ 0) THEN BEGIN
           ;; no MODIS data available
           laimin[*] = 0.
           laimax[*] = 0.
        ENDIF ELSE BEGIN
           ;; MODIS data available, check for gaps
           IF (count GT 0) THEN BEGIN
              ;; gaps available, fill gaps
              laimin[idx]=nodata
              laimax[idx]=nodata
              
              gaps = where(laimin EQ nodata, gapcount)
              dia = 4
              WHILE (gapcount GT 0L) AND (dia LT 32) DO BEGIN
                 
                 print,'Trying to fill ',gapcount, $
                       ' holes in LAImin/LAImax with radius: ',dia
                 
                 laimin[gaps] = (fill_dist(laimin,dia,nodata))[gaps]
                 laimax[gaps] = (fill_dist(laimax,dia,nodata))[gaps]
                 
                 dia=dia*2
                 gaps = where(laimin EQ nodata, gapcount)
                 
              ENDWHILE
              
              ;; remaining gaps: probably water, fill with 0.
              IF gapcount GT 0L THEN laimin[gaps] = 0.
              IF gapcount GT 0L THEN laimax[gaps] = 0.
              
           ENDIF
        ENDELSE

        ;; scale LAI with max/min values
        FOR m = 0,nmonth-1 DO lai[*,*,m] = (((lai[*,*,m] - laimin)/((laimax - laimin)>0.01))>0.0)<1.0
     ENDIF

     
     ;; C4 is the sum of LAI for those months that satisfy the
     ;; C4 growth criteria (T>22C and p>25mm) over the sum of
     ;; LAI for all months
     
     lai = (lai-0.25)>0.0

     c4frac = (total(lai*((thighres GT 22.0) AND (phighres GT 25.0)),3,/nan) / $
               (total(lai,3,/nan)>0.1))<1.0
     
     ;; another potential fix for more consistency w. Still
     ;; et al. 2003 C4/C3 distribution: they used NDVI
     ;; instead of LAI
     c4frac = sqrt(c4frac)

     lai = 0
     
     time_laicalc += systime(1)
     
  ENDIF ELSE BEGIN
     ;; no C4 present in current patch
     c4frac = make_array(nlon,nlat,/float,value=0.)
  ENDELSE
  
  thighres = 0
  phighres = 0

  time_cropread -= systime(1)

;; read crop data (Leff et al.)
  crop_prefix = 'crops'
  crop_names = ['barley','cassava','cotton',$
                'grnuts','maize','millet', $
                'oilpalm','others','potat','pulses', $
                'rape','rice','rye','sorghum',$
                'soy','sugarb','sugarc','sunfl','wheat']
  ncrops = n_elements(crop_names)

  read_cropdata, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                 crop_prefix, crop_names, crops

  allcrops = make_array(nlon,nlat,ncrops,/float,value=!values.f_nan)
  FOR c=0,ncrops-1 DO BEGIN
     read_cropdata, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                    crop_prefix, crop_names[c], tempcrop
     allcrops[*,*,c] = tempcrop
  ENDFOR
  tempcrop = 0
  
  time_cropread += systime(1)
  ;;window,0,xsize=nlon*wsc,ysize=nlat*wsc
  ;;tv,bytscl(congrid(crops,nlon*wsc,nlat*wsc),min=0.,max=100)

  time_avhrrread -= systime(1)

;; read AVHRR VCF
  avhrr_prefix = 'avhrr/vcf'
  avhrr_names = ['gl-latlong-treecover.bin','gl-latlong-broadleaf.bin', $
                 'gl-latlong-needleleaf.bin', $
                 'gl-latlong-deciduous.bin','gl-latlong-evergreen.bin']

  read_avhrrvcf, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                 avhrr_prefix, avhrr_names[0], treecover
  read_avhrrvcf, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                 avhrr_prefix, avhrr_names[1], broadleaf
  read_avhrrvcf, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                 avhrr_prefix, avhrr_names[2], needleleaf
  read_avhrrvcf, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                 avhrr_prefix, avhrr_names[3], deciduous
  read_avhrrvcf, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                 avhrr_prefix, avhrr_names[4], evergreen
  
  time_avhrrread += systime(1)

;  window,1,xsize=nlon*wsc,ysize=nlat*wsc
;  tv,bytscl(congrid(treecover,nlon*wsc,nlat*wsc),min=0.,max=100)
;  tv,bytscl(congrid(deciduous+evergreen,nlon*wsc,nlat*wsc),min=0.,max=100)
;  tv,bytscl(congrid(treecover,nlon*wsc,nlat*wsc),min=0.,max=100)

  time_lcread -= systime(1)

;; read MOD12Q1 to determine grass/shrub fractions (from IGBP class)

  IGBP_version = 5

  read_modistile, basedir, 2004, 1, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                  'IGBP', obsvarunits, obsvarlongnames, IGBP_version, lc_data, lc_error, lc_mask, $
                  pollon=pollon, pollat=pollat,/initialize

;; IGBP classes of MOD12Q1:
;;ID Name                                  Shrub     Grass  	
;; 0 Water                                   -          -
;; 1 Evergreen Needleleaf Forest            90%        10%
;; 2 Evergreen Broadleaf Forest             90%        10%
;; 3 Deciduous Needleleaf Forest            90%        10%
;; 4 Deciduous Broadleaf Forest             90%        10%
;; 5 Mixed Forests                          90%        10%
;; 6 Closed Shrublands                      80%        20%
;; 7 Open Shrublands                        80%        20%
;; 8 Woody Savannas                         50%        50%
;; 9 Savannas                               10%        90%
;;10 Grasslands                              0%       100%
;;11 Permanent Wetlands                     50%        50%
;;12 Croplands                               0%       100%
;;13 Urban and Built-Up                      0%       100%
;;14 Cropland/Natural Vegetation Mosaic      0%       100%
;;15 Snow and Ice                            0%       100%
;;16 Barren or Sparsely Vegetated           10%        90%

;; biome arrays
  water12 = make_array(nlon,nlat,/float,value=!values.f_nan)
  shrub = make_array(nlon,nlat,/float,value=!values.f_nan)
  grass = make_array(nlon,nlat,/float,value=!values.f_nan)
;  urban = make_array(nlon,nlat,/float,value=!values.f_nan)
;  snowice = make_array(nlon,nlat,/float,value=!values.f_nan)

  shrubidx=make_array(256)
  grassidx=make_array(256)

  ;; shrub/grass fraction based on landcover: always amount to 100%:

;; new suggestion (2008) by P. Lawrence
;                    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
;  shrubidx[0:16]=[ 50., 90., 90., 90., 90., 90., 80., 80., 50., 10.,  0., 50.,  0.,  0.,  0.,  0., 10.]
;  grassidx[0:16]=[ 50., 10., 10., 10., 10., 10., 20., 20., 50., 90.,100., 50.,100.,100.,100.,100., 90.]
  
;; my suggestion:
  shrubidx[0:16]=[ 50., 50., 50., 50., 50., 50., 80., 50., 50., 10.,  0., 50.,  0.,  0.,  0.,  0., 10.]
  grassidx[0:16]=[ 50., 50., 50., 50., 50., 50., 20., 50., 50., 90.,100., 50.,100.,100.,100.,100., 90.]

  ;; LC 255 or NAN is interrupted space, set it to water
  baddata = where(finite(lc_data,/nan),badcount)
  IF badcount GT 0 THEN lc_data[baddata] = 0.       
  lc_data = byte(lc_data)
  baddata = where(lc_data EQ 255B,badcount)
  IF badcount GT 0 THEN lc_data[baddata] = 0B

  ;; regridding
  xyindex = yidmodis*nlon + xidmodis
  hist = histogram(xyindex,omin=omin,reverse_indices = R)

  FOR i=0L,n_elements(hist)-1L DO BEGIN
     IF R[i] NE R[i+1] THEN BEGIN
        ;; open and closed shrublands (IGBP 6/7)
        shrub[i+omin] = mean(shrubidx[lc_data[R[R[i]:R[i+1]-1]]])
        ;; from grasslands (IGBP 10)
        grass[i+omin] = mean(grassidx[lc_data[R[R[i]:R[i+1]-1]]])

        ;; water mask (IGBP = 0)
        water12[i+omin] = mean(lc_data[R[R[i]:R[i+1]-1]] EQ 0B)
        
        ;; from urban (IGBP 13) --> maybe exclude urban pixels
        ;; later (e.g. if more than x% urban then do not assimilate
;        urban[i+omin] = mean(lc_data[R[R[i]:R[i+1]-1]] EQ 13B)*100.

        ;; what to do with snow+ice and wetlands? Maybe exclude
        ;; them later
;        snowice[i+omin] = mean(lc_data[R[R[i]:R[i+1]-1]] EQ 15B)*100.
        ;;wetland[i+omin] = mean(lc_data[R[R[i]:R[i+1]-1]] EQ 11B)*100.
     ENDIF
  ENDFOR

  lc_data = 0
  lc_mask = 0
  lc_error = 0
  xmodis = 0
  ymodis = 0
  hmodis = 0
  vmodis = 0
  xidmodis = 0
  yidmodis = 0
  xyindex = 0
  hist = 0
  R = 0

  time_lcread += systime(1)

;; read MOD44 to separate tree, herbaceous and bare ground

;; Info from: http://edcdaac.usgs.gov/modis/mod44b.asp

;; Version 3:
;; 0-100 % cover
;; 251 bad data
;; 253 water
;; 255 interrupted space

;; Version 5:
;; 0-100 % cover
;; 200 water
;; 253 fill value

  time_vcfread -= systime(1)

  VCF_version = 5

  tree = make_array(nlon,nlat,/float,value=!values.f_nan)
  herb = make_array(nlon,nlat,/float,value=!values.f_nan)
  bare = make_array(nlon,nlat,/float,value=!values.f_nan)
  water44 = make_array(nlon,nlat,/float,value=!values.f_nan)

  ;; read tree cover fraction
  IF VCF_version EQ 3 THEN BEGIN
     read_modistile, basedir, 2000, 305, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'TREE', obsvarunits, obsvarlongnames, VCF_version, tree_data, tree_error, tree_mask, $
                     pollon=pollon, pollat=pollat,/initialize
     
     ;; read herbaceous cover fraction
     read_modistile, basedir, 2000, 305, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'HERB', obsvarunits, obsvarlongnames, VCF_version, herb_data, herb_error, herb_mask, $
                     pollon=pollon, pollat=pollat,/initialize
     
     ;; read bare soil fraction
     read_modistile, basedir, 2000, 305, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'BARE', obsvarunits, obsvarlongnames, VCF_version, bare_data, bare_error, bare_mask, $
                     pollon=pollon, pollat=pollat,/initialize

     ;; extract water mask (bare fraction values > 100)
     water_data = 100. * float(tree_data EQ 253.)

     ;; set missing values to NAN
     badvalues = where((tree_data GT 100.) OR (herb_data GT 100.) OR (bare_data GT 100.),badcount)
     IF badcount GT 0L THEN BEGIN
        tree_data[badvalues]=!values.f_nan
        herb_data[badvalues]=!values.f_nan
        bare_data[badvalues]=!values.f_nan
     ENDIF
     
     ;; regridding
     xyindex = yidmodis*nlon + xidmodis
     hist = histogram(xyindex,omin=omin,reverse_indices = R)
     FOR i=0L,n_elements(hist)-1L DO BEGIN
        IF R[i] NE R[i+1] THEN BEGIN
           IF total(finite(tree_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              tree[i+omin] = mean(tree_data[R[R[i]:R[i+1]-1]],/nan)
           IF total(finite(herb_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              herb[i+omin] = mean(herb_data[R[R[i]:R[i+1]-1]],/nan)
           IF total(finite(bare_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              bare[i+omin] = mean(bare_data[R[R[i]:R[i+1]-1]],/nan)
           IF total(finite(water_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              water44[i+omin] = mean(water_data[R[R[i]:R[i+1]-1]],/nan)
        ENDIF
     ENDFOR

  ENDIF ELSE BEGIN

     ;; collectiion 5 puts values outside of projection into the same category as
     ;; fill values, so read water mask to differentiate
     read_modistile, basedir, 2000, 55, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'WATER', obsvarunits, obsvarlongnames, 5, water_data, water_error, water_mask, $
                     pollon=pollon, pollat=pollat,/initialize

     ;; fill missing area (with no MODIS tiles) with water
     badvalues = where(finite(water_data,/nan),badcount)
     IF badcount GT 0L THEN water_data[badvalues] = 1.0
     ;; generate IEEE NAN missing values from missing MODIS flag
     badvalues = where(water_data EQ 253.,badcount)
     IF badcount GT 0L THEN water_data[badvalues] = !values.f_nan
     ;; set water to 100%
     goodvalues = where(finite(water_data),goodcount)
     IF goodcount GT 0L THEN water_data[goodvalues] = 100.*float(water_data[goodvalues] EQ 1.)

     ;; collection 5 only contains tree fraction, so still read
     ;; collection 3 herbal and bare fraction
     read_modistile, basedir, 2004, 65, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'TREE', obsvarunits, obsvarlongnames, 5, tree_data, tree_error, tree_mask, $
                     pollon=pollon, pollat=pollat,/initialize

     badvalues = where(tree_data GT 100.,badcount)
     IF badcount GT 0L THEN tree_data[badvalues] = !values.f_nan

     ;; regridding
     xyindex = yidmodis*nlon + xidmodis
     hist = histogram(xyindex,omin=omin,reverse_indices = R)
     FOR i=0L,n_elements(hist)-1L DO BEGIN
        IF R[i] NE R[i+1] THEN BEGIN
           IF total(finite(tree_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              tree[i+omin] = mean(tree_data[R[R[i]:R[i+1]-1]],/nan)
           IF total(finite(water_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              water44[i+omin] = mean(water_data[R[R[i]:R[i+1]-1]],/nan)
        ENDIF
     ENDFOR

     ;; fix for high latitude water mask reprojection problems (no valid pixels)
     badidx = where(finite(water44,/nan),badcount)
     dia = 4
     WHILE (badcount GT 0L) AND (dia LT 32) DO BEGIN
        water44[badidx] = nodata
        water44 = fill_dist(water44,dia,nodata)
        badidx = where(water44 EQ nodata,badcount)
        IF badcount GT 0L THEN water44[badidx] = !values.f_nan
        dia *= 2
     ENDWHILE

     ;; check for missing water mask values (none should remain)
     badidx = where(finite(water44,/nan),badcount)
     IF badcount GT 0L THEN stop,"Missing values in reprojected water mask found"
      
     ;; read herbaceous cover fraction
     read_modistile, basedir, 2000, 305, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'HERB', obsvarunits, obsvarlongnames, 3, herb_data, herb_error, herb_mask, $
                     pollon=pollon, pollat=pollat,/initialize

     badvalues = where(herb_data GT 100.,badcount)
     IF badcount GT 0L THEN herb_data[badvalues] = !values.f_nan
     
     ;; read bare soil fraction
     read_modistile, basedir, 2000, 305, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     'BARE', obsvarunits, obsvarlongnames, 3, bare_data, bare_error, bare_mask, $
                     pollon=pollon, pollat=pollat,/initialize

     badvalues = where(bare_data GT 100.,badcount)
     IF badcount GT 0L THEN bare_data[badvalues] = !values.f_nan

     ;; regridding
     xyindex = yidmodis*nlon + xidmodis
     hist = histogram(xyindex,omin=omin,reverse_indices = R)
     FOR i=0L,n_elements(hist)-1L DO BEGIN
        IF R[i] NE R[i+1] THEN BEGIN
           IF total(finite(herb_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              herb[i+omin] = mean(herb_data[R[R[i]:R[i+1]-1]],/nan)
           IF total(finite(bare_data[R[R[i]:R[i+1]-1]])) GT 0 THEN $
              bare[i+omin] = mean(bare_data[R[R[i]:R[i+1]-1]],/nan)
        ENDIF
     ENDFOR

  ENDELSE

  tree_data = 0
  herb_data = 0
  bare_data = 0
  water_data = 0
  tree_error = 0
  herb_error = 0
  bare_error = 0
  water_error = 0
  tree_mask = 0
  herb_mask = 0
  bare_mask = 0
  water_mask = 0
  xyindex = 0
  xmodis = 0
  ymodis = 0
  hmodis = 0
  vmodis = 0
  xidmodis = 0
  yidmodis = 0
  hist = 0
  R = 0
 
  ;; set missing data to bare soil (happens mostly in polar and mountain areas)
  badidx = where(finite(tree,/nan) OR finite(herb,/nan) OR finite(bare,/nan), badcount)
  IF badcount GT 0L THEN BEGIN
     tree[badidx] = 0.
     herb[badidx] = 0.
     bare[badidx] = 100.-water44[badidx]
  ENDIF

  ;; bare+herb+tree+water can be more than 100%
  ;; first, subtract from bare, then herb, then tree, then water
  sum = tree+herb+bare+water44
  more = where(sum GT 100., morecount)
  IF morecount NE 0L THEN bare[more] = (bare[more] - (sum[more]-100.))>0.
  sum = tree+herb+bare+water44
  more = where(sum GT 100., morecount)
  IF morecount NE 0L THEN herb[more] = (herb[more] - (sum[more]-100.))>0.
  sum = tree+herb+bare+water44
  more = where(sum GT 100., morecount)
  IF morecount NE 0L THEN tree[more] = (tree[more] - (sum[more]-100.))>0.
  more = where(sum GT 100., morecount)
  IF morecount NE 0L THEN water44[more] = (water44[more] - (sum[more]-100.))>0.
  sum = tree+herb+bare+water44
  more = where(sum GT 100., morecount)
  IF morecount NE 0L THEN stop,"tree+herb+bare+water greater than 100%"
  sum = 0
  
  ;; bare+herb+tree+water can be less than 100%
  ;; first, add bare, then herb, then tree, then water
  sum = tree+herb+bare+water44
  less = where(sum LT 100., lesscount)
  IF lesscount NE 0L THEN bare[less] = (bare[less] + (100.-sum[less]))<100.
  sum = tree+herb+bare+water44
  less = where(sum LT 100., lesscount)
  IF lesscount NE 0L THEN herb[less] = (herb[less] + (100.-sum[less]))<100.
  sum = tree+herb+bare+water44
  less = where(sum LT 100., lesscount)
  IF lesscount NE 0L THEN tree[less] = (tree[less] + (100.-sum[less]))<100.
  less = where(sum LT 100., lesscount)
  IF lesscount NE 0L THEN water44[less] = (water44[less] + (100.-sum[less]))<100.
  sum = tree+herb+bare+water44
  less = where(sum LT 100., lesscount)
  IF lesscount NE 0L THEN stop,"tree+herb+bare+water less than 100%"
  sum = 0

  time_vcfread += systime(1)

  ;; now determine cover fractions of PFT's

  ;; tree fractions: need to fill gaps where there is tree % but no
  ;; fraction of type broadleaf, deciduous, needleleaf, evergreen -->
  ;; create fill type map

  time_calc -= systime(1)

  gaps = where((tree GT 0.) AND (treecover EQ 0.), gapcount)

  ;; we only need to fill leaf types where treecover is 0 and tree
  ;; is larger than 0. Set all treecover of 0 to nodata temporarily
  IF gapcount GT 0L THEN BEGIN
     treecovert  = treecover
     broadleaft  = broadleaf
     deciduoust  = deciduous

     idx = where(treecovert EQ 0.,count)

     treecovert[idx]  = nodata
     broadleaft[idx]  = nodata
     deciduoust[idx]  = nodata
  ENDIF

  dia = 2
  WHILE (gapcount GT 0L) AND (dia LT 32) DO BEGIN

     print,'Trying to fill ',gapcount,' holes in treecover with radius: ',dia

     treecovert  = fill_dist(treecovert,dia,nodata)

     broadleaft  = fill_dist(broadleaft,dia,nodata)
     deciduoust  = fill_dist(deciduoust,dia,nodata)

     treecover[gaps]  = treecovert[gaps]

     broadleaf[gaps]  = broadleaft[gaps]
     deciduous[gaps]  = deciduoust[gaps]

     dia=dia*2
     gaps = where((tree GT 0.) AND ((treecover EQ 0.) OR (treecover EQ nodata)), gapcount)

  ENDWHILE

  treecovert = 0
  broadleaft = 0
  deciduoust = 0

  IF gapcount GT 0L THEN BEGIN
     ;; assign remaining (faulty?) tree cover to bare cover
     bare[gaps] = bare[gaps] + tree[gaps]
     tree[gaps] = 0.

     treecover[gaps] = 0.
     broadleaf[gaps] = 0.
     deciduous[gaps] = 0.
  ENDIF

  broadleafpct = broadleaf / (treecover>1.)
  deciduouspct = deciduous / (treecover>1.)
  needleleafpct = ((treecover - broadleaf) / (treecover>1.))>0.
  evergreenpct  = ((treecover - deciduous) / (treecover>1.))>0.

  treecover = 0
  broadleaf = 0
  needleleaf = 0
  deciduous = 0
  evergreen = 0

  dbf = make_array(nlon,nlat,/float,value=!values.f_nan)
  ebf = make_array(nlon,nlat,/float,value=!values.f_nan)
  dnf = make_array(nlon,nlat,/float,value=!values.f_nan)
  enf = make_array(nlon,nlat,/float,value=!values.f_nan)

  ;; deciduous
  di = where(deciduouspct EQ 0.,dc)
  IF dc GT 0 THEN BEGIN
     dbf[di] = 0.
     dnf[di] = 0.
  ENDIF

  ;; evergreen
  di = where(evergreenpct EQ 0.,dc)
  IF dc GT 0 THEN BEGIN
     ebf[di] = 0.
     enf[di] = 0.
  ENDIF

  ;; broadleaf deciduous and needleleaf deciduous
  di = where((deciduouspct GT 0.) AND (deciduouspct GE broadleafpct),dc)
  IF dc GT 0 THEN BEGIN
     dbf[di] = broadleafpct[di]*tree[di]
     dnf[di] = (deciduouspct[di] - broadleafpct[di])*tree[di]
  ENDIF
  
  di = where((deciduouspct GT 0.) AND (deciduouspct LT broadleafpct),dc)
  IF dc GT 0 THEN BEGIN
     dbf[di] = deciduouspct[di]*tree[di]
     dnf[di] = 0.
  ENDIF

  ;; broadleaf evergreen and needleleaf evergreen
  di = where((evergreenpct GT 0.) AND (evergreenpct GE needleleafpct),dc)
  IF dc GT 0 THEN BEGIN
     ebf[di] = (evergreenpct[di]-needleleafpct[di])*tree[di]
     enf[di] = needleleafpct[di]*tree[di]
  ENDIF

  di = where((evergreenpct GT 0.) AND (evergreenpct LT needleleafpct),dc)
  IF dc GT 0 THEN BEGIN
     ebf[di] = 0.
     enf[di] = evergreenpct[di]*tree[di]
  ENDIF
  
  evergreenpct = 0
  deciduouspct = 0
  broadleafpct = 0
  needleleafpct = 0

  ;; check for consistency, and add trees if needed in
  ;; places where the leaf type does not sum up to total
  ;; tree cover.
  di = where((round(dbf+dnf+ebf+enf) NE round(tree)) AND (dbf GE ebf) AND (dbf GE dnf) AND (dbf GE enf),dc)
  IF dc GT 0 THEN dbf[di] = tree[di] - (ebf[di]+dnf[di]+enf[di])

  di = where((round(dbf+dnf+ebf+enf) NE round(tree)) AND (ebf GE dbf) AND (ebf GE dnf) AND (ebf GE enf),dc)
  IF dc GT 0 THEN ebf[di] = tree[di] - (dbf[di]+dnf[di]+enf[di])

  di = where((round(dbf+dnf+ebf+enf) NE round(tree)) AND (dnf GE dbf) AND (dnf GE ebf) AND (dnf GE enf),dc)
  IF dc GT 0 THEN dnf[di] = tree[di] - (dbf[di]+ebf[di]+enf[di])

  di = where((round(dbf+dnf+ebf+enf) NE round(tree)) AND (enf GE dbf) AND (enf GE ebf) AND (enf GE dnf),dc)
  IF dc GT 0 THEN enf[di] = tree[di] - (dbf[di]+ebf[di]+dnf[di])

  ;; calculate herbaceous fractions: crop + grass + shrub
  
  cro = make_array(nlon,nlat,/float,value=0.)
  
  ;; limit crop fraction and herb / bare fraction
  ci = where((herb GT 0.) AND (crops GT (herb+bare)),count)
  IF count GT 0 THEN BEGIN
     cro[ci] = herb[ci]+bare[ci]
     herb[ci] = 0.
     bare[ci] = 0.
  ENDIF

  ci = where((herb GT 0.) AND (crops GT herb) AND (crops LE (herb+bare)),count)
  IF count GT 0 THEN BEGIN
     bare[ci] = herb[ci] + bare[ci] - crops[ci]
     herb[ci] = 0.
     cro[ci] = crops[ci]
  ENDIF

  ci = where((herb GT 0.) AND (crops LE herb),count)
  IF count GT 0 THEN BEGIN
     herb[ci] = herb[ci] - crops[ci]
     cro[ci] = crops[ci]
  ENDIF

  crops = 0

  grassshrub = grass+shrub
  gaps = where((herb GT 0.) AND (grassshrub EQ 0.), gapcount)

  IF gapcount GT 0L THEN BEGIN
     
     grassshrubt = grassshrub
     grasst = grass

     idx = where(grassshrubt EQ 0.,count)

     grassshrubt[idx] = nodata
     grasst[idx] = nodata
  ENDIF

  dia = 2
  WHILE (gapcount GT 0L) AND (dia LT 32) DO BEGIN

     print,'trying to fill holes in grass + shrub cover with radius ',dia,gapcount

     grassshrubt = fill_dist(grassshrubt,dia,nodata)
     grasst = fill_dist(grasst,dia,nodata)

     grassshrub[gaps] = grassshrubt[gaps]
     grass[gaps] = grasst[gaps]

     gaps = where((herb GT 0.) AND ((grassshrub EQ 0.) OR (grassshrub EQ nodata)), gapcount)
     dia = dia*2

  ENDWHILE
  
  grassshrubt = 0
  grasst = 0

  IF gapcount GT 0L THEN BEGIN
     ;; assign remaining herbaceous gaps to bare soil
     bare[gaps] = bare[gaps] + herb[gaps]
     herb[gaps] = 0.

     grass[gaps] = 0.
     grassshrub[gaps] = 0.
  ENDIF

  ;; divide herbaceous (non-crop) into grass and shrub
  shr = make_array(nlon,nlat,/float,value=0.)
  gra = make_array(nlon,nlat,/float,value=0.)

  grasspct = grass/((grassshrub)>1.)
  shrubpct = ((grassshrub-grass)/((grassshrub)>1.))>0.

  idx = where(shrubpct GT 0.,count,complement=idx2,ncomplement=count2)
  IF count GT 0 THEN BEGIN
     shr[idx] = herb[idx] * shrubpct[idx]
     gra[idx] = (1.-shrubpct[idx]) * herb[idx]
  ENDIF
  IF count2 GT 0 THEN BEGIN
     shr[idx2] = 0.
     gra[idx2] = herb[idx2]
  ENDIF
  
  grassshrub = 0
  grass = 0
  shrub = 0
  shrubpct = 0
  grasspct = 0

  IF VCF_version EQ 5 THEN BEGIN
     water = water44
  ENDIF ELSE BEGIN
     ;; now we have areas that are water in MODIS but land in AVHRR or
     ;; vice versa: water44 is taken as true land/ocean %area here
     ;; take MOD44 landmask for most areas, except for antarctica and for
     ;; >80N, where only MOD12 has data
     water = water44 * ((lat GE -60.) AND (lat LT 80.)) + water12 * ((lat LT -60.) OR (lat GE 80.))
     ;; 80S-90S all is land (bare+snow and ice) (MOD12Q1 also has errors)
     water = water * (lat GE -80.) + 1.*(lat LT -80.)
  ENDELSE

  water44 = 0
  water12 = 0

  ncvar = bytarr(nlon,nlat,nvar)

  ;; now we apply the climate rules to distribute between
  ;; arctic, boreal, temperate and tropical PFT's

  ;; set bare fraction to 100% for antarctica
  ncvar[*,*,0]= round(bare)

  ncvar[*,*,1] = round(enf) * ((tmin GT -19.0) AND (gdd GT 1200.0))
  ncvar[*,*,2] = round(enf) * ((tmin LE -19.0) OR (gdd LE 1200.0))

  ncvar[*,*,3] = round(dnf)

  ncvar[*,*,4] = round(ebf) * (tmin GT 15.5)
  ncvar[*,*,5] = round(ebf) * (tmin LE 15.5)

  ncvar[*,*,6] = round(dbf) * (tmin GT 15.5)
  ncvar[*,*,7] = round(dbf) * ((tmin GT -15.0) AND (tmin LE 15.5) AND (gdd GT 1200.0))
  ncvar[*,*,8] = round(dbf) * ((tmin LE -15.0) OR (gdd LE 1200.0))

  ncvar[*,*,9] = round(shr) * ((tmin GT -19.0) AND (gdd GT 1200.0) AND $
                               (pann GT 520.0) AND (pwin GT 2./3.*pann))
  ncvar[*,*,10] = round(shr) * ((tmin GT -19.0) AND (gdd GT 1200.0) AND $
                                ((pann LE 520.0) OR (pwin LE 2./3.*pann)))
  ncvar[*,*,11]= round(shr) * ((tmin LE -19.0) OR (gdd LE 1200.0))

  ncvar[*,*,12]= round(gra) * (gdd LE 1000.0)

  ncvar[*,*,13]= round(gra) * ((gdd GT 1000.0) AND ((tmax LE 22.0) OR (pmaxtg22 LE 25.0)) $
                               AND ((tmin LE 22.0) OR (pmin LE 25.0)))
  ncvar[*,*,14]= round(gra) * ((gdd GT 1000.0) AND (tmin GT 22.0) AND (pmin GT 25.0))
  
  ;; distribute grasses according to LAI-derived c4-fraction if neither of the above criteria is met
  leftover = round(gra) - round(total(float(ncvar[*,*,12:14]),3))
  ncvar[*,*,13]= ncvar[*,*,13] + (1.-c4frac)*leftover
  ncvar[*,*,14]= ncvar[*,*,14] + c4frac*leftover
  
  c1 = 15
  c2 = ncrops+15-1
  
  totalcrops = total(allcrops,3)>0.01
  FOR v=c1,c2 DO BEGIN
     ncvar[*,*,v] = round(allcrops[*,*,v-c1]/totalcrops * cro)
  ENDFOR
  
  ncvar[*,*,34]= round(water)

  totalcrops = 0
  allcrops = 0
  c4frac = 0

  bare = 0
  tree = 0
  herb = 0

  enf = 0
  dnf = 0
  ebf = 0
  dbf = 0
  shr = 0
  gra = 0
  cro = 0

  water = 0
  
  IF verbose THEN BEGIN
     FOR v=0,nvar-1 DO print,ncvarnames[v],'% :',mean(ncvar[*,*,v],/nan)
     print,'Minimum Total land+water %: ',min(total(ncvar,3))
     print,'Maximum Total land+water %: ',max(total(ncvar,3))
  ENDIF

  ;; check for >100% total everywhere
  more = where(total(ncvar,3) GT 100,morecount)
  WHILE (morecount GT 0L) DO BEGIN
     v = fix(randomu(seed)*nvar)<(nvar-1)
     temp = fix(ncvar[*,*,v])
     temp[more] = (temp[more] - 1*(temp[more] NE 100))>0
     ncvar[*,*,v] = byte(temp)
     more = where(total(ncvar,3) GT 100,morecount)
  ENDWHILE

  ;; check for <100% total everywhere
  less = where(total(ncvar,3) LT 100,lesscount)
  WHILE (lesscount GT 0L) DO BEGIN
     v = fix(randomu(seed)*nvar)<(nvar-1)
     temp = fix(ncvar[*,*,v])
     temp[less] = (temp[less] + 1*(temp[less] NE 0))<100
     ncvar[*,*,v] = byte(temp)
     less = where(total(ncvar,3) LT 100,lesscount)
  ENDWHILE
  
  IF verbose THEN BEGIN
     FOR v=0,nvar-1 DO print,ncvarnames[v],'% :',mean(ncvar[*,*,v],/nan)
     print,'Minimum Total land+water %: ',min(total(ncvar,3))
     print,'Maximum Total land+water %: ',max(total(ncvar,3))
  ENDIF
  
  time_calc += systime(1)

  ;; --------- NETCDF WRITE DATA SECTION ---------------

  time_netcdf -= systime(1)

  ;; Now write the data to the NetCDF file(s)

  FOR file = 0, nfiles-1 DO BEGIN

     ;; only open NetCDF file[s] if they are not in use by another process
     file_lock,pftdir,pftfile[file],/lock,id=100L*xproc+yproc
     
     ;; inquire id's in existing file and check if healthy
     ncid[file] = ncdf_open(pftdir+pftfile[file],/write)
     
     latid[file] = ncdf_varid(ncid[file],'lat')
     lonid[file] = ncdf_varid(ncid[file],'lon')
     IF manyfiles THEN BEGIN
        v = file
        ncvarid[v] = ncdf_varid(ncid[file],ncvarnames[v])
     ENDIF ELSE BEGIN
        FOR v = 0,nvar-1 DO ncvarid[v] = ncdf_varid(ncid[file],ncvarnames[v])
     ENDELSE
     
     ;; write longitude, latitude
     NCDF_VARPUT, ncid[file], lonid[file], reform(lon[*,0]), offset=[tx0], count=[nlon]
     NCDF_VARPUT, ncid[file], latid[file], reform(lat[0,*]), offset=[ty0], count=[nlat]
     
     ;; write monthly data
     IF manyfiles THEN BEGIN
        v = file
        NCDF_VARPUT, ncid[file], ncvarid[v],reform(ncvar[*,*,v]), offset=[tx0,ty0],count=[nlon,nlat]
     ENDIF ELSE BEGIN
        FOR v = 0, nvar-1 DO BEGIN
           NCDF_VARPUT, ncid[file], ncvarid[v],reform(ncvar[*,*,v]), offset=[tx0,ty0],count=[nlon,nlat]
        ENDFOR
     ENDELSE
     
     ;; close NetCDF file
     NCDF_CLOSE, ncid[file] 
     
     ;; remove lock
     file_lock,pftdir,pftfile[file],/unlock,id=100L*xproc+yproc
     
  ENDFOR

  time_netcdf += systime(1)

  time_total += systime(1)

  print
  print,"Total wall clock by process: "
  print,'netcdf r/w   : ',time_netcdf
  print,'climate read : ',time_climread
  print,'climate calc : ',time_climcalc
  print,'lai read     : ',time_lairead
  print,'lai calc     : ',time_laicalc
  print,'crop read    : ',time_cropread
  print,'avhrr read   : ',time_avhrrread
  print,'lc read      : ',time_lcread
  print,'vcf read     : ',time_vcfread
  print,'final calc   : ',time_calc
  print,'TOTAL        : ',time_total

END
