PRO regridpft

;; IDL Code by Reto Stockli, Blue Marble Research, 2011/08/16

;; This script reprojects the global pft dataset from a regular
;; cylindrical coordinate system to
;; a) a rotated pole coordinate
;; b) smaller dimensions
;; system. It currently uses bilinear interpolation. 

  !except = 2

;; flags
  manyfilesin = 1
  manyfilesout = 0
  climatepft = 1
  readallcrops = 1
  rotated = 0
  
;; version
  version = "v1.67"


;; directories
  basedir = '/project/msclim/stockli'
;;  basedir = '/Users/stockli'
  indir = basedir + '/pft/'
  outdir = indir
  inregion = 'global'
  outregion = 'global_0.5'

;; output grid 
  IF rotated THEN BEGIN
     ;; rotated pole grid
     olonmin = -17.0
     olonmax = 8.0
     olatmin = -10.0
     olatmax = 11.0

     pollon = -170.0
     pollat = 43.0
     olonlen = 120
     olatlen = 120
  ENDIF ELSE BEGIN
     ;; regular cylindrical
     olonmin = -180.0
     olonmax = 180.0
     olatmin = -90.0
     olatmax = 90.0
     dolon = 0.5
     dolat = 0.5
  ENDELSE

;; netcdf names

  IF climatepft THEN BEGIN

     IF readallcrops THEN BEGIN

        nvar = 17+19-1
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
                      'cro_sug',$
                      'cro_sgb',$
                      'cro_sgc',$
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

     ENDIF ELSE BEGIN

        nvar = 17
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
                      'cro_all',$
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
                          'Evergreen Broadleaf Shrubs (-)',$
                          'Deciduous Broadleaf Shrubs (Temperate)',$
                          'Deciduous Broadleaf Shrubs (Boreal)',$
                          'C3 Grass (Arctic)',$
                          'C3 Grass (Non-Arctic)',$
                          'C4 Grass (-)',$
                          'Croplands (-)',$
                          'Water (-)']
        
        ncvarunits = replicate('%',nvar)

     ENDELSE

  ENDIF ELSE BEGIN
     
     nvar = 9
     ncvarnames = ['bar','enf','ebf','dnf','dbf','shr','gra','cro','wat']
     ncvarlongnames = ['Bare Soil', $
                       'Evergreen Needleleaf Forest', $
                       'Evergreen Broadleaf Forest', $
                       'Deciduous Needleleaf Forest', $
                       'Deciduous Broadleaf Forest', $
                       'Shrubs', $
                       'Grasslands', $
                       'Croplands', $
                       'Water']
     ncvarunits = replicate('%',nvar)

  ENDELSE

;; -------- START MAIN CODE ---------

  ;; calculate current time, format for NetCDF output
  caldat,systime(/julian,/utc),prod_month,prod_day,prod_year,prod_hour,prod_minute,prod_second
  
  prod_date = string(prod_year,format='(I4)')+'-'+string(prod_month,format='(I2.2)')+'-'+ $
              string(prod_day,format='(I2.2)')+' '+string(prod_hour,format='(I2.2)')+':'+ $
              string(prod_minute,format='(I2.2)')+':'+string(round(prod_second),format='(I2.2)')

  IF rotated THEN BEGIN
     ;; define output resolution (degrees per pixel)
     dolon = 1d0/double(olonlen)
     dolat = 1d0/double(olatlen)
     
     nolon = round((olonmax - olonmin)*double(olonlen))
     nolat = round((olatmax - olatmin)*double(olatlen))
     
     ;; reproject these arrays to geographical olon/olat pairs
     RPOL2LL,OLON,OLAT,LON,LAT,pollon=pollon,pollat=pollat

     ;; define boundaries which we will read from input topography dataset
     olonmin = min(lon)-0.5
     olonmax = max(lon)+0.5
     olatmin = min(lat)-0.5
     olatmax = max(lat)+0.5

  ENDIF ELSE BEGIN

     nolon = round((olonmax - olonmin)/dolon)
     nolat = round((olatmax - olatmin)/dolat)

     ncvar = bytarr(nolon,nolat,nvar)

  ENDELSE

  ;; output lon/lat arrays (centered on pixel)
  olon = ((dindgen(nolon) + 0.5d0)*dolon + double(olonmin)) # replicate(1.d0,nolat)    ; W->E
  olat = replicate(1.d0,nolon) # ((dindgen(nolat) + 0.5d0)*dolat + double(olatmin))    ; S->N
  
  ;; define output array
  ncvar = bytarr(nolon,nolat,nvar)

  ;; define file names
  IF manyfilesin THEN BEGIN
     nfilesin = nvar
     infile=ncvarnames+'.'+inregion+'.nc' 
  ENDIF ELSE BEGIN 
     nfilesin = 1
     infile='pft.'+inregion+'.nc'
  ENDELSE

  IF manyfilesout THEN BEGIN
     nfilesout = nvar
     outfile=ncvarnames+'.'+outregion+'.nc' 
  ENDIF ELSE BEGIN 
     nfilesout = 1
     outfile='pft.'+outregion+'.nc'
  ENDELSE

  ;; read original NetCDF file

  ;; inquire geographic domain of input file
  ncid = ncdf_open(indir+infile[0])
  ncdf_varget,ncid,'lon',lons_pft
  ncdf_varget,ncid,'lat',lats_pft

  dlon_pft = abs(lons_pft[1]-lons_pft[0])
  dlat_pft = abs(lats_pft[1]-lats_pft[0])
  
  latmin_pft = min(lats_pft) - 0.5d0*dlat_pft
  lonmin_pft = min(lons_pft) - 0.5d0*dlon_pft
  latmax_pft = max(lats_pft) + 0.5d0*dlat_pft
  lonmax_pft = max(lons_pft) + 0.5d0*dlon_pft
  
  nlon_pft = n_elements(lons_pft)
  nlat_pft = n_elements(lats_pft)
  
  xmin_pft = -1L
  xmax_pft = -1L
  ymin_pft = -1L
  ymax_pft = -1L
  
  ;; define precision for getting rid of rounding errors below
  pr = 1.d-6

  FOR i=nlon_pft-1L,0L,-1L DO IF (round((lons_pft[i]-0.5d0*dlon_pft)/pr,/L64) LE round(olonmin/pr,/L64)) $
     AND (xmin_pft EQ -1L) THEN xmin_pft = i
  FOR i=nlat_pft-1L,0L,-1L DO IF (round((lats_pft[i]-0.5d0*dlat_pft)/pr,/L64) LE round(olatmin/pr,/L64)) $
     AND (ymin_pft EQ -1L) THEN ymin_pft = i
  FOR i=0L,nlon_pft-1L,1L DO IF (round((lons_pft[i]+0.5d0*dlon_pft)/pr,/L64) GE round(olonmax/pr,/L64)) $
     AND (xmax_pft EQ -1L) THEN xmax_pft = i
  FOR i=0L,nlat_pft-1L,1L DO IF (round((lats_pft[i]+0.5d0*dlat_pft)/pr,/L64) GE round(olatmax/pr,/L64)) $
     AND (ymax_pft EQ -1L) THEN ymax_pft = i
  

  IF ((xmin_pft EQ (-1L)) OR (ymin_pft EQ (-1L)) OR (xmax_pft EQ (-1L)) OR (ymax_pft EQ (-1L))) THEN BEGIN
     print,'pft dataset does not cover requested area. stopping.'
     print,'We need: ',olatmin,olatmax,olonmin,olonmax
     print,'We get:  ',latmin_pft,latmax_pft,lonmin_pft,lonmax_pft
     stop
  ENDIF

  nx_pft = xmax_pft - xmin_pft + 1L
  ny_pft = ymax_pft - ymin_pft + 1L

  ncdf_close,ncid

  ;; reprojection indices
  IF rotated THEN BEGIN
     lons_data = lons_pft[xmin_pft:xmax_pft]
     lats_data = lats_pft[ymin_pft:ymax_pft]
     
     lonidx = (lon - lons_data[0])/dlon_pft
     latidx = (lat - lats_data[0])/dlat_pft
  ENDIF

  data = bytarr(nx_pft,ny_pft)
 
  ;; now read pft's
  FOR v=0,nvar-1 DO BEGIN
;  FOR v=0,0 DO BEGIN

     print,'Reprojecting variable ',v
        
     IF manyfilesin THEN BEGIN
        ncid = ncdf_open(indir+infile[v])
     ENDIF ELSE BEGIN
        ncid = ncdf_open(indir+infile[0])
     ENDELSE
     
     ncdf_varget,ncid, ncvarnames[v], data, offset = [xmin_pft,ymin_pft], count = [nx_pft,ny_pft]      
     ncdf_close,ncid
     
     IF rotated THEN BEGIN
        ncvar[*,*,v] = round(bilinear(float(data),lonidx,latidx))
     ENDIF ELSE BEGIN
        ncvar[*,*,v] = round(rebin(float(data),nolon,nolat))
     ENDELSE
     
  ENDFOR

;  goto,bla

  ;; check for >100% total everywhere
  more = where(total(ncvar,3) GT 100,morecount)
  WHILE (morecount GT 0) DO BEGIN
     print,'Reducing total pft to 100 for ',morecount,' points'
     v = fix(randomu(seed)*nvar)
     temp = fix(ncvar[*,*,v])
     temp[more] = (temp[more] - 1*(temp[more] NE 100.))>0
     ncvar[*,*,v] = byte(temp)
     more = where(total(ncvar,3) GT 100,morecount)
  ENDWHILE
  
  ;; check for <100% total everywhere
  less = where(total(ncvar,3) LT 100,lesscount)
  WHILE (lesscount GT 0) DO BEGIN
     print,'Augmenting total pft to 100 for ',lesscount,' points'
     v = fix(randomu(seed)*nvar)
     temp = fix(ncvar[*,*,v])
     temp[less] = (temp[less] + 1*(temp[less] NE 0.))<100
     ncvar[*,*,v] = byte(temp)
     less = where(total(ncvar,3) LT 100,lesscount)
  ENDWHILE
  
;  bla:

  FOR v=0,nvar-1 DO print,ncvarnames[v],'% :',mean(ncvar[*,*,v],/nan)
  print,'Minimum Total land+water (final map) %: ',min(total(ncvar,3))

  ;; Write NetCDF file(s)

  ;; ---------- BEGIN NETCDF DEFINE ----------

  ncid = lonarr(nfilesout)
  olonid = lonarr(nfilesout)
  olatid = lonarr(nfilesout)
  ncvarid = lonarr(nvar)

  FOR file = 0, nfilesout-1 DO BEGIN
     ;; create NetCDF file (header and dimensions etc.)
     ;; if xchunk or ychunk are 0 independently if it exists.

     ncid[file] = NCDF_CREATE(outdir+outfile[file],/CLOBBER)
     NCDF_CONTROL,ncid[file],/NOFILL
     IF rotated THEN BEGIN
        olondim = NCDF_DIMDEF(ncid[file], 'rlon',nolon)
        olatdim = NCDF_DIMDEF(ncid[file], 'rlat',nolat)

        olonid[file] = NCDF_VARDEF(ncid[file], 'rlon', [olondim], /DOUBLE)
        NCDF_ATTPUT,ncid[file],olonid[file],'long_name','Longitude in rotated pole grid'
        NCDF_ATTPUT,ncid[file],olonid[file],'units','degrees'
        
        olatid[file] = NCDF_VARDEF(ncid[file], 'rlat', [olatdim], /DOUBLE)
        NCDF_ATTPUT,ncid[file],olatid[file],'long_name','Latitude in rotated pole grid'
        NCDF_ATTPUT,ncid[file],olatid[file],'units','degrees'
     
        poleid = NCDF_VARDEF(ncid[file], 'rotated_pole', /CHAR)
        NCDF_ATTPUT,ncid[file],poleid,'grid_mapping_name','rotated_latitude_longitude'
        NCDF_ATTPUT,ncid[file],poleid,'grid_north_pole_latitude',pollat
        NCDF_ATTPUT,ncid[file],poleid,'grid_north_pole_longitude',pollon
     ENDIF ELSE BEGIN
        olondim = NCDF_DIMDEF(ncid[file], 'lon',nolon)
        olatdim = NCDF_DIMDEF(ncid[file], 'lat',nolat)

        olonid[file] = NCDF_VARDEF(ncid[file], 'lon', [olondim], /DOUBLE)
        NCDF_ATTPUT,ncid[file],olonid[file],'long_name','Longitude'
        NCDF_ATTPUT,ncid[file],olonid[file],'units','degrees_east'
        
        olatid[file] = NCDF_VARDEF(ncid[file], 'lat', [olatdim], /DOUBLE)
        NCDF_ATTPUT,ncid[file],olatid[file],'long_name','Latitude'
        NCDF_ATTPUT,ncid[file],olatid[file],'units','degrees_north'
     ENDELSE
     

;;     missing = 255B ;; not used yet
     missing = 128B ;; not used yet
     scale = 1.0
     offset = 0.0
     
     ;; create data variables
     IF manyfilesout THEN BEGIN
        v=file
        
        ;; 2D data variable
        ncvarid[v] = NCDF_VARDEF(ncid[file],ncvarnames[v], [olondim,olatdim], /BYTE)

        ;; Attributes for data variables
        NCDF_ATTPUT, ncid[file], ncvarid[v], '_FillValue',missing
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'valid_range',[0,100]           
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'axis','YX'                              
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'long_name',ncvarlongnames[v]
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'units', ncvarunits[v]
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'scale_factor',scale
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'add_offset',offset
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'version', version
        NCDF_ATTPUT, ncid[file], ncvarid[v], 'prod_date', prod_date
         
     ENDIF ELSE BEGIN
        FOR v = 0,nvar-1 DO BEGIN
           
           ;; 2D data variable
           ncvarid[v] = NCDF_VARDEF(ncid[file],ncvarnames[v], [olondim,olatdim], /BYTE)

           ;; Attributes for data variables
           NCDF_ATTPUT, ncid[file], ncvarid[v], '_FillValue',missing
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'valid_range',[0,100]
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'axis','YX'                              
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'long_name',ncvarlongnames[v]
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'units', ncvarunits[v]
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'scale_factor',scale
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'add_offset',offset
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'version', version
           NCDF_ATTPUT, ncid[file], ncvarid[v], 'prod_date', prod_date
           
        ENDFOR
     ENDELSE
     
     ;; create global attributes
     IF rotated THEN BEGIN
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Title','Fractional biome distribution on rotated pole grid'
     ENDIF ELSE BEGIN
        NCDF_ATTPUT, ncid[file], /GLOBAL, 'Title','Fractional biome distribution'
     ENDELSE

     NCDF_ATTPUT, ncid[file], /GLOBAL, 'Institution','Blue Marble Research'
     NCDF_ATTPUT, ncid[file], /GLOBAL, 'Source','MOD44B, MOD12Q, AVHRR VCF, Crop database'
     NCDF_ATTPUT, ncid[file], /GLOBAL, 'Author','Reto Stockli'
     NCDF_ATTPUT, ncid[file], /GLOBAL, 'References','Reference Publication: P.J. Lawrence et al. (2007). Representing a new MODIS consistent land surface in the Community Land Model (CLM 3.0). JGR, Vol. 112, G01023'
     NCDF_ATTPUT, ncid[file], /GLOBAL, 'Conventions','CF-1.4'
   
     NCDF_CONTROL, ncid[file], /ENDEF ; Put file in data mode. 
     
     ;; close NetCDF file
     NCDF_CLOSE, ncid[file] 
     
     ;; inquire id's in existing file and check if healthy
     ncid[file] = ncdf_open(outdir+outfile[file],/write)
     
     IF rotated THEN BEGIN
        olatid[file] = ncdf_varid(ncid[file],'rlat')
        olonid[file] = ncdf_varid(ncid[file],'rlon')
     ENDIF ELSE BEGIN
        olatid[file] = ncdf_varid(ncid[file],'lat')
        olonid[file] = ncdf_varid(ncid[file],'lon')
     ENDELSE

     IF manyfilesout THEN BEGIN
        v = file
        ncvarid[v] = ncdf_varid(ncid[file],ncvarnames[v])
     ENDIF ELSE BEGIN
        FOR v = 0,nvar-1 DO ncvarid[v] = ncdf_varid(ncid[file],ncvarnames[v])
     ENDELSE
     
     ;; write longitude, latitude
     IF rotated THEN BEGIN
        NCDF_VARPUT, ncid[file], olonid[file], reform(olon[*,0]), offset=[0], count=[nolon]
        NCDF_VARPUT, ncid[file], olatid[file], reform(olat[0,*]), offset=[0], count=[nolat]
     ENDIF ELSE BEGIN
        NCDF_VARPUT, ncid[file], olonid[file], reform(olon[*,0]), offset=[0], count=[nolon]
        NCDF_VARPUT, ncid[file], olatid[file], reform(olat[0,*]), offset=[0], count=[nolat]
     ENDELSE
     
     ;; write monthly data
     IF manyfilesout THEN BEGIN
        v = file
        NCDF_VARPUT, ncid[file], ncvarid[v],reform(ncvar[*,*,v]), offset=[0,0],count=[nolon,nolat]
     ENDIF ELSE BEGIN
        FOR v = 0, nvar-1 DO BEGIN
           NCDF_VARPUT, ncid[file], ncvarid[v],reform(ncvar[*,*,v]), offset=[0,0],count=[nolon,nolat]
        ENDFOR
     ENDELSE
     
     ;; close NetCDF file
     NCDF_CLOSE, ncid[file] 
        
  ENDFOR

  stop

END
