PRO read_pftdata, basedir, pftdir, pftregion, lon, lat, pft, $
                   dlon=dlon, dlat=dlat, pftnames=pftnames
;; Read NetCDF files containing PFT maps
;; Reto Stockli (Blue Marble Research), 2009/05/01

;; define input files

;; 35 PFT's
  ncvarnames = [$
               'bar_all',$
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

  nvar = n_elements(ncvarnames)

  IF keyword_set(dlon) AND keyword_set(dlat) THEN regular=1 ELSE regular=0

  ;; check if we have a single file or multiple files
  manyfiles = 0
  pattern = basedir[0] + '/'+pftdir[0]+'/pft.'+pftregion+'.nc'
  infiles = file_search(pattern)
  IF infiles[0] EQ '' THEN BEGIN
     manyfiles = 1

     ;; and find number of pft's
     pattern = basedir[0] + '/'+pftdir[0]+'/*.'+pftregion+'.nc'
     infiles = file_search(pattern)

     IF n_elements(infiles) NE nvar THEN stop

     ;; rearrange infiles to account for pft name sequence
     infiles =  basedir[0] + '/'+pftdir[0]+'/'+ncvarnames+'.'+pftregion+'.nc'
  ENDIF ELSE BEGIN
     ncid = ncdf_open(infiles[0])
     ;; derive number of pft's
     info = ncdf_inquire(ncid)
     ncdf_close,ncid  
  
     IF (info.nvars - info.ndims) NE nvar THEN stop

  ENDELSE

  pftnames=ncvarnames

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

;; inquire geographic domain of input file
  ncid = ncdf_open(infiles[0])

  varid = ncdf_varid(ncid,'lon')
  IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlon')
  ncdf_varget,ncid,varid,lons_pft

  varid = ncdf_varid(ncid,'lat')
  IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlat')
  ncdf_varget,ncid,varid,lats_pft

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

  FOR i=nlon_pft-1L,0L,-1L DO IF (round((lons_pft[i]-0.5d0*dlon_pft)/pr,/L64) LE round(lonmin/pr,/L64)) $
     AND (xmin_pft EQ -1L) THEN xmin_pft = i
  FOR i=nlat_pft-1L,0L,-1L DO IF (round((lats_pft[i]-0.5d0*dlat_pft)/pr,/L64) LE round(latmin/pr,/L64)) $
     AND (ymin_pft EQ -1L) THEN ymin_pft = i
  FOR i=0L,nlon_pft-1L,1L DO IF (round((lons_pft[i]+0.5d0*dlon_pft)/pr,/L64) GE round(lonmax/pr,/L64)) $
     AND (xmax_pft EQ -1L) THEN xmax_pft = i
  FOR i=0L,nlat_pft-1L,1L DO IF (round((lats_pft[i]+0.5d0*dlat_pft)/pr,/L64) GE round(latmax/pr,/L64)) $
     AND (ymax_pft EQ -1L) THEN ymax_pft = i

  IF ((xmin_pft EQ (-1L)) OR (ymin_pft EQ (-1L)) OR (xmax_pft EQ (-1L)) OR (ymax_pft EQ (-1L))) THEN BEGIN
     print,'pft dataset does not cover requested area. stopping.'
     print,'We need: ',latmin,latmax,lonmin,lonmax
     print,'We get:  ',latmin_pft,latmax_pft,lonmin_pft,lonmax_pft
     print,xmin_pft,xmax_pft,ymin_pft,ymax_pft
     stop
  ENDIF

  nx_pft = xmax_pft - xmin_pft + 1L
  ny_pft = ymax_pft - ymin_pft + 1L

  IF (NOT regular) THEN BEGIN
     ;; adjust for boundary conditions in bilinear interpolation
     lontmp = (lon>(lonmin_pft+0.5d0*dlon_pft))<(lonmax_pft-0.5d0*dlon_pft)
     lattmp = (lat>(latmin_pft+0.5d0*dlat_pft))<(latmax_pft-0.5d0*dlat_pft)

     ;; derive bilinear interpolation indices
     x_pft = (lontmp - lons_pft[xmin_pft])/dlon_pft
     y_pft = (lattmp - lats_pft[ymin_pft])/dlat_pft

  ENDIF

  ncdf_close,ncid
;; read pft's

  pft = fltarr(nlon,nlat,nvar)

  IF manyfiles THEN BEGIN

     FOR v=0,nvar-1 DO BEGIN

        ncid = ncdf_open(infiles[v])
        ncdf_varget,ncid, ncvarnames[v], data, offset = [xmin_pft,ymin_pft], count = [nx_pft,ny_pft]      
        ncdf_close,ncid

        IF regular THEN BEGIN
           ;; scale pft's to output dimension if needed
           IF (nx_pft ne nlon) OR (ny_pft ne nlat) THEN BEGIN
              
              IF (nx_pft EQ (nx_pft/nlon)*nlon) AND (ny_pft EQ (ny_pft/nlat)*nlat) THEN BEGIN
                 ;; simple rebinning
                 pft[*,*,v] = rebin(float(data),nlon,nlat)
              ENDIF ELSE BEGIN
                 IF (nlon*nx_pft GT 20000) OR (nlat*ny_pft GT 20000) THEN BEGIN
                    ;; NN rebinning
                    pft[*,*,v] = congrid(float(data),nlon,nlat)
                 ENDIF ELSE BEGIN
                    ;; find common denominator, then do magnifying
                    pft[*,*,v] = rebin(congrid(float(data),nlon*nx_pft,nlat*ny_pft),nlon,nlat)
                 ENDELSE
              ENDELSE
              
           ENDIF ELSE BEGIN
              
              pft[*,*,v] = float(data)
              
           ENDELSE

        ENDIF ELSE BEGIN
           ;; use bilinear interpolation for irregular grid
           pft[*,*,v] = interpolate(float(data),x_pft,y_pft)
        ENDELSE

     ENDFOR

  ENDIF ELSE BEGIN

     ncid = ncdf_open(infiles[0])
     
     FOR v=0,nvar-1 DO BEGIN
        
        ncdf_varget,ncid, ncvarnames[v], data, offset = [xmin_pft,ymin_pft], count = [nx_pft,ny_pft]      
        
        IF regular THEN BEGIN
           ;; scale pft's to output dimension if needed
           IF (nx_pft ne nlon) OR (ny_pft ne nlat) THEN BEGIN
              
              IF (nx_pft EQ (nx_pft/nlon)*nlon) AND (ny_pft EQ (ny_pft/nlat)*nlat) THEN BEGIN
                 ;; simple rebinning
                 pft[*,*,v] = rebin(float(data),nlon,nlat)
              ENDIF ELSE BEGIN
                 IF (nlon*nx_pft GT 20000) OR (nlat*ny_pft GT 20000) THEN BEGIN
                    ;; NN rebinning
                    pft[*,*,v] = congrid(float(data),nlon,nlat)
                 ENDIF ELSE BEGIN
                    ;; find common denominator, then do magnifying
                    pft[*,*,v] = rebin(congrid(float(data),nlon*nx_pft,nlat*ny_pft),nlon,nlat)
                 ENDELSE
              ENDELSE
              
           ENDIF ELSE BEGIN
              
              pft[*,*,v] = float(data)
              
           ENDELSE
        
        ENDIF ELSE BEGIN
           ;; use bilinear interpolation for irregular grid
           pft[*,*,v] = interpolate(float(data),x_pft,y_pft)
        ENDELSE

     ENDFOR
     
     ncdf_close,ncid

  ENDELSE

END
