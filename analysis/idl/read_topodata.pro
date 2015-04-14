PRO read_topodata, basedir, topodir, toporegion, lon, lat, topo, dlon=dlon, dlat=dlat
;; read topography file (is a composite of GTOPO30 > 60N, CGIAR SRTM 60N-60S)
;; Reto Stockli (Blue Marble Research), February 2009

  ;; check whether we have a regular grid or not
  IF keyword_set(dlon) AND keyword_set(dlat) THEN regular=1 ELSE regular=0

;; define input file
  infile = (basedir+'/'+topodir+'/topo.'+toporegion+'.nc')[0]
  temp = file_search(infile)
  IF temp EQ '' THEN BEGIN
     print,'TOPO file: ',infile,' not found. stopping.'
     stop
  ENDIF

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
  ncid = ncdf_open(infile)

  varid = ncdf_varid(ncid,'lon')
  IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlon')
  ncdf_varget,ncid,varid,lons_topo

  varid = ncdf_varid(ncid,'lat')
  IF varid EQ -1 THEN varid = ncdf_varid(ncid,'rlat')
  ncdf_varget,ncid,varid,lats_topo

  IF n_elements(lons_topo) EQ 1L THEN dlon_topo = dlon ELSE $
     dlon_topo = abs(lons_topo[1]-lons_topo[0])
  IF n_elements(lats_topo) EQ 1L THEN dlat_topo = dlat ELSE $
     dlat_topo = abs(lats_topo[1]-lats_topo[0])

  latmin_topo = min(lats_topo) - 0.5d0*dlat_topo
  lonmin_topo = min(lons_topo) - 0.5d0*dlon_topo
  latmax_topo = max(lats_topo) + 0.5d0*dlat_topo
  lonmax_topo = max(lons_topo) + 0.5d0*dlon_topo

  nlon_topo = n_elements(lons_topo)
  nlat_topo = n_elements(lats_topo)

  xmin_topo = -1L
  xmax_topo = -1L
  ymin_topo = -1L
  ymax_topo = -1L

;; define precision for getting rid of rounding errors below
  pr = 1.d-6

  FOR i=nlon_topo-1L,0L,-1L DO IF (round((lons_topo[i]-0.5d0*dlon_topo)/pr,/L64) LE round(lonmin/pr,/L64)) $
     AND (xmin_topo EQ -1L) THEN xmin_topo = i
  FOR i=nlat_topo-1L,0L,-1L DO IF (round((lats_topo[i]-0.5d0*dlat_topo)/pr,/L64) LE round(latmin/pr,/L64)) $
     AND (ymin_topo EQ -1L) THEN ymin_topo = i
  FOR i=0L,nlon_topo-1L,1L DO IF (round((lons_topo[i]+0.5d0*dlon_topo)/pr,/L64) GE round(lonmax/pr,/L64)) $
     AND (xmax_topo EQ -1L) THEN xmax_topo = i
  FOR i=0L,nlat_topo-1L,1L DO IF (round((lats_topo[i]+0.5d0*dlat_topo)/pr,/L64) GE round(latmax/pr,/L64)) $
     AND (ymax_topo EQ -1L) THEN ymax_topo = i
  
  IF ((xmin_topo EQ (-1L)) OR (ymin_topo EQ (-1L)) OR (xmax_topo EQ (-1L)) OR (ymax_topo EQ (-1L))) THEN BEGIN
     print,'topography dataset does not cover requested area. stopping.'
     print,'We need: ',latmin,latmax,lonmin,lonmax
     print,'We get:  ',latmin_topo,latmax_topo,lonmin_topo,lonmax_topo
     print,xmin_topo,xmax_topo,ymin_topo,ymax_topo
     stop
  ENDIF

  nx_topo = xmax_topo - xmin_topo + 1L
  ny_topo = ymax_topo - ymin_topo + 1L

  IF (NOT regular) THEN BEGIN
     ;; adjust for boundary conditions in bilinear interpolation
     lontmp = (lon>(lonmin_topo+0.5d0*dlon_topo))<(lonmax_topo-0.5d0*dlon_topo)
     lattmp = (lat>(latmin_topo+0.5d0*dlat_topo))<(latmax_topo-0.5d0*dlat_topo)

     ;; derive bilinear interpolation indices
     x_topo = (lontmp - lons_topo[xmin_topo])/dlon_topo
     y_topo = (lattmp - lats_topo[ymin_topo])/dlat_topo
  ENDIF

;; read topography
  ncdf_varget,ncid, 'Z', data, offset = [xmin_topo,ymin_topo], count = [nx_topo,ny_topo]

  ncdf_close,ncid

  data = float(data)

  IF regular THEN BEGIN

;; scale topography to output dimension if needed
     IF (nx_topo ne nlon) OR (ny_topo ne nlat) THEN BEGIN

        IF (nlon GT nx_topo) OR (nlat GT ny_topo) THEN BEGIN
           ;; magnifying: bilinear interpolation
           topo = congrid(data,nlon,nlat,/interp,/center)

        ENDIF ELSE BEGIN
           ;; minifying: two methods
           IF ((nx_topo EQ (nx_topo/nlon)*nlon) AND (ny_topo EQ (ny_topo/nlat)*nlat)) OR $
              ((float(nlon)/float(nx_topo) EQ float(nlon/nx_topo)) AND $
               (float(nlat)/float(ny_topo) EQ float(nlat/ny_topo))) THEN BEGIN
              ;; simple rebinning -- fast
              topo = rebin(data,nlon,nlat) 
           ENDIF ELSE BEGIN
              ;; new dimensions are not integer factor of old dimensions
              topo = fix(rebin_new(data,nlon,nlat))
           ENDELSE

        ENDELSE

     ENDIF ELSE BEGIN
        ;; no regridding
        topo = data
     ENDELSE
     
  ENDIF ELSE BEGIN
     ;; use bilinear interpolation for irregular grid
     topo = interpolate(data,x_topo,y_topo)
  ENDELSE


END
