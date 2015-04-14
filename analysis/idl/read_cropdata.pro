PRO read_cropdata, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                  crop_prefix, filenames, outdata

  ;; directories
  infiles = basedir + '/'+crop_prefix + '/' + filenames + '_5min.nc'

  ;; size of global dataset
  wx = 4320L
  wy = 2160L
  res = 0.083333333D
  
  ;; size and position of tile to be read
  ;; read safe border if possible
  x0 = (round((lonmin + 180.)/res)-1)>0
  y0 = (round((90.-latmax)/res)-1)>0
  x1 = (round((lonmax + 180.)/res))<(wx-1)
  y1 = (round((90.-latmin)/res))<(wy-1)

  ;; size of safe border
  dx0 = round((lonmin + 180.)/res)   - x0
  dx1 = x1 - round((lonmax + 180.)/res-1)
  dy0 = round((90.-latmax)/res)      - y0
  dy1 = y1 -    round((90.-latmin)/res-1)

  ;; total size including safe border
  nx1 = x1 - x0 + 1
  ny1 = y1 - y0 + 1

  ;; actual size of read dataset
  nx = x1 - x0 + 1 - dx0 - dx1
  ny = y1 - y0 + 1 - dy0 - dy1

;  print,x0,x1,y0,y1
;  print,dx0,dx1,dy0,dy1

  ;; arrays
  ndata = n_elements(filenames)
  data = fltarr(nx1,ny1)
  data1 = fltarr(nx1,ny1)

  ;; read areas from netcdf files
  FOR f=0,ndata-1 DO BEGIN
     ncid = ncdf_open(infiles[f])
     ncdf_varget,ncid,'farea',data1,offset=[x0,y0,0,0],count=[nx1,ny1,1,1]
     NCDF_ATTGET,ncid,'farea','missing_value',missing_value
     badvalues = where(data1 EQ missing_value,badcount)
     IF badcount GT 0 THEN data1[badvalues]=0. ; set missing values to 0
     data = data + data1 * 100.
     ncdf_close,ncid
  ENDFOR

  ;; regrid to desired output grid if needed

  ;; number of lon/lat including safe border
  nlon1 = round((lonmax - lonmin + (dx0+dx1)*res )/dlon)
  nlat1 = round((latmax - latmin + (dy0+dy1)*res )/dlat)
  ;; number of lon/lat on output
  nlon = round((lonmax - lonmin)/dlon)
  nlat = round((latmax - latmin)/dlat)

  IF (nlon NE nx) OR (nlat NE ny) THEN BEGIN
      IF (nlon LE nx) AND (nlat LE ny) THEN BEGIN
;;          print,'crop averaging'
          ;; find common denominator for scaling, then do minifying
          IF nlon*nlat GT 1 THEN $
            outdata = reverse(rebin(congrid(data[dx0:nx+dx0-1,dy0:ny+dy0-1],nlon*nx,nlat*ny),nlon, nlat),2) $
          ELSE $
            outdata = rebin(congrid(data[dx0:nx+dx0-1,dy0:ny+dy0-1],nlon*nx,nlat*ny),nlon, nlat)
      ENDIF ELSE BEGIN
;;          print,'crop magnifying'
          dx2 = round(dx0*res/dlon)
          dy2 = round(dy0*res/dlat)
          outdata = (reverse(congrid(data,nlon1, nlat1,/interp,/center),2))[dx2:dx2+nlon-1,dy2:dy2+nlat-1]
;;          outdata = (reverse(congrid(data,nlon1, nlat1),2))[dx2:dx2+nlon-1,dy2:dy2+nlat-1]
      ENDELSE
  ENDIF ELSE BEGIN
     outdata = reverse(data[dx0:nx+dx0-1,dy0:ny+dy0-1],2)
  ENDELSE

  ;print,'crops: ',mean(outdata)

END
