PRO read_avhrrvcf, basedir, lonmin, latmin, lonmax, latmax, dlon, dlat, $
                  avhrr_prefix, filename, outdata

  ;; directories
  infile = basedir + '/'+avhrr_prefix + '/' + filename

  ;; size of global dataset
  wx = 43200L
  wy = 21600L
  res = 0.0083333333D
  
  ;; size and position of tile to be read
  nx = round((lonmax - lonmin)/res)
  ny = round((latmax - latmin)/res)
  x0 = round((lonmin + 180.)/res)
  y0 = round((90.-latmax)/res)

  ;; arrays
  data = bytarr(nx,ny)
  line = bytarr(nx)

  ;; read area line by line
  get_lun,lun
  openr,lun,infile
  FOR y = 0L,ny-1L DO BEGIN
     point_lun,lun,(y+y0)*wx + x0
     readu,lun,line
     data[*,y] = line
  ENDFOR
  close,lun
  free_lun,lun

;;For treecover:
;;10 - 80 percent tree cover
;;254     non-vegetated
;;255     tree cover less than 10%
;;For evergreen, deciduous, broadleaf, and needleleaf:
;;10 - 80 percent cover for indicated leaf longevity and type
;;        (% evergreen + % deciduous = % tree cover; and
;;         % broadleaf + % needleleaf = % tree cover)

  badvalues = where (data EQ 255B,badcount)
  IF badcount GT 0 THEN data[badvalues]=0B
  badvalues = where (data EQ 254B,badcount)
  IF badcount GT 0 THEN data[badvalues]=0B

  ;; regrid to desired output grid if needed
  nlon = round((lonmax - lonmin)/dlon)
  nlat = round((latmax - latmin)/dlat)

  data = float(reverse(data,2))

  ;; rescale DATA if needed
  IF (nlon NE nx) OR (nlat NE ny) THEN BEGIN

     IF (nlon GT nx) OR (nlat GT ny) THEN BEGIN
        ;; magnifying: bilinear interpolation
        outdata = congrid(data,nlon,nlat,/interp,/center)

     ENDIF ELSE BEGIN
        ;; minifying

        IF ((nx EQ (nx/nlon)*nlon) AND (ny EQ (ny/nlat)*nlat)) OR $
           ((float(nlon)/float(nx) EQ float(nlon/nx)) AND $
            (float(nlat)/float(ny) EQ float(nlat/ny))) THEN BEGIN
           ;; simple rebinning -- fast
           outdata = rebin(data,nlon,nlat) 
        ENDIF ELSE BEGIN
           ;; use new rebinning for non-integer dimensional scaling
           outdata = fix(rebin_new(data,nlon,nlat))
        ENDELSE

     ENDELSE

  ENDIF ELSE BEGIN
     ;; no regridding
     outdata = data

  ENDELSE

END
