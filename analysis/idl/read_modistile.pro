;; IDL Module to handle MODIS Land Products on sinusoidal grid
;; includes QA screening

;; 2009/06/07 Reto Stockli (Blue Marble Research)

FUNCTION btest,number,j
  ;; tests whether bit j is set in number and returns 0 (not set) or 1 (set)
  powerOfTwo = 2L^j
  IF (LONG(number) AND powerOfTwo) EQ powerOfTwo THEN bit = 1 ELSE bit = 0
  RETURN, bit
END

FUNCTION btestarr,array,j
  ;; tests whether bit j is set in each
  ;; element of array and returns array
  ;; indices for elements where bit is set
  powerOfTwo = 2L^j
  n = n_elements(array)
  bitindex = where((long(array) AND powerOfTwo) EQ powerOfTwo,nset)
  RETURN, bitindex
END

PRO read_file, ftemp, varname,h,v, hmodis,vmodis,xidmodis,yidmodis,xmodis,ymodis, $
               valout, units, long_name, modres
  
  nx = 1200L*modres
  ny = 1200L*modres
  
  npmodis=n_elements(hmodis)

  sd_id = hdf_sd_start(ftemp,/read)
  index = hdf_sd_nametoindex(sd_id,varname)
  sds_id = hdf_sd_select(sd_id,index)

  units=''
  longname=''
  hdf_fill=!values.f_nan
  scale_factor=1
  add_offset=0
  valid_range=[!values.f_nan,!values.f_nan]

  IF HDF_SD_ATTRFIND(sds_id,"units") GE 0 THEN $
     HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"units"),DATA=units
  IF HDF_SD_ATTRFIND(sds_id,"long_name") GE 0 THEN $
     HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"long_name"),DATA=long_name
  IF HDF_SD_ATTRFIND(sds_id,"_FillValue") GE 0 THEN $
     HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"_FillValue"),DATA=hdf_fill

  ;; do not get scale/offset in QC flags
  IF (strpos(varname,'QC') LT 0) AND (strpos(varname,'Quality') LT 0) $
     AND (strpos(varname,'reliability') LT 0) THEN BEGIN
     IF HDF_SD_ATTRFIND(sds_id,"scale_factor") GE 0 THEN  $
        HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"scale_factor"),DATA=scale_factor
     IF HDF_SD_ATTRFIND(sds_id,"add_offset") GE 0 THEN  $
        HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"add_offset"),DATA=add_offset
  ENDIF

  IF HDF_SD_ATTRFIND(sds_id,"valid_range") GE 0 THEN  $
     HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"valid_range"),DATA=valid_range

  add_offset = add_offset[0]
  scale_factor = scale_factor[0]
  hdf_fill = hdf_fill[0]

  ;; ATTENTION: NDVI/EVI have reversed scale factor!!!
  IF (varname EQ '1 km 16 days EVI') OR (varname EQ '1 km 16 days NDVI') THEN scale_factor = 1./scale_factor
  IF (varname EQ '250m 16 days EVI') OR (varname EQ '250m 16 days NDVI') THEN scale_factor = 1./scale_factor
  IF (varname EQ '500m 16 days EVI') OR (varname EQ '500m 16 days NDVI') THEN scale_factor = 1./scale_factor

  ;; ATTENTION: GPP/NPP are in kgC/m2/8days, we want gC/m2/day
  IF (varname EQ 'PsnNet_1km') OR (varname EQ 'Gpp_1km') THEN scale_factor=scale_factor*1000./8.

  index = where ((hmodis eq h) and (vmodis eq v),numpix)

  IF (numpix GT 500L) THEN BEGIN
     ;; read whole tile into memory and then distribute pixels (faster for many pixels)

     hdf_sd_getdata,sds_id,val

     xymodis = ymodis[index]*nx + xmodis[index]

     IF (scale_factor NE 1) OR (add_offset NE 0) THEN BEGIN

        ;; convert to floating point
        val = float(val)
        
        ;; check for out-of-range
        IF finite(valid_range[0]) AND finite(valid_range[1]) THEN BEGIN
           baddata = where ((val LT valid_range[0]) OR (val GT valid_range[1]),badcount)
           IF badcount GT 0 THEN val[baddata]=!values.f_nan
        ENDIF
        
        ;; check for fill values
        IF finite(hdf_fill) THEN BEGIN
           baddata = where(val EQ hdf_fill, badcount)
           IF badcount GT 0 THEN val[baddata]=!values.f_nan
        ENDIF

        ;; re-scale data
        gooddata = where(finite(val),goodcount)
        IF goodcount GT 0L THEN BEGIN
           val[gooddata] = ( val[gooddata] - add_offset)*scale_factor
        ENDIF

        valout[index] = val[xymodis]
     ENDIF ELSE BEGIN
        ;; handle byte/integer/QC values w/o out-of-range checking
        valout[index] = val[xymodis]
     ENDELSE

     val = 0
     xymodis = 0

  ENDIF ELSE BEGIN

     IF numpix GT 0L THEN BEGIN
              
        FOR p = 0,numpix-1 DO BEGIN
           
           hdf_sd_getdata,sds_id,val,start=[xmodis[index[p]],ymodis[index[p]]], count=[1,1]

           IF (scale_factor NE 1) OR (add_offset NE 0) THEN BEGIN

              ;; convert to floating point
              val = float(val[0])
              
              ;; check for out-of-range
              IF finite(valid_range[0]) AND finite(valid_range[1]) THEN BEGIN
                 IF (val LT valid_range[0]) OR (val GT valid_range[1]) THEN val=!values.f_nan
              ENDIF
              
              ;; check for fill values
              IF finite(hdf_fill) THEN BEGIN
                 IF val EQ hdf_fill THEN val=!values.f_nan
              ENDIF
              
              ;; re-scale data
              IF finite(val) THEN BEGIN
                 valout[index[p]] = ( val - add_offset ) * scale_factor
              ENDIF
           ENDIF ELSE BEGIN
              ;; handle byte/integer/QC values w/o out-of-range checking
              valout[index[p]] = val
           ENDELSE
           
        ENDFOR
        
     ENDIF

  ENDELSE
  
  hdf_sd_endaccess,sds_id
  hdf_sd_end,sd_id

END  

PRO read_modis_static,datadir,year,day,modis_prefix,modis_name,hmin,hmax,vmin,vmax, $
                      hmodis,vmodis,xidmodis,yidmodis,xmodis,ymodis, $
                      modis_data,modis_units, modis_longname, npmodis, modres

  syear = string(year,format='(I4)')
  sday = string(day,format='(I3.3)')

  sdate = syear+sday

  FOR h=hmin,hmax DO BEGIN
     FOR v=vmin,vmax DO BEGIN

        sv=string(format='(I2.2)',v)
        sh=string(format='(I2.2)',h)

        infile=(file_search(datadir+'/'+sdate+'/'+modis_prefix+'.A'+sdate+'.h'+sh+'v'+sv+'*.hdf'))[0]

        IF (infile NE '') THEN BEGIN
           read_file, infile,modis_name,h,v,hmodis,vmodis,xidmodis,yidmodis, $
                      xmodis,ymodis, modis_data, modis_units, modis_longname, modres
        ENDIF ELSE BEGIN
           print,'File ',datadir+'/'+sdate+'/'+modis_prefix+'*.h'+sh+'v'+sv+'*.hdf not available !!!'
        ENDELSE
        
     ENDFOR
  ENDFOR

END

PRO ll2tile, rlon,rlat,x,y,h,v, modres, pollon, pollat

  ;; calculates pixel-centered MODIS h/v and x/y from lon/lat
  ;; new: rotated pole supported

  IF (abs(pollat-90.d0) lt 1.d-3) THEN BEGIN
      ;; no rotated pole
      lon = rlon
      lat = rlat
  ENDIF ELSE BEGIN
      ;; rotated pole
      rpol2ll,rlon,rlat,lon,lat,pollon=pollon,pollat=pollat
  ENDELSE

  nx = 1200L*modres  ; indices go from 0 .. nx - 1
  ny = 1200L*modres
  nh = 36L
  nv = 18L
    
  npt = n_elements(lon)

  x = lonarr(npt)
  h = lonarr(npt)

  h[*] = -1

  v = (long(( 90.d0 - double(lat) ) * double(nv) / 180.d0)<(nv-1L))>0

  y = ((round(( ( 90.d0 - double(lat) ) * double(nv) / 180.d0 - double(v) ) $
              * double(ny) - 0.5d0))<(ny-1L))>0L

  FOR htmp=0L,nh-1L DO BEGIN
     xtmp = round( ( (lon * (sin( (double(y) + 0.5d0 + double(v*ny))*!dpi/double(ny*nv) )) $
                  + 180.d0) / 360.d0 * double(nh) - double(htmp) ) * double(nx) - 0.5d0 )
     goodx = where((xtmp GE 0L) AND (xtmp LT nx), goodcount)

     IF goodcount GT 0 THEN BEGIN
        x[goodx]=xtmp[goodx]
        h[goodx]=htmp
     ENDIF
  ENDFOR

  a=where(h LT 0,acount)
  IF acount GT 0 THEN BEGIN
     print,'WRONG H INDEX IN LL2TILE: ', h
     stop
  ENDIF

END

PRO tile2ll, x,y,h,v,rlon,rlat, modres, pollon, pollat

  ;; calculates pixel-centered lon/lat from MODIS x/y + h/v
  ;; new: rotated pole supported

  nx = 1200L*modres    ; indices go from 0 .. nx - 1
  ny = 1200L*modres
  nh = 36L
  nv = 18L

  lat= 90.d0 - ( double(v) + (double(y) + 0.5d0) / double(ny) ) * 180.d0 / double(nv)
  lon= ( -180.d0 + ( double(h) + (double(x) + 0.5d0) / double(nx) ) * 360.d0 / double(nh) ) / $
       (sin( (double(y) + 0.5d0 + double(v*ny))*!dpi/double(ny*nv) )>1.d-8)

  IF (abs(pollat-90.d0) lt 1.d-3) THEN BEGIN
      ;; no rotated pole
      rlon = lon
      rlat = lat
  ENDIF ELSE BEGIN
      ;; rotated pole
      ll2rpol,lon,lat,rlon,rlat,pollon=pollon,pollat=pollat
  ENDELSE

END 

PRO read_modistile, basedir, year, day, lonmin, lonmax, latmin, latmax, dlon, dlat,  $
                    obsvarname, obsunits, obstitle, version, obsvar, obsvarerror, obsmask, $
                    pollon=pollon, pollat=pollat, initialize=initialize

  COMMON MODISTILE,vmodis,hmodis,xmodis,ymodis,xidmodis,yidmodis,modis_lon,modis_lat, $
     npmodis,hmin,hmax,vmin,vmax

;; reads specified MODIS product from sinusoidal tiles for
;; specified area, time and projection

;; check whether we have a rotated pole (needs special treatment for
;; irregular MODIS grid)
  rotatedpole = 0
  IF keyword_set(pollon) AND keyword_set(pollat) THEN BEGIN
     IF (abs(pollat-90.d0) GE 1.d-3) THEN rotatedpole = 1
  ENDIF

;; define MODIS product metadata and QA weighting
;; qc bit arrays go from left (MSB) to right (LSB)
  CASE obsvarname OF
     'LSTDAY' : BEGIN
        modis_prefix='MOD11A1'
        modis_name='LST_Day_1km'
        qc1_name='QC_Day'
        qc2_name=''
        obsunits='K'
        obstitle='Daytime Land Surface Temperature'
        minerror = 0.5 ;; [K]
        qc1flags = [4,2,0,0,0,0,9,0] 
        modres = 1L  
     END
     'LSTNIGHT' : BEGIN
        modis_prefix='MOD11A1'
        modis_name='LST_Night_1km'
        qc1_name='QC_Night'
        qc2_name=''
        obsunits='K'
        obstitle='Nighttime Land Surface Temperature'
        minerror = 0.5 ;; [K]
        qc1flags = [4,2,0,0,0,0,9,0] 
        modres = 1L  
     END

     'LAI': BEGIN
        modis_prefix='MOD15A2'
        modis_name='Lai_1km'
        qc1_name='FparLai_QC'
        qc2_name='FparExtra_QC'
        obsunits='m2/m2'
        obstitle='Leaf Area Index'
        minerror = 0.1
        qc1flags = [9,4,2,9,9,9,0,1] ;fparlaiqc (MOD15 QA bits 7 - 0)
        qc2flags = [0,9,9,9,5,9,9,9] ;fparextraqc (MOD15 QA bits 7 - 0)
        modres = 1L  
     END
     'FPAR': BEGIN
        modis_prefix='MOD15A2'
        modis_name='Fpar_1km'
        qc1_name='FparLai_QC'
        qc2_name='FparExtra_QC'
        obsunits='-'
        obstitle='Fraction of PAR absorbed by plants'
        minerror = 0.02
        qc1flags = [9,4,2,9,9,9,0,1] ;fparlaiqc (MOD15 QA bits 7 - 0)
        qc2flags = [0,9,9,9,5,9,9,9] ;fparextraqc (MOD15 QA bits 7 - 0)
        modres = 1L 
     END
     'NDVI': BEGIN
        modis_prefix='MOD13A2'
        modis_name='1 km 16 days NDVI'
        qc1_name='1 km 16 days VI Quality'
        qc2_name='1 km 16 days pixel reliability'
        modres = 1L 
        obsunits='-'
        obstitle='Normalized Difference Vegetation Index'
        minerror = 0.01
        qc1flags = [9,9,9,9,0,7,7,6,4,2,9,8,4,1,9,1] 
        qc2flags = [9,9,9,9,9,9,9,0]
     END
     'EVI': BEGIN
        modis_prefix='MOD13A2'  ;; 1km
        modis_name='1 km 16 days EVI'
        qc1_name='1 km 16 days VI Quality'
        qc2_name='1 km 16 days pixel reliability'
        modres = 1L
        obsunits='-'
        obstitle='Enhanced Vegetation Index'
        minerror = 0.01
        qc1flags = [9,9,9,9,0,7,7,6,4,2,9,8,4,1,9,1] 
        qc2flags = [9,9,9,9,9,9,9,0]
     END
     'GPP': BEGIN
        modis_prefix='MOD17A2'
        modis_name='Gpp_1km'
        qc1_name='Psn_QC_1km'
        qc2_name=''
        obsunits='gC/m2/day'
        obstitle='Gross Primary Production'
        minerror = 0.05
        modres = 1L 
        qc1flags = [9,4,1,2,1,2,0,1] 
     END
     'NPP': BEGIN
        modis_prefix='MOD17A2'
        modis_name='PsnNet_1km'
        qc1_name='Psn_QC_1km'
        qc2_name=''
        obsunits='gC/m2/day'
        obstitle='Net Primary Production'
        minerror = 0.05
        modres = 1L   
        qc1flags = [9,4,1,2,1,2,0,1] 
     END
     'WSA': BEGIN
        modis_prefix='MCD43B3'
        modis_name='Albedo_WSA_shortwave'
        qc1_name=''
        qc2_name=''
        obsunits='-'
        obstitle='White Sky Albedo'
        minerror = 0.01
        modres = 1L   
     END
     'BSA': BEGIN
        modis_prefix='MCD43B3'
        modis_name='Albedo_BSA_shortwave'
        qc1_name=''
        qc2_name=''
        obsunits='-'
        obstitle='Black Sky Albedo'
        minerror = 0.01
        modres = 1L   
     END
     'SNOW': BEGIN
        modis_prefix='MOD10A1'
        modis_name='Fractional_Snow_Cover'
        obstitle='Fractional Snow Cover'
        qc1_name=''
        qc2_name=''
        obsunits='%'
        minerror = 1.0
        modres = 2L   
     END
     'FIRE': BEGIN
        modis_prefix='MOD14A2'
        modis_name='FireMask'
        qc1_name=''
        qc2_name=''
        obsunits='Number of 1 km^2 Pixels'
        obstitle='Fire Count'
        minerror = 0.01
        modres = 1L   
     END
     'IGBP': BEGIN
        IF version EQ 4 THEN modis_prefix='MOD12Q1'
        IF version EQ 5 THEN modis_prefix='MCD12Q1'
        modis_name='Land_Cover_Type_1'
        qc1_name=''
        qc2_name=''
        obsunits=''
        obstitle='IGBP Land Cover Classification'
        minerror = 0.0
        IF version EQ 4 THEN modres = 1L   
        IF version EQ 5 THEN modres = 2L   
     END
     'UMD': BEGIN
        modis_prefix='MOD12Q1'
        modis_name='Land_Cover_Type_2'
        qc1_name=''
        qc2_name=''
        obsunits=''
        obstitle='UMD Land Cover Classification'
        minerror = 0.0
        modres = 1L   
     END
     'WATER': BEGIN
        modis_prefix='MOD44W'
        modis_name='water_mask'
        qc1_name=''
        qc2_name=''
        obsunits=''
        obstitle='Water Mask'
        minerror = 0.0
        modres = 4L
     END
     'TREE': BEGIN
        modis_prefix='MOD44B'
        IF (version EQ 3) THEN modis_name='Percent Tree Cover'
        IF (version EQ 5) THEN modis_name='Percent_Tree_Cover'
        qc1_name=''
        qc2_name=''
        obsunits=''
        obstitle='Percent Tree Cover'
        minerror = 0.0
        IF (version EQ 3) THEN modres = 2L   
        IF (version EQ 5) THEN modres = 4L
     END
     'HERB': BEGIN
        modis_prefix='MOD44B'
        modis_name='Percent Nontree Vegetation'
        qc1_name=''
        qc2_name=''
        obsunits=''
        obstitle='Percent Herbaceous Cover'
        minerror = 0.0
        modres = 2L   
     END
     'BARE': BEGIN
        modis_prefix='MOD44B'
        modis_name='Percent Bare'
        qc1_name=''
        qc2_name=''
        obsunits=''
        obstitle='Percent Bare Soil'
        minerror = 0.0
        modres = 2L   
     END
     ELSE: BEGIN
        print,"MODIS product not defined in read_modistile.pro"
        stop
     END
  ENDCASE


  ;; initializing geographic boundary conditions
  IF ((day EQ 0) OR (year EQ 0)) OR keyword_set(initialize) THEN BEGIN

     ;; define modis spatial array dimensions
     nx = 1200L*modres
     ny = 1200L*modres
     res = 360.d/36.d/double(nx)

     ;; define number of elements in lon/lat
     nlon = round((lonmax-lonmin)/dlon)
     nlat = round((latmax-latmin)/dlat)
     
     ;; create lon/lat arrays (centered on grid point)
     lon = (dindgen(nlon)*dlon + lonmin + 0.5d*dlon) # replicate(1.d,nlat)    ; W->E
     lat = replicate(1.d,nlon) # (dindgen(nlat)*dlat + latmin + 0.5d*dlat)    ; S->N
     
     ;; preliminary check for number of expected observations in area
     npmodis = 0L

     IF rotatedpole THEN BEGIN

        FOR j=0L,nlat-1L DO BEGIN
           FOR k=0L,((round(dlat/res)-1L)>0L) DO BEGIN
              FOR i=0L,nlon-1L DO BEGIN
                 FOR l=0L,((round(dlon/res)-1L)>0L) DO BEGIN
                    
                    ;; we feed pixel-centered lon/lat to ll2tile
                    ll2tile,lon[i,j]-0.5d*dlon+(double(l)+0.5d)*res,$
                            lat[i,j]-0.5d*dlat+(double(k)+0.5d)*res, $
                            x0,y0,h0,v0, modres,pollon,pollat
                    
                    npmodis = npmodis + 1L
                    
                 ENDFOR
              ENDFOR
           ENDFOR
        ENDFOR

     ENDIF ELSE BEGIN

        FOR j=0L,nlat-1L DO BEGIN

           FOR k=0L,((round(dlat/res)-1L)>0L) DO BEGIN
              ;; we feed pixel-centered lon/lat to ll2tile
              ll2tile,lon[*,j]-0.5d*dlon+0.5d*res,$
                      lat[*,j]-0.5d*dlat+(double(k)+0.5d)*res,$
                      x0,y0,h0,v0, modres,pollon,pollat
              
              ll2tile,lon[*,j]+0.5d*dlon-0.5d*res,$
                      lat[*,j]-0.5d*dlat+(double(k)+0.5d)*res,$
                      x1,y1,h1,v1, modres,pollon,pollat

              FOR i=0L,nlon-1L DO BEGIN

                 xmin0 = x0[i]
                 hmin0 = h0[i]
                 
                 xmax0 = x1[i]
                 hmax0 = h1[i]
                 
                 IF (xmax0 LT xmin0) AND (hmax0 EQ hmin0) THEN xmax0=xmin0
                 
                 IF (xmax0 LT 0) AND (hmax0 GT hmin0) THEN BEGIN
                    xmax0 = nx - 1L
                    hmax0 = hmax0 - 1L
                 ENDIF
                 
                 npmodis = npmodis + (hmax0-hmin0)*nx + (xmax0-xmin0+1L)
                 
              ENDFOR
           ENDFOR
        ENDFOR

     ENDELSE

     print,'We expect ',npmodis,' MODIS pixels in the selected area (maximum).'

     ;; allocate static arrays based on number of expected observations in area
     vmodis=intarr(npmodis)
     hmodis=intarr(npmodis)
     xmodis=intarr(npmodis)
     ymodis=intarr(npmodis)
     xidmodis=intarr(npmodis)
     yidmodis=intarr(npmodis)

     modis_lon=dblarr(npmodis)
     modis_lat=dblarr(npmodis)

;; now derive modis indices and store them into above arrays
     npmodis = 0L

     IF rotatedpole THEN BEGIN

        FOR j=0L,nlat-1L DO BEGIN
           FOR k=0L,((round(dlat/res)-1L)>0L) DO BEGIN
              FOR i=0L,nlon-1L DO BEGIN
                 FOR l=0L,((round(dlon/res)-1L)>0L) DO BEGIN
                    
                    ;; we feed pixel-centered lon/lat to ll2tile
                    ll2tile,lon[i,j]-0.5d*dlon+(double(l)+0.5d)*res,$
                            lat[i,j]-0.5d*dlat+(double(k)+0.5d)*res, $
                            x0,y0,h0,v0, modres,pollon,pollat
                    
                    n = 1
                    
                    vmodis[npmodis:npmodis+n-1L] = v0[0]
                    ymodis[npmodis:npmodis+n-1L] = y0[0]

                    xidmodis[npmodis:npmodis+n-1L] = i
                    yidmodis[npmodis:npmodis+n-1L] = j
                    
                    hmodis[npmodis:npmodis+n-1L] = h0[0]
                    xmodis[npmodis:npmodis+n-1L] = x0[0]
                    
                    npmodis = npmodis + n

                 ENDFOR
              ENDFOR
           ENDFOR
        ENDFOR

     ENDIF ELSE BEGIN

        FOR j=0L,nlat-1L DO BEGIN
           
           FOR k=0L,((round(dlat/res)-1L)>0L) DO BEGIN
              ;; we feed pixel-centered lon/lat to ll2tile
              ll2tile,lon[*,j]-0.5d*dlon+0.5d*res,$
                      lat[*,j]-0.5d*dlat+(double(k)+0.5d)*res,$
                      x0,y0,h0,v0, modres,pollon,pollat
              ll2tile,lon[*,j]+0.5d*dlon-0.5d*res,$
                      lat[*,j]-0.5d*dlat+(double(k)+0.5d)*res,$
                      x1,y1,h1,v1, modres,pollon,pollat
              
              FOR i=0L,nlon-1L DO BEGIN
                 xmin0 = x0[i]
                 hmin0 = h0[i]
                 
                 xmax0 = x1[i]
                 hmax0 = h1[i]
                 
                 IF (xmax0 LT xmin0) AND (hmax0 EQ hmin0) THEN xmax0=xmin0
                 
                 IF (xmax0 LT 0) AND (hmax0 GT hmin0) THEN BEGIN
                    xmax0 = nx - 1L
                    hmax0 = hmax0 - 1L
                 ENDIF
                 
                 n = (hmax0-hmin0)*nx + (xmax0-xmin0+1L)
                 
                 vmodis[npmodis:npmodis+n-1L] = v0[0]
                 ymodis[npmodis:npmodis+n-1L] = y0[0]
                 xidmodis[npmodis:npmodis+n-1L] = i
                 yidmodis[npmodis:npmodis+n-1L] = j
                 
                 xi = indgen(n)
                 
                 hmodis[npmodis:npmodis+n-1L] = hmin0 + (xmin0 + xi)/nx
                 xmodis[npmodis:npmodis+n-1L] = (xmin0 + xi) MOD nx
                 
                 npmodis = npmodis + n
                 
              ENDFOR

           ENDFOR
        ENDFOR

     ENDELSE    

     ;; get MODIS lon/lat to each x/y/h/v value
     tile2ll,xmodis,ymodis, hmodis,vmodis, $
             modis_lon, modis_lat, modres, pollon, pollat        
     
     ;; limit valid MODIS lat/lon to regional bounds
     modis_lon = (modis_lon>lonmin)<lonmax
     modis_lat = (modis_lat>latmin)<latmax

     hmax = max(hmodis)
     hmin = min(hmodis)
     vmax = max(vmodis)
     vmin = min(vmodis)
     
     print,'Check for requested tile bounds: '
     print,'MODIS lon/lat bounds (lonmin, lonmax, latmin, latmax): ',lonmin,lonmax,latmin,latmax
     print,'MODIS tile bounds (hmin, hmax, vmin, vmax): ',hmin,hmax,vmin,vmax

;     print,xmodis
;     print,ymodis

  ENDIF

  ;; reading data
  IF (day NE 0) AND (year NE 0) THEN BEGIN

     ;; define MODIS collection (version)
     sversion = string(version,format='(I3.3)')
     
     ;; define directories
     datadir=basedir+'/modis/'+modis_prefix+'.'+sversion
     datafile=modis_prefix

     ;; set up variables
     obsvar=make_array(npmodis,/float,value=!values.f_nan)
     obsvarerror=make_array(npmodis,/float,value=!values.f_nan)
     obsmask=make_array(npmodis,/byte,value=0B)

     syear = string(year,format='(I4)')
     sday = string(day,format='(I3.3)')
     sdate = syear+sday

;; read the MODIS data for given date, and tile by tile
;; fill observation and quality flag vectors with up to npmodis
;; observations
;; then quality screen the observation vector based on above quality
;; flags filters and return observation + error vector together with
;; spatial indices

     IF qc1_name NE '' THEN $
        IF (n_elements(qc1flags) EQ 8) THEN qc1 = make_array(npmodis,/byte,value=0B) $
        ELSE qc1 = make_array(npmodis,/uint,value=0)
     IF qc2_name NE '' THEN $
        IF (n_elements(qc2flags) EQ 8) THEN qc2 = make_array(npmodis,/byte,value=0B) $
        ELSE qc2 = make_array(npmodis,/uint,value=0)

     ;; cycle through L4 tiles and read data+QA
     FOR h=hmin,hmax DO BEGIN
        FOR v=vmin,vmax DO BEGIN

           sv=string(format='(I2.2)',v)
           sh=string(format='(I2.2)',h)

           infile=(file_search(datadir+'/'+sdate+'/'+datafile+'.A'+sdate+'.h'+sh+'v'+sv+'*.hdf'))[0]

           IF (infile NE '') THEN BEGIN
              read_file, infile,modis_name,h,v, hmodis,vmodis,xidmodis,yidmodis, $
                         xmodis,ymodis, obsvar, units, long_name, modres
              IF qc1_name NE '' THEN read_file, infile,qc1_name,h,v, hmodis,vmodis,xidmodis,yidmodis, $
                                                xmodis,ymodis, qc1, units, long_name, modres
              IF qc2_name NE '' THEN read_file, infile,qc2_name,h,v, hmodis,vmodis,xidmodis,yidmodis, $
                                                xmodis,ymodis, qc2, units, long_name, modres
           ENDIF ELSE BEGIN
;;              print,"MODIS Tile not found: ",datafile+'.A'+sdate+'.h'+sh+'v'+sv+'*.hdf'
           ENDELSE
           
        ENDFOR
     ENDFOR

     ;; evaluate data quality and screen w. QA information from
     ;; MODIS Products
     goodobs = where(finite(obsvar),goodcount)
     
     IF (goodcount GT 0L) THEN BEGIN

        CASE modis_name OF

           'Fractional_Snow_Cover' : BEGIN
              ;; daily snow: no QA, discrete distribution: 
              ;; value 1-100 : snow cover %  ok
              ;; value 200 : missing data
              ;; value 201 : no decision
              ;; value 211 : night
              ;; value 225 : land   ok
              ;; value 237 : inland water ok
              ;; value 239 : ocean
              ;; value 250 : cloud
              ;; value 254 : detector saturated
              ;; value 255 : missing
              
              bval = byte(round(obsvar[goodobs]))
              
              goodobs2 = where((bval LE 100B) OR (bval EQ 225B) OR $
                               (bval EQ 237B),goodcount2)     

              IF (goodcount2 GT 0L) THEN BEGIN
                 
                 obsvar[goodobs[goodobs2]] = (bval[goodobs2] LE 100B)*float(bval[goodobs2]) + $
                                             (bval[goodobs2] GT 100B)*0.
                 obsvarerror[goodobs[goodobs2]] = minerror
                 obsmask[goodobs[goodobs2]] = 1B
                 
              ENDIF
              
           END

           'Maximum_Snow_Extent' : BEGIN
              ;; snow 8 day composite: no QA, discrete distribution: 
              ;; value 25  : land,  no snow
              ;; value 37  : lake,  no snow
              ;; value 39  : ocean, no snow
              ;; value 100 : lake,  ice
              ;; value 200 : land,  snow cover
              
              bval = byte(round(obsvar[goodobs]))
              
              goodobs2 = where((bval EQ 25B) OR (bval EQ 37B) OR $
                               (bval EQ 39B) OR (bval EQ 100B) OR $
                               (bval EQ 200B),goodcount2)     

              IF (goodcount2 GT 0L) THEN BEGIN
                 
                 obsvar[goodobs[goodobs2]] = (bval[goodobs2] EQ 100B)*100. + $
                                             (bval[goodobs2] EQ 200B)*100.
                 obsvarerror[goodobs[goodobs2]] = minerror
                 obsmask[goodobs[goodobs2]] = 1B
                 
              ENDIF
              
           END

           'FireMask' : BEGIN
              ;; FIRE: QA part of fire mask
              
              bval = byte(round(obsvar[goodobs]))
              
              goodobs2 = where((bval EQ 7B) OR (bval EQ 8B) OR (bval EQ 9B),goodcount2)     
              IF (goodcount2 GT 0L) THEN BEGIN
                 
                 obsvar[goodobs[goodobs2]] = ((bval[goodobs2] GE 7B) AND (bval[goodobs2] LE 9B))
                 obsvarerror[goodobs[goodobs2]] = minerror + minerror * $
                                                  ((bval[goodobs2] EQ 8B)*4. + (bval[goodobs2] EQ 7B)*9.) 
                 obsmask[goodobs[goodobs2]] = 1B
                 
              ENDIF

           END

           'Land_Cover_Type_1' : BEGIN
              ;; Land Cover, no QA Screening
              
              obsvarerror[goodobs] = minerror
              obsmask[goodobs] = 1B
              
           END

           'Land_Cover_Type_2' : BEGIN
              ;; Land Cover, no QA Screening
              
              obsvarerror[goodobs] = minerror
              obsmask[goodobs] = 1B
              
           END
           
           'Percent Tree Cover' : BEGIN
              ;; MODIS VCF, no QA Screening
              
              obsvarerror[goodobs] = minerror
              obsmask[goodobs] = 1B
              
           END
           
           'Percent Nontree Vegetation' : BEGIN
              ;; MODIS VCF, no QA Screening
              
              obsvarerror[goodobs] = minerror
              obsmask[goodobs] = 1B
              
           END

           'Percent Bare' : BEGIN
              ;; MODIS VCF, no QA Screening
              
              obsvarerror[goodobs] = minerror
              obsmask[goodobs] = 1B
              
           END
           
           ELSE : BEGIN
              ;; other variables: regular QA filtering

              obsvarerror[*] = minerror
              obsmask[goodobs] = 1B
              
              ;; bitwise screening QA flags (btests starts from LSB -> MSB with argument 0 .. nbits-1
              ;; above in the bitmask, MSB is left and LSB is right
              
              IF qc1_name NE '' THEN BEGIN
                 nbits=n_elements(qc1flags)
                 FOR b=0,nbits-1 DO BEGIN 
                    ;; we begin left (MSB) and move to the right (LSB) for each bitarray
                    idx = btestarr(qc1,nbits-b-1)
                    IF idx[0] GE 0 THEN BEGIN
                       IF (qc1flags[b] eq 9) THEN obsmask[idx]=0B
                       obsvarerror[idx]=obsvarerror[idx]+float(qc1flags[b])*minerror
                    ENDIF
                 ENDFOR
                 idx = where(qc1 EQ 65B,count)
                 IF count GT 0L THEN obsmask[idx] = 0B
              ENDIF
              
              IF qc2_name NE '' THEN BEGIN
                 nbits=n_elements(qc2flags)
                 FOR b=0,nbits-1 DO BEGIN 
                    ;; we begin left (MSB) and move to the right (LSB) for each bitarray
                    idx = btestarr(qc2,nbits-b-1)
                    IF idx[0] GE 0 THEN BEGIN
                       IF (qc2flags[b] eq 9) THEN obsmask[idx]=0B
                       obsvarerror[idx]=obsvarerror[idx]+float(qc2flags[b])*minerror
                    ENDIF
                 ENDFOR
              ENDIF

              badidx = where(obsmask EQ 0B,badcount)
              IF badcount GT 0L THEN BEGIN
                 obsvar[badidx] = !values.f_nan
                 obsvarerror[badidx] = !values.f_nan
              ENDIF

;              print,day
;              FOR p=0,npmodis-1 DO BEGIN
;                 print,xmodis[p],ymodis[p],obsvar[p],obsvarerror[p]
;              ENDFOR

;              IF day EQ 1 THEN stop

           END

        ENDCASE ;; modis name

     ENDIF ;; any valid observations

  ENDIF ;; if reading

END


