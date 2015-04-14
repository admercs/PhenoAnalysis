PRO regrid_modis,datadir, year, lonmin, lonmax, latmin, latmax, dlon, dlat, pollon, pollat, $
                 obsvar, obsvarerror, obsvarname, obsvarunits, obsvarlongname, modis_version

  COMMON MODISTILE,vmodis,hmodis,xmodis,ymodis,xidmodis,yidmodis,modis_lon,modis_lat, $
     npmodis,hmin,hmax,vmin,vmax


  ;; set up diagnostic timers
  time_init = 0.
  time_read = 0.
  time_grid = 0.

  ;; flags
  errscale = 1.0 ;; scale qualitative error bars
  usewgt = 0     ;; use weighted regridding 

  ;; diagnose leap year
  leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
  if leapyear THEN ndayofyear = 366 ELSE ndayofyear = 365

  ;; initalize MODIS grid
  read_modistile, datadir, 0, 0, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                  obsvarname, obsvarunits, obsvarlongname, modis_version, $
                  obsvartemp, obsvarerrortemp, obsmasktemp, $
                  pollon=pollon, pollat=pollat,/initialize
  
  ;; rebin observation grid to plotting grid
  nlon = round((lonmax-lonmin)/dlon)
  nlat = round((latmax-latmin)/dlat)
  xyindex = long(yidmodis)*nlon + long(xidmodis)
  hist = histogram(xyindex,omin=omin,reverse_indices = R)

  npobs=npmodis
  ntobs=ndayofyear           
  dtobs=24.*3600.               ; observations time step [s]

  ;; initialize yearly data arrays
  obsvar=make_array(nlon,nlat,ntobs,/float,value=!values.f_nan)
  obsvarerror=make_array(2,nlon,nlat,ntobs,/float,value=!values.f_nan)

  FOR t=0,ntobs-1 DO BEGIN

     ;; read MODIS data day by day (slow, for large grids)
     read_modistile, datadir, year, t+1, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                     obsvarname, obsvarunits, obsvarlongname, modis_version, $
                     obsvartemp, obsvarerrortemp, obsmasktemp, $
                     pollon=pollon, pollat=pollat

     IF total(finite(obsvartemp)) GT 0L THEN BEGIN

        tmp = make_array(nlon,nlat,/float,value=!values.f_nan)
        tmperr = make_array(nlon,nlat,/float,value=!values.f_nan)

        FOR i=0L,n_elements(hist)-1L DO BEGIN
           IF R[i] NE R[i+1] THEN BEGIN

              idx = where(obsmasktemp[R[R[i]:R[i+1L]-1L]],validcount)
              
              IF validcount GE 1 THEN BEGIN

                 IF usewgt THEN BEGIN
                    ;; variables are averaged weighted by observation error
                    tmp[i+omin] = total(obsvartemp[R[R[i]:R[i+1]-1]] $
                                        * float(obsmasktemp[R[R[i]:R[i+1]-1]]) / $
                                        obsvarerrortemp[R[R[i]:R[i+1]-1]],/nan) $
                                  / total(obsmasktemp[R[R[i]:R[i+1]-1]]/$
                                          obsvarerrortemp[R[R[i]:R[i+1]-1]],/nan)
                    
                    tmperr[i+omin] = total(1./obsvarerrortemp[R[R[i]:R[i+1]-1]] $
                                           * float(obsmasktemp[R[R[i]:R[i+1]-1]]),/nan) $
                                     / total(1./obsvarerrortemp[R[R[i]:R[i+1]-1]]^2 * $
                                             float(obsmasktemp[R[R[i]:R[i+1]-1]]),/nan) * errscale
                 ENDIF ELSE BEGIN
                    
                    ;; variables are averaged linearly
                    tmp[i+omin] = total(obsvartemp[R[R[i]:R[i+1]-1]] $
                                        * float(obsmasktemp[R[R[i]:R[i+1]-1]]),/nan) $
                                  / total(float(obsmasktemp[R[R[i]:R[i+1]-1]]))
                    
                    tmperr[i+omin] = total(obsvarerrortemp[R[R[i]:R[i+1]-1]] $
                                           * float(obsmasktemp[R[R[i]:R[i+1]-1]]),/nan) $
                                     / total(float(obsmasktemp[R[R[i]:R[i+1]-1]])) * errscale

                 ENDELSE
              ENDIF
              
           ENDIF
        ENDFOR ;; 2D histogram loop

        obsvar[*,*,t] = tmp
        obsvarerror[0,*,*,t] = tmp + tmperr
        obsvarerror[1,*,*,t] = tmp - tmperr

     ENDIF

  ENDFOR  ;; DoY loop

END
