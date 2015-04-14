PRO phenology_global

  ;; this script plots a global analysis of temperature versus light
  ;; limitation factors

  COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, range, $
     xrange, yrange, length, nyears, legenddx, legenddy, aspect, pollon, pollat, psout, $
     ysize_ps_default, ysize_x_default, xsize, ysize
  
  !except=2

  ;; this code read and plots the climate sensitivity experiments for
  ;; evaluating the light-limitation of temperate vegetation over Switzerland

;  basedir = '/Users/stockli/pheno_analysis-output/'
;  plotdir = '/Users/stockli/publications/pheno_light/maps-test/'
;  datadir = '/Users/stockli/'
;  writedir = '/Users/stockli/publications/pheno_light/data/'

  basedir = '/store/users/stockli/pheno_analysis-output/'
  plotdir = '/scratch/eiger/julier/stockli/temp/'
  datadir = '/store/users/stockli/'
  writedir = '/scratch/eiger/julier/stockli/temp/'

  sitename = 'Global'
  region = 'global_0.1'
  model = 'analysis'

  experiments = ['Global-Prediction-0K','Global-Prediction-2K','Global-Prediction-4K','Global-Prediction-6K']
  plotnames = ['+0K','+2K','+4K','+6K']

;  experiments = ['Global-Prediction-6K']
;  plotnames = ['+6K']

;  experiments = ['Global-Prediction-0K','Global-Prediction-6K']
;  plotnames = ['+0K','+6K']

  years=1959 + indgen(2010-1959+1)
;  years=[2007,2008,2009]
;  years=[2009]
;  years=[2002,2003,2004,2005,2006,2007,2008,2009]

  psout = 1 ;; plot to EPS/PDF/PNG instead of X11
  usecontour = 1 ;; use contour map instead of pixel plot
  plotproj = 1  ;; reproject regular grid to Mercator projection [partly INOP]

  plotlight = 0 ;; plot relative photoperiod limitation
  plotsos = 0   ;; plot earliest SOS dates 
  plotvar = 1   ;; plot SOS date variability
  plotrelvar = 0 ;; plot relative SOS date variability
  screenvarlight = 0 ;; screen above plots by 50% relative limitation and 50% relative SOS date variability
  contourvarlight = 0 ;; plot contour line of area with at least 50% relative limitation and less than 50% relative SOS date variability


  subregion = "Northern_Hemisphere"
  lonmin = -130.0
  lonmax = 150.0
  latmin = 25.0
  latmax = 65.0
  dlon = 1.0
  dlat = 1.0

;  subregion = "Swiss"
;  lonmin = 8.5
;  lonmax = 8.5
;  latmin = 47.5
;  latmax = 47.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "Morgan_Monroe"
;  lonmin = -86.5
;  lonmax = -86.5
;  latmin = 39.5
;  latmax = 39.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "Harvard_Forest"
;  lonmin = -72.5
;  lonmax = -72.5
;  latmin = 42.5
;  latmax = 42.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "San_Marino"
;  lonmin = 12.5
;  lonmax = 12.5
;  latmin = 43.5
;  latmax = 43.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "Stockholm"
;  lonmin = 18.5
;  lonmax = 18.5
;  latmin = 60.5
;  latmax = 60.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "Tokyo"
;  lonmin = 139.5
;  lonmax = 139.5
;  latmin = 36.5
;  latmax = 36.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "Bordeaux"
;  lonmin = -0.5
;  lonmax = -0.5
;  latmin = 44.5
;  latmax = 44.5
;  dlon = 1.0
;  dlat = 1.0
  
;  subregion = "Georgia"
;  lonmin = -83.5
;  lonmax = -83.5
;  latmin = 34.5
;  latmax = 34.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "London"
;  lonmin = -0.5
;  lonmax = -0.5
;  latmin = 51.5
;  latmax = 51.5
;  dlon = 1.0
;  dlat = 1.0

;  subregion = "Frankfurt"
;  lonmin = 8.5
;  lonmax = 8.5
;  latmin = 50.5
;  latmax = 50.5
;  dlon = 1.0
;  dlat = 1.0

  ;; parameters
  FNAN = !values.f_nan
  nodata = -9999.0
  thval = 0.35
  laierr = 0.34


  ;; start processing here
  
  nyears = n_elements(years)
  nexp = n_elements(experiments)

  ;; check whether we have a single gridpoint
  IF (lonmin EQ lonmax) AND (latmin EQ latmax) THEN BEGIN
     is_site = 1
     nlon = 1
     nlat = 1s
     lons = lonmin
     lats = latmin
  ENDIF ELSE BEGIN
     is_site = 0
     nlon = round((lonmax-lonmin)/dlon)
     nlat = round((latmax-latmin)/dlat)
     lons = findgen(nlon)*dlon + lonmin + 0.5*dlon ; W->E
     lats = findgen(nlat)*dlat + latmin + 0.5*dlat ; S->N
  ENDELSE


  ;; read pft distribution of plotting grid           
  read_pftdata, datadir, 'pft', region, lons, lats, pft, dlon=dlon, dlat=dlat, pftnames=pftnames

  ;; only plot areas where temperate forests cover more than 5% of the area

  ;; temperate
;  mask = (pft[*,*,1]+pft[*,*,5]+pft[*,*,7]+pft[*,*,10]+pft[*,*,13]) GT 10.0

  ;; temperate + boreal
;  mask = (pft[*,*,1]+pft[*,*,2]+pft[*,*,3]+pft[*,*,5]+pft[*,*,7]+ $
;          pft[*,*,8]+pft[*,*,10]+pft[*,*,11]+pft[*,*,13]) GT 5.0

  ;; temperate + boreal + arctic
;  mask = (pft[*,*,1]+pft[*,*,2]+pft[*,*,3]+pft[*,*,5]+pft[*,*,7]+ $
;          pft[*,*,8]+pft[*,*,10]+pft[*,*,11]+pft[*,*,12]) GT 5.0

  ;; tall trees
  mask = (pft[*,*,1]+pft[*,*,2]+pft[*,*,3]+pft[*,*,4]+pft[*,*,5]+ $
          pft[*,*,6]+pft[*,*,7]+pft[*,*,8]+pft[*,*,9]) GT 5.0

;  mask[*] = 1

  dtmod=!values.f_nan           ; we do not know model time step yet
  ndays = 366

  FOR p=0,n_elements(pftnames)-1 DO BEGIN
     print,pftnames[p],mean(pft[*,*,p])
  ENDFOR

  ;; set up global arrays
  sos = make_array(nlon,nlat,nyears,nexp,/float,value=!values.f_nan)

  IF is_site THEN BEGIN
     light = make_array(nlon,nlat,ndays,nyears,nexp,/float,value=!values.f_nan)
     temperature = make_array(nlon,nlat,ndays,nyears,nexp,/float,value=!values.f_nan)
     moisture = make_array(nlon,nlat,ndays,nyears,nexp,/float,value=!values.f_nan)
  ENDIF ELSE BEGIN
     light = make_array(nlon,nlat,nyears,nexp,/float,value=!values.f_nan)
     temperature = make_array(nlon,nlat,nyears,nexp,/float,value=!values.f_nan)
     moisture = make_array(nlon,nlat,nyears,nexp,/float,value=!values.f_nan)
     percent = make_array(nlon,nlat,nyears,nexp,/float,value=!values.f_nan)
  ENDELSE

  laimin = make_array(nlon,nlat,nyears,/float,value=!values.f_nan)
  laimax = make_array(nlon,nlat,nyears,/float,value=!values.f_nan)
  laijan = make_array(nlon,nlat,nyears,/float,value=!values.f_nan)
  laiaug = make_array(nlon,nlat,nyears,/float,value=!values.f_nan)

  ;; experiments
  FOR e = 0,nexp-1 DO BEGIN
     
     ;; loop for laimin/laimax calculation
     FOR l = 0, 1 DO BEGIN

        ;; year loop
        FOR y = 0,nyears-1 DO BEGIN
           
           read_pheno_regional,basedir,1,experiments[e],sitename,years[y],'LAI',dtmod, $
                               lai,error,varname,varunits,varlongname, $
                               modlons, modlats, modelname,modellongname, nodata, has_model, $
                               nlonmod, nlatmod, 1, 1, 0, 0, $
                               lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax
           
           IF l EQ 1 THEN BEGIN
              read_pheno_regional,basedir,1,experiments[e],sitename,years[y],'TEMP_FAC',dtmod, $
                                  temp_fac,error,varname,varunits,varlongname, $
                                  modlons, modlats, modelname,modellongname, nodata, has_model, $
                                  nlonmod, nlatmod, 1, 1, 0, 0, $
                                  lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax
              
              read_pheno_regional,basedir,1,experiments[e],sitename,years[y],'LIGHT_FAC',dtmod, $
                                  light_fac,error,varname,varunits,varlongname, $
                                  modlons, modlats, modelname,modellongname, nodata, has_model, $
                                  nlonmod, nlatmod, 1, 1, 0, 0, $
                                  lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax
              
              read_pheno_regional,basedir,1,experiments[e],sitename,years[y],'MOIST_FAC',dtmod, $
                                  moist_fac,error,varname,varunits,varlongname, $
                                  modlons, modlats, modelname,modellongname, nodata, has_model, $
                                  nlonmod, nlatmod, 1, 1, 0, 0, $
                                  lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax
           ENDIF
           
           IF l EQ 0 THEN BEGIN
              ;; get minimum and maximum LAI of each year
              laimin[*,*,y] = min(lai,dimension=3)
              laimax[*,*,y] = max(lai,dimension=3)
              laijan[*,*,y] = lai[*,*,30]
              laiaug[*,*,y] = lai[*,*,220]

           ENDIF ELSE BEGIN
              ;; derive SOS
              IF nyears EQ 1 THEN BEGIN
                 mmin = laimin
                 mmax = laimax
                 mjan = laijan
                 maug = laiaug
              ENDIF ELSE BEGIN
                 mmin = mean(laimin,dimension=3)
                 mmax = mean(laimax,dimension=3)
                 mjan = mean(laijan,dimension=3)
                 maug = mean(laiaug,dimension=3)
              ENDELSE

              idx = where(finite(mmin) AND finite(mmax) AND finite(mjan) AND finite(maug) AND $
                          ((mmax - mmin) GT 1.0) AND (mjan LT 2.0) AND (maug GT 2.0) ,count)

              threshold = make_array(nlon,nlat,/float,value=9999.0)
              IF count GT 0 THEN threshold[idx]=(mmax[idx]-mmin[idx])*thval + mmin[idx]
              
              ndaystemp = (size(lai,/dimensions))[2]

              temp = make_array(nlon,nlat,/float,value=0.0)
              FOR d=30,ndaystemp-1 DO BEGIN
                 idx=where(lai[*,*,d] GT threshold,count)
                 if count GT 0 THEN temp[idx] = (temp[idx] EQ 0.0)*float(d) + (temp[idx] NE 0.0)*temp[idx]
              ENDFOR
              idx = where((temp EQ 0.0) OR (temp LT 35),count)
              IF count GT 0 THEN temp[idx] = !values.f_nan
              sos[*,*,y,e] = temp

              IF is_site THEN BEGIN
                 ;; save daily limitation factors
                 temperature[*,*,0:ndaystemp-1,y,e] = temp_fac
                 light[*,*,0:ndaystemp-1,y,e] = light_fac
                 moisture[*,*,0:ndaystemp-1,y,e] = moist_fac
              ENDIF ELSE BEGIN
                 ;; derive temp/light/moist limitation at SOS
                 idx = where(finite(temp),count)
                 idx2d = array_indices([nlon,nlat],idx,/dimensions)
                 IF count GT 0 THEN BEGIN
                    FOR i=0,count-1 DO BEGIN
                       temperature[idx2d[0,i],idx2d[1,i],y,e] = temp_fac[idx2d[0,i],idx2d[1,i],fix(temp[idx[i]])]
                       light[idx2d[0,i],idx2d[1,i],y,e] = light_fac[idx2d[0,i],idx2d[1,i],fix(temp[idx[i]])]
                       moisture[idx2d[0,i],idx2d[1,i],y,e] = moist_fac[idx2d[0,i],idx2d[1,i],fix(temp[idx[i]])]
                    ENDFOR
                 ENDIF
              ENDELSE
           ENDELSE

        ENDFOR

     ENDFOR

  ENDFOR

  IF ~is_site THEN BEGIN
     
     ;; derive % light limitation at SOS date
     idx = where(finite(light) AND finite(temperature) AND finite(moisture),count)
     IF count GT 0 THEN BEGIN
        percent[idx] = 100.*(1.-light[idx]) / (((1.-light[idx]) + (1.-temperature[idx]) + (1.-moisture[idx]))>0.01)
     ENDIF

     FOR e = 0,nexp-1 DO BEGIN

        ;; prepare for graphics output
        IF psout THEN BEGIN
           
           ;; setting up plotting to PS device
           set_plot,'ps'

           xsize = 16.0

           IF plotproj THEN BEGIN
              ysize = xsize * nlat / nlon * 1.5
           ENDIF ELSE BEGIN
              ysize = xsize * nlat / nlon
           ENDELSE

           thick=ysize/4.0      ; line thickness
           charthick=2.0        ; char thickness
           charsize=ysize/3.0   ; char size
           
           !P.FONT=0
           
        ENDIF ELSE BEGIN
           
           ;; set up plotting window in X
           set_plot,'x'

           xsize = 1000L

           IF plotproj THEN BEGIN
              ysize = round(xsize * nlat / nlon * 1.5)
           ENDIF ELSE BEGIN
              ysize = xsize * nlat / nlon
           ENDELSE
           
           thick=1.0            ; line thickness
           charthick=1.5        ; char thickness 
           charsize=1.0         ; char size

           !P.FONT=1
           
           device,decomposed=0,set_font='HELVETICA',/tt_font
           
        ENDELSE

        IF plotsos THEN ctnum = 0 ELSE ctnum = 26

        loadct,ctnum,file='brewer.tbl' 

        tvlct,red,green,blue,/get
        ntable=!d.table_size
     
        IF ~plotsos THEN BEGIN
           red = reverse(red)
           green = reverse(green)
           blue = reverse(blue)
        ENDIF

        ;; red
        red[ntable-2]=255B
        green[ntable-2]=50B
        blue[ntable-2]=50B

        ;; grey
        red[ntable-3]=128B
        green[ntable-3]=128B
        blue[ntable-3]=128B
        
        ;; load color bar
        IF psout THEN BEGIN
           red[0]=255B
           green[0]=255B
           blue[0]=255B
           red[ntable-1]=0B
           green[ntable-1]=0B
           blue[ntable-1]=0B
           tvlct,red,green,blue
        ENDIF ELSE BEGIN
           red[0]=0
           green[0]=0
           blue[0]=0
           red[ntable-1]=255B
           green[ntable-1]=255B
           blue[ntable-1]=255B
           tvlct,red,green,blue
        ENDELSE

;  tvscl,congrid(pft[*,*,7] GT 5.0,nlon*2,nlat*2)
;  tv,congrid(bytscl(sos[*,*,0,0],30,200,top=252),nlon*4,nlat*4)
;  tv,congrid(bytscl(temperature[*,*,0,0],0,1,top=252),nlon*4,nlat*4)
;  tv,congrid(bytscl(light[*,*,0,0],0,1,top=252),nlon*4,nlat*4)
;  tv,congrid(bytscl(moisture[*,*,0,0],0,1,top=252),nlon*4,nlat*4)

;  tv,congrid(bytscl(moist_fac[*,*,90],0,1,top=252),nlon*4,nlat*4)
;  tv,congrid(bytscl(light_fac[*,*,90],0,1,top=252),nlon*4,nlat*4)
;  tv,congrid(bytscl(temp_fac[*,*,90],0,1,top=252),nlon*4,nlat*4)
;  tv,congrid(bytscl(percent[*,*,0,0],0,100,top=252),nlon*4,nlat*4)

;;  tv,congrid(bytscl(laimin,0,6,top=252),nlon*4,nlat*4)
;;  tv,congrid(bytscl(laimax,0,6,top=252),nlon*4,nlat*4)
;;  tv,congrid(bytscl(threshold,0,6,top=252),nlon*4,nlat*4)

        IF psout THEN BEGIN
           IF plotvar THEN filename = plotdir + 'sosvarmap_'+experiments[e]
           IF plotsos THEN filename = plotdir + 'sosmap_'+experiments[e]
           IF plotlight THEN filename = plotdir + 'lightmap_'+experiments[e]
           IF plotrelvar THEN filename = plotdir + 'sosrelvarmap_'+experiments[e]
           IF plotproj THEN filename += '_proj'
           device, /encapsul, bits_per_pixel=8,/color, font_size=5, $
                   filename=filename+'.eps', preview=1,/isolatin1, /helvetica, $
                   xsize=xsize, ysize=ysize, decomposed=0, language_level=2
        ENDIF ELSE BEGIN
           window,0,xsize=xsize,ysize=ysize        
        ENDELSE

        IF nyears GT 1 THEN BEGIN
           ;; standard deviation of SOS dates, including NAN treatment
           outvar = sqrt(total(sos[*,*,*,e]^2,3,/nan)/(total(finite(sos[*,*,*,e]),3)>1.0) - $
                         (total(sos[*,*,*,e],3,/nan)/(total(finite(sos[*,*,*,e]),3)>1.0))^2)
           idx = where(total(finite(sos[*,*,*,e]),3) EQ 0,count)
           IF count GT 0 THEN outvar[idx] = !values.f_nan

           ;; get minimum SOS date
;;         outsos = min(sos[*,*,*,e],dimension=3,/nan)
           
           ;; get 5% quantile of SOS dates
           outsos = make_array(nlon,nlat,/float,value=!values.f_nan)
           FOR b = 0L,nlat-1L DO BEGIN
              FOR a = 0L,nlon-1 DO BEGIN
                 outsos[a,b] = sos[a,b,(sort(reform(sos[a,b,*,e])))[0.05*nyears],e]
              ENDFOR
           ENDFOR

           outlight = total(percent[*,*,*,e],3,/nan)/(total(finite(percent[*,*,*,e]),3)>1.0)
           idx = where(total(finite(percent[*,*,*,e]),3) EQ 0,count)
           IF count GT 0 THEN outlight[idx] = !values.f_nan

           outvar0 = sqrt(total(sos[*,*,*,0]^2,3,/nan)/(total(finite(sos[*,*,*,0]),3)>1.0) - $
                          (total(sos[*,*,*,0],3,/nan)/(total(finite(sos[*,*,*,0]),3)>1.0))^2)
           idx = where(total(finite(sos[*,*,*,0]),3) EQ 0,count)
           IF count GT 0 THEN outvar0[idx] = !values.f_nan

           outrelvar = make_array(nlon,nlat,/float,value=!values.f_nan)
           idx = where(finite(outvar) AND finite(outvar0) AND (outvar0 NE 0.0),count)
           IF count NE 0 THEN outrelvar[idx] = outvar[idx] / outvar0[idx] * 100.

        ENDIF ELSE BEGIN
           outvar = make_array(nlon,nlat,/float,value=!values.f_nan)
           outsos = sos[*,*,*,e]
           outlight = percent[*,*,*,e]
           outrelvar = make_array(nlon,nlat,/float,value=!values.f_nan)
        ENDELSE

        IF plotsos THEN outdata = outsos
        IF plotvar THEN outdata = outvar
        IF plotrelvar THEN outdata = outrelvar
        IF plotlight THEN outdata = outlight

        dl = 0.05

        x0 = dl * ysize/xsize
        x1 = 1.0 - dl * ysize/xsize
        y0 = dl
        y1 = 1.0 - dl
        nx = (x1 - x0) * xsize
        ny = (y1 - y0) * ysize


        IF plotvar THEN BEGIN
           lmin = 0
           lmax = 9.
           nlev = 10.
        ENDIF        
        IF plotsos THEN BEGIN
           lmin = 70.
           lmax = 110.
           nlev = 10
        ENDIF
        IF plotlight THEN BEGIN
           lmin = 10.
           lmax = 90.
           nlev = 9
        ENDIF
        IF plotrelvar THEN BEGIN
           lmin = 10.
           lmax = 90.
           nlev = 9
        ENDIF

        ;; set up map projection
        IF plotproj THEN BEGIN
           map = map_proj_init(8,datum=8,limit=[latmin,lonmin,latmax,lonmax])
        ENDIF ELSE BEGIN
           map = map_proj_init(8,datum=8,limit=[latmin,lonmin,latmax,lonmax])
        ENDELSE

        IF usecontour THEN BEGIN

           dlev = (lmax-lmin)/float(nlev)        
           levels = findgen(nlev)*dlev + lmin

           idx = where(finite(outdata),count)
           IF count GT 0 THEN outdata[idx] = (outdata[idx]>(lmin+0.001))<(lmax-dlev)
           
           idx = where(mask EQ 0,count)
           IF count GT 0 THEN outdata[idx] = !values.f_nan

           IF screenvarlight THEN BEGIN
              ;; remove pixels where relative light limitation is too low or
              ;; relative SOS date variability is too high
              idx = where((outlight LT 50.) OR (outrelvar GT 50.),count)
              IF count GT 0 THEN outdata[idx] = !values.f_nan
              
           ENDIF

           IF plotproj THEN BEGIN

;              xy = map_proj_forward(lons#replicate(1.0,nlat),replicate(1.0,nlon)#lats,map_structure=map)
;              triangulate,reform(xy[0,*]),reform(xy[1,*]),tr
;              projdata=trigrid(reform(xy[0,*]),reform(xy[1,*]),outdata[*],tr,[0,0],map.uv_box,nx=2*nlon,ny=4*nlat, $
;                               missing=!values.f_nan)

              contour,outdata, levels=levels, $
                      c_colors = round((findgen(nlev) + 0.5)/float(nlev) * (ntable-5))+1B, /cell_fill, $
                      color=ntable-1, position=[x0, y0, x1, y1], /normal, $
                      xstyle=5,ystyle=5, c_charsize=charsize, c_thick=thick*1.5

           ENDIF ELSE BEGIN

              contour,outdata,levels=levels, $
                      c_colors = round((findgen(nlev) + 0.5)/float(nlev) * (ntable-5))+1B, /cell_fill, $
                      color=ntable-1, position=[x0, y0, x1, y1], /normal, $
                      xstyle=5,ystyle=5, c_charsize=charsize, c_thick=thick*1.5
           ENDELSE

        ENDIF ELSE BEGIN
           outbyte = make_array(nlon,nlat,/byte,value=0B)
           idx = where(finite(outdata),count)        
           IF count GT 0 THEN outbyte[idx] = (bytscl(outdata[idx],lmin,lmax,top=ntable-4)>1B) * mask[idx]
           IF psout THEN tv,congrid(outbyte,nlon*5,nlat*5),x0*xsize,y0*ysize,xsize=nx,ysize=ny,/centimeters $
           ELSE tv,congrid(outbyte,nx,ny),x0*xsize,y0*ysize,xsize=nx,ysize=ny,/device
        ENDELSE

        IF contourvarlight THEN BEGIN
           ;; highlight pixels where relative light limitation is high
           ;; enough and relative SOS date variability is low
           a = make_array(nlon,nlat,/float,value=0.)
           idx = where((outlight GE 50.) AND (outrelvar LE 50.),count)
           IF count GT 0 THEN a[idx] = 1.0

           contour,a,levels=[0.5],c_colors=[ntable-2],/overplot,c_thick=thick*1.5
           
        ENDIF

        xyouts,1.-2*x0,2*y0,plotnames[e],charsize=charsize*2.5,charthick=charthick*2, $
               align=1.0,color=ntable-1,/normal

        plot, map.uv_box[[0, 2]], map.uv_box[[1, 3]], Position=[x0, y0, x1, y1], $
              /nodata, XStyle=5, YStyle=5, /NoErase,/normal
        map_continents,color=ntable-3,mlinethick=thick*1.5,fill=0, map_structure=map

        plots,[x0,x0,x1,x1,x0],[y0,y1,y1,y0,y0],color=ntable-1,thick=thick*2.5,/normal

        IF psout THEN BEGIN
           device,/close      
           set_plot,'x'
           
           spawn,'convert -density 300 -resize 1600 -flatten +antialias '+filename+'.eps '+filename+'.png'
           spawn,'epstopdf '+filename+'.eps'

        ENDIF 
        
     ENDFOR

     IF plotvar THEN BEGIN
        colbar_standalone,divisions=nlev-1, units='Days', varname='SOS Date Variability', /ctrev, $
                          filename=plotdir+'sosvarmap.cb.eps',range=[lmin,lmax],psout=psout,ctnum=ctnum,/horizontal
     ENDIF
     IF plotsos THEN BEGIN
        colbar_standalone,divisions=nlev-1, units='Day Month', varname='Earliest SOS Date', $
                          filename=plotdir+'sosmap.cb.eps',range=[lmin,lmax],psout=psout,ctnum=ctnum,/doy,/horizontal
     ENDIF
     IF plotlight THEN BEGIN
        colbar_standalone,divisions=nlev-1, units='%', varname='Photoperiod limitation on green-up', /ctrev, $
                          filename=plotdir+'lightmap.cb.eps',range=[lmin,lmax],psout=psout,ctnum=ctnum,/horizontal
     ENDIF
     IF plotrelvar THEN BEGIN
        colbar_standalone,divisions=nlev-1, units='%', varname='Relative SOS Date Variability', /ctrev, $
                          filename=plotdir+'sosrelvarmap.cb.eps',range=[lmin,lmax],psout=psout,ctnum=ctnum,/horizontal
     ENDIF

  ENDIF ELSE BEGIN
     ;; print site time series and parameters to file

     openw,lun,writedir+'sos_'+subregion+'.dat', width=120,/get_lun
     printf,lun,"Year ",experiments
     FOR z=0,nyears-1 DO BEGIN
        printf,lun,years[z],reform(sos[0,0,z,*])
     ENDFOR
     free_lun,lun

     openw,lun,writedir+'temp_fac_'+subregion+'.dat', width=120,/get_lun
     printf,lun,"Year Doy ",experiments
     FOR z=0,nyears-1 DO BEGIN
        FOR d=0,ndays-1 DO BEGIN
           printf,lun,years[z],d+1,reform(temperature[0,0,d,z,*])
        ENDFOR
     ENDFOR
     free_lun,lun

     openw,lun,writedir+'light_fac_'+subregion+'.dat', width=120,/get_lun
     printf,lun,"Year Doy ",experiments
     FOR z=0,nyears-1 DO BEGIN
        FOR d=0,ndays-1 DO BEGIN
           printf,lun,years[z],d+1,reform(light[0,0,d,z,*])
        ENDFOR
     ENDFOR
     free_lun,lun

     openw,lun,writedir+'moist_fac_'+subregion+'.dat', width=120,/get_lun
     printf,lun,"Year Doy ",experiments
     FOR z=0,nyears-1 DO BEGIN
        FOR d=0,ndays-1 DO BEGIN
           printf,lun,years[z],d+1,reform(moisture[0,0,d,z,*])
        ENDFOR
     ENDFOR
     free_lun,lun

  ENDELSE

END
