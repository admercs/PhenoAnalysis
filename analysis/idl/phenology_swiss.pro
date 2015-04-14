PRO phenology_swiss

  COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, $
     range, xrange, yrange, length, nyears, legenddx, legenddy, aspect


  ;; this code read and plots the climate sensitivity experiments for
  ;; evaluating the light-limitation of temperate vegetation over Switzerland

  ;; flags
  !except=2 ;; print fpe's or not
  write=0   ;; write ascii files with statistics for plotting in R

  ;; directories
  basedir = '/Users/stockli/pheno_analysis-output/'
  datadir = '/Users/stockli/publications/pheno_light/data/'

  ;; case names
  sitename = 'Switzerland'
  model = 'analysis'
  experiments = ['Swiss-Prediction-0K','Swiss-Prediction-2K','Swiss-Prediction-4K','Swiss-Prediction-6K']

  ;; time extent
  years=indgen(53)+1958

  ;; parameters
  FNAN = !values.f_nan
  nodata = -9999.0
  thval = 0.35
  laierr = 0.34

  nyears = n_elements(years)
  nexp = n_elements(experiments)
  ndays = 366

  lai = make_array(ndays,nyears,nexp,/float,value=FNAN)
  temp_fac = make_array(ndays,nyears,nexp,/float,value=FNAN)
  moist_fac = make_array(ndays,nyears,nexp,/float,value=FNAN)
  light_fac = make_array(ndays,nyears,nexp,/float,value=FNAN)

  ;; read output data
  FOR y = 0,nyears-1 DO BEGIN

     FOR e = 0,nexp-1 DO BEGIN
        year = years[y]
        syear = string(year,format='(I4.4)')

        read_pheno_local,basedir,1,experiments[e],sitename,year,'LAI',86400, $
                         data,error,varlongname,varunits, 0, 0, $
                         modelshort,modelname, nodata, modelok

        lai[0:n_elements(data)-1,y,e] = reform(data)

        read_pheno_local,basedir,1,experiments[e],sitename,year,'TEMP_FAC',86400, $
                         data,error,varlongname,varunits, 0, 0, $
                         modelshort,modelname, nodata, modelok

        temp_fac[0:n_elements(data)-1,y,e] = reform(data)

        read_pheno_local,basedir,1,experiments[e],sitename,year,'MOIST_FAC',86400, $
                         data,error,varlongname,varunits, 0, 0, $
                         modelshort,modelname, nodata, modelok

        moist_fac[0:n_elements(data)-1,y,e] = reform(data)

        read_pheno_local,basedir,1,experiments[e],sitename,year,'LIGHT_FAC',86400, $
                         data,error,varlongname,varunits, 0, 0, $
                         modelshort,modelname, nodata, modelok

        light_fac[0:n_elements(data)-1,y,e] = reform(data)

     ENDFOR

  ENDFOR

  ;; derive modeled SOS 
  sos_m = make_array(nyears,nexp,/float,value=FNAN)
  sos_m_e = make_array(nyears,nexp,/float,value=FNAN)
  
  FOR e=0,nexp-1 DO BEGIN

     ;; version 1
     threshold=(max(lai[*,*,e],/nan)-min(lai[*,*,e],/nan))*thval + min(lai[*,*,e],/nan)

     ;; version 2
     threshold=(max(mean(lai[*,*,e],dimension=2,/nan),/nan)-min(mean(lai[*,*,e],dimension=2,/nan),/nan))*thval $
               + min(mean(lai[*,*,e],dimension=2,/nan),/nan)

     FOR y=0,nyears-1 DO BEGIN           
        ;; mean
        isgrowing=where(lai[30:ndays-1,y,e] GT threshold,growcount)
        if growcount GT 1 THEN sos_m[y,e]=isgrowing[0] + 30

        ;; error
        isgrowing=where(lai[30:ndays-1,y,e]+laierr GT threshold,growcount)
        if growcount GT 1 THEN sos_m_e[y,e]=(isgrowing[0] + 30 - sos_m[y,e]) * 0.5
        isgrowing=where(lai[30:ndays-1,y,e]-laierr GT threshold,growcount)
        if growcount GT 1 THEN sos_m_e[y,e]=sos_m_e[y,e] + (isgrowing[0] + 30 - sos_m[y,e]) * 0.5
     ENDFOR
  ENDFOR

  ;; the lai error translates to 0-5 day SOS error
  sos_m_e[*] = 5.0

  ;; read observed SOS
  sos_o = make_array(nyears,/float,value=FNAN)
  sos_o_e = make_array(nyears,/float,value=0.)

  temp=''
  y=0

  IF 1 EQ 0 THEN BEGIN
     dates=fltarr(4)
     openr,lun,'/Users/stockli/modelfarm/tower_data/selected/phenology/Swiss_Lowland/INDEX.2006.uncert.dat',/get_lun
     readf,lun,temp
     
     WHILE EOF(lun) EQ 0 DO BEGIN
        readf,lun,temp
        dates[*] = strsplit(temp,';',/extract,/preserve_null)
        yi = where(years eq round(dates[0]),ycount)
        IF ycount GT 0 THEN BEGIN
           sos_o[yi] = dates[1]
           sos_o_e[yi] = (dates[2]-dates[1])/2.
        ENDIF
     END
     free_lun,lun
  ENDIF ELSE BEGIN
     temp=''
     year=0
     sos=0
     openr,lun,'/Users/stockli/modelfarm/tower_data/selected/phenology/Swiss_Lowland/RutishauserEtAl_2007_StatSpringPlantReconstruction_JGR_updated_2010.dat',/get_lun
     readf,lun,temp

     WHILE EOF(lun) EQ 0 DO BEGIN
        readf,lun,year,sos
        yi = where(years eq year,ycount)
        IF ycount GT 0 THEN BEGIN
           sos_o[yi] = sos
           sos_o_e[yi] = 5.4
        ENDIF
     END
     free_lun,lun
  ENDELSE

  ;; plot stuff
  
  ;; set up display
  !P.FONT=1
  device,decomposed=0,set_font='HELVETICA',/tt_font

  xsize = 1000
  ysize = 500
  aspect = 2.0
  charthick = 1.0
  charsize = 1.0
  thick = 1.0
  colored=1
  legend=0

  loadct,50,file='brewer.tbl'
  tvlct,red,green,blue,/get
  green=shift(green,1)
  blue=shift(blue,1)
  red=shift(red,1)
  red[0]=0
  green[0]=0
  blue[0]=0
  ntable=!d.table_size
  red[ntable-2]=150B
  green[ntable-2]=150B
  blue[ntable-2]=150B
  red[ntable-1]=255B
  green[ntable-1]=255B
  blue[ntable-1]=255B
  tvlct,red,green,blue


  ;; plot time series
  t = julday(1,1,years,0,0,0)
  x = sos_o
  x_e = sos_o_e
  y = sos_m
  y_e = sos_m_e
  ytitle = 'Start of Season (days)'
  xtitle = 'Years'
  xrange=[min(t)-365.,max(t)+365.]
  yrange=[min([min(x-x_e,/nan),min(y,/nan)],/nan)-5,$
          max([max(x+x_e,/nan),max(y,/nan)],/nan)+5]

  result = LABEL_DATE(date_format="%Y")

  window,0,xsize=xsize,ysize=ysize
  plot,t,replicate(0,nyears),xstyle=9,ystyle=9, xrange=xrange, yrange=yrange, /nodata, xthick=2.0, ythick=2.0, $
       thick=thick,charsize=2*charsize,charthick=charthick,xminor=1,yminor=1, xticks = 15, $
       xtitle=xtitle, ytitle = ytitle, color = ntable-1, xtickformat = 'LABEL_DATE', $
       position=[0.18,0.10,0.95,0.92], title = sitename

;  x = sos_o[52-30:52]

  good = where(finite(x),goodcount)
  a=regress(indgen(goodcount),x[good],const=b,mcorrelation=r, $
            measure_errors=replicate(0.5,goodcount),/double, ftest=fvalue)
  print,'OBS trend '+string(a)

;  y = sos_m[52-30:52,*]

  FOR e=0,nexp-1 DO BEGIN 
     overplot,t,y[*,e],experiments[e],e+1,1
     overplot,t,transpose([[y[*,e]-y_e[*,e]],[y[*,e]+y_e[*,e]]]),experiments[e],e+1,1

     good = where(finite(y[*,e]),goodcount)
     a=regress(indgen(goodcount),y[good,e],const=b,mcorrelation=r, $
               measure_errors=replicate(0.5,goodcount),/double, ftest=fvalue)
     print,experiments[e]+' trend '+string(a)

     index = where(finite(x) AND finite(y[*,e]),count)
     
     theta=acos(correlate(x[index],y[index,e]))
     rms=sqrt(total((y[index,e] - x[index])^2.0)/float(count))
     bias = total(y[index,e] - x[index])/float(count)
     
     r=correlate(x[index],y[index,e])
     var=1/(count-3.0)
     zvalue = 0.5*alog((1+r)/(1-r))/sqrt(var)
     if (zvalue lt 0) then pvalue = 2*(gauss_pdf(zvalue)) $
     else pvalue = 2*(1-gauss_pdf(zvalue))

     print,experiments[e]+' correlation: ',cos(theta), r,' bias: ',bias,' rms: ',rms, ' pvalue: ',pvalue
  ENDFOR

  overplot,t,transpose([[x-x_e],[x+x_e]]),'OBS',0,1
  overplot,t,x,'OBS',0,1


  ;; save data
  IF write THEN BEGIN

     openw,lun,datadir+'sos_'+sitename+'.dat', width=120,/get_lun
     printf,lun,"Year OBS ",experiments
     FOR z=0,nyears-1 DO BEGIN
        printf,lun,years[z],sos_o[z],reform(sos_m[z,*])
     ENDFOR
     free_lun,lun

     openw,lun,datadir+'temp_fac_'+sitename+'.dat', width=120,/get_lun
     printf,lun,"Year Doy ",experiments
     FOR z=0,nyears-1 DO BEGIN
        FOR d=0,ndays-1 DO BEGIN
           printf,lun,years[z],d+1,reform(temp_fac[d,z,*])
        ENDFOR
     ENDFOR
     free_lun,lun

     openw,lun,datadir+'light_fac_'+sitename+'.dat', width=120,/get_lun
     printf,lun,"Year Doy ",experiments
     FOR z=0,nyears-1 DO BEGIN
        FOR d=0,ndays-1 DO BEGIN
           printf,lun,years[z],d+1,reform(light_fac[d,z,*])
        ENDFOR
     ENDFOR
     free_lun,lun

     openw,lun,datadir+'moist_fac_'+sitename+'.dat', width=120,/get_lun
     printf,lun,"Year Doy ",experiments
     FOR z=0,nyears-1 DO BEGIN
        FOR d=0,ndays-1 DO BEGIN
           printf,lun,years[z],d+1,reform(moist_fac[d,z,*])
        ENDFOR
     ENDFOR
     free_lun,lun

  ENDIF

  good = where(finite(x),goodcount)
  a=regress(indgen(goodcount),x[good],const=b,mcorrelation=r, $
            measure_errors=replicate(0.5,goodcount),/double)
  print,'OBS trend '+string(a)

  ;; plot SOS variance
  histrange=[min([min(x,/nan),min(y,/nan)],/nan)-5,$
             max([max(x,/nan),max(y,/nan)],/nan)+5]

  binsize = 1.
  h_o = histogram(x,min=histrange[0],max=histrange[1],binsize = binsize,locations=locations)
  
  m_o = mean(x,/nan)
  s_o = stddev(x,/nan)


  nbins = n_elements(h_o)
  h_m = make_array(nbins,nexp,/float,value=0)
  m_m = make_array(nexp,/float,value=0)
  s_m = make_array(nexp,/float,value=0)
  FOR e=0,nexp-1 DO BEGIN
     h_m[*,e] = histogram(y[*,e],min=histrange[0],max=histrange[1],nbins=nbins,locations=locations)
     m_m[e] = mean(y[*,e],/nan)
     s_m[e] = stddev(y[*,e],/nan)
  ENDFOR

  t = julday(1,1,2009,0,0,0) + (findgen(nbins)*float(binsize) + histrange[0])

  xrange=[t[0]-2.0*binsize,t[nbins-1]+2.0*binsize]
  yrange=[min(h_m)-1.0,max(h_m)+1.0]
  ytitle = 'Number of Occurences'
  xtitle = 'Time'

  window,1,xsize=xsize,ysize=ysize
  result = LABEL_DATE(date_format="%D %M")
  plot,t,replicate(0,nbins),xstyle=9,ystyle=9, xrange=xrange, $
       yrange=yrange, /nodata, xthick=2.0, ythick=2.0, $
       thick=thick,charsize=2*charsize,charthick=charthick,xminor=5*binsize, yminor=1, xticks = (nbins+4)/5, $
       xtitle='Time', ytitle = ytitle, color = ntable-1, xtickformat = 'LABEL_DATE', $
       position=[0.18,0.10,0.95,0.92], title = sitename

  ;; plot bars
  FOR b=0,nbins-1 DO plots,[t[b]-0.2*binsize,t[b]-0.2*binsize,t[b]+0.2*binsize,t[b]+0.2*binsize], $
                           [0,h_o[b],h_o[b],0], thick=thick,color=ntable-1,/data
  
  FOR e=0,nexp-1 DO BEGIN
     FOR b=0,nbins-1 DO plots,[t[b]-0.2*binsize,t[b]-0.2*binsize,t[b]+0.2*binsize,t[b]+0.2*binsize]+0.1*float(1+e), $
                              [0,h_m[b,e],h_m[b,e],0], thick=thick,color=e+1,/data
  ENDFOR

  ;; plot gaussian distributions
  dt = 0.1
  nt = round((histrange[1]-histrange[0]+1)/dt)
  t_g = findgen(nt)*dt+histrange[0]

  start_params = [max(h_o),m_o,s_o]
  nparams = n_elements(start_params)
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, nparams)
  parinfo[0].limited = [1,1]
  parinfo[0].limits  = [start_params[0]*0.1,start_params[0]*10.0]
  parinfo[0].fixed = 0
  parinfo[1].limited = [1,1]
  parinfo[1].limits  = [start_params[1]-50.,start_params[1]+50.0]
  parinfo[1].fixed = 1
  parinfo[2].limited = [1,1]
  parinfo[2].limits  = [start_params[2]*0.1,start_params[2]*10.0]
  parinfo[2].fixed = 1
  
  parinfo.value = start_params
  params = MPFITFUN("gauss_fun",locations, h_o, parinfo=parinfo)
  
  g_o = gauss_fun(t_g,params)  
  plots,t_g+julday(1,1,2009,0,0,0),g_o, thick=thick*2.0,color=ntable-1,/data

  FOR e=0,nexp-1 DO BEGIN
     start_params = [max(h_m[*,e]),m_m[e],s_m[e]]
     nparams = n_elements(start_params)
     parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                          limits:[0.D,0]}, nparams)
     parinfo[0].limited = [1,1]
     parinfo[0].limits  = [start_params[0]*0.1,start_params[0]*10.0]
     parinfo[0].fixed = 0
     parinfo[1].limited = [1,1]
     parinfo[1].limits  = [start_params[1]-50.,start_params[1]+50.0]
     parinfo[1].fixed = 1
     parinfo[2].limited = [1,1]
     parinfo[2].limits  = [start_params[2]*0.1,start_params[2]*10.0]
     parinfo[2].fixed = 1
     
     parinfo.value = start_params
     params = MPFITFUN("gauss_fun",locations, h_m[*,e], parinfo=parinfo)
     
     g_m = gauss_fun(t_g,params)
     plots,t_g+julday(1,1,2009,0,0,0),g_m, thick=thick*2.0,color=1+e,/data
  ENDFOR

  ;; check the average influence of light limitation on the total
  ;; growth constraint
  FOR e=0,nexp-1 DO BEGIN

     amount = make_array(nyears,/float,value=!values.f_nan)
     FOR z=0,nyears-1 DO BEGIN
        d = round(y[z,e])
        amount[z] = 100.*(1.-light_fac[d,z,e]) / ((1.-light_fac[d,z,e]) + (1.-temp_fac[d,z,e]) + (1.-moist_fac[d,z,e]))
     ENDFOR

     print,experiments[e],":",mean(amount,/nan),s_m[e]

  ENDFOR
  
  ;; plot evolution of temperature, moisture and light limitation

  t = julday(1,1,2009,0,0,0) + findgen(ndays)

  xrange=[t[0]-0.5,t[ndays-2]]
  yrange=[-0.1,1.1]
  ytitle = 'Growth Constraint'
  xtitle = 'Time'

  window,2,xsize=xsize,ysize=ysize
  result = LABEL_DATE(date_format="%M")
  plot,t,replicate(0,ndays),xstyle=9,ystyle=9, xrange=xrange, $
       yrange=yrange, /nodata, xthick=2.0, ythick=2.0, $
       thick=thick,charsize=2*charsize,charthick=charthick, xminor=1, yminor=1, xtickinterval = 30.5, $
       xtitle='Time', ytitle = ytitle, color = ntable-1, xtickformat = 'LABEL_DATE', $
       position=[0.18,0.10,0.95,0.92], title = sitename

  nyears_save = nyears
  nyears = 1

  t_m = mean(light_fac[*,*,0],dimension=2,/nan)
  
  overplot,t,t_m,experiments[0],0,1

  FOR e=0,nexp-1 DO BEGIN

     t_m = mean(temp_fac[*,*,e],dimension=2,/nan)
     t_s = stddev(temp_fac[*,*,e],dimension=2,/nan)
     
     overplot,t,t_m,experiments[e],e+1,1
     overplot,t,transpose([[t_m-t_s],[t_m+t_s]]),experiments[e],e+1,1


  ENDFOR
  nyears = nyears_save


  stop
  
END
