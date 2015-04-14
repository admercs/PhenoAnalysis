PRO enkf

  !except = 2

  !p.thick=2
  !p.charsize=2
  !p.charthick=2
  xsize=800
  ysize=600
  device,decomposed=0,set_font='HELVETICA',/tt_font
  window,0,xsize=xsize,ysize=ysize
  loadct,50,file='../idl/brewer.tbl' 
  tvlct,red,green,blue,/get
  green=shift(green,1)
  blue=shift(blue,1)
  red=shift(red,1)
  red[0]=0
  green[0]=0
  blue[0]=0
  ntable=!d.table_size
  red[ntable-1]=255B
  green[ntable-1]=255B
  blue[ntable-1]=255B
  tvlct,red,green,blue

  libpath = '../idl_lib/'

  ndim = 3L
  nrens = 200L
  nrobs = 500L
  ncall = 1

  obsagg = 2

  single = 1

  w = fltarr(ndim)
  w[*] = [2.0,1.0,0.5]
  offset = [100.0,500.0,900.0]
  scale  = [100.0,100.0,100.0]
  obsmean = 700.0
  obsvar = 100.0
  obserror = 100.0

  w /= total(w)

  
  seed = systime(1)

  seed = 234

  states = make_array(ndim,nrens,/float,value=0.0)
  FOR d=0,ndim-1 DO states[d,*] = randomn(seed,nrens,/normal) * scale[d] + offset[d]

  states_prior = states

  FOR c=1,ncall DO BEGIN

     CASE obsagg OF

        0: BEGIN
           ;; disaggregation of observations to individual states by
           ;; weighting observations
           
           HA = make_array(nrobs*ndim,nrens,/float,value=0.0)
           obsdat = make_array(nrobs*ndim,/float,value=0.0)
           obserr = make_array(nrobs*ndim,/float,value=0.0)
           p = 0
           FOR o=0,nrobs-1 DO BEGIN
              FOR d=0,ndim-1 DO BEGIN
                 
                 obsdat[p] = (randomn(seed,/normal) * obsvar + obsmean) * w[d]
                 obserr[p] = obserror * sqrt(w[d])
                 
                 val_mean = mean(states[d,*])
;           val_var = variance(states[d,*])
                 
;           HA[p,*] = states[d,*]
                 
                 HA[p,*] = (states[d,*] - val_mean) * sqrt(w[d])
                 
;;           HA[p,*] *= sqrt(val_var * w[d] / variance(HA[p,*]))
                 
                 HA[p,*] += (val_mean * w[d])
                 
                 p += 1
              ENDFOR
           ENDFOR
           
           nrobs *= ndim
           
        END

        1: BEGIN

           ;; linear aggregation of centered ensemble means and variances
           
           HA = make_array(nrobs,nrens,/float,value=0.0)
           obsdat = make_array(nrobs,/float,value=0.0)
           obserr = make_array(nrobs,/float,value=0.0)
           
           val_mean = fltarr(ndim)
;;     val_var = fltarr(ndim)
           
           FOR o=0,nrobs-1 DO BEGIN
              obsdat[o] = randomn(seed,/normal) * obsvar + obsmean
              obserr[o] = obserror
              
              FOR d=0,ndim-1 DO BEGIN
                 val_mean[d] = mean(states[d,*])
;;           val_var[d] = variance(states[d,*])
              ENDFOR
              
              FOR r=0,nrens-1 DO BEGIN
                 HA[o,r] = total((states[*,r] - val_mean)*sqrt(w))
              ENDFOR
              
;;        HA[o,*] *= sqrt(total(val_var * w ) / variance(HA[o,*]))
              
              HA[o,*] += total(val_mean * w)
              
           ENDFOR
        END

        2: BEGIN
           ;; linear aggregation of individual state ensemble members
           
           HA = make_array(nrobs,nrens,/float,value=0.0)
           obsdat = make_array(nrobs,/float,value=0.0)
           obserr = make_array(nrobs,/float,value=0.0)

           nsample = 1000

           FOR o=0,nrobs-1 DO BEGIN
              obsdat[o] = randomn(seed,/normal) * obsvar + obsmean
              obserr[o] = obserror

              idx = fix(randomu(seed,nsample,ndim) * nrens)

              FOR r=0,nrens-1 DO BEGIN
                 HA[o,r] += total(states[*,r] * w)
              ENDFOR

           ENDFOR

        END

     ENDCASE

     IF c EQ 1 THEN HA_save = HA

     IF single THEN BEGIN

        status = CALL_EXTERNAL(libpath+'enkf_wrapper_r4.so', 'enkf_wrapper_r4', $
                               ndim, nrens, nrobs, states, HA, obsdat, obserr, $
                               VALUE=[0,0,0,0,0,0,0],/i_value,/unload)

     ENDIF ELSE BEGIN

        dstates = double(states)
        dHA = double(HA)
        dobsdat = double(obsdat)
        dobserr = double(obserr)
        
        status = CALL_EXTERNAL(libpath+'enkf_wrapper_r8.so', 'enkf_wrapper_r8', $
                               ndim, nrens, nrobs, dstates, dHA, dobsdat, dobserr, $
                               VALUE=[0,0,0,0,0,0,0],/i_value,/unload)
        
        states = float(dstates)

     ENDELSE

  ENDFOR

  nstep = 100
  xrange = [min([states[*],states_prior[*],obsdat]),max([states[*],states_prior[*],obsdat])]
  step = (xrange[1] - xrange[0])/float(nstep)
  xval = findgen(nstep+1)*step + xrange[0]
  
  yrange = [0,nrens/5]

  plot,xval,xval, xrange=xrange,yrange=yrange,color=ntable-1,/nodata
  ly = 0.85
  FOR d=0,ndim-1 DO BEGIN
     oplot,xval,histogram(states_prior[d,*],binsize=step,min=xrange[0],max=xrange[1]),$
           color=2,linestyle=d*2
     oplot,xval,histogram(states[d,*],binsize=step,min=xrange[0],max=xrange[1]),$
           color=1,linestyle=d*2

     plots,[0.2,0.25],[ly,ly],color=2,linestyle=d*2,/normal
     xyouts,0.30,ly,"Prior "+string(d+1,format='(I1)'),$
            color=ntable-1,charsize=2.0,charthick=1.0,/normal
     ly -= 0.05
     plots,[0.2,0.25],[ly,ly],color=1,linestyle=d*2,/normal
     xyouts,0.30,ly,"Posterior "+string(d+1,format='(I1)'),$
            color=ntable-1,charsize=2.0,charthick=1.0,/normal
     ly -= 0.05
     
  ENDFOR
;;  plots,obsdat,1,color=3,psym=4,symsize=3
  oplot,xval,histogram(obsdat,binsize=step,min=xrange[0],max=xrange[1])*nrens/nrobs,$
        color=3,linestyle=0
  plots,[0.2,0.25],[ly,ly],color=3,linestyle=0,/normal
  xyouts,0.30,ly,"OBS", color=ntable-1,charsize=2.0,charthick=1.0,/normal
  ly -= 0.05

  oplot,xval,histogram(HA_save[*],binsize=step,min=xrange[0],max=xrange[1]) / nrobs,$
        color=4,linestyle=0
  plots,[0.2,0.25],[ly,ly],color=4,linestyle=0,/normal
  xyouts,0.30,ly,"HA", color=ntable-1,charsize=2.0,charthick=1.0,/normal
  ly -= 0.05


  print,status

  stop
  
END
