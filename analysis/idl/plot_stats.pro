PRO plot_stats

  title = 'Assimilation Performance'

  lai_mad  = [2.2860933,0.91795359,0.46013505,0.40656871,0.33916844]
  fpar_mad = [0.3097218,0.14264538,0.07383467,0.06965413,0.06492696]

  nobs = [0,13943482,50768668,213642410,869605738]

  lai_var = [2.7932509,0.44688009,0.17191457,0.13850928,0.11279419]
  fpar_var = [0.32614888,0.057623868,0.026195720,0.014948076,0.011999570]

  xtitle = 'Assimilation Experiment'
  xrange = [-0.1,4.1]
  xtickvalue = [0,1,2,3,4]
  xtickname = ['PRIOR','4','16','64','256']

  nobsrange = [min(nobs)-0.05*(max(nobs)-min(nobs)),max(nobs)+0.05*(max(nobs)-min(nobs))]
  nobsrange = [-0.05*(max(nobs)-min(nobs)),max(nobs)+0.05*(max(nobs)-min(nobs))]
  nobstitle = 'Observation Count'


  laititle = 'LAI Error & Uncertainty'

  lairange = [-0.05*(max([lai_mad,lai_var])-min([lai_mad,lai_var])),$
              max([lai_mad,lai_var])+0.1*(max([lai_mad,lai_var])-min([lai_mad,lai_var]))]

  fpartitle = 'FPAR Error & Uncertainty'

  fparrange = [-0.05*(max([fpar_mad,fpar_var])-min([fpar_mad,fpar_var])),$
               max([fpar_mad,fpar_var])+0.05*(max([fpar_mad,fpar_var])-min([fpar_mad,fpar_var]))]

  
  filename = '../plots/plots/stats/assimilation_stats.eps'
  ysize_ps_default=7.5          ; default height of the plot [cm]
  aspect = 1.25

  set_plot,'ps' 

  ysize = ysize_ps_default

  thick=ysize/8.0*5.0           ; plot line thickness
  charthick=2.0                 ; charthickness relative to 1.0
  charsize=ysize/11.0           ; charsize relative to 10pt.

  !P.FONT=0
  print,'Writing EPS: ',filename
  DEVICE, /encapsul, bits_per_pixel=8,/color, font_size=5, $
          Filename=filename, preview=0,/isolatin1, $
          xsize=ysize*aspect, ysize=ysize,/HELVETICA

  loadct,50,file='brewer.tbl' 
  tvlct,red,green,blue,/get
  green=shift(green,1)
  blue=shift(blue,1)
  red=shift(red,1)

  red[0]=255B
  green[0]=255B
  blue[0]=255B
  ntable=!d.table_size
  red[ntable-1]=0B
  green[ntable-1]=0B
  blue[ntable-1]=0B
  tvlct,red,green,blue

  plot,[0],[0],yrange=nobsrange,xrange=xrange,xstyle=9,ystyle=9, $
       xtickv=xtickvalue,xtickname=xtickname,color=ntable-1, /nodata,$
       xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
       xthick=thick, ythick=thick, xtitle=xtitle, title=title, $
       /normal, ytitle=nobstitle, xticklen = -0.01*aspect, yticklen=-0.01, $
       position = [0.15,0.10,0.75,0.92]
                    

  xindex = indgen(5)

  oplot,xindex,nobs,thick=thick,color=ntable-1

  axis,0.8,yaxis=1,yrange=lairange,ystyle=1, color=1, /save, $
       yminor=1,charsize=2.0*charsize,charthick=charthick, $
       ythick=thick,/normal,ytitle=laititle,yticklen=-0.01,/data

  oplot,xindex,lai_mad,thick=thick,color=1
  oplot,xindex,lai_var,thick=thick,color=1,linestyle=2

  axis,0.9,yaxis=1,yrange=fparrange,ystyle=1, color=2, /save, $
       yminor=1,charsize=2.0*charsize,charthick=charthick, $
       ythick=thick,/normal,ytitle=fpartitle,yticklen=-0.01,/data

  oplot,xindex,fpar_mad,thick=thick,color=2
  oplot,xindex,fpar_var,thick=thick,color=2,linestyle=2

  DEVICE,/CLOSE      
  set_plot,'x'
  
  pixelsize=1024
  spixelsize=''
  reads,pixelsize*aspect,spixelsize
  spawn,'convert -density 300 -flatten -resize '+spixelsize+' '+filename+' '+filename+'.png'
  spawn,'epstopdf '+filename
  

  stop

END
