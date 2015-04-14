PRO colbar_standalone, DIVISIONS=divisions, NCOLORS=ncolors, range=range, horizontal=horizontal, $
                       UNITS=units, varname=varname, CRANGE=crange, ctnum=ctnum, ctrev=ctrev, $
                       FILENAME=filename, psout=psout, monthly=monthly, doy=doy
                       
  tvbar = 0 ;; use TV instead of Polyfill

  monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

  IF psout AND NOT keyword_set(filename) THEN BEGIN
     print,'Please supply a filename for the colorbar'
     stop
  ENDIF

  IF keyword_set(doy) AND ~keyword_set(horizontal) THEN fac=0.25 ELSE fac=0.15 

  IF psout THEN BEGIN
     ;; setting up plotting to PS device
     set_plot,'ps'

     size_default=7.5        ; default height/width of the plot [cm]

     IF keyword_set(horizontal) THEN BEGIN
        ysize = size_default * fac  
        xsize = size_default
     ENDIF ELSE BEGIN
        ysize = size_default     
        xsize = size_default * fac
     ENDELSE

     thick=size_default/8.0*3.0 ;; plot line thickness
     charthick=2.0              ;; charthickness relative to 1.0
     charsize=size_default/11.0 ;; charsize relative to 10pt.

     print,'Writing EPS: ',filename
     DEVICE, /encapsul, bits_per_pixel=8,/color, font_size=5, $
             Filename=filename, preview=0,/isolatin1, $
             xsize=xsize, ysize=ysize,/HELVETICA

     !P.FONT=0

  ENDIF ELSE BEGIN

     ;; set up plotting window in X
     set_plot,'x'
     
     size_default=400.0        ; default height of the plot [pixels]

     IF keyword_set(horizontal) THEN BEGIN
        ysize = round(size_default * fac * 1.25)
        xsize = round(size_default)
     ENDIF ELSE BEGIN
        ysize = round(size_default)  
        xsize = round(size_default * fac * 1.25)
     ENDELSE
     
     thick=2.0                                ;; plot line thickness
     charthick=1.5                            ;; charthickness relative to 1.0
     charsize=1.0 * sqrt(size_default/400.0)  ;; charsize relative to 10pt.
     window,30,xsize=xsize,ysize=ysize
 
     device,decomposed=0,set_font='Helvetica',/tt_font

     !P.FONT=1
    
  ENDELSE

  loadct,ctnum,file='brewer.tbl' 

  tvlct,red,green,blue,/get
  ntable=!d.table_size
  
  IF keyword_set(ctrev) THEN BEGIN
     red = reverse(red)
     green = reverse(green)
     blue = reverse(blue)
  ENDIF

  red[ntable-2]=20B
  green[ntable-2]=20B
  blue[ntable-2]=150B
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
  
  IF keyword_set(monthly) THEN BEGIN
     divisions = 12
     ncolors = 12
     range = [0,12]
  ENDIF
  IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 12
  IF N_ELEMENTS(ncolors) EQ 0 THEN ncolors = divisions
  IF N_ELEMENTS(varname) EQ 0 THEN varname = 'test'
  IF N_ELEMENTS(units) EQ 0 THEN units = '-'
  IF N_ELEMENTS(crange) EQ 0 THEN crange=[1,ntable-4]
  IF N_ELEMENTS(range) EQ 0 THEN range=[0,1]
  
  ;; generate colors for color bar
  colors=byte(((findgen(ncolors) + 0.5)/float(ncolors) * float(crange[1]-crange[0]))+crange[0])

  title=varname+' ('+units+')'
  
  format = '(f6.1)'
  ticklen = 0.0

  IF keyword_set(horizontal) THEN BEGIN
     barx0 = 0.05
     bary0 = 0.75
     bardx = 0.9
     bardy = 0.2
     yticks=1
     xticks=divisions
     yrange=[0,1]
     xrange=range
  ENDIF ELSE BEGIN
     barx0 = 0.05
     bary0 = 0.05
     bardx = 0.2
     bardy = 0.9
     yticks=divisions
     xticks=1
     yrange=range
     xrange=[0,1]
  ENDELSE

  ;; build plot area
  PLOT, xrange, yrange, /NODATA, XTICKS=xticks, YTICKS=yticks, XTHICK=2, YTHICK=2, $
        XSTYLE=5, YSTYLE=5, charthick=charthick, XMINOR=1, YMINOR=1, $
        COLOR=ntable-1, CHARSIZE=2.0*charsize, /NOERASE,/normal, $
        YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=ticklen , XTICKLEN=ticklen, $
        POSITION=[barx0,bary0,barx0+bardx,bary0+bardy]

  ;; build and draw color bar
  IF tvbar THEN BEGIN
     IF keyword_set(horizontal) THEN BEGIN
        bar = colors # REPLICATE(1B,20)
     ENDIF ELSE BEGIN
        bar = REPLICATE(1B,20) # colors
     ENDELSE
        
     IF ~keyword_set(psout) THEN bar = CONGRID(bar, xsize*bardx, ysize*bardy)
     TV, bar, barx0,bary0,xsize=bardx,ysize=bardy,/normal      
  ENDIF ELSE BEGIN
     FOR c=0,ncolors-1 DO BEGIN
        IF keyword_set(horizontal) THEN BEGIN
           dx = bardx / ncolors
           dy = bardy
           x0=barx0+float(c)*dx
           x1=barx0+float(c+1)*dx
           y0=bary0
           y1=bary0+bardy
        ENDIF ELSE BEGIN
           dx = bardx 
           dy = bardy / ncolors
           x0=barx0
           x1=barx0+bardx
           y0=bary0+float(c)*dy
           y1=bary0+float(c+1)*dy
        ENDELSE
        polyfill,[x0,x0,x1,x1],[y0,y1,y1,y0],/normal,color=colors[c]
    ENDFOR
  ENDELSE

  ;; plot axes
  IF keyword_set(monthly) THEN BEGIN
     ;; month names
     IF keyword_set(horizontal) THEN BEGIN
        AXIS, XAXIS=0, XRANGE=xrange, XTICKS=xticks, XTHICK=2, YTHICK=2, $
              XTICKLEN=ticklen, XSTYLE=0, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
              XTITLE=title, XMINOR=1, charthick=charthick, $
              XTICKNAME = ['      ','J','F','M','A','M','J','J','A','S','O','N','D','       '], $
              XTICKV = findgen(14)-0.5  
     ENDIF ELSE BEGIN
        AXIS, YAXIS=1, YRANGE=yrange, YTICKS=yticks, XTHICK=2, YTHICK=2, $
              YTICKLEN=ticklen, YSTYLE=0, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
              YTITLE=title, YMINOR=1, charthick=charthick, $
              YTICKNAME = ['      ','J','F','M','A','M','J','J','A','S','O','N','D','       '], $
              YTICKV = findgen(14)-0.5  
     ENDELSE

  ENDIF ELSE BEGIN
     IF keyword_set(doy) THEN BEGIN
        ;; day and month names
        time0 = (timegen(1,start=julday(1,1,2001)+range[0]-1))[0]
        time1 = (timegen(1,start=julday(1,1,2001)+range[1]-1))[0]
        time = time0 + findgen(divisions+1)*(time1-time0)/float(divisions)
        caldat,time,month,day
        IF keyword_set(horizontal) THEN BEGIN
           ticknames = string(day,format='(I2)')+'!C'+monthnames[month-1]+' '
           AXIS, XAXIS=0, XRANGE=xrange, XTICKS=xticks, XTHICK=2, YTHICK=2, $
                 XTICKLEN=ticklen, XSTYLE=0, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
                 XTITLE='!C'+title, XMINOR=1, charthick=charthick, XTICKNAME = ticknames  
         ENDIF ELSE BEGIN
            ticknames = string(day,format='(I2)')+' '+monthnames[month-1]+' '
            AXIS, YAXIS=1, YRANGE=yrange, YTICKS=yticks, XTHICK=2, YTHICK=2, $
                  YTICKLEN=ticklen, YSTYLE=0, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
                  YTITLE=title, YMINOR=1, charthick=charthick, YTICKNAME = ticknames  
        ENDELSE

     ENDIF ELSE BEGIN
        ;; all other cases
        IF keyword_set(horizontal) THEN BEGIN
           AXIS, XAXIS=0, XRANGE=xrange, XTICKS=xticks, XTHICK=2, YTHICK=2, $
                 XTICKLEN=ticklen, XSTYLE=0, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
                 XTITLE=title, XMINOR=1, charthick=charthick  
        ENDIF ELSE BEGIN
           AXIS, YAXIS=1, YRANGE=yrange, YTICKS=yticks, XTHICK=2, YTHICK=2, $
                 YTICKLEN=ticklen, YSTYLE=0, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
                 YTITLE=title, YMINOR=1, charthick=charthick  
        ENDELSE
     ENDELSE
  ENDELSE

  IF keyword_set(horizontal) THEN BEGIN
     plots,[barx0,barx0,barx0+bardx,barx0+bardx],[bary0,bary0+bardy,bary0+bardy,bary0], $
           thick=2,color=ntable-1,/normal
  ENDIF ELSE BEGIN
     plots,[barx0+bardx,barx0,barx0,barx0+bardx],[bary0,bary0,bary0+bardy,bary0+bardy], $
           thick=2,color=ntable-1,/normal
  ENDELSE

  IF psout THEN BEGIN
     DEVICE,/CLOSE      
     set_plot,'x'
     
     IF keyword_set(horizontal) THEN BEGIN
        spawn,'convert -density 300 -resize x300 -flatten '+filename+' '+filename+'.png'
     ENDIF ELSE BEGIN
        spawn,'convert -density 300 -resize 300 -flatten '+filename+' '+filename+'.png'
     ENDELSE
     spawn,'epstopdf '+filename     
  ENDIF 
 
END
