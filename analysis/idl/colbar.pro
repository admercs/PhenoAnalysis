PRO colbar, DIVISIONS=divisions, NCOLORS=ncolors, $
            UNITS=units, varname=varname, CRANGE=crange, monthly=monthly
  COMMON share
  
  IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 12
  IF N_ELEMENTS(ncolors) EQ 0 THEN ncolors = 120
  IF N_ELEMENTS(varname) EQ 0 THEN varname = 'test'
  IF N_ELEMENTS(units) EQ 0 THEN units = '-'
  IF N_ELEMENTS(crange) EQ 0 THEN crange=[1,252]
  
  title=varname+' ('+units+')'
  
  format = '(f6.1)'
  ticklen = 0.5
  barx0 = 0.01
  bary0 = 0.05
  bardx = 0.03
  bardy = 0.9

  bar = REPLICATE(1B,20) # (BINDGEN(crange[1]-crange[0]+1)+crange[0])
  bar = BYTE(congrid(bar,20,ncolors))
  IF psout EQ 0 THEN bar = CONGRID(bar, xsize*bardx, ysize*bardy)
  
  TV, bar, barx0,bary0,xsize=bardx,ysize=bardy,/normal      
  
  PLOT, range, crange, /NODATA, XTICKS=1, XTHICK=2, YTHICK=2, $
        YTICKS=divisions, XSTYLE=1, YSTYLE=9, charthick=charthick, $
        COLOR=ntable-1, CHARSIZE=2.0*charsize, /NOERASE,/normal, $
        YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
        YMINOR=2, POSITION=[barx0,bary0,barx0+bardx,bary0+bardy],/DEVICE
     
  IF keyword_set(monthly) THEN BEGIN
     ;; month names
     AXIS, YAXIS=1, YRANGE=range, YTICKS=divisions, XTHICK=2, YTHICK=2, $
           YTICKLEN=ticklen, YSTYLE=1, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
           YMINOR=2, charthick=charthick, YTICKNAME = strarr(13) + '      '
     
     AXIS, YAXIS=1, YRANGE=range, YTICKS=divisions, XTHICK=2, YTHICK=2, $
           YTICKLEN=0, YSTYLE=1, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
           YTITLE=title, YMINOR=2, charthick=charthick, $
           YTICKNAME = ['      ','J','F','M','A','M','J','J','A','S','O','N','D','       '], $
           YTICKV = findgen(14)-0.5   
  ENDIF ELSE BEGIN
     ;; regular color bar
     AXIS, YAXIS=1, YRANGE=range, YTICKS=divisions, XTHICK=2, YTHICK=2, $
           YTICKLEN=ticklen, YSTYLE=1, COLOR=ntable-1, CHARSIZE=2.0*charsize, $
           YTITLE=title, YMINOR=2, charthick=charthick
  ENDELSE
END
