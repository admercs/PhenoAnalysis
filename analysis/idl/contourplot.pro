PRO contourplot,x,y,plotname,p,dt
   COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, range, xrange, yrange, length, nyears, legenddx, legenddy
   print,'plotting: ',plotname
   
   IF colored THEN BEGIN
      color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>1)
      linestyle=replicate(0,10)
   ENDIF ELSE BEGIN
      color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>50)
      linestyle=indgen(10)
   ENDELSE
   
   badvalues=where((y LT range[0]),count)
   IF count GT 0 THEN y[badvalues]=!values.f_nan
   
   ntime=round(86400./dt)
   IF ntime GT 3 THEN BEGIN
       nday=(size(y,/dimensions))[0]/ntime
       vartemp=fltarr(ntime,nday)
       vartemp[*,*]=y

       nlevels=20
       contour,transpose(vartemp),(indgen(nday)+0.5)*ntime,indgen(ntime)*24./ntime+0.5*24./ntime,/cell_fill,/overplot, $
         levels=indgen(nlevels)*(range[1]-range[0])/(nlevels-1)+range[0], /follow, $
         c_colors=(indgen(nlevels)*(ntable)/(nlevels-1)-2)>1B

   ENDIF

;   IF legend THEN $
;      contour,transpose(y),indgen(365),indgen(nday)*24./nday,/overplot, $
;      charsize=charsize,charthick=charthick,thick=1.0,$
;      nlevels=1,/follow, levels=[0,100,200]
   
   RETURN
END
