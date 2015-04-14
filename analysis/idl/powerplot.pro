PRO powerplot,dt,x,y,plotname,p
   COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, range, xrange, yrange, length, nyears, legenddx, legenddy
   print,'plotting: ',plotname
   
   IF colored THEN BEGIN
      ;;color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>1)
      color=[ntable-1,indgen(9)+1]
      linestyle=replicate(0,10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDIF ELSE BEGIN
      color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>50)
      linestyle=indgen(10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDELSE
;   IF (colored) AND (p NE 0) THEN color=p*(ntable-2)/nplot ELSE linestyle=p
;   IF (p NE 0) THEN color=p*(ntable-2)/nplot
   
   nt=n_elements(x)

   n21=nt/2.+1
   freq=x
   freq[n21]=n21-nt+findgen(n21-2)
   freq=(freq/(2.*365.*(dt/3600.)))[0:n21-1] ; frequency in day-1
   f=(abs(fft(y,-1)))[0:n21-1]

   oplot,congrid(freq,20),congrid(f,20),thick=thick,color=color[p],$
     linestyle=linestyle[p],psym=symbols[p],symsize=2
   oplot,freq,f,thick=thick,color=color[p],linestyle=linestyle[p]

   IF legend THEN BEGIN
       x0=0.20+legenddx
       y0=0.85+legenddy
       dx=0.1
       dy=0.05
       plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
         linestyle=linestyle[p], thick=thick,/normal
       plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
         linestyle=linestyle[p], thick=thick,/normal,psym=symbols[p],symsize=2
       xyouts,x0+dx+0.01,y0-p*dy-0.01,plotname,charsize=2.0*charsize, $
         charthick=charthick,/normal,color=ntable-1
   ENDIF
      
   RETURN
END
