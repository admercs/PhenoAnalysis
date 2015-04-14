PRO overplot,x,y,plotname,p,connect
   COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, $
     range, xrange, yrange, length, nyears, legenddx, legenddy, aspect
 
   ;; polyfill error bars or not
   use_poly = 0
  
   print,'plotting: ',plotname
   dims=size(y,/dimensions)
   IF size(dims,/dimensions) GT 1 THEN BEGIN
       nt=dims[1]
       np=dims[0]
       means=fltarr(nt)
       stddevs=fltarr(nt)
       mins=fltarr(nt)
       maxs=fltarr(nt)
       IF np EQ 2 THEN BEGIN
           mins=reform(y[1,*])
           maxs=reform(y[0,*])
       ENDIF ELSE BEGIN
           FOR t=0L,nt-1 DO BEGIN
               means[t]=mean(y[*,t],/nan)
               ;stddevs[t]=meanabsdev(y[*,t],/nan)
               stddevs[t]=stddev(y[*,t],/nan)
           ENDFOR
           mins=means-stddevs
           maxs=means+stddevs
       ENDELSE
   ENDIF ELSE BEGIN
       nt=dims[0]
       np=1
   ENDELSE

   IF colored THEN BEGIN
      ;; for rainbow: 
      ;;color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>1)
      ;; for colorbrewer 9 color qualitative cb
      color=[ntable-1,indgen(9)+1]
      ;color=[255,254,253,50]
      linestyle=replicate(0,10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDIF ELSE BEGIN
      color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>50)
      linestyle=indgen(10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDELSE

   IF np GT 1 THEN BEGIN
       vi=where(finite(mins) AND finite(maxs),vn)
       dx=0.01*(max(x)-min(x))/float(aspect)
       IF (vn gt 300L*nyears) THEN BEGIN
          FOR l=0L,20L*nyears DO BEGIN
             t=(l*nt/(20L*nyears))<(nt-1L)
             t1=((l+1)*nt/(20L*nyears))<(nt-2L)  ;; -2 because of most years only having 365 valid days ;-)
;             print,t,t1
             IF use_poly THEN BEGIN
                IF finite(mins[t]) AND finite(mins[t1]) AND finite(maxs[t]) AND finite(maxs[t1]) THEN $
                   polyfill,[x[t],x[t1],x[t1],x[t]],[mins[t],mins[t1],maxs[t1],maxs[t]],color=color[p]
             ENDIF ELSE BEGIN
                plots,[x[t],x[t]],[mins[t],maxs[t]],color=color[p], linestyle=linestyle[p], thick=1
                plots,[x[t]-dx,x[t]+dx],[mins[t],mins[t]],color=color[p], linestyle=linestyle[p], thick=1
                plots,[x[t]-dx,x[t]+dx],[maxs[t],maxs[t]],color=color[p], linestyle=linestyle[p], thick=1
             ENDELSE
          ENDFOR
       ENDIF ELSE BEGIN
          FOR t=0L,nt-1L DO BEGIN
             IF use_poly AND (t LT (nt-1L)) THEN BEGIN
                IF finite(mins[t]) AND finite(mins[t+1]) AND finite(maxs[t]) AND finite(maxs[t+1]) THEN $
                   polyfill,[x[t],x[t+1],x[t+1],x[t]],[mins[t],mins[t+1],maxs[t+1],maxs[t]],color=color[p]
             ENDIF ELSE BEGIN
                plots,[x[t],x[t]],[mins[t],maxs[t]],color=color[p], linestyle=linestyle[p], thick=1
                plots,[x[t]-dx,x[t]+dx],[mins[t],mins[t]],color=color[p], linestyle=linestyle[p], thick=1
                plots,[x[t]-dx,x[t]+dx],[maxs[t],maxs[t]],color=color[p], linestyle=linestyle[p], thick=1
             ENDELSE
          ENDFOR
       ENDELSE

   ENDIF ELSE BEGIN

       IF connect THEN oplot,x,y,thick=thick,color=color[p],linestyle=linestyle[p]
       IF connect THEN BEGIN
          vi=where(finite(y),vn)
          IF vn GT 50L*nyears THEN BEGIN
             oplot,congrid(x,20L*nyears),congrid(y,20L*nyears),thick=thick,color=color[p],$
                   linestyle=linestyle[p],psym=symbols[p],symsize=2 
          ENDIF ELSE BEGIN
             oplot,x,y,thick=thick,color=color[p],$
                   linestyle=linestyle[p],psym=symbols[p],symsize=2
          ENDELSE
       ENDIF ELSE BEGIN
           oplot,x,y,thick=thick,color=color[p],$
             linestyle=linestyle[p],psym=symbols[p],symsize=2
       ENDELSE

       IF legend THEN BEGIN
           x0=0.20/aspect+legenddx
           y0=0.85+legenddy
           dx=0.1/aspect
           dy=0.05
           plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal
           plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal,psym=symbols[p],symsize=2
           xyouts,x0+dx+0.002,y0-p*dy-0.01,plotname,charsize=2.0*charsize, $
             charthick=charthick,/normal,color=ntable-1
       ENDIF
   ENDELSE

   RETURN
END
