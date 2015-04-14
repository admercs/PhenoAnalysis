PRO histoplot,x,y,plotname,p
   COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, $
     range, xrange, yrange, length, nyears, legenddx, legenddy
   
   nh=100 ; number of histogram bins
   barplot=0 ; plot histogram as bars

   print,'plotting: ',plotname
   dims=size(y,/dimensions)
   IF size(dims,/dimensions) GT 1 THEN BEGIN
       nt=dims[1]
       np=dims[0]

       nplottemp=nplot*2
       BINSIZE = (MAX(y,/nan) - MIN(y,/nan)) / (nh - 1)
       goodt=where(finite(y[0,*]),goodtcount)

       IF goodtcount GT 2 THEN BEGIN
           startt=goodt[0]
           endt=goodt[goodtcount-1]
           
           histystart=histogram(y[*,startt],max=max(y,/nan),min=min(y,/nan),nbins=nh,locations=indexystart,/nan)
           histyend=histogram(y[*,endt],max=max(y,/nan),min=min(y,/nan),nbins=nh,locations=indexyend,/nan)
           
       ENDIF ELSE BEGIN
           histystart=!values.f_nan
           histyend=!values.f_nan
           indexystart=!values.f_nan
           indexyend=!values.f_nan
       ENDELSE
   ENDIF ELSE BEGIN
       nt=dims[0]
       np=1

       nplottemp=nplot

       BINSIZE = (MAX(y,/nan) - MIN(y,/nan)) / (nh - 1)

       histy=histogram(y,max=max(y,/nan),min=min(y,/nan),nbins=nh,locations=indexy,/nan)
   ENDELSE

   IF colored THEN BEGIN
      ;;color=fix([ntable-1,[reverse(findgen((nplottemp-1)>1)/((nplottemp-2)>1))*(ntable-2)]]>1)
      ;; for colorbrewer 9 color qualitative cb
      color=[ntable-1,indgen(9)+1]
      linestyle=replicate(0,10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDIF ELSE BEGIN
      color=fix([ntable-1,[reverse(findgen((nplottemp-1)>1)/((nplottemp-2)>1))*(ntable-2)]]>50)
      linestyle=indgen(10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDELSE
;   IF (colored) AND (p NE 0) THEN color=p*(ntable-2)/nplottemp ELSE linestyle=p
;   IF (p NE 0) THEN color=p*(ntable-2)/nplottemp
   

   IF np EQ 1 THEN BEGIN

       IF barplot THEN BEGIN
           FOR t=0,nh-1 DO BEGIN
               x0=indexy[t]
               x1=indexy[t]+binsize
               y0=0
               y1=float(histy[t])/max(histy,/nan)
               polyfill,[x0,x0,x1,x1],[y0,y1,y1,y0],color=color[p]
           ENDFOR
       ENDIF ELSE BEGIN
           oplot,indexy,float(histy)/max(histy,/nan),thick=thick,color=color[p],$
             linestyle=linestyle[p] ;,psym=symbols[p],symsize=2
           oplot,indexy,float(histy)/max(histy,/nan),thick=thick,color=color[p],$
             linestyle=linestyle[p],psym=symbols[p],symsize=2
       ENDELSE
   ENDIF ELSE BEGIN
       IF barplot THEN BEGIN
           FOR t=0,nh-1 DO BEGIN
               x0=indexystart[t]
               x1=indexystart[t]+binsize
               y0=0
               y1=float(histystart[t])/max([histystart,histyend],/nan)
               polyfill,[x0,x0,x1,x1],[y0,y1,y1,y0],color=color[p]
           ENDFOR
           FOR t=0,nh-1 DO BEGIN
               x0=indexyend[t]
               x1=indexyend[t]+binsize
               y0=0
               y1=float(histyend[t])/max([histystart,histyend],/nan)
               polyfill,[x0,x0,x1,x1],[y0,y1,y1,y0],color=color[p+nplot]
           ENDFOR
       ENDIF ELSE BEGIN
           oplot,indexystart,float(histystart)/max([histystart,histyend],/nan),thick=thick,color=color[p],$
             linestyle=linestyle[p] ;,psym=symbols[p],symsize=2
           oplot,indexyend,float(histyend)/max([histystart,histyend],/nan),thick=thick,color=color[p+nplot],$
             linestyle=linestyle[p+nplot] ;,psym=symbols[p],symsize=2
           oplot,indexystart,float(histystart)/max([histystart,histyend],/nan),thick=thick,color=color[p],$
             linestyle=linestyle[p],psym=symbols[p],symsize=2
           oplot,indexyend,float(histyend)/max([histystart,histyend],/nan),thick=thick,color=color[p+nplot],$
             linestyle=linestyle[p+nplot],psym=symbols[p+nplot],symsize=2
       ENDELSE
   ENDELSE

   IF legend THEN BEGIN
       x0=0.20+legenddx
       y0=0.85+legenddy
       dx=0.1
       dy=0.05

       IF np EQ 1 THEN BEGIN
           plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal
           plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal,psym=symbols[p],symsize=2
           xyouts,x0+dx+0.01,y0-p*dy-0.01,plotname,charsize=2.0*charsize, $
             charthick=charthick,/normal,color=ntable-1
       ENDIF ELSE BEGIN
           plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal
           plots,x0+findgen(3)/3.0*dx,replicate(y0-p*dy,4),color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal,psym=symbols[p],symsize=2
           xyouts,x0+dx+0.01,y0-p*dy-0.01,plotname+': Start',charsize=2.0*charsize, $
             charthick=charthick,/normal,color=ntable-1

           plots,x0+findgen(3)/3.0*dx,replicate(y0-(p+nplot)*dy,4),color=color[p+nplot], $
             linestyle=linestyle[p+nplot], thick=thick,/normal
           plots,x0+findgen(3)/3.0*dx,replicate(y0-(p+nplot)*dy,4),color=color[p+nplot], $
             linestyle=linestyle[p+nplot], thick=thick,/normal,psym=symbols[p+nplot],symsize=2
           xyouts,x0+dx+0.01,y0-(p+nplot)*dy-0.01,plotname+': End',charsize=2.0*charsize, $
             charthick=charthick,/normal,color=ntable-1
       ENDELSE
   ENDIF

   RETURN
END
