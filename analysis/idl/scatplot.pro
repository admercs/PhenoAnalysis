PRO scatplot,x,y,plotname,p
COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, range, xrange, yrange, length, nyears, legenddx, legenddy
   print,'scatter-plotting: ',plotname


   binning=0   ; do not plot all scatter points but bin them
   nb = 50 ; n=100 bins
   
   dims=size(x,/dimensions)
   IF size(dims,/dimensions) GT 1 THEN BEGIN
       nt=dims[1]
       nx=dims[0]
       meanx=fltarr(nt)
       stddevx=fltarr(nt)
       minx=fltarr(nt)
       maxx=fltarr(nt)
       IF nx EQ 2 THEN BEGIN
          meanx = (reform(x[1,*]) + reform(x[0,*]))*0.5
          minx=reform(x[1,*])
          maxx=reform(x[0,*])
       ENDIF ELSE BEGIN
           FOR t=0L,nt-1 DO BEGIN
               meanx[t]=mean(x[*,t],/nan)
               ;stddevx[t]=meanabsdev(x[*,t],/nan)
               good=where(finite(x[*,t]),count)
               IF count GE 2 THEN stddevx[t]=stddev(x[*,t],/nan) ELSE stddevx[t]=!values.f_nan
           ENDFOR
           minx=meanx-stddevx
           maxx=meanx+stddevx
       ENDELSE
   ENDIF ELSE BEGIN
       nt=dims[0]
       nx=1
       meanx=x
   ENDELSE

 

   dims=size(y,/dimensions)
   IF size(dims,/dimensions) GT 1 THEN BEGIN
       nt=dims[1]
       ny=dims[0]
       meany=fltarr(nt)
       stddevy=fltarr(nt)
       miny=fltarr(nt)
       maxy=fltarr(nt)
       IF ny EQ 2 THEN BEGIN
          meany = (reform(y[1,*]) + reform(y[0,*]))*0.5
          miny=reform(y[1,*])
          maxy=reform(y[0,*])
       ENDIF ELSE BEGIN
           FOR t=0L,nt-1 DO BEGIN
               meany[t]=mean(y[*,t],/nan)
               ;stddevy[t]=meanabsdev(y[*,t],/nan)
               good=where(finite(y[*,t]),count)
               IF count GE 2 THEN stddevy[t]=stddev(y[*,t],/nan) ELSE stddevy[t]=!values.f_nan
           ENDFOR
           miny=meany-stddevy
           maxy=meany+stddevy
       ENDELSE
   ENDIF ELSE BEGIN
       nt=dims[0]
       ny=1
       meany=y
   ENDELSE

   print,'nx: ',nx
   print,'ny: ',ny

   xsize=float(xrange[1]-xrange[0])
   ysize=float(yrange[1]-yrange[0])

   IF colored THEN BEGIN
      ;; for rainbow: 
      ;;color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>1)
      ;; for colorbrewer 9 color qualitative cb
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


   ;; plot x variable and errors
   IF nx GT 1L THEN BEGIN
      IF binning THEN BEGIN
         xmin=min(meanx,/nan)
         xmax=max(meanx,/nan)
         dx=(xmax-xmin)/float(nb)
         meanyb=fltarr(nb)
         minxb=fltarr(nb)
         maxxb=fltarr(nb)
         
         FOR t=0L,nb-1L DO BEGIN
            a=where((meanx GE xmin+dx*t) AND (meanx LT xmin+dx*(t+1L)),count)
            IF count GT 0 THEN BEGIN
               meanyb[t]=mean(meany[a],/nan)
               minxb[t]=mean(minx[a],/nan)
               maxxb[t]=mean(maxx[a],/nan)


               plots,[minxb[t],maxxb[t]],[meanyb[t],meanyb[t]],color=color[p], linestyle=linestyle[p], thick=1
               plots,[minxb[t],minxb[t]],[meanyb[t]-0.01*ysize,meanyb[t]+0.01*ysize],color=color[p], linestyle=linestyle[p], thick=1
               plots,[maxxb[t],maxxb[t]],[meanyb[t]-0.01*ysize,meanyb[t]+0.01*ysize],color=color[p], linestyle=linestyle[p], thick=1
            ENDIF
         ENDFOR
      ENDIF ELSE BEGIN
         FOR t=0L,nt-1L DO BEGIN
            plots,[minx[t],maxx[t]],[meany[t],meany[t]],color=color[p], linestyle=linestyle[p], thick=1
            plots,[minx[t],minx[t]],[meany[t]-0.01*ysize,meany[t]+0.01*ysize],color=color[p], linestyle=linestyle[p], thick=1
            plots,[maxx[t],maxx[t]],[meany[t]-0.01*ysize,meany[t]+0.01*ysize],color=color[p], linestyle=linestyle[p], thick=1
            
         ENDFOR
      ENDELSE
      
   ENDIF

   ;; plot y variable and errors
   IF ny GT 1L THEN BEGIN
      
      IF binning THEN BEGIN
         xmin=min(meanx,/nan)
         xmax=max(meanx,/nan)
         dx=(xmax-xmin)/float(nb)
         meanxb=fltarr(nb)
         minyb=fltarr(nb)
         maxyb=fltarr(nb)
         
         FOR t=0L,nb-1L DO BEGIN
            a=where((meanx GE xmin+dx*t) AND (meanx LT xmin+dx*(t+1L)),count)
            IF count GT 0 THEN BEGIN
               meanxb[t]=mean(meanx[a],/nan)
               minyb[t]=mean(miny[a],/nan)
               maxyb[t]=mean(maxy[a],/nan)
            
               plots,[meanxb[t],meanxb[t]],[minyb[t],maxyb[t]],color=color[p], linestyle=linestyle[p], thick=1
               plots,[meanxb[t]-0.01*xsize,meanxb[t]+0.01*xsize],[minyb[t],minyb[t]],color=color[p], linestyle=linestyle[p], thick=1
               plots,[meanxb[t]-0.01*xsize,meanxb[t]+0.01*xsize],[maxyb[t],maxyb[t]],color=color[p], linestyle=linestyle[p], thick=1
            ENDIF
         ENDFOR
      ENDIF ELSE BEGIN
         FOR t=0L,nt-1L DO BEGIN
            plots,[meanx[t],meanx[t]],[miny[t],maxy[t]],color=color[p], linestyle=linestyle[p], thick=1
            plots,[meanx[t]-0.01*xsize,meanx[t]+0.01*xsize],[miny[t],miny[t]],color=color[p], linestyle=linestyle[p], thick=1
            plots,[meanx[t]-0.01*xsize,meanx[t]+0.01*xsize],[maxy[t],maxy[t]],color=color[p], linestyle=linestyle[p], thick=1
            
         ENDFOR
      ENDELSE
      
   ENDIF

   IF (nx EQ 1) AND (ny EQ 1) THEN BEGIN

       IF binning THEN BEGIN
           xmin=min(x,/nan)
           xmax=max(x,/nan)
           dx=(xmax-xmin)/float(nb)
           xb=fltarr(nb)
           yb=fltarr(nb)

           FOR i=0L,nb-1L DO BEGIN
               a=where((x GE xmin+dx*i) AND (x LT xmin+dx*(i+1)),count)
               IF count GT 0 THEN BEGIN
                   xb[i]=mean(x[a],/nan)
                   yb[i]=mean(y[a],/nan)
               ENDIF
           ENDFOR

           oplot,xb,yb,thick=thick,color=color[p],linestyle=linestyle[p], $
             psym=symbols[p],symsize=1.0
       ENDIF ELSE BEGIN
           oplot,x,y,thick=thick,color=color[p],linestyle=linestyle[p], $
             psym=symbols[p],symsize=1.0
       ENDELSE
   ENDIF

   print,'scatter-plotting legend: ',plotname
   
   IF colored THEN BEGIN
      ;; for rainbow: 
      ;;color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>1)
      ;; for colorbrewer 9 color qualitative cb
      color=[ntable-1,indgen(9)+1]
      linestyle=replicate(0,10)
      symbols=indgen(10)+1
      symbols[2]=7
   ENDIF ELSE BEGIN
      color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>50)
      linestyle=indgen(10)
      symbols=indgen(10)+1
      symbols[2]=7
   ENDELSE
   

   IF regline AND (nx EQ 1) AND (ny EQ 1) THEN BEGIN
      index=where(finite(x) AND finite(y),count)
      IF count GT 0 THEN BEGIN
                
         ; ordinary least squares
         a=regress(x[index],y[index],const=b,mcorrelation=r, $
                   measure_errors=replicate(10.,n_elements(index)),/double)

         rms=sqrt(total((y[index] - x[index])^2.0)/float(count))
         
;         plots,[max([range[0],(range[0]-b)/a]),min([range[1],(range[1]-b)/a])], $
;            [max([range[0],(range[0]-b)/a])*a+b,min([range[1],(range[1]-b)/a])*a+b], $
;            color=color[p],linestyle=linestyle[p], thick=thick, /data
         
;         plotname=plotname+': R!E2!N='+strtrim(string(r,format='(f8.2)'),2)+ $
;            ', a='+strtrim(string(a,format='(f8.2)'),2)
         
         plots,[range[0],range[1]], [range[0],range[1]], $
            color=ntable-1,linestyle=linestyle[0], thick=thick, /data
         
         plotname=plotname+': r='+strtrim(string(r,format='(f8.2)'),2)+' rmse='+strtrim(string(rms,format='(f8.2)'),2)
         
         print,"R: ",r
         print,"RMSE: ",rms
         print,"scale:  ",a
         print,"offset:  ",b
      ENDIF
   ENDIF
   
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
