PRO taylorplot,x,y,plotname,p, varname
   COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline,range, xrange, yrange, length, nyears, legenddx, legenddy
   
   printstat = 1 ;; print statistics
   plotvarname = 1  ;; plot variable name

   dims=size(y,/dimensions)
   nt=dims[0]

; Make a vector of 16 points, A[i] = 2pi/16: 
   A = FINDGEN(17) * (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
;;   USERSYM, COS(A), SIN(A), /FILL 
   USERSYM, COS(A), SIN(A)

   IF colored THEN BEGIN
      ;;color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>1)
      color=[ntable-1,indgen(8)+1,indgen(8)+1,indgen(8)+1,indgen(8)+1]
      linestyle=replicate(0,30)
;      symbols=(indgen(10)+1)<7
;      symbols[2]=7
      symbols = replicate(5,30)
   ENDIF ELSE BEGIN
      color=fix([ntable-1,[reverse(findgen((nplot-1)>1)/((nplot-2)>1))*(ntable-2)]]>50)
      linestyle=indgen(10)
      symbols=(indgen(10)+1)<7
      symbols[2]=7
   ENDELSE


   index=where(finite(x) AND finite(y),count)

   IF count GT 0 THEN BEGIN

      ;; taylor variables
       stdx=stddev(x[index])
       stdy=stddev(y[index])
       varx=stdx^2
       vary=stdy^2
       stdnorm=stdy/stdx
       rval = (correlate(x[index],y[index]))>0.0
       theta=acos(rval)

;       sig_diff_mean = (TM_TEST(x[index],y[index],/unequal))[1] LT 0.05
;       sig_diff_var = (FV_TEST(x[index],y[index]))[1] LT 0.05

       ;; new Taylor plot
;;       varx=stdx^2
;;       vary=stdy^2
;;       normbias=(vary-varx)/(range[1]-range[0])/0.5 + 1.0
;;       normbias=(stdy-stdx)/(stdy+stdx) + 1.0
;;       normbias=bias/(range[1]-range[0])/0.5 + 1.0
;;       stdnorm=normbias
       
       IF printstat THEN BEGIN
          meany=mean(y[index])
          meanx=mean(x[index])
          bias=meany-meanx
          mab=mean(abs(y[index]-x[index]))
          rms=sqrt(total((y[index] - x[index])^2 )/float(count-1) )       
          cp_rms=sqrt( total(((y[index]-meany) - (x[index]-meanx))^2 )/float(count) )       
;;          cp_rms2=sqrt(stdy^2 + stdx^2 - 2.*stdx*stdy*cos(theta))


          sig_diff_bias = abs(bias) GT 0.05*(range[1]-range[0])
          sig_diff_rms = rms GT 0.05*(range[1]-range[0])
          sig_diff_corr = 1 - corr_ttest(corr=rval,N=count,sl=0.05)
          
          rms=sqrt(total((y[index] - x[index])^2.0)/float(count))
          rms2=sqrt(bias^2 + cp_rms^2)
          
          print,plotname, ': ',varname
          print,'N: ',count,' sd/sd: ',stdnorm,' r: ',rval,' Bias: ',bias,' cpRMS: ',cp_rms, ' RMS: ',rms
          print,'MAB: ',mab
          print,'R: ',cos(theta),' RMS: ',rms,format='(A,f5.2,A,f5.1)'
          print,'stdx: ',stdx,' stdy: ',stdy
          print,'varx: ',varx,' vary: ',vary
          print,'meanx: ',meanx,' meany: ',meany
          print,'significant bias: ',sig_diff_bias
          print,'significant rms:  ',sig_diff_rms
          print,'significant corr: ',sig_diff_corr
          print,' '
          
       ENDIF

;       IF (stdnorm LE 2.0) AND (stdnorm GE 0.0) THEN BEGIN
          stdnorm = (stdnorm > 0.05) < 1.95

          IF sig_diff_bias OR sig_diff_rms OR sig_diff_corr THEN BEGIN
             symsize = 2.5
             symthick = 2.0
             symbols = replicate(4,30)
          ENDIF ELSE BEGIN 
             symsize = 2.5 
             symthick = 3.0
          ENDELSE
           oplot,[stdnorm],[theta],/polar, thick=symthick,color=color[p],$
             linestyle=linestyle[p],psym=symbols[p],symsize=symsize,noclip=1

           IF plotvarname THEN BEGIN
              xc=cos(theta)*stdnorm
              yc=sin(theta)*stdnorm
              xyouts,xc*1.03,yc*1.03,varname,charsize=2.0*charsize, $
                     charthick=1.5*charthick,color=color[p],alignment=0.0
           ENDIF
;       ENDIF

       IF legend THEN BEGIN
           x0=0.65+legenddx
           y0=0.95+legenddy
           dy=0.05
           plots,[x0],[y0-p*dy],color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal
           plots,[x0],[y0-p*dy],color=color[p], $
             linestyle=linestyle[p], thick=thick,/normal,psym=symbols[p],symsize=2
           xyouts,x0+0.03,y0-p*dy-0.01,plotname,charsize=2.0*charsize, $
             charthick=charthick,/normal,color=ntable-1
       ENDIF

   ENDIF

   RETURN
END
