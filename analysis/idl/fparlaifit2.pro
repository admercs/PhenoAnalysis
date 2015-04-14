PRO fparlaifit2,x,a,f,pder

;; This function fits a double logarithmic function between fpar and
;; lai time series and returns three coefficients that determine this
;; function

;; To be used with IDL's CURVEFIT procedure

;; 2009/05/31 Reto Stockli (Blue Marble Research)

; f(x) = ln(1-x/a)/ln(1-b) * c * a
; df/da = 1/(1-x/a) *(-x/a^2) * c * a / ln(1-b) + ln(1-x/a)/ln(1-b) * c
; df/db = ln(1-x/a)/(ln(1-b))^2 *(1./(1-b)) * (-1) * c * a
; df/dc = ln(1-x/a)/ln(1-b) * a

; a[0] = vcover
; a[1] = fparsat
; a[2] = laimax

f = alog( (1.-x/(a[0]>0.001))>0.001 ) / alog( (1.- a[1])>0.001 ) * a[2] * a[0]

pder = [[1./(1.-x/(a[0]>0.001))*(-x/((a[0]^2)>0.001))*a[2]*a[0]/alog((1.-a[1])>0.001) + $
         alog((1.-x/(a[0]>0.001))>0.001) / alog(( 1.- a[1])>0.001 ) * a[2] ], $
        [alog((1.-x/(a[0]>0.001))>0.001 ) / (alog(( 1.- a[1])>0.001 ))^2.0 * a[2] * a[0] / (( 1.- a[1])>0.001 ) * (-1.)], $
        [alog((1.-x/(a[0]>0.001))>0.001 ) / alog(( 1.- a[1])>0.001 ) * a[0]] ]
END
