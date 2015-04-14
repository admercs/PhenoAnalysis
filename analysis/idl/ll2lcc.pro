PRO ll2lcc, lon, lat, x, y, lambda0, phi0, phi1, phi2, $
            dx=dx, dy=dy, x0=x0, y0=y0, a=a, b=b, f=f

;; projects latitude/longitude into a lambert conic conformal
;; projection

;; Input Arguments
;; lon : longitude [degrees]
;; lat : latitude [degrees]
;; lambda0 : reference longitude [degrees]
;; phi0 : reference latitude [degrees]
;; phi1 : standard parallel 1 [degrees]
;; phi2 : standard parallel 2 [degrees]

;; Optional Input Arguments
;; dx : x grid spacing [m]
;; dy : y grid spacing [m]
;; x0 : x grid offset [m]
;; y0 : y grid offset [m]
;; a  : ellipse major axis [m]
;; b  : ellipse minor axis [m]
;; f  : ellipse flattening [-] 

;; Output Arguments
;; x : lcc x coordinate [m or grid points]
;; y : lcc y coordinate [m or grid points]

  ;; No elliptical spheriod allowed in this code yet
  ;; NARR sphere definition (a=b, f=0)
  IF keyword_set
  a=6371200.0d0
  b=6371200.0d0
  f=0.0d0
  

  ;; corner points:
;; (1,1)   1.000N, 145.500W
;; (1,277) 46.635N, 148.639E
;; (349,277)       46.352N, 2.566W
;; (349,1) 0.897N, 68.318W


;; projection constants
  radeg=360.d0/(2.d0*!dpi)
  dtor=2d0*!dpi/360.d0

  ;; secand cone (1st and 2nd standard parallel different)
  secant = abs(phi1-phi2) GE 1.d-10

  rphi1 = phi1 * dtor
  rphi2 = phi2 * dtor
  rphi0 = phi0 * dtor
  rlambda0 = lambda0 * dtor

  IF (secant) THEN BEGIN
     n = ( alog(cos(rphi1)) - alog(cos(rphi2)) ) / $
         ( alog(tan(0.25d0*!dpi + 0.5d0*rphi2)) - alog(tan(0.25d0*!dpi + 0.5d0*rphi1)) )
  ENDIF ELSE BEGIN
     n = sin(rphi1)
  ENDELSE

  F = cos(rphi1) * tan(0.25d0*!dpi + 0.5d0*rphi1)^n / n
  rrho0 = F * tan(0.25d0*!dpi + 0.5d0*rphi0)^(-n)

  IF forward THEN BEGIN
     ;; calculate x/y from phi/lambda: forward projection

     rphi = phi*dtor
     rlambda = lambda*dtor

     rdlambda = rlambda - rlambda0

     ;; check if dlambda is outside of +/- PI
     IF (abs(rdlambda) GT !dpi) THEN BEGIN
        rdlambda += !dpi  ;; adjust to 0..2pi rad
        rdlambda -= 2.d0*!dpi * float(floor(rdlambda / (2.d0* !dpi))) ;; remove integral # of revolutions
        rdlambda -= !dpi ;; adjust back to -pi..pi rad
     ENDIF

     rrho = F * tan(0.25d0*!dpi + 0.5d0*rphi)^(-n)

     x = rrho * sin(n*rdlambda)
     y = rrho0 - rrho*cos(n*rdlambda)
     
     x = a*x/dx - x0
     y = a*y/dy - y0

     print,x,y

END
