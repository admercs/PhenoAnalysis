FUNCTION fill_dist, data, rad, nodata
;; inverse distance weighted spatial smoothing of dataset
;; which removes holes no larger than the circular area with radius rad

;; 2009/06/02 Reto Stockli (Blue Marble Research)
;; 2009/11/05 added minimum # of valid pixels dependent on radius (10%
;; valid pixels needed)

;; arguments:
;; data: input dataset (2D)
;; rad: fill radius (pixels) ---> this is actually the fill rad
;; nodata: nodata value of input dataset

;; output: smoothed version with holes removed

  dims=size(data,/dimensions)
  nx=dims[0]
  ny=dims[1]
  warr=fltarr(nx,ny)
  wtot=fltarr(nx,ny)
  wcnt=intarr(nx,ny)

  ;; first put the array on a larger grid with nodata values covering
  ;; the excess radius
  data2 = fltarr(nx+2*rad,ny+2*rad)
  data2[*] = nodata
  data2[rad:rad+nx-1,rad:rad+ny-1] = data

  maxcnt = 0
  FOR x=-rad,rad DO BEGIN
     FOR y=-rad,rad DO BEGIN
        arr=(shift(data2,x,y))[rad:rad+nx-1,rad:rad+ny-1]
        dist=sqrt(x^2+y^2)>0.5
        circle = dist LE rad
        maxcnt += circle
        warr=warr+float(arr NE nodata)*arr/dist*float(circle)
        wtot=wtot+float(arr NE nodata)/dist*float(circle)
        wcnt=wcnt+(arr NE nodata)*circle
     ENDFOR
  ENDFOR

  mincnt = round(0.1 * maxcnt)>1

  newdata=float(data NE nodata)*data + $
           float(data EQ nodata)*(warr/(wtot>0.001)*float(wcnt GE mincnt) +  $
                                  float(wcnt LT mincnt)*nodata)

  RETURN, newdata
END
