FUNCTION fill_gauss, data, dia, nodata
;; runs a spatial gaussian smoothing on the data to fill small holes
;; TODO: the sigma calculation is wrong. Please update!


;; arguments:
;; data: input dataset (2D)
;; dia: fill diameter (pixels)
;; nodata: nodata value of input dataset

;; output: smoothed version with holes removed

  dims=size(data,/dimensions)
  nx=dims[0]
  ny=dims[1]
  warr=fltarr(nx,ny)
  wtot=fltarr(nx,ny)
  wcnt=intarr(nx,ny)

  ;; derive gaussian sigma parameter from gaussian FWHM (=0.5 times the diameter)
  sigma= dia / 2.35482 

  ;; first put the array on a larger grid with nodata values covering
  ;; the excess radius
  data2 = fltarr(nx+2*dia,ny+2*dia)
  data2[*] = nodata
  data2[dia:dia+nx-1,dia:dia+ny-1] = data

  maxcnt = 0
  FOR x=-dia,dia DO BEGIN
     FOR y=-dia,dia DO BEGIN
        arr=(shift(data2,x,y))[dia:dia+nx-1,dia:dia+ny-1]
        g=exp(-(x^2+y^2)/(2.*sigma^2))
        dist=sqrt(x^2+y^2)>0.5
        circle=dist LE dia
        maxcnt += circle
        warr=warr+float(arr NE nodata)*arr*g*float(circle)
        wtot=wtot+float(arr NE nodata)*g*float(circle)
        wcnt=wcnt+(arr NE nodata)*circle
     ENDFOR
  ENDFOR

  mincnt = round(0.1 * maxcnt)>1

  newdata=float(data NE nodata)*data + $
          float(data EQ nodata)*(warr/(wtot>0.001)*float(wcnt GE mincnt) + $
                                 float(wcnt LT mincnt)*nodata)

  RETURN, newdata
END

