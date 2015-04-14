;+
; :Description: 
;   performs averaging to output dimensions that are not necessarily
;   integer factors of the input dimensions. This is an extension of
;   the IDL-internal rebin function, but is much slower than rebin.
;
; :Categories:
;   Statistics
;
; :Params:
;   data : in, required, type=fltarr
;     2D data array (input)
;   nx_out : in, required, type=integer 
;     first output dimension length
;   ny_out : in, required, type=integer 
;     second output dimension length
;
; :Keywords:
;  missing : in, optional, type="same type as data" 
;    missing value in data and outdata
;
; :Returns: 
;   rebinned 2D float array of nx_out and ny_out dimensions
;
; :Todo:
;   modify for working with 1D and 3D arrays
;
;   make work for magnifying as well
;
;   make work for non-floating point data as well
;
; :Requires:
;   IDL 7.1
;
; :Examples:
;
; :History:
;   2010/08/03 Reto Stockli
;
;   2011/03/18 : added standard idldoc-compatible rst-style header
;
;   2013/07/15 : improved handing of NAN missing value
;
; :Author:
;   Reto Stockli (Blue Marble Research)
;    
; :Copyright:
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;
;   You should have received a copy of the GNU General Public License
;   along with this program.  If not, see http://www.gnu.org/licenses/.
;
;-
FUNCTION rebin_new, data, nx_out, ny_out, missing=missing

  ;; temporary storage and conversion to float
  indata = float(data)

  dims = size(indata,/dimensions)

  IF n_elements(dims) EQ 2 THEN BEGIN
     nx = dims[0]
     ny = dims[1]
  ENDIF ELSE BEGIN
     nx = dims[0]
     ny = 1L
  ENDELSE

  ;; set NaN
  IF keyword_set(missing) THEN BEGIN
     IF finite(missing) THEN BEGIN
        idx = where(indata EQ missing,count)
        IF count GT 0L THEN indata[idx] = !values.f_nan
     ENDIF
  ENDIF

  ;; create output array
  outdata = fltarr(nx_out,ny_out)

  ;; create indices
  x0 = lindgen(nx_out)*nx/nx_out
  x1 = (lindgen(nx_out)+1)*nx/nx_out-1
  y0 = lindgen(ny_out)*ny/ny_out
  y1 = (lindgen(ny_out)+1)*ny/ny_out-1

  ;; rebinning
  FOR y=0,ny_out-1 DO BEGIN
     FOR x=0,nx_out-1 DO BEGIN
        outdata[x,y] = mean(indata[x0[x]:x1[x],y0[y]:y1[y]],/NAN)
     ENDFOR
  ENDFOR

  ;; set missing value for output
  IF keyword_set(missing) THEN BEGIN
     IF finite(missing) THEN BEGIN
        idx = where(finite(outdata,/NAN),count)
        IF count GT 0L THEN outdata[idx] = missing
     ENDIF
  ENDIF

  RETURN, outdata

END
