;+
; NAME:
;	ll2rpol
;
; PURPOSE:
;	convert lat/lon to physical grid in
;	rotated pole coord.
;
;	Adaption of the DWD-Functions to convert real geographical 
;	coordinates (LAT,LON) into coordinates in the roteted system 
;	(RLAT,RLON). All angles are in degrees (north>0, east>0)
;   
; CATEGORY:
;
; CALLING SEQUENCE:
;	LL2PH,LON,LAT,RLON,RLAT,pollon=pollon,pollat=pollat
;
; EXAMPLE:
;
; INPUTS:
;	lon,lat:	geogr. coord.
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;	pollon,pollat:	coord. of the rotated pole
;			(default is pole for SM model)
;
; OUTPUTS:
;	rlon,rlat:	physical (rot.) coord.
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	
; MODIFICATION HISTORY:
;     05/90   D.MAJEWSKI (DWD), G. DE MORSIER (SMA) 
;     03/93   D.BRESCH (ETHZ)
; David N. Bresch, 970701
; Reto Stockli 2009/04/25 (Blue Marble Research), vectorized
;-
pro LL2RPOL,LON,LAT,RLON,RLAT,pollon=pollon,pollat=pollat

;; define pole-position:
if n_elements(pollat) eq 0 then pollat=45.d0
if n_elements(pollon) eq 0 then pollon=-175.d0

zrpi18=360.d0/(2.d0*!dpi)
zpir18=2d0*!dpi/360.d0
zsinpol = sin(zpir18*pollat)
zcospol = cos(zpir18*pollat)
zlonpol = zpir18*pollon

;; do case without rotated pole
if (abs(pollat-90.d0) lt 1.d-3) then begin
    rlat=lat
    rlon=lon
    return
endif

;; first, the conversion of LAT to RLAT:
ZLAT    = ZPIR18*LAT
ZLON    = LON
idx = where(ZLON GT 180.d0,count)
if count gt 0 then ZLON[idx] = ZLON[idx] - 360.d0
ZLON    = ZPIR18*ZLON
ZARG    = ZCOSPOL*COS(ZLAT)*COS(ZLON-ZLONPOL) + ZSINPOL*SIN(ZLAT)
RLAT = ZRPI18*ASIN(ZARG)
 
;; now, the conversion for RLON follws: 
ZLAT    = ZPIR18*LAT
ZLON    = LON
idx = where(ZLON GT 180.d0,count)
if count gt 0 then ZLON[idx] = ZLON[idx] - 360.d0
ZLON    = ZPIR18*ZLON
ZARG1   = - SIN(ZLON-ZLONPOL)*COS(ZLAT)
ZARG2   = - ZSINPOL*COS(ZLAT)*COS(ZLON-ZLONPOL)+ZCOSPOL*SIN(ZLAT)

RLON = ZRPI18*ATAN(ZARG1,ZARG2)

END
