;+
; NAME:
;	rpol2ll
;
; PURPOSE:
;	convert physical grid to lat/lon in
;	rotated ploe coord.
;
;     Adaption of the DWD-Functions to convert rotated pole coordinates
;     (RLON,RLAT) into geogrphic coordinates (LON,LAT). The location of 
;     the rotated pole is passed trough the common block  /rotpol/. The 
;     first two arguments are input, the second two are output. All angles 
;     are in degrees (north>0, east>0)
;   
; CATEGORY:
;
; CALLING SEQUENCE:
;       rpol2ll,rlon,rlat,lon,lat,pollon=pollon,pollat=pollat
;
; EXAMPLE:
;
; INPUTS:
;	rlon,rlat:	physical (rot.) coord.
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;	pollon,pollat:	coord. of the rotated pole
;			(default is pole for SM model)
;
; OUTPUTS:
;	lon,lat:	geogr. coord.
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
;     05/90   D.MAJEWSKI (DWD)
;     03/93   D.BRESCH (ETHZ)
; David N. Bresch, 970701
; Reto Stockli, 2009/04/24 (Blue Marble Research), vectorized
;-
PRO RPOL2LL,RLON,RLAT,LON,LAT,pollon=pollon,pollat=pollat

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
    lat=rlat
    lon=rlon
    return
endif

;; first, the conversion of RLAT to LAT:
ZRLAT  = ZPIR18*RLAT
ZRLON  = RLON

idx = where(ZRLON GT 180.d0,count)
IF count GT 0 THEN ZRLON[idx] = ZRLON[idx] - 360.d0
ZRLON = ZPIR18*ZRLON
ARG   = ZCOSPOL*COS(ZRLAT)*COS(ZRLON) + ZSINPOL*SIN(ZRLAT)
LAT   = ZRPI18*ASIN(ARG)

;; follows conversion of RLON to LON:
ZRLAT = ZPIR18*RLAT
ZRLON = RLON
idx = where(ZRLON GT 180.d0,count)
IF count GT 0 THEN ZRLON[idx] = ZRLON[idx] - 360.d0

ZRLON = ZPIR18*ZRLON

ZARG1 = SIN(ZLONPOL)*(-ZSINPOL*COS(ZRLON)*COS(ZRLAT) + $
                      ZCOSPOL*SIN(ZRLAT)) - COS(ZLONPOL)*SIN(ZRLON)*COS(ZRLAT)
ZARG2   = COS(ZLONPOL)*(-ZSINPOL*COS(ZRLON)*COS(ZRLAT) + $
                        ZCOSPOL*SIN(ZRLAT)) + SIN(ZLONPOL)*SIN(ZRLON)*COS(ZRLAT)

LON = ZRPI18*ATAN(ZARG1,ZARG2)

END
