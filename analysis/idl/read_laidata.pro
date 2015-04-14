PRO read_laidata, obsvarname, dtobs, year,obsdir,sitename, obsvar, obsvarerror, obstitle, obsunits, nodata, sitelon, sitelat, has_obs

syear=''
reads,year,syear
syear=strtrim(syear,2)
 
leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
IF leapyear THEN ndays = 366 ELSE ndays = 365

porosity=55. ; for soil moisture observation plots
;porosity=0.6                     ; for soil moisture observation plots
lv=2.52*1.e6                    ; latent heat of vaporization (J/kg)

has_obs=0
noobs=-9999.
infiles=(file_search(obsdir+'phenology/'+sitename+'/lai.'+syear+'.dat'))[0]

dtobs = 24.*3600.
ntobs = 365
obsvar = replicate(nodata,ntobs)
obsvarerror = replicate(nodata,ntobs)

IF infiles NE '' THEN BEGIN
   has_obs = 1
   openr,1,infiles

   day=0
   lai=0.
   laistddev=0.
   WHILE (EOF(1) EQ 0) DO BEGIN
      ;; implement shift by one day after FEB29 of a leap year!!!
      readf,1,day,lai,laistddev
      obsvar[day-1]=lai
      obsvarerror[day-1]=laistddev
   END
   close,1

   obstitle='Leaf Area Index'
   obsunits='m2/m2'

ENDIF



END
