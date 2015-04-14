COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, range, $
   xrange, yrange, length, nyears, legenddx, legenddy, aspect, pollon, pollat, psout, $
   ysize_ps_default, ysize_x_default, xsize, ysize

  COMMON MODISTILE,vmodis,hmodis,xmodis,ymodis,xidmodis,yidmodis,modis_lon,modis_lat, $
     npmodis,hmin,hmax,vmin,vmax
  
  DEVICE,retain=2


;FOR sites = 1,30 DO BEGIN

;; TODO: re-implement land-screening for MODIS read!

;; Detect FPE and print code lines where FPE occurs during runtime (res)
  !EXCEPT = 2

;; define nodata value
  nodata = -9999.               ; missing data value used throughout the code (if not !values.f_nan)

;; Directories

;  openu,35,'test.dat'

;; the pheno analysis output directory (where the output of the
;; experiments are stored          print*,"SAVE: ",p

  basedir='/Users/stockli/phenoanalysis-output/'
;  basedir='/glade/p/cgd/tss/people/stockli/phenoanalysis-output/'

;; MODIS and pft data is stored here
;  datadir='/Users/stockli/phenoanalysis-datasets'
  datadir='/Users/stockli'
;  datadir='/glade/p/cgd/tss/people/stockli'

;; ground-measured LAI time-series *not supported right now
  laidir='' 

;; file for site/regional information
  sitedir='../../data/sites/' 
  sitefile = 'example.dat'
;  sitefile = 'global4.dat'
;  sitefile = 'garonna.dat'
;  sitefile = 'swiss.dat'
;  sitefile = 'harvard.dat'

;; experiment information
  sites=[4]

  experiments = ['Example1']
;  experiments = ['Example3']
;  experiments = ['Prediction']
;  experiments = ['Assimilation']

;; name of variables to plot
  varnames=['LAI']
;  varnames=['FPAR']

;varnames=['GSI']
;varnames=['LAI','TEMP_FAC','MOIST_FAC','LIGHT_FAC'] ;,'FORCE_FAC','CHILL_FAC']
varnames=['FPAR','TEMP_FAC','MOIST_FAC','LIGHT_FAC'] ;,'FORCE_FAC','CHILL_FAC']
;varnames=['TMIN_AVE','TMIN_MAX','TMIN_MIN']
;varnames=['VPD_AVE','VPD_MAX','VPD_MIN']
;varnames=['PHOTO_AVE','PHOTO_MAX','PHOTO_MIN']
;varnames=['RG_AVE','RG_MAX','RG_MIN']
;varnames=['RAIN_AVE','RAIN_MAX','RAIN_MIN']
;varnames=['CHILL_AVE','CHILL_MAX','CHILL_MIN']
;varnames=['LAI','LAI_MIN','LAI_MAX','LAI_SAT']
;varnames=['FPAR','FPAR_SAT']
;varnames=['TMIN','TMEAN']
;varnames=['GROWTHRATE','DECAYRATE']
;varnames=['GROWTHRATE']
;varnames=['DECAYRATE']
;varnames=['LAI_MIN']
;varnames=['LAI_MAX']
;varnames=['TMIN']
;varnames=['TMIN_MIN']
;varnames=['TMIN_MAX']
;varnames=['PHOTO']
;varnames=['PHOTO_MIN']
;varnames=['PHOTO_MAX']
;varnames=['VPD']
;varnames=['VPD_MIN']
;varnames=['VPD_MAX']
;varnames=['FVCOVER']
;varnames=['FPAR_SAT']
;varnames=['LAI_SAT']
;varnames=['TMIN_TAVE']
;varnames=['PHOTO_TAVE']
;varnames=['VPD_TAVE']
;varnames=['TMIN_TAVE','PHOTO_TAVE','VPD_TAVE','CHILL_TAVE']
;varnames=['TMIN_AVE']
;varnames=['PHOTO_AVE']
;varnames=['VPD_AVE']
;varnames=['MOIST_FAC']
;varnames=['LIGHT_FAC']
;varnames=['TEMP_FAC']

  obsvarnames=varnames
  testvarname=''                ; varname of the variable to plot against with scatter plots (shown on x-axis)

;; range of primary axis (if not autorange is chosen)
  range=[0.0,7.5]
;  range=[0.0,1.0]
;  range=[0,300]

;; range of second axis (if multiple variables with different units
;; are plotted and autorange is not chosen)
  range2=[0,1]

;; range of axis to plot agains with scatter plots (if not autorange
;; is chosen)
  testrange=[0,1]

;; range of observation axis in scatter plots (if not autorange is chosen)
  obsrange=[0.,1.2]

; User definable plot arguments
  years=[2002]

;  years=[2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010]
;years=indgen(51)+1958

  modis_version = 5

;; rotated pole (set to 180./90. for default)
;pollon = -170.0
;pollat = 43.0
  pollon = 180.0
  pollat = 90.0

;; flags
  colored=1                     ; produce color/monochrome output
  psout=0                       ; write post script files 
  integrate=0                   ; temporally integrate data on plots
  taylor=0                      ; make a taylor plot
  histo=0                       ; make a histogram plot of the data
  scatter=0                     ; make scatter plots
  scatterall=0                  ; plot individual observations on scatter plot
  contour=0                     ; plot contour fish plots
  power=0                       ; power-spectrum plot
  regline=1                     ; plot regression line or stddev shade in diurplot
  legend=1                      ; plot a legend
  autorange=1                   ; automatically guess plot range
  plotobs=1                     ; plot observations also (mandatory with scatter plot)
  saveobs=0                     ; read saved observations instead of MODIS tiles
  simulator=0                   ; plot simulated observations from data assimilation
  plotmodels=1                  ; plot observations also (mandatory with scatter plot)
  dtplot=24.*3600.              ; seconds per time step

  plotexperimentname=1  ; plot experimentname in legend
  plotmodelname=0       ; plot model names in legend
  plotvarname=1         ; plot variable names in legend
  ignoregaps=1          ; do not sync obs and model gaps
  selpft = 1            ; selected pft class in case pft's are not averaged in model output
  selhgt = 1            ; select hgt class in case hgt's are not averaged in model output
  landscreen = 0        ; screen observations by land cover class of the site (1) or not (0) 
  errscale = 0.5        ; for plotting purposes: scale the QA-based Error bars (arbitrary units) by this factor

; final data output options
  interannual=0
  printascii=0     
  printmonthlylai=0
  printlai=0

; Design
;;  aspectset = 1                 ; force aspect ratio
  ysize_ps_default=7.5          ; default height of the plot [cm]
  ysize_x_default=450           ; default height of the plot [pixels]

  legenddx=0.02
  legenddy=-0.0

  IF simulator THEN BEGIN
     plotobs = 1
  ENDIF

;; plots go here
  plotdir=basedir+experiments[0]+'/plots/'

;; Time Attributes
  monstr_short=['J','F','M','A','M','J','J','A','S','O','N','D']
  monstr_long=['JAN','FEB','MAR','APR','MAY','JUN', $
               'JUL','AUG','SEP','OCT','NOV','DEC']
  monday_noleap=[31,28,31,30,31,30,31,31,30,31,30,31]
  monday_leap=[31,29,31,30,31,30,31,31,30,31,30,31]
  
  monsum_noleap=[0,31,59,90,120,151,181,212,243,273,304,334,365]
  monsum_leap=[0,31,60,91,121,152,182,213,244,274,305,335,366]
  
  
; Other arrays
  nyears=long(n_elements(years))

  IF nyears GT 1 THEN aspectset = 3 ELSE aspectset = 1
  
;; calculate number of days in time axis
  ndays = 0L
  FOR y=0,nyears-1 DO BEGIN
     year = years[y]
     leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
     
     IF leapyear THEN ndays += 366L ELSE ndays += 365L
     
  ENDFOR
  
  models = 1  ;; we only have pheno_analysis as the model right now
  
   nsites=n_elements(sites)
   nobs=n_elements(obsvarnames)
   nexperiments=n_elements(experiments)
   nvarnames=n_elements(varnames)
   nobsvarnames=n_elements(obsvarnames)
   nmodels=n_elements(models)

   IF taylor THEN BEGIN
       ntaylor=nexperiments  ; Taylor plot legend will only contain number of models

       modelstemp=models
       experimentstemp=experiments
       varnamestemp=varnames
       sitestemp=sites

       ; varnames x models
       varnames=replicate(varnamestemp[0],nmodels)
       FOR v=1,nvarnames-1 DO BEGIN
           experiments=[experiments,experimentstemp]
           models=[models,modelstemp]
           varnames=[varnames,replicate(varnamestemp[v],nmodels)]
       ENDFOR

       modelstemp=models
       experimentstemp=experiments
       varnamestemp=varnames
       sitestemp=sites

       ; varnames x models x sites
       sites=replicate(sitestemp[0],nmodels*nvarnames)
       FOR s=1,nsites-1 DO BEGIN
           experiments=[experiments,experimentstemp]
           models=[models,modelstemp]
           varnames=[varnames,varnamestemp]
           sites=[sites,replicate(sitestemp[s],nmodels*nvarnames)]
       ENDFOR
       
       nsites=n_elements(sites)
       nvarnames=n_elements(varnames)
       nexperiments=n_elements(experiments)
       nmodels=n_elements(models)
       nobs=n_elements(obsvarnames)
       nplot=n_elements(varnames)+n_elements(obsvarnames)
   ENDIF ELSE BEGIN

      IF (nvarnames GT 1) THEN BEGIN
         
         IF (nsites NE nvarnames) AND plotmodels THEN BEGIN
            nsites=nvarnames
            sites=replicate(sites[0],nsites)
         ENDIF
         
         IF nexperiments NE nvarnames THEN BEGIN 
            nexperiments=nvarnames
            experiments=replicate(experiments[0],nexperiments)
         ENDIF
         
         IF (nmodels NE nvarnames) AND plotmodels THEN BEGIN
            nmodels=nvarnames
            models=replicate(models[0],nmodels)
         ENDIF
         
      ENDIF ELSE BEGIN
         
         IF nexperiments NE nsites THEN BEGIN
            nexperiments = nsites
            experiments=replicate(experiments[0],nexperiments)
         ENDIF

         IF nmodels NE nexperiments THEN BEGIN 
            nmodels=nexperiments
            models=replicate(models[0],nmodels)
         ENDIF
                  
         IF (nvarnames NE nmodels) AND plotmodels THEN BEGIN
            nvarnames=nmodels
            varnames=replicate(varnames[0],nmodels)
         ENDIF

         IF (nsites NE nvarnames) AND plotmodels THEN BEGIN
            nsites=nvarnames
            sites=replicate(sites[0],nsites)
         ENDIF
         
      ENDELSE
      
      IF plotobs THEN nobs=n_elements(obsvarnames) ELSE nobs=1
      
      nplot=nobs+nmodels
      
   ENDELSE

   ssites=strarr(nsites)
   reads,sites,ssites
   ssites=strtrim(ssites,2)       

   IF nsites GT 1 THEN BEGIN
       test=where(sites EQ sites[0],testcount)
       IF testcount NE nsites THEN manysites=1 ELSE manysites=0
   ENDIF ELSE manysites=0

   IF nmodels GT 1 THEN BEGIN
       test=where(models EQ models[0],testcount)
       IF testcount NE nmodels THEN manymodels=1 ELSE manymodels=0
   ENDIF ELSE manymodels=0

   IF varnames[0] EQ '' THEN varnames=obsvarnames
   varname=varnames[0]          ; varname for the observation data

   ; initial time set-up
   dtobs=!values.f_nan
   dtmod=!values.f_nan
   ntplot=ndays*24L*3600L/long(dtplot) ;; number of plotting time steps
     
   ; Data arrays
   modtitle=strarr(nmodels)
   modelshort=strarr(nmodels)
   modelnames=strarr(nmodels)
   modvartitle=strarr(nmodels)
   modvarunits=strarr(nmodels)
   

   npobs=1   ; number of observation points per time step
   obsvartitle=strarr(nobs)
   obsvarunits=strarr(nobs)
   testvartitle=strarr(nobs)
   testvarunits=strarr(nobs)
   obstitle=strarr(nobs)
   
   ; STRINGS
   syear0=''
   syear1=''
   reads,years[0],syear0
   reads,years[nyears-1],syear1
   syear0=strtrim(syear0,2)
   syear1=strtrim(syear1,2)

   ; test flags
   IF (integrate+scatter+power+contour+taylor+histo) GT 1 THEN BEGIN
       print,'Please only specify one plot method!'
       stop
   ENDIF

   IF (nyears GT 1) AND contour THEN BEGIN
       print,'No multiple years plot possible in contour mode'
       stop
   ENDIF
   

   ; open Site file
   stemp=''
   siteshort=''
   sitename=''
   sitelat=0.
   sitelon=0.

   nsitemax = 100
   siteshorts=strarr(nsitemax)
   sitenames=strarr(nsitemax)
   sitelonmins=fltarr(nsitemax)
   sitelonmaxs=fltarr(nsitemax)
   sitelatmins=fltarr(nsitemax)
   sitelatmaxs=fltarr(nsitemax)
   sitedlons=fltarr(nsitemax)
   sitedlats=fltarr(nsitemax)

   openr,1,sitedir+sitefile
   readf,1,stemp
   readf,1,stemp
   WHILE (NOT EOF(1)) DO BEGIN
       readf,1,stemp
       arr=strsplit(stemp,',',/extract)
       s = arr[0]-1
       siteshorts[s]=strtrim(arr[1],2)
       sitenames[s]=strtrim(arr[2],2)
       sitelonmins[s]=float(arr[3])
       sitelonmaxs[s]=float(arr[4])
       sitelatmins[s]=float(arr[5])
       sitelatmaxs[s]=float(arr[6])
       IF strtrim(arr[7],2) EQ "NA" THEN sitedlons[s] = 0.05 ELSE sitedlons[s]=float(arr[7])
       IF strtrim(arr[8],2) EQ "NA" THEN sitedlats[s] = 0.05 ELSE sitedlats[s]=float(arr[8])
   ENDWHILE
   close,1

   sitename=sitenames[sites[0]-1]

   ;; geographic bounds for plotting and observation array       
   lonmin = sitelonmins[sites[0]-1] 
   lonmax = sitelonmaxs[sites[0]-1]
   latmin = sitelatmins[sites[0]-1]
   latmax = sitelatmaxs[sites[0]-1]
   dlon = sitedlons[sites[0]-1]
   dlat = sitedlats[sites[0]-1]

   IF lonmin EQ lonmax THEN BEGIN
      dlon = 0.05
      lonmin = lonmin - 0.5*dlon
      lonmax = lonmax + 0.5*dlon
   ENDIF

   IF latmin EQ latmax THEN BEGIN
      dlat = 0.05
      latmin = latmin - 0.5*dlat
      latmax = latmax + 0.5*dlat
   ENDIF

   nlon = round((lonmax-lonmin)/dlon)
   nlat = round((latmax-latmin)/dlat)

   IF (nlon GT 1) OR (nlat GT 1) THEN stop,"Only single grid points allowed"

   IF nyears GT 1 THEN stime=syear0+'-'+syear1 ELSE stime=syear0

   CASE 1 OF
      contour: filename=plotdir+'/'+sitename+'.contour.'+stime+'.'+varname
      scatter: IF testvarname EQ '' THEN $
        filename=plotdir+'/'+sitename+'.scatter.'+stime+'.'+varname ELSE $
        filename=plotdir+'/'+sitename+'.scatter.'+stime+'.'+varname+'_'+testvarname
      power: filename=plotdir+'/'+sitename+'.power.'+stime+'.'+varname
      taylor: filename=plotdir+'/'+sitename+'.taylor.'+stime+'.'+varname
      integrate: filename=plotdir+'/'+sitename+'.integrate.'+stime+'.'+varname
      histo: filename=plotdir+'/'+sitename+'.histogram.'+stime+'.'+varname
      ELSE: filename=plotdir+'/'+sitename+'.dayplot.'+stime+'.'+varname
   ENDCASE
   IF plotexperimentname THEN filename=filename+'.'+experiments[0]
   filename=filename+'.eps'

   monstr=monstr_short

   IF ignoregaps THEN gaps=0 ELSE gaps=1
  
   firstobs=1
   firstmodel=1

   has_obs = bytarr(nobs)
   has_model = bytarr(nmodels)

   FOR y=0L,nyears-1 DO BEGIN
       year=years[y]

       leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))       
       IF leapyear THEN ndayofyear = 366L ELSE ndayofyear = 365L

       ;; read OBSERVATIONS
       IF plotobs THEN BEGIN
           obsvartitletemp=''
           obsvarunitstemp=''
           FOR o=0,nobs-1 DO BEGIN
              obsvarname=obsvarnames[o]

              ;; regular MODIS tiles
              IF (simulator EQ 0) AND (saveobs EQ 0) AND ((obsvarname EQ 'LAI') OR (obsvarname EQ 'FPAR')) THEN BEGIN 
                 
                 ;; initalize MODIS grid
                 read_modistile, datadir, 0, 0, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                                 obsvarnames[o], obsvarunitstemp, obsvarlongnamestemp, modis_version, $
                                 obsvartemp, obsvarerrortemp, obsmasktemp, $
                                 pollon=pollon, pollat=pollat,/initialize
                 
                 npobs=npmodis  ; number of observations for whole grid
                 
                 FOR day=1,ndayofyear DO BEGIN
                    
                    ;; read MODIS data
                    read_modistile, datadir, year, day, lonmin, lonmax, latmin, latmax, dlon, dlat, $
                                    obsvarnames[o], obsvarunitstemp, obsvarlongnamestemp, modis_version, $
                                    obsvartemp, obsvarerrortemp, obsmasktemp, $
                                    pollon=pollon, pollat=pollat
                    
                    obstitle[o]='MODIS'
                    IF plotvarname THEN BEGIN
                       obstitle[o]=obstitle[o]+' '+obsvarname
                    ENDIF
                    
                    IF firstobs THEN BEGIN
                       IF finite(dtobs,/nan) THEN dtobs=24.*3600.
                       ntobs=ndays*24L*3600L/long(dtobs)                       
                       obsvar=make_array(npobs,ntobs,nobs,/float,value=!values.f_nan)
                       obsvarerror=make_array(npobs,ntobs,nobs,/float,value=!values.f_nan)
                       obsmask=make_array(npobs,ntobs,nobs,/byte)
                       firstobs=0
                    ENDIF
                    
                    badindex=where(obsvartemp EQ nodata,badcount)
                    IF badcount GT 0 THEN obsvartemp[badindex]=!values.f_nan
                    
                    obsdays0 = 0L
                    FOR y0=0,y-1 DO BEGIN
                       year0 = years[y0]
                       leapyear = ((((year0 mod 4) eq 0) and ((year0 mod 100) ne 0)) or ((year0 mod 400) eq 0))
                       
                       IF leapyear THEN obsdays0 += 366L ELSE obsdays0 += 365L
                    ENDFOR
                    obsindex=(obsdays0 + day - 1L)*24L*3600L/long(dtobs)
                    
                    obsvar[*,obsindex,o]=obsvartemp
                    obsvarerror[*,obsindex,o]=obsvarerrortemp
                    obsmask[*,obsindex,o]=obsmasktemp
                    
                    obsvartitle[o]=obsvartitletemp
                    obsvarunits[o]=obsvarunitstemp
                    
                 ENDFOR
              ENDIF ;; regular MODIS HDF tiles
              
              ;; simulated FPAR&LAI from DA model
              IF simulator OR saveobs THEN BEGIN

                  obsvartitletemp=''
                  obsvarunitstemp=''
                  read_pheno_local,basedir,2+saveobs,experiments[o],sitenames[sites[o]-1],year,obsvarname,dtobs, $
                    obsvartemp,obsvarerrortemp,obsvartitletemp,obsvarunitstemp, selpft, selhgt, $
                    obsshorttemp,obsnamestemp, nodata, obsok
                  
                  npobs = 1

                  has_obs[o] = obsok
                  
; LAI time series from sfc observations
;                    read_laidata, obsvarname, dtobs, year,laidir,sitenames[sites[o]-1], obsvartemp, $
;                                  obsvarerrortemp, obsvartitletemp, obsvarunitstemp, $
;                                  nodata, sitelons[sites[o]-1], sitelats[sites[o]-1], has_obs[o]
                    

                  IF has_obs[o] THEN BEGIN
                      IF (nobs GT 1) THEN $
                        obstitle[o]='OBS: '+obsvarname $
                      ELSE $
                        obstitle[o]='OBS'
                      
                      IF finite(dtobs) THEN BEGIN
                          IF firstobs THEN BEGIN
                              ntobs=ndays*24L*3600L/long(dtobs)                       
                              obsvar=make_array(npobs,ntobs,nobs,/float,value=!VALUES.F_NAN)
                              obsvarerror=make_array(npobs,ntobs,nobs,/float,value=!values.f_nan)
                              obsmask=make_array(npobs,ntobs,nobs,/byte)
                              firstobs=0
                          ENDIF
                          
                          obsdays0 = 0L
                          FOR y0=0,y-1 DO BEGIN
                             year0 = years[y0]
                             leapyear = ((((year0 mod 4) eq 0) and ((year0 mod 100) ne 0)) or ((year0 mod 400) eq 0))
                               
                             IF leapyear THEN obsdays0 += 366L ELSE obsdays0 += 365L
                          ENDFOR
                          obsindex=obsdays0*24L*3600L/long(dtobs) + lindgen(n_elements(obsvartemp))

                          ;; TO BE TESTED
                          obsvar[0,obsindex,o]=obsvartemp 
                          obsvarerror[0,obsindex,o]=obsvarerrortemp
                          obsmask[0,obsindex,o] = finite(obsvartemp) AND finite(obsvarerrortemp)
                          
                          obsvartitle[o]=obsvartitletemp
                          obsvarunits[o]=obsvarunitstemp
                          
                       ENDIF

                  ENDIF
              ENDIF ;; simulator

          ENDFOR
      ENDIF

   ; read MODEL DATA   
       IF plotmodels THEN BEGIN
           FOR m=0,nmodels-1 DO BEGIN
               
               modvartitletemp=''
               modvarunitstemp=''
               read_pheno_local,basedir,models[m],experiments[m],sitenames[sites[m]-1],year,varnames[m],dtmod, $
                 modvartemp,modvarerrortemp,modvartitletemp,modvarunitstemp, selpft, selhgt, $
                 modelshorttemp,modelnamestemp, nodata, modelok

               has_model[m] = modelok

               IF has_model[m] THEN BEGIN

                  IF finite(dtmod) THEN BEGIN
                     IF firstmodel THEN BEGIN
                        ntmod=ndays*24L*3600L/long(dtmod)      
                        modvar=make_array(ntmod,nmodels,/float,value=!VALUES.F_NAN) 
                        modvarerror=make_array(ntmod,nmodels,/float,value=!VALUES.F_NAN)
                        firstmodel=0
                     ENDIF
                     
                     modvartitle[m]=modvartitletemp
                     modvarunits[m]=modvarunitstemp
                     modelshort[m]=modelshorttemp
                     modelnames[m]=modelnamestemp

                     modtitle[m]=''
                     IF plotmodelname THEN BEGIN
                        modtitle[m]=modtitle[m]+' '+modelshort[m]
                     ENDIF
                     IF plotexperimentname THEN BEGIN
                         modtitle[m]=modtitle[m]+' '+experiments[m]
                     ENDIF                     
                     IF plotvarname THEN BEGIN
                        modtitle[m]=modtitle[m]+' '+varnames[m]
                     ENDIF

                     moddays0 = 0L
                     FOR y0=0,y-1 DO BEGIN
                        year0 = years[y0]
                        leapyear = ((((year0 mod 4) eq 0) and ((year0 mod 100) ne 0)) or ((year0 mod 400) eq 0))
                        
                        IF leapyear THEN moddays0 += 366L ELSE moddays0 += 365L
                     ENDFOR

                     dims = size(modvartemp,/dimensions)
                     modindex=moddays0*24L*3600L/long(dtmod) + lindgen(dims[0])
                     modvar[modindex,m]=modvartemp 
                     modvarerror[modindex,m]=modvarerrortemp 
                  ENDIF
               ENDIF ;; if model has data
            ENDFOR ;; for models
        ENDIF ;; plot model output

    ENDFOR ; number of years

   IF firstobs THEN plotobs=0
   IF firstmodel THEN plotmodel=0
   IF plotobs THEN IF finite(dtobs,/nan) OR finite(max(obsvar,/nan),/nan) THEN plotobs=0
   IF plotmodels THEN IF finite(dtmod,/nan) OR finite(max(modvar,/nan),/nan) THEN plotmodels=0

   ; plot data 2D or 1D
   IF plotobs THEN BEGIN
       start1=0L
       end1=ndays*24L*3600L/long(dtobs)
       mtotobs=end1-start1
       index1=lindgen(end1-start1)+start1
   ENDIF
   IF plotmodels THEN BEGIN
       start2=0L
       end2=ndays*24L*3600L/long(dtmod)
       mtotmod=end2-start2
       index2=lindgen(end2-start2)+start2
   ENDIF

   startp=0L
   endp=ndays*24L*3600L/long(dtplot)
   mtot=endp-startp

   axisindex=bytarr(nplot)
   secondaxis=0
   IF taylor EQ 0 THEN BEGIN
      IF plotobs THEN FOR o=0,nobs-1 DO IF (obsvarunits[o] NE obsvarunits[0]) AND (has_obs[o]) THEN axisindex[o]=1
      IF plotmodels AND plotobs THEN FOR m=0,nmodels-1 DO IF modvarunits[m] NE obsvarunits[0] THEN axisindex[m+nobs]=1
      IF plotmodels AND (plotobs EQ 0) THEN FOR m=0,nmodels-1 DO IF modvarunits[m] NE modvarunits[0] THEN axisindex[m+nobs]=1
      IF (max(axisindex) GT 0) THEN secondaxis=1 
   ENDIF

  IF plotobs THEN BEGIN
      FOR o=0,nobs-1 DO BEGIN
          IF axisindex[o] EQ 1 THEN BEGIN
              secondaxistitle=obsvartitle[o]
              secondaxisunits=obsvarunits[o]
          ENDIF
      ENDFOR
  ENDIF
                         
  IF plotmodels THEN BEGIN
      FOR m=0,nmodels-1 DO BEGIN
          IF axisindex[m+nobs] EQ 1 THEN BEGIN
              secondaxistitle=modvartitle[m]
              secondaxisunits=modvarunits[m]
          ENDIF
      ENDFOR
  ENDIF

   ; Check for missing obs and create a mask with NAN
   IF plotmodels THEN BEGIN
      IF manysites THEN badvar=bytarr(ntmod,nobs) ELSE badvar=bytarr(ntmod,1+secondaxis) 
   ENDIF

   IF power THEN BEGIN
       FOR m=0,nmodels-1 DO BEGIN
          badindex=where(finite(modvar[*,m],/nan),count)
          IF count GT 0 THEN BEGIN
             modvar[badindex,m]=mean(modvar[*,m],/nan)
             modvarerror[badindex,m]=mean(modvarerror[*,m],/nan)
          ENDIF
       ENDFOR
   ENDIF

; mask model data with NAN
   IF plotmodels THEN BEGIN
       FOR m=0,nmodels-1 DO BEGIN
          IF gaps THEN BEGIN
             IF manysites THEN badindex=where(badvar[*,m] OR finite(modvar[*,m],/nan),count) ELSE $
                badindex=where(badvar[*,axisindex[nobs+m]] OR finite(modvar[*,m],/nan),count)
          ENDIF ELSE badindex=where(finite(modvar[*,m],/nan),count)         
          IF count GT 0 THEN BEGIN
             modvar[badindex,m]=!VALUES.F_NAN
             modvarerror[badindex,m]=!VALUES.F_NAN
          ENDIF
       ENDFOR
   ENDIF

; rescale data to the plotting time interval (averaging)
   IF plotobs THEN BEGIN
       IF (dtobs NE dtplot) THEN BEGIN
           IF mtotobs GT mtot THEN mtotobs=mtotobs/mtot*mtot
           IF ntobs GT ntplot THEN ntobs=ntobs/ntplot*ntplot
           indexp=rebin(index1[0:mtotobs-1],mtot)*dtobs/dtplot
           
           obsvarnew=make_array(npobs,ntplot,nobs,/float)
           obsmasknew=make_array(npobs,ntplot,nobs,/float)
           obsvarerrornew=make_array(npobs,ntplot,nobs,/float)
           testvarnew=make_array(npobs,ntplot,nobs,/float)
           testmasknew=make_array(npobs,ntplot,nobs,/float)
           testvarerrornew=make_array(npobs,ntplot,nobs,/float)
           FOR o=0,nobs-1 DO BEGIN
              FOR p=0,npobs-1 DO BEGIN
                 obsvarnew[p,*,o]=nan_rebin(obsvar[p,0:ntobs-1,o],ntplot)
                 obsmasknew[p,*,o]=nan_rebin(obsmask[p,0:ntobs-1,o],ntplot)
                 obsvarerrornew[p,*,o]=nan_rebin(obsvarerror[p,0:ntobs-1,o],ntplot)
                 testvarnew[p,*,o]=nan_rebin(testvar[p,0:ntobs-1,o],ntplot)          
                 testmasknew[p,*,o]=nan_rebin(testmask[p,0:ntobs-1,o],ntplot)
                 testvarerrornew[p,*,o]=nan_rebin(testvarerror[p,0:ntobs-1,o],ntplot)
              ENDFOR
           ENDFOR
           obsmask=obsmasknew
           obsvarerror=obsvarerrornew
           testvar=testvarnew
           testvarerror=testvarerrornew
           testmask=testmasknew
           obsmask=obsmasknew
           obsvarerror=obsvarerrornew
                 
           obsvar=obsvarnew

           obsvarnew=0
           obsvarerrornew=0
           obsmasknew=0
           testvarnew=0
           testvarerrornew=0
           testmasknew=0
           obsmasknew=0
           ntobs=ntplot
           dtobs=dtplot
       ENDIF ELSE indexp=index1
   ENDIF

   IF plotmodels THEN BEGIN
       IF dtmod NE dtplot THEN BEGIN
           mtotmod=mtotmod/mtot*mtot
           ntmod=ntmod/ntplot*ntplot
           indexp=rebin(index2[0:mtotmod-1],mtot)*dtmod/dtplot
           
           modvarnew=make_array(ntplot,nmodels,/float)
           modvarerrornew=make_array(ntplot,nmodels,/float)
           FOR m=0,nmodels-1 DO BEGIN
              modvarnew[*,m]=nan_rebin(modvar[0:ntmod-1,m],ntplot)
              modvarerrornew[*,m]=nan_rebin(modvarerror[0:ntmod-1,m],ntplot)
           ENDFOR
           modvar=modvarnew
           modvarnew=0
           modvarerror=modvarerrornew
           modvarerrornew=0
           
           ntmod=ntplot
           dtmod=dtplot
       ENDIF ELSE indexp=index2
   ENDIF

   ;; integrative/cumulative plots
   IF integrate THEN BEGIN
       IF plotobs THEN BEGIN
           obsvarsqrt=obsvar
           FOR o=0,nobs-1 DO BEGIN
               badindex=where(finite(obsvar[*,o],/nan),count)
               IF count GT 0 THEN obsvar[badindex,o]=0.0
               IF startp NE 0 THEN obsvar[0:startp-1,o]=0.
               IF water THEN obsvar[*,o]=total(obsvar[*,o],/cumulative) ELSE  $
                 obsvar[*,o]=total(obsvar[*,o]*dtplot,/cumulative)
               IF count GT 0 THEN obsvar[badindex,o]=!VALUES.F_NAN   
           ENDFOR
       ENDIF
       IF plotmodels THEN BEGIN
           FOR m=0,nmodels-1 DO BEGIN
              badindex=where(finite(modvar[*,m],/nan),count)
              IF count GT 0 THEN modvar[badindex,m]=0.0
              IF startp NE 0 THEN modvar[0:startp-1,m]=0.
              IF water THEN modvar[*,m]=total(modvar[*,m],/cumulative) ELSE $
                 modvar[*,m]=total(modvar[*,m]*dtplot,/cumulative)
              IF count GT 0 THEN modvar[badindex,m]=!VALUES.F_NAN  
           ENDFOR
       ENDIF
   ENDIF

   ;; create the plot window
   aspect=1
   IF contour THEN aspectset=2.0
   IF scatter THEN aspectset=1.0

   IF (scatter EQ 0) AND (power EQ 0) AND (contour EQ 0) AND (taylor EQ 0) AND (histo EQ 0) THEN aspect=nyears

   IF aspectset GT 0 THEN aspect=aspectset

   ; create yaxis and plot titles
   !y.ticklen=0.01/aspect

   ;; create x axis (e.g. time) 
   IF nyears EQ 1 THEN BEGIN
       !x.ticklen=0.01
       xticks=12
       xminor=1
       xtickname=['      '+monstr_short,' ']
       xtitle=stime
   ENDIF ELSE BEGIN
       xticks=nyears
       xminor=12
       xtickname=replicate(' ',nyears+1)
       xtitle='Years'
       IF nyears LE 5 THEN !x.ticklen=0.5 ELSE !x.ticklen=0.01

   ENDELSE
   
   
   IF (plotmodels EQ 0) OR (modvarunits[0] EQ '') THEN modvarunits[*]=obsvarunits[0]
   IF (plotmodels EQ 0) OR (modvartitle[0] EQ '') THEN modvartitle[*]=obsvartitle[0]

   ; remove underscores from the sitename
   FOR s=0,nsites-1 DO BEGIN
       REPEAT BEGIN
           sn=sitenames[sites[s]-1]
           pos=strpos(sn,'_')
           IF pos GE 0 THEN strput,sn,' ',pos
           sitenames[sites[s]-1]=sn
       ENDREP UNTIL (pos LT 0)
   ENDFOR
   REPEAT BEGIN
       sn=sitename
       pos=strpos(sn,'_')
       IF pos GE 0 THEN strput,sn,' ',pos
       sitename=sn
   ENDREP UNTIL (pos LT 0)



   CASE 1 OF
      integrate: BEGIN
          IF water THEN BEGIN
              IF plotobs THEN ytitle=obsvartitle[0]+' ('+obsvarunits[0]+')' ELSE ytitle=modvartitle[0]+' ('+modvarunits[0]+')'
              ytitle='Integrated '+ytitle
              IF secondaxis THEN secondaxistitle='Integrated '+secondaxistitle
          ENDIF ELSE BEGIN
              IF plotobs THEN BEGIN
                  IF strmid(obsvarunits[0],strlen(obsvarunits[0])-1,1) EQ 's' THEN  $
                    obsvarunits[0]=strmid(obsvarunits[0],0,strlen(obsvarunits[0])-2)
                  IF strmid(obsvarunits[0],strlen(obsvarunits[0])-3,3) EQ 's-1' THEN  $
                    obsvarunits[0]=strmid(obsvarunits[0],0,strlen(obsvarunits[0])-3)
                  IF strmid(obsvarunits[0],0,4) EQ 'W/m2' THEN obsvarunits[0]='J/m2'
                  IF strmid(obsvarunits[0],0,5) EQ 'W m-2' THEN obsvarunits[0]='J/m2'
                  ytitle='Integrated '+ obsvartitle[0]+' ('+obsvarunits[0]+')'
              ENDIF ELSE BEGIN
                  IF strmid(modvarunits[0],strlen(modvarunits[0])-1,1) EQ 's' THEN  $
                    modvarunits[0]=strmid(modvarunits[0],0,strlen(modvarunits[0])-2)
                  IF strmid(modvarunits[0],strlen(modvarunits[0])-3,3) EQ 's-1' THEN  $
                    modvarunits[0]=strmid(modvarunits[0],0,strlen(modvarunits[0])-3)
                  IF strmid(modvarunits[0],0,4) EQ 'W/m2' THEN modvarunits[0]='J/m2'
                  IF strmid(modvarunits[0],0,5) EQ 'W m-2' THEN modvarunits[0]='J/m2'
                  ytitle='Integrated '+ modvartitle[0]+' ('+modvarunits[0]+')'
              ENDELSE
              IF secondaxis THEN BEGIN
                 IF strmid(secondaxisunits[0],strlen(secondaxisunits[0])-1,1) EQ 's' THEN $
                    secondaxisunits[0]=strmid(secondaxisunits[0],0,strlen(secondaxisunits[0])-2)
                 secondaxistitle='Integrated '+secondaxistitle
              ENDIF
          ENDELSE
          IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
      END
      contour: BEGIN
         ytitle='hours'
         IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
         IF plotobs THEN title=title+'OBS' ELSE BEGIN
             IF plotexperimentname THEN title=title+': '+experiments[0]
         ENDELSE
      END
      power: BEGIN
          IF plotobs THEN ytitle=obsvartitle[0]+' (power)' ELSE ytitle=modvartitle[0]+' (power)'
          IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
          xtitle='Frequency (per day) '+stime
      END
      taylor: BEGIN
          ytitle='standard deviation (normalized)'
          IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
          IF plotexperimentname THEN title=title ;+': '+experiments[0]   
          IF manysites THEN title = title + ': '+varnames[0]
          xtitle='standard deviation (normalized)'
      END
      scatter: BEGIN
          IF testvarname EQ '' THEN BEGIN
              IF plotobs THEN ytitle='MODEL: '+obsvartitle[0]+' ('+obsvarunits[0]+')' ELSE $
                ytitle='MODEL: '+modvartitle[0]+' ('+modvarunits[0]+')'
              IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
              title=sitename
          ENDIF ELSE BEGIN
              IF plotobs THEN ytitle=obsvartitle[0]+' ('+obsvarunits[0]+')' ELSE $
                ytitle=modvartitle[0]+' ('+modvarunits[0]+')'
              IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
              title=sitename
          ENDELSE
      END
      histo: BEGIN
          !x.ticklen=0.01
          ytitle='Normalized Density'
            IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
          ;title=title+': '+xtitle
      END
      ELSE: BEGIN
          ;IF nobs GT 1 THEN BEGIN
          ;    IF plotobs THEN ytitle='Mean Flux ('+obsvarunits[0]+')' ELSE ytitle='Mean ('+modvarunits[0]+')'
          ;ENDIF ELSE BEGIN
              IF plotobs THEN ytitle=obsvartitle[0]+' ('+obsvarunits[0]+')' ELSE ytitle=modvartitle[0]+' ('+modvarunits[0]+')'
          ;ENDELSE
          ;IF manymodels OR (plotmodels EQ 0) THEN title=sitename ELSE title=modelshort[0]+': '+sitename
          title=sitename
      END
   ENDCASE
     
   ; auto-guess the plot range and set y axis naming
   IF autorange THEN BEGIN
       minmod=make_array(nmodels,/float,value=!values.f_nan)
       maxmod=make_array(nmodels,/float,value=!values.f_nan)
       minobs=make_array(nobs,/float,value=!values.f_nan)
       maxobs=make_array(nobs,/float,value=!values.f_nan)
       IF plotobs THEN BEGIN
           FOR o=0,nobs-1 DO BEGIN
              automin=min(obsvar[*,indexp,o],max=automax,/NAN) 
              
              minobs[o]=automin
              maxobs[o]=automax
           ENDFOR
       ENDIF
       IF plotmodels THEN BEGIN
           FOR m=0,nmodels-1 DO BEGIN
               automin=min(modvar[indexp,m],max=automax,/NAN)
               minmod[m]=automin
               maxmod[m]=automax
           ENDFOR
       ENDIF
       automin=min(([minobs,minmod])[where(axisindex EQ 0)],/NAN)
       automax=max(([maxobs,maxmod])[where(axisindex EQ 0)],/NAN)
       range=[automin-0.05*abs(automin),automax+0.05*automax]
       IF max(axisindex) EQ 1 THEN BEGIN
           automin=min(([minobs,minmod])[where(axisindex EQ 1)],/NAN)
           automax=max(([maxobs,maxmod])[where(axisindex EQ 1)],/NAN)
           secondaxisrange=[automin-0.05*abs(automin),automax+0.05*automax]
       ENDIF
   ENDIF ELSE secondaxisrange=range2
   

   IF psout THEN BEGIN
      set_plot,'ps' 

      ysize = ysize_ps_default

      thick=ysize/8.0*5.0       ; plot line thickness
      charthick=2.0             ; charthickness relative to 1.0
      charsize=ysize/11.0     ; charsize relative to 10pt.

      !P.FONT=0
      print,'Writing EPS: ',filename
      DEVICE, /encapsul, bits_per_pixel=8,/color, font_size=5, $
         Filename=filename, preview=0,/isolatin1, $
         xsize=ysize*aspect, ysize=ysize,/HELVETICA

      IF colored THEN loadct,50,file='brewer.tbl' ELSE loadct,17,file='brewer.tbl'
      tvlct,red,green,blue,/get
      green=shift(green,1)
      blue=shift(blue,1)
      red=shift(red,1)
 
;      IF colored THEN loadct,25 ELSE loadct,0  
;      tvlct,red,green,blue,/get
      IF colored EQ 0 THEN BEGIN
         red=reverse(red)
         green=reverse(green)
         blue=reverse(blue)
      ENDIF
      red[0]=255B
      green[0]=255B
      blue[0]=255B
      ntable=!d.table_size
      ;red[ntable-2]=100B
      ;green[ntable-2]=100B
      ;blue[ntable-2]=100B
      red[ntable-1]=0B
      green[ntable-1]=0B
      blue[ntable-1]=0B
      tvlct,red,green,blue
   ENDIF ELSE BEGIN
      set_plot,'x'

      ysize = ysize_x_default
      xsize = ysize * aspect

      device,get_screen_size = screen_size
      ;IF xsize GT screen_size[0] THEN BEGIN
      ;   xsize = screen_size[0]
      ;   ysize = xsize / aspect
      ;ENDIF

      thick=2.0      ; plot line thickness
      charthick=1.5             ; charthickness relative to 1.0
      charsize=1.0 * sqrt(ysize/float(ysize_x_default))     ; charsize relative to 10pt.

      !P.FONT=1
      device,decomposed=0,set_font='HELVETICA',/tt_font
      window,0,xsize=xsize,ysize=ysize

      IF colored THEN loadct,50,file='brewer.tbl' ELSE loadct,17,file='brewer.tbl'
      tvlct,red,green,blue,/get
      green=shift(green,1)
      blue=shift(blue,1)
      red=shift(red,1)
      ;IF colored THEN loadct,25 ELSE loadct,0
      ;tvlct,red,green,blue,/get
      red[0]=0
      green[0]=0
      blue[0]=0
      ntable=!d.table_size
      ;red[ntable-2]=150B
      ;green[ntable-2]=150B
      ;blue[ntable-2]=150B
      red[ntable-1]=255B
      green[ntable-1]=255B
      blue[ntable-1]=255B
      tvlct,red,green,blue
 
         
   ENDELSE

   q = 0
   FOR ax=0,secondaxis DO BEGIN

       IF ax EQ 0 THEN BEGIN
          ;; plot axes and titles
           CASE 1 OF
               scatter: BEGIN
                   !x.ticklen=0.01
                   IF (obsvarnames[0] EQ varnames[0]) OR (plotmodels EQ 1) THEN BEGIN
                       IF testvarname NE '' THEN BEGIN
                           yrange=range
                           xrange=testrange
                           plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1, $
                             color=ntable-1, /nodata,$
                             xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                             xthick=2.0, ythick=2.0, xtitle=testvartitle[0]+ $
                             ' ('+testvarunits[0]+')', ytitle=ytitle, title=title, $
                             /normal,position=[0.18,0.10,1.-(0.05+0.1*secondaxis),0.92]
                       ENDIF ELSE BEGIN
                           xrange=range
                           yrange=range
                           plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1, $
                             color=ntable-1, /nodata,$
                             xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                             xthick=2.0, ythick=2.0, xtitle=obstitle[0]+': '+modvartitle[0]+ $
                             ' ('+modvarunits[0]+')', ytitle=ytitle, title=title, $
                             /normal,position=[0.18,0.10,1.-(0.05+0.1*secondaxis),0.92]
                       ENDELSE
                   ENDIF ELSE BEGIN
                       xrange=obsrange
                       yrange=range
                       plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1, $
                         color=ntable-1, /nodata,$
                         xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                         xthick=2.0, ythick=2.0, xtitle=obstitle+': '+obsvartitle[0]+ $
                         ' ('+obsvarunits[0]+')', ytitle=ytitle, title=title, $
                         /normal,position=[0.18,0.10,1.-(0.05+0.1*secondaxis),0.92]
                   ENDELSE
               END
               histo: BEGIN
                   yrange=[0,1]
                   xrange=range
                   plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1, $
                     color=ntable-1, /nodata,$
                     xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                     xthick=2.0, ythick=2.0, xtitle=modvartitle[0]+ $
                     ' ('+modvarunits[0]+')', ytitle=ytitle, title=title, $
                     /normal,position=[0.18,0.10,1.-(0.05+0.1*secondaxis),0.92]
               END
               contour: BEGIN
                   yrange=[0,24]
                   xrange=[startp,endp]
                   plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1, $
                     xticks=xticks,xtickname=xtickname, color=ntable-1, /nodata, $
                     ytickname=['00','06','12','18','24'] , yticks=4, $
                     xminor=xminor,yminor=3,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                     xthick=2.0, ythick=2.0, xtitle=xtitle, ytitle=ytitle, title=title, $
                     /normal, ticklen=0.02,position=[0.18,0.10,0.98,0.92],xticklen=-0.005,yticklen=-0.01
               END
               
               power: BEGIN          
                   !x.ticklen=0.01
                   n21=ntobs/2.+1
                   freq=indexp
                   freq[n21]=n21-ntobs+findgen(n21-2)
                   freq=freq/(ntobs*dtobs/(24.*3600.)) ; frequency in day-1
                   
                  
                   yrange=range
                   xrange=[1./365.,12.]
                   plot,freq[0:n21-1],freq[0:n21-1],xstyle=1,ystyle=1,yrange=yrange, $
                     color=ntable-1,xrange=xrange,/nodata,$
                     charsize=2.0*charsize,charthick=charthick, thick=thick, $
                     xthick=2.0, ythick=2.0, xtitle=xtitle, ytitle=ytitle, title=title, $
                     /normal,position=[0.18,0.10,0.95,0.92],/xlog,/ylog
                   
               END
               taylor: BEGIN
                   !x.ticklen=0.01
                   !y.ticklen=0.01
                   maxdev=2.0
                   xrange=[0.,maxdev]
                   yrange=[0.,maxdev]
                   plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=9,ystyle=9, $
                     color=ntable-1, /nodata,$
                     xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                     xthick=2.0, ythick=2.0, xtitle=xtitle, ytitle=ytitle, $
                     /normal,position=[0.15,0.125,0.85,0.825],/polar,noclip=1

                   npoints=100
                   oplot,replicate(1.0,npoints+1),findgen(npoints+1)/npoints*!pi/2.,/polar,noclip=1,color=ntable-1
                   oplot,replicate(maxdev,npoints+1),findgen(npoints+1)/npoints*!pi/2.,/polar,thick=thick,noclip=1,color=ntable-1

                   npoints=13
                   points=[findgen(11)/10.,0.95, 0.99]
                   FOR i=0,npoints-1 DO BEGIN
                       theta=acos(points[i])
                       r=maxdev*1.05
                       x=cos(theta)*r
                       y=sin(theta)*r
                       sint=''
                       sflt=''
                       corr=points[i]
                       reads,fix(corr),sint
                       sint=strtrim(sint,2)
                       reads,fix(100.*(corr-fix(corr))),sflt
                       sflt=strtrim(sflt,2)
                       IF (strlen(sflt) EQ 2) AND (strmid(sflt,1,1) EQ '0') THEN sflt=strmid(sflt,0,1)
                       oplot,[maxdev*0.95,maxdev],[theta,theta],/polar,thick=thick,noclip=1,color=ntable-1
                       xyouts,x,y,sint+'.'+sflt,charsize=2.0*charsize,charthick=charthick, orientation=theta/!pi*180.,color=ntable-1
                   ENDFOR
                   xyouts,maxdev*0.8,maxdev*0.8,'Correlation',charsize=2.0*charsize,charthick=charthick, orientation=-45., align=0.5,color=ntable-1
                   xyouts,0.5,0.94,title,charsize=2.5*charsize,charthick=charthick, align=0.5,/normal,color=ntable-1

               END
               ELSE: BEGIN
                   xrange=[startp,endp]
                   yrange=range
                   plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1+8*secondaxis, $
                     xticks=xticks,xtickname=xtickname,color=ntable-1, /nodata,$
                     xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                     xthick=2.0, ythick=2.0, xtitle=xtitle, title=title, $
                     /normal,position=[0.18/aspect,0.10,1.-(0.05+0.1*secondaxis)/aspect,0.92], ytitle=ytitle
                   
                   ;; create centered labels with year
                   IF nyears GT 1 THEN BEGIN
                      FOR y=0,nyears-1 DO BEGIN
                         dx = (1.-(0.05+0.1*secondaxis)/aspect - 0.18/aspect)/nyears
                         xpos = 0.18/aspect + dx/2. + y*dx
                         xyouts,xpos,.06,strtrim(years[y],2),color=ntable-1,charsize=2.0*charsize,charthick=charthick, $
                                Alignment=0.5,Orientation=0.0, /Normal
                      ENDFOR
                         

                   ENDIF

                   ;; different color y axis
                   ;xyouts,0.07/aspect,(0.92-0.1)/2.+0.1,ytitle,color=253,charsize=2.0*charsize,charthick=charthick, $
                   ;      Alignment=0.5,Orientation=90.0, /Normal
               END
           ENDCASE
           
       ENDIF ELSE BEGIN
           yrange=secondaxisrange
           axis,yaxis=1,yrange=yrange,ystyle=1, $
              color=ntable-1, /save,$
             yminor=1,charsize=2.0*charsize,charthick=charthick, $
             ythick=2.0,/normal, ytitle=secondaxistitle+' ('+secondaxisunits+')'

           ;; different color y axis
           ;xyouts,1.-0.07/aspect,(0.92-0.1)/2.+0.1,secondaxistitle+' ('+secondaxisunits+')', $
           ;       color=50,charsize=2.0*charsize,charthick=charthick, $
           ;       Alignment=0.5,Orientation=90.0, /Normal
       ENDELSE
           
   ; plot data
       CASE 1 OF
          integrate: BEGIN
             IF plotobs THEN BEGIN
                FOR o=0,nobs-1 DO BEGIN
                   IF axisindex[o] EQ ax THEN BEGIN
                      
                      overplot,indexp,obsvar[indexp,o],obstitle[o],o,1    
                      
                   ENDIF
                ENDFOR
             ENDIF
             IF plotmodels THEN BEGIN
                FOR m=0,nmodels-1 DO IF axisindex[nobs+m] EQ ax THEN $
                   overplot,indexp,modvar[indexp,m],modtitle[m],m+nobs,1
             ENDIF
          END
          power: BEGIN
             IF plotobs THEN FOR o=0,nobs-1 DO powerplot,dtplot,indexp,obsvar[indexp,o],obstitle,0         
             IF plotmodels THEN FOR m=0,nmodels-1 DO powerplot,dtplot,indexp,modvar[indexp,m],modtitle[m],m+1
          END
          scatter: BEGIN
             IF testvarname EQ '' THEN BEGIN

                ;; plot obs with their spatial variability
                obsvartemp=obsvar[*,*,0]
                obserrortemp=obsvarerror[*,*,0]
                
                obsmasktemp=round(obsmask[*,*,0])
                badobs=where(obsmasktemp EQ 0,count)
                IF count GT 0 THEN BEGIN
                   obserrortemp[badobs]=!values.f_nan
                   obsvartemp[badobs]=!values.f_nan
                ENDIF
                
                goodobsall=where(finite(obsvartemp),goodobsallcount)
                
                ;; created weigthed mean of best obs
                IF scatterall THEN BEGIN
                   obstemp=make_array(goodobsallcount,/float)
                   modtemp=make_array(goodobsallcount,nmodels,/float)
                ENDIF ELSE BEGIN
                   obstemp=make_array(ntplot,/float,value=!values.f_nan)
                   obsminmaxtemp = make_array(2,ntplot,/float,value=!values.f_nan)
                ENDELSE
                
                idx = 0
                
                FOR t=0L,ntplot-1 DO BEGIN
                   
                   goodobs=where(finite(obsvartemp[*,t]),count)
                   
                   IF scatterall THEN BEGIN
                      IF count GT 0 THEN BEGIN
                         obstemp[idx:idx+count-1] = obsvartemp[goodobs,t]
                         FOR m=0,nmodels-1 DO BEGIN
                            modtemp[idx:idx+count-1,m] = modvar[t,m]
                         ENDFOR
                         idx=idx+count
                      ENDIF
                   ENDIF ELSE BEGIN
                      ;;print,t,count,mean(obserrortemp[*,t],/nan)
                      IF (count GE 2) THEN BEGIN   
                         obstemp[t]=mean(obsvartemp[goodobs,t],/nan)
;                                    obstemp[t]=total(obsvartemp[goodobs,t]/obserrortemp[goodobs,t],/nan) $
;                                               /total(1./obserrortemp[goodobs,t],/nan)
                         obsminmaxtemp[0,t] = mean(obsvartemp[goodobs,t]+obserrortemp[goodobs,t]*errscale,/nan)
                         obsminmaxtemp[1,t] = mean(obsvartemp[goodobs,t]-obserrortemp[goodobs,t]*errscale,/nan)
                      ENDIF
                   ENDELSE
                   
                ENDFOR
                
                FOR m=0,nmodels-1 DO BEGIN
                   modminmaxtemp = make_array(2,ntplot,/float,value=!values.f_nan)
                   modminmaxtemp[0,*]=modvar[*,m]+modvarerror[*,m]
                   modminmaxtemp[1,*]=modvar[*,m]-modvarerror[*,m]
                   IF npobs EQ 1 THEN BEGIN
                      scatplot,obsvar[indexp],modminmaxtemp[*,indexp],modtitle[m],m+nobs
                      scatplot,obsvar[indexp],reform(modvar[indexp,m]),modtitle[m],m+nobs  
                   ENDIF ELSE BEGIN
                      IF scatterall THEN BEGIN
                         scatplot,obstemp,modtemp[*,m],modtitle[m],m+nobs  
                      ENDIF ELSE BEGIN
                         scatplot,obsminmaxtemp[*,indexp],modminmaxtemp[*,indexp],modtitle[m],m+nobs
                         scatplot,obstemp[indexp],reform(modvar[indexp,m]),modtitle[m],m+nobs  
                      ENDELSE
                   ENDELSE
                ENDFOR
                
             ENDIF ELSE BEGIN
                ;; plotting with testvar

                FOR o=0,nobs-1 DO BEGIN
                   ;; plot obs with their spatial variability
                   obsvartemp=obsvar[*,*,o]
                   obserrortemp=obsvarerror[*,*,o]
                   obsmasktemp=round(obsmask[*,*,o])
                   
                   badobs=where(obsmasktemp EQ 0,count)
                   IF count GT 0 THEN BEGIN
                      obserrortemp[badobs]=!values.f_nan
                      obsvartemp[badobs]=!values.f_nan
                   ENDIF
                   
                   obstemp=make_array(ntplot,/float,value=!values.f_nan)
                   ;; created weigthed mean of best obs
                   FOR t=0L,ntplot-1 DO BEGIN
                      goodobs=where(finite(obsvartemp[*,t]),count)
                      ;; print,t,count,mean(obserrortemp[*,t],/nan)
                      IF (count GE 2) THEN BEGIN
                         obstemp[t]=mean(obsvartemp[goodobs,t],/nan) 
;                                     obstemp[t]=total(obsvartemp[goodobs,t]/obserrortemp[goodobs,t],/nan)/ $
;                                                total(1./obserrortemp[goodobs,t],/nan)
                      ENDIF
                   ENDFOR
                   
                   testvartemp=testvar[*,*,o]
                   testerrortemp=testvarerror[*,*,o]
                   testmasktemp=round(testmask[*,*,o])
                   
                   badtest=where(testmasktemp EQ 0,count)
                   IF count GT 0 THEN BEGIN
                      testerrortemp[badobs]=!values.f_nan
                      testvartemp[badobs]=!values.f_nan
                   ENDIF
                   
                   testtemp=make_array(ntplot,/float,value=!values.f_nan)
                   ;; created weigthed mean of best obs
                   FOR t=0L,ntplot-1 DO BEGIN
                      goodtest=where(finite(testvartemp[*,t]),count)
                      IF (count GE 2) THEN BEGIN
                         testtemp[t]=mean(testvartemp[goodtest,t],/nan) 
;                                     testtemp[t]=total(testvartemp[goodtest,t]/testerrortemp[goodtest,t],/nan)/ $
;                                                 total(1./testerrortemp[goodtest,t],/nan) 
                      ENDIF
                   ENDFOR
                   
                   scatplot,testvartemp[*,indexp],obsvartemp[*,indexp],obstitle,o
                   scatplot,testtemp[indexp],obstemp[indexp],obstitle,o
                ENDFOR
                
                IF plotmodels THEN BEGIN
                   FOR m=0,nmodels-1 DO BEGIN
                      modminmaxtemp = make_array(2,ntplot,/float,value=!values.f_nan)
                      modminmaxtemp[0,*]=modvar[*,m]+modvarerror[*,m]
                      modminmaxtemp[1,*]=modvar[*,m]-modvarerror[*,m]
                      scatplot,testvartemp[*,indexp],modminmaxtemp[*,indexp],modtitle[m],m+nobs
                      scatplot,testtemp[indexp],reform(modvar[indexp,m]),modtitle[m],m+nobs          
                   ENDFOR
                ENDIF
                   
             ENDELSE
          END
          taylor: BEGIN
             nplot=ntaylor+1
             FOR m=0,nmodels-1 DO BEGIN
                
                ;; plot obs with their spatial variability
                obsvartemp=obsvar[*,*,m]
                obserrortemp=obsvarerror[*,*,m]
                
                ;; created weigthed mean of best obs
                obstemp=make_array(ntplot,/float,value=!values.f_nan)
                obsminmaxtemp = make_array(2,ntplot,/float,value=!values.f_nan)
                FOR t=0L,ntplot-1 DO BEGIN
                   goodobs=where(finite(obsvartemp[*,t]) AND finite(obserrortemp[*,t]),count)
                                ;print,t,count,mean(obserrortemp[*,t],/nan)
                   IF (count GE 2) THEN BEGIN   
                      obstemp[t]=mean(obsvartemp[goodobs,t],/nan)
;                           obstemp[t]=total(obsvartemp[goodobs,t]/obserrortemp[goodobs,t],/nan) $
;                                      /total(1./obserrortemp[goodobs,t],/nan)
                      obsminmaxtemp[0,t] = mean(obsvartemp[goodobs,t]+obserrortemp[goodobs,t]*errscale,/nan)
                      obsminmaxtemp[1,t] = mean(obsvartemp[goodobs,t]-obserrortemp[goodobs,t]*errscale,/nan)
                   ENDIF
                   
                ENDFOR
                
                modtemp = modvar[*,m]
                
                IF m GT ntaylor THEN legend=0
                IF manysites THEN $
                   taylorplot,obstemp[indexp],modtemp[indexp],modtitle[m],(m-m/ntaylor*ntaylor)+1,ssites[m] $
                ELSE $
                   taylorplot,obstemp[indexp],modtemp[indexp],modtitle[m],(m-m/ntaylor*ntaylor)+1,varnames[m]
             ENDFOR
             
          END
          contour: BEGIN
             IF plotobs THEN BEGIN
                contourplot,indexp,obsvar[indexp,0],obstitle,0,dtplot           
             ENDIF ELSE BEGIN
                contourplot,indexp,modvar[indexp,0],modtitle[0],1,dtplot
             ENDELSE
          END
          histo: BEGIN
             IF plotobs THEN BEGIN
                FOR o=0,nobs-1 DO histoplot,indexp,obsvar[*,indexp,o],obstitle[o],o
             ENDIF
             IF plotmodels THEN BEGIN
                FOR m=0,nmodels-1 DO histoplot,indexp,modvar[indexp,m],modtitle[m],nobs+m
             ENDIF
          END
          ELSE: BEGIN

             IF plotobs THEN BEGIN   
                FOR o=0,nobs-1 DO IF axisindex[o] EQ ax THEN BEGIN
                   IF obstitle[o] NE '' THEN BEGIN
                      ;; plot obs with error
                      obsvartemp=obsvar[*,*,o]
                      obsdevtemp=make_array(ntplot,/float,value=!values.f_nan)
                      obserrortemp=obsvarerror[*,*,o]
                      obstemp=make_array(ntplot,/float,value=!values.f_nan)
                      obsminmaxtemp = make_array(2,ntplot,/float,value=!values.f_nan)
                      
                      ;; created weigthed mean of best obs
                      FOR t=0L,ntplot-1 DO BEGIN

                         goodobs=where(finite(obsvartemp[*,t]) AND finite(obserrortemp[*,t]),count)
                         IF (count GE 1) THEN BEGIN   
                            obstemp[t]=mean(obsvartemp[goodobs,t],/nan)
                            obsminmaxtemp[0,t] = obstemp[t] + mean(obserrortemp[goodobs,t],/nan)
                            obsminmaxtemp[1,t] = obstemp[t] - mean(obserrortemp[goodobs,t],/nan)
                         ENDIF
                         
                      ENDFOR
                      IF total(finite(obstemp[indexp])) GT 0.0 THEN BEGIN
                         overplot,indexp,obsminmaxtemp[*,indexp],obstitle[o],q,0
                         overplot,indexp,obstemp[indexp],obstitle[o],q,0
                         q=q+1
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
             
             IF plotmodels THEN BEGIN
                FOR m=0,nmodels-1 DO BEGIN
                   IF axisindex[nobs+m] EQ ax THEN BEGIN
                      modminmaxtemp = make_array(2,ntplot,/float,value=!values.f_nan)
                      modminmaxtemp[0,*]=modvar[*,m]+modvarerror[*,m]
                      modminmaxtemp[1,*]=modvar[*,m]-modvarerror[*,m]
                      overplot,indexp,modminmaxtemp[*,indexp],modtitle[m],q,0
                      overplot,indexp,reform(modvar[indexp,m]),modtitle[m],q,1
                      q=q+1
                   ENDIF
                ENDFOR
             ENDIF

          END
       ENDCASE
    ENDFOR
   
   ;; plot the color bar
   IF contour THEN BEGIN
       IF plotmodels THEN colbar,divisions=6,ncolors=60, units=modvarunits[0], varname=modvartitle[0] $
       ELSE colbar,divisions=6,ncolors=60, units=obsvarunits[0], varname=obsvartitle[0]
   ENDIF

   IF psout THEN BEGIN
       DEVICE,/CLOSE      
       set_plot,'x'

       pixelsize=1024
       spixelsize=''
       reads,pixelsize*aspect,spixelsize
       spawn,'convert -density 300 -resize '+spixelsize+' -flatten '+filename+' '+filename+'.png'
;       spawn,'convert -quality 90 '+filename+'.png '+filename+'.jpg'
       spawn,'epstopdf '+filename
       
   ENDIF


   ; latex table output

   latextable = 0

   IF latextable THEN BEGIN

      ;text = strtrim(sitenames[sites[0]-1],2)
      text = string(sites[0],format='(I3)')

      FOR m=0,nmodels-1 DO BEGIN

         y = modvar[*,m]

         ;; plot obs with their spatial variability
         obsvartemp=obsvar[*,*,m]
         obsdevtemp=make_array(ntplot,/float,value=!values.f_nan)
         obserrortemp=obsvarerror[*,*,m]
         obsmasktemp=round(obsmask[*,*,m])
         obstemp=make_array(ntplot,/float)
         badobs=where(obsmasktemp EQ 0,count)
         IF count GT 0 THEN BEGIN
            obserrortemp[badobs]=!values.f_nan
            obsvartemp[badobs]=!values.f_nan
         ENDIF
         ;; created weigthed mean of best obs
         FOR t=0L,ntplot-1 DO BEGIN
            goodobs=where(finite(obsvartemp[*,t]) AND finite(obserrortemp[*,t]),count)
            ;;print,t,count,mean(obserrortemp[*,t],/nan)
            IF (count GE 2) THEN BEGIN   
               ;;obstemp[t]=mean(obsvartemp[goodobs,t],/nan)
               obstemp[t]=total(obsvartemp[goodobs,t]/obserrortemp[goodobs,t],/nan) $
                          /total(1./obserrortemp[goodobs,t],/nan)
            ENDIF ELSE obstemp[t]=!values.f_nan 
         ENDFOR

         x = obstemp

         index=where(finite(x) AND finite(y),count)

         IF count GT 0 THEN BEGIN
            stdx=stddev(x[index],/nan)
            stdy=stddev(y[index],/nan)
            stdnorm=stdy/stdx
            
            meany=mean(y[index],/nan)
            meanx=mean(x[index],/nan)
            bias=meany-meanx
            
            cp_rms=sqrt( total(((y[index]-meany) - (x[index]-meanx))^2 )/float(count) )
                        
            theta=acos(correlate(x[index],y[index]))
            xc=cos(theta)*stdnorm
            yc=sin(theta)*stdnorm
            
            ;r=correlate(x[index],y[index])
            ;var=1./(float(count)-3.0)
            ;zvalue = 0.5*alog((1.+r)/(1.-r))/sqrt(var)
            ;if (zvalue lt 0.) then pvalue = 2.*(gauss_pdf(zvalue)) $
            ;else pvalue = 2.*(1.-gauss_pdf(zvalue))

            ; pvalue = (tm_test(x[index],y[index]))[1]

            ;; works:
            is_significant = CORR_TTEST(corr=correlate(x[index],y[index]),N=count,sl=0.0001, ro=0.)

            ;; cp_rms2=sqrt(stdy^2 + stdx^2 - 2.*stdx*stdy*cos(theta))

            rms=sqrt(total((y[index] - x[index])^2.0)/float(count))
            rms2=sqrt(bias^2 + cp_rms^2)

            text=text+' & '
            IF is_significant THEN text = text + '\textbf{'
            text=text+strtrim(string(cos(theta),format='(f5.2)'),2)+' ('+ strtrim(string(rms,format='(f5.2)'),2)+') '
            IF is_significant THEN text = text + '} '
            
;            IF (pvalue LE 0.01) AND (pvalue GT 0.001) THEN text = text + '\tablenotemark{a}'
;            IF (pvalue LE 0.001) AND (pvalue GT 0.0001) THEN text = text + '\tablenotemark{b}'
;            IF (pvalue LE 0.0001) THEN text = text + '\tablenotemark{c}'

         ENDIF
      ENDFOR
      
      text = text + ' \\'

      print,text
      ;openw,25,plotdir+'/'+'lai_statistics.txt'
      printf,25,text
      ;close,25


   ENDIF

   latextable2 = 0

   IF latextable2 THEN BEGIN

      ;text = strtrim(sitenames[sites[0]-1],2)
      text = string(sites[0],format='(I3)')

      FOR m=0,nmodels-1 DO BEGIN

         y = modvar[*,*,m]

         a=where(finite(y[0,*]),count)
         lastvalid = a[count-1]

         mval = mean(y[*,lastvalid],/nan)
         stdval = stddev(y[*,lastvalid],/nan)

         text=text +' & '+strtrim(string(mval,format='(f8.1)'),2)+' $\pm$ '+strtrim(string(stdval,format='(f8.1)'),2)
         ;text=text +' & '+strtrim(string(mval,format='(f8.2)'),2)+' $\pm$ '+strtrim(string(stdval,format='(f8.2)'),2)

      ENDFOR

      text = text + ' \\'
      print,text
      ;openw,25,plotdir+'/'+'climate_param_statistics.txt'
      ;openw,25,plotdir+'/'+'structural_param_statistics.txt'
      ;printf,25,text
      ;close,25

   ENDIF

   ; curve fitting

   fit=0

   IF fit THEN BEGIN
      ; fpar/lai fitting
      goodobs = where (finite(testvar) AND finite(obsvar))
      x=testvar[goodobs]
      y=obsvar[goodobs]
      ; vcover, fparsat, laisat
      init_coefs=[0.9,0.99,5.0]
      ;init_coefs=[0.745,0.95,3.7]
      fita=[1,1,1]
      weights=1/y
      ;result = lmfit(x,y,init_coefs,measure_errors=0.1*y,fita=fita,function_name='fparlaifit',/double)
      ;result = curvefit(x,y,weights,init_coefs,sigma,function_name='fparlaifit2')
      lai = alog( (1.-testvar/(init_coefs[0]>1.e-3))>1.e-3 ) / alog( (1.- init_coefs[1])>0.001 ) * init_coefs[2] * init_coefs[0]

      plot,testvar,obsvar,psym=4,color=200
      oplot,testvar,lai,psym=2,color=250

      stop

   ENDIF

   
   ; ASCII file output
 
   IF printmonthlylai THEN BEGIN
       ; ASCII monthly phenology output
       
       monthvar=make_array(12,/float,value=!values.f_nan)
       montherr=make_array(12,/float,value=!values.f_nan)
       
       ;badobs=where((obsmask EQ 0) OR finite(obsvar,/nan))
       ;obsvar[badobs]=!values.f_nan
       ;obsvarerror[badobs]=!values.f_nan
              
       badmod=where(finite(modvar,/nan),badcount)
       IF badcount GT 0 THEN modvar[badmod]=!values.f_nan
              
       FOR m=0,11 DO BEGIN
           dt=ntplot/float(nyears)/12.
           nt=ntplot/float(nyears)
           
           index=round((m*dt+findgen(dt))#replicate(1.0,nyears)+replicate(nt,dt)#findgen(nyears))
           
           ;goodobs=where(finite(obsvar[*,index]),goodcount)
           goodmod=where(finite(modvar[*,index]),goodcount)
           IF goodcount GT nyears THEN BEGIN
               ;monthvar[m]=total(obsvar[*,index]/obsvarerror[*,index],/nan)/total(1./obsvarerror[*,index],/nan)
               ;montherr[m]=stddev(obsvar[*,index],/nan)
               monthvar[m]=mean(modvar[*,index])
               montherr[m]=stddev(modvar[*,index],/nan)
           ENDIF
           
           
       ENDFOR
       
       
       badobs=where(finite(monthvar,/nan),badcount)
       goodobs=where(finite(monthvar),goodcount)

       IF badcount GT 0 THEN BEGIN
           monthvar1=interpol(monthvar[goodobs],goodobs,indgen(12))
           monthvar=shift(monthvar,6)
           goodobs=where(finite(monthvar),goodcount)
           monthvar2=shift(interpol(monthvar[goodobs],goodobs,indgen(12)),6)
           monthvar=(monthvar1+monthvar2)/2.0
       ENDIF

       ;badobs=where(finite(montherr,/nan),badcount)
       ;IF badcount GT 0 THEN montherr[badobs]=nodata

       ;filedir=laidir+'modisdata/'+sitename+'/'
       filedir=laidir+'modeldata/'+sitename+'/'
       if file_test(filedir,/directory) eq 0 THEN file_mkdir,filedir
       filename=filedir+strlowcase(varnames[0])+'.dat'
       openw,1,filename
       ;printf,1,sitename,',',varnames[0],' [',obsvarunits[0],']',', 12 satellite estimates'
       printf,1,sitename,',',varnames[0],' [',obsvarunits[0],']',', 12 model predictions'
       printf,1,'DoY     Mean    Stddev'
       FOR m=0,11 DO BEGIN
           printf,1,(monsum[m]+monsum[m+1])/2,' ',monthvar[m],' ',montherr[m],format='(i4,a1,f7.2,a1,f7.2)'
       ENDFOR
       close,1
   ENDIF

   IF printlai THEN BEGIN
       ; ASCII phenology output
       
       nt=ntplot/nyears

       var=make_array(nt,/float,value=!values.f_nan)
       err=make_array(nt,/float,value=!values.f_nan)
       
       badobs=where((obsmask EQ 0) OR finite(obsvar,/nan))
       obsvar[badobs]=!values.f_nan
       obsvarerror[badobs]=!values.f_nan
       
       
       FOR m=0,nt-1 DO BEGIN
           dt=1
           
           index=round((m*dt+findgen(dt))#replicate(1.0,nyears)+replicate(nt,dt)#findgen(nyears))
           
           goodobs=where(finite(obsvar[*,index]),goodcount)
           IF goodcount GT nyears THEN BEGIN
               var[m]=total(obsvar[*,index]/obsvarerror[*,index],/nan)/total(1./obsvarerror[*,index],/nan)
               err[m]=stddev(obsvar[*,index],/nan)
           ENDIF
           
           
       ENDFOR
       
       badobs=where(finite(var,/nan),badcount)
       IF badcount GT 0 THEN var[badobs]=nodata

       badobs=where(finite(err,/nan),badcount)
       IF badcount GT 0 THEN err[badobs]=nodata

       filedir=laidir+'modisdata/'+sitename+'/'
       if file_test(filedir,/directory) eq 0 THEN file_mkdir,filedir
       filename=filedir+strlowcase(varnames[0])+'.'+stime+'.dat'
       openw,1,filename
       printf,1,sitename,',',varnames[0],' [',obsvarunits[0],']',', ',nt,' measurements'
       printf,1,'DOY, Mean, Stddev'
       FOR m=0,nt-1 DO BEGIN
           printf,1,round((m+0.5)*dtplot/3600./24.),' ',var[m],' ',err[m]
       ENDFOR
       close,1
   ENDIF

   IF printascii THEN BEGIN
       ; ASCII file output
              
       nt=ntplot/nyears
       badmod=where(finite(modvar,/nan),badcount)
       IF badcount GT 0 THEN modvar[badmod]=nodata


       filename=plotdir+'/'+sitename+'.'+strlowcase(varnames[0])+'.ascii.'+stime+'.dat'
       openw,lun,filename,/get_lun
       printf,lun,sitename,': Variable: ',varnames[0], ' (',modvartitle[0],') Units: [',modvarunits[0],']'
       printf,lun,'Year ','DOY ','        Mean','       Stddev'
       FOR y=0,nyears-1 DO BEGIN
          FOR m=0,nt-1 DO BEGIN
             printf,lun,years[y],' ',round((m+0.5)*dtplot/3600./24.),' ',reform(modvar[0,m+y*nt]), $
                    ' ',reform(modvar[1,m+y*nt]), format='(I4,A1,I3,A1,F12.4,A1,F12.4)'
          ENDFOR
       ENDFOR
       free_lun,lun
   ENDIF

   IF interannual THEN BEGIN

      harvard = 0
      mmsf = 0
      swiss = 1
 
      thval=0.5
      
      sos_o = replicate(!values.f_nan,nyears)
      eos_o = replicate(!values.f_nan,nyears)
      sos_o_e = replicate(0.,nyears)
      eos_o_e = replicate(0.,nyears)

      IF swiss THEN BEGIN

         temp=''
         y=0
         dates=make_array(4)
         openr,1,datadir+'/modelfarm/tower_data/selected/phenology/Swiss_Lowland/INDEX.2006.uncert.dat'
         readf,1,temp
         
         WHILE EOF(1) EQ 0 DO BEGIN
            readf,1,temp
            dates[*] = strsplit(temp,';',/extract,/preserve_null)
            yi = where(years eq round(dates[0]),ycount)
            IF ycount GT 0 THEN BEGIN
               sos_o[yi] = dates[1]
               sos_o_e[yi] = (dates[2]-dates[1])/2.
            ENDIF
         END
         close,1

      ENDIF

      IF mmsf THEN BEGIN
         threshold=(max(obsvar,/nan)-min(obsvar,/nan))*thval + min(obsvar,/nan)

         nt=ntplot/float(nyears)
         FOR y=0,nyears-1 DO BEGIN          
            index = indgen(nt-50) + y*nt + 50
            isgrowing=where(obsvar[index] GT threshold,growcount)
            isdormant=where(obsvar[index] LE threshold,dormcount)
            lastday = where(isdormant LT isgrowing[0],lastdaycount)          ; last dormant day (spring)
            firstday = where(isdormant GT isgrowing[growcount-1],firstdaycount) ; first dormant (autumn)
            IF (growcount GE 1) AND (lastdaycount GE 1) THEN BEGIN
               val1=obsvar[index[isgrowing[0]]]
               val0=obsvar[index[isdormant[lastday[lastdaycount-1]]]]
               s0 = (threshold - val0)/(val1-val0)
               sos_o[y]=(s0*isgrowing[0]+(1.-s0)*isdormant[lastday[lastdaycount-1]])+1+50
               sos_o_e[y] = (isgrowing[0] - isdormant[lastday[lastdaycount-1]])/2.
            ENDIF
            IF (growcount GE 1) AND (firstdaycount GE 1) THEN BEGIN
               val0=obsvar[index[isgrowing[growcount-1]]]
               val1=obsvar[index[isdormant[firstday[0]]]]
               s0 = (threshold - val0)/(val1-val0)
               eos_o[y]=((1.-s0)*isgrowing[growcount-1]+s0*isdormant[firstday[0]])+1+50
               eos_o_e[y] = (isdormant[firstday[0]]-isgrowing[growcount-1])/2.
            ENDIF
         ENDFOR

         ; CO2
         sos_o=[117,110,109,116,108,108,109]
         eos_o=[293,296,297,296,296,306,297]

      ENDIF
 
      ; read harvard forest sos and eos data
      IF harvard THEN BEGIN
         openr,1,datadir+'/modelfarm/tower_data/selected/phenology/Harvard_Forest/Harvard_Forest.spr.dat'
         openr,2,datadir+'/modelfarm/tower_data/selected/phenology/Harvard_Forest/Harvard_Forest.fall.dat'
         fval = 0.
         ferr = 0.
         ytest = 0
         FOR y=0,nyears-1 DO BEGIN
            readf,1,ytest,fval,ferr
            sos_o[y]=fval
            sos_o_e[y]=ferr
            readf,2,ytest,fval,ferr
            eos_o[y]=fval
            eos_o_e[y]=ferr
         ENDFOR
         close,1
         close,2
      ENDIF

      sos_m = make_array(nyears, nmodels,/float,value=!values.f_nan)
      eos_m = make_array(nyears, nmodels,/float,value=!values.f_nan)

      
      FOR m=0,nmodels-1 DO BEGIN

        ; if m EQ 1 THEN thval=0.25
        ; if m EQ 2 THEN thval=0.75

         threshold=(max(modvar[*,m],/nan)-min(modvar[*,m],/nan))*thval + min(modvar[*,m],/nan)
         ;threshold=median(modvar[*,m])
         nt=ntplot/float(nyears)
         FOR y=0,nyears-1 DO BEGIN           
            index = indgen(nt-50) + y*nt + 50
            isgrowing=where(modvar[index,m] GT threshold,growcount)
            if growcount GT 1 THEN sos_m[y,m]=isgrowing[0]+1+50
            if growcount GT 1 THEN eos_m[y,m]=isgrowing[growcount-1]+1+50

         ENDFOR
      ENDFOR

      secondaxis=0

      IF psout THEN BEGIN
         set_plot,'ps' 
         filename = plotdir+'/'+sitename+'.'+strlowcase(varnames[0])+'.sos.'+stime+'.eps'
         print,'Writing EPS: ',filename
         DEVICE, /encapsul, bits_per_pixel=8,/color, font_size=5, $
                 Filename=filename, preview=0,/isolatin1, $
                 xsize=ysize*aspect, ysize=ysize,/HELVETICA
         tvlct,red,green,blue

      ENDIF ELSE BEGIN
         window,1,xsize=xsize,ysize=ysize
         tvlct,red,green,blue
      ENDELSE

      ;experiments=['PROGNOSTIC (50%)','PROGNOSTIC (25%)','PROGNOSTIC (75%)']
      ;experiments=['PROGNOSTIC LAI','TEMP_FAC','LIGHT_FAC']

      x = sos_o
      x_e = sos_o_e
      y = sos_m
      ytitle = 'Start of Season (days)'


;      x = eos_o
;      x_e = eos_o_e
;      y = eos_m
;      ytitle = 'End of Season (days)'

      xrange = [-0.5,nyears-0.5]
      yrange=[min([min(x-x_e,/nan),min(y,/nan)],/nan)-5,$
              max([max(x+x_e,/nan),max(y,/nan)],/nan)+5]

      plot,indgen(nyears),replicate(0,nyears),xstyle=1,ystyle=1,xrange = xrange, yrange=yrange, /nodata, xthick=2.0, ythick=2.0, $
           thick=thick,charsize=2*charsize,charthick=charthick,xminor=1,yminor=1, $
           xtitle='Years', ytitle = ytitle, color = ntable-1,xticks=nyears, xtickname=replicate(' ',nyears+1), $
           position=[0.18/aspect,0.10,1.-(0.05+0.1*secondaxis)/aspect,0.92], title = sitename

      FOR year=0,nyears-1,3 DO BEGIN
         dx = (1.-(0.05)/aspect - 0.18/aspect)/nyears
         xpos = 0.18/aspect + dx/2. + year*dx
         xyouts,xpos,.06,strtrim(years[year],2),color=ntable-1,charsize=1.5*charsize,charthick=charthick, $
                Alignment=0.5,Orientation=0.0, /Normal 
      ENDFOR

      FOR m=0,nmodels-1 DO BEGIN 
         overplot,indgen(nyears),y[*,m],experiments[m],m+1,1
         good = where(finite(y[*,m]),goodcount)
         a=regress(indgen(goodcount),y[good,m],const=b,mcorrelation=r, $
                   measure_errors=replicate(0.5,goodcount),/double, ftest=fvalue)
         print,experiments[m]+' trend '+string(a)

         index = where(finite(x) AND finite(y[*,m]),count)
         
         theta=acos(correlate(x[index],y[index,m]))
         rms=sqrt(total((y[index,m] - x[index])^2.0)/float(count))
         bias = total(y[index,m] - x[index])/float(count)
         
         r=correlate(x[index],y[index,m])
         var=1/(count-3.0)
         zvalue = 0.5*alog((1+r)/(1-r))/sqrt(var)
         if (zvalue lt 0) then pvalue = 2*(gauss_pdf(zvalue)) $
         else pvalue = 2*(1-gauss_pdf(zvalue))

         print,experiments[m]+' correlation: ',cos(theta), r,' bias: ',bias,' rms: ',rms, ' pvalue: ',pvalue
      ENDFOR

      overplot,indgen(nyears),transpose([[x-x_e],[x+x_e]]),'OBS',0,1
      overplot,indgen(nyears),x,'OBS',0,1

      ;oplot,x,thick=thick,color=255
      good = where(finite(x),goodcount)
      a=regress(indgen(goodcount),x[good],const=b,mcorrelation=r, $
                measure_errors=replicate(0.5,goodcount),/double)
      print,'OBS trend '+string(a)


      IF psout THEN BEGIN
         DEVICE,/CLOSE      
         set_plot,'x'
         
         pixelsize=800
         spixelsize=''
         reads,pixelsize*aspect,spixelsize
         spawn,'convert -density 300 -resize '+spixelsize+' -flatten '+filename+' '+filename+'.png'
         spawn,'convert -quality 90 '+filename+'.png '+filename+'.jpg'
         spawn,'epstopdf '+filename
      ENDIF

      sos_m = string(sos_m)
      nothing = where(variability EQ -9999,nothingcount)
      IF (nothingcount GT 0) THEN svariability[nothing] = 'NA'
      
      filedir=plotdir+'/'
      filename=filedir+sitename+'.'+strlowcase(varnames[0])+'.interannual.'+stime+'.dat'
      openw,1,filename, width=120
      printf,1,sitename,': SOS (DOY) from ',thval,' x ',varnames[0],' magnitude threshold, ',nyears,' years'
      printf,1,'Year  ',experiments
      FOR y=0,nyears-1 DO BEGIN
         printf,1,years[y],reform(svariability[y,*])
      ENDFOR
      close,1

   ENDIF

   stop


END

