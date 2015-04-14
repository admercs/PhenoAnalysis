COMMON SHARE,ntable,colored,charsize,charthick,legend,nplot,thick, regline, range, $
   xrange, yrange, length, nyears, legenddx, legenddy, aspect, pollon, pollat, psout, $
   ysize_ps_default, ysize_x_default, xsize, ysize

  COMMON MODISTILE,vmodis,hmodis,xmodis,ymodis,xidmodis,yidmodis,modis_lon,modis_lat, $
     npmodis,hmin,hmax,vmin,vmax
  
;; this code analyzes regional prognostic/assimilation phenology
;; from the pheno_analysis model and compares it with MODIS
;; phenology

  DEVICE,retain=2

  ;; Detect FPE and print code lines where FPE occurs during runtime (res)
  !EXCEPT = 2

  ;; define nodata value
  nodata = -9999.               ; missing data value used throughout the code (if not !values.f_nan)

  ;; Directories

;; the pheno analysis output directory (where the output of the
;; experiments are stored
  basedir='/Users/stockli/phenoanalysis-output/'
;  basedir='/scratch/stockli/phenoanalysis-output/'
  basedir='/glade/p/cgd/tss/people/stockli/phenoanalysis-output/'
  
  ;; MODIS data is stored here
;  datadir='/Users/stockli/phenoanalysis-datasets'
  datadir='/Users/stockli'
;  datadir='/scratch/stockli'
  datadir='/glade/p/cgd/tss/people/stockli/phenoanalysis-datasets'
;  datadir='/glade/p/cgd/tss/people/stockli'

  ;; file for site/regional information
  sitedir='../../data/sites/' 
;  sitefile = 'example.dat'
;  sitefile = 'garonna.dat'
  sitefile = 'global.dat'
;  sitefile = 'global4.dat'
;  sitefile = 'swiss.dat'
;  sitefile = 'cosmo.dat'

  ;; define time extent
  years = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012]
  years = [2008]

  
;  years = indgen(50)+1960

  ;; selection of boundary file domain (pft/topo)
;;     region = 'Colorado_Front_Range-FULL'
  region = 'Global-FULL'

  ;; selection of expemiment
;     experiments = ['Example2']
;  experiments = ['Example3']
;  experiments = ['Garonna']
     experiments = ['Global']
;     experiments = ['Swiss-Assimilation']
;     experiments = ['Swiss-Prediction']
;     experiments = ['COSMO']

  ;; select variables
;                              
;     varnames = 'PFT'
;     varnames = 'TOPO'
     varnames = ['FPAR']
;     varnames = ['LAI']
;  varnames = ['LAI','TEMP_FAC','MOIST_FAC','LIGHT_FAC'] ;,'FORCE_FAC','CHILL_FAC']
;     varnames = ['FPAR','TEMP_FAC','MOIST_FAC','LIGHT_FAC'] ;,'FORCE_FAC','CHILL_FAC']
;     varnames  = ['LAI','TMIN_AVE']
;     varnames = ['MOIST_FAC']
;     varnames = ['TEMP_FAC']
;     varnames = ['LIGHT_FAC']
;     varnames = ['TMIN_AVE']
;     varnames = ['PHOTO_AVE']
;     varnames = ['VPD_AVE']
;     varnames = ['RG_AVE']
;     varnames = ['GSI']

  ;; Specify Plotting grid independent of Model or Observation grid
  ;; either (centered on site lon/lat, also suitable for choice of many sites):
  pollon = -180.0
  pollat = 90.0
;  pollon = -170.0
;  pollat = 43.0

  ;; version of MODIS L3 tiles (3,4,5 ..)
  modis_version = 5

  ;; plotting time step (daily)
  dtplot = 24.*3600.            ; plot interval [s]

  ;; first plotting time step (0 = 1 January)
  t0 = 176

  sites0 = 1
  sites1 = sites0

;  sites0 = (7-1)*4 + 1
;  sites1 = sites0

  FOR allsites=sites0,sites1 DO BEGIN

     ;; selection of area
     sites = allsites

     range  = [0,1]             ; plot-range magnitude
     range2 = [0,1]
     
     obsvarnames = varnames

     ;; Flags
     series = 1                 ; create regional time series plot instead of map
     scatter = 0                ; create regional scatter plot instead of map
     taylor  = 0                ; create regional taylor plot instead of map
     regline = 1

     autorange = 0              ; auto-guess plot range
     psout = 0                  ; output EPS of map instead of timeseries
     psseries = 0               ; output EPS of timeseries instead of map
     simulator = 0              ; plot simulator output or not
     saveobs = 1                ; plot observations saved from model run
     plotmodels = 1             ; plot model output or not
     plotobs = 1                ; plot observations or not
     ploterror = 1              ; plot error bars on scatter plots
     plotpft = 0                ; plot PFT map
     plottopo = 0               ; plot topography map
     
     use_grid = 1               ; use grid_xxx.dat grid information (1) or above regular grid info (0)
     plottitle = 1              ; plot a title on each window
     plotdate = 1               ; plot the date on each window
     plotexperiment = 0         ; plot experiment name on each window
     ignoregaps = 1             ; plot model data when observations are missing
     regrid = 0                 ; regrid output by triangulation or simply scale to output
     plot_map = 0               ; also plot continents and rivers etc.
     hammeraitoff = 0           ; create plot in hammer aitoff projection (only useful for global plot)
     interactive = 1            ; interactively select area on the map or simply plot full area    
     makelatextable = 0         ; save correlation statistics in LaTeX table (only if scatter is selected)
     makedifference = 1         ; plot difference map
     makeabsdifference = 1      ; plot absolute difference map (only valid with makedifference set to 1 also)
     makemean = 0               ; plot annual mean + seasonal variability (mandatory for all flags below)
     makevar = 0                ; plot interannual variability and decadal trends (only if makemean is also selected)
     makeseas = 0               ; plot mean seasonal cycle
     makeannvar = 0             ; plot interannual variability
     makecontour = 0            ; plot contours of variability as second layer
     maketiming = 0             ; plot mean spring date and intereannual variability of mean spring date

     dtmodis = 8                ; modis compositing period [days]

     ;; Design
     colored = 1
     ysize_ps_default=8         ; default height of the plot [cm]
     ysize_screen_default=400   ; default height of the plot [pixels]
     xsize_screen_max=1280      ; maximum horizontal screen size [pixels]
     
     xylim = 0.075              ; xy margins
     cbsize = 0.2               ; cb x size
     aspect = 0
     legend = 1                 ; plot legend (for the time series plots)
     legenddx = 0.05
     legenddy = -0.0
     colorbar = 1               ; 0: none, 1: within window, 2: separate window

     ;; pft-specific selectors
     selpft = 34 ;; plot this pft only (if available from model output, 1..npft) 5,8,3,15

     ;; for averaged pft classes: only plot pixels (obs and model)
     ;; that satisfy these conditions
     minsumpft = 75. 
     minselpft = 0.1 * (100.-minsumpft)

     ;; elevation specific selectors (not implemented yet)
     selhgt = 0 ;; plot this hgt only (if available from model output, 1..nhgt)

     ;; check for correct plot selection
     IF (taylor+scatter+series) GT 1 THEN BEGIN
        print,'Please select either Series or Taylor or Scatter plot.'
        stop
     ENDIF

     IF makelatextable AND (~scatter AND ~taylor) THEN BEGIN
        print,'Please select SCATTER OR TAYLOR (scatter/taylor plots) in order to print a latex table with statistics'
        stop
     ENDIF

     IF psout AND psseries THEN BEGIN
        print,'Please select either psout or psmap.'
     ENDIF
     
     IF (interactive+psout+psseries) EQ 0 THEN BEGIN
        print,'Please select interactive if psout and psseries are 0.'
        stop
     ENDIF

     IF (interactive+psout) GT 1 THEN BEGIN
        print,'Interactive and psout cannot be selected simultaneously.'
        stop
     ENDIF

     IF plotpft AND plottopo THEN BEGIN
        print,'Please select either pft or topography for plotting.'
        stop
     ENDIF

     IF plotpft AND (selpft EQ 0) THEN BEGIN
        print,'Please set selpft to 1-35 when plotting pft cover.'
        stop
     ENDIF

     pftnames = [$
                'bar_all',$
                'enf_tem',$
                'enf_bor',$
                'dnf_bor',$
                'ebf_tro',$
                'ebf_tem',$
                'dbf_tro',$
                'dbf_tem',$
                'dbf_bor',$
                'ebs_all',$
                'dbs_tem',$
                'dbs_bor',$
                'c3g_arc',$
                'c3g_nar',$
                'c4g_all',$
                'cro_brl',$
                'cro_cas',$
                'cro_cot',$
                'cro_grn',$
                'cro_mze',$
                'cro_mil',$
                'cro_oil',$
                'cro_oth',$
                'cro_pot',$
                'cro_pul',$
                'cro_rap',$
                'cro_ric',$
                'cro_rye',$
                'cro_sor',$
                'cro_soy',$
                'cro_sgb',$
                'cro_sgc',$
                'cro_sun',$
                'cro_wht',$
                'wat_all']

     pftlongnames = [$
       'Bare Soil (-)',$
       'Evergreen Needleleaf Trees (Temperate)',$
       'Evergreen Needleleaf Trees (Boreal)',$
       'Deciduous Needleleaf Trees (Boreal)',$
       'Evergreen Broadleaf Trees  (Tropical)',$
       'Evergreen Broadleaf Trees  (Temperate)',$
       'Deciduous Broadleaf Trees  (Tropical)',$
       'Deciduous Broadleaf Trees  (Temperate)',$
       'Deciduous Broadleaf Trees  (Boreal)',$
       'Evergreen Broadleaf Shrubs (All)',$
       'Deciduous Broadleaf Shrubs (Temperate)',$
       'Deciduous Broadleaf Shrubs (Boreal)',$
       'C3 Grass (Arctic)',$
       'C3 Grass (Non-Arctic)',$
       'C4 Grass (All)',$
       'Crops (Barley)',$
       'Crops (Cassava)',$
       'Crops (Cotton)',$
       'Crops (Groundnuts)',$
       'Crops (Maize)',$
       'Crops (Millet)',$
       'Crops (Oilpalm)',$
       'Crops (Other)',$
       'Crops (Potato)',$
       'Crops (Pulses)',$
       'Crops (Rape)',$
       'Crops (Rice)',$
       'Crops (Rye)',$
       'Crops (Sorghum)',$
       'Crops (Soy)',$
       'Crops (Sugarbeets)',$
       'Crops (Sugarcane)',$
       'Crops (Sunflower)',$
       'Crops (Wheat)',$
       'Water (-)']

     IF plotpft OR plottopo THEN BEGIN
        IF plottopo THEN BEGIN
           plottitle = 0
           varnames = "Z"
           range = [0,2000]
        ENDIF
        IF plotpft THEN BEGIN
           plottitle = 1
           modtitle = pftlongnames[selpft-1]
           varnames = "PFT_PCT"
           range = [0,50]
        ENDIF
        plotmodels = 1
        plotobs = 0
        plotdate = 0
     ENDIF

     ;; maps and plots go here
     IF plottopo OR plotpft OR makedifference THEN BEGIN        
        mapdir=basedir+experiments[0]+'/static/'
     ENDIF ELSE BEGIN
        mapdir=basedir+experiments[0]+'/maps/'
     ENDELSE
     plotdir=basedir+experiments[0]+'/plots/'
     tabledir=basedir+experiments[0]+'/tables/' 

     ;; ---- START OF MAIN CODE -----
     
     ;; prescribe the use of model number 1 (phenology model)
     models = 1

     ;; Time Attributes
     monstr_short=['J','F','M','A','M','J','J','A','S','O','N','D']
     monstr_long=['JAN','FEB','MAR','APR','MAY','JUN', $
                  'JUL','AUG','SEP','OCT','NOV','DEC']
     monday_noleap=[31,28,31,30,31,30,31,31,30,31,30,31]
     monday_leap=[31,29,31,30,31,30,31,31,30,31,30,31]     
     monsum_noleap=[0,31,59,90,120,151,181,212,243,273,304,334,365]
     monsum_leap=[0,31,60,91,121,152,182,213,244,274,305,335,366]

     ;; always allocate 366 days per year 
     nyears=long(n_elements(years))
     dpy=366L
     ntplot = dpy * 24L * 3600L / long(dtplot)

     dtobs=!values.f_nan        ; we do not know observation time step yet
     dtmod=!values.f_nan        ; we do not know model time step yet
     dtsim=!values.f_nan        ; we do not know model time step yet

     ntmod = 0
     ntobs = 0
     ntsim = 0

     firstobs=1
     firstmodel=1

     has_obs=0
     has_model=0
     has_sim=0
     nlon=0
     nlat=0
     dlon=0.0
     dlat=0.0
     

     IF simulator THEN plotobs=1

     nexperiments = n_elements(experiments)
     
     nsites = n_elements(sites)

     nvarnames = n_elements(varnames)
     IF nexperiments NE nvarnames THEN BEGIN
        nexperiments=nvarnames
        experiments=replicate(experiments[0],nvarnames)
     ENDIF
     nmodels = n_elements(models)
     IF nmodels NE nexperiments THEN BEGIN
        nmodels=nexperiments
        models=replicate(models[0],nmodels)
     ENDIF
     nobs =  n_elements(obsvarnames)
     IF nobs NE nexperiments THEN BEGIN
        nobs=nexperiments
        obsvarnames=replicate(obsvarnames[0],nobs)
     ENDIF

     IF nsites GT 1 THEN BEGIN
        manysites = 1 
        IF nsites NE nexperiments THEN BEGIN
           nexperiments=nsites
           experiments = replicate(experiments[0],nexperiments)
        ENDIF
     ENDIF ELSE BEGIN 
        manysites = 0
        IF nsites NE nexperiments THEN BEGIN
           nsites=nexperiments
           sites = replicate(sites[0],nsites)
        ENDIF
     ENDELSE

     
     nplot=(plotobs+plotmodels)>1

     IF (sites[0] EQ sites0) THEN BEGIN
        corr_all = make_array(nmodels,sites1-sites0+2,/float)
        rmse_all = make_array(nmodels,sites1-sites0+2,/float)
        name_all = make_array(sites1-sites0+2,/string)
     ENDIF

     IF aspect EQ 0 THEN aspect = nyears<3

     ;; Time arrays
     syear0=''
     syear1=''
     reads,years[0],syear0
     reads,years[nyears-1],syear1
     syear0=strtrim(syear0,2)
     syear1=strtrim(syear1,2)
     IF nyears GT 1 THEN stime=syear0+'-'+syear1 ELSE stime=syear0

     ;; Data arrays
     modelnames=strarr(nmodels)
     modellongnames=strarr(nmodels)

     modvarnames=strarr(nmodels)
     modvarunits=strarr(nmodels)
     modvarlongnames=strarr(nmodels)

     obsvarunits=strarr(nobs)
     obsvarlongnames=strarr(nobs)
     
     ;; open Site file

     stemp=''
     siteshort=''
     sitename=''
     sitelat=0.
     sitelon=0.

     siteshorts=strarr(nsites)
     sitenames=strarr(nsites)
     sitelonmins=dblarr(nsites)
     sitelonmaxs=dblarr(nsites)
     sitelatmins=dblarr(nsites)
     sitelatmaxs=dblarr(nsites)
     sitedlons=dblarr(nsites)
     sitedlats=dblarr(nsites)
     
     FOR s = 0,nsites-1 DO BEGIN
        openr,1,sitedir+"/"+sitefile
        readf,1,stemp
        readf,1,stemp

        WHILE (NOT EOF(1)) DO BEGIN
           readf,1,stemp
           
           arr=strsplit(stemp,',',/extract)
           snum = arr[0]
           IF snum EQ sites[s] THEN BEGIN
              siteshorts[s]=strtrim(arr[1],2)
              sitenames[s]=strtrim(arr[2],2)
              sitelonmins[s]=double(strtrim(arr[3],2)+"d0")
              sitelonmaxs[s]=double(strtrim(arr[4],2)+"d0")
              sitelatmins[s]=double(strtrim(arr[5],2)+"d0")
              sitelatmaxs[s]=double(strtrim(arr[6],2)+"d0")
              IF strtrim(arr[7],2) EQ "NA" THEN sitedlons[s] = 0.05d0 ELSE $
                 sitedlons[s]=double(strtrim(arr[7],2)+"d0")
              IF strtrim(arr[8],2) EQ "NA" THEN sitedlats[s] = 0.05d0 ELSE $
                 sitedlats[s]=double(strtrim(arr[8],2)+"d0")
           ENDIF
        ENDWHILE
        close,1

        IF sitenames[s] EQ '' THEN BEGIN
           print,'Site Number ',sites[s],' not found in ',sitedir,'/',sitefile
           stop
        ENDIF
     ENDFOR

     ;; geographic bounds for plotting and observation array       
     pollon = 180.d0
     pollat = 90.d0
     lonmin = sitelonmins[0] 
     lonmax = sitelonmaxs[0]
     latmin = sitelatmins[0]
     latmax = sitelatmaxs[0]
     dlon = sitedlons[0]
     dlat = sitedlats[0]

     print,'Grid detected with (lonmin/lonmax/latmin/latmax) : ',lonmin,lonmax,latmin,latmax
     ;; geographic bounds for plotting and observation array (single
     ;; site, known lon/lat boundaries)   
     nlon = round((lonmax-lonmin)/dlon)
     nlat = round((latmax-latmin)/dlat)
     
     print,'Number of lon/lat points : ',nlon,nlat
     
     ;; create lon/lat arrays (centered on grid point)
;        lons = (findgen(nlon)*dlon + lonmin + 0.5*dlon) # replicate(1.0,nlat) ; W->E
;        lats = replicate(1.0,nlon) # (findgen(nlat)*dlat + latmin + 0.5*dlat) ; S->N
     lons = dindgen(nlon)*dlon + lonmin + 0.5d0*dlon ; W->E
     lats = dindgen(nlat)*dlat + latmin + 0.5d0*dlat ; S->N
     
     ;; check whether we have a rotated pole (needs special treatment for map grid)
     IF (abs(pollat-90.d0) lt 1.d-3) THEN BEGIN
        rotatedpole = 0
     ENDIF ELSE BEGIN
        rotatedpole = 1
        print,'Rotated Pole detected with North Pole at (lon/lat) : ',pollon,pollat
     ENDELSE
     
     
     IF plotmodels THEN BEGIN
        ;; read model data

        FOR m=0,nmodels-1 DO BEGIN

           FOR y=0L,nyears-1 DO BEGIN
              year=years[y]

              pftave = 1        ; average pfts
              hgtave = 1        ; average height classes

              read_pheno_regional,basedir,1,experiments[m],sitenames[m],year,varnames[m],dtmod, $
                                  modvartemp,modvarerrortemp,modvarnamestemp,modvarunitstemp,modvarlongnamestemp, $
                                  modlons, modlats, modelnamestemp,modellongnamestemp, nodata, has_model, $
                                  nlonmod, nlatmod, pftave, hgtave, selpft, selhgt, pft, topo, $
                                  lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax

              IF has_model THEN BEGIN              
                 IF finite(dtmod) THEN BEGIN
                    IF firstmodel THEN BEGIN
                       ;; first read: define output grid
                       
                       ;; do not use topography masking right now
                       topomask = bytarr(nlon,nlat)
                       topomask[*] = 1B

                       dims = (size(pft,/dimensions))
                       IF n_elements(dims) EQ 3 THEN BEGIN
                          npft = dims[2]
                          pft_ave = fltarr(npft)
                          
                          landidx = where(pft[*,*,npft-1] LT 50.,count)

                          IF count EQ 0 THEN BEGIN
                             print,'Only water pixels selected. Please choose another grid area'
                             stop
                          ENDIF

                          FOR p = 0,npft-1 DO pft_ave[p] = mean((pft[*,*,p])[landidx])
                          
                          ;; screen by minimum pft appearance (do never
                          ;; select water pixels (npft-1) in observations
                          pftidx = where(pft_ave[0:npft-2] GT minselpft,count)
                          
                          ;; overrride for plotting: screen only water pixels (npft-1)
                          count2=npft-1
                          pftidx2 = indgen(count2)
                          
                          IF count2 EQ 1 THEN pftmask = reform(pft[*,*,pftidx2]) GE minsumpft ELSE $
                             pftmask = total(pft[*,*,pftidx2],3) GE minsumpft
                          
                          ;; set up plot mask based on pft/topography
                          ;; screening (only pft employed right now)
                          gridmask = pftmask AND topomask

                       ENDIF ELSE BEGIN
                          gridmask = make_array(nlon,nlat,/byte,value=1B)
                          npft = 1
                       ENDELSE

                       noplotidx = where(gridmask EQ 0,noplotcount)

                       ;; allocate model arrays
                       ntmod=dpy*24L*3600L/long(dtmod)     
                       modvar=make_array(nlon, nlat,ntmod,nyears,nmodels,/float,value=!VALUES.F_NAN)
                       IF ploterror THEN modvarerror=make_array(nlon, nlat,ntmod,nyears,nmodels,/float,value=!VALUES.F_NAN)
                       
                       ;; not new read anymore, proceed with standard read
                       firstmodel=0
                    ENDIF

                    modvarnames[m]=modvarnamestemp
                    modvarunits[m]=modvarunitstemp
                    modvarlongnames[m]=modvarlongnamestemp

                    modelnames[m]=modelnamestemp
                    modellongnames[m]=modellongnamestemp

                    badmod = where(modvartemp EQ nodata,badcount)
                    IF badcount GT 0 THEN modvartemp[badmod] = !values.f_nan

                    ;; rebin model grid to plotting grid
                    leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
                    IF leapyear THEN ntend = dpy*24L*3600L/long(dtmod) ELSE ntend = (dpy-1)*24L*3600L/long(dtmod)
                    FOR t = 0,ntend-1 DO BEGIN 
                       modvar[*,*,t,y,m]=rebin(modvartemp[*,*,t],nlon,nlat)
                       IF ploterror THEN modvarerror[*,*,t,y,m]=rebin(modvarerrortemp[*,*,t],nlon,nlat)
                    ENDFOR

                    modvartemp = 0
                    modvarerrortemp = 0

                 ENDIF ;; has data
              ENDIF    ;; has model

           ENDFOR ;; years

           ;; select PFT
           IF plotpft THEN BEGIN
;              modvarlongnames[0] = strupcase(pftnames[selpft-1])+" Plant Functional Type"
           ENDIF

           IF noplotcount GT 0 THEN BEGIN

              FOR ly=long64(0),nyears-1 DO BEGIN
                 FOR ld=long64(0),ntmod-1 DO BEGIN
                    idx = long64(noplotidx) + $
                          long64(nlon)*long64(nlat)*ld + $
                          long64(nlon)*long64(nlat)*long64(ntmod)*ly + $
                          long64(nlon)*long64(nlat)*long64(ntmod)*long64(nyears)*long64(m)
                    modvar[idx] = !values.f_nan
                 ENDFOR
              ENDFOR

              IF ploterror THEN BEGIN
                 FOR ly=long64(0),nyears-1 DO BEGIN
                    FOR ld=long64(0),ntmod-1 DO BEGIN
                       idx = long64(noplotidx) + $
                             long64(nlon)*long64(nlat)*ld + $
                             long64(nlon)*long64(nlat)*long64(ntmod)*ly + $
                             long64(nlon)*long64(nlat)*long64(ntmod)*long64(nyears)*long64(m)
                       modvarerror[idx] = !values.f_nan
                    ENDFOR
                 ENDFOR
              ENDIF

           ENDIF

        ENDFOR ;; models

     ENDIF ;; plotmodels

     IF plotobs THEN BEGIN
        ;; read observation data

        FOR o=0,nobs-1 DO BEGIN
           
           FOR y=0L,nyears-1 DO BEGIN
              year=years[y]

              IF simulator OR saveobs THEN BEGIN
                 pftave = 1
                 hgtave = 1
                 read_pheno_regional,basedir,2+saveobs,experiments[o],sitenames[o],year,varnames[o],dtsim, $
                                     simvartemp,simvarerrortemp,simvarnamestemp,simvarunitstemp,simvarlongnamestemp, $
                                     simlons, simlats, simnamestemp,simlongnamestemp, nodata, has_sim, $
                                     nlonsim, nlatsim, pftave, hgtave, selpft, selhgt, pft, topo, $
                                     lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax

                 IF ~keyword_set(noplotids) THEN BEGIN
                    ;; do not use topography masking right now
                    topomask = bytarr(nlon,nlat)
                    topomask[*] = 1B

                    dims = (size(pft,/dimensions))
                    IF n_elements(dims) EQ 3 THEN BEGIN
                       npft = dims[2]
                       pft_ave = fltarr(npft)
                       
                       landidx = where(pft[*,*,npft-1] LT 50.,count)
                       
                       IF count EQ 0 THEN BEGIN
                          print,'Only water pixels selected. Please choose another grid area'
                          stop
                       ENDIF
                       
                       FOR p = 0,npft-1 DO pft_ave[p] = mean((pft[*,*,p])[landidx])
                       
                       ;; screen by minimum pft appearance (do never
                       ;; select water pixels (npft-1) in observations
                       pftidx = where(pft_ave[0:npft-2] GT minselpft,count)
                       
                       ;; overrride for plotting: screen only water pixels (npft-1)
                       count2=npft-1
                       pftidx2 = indgen(count2)
                       
                       IF count2 EQ 1 THEN pftmask = reform(pft[*,*,pftidx2]) GE minsumpft ELSE $
                          pftmask = total(pft[*,*,pftidx2],3) GE minsumpft
                       
                       ;; set up plot mask based on pft/topography
                       ;; screening (only pft employed right now)
                       gridmask = pftmask AND topomask
                       
                    ENDIF ELSE BEGIN
                       gridmask = make_array(nlon,nlat,/byte,value=1B)
                       npft = 1
                    ENDELSE
                    
                    noplotidx = where(gridmask EQ 0,noplotcount)

                 ENDIF

                 IF has_sim THEN BEGIN

                    has_obs = 1
                    IF firstobs THEN BEGIN
                       ntobs=dpy*24L*3600L/long(dtsim)     
                       obsvar=make_array(nlon, nlat,ntobs,nyears,nobs,/float,value=!VALUES.F_NAN)
                       IF ploterror THEN obsvarerror=make_array(2,nlon, nlat,ntobs,nyears,nobs,/float,value=!VALUES.F_NAN)
                       obsdays=intarr(ntobs,nyears,nobs)
                       firstsim=0
                       dtobs = dtsim
                       firstobs = 0
                    ENDIF

                    leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
                    IF leapyear THEN ntend = dpy*24L*3600L/long(dtsim) ELSE ntend = (dpy-1)*24L*3600L/long(dtsim)
                    FOR t = 0,ntend-1 DO BEGIN 
                       obsvar[*,*,t,y,o]=rebin(simvartemp[*,*,t],nlon,nlat)
                       IF ploterror THEN BEGIN
                          obsvarerror[0,*,*,t,y,o]=obsvar[*,*,t,y,o]+rebin(simvarerrortemp[*,*,t],nlon,nlat)
                          obsvarerror[1,*,*,t,y,o]=obsvar[*,*,t,y,o]-rebin(simvarerrortemp[*,*,t],nlon,nlat)
                       ENDIF
                    ENDFOR

                    obsdays[*,*,o] = indgen(ntobs*nyears) + 1                    
                    
                    obsvarlongnames[o]=simvarlongnamestemp
                    obsvarunits[o]=simvarunitstemp
                    
                 ENDIF

              ENDIF ELSE BEGIN
                 
                 IF (obsvarnames[o] EQ 'LAI') OR (obsvarnames[o] EQ 'FPAR') THEN BEGIN

                    ;; read and composite L3 MOD15A2 tiles
                    regrid_modis,datadir, year, lonmin, lonmax, latmin, latmax, dlon, dlat, pollon, pollat, $
                                 obsvartemp, obsvarerrortemp, obsvarnames[o], obsvarunitstemp, obsvarlongnamestemp, $
                                 modis_version
                    
                    IF firstobs THEN BEGIN
                       dims = size(obsvartemp,/dimensions)
                       ntobs=dpy     
                       dtobs=3600L*24L
                       obsvar=make_array(nlon,nlat,ntobs,nyears,nobs,/float,value=!values.f_nan)
                       obsvariance=make_array(nlon,nlat,ntobs,nyears,nobs,/float,value=!values.f_nan)
                       obsvarabsdev=make_array(nlon,nlat,ntobs,nyears,nobs,/float,value=!values.f_nan)
                       IF ploterror THEN obsvarerror=make_array(2,nlon,nlat,ntobs,nyears,nobs,/float,value=!values.f_nan)
                       obsdays=intarr(ntobs,nyears,nobs)  
                       obsvarlongnames=strarr(nobs)
                       obsvarunits=strarr(nobs)
                       firstobs = 0
                    ENDIF
                    
                    leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
                    IF leapyear THEN ntend = dpy ELSE ntend = dpy-1
                    FOR t = 0,ntend-1 DO BEGIN 
                       obsvar[*,*,t,y,o]=rebin(obsvartemp[*,*,t],nlon,nlat)
                       IF ploterror THEN BEGIN
                          obsvarerror[0,*,*,t,y,o]=rebin(reform(obsvarerrortemp[0,*,*,t]),nlon,nlat)
                          obsvarerror[1,*,*,t,y,o]=rebin(reform(obsvarerrortemp[1,*,*,t]),nlon,nlat)
                       ENDIF
                    ENDFOR

                    obsdays[*,*,o] = indgen(ntobs*nyears) + 1                    
                    obsvarlongnames[o] = obsvarlongnamestemp
                    obsvarunits[o]=obsvarunitstemp

                 ENDIF

              ENDELSE
              
           ENDFOR ;; years

        ENDFOR ;; sites

        ;; simplify obsdays
        obsdays = obsdays[*,0,0]
        
        IF noplotcount GT 0L THEN BEGIN

           FOR lo=long64(0),nobs-1 DO BEGIN
              FOR ly=long64(0),nyears-1 DO BEGIN
                 FOR ld=long64(0),ntobs-1 DO BEGIN
                    idx = long64(noplotidx) + $
                          long64(nlon)*long64(nlat)*ld + $
                          long64(nlon)*long64(nlat)*long64(ntobs)*ly + $
                          long64(nlon)*long64(nlat)*long64(ntobs)*long64(nyears)*lo
                    obsvar[idx] = !values.f_nan                    
                 ENDFOR
              ENDFOR
           ENDFOR

        ENDIF
        
     ENDIF ;; plotobs

     axisindex=bytarr(nmodels+nobs)
     secondaxis=0
     IF plotobs THEN FOR o=0,nobs-1 DO IF (obsvarunits[o] NE obsvarunits[0]) THEN axisindex[o]=1
     IF plotmodels AND plotobs THEN FOR m=0,nmodels-1 DO IF modvarunits[m] NE obsvarunits[0] THEN axisindex[m+nobs]=1
     IF plotmodels AND (plotobs EQ 0) THEN FOR m=0,nmodels-1 DO IF modvarunits[m] NE modvarunits[0] THEN axisindex[m+nobs]=1
     IF (max(axisindex) GT 0) THEN secondaxis=1 
     
     IF firstobs THEN plotobs = 0
     IF firstmodel THEN plotmodel = 0
     IF plotpft OR plottopo THEN plotmodels = 1
     
     ;; auto-guess the plot range and set y axis naming
     IF autorange THEN BEGIN
        minmod=make_array(nmodels,/float,value=!values.f_nan)
        maxmod=make_array(nmodels,/float,value=!values.f_nan)
        minobs=make_array(nobs,/float,value=!values.f_nan)
        maxobs=make_array(nobs,/float,value=!values.f_nan)
        IF plotobs THEN BEGIN
           FOR o=0,nobs-1 DO BEGIN
              automin=min(obsvar[*,*,*,*,o],max=automax,/NAN)             
              minobs[o]=automin
              maxobs[o]=automax
           ENDFOR
        ENDIF
        IF plotmodels THEN BEGIN
           FOR m=0,nmodels-1 DO BEGIN
              automin=min(modvar[*,*,*,*,m],max=automax,/NAN)
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

     ;; calculate difference btw obs and mod if needed
     IF makedifference THEN BEGIN
        nplot = 1

        IF plotobs THEN BEGIN

           plotmodels = 0
           diffvar=make_array(nlon,nlat,ntobs,nyears,/float)
           idx=where(finite(reform(modvar[*,*,obsdays[*]-1,*,0]),/nan) OR finite(obsvar[*,*,*,*,0],/nan),count)
           diffvar[*] = reform(modvar[*,*,obsdays[*]-1,*,0]) - obsvar[*,*,*,*,0]
           IF makeabsdifference THEN diffvar=abs(diffvar)
           IF count GT 0 THEN diffvar[idx] = !values.f_nan
           
           IF nyears EQ 1 THEN diffcount = total(finite(diffvar),3) ELSE $ 
              diffcount = total(total(finite(diffvar),4),3)
           IF nyears EQ 1 THEN diffmean = total(diffvar,3,/nan) / float(diffcount>1) ELSE $
              diffmean = total(total(diffvar,4,/nan),3,/nan) / float(diffcount>1)
           idx = where(diffcount EQ 0,count)
           IF count GT 0 THEN diffmean[idx] = !values.f_nan

           obsvar[*] = !values.f_nan
           obsvarerror[*] = !values.f_nan
           obsvar[*,*,*,*,0] = diffvar

        ENDIF ELSE BEGIN

           diffvar = make_array(nlon,nlat,ntmod*nyears,/float)
           diffvar[*] = modvarerror[*,*,*,*,0]
           diffcount =  total(finite(diffvar),3)
           difftotal = total(diffvar,3,/nan)
           diffmean = difftotal / float(diffcount>1)
           idx = where(diffcount EQ 0,count)
           IF count GT 0 THEN diffmean[idx] = !values.f_nan

           modvar[*] = !values.f_nan
           modvarerror[*] = !values.f_nan

           modvar[*,*,*,0,0] = diffvar

        ENDELSE

        ;; adjust ranges
        IF makeabsdifference THEN BEGIN
           range[0] = 0.0
           range[1] = range[1]*0.5
        ENDIF ELSE BEGIN
           range[0] = -range[1]*0.5
           range[1] = range[1]*0.5
        ENDELSE

        diffwgts = replicate(1.0,nlon) # cos(lats*!dtor)

        IF count GT 0 THEN diffmean[idx] = 0.
        diffall = total(diffmean * diffwgts,/nan) / total(diffwgts * (diffcount GT 0))

        IF makeabsdifference THEN $
           print,'Mean Absolute Difference: ',diffall ELSE $
              print,'Mean Difference: ',diffall

     ENDIF


     ;; calculate observation statistics
     IF plotobs THEN BEGIN
        FOR o=0,nobs-1 DO BEGIN

           ;; calculate annual mean and seasonal variability
           IF makemean THEN BEGIN
              ;; calculate mean seasonal cycle
              IF nyears GT 1 THEN BEGIN
                 seas = total(reform(obsvar[*,*,*,*,o]),4,/nan)
                 seascount = total(finite(reform(obsvar[*,*,*,*,o])),4)
                 idx = where(seascount EQ 0,count)
                 seas = seas / float(seascount>1)
                 IF count GT 0 THEN seas[idx] = !values.f_nan                 
              ENDIF ELSE seas = reform(obsvar[*,*,*,0,o])

              IF makeseas THEN BEGIN

                 obsvar[*,*,*,0,o] = seas
                 IF nyears GT 1 THEN obsvar[*,*,*,1:nyears-1,o] = !values.f_nan

                 IF ploterror THEN BEGIN
                    FOR e=0,1 DO BEGIN
                       IF nyears GT 1 THEN BEGIN
                          seas = total(reform(obsvarerror[e,*,*,*,*,o]),4,/nan)
                          seascount = total(finite(reform(obsvarerror[e,*,*,*,*,o])),4)
                          idx = where(seascount EQ 0,count)
                          seas = seas / float(seascount>1)
                          IF count GT 0 THEN seas[idx] = !values.f_nan                 
                       ENDIF ELSE seas = reform(obsvarerror[e,*,*,*,0,o])

                       obsvarerror[e,*,*,*,0,o] = seas
                       IF nyears GT 1 THEN obsvarerror[e,*,*,*,1:nyears-1,o] = !values.f_nan

                    ENDFOR
                 ENDIF

              ENDIF ELSE BEGIN

                 ;; calculate seasonal mean
                 seasmean =  total(seas,3,/nan)
                 seascount = total(finite(seas),3)
                 idx = where(seascount EQ 0,count)
                 seasmean = seasmean / float(seascount>1)
                 IF count GT 0 THEN seasmean[idx] = !values.f_nan

;                 obsvar[*,*,*,*,o] = !values.f_nan
;                 obsvar[*,*,ntobs/2,0,o] = seasmean

                 ;; calculate seasonal variance
                 seasmean2 = make_array(nlon,nlat,ntobs,/float)
                 FOR i=0,ntobs-1 DO seasmean2[*,*,i] = seasmean
                 seasvar = total((seas-seasmean2)^2,3,/nan)
                 seascount = total(finite(seas),3)
                 idx = where(seascount LE 1,count)
                 seasvar = seasvar / ((seascount-1)>1)
                 IF count GT 0 THEN seasvar[idx] = !values.f_nan

                 ;; calculate seasonal amplitude
                 seasamp = max(seas,dimension=3,/nan) - min(seas,dimension=3,/nan)

;                 obsvar[*,*,*,*,o] = !values.f_nan
;                 obsvar[*,*,ntobs/2,0,o] = seasamp

                 IF makevar AND (NOT maketiming) THEN BEGIN
                    ;; calculate annual mean
                    ann =  total(reform(obsvar[*,*,*,*,o]),3,/nan)
                    anncount = total(finite(reform(obsvar[*,*,*,*,o])),3)
                    idx = where(anncount EQ 0,count)
                    ann = ann / float(anncount>1)
                    IF count GT 0 THEN ann[idx] = !values.f_nan

                    ;; calculate all-mean
                    IF nyears GT 1 THEN BEGIN
                       annmean = total(ann,3,/nan)
                       anncount = total(finite(ann),3)
                       idx = where(anncount EQ 0,count)
                       annmean = annmean / float(anncount>1)
                       IF count GT 0 THEN annmean[idx] = !values.f_nan
                    ENDIF ELSE annmean = ann

                    IF (nyears GT 1) THEN BEGIN
                       ;; calculate interannual variance
                       ann2 = make_array(nlon,nlat,nyears,/float)
                       FOR i=0,nyears-1 DO ann2[*,*,i] = annmean
                       annvar = total((ann-ann2)^2,3,/nan)
                       anncount = total(finite(ann),3)
                       idx = where(anncount LE 1,count)
                       annvar = annvar / ((anncount-1)>1)
                       IF count GT 0 THEN annvar[idx] = !values.f_nan
                       
                       ;; calculate trend
;;                        anntrend = make_array(nlon,nlat,/float,value=!values.f_nan)
;;                        idx = where(finite(annmean),count)
;;                        IF count GT 0 THEN BEGIN
;;                           xval = findgen(nyears)
;;                           ind = array_indices(annmean,idx)
;;                           FOR i=0,count-1 DO BEGIN
;;                              yval = reform(ann[ind[0,i],ind[1,i],*])
;;                              IF stddev(yval) GT 1.e-3 THEN $
;;                                 anntrend[ind[0,i],ind[1,i]] = regress(xval,yval) ELSE $
;;                                    anntrend[ind[0,i],ind[1,i]] = !values.f_nan
;;                           ENDFOR
;;                        ENDIF
                       

                       IF makeannvar THEN BEGIN
                          ;; calculate interannual variability
                          
                          obsvar[*,*,*,*,o] = !values.f_nan
                          obsvar[*,*,ntobs/2,*,o] = ann-ann2

                       ENDIF ELSE BEGIN
                          obsvar[*,*,*,*,o] = !values.f_nan
                          obsvar[*,*,ntobs/2,*,o] = ann
;;                          obsvar[*,*,ntobs/2,*,o] = anntrend * 10.0
                       ENDELSE

                    ENDIF

                 ENDIF

              ENDELSE

              IF maketiming THEN BEGIN
                 ;; calculate annual spring date where the seasonal
                 ;; variance is large enough              
                 badidx = where((seasamp LT 1.0) OR finite(seasmean,/nan),badcount)

                 spring = make_array(nlon,nlat,nyears,/float,value=!values.f_nan)
                 FOR i=0,nyears-1 DO BEGIN
                    nshift = 5
                    diff = reform(obsvar[*,*,*,i,o]) - shift(reform(obsvar[*,*,*,i,o]),0,0,nshift)
                    diff[*,*,0:nshift+1] = !values.f_nan
                    a = max(diff,idx,dimension=3,/nan)
                    idx3d = float(reform((array_indices(diff,idx))[2,*]))
                    IF badcount GT 0 THEN idx3d[badidx] = !values.f_nan
                    spring[*,*,i] = idx3d                 
                 ENDFOR

                 ;; mean spring date
                 IF nyears GT 1 THEN BEGIN
                    springmean = total(spring,3,/nan)
                    springcount = total(finite(spring),3)
                    idx = where(springcount EQ 0,count)
                    springmean = springmean / float(springcount>1)
                    IF count GT 0 THEN springmean[idx] = !values.f_nan
                 ENDIF ELSE springmean = spring

                 IF makevar AND (nyears GT 1) THEN BEGIN
                    ;; calculate spring date interannual variance
                    spring2 = make_array(nlon,nlat,nyears,/float)
                    FOR i=0,nyears-1 DO spring2[*,*,i] = springmean
                    springvar = total((spring-spring2)^2,3,/nan)
                    springcount = total(finite(spring),3)
                    idx = where(springcount LE 1,count)
                    springvar = springvar / ((springcount-1)>1)
                    IF count GT 0 THEN springvar[idx] = !values.f_nan
                    
                    ;; calculate spring date trend
                    springtrend = make_array(nlon,nlat,/float,value=!values.f_nan)
                    idx = where(finite(springmean),count)
                    IF count GT 0 THEN BEGIN
                       xval = findgen(nyears)
                       ind = array_indices(springmean,idx)
                       FOR i=0,count-1 DO BEGIN
                          yval = reform(spring[ind[0,i],ind[1,i],*])
                          IF stddev(yval) gt 1.e-3 THEN $
                             springtrend[ind[0,i],ind[1,i]] = regress(xval,yval) ELSE $
                                springtrend[ind[0,i],ind[1,i]] = !values.f_nan
                       ENDFOR
                    ENDIF

;                    obsvar[*,*,*,*,o] = !values.f_nan
;                    obsvar[*,*,0,0,o] = sqrt(springvar)
;                    obsvar[*,*,1,0,o] = springtrend * 10.0

                 ENDIF

              ENDIF

           ENDIF

        ENDFOR
     ENDIF ;; plotobs

     ;; calculate Model statistics
     IF plotmodels THEN BEGIN
        FOR m=0,nmodels-1 DO BEGIN

           ;; calculate annual mean and seasonal variability
           IF makemean THEN BEGIN
              ;; calculate mean seasonal cycle
              IF nyears GT 1 THEN BEGIN
                 seas = total(reform(modvar[*,*,*,*,m]),4,/nan)
                 seascount = total(finite(reform(modvar[*,*,*,*,m])),4)
                 idx = where(seascount EQ 0,count)
                 seas = seas / float(seascount>1)
                 IF count GT 0 THEN seas[idx] = !values.f_nan                 
              ENDIF ELSE seas = reform(modvar[*,*,*,0,m])

              IF makeseas THEN BEGIN

                 IF ploterror THEN BEGIN
                    IF nyears GT 1 THEN BEGIN
                       seas = total(reform(modvarerror[*,*,*,*,m]),4,/nan)
                       seascount = total(finite(reform(modvarerror[*,*,*,*,m])),4)
                       idx = where(seascount EQ 0,count)
                       seas = seas / float(seascount>1)
                       IF count GT 0 THEN seas[idx] = !values.f_nan                 
                    ENDIF ELSE seas = reform(modvarerror[*,*,*,0,m])
                    
                    modvarerror[*,*,*,0,m] = seas
                    IF nyears GT 1 THEN modvarerror[*,*,*,1:nyears-1,m] = !values.f_nan
                    
                 ENDIF

              ENDIF ELSE BEGIN

                 ;; calculate seasonal mean
                 seasmean =  total(seas,3,/nan)
                 seascount = total(finite(seas),3)
                 idx = where(seascount EQ 0,count)
                 seasmean = seasmean / float(seascount>1)
                 IF count GT 0 THEN seasmean[idx] = !values.f_nan

                 ;; calculate seasonal variance
                 seasmean2 = make_array(nlon,nlat,ntmod,/float)
                 FOR i=0,ntmod-1 DO seasmean2[*,*,i] = seasmean
                 seasvar = total((seas-seasmean2)^2,3,/nan)
                 seascount = total(finite(seas),3)
                 idx = where(seascount LE 1,count)
                 seasvar = seasvar / ((seascount-1)>1)
                 IF count GT 0 THEN seasvar[idx] = !values.f_nan

                 ;; calculate seasonal amplitude
                 seasamp = max(seas,dimension=3,/nan) - min(seas,dimension=3,/nan)

                 IF makevar AND (NOT maketiming) THEN BEGIN
                    ;; calculate annual mean
                    ann =  total(reform(modvar[*,*,*,*,m]),3,/nan)
                    anncount = total(finite(reform(modvar[*,*,*,*,m])),3)
                    idx = where(anncount EQ 0,count)
                    ann = ann / float(anncount>1)
                    IF count GT 0 THEN ann[idx] = !values.f_nan

                    ;; calculate all-mean
                    IF nyears GT 1 THEN BEGIN
                       annmean = total(ann,3,/nan)
                       anncount = total(finite(ann),3)
                       idx = where(anncount EQ 0,count)
                       annmean = annmean / float(anncount>1)
                       IF count GT 0 THEN annmean[idx] = !values.f_nan
                    ENDIF ELSE annmean = ann

                    IF (nyears GT 1) THEN BEGIN
                       ;; calculate interannual variance
                       ann2 = make_array(nlon,nlat,nyears,/float)
                       FOR i=0,nyears-1 DO ann2[*,*,i] = annmean
                       annvar = total((ann-ann2)^2,3,/nan)
                       anncount = total(finite(ann),3)
                       idx = where(anncount LE 1,count)
                       annvar = annvar / ((anncount-1)>1)
                       IF count GT 0 THEN annvar[idx] = !values.f_nan
                       
                       ;; calculate trend
;;                        anntrend = make_array(nlon,nlat,/float,value=!values.f_nan)
;;                        idx = where(finite(annmean),count)
;;                        IF count GT 0 THEN BEGIN
;;                           xval = findgen(nyears)
;;                           ind = array_indices(annmean,idx)
;;                           FOR i=0,count-1 DO BEGIN
;;                              yval = reform(ann[ind[0,i],ind[1,i],*])
;;                              IF stddev(yval) ge 1.e-3 THEN $
;;                                 anntrend[ind[0,i],ind[1,i]] = regress(xval,yval) $
;;                              ELSE anntrand[ind[0,i],ind[1,i]] = !values.f_nan
;;                           ENDFOR
;;                        ENDIF
                       

                       IF makeannvar THEN BEGIN
                          ;; calculate interannual variability
                          
                          modvar[*,*,*,*,m] = !values.f_nan
                          modvar[*,*,obsdays[ntobs/2],*,m] = ann-ann2

                       ENDIF ELSE BEGIN
                          modvar[*,*,*,*,m] = !values.f_nan
                          modvar[*,*,obsdays[ntobs/2],*,m] = ann
;                          modvar[*,*,obsdays[ntobs/2],*,m] = anntrend * 10.0
                       ENDELSE

                    ENDIF

                 ENDIF

              ENDELSE

              IF maketiming THEN BEGIN
                 ;; calculate annual spring date where the seasonal
                 ;; variance is large enough              
                 badidx = where((seasamp LT 1.0) OR finite(seasmean,/nan),badcount)

                 spring = make_array(nlon,nlat,nyears,/float,value=!values.f_nan)
                 FOR i=0,nyears-1 DO BEGIN
                    nshift = 5
                    diff = reform(modvar[*,*,*,i,m]) - shift(reform(modvar[*,*,*,i,m]),0,0,nshift)
                    diff[*,*,0:nshift+1] = !values.f_nan
                    a = max(diff,idx,dimension=3,/nan)
                    idx3d = float(reform((array_indices(diff,idx))[2,*]))
                    IF badcount GT 0 THEN idx3d[badidx] = !values.f_nan
                    spring[*,*,i] = idx3d                 
                 ENDFOR

                 ;; mean spring date
                 IF nyears GT 1 THEN BEGIN
                    springmean = total(spring,3,/nan)
                    springcount = total(finite(spring),3)
                    idx = where(springcount EQ 0,count)
                    springmean = springmean / float(springcount>1)
                    IF count GT 0 THEN springmean[idx] = !values.f_nan
                 ENDIF ELSE springmean = spring

                 modvar[*,*,*,*,m] = !values.f_nan
                 modvar[*,*,0,0,m] = springmean/365.*12.
                 modvar[*,*,1,0,m] = !values.f_nan

                 IF makevar AND (nyears GT 1) THEN BEGIN
                    ;; calculate spring date interannual variance
                    spring2 = make_array(nlon,nlat,nyears,/float)
                    FOR i=0,nyears-1 DO spring2[*,*,i] = springmean
                    springvar = total((spring-spring2)^2,3,/nan)
                    springcount = total(finite(spring),3)
                    idx = where(springcount LE 1,count)
                    springvar = springvar / ((springcount-1)>1)
                    IF count GT 0 THEN springvar[idx] = !values.f_nan
                    
                    ;; calculate spring date trend
                    springtrend = make_array(nlon,nlat,/float,value=!values.f_nan)
                    idx = where(finite(springmean),count)
                    IF count GT 0 THEN BEGIN
                       xval = findgen(nyears)
                       ind = array_indices(springmean,idx)
                       FOR i=0,count-1 DO BEGIN
                          yval = reform(spring[ind[0,i],ind[1,i],*])
                          springtrend[ind[0,i],ind[1,i]] = regress(xval,yval)
                       ENDFOR
                    ENDIF

                    modvar[*,*,*,*,m] = !values.f_nan
                    modvar[*,*,0,0,m] = sqrt(springvar)
                    modvar[*,*,1,0,m] = springtrend * 10.0

                 ENDIF

              ENDIF

           ENDIF

        ENDFOR
     ENDIF ;; plotmodels

     IF makeseas THEN nyears = 1
     
     ;; define first time step (year/day)
     t = t0
     y = 0
     
     ;; static maps are stored in first time step
     IF plotpft OR plottopo THEN t=0

     ;; play forward or backward in time
     first = 1
     keyname = ''
     WHILE (keyname NE 'x') DO BEGIN

        ;; create filename

        filename = sitenames[0]+'.'+varnames[0]+'.'+experiments[0]
        IF plotpft THEN filename += '.'+pftnames[selpft-1]
        IF makedifference THEN BEGIN
           IF makeabsdifference THEN BEGIN
              filename += '.mad'
           ENDIF ELSE BEGIN
              filename += '.mbe'
           ENDELSE
        ENDIF
        IF makemean THEN BEGIN
           IF maketiming THEN BEGIN
              IF makevar THEN BEGIN
                 filename += '.datevar'
              ENDIF ELSE BEGIN
                 filename += '.datemean'
              ENDELSE
           ENDIF ELSE BEGIN
              IF makevar THEN BEGIN
                 IF makeannvar THEN BEGIN
                    filename += '.annvar'
                 ENDIF ELSE BEGIN
                    filename += '.var'
                 ENDELSE
              ENDIF ELSE BEGIN
                 IF makeseas THEN BEGIN
                    filename += '.seas'
                 ENDIF ELSE BEGIN
                    filename += '.mean'
                 ENDELSE
              ENDELSE
           ENDELSE
        ENDIF

        IF ~(plottopo OR plotpft) THEN BEGIN
           IF scatter THEN filename += '.scatter'
           IF taylor THEN filename += '.taylor'
           IF series THEN filename += '.series'
        ENDIF

        ;; create time for plot
        year = years[y]
        
        leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
        IF leapyear THEN tend = dpy*24L*3600L/long(dtplot) ELSE tend = (dpy-1)*24L*3600L/long(dtplot)
        
        IF plotobs THEN temp = min(abs((t/dtmodis)*dtmodis+1-obsdays),obst)

        month=1
        IF leapyear THEN BEGIN
           WHILE monsum_leap[month] LT (t+1) DO month=month+1
           smonth = monstr_long[month-1]
           day = t+1 - monsum_leap[month-1]
        ENDIF ELSE BEGIN
           WHILE monsum_noleap[month] LT (t+1) DO month=month+1
           smonth = monstr_long[month-1]
           day = t+1 - monsum_noleap[month-1]
        ENDELSE
        
        sday = string(day,format='(I2.2)')
        syear = string(year,format='(I4)')
        smo = string(month,format='(I2.2)')

        ;; map is plotted by day
        IF psout AND NOT makedifference THEN stime = syear+smo+sday

        IF ~(plottopo OR plotpft) THEN filename += '.'+stime

        IF first THEN print,'Plotting : ',stime
        
        ;; create correct labels
        IF obsvarunits[0] EQ '' THEN units = modvarunits[0] ELSE units = obsvarunits[0]
        IF obsvarlongnames[0] EQ '' THEN longname = modvarlongnames[0] ELSE longname = obsvarlongnames[0]
        IF makemean THEN BEGIN
           IF maketiming THEN BEGIN
              IF makevar THEN BEGIN
                 longname = 'Mean Spring Date'
                 units = 'Months'
              ENDIF ELSE BEGIN
                 longname = 'Spring Date Interannual Variability'
                 units = 'Months'
              ENDELSE
           ENDIF ELSE BEGIN
              IF makevar THEN longname = 'LAI Variability' ELSE $
                 longname = 'LAI Annual Mean'
           ENDELSE
        ENDIF

        ;; plot map if interactive plotting is chosen
        IF ((interactive AND first) OR psout) THEN BEGIN

           ;; plot map of first observation and model only
           m = 0
           o = 0

           IF (colorbar EQ 0) OR (colorbar EQ 2) THEN cbsize = 0.0

           IF psout THEN BEGIN

              ;; setting up plotting to PS device
              set_plot,'ps'

              ;; set up EPS sizes
              ysize = ysize_ps_default     
              xsize = float(ysize) * (cbsize + float(nlon) * dlon / (float(nlat) * dlat) * float(nplot))

              thick=ysize/8.0*3.0 ; plot line thickness
              charthick=2.0       ; charthickness relative to 1.0
              charsize=ysize/11.0 ; charsize relative to 10pt.

              !P.FONT=0

           ENDIF ELSE BEGIN
              ;; set up plotting window in X
              set_plot,'x'
              
              ysize = round(ysize_screen_default)
              xsize = round(float(ysize) * (cbsize + nlon * dlon / (nlat * dlat) * nplot))
              IF xsize GT xsize_screen_max THEN BEGIN
                 ysize = fix(long(ysize) * xsize_screen_max / xsize)
                 xsize = xsize_screen_max
              ENDIF

              device,get_screen_size = screen_size
              
              thick=2.0                               ; plot line thickness
              charthick=1.5                           ; charthickness relative to 1.0
              charsize=1.0 * sqrt(ysize/float(ysize)) ; charsize relative to 10pt.

              !P.FONT=1

              device,decomposed=0,set_font='HELVETICA',/tt_font
              ;; 2D map window
              window,0,xsize=xsize,ysize=ysize
              
           ENDELSE

           IF makedifference AND NOT makeabsdifference THEN BEGIN
              loadct,26,file='brewer.tbl' 
           ENDIF ELSE BEGIN
              IF (varnames[0] EQ 'LAI') OR plotpft THEN BEGIN
                 loadct,14,file='brewer.tbl' 
              ENDIF ELSE BEGIN
                 loadct,16,file='brewer.tbl' 
              ENDELSE
           ENDELSE
           
           IF maketiming THEN loadct,36,file='brewer.tbl'

           tvlct,red,green,blue,/get
           ntable=!d.table_size

           IF maketiming THEN BEGIN
              ntable = 12+4
              red=shift(red,1)
              green=shift(green,1)
              blue=shift(blue,1)
           ENDIF
           
           red[ntable-2]=20B
           green[ntable-2]=20B
           blue[ntable-2]=150B
           red[ntable-3]=128B
           green[ntable-3]=128B
           blue[ntable-3]=128B

           ;; load color bar
           IF psout THEN BEGIN
              red[0]=255B
              green[0]=255B
              blue[0]=255B
              red[ntable-1]=0B
              green[ntable-1]=0B
              blue[ntable-1]=0B
              tvlct,red,green,blue
           ENDIF ELSE BEGIN
              red[0]=0
              green[0]=0
              blue[0]=0
              red[ntable-1]=255B
              green[ntable-1]=255B
              blue[ntable-1]=255B
              tvlct,red,green,blue
           ENDELSE

           IF colorbar EQ 2 THEN BEGIN
              ;; plot the color bar (standalone EPS file)               
              
              IF varnames[0] EQ 'LAI' THEN BEGIN
                 IF maketiming THEN BEGIN
                    colbar_standalone,divisions=12, units=units, varname=longname, crange=[1,12], $
                                      filename=mapdir+filename+'.cb.eps', monthly=1
                 ENDIF ELSE BEGIN
                    colbar_standalone,divisions=6, units=units, varname=longname, filename=mapdir+filename+'.cb.eps'
                 ENDELSE
              ENDIF ELSE BEGIN
                 colbar_standalone,divisions=5, units=units, varname=longname, filename=mapdir+filename+'.cb.eps'
              ENDELSE
           ENDIF

           IF psout THEN BEGIN
              
              set_plot,'ps'

              print,'Writing EPS: ',mapdir+filename+'.eps'
              DEVICE, /encapsul, bits_per_pixel=8,/color, font_size=5, $
                      Filename=mapdir+filename+'.eps', preview=0,/isolatin1, $
                      xsize=xsize, ysize=ysize,/HELVETICA
           ENDIF ELSE BEGIN
              wset,0
           ENDELSE

           ;; plot obs+model 2D maps
           FOR p = 0, nplot-1 DO BEGIN

              ;; window margins
              x0 = (cbsize + xylim + float(p) * nlon * dlon / (nlat * dlat) ) / (cbsize + float(nplot) * nlon * dlon / (nlat * dlat) ) * xsize
              x1 = (cbsize - xylim + float(p+1) * nlon * dlon / (nlat * dlat) ) / (cbsize + float(nplot) * nlon * dlon / (nlat * dlat) ) *xsize
              y0 = xylim * ysize
              y1 = ysize - xylim * ysize

              xn0 = x0 / xsize
              xn1 = x1 / xsize
              yn0 = y0 / ysize
              yn1 = y1 / ysize

              IF regrid THEN BEGIN

                 ;; --> not working right now!!!
                 IF rotatedpole THEN BEGIN
                    map = map_proj_init(cylindrical,limit=[latmin,lonmin,latmax,lonmax],center_longitude = 180.+pollon,center_latitude=90.-pollat)
                 ENDIF ELSE BEGIN
                    map = map_proj_init(cylindrical,limit=[latmin,lonmin,latmax,lonmax])
                 ENDELSE

                 IF (plotmodels AND (NOT plotobs)) OR (plotmodels AND plotobs AND (p EQ 1)) THEN $
                    outmod = map_image(reform(modvar[*,*,t,y,m]), modlons, modlats, $
                                       xstart=xp0, ystart=yp0, xsize = npx, ysize = npy, missing = !values.f_nan)
                 
                 IF (plotobs AND (NOT plotmodels)) OR (plotmodels AND plotobs AND (p EQ 0)) THEN $
                    outobs = map_image(reform(obsvar[*,*,obst,y,o]), reform(lons[*,0]), reform(lats[0,*]), $
                                       xstart=xp0, ystart=yp0, xsize = npx, ysize = npy, missing = !values.f_nan)
                 
              ENDIF ELSE BEGIN
                 ;; --> use this one:
                 IF psout THEN BEGIN
                    nx=nlon*10
                    ny=nlat*10
                 ENDIF ELSE BEGIN
                    nx = round(x1 - x0)
                    ny = round(y1 - y0)
                 ENDELSE

                 IF (plotmodels AND (NOT plotobs)) OR (plotmodels AND plotobs AND (p EQ 1)) THEN BEGIN
                    outmod =congrid(reform(modvar[*,*,t,y,m],nlon,nlat),nx,ny) 
                 ENDIF
                 IF (plotobs AND (NOT plotmodels)) OR (plotmodels AND plotobs AND (p EQ 0)) THEN BEGIN
                    IF hammeraitoff THEN BEGIN
                       set_plot,'z'
                       device,set_resolution=[nx,ny]
                       MAP_SET,0,0,0,/HAMMER,XMARGIN=[0,0],YMARGIN=[0,0],/noborder
;;                             !map = map_proj_init(13,datum=8,limit=[latmin,lonmin,latmax,lonmax])
                       outobs =map_image(reform(obsvar[*,*,obst,y],nlon,nlat),/bilinear,compress=1,missing=!values.f_nan)
                       IF psout AND NOT psseries THEN set_plot,'ps' ELSE set_plot,'x'
                    ENDIF ELSE BEGIN
                       outobs =congrid(reform(obsvar[*,*,obst,y],nlon,nlat),nx,ny)
                    ENDELSE
                 ENDIF

              ENDELSE

              IF (plotmodels AND (NOT plotobs)) OR (plotmodels AND plotobs AND (p EQ 1)) THEN BEGIN
                 badidx = where(finite(outmod,/nan),badcount)
                 IF badcount GT 0 THEN outmod[badidx] = nodata
                 outbyte = (bytscl(outmod,min=range[0],max=range[1],top=ntable-4)>1B)*(outmod NE nodata) + $
                           (ntable-3)*(outmod EQ nodata)
                 IF psout AND NOT psseries THEN $
                    tv,outbyte,x0,y0,xsize=x1 - x0,ysize=y1 - y0,/centimeters $
                 ELSE $
                    tv,outbyte,x0,y0,xsize=nx,ysize=ny,/device

                 
              ENDIF
              IF (plotobs AND (NOT plotmodels)) OR (plotmodels AND plotobs AND (p EQ 0)) THEN BEGIN
                 badidx = where(finite(outobs,/nan),badcount)
                 IF badcount GT 0 THEN outobs[badidx] = nodata
                 outbyte = (bytscl(outobs,min=range[0],max=range[1],top=ntable-4)>1B)*(outobs NE nodata) + $
                           (ntable-3)*(outobs EQ nodata)
                 IF psout AND NOT psseries THEN $
                    tv,outbyte,x0,y0,xsize=x1-x0,ysize=y1-y0,/centimeters $
                 ELSE $
                    tv,outbyte,x0,y0,xsize=nx,ysize=ny,/device
              ENDIF
              
              IF plottitle THEN BEGIN

                 IF ~plotpft THEN BEGIN
                    IF simulator THEN obstitle = 'SATELLITE SIMULATION' ELSE obstitle = 'SATELLITE OBSERVATION'
                    modtitle = 'MODEL PREDICTION'
                    IF makedifference AND NOT makeabsdifference THEN BEGIN
                       obstitle = 'Mean Bias Error MODEL-SATELLITE'
                       modtitle = 'MODEL UNCERTAINTY'
                    ENDIF
                    IF makeabsdifference AND makedifference THEN BEGIN
                       obstitle = 'Mean Absolute Bias MODEL-SATELLITE'
                       modtitle = 'MODEL UNCERTAINTY'
                    ENDIF
                 ENDIF
                 
                 IF (p EQ 0) AND (nplot GE 2) THEN $
                    xyouts,x0/xsize,1.-xylim*0.75,obstitle,charsize=charsize*2.5,charthick=charthick*2, $
                           align=0.0,color=ntable-1,/normal $
                 ELSE BEGIN
                    IF (p EQ 1) OR plotmodels THEN BEGIN
                       xyouts,x1/xsize,1.-xylim*0.75,modtitle,$
                              charsize=charsize*2.5,charthick=charthick*2, $
                              align=1.0,color=ntable-1,/normal
                    ENDIF ELSE BEGIN
                       xyouts,x1/xsize,1.-xylim*0.75,obstitle,$
                              charsize=charsize*2.5,charthick=charthick*2, $
                              align=1.0,color=ntable-1,/normal
                    ENDELSE
                 ENDELSE
              ENDIF

              IF makecontour THEN BEGIN
                 ;; contour plot of variability or trends
                 IF maketiming THEN BEGIN
                    contour,reform(modvar[*,*,1,0]),levels=[10,20], c_labels=[1,1], c_linestyle=[0,0], $
                            color=ntable-1, Position=[x0/xsize, y0/ysize, x1/xsize, y1/ysize],/normal,/noerase, $
                            xstyle=5,ystyle=5,/zstyle,c_charsize=charsize,c_thick=thick*0.5
                    
                 ENDIF ELSE BEGIN
                    IF makevar THEN BEGIN
                       contour,reform(modvar[*,*,1,0]),levels=[-0.1,0.1], c_labels=[1,1], c_linestyle=[0,0], $
                               color=ntable-1, Position=[x0/xsize, y0/ysize, x1/xsize, y1/ysize],/normal,/noerase, $
                               xstyle=5,ystyle=5,/zstyle,c_charsize=charsize,c_thick=thick*0.5
                    ENDIF ELSE BEGIN
                       contour,reform(modvar[*,*,1,0]),levels=[1.0,2.0], c_labels=[1,1], c_linestyle=[0,0], $
                               color=ntable-1, Position=[x0/xsize, y0/ysize, x1/xsize, y1/ysize],/normal,/noerase, $
                               xstyle=5,ystyle=5,/zstyle,c_charsize=charsize,c_thick=thick*0.5
                    ENDELSE
                 ENDELSE
              ENDIF

              IF plot_map THEN BEGIN
                 IF rotatedpole THEN BEGIN
                    rpol2ll,lonmin,latmax,lon0,lat0,pollon=pollon,pollat=pollat
                    rpol2ll,lonmax,latmax,lon1,lat1,pollon=pollon,pollat=pollat
                    rpol2ll,lonmin,latmin,lon3,lat3,pollon=pollon,pollat=pollat
                    rpol2ll,lonmax,latmin,lon2,lat2,pollon=pollon,pollat=pollat
                    
                    map_set,90.-pollat,180.+pollon,limit=[lat0,lon0,lat1,lon1,lat2,lon2,lat3,lon3],/cylindrical,/noerase,$
                            xmargin=[0,0],ymargin=[0,0],/noborder
                    map = !map
                 ENDIF ELSE BEGIN
                    IF hammeraitoff THEN BEGIN
                       map = map_proj_init(13,datum=8,limit=[latmin,lonmin,latmax,lonmax])
                    ENDIF ELSE BEGIN
                       map = map_proj_init(8,datum=8,limit=[latmin,lonmin,latmax,lonmax])
                    ENDELSE
                 ENDELSE
                 plot, map.uv_box[[0, 2]], map.uv_box[[1, 3]], Position=[x0/xsize, y0/ysize, x1/xsize, y1/ysize], $
                       /NoData, XStyle=5, YStyle=5, /NoErase,/normal

                 IF (lonmax-lonmin) LT 40. THEN BEGIN
;                    map_continents,/coasts,/rivers,color=ntable-2,mlinethick=2,/hires,fill=0, map_structure=map
                    map_continents,/coasts,color=ntable-2,mlinethick=2,/hires,fill=0, map_structure=map
                    map_continents,/countries,color=ntable-1,mlinethick=1,/hires,fill=0, map_structure=map
                 ENDIF ELSE BEGIN
                    map_continents,color=ntable-1,mlinethick=2,fill=0, map_structure=map
                 ENDELSE

                 IF (NOT hammeraitoff) AND ((max(lons)-min(lons)) LT 300.) THEN plot,[mean(lons)],[mean(lats)], Position=[x0/xsize, y0/ysize, x1/xsize, y1/ysize], $
                    /NoData, XStyle=1, YStyle=1, xrange=[lonmin,lonmax],yrange=[latmin,latmax], $
                    /NoErase,/normal, xticks=4,yticks=4,color=ntable-1,charsize=charsize*1.5,charthick=charthick, $
                    thick=thick, xthick=thick, ythick=thick,xtitle='LON',ytitle='LAT'
              ENDIF

           ENDFOR

           IF plotdate THEN BEGIN

              IF nplot EQ 1 THEN plotdatex = x0/xsize ELSE plotdatex = 0.47

              IF makedifference OR makemean OR makevar THEN BEGIN
                 xyouts,plotdatex,1.-xylim*0.75,stime,charsize=charsize*2.5,charthick=charthick*2, $
                        align=0.0,color=ntable-1,/normal
              ENDIF ELSE BEGIN
                 
                 xyouts,plotdatex,1.-xylim*0.75,sday,charsize=charsize*2.5,charthick=charthick*2, $
                        align=0.0,color=ntable-1,/normal
                 xyouts,plotdatex+0.05 / (float(nplot) * nlon * dlon / (nlat * dlat)), $
                        1.-xylim*0.75,smonth,charsize=charsize*2.5,charthick=charthick*2, $
                        align=0.0,color=ntable-1,/normal
                 xyouts,plotdatex+0.15 / (float(nplot) * nlon * dlon / (nlat * dlat)), $
                        1.-xylim*0.75,syear,charsize=charsize*2.5,charthick=charthick*2, $
                        align=0.0,color=ntable-1,/normal
              ENDELSE

              IF plotexperiment THEN xyouts,(0.5*(x1-x0)+x0)/xsize,1.-xylim*0.75,experiments[0],charsize=charsize*2.5,charthick=charthick*2, $
                                            align=0.5,color=ntable-1,/normal

           ENDIF

           IF colorbar EQ 1 THEN BEGIN
              ;; plot the color bar                
              IF varnames[0] EQ 'LAI' THEN BEGIN
                 IF maketiming THEN BEGIN
                    colbar,divisions=12, units=units, varname=longname, crange=[1,12],monthly=1
                 ENDIF ELSE BEGIN
                    colbar,divisions=6, units=units, varname=longname
                 ENDELSE
              ENDIF ELSE BEGIN
                 colbar,divisions=5, units=units, varname=longname
              ENDELSE
           ENDIF
           
           IF psout THEN BEGIN
              DEVICE,/CLOSE      
              set_plot,'x'
              
              spawn,'convert -depth 8 -density 300 -resize 1280 -flatten '+mapdir+filename+'.eps '+mapdir+filename+'.png'
;;              spawn,'epstopdf '+mapdir+filename+'.eps'
              spawn,'rm '+mapdir+filename+'.eps'

              IF plottopo OR plotpft THEN stop
              IF makedifference THEN stop

           ENDIF 

        ENDIF ;; (interactive or psout) AND first

        IF interactive OR psseries THEN BEGIN

           wset,0
           
           keyname = GET_KBRD(0,/key_name)

           cursor, cx, cy, 0, /device

           IF (!mouse.button eq 1) THEN BEGIN


              IF makeseas OR scatter OR taylor THEN aspect = 1

              ;; time series window
              IF psseries THEN BEGIN
                 
                 set_plot,'ps'

                 ;; set up EPS sizes
                 tsysize = ysize_ps_default     
                 tsxsize = aspect*tsysize
                 
                 thick=tsysize/8.0*3.0 ; plot line thickness
                 charthick=2.0         ; charthickness relative to 1.0
                 charsize=tsysize/11.0 ; charsize relative to 10pt.
                 
                 !P.FONT=0
                 
                 print,'Writing EPS: ',plotdir+filename+'.eps'
                 DEVICE, /encapsul, bits_per_pixel=8,/color, font_size=5, $
                         Filename=plotdir+filename+'.eps', preview=0,/isolatin1, $
                         xsize=tsxsize, ysize=tsysize,/HELVETICA

                 IF colored THEN loadct,50,file='brewer.tbl' ELSE loadct,17,file='brewer.tbl'
                 tvlct,red,green,blue,/get
                 ntable=!d.table_size
                 green=shift(green,1)
                 blue=shift(blue,1)
                 red=shift(red,1)

                 red[0]=255B
                 green[0]=255B
                 blue[0]=255B
                 red[ntable-1]=0B
                 green[ntable-1]=0B
                 blue[ntable-1]=0B
              ENDIF ELSE BEGIN
                 tsysize = ysize_screen_default
                 tsxsize = aspect*tsysize

                 IF tsxsize GT xsize_screen_max THEN BEGIN
                    tsysize = fix(long(tsysize) * xsize_screen_max / tsxsize)
                    tsxsize = xsize_screen_max
                 ENDIF

                 window,1,xsize=tsxsize,ysize=tsysize
                 IF colored THEN loadct,50,file='brewer.tbl' ELSE loadct,17,file='brewer.tbl'
                 tvlct,red,green,blue,/get
                 ntable=!d.table_size
                 green=shift(green,1)
                 blue=shift(blue,1)
                 red=shift(red,1)

                 red[0]=0B
                 green[0]=0B
                 blue[0]=0B
                 red[ntable-1]=255B
                 green[ntable-1]=255B
                 blue[ntable-1]=255B
              ENDELSE

              tvlct,red,green,blue
              
              if (nplot EQ 2) THEN p=1 ELSE p=0
              ;; window margins
              x0 = (cbsize + xylim + float(p) * nlon * dlon / (nlat * dlat) ) / $
                   (cbsize + float(nplot) * nlon * dlon / (nlat * dlat) ) * xsize
              x1 = (cbsize - xylim + float(p+1) * nlon * dlon / (nlat * dlat) ) / $
                   (cbsize + float(nplot) * nlon * dlon / (nlat * dlat) ) *xsize
              y0 = xylim * ysize
              y1 = ysize - xylim * ysize

              cx = fix( ( ( ( (!mouse.x - x0)/(x1-x0) )>0.)<0.99999)*(nlon) )
              cy = fix( ( ( ( (!mouse.y - y0)/(y1-y0) )>0.)<0.99999)*(nlat) )
              
              print,'Cursor: ', cx, cy, 'lon/lat: ',lons[cx],lats[cy]
              
              print,'PFT distribution: '
              FOR p = 0,npft-1 DO IF pft[cx,cy,p] GT 5.0 THEN $
                 print,p+1,' ',strtrim(pftnames[p]),': ',fix(pft[cx,cy,p]),' %'
              print

              IF series THEN BEGIN
                 ;; plot time series
                 
                 ;; set up plot labeling
                 IF makedifference THEN BEGIN
                    IF plotobs THEN ytitle = obsvarlongnames[0] + ' Error ('+obsvarunits[0]+')'
                    IF plotmodels THEN ytitle = modvarlongnames[0] + ' Error ('+modvarunits[0]+')'
                 ENDIF ELSE BEGIN
                    IF plotobs THEN ytitle = obsvarlongnames[0] + ' ('+obsvarunits[0]+')'
                    IF plotmodels THEN ytitle = modvarlongnames[0] + ' ('+modvarunits[0]+')'
                 ENDELSE

                 secondaxistitle = ""
                 secondaxisunits = ""
                 IF plotobs THEN BEGIN
                    FOR o=0,nobs-1 DO BEGIN
                       IF axisindex[o] EQ 1 THEN BEGIN
                          secondaxistitle=obsvarlongnames[o]
                          secondaxisunits=obsvarunits[o]
                       ENDIF
                    ENDFOR
                 ENDIF
                 
                 IF plotmodels THEN BEGIN
                    FOR m=0,nmodels-1 DO BEGIN
                       IF axisindex[m+nobs] EQ 1 THEN BEGIN
                          secondaxistitle=modvarlongnames[m]
                          secondaxisunits=modvarunits[m]
                       ENDIF
                    ENDFOR
                 ENDIF
                 
                 IF secondaxistitle EQ "Light Control Factor" THEN secondaxistitle = "Control Factor"

                 IF nyears EQ 1 THEN BEGIN
                    xticks=round(dpy/30.5)
                    xminor=1
                    xtickname='       '+[strtrim(monstr_short,2),'        ']
                    xtitle=stime
                    !x.ticklen=0.01
                 ENDIF ELSE BEGIN
                    xticks=nyears
                    xminor=12
                    xtickname=replicate(' ',nyears+1)
                    xtitle='Years'
                    IF nyears LE 5 THEN !x.ticklen=0.5 ELSE !x.ticklen=0.01              
                 ENDELSE


                 title = sitenames[0]

                 ;; create plot axis
                 xrange=[0,nyears*dpy-1]
                 yrange=range
                 plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1+8*secondaxis, $
                      xticks=xticks,xtickname=xtickname,color=ntable-1, /nodata,$
                      xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                      xthick=2.0, ythick=2.0, xtitle=xtitle, title=title, $
                      /normal,position=[0.18/aspect,0.10,1.-(0.05+0.1*secondaxis)/aspect,0.92], ytitle=ytitle, $
                      xticklen = 0.01, yticklen = 0.01 / aspect
                 
                 ;; create centered labels with year
                 IF nyears GT 1 THEN BEGIN
                    FOR yi=0,nyears-1 DO BEGIN
                       dx = (1.-(0.05+0.1*secondaxis)/aspect - 0.18/aspect)/nyears
                       xpos = 0.18/aspect + dx/2. + yi*dx
                       xyouts,xpos,.06,strtrim(years[yi],2),color=ntable-1,charsize=2.0*charsize,charthick=charthick, $
                              Alignment=0.5,Orientation=0.0, /Normal
                    ENDFOR
                 ENDIF
                 
                 q = 0

                 ;; plot data and models
                 FOR ax=0,secondaxis DO BEGIN

                    IF ax EQ 1 THEN BEGIN
                       yrange=secondaxisrange
                       axis,yaxis=1,yrange=yrange,ystyle=1, $
                            color=ntable-1, /save, $
                            yminor=1,charsize=2.0*charsize,charthick=charthick, $
                            ythick=2.0,/normal, ytitle=secondaxistitle+' ('+secondaxisunits+')'
                    ENDIF

                    ;; plot data
                    IF plotobs THEN BEGIN
                       indexobs = obsdays
                       FOR yi = 1,nyears-1 DO indexobs = [indexobs,obsdays+yi*dpy]
                       
                       FOR o=0,nobs-1 DO BEGIN
                          IF axisindex[o] EQ ax THEN BEGIN
                             
                             obstitle='OBS: '+varnames[o]

                             ;; plot observation spread
                             IF ploterror THEN BEGIN
                                obstemp = make_array(2,ntobs*nyears,/float)
                                nobstemp = make_array(2,ntobs*nyears,/long)
                                wgttemp = make_array(2,ntobs*nyears,/float)
                                FOR i=0L,n_elements(cx)-1L DO BEGIN
                                   temp = obsvarerror[*,cx[i],cy[i],*,*,o]
                                   badidx = where(finite(temp,/nan),badcount,complement=goodidx,ncomplement=goodcount)
                                   IF goodcount GT 0 THEN BEGIN
                                      obstemp[goodidx] += temp[goodidx]*cos(lats[cy[i]]*!dtor)
                                      nobstemp[goodidx] += 1L
                                      wgttemp[goodidx] += cos(lats[cy[i]]*!dtor)
                                   ENDIF
                                ENDFOR
                                badidx = where(nobstemp LE 0.25*n_elements(cx),badcount,complement=goodidx,ncomplement=goodcount)
                                IF badcount GT 0 THEN obstemp[badidx] = !values.f_nan
                                IF goodcount GT 0 THEN obstemp[goodidx] /= wgttemp[goodidx]

                                IF total(finite(obstemp)) GT 0.0 THEN BEGIN                         
                                   overplot,indexobs,obstemp,obstitle,q,0
                                ENDIF
                                
                             ENDIF
                             
                             ;; plot observation mean
                             obstemp = make_array(ntobs*nyears,/float)
                             nobstemp = make_array(ntobs*nyears,/long)
                             wgttemp = make_array(ntobs*nyears,/float)
                             FOR i=0L,n_elements(cx)-1L DO BEGIN
                                temp = obsvar[cx[i],cy[i],*,*,o]
                                badidx = where(finite(temp,/nan),badcount,complement=goodidx,ncomplement=goodcount)
                                IF goodcount GT 0 THEN BEGIN
                                   obstemp[goodidx] += temp[goodidx]*cos(lats[cy[i]]*!dtor)
                                   nobstemp[goodidx] += 1L
                                   wgttemp[goodidx] += cos(lats[cy[i]]*!dtor)
                                ENDIF
                             ENDFOR
                             badidx = where(nobstemp LE 0.25*n_elements(cx),badcount,complement=goodidx,ncomplement=goodcount)
                             IF badcount GT 0 THEN obstemp[badidx] = !values.f_nan
                             IF goodcount GT 0 THEN obstemp[goodidx] /= wgttemp[goodidx]
                             IF total(finite(obstemp)) GT 0.0 THEN BEGIN                         
                                overplot,indexobs,obstemp,obstitle,q,0
                                q=q+1
                             ENDIF   
                          ENDIF
                       ENDFOR
                    ENDIF

                    IF plotmodels THEN BEGIN
                       FOR m = 0,nmodels-1 DO BEGIN
                          IF axisindex[nobs+m] EQ ax THEN BEGIN
                             indexmod = indgen(ntmod*nyears)

                             modtitle='MOD: '+varnames[m]
                             
                             ;; plot ensemble spread
                             IF ploterror THEN BEGIN

                                modtemp = make_array(2,ntmod*nyears,/float)
                                nmodtemp = make_array(2,ntmod*nyears,/long)
                                wgttemp = make_array(2,ntmod*nyears,/float)
                                FOR i=0L,n_elements(cx)-1L DO BEGIN
                                   temp = make_array(2,1,1,ntmod,nyears,/float)
                                   temp[0,0,0,*,*] = modvar[cx[i],cy[i],*,*,m] + modvarerror[cx[i],cy[i],*,*,m]
                                   temp[1,0,0,*,*] = modvar[cx[i],cy[i],*,*,m] - modvarerror[cx[i],cy[i],*,*,m]

                                   IF plotobs AND NOT ignoregaps THEN BEGIN
                                      obstemp = obsvarerror[*,cx[i],cy[i],congrid(indgen(ntobs),ntmod),*,m]
                                      badidx = where(finite(temp,/nan) OR finite(obstemp,/nan),badcount,$
                                                     complement=goodidx,ncomplement=goodcount)
                                   ENDIF ELSE BEGIN
                                      badidx = where(finite(temp,/nan),badcount,complement=goodidx,ncomplement=goodcount)
                                   ENDELSE
                                   IF goodcount GT 0 THEN BEGIN
                                      modtemp[goodidx] += temp[goodidx]*cos(lats[cy[i]]*!dtor)
                                      nmodtemp[goodidx] += 1L
                                      wgttemp[goodidx] += cos(lats[cy[i]]*!dtor)
                                   ENDIF
                                ENDFOR
                                badidx = where(nmodtemp LE 0.25*n_elements(cx),badcount,complement=goodidx,ncomplement=goodcount)
                                IF badcount GT 0 THEN modtemp[badidx] = !values.f_nan
                                IF goodcount GT 0 THEN modtemp[goodidx] /= wgttemp[goodidx]
                                
                                overplot,indexmod,modtemp,modtitle,q,1
                             ENDIF
                             
                             ;; plot ensemble mean
                             modtemp = make_array(ntmod*nyears,/float)
                             nmodtemp = make_array(ntmod*nyears,/long)
                             wgttemp = make_array(ntmod*nyears,/float)
                             FOR i=0L,n_elements(cx)-1L DO BEGIN
                                temp = modvar[cx[i],cy[i],*,*,m]
                                IF plotobs AND NOT ignoregaps THEN BEGIN
                                   obstemp = obsvar[cx[i],cy[i],congrid(indgen(ntobs),ntmod),*,m]
                                   badidx = where(finite(temp,/nan) OR finite(obstemp,/nan),badcount,$
                                                  complement=goodidx,ncomplement=goodcount)
                                ENDIF ELSE BEGIN
                                   badidx = where(finite(temp,/nan),badcount,complement=goodidx,ncomplement=goodcount)
                                ENDELSE
                                IF goodcount GT 0 THEN BEGIN
                                   modtemp[goodidx] += temp[goodidx]*cos(lats[cy[i]]*!dtor)
                                   nmodtemp[goodidx] += 1L
                                   wgttemp[goodidx] += cos(lats[cy[i]]*!dtor)
                                ENDIF
                                
                             ENDFOR

                             badidx = where(nmodtemp LE 0.25*n_elements(cx),badcount,complement=goodidx,ncomplement=goodcount)
                             IF badcount GT 0 THEN modtemp[badidx] = !values.f_nan
                             IF goodcount GT 0 THEN modtemp[goodidx] /= wgttemp[goodidx]

                             overplot,indexmod,modtemp,modtitle,q,1
                             q=q+1
                          ENDIF
                       ENDFOR
                    ENDIF

                 ENDFOR ;; axis

              ENDIF ;; series

              IF scatter THEN BEGIN

                 secondaxis = 0
                 
                 IF plotobs THEN xtitle = 'OBS: ' + obsvarnames[0] + ' ('+obsvarunits[0]+')'
                 IF plotmodels THEN ytitle = 'MOD: ' + modvarnames[0] + ' ('+modvarunits[0]+')'
                 
                 IF plotobs THEN xtitle = 'OBS: ' + longname + ' ('+obsvarunits[0]+')'
                 IF plotmodels THEN ytitle = 'MOD: ' + longname + ' ('+modvarunits[0]+')'
                 
                 !x.ticklen=0.01
                 !y.ticklen=0.01
                 
                 modtitle = experiments           
                 title = sitenames[0]
                 
                 ;; create plot axis

                 xrange=range
                 yrange=range
                 plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=1,ystyle=1, $
                      color=ntable-1, /nodata,$
                      xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                      xthick=2.0, ythick=2.0, xtitle=xtitle, ytitle=ytitle, title=title, $
                      /normal,position=[0.18,0.10,1.-(0.05+0.1*secondaxis),0.92]

                 FOR m = 0,nmodels-1 DO BEGIN
                    obstemp = make_array(n_elements(cx),ntobs*nyears,/float)
                    modtemp = make_array(n_elements(cx),ntobs*nyears,/float)

                    FOR i=0L,n_elements(cx)-1L DO BEGIN
                       obstemp[i,*] = obsvar[cx[i],cy[i],*,*,m]
                       modtemp[i,*] = modvar[cx[i],cy[i],obsdays[*]-1,*,m]
                    ENDFOR

                    obstemp = obstemp[*]
                    modtemp = modtemp[*]
                    
                    IF ploterror THEN BEGIN

                       modmintemp = make_array(n_elements(cx),ntobs*nyears,/float)
                       modmaxtemp = make_array(n_elements(cx),ntobs*nyears,/float)
                       obsmintemp = make_array(n_elements(cx),ntobs*nyears,/float)
                       obsmaxtemp = make_array(n_elements(cx),ntobs*nyears,/float)

                       FOR i=0L,n_elements(cx)-1L DO BEGIN
                          obsmaxtemp[i,*] = obsvarerror[0,cx[i],cy[i],*,*,m]
                          obsmintemp[i,*] = obsvarerror[1,cx[i],cy[i],*,*,m]
                          modmaxtemp[i,*] = modtemp[i,*] + modvarerror[cx[i],cy[i],obsdays[*]-1,*,m]
                          modmintemp[i,*] = modtemp[i,*] - modvarerror[cx[i],cy[i],obsdays[*]-1,*,m]
                       ENDFOR

                       ntemp = n_elements(modmaxtemp)
                       moderrtemp = make_array(2,ntemp,/float)
                       obserrtemp = make_array(2,ntemp,/float)
                       moderrtemp[0,*] = modmaxtemp
                       moderrtemp[1,*] = modmintemp
                       obserrtemp[0,*] = obsmaxtemp
                       obserrtemp[1,*] = obsmintemp
                    ENDIF

                    indexp = where(finite(obstemp),count)
                    IF count GT 0 THEN BEGIN
                       scatplot,obstemp[indexp],modtemp[indexp],modtitle[m],m+1
                       IF ploterror THEN scatplot,obserrtemp[*,indexp],moderrtemp[*,indexp],modtitle[m],m+1

                       ;; ordinary least squares
                       a=regress(obstemp[indexp],modtemp[indexp],const=b,mcorrelation=corr, $
                                 measure_errors=replicate(10.,n_elements(indexp)),/double)
                       
                       rms=sqrt(total((obstemp[indexp] - modtemp[indexp])^2.0)/float(count))

                       corr_all[m,sites[0]-sites0] = corr
                       rmse_all[m,sites[0]-sites0] = corr

                    ENDIF ELSE BEGIN

                       corr_all[m,sites[0]-sites0] = -9999.
                       rmse_all[m,sites[0]-sites0] = -9999.

                    ENDELSE

                 ENDFOR

              ENDIF ;; scatter

              IF taylor THEN BEGIN

;                 modtitle = experiments 
;                 IF manysites THEN BEGIN
                 title = varnames[0] + ': '+experiments[0] 
;                 ENDIF ELSE BEGIN
;                    title = sitenames[0]
;                 ENDELSE

                 xtitle='standard deviation (normalized)'
                 ytitle='standard deviation (normalized)'
;                 xtitle='bias (normalized)'
;                 ytitle='bias (normalized)'
                 
                 !x.ticklen=0.01
                 !y.ticklen=0.01

                 ;; regular Taylor plot: normdev,R
                 maxdev=2.0
                 xrange=[0.,maxdev]
                 yrange=xrange

                 ;; create plot axis

                 plot,[0],[0],yrange=yrange,xrange=xrange,xstyle=9,ystyle=9, $
                      color=ntable-1, /nodata,$
                      xminor=1,yminor=1,charsize=2.0*charsize,charthick=charthick, thick=thick, $
                      xthick=thick, ythick=thick, xtitle=xtitle, ytitle=ytitle, $
                      /normal,position=[0.15,0.10,0.85,0.825],/polar,noclip=1
                 
                 npoints=100
                 oplot,replicate(1.0,npoints+1),findgen(npoints+1)/npoints*!pi/2.,/polar,noclip=1,color=ntable-1
                 oplot,replicate(maxdev,npoints+1),findgen(npoints+1)/npoints*!pi/2.,/polar,thick=thick,noclip=1,color=ntable-1

                 npoints=13
                 points=[findgen(11)/10.,0.95, 0.99]
                 FOR i=0,npoints-1 DO BEGIN
                    theta=acos(points[i])
                    r=maxdev*1.05
                    xi=cos(theta)*r
                    yi=sin(theta)*r
                    sint=''
                    sflt=''
                    corr=points[i]
                    reads,fix(corr),sint
                    sint=strtrim(sint,2)
                    reads,fix(100.*(corr-fix(corr))),sflt
                    sflt=strtrim(sflt,2)
                    IF (strlen(sflt) EQ 2) AND (strmid(sflt,1,1) EQ '0') THEN sflt=strmid(sflt,0,1)
                    oplot,[maxdev*0.95,maxdev],[theta,theta],/polar,thick=thick,noclip=1,color=ntable-1
                    xyouts,xi,yi,sint+'.'+sflt,charsize=2.0*charsize,charthick=charthick, orientation=theta/!pi*180.,color=ntable-1
                 ENDFOR
                 xyouts,maxdev*0.8,maxdev*0.8,'Correlation',charsize=2.0*charsize,charthick=charthick, orientation=-45., align=0.5,color=ntable-1
                 xyouts,0.5,0.94,title,charsize=2.5*charsize,charthick=charthick, align=0.5,/normal,color=ntable-1
                 
                 ;; plot data and models

                 FOR m = 0,nmodels-1 DO BEGIN

                    IF selpft GT 0 THEN BEGIN

                       ;; latitude-weighted spatial compositing by pft
                       obstemp = make_array(npft,ntobs*nyears,/float)
                       modtemp = make_array(npft,ntobs*nyears,/float)
                       wgttemp = make_array(npft,ntobs*nyears,/float)
                       ntemp = make_array(npft,ntobs*nyears,/float)
                       
                       FOR p=0,npft-1 DO BEGIN
                          FOR i=0L,n_elements(cx)-1L DO BEGIN
                             IF pft[cx[i],cy[i],p] GT 25 THEN BEGIN
                                temp = reform(obsvar[cx[i],cy[i],*,*,m])
                                temp2 = reform(modvar[cx[i],cy[i],obsdays[*]-1,*,m])
                                idx = where(finite(temp) AND finite(temp2),count)
                                IF count GT 0L THEN BEGIN
                                   obstemp[p,idx] += temp[idx]*cos(lats[cy[i]]*!dtor)
                                   modtemp[p,idx] += temp2[idx]*cos(lats[cy[i]]*!dtor)
                                   wgttemp[p,idx] += cos(lats[cy[i]]*!dtor)
                                   ntemp[p,idx] += 1L
                                ENDIF
                             ENDIF
                          ENDFOR
                       ENDFOR
                       
                       idx = where(ntemp GT 0L,count)
                       IF count GT 0L THEN BEGIN
                          obstemp[idx] /= wgttemp[idx]
                          modtemp[idx] /= wgttemp[idx]
                       ENDIF
                       
                       idx = where(ntemp EQ 0L,count)
                       IF count GT 0L THEN BEGIN
                          obstemp[idx] = !values.f_nan
                          modtemp[idx] = !values.f_nan
                       ENDIF
                       
                       FOR p = 0,npft-3 DO BEGIN
                          idx = where(finite(obstemp[p,*]),count)
                          IF count GT 0L THEN $
                             taylorplot,obstemp[p,idx],modtemp[p,idx],modtitle[m],p,string(p+1,format='(I2)')
                       ENDFOR

                    ENDIF ELSE BEGIN
                       obstemp = make_array(n_elements(cx),ntobs*nyears,/float)
                       modtemp = make_array(n_elements(cx),ntobs*nyears,/float)
                       
                       FOR i=0,n_elements(cx)-1 DO BEGIN
                          obstemp[i,*] = obsvar[cx[i],cy[i],*,*,m]
                          modtemp[i,*] = modvar[cx[i],cy[i],obsdays[*]-1,*,m]
                       ENDFOR
                       
                       obstemp = obstemp[*]
                       modtemp = modtemp[*]
                       
                       indexp = where(finite(obstemp),count)
                       
                       IF count GT 0 THEN BEGIN
                          IF manysites THEN $
                             taylorplot,obstemp[indexp],modtemp[indexp],modtitle[m],4,sitenames[m] $
                          ELSE $
                             taylorplot,obstemp[indexp],modtemp[indexp],modtitle[m],m+1,varnames[m]
                       ENDIF

                    ENDELSE

                 ENDFOR

              ENDIF ;; taylor
              
              IF psseries THEN BEGIN
                 DEVICE,/CLOSE      
                 set_plot,'x'
                 
                 spawn,'convert -depth 8 -density 300 -resize 1024 -flatten '+plotdir+filename+'.eps '+plotdir+filename+'.png'
                 spawn,'epstopdf '+plotdir+filename+'.eps'
                 spawn,'rm '+plotdir+filename+'.eps'

              ENDIF 

;              IF psseries THEN BEGIN
;                 t=tend
;                 y=nyears
;              ENDIF

           ENDIF ;; mouse button
           
           wait, 0.1
           
        END ;; interactive or psseries

        first = 0

        IF psout THEN t=t+1
        
        IF keyname EQ 'LEFT' THEN BEGIN
           t=t-1  
           first = 1
        ENDIF
        IF keyname EQ 'RIGHT' THEN BEGIN
           t=t+1
           first = 1
        ENDIF
        IF keyname EQ 'UP' THEN BEGIN
           y=y-1
           first = 1
        ENDIF
        IF keyname EQ 'DOWN' THEN BEGIN
           y=y+1
           first = 1
        ENDIF
        
        IF t EQ tend THEN BEGIN
           t=0
           y=y+1
        ENDIF
        IF y GE nyears THEN keyname='x'
        
        IF t LT 0 THEN BEGIN
           t=tend
           y=y-1
        ENDIF
        IF y LT 0 THEN keyname='x'

     ENDWHILE ;; while loop

  ENDFOR

  ;; print a table with statistical values in LaTeX format
  IF makelatextable THEN BEGIN
     
     filename='table.'+experiments[0]+'.'+varnames[0]+'.tex'
     openw,lun,tabledir+filename,/get_lun
     
     ;; create table of Correlation

      ;text = strtrim(sitenames[sites[0]-1],2)

      FOR p=0,npft-3 DO BEGIN

;         text = string(p+1,format='(I3)')
         text = ''

         x = reform(obstemp[p,*])
         y = reform(modtemp[p,*])

         index=where(finite(x) AND finite(y),count)

         IF count GT 0 THEN BEGIN
            stdx=stddev(x[index])
            stdy=stddev(y[index])
            varx=stdx^2
            vary=stdy^2
            meany=mean(y[index])
            meanx=mean(x[index])
            bias=meany-meanx
            mab=mean(abs(y[index]-x[index]))
            cp_rms=sqrt( total(((y[index]-meany) - (x[index]-meanx))^2 )/float(count-1) )       
            rms=sqrt(total((y[index] - x[index])^2 )/float(count-1) )       
            stdnorm=stdy/stdx
            rval = correlate(x[index],y[index])

;            tm = TM_TEST(x[index],y[index])
;            fv = FV_TEST(x[index],y[index]) 
;            rs = RS_TEST(x[index],y[index]) 
          
            sig_diff_bias = abs(bias) GT 0.05*(range[1]-range[0])
            sig_diff_rms = rms GT 0.05*(range[1]-range[0])
            sig_diff_corr = 1 - corr_ttest(corr=rval,N=count,sl=0.05)

            print,p+1, ': ',experiments[0],': ',varnames[0]
            print,'N: ',count,' sd/sd: ',stdnorm,' r: ',rval,' Bias: ',bias,' cpRMS: ',cp_rms, ' RMS: ',rms
            print,'MAB: ',mab
            print,'R: ',rval,' RMS: ',rms,format='(A,f5.2,A,f5.1)'
            print,'stdx: ',stdx,' stdy: ',stdy
            print,'varx: ',varx,' vary: ',vary
            print,'meanx: ',meanx,' meany: ',meany
            print,'significant bias: ',sig_diff_bias
            print,'significant rms:  ',sig_diff_rms
            print,'significant corr: ',sig_diff_corr
            print,' '

            text=text+' & '
            IF ~sig_diff_bias THEN text = text + '\textbf{'
            text=text+strtrim(string(bias,format='(f5.2)'),2)
            IF ~sig_diff_bias THEN text = text + '}'

            text=text+' ('

            IF ~sig_diff_rms THEN text = text + '\textbf{'
            text=text+strtrim(string(rms,format='(f5.2)'),2)
            IF ~sig_diff_rms THEN text = text + '}'

;            text=text+', '
;
;            IF ~sig_diff_corr THEN text = text + '\textbf{'
;            text=text+strtrim(string(rval,format='(f5.2)'),2)
;            IF ~sig_diff_corr THEN text = text + '}'
      
            text=text+')'

;;            text = text + ' \\'

            print,text
            printf,lun,text
 
         ENDIF
      ENDFOR

      free_lun,lun

  ENDIF
  
  stop
  
  END
