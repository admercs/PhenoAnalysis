PRO read_pheno_regional,outputdir,model,experiment,sitename,year,varname,dtmod, $
                        modvar,modvarerror,modvarname,modvarunit,modvarlongname, modlons, modlats, $
                        modelshort,modelname, nodata, has_model, nlonmod, nlatmod, $
                        pftave, hgtave, selpft, selhgt, pft, topo, $
                        lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax

  syear=''
  reads,year,syear
  syear=strtrim(syear,2)
  
  IF keyword_set(lonmin) OR keyword_set(lonmax) OR keyword_set(latmin) OR keyword_set(latmax) THEN $
     has_grid=1 ELSE has_grid=0


  monsum_noleap=[0,31,59,90,120,151,181,212,243,273,304,334,365]
  monsum_leap=[0,31,60,91,121,152,182,213,244,274,305,335,366]
  
  leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
  
  ;; read gridded observations (saveobs):
  ;; model = 3

  ;; this is the filename syntax of the individual models
  modelprefix=['analysis','simulator','observations']
  modelsuffix=['nc','nc','nc']
  modelnames =['Analysis','Simulator','Observations']

  modelname = modelnames[model-1]
  modelshort = modelname

  infiles=file_search(outputdir+experiment+'/output/'+ sitename + '.' + $
                      modelprefix[model-1] + '.' + syear + '*.' + modelsuffix[model-1])

  has_model=0
  IF infiles[0] NE '' THEN BEGIN
     first=1
     FOR f=0,n_elements(infiles)-1 DO BEGIN
        infile=infiles[f]
        ncid=NCDF_OPEN(infile)
        print,'Reading Model Data: ',infile
        ;; get time axis
        varid=NCDF_VARID(ncid,'time')
        NCDF_VARGET,ncid,varid,modtime
        NCDF_ATTGET,ncid,varid,'units',timeunits
        NCDF_ATTGET,ncid,varid,'calendar',calendar
        timeunits=strtrim(string(timeunits),2)
        calendar=strtrim(string(calendar),2)
        timeunits=strsplit(timeunits,' ',/extract)
        day0=0
        month0=0
        year0=0
        date0=strsplit(timeunits[2],'-',/extract)
        reads,date0[0],year0
        reads,date0[1],month0
        reads,date0[2],day0
        
        ;; check for calendar in model output
        ;; most models have only 365 day years
        IF (calendar EQ 'gregorian') AND leapyear THEN BEGIN
           monsum = monsum_leap 
           ndays = 366
        ENDIF ELSE BEGIN
           monsum = monsum_noleap
           ndays = 365
        ENDELSE

        CASE timeunits[0] OF
           'hours': timestepmod=3600.*(modtime[1]-modtime[0])
           'days': timestepmod=24.*3600.*(modtime[1]-modtime[0])
           ELSE:
        ENDCASE
        
        IF finite(dtmod,/nan) THEN dtmod=float(round(timestepmod))
        
        ;; get lons/lats
        varid=NCDF_VARID(ncid,'x')
        IF (varid LT 0) THEN varid=NCDF_VARID(ncid,'lon')
        NCDF_VARGET,ncid,varid,modlons
        varid=NCDF_VARID(ncid,'y')
        IF (varid LT 0) THEN varid=NCDF_VARID(ncid,'lat')
        NCDF_VARGET,ncid,varid,modlats
        
        ;; inquire file
        info=ncdf_inquire(ncid)
        ndims=info.ndims
        dimnames = strarr(info.ndims)
        dimsizes = intarr(info.ndims)
        FOR d=0,info.ndims-1 DO BEGIN
           ncdf_diminq,ncid,d,n,s
           dimnames[d]=n
           dimsizes[d]=s
        ENDFOR  

        modvarname = varname
        
        ;; inquire actual data
        varid=NCDF_VARID(ncid,varname)

        IF varid GE 0L THEN BEGIN

           has_model=1

           varinfo = ncdf_varinq(ncid,varid)

           ;; set up dimensions on first read (for model
           ;; simulations with many files per year)
           IF first THEN BEGIN

              ;; read pft distribution
              ncdf_varget,ncid,"PFT_PCT",pft

              ;; read topography distribution
              ncdf_varget,ncid,"Z",topo

              ;; inquire attributes
              has_fillvalue = 0B
              has_add_offset = 0B
              has_scale_factor = 0B
              FOR attnum=0,varinfo.natts-1 DO BEGIN
                 attname = ncdf_attname(ncid,varid,attnum)
                 IF attname EQ '_FillValue' THEN BEGIN
                    has_fillvalue = 1B
                    ncdf_attget,ncid,varid,'_FillValue',fillvalue
                 ENDIF
                 IF attname EQ 'add_offset' THEN BEGIN
                    has_add_offset = 1B
                    ncdf_attget,ncid,varid,'add_offset',add_offset
                 ENDIF 
                 IF attname EQ 'scale_factor' THEN BEGIN
                    has_scale_factor = 1B
                    ncdf_attget,ncid,varid,'scale_factor',scale_factor
                 ENDIF
              ENDFOR
              
              IF varname EQ "PFT_PCT" THEN BEGIN
                 has_scale_factor = 0B
                 has_add_offset = 0B
              ENDIF

              ;; inquire data names and units
              NCDF_ATTGET,ncid,varid,'long_name',tmp
              modvarlongname=strtrim(string(tmp),2)
              NCDF_ATTGET,ncid,varid,'units',tmp
              modvarunit=strtrim(string(tmp),2)

              ;; predefine space/time dimensions
              nlonmod = 1
              nlatmod = 1
              nensmod = 1
              ntmod = 1
              
              ;; check if this is a 1D model output
              idx =  where(dimnames[varinfo.dim] EQ 'lon',count)
              idx2 =  where(dimnames[varinfo.dim] EQ 'lat',count2)
              IF (count GT 0) AND (count2 GT 0) THEN BEGIN
                 IF (dimsizes[varinfo.dim[idx]] EQ 1) AND (dimsizes[varinfo.dim[idx2]] EQ 1) THEN BEGIN
                    print,"WARNING: This code is not suitable to display 1D model output."
                    print,"Please set $nlon and $nlat > 1 and rerun the data assimilation"
                    stop
                 ENDIF
              ENDIF
              
              ;; check if the data is stored per hgt
              ;; then average or select specific hgt's
              idx = where(dimnames[varinfo.dim] EQ 'Z',count)
              IF count GT 0 THEN BEGIN
                 print,"The code is not yet ready to analyze individual subgrid-elevation patches. Stopping."
                 print,"Solution: please run the data assimilation with $averagehgt = 1; "
                 stop
              ENDIF

              ;; inquire space/time dimensions
              FOR d=0,n_elements(varinfo.dim)-1 DO BEGIN

                 CASE dimnames[varinfo.dim[d]] OF
                    'ens' : nensmod = dimsizes[varinfo.dim[d]]
                    'pft' : npftmod = dimsizes[varinfo.dim[d]]
                    'Z'   : nhgtmod = dimsizes[varinfo.dim[d]]
                    'lon' : nlonmod = dimsizes[varinfo.dim[d]]
                    'lat' : nlatmod = dimsizes[varinfo.dim[d]]
                    'x'   : nlonmod = dimsizes[varinfo.dim[d]]
                    'y'   : nlatmod = dimsizes[varinfo.dim[d]]
                    'time': ntmod = dimsizes[varinfo.dim[d]]
                    ELSE: BEGIN
                       print,'Unknown dimension found in NetCDF: ',dimnames[varinfo.dim[d]],dimsizes[varinfo.dim[d]]
                    END
                 ENDCASE

              ENDFOR

              ;; adjust lon/lat extent according to required grid
              IF has_grid THEN BEGIN

                 lonmin_mod = min(modlons,max=lonmax_mod)
                 latmin_mod = min(modlats,max=latmax_mod)
                 dlon = modlons[1] - modlons[0]
                 dlat = modlats[1] - modlats[0]

                 lonmin_mod -= 0.5d0*dlon
                 lonmax_mod += 0.5d0*dlon
                 latmin_mod -= 0.5d0*dlat
                 latmax_mod += 0.5d0*dlat

                 xmin = -1L
                 xmax = -1L
                 ymin = -1L
                 ymax = -1L

;; define precision for getting rid of rounding errors below
                 pr = 1.d-4

                 FOR i=nlonmod-1L,0L,-1L DO IF (round((modlons[i]-0.5d0*dlon)/pr,/L64) LE round(lonmin/pr,/L64)) $
                    AND (xmin EQ -1L) THEN xmin = i
                 FOR i=nlatmod-1L,0L,-1L DO IF (round((modlats[i]-0.5d0*dlat)/pr,/L64) LE round(latmin/pr,/L64)) $
                    AND (ymin EQ -1L) THEN ymin = i
                 FOR i=0L,nlonmod-1L,1L DO IF (round((modlons[i]+0.5d0*dlon)/pr,/L64) GE round(lonmax/pr,/L64)) $
                    AND (xmax EQ -1L) THEN xmax = i
                 FOR i=0L,nlatmod-1L,1L DO IF (round((modlats[i]+0.5d0*dlat)/pr,/L64) GE round(latmax/pr,/L64)) $
                    AND (ymax EQ -1L) THEN ymax = i

                 IF ((xmin EQ (-1L)) OR (ymin EQ (-1L)) OR (xmax EQ (-1L)) OR (ymax EQ (-1L))) THEN BEGIN
                    print,'model phenology does not cover requested area. stopping.'
                    print,'We need: ',latmin,latmax,lonmin,lonmax
                    print,'We get:  ',latmin_mod,latmax_mod,lonmin_mod,lonmax_mod
                    print,xmin,xmax,ymin,ymax
                    stop
                 ENDIF

                 nlonmod = xmax - xmin + 1L
                 nlatmod = ymax - ymin + 1L
                 
              ENDIF ELSE BEGIN
                 xmin = 0
                 xmax = nlonmod-1
                 ymin = 0
                 ymax = nlatmod-1
              ENDELSE

              ntmod=round(float(ndays)*3600./dtmod*24.)

              ;; check if a separate model error (standard deviation) is
              ;; also present in the file
              varid2=NCDF_VARID(ncid,varname+'_stddev')
              
              IF varid2 GE 0 THEN has_error = 1B ELSE has_error = 0B
              
              ;; initialize variables
              modvar=make_array(nlonmod,nlatmod,ntmod,/float,value=!values.f_nan)
              modvarerror=make_array(nlonmod,nlatmod,ntmod,/float,value=!values.f_nan)
              
           ENDIF

           ;; now read variable time by time step (more memory efficient
           ;; with large ensembles)
           
           ncoffset=replicate(0,n_elements(varinfo.dim))
           nccount=dimsizes[varinfo.dim]

           IF has_grid THEN BEGIN
              lonidx = where(dimnames[varinfo.dim] EQ 'lon',haslon)
              ncoffset[lonidx] = xmin
              nccount[lonidx] = nlonmod
              latidx = where(dimnames[varinfo.dim] EQ 'lat',haslat)
              ncoffset[latidx] = ymin
              nccount[latidx] = nlatmod
           ENDIF

           timeidx = where(dimnames[varinfo.dim] EQ 'time',hastime)
           IF ~hastime THEN BEGIN
             ntmod = 1
           ENDIF

           FOR t = 0, ntmod-1 DO BEGIN

              IF hastime THEN BEGIN
                 ncoffset[timeidx] = t
                 nccount[timeidx] = 1
              ENDIF 
 
              pftidx = where(dimnames[varinfo.dim] EQ 'pft',haspft)
              IF haspft THEN BEGIN
                 IF selpft GT 0 THEN BEGIN
                    ncoffset[pftidx] = selpft-1
                    nccount[pftidx] = 1
                 ENDIF ELSE BEGIN
                    ncoffset[pftidx] = 0
                    nccount[pftidx] = npftmod
                 ENDELSE
              ENDIF

              NCDF_VARGET,ncid,varid,vart, offset=ncoffset, count=nccount

              IF size(vart,/type) EQ 1L THEN is_byte = 1B ELSE is_byte = 0B

              ;; we need both scale_factor and add_offset to be present for
              ;; decompressing byte/integer values to floating point values
              IF has_scale_factor AND has_add_offset THEN BEGIN
                 ;; scaled output defaults to floating point precision (ok?)
                 IF has_fillvalue THEN BEGIN
                    ;; only scale valid data
                    idx = where(vart NE fillvalue,count,complement=idx2,ncomplement=count2)
                    IF (is_byte) THEN BEGIN
                       vart = (vart GE 128B)*(float(vart)-256.0) + (vart LE 127B)*float(vart)
                    ENDIF ELSE vart = float(vart)                 
                    IF count GT 0L THEN vart[idx] *= scale_factor
                    IF count GT 0L THEN vart[idx] += add_offset
                    IF count2 GT 0L THEN vart[idx2] = !values.f_nan
                 ENDIF ELSE BEGIN
                    ;; scale all data
                    IF (is_byte) THEN BEGIN
                       vart = (vart GE 128B)*(float(vart)-256.0) + (vart LE 127B)*float(vart)
                    ENDIF ELSE vart = float(vart)                 
                    vart *= scale_factor
                    vart += add_offset
                 ENDELSE
              ENDIF ELSE BEGIN
                 IF has_fillvalue THEN BEGIN
                    badvar=where(vart EQ fillvalue,badcount)
                    vart = float(vart)
                    IF badcount GT 0 THEN vart[badvar]=!values.f_nan
                 ENDIF
              ENDELSE

              IF has_error THEN BEGIN
                 NCDF_VARGET,ncid,varid2,varerrt, offset=ncoffset, count=nccount

                 IF size(varerrt,/type) EQ 1L THEN is_byte = 1B ELSE is_byte = 0B

                 ;; we need both scale_factor and add_offset to be present for
                 ;; decompressing byte/integer values to floating point values
                 IF has_scale_factor AND has_add_offset THEN BEGIN
                    ;; scaled output defaults to floating point precision (ok?)
                    IF has_fillvalue THEN BEGIN
                       ;; only scale valid data
                       idx = where(varerrt NE fillvalue,count,complement=idx2,ncomplement=count2)
                       IF (is_byte) THEN BEGIN
                          varerrt = (varerrt GE 128B)*(float(varerrt)-256.0) + (varerrt LE 127B)*float(varerrt)
                       ENDIF ELSE varerrt = float(varerrt)                 
                       IF count GT 0L THEN varerrt[idx] *= scale_factor
                       IF count GT 0L THEN varerrt[idx] += add_offset
                       IF count2 GT 0L THEN varerrt[idx2] = !values.f_nan
                    ENDIF ELSE BEGIN
                       ;; scale all data
                       IF (is_byte) THEN BEGIN
                          varerrt = (varerrt GE 128B)*(float(varerrt)-256.0) + (varerrt LE 127B)*float(varerrt)
                       ENDIF ELSE varerrt = float(varerrt)                 
                       varerrt *= scale_factor
                       varerrt += add_offset
                    ENDELSE
                 ENDIF ELSE BEGIN
                    IF has_fillvalue THEN BEGIN
                       badvar=where(varerrt EQ fillvalue,badcount)
                       varerrt = float(varerrt)
                       IF badcount GT 0 THEN varerrt[badvar]=!values.f_nan
                    ENDIF
                 ENDELSE
              ENDIF

              ;; check if the data is stored per pft
              ;; then average or select specific pft's
              idx = where(dimnames[varinfo.dim] EQ 'pft',count)
              IF (count GT 0) AND (selpft EQ 0) THEN BEGIN
                 idx=idx[0]
                 vart=total(vart,idx+1,/nan)/float(dimsizes[varinfo.dim[idx]])
              ENDIF ELSE BEGIN
                 goto,bla
                 dimt=dimsizes[varinfo.dim]
                 ndimt=(size(dimt,/dimensions))[0]
                 
                 IF ndimt EQ 5 THEN BEGIN
                    CASE idx OF
                       0: vart=total(vart[selpft-1,*,*,*,*],idx+1,/nan)
                       1: vart=total(vart[*,selpft-1,*,*,*],idx+1,/nan)
                       2: vart=total(vart[*,*,selpft-1,*,*],idx+1,/nan)
                       3: vart=total(vart[*,*,*,selpft-1,*],idx+1,/nan)
                       4: vart=vart[*,*,*,*,selpft-1]
                       ELSE:
                    END
                    IF has_error THEN BEGIN
                       CASE idx OF
                          0: varerrt=total(varerrt[selpft-1,*,*,*,*],idx+1,/nan)
                          1: varerrt=total(varerrt[*,selpft-1,*,*,*],idx+1,/nan)
                          2: varerrt=total(varerrt[*,*,selpft-1,*,*],idx+1,/nan)
                          3: varerrt=total(varerrt[*,*,*,selpft-1,*],idx+1,/nan)
                          4: varerrt=varerrt[*,*,*,*,selpft-1]
                          ELSE:
                       END
                    ENDIF
                 ENDIF

                 IF ndimt EQ 4 THEN BEGIN
                    CASE idx OF
                       0: vart=total(vart[selpft-1,*,*,*],idx+1,/nan)
                       1: vart=total(vart[*,selpft-1,*,*],idx+1,/nan)
                       2: vart=total(vart[*,*,selpft-1,*],idx+1,/nan)
                       3: vart=vart[*,*,*,selpft-1]
                       ELSE:
                    END
                    IF has_error THEN BEGIN
                       CASE idx OF
                          0: varerrt=total(varerrt[selpft-1,*,*,*],idx+1,/nan)
                          1: varerrt=total(varerrt[*,selpft-1,*,*],idx+1,/nan)
                          2: varerrt=total(varerrt[*,*,selpft-1,*],idx+1,/nan)
                          3: varerrt=varerrt[*,*,*,selpft-1]
                          ELSE:
                       END
                    ENDIF
                 ENDIF

                 IF ndimt EQ 3 THEN BEGIN
                    CASE idx OF
                       0: vart=total(vart[selpft-1,*,*],idx+1,/nan)
                       1: vart=total(vart[*,selpft-1,*],idx+1,/nan)
                       2: vart=vart[*,*,selpft-1]
                       ELSE:
                    END
                    IF has_error THEN BEGIN
                       CASE idx OF
                          0: varerrt=total(varerrt[selpft-1,*,*],idx+1,/nan)
                          1: varerrt=total(varerrt[*,selpft-1,*],idx+1,/nan)
                          2: varerrt=varerrt[*,*,selpft-1]
                          ELSE:
                       END
                    ENDIF
                 ENDIF

                 IF ndimt EQ 2 THEN BEGIN
                    CASE idx OF
                       0: vart=total(vart[selpft-1,*],idx+1,/nan)
                       1: vart=vart[*,selpft-1]
                       ELSE:
                    END
                    IF has_error THEN BEGIN
                       CASE idx OF
                          0: varerrt=total(varerrt[selpft-1,*],idx+1,/nan)
                          1: varerrt=varerrt[*,selpft-1]
                          ELSE:
                       END
                    ENDIF
                 ENDIF
                 bla:
              ENDELSE

              vart = reform(vart)
              IF has_error THEN varerrt = reform(varerrt)

              ;; create ensemble standard deviation
              IF (nensmod GT 2) THEN BEGIN
                 varerrt = fltarr(nlonmod,nlatmod)
                 FOR y = 0, nlatmod-1 DO BEGIN
                    FOR x = 0, nlonmod-1 DO BEGIN
                       varerrt[x,y] = stddev(vart[*,x,y],/nan)
                    ENDFOR
                 ENDFOR
                 ;; average ensembles
                 vart = total(vart,1,/nan)/float(nensmod)
              ENDIF

              ;; create mean + stddev from compressed ensemble
              IF (nensmod EQ 2) THEN BEGIN
                 varerrt = reform(vart[1,*,*])
                 vart = reform(vart[0,*,*])
              ENDIF

              ;; put model values on correct time axis      
              time0=float((day0-1)+monsum(month0-1)+float(ndays)*(year0-year))*24.*3600./dtmod
              
              IF (experiment EQ 'regional_orig') OR (experiment EQ 'local_orig') THEN vartdev = vartdev*0.25

              CASE timeunits[0] OF
                 'seconds': modtime2=round(modtime[t]/dtmod + time0)
                 'hours': modtime2=round(modtime[t]*3600./dtmod + time0)
                 'days':  modtime2=round(modtime[t]*24.*3600./dtmod + time0)
                 ELSE: BEGIN
                    print,"Unknown model time units. stopped."
                    stop
                 END
              ENDCASE

              IF (modtime2 GE 0) AND (modtime2 LT ntmod) THEN BEGIN
                 modvar[*,*,modtime2]=vart
                 IF has_error THEN modvarerror[*,*,modtime2]=varerrt
              ENDIF

              vart = 0
              varerrt = 0
              
              first=0           ; not first file
           ENDFOR               ;; number of time steps in file

        ENDIF ELSE has_model = 0

        NCDF_CLOSE,ncid ;; close file
     ENDFOR             ; number of files
  ENDIF                 ; if there is a file to read 

END
