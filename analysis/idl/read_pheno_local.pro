PRO read_pheno_local, outputdir, model, experiment, sitename, year, varname, dtmod, modvar, modvarerror, $
                      modvartitle, modvarunits, selpft, selhgt, modelshort,modelname, nodata, has_model

  syear=''
  reads,year,syear
  syear=strtrim(syear,2)

  lv       = 2.52*1.e6          ; latent heat of vaporization (J/kg)
  porosity = 0.5                ; assumed porosity if there is no information (-)

  monsum_noleap=[0,31,59,90,120,151,181,212,243,273,304,334,365]
  monsum_leap=[0,31,60,91,121,152,182,213,244,274,305,335,366]
  
  leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))

;  model = 3

; this is the filename syntax of the individual models
  modelprefix=['analysis','simulator','observations']
  modelsuffix=['nc','nc','nc']
  modelnames =['Analysis','Simulator','Observations']

  modelname = modelnames[model-1]
  modelshort = modelname

  has_model = 0B
  has_error = 0B
  first = 1B

  infiles=file_search(outputdir+experiment+'/output/'+ sitename + '.' + $
                      modelprefix[model-1] + '.' + syear + '*.' + modelsuffix[model-1])

  IF infiles[0] NE '' THEN BEGIN
     FOR f=0,n_elements(infiles)-1 DO BEGIN
        infile=infiles[f]
        ncid=NCDF_OPEN(infile)
        print,'Reading Model Data: ',infile

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
           'seconds': timestepmod=(modtime[1]-modtime[0])
           'hours': timestepmod=3600.*(modtime[1]-modtime[0])
           'days': timestepmod=24.*3600.*(modtime[1]-modtime[0])
           ELSE:
        ENDCASE

        IF finite(dtmod,/nan) THEN dtmod=float(round(timestepmod))

        varid=NCDF_VARID(ncid,varname)

        IF varid GE 0L THEN BEGIN

           has_model=1

           NCDF_VARGET,ncid,varid,vart

           IF size(vart,/type) EQ 1L THEN is_byte = 1B ELSE is_byte = 0B

           ;; inquire attributes
           varinfo = ncdf_varinq(ncid,varid)
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

           NCDF_ATTGET,ncid,varid,'long_name',tmp
           modvartitle=strtrim(string(tmp),2)
           NCDF_ATTGET,ncid,varid,'units',tmp
           modvarunits=strtrim(string(tmp),2)

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

           ;; check if this is a 1D model output
           idx =  where(dimnames[varinfo.dim] EQ 'lon',count)
           IF count GT 0 THEN BEGIN
              IF dimsizes[varinfo.dim[idx]] NE 1 THEN BEGIN
                 print,"This code is not suitable to display 2D model output. Stopping"
                 print,"Please set $nlon to 1 and rerun the data assimilation"
                 stop
              ENDIF
           ENDIF

           idx =  where(dimnames[varinfo.dim] EQ 'lat',count)
           IF count GT 0 THEN BEGIN
              IF dimsizes[varinfo.dim[idx]] NE 1 THEN BEGIN
                 print,"This code is not suitable to display 2D model output. Stopping"
                 print,"Please set $nlat to 1 and rerun the data assimilation"
                 stop
              ENDIF
           ENDIF

           ;; check if the data is stored per hgt
           ;; then average or select specific hgt's
           idx = where(dimnames[varinfo.dim] EQ 'Z',count)
           IF count GT 0 THEN BEGIN
              print,"The code is not ready to analyze individual subgrid-elevation patches. Stopping."
              print,"Solution: please run the data assimilation with $averagehgt = 1; "
              stop
           ENDIF

           ;; check if a separate model error (standard deviation) is
           ;; also present in the file
           varid=NCDF_VARID(ncid,varname+'_stddev')

           IF varid GE 0 THEN BEGIN

              NCDF_VARGET,ncid,varid,varerrt
              
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
              
              has_error = 1B
           ENDIF

           ;; check if the data is stored per pft
           ;; then average or select specific pft's
           idx = where(dimnames[varinfo.dim] EQ 'pft',count)
           IF count GT 0 THEN BEGIN
              idx=idx[0]
              IF (selpft EQ 0) OR (dimsizes[varinfo.dim[idx]] EQ 1) THEN BEGIN
                 vart=total(vart,idx+1,/nan)/float(dimsizes[varinfo.dim[idx]])
              ENDIF ELSE BEGIN
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

              ENDELSE
           ENDIF

           vart = reform(vart)
           
           IF has_error THEN varerrt = reform(varerrt)

           ndim = n_elements(size(vart,/dimensions))

           IF first THEN BEGIN
              ntmod=round(ndays*3600./dtmod*24.)

              ndim = n_elements(size(vart,/dimensions))
              IF (ndim GT 1) AND ((size(vart,/dimensions))[0] GT 1) THEN BEGIN
                 npmod=(size(vart,/dimensions))[0]
              ENDIF ELSE  BEGIN
                 npmod = 1
              ENDELSE

              modvar=make_array(ntmod,/float,value=!values.f_nan)
              modvarerror=make_array(ntmod,/float,value=!values.f_nan)
           ENDIF

           ;; put model values on correct time axis      
           time0=float((day0-1)+monsum(month0-1)+ndays*(year0-year))*24.*3600./dtmod

           CASE timeunits[0] OF
              'seconds': modtime2=long(modtime/dtmod + time0 - 0.1)
              'hours': modtime2=long(modtime*3600./dtmod + time0 - 0.1)
              'days':  modtime2=long(modtime*24.*3600./dtmod + time0 - 0.1)
              ELSE: BEGIN
                 print,"Unknown model time units. stopped."
                 stop
              END
           ENDCASE

           index=where((modtime2 GE 0) AND (modtime2 LT ntmod),timecount)
           IF timecount GT 0 THEN BEGIN
              IF npmod EQ 1 THEN modvar[modtime2[index]]=vart[index] ELSE $ 
                 modvar[modtime2[index]]=total(vart[*,index],1)/float(npmod)

              IF has_error THEN  modvarerror[modtime2[index]]=varerrt[index]
           ENDIF

        ENDIF ELSE BEGIN
           IF first THEN BEGIN
              ntmod=round(ndays*3600./dtmod*24.)
              npmod=1
              modvar=make_array(ntmod,/float,value=!values.f_nan)
              modvarerror=make_array(ntmod,/float,value=!values.f_nan)
           ENDIF
        ENDELSE

        NCDF_CLOSE,ncid
        first=0B                 ; not first file
     ENDFOR                     ; number of files          

  ENDIF                         ; if there is a file to read 
END
