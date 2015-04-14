PRO netcdf_init,filename,version,dims,vars,global,ncid,dimid,dimvarid,varid

  ;; Code which initializes a NetCDF file when provided with the
  ;; set-up of dimensions, units and variable names. Returns NetCDF
  ;; file ID, dimension id's and variable id's for later writing data
  ;; to this file
  
  ;; todo: test with multiple lon/lat sizes (e.g. multiple resolutions
  ;; in one file

  ;; 27 August 2008 Reto Stockli
  ;; 02 September 2009 RES: added char+short+long integer support

  ndim = (size(dims,/dimensions))[0]
  nvar = (size(vars,/dimensions))[0]
  nglobal = (size(global,/dimensions))[0]

  ndimvar = 0
  FOR d = 0,ndim-1 DO BEGIN
     seldimid = where(dims[d].dims NE 0,dimcount)
     ndimvar = ndimvar+(dimcount GT 0)
  ENDFOR

  ncid = -1
  dimid = replicate(-1,ndim)
  dimvarid = replicate(-1,ndim)
  varid = replicate(-1,nvar)

  ;; calculate current time, format for NetCDF output
  caldat,systime(/julian,/utc),month,day,year,hour,minute,second
  
  prod_date = string(year,format='(I4)')+'-'+string(month,format='(I2.2)')+'-'+ $
              string(day,format='(I2.2)')+' '+string(hour,format='(I2.2)')+':'+ $
              string(minute,format='(I2.2)')+':'+string(round(second),format='(I2.2)')

  ;; check for the file, if it exists, open it for appending and first read dimension and
  ;; variable id's else create one

  exist = 1B
  temp = file_search(filename)

  IF strtrim(temp) EQ '' THEN BEGIN
     print,filename+' does not exist; creating new one.'
     exist = 0B
  ENDIF
 
  IF exist NE 0B THEN BEGIN

     ncid = ncdf_open(filename,/write)

     info = ncdf_inquire(ncid)

     ;; check number of dimensions and variables
     IF ((info.ndims NE ndim) OR (info.nvars NE (nvar+ndimvar)) OR $
         (info.ngatts NE nglobal)) THEN exist = 0B

     IF exist NE 0B THEN BEGIN

        ;; get dimension id's
        FOR d=0,ndim-1 DO BEGIN
           dimid[d] = ncdf_dimid(ncid,dims[d].name)

           IF dimid[d] LT 0 THEN BEGIN
               exist = 0B
           ENDIF ELSE BEGIN
               ncdf_diminq,ncid,dimid[d],name,size
               IF size NE dims[d].size THEN exist = 0B
           ENDELSE

           seldimid = where(dims[d].dims NE 0,dimcount)
           IF dimcount GT 0 THEN BEGIN
              dimvarid[d] = ncdf_varid(ncid,dims[d].name)
              IF dimvarid[d] LT 0 THEN exist = 0B
           ENDIF

        ENDFOR

        ;; get variable id's
        FOR v=0,nvar-1 DO BEGIN
           varid[v] = ncdf_varid(ncid,vars[v].name)
           IF varid[v] LT 0 THEN exist = 0B
        ENDFOR

     ENDIF

     ;; check for completeness of existing file 
     IF exist EQ 0B THEN BEGIN
        print,filename+' exists but is incomplete; creating new one.'
        ncdf_close,ncid
     ENDIF ELSE BEGIN
        print,filename+' exists and is healthy; opening for appending/updating.'
     ENDELSE

  ENDIF

  IF exist EQ 0 THEN BEGIN
     ;; file does not exist or is incomplete, so create new one

     ;; create NetCDF file (header and dims etc.)
     ncid = NCDF_CREATE(filename,/CLOBBER)
     
     ;; enter NetCDF control mode
     NCDF_CONTROL,ncid,/FILL
     ;NCDF_CONTROL,ncid,/NOFILL
     ;NCDF_CONTROL,ncid,/VERBOSE
     
     ;; define dims
     FOR d=0,ndim-1 DO BEGIN
        IF (dims[d].size EQ -1) THEN $
           dimid[d] = NCDF_DIMDEF(ncid,dims[d].name,/UNLIMITED) $
        ELSE $
           dimid[d] = NCDF_DIMDEF(ncid,dims[d].name,dims[d].size)

        IF dimid[d] LT 0 THEN BEGIN
           print,'Error with dimension ',dims[d].name,'. Aborting ...'
           stop
        ENDIF
     ENDFOR
     
     ;; define dimension variables
     FOR d=0,ndim-1 DO BEGIN
        
        seldimid = where(dims[d].dims NE 0,dimcount)
        
        IF dimcount GT 0 THEN BEGIN
           ;; this dimension has a variable associated with it
           
           CASE dims[d].type OF
              'float' : dimvarid[d] = NCDF_VARDEF(ncid, dims[d].name,dimid[seldimid],/FLOAT)
              'double' : dimvarid[d] = NCDF_VARDEF(ncid, dims[d].name,dimid[seldimid],/DOUBLE)
              'long' : dimvarid[d] = NCDF_VARDEF(ncid, dims[d].name,dimid[seldimid],/LONG)
              'integer' : dimvarid[d] = NCDF_VARDEF(ncid, dims[d].name,dimid[seldimid],/SHORT)
              'byte' : dimvarid[d] = NCDF_VARDEF(ncid, dims[d].name,dimid[seldimid],/BYTE)
              'char' : dimvarid[d] = NCDF_VARDEF(ncid, dims[d].name,dimid[seldimid],/CHAR)
               ELSE : BEGIN
                 print, "Error 400: NetCDF type ', dims[d].type, 'not yet implemented"
                 stop
              END
           ENDCASE
           
           IF dimvarid[d] LT 0 THEN BEGIN
              print,'Error with dimension variable ',dims[d].name,'. Aborting ...'
              stop
           ENDIF
           
           NCDF_ATTPUT,ncid,dimvarid[d],'long_name',dims[d].longname
           NCDF_ATTPUT,ncid,dimvarid[d],'units',dims[d].units
           
           IF dims[d].name EQ 'time' THEN NCDF_ATTPUT, ncid,dimvarid[d], 'calendar','gregorian'
           
        ENDIF
        
     ENDFOR
     
     ;; define data variables
     FOR v=0,nvar-1 DO BEGIN
        
        seldimid = where(vars[v].dims NE 0,dimcount)
        
        IF dimcount GT 0 THEN BEGIN
           ;; this variable has one or several dimensions
           
           CASE vars[v].type OF
              'float' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,dimid[seldimid],/FLOAT)
              'double' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,dimid[seldimid],/DOUBLE)
              'long' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,dimid[seldimid],/LONG)
              'integer' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,dimid[seldimid],/SHORT)
              'byte' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,dimid[seldimid],/BYTE)
              'char' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,dimid[seldimid],/CHAR)
              ELSE : BEGIN
                 print, "Error 400: NetCDF type ', vars[v].type, 'not yet implemented"
                 stop
              END
           ENDCASE
        ENDIF ELSE BEGIN
           ;; this variable is a simple scalar
        
            CASE vars[v].type OF
              'float' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,/FLOAT)
              'double' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,/DOUBLE)
              'long' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,/LONG)
              'integer' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,/SHORT)
              'byte' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,/BYTE)
              'char' : varid[v] = NCDF_VARDEF(ncid, vars[v].name,/CHAR)
              ELSE : BEGIN
                 print, "Error 400: NetCDF type ', vars[v].type, 'not yet implemented"
                 stop
              END
           ENDCASE
        ENDELSE
           
        IF varid[v] LT 0 THEN BEGIN
           print,'Error with variable ',vars[v].name,'. Aborting ...'
           stop
        ENDIF
        
        NCDF_ATTPUT,ncid,varid[v],'long_name',vars[v].longname
        NCDF_ATTPUT,ncid,varid[v],'units',vars[v].units
        CASE vars[v].type OF
           'float':   NCDF_ATTPUT,ncid,varid[v],'missing_value',float(vars[v].missing)
           'double':  NCDF_ATTPUT,ncid,varid[v],'missing_value',double(vars[v].missing)
           'long':    NCDF_ATTPUT,ncid,varid[v],'missing_value',long(vars[v].missing)
           'integer': NCDF_ATTPUT,ncid,varid[v],'missing_value',fix(vars[v].missing)
;           'byte':
;           NCDF_ATTPUT,ncid,varid[v],'missing_value',byte(vars[v].missing)
           ;; fix for now:
           'byte':    NCDF_ATTPUT,ncid,varid[v],'missing_value',-127
           'char':    NCDF_ATTPUT,ncid,varid[v],'missing_value',string(vars[v].missing)
           ELSE : BEGIN
              print, "Error 403: NetCDF type ', vars[v].type, 'not yet implemented"
              stop
           END
        ENDCASE

;        NCDF_ATTPUT,ncid,varid[v],'_FillValue',vars[v].missing ;;
;        does not seem to work in IDL?
        NCDF_ATTPUT,ncid,varid[v],'version',version
        NCDF_ATTPUT,ncid,varid[v],'prod_date',prod_date
        ;; NCDF_ATTPUT,ncid,varid[v],'axis',vars[v].axis
        ;; NCDF_ATTPUT,ncid,varid[v],'scale_factor',vars[v].scale
        ;; NCDF_ATTPUT,ncid,varid[v],'add_offset',vars[v].offset
        
     ENDFOR
     
     ;; define global attributes
     FOR g=0,nglobal-1 DO BEGIN
        
        NCDF_ATTPUT, ncid, /GLOBAL, global[g].name,global[g].text
        
     ENDFOR
     
     NCDF_CONTROL, ncid, /ENDEF ; Put file in data mode. 

  ENDIF

END
