PRO netcdf_write,ncid,varid,data,offset=OFFSET,count=COUNT

  IF (n_elements(offset) GT 0) AND (n_elements(count) GT 0) THEN BEGIN

     NCDF_VARPUT,ncid,varid,data,count=count,offset=offset

  ENDIF ELSE BEGIN

     NCDF_VARPUT,ncid,varid,data

  ENDELSE

END
