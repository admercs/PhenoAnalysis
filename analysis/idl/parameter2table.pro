PRO parameter2table

  basedir = '../../'

  experiment = 'Global-Prediction-0K'

  selpar = [0,1,2,3,4,5]
;  selpar = [10,11,12,13,14,15]
;  selpar = [17,18,19]

;  selpar = [0,1,2,3,4,5,10,11,12,13,14,15,17,18,19]
  
  sitefile = 'sites_global4.dat'

  csv = 1

  npft = 35
  nparam = 22

  paramname = strarr(nparam)

  paramval = fltarr(npft,nparam)
  paramvariance = fltarr(npft,nparam)

  pftname = strarr(npft)
  pftobs = lonarr(npft)

  temp = ''

  parameterdir = basedir + 'data/parameters/'
  parameterfile = 'prediction_parameters_regional_'+experiment+'.dat'
  
  openr,lun,parameterdir+parameterfile,/get_lun
  readf,lun,temp
  
  val1=0.
  val2=0.
  str1=''
  str2=''
  
  FOR pft = 0,npft-1 DO BEGIN
     readf,lun,temp
     a = strtrim(strmid(temp,5,13),2)
     pos=strpos(a,'_')
     IF pos GE 0 THEN strput,a,' ',pos
     pftname[pft] = a
     pftobs[pft] = pftobs[pft] + long(strmid(temp,50,20))
     FOR param = 0,nparam-1 DO BEGIN
        
        readf,lun,temp
        paramval[pft,param] = float(strmid(temp,0,12))
        paramvariance[pft,param] = float(strmid(temp,19,14))
        a = strtrim(strmid(temp,44,20))
        pos=strpos(a,'_')
        IF pos GE 0 THEN strput,a,' ',pos
        paramname[param] = a
        
        print,paramname[param],paramval[pft,param],paramvariance[pft,param]
        
     ENDFOR
  ENDFOR
  
  free_lun,lun

  IF csv THEN BEGIN
     ;; excel table
     suffix = 'csv'
     separator = ', '
     pm = ' +/-'
     eol = ''

  ENDIF ELSE BEGIN
     ;; latex table
     suffix = 'tex'
     separator = ' & '
     pm = ' $\pm$ '
     eol = ' \\'

  ENDELSE

  openw,lun,'parameters_'+experiment+'.' + suffix,/get_lun

  text = "PFT"
  FOR p=0,n_elements(selpar)-1 DO BEGIN
         param = selpar[p]

      text = text + separator + strtrim(paramname[param],2)

  ENDFOR

  text = text + eol
  printf,lun,text

  paramstddev = sqrt(paramvariance)

  FOR pft=0,npft-1 DO BEGIN
;;     print,pftname[pft],pftobs[pft]
     text = pftname[pft]
     
     FOR p=0,n_elements(selpar)-1 DO BEGIN
         param = selpar[p]
         text=text+separator+strtrim(string(paramval[pft,param],format='(f8.1)'),2)
;         text=text+separator+strtrim(string(paramval[pft,param],format='(f8.1)'),2)+$
;              pm+strtrim(string(paramstddev[pft,param],format='(f8.1)'),2)
;         text=text+separator+strtrim(string(paramval[pft,param],format='(f8.2)'),2)+$
;              pm+strtrim(string(paramstddev[pft,param],format='(f8.2)'),2)
;         text=text+separator+strtrim(string(paramval[pft,param],format='(f8.3)'),2)+$
;              pm+strtrim(string(paramvariance[pft,param],format='(f8.3)'),2)
     ENDFOR

     text = text + eol
     printf,lun,text
  ENDFOR

  free_lun,lun

  stop

END
