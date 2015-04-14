PRO convert_sites
  
  
  infile='../../sites/sites_selected.dat'
  outfile='../../sites/selected.dat'
  nsites=22
  dlon = 0.0
  dlat = 0.0

  sitetype = "grid"
  projection = "lonlat"
  projparam = "NA"


  ;; read old site file
  stemp=''

  names=strarr(nsites)
  shorts=strarr(nsites)
  lats=fltarr(nsites)
  lons=fltarr(nsites)
  elevations=intarr(nsites)

  openr,lun,infile,/get_lun
  readf,lun,stemp
  readf,lun,stemp
  FOR s=0,nsites-1 DO BEGIN
     readf,lun,stemp

     arr=strsplit(stemp,',',/extract)
     shorts[s]=strtrim(arr[1],2)
     names[s]=strtrim(arr[2],2)
     lons[s]=float(arr[3])
     lats[s]=float(arr[4])
     elevations[s]=fix(arr[5])
  ENDFOR
  free_lun,lun

  ;; write new site file
  openw,lun,outfile,/get_lun
  stemp="# This is the site / region geographic boundary conditions file. Columns are separated by commas. For site locations xrange/yrange are scalars. "
  printf,lun,stemp
  stemp="No.,    Abbr.,  Name,                           xmin,           xmax,           ymin,           ymax,           deltax,         deltay,         Sitetype,       Projection,     Projparam,      Elevation"
  printf,lun,stemp
  FOR s=0,nsites-1 DO BEGIN
     printf,lun,s+1,shorts[s],names[s],lons[s]-0.5*dlon,lons[s]+0.5*dlon,lats[s]-0.5*dlat,lats[s]+0.5*dlat,$
            dlon,dlat,sitetype,projection,projparam,elevations[s], $
            format='(I3,",",4X,A5,",",2X,A30,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",A8,",",A8,",",A8,",",I8)'
  ENDFOR

  free_lun,lun

  stop
END

