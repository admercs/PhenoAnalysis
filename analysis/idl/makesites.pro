PRO makesites
  
;; this script creates a sites.dat file by pointing on a world map and
;; selecting lon/lat center points

  ;; flags
  writeall = 1   ;; write all cursor position, don't ask

  ;; directories
;  datadir = '/project/msclim/stockli'
  datadir = '/Users/stockli'

  pftdir = datadir+'/pft_sumcrops/'
  pftfile = 'pft.global_0.1.nc'

  biomesoildir = datadir+'/biome_soil/'
  biomesoilfile = 'biome_soil.global.nc'

  topodir = datadir+'/topo/'
  topofile = 'topo.global_0.1.nc'

  sitedir  = '../../sites/'
  sitefile = 'sites_xxx.dat'

  ;; average pft and elevation distribution over grid size
  dlon = 0.5
  dlat = 0.5

  ;; maximum number of sites
  maxsites = 2048

  ;; allocate arrays
  sitename = strarr(maxsites)
  siteshort = strarr(maxsites)
  plon = fltarr(maxsites)
  plat = fltarr(maxsites)
  phgt = intarr(maxsites)
  pbiome = intarr(maxsites)
  psoil = intarr(maxsites)

  ;; initialize constants
  thgt = 0
  psoilname = 'NA'
  pbiomename = 'NA'
  dtsite = 360
  

  ;; reset site counter
  ns = -1

  ;; open sitefile or append to site file
  ret=file_search(sitedir+sitefile,count=count)
  IF count NE 0 THEN BEGIN
     temp=''
     print,'Site file already exists. Appending'

     ;; read current file content
     openr,lun,sitedir+sitefile,/get_lun,width=300
     readf,lun,temp
     readf,lun,temp
     WHILE NOT EOF(lun) DO BEGIN
        ns += 1
        readf,lun,temp
        temp2=strsplit(temp,",",/extract)
        siteshort[ns] = temp2[1]
        sitename[ns] = temp2[2]
        plon[ns] = temp2[3]
        plat[ns] = temp2[4]
        phgt[ns] = temp2[5]
        pbiome[ns] = temp2[8]
        psoil[ns] = temp2[10]        
     ENDWHILE
     free_lun,lun
  ENDIF ELSE BEGIN
     print,'New site file.'
  ENDELSE

  openw,lun,sitedir+sitefile,/get_lun,width=300  
  printf,lun,'ID   Short  Site      Longitude Latitude Elev   Height    Biome Type   SiB Nr.    Soil Type    FAO Nr.    dt'
  printf,lun,'                      [deg]     [deg]    [m]    [m]                                                       [min]'

  FOR s = 0,ns DO BEGIN
     printf,lun,s+1,siteshort[s],sitename[s],plon[s],plat[s],phgt[s], $
            thgt,pbiomename,pbiome[s],psoilname,psoil[s],dtsite, $
            format='(I4,",",A6,",",A12,",",F12.4,",",F12.4,",",I8,",",I8,",",A10,",",I4,",",A10,",",I4,",",I4)'
  ENDFOR

  ;; advance to new position
  ns += 1

  ;; display xy extent
  nx = 1200
  ny = 600

  ;; set up graphics
  device,retain=2
  device,decomposed=0
  window,0,xsize=nx,ysize=ny
  wset,0
  loadct,10,file='brewer.tbl'
  tvlct,red,green,blue,/get
  red[0]=0
  green[0]=0
  blue[0]=0
  red[255]=20
  green[255]=200
  blue[255]=255
  tvlct,red,green,blue
  
  ;; read pft's
  ncid=ncdf_open(pftdir+pftfile)
  
  dimid = ncdf_dimid(ncid,'lon')
  ncdf_diminq,ncid,dimid,name,nlon
  dimid = ncdf_dimid(ncid,'lat')
  ncdf_diminq,ncid,dimid,name,nlat

  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat

  ret = ncdf_inquire(ncid)

  npft = ret.nvars - ret.ndims

  pft = bytarr(nlon,nlat,npft)
  pftnames = strarr(npft)
  pftlongnames = strarr(npft)

  FOR p = 0, npft-1 DO BEGIN

     ncdf_varget,ncid,p+ret.ndims,temp
     pft[*,*,p] = temp

     temp = ncdf_varinq(ncid,p+ret.ndims)
     pftnames[p] = temp.name

     ncdf_attget,ncid,p+ret.ndims,'long_name',temp
     pftlongnames[p] = string(temp)


  ENDFOR

  ncdf_close,ncid

  ; read biome and soil maps (assume same resolution)
  ncid = ncdf_open(biomesoildir+biomesoilfile)
  ncdf_varget,ncid,'lc',biome
  ncdf_varget,ncid,'tex_sub',soil
  ncdf_close,ncid

  ; read elevation map (assume same resolution)
  ncid = ncdf_open(topodir+topofile)
  ncdf_varget,ncid,'Z',hgt
  ncdf_close,ncid

  keyname = ''
  quit = 0
  p = 0

  WHILE quit EQ 0 DO BEGIN

     print,'displaying pft: ',p+1,' : ',pftlongnames[p]

     ;; display pft map
     tv,bytscl(congrid(reform(pft[*,*,p]),nx,ny),min=0,max=100,top=254)

;;     sitesel = [4,6,9,13,17,21,25,29,33,37,43,47,49,55,57,61]

     ;; display current site locations
     FOR s=0,ns-1 DO BEGIN

;        idx = where(sitesel EQ (s+1),count)
;        IF (count GT 0) THEN BEGIN
           x = round((plon[s] + 180.)/360. * float(nx))
           y = round((plat[s] + 90.)/180.  * float(ny))
;           print,s,plon[s],plat[s],x,y
           plots,x,y,psym=6,symsize=2,thick=2,color=255,/device
;        ENDIF
     ENDFOR

     WHILE keyname EQ '' DO BEGIN
     
        keyname = GET_KBRD(0,/key_name)
        
        cursor, cx, cy, 0, /device
        
        IF (!mouse.button eq 1) THEN BEGIN
           
           lon = (((!mouse.x/float(nx) )>0.)<0.99999) * 360. - 180.
           lat = (((!mouse.y/float(ny) )>0.)<0.99999) * 180. - 90.

           ;; ROUND to 0.5 degrees
           plon[ns] = float(round(lon*4))/4.
           plat[ns] = float(round(lat*4))/4.
           
           ppft = bytarr(npft)
           FOR ipft = 0, npft-1 DO BEGIN
              ppft[ipft] = mean(pft[(plon[ns]-0.5*dlon+180.)/360.*nlon:(plon[ns]+0.5*dlon+180.)/360.*nlon, $
                              (plat[ns]-0.5*dlat+90.)/180.*nlat:(plat[ns]+0.5*dlat+90.)/180.*nlat,ipft])
           ENDFOR

           maxpft = max(ppft,pftindex)

           pftname = pftnames[pftindex]
           
           sitename[ns] = pftname + string(ns+1,format='(I4.4)')
           siteshort[ns] = 'G'+ string(ns+1,format='(I4.4)')

           pbiome[ns] = biome[(plon[ns]+180.)/360.*nlon,(plat[ns]+90.)/180.*nlat]
           psoil[ns] = soil[(plon[ns]+180.)/360.*nlon,(plat[ns]+90.)/180.*nlat]
           phgt[ns] = round(mean(hgt[(plon[ns]-0.5*dlon+180.)/360.*nlon:(plon[ns]+0.5*dlon+180.)/360.*nlon, $
                           (plat[ns]-0.5*dlat+90.)/180.*nlat:(plat[ns]+0.5*dlat+90.)/180.*nlat]))

           print,ns+1, plon[ns], plat[ns],' Z: ',phgt[ns],' m maxpft: ',pftname,' pfts %: ', ppft


           IF writeall THEN BEGIN
              ;; write each site except those that have too much water
              ;; area

              IF ppft[npft-1] LT 50. THEN BEGIN

                 printf,lun,ns+1,siteshort[ns],sitename[ns],plon[ns],plat[ns],phgt[ns], $
                        thgt,pbiomename,pbiome[ns],psoilname,psoil[ns],dtsite, $
                        format='(I4,",",A6,",",A12,",",F12.4,",",F12.4,",",I8,",",I8,",",A10,",",I4,",",A10,",",I4,",",I4)'
                 
                 ;; display new site location
                 x = round((plon[ns] + 180.)/360. * float(nx))
                 y = round((plat[ns] + 90.)/180.  * float(ny))
                 plots,x,y,psym=6,symsize=2,thick=2,color=255,/device

                 ns +=1
              ENDIF ELSE BEGIN
                print,'Too much water in this area, not saving it ...'
             ENDELSE
           ENDIF
           
        ENDIF

        wait, 0.25

     END

     IF keyname EQ 'LEFT' THEN p=p-1
     IF keyname EQ 'RIGHT' THEN p=p+1

     IF NOT writeall THEN BEGIN
        IF keyname EQ 'DOWN' THEN BEGIN
           ;; write new entry in site file
           
           IF pftname EQ pftnames[p] THEN BEGIN
              
              printf,lun,ns+1,siteshort[ns],sitename[ns],plon[ns],plat[ns],phgt[ns], $
                     thgt,pbiomename,pbiome[ns],psoilname,psoil[ns],dtsite, $
                     format='(I4,",",A6,",",A12,",",F12.4,",",F12.4,",",I8,",",I8,",",A10,",",I4,",",A10,",",I4,",",I4)'
              
              ns +=1
              
           ENDIF ELSE BEGIN
              
              print,"Maximum %PFT of this site does not correspond to displayed PFT map."
              
           ENDELSE
        ENDIF
     ENDIF

     IF keyname EQ 'UP' THEN quit=1

     keyname = ''

     IF p LT 0 THEN p=0
     IF p GE npft THEN p=npft-1

  END

  ;; close site file
  free_lun,lun



  stop

END
