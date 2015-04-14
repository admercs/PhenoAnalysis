PRO climate_sites
    
   basedir='/Users/stockli'
;;   basedir='/project/msclim/stockli'
   plotdir='../plots/plots/maps/'   
   sitefile='../sites/sites_global256.dat'    
   plotfile='climate_global_sites_256.eps'

   ;; parameters and flags
   nsites = 256
   landonly = 1

   sites64 = indgen(64)
   sites16 = [4,6,9,13,17,21,25,29,33,37,43,47,49,55,57,61] - 1
   sites4 =  [9,17,29,55] - 1

;; read climate fields (few data, ascii format on a 0.5x0.5 degree
;; global grid). They go from north to south and are centered on 0.25
;; degrees
;; first point is therefore: -179.75W, 89.75N
     prec_file = basedir + '/wilmott_matsuura/precip_corrected.clim'
     temp_file = basedir + '/wilmott_matsuura/lw_temp_dem.clim'

     nmonth = 12
     cnlon = 720
     cnlat = 360
     cdlon = 360./cnlon
     cdlat = 180./cnlat

   lonmin=-180.0
   lonmax= 180.0
   latmin= -90.0
   latmax= 90.0
   map_limit=[latmin,lonmin,latmax,lonmax]

;; create lon/lat arrays (centered on pixel)

   clon = ((findgen(cnlon) + 0.5)*cdlon + lonmin) # replicate(1.0,cnlat) ; W->E
   clat = replicate(1.0,cnlon) # ((findgen(cnlat) + 0.5)*cdlat + latmin) ; S->N
   area = cos(reform(clat[0,*])/180.*!pi) ;; * (2.*6371.0*!pi/float(cnlon))^2

     prec_glob = fltarr(cnlon,cnlat,nmonth)
     temp_glob = fltarr(cnlon,cnlat,nmonth)
     ctemp = fltarr(nmonth+3)

     openr,lun,prec_file,/get_lun
     FOR x=0,cnlon-1 DO BEGIN
         FOR y=0,cnlat-1 DO BEGIN
             readf,lun,ctemp, format = '(2F8.3,13F7.1)'
             prec_glob[x,y,*]=ctemp[2:2+nmonth-1]
         ENDFOR
     ENDFOR
     close,lun
     free_lun,lun

     openr,lun,temp_file,/get_lun
     FOR x=0,cnlon-1 DO BEGIN
         FOR y=0,cnlat-1 DO BEGIN
             readf,lun,ctemp, format = '(2F8.3,13F8.1)'
             temp_glob[x,y,*]=ctemp[2:2+nmonth-1]
         ENDFOR
     ENDFOR
     close,lun
     free_lun,lun

;; arrange climate fields from south to north
     temp_glob = reverse(temp_glob,2)
     prec_glob = reverse(prec_glob,2)


     IF landonly THEN BEGIN
        ;; read pft dataset to screen water areas. water areas receive
        ;; 0% weight in averaging

        region = 'global_0.1'
        read_pftdata, basedir, 'pft', region, clon, clat, pft, dlon=cdlon, dlat=cdlat

        npft = (size(pft,/dimensions))[2]
        landfrac = 1. - reform(pft[*,*,npft-1])/100.

     ENDIF ELSE BEGIN
        landfrac = replicate(1.0,cnlon,cnlat)
     ENDELSE

   ;; open Site file
   stemp=''
   siteshort=''
   sitename=''
   sitelat=0.
   sitelon=0.

   names=strarr(nsites)
   lats=fltarr(nsites)
   lons=fltarr(nsites)
   chosen=replicate(1,nsites)
   dlats=fltarr(nsites)
   dlons=fltarr(nsites)
   align=fltarr(nsites)
   precipitation=fltarr(nsites)
   temperature=fltarr(nsites)

   openr,1,sitefile
   readf,1,stemp
   readf,1,stemp
   FOR s=0,nsites-1 DO BEGIN
       readf,1,stemp

       arr=strsplit(stemp,',',/extract)
       siteshort=strtrim(arr[1],2)
       sitename=strtrim(arr[2],2)
       sitelon=arr[3]
       sitelat=arr[4]
       names[s]=sitename
       lats[s]=sitelat
       lons[s]=sitelon
   ENDFOR
   close,1

   ;; assign mean temperature and mean precipitation to each site
   FOR s=0,nsites-1 DO BEGIN
      x=round((lons[s]+180.)/360. * float(cnlon))
      y=round((lats[s]+90.)/180. * float(cnlat))
      temperature[s] = mean(temp_glob[x,y,*])
      precipitation[s] = mean(prec_glob[x,y,*])
   ENDFOR
   
   FOR s=0,nsites-1 DO BEGIN
      print,s+1,' ',names[s],lons[s],lats[s],temperature[s],precipitation[s]
   ENDFOR
      

   ;; create 2D historgram of global temp/precip weighted by area
   meanprec = total(prec_glob,3)/nmonth
   meantemp = total(temp_glob,3)/nmonth

   numlon = round(cnlon*area)
   idx = long((findgen(numlon[0])+0.5)/float(numlon[0])*cnlon)
   idx2 = where(landfrac[idx,0] GT 0.5,count)
   IF count GT 0 THEN idx3 = idx[idx2]+0L*cnlon
   FOR y=1L,cnlat-1L DO BEGIN
      idx = fix((findgen(numlon[y])+0.5)/float(numlon[y])*cnlon)
      idx2 = where(landfrac[idx,y] GT 0.5,count)
      IF count GT 0 THEN idx3 = [idx3,idx[idx2]+y*cnlon]    
   ENDFOR

   meanprec1d = meanprec[idx3]
   meantemp1d = meantemp[idx3]

   area2d = replicate(1.0,cnlon) # area

   landarea = total(landfrac * area2d)

   sitex=round((lons+180.)/360. * float(cnlon))
   sitey=round((lats+90.)/180. * float(cnlat))
   sitearea = total(landfrac[sitex[sites4],sitey[sites4]]*area2d[sitex[sites4],sitey[sites4]])

   print,sitearea/landarea

   stop

   set_plot,'ps'
   

   DEVICE, /encapsul, bits_per_pixel=8, /COLOR,$
      Filename=plotdir+plotfile, $
      xsize=35, ysize=20,/Helvetica
   !p.font = 0

; Make a vector of 16 points, A[i] = 2pi/16: 
   A = FINDGEN(17) * (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
   USERSYM, COS(A), SIN(A), /FILL 

   yrange=[-20,30]
   xrange=[0,250]
   xstep=5
   ystep=1

   h2d = hist_2d(meanprec1d,meantemp1d,bin1=xstep,bin2=ystep,min1=xrange[0],max1=xrange[1],min2=yrange[0],max2=yrange[1])

   dims= size(h2d,/dimensions)
   nx=dims[0]
   ny=dims[1]

   loadct,45,file='brewer.tbl'
   tvlct,red1,green1,blue1,/get

   loadct,11,file='brewer.tbl'
   tvlct,red,green,blue,/get
 ;  green=shift(green,1)
 ;  blue=shift(blue,1)
 ;  red=shift(red,1)
   red[0]=255B
   green[0]=255B
   blue[0]=255B
   ntable=!d.table_size
   red[ntable-1]=0B
   green[ntable-1]=0B
   blue[ntable-1]=0B
   red[1:4] = red1[0:3]
   green[1:4] = green1[0:3]
   blue[1:4] = blue1[0:3]
   tvlct,red,green,blue

   plot,[0],[0],/nodata,xstyle=1,ystyle=1,xrange=xrange,yrange=yrange, color=ntable-1,$
        xtitle='mean precipitation [mm/month]',ytitle='mean temperature [deg C]', $
        xminor=1,yminor=1,xthick=5.0,ythick=5.0,charsize=1.5,charthick=2.0, $
        xticklen=-0.01*35./20.,yticklen=-0.01

   nlevels = 20
   hrange=[0,150]

   contour,h2d>hrange[0],findgen(nx)/float(nx-1)*(xrange[1]-xrange[0])+xrange[0],$
           findgen(ny)/float(ny-1)*(yrange[1]-yrange[0])+yrange[0], $
           levels=indgen(nlevels)*(hrange[1]-hrange[0])/(nlevels-1)+hrange[0], $
           c_colors=(indgen(nlevels)*(ntable)/(nlevels-1)-2)>5B,/follow,/overplot,/cell_fill
   

   FOR l=0,3 DO BEGIN
      FOR i=0,nsites-1 DO BEGIN

         ;; guarantee that main sites are on top
         ok = 0

         idx = where(sites4 EQ i,count)
         IF count GT 0 THEN BEGIN
            IF l EQ 3 THEN ok = 1
            psym = 6
            color = 1
            symsize = 1.75
            symthick = 8.0
         ENDIF ELSE BEGIN
            symsize = 1.25
            symthick = 5.0
            idx = where(sites16 EQ i,count)
            IF count GT 0 THEN BEGIN
               IF l EQ 2 THEN ok = 1
               psym = 5
               color = 2
            ENDIF ELSE BEGIN
               idx = where(sites64 EQ i,count)
               IF count GT 0 THEN BEGIN
                  IF l EQ 1 THEN ok = 1
                  psym = 4
                  color = 3
               ENDIF ELSE BEGIN
                  IF l EQ 0 THEN ok = 1
                  psym = 8
                  color = 4
               ENDELSE
            ENDELSE
         ENDELSE
         
         IF ok THEN $
            plots,precipitation[i],temperature[i],symsize=symsize,psym=psym,thick=symthick,color=color

      ENDFOR

   ENDFOR
   
   DEVICE,/close
   set_plot,'x'

  ;; spawn,"convert -density 300 -resize 1000 -quality 90 "+plotdir+plotfile+" "+plotdir+plotfile+".jpg"
   spawn,"convert -density 300 -resize 1000 "+plotdir+plotfile+" "+plotdir+plotfile+".png"
   spawn,"epstopdf "+plotdir+plotfile

   stop
   
END

