PRO map_sites
  
  set_plot,'x'
  device,decomposed=0,pseudo_color=8,set_font='HELVETICA',/tt_font   
  !p.font=1
  loadct,13                     ; rainbow lut
  tvlct,R,G,B,/get
  n_table=!D.table_size
  R[0]=150
  G[0]=150
  B[0]=250
  R[n_table-1]=0                ; black
  G[n_table-1]=0
  B[n_table-1]=0
  R[n_table-2]=255              ; white
  G[n_table-2]=255
  B[n_table-2]=255
  R[n_table-3]=255              ; red
  G[n_table-3]=0
  B[n_table-3]=0
  R[n_table-4]=150              ; light red
  G[n_table-4]=150
  B[n_table-4]=100
  tvlct,R,G,B
  
  plotdir='../plots/plots/maps/'
  
  sitefile='../sites/sites_global256.dat'

  sites64 = indgen(64)
  sites16 = [4,6,9,13,17,21,25,29,33,37,43,47,49,55,57,61] - 1
  sites4 =  [9,17,29,55] - 1

                                ; open Site file
  nsites=256
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
     names[s]=siteshort
     lats[s]=sitelat
     lons[s]=sitelon
  ENDFOR
  close,1

  align[7]=1.0
  
  FOR i=0,nsites-1 DO BEGIN
     print,i+1,' ',names[i],lons[i],lats[i]
  ENDFOR
  
  lonmin=-180.0
  lonmax= 180.0
  latmin= -90.0
  latmax= 90.0
  map_limit=[latmin,lonmin,latmax,lonmax]
  
  set_plot,'ps'
  
  loadct,45,file='brewer.tbl'
  tvlct,R,G,B,/get
  R=shift(R,1)
  G=shift(G,1)
  B=shift(B,1)
  n_table=!D.table_size
  R[0]=150
  G[0]=175
  B[0]=255
  R[n_table-1]=0                ; black
  G[n_table-1]=0
  B[n_table-1]=0
  R[n_table-2]=150              ; white
  G[n_table-2]=220
  B[n_table-2]=150
  R[n_table-3]=200              ; red
  G[n_table-3]=0
  B[n_table-3]=0
  R[n_table-4]=0                ; light blue
  G[n_table-4]=0
  B[n_table-4]=200
  tvlct,R,G,B
  
  plotfile='global_sites_256.eps'
  DEVICE, /encapsul, bits_per_pixel=8, /COLOR,$
          Filename=plotdir+plotfile, $
          xsize=35, ysize=20,/Helvetica
  !p.font = 0
  map_set,0,0,/cylindrical,limit=map_limit,xmargin=[1,1],ymargin=[1,1],/noborder
;   map_grid,charsize=1.3,color=n_table-1,londel=10,latdel=10,/horizon,/fill_horizon
  map_continents,/fill_continents,color=n_table-2,/hires
  map_continents,/countries,colors=n_table-1,thick=2
                                ;map_grid,charsize=1.3,color=n_table-1,londel=10,latdel=10,label=2
                                ;xyouts,0.5,0.95,'Carboeurope Sites',/normal,color=n_table-1, $
                                ;   charsize=1.5,alignment=0.5

; Make a vector of 16 points, A[i] = 2pi/16: 
  A = FINDGEN(17) * (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
  USERSYM, COS(A), SIN(A),/fill

  FOR l=0,3 DO BEGIN
     ;; perform two loops to plot 4 region symbols on top of all others
     FOR i=0,nsites-1 DO BEGIN

        ;; guarantee that main sites are on top
        ok = 0

        idx = where(sites4 EQ i,count)
        IF count GT 0 THEN BEGIN

           IF l EQ 3 THEN ok=1

           psym = 6
           color = 1
                      
           symsize = 1.75
           symthick = 8.0

        ENDIF ELSE BEGIN
           
           symsize = 1.25
           symthick = 5.0

           idx = where(sites16 EQ i,count)
           IF count GT 0 THEN BEGIN
              IF l EQ 2 THEN ok=1
              psym = 5
              color = 2
           ENDIF ELSE BEGIN
              idx = where(sites64 EQ i,count)
              IF count GT 0 THEN BEGIN
                 IF l EQ 1 THEN ok=1
                 psym = 4
                 color = 3
              ENDIF ELSE BEGIN
                 IF l EQ 0 THEN ok=1
                 psym = 8
                 color = 4
              ENDELSE
           ENDELSE
           
        ENDELSE
        
        IF ok THEN plots,lons[i],lats[i],symsize=symsize,psym=psym,thick=symthick,color=color

     ENDFOR
     
  ENDFOR
  
  DEVICE,/close
  set_plot,'x'

;   spawn,"convert -density 300 -resize 1000 -quality 90 "+plotdir+plotfile+" "+plotdir+plotfile+".jpg"
  spawn,"convert -density 300 -resize 800 "+plotdir+plotfile+" "+plotdir+plotfile+".png"
  spawn,"epstopdf "+plotdir+plotfile

  stop
  
END

