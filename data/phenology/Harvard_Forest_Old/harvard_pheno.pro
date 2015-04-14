PRO harvard_pheno

indir=''

year0 = 1990 ; start year
year1 = 2006 ; end year

ntrees = 6  ; maximun number of trees

species=['QURU','ACRU','BELE','TSCA','QUAL','QUVE']
species=['QURU','ACRU','QUAL']

season='spr'
phases=['BBRK','L75']

season='fall'
phases=['LCOLOR','LFALL']

nyears = year1-year0+1

threshold = 50 ; % threshold for diagnosing phase entry date

nspecies=n_elements(species)
nphases=n_elements(phases)
phasepos=intarr(nphases)
phaseval=intarr(nphases)
temp=''
date=intarr(3)
tree=0


openw,2,'Harvard_Forest.'+season+'.dat'


FOR year=year0,year1 DO BEGIN

   datematrix = fltarr(nspecies,ntrees,nphases)*!values.f_nan
   sp_mean=fltarr(nspecies)*!values.f_nan
   sp_dev=fltarr(nspecies)*!values.f_nan

   leapyear = ((((year mod 4) eq 0) and ((year mod 100) ne 0)) or ((year mod 400) eq 0))
   IF leapyear THEN BEGIN 
      ndays = 366 
      monthdays=[0,31,29,31,30,31,30,31,31,30,31,30,31]
   ENDIF ELSE BEGIN 
      ndays = 365
      monthdays=[0,31,28,31,30,31,30,31,31,30,31,30,31]
   ENDELSE

   syear=string(year,format='(I4)')

   IF file_search(indir+season+syear+'.txt') NE '' THEN BEGIN
      openr,1,indir+season+syear+'.txt'
      ;; read header
      readf,1,temp
      temp2=strtrim(strsplit(temp,',',/extract,/preserve_null),2)
      npos=n_elements(temp2)
      FOR p=0,nphases-1 DO phasepos[p]=where(temp2 EQ phases[p])
      
      treepos = where(temp2 EQ 'TREEID')
      
      WHILE (EOF(1) EQ 0) DO BEGIN
         readf,1,temp
         temp2=strtrim(strsplit(temp,',',/extract,/preserve_null),2)
         
         reads,strsplit(temp2[0],'/',/extract,/preserve_null),date
         doy = total(monthdays[0:date[0]])+date[1]
         
         treeid = strsplit(temp2[treepos],'-',/extract,/preserve_null)
         reads,treeid[1],tree
         s=(where(treeid[0] EQ species,scount))[0]
         
         IF scount GT 0 THEN BEGIN
            reads,temp2[phasepos],phaseval
            t = tree-1
            FOR p=0,nphases-1 DO BEGIN
               IF finite(datematrix[s,t,p],/nan) AND (phaseval[p] GE threshold) THEN BEGIN
                  datematrix[s,t,p]=doy
                                ;print,treeid[0],'-',treeid[1],s,t,p,doy
               ENDIF
            ENDFOR
         ENDIF
         
      END
      close,1
   ENDIf

   FOR s = 0, nspecies-1 DO BEGIN
      sp_mean[s] = mean(datematrix[s,*,*],/nan)
      sp_dev[s] = meanabsdev(datematrix[s,*,*],/nan)
      
      print,year,s,sp_mean[s],sp_dev[s]
   ENDFOR
   printf,2,year,mean(sp_mean,/nan),mean(sp_dev,/nan)


ENDFOR

close,2

stop
         
         


stop

END
