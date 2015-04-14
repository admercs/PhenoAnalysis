PRO hdfstructure,sd_id

  HDF_SD_FILEINFO, sd_id, datasets, attributes
  FOR I=0,datasets-1 DO begin   ;start cycle SDS's
     sds_id=HDF_SD_SELECT(sd_id,I)
     HDF_SD_GETINFO,sds_id,DIMS=dims,NDIMS=ndims,NAME=name,NATTS=natts
     PRINT,FORMAT='(I0,".",3A0,4(I0,:,"x"))', I, "Short name:  ", name , ", size: " , dims
     FOR J=0,natts-1	DO BEGIN ;start cycle attributes for every SDS
        HDF_SD_ATTRINFO,sds_id,J,NAME=name,DATA=attr_dat
        PRINT,FORMAT='(A0,":  ",5(A0))',name,attr_dat
        attr_dat=""
        name=""
     ENDFOR                     ;end cycle attributes
  ENDFOR

END

PRO read_avhrrltdr, basedir, year, sitename, sitelon, sitelat, obsvarname, nodata, obsvar, obserror, obsmask, obsdays, obsunits, obstitle, has_obs

; parameters
dx=1 ; number of x pixels +/- (0..x)
dy=1 ; number of y pixels +/- (0..y)

; vars
temp=''
avhrryear=0
ahvrrday=0
has_obs=0

syear=''
sday=''

; MSB -> LSB: 15 -> 0
bits=[0,0,0,9,9,9,9,9,0,9,0,9,9,9,9,9]
;bits=[0,0,0,9,9,9,9,9,0,9,0,0,9,0,0,0]
nbits=16

np=(2*dx+1)*(dy*2+1)


CASE obsvarname OF
    'LAI': BEGIN
        indir=basedir+'ndvi/ltdr_ndvi/'
        ltdr_prefix='AVH15C2'
        infile_prefix=ltdr_prefix+'.A'
        dt=8
        nt=46
        obsunits='m2/m2'
        obstitle='Leaf Area Index'
        hdfname='LAI'
        QAname='QA'
    END
    'FPAR': BEGIN
        indir=basedir+'ndvi/ltdr_ndvi/'
        ltdr_prefix='AVH15C2'
        infile_prefix=ltdr_prefix+'.A'
        dt=8
        nt=46
        obsunits='-'
        obstitle='Fraction of PAR absorbed by plants'
        hdfname='FPAR'
        QAname='QA'
    END
    'NDVI': BEGIN
        indir=basedir+'ndvi/ltdr_ndvi/'
        ltdr_prefix='AVH13C1'
        infile_prefix=ltdr_prefix+'.A'
        dt=1
        nt=365
        obsunits='-'
        obstitle='Normalized Difference Vegetation Index'
        hdfname='NDVI'
        QAname='QA'
    END
ENDCASE

obsvar=fltarr(np,nt)
obsvar[*]=nodata
obserror=fltarr(np,nt)
obsmask=fltarr(np,nt)

biomemask=bytarr(np,nt)
obsdays=intarr(nt)
obsvartemp=fltarr(np)

reads,year,syear
syear=strtrim(syear,2)

FOR avhrrday=1,365 DO BEGIN

    obsdays[(avhrrday-1)/dt]=avhrrday-1

    reads,avhrrday,sday
    sday=strtrim(sday,2)
    IF (avhrrday LT 10) THEN sday='00'+sday
    IF (avhrrday GE 10) AND (avhrrday LT 100) THEN sday='0'+sday

    infile=(file_search(indir+infile_prefix+syear+sday+'*.hdf'))[0]

    IF infile NE '' THEN BEGIN
; open HDF file

        ;fid=EOS_GD_OPEN(infile,/READ)
        ;gid=EOS_GD_ATTACH(fid,'Grid')

        ;result = EOS_GD_READATTR(gid,'_FillValue',hdffill)

        FileHandle=HDF_OPEN(infile,/READ)
        sd_id=HDF_SD_START(infile,/READ)

        ;hdfstructure,sd_id

        sds_no = HDF_SD_NAMETOINDEX(sd_id, hdfname)
        sds_id = HDF_SD_SELECT(SD_ID,sds_no)
        HDF_SD_GETINFO, sds_id, DIMS=dims
        hdffill=-9999.
        IF HDF_SD_ATTRFIND(sds_id,"_FillValue") GE 0 THEN $
          HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"_FillValue"),DATA=hdffill
        
        hdfrange=[-1000,10000]
        ; why does it not find the valid_range attribute?
                                ;IF HDF_SD_ATTRFIND(sds_id,"valid_range") GE 0 THEN $
                                ;  HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"valid_range"),DATA=hdfrange

        scale=1.
        IF HDF_SD_ATTRFIND(sds_id,"scale_factor") GE 0 THEN  $
          HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"scale_factor"),DATA=scale
        offset=0.
        IF HDF_SD_ATTRFIND(sds_id,"add_offset") GE 0 THEN  $
          HDF_SD_ATTRINFO,sds_id,HDF_SD_ATTRFIND(sds_id,"add_offset"),DATA=offset
        
        scale=float(scale[0])
        offset=float(offset[0])
        hdffill=hdffill[0]
        nx=dims[0]
        ny=dims[1]

        x0 = round((180.+sitelon)/360.*float(nx))
        y0 = round((90.-sitelat)/180.*float(ny))


        ; get data
        HDF_SD_GETDATA, sds_id, data,start=[x0-dx,y0-dy],count=[2*dx+1,2*dy+1]

        ;HDF_SD_GETDATA, sds_id, array
        ;write_png,'image.png',reverse(bytscl(array,min=0,max=6000),2)

        ; get QA bits
        sds_no = HDF_SD_NAMETOINDEX(sd_id, QAname)
        sds_id = HDF_SD_SELECT(SD_ID,sds_no)
        HDF_SD_GETDATA, sds_id, qa,start=[x0-dx,y0-dy],count=[2*dx+1,2*dy+1]

        ;HDF_SD_GETDATA, sds_id, array
        ;array=((array AND 1) OR (ishft(array,-1) AND 1))*255B
        ;array=((ishft(array,-6) AND 1))*255B
        ;write_png,'qa.png',reverse(byte(array),2)

;        stop

; close HDF file        
        HDF_SD_END,sd_id
        HDF_CLOSE,FileHandle
        
        baddata=where((data LT hdfrange[0]) OR (data GT hdfrange[1]) OR (data EQ hdffill),badcount)
        obsvartemp=scale*(float(data) - offset)
        IF badcount GT 0 THEN obsvartemp[baddata] = nodata

        ; estimate 100% error range
        maxerror=float(hdfrange(1)-hdfrange(0))*scale ; 100% error is full variable range
        
        ; make variables one-dimensional
        obsvartemp=(obsvartemp)[*]
        qa=(qa)[*]

        FOR p=0,np-1 DO BEGIN
            qflag=1
            verybad=0
            error=0.01*maxerror ; assume 1% error for best quality
                                ; and 10% error for every quality flag
                                ; which is weighted as above 
            FOR b=0,nbits-1 DO BEGIN
                IF (b EQ 3) AND (ishft(qa[p],-b) AND 1) THEN print,'Pixel is over Water'
                IF (ishft(qa[p],-b) AND 1) THEN error=error+float(bits[nbits-1-b])*0.1*maxerror
                IF (ishft(qa[p],-b) AND 1) AND (bits[nbits-1-b] EQ 9) THEN verybad=1
            ENDFOR
            
            obserror[p,(avhrrday-1)/dt]=error
            obsmask[p,(avhrrday-1)/dt]=1-verybad
            obsvar[p,(avhrrday-1)/dt]= obsvartemp[p]
            
        ENDFOR
    ENDIF


ENDFOR

gooddata = where(obsvar NE nodata,goodcount)
IF goodcount GT 0 THEN has_obs=1

END
