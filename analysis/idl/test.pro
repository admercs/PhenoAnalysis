PRO test

  device,decomposed=0
  loadct,13

  nx = 360L
  ny = 180L

  dx = 360./float(nx)
  dy = 180./float(ny)

  lon = (findgen(nx)/float(nx) * 360.0 - 180.0 + 0.5 * dx) # replicate(1.0,ny)
  lat = replicate(1.0,nx) # (findgen(ny)/float(ny) * 180.0 - 90.0 + 0.5 * dy)

  indir='/Users/stockli/pft/'

  varname = 'c4g_all'

  infile1='pft.global_0.1.nc'
  infile2=varname+'.Global-FULL.nc'

  ncid = ncdf_open(indir+infile1)
  ncdf_varget,ncid,varname,data1
  ncdf_close,ncid

  ncid = ncdf_open(indir+infile1)
  ncdf_varget,ncid,'wat_all',data3
  ncdf_close,ncid

  ncid = ncdf_open(indir+infile2)
  ncdf_varget,ncid,varname,data2
  ncdf_close,ncid


  data1 = rebin(data1,nx,ny)
  data2 = rebin(data2,nx,ny)
  data3 = rebin(data3,nx,ny)

  wgt = cos(lat*!dtor)
  print,"C4G (OK): ",mean(data1*wgt)/mean(wgt)
  print,"C4G (NEW): ",mean(data2*wgt)/mean(wgt)
  print,"Water: ",mean(data3*wgt)/mean(wgt)



  window,0,xsize=nx,ysize=ny
  tv,bytscl(data1,min=0,max=100)

  window,1,xsize=nx,ysize=ny
  tv,bytscl(data2,min=0,max=100)


  stop

END
