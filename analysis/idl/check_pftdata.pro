PRO check_pftdata

select = 0  ;; 0: 17 pft's, 1: 35 pft's
basedir = '/Users/stockli'
pftdir = 'pft'
pftregion = 'global'

missing = 129B

;; 17 PFT's
  ncvarnames17 = ['bar_all',$
                'enf_tem',$
                'enf_bor',$
                'dnf_bor',$
                'ebf_tro',$
                'ebf_tem',$
                'dbf_tro',$
                'dbf_tem',$
                'dbf_bor',$
                'ebs_all',$
                'dbs_tem',$
                'dbs_bor',$
                'c3g_arc',$
                'c3g_nar',$
                'c4g_all',$
                'cro_all',$
                'wat_all']

;; 35 PFT's
  ncvarnames35 = ['bar_all',$
                'enf_tem',$
                'enf_bor',$
                'dnf_bor',$
                'ebf_tro',$
                'ebf_tem',$
                'dbf_tro',$
                'dbf_tem',$
                'dbf_bor',$
                'ebs_all',$
                'dbs_tem',$
                'dbs_bor',$
                'c3g_arc',$
                'c3g_nar',$
                'c4g_all',$
                'cro_brl',$
                'cro_cas',$
                'cro_cot',$
                'cro_grn',$
                'cro_mze',$
                'cro_mil',$
                'cro_oil',$
                'cro_oth',$
                'cro_pot',$
                'cro_pul',$
                'cro_rap',$
                'cro_ric',$
                'cro_rye',$
                'cro_sor',$
                'cro_soy',$
                'cro_sug',$
                'cro_sgb',$
                'cro_sgc',$
                'cro_wht',$
                'wat_all']


IF select THEN ncvarnames = ncvarnames35 ELSE ncvarnames = ncvarnames17

nvar = n_elements(ncvarnames)

infiles =  basedir + '/'+pftdir+'/'+ncvarnames+'.'+pftregion+'.nc'

FOR v=0,nvar-1 DO BEGIN

   print,'Reading: ',ncvarnames[v]

   ncid = ncdf_open(infiles[v])
   ncdf_varget,ncid, ncvarnames[v], data      
   ncdf_close,ncid

   a = where(data EQ missing,count)

   IF count GT 0 THEN print,'missing data points: ',count


ENDFOR




END
