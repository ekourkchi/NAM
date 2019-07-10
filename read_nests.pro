PRO read_nests,catfile,cat_str,np,Vsun

; Read galaxy catalogs from Brent Tully that begin with 'nests'

IF (N_PARAMS() EQ 0) THEN BEGIN
	PRINT,'read_nests,catfile,cat_str,np,Vsun'
	RETURN
ENDIF

IF (Vsun lt 100) then begin
	PRINT, 'read_nests: Vsun is too lo'
	RETURN
ENDIF
	
o_err_max = 0.10d0
dmodlimit = 32.5999d0 ; =33Mpc, distance where lumlimit2 starts
dmodlimit2 = 29.2201d0 ; =6.98 Mpc, distance mod where lumlimit starts
lumlimit = 9.5d0  ; Luminosity Limit
lumlimit2 = 10.3d0

IF ~keyword_set(Vsun) then Vsun = 239.  
u= 9.8 & v = 11.6 & w = 5.9 


form='(i6,1x,f8.4,1x,f8.4,1x,i3,1x,f6.2,1x,f5.2,1x,f4.2,1x,i3,'
form = form + '1x,f5.0,1x,f5.0,1x,f5.0,1x,i4,1x,f5.2,1x,f6.3,1x,'
form = form + 'f6.3,2(1x,f8.4),8(1x,f6.0),1x,a22,1x,i7)'

  K_lum_cor = double(1.)
;  catdir = '~/Documents/nam/Catalogs/'   
  CLOSE,/all
  OPENR,/get_lun,rdunit,catfile
  name1 = ''
  pgc = 1L
  mtype = 0
  ng = 0
  sgl = 1d0
  sgb = 1d0
  gl = 1d0
  gb = 1d0
  nbg_id = ''

  ; count number of mass tracers
  cat_str = 	{ pgc : 0L,$
   		name : '',$
		pgc1 : 0L,$
		sgl : 0d0,$
  		sgb : 0d0,$
  		nd : 0L,$
		x : 0d0, $
		y : 0d0, $
		z : 0d0, $
		dist_catalog : 0d0,$
		obs_err : 0d0,$
  		cz_helio : 0d0,$
		cz_err : 0d0,$
  		cz_catalog : 0d0,$
  		sigma_zi : 0d0,$
  		luminosity : 0d0,$
  		mtype : 0 $
  		}
		

  READF,rdunit,name1
  np = 0
  WHILE (~EOF(rdunit)) DO BEGIN
	   READF,rdunit,name1
	    np++
  ENDWHILE
  CLOSE, rdunit

  cat_str = REPLICATE(cat_str,np+1)

  OPENR,rdunit,catfile
  OPENR,/get_lun,disunit,'distances.dat'
  READF,rdunit,name1
  np2=0 
  FOR nn = 1, np DO BEGIN
 	READF,rdunit,format=form, $
 	 pgc,sgl,sgb,Nd,dMpc,dmod,d_e,Nv,Vls1,v_rms,cz_e,NK,Kmag,logLK,$
	 logMass,gl,gb,vh,vgsre,vgsrr,vgsrm,vlgkm,vlgcv,vls,vcmb,name1,pgc1

 	IF (d_e eq 0.0) then d_e = !VALUES.D_NAN
 	IF (dMpc eq 0.0) then dMpc = !VALUES.D_NAN
 	IF (v_rms eq 0.0) then rms = !VALUES.D_NAN
        IF (pgc GT 900000) then logLK = logLK + ALOG10(1.5)

	; Decide which objects in Nests to skip
	; If in 80000 series it is a true test particle, so do not skip
	IF (pgc LT 800000 OR PGC GT 800100) THEN BEGIN
		; If beyond 33 Mpc an lumlimit2 then skip
	IF (dmod GT dmodlimit AND logLK LT lumlimit2) THEN CONTINUE
		; If beyond 7 Mpc an lumlimit then skip if no 'good' distance measure
	IF (dmod GT dmodlimit2 AND logLK LT lumlimit)  THEN $ 
	       IF (d_e gt o_err_max OR finite(d_e,/NAN)) THEN CONTINUE
        ENDIF

   	np2++
   	cat_str[np2].pgc = pgc
	if (strtrim(name1,2) eq 'Galaxy') then name1 ='MW'
;;;	if (strmid(name1,0,3) eq 'PGC' ) then name1='P'+strmid(name1,3)
;;;;	if (strmid(name1,0,3) eq 'ESO' ) then name1='E'+strmid(name1,3)
;;	if (strmid(name1,0,3) eq 'NGC' ) then name1='N'+strmid(name1,3)
;;	if (strmid(name1,0,3) eq 'UGC' ) then name1='U'+strmid(name1,3)
;;	if (strmid(name1,0,2) eq 'IC' ) then name1='I'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'E0' ) then name1='E'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'E0' ) then name1='E'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'U0' ) then name1='U'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'U0' ) then name1='U'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'N0' ) then name1='N'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'N0' ) then name1='N'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'I0' ) then name1='I'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'I0' ) then name1='I'+strmid(name1,2)
;;	if (strmid(name1,0,2) eq 'C0' ) then name1='C'+strmid(name1,2)
;;	pos = strpos(name1,'-')
;;	if (pos ne -1) then name1 = strmid(name1,0,pos)+strmid(name1,pos+1)
	cat_str[np2].name = strtrim(name1,2)
        cat_str[np2].sgl = sgl
        cat_str[np2].sgb = sgb
        cat_str[np2].nd = nd
        cat_str[np2].obs_err = d_e
        cat_str[np2].cz_helio = Vh
        cat_str[np2].cz_err = cz_e
;       cat_str[np2].sigma_zi = 10.
	cat_str[np2].luminosity=10.^logLK
	cat_str[np2].mtype = mtype
	cat_str[np2].pgc1 = pgc1
        cat_str[np2].cz_catalog = VhToVg(cat_str[np2],u,v,w,Vsun)
        READF, disunit, icnt3, dist_adjust
	if icnt3 ne pgc then begin
		PRINT,'read_nests: distances.dat catalog not in sync'
		stop
	endif
	if (finite(d_e)) then $
        	cat_str[np2].dist_catalog =dMpc $
	else $
        	cat_str[np2].dist_catalog =dist_adjust
	if ~finite(dMpc) then $
 		cat_str[np2].dist_catalog = cat_str[np2].cz_catalog/75.
;	dMpc2=10d0^((dmod-25d0)/5.d0) 
	sphereToRec,cat_str[np2].dist_catalog,sgl*!dtor,sgb*!dtor,x,y,z
	cat_str[np2].x = x
	cat_str[np2].y = y
        cat_str[np2].z = z
 ENDFOR
 np = np2
 cat_str = cat_str[0:np]
 CLOSE,rdunit
 FREE_LUN, rdunit
 CLOSE,disunit
 RETURN
END
