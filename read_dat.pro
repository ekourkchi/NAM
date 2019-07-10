PRO READ_DAT,suffix_idtag,np,cat_str,dat_str,h0

; Read .dat output file from NAM into structure

IF (N_PARAMS() EQ 0) THEN BEGIN
	PRINT, 'Usage: READ_DAT,suffix_idtag,np,cat_str,dat_str,h0'
	RETURN
ENDIF
datfile = 'test'+suffix_idtag+'.dat'
;FORM160 = "( I3, A9, F11.5, 2F9.2, 2F7.2, E11.3)"
FORM160 = "( I3, I9, F11.3, F9.2, F10.4, F6.2, F7.2, E11.3)"
namei = ''
j=1 & trl=1
dat_str = {icnt:0L,name:'',mass_cat:0d0,mass_calc:0d0,cz_cat:0d0,cz_calc:0d0,$
	d_cat:0d0,d_calc:0d0,o_err:0d0,cz_err:0d0,chi:0d0,sos:0d0,bad:0d0,$
	Rcut:0d0,nsos:0d0,x:0d0,y:0d0,z:0d0,sgl:0d0,sgb:0d0}
dat_str = replicate(dat_str,np+1)
dat_str.icnt = lindgen(np+1)
header = ''
openr,/get_lun,rdunit,datfile
readf,rdunit,header
line = ''
for i= 1,np do begin
	readf,rdunit,line
;	print,line
	reads,form=FORM160,line,$
	       	trl,icnti,mass_cati,cz_cati,cz_calci,d_cati,o_erri,sosi 
	;namei = strtrim(namei,2)
	dat_str[i].icnt = icnti
	dat_str[i].mass_cat = mass_cati
	dat_str[i].mass_calc = mass_cati
	dat_str[i].cz_calc = cz_calci
	dat_str[i].d_calc = d_cati
	dat_str[i].o_err = o_erri
	dat_str[i].sos = sosi
	delcz = cz_calci - cz_cati
	dat_str[i].chi = SQRT(delcz^2/((o_erri*d_cati*h0*100.)^2+cat_str[i].cz_err^2 + 20d0^2))
endfor
dat_str.cz_cat = cat_str.cz_catalog
dat_str.name = cat_str.name
dat_str.x = cat_str.x
dat_str.y = cat_str.y
dat_str.z = cat_str.z
dat_str.sgl = cat_str.sgl
dat_str.sgb = cat_str.sgb
dat_str.cz_err = cat_str.cz_err
dat_str.d_cat = cat_str.dist_catalog

free_lun,rdunit
return
end
