function interpgrid,sgl,sgb,din,grid,cone=cone,full=full

h0 = .75
if ~keyword_set(full) then full = 0
if ~keyword_set(cone) then cone = 20
if ~keyword_set(grid) then begin
	tnstp = 3581577L
	grid = fltarr(tnstp,6)
	openr,/get_lun,runit,'grid_ascii.dat'
	for i=0L,tnstp-1 do begin
		readf,runit,format='(3(f5.1,1x),f7.3,1x,f7.3,1x,f7.3)',x,y,z,vx,vy,vz
		grid[i,*] = [x,y,z,vx,vy,vz]
	endfor
	close,runit
	free_lun,runit

endif
if n_elements(din) EQ 1 then din1 = din else din1=0
if full then din=findgen(38)+1
xhat= cos(sgl*!dtor)*cos(sgb*!dtor)
yhat= sin(sgl*!dtor)*cos(sgb*!dtor)
zhat= sin(sgb*!dtor)

x = xhat*din
y = yhat*din
z = zhat*din


sz = size(grid,/dim)
nstp = LONG(sz[0]^(1./3.))
; Get 8 points nearest
ndis = n_elements(din)
cz = fltarr(ndis)
for i = 0, ndis-1 do begin
	dis2 = (grid[*,0] - x[i])^2 + (grid[*,1] - y[i])^2 + (grid[*,2] - z[i])^2
	; Reduce number of distances to check
	nearby = where(dis2 lt 10.,count)
	dis3 = dis2[nearby]
	; Sort distances
	srt = sort(dis3)
	nearest = nearby[srt[0:7]]
	best = grid[nearest[0],*]
	cz[i] = (best[3]*xhat + best[4]*yhat + best[5]*zhat + h0*din[i])*100.
endfor
vsun = 240.
read_nests,'nests_VI38_sort.cat',cat_str,np,Vsun
READ_DAT,'V6',1381,cat_str,dat_str,h0
if ndis gt 1 then begin
	p=plot(din,cz,symbol='dot',linestyle='-',xtitle='D_m [Mpc]',ytitle='cz_m [km/s]',$
		title='SGL = '+strtrim(string(sgl),2)+', SGB = '+strtrim(string(sgb),2),$
		thick=2)
        ind = where(abs(cat_str.sgl -sgl) lt cone/cos(sgb*!dtor)$ 
		and abs(cat_str.sgb - sgb) lt cone  $
		and cat_str.dist_catalog gt 5,count) 
	if (count ge 1) then $
		p=plot(/over,dat_str[ind].d_calc,dat_str[ind].cz_calc,symbol='+',line='')
		p=plot(/over,findgen(40),h0*findgen(40)*100.,line='dash')
		if din1[0] gt 0 then p=plot(/over,[din1,din1],[-1000,4000],color='blue')
endif
return,cz
end

 

