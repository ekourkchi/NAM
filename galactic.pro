PRO GALACTIC,ra,dec,year,gl,gb,j,fk4=fk4,supergalactic=supergalactic
if (j eq 1) then begin
	glactc,ra/!dtor,dec/!dtor,year,gl,gb,j,/degree,fk4=fk4,$
		supergalactic=supergalactic
	gl *= !dtor
	gb *= !dtor
endif
if (j eq 2) then begin
	glactc,ra,dec,year,gl/!dtor,gb/!dtor,j,/degree,fk4=fk4,$
		supergalactic=supergalactic
	ra *= !dtor
	dec *= !dtor
endif
return
end

