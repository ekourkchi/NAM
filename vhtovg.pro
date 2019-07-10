FUNCTION VhToVg,cat_str,u,v,w,Vsun
; Converts from Vheliocentric to VGalactic
; Adds projection of Vsunv in direction of galaxies to their heliocentric velocity
; Vsunv is input as u, v, w, Vsun

IF (N_PARAMS() EQ 0) THEN BEGIN
	PRINT,'Usage: VhToVg,cat_str,u,v,w,Vsun'
	RETURN,0
ENDIF

	; Convert degrees to radians
	sglr = cat_str.sgl/!radeg
	sgbr = cat_str.sgb/!radeg
	; Get unit vectors to all galaxies in SG coords
	xhat = [[cos(sglr)*cos(sgbr)],[sin(sglr)*cos(sgbr)],[sin(sgbr)]]
	Vsunv = [u, Vsun+v, w]
	; Convert from Galactic Rectilinear to Galactic Coordinates
	recToSphere,vsunv[0], Vsunv[1], Vsunv[2], Vsol, Vsun_gl, Vsun_gb
	; Transform velocity from galactic to equatorial
	galactic, Vsun_RA, Vsun_DE, 2000, Vsun_gl, Vsun_gb, 2
	; Transform from equatorial to supergalactic
	galactic, Vsun_RA, Vsun_DE, 2000, Vsun_sgl, Vsun_sgb, 1, /supergalactic
	; Transform from spherical to rectilinear
	sphereToRec, Vsol, Vsun_sgl, Vsun_sgb, Vsun_sgx, Vsun_sgy, Vsun_sgz
	vsun_sg= [vsun_sgx, vsun_sgy, vsun_sgz]
	return, cat_str.cz_helio + xhat#vsun_sg
END

