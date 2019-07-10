pro spheretorec, r, lon, lat, x, y, z
; convert spherical r,lon,lat positions to Rectilinear x,y,z
; lon - longitude in radians
; lat - latitude in radian

; Right Handed Cartesian system
; +x towards lat=0, lon=0
; +y toward lat=0, lon= pi/2 
; +z toward lat=pi/2

dxy = r*cos(lat)
x = dxy*cos(lon)
y = dxy*sin(lon)
z = r*sin(lat)
return
end
