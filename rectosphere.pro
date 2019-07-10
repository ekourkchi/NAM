pro rectosphere, x, y, z, r, lon, lat
; convert Rectilinear x,y,z to spherical coordinates

; Inputs
; x,y,z cartesian positions in arbitrary units

; Outputs
; r - radius in same arbitrary units
; lat - latitude in radians
; lon - longitude in radians

; Cartesian system
; +x towards lon=0, lat=0
; +y toward lon=pi/2, lat=0
; +z toward lon=0, lat=pi/2

dxy2 = (x^2 + y^2)
dxy = sqrt(dxy2)
lat = atan(z/dxy)
lon = atan(y,x)
r = sqrt(dxy2 + z^2)
return
end
