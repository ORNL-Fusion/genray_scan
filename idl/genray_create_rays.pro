function genray_create_rays, $
    xOffSet = _xOffSet, $
    zOffSet = _zOffSet, $
    angle_deg = _angle_deg, $
    width_m = width_m, $
    spread_deg = _spread_deg, $
    rayDensity = _rayDensity, $
    beam = _beam, $
    scan = _scan 

	if keyword_set(_beam) then beam = _beam else beam = 1	
	if keyword_set(_scan) then scan = _scan else scan = 0	
	if keyword_set(_xOffSet) then xOffSet = _xOffSet else xOffSet = 0.07	
	if keyword_set(_zOffSet) then zOffSet = _zOffSet else zOffSet = 0.005	
	if keyword_set(_angle_deg) then angle_deg = _angle_deg else angle_deg = 20	
	if keyword_set(_width_m) then width_m = _width_m else width_m = 0.05 	
	if keyword_set(_spread_deg) then spread_deg = _spread_deg else spread_deg = 10.0	
	if keyword_set(_rayDensity) then rayDensity = _rayDensity else rayDensity = 3	

if scan then begin

    ; Scan ray file

    nRayZ = 10  
    nRayA = 20 

    zMinRay = 2.55
    zMaxRay = 2.65

    aMinRay = -35
    aMaxRay = +35

    rayR = 7d-2
 
    zRay = fIndGen(nRayZ)*(zMaxRay-zMinRay)/(nRayZ-1)+zMinRay - g_z0
    aRay = fIndGen(nRayA)*(aMaxRay-aMinRay)/(nRayA-1)+aMinRay

    zRay2D = rebin(zRay,nRayZ,nRayA)
    aRay2D = transpose(rebin(aRay,nRayA,nRayZ))

    openw, lun, "rays.txt", /get_lun
    
    printf, lun, 'ncone=', nRayZ*nRayA
    printf, lun, 'powtot=', zRay2D[*]*0+0.02d6 
    printf, lun, 'zst=', zRay2D[*]
    printf, lun, 'rst=', zRay2D[*]*0+rayR
    printf, lun, 'phist=', zRay2D[*]*0
    printf, lun, 'betast=', aRay2D[*]
    printf, lun, 'alfast=', zRay2D[*]*0+180.0
    printf, lun, 'alpha1=', zRay2D[*]*0
    printf, lun, 'alpha2=', zRay2D[*]*0

    close, lun

endif


if beam eq 1 then begin

    ; Square beam
 
    c1 = fIndGen(rayDensity)/(rayDensity-1)*width_m-width_m/2
    c2 = fIndGen(rayDensity)/(rayDensity-1)*width_m-width_m/2

    c1_2d = rebin(c1,rayDensity,rayDensity)
    c2_2d = transpose(rebin(c2,rayDensity,rayDensity))

    ; Extract circular beam 

    iiKeep = where( sqrt( c1_2d^2 + c2_2d^2 ) le width_m/2, iiKeepCnt )

    xArr = c1_2d[iiKeep] * tan( angle_deg * !dtor ) + xOffSet 
    yArr = c2_2d[iiKeep]
    zArr = c1_2d[iiKeep] + zOffSet

    rArr = sqrt( xArr^2 + yArr^2 )
    pArr = atan( yArr, xArr ) * !radeg

    ; This seems to be the angle away from the vector from the origin to the point in the x-y plane
    ; i.e., for parallel rays in the x-y plane we need to subtract an angle 
    
    alp = sin ( atan( c1_2d[iiKeep], c2_2d[iiKeep] ) ) * spread_deg + 0 - atan(c1_2d[iiKeep],xArr)*!radeg 
    alp = alp*0 + 180

    ; Angle in the x-z plane
    bet = cos ( atan( c1_2d[iiKeep], c2_2d[iiKeep] ) ) * spread_deg + angle_deg
    bet = bet * 0

    pow = xArr*0 + 0.02d6

    txt = ['ncone='+string(iiKeepCnt,format='(i3.3)') ]
    txt = [txt, 'powtot='+ StrJoin(string(pow)) ]
    txt = [txt, 'zst='+    StrJoin(string(zArr)) ]
    txt = [txt, 'rst='+    StrJoin(string(rArr)) ] 
    txt = [txt, 'phist='+  StrJoin(string(pArr)) ] 
    txt = [txt, 'betast='+ StrJoin(string(bet)) ] 
    txt = [txt, 'alfast='+ StrJoin(string(alp)) ] 
    txt = [txt, 'alpha1='+ StrJoin(string(pow*0)) ]
    txt = [txt, 'alpha2='+ StrJoin(string(pow*0)) ]

endif

return, txt

end