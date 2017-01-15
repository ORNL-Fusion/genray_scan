function genray_create_rays, $
    xOffSet = _xOffSet, $
    zOffSet = _zOffSet, $
    angle_deg = _angle_deg, $
    width_m = width_m, $
    spread_deg_x = _spread_deg_x, $
    spread_deg_z = _spread_deg_z, $
    rayDensity = _rayDensity, $
    beam = _beam, $
    scan = _scan, $
    fwhm = _fwhm 

	if keyword_set(_beam) then beam = _beam else beam = 1	
	if keyword_set(_scan) then scan = _scan else scan = 0	
	if keyword_set(_xOffSet) then xOffSet = _xOffSet else xOffSet = 0.	
	if keyword_set(_zOffSet) then zOffSet = _zOffSet else zOffSet = 0.0	
	if keyword_set(_angle_deg) then angle_deg = _angle_deg else angle_deg = 0.0
	if keyword_set(_width_m) then width_m = _width_m else width_m = 0.05 	
	if keyword_set(_spread_deg_x) then spread_deg_x = _spread_deg_x else spread_deg_x = 0
    if keyword_set(_spread_deg_z) then spread_deg_z = _spread_deg_z else spread_deg_z = 0
	if keyword_set(_rayDensity) then rayDensity = _rayDensity else rayDensity = 3	
	if keyword_set(_fwhm) then fwhm = _fwhm else fwhm = width_m / 2.0	

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

    r_2d = sqrt( c1_2d^2 + c2_2d^2 )

    iiKeep = where( r_2d le width_m/2, iiKeepCnt )

    if iiKeepCnt gt 199 then begin

        print, 'ERROR: Too many rays'
        print, 'Adjust "nconea" first'
        stop

    endif

    xArr = c1_2d[iiKeep] * tan( angle_deg * !dtor ) + xOffSet 
    yArr = c2_2d[iiKeep]
    zArr = c1_2d[iiKeep] + zOffSet

    rArr = sqrt( xArr^2 + yArr^2 )
    pArr = atan( yArr, xArr ) * !radeg

    ; Angle in the x-y plane
    alp = 180 + $
        cos ( atan( c1_2d[iiKeep], c2_2d[iiKeep] ) + !pi ) * spread_deg_x $
        * sqrt ( c1_2d[iiKeep]^2 + c2_2d[iiKeep]^2 ) / ( width_m / 2 )

    ; Angle in the x-z plane
    bet = angle_deg + $
        sin ( atan( c1_2d[iiKeep], c2_2d[iiKeep] ) ) * spread_deg_z $ 
        * sqrt ( c1_2d[iiKeep]^2 + c2_2d[iiKeep]^2 ) / ( width_m / 2 )

    ; Gaussian HE11 beam profile
    c = fwhm / 2.35482 ; See https://en.wikipedia.org/wiki/Gaussian_function 
    
    pow = 100e3 * exp ( -r_2d[iiKeep]^2 / (2*c^2) )

    txt = ['ncone='+string(iiKeepCnt,format='(i3.3)') ]
    txt = [txt, ['powtot=', string(pow)] ]
    txt = [txt, ['zst='   , string(zArr)] ]
    txt = [txt, ['rst='   , string(rArr)] ] 
    txt = [txt, ['phist=' , string(pArr)] ] 
    txt = [txt, ['betast=', string(bet)] ] 
    txt = [txt, ['alfast=', string(alp)] ] 
    txt = [txt, ['alpha1=', string(pow*0)] ]
    txt = [txt, ['alpha2=', string(pow*0)] ]

endif

return, txt

end
