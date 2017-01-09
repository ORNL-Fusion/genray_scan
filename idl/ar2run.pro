
    g_z0 = 2.6

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

    ; Beam ray file

    rOffSet = 0.07
    zOffSet = 0.005
    angle = 0 
    width = 0.062/4
    spread = 10.0
    nPts = 3 
    nAngles = 3
   
    dY = width/2*cos((90-angle)*!dtor)
    y1 = rOffset-dy
    y2 = rOffset+dy

    dX = width/2*sin((90-angle)*!dtor)
    x1 = zOffset-dX
    x2 = zOffset+dX

    xArr = fIndGen(nPts)/(nPts-1)*(x2-x1)+x1
    yArr = fIndGen(nPts)/(nPts-1)*(y2-y1)+y1

    rRay2D = fltArr(nPts,nAngles) 
    zRay2D = fltArr(nPts,nAngles) 
    aRay2D = fltArr(nPts,nAngles) 

    da = spread*2/(nAngles-1)
    for a=0,nAngles-1 do begin
        rRay2D[*,a] = yArr
        zRay2D[*,a] = xArr
        aRay2D[*,a] = a*dA-spread+angle
    endfor 

    openw, lun, "beam.txt", /get_lun
    
    printf, lun, 'ncone=', n_elements(rRay2D[*])
    printf, lun, 'powtot=', zRay2D[*]*0+0.02d6 
    printf, lun, 'zst=', zRay2D[*]
    printf, lun, 'rst=', rRay2D[*]
    printf, lun, 'phist=', zRay2D[*]*0
    printf, lun, 'betast=', aRay2D[*]
    printf, lun, 'alfast=', zRay2D[*]*0+180.0
    printf, lun, 'alpha1=', zRay2D[*]*0
    printf, lun, 'alpha2=', zRay2D[*]*0


    close, lun


	freq = 28d9
    wrf = 2*!pi*freq

    dispersionAxialWaveLength = 0.00
    r0 = 100
    radius = 0.06

	nphi = 0;2*!pi*r0/dispersionAxialWaveLength

    bField_flat = 0
	bField_eqdsk = 0
    bField_linear = 0
    bField_numeric = 1

    bFieldFileName = 'input/bFieldFile.txt'
    nHeader = 9
    nLines = file_lines(bFieldFileName)

    bFieldRead = replicate({r:0.0, z:0.0, br:0.0, bt:0.0, bz:0.0, bmag:0.0, psi:0.0},nLines-nHeader)
    openr, lun, bFieldFileName, /get_lun
    skip_lun, lun, nHeader, /lines     
    readf, lun, bFieldRead
    free_lun, lun

    bField_r1 = bFieldRead.z+r0
    bField_z1 = bFieldRead.r

    bField_r2 = +bField_r1 
    bField_z2 = -bField_z1 

    bField_br = -bFieldRead.bz
    bField_bt = +bFieldRead.bt
    bField_bz = +bFieldRead.br

    bField_psi = bFieldRead.psi

    bField_r = [bField_r1,bField_r2]
    bField_z = [bField_z1,bField_z2]

    bField_psi = [bField_psi,bField_psi]
    bField_br = -[bField_br,bField_br]
    bField_bt = [bField_bt,-bField_bt]
    bField_bz = [bField_bz,-bField_bz]

	gaussian_profiles = 0
	flux_profiles = 0
   	numeric_flux_profiles = 0
	fred_namelist_input = 0
    mpex_profiles = 1

    x0 = r0
    y0 = 0
    rCenter = r0+2.6
    zCenter = y0

	atomicZ	= [-1,1]
	amu = [_me_amu,2]

	nn = [	[0.0,		0.0,		9.0,	10.0],$ ; electrons are spec 0
			[0.0,	    5d19,		9.0,	10.0] ]

	tt = [	[1e3,	1e3,		6.0,	4.0],$
			[1e3,	1e3,		6.0,	4.0] ]

	; Grid
	nR = 257 
	nZ = 128 
	rMin = r0+2.3
	rMax = r0+2.9
	zMin = -radius*2
	zMax = +radius*2

    rDomainBox = [rMin,rMax,rMax,rMin,rMin]
    zDomainBox = [zMin,zMin,zMax,zMax,zMin]

    rLim = rDomainBox
    zLim = zDomainBox

    angle = fIndGen(360)/359*2*!pi
    rlcfs = (0.1*cos(angle))+rCenter
    zlcfs = (0.06*sin(angle)) 

    triangulate, bField_r, bField_z, tr, b

    ; These should be moved inside of ar2_create_input

    rTmp=fIndGen(nR)/(nR-1)*(rMax-rMin)+rMin
    zTmp=fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

    rMPEX = rTmp
    zMPEX = zTmp

    psi2d = trigrid(bField_r, bField_z, bField_psi, tr, nx=nR, ny=nZ, xOut=rTmp, yOut=zTmp)
    psi2d = psi2d/max(psi2d)
    psi2d = 1-psi2d

    br2d_numeric = trigrid(bField_r, bField_z, bField_br, tr, nx=nR, ny=nZ, xOut=rTmp, yOut=zTmp)
    bt2d_numeric = trigrid(bField_r, bField_z, bField_bt, tr, nx=nR, ny=nZ, xOut=rTmp, yOut=zTmp)
    bz2d_numeric = trigrid(bField_r, bField_z, bField_bz, tr, nx=nR, ny=nZ, xOut=rTmp, yOut=zTmp)

    bMag2d_numeric = sqrt(br2d_numeric^2+bt2d_numeric^2+bz2d_numeric^2)
    bMagMPEX = bMag2d_numeric

    neSig = 0.05
    neMax = 6e19
    neMin = 1e17
    neMPEX = (neMax * exp(-(1-psi2d)^1.5/neSig^2))>neMin

    teSig = 0.05 
    teMax = 4 
    teMin = 4
    teMPEX = (teMax * exp(-(1-psi2d)^1.5/teSig^2))>teMin

    ; Read in GENRAY data


    genrayFileName = 'genray.nc'

	cdfId = ncdf_open ( genrayFileName, /noWrite ) 
        nCdf_varGet, cdfId, 'densprofxz', g_densprofxz 
        nCdf_varGet, cdfId, 'temprof', g_temprof 
        nCdf_varGet, cdfId, 'wx', g_wx 
        nCdf_varGet, cdfId, 'wy', g_wy 
        nCdf_varGet, cdfId, 'wz', g_wz 
        nCdf_varGet, cdfId, 'nrayelt', g_nrayelt 
        nCdf_varGet, cdfId, 'w_x_densprof_nc', g_w_x_densprof_nc 
        nCdf_varGet, cdfId, 'eqdsk_x', g_eqdsk_x 
        nCdf_varGet, cdfId, 'eqdsk_z', g_eqdsk_z 
        nCdf_varGet, cdfId, 'eqdsk_psi', g_eqdsk_psi
        nCdf_varGet, cdfId, 'eqdsk_beq', g_eqdsk_beq
        nCdf_varGet, cdfId, 'bmodprofxz', g_bmodprofxz
        nCdf_varGet, cdfId, 'wnper', g_wnper
        nCdf_varGet, cdfId, 'wnpar', g_wnpar
        nCdf_varGet, cdfId, 'vgr_r', g_vgr_r
        nCdf_varGet, cdfId, 'wn_r', g_wn_r
        nCdf_varGet, cdfId, 'cwexde', g_cwexde
        nCdf_varGet, cdfId, 'cweyde', g_cweyde
        nCdf_varGet, cdfId, 'cwezde', g_cwezde
        nCdf_varGet, cdfId, 'sb_x', g_sb_x
        nCdf_varGet, cdfId, 'sb_y', g_sb_y
        nCdf_varGet, cdfId, 'sb_z', g_sb_z
        nCdf_varGet, cdfId, 'delpwr', g_delpwr
        nCdf_varGet, cdfId, 'ste', g_ste
        nCdf_varGet, cdfId, 'i_ox_conversion', g_i_ox
        nCdf_varGet, cdfId, 'Rgrid', g_pwr_r
        nCdf_varGet, cdfId, 'Zgrid', g_pwr_z
        nCdf_varGet, cdfId, 'spwr_rz_e', g_pwr_e
        nCdf_varGet, cdfId, 'spwr_rz_i', g_pwr_i
        nCdf_varGet, cdfId, 'spwr_rz_cl', g_pwr_cl
	ncdf_close, cdfId

    gEx = complex(g_cwexde[*,*,0],g_cwexde[*,*,1])
    gEy = complex(g_cweyde[*,*,0],g_cweyde[*,*,1])
    gEz = complex(g_cwezde[*,*,0],g_cwezde[*,*,1])
    gEVecMag = sqrt(gEx^2+gEy^2+gEz^2)

    gBVecMag = sqrt(g_sb_x^2+g_sb_y^2+g_sb_z^2)
    gBx = g_sb_x/gBVecMag
    gBy = g_sb_y/gBVecMag
    gBz = g_sb_z/gBVecMag

    gEDotB = gEx*gBx+gEy*gBy+gEz*gBz
    gAngle = acos(gEDotB/(gBVecMag*gEVecMag))*!radeg

    rGENRAY = g_eqdsk_z+g_z0+r0
    zGENRAY = g_eqdsk_x

    g_ray_r = g_wz*1e-2+g_z0
    g_ray_z = g_wx*1e-2

    g_wkPer = g_wnper*wrf/_c
    g_wkPar = g_wnpar*wrf/_c
    g_wkMag = sqrt(g_wkPer^2+g_wkPar^2)
    g_kr = g_wn_r*wrf/_c
    g_vpr = wrf/g_kr

    neGENRAY = transpose(g_densprofxz)*1e6
    teGENRAY = interpolate(teMPEX,(rGENRAY-rMin)/(rMax-rMin)*(nR-1),(zGENRAY-zMin)/(zMax-zMin)*(nZ-1),/grid)
    bMagGENRAY = transpose(g_bmodprofxz)
    
    plotRays = 0

    ; Plot up a single ray characteristics

    if plotRays then begin
    rayN = 8 
    rN = g_nrayelt[rayN]

    p=plot(g_wkPer[0:rN-1,rayN],layout=[2,4,1],/ylog,title='kPer,kPar')
    p=plot(g_wkPar[0:rN-1,rayN],/over,color='r')
    p=plot(-g_wkPar[0:rN-1,rayN],/over,color='b')
    p=plot(g_vpr[0:rN-1,rayN],layout=[2,4,2],title='vPhase_r',/current)
    p=plot(g_delpwr[0:rN-1,rayN]/max(g_delpwr[0:rN-1,rayN]),layout=[2,4,3],title='Power',/current)
    p=plot(gEx[0:rN-1,rayN],layout=[2,4,4],title='E*/E',/current,thick=2)
    p=plot(imaginary(gEx[0:rN-1,rayN]),/over)
    p=plot(gEy[0:rN-1,rayN],/over,thick=2,color='b')
    p=plot(imaginary(gEx[0:rN-1,rayN]),/over,color='b')
    p=plot(gEz[0:rN-1,rayN],/over,thick=2,color='r')
    p=plot(imaginary(gEz[0:rN-1,rayN]),/over,color='r')
    p=plot(g_ste[0:rN-1,rayN]*1e3,layout=[2,4,5],title='Temp [eV]',/current,thick=2)

    endif

    nRays = n_elements(g_wx[0,*])

    nL=25

    bMagUse = bMagMPEX
    neUse = neMPEX
    teUse = teMPEX
    rUse = rMPEX
    zUse = zMPEX

    ;bMagUse = bMagGENRAY
    ;neUse = neGENRAY
    ;teUse = teGENRAY
    ;rUse = rGENRAY
    ;zUse = zGENRAY

    wce = _e*bMagUse/_me
    wpe = sqrt(neUse*_e^2/(_me*_e0))
    L = 1-wpe^2/wrf^2 * wrf/(wrf+wce)
    R = 1-wpe^2/wrf^2 * wrf/(wrf-wce)
    P = 1-wpe^2/wrf^2
    wUH = sqrt(wpe^2+wce^2)
    n_X = sqrt(-(wpe^2-wrf*(wrf+wce))*(wpe^2-wrf*(wrf-wce)) / (wrf^2*(wpe^2+wce^2-wrf^2)))
    k_X = (n_X*wrf/_c)>0
    nPar_crit = (wce/wrf) / (1+wce/wrf)
    kPar_crit = nPar_crit*wrf/_c

    DB_kPar = 10000.0
    DB_nPar = DB_kPar*_c/wrf

    vTh = sqrt(2*teUse*_e/_me)

    rho_L = vTh / wce

    layout=[3,2]
    levels = 10.0^[15,16,17,18,19,20]
    colors = 100-bytscl(alog10(levels)-min(alog10(levels)))/3.+150
    c=contour(neUse,rUse-r0,zUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,1],xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])
    c=contour(wrf/wce,rUse-r0,zUse,c_value=[1,2,3,4,5],/over, rgb_table=0,c_color=0,c_label_show=1)
    c=contour(wrf/wUH,rUse-r0,zUse,c_value=[1],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(neUse,rUse-r0,zUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,3],/current,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])
    c=contour(wrf/wpe,rUse-r0,zUse,c_value=[0.5,1,2],/over, rgb_table=1,rgb_indices=[150,150,150],c_label_show=1)
    c=contour(wpe/wce,rUse-r0,zUse,c_value=[1,2,3,4,5,6],/over, rgb_table=8,c_color=[150],c_thick=2)

    c=contour(neUse,rUse-r0,zUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,4],/current,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])
    c=contour(L,rUse-r0,zUse,c_value=[0],/over, rgb_table=0,c_color=0,c_thick=2)
    c=contour(R,rUse-r0,zUse,c_value=[0],/over, rgb_table=1,c_color=[150],c_thick=2)
    c=contour(P,rUse-r0,zUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(k_X,rUse-r0,zUse,c_value=[0.1,1,10,100,1000],rgb_table=1,layout=[layout,2],/current,title='X Mode',xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])
    c=contour(wrf/wUH,rUse-r0,zUse,c_value=[1],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(neUse,rUse-r0,zUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,5],/current,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])
    c=contour(kPar_crit,rUse-r0,zUse,c_value=[0,10,50,100,200,300,400],/over, rgb_table=0,c_color=0,c_thick=1,c_label_show=1)
    c=contour(P,rUse-r0,zUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)

    ;v=vector(br2d_numeric,bz2d_numeric,rUse-r0,zUse,/auto_subsample,/auto_color,/auto_range,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])

    c=contour(neUse,rUse-r0,zUse,n_lev=15,c_value=10d0^[15,16,17,18,19,20],c_color=0,c_label_show=1,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])
    pp=plot(zUse,neUse[nR/2,*],/ylog)
    pp=plot(zUse,neUse[nR/2,*])

    c=contour(neUse,rUse-r0,zUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax])

    if plotRays then begin
    c=contour(transpose(alog10(g_pwr_e)),g_pwr_z+g_z0,g_pwr_r,/over)
    kLevels = 10^fIndGen(5)
    kColors = bytScl(kLevels,top=254)+1
    for ray=0,nRays-1 do begin
        pp=plot(g_ray_r[0:g_nrayelt[ray]-1,ray],g_ray_z[0:g_nrayelt[ray]-1,ray],$
                /over,xrange=[rMin,rMax]-r0,yRange=[zMin,zMax],$
                vert_colors=255-(bytSCl(alog10(((abs(g_wkPar[0:g_nrayelt[ray]-1,ray]))>100)<10.^5),top=254,min=2,max=5)+1),$
                rgb_table=1)
    endfor
    endif
    c=contour(L,rUse-r0,zUse,c_value=[0],/over, rgb_table=0,c_color=0,c_thick=2)
    c=contour(R,rUse-r0,zUse,c_value=[0],/over, rgb_table=1,c_color=[150],c_thick=2)
    c=contour(P,rUse-r0,zUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)
    c=contour(wrf/wUH,rUse-r0,zUse,c_value=[1],/over, rgb_table=7,c_color=[254],c_thick=2)
    c=contour(wrf/wce,rUse-r0,zUse,c_value=[1,2,3,4,5,6,7],/over, rgb_table=0,c_color=0,c_label_show=1,c_thick=2)

    for harm=1,6 do begin 
        wceDB_m = harm*wce/(1.0-3*DB_nPar*vTh/_c)
        wceDB_p = harm*wce/(1.0+3*DB_nPar*vTh/_c)
        c=contour(wrf/wceDB_m,rUse-r0,zUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
        c=contour(wrf/wceDB_p,rUse-r0,zUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
    endfor     

    if plotRays then begin
    minAbsZ = fltArr(nRays) 
    maxAbsKPar = fltArr(nRays)
    nBins=41
    hist = fltArr(nBins)
    minH = -0.06
    maxH = +0.06
    bins = fIndGen(nBins)/(nBins-1)*(maxH-minH)+minH
    for ray=0,nRays-1 do begin
        minAbsZ[ray] = min(abs(g_ray_z[0:g_nrayelt[ray]-1,ray]))
        maxAbsKPar[ray] = min(abs(g_wkPar[0:g_nrayelt[ray]-1,ray]))
        for i=1,g_nrayelt[ray]-1 do begin
            value = g_delpwr[i-1,ray]-g_delpwr[i,ray]
            hh = (g_ray_z[i,ray]+minH)/(maxH-minH)*(nBins-1)
            hist[hh] = hist[hh] + value
        endfor
    endfor

    endif


    ;   Extract a line domain

    r1 = 2.57+r0
    z1 = 0.12
    angle = 70
    len = 0.14
    nPts = nR

    dS = len/nPts

    r2 = +len * cos(angle*!dtor) + r1
    z2 = -len * sin(angle*!dtor) + z1

    dR = (r2-r1)/nPts
    dZ = (z2-z1)/nPts

    rPts = fIndGen(nPts)*dR + r1
    zPts = fIndGen(nPts)*dZ + z1

    p = plot([r1,r2]-r0,[z1,z2],thick=2,/over)
    p = plot(rPts-r0,zPts,thick=2,/over)

    rNorm = (rPts-rUse[0])/(rUse[-1]-rUse[0])*(n_elements(rUse)-1)
    zNorm = (zPts-zUse[0])/(zUse[-1]-zUse[0])*(n_elements(zUse)-1)

    neLine = interpolate(neUse,rNorm,zNorm)
    teLine = interpolate(teUse,rNorm,zNorm)

    brLine = interpolate(br2D_numeric,rNorm,zNorm)
    btLine = interpolate(bt2D_numeric,rNorm,zNorm)
    bzLine = interpolate(bz2D_numeric,rNorm,zNorm)

    lineUnitVec = [(r2-r1)/len,(z2-z1)/len]
    acrossUnitVec = [-lineUnitVec[1],lineUnitVec[0]]
    bAlongLine = fltArr(nPts)
    bAcrossLine = fltArr(nPts)
    for i=0,nPts-1 do begin
        bAlongLine[i] = lineUnitVec[0] * brLine[i] + lineUnitVec[1] * bzLine[i]
        bAcrossLine[i] = acrossUnitVec[0] * brLine[i] + acrossUnitVec[1] * bzLine[i]
    endfor

    ; Overwrite 2d arrays with new 1-D slice profiles

    rMin = r0
    rMax = r0+len
    zMin = -0.1
    zMax = +0.1

    print, 'rMin: ', rMin
    print, 'rMax: ', rMax
    print, 'zMin: ', zMin
    print, 'zMax: ', zMax

    nempex = rebin(neline,npts,nZ)
    teMPEX = rebin(teLine,nPts,nZ)

    br2d_numeric = rebin(bAlongLine,nPts,nZ)
    bt2d_numeric = rebin(btLine,nPts,nZ)
    bz2d_numeric = rebin(bAcrossLine,nPts,nZ)

stop
