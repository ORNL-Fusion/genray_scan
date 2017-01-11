pro genray_scan, runGENRAY = _runGENRAY

    if keyword_set( _runGENRAY ) then runGENRAY = _runGENRAY else runGENRAY = 0

    @constants

    freq = 28e9
    t_eV = 4.0
    wrf = 2*!pi*freq
    
    T_eV = 10.0

	templateDir = 'template3'    	
    genrayBinary = expand_path( '~/code/genray-c/genray-c_160826.1/xgenray' )

    nC = 25
    curMin = 175
    curMax = 285
    curc = fIndGen(nC)/(nC-1)*(curMax-curMin)+curMin 

	nC = n_elements(curc)

    ; Beam parameters

    xOffSet = 0.10
    zOffSet = -0.02
    angle_deg = 20 
    width_m = 0.02
    spread_deg = 1.0
    rayDensity = 5

    run = 0

    cd, current = rootDir

    for c=0,nC-1 do begin

        ; Stage run

        thisDir = 'run'+string(run,format='(i3.3)')
        print, thisDir

        setParms = 0
        if runGENRAY then begin
        if file_test(thisDir,/directory) eq 0 

            file_delete, thisDir, /recursive, /allow_nonexistent
            file_copy, templateDir, thisDir, /recursive, /overwrite
            setParams = 1

        endif
        endif
        cd, thisDir

        ; Update run parameters 

        if runGENRAY then begin
        if setParams then begin
            
            genray_set_params, current = curc[c]

            rayTxt = genray_create_rays( xOffSet=xOffSet, zOffSet=zOffSet, $
                angle_deg=angle_deg, width_m=width_m, spread_deg=spread_deg, $
                rayDensity=rayDensity )  

            genray_set_params, rayTxt = rayTxt

            genray_set_params, /density

            genray_set_params, T_eV = T_eV

        endif
        endif

        ; Run genray

        if runGENRAY then begin
        if file_test( 'genray.nc') eq 0 then begin
            print, 'Running GENRAY ...'
            spawn, genrayBinary+' > genray.log', stdOut, stdErr
            print, 'DONE'
        endif else begin
            print, 'Nothing to do.'
        endelse
        endif

        ; Get ouput

        thisGR = genray_read_output()

        if size(pwr_e_percent,/type) eq 0 then begin
            nRPwr = n_elements(thisGR.spwr_rz_e[*,0])
            nRays = n_elements(thisGR.transm_ox)
            pwr_e_percent = fltArr(nRPwr,nC)
            pwr_i_percent = fltArr(nRPwr,nC)
            pwr_c_percent = fltArr(nRPwr,nC)
        endif

        launchedPower_kW = total( thisgr.delpwr[0,*] ) * 1e-7 / 1e3; covert from erg/s to J/s=W to kW

        pwr_e_percent[*,run] = total(thisgr.spwr_rz_e,2) * 1e-7 / 1e3 / launchedPower_kW * 100 
        pwr_i_percent[*,run] = total(thisgr.spwr_rz_i,2) * 1e-7 / 1e3 / launchedPower_kW * 100
        pwr_c_percent[*,run] = total(thisgr.spwr_rz_cl,2) * 1e-7 / 1e3 / launchedPower_kW * 100

        cd, rootDir
        ++run

    endfor

    r = thisGR.pwr_r
    z = thisGR.pwr_z

    maxLevel = 5.5e11
    scaling = 1

    levels = [-4,-3,-2,-1,0,1,2] 
    colors = 255-bytscl(levels,top=244)

    c=contour(alog10(pwr_e_percent),r,curc/40*1e3,c_value=levels,c_color=colors,rgb_table=3,xRange=[0.0,.08],$
            xtitle='r [m]', ytitle='Coil Current [A]', $
            title='Log10(Deposited Power) (e), Max Value: '+string(max(pwr_e_percent),format='(f4.1)')+'%', $
            layout=[3,1,1],/fill)
    c=contour(alog10(pwr_i_percent),r,curc/40*1e3,c_value=levels,c_color=colors,rgb_table=3,xRange=[0.0,.08],$
            xtitle='r [m]', ytitle='Coil Current [A]', title='Deposited Power (i)', layout=[3,1,2], /current)
    c=contour(alog10(pwr_c_percent),r,curc/40*1e3,c_value=levels,c_color=colors,rgb_table=3,xRange=[0.0,.08],$
            xtitle='r [m]', ytitle='Coil Current [A]', $
            title='Log10(Deposited Power) (e-i), Max Value: '+string(max(pwr_c_percent),format='(f4.1)')+'%', $
            layout=[3,1,3], /current,/fill)
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
    cdfId = ncdf_open ( fileList[0], /noWrite ) 
        nCdf_varGet, cdfId, 'spwr_rz_e', g_pwr_e
        nCdf_varGet, cdfId, 'transm_ox', g_transm_ox
    ncdf_close, cdfId

    nRPwr = n_elements(g_pwr_e[*,0])
    nRays = n_elements(g_transm_ox)

    resCurC5 = [175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,255,260,265,285]
    res5 = [0,0.02,0.035,0.055,0.065,0.075,0.085,0.09]

    resCurC4 = [220,225,230,235,240,245,255,260,265,285]
    res4 = [0,0.02,0.035,0.045,0.06,0.07,0.085,0.09,0.095,0.115]

    res3 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    over=0
    over2=0

    curc = fltArr(nFiles)

    pwrall = fltArr(nRPwr,nFiles)
    transAll = fltArr(nRays,nFiles)

    nDensPwr = 30
    dMin = 1e17
    dMax = 6e19
    dGrid = fIndGen(nDensPwr)/(nDensPwr-1)*(dMax-dMin)+dMin
    densityPwrAll = fltArr(nDensPwr,n_elements(curc))

    for f=0,nFiles-1 do begin

        curc[f] = float(thisCurc)

    	cdfId = ncdf_open ( fileList[f], /noWrite ) 
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

        rGENRAY = g_eqdsk_x
        zGENRAY = g_eqdsk_z

        g_ray_r = g_wz*1e-2
        g_ray_z = g_wx*1e-2

        g_wkPer = g_wnper*wrf/_c
        g_wkPar = g_wnpar*wrf/_c
        g_wkMag = sqrt(g_wkPer^2+g_wkPar^2)
        g_kr = g_wn_r*wrf/_c
        g_vpr = wrf/g_kr

        neGENRAY = transpose(g_densprofxz)*1e6
        teGENRAY = neGENRAY*0+t_eV 
        bMagGENRAY = transpose(g_bmodprofxz)
            
        print, fileList[f], size(g_pwr_e)

        r = g_pwr_r
        z = g_pwr_z

        thisi = (r-rGENRAY[0])/(rGENRAY[-1]-rGENRAY[0])*(n_elements(rGENRAY)-1)
        thisj = (z-zGENRAY[0])/(zGENRAY[-1]-zGENRAY[0])*(n_elements(zGENRAY)-1)
        thisDensity = interpolate(g_densprofxz*1e6,thisi,thisj,/grid)

        for i=0,n_elements(r)-1 do begin
            for j=0,n_elements(z)-1 do begin
                thisPwr = g_pwr_e[i,j]
                dIndex = (thisDensity[i,j]-dMin)/(dMax-dMin)*(nDensPwr-1)
                densityPwrAll[dIndex,f] = densityPwrAll[dIndex,f] + thisPwr
                ;print, thisDensity, dIndex, nDensPwr
                ;print, r[i], z[j], thisi, thisj, n_elements(rGENRAY),n_elements(zGENRAY),thisDensity
            endfor
        endfor

        print, min(r), max(r)
        thisPwr = total(g_pwr_e,2)

        pwrAll[*,f] = thisPwr

        p=plot(r,thisPwr,over=over,layout=[2,1,1],title='r Grid',/current)    
        over = (over + 1)<1
        p=plot(dGrid,densityPwrAll[*,f],over=over2,layout=[2,1,2],/current,title='Density Grid')    
        over2 = (over2 + 1)<1

        transAll[*,f] = g_transm_ox

    endfor

    maxLevel = 1.5e11
    scaling = 2

    nLevs = 20
    levelsA = (fIndGen(nLevs)/(nLevs-1))
    levels = levelsA^scaling
    levels = levels/max(levels)*maxLevel
    colors = 255-(bytScl(levelsA,top=254)+1)

    c=contour(pwrAll,r,curc,c_value=levels,c_color=colors,/fill,rgb_table=3,xRange=[0.0,.08],$
            xtitle='r [m]', ytitle='Coil Current [A]', title='(Unmapped)')
    p=plot(res4,resCurC4,/over,thick=2)
    p=plot(res5,resCurC5,/over,thick=2)

    c.save, 'pscan.png', resolution=300, /transpar

    nR = n_elements(rGENRAY) 
    nZ = n_elements(zGENRAY)
    iiP = where(rGENRAY ge 0)
    rDensityGrid = interpol(rGENRAY[iiP],g_densprofxz[iiP,nZ/2]*1e6,dGrid)

    c=contour(densityPwrAll,rDensityGrid,curc,c_value=levels,c_color=colors,/fill,rgb_table=3,xRange=[0.0,.08],$
            xtitle='r [m]', ytitle='Coil Current [A]', title= '(Mapped)')
    p=plot(res4,resCurC4,/over,thick=2)
    p=plot(res5,resCurC5,/over,thick=2)

    c.save, 'pscanMapped.png', resolution=300, /transpar

    nLevs = 10
    levels = fIndGen(nLevs)/(nLevs-1)
    colors = 255-(bytScl(levels,top=254)+1)
    c=contour(transAll,fIndGen(nRays),curc,xTitle='Ray Number',yTitle='Coil Current [A]',c_value=levels,c_color=colors,/fill,$
            rgb_table=3,position=[0.15,0.15,0.7,0.9])
    cb = colorbar(TITLE='Transmission P(O)/P(X)',orientation=1,taper=0,tickInterval=0.2,textPos=1,tickFormat='(f3.1)',border=1)

stop 
end
