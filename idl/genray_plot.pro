pro genray_plot

    @constants

    freq = 28d9
    wrf = 2*!pi*freq
    teMPEX = 20.0 ; eV

    r0 = 0
    g_z0 = 2.6
    radius = 0.06

	nX = 128 
    nY = 128 
	nZ = 128 

	xMin = -radius*2
	xMax = +radius*2

   	yMin = -radius*2
	yMax = +radius*2

	zMin = r0+2.3
	zMax = r0+2.9

    xTmp=fIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
    yTmp=fIndGen(nY)/(nY-1)*(yMax-yMin)+yMin
    zTmp=fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

    xMPEX = xTmp
    yMPEX = yTmp
    zMPEX = zTmp

    x2D = rebin(xMPEX,nX,nZ)
    z2D = transpose(rebin(zMPEX,nZ,nX))

    x3D = rebin(xMPEX,nX,nY,nZ)
    y3D = transpose(rebin(yMPEX,nY,nX,nZ),[1,0,2])
    z3D = transpose(rebin(zMPEX,nZ,nX,nY),[1,2,0])

    gr = genray_read_output()

    ; Add offset

    gr.z = gr.z + g_z0
    gr.pwr_z = gr.pwr_z + g_z0  
    gr.ray_z = gr.ray_z + g_z0  

    ; Plot up a single ray characteristics

    rayN = 0 
    rN = gr.nrayelt[rayN]

    layout = [2,2]
    p=plot(gr.wkPer[0:rN-1,rayN],layout=[layout,1],/ylog,title='kPer,kPar')
    p=plot(gr.wkPar[0:rN-1,rayN],/over,color='r')
    p=plot(-gr.wkPar[0:rN-1,rayN],/over,color='b')
    p=plot(gr.vpr[0:rN-1,rayN],layout=[layout,2],title='vPhase_r',/current)
    p=plot(gr.delpwr[0:rN-1,rayN]/max(gr.delpwr[0:rN-1,rayN]),layout=[layout,3],title='Power',/current)
    p=plot(gr.Ex[0:rN-1,rayN],layout=[layout,4],title='E*/E',/current,thick=2)
    p=plot(imaginary(gr.Ex[0:rN-1,rayN]),/over)
    p=plot(gr.Ey[0:rN-1,rayN],/over,thick=2,color='b')
    p=plot(imaginary(gr.Ex[0:rN-1,rayN]),/over,color='b')
    p=plot(gr.Ez[0:rN-1,rayN],/over,thick=2,color='r')
    p=plot(imaginary(gr.Ez[0:rN-1,rayN]),/over,color='r')

    nRays = n_elements(gr.ray_x[0,*])

    nL=25

    nXG = n_elements(gr.x)
    xMinG = gr.x[0]
    xMaxG = gr.x[-1]

    nZG = n_elements(gr.z)
    zMinG = gr.z[0]
    zMaxG = gr.z[-1]

    iix = (xMPEX-xMinG)/(xMaxG-xMinG)*(nXG-1)
    iiz = (zMPEX-zMinG)/(zMaxG-zMinG)*(nZG-1)

    bMag_xz = transpose(interpolate(gr.bMag_xz,iix,iiz,/grid))
    dens_xz = transpose(interpolate(gr.dens_xz,iix,iiz,/grid))

    teUse = teMPEX
    xUse = xMPEX
    zUse = zMPEX

    wce = _e*bMag_xz/_me
    wpe = sqrt(dens_xz*_e^2/(_me*_e0))
    L = 1-wpe^2/wrf^2 * wrf/(wrf+wce)
    R = 1-wpe^2/wrf^2 * wrf/(wrf-wce)
    P = 1-wpe^2/wrf^2
    wUH = sqrt(wpe^2+wce^2)
    n_X = sqrt(-(wpe^2-wrf*(wrf+wce))*(wpe^2-wrf*(wrf-wce)) / (wrf^2*(wpe^2+wce^2-wrf^2)))
    k_X = (n_X*wrf/_c)>0
    nPar_crit = (wce/wrf) / (1+wce/wrf)
    kPar_crit = nPar_crit*wrf/_c

    DB_kPar = 2000.0
    DB_nPar = DB_kPar*_c/wrf

    vTh = sqrt(2*teUse*_e/_me)

    rho_L = vTh / wce

    launchedPower_kW = total( gr.delpwr[0,*] ) * 1e-7 / 1e3; covert from erg/s to J/s=W to kW

    pow_e_percent_2D = gr.spwr_rz_e * 1e-7 / 1e3 / launchedPower_kW * 100 
    pow_i_percent_2D = gr.spwr_rz_i * 1e-7 / 1e3 / launchedPower_kW * 100
    pow_c_percent_2D = gr.spwr_rz_cl * 1e-7 / 1e3 / launchedPower_kW * 100

    pow_e_percent = total(gr.spwr_rz_e,2) * 1e-7 / 1e3 / launchedPower_kW * 100 
    pow_i_percent = total(gr.spwr_rz_i,2) * 1e-7 / 1e3 / launchedPower_kW * 100
    pow_c_percent = total(gr.spwr_rz_cl,2) * 1e-7 / 1e3 / launchedPower_kW * 100

    totalAbsorbedPower_e = total(pow_e_percent) 
    totalAbsorbedPower_i = total(pow_i_percent) 
    totalAbsorbedPower_c = total(pow_c_percent) 

    totalAbsorbedPower = totalAbsorbedPower_e + totalAbsorbedPower_i + totalAbsorbedPower_c  
    totalReflectedPower = 100-totalAbsorbedPower
 
    xRange = [zMin,zMax]
    yRange = [xMin,xMax]

    layout=[3,2]
    levels = 10.0^[17.1,18,19,20]
    colors = 100-bytscl(alog10(levels)-min(alog10(levels)))/3.+150

    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,$
            rgb_table=3,layout=[layout,1],xrange=xrange,yRange=yrange)
    c=contour(wrf/wce,zUse,xUse,c_value=[1,2,3,4,5],/over, rgb_table=0,c_color=0,c_label_show=1)
    c=contour(wrf/wUH,zUse,xUse,c_value=[1],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,3],/current,xrange=xrange,yRange=yrange)
    c=contour(wrf/wpe,zUse,xUse,c_value=[0.5,1,2],/over, rgb_table=1,rgb_indices=[150,150,150],c_label_show=1)
    c=contour(wpe/wce,zUse,xUse,c_value=[1,2,3,4,5,6],/over, rgb_table=8,c_color=[150],c_thick=2)

    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,4],/current,xrange=xrange,yRange=yrange)
    c=contour(L,zUse,xUse,c_value=[0],/over, rgb_table=0,c_color=0,c_thick=2)
    c=contour(R,zUse,xUse,c_value=[0],/over, rgb_table=1,c_color=[150],c_thick=2)
    c=contour(P,zUse,xUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(k_X,zUse,xUse,c_value=[0.1,1,10,100,1000],rgb_table=1,layout=[layout,2],/current,title='X Mode',xrange=xrange,yRange=yrange)
    c=contour(wrf/wUH,zUse,xUse,c_value=[1],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,5],/current,xrange=xrange,yRange=yrange)
    c=contour(kPar_crit,zUse,xUse,c_value=[0,10,50,100,200,300,400],/over, rgb_table=0,c_color=0,c_thick=1,c_label_show=1)
    c=contour(P,zUse,xUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)

    c=contour(dens_xz,zUse,xUse,n_lev=15,c_value=10d0^[15,16,17,18,19,20],c_color=0,c_label_show=1,xrange=xrange,yRange=yrange)
    pp=plot(xUse,dens_xz[*,nX/2],/ylog)
    pp=plot(xUse,dens_xz[*,nX/2])


    ; Single panel plot with density
    ; ------------------------------

    _fs = 18
    _dim=[1500,800]

    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,$
            xrange=xrange,yRange=yrange,aspect_ratio=1,title='Ray Trajectories and Cutoffs',$
            xTitle='z [m]', yTitle='r [m]',dim=_dim, font_size=_fs)

    kLevels = 10^fIndGen(5)
    kColors = bytScl(kLevels,top=254)+1
    for ray=0,nRays-1 do begin
        thisRay_r = sqrt( (gr.ray_x[0:gr.nrayelt[ray]-1,ray])^2 + (gr.ray_y[0:gr.nrayelt[ray]-1,ray])^2 )
        pp=plot(gr.ray_z[0:gr.nrayelt[ray]-1,ray], thisRay_r,$
                /over,xrange=xrange,yRange=yrange,$
                vert_colors=255-(bytSCl(alog10(((abs(gr.wkPar[0:gr.nrayelt[ray]-1,ray]))>100)<10.^5),top=254,min=2,max=5)+1),$
                rgb_table=1)
    endfor
    c=contour(L,zUse,xUse,c_value=[0],/over, rgb_table=0,c_color=0,c_thick=2)
    c=contour(R,zUse,xUse,c_value=[0],/over, rgb_table=1,c_color=[150],c_thick=2)
    cp=contour(P,zUse,xUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)
    contour, p, zUse, xUse, levels=[0],path_xy=cp_path, /path_data_coords
    c=contour(wrf/wUH,zUse,xUse,c_value=[1],/over, rgb_table=7,c_color=[254],c_thick=2)
    c=contour(wrf/wce,zUse,xUse,c_value=[1,2,3,4,5,6,7],/over, rgb_table=0,c_color=0,c_label_show=1,c_thick=2)

    for harm=1,6 do begin 
        wceDB_m = harm*wce/(1.0-3*DB_nPar*vTh/_c)
        wceDB_p = harm*wce/(1.0+3*DB_nPar*vTh/_c)
        c=contour(wrf/wceDB_m,zUse,xUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
        c=contour(wrf/wceDB_p,zUse,xUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
    endfor     

    c.save, 'density-cutoffs.png', resolution=300

    ; --------------------------


    ; Single panel plot with no density
    ; ---------------------------------
    cB=contour(dens_xz,zUse,xUse,/nodata,/fill,c_value=levels,rgb_indices=colors,$
            rgb_table=3,xrange=xrange,yRange=yrange,aspect_ratio=1,title='log10( Power Absorbtion (e) )',$
            xTitle='z [m]', yTitle='r [m]',dim=_dim, font_size=_fs)

    _nL = 30
    _min = -15
    _max = 2
    _levels = fIndGen(_nL)/(_nL-1)*(_max-_min)+_min
    _colors = 255-bytscl(_levels,top=244)
 
    cB=contour(transpose(alog10(pow_e_percent_2D)),gr.pwr_z,gr.pwr_r,/over,/fill,c_value=_levels,c_color=_colors,rgb_table=3)

    kLevels = 10^fIndGen(5)
    kColors = bytScl(kLevels,top=254)+1
    for ray=0,nRays-1 do begin
        thisRay_r = sqrt( (gr.ray_x[0:gr.nrayelt[ray]-1,ray])^2 + (gr.ray_y[0:gr.nrayelt[ray]-1,ray])^2 )
        pp=plot(gr.ray_z[0:gr.nrayelt[ray]-1,ray], thisRay_r,$
                /over,xrange=xrange,yRange=yrange,$
                vert_colors=255-(bytSCl(alog10(((abs(gr.wkPar[0:gr.nrayelt[ray]-1,ray]))>100)<10.^5),top=254,min=2,max=5)+1),$
                rgb_table=1,thick=2)
    endfor
    cB=contour(L,zUse,xUse,c_value=[0],/over, rgb_table=0,c_color=0,c_thick=2)
    cB=contour(R,zUse,xUse,c_value=[0],/over, rgb_table=1,c_color=[150],c_thick=2)
    cp=contour(P,zUse,xUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)
    contour, p, zUse, xUse, levels=[0],path_xy=cp_path, /path_data_coords
    c=contour(wrf/wUH,zUse,xUse,c_value=[1],/over, rgb_table=7,c_color=[254],c_thick=2)
    c=contour(wrf/wce,zUse,xUse,c_value=[1,2,3,4,5,6,7],/over, rgb_table=0,c_color=0,c_label_show=1,c_thick=2)

    for harm=1,6 do begin 
        wceDB_m = harm*wce/(1.0-3*DB_nPar*vTh/_c)
        wceDB_p = harm*wce/(1.0+3*DB_nPar*vTh/_c)
        c=contour(wrf/wceDB_m,zUse,xUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
        c=contour(wrf/wceDB_p,zUse,xUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
    endfor     

    cB.save, 'ray-absorbtion.png', resolution=300

    ; --------------------------


    iix = ( sqrt(x3D^2+y3D^2) -xMinG)/(xMaxG-xMinG)*(nXG-1)
    iiz = ( z3D-zMinG)/(zMaxG-zMinG)*(nZG-1)

    range = [-1,1]*0.15 
    p=plot3d( (gr.ray_x[0,*])[*], (gr.ray_y[0,*])[*], (gr.ray_z[0,*])[*], $
        lineStyle='none', symbol='s', sym_size=1.0, $
        xRange = range, yRange = range, aspect_ratio=1.0, aspect_z = 1.0 )
    for n=0,n_elements(gr.ray_x[0,*])-1 do begin
        p=plot3d( gr.ray_x[0:gr.nrayelt[n]-1,n], $
            gr.ray_y[0:gr.nrayelt[n]-1,n], gr.ray_z[0:gr.nrayelt[n]-1,n]$
            , /over , zTitle='Z', thick=2 )
    endfor

    nPlasma = 12 
    theta = fIndGen(nPlasma)/nPlasma*360
    plasma_outline_r = cp_path[1,*]
    plasma_outline_z = cp_path[0,*]
    iiKeep = where( abs(plasma_outline_r) lt 0.05 and plasma_outline_r gt 0.001 )
    _r = plasma_outline_r[iiKeep]
    _z = plasma_outline_z[iiKeep]
    for n=0,nPlasma-1 do begin
        _x = _r * cos(theta[n])
        _y = _r * sin(theta[n])
        p=plot3d(_x,_y,_z,/over,thick=3,color='r',transparency=60)
    endfor

    c = contour(transpose(gr.spwr_rz_e), gr.pwr_z, gr.pwr_r, layout=[1,3,1],/fill,title='power absorped (e)')
    c = contour(transpose(gr.spwr_rz_i), gr.pwr_z, gr.pwr_r, layout=[1,3,2],/current,/fill,title='power absorped (i)')
    c = contour(transpose(gr.spwr_rz_cl), gr.pwr_z, gr.pwr_r, layout=[1,3,3],/current,/fill,title='power absorped (cl)')

    nL = string(10B)
    p=plot(gr.pwr_r,total(pow_e_percent,/cum),$
            title='Absorped Power: '+string(totalAbsorbedPower,format='(f4.1)')+'%, ' $
            +'Reflected Power: '+string(totalReflectedPower,format='(f4.1)')+'%'+nL+nL $
            +' (e: '+string(totalAbsorbedPower_e,format='(f4.1)')+'%, ' $
            +' i: '+string(totalAbsorbedPower_i,format='(f4.1)')+'%, ' $
            +' c: '+string(totalAbsorbedPower_c,format='(f4.1)')+'%)', $
            yTitle = 'Absorbed Power', xTitle='r [m]',thick=2, color='b')
    p=plot(gr.pwr_r,total(pow_i_percent,/cum),color='r',/over,thick=2)
    p=plot(gr.pwr_r,total(pow_c_percent,/cum),/over,thick=2)
    thisDensity = interpol( gr.dens_xy[*,n_elements(gr.x)/2]/max(gr.dens_xy[*,n_elements(gr.x)/2])*100, gr.x, gr.pwr_r )
    pp=plot(gr.pwr_r,thisDensity,/over)

    ; Plot kPar 

    p=plot( gr.wkpar[0:gr.nrayelt[0]-1,0], zTitle='kPar')
    for n=1,n_elements(gr.ray_x[0,*])-1 do begin
        p=plot( gr.wkpar[0:gr.nrayelt[n]-1,n], /over)
    endfor

   

stop
end
