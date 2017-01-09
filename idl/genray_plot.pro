pro genray_plot

    @constants

    freq = 28d9
    wrf = 2*!pi*freq
    teMPEX = 2.0 ; eV

    r0 = 0
    g_z0 = 2.6
    radius = 0.06

	nX = 64 
    nY = 64 
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
stop
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

    DB_kPar = 10000.0
    DB_nPar = DB_kPar*_c/wrf

    vTh = sqrt(2*teUse*_e/_me)

    rho_L = vTh / wce

    xRange = [zMin,zMax]
    yRange = [xMin,xMax]

    layout=[3,2]
    levels = 10.0^[15,16,17,18,19,20]
    colors = 100-bytscl(alog10(levels)-min(alog10(levels)))/3.+150
    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,layout=[layout,1],xrange=xrange,yRange=yrange)
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

    ;v=vector(br2d_numeric,bz2d_numeric,xUse,zUse,/auto_subsample,/auto_color,/auto_range,xrange=[rMin,rMax],yRange=[zMin,zMax])

    c=contour(dens_xz,zUse,xUse,n_lev=15,c_value=10d0^[15,16,17,18,19,20],c_color=0,c_label_show=1,xrange=xrange,yRange=yrange)
    pp=plot(xUse,dens_xz[*,nX/2],/ylog)
    pp=plot(xUse,dens_xz[*,nX/2])

    c=contour(dens_xz,zUse,xUse,/fill,c_value=levels,rgb_indices=colors,rgb_table=3,xrange=xrange,yRange=yrange)

    c=contour(transpose(alog10(gr.pwr_e)),gr.pwr_z,gr.pwr_r,/over)

    kLevels = 10^fIndGen(5)
    kColors = bytScl(kLevels,top=254)+1
    for ray=0,nRays-1 do begin
        pp=plot(gr.ray_z[0:gr.nrayelt[ray]-1,ray],gr.ray_x[0:gr.nrayelt[ray]-1,ray],$
                /over,xrange=xrange,yRange=yrange,$
                vert_colors=255-(bytSCl(alog10(((abs(gr.wkPar[0:gr.nrayelt[ray]-1,ray]))>100)<10.^5),top=254,min=2,max=5)+1),$
                rgb_table=1)
    endfor

    c=contour(L,zUse,xUse,c_value=[0],/over, rgb_table=0,c_color=0,c_thick=2)
    c=contour(R,zUse,xUse,c_value=[0],/over, rgb_table=1,c_color=[150],c_thick=2)
    c=contour(P,zUse,xUse,c_value=[0],/over, rgb_table=7,c_color=[150],c_thick=2)
    c=contour(wrf/wUH,zUse,xUse,c_value=[1],/over, rgb_table=7,c_color=[254],c_thick=2)
    c=contour(wrf/wce,zUse,xUse,c_value=[1,2,3,4,5,6,7],/over, rgb_table=0,c_color=0,c_label_show=1,c_thick=2)

    for harm=1,6 do begin 
        wceDB_m = harm*wce/(1.0-3*DB_nPar*vTh/_c)
        wceDB_p = harm*wce/(1.0+3*DB_nPar*vTh/_c)
        c=contour(wrf/wceDB_m,zUse,xUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
        c=contour(wrf/wceDB_p,zUse,xUse,c_value=[1],/over, rgb_table=0,c_color=0,c_label_show=0)
    endfor     

    
    iix = ( sqrt(x3D^2+y3D^2) -xMinG)/(xMaxG-xMinG)*(nXG-1)
    iiz = ( z3D-zMinG)/(zMaxG-zMinG)*(nZG-1)

    ;dens3D = interpolate(gr.dens_xz,iix,iiz)
    ;v=volume(dens3D, volume_dimensions=[xmax-xmin,ymax-ymin,zmax-zmin], $
    ;    volume_location=[xmin,ymin,zmin], opacity_table0=fIndGen(256)^2/(256.0^2)*255 )
 
    range = [-1,1]*0.15 
    p=plot3d( (gr.ray_x[0,*])[*], (gr.ray_y[0,*])[*], (gr.ray_z[0,*])[*], $
        lineStyle='none', symbol='s', sym_size=1.0, $
        xRange = range, yRange = range, aspect_ratio=1.0, aspect_z = 1.0 )
    for n=0,n_elements(gr.ray_x[0,*])-1 do begin
        p=plot3d( gr.ray_x[*,n], gr.ray_y[*,n], gr.ray_z[*,n], /over , zTitle='Z' )
    endfor

stop
end