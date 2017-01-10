function genray_read_output

    @constants

    genrayFileName = 'genray.nc' 
    cdfId = ncdf_open ( genrayFileName, /noWrite ) 
        nCdf_varGet, cdfId, 'freqcy', g_freqcy 
        nCdf_varGet, cdfId, 'densprofxz', g_densprofxz 
        nCdf_varGet, cdfId, 'densprofyz', g_densprofyz 
        nCdf_varGet, cdfId, 'densprofxy', g_densprofxy 
        nCdf_varGet, cdfId, 'temprof', g_temprof 
        nCdf_varGet, cdfId, 'wx', g_wx 
        nCdf_varGet, cdfId, 'wy', g_wy 
        nCdf_varGet, cdfId, 'wz', g_wz 
        nCdf_varGet, cdfId, 'nrayelt', g_nrayelt 
        nCdf_varGet, cdfId, 'w_x_densprof_nc', g_w_x_densprof_nc 
        nCdf_varGet, cdfId, 'eqdsk_x', g_eqdsk_x 
        nCdf_varGet, cdfId, 'eqdsk_y', g_eqdsk_y 
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
        nCdf_varGet, cdfId, 'transm_ox', g_transm_ox

	ncdf_close, cdfId

    wrf = 2 * !pi * g_freqcy

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

    g_ray_x = g_wx*1e-2
    g_ray_y = g_wy*1e-2
    g_ray_z = g_wz*1e-2

    g_wkPer = g_wnper*wrf/_c
    g_wkPar = g_wnpar*wrf/_c
    g_wkMag = sqrt(g_wkPer^2+g_wkPar^2)
    g_kr = g_wn_r*wrf/_c
    g_vpr = wrf/g_kr

    bMag_xz = g_bmodprofxz

    ne_xz = g_densprofxz*1e6
    ne_yz = g_densprofyz*1e6
    ne_xy = g_densprofxy*1e6

    return, { $
        freq : g_freqcy, $
        Ex : gEx, $
        Ey : gEy, $
        Ez : gEz, $
        EVecMag : gEVecMag, $
        BVecMag : gBVecMag, $ 
        Bx : gBx, $ 
        By : gBy, $ 
        Bz : gBz, $ 
        EDotB : gEDotB, $ 
        Angle : gAngle, $
        x : g_eqdsk_x, $ 
        y : g_eqdsk_y, $ 
        z : g_eqdsk_z, $
        wkPer : g_wkPer, $ 
        wkPar : g_wkPar, $ 
        wkMag : g_wkMag, $ 
        kr : g_kr, $
        vpr : g_vpr, $
        dens_xz : ne_xz, $ 
        dens_yz : ne_yz, $ 
        dens_xy : ne_xy, $ 
        bMag_xz :bMag_xz, $
        nrayelt : g_nrayelt, $ 
        pwr_e : g_pwr_e, $
        pwr_i : g_pwr_i, $
        transm_ox : g_transm_ox, $
        pwr_r : g_pwr_r, $ 
        pwr_z : g_pwr_z, $
        delpwr : g_delpwr, $
        wx : g_wx, $
        wy : g_wy, $
        wz : g_wz, $
        ray_x : g_ray_x, $
        ray_y : g_ray_y, $
        ray_z : g_ray_z, $
        spwr_rz_e : g_pwr_e, $
        spwr_rz_i : g_pwr_i, $
        spwr_rz_cl : g_pwr_cl }

end
