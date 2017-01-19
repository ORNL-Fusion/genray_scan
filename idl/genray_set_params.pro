pro genray_set_params, $
        current = current, $
        rayTxt = rayTxt, $
        density = density, $ 
        T_eV = T_eV, $
        curr_multiplier = curr_multiplier 

    fileName = 'genray.in'
    backupFileName = 'genray.in.template'

    nLines = file_lines( fileName ) 

    data = strArr(nLines)

    ; Read genray.in

    openr, lun, fileName, /get_lun
    thisLine = ''
    for n=0,nLines-1 do begin
        readf, lun, thisLine
        data[n] = thisLine
    endfor
    free_lun, lun

    ; Template ray block is from line 1119:1128

    aboveRays = data[0:1119]
    rayBlock = data[1120:1128]
    belowRays = data[1129:-1] 

    ; Overwrite coil current
    ; ----------------------

    if keyword_set(current) then begin 

        lineNo = 180 

        nCoils = n_elements(curr_multiplier)

        print, 'UPDATING CURRENTS BLOCK'
        newStr = 'curc= '
        for _c=0,nCoils-1 do begin
            cStr = string( current*curr_multiplier[_c], format='(i3.3)')
            if _c eq nCoils-1 then begin
                newStr = newStr + cStr+'.d3'
            endif else begin
                newStr = newStr + cStr+'.d3,  '
            endelse
        endfor
        aboveRays[lineNo] = newStr

    endif


    ; Overwrite the rays with the rayTxt block
    ; ----------------------------------------

    if keyword_set(rayTxt) then begin

        rayBlock = rayTxt

        print, 'UPDATING RAYS BLOCK'

    endif

    template_wall_rmax = 0.09 ; I had to calibrate this number by hand. WTF. Where is it set? 

    if keyword_set(density) or keyword_set(T_eV) then begin

        nProfile = 64
        lineNo = 697

        newStr = ' ndens='+string(nProfile,format='(i3.3)')
        aboveRays[lineNo] = newStr

    endif

    ; Overwrite the density profile
    ; -----------------------------

    if keyword_set(density) then begin

        densStart = 143 ; Within the belowRays block 

        floor_ = density.floor_
        mag = density.mag 
        offset = density.offset 
        width = density.width 

        xMin = 0
        xMax = template_wall_rmax
        nX = nProfile
        x = fIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
        y = ( mag * ( 1.0 - tanh( (x*1e2-offset)/width ) ) ) > floor_ 

        densBlock = 'prof = '+string(y[0])
        for i=1,nX-1 do begin
            densBlock = [ densBlock, string(y[i]) ]
        endfor

        belowRays[densStart:densStart+nProfile-1] = densBlock

        save, x, y, fileName='density.sav'

    endif    


    ; Overwrite the Temp [keV] profile
    ; -----------------------------

    if keyword_set(T_eV) then begin

        tempStart = 237 ; Within the belowRays block 

        nTemp = nProfile

        T = fltArr(nTemp) + T_eV * 1e-3

        tempBlock = 'prof = '+string(T[0])
        for i=1,nTemp-1 do begin
            tempBlock = [ tempBlock, string(T[i]) ]
        endfor

        belowRays[tempStart:tempStart+nTemp-1] = tempBlock

    endif    

    data2 = [ aboveRays, rayBlock, belowRays ]

    nLines = n_elements(data2)

    ; Write new genray.in

    while file_test(backupFileName) do begin
        backupFileName = backupFileName+'_'
    endwhile

    file_move, fileName, backupFileName

    openw, lun, fileName, /get_lun
    for n=0,nLines-1 do begin
        printf, lun, data2[n]
    endfor
    free_lun, lun

end
