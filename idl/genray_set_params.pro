pro genray_set_params, current = current, rayTxt = rayTxt, density = density, T_eV = T_eV 

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

        print, 'UPDATING CURRENTS BLOCK'
        cStr = string( current, format='(i3.3)')
        cStrMod = string( current+1, format='(i3.3)')
        newStr =  ' curc= 0.d3,  '+cStrMod+'.d3,  '+cStr+'.d3,  '+cStr+'.d3'
        aboveRays[lineNo] = newStr

    endif


    ; Overwrite the rays with the rayTxt block
    ; ----------------------------------------

    if keyword_set(rayTxt) then begin

        rayBlock = rayTxt

        print, 'UPDATING RAYS BLOCK'

    endif

    template_wall_rmax = 0.12  

    ; Overwrite the density profile
    ; -----------------------------

    if keyword_set(density) then begin

        densStart = 143 ; Within the belowRays block 

        nDens = 64

        floor_ = 1e17
        mag = 9e18
        offset = 3.2 ; cm
        width = 0.38 ; cm
        
        xMin = 0
        xMax = template_wall_rmax
        nX = nDens
        x = fIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
        y = ( mag * ( 1.0 - tanh( (x*1e2-offset)/width ) ) ) > floor_ 

        densBlock = 'prof = '+string(y[0])
        for i=1,nX-1 do begin
            densBlock = [ densBlock, string(y[i]) ]
        endfor

        belowRays[densStart:densStart+nDens-1] = densBlock

    endif    


    ; Overwrite the Temp [keV] profile
    ; -----------------------------

    if keyword_set(T_eV) then begin

        tempStart = 237 ; Within the belowRays block 

        nTemp = 64

        template_wall_rmax = 0.12  

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
