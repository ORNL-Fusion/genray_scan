pro genray_set_params, current = current, rayTxt = rayTxt

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

    ; Parse genray.in and modify

    for n=0,nLines-1 do begin

        ; Overwrite coil current

        if keyword_set(current) then begin 

            templateStr = ' curc=   264.d3,  264.d3,  264.d3,  264.d3 ! NEW ! [A] ! For model_b=2 only'

            if strCmp( data[n], templateStr ) then begin

                print, 'UPDATING CURRENTS BLOCK'
                cStr = string( current, format='(i3.3)')
                newStr =  ' curc=   '+cStr+'.d3,  '+cStr+'.d3,  '+cStr+'.d3,  '+cStr+'.d3'
                data[n] = newStr

            endif 

        endif

    endfor

    ; Overwrite the rays with the rayTxt block

    if keyword_set(rayTxt) then begin

        ; Template ray block is from line 1119:1128

        top = data[0:1118]
        bottom = data[1129:-1] 

        data = [top,rayTxt,bottom]

        print, 'UPDATING RAYS BLOCK'
stop
    endif

    nLines = n_elements(data)

    ; Write new genray.in

    while file_test(backupFileName) do begin
        backupFileName = backupFileName+'_'
    endwhile

    file_move, fileName, backupFileName

    openw, lun, fileName, /get_lun
    for n=0,nLines-1 do begin
        printf, lun, data[n]
    endfor
    free_lun, lun

end
