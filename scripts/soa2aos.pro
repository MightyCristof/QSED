FUNCTION soa2aos, table


;; check for input table
if (table eq !null) then return, !null

;; rip tag info
tags = tag_names(table)
ntags = n_elements(tags)

;; first tag array data shape for initialization
sz = size(table.(0))
;; length of output structure (always the last dimension for 2D and higher arrays)
;; if 1D array, histogram 1D arrays and pick highest value as length
if (sz[0] eq 1) then begin
    len = []
    for i = 0,ntags-1 do if ((size(table.(i)))[0] eq 1) then len = [len,(size(table.(i)))[1]]
    hr = histogram(len,locations=xr)
    hmax = max(hr,imax)
    len = xr[imax]
endif else len = sz[sz[0]]
             
for t = 0L,ntags-1L do begin
    sz = size(table.(t))
    case sz[0] of
        1:  begin
            if (t eq 0) then begin
                if (sz[sz[0]] eq len ) then new_table = create_struct(tags[0],(table.(0)[0])) else $
                                            new_table = create_struct(tags[0],table.(0))
            endif else begin
                if (sz[sz[0]] eq len) then new_table = create_struct(new_table,tags[t],(table.(t))[0]) else $
                                           new_table = create_struct(new_table,tags[t],table.(t))
            endelse
            end
        else: begin
            str = replicate('*',sz[0])
            str[-1] = '0'
            str = strjoin(str,',')
            if (t eq 0) then re = execute('new_table = create_struct(tags[0],(table.(0)['+str+']))') else $
                             re = execute('new_table = create_struct(new_table,tags[t],(table.(t))['+str+'])')
            endelse
    endcase
endfor

;; replicate full structure length
new_table = replicate(new_table,len)

;; procede through tags and fill data
for t = 0L,ntags-1L do new_table.(t) = table.(t)

return, new_table


END





;for t = 1L,ntags-1L do begin
;    sz = size(table.(t))
;    case sz[0] of
;        1:  begin
;            if (t eq 0) then new_table = create_struct(tags[0],(table.(0)[0])) else $
;                             if (sz[sz[0]] eq len) then new_table = create_struct(new_table,tags[t],(table.(t))[0]) else $
;                                                        new_table = create_struct(new_table,tags[t],table.(t))
;            end
;        else: begin
;            str = strarr(sz[0])+'*'
;            str[-1] = '0'
;            str = strjoin(str,',')
;            if (t eq 0) then re = execute('new_table = create_struct(tags[0],(table.(0)['+str+']))') else $
;                             re = execute('new_table = create_struct(new_table,tags[t],(table.(t))['+str+'])')
;            endelse
;    endcase
;endfor


















