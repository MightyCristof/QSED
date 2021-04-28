FUNCTION medabsdev, data, $
                    DIM = dim
                    

if n_elements(dim) eq 0 then dim = 0

medn = median(data,dim=dim)
sz = n_elements(medn)
mad = dblarr(sz)

case dim of
    2: for i = 0,sz-1 do mad[i] = median(abs(data[i,*]-medn[i]))
    1: for i = 0,sz-1 do mad[i] = median(abs(data[*,i]-medn[i]))
    0: mad[0] = median(abs(data-medn))
endcase

return, mad


END






