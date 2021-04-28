FUNCTION dlum, z, $
               SQ = sq


dl = dluminosity(z)         ;; dL in pc
dl *= !const.parsec * 1d2   ;; convert dL to cm

if keyword_set(sq) then dl = dl^2.

return, dl


END



