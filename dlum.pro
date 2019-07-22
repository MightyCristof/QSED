FUNCTION dlum, z, $
               SQUARED = squared


testz = 10.^(dindgen(150)/100.)-(1.-min(z)>0.)      ;; range of z values to calculate dL
dl = lumdist(testz,h0=70,omega_m=0.3,lambda0=0.7)   ;; dL in Mpc => dl = (1+z)c ¼dz/H(z)
dl = interpol(dl,testz,z)
dl *= 1e6 * !const.parsec * 1e2                     ;; dL converted from Mpc to cm
if keyword_set(squared) then dl = dl^2.             ;; dL^2

return, dl


END



