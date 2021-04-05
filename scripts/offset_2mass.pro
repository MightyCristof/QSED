PRO offset_2mass, filt, $
                  flx, $
                  err, $
                  gobs


;; results of POLY_FIT.PRO against output UK-2MASS fluxes
rel = [[-0.85019920,1.6276533,-0.095407501], $
       [-0.86744305,1.6064217,-0.084541737], $
       [-1.1354586,1.7149793,-0.094379308]]
       
;; normalize UKIDSS and 2MASS data
i2m = where(strmatch(filt,'TWOM*'),ct2m)
iuk = where(strmatch(filt,'UK*'),ctuk)
if (ct2m eq 3 and ctuk eq 4) then begin
    iuk = iuk[1:-1]
    ;iboth = where(total(obs.bin[[i2m,iuk]],1) eq 6,ct)
    ;if (ct eq 0) then stop
    for i = 0,ct2m-1 do begin
        iboth = where(total(gobs[[i2m[i],iuk[i]],*],1) eq 2,ct)
        flx_2m = reform(alog10(flx[i2m[i],iboth]))
        var_2m = reform((err[i2m[i],iboth]/(flx[i2m[i],iboth]*alog(10)))^2)
        flx_uk = reform(alog10(flx[iuk[i],iboth]))
        var_uk = reform((err[iuk[i],iboth]/(flx[iuk[i],iboth]*alog(10)))^2)
        ;; UKIDSS-2MASS relation
        ;rel2m = poly_fit(flx_2m,flx_uk,2)
        ;fit2m = rel2m[0] + rel2m[1]*flx_2m + rel2m[2]*flx_2m^2
        fitx = flx_2m[sort(flx_2m)]
        fity = rel[0,i] + rel[1,i]*fitx + rel[2,i]*fitx^2
        ;; one-to-one line
        onex = flx_uk[sort(flx_uk)]
        oney = onex
        ;; match relation to one-to-one line
        loc = value_locate(oney,fity)
        if (n_elements(loc) eq 1 and (loc)[0] eq -1) then stop
        ;; difference in 2MASS flux
        diff = onex[loc] - fitx
        ;; index of relation turnover
        rmin = min(abs(onex[loc]-fitx),imin)
        ;; average offsets
        stat = moment(diff)
        stat_lo = moment(diff[0:imin])
        stat_hi = moment(diff[imin:-1])
        temp_flx = flx_2m
        temp_var = var_2m
        ilo = where(flx_2m le (4.6<fitx[imin]),ct)
        if (ct gt 0.) then begin
            temp_flx[ilo] += stat_lo[0]
            temp_var[ilo] += stat_lo[1]
        endif
        ihi = where(flx_2m gt fitx[imin],ct)
        if (ct gt 0.) then begin
            temp_flx[ihi] += stat_hi[0]
            temp_var[ihi] += stat_hi[1]
        endif
        corr_flx = 10.^temp_flx
        corr_err = sqrt(corr_flx^2*alog(10)^2*temp_var)
        flx[i2m[i],iboth] = corr_flx
        err[i2m[i],iboth] = corr_err
    endfor
endif
   
   
END



     
;        ;; difference between 2MASS and relation
;        diff = fit - flx_2m
;        ;; difference between relation and one-to-one line
;        reld = fit - flx_uk
;        ;; only adjust sources close to one-to-one line
;        resistant_mean,diff,2.0,mn,sigmn,nrej,goodvec=iadj
;        ;; first flux correction
;        corr = flx_2m
;        corr[iadj] += reld[iadj]
;        p = plot(flx_uk,flx_2m,'.',transparency=95,xra=[1,7],yra=[1,7],aspect_ratio=1,title=filt[i2m[i]])
;        p = plot(flx_uk,corr,'.r',/ov)
;        p = plot(p.xra,p.yra,'--r',/ov)        
;        ;; difference between unadjusted one-to-one line and 2MASS corrected fluxes
;        ileft = exclude(flx_2m,iadj)
;        diff = flx_uk[ileft] - flx_2m[ileft]
;        resistant_mean,diff,1.0,mn,sigmn,nrej,goodvec=ivalid
;
;
;
;        ileft = exclude(flx_2m,iadj)
;        ;; STDDEV of difference from relation
;        sigm = stddev(diff)
;        ;; scatter 2MASS around one-to-one line
;        left = flx_uk[ileft] - randomn(seed,n_elements(ileft))*sigm
;        corr = flx_uk - randomn(seed,n_elements(flx_uk))*sigm


        


