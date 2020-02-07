;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   read_sed_phot
;   
; PURPOSE:
;   Read photometry from .FITS file and format for SED modeling.
;   
; CALLING SEQUENCE:
;   read_sed_phot, file, [, /MASK, /NIR, /DERED, /FORCED_PHOT, /CORR_2MASS, /MIN_ERR ]
;
; INPUTS:
;   file			- String containing the name of the photometry data file.
;	
; OPTIONAL INPUTS:
;   /MASK           - Use mangle mask to remove sources in WISE near bright stars.
;   /NIR            - Only keep sources with at least 1 NIR detection.
;   /DERED          - Correct fluxes for Galactic extinction.
;   /FORCED_PHOT    - Replace AllWISE with unWISE forced photometry where S/N is greater.
;   /CORR_2MASS     - If using UKIDSS & 2MASS, add a correction fator to 2MASS data.
;   /MIN_ERR        - Require a minimum error of ±0.05 mag (5% flux)
;   
; OUTPUTS:
;   obs				- Structure containing all data for SED modeling. Subject to change!
;                     Common contents include: Object ID, position+error, mag+error, 
;                     flux+error, good photometry flag (bin), redshift+error, string
;                     containing all redshift+error from multiple sources (zall+zall_err).
;   band			- Array of central wavelength for photometric data. 
;   
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;   WISE reject mask directory and mask name must be changed for personal use.
;   
;   At current, sources need to pass the following criteria:
;       1. Have redshift information
;       2. Photometry S/N≥3 (≥1 for forced photometry)
;       3. Photometry must be "clean" (SDSS clean=1)
;       4. At least 7 photometric bands to (hopefully) avoid overfitting
;   
;   If the DERED keyword is set, the path to dust map directory must be set prior 
;   to running script.
;   
; EXAMPLES:
;
; REVISION HISTORY:
;   2015-Jul-29  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO read_sed_phot, file, $
	               MASK = mask, $
                   DERED = dered, $
                   FORCED_PHOT = forced_phot, $
                   MIN_ERR = min_err

nfiles = n_elements(file)
if (nfiles gt 1) then begin
    outfile = strarr(nfiles)
    for i = 0,n_elements(file)-1 do begin
        temp = strsplit(file[i],'/.',/extract)
        outfile[i] = temp[-2]+'_flux.sav';temp[where(strmatch(temp,'*part*'),outlen)]+'_flux.sav'
        outlen=1
        if (outlen eq 0) then stop
    endfor
endif else begin
    outfile = (strsplit(file,'/.',/extract))[-2]+'_flux.sav'
endelse

if keyword_set(mask) then begin
    ;; load reject mask
    mask_dir = '~/IDLWorkspace/libraries/cmc/mangle_masks/'
    wise_mask = mask_dir+'wise_mask_allwise_stars_pix.ply'
    read_mangle_polygons,wise_mask,allwise
endif

for f = 0,nfiles-1 do begin
    ;; read data
    data = mrdfits(file[f],1)
    
    ;; SDSS and XDQSOz indices
    iisdss = data.ra_sdss ne -9999. and data.dec_sdss ne -9999.
    isdss = where(iisdss,sdsslen)
    ixdqso = where(~iisdss,xdqsolen,/null)
    
    ;; number of catalog sources
    ndata = n_elements(data)
    ;; source data vectors
    objid = lon64arr(ndata)
    ra = dblarr(ndata) & e_ra = dblarr(ndata)
    dec = dblarr(ndata) & e_dec = dblarr(ndata)
    ;; source data fill
    objid[isdss] = data[isdss].objid
    objid[ixdqso] = data[ixdqso].objid_xdqso
    ra[*] = data.ra
    dec[*] = data.dec
    ;; positional errors in arcsec
    e_ra[isdss] = data[isdss].raerr
    e_dec[isdss] = data[isdss].decerr
    e_ra[ixdqso] = data[ixdqso].sigra
    e_dec[ixdqso] = data[ixdqso].sigdec

    ;; photometry variables (mag, e_mag, flux, e_flux)
    band = ['SDSS1','SDSS2','SDSS3','SDSS4','SDSS5', $
            'WISE1','WISE2','WISE3','WISE4', $
            'UK1','UK2','UK3','UK4', $
            'TWOM1','TWOM2','TWOM3', $
            'GALEX1','GALEX2']
    mag_vars = ['DERED_U','DERED_G','DERED_R','DERED_I','DERED_Z', $
                'W1','W2','W3','W4', $
                'YPETROMAG','J_1PETROMAG','HPETROMAG','KPETROMAG', $
                'J_M_STDAP','H_M_STDAP','K_M_STDAP', $
                'FUV_MAG','NUV_MAG']
    e_mag_vars = ['MODELMAGERR_u','MODELMAGERR_g','MODELMAGERR_r','MODELMAGERR_i','MODELMAGERR_z', $
                  'W1ERR','W2ERR','W3ERR','W4ERR', $
                  'YPETROMAGERR','J_1PETROMAGERR','HPETROMAGERR','KPETROMAGERR', $
                  'J_MSIG_STDAP','H_MSIG_STDAP','K_MSIG_STDAP', $
                  'FUV_MAGERR','NUV_MAGERR']
    flux_vars = mag_vars+'_FLUX'
    dered_vars = flux_vars+'_UNRED'
    e_flux_vars = e_mag_vars+'_FLUX'
    nbands = n_elements(band)

    ;; calculate flux & errors in microjanskys
    for i = 0,nbands-1 do begin
        re = execute(mag_vars[i]+'=data.'+mag_vars[i])
        re = execute(e_mag_vars[i]+'=data.'+e_mag_vars[i])
        re = execute(flux_vars[i]+'=magflux('+mag_vars[i]+','+e_mag_vars[i]+',band[i],err='+e_flux_vars[i]+')')
    endfor
	
	;; correct for Galactic extinction
    if keyword_set(dered) then $
        for i = 0,nbands-1 do $
            ;; if using SDSS photometry that has not been deredenned, comment out the line below
            if (strmatch(band[i],'SDSS*') eq 0) then $
                re = execute(flux_vars[i]+'=mw_ext_corr(ra,dec,'+flux_vars[i]+',band[i])')
	
    ;; conversions needed for nanomaggies & Vega2AB flux conversions
    ;; XDQSOz and unWISE in nanomaggies
    ;; XDQSOz and WISE photometry zero pt also in Vega
    nmgy2mujy = 3.631												;; 1 nMgy = 3.631 uJy
    v2ab_flux = 10.^([-2.699d,-3.339d,-5.174d,-6.620d]/2.5)			;; Vega2AB flux

    if (xdqsolen gt 0.) then begin
        ;; add XDQSO flux and calculate magnitudes
        ;;             SDSS     GALEX   UKIDSS   WISE
        ;; psfflux==[u,g,r,i,z,NUV,FUV,Y,J,H,Ks,W1,W2]
        psfband = ['SDSS1','SDSS2','SDSS3','SDSS4','SDSS5','GALEX2','GALEX1','UK1','UK2','UK3','UK4','WISE1','WISE2']
        match2,band,psfband,ix2f,if2x									;; EXTREMELY IMPORTANT!!!
        xdb = psfband[ix2f]
        xdf = nmgy2mujy * data[ixdqso].psfflux[ix2f]					;; here there is no WISE3/WISE4 match... so these are
        xde = nmgy2mujy * 1./sqrt(data[ixdqso].psfflux_ivar[ix2f])		;; replaced with SDSS1 (where psfband[ixd2f=-1]).
        																;; note: not indexing back of array w/ -1 (would be WISE2)
        ;; remove non-finite and null values														
        ifin = where(~finite(xdf) or ~finite(xde) or xdf eq 0. or xde eq 0.,finct)
        if (finct gt 0.) then begin
            xdf[ifin] = 0.
            xde[ifin] = 0.
        endif
        ;; remove observations where flux+err is not greater than background (flux == 0).
        ;; observations where flux+err > 0, set flux = (flux+error) / 2
        ineg = where(xdf lt 0.,negct)
        if (negct gt 0.) then begin
            pos_val = ((xdf[ineg]+xde[ineg])>0.)/2.
            xdf[ineg] = pos_val
        endif
        ;; convert WISE to AB 
        iwise = where(strmatch(xdb,'WISE*'),ct)
        if (ct eq 0) then stop
        v2ab = rebin(v2ab_flux[0:1],2,n_elements(xdf[0,*]))
        xdf[iwise,*] *= v2ab
        xde[iwise,*] *= v2ab
        
        ;; add XDQSOz photometry to full photometry
        for i = 0,n_elements(band)-1 do begin
            if (ix2f[i] eq -1) then continue							;; skip the unmatched WISE3/WISE4
            re = execute(flux_vars[i]+'[ixdqso] = xdf[i,*]')
            re = execute(e_flux_vars[i]+'[ixdqso] = xde[i,*]')
            ;; compute XDQSOz mags and errors
            re = execute(mag_vars[i]+'[ixdqso] = magflux(xdf[i,*],xde[i,*],xdb[i],err='+e_mag_vars[i]+'[ixdqso],/flux_in)')
        endfor
    endif
    
    ;; add unWISE flux and calculate magnitudes
    if keyword_set(forced_phot) then begin
        unwb = band[where(strmatch(band,'WISE*'),ct)]			;; unWISE band
        if (ct eq 0) then stop
        unwf = ['w1','w2','w3','w4']+'_nanomaggies'					;; unWISE flux
        e_unwf = unwf + '_ivar'										;; unWISE inverse variance
        unw = unwf+'_mag'											;; unWISE mags
        e_unw = e_unwf+'_mag'										;; unWISE mag errors    
    
        for i = 0,n_elements(unwb)-1 do begin
            ;; convert from nanomaggies to microjansky, from Vega to AB
            re = execute(unwf[i]+'= nmgy2mujy * v2ab_flux[i] * data.'+unwf[i])
            re = execute(e_unwf[i]+'= nmgy2mujy * v2ab_flux[i] * 1./sqrt(data.'+e_unwf[i]+')')
            ;; remove -9999 detections
            re = execute('inin = where(data.'+unwf[i]+' eq -9999. or data.'+e_unwf[i]+' eq -9999.,ninct)')
            if (ninct ne 0) then re = execute(unwf[i]+'[inin] = 0. & '+e_unwf[i]+'[inin] = 0.')
            ;; remove non-finite values
            re = execute('ifin = where(~finite('+unwf[i]+') or ~finite('+e_unwf[i]+'),finct)')
            if (finct ne 0) then re = execute(unwf[i]+'[ifin] = 0. & '+e_unwf[i]+'[ifin] = 0.')
            ;; indices of negative fluxes
            re = execute('ineg = where('+unwf[i]+' lt 0.,negct)')
            if (negct ne 0) then begin
                ;; where flux+error > 0, set flux = (flux+error) / 2 
                re = execute('pos_val = (('+unwf[i]+'[ineg] + '+e_unwf[i]+'[ineg])>0.)/2.')
                re = execute(unwf[i]+'[ineg] = pos_val')
            endif
            ;; compute unWISE mags and errors
            re = execute(unw[i]+' = magflux('+unwf[i]+','+e_unwf[i]+',unwb[i],err='+e_unw[i]+',/flux_in)')
        endfor		
    
        ;; add unWISE photometry to full photometry
        for i = 0,n_elements(unwb)-1 do begin
            iw = where(strmatch(band,unwb[i]),wct)				;; match unWISE band
            if (wct eq 0) then stop
            re = execute('sn_wise = '+flux_vars[iw[0]]+'/'+e_flux_vars[iw[0]])
            re = execute('sn_unw = '+unwf[i]+'/'+e_unwf[i])
            irep = where(finite(sn_unw) and sn_unw gt sn_wise)
            re = execute(mag_vars[iw[0]]+'[irep] = '+unw[i]+'[irep]')
            re = execute(e_mag_vars[iw[0]]+'[irep] = '+e_unw[i]+'[irep]')
            re = execute(flux_vars[iw[0]]+'[irep] = '+unwf[i]+'[irep]')
            re = execute(e_flux_vars[iw[0]]+'[irep] = '+e_unwf[i]+'[irep]')
        endfor
    endif
    
    ;; combine photometry
    mag = dblarr(nbands,ndata)
    e_mag = dblarr(nbands,ndata)
    flux = dblarr(nbands,ndata)
    e_flux = dblarr(nbands,ndata)
    for i = 0,nbands-1 do begin
        re = execute('mag[i,*]='+mag_vars[i])
        re = execute('e_mag[i,*]='+e_mag_vars[i])
        re = execute('flux[i,*]='+flux_vars[i])
        re = execute('e_flux[i,*]='+e_flux_vars[i])
    endfor
    ;; minimum photometric errors of ±0.05 mag (5% flux)
    if keyword_set(min_err) then e_flux = e_flux > sqrt((-0.4*alog(10)*flux*0.05)^2)
    
    ;; remove observations that fail S/N cut
    sn = flux/e_flux
    isn = where(finite(sn) and sn ge 3.,complement=ibadsn,ncomplement=nbad)
    if (nbad gt 0) then begin
        flux[ibadsn] = 0.
        e_flux[ibadsn] = 0.
        mag[ibadsn] = -9999.
        e_mag[ibadsn] = -9999.
    endif
    ;; byte index for good photometry
    bin = flux gt 0. and e_flux gt 0.

    ;; remove 2MASS where UKIDSS is available
    i2m = where(strmatch(band,'TWOM?'),n2m)
    iuk = where(strmatch(band,'UK?'),nuk)
    if (n2m eq 3 and nuk eq 4) then begin
        iuk = iuk[1:-1]
        jhk_uk = bin[iuk,*]
        jhk_2m = bin[i2m,*]
        irem = where(total(jhk_uk,1) gt 0,nrem)
        if (nrem gt 0) then jhk_2m[*,irem] = 0
        bin[i2m,*] = jhk_2m
    endif
    
    ;; SOURCE DATA FILL
    ;; output data structure	
    obs = {objid: long64(0), $
           ra: 0d, $
           e_ra: 0d, $
           dec: 0d, $
           e_dec: 0d, $
           mag: dblarr(nbands), $
           e_mag: dblarr(nbands), $
           flux: dblarr(nbands), $
           e_flux: dblarr(nbands), $
           bin: bytarr(nbands), $
           z: 0d, $
           zerr: 0d, $
           ztype:'', $
           zall: '', $
           zallerr: '', $
           class: '' $
           }
    obs = replicate(obs,ndata)
    
    ;; object ID and positions
    obs[isdss].objid = data[isdss].objid
    if (xdqsolen gt 0.) then obs[ixdqso].objid = long64(strtrim(data[ixdqso].objid_xdqso,2))
    obs.ra = ra
    obs.dec = dec
    obs.e_ra = e_ra
    obs.e_dec = e_dec
    ;; photometry data fill
    obs.mag = mag
    obs.e_mag = e_mag
    obs.flux = flux
    obs.e_flux = e_flux
    ;; good photometry flag
    obs.bin = bin
    
    ;; 2MASS photometry quality flags 'ABC' for reliable data
    if (n2m eq 3) then begin
        qual = strtrim(data.ph_qual_2mass,2)
        qual = transpose([[strmatch(strmid(qual,0,1),'[ABC]')],[strmatch(strmid(qual,1,1),'[ABC]')],[strmatch(strmid(qual,2,1),'[ABC]')]])
        bin2m = bin[i2m,*]
        bin2m[where(qual eq 0,/null)] = 0
        bin[i2m,*] = bin2m
    endif
    
    ;; full redshift data set
    ;; (1) ZP     == SDSS DR14 phot-z
    ;; (2) PEAKZ  == XDQSOz (DiPompeo+15)
    ;; (3) ZS     == SDSS DR14 spec-z
    ;; (4) Z_SUPP == Reyes+08, Lacy+13, Hainline+14, Yuan+16
    zstr = ['zp','peakz','zs','z_zsupp']
    e_zstr = ['zperr','peakfwhm','zserr','cat_zsupp']	;; on (4) returns source catalog
    zall = strarr(ndata)			;; all redshift data
    zallerr = strarr(ndata)			;; all redshift error data
    z = dblarr(ndata)				;; "best" redshift value
    zerr = dblarr(ndata)			;; "best" redshift error value
    ;; trust SDSS photometric redshifts only where reliable; -1 <= photoerrorclass <= 3 (photoerrorclass == 1 is best match)
    iphotz = where(data.photoerrorclass ge -1 and data.photoerrorclass lt 3,complement=badz,ncomplement=nbadz)
    if (nbadz gt 0) then data[badz].zp = -9999.
    ;; sort redshift data
    for i = 0,n_elements(zstr)-1 do begin
        re = execute('iz = where(finite(data.'+zstr[i]+') and data.'+zstr[i]+' gt 0.,zlen)')
        if (zlen gt 0.) then begin
            re = execute('zall[iz] += strtrim(data[iz].'+zstr[i]+',2)')
            re = execute('z[iz] = data[iz].'+zstr[i])
            re = execute('zallerr[iz] += strtrim(data[iz].'+e_zstr[i]+',2)')
            if (e_zstr[i] ne 'cat_zsupp') then re = execute('zerr[iz] = data[iz].'+e_zstr[i]) else $
                                               zerr[iz] = z[iz]*0.05
        endif
        zall += ','
        zallerr += ','
    endfor
    ztype = zorig(zall)
    
    ;; finalized source redshift information
    obs.z = z
    obs.zerr = zerr
    obs.ztype = ztype
    obs.zall = zall
    obs.zallerr = zallerr
    
    ;; SDSS recognized "CLASS"
    obs.class = data.class
    
    ;; QUALITY CUTS (goodbye, sweet sources)
    ;; keep sources that...
    
    ;; ..have clean photometry
    ;; NOTE: After this step OBS and DATA are no longer the same length!
    ikeep = where(data.clean eq 1,ct)
    if (ct eq 0) then continue
    obs = obs[ikeep]
    		
    ;; ...have redshift data!; loss of photo-z from photoerrorclass, etc. (sanity check)
    ikeep = where(obs.z gt 0.,ct)
    if (ct eq 0) then continue
    obs = obs[ikeep]
    
    ;; ...have redshifts within specified range (phot-z ≤ 0.6; spec-z ≤ 1.0)
    iizs = strmatch(obs.ztype,'ZS*') and obs.z le 1.0
    iizp = strmatch(obs.ztype,'ZP') and obs.z le 0.6
    ikeep = where(iizs or iizp,ct)
    if (ct gt 0) then obs = obs[ikeep]
    
    ;; ...have constrained redshift errors
    ikeep = where(obs.zerr gt 0.,ct)
    if (ct eq 0) then continue
    obs = obs[ikeep]
    
    ;; ...can be resampled in redshift space (redshift cuttoff z == 3)
    ikeep = where(obs.z+4.*obs.zerr lt 3.,ct)
    if (ct eq 0.) then continue
    obs = obs[ikeep]
    
    ;; ...have a minimum number of 7 photometric bands
    ikeep = where(total(obs.bin,1) ge 7,ct)
    if (ct eq 0) then continue 
    obs = obs[ikeep]
    
    ;; ...have detections in all four WISE bands
    iwise = where(strmatch(band,'WISE?'),ct)
    if (ct eq 0) then continue
    ikeep = where(total(obs.bin[iwise],1) eq 4,ct)
    if (ct gt 0) then obs = obs[ikeep]
    
    ;; ...aren't in the mask
    if keyword_set(mask) then begin
        euler,obs.ra,obs.dec,gal_l,gal_b,1
        in_allwise=is_in_window_pix(ra=gal_l,dec=gal_b,allwise,scheme='6s')
        ikeep = where(~in_allwise,ct)				;; sources not in reject mask
        if (ct eq 0.) then continue
        obs = obs[ikeep]
    endif
    
    ;; finally, no duplicate objects in data set!
    if (n_elements(uniq(obs.objid,sort(obs.objid))) ne n_elements(obs)) then begin
        print, 'DUPLICATE SOURCE DETECTED!'
        obsid = obs.objid
        totbin = total(obs.bin,1)
        ;; remove duplicate objects, keep sources with most photometry
        ikeep = rem_dup(obsid,totbin)
        obs = obs[ikeep]
    endif
    
    save,obs,band,/compress,file=outfile[f]
endfor


END



