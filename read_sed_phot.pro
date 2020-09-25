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
;   /MIN_ERR        - Require a minimum error of ±0.05 mag (5% flux).
;   /ACCEPT         - If set, keep only sources which pass all quality cuts.
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
                   MIN_ERR = min_err, $
                   ACCEPT = accept


;on_error,2
;npar = n_params()

;; load mask only once
if keyword_set(mask) then begin
    ;; load reject mask
    mask_dir = '~/IDLWorkspace/libraries/cmc/mangle_masks/'
    wise_mask = mask_dir+'wise_mask_allwise_stars_pix.ply'
    read_mangle_polygons,wise_mask,allwise
endif

nfiles = n_elements(file)
outfile = strarr(nfiles)
for i = 0,nfiles-1 do outfile[i] = strsplit((strsplit(file[i],'/',/extract))[-1],'fits',/regex,/extract)+'sav'

for f = 0,nfiles-1 do begin
    ;; read data
    data = mrdfits(file[f],1)
    
    ;; instrument indices
    iisdss = data.ra_sdss ne -9999. and data.dec_sdss ne -9999.
    isdss = where(iisdss,sdsslen)
    iixdqso = data.ra_xdqso ne -9999. and data.dec_xdqso ne -9999.
    ixdqso = where(iixdqso and ~iisdss,xdqsolen,/null)
    iizsupp = data.ra_zsupp ne -9999. and data.dec_zsupp ne -9999.
    iiwise = data.ra_wise ne -9999. and data.dec_wise ne -9999.
    iiunwise = data.ra_unwise ne -9999. and data.dec_unwise ne -9999.
    iiukidss = data.ra_ukidss ne -9999. and data.dec_ukidss ne -9999.
    iitwom = data.ra_2mass ne -9999. and data.dec_2mass ne -9999.
    iigalex = data.ra_galex ne -9999. and data.dec_galex ne -9999.
    
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
        match2,band,psfband,iix2f,if2x									;; EXTREMELY IMPORTANT!!!
        ix2f = iix2f[where(iix2f ne -1,/null)]
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
        if (negct gt 0.) then xdf[ineg] = ((xdf[ineg]+xde[ineg])>0.)/2.
        ;; convert WISE to AB 
        iw12 = where(strmatch(xdb,'WISE*'),ct)
        if (ct eq 0) then stop
        v2ab = rebin(v2ab_flux[0:1],2,n_elements(xdf[0,*]))
        xdf[iw12,*] *= v2ab
        xde[iw12,*] *= v2ab
        
        ;; add XDQSOz photometry to full photometry
        for i = 0,n_elements(band)-1 do begin
            ;; check for band match
            im = (where(strmatch(xdb,band[i]),ct))[0]
            if (ct eq 0) then continue
            ;; check for valid photometry
            ipos = where(xdf[im,*] gt 0. and xde[im,*] gt 0.,poslen)
            if (poslen eq 0) then continue
            re = execute(flux_vars[i]+'[ixdqso[ipos]] = xdf[im,ipos]')
            re = execute(e_flux_vars[i]+'[ixdqso[ipos]] = xde[im,ipos]')
            re = execute(mag_vars[i]+'[ixdqso[ipos]] = magflux(xdf[im,ipos],xde[im,ipos],xdb[im],err=xdmerr,/flux_in)')
            re = execute(e_mag_vars[i]+'[ixdqso[ipos]] = xdmerr')
        endfor
    endif

    ;; add unWISE flux and calculate magnitudes
    if keyword_set(forced_phot) then begin
        iirepir = bytarr(4,ndata)
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
            ;; where flux+error > 0, set flux = (flux+error) / 2 
            re = execute('ineg = where('+unwf[i]+' lt 0.,negct)')
            if (negct ne 0) then re = execute(unwf[i]+'[ineg] = (('+unwf[i]+'[ineg] + '+e_unwf[i]+'[ineg])>0.)/2.')
            ;; compute unWISE mags and errors
            ;; keyword error out works here as it is creating the array
            re = execute(unw[i]+' = magflux('+unwf[i]+','+e_unwf[i]+',unwb[i],err='+e_unw[i]+',/flux_in)')
        endfor		
    
        ;; add unWISE photometry to full photometry
        for i = 0,n_elements(unwb)-1 do begin
            iw = (where(strmatch(band,unwb[i]),wct))[0]				;; match unWISE band
            if (wct eq 0) then stop
            ;; use best S/N
            re = execute('sn_sampl = '+flux_vars[iw]+'/'+e_flux_vars[iw])
            re = execute('sn_unw = '+unwf[i]+'/'+e_unwf[i])
            ;; first pass, S/N unWISE greater than sample (disregards non-finite values)
            iigt = sn_unw gt sn_sampl
            ;; replace non-finite sample with real-valued unWISE
            iifi = finite(sn_unw) and ~finite(sn_sampl)
            ;; combine
            irep = where(iigt or iifi,repct,complement=iww)
            if (repct eq 0) then continue
            iirepir[i,irep] = 1
            re = execute(mag_vars[iw]+'[irep] = '+unw[i]+'[irep]')
            re = execute(e_mag_vars[iw]+'[irep] = '+e_unw[i]+'[irep]')
            re = execute(flux_vars[iw]+'[irep] = '+unwf[i]+'[irep]')
            re = execute(e_flux_vars[iw]+'[irep] = '+e_unwf[i]+'[irep]')
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
    
    ;; remove 2MASS band where UKIDSS band is available
    if (total(where([strmatch(band,'TWOM?'),strmatch(band,'UK?')])) gt 1) then begin
        jhk =['1','2','3','4']
        for i = 0,2 do begin
            ijhk = where([strmatch(band,'TWOM'+jhk[i]),strmatch(band,'UK'+jhk[i+1])],njhk)
            if (njhk eq 2) then bin[ijhk[0],where(bin[ijhk[1]-nbands,*],/null)] = 0
        endfor
    endif
    
    ;; 2MASS photometry quality flags 'ABC' for reliable data
    i2m = where(strmatch(band,'TWOM?'),n2m)
    if (n2m gt 0) then begin
        qual = strtrim(data.ph_qual_2mass,2)
        iiqual = transpose([[strmatch(strmid(qual,0,1),'[ABC]')],[strmatch(strmid(qual,1,1),'[ABC]')],[strmatch(strmid(qual,2,1),'[ABC]')]])
        for i = 0,2 do begin
            ijhk = where(strmatch(band,'TWOM'+jhk[i]),njhk)
            if (njhk eq 1) then bin[ijhk[0],where(iiqual[i,*] eq 0,/null)] = 0
        endfor
    endif
    ;; 2MASS indices
    ;i2mfilt = where(strmatch(band,'TWOM?'))
    ;bin2m = bin[i2mfilt,*]
    ;iitwom = total(bin2m,1) gt 0
    
    ;; full redshift data set in ascending order of use
    ;; ZP     == SDSS DR14 phot-z
    ;; PEAKZ  == XDQSOz (DiPompeo+15)
    ;; ZS     == SDSS DR14 spec-z
    ;; Z_SUPP == Reyes+08, Lacy+13, Hainline+14, Yuan+16
    zdat = ['ZP','PEAKZ','ZS','ZSUPP']                ;; name of redshift tag
    e_zdat = ['ZPERR','PEAKFWHM','ZSERR','ZSUPP']	;; name of redshift error tag; N/A for Z_SUPP so placeholder
    zstr = strarr(ndata)			;; all redshift data
    zerrstr = strarr(ndata)			;; all redshift error data
    zbin = strarr(ndata)            ;; valid redshift flag
    z = dblarr(ndata)				;; "best" redshift value
    zerr = dblarr(ndata)			;; "best" redshift error value
    ;; trust SDSS photometric redshifts only where reliable; -1 <= photoerrorclass <= 3 (photoerrorclass == 1 is best match)
    iigdzp = finite(data.photoerrorclass) and data.photoerrorclass ne -9999
    ;if (nbadz gt 0) then data[badz].zp = -9999.
    ;; sort redshift data
    for i = 0,n_elements(zdat)-1 do begin
        re = execute('z_value = data.'+zdat[i]+' & z_error = data.'+e_zdat[i])    
        if (zdat[i] eq 'ZSUPP') then begin
            isupp = where(strtrim(z_error,2),nsupp,complement=invalid)
            if (nsupp gt 0) then z_error[isupp] = z_value[isupp]*0.05
            z_error[invalid] = 0.       
        endif
        iiz = finite(z_value) and finite(z_error) and z_value gt 0. and z_error gt 0.
        iz = where(finite(z_value) and finite(z_error) and z_value gt 0. and z_error gt 0.,zlen)
        if (zlen gt 0.) then begin
            zstr[iz] += strtrim(z_value[iz],2)
            zerrstr[iz] += strtrim(z_error[iz],2)
        endif
        ;; store photometric redshift but do not flag it for use
        if (zdat[i] eq 'ZP') then iz = where(iiz and iigdzp,/null)        
        zbin[iz] += '1'
        ;; fill string for next iteration
        zstr += ','
        zerrstr += ','
        zbin += ','
    endfor
    ;; determine best redshift
    ztype = zorigin(zbin,zstr,zpref=zpref)
    !null = zorigin(zbin,zerrstr,zpref=zerrpref)
    
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
           repir: bytarr(4), $
           z: 0d, $
           zerr: 0d, $
           zstr: '', $
           zerrstr: '', $
           zbin: '', $
           ztype: '', $
           photoerrorclass: -9999, $
           class: '', $
           iiaccept: 1b $
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
    ;; IR replacement flag
    if keyword_set(forced_phot) then obs.repir = iirepir
    ;; finalized source redshift information
    obs.z = zpref
    obs.zerr = zerrpref
    obs.zstr = zstr
    obs.zerrstr = zerrstr
    obs.zbin = zbin
    obs.ztype = ztype
    obs.photoerrorclass = data.photoerrorclass > (-9999)
    ;; SDSS recognized "CLASS"
    obs.class = data.class
    
    ;; QUALITY CUTS (goodbye, sweet sources)
    ;; keep sources that...
    
    ;; NOTE: If ACCEPT keyword set, after this step OBS and DATA are no longer the same length!
    
    ;; ...have redshifts within specified range (0 < phot-z ≤ 0.6; 0 < spec-z ≤ 1.0)
    iizs = obs.z gt 0. and strmatch(obs.ztype,'ZS*') and obs.z le 1.0
    iizp = obs.z gt 0. and (strmatch(obs.ztype,'ZP') or strmatch(obs.ztype,'PEAKZ')) and obs.z le 0.6
    iizrang = iizs or iizp
    iaccept = where(iizrang,ct,complement=irem)
    if (ct eq 0) then continue
    obs[irem].iiaccept = 0
    
    ;; ...have constrained redshift errors
    iizerr = obs.zerr ne 0
    iaccept = where(iizerr ne 0.,ct,complement=irem)
    if (ct eq 0) then continue
    obs[irem].iiaccept = 0
    iizgood = iizrang and iizerr

    ;; ..have clean photometry
    iiclean = data.clean eq 1
    if (xdqsolen gt 0) then iiclean[ixdqso] = 1               ;; XDQSO photometry already verified
    iaccept = where(iiclean eq 1,ct,complement=irem)
    if (ct eq 0) then continue
    obs[irem].iiaccept = 0

    ;; ...can be resampled in redshift space (redshift cuttoff z == 3)
    ;iiaccept = where(obs.z+4.*obs.zerr lt 3.,ct,complement=irem)
    ;if (ct eq 0.) then continue
    ;obs[irem].iiaccept = 0
    
    ;; ...have a minimum number of 7 photometric bands
    iibands = total(obs.bin,1) ge 7
    iaccept = where(iibands,ct,complement=irem)
    if (ct eq 0) then continue 
    obs[irem].iiaccept = 0
    
    ;; ...have detections in all four WISE bands
    iwise = where(strmatch(band,'WISE?'),ct)
    if (ct eq 0) then continue
    iifourw = total(obs.bin[iwise],1) eq 4
    iaccept = where(iifourw,ct,complement=irem)
    if (ct eq 0) then continue
    obs[irem].iiaccept = 0

    ;; ...aren't in the mask
    iinomsk = 0b
    if keyword_set(mask) then begin
        euler,obs.ra,obs.dec,gal_l,gal_b,1
        in_allwise=is_in_window_pix(ra=gal_l,dec=gal_b,allwise,scheme='6s')
        iinomsk = ~in_allwise
        iaccept = where(iinomsk,ct,complement=irem)		 ;; sources not in reject mask
        if (ct eq 0.) then continue
        obs[irem].iiaccept = 0
    endif
    
    ;; finally, no duplicate objects in data set!
    iinodup = bytarr(ndata)
    if (n_elements(uniq(obs.objid,sort(obs.objid))) ne n_elements(obs)) then begin
        print, 'DUPLICATE SOURCE DETECTED!'
        obsid = obs.objid
        totbin = total(obs.bin,1)
        ;; remove duplicate objects, keep sources with most photometry
        ikeep = rem_dup(obsid,totbin)
        iinodup[ikeep] = 1
        iaccept = where(iinodup,ct,complement=irem)     ;; remove duplicates
        if (ct eq 0.) then continue
        obs[irem].iiaccept = 0
    endif

    if keyword_set(accept) then obs = obs[where(obs.iiaccept,/null)]
    
    save,obs,band,iisdss,iixdqso,iizsupp,iiwise,iiunwise,iiukidss,iitwom,iigalex,iiclean,iizrang,iizerr,iizgood,iibands,iifourw,iinomsk,iinodup,/compress,file=outfile[f]
endfor


END



