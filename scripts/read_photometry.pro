;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   read_photometry
;   
; PURPOSE:
;   Read photometry from .FITS file and format for SED modeling.
;   
; CALLING SEQUENCE:
;   read_sed_phot, file, [, /DERED, /MIN_ERR, /ACCEPT ]
;
; INPUTS:
;   file			- String containing the name of the photometry data file.
;	
; OPTIONAL INPUTS:
;   /MASK           - Use mangle mask to remove sources in WISE near bright stars.
;   /DERED          - Correct fluxes for Galactic extinction.
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
;       2. Photometry S/N≥3
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
PRO read_photometry, file, $
                     MIN_ERR = min_err, $
                     ACCEPT = accept


;; number of files/output file names (for large datasets spanning multiple files)
nfiles = n_elements('data/'+file)
outfile = strarr(nfiles)
for i = 0,nfiles-1 do outfile[i] = strsplit((strsplit(file[i],'/',/extract))[-1],'.fits',/regex,/extract)+'_READY.sav'

;; for each data file
for f = 0,nfiles-1 do begin
    ;; read data
    data = mrdfits('data/'+file[f],1)
    
    ;; valid data indices
    idata = where(data.ra ne -9999. and data.dec ne -9999.,ndata)
    
    ;; source data vectors
    objid = lon64arr(ndata)
    ra = dblarr(ndata) & e_ra = dblarr(ndata)
    dec = dblarr(ndata) & e_dec = dblarr(ndata)
    ;; source data fill
    objid[idata] = data[idata].objid
    ra[idata] = data[idata].ra
    dec[idata] = data[idata].dec
    ;; positional errors in arcsec
    e_ra[idata] = data[idata].raerr
    e_dec[idata] = data[idata].decerr

    ;; photometry variables (mag, e_mag, flux, e_flux)
    band = ['SDSS1','SDSS2','SDSS3','SDSS4','SDSS5', $
            'WISE1','WISE2','WISE3','WISE4', $
            'UK1','UK2','UK3','UK4', $
            'TWOM1','TWOM2','TWOM3', $
            'GALEX1','GALEX2']
    mag_vars = ['CMODELMAG_U','CMODELMAG_G','CMODELMAG_R','CMODELMAG_I','CMODELMAG_Z', $
                'W1','W2','W3','W4', $
                'YPETROMAG','J_1PETROMAG','HPETROMAG','KPETROMAG', $
                'J_M_STDAP','H_M_STDAP','K_M_STDAP', $
                'FUV_MAG','NUV_MAG']
    e_mag_vars = ['CMODELMAGERR_U','CMODELMAGERR_G','CMODELMAGERR_R','CMODELMAGERR_I','CMODELMAGERR_Z', $
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
        
    ;; add redshift data
    z = dblarr(ndata)
    zerr = dblarr(ndata)
    ;; photometric redshift
    zp = data[idata].zp
    zperr = data[idata].zperr
    ;; spectroscopic redshift
    zs = data[idata].zs
    zserr = data[idata].zserr
    ;; combine redshifts, choose spectroscopic where available
    z[*] = zp
    z[where(zs gt 0.,/null)] = zs[where(zs gt 0.,/null)]
    zerr[*] = zperr
    zerr[where(zs gt 0.,/null)] = zserr[where(zs gt 0.,/null)]
    
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
           iiaccept: 1b}
    obs = replicate(obs,ndata)

    ;; object ID and positions
    obs.objid = objid
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
    ;; finalized source redshift information
    obs.z = z
    obs.zerr = zerr
    
    ;; QUALITY CUTS (goodbye, sweet sources)
    ;; keep sources that...
        
    ;; ...have constrained redshift errors
    iizerr = obs.zerr ne 0 and obs.zerr lt 0.5
    iaccept = where(iizerr ne 0.,ct,complement=irem,ncomplement=nrem)
    if (ct eq 0) then continue
    if (nrem gt 0) then obs[irem].iiaccept = 0

    ;; ...have a minimum number of 7 photometric bands
    iibands = total(obs.bin,1) ge 7
    iaccept = where(iibands,ct,complement=irem,ncomplement=nrem)
    if (ct eq 0) then continue 
    if (nrem gt 0) then obs[irem].iiaccept = 0
        
    if keyword_set(accept) then obs = obs[where(obs.iiaccept,/null)]
    
    save,obs,band,/compress,file='data/'+outfile[f]
endfor


END



