;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   read_sed_phot
;   
; PURPOSE:
;   Read photometry from .FITS file and format for SED modeling.
;   
; CALLING SEQUENCE:
;   read_sed_phot, file, [, /MASK, /NIR, /FORCED_SN ]
;
; INPUTS:
;   file			- String containing the name of the photometry data file.
;	
; OPTIONAL INPUTS:
;   /MASK			- Use mangle mask to remove sources in WISE near bright stars.
;   /NIR			- Only keep sources with at least 1 NIR detection.
;	/FORCED_SN		- Require all WISE forced photometry to have S/N≥1.
;	
; OUTPUTS:
;	obs				- Structure containing all data for SED modeling. Subject to change!
;                     Common contents include: Object ID, RA, Dec, mag+error, flux+error,
;                     good photometry flag (bin), redshift+error, string containing all
;                     redshift+error from multiple sources (zarr+e_zarr).
;	band			- Array of central wavelength for photometric data. 
;   
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	WISE reject mask directory and mask name must be changed for personal use.
;
;	At current, sources need to pass the following criteria:
;		1. Have redshift information
;		2. Photometry S/N≥3 (≥1 for forced photometry)
;		3. Photometry must be "clean" (SDSS clean=1)
;		4. At least 7 photometric bands to (hopefully) avoid overfitting
;   
; EXAMPLES:
;
; REVISION HISTORY:
;   2015-Jul-29  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO read_sed_phot, file, $
	               MASK = mask, $
	               NIR = nir, $
	               FORCED_SN = forced_sn


outfile = file
for i = 0,n_elements(file)-1 do begin
	temp = strsplit(file[i],'/.',/extract)
	outfile[i] = temp[-3]+'_flux.sav'
endfor

if keyword_set(mask) then begin
	;; load reject mask
	mask_dir = '~/IDLWorkspace/libraries/cmc/mangle_masks/'
	wise_mask = mask_dir+'wise_mask_allwise_stars_pix.ply'
	read_mangle_polygons,wise_mask,allwise
endif

for f = 0,n_elements(file)-1 do begin
	;; read data
	data = mrdfits(file[f],1)
	
	;; restrict to sources with NIR photometry 
	if keyword_set(nir) then begin
		inir = where(data.ra_ukidss gt -9999. and data.dec_ukidss gt -9999.,nirlen)
		if (nirlen gt 0) then data = data[inir] else continue
	endif

	;; keep only sources with unWISE data
	iunwise = where(data.ra_unwise ne -9999. and data.dec_unwise ne -9999.,unwiselen)
	if (unwiselen eq 0) then continue
	data = data[iunwise]
	;; SDSS and XDQSOz indices
	iisdss = data.ra_sdss ne -9999. and data.dec_sdss ne -9999.
	isdss = where(iisdss,sdsslen)
	ixdqso = where(~iisdss,xdqsolen)
    
    ;; photometry variables (mag, e_mag, flux, e_flux)
	filt = ['SDSS1','SDSS2','SDSS3','SDSS4','SDSS5', $
			'WISE1','WISE2','WISE3','WISE4', $
			'UK1','UK2','UK3','UK4' $
			]
	mag_vars = ['DERED_U','DERED_G','DERED_R','DERED_I','DERED_Z', $
				'W1','W2','W3','W4', $
				'YPETROMAG','J_1PETROMAG','HPETROMAG','KPETROMAG' $
				]
	e_mag_vars = ['MODELMAGERR_u','MODELMAGERR_g','MODELMAGERR_r','MODELMAGERR_i','MODELMAGERR_z', $
				  'W1ERR','W2ERR','W3ERR','W4ERR', $
				  'YPETROMAGERR','J_1PETROMAGERR','HPETROMAGERR','KPETROMAGERR' $
				  ]
	flux_vars = mag_vars+'_FLUX'
	e_flux_vars = e_mag_vars+'_FLUX'
	nfilts = n_elements(filt)

	;; calculate flux & errors in microjanskys
	for i = 0,nfilts-1 do begin
		re = execute(mag_vars[i]+'=data.'+mag_vars[i])
		re = execute(e_mag_vars[i]+'=data.'+e_mag_vars[i])
		re = execute(flux_vars[i]+'=magflux('+mag_vars[i]+','+e_mag_vars[i]+',filt[i],err='+e_flux_vars[i]+')')
	endfor
	
	;; conversions needed for nanomaggies & Vega2AB flux conversions
	;; XDQSOz and unWISE in nanomaggies
	;; XDQSOz and WISE photometry zero pt also in Vega
	nmgy2mujy = 3.631												;; 1 nMgy = 3.631x10^-6 Jy
	v2ab_flux = 10.^([-2.699d,-3.339d,-5.174d,-6.620d]/2.5)			;; Vega2AB flux

	if (xdqsolen gt 0.) then begin
		;; add XDQSO flux and calculate magnitudes
		;;             SDSS     GALEX   UKIDSS   WISE
		;; psfflux==[u,g,r,i,z,NUV,FUV,Y,J,H,Ks,W1,W2]
		psfband = ['SDSS1','SDSS2','SDSS3','SDSS4','SDSS5','GALEX2','GALEX1','UK1','UK2','UK3','UK4','WISE1','WISE2']
		match2,filt,psfband,ix2f,if2x									;; EXTREMELY IMPORTANT!!!
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
		iwise = where(strmatch(xdb,'WISE*'))
		v2ab = rebin(v2ab_flux[0:1],2,n_elements(xdf[0,*]))
		xdf[iwise,*] *= v2ab
		xde[iwise,*] *= v2ab

		;; add XDQSOz photometry to full photometry
		for i = 0,n_elements(filt)-1 do begin
			if (ix2f[i] eq -1) then continue							;; skip the unmatched WISE3/WISE4
			re = execute(flux_vars[i]+'[ixdqso] = xdf[i,*]')
			re = execute(e_flux_vars[i]+'[ixdqso] = xde[i,*]')
			;; compute XDQSOz mags and errors
			re = execute(mag_vars[i]+'[ixdqso] = magflux(xdf[i,*],xde[i,*],xdb[i],err='+e_mag_vars[i]+'[ixdqso],/flux_in)')
		endfor
	endif

	;; add unWISE flux and calculate magnitudes
	unwb = filt[where(strmatch(filt,'WISE*'))]					;; unWISE band
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
		iw = where(strmatch(filt,unwb[i]),wct)				;; match unWISE band
		if (wct eq 0) then stop
		re = execute(mag_vars[iw[0]]+' = '+unw[i])
		re = execute(e_mag_vars[iw[0]]+' = '+e_unw[i])
		re = execute(flux_vars[iw[0]]+' = '+unwf[i])
		re = execute(e_flux_vars[iw[0]]+' = '+e_unwf[i])
	endfor
    
	;; output data structure	
    ndata = n_elements(data)
	obs = {objid: long64(0), $
		   ra: 0d, $
		   dec: 0d, $
		   mag: dblarr(nfilts), $
		   e_mag: dblarr(nfilts), $
		   flux: dblarr(nfilts), $
		   e_flux: dblarr(nfilts), $
		   ;raw_mag: dblarr(nfilts), $
		   ;raw_e_mag: dblarr(nfilts), $
		   ;raw_flux: dblarr(nfilts), $
		   ;raw_e_flux: dblarr(nfilts), $		   
		   bin: bytarr(nfilts), $
		   z: 0d, $
		   e_z: 0d, $
		   zarr: '', $
		   e_zarr: '' $
		   }
	obs = replicate(obs,ndata)
	
	;; source data fill
	obs[isdss].objid = data[isdss].objid
	obs[ixdqso].objid = long64(strtrim(data[ixdqso].objid_xdqso,2))
	obs.ra = data.ra
	obs.dec = data.dec
	
	for i = 0,n_elements(filt)-1 do begin
		re = execute('obs.mag[i]='+mag_vars[i])
		re = execute('obs.e_mag[i]='+e_mag_vars[i])
		re = execute('obs.flux[i]='+flux_vars[i])
		re = execute('obs.e_flux[i]='+e_flux_vars[i])
		;re = execute('obs.raw_e_flux[i]='+e_flux_vars[i])
	endfor
	;; minimum photometric errors of ±0.05 mag (5% flux)
	obs.e_flux = obs.e_flux > sqrt((-0.4d*(0.05*obs.flux*alog(10.)))^2)
	
	;; S/N > 3 for optical/NIR photometry
	inotir = where(~strmatch(filt,'WISE*'))
	flux = obs.flux[inotir]
	e_flux = obs.e_flux[inotir]
	mag = obs.mag[inotir]
	e_mag = obs.e_mag[inotir]
	sn = flux/e_flux
	isn = where(~finite(sn) or sn lt 3.,snct)
	;; remove observations that fail S/N cut
	if (snct ne 0) then begin
		flux[isn] = 0.
		e_flux[isn] = 0.
		mag[isn] = -9999.
		e_mag[isn] = -9999.
	endif
	obs.flux[inotir] = flux
	obs.e_flux[inotir] = e_flux
	obs.mag[inotir] = mag
	obs.e_mag[inotir] = e_mag

	;; S/N > 1 for forced photometry, otherwise trust it
	if keyword_set(forced_sn) then begin
		iir = where(strmatch(filt,'WISE*'))
		flux = obs.flux[iir]
		e_flux = obs.e_flux[iir]
		mag = obs.mag[iir]
		e_mag = obs.e_mag[iir]
		sn = flux/e_flux
		isn = where(~finite(sn) or sn lt 1.,snct)
		;; remove observations that fail S/N cut
		if (snct ne 0) then begin
			flux[isn] = 0.
			e_flux[isn] = 0.
			mag[isn] = -9999.
			e_mag[isn] = -9999.
		endif
		obs.flux[iir] = flux
		obs.e_flux[iir] = e_flux
		obs.mag[iir] = mag
		obs.e_mag[iir] = e_mag
	endif

	;; byte index for good photometry
	obs.bin = obs.flux gt 0. and obs.e_flux gt 0.				

	;; full redshift data set
	;; (1) ZP     == SDSS DR14 phot-z
	;; (2) PEAKZ  == XDQSOz (DiPompeo+15)
	;; (3) ZS     == SDSS DR14 spec-z
	;; (4) Z_SUPP == Reyes+08, Lacy+13, Hainline+14, Yuan+16
	zstr = ['zp','peakz','zs','z_zsupp']
	e_zstr = ['zperr','peakfwhm','zserr','cat_zsupp']	;; on (4) returns source catalog
	zarr = strarr(ndata)			;; all redshift data
	e_zarr = strarr(ndata)			;; all redshift error data
	z = dblarr(ndata)				;; "best" redshift value
	e_z = dblarr(ndata)				;; "best" redshift error value
	;; trust SDSS photometric redshifts only when photoerrorclass == 1
	iphotz = where(data.clean eq 1,complement=inphotz,ncomplement=nphotz)
	if (nphotz gt 0) then data[inphotz].zp = -9999.
	;; sort redshift data
	for i = 0,n_elements(zstr)-1 do begin
		re = execute('iz = where(finite(data.'+zstr[i]+') and data.'+zstr[i]+' gt 0.,zlen)')
		if (zlen gt 0.) then begin
			re = execute('zarr[iz] += strtrim(data[iz].'+zstr[i]+',2)')
			re = execute('z[iz] = data[iz].'+zstr[i])
			re = execute('e_zarr[iz] += strtrim(data[iz].'+e_zstr[i]+',2)')
			if (e_zstr[i] ne 'cat_zsupp') then re = execute('e_z[iz] = data[iz].'+e_zstr[i]) else $
			                                   e_z[iz] = -9999.
		endif
		zarr += ','
		e_zarr += ','
	endfor
	
	;; save all redshift information
	obs.z = z
	obs.e_z = e_z
	obs.zarr = zarr
	obs.e_zarr = e_zarr
	
	;; QUALITY CUTS (goodbye, sweet sources)
	;; keep sources that...

	;; ..have clean photometry
	;; NOTE: After this step OBS and DATA are no longer the same length!
	ikeep = where(data.clean,ct)
	if (ct eq 0) then continue
	obs = obs[ikeep]
			
	;; ...have redshift data!; loss of photo-z from photoerrorclass, etc. (sanity check)
	ikeep = where(obs.z gt 0.,ct)
	if (ct eq 0) then continue
	obs = obs[ikeep]
	
	;; ...have a minimum number of 7 photometric bands
	ikeep = where(total(obs.bin,1) ge 7,ct)
	if (ct eq 0) then continue 
	obs = obs[ikeep]
		
	;; ...aren't in the mask
	if keyword_set(mask) then begin
		euler,obs.ra,obs.dec,gal_l,gal_b,1
		in_allwise=is_in_window_pix(ra=gal_l,dec=gal_b,allwise,scheme='6s')
		ikeep = where(~in_allwise,ct)				;; sources not in reject mask
		if (ct eq 0.) then continue
		obs = obs[ikeep]
	endif

	;; ...have NIR photometry; loss from S/N
	if keyword_set(nir) then begin
		inir = where(strmatch(filt,'UK*'),nirlen)
		if (nirlen eq 0) then stop
		ikeep = where(total(obs.bin[inir],1) ge 1,ct)
		if (ct eq 0) then stop
		obs = obs[ikeep]
	endif

	;; finally, no duplicate objects in data set!
	if (n_elements(uniq(obs.objid,sort(obs.objid))) ne n_elements(obs)) then begin
		print, 'DUPLICATE SOURCE DETECTED!'
		stop
	endif

	band = filt
	save,obs,band,/compress,file=outfile[f]
endfor


END







