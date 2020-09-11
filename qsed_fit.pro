;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	qsed_fit
;
; PURPOSE:
;	
; CALLING SEQUENCE:
;   qsed_fit, phot, filts
;	
; INPUTS:
;	phot			- Input IDL structure of source photometry and data.
;   filts			- String array of photometric filters, matched to phot.
;	
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   param			- Array of best fit parameters (per object):
;					  [E(B-v), z, C_AGN, C_ELL, C_SFG, C_IRR, chi-square, DoF]
;	band			- String array of photometric filters, matched phot to template.
;	wave			- Array of central wavelength values, matched phot to template.
;
; OPTIONAL OUTPUTS:
;   FLAT			- Return chi-square closest to 1, instead of absolute minimum.
;   
; COMMENTS:
;   It is important to note that this procedure matches the input redshift of a source
;	to the template redshift (to the third decimal). This matched redshift is saved
;	as the best-fit parameter redshift in the second column of param. While this redshift
;	is generally acceptable, it can cause issues as it is not the measured redshift of 
;	the source (e.g., For source redshifts between zero and the first bin in the template
;	redshift grid, best-fit redshift can be set to zero. This can lead to a luminosity 
;	distance of zero, no calculated luminosity, etc.) It is always best to use the source 
;	redshift in all post processing!
;	
; EXAMPLES:
;
; PROCEDURES CALLED:
;	NDX2.PRO
;	
; REVISION HISTORY:
;   2016-Aug-24  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION qsed_fit, phot, $
				   filts, $
				   FLAT = flat


;; load object data
obj_vars = tag_names(phot)
;; extract photometry variables
for i = 0,n_elements(obj_vars)-1 do re = execute(obj_vars[i]+' = phot.'+obj_vars[i])

;; load template data
common _galtemp
temp_vars = scope_varname(common='_galtemp')

;; modeling variables
fit_vars = ['PARAM','BAND','WAVE']
;; degrees of freedom - color and redshift
cdof = 1
zdof = 0

;; match source redshift with template redshift
loc = value_locate(ztemp,round(z*1000.d)/1000.d)		;; must use double to round correctly!
il = where(loc eq -1,ct)
if (ct ne 0) then print, 'Template redshift range insufficient to model data', stop

;; match template bandpass with observed photometry
ifilt = []
for i = 0,n_elements(filts)-1 do ifilt = [ifilt,where(filts[i] eq wavband,/null)]
tmp = temp[ifilt,*,loc,*]
band = wavband[ifilt]
wave = obswav[ifilt]				;; needed to match properly for plotting!
ztmp = ztemp[loc]

;;commonly used array sizes
numt = n_elements(tmp[0,0,0,*])		;; number of galaxy tmplates
clen = n_elements(ebv_agn)        	;; length of color vector
zlen = n_elements(ztmp)          	;; length of redshift vector
nobj = n_elements(ztmp)		      	;; number of sources

;; create best fit array
param = dblarr(4+numt,nobj)

;tic
;; loop through all objects
for i = 0,nobj-1 do begin
	;; temporary variables
	iband = where(bin[*,i] eq 1,blen)					;; choose filters with data
	obj_flux = rebin(flux[iband,i],blen,clen)			;; match flux to temp dimensions
	obj_e_flux = rebin(e_flux[iband,i],blen,clen)		;; match error to temp dimensions
	;; solve for template contribution/coefficients
	coeff = nnlls(obj_flux,obj_e_flux,tmp[iband,*,i,*],/nneg)
	csz = size(coeff,/dim)
	dof = (blen-1)-total(coeff ne 0.,3)-zdof-cdof		;; DoF array
	;; calculate chi-square
	chi = ndx2(rebin(obj_flux,blen,clen,1,csz[-1]),rebin(obj_e_flux,blen,clen,1,csz[-1]),total(rebin(tmp[iband,*,i,*],blen,clen,1,numt,csz[-1])*rebin(reform(coeff,1,clen,1,numt,csz[-1]),blen,clen,1,numt,csz[-1]),4),d=1)
	;; ensure there is at least one positive template coefficient
	ipos = where(total(coeff gt 0.,3),poslen)
	if (poslen gt 0) then begin
	    if keyword_set(flat) then !NULL = min(abs(1.-chi[ipos]/dof[ipos]),imin) else $
	                              !NULL = min(chi[ipos]/dof[ipos],imin)
	endif else stop
	;; find minimum chi-square
	ind = array_indices(chi,ipos[imin])
	;; fill best-fit array
	param[*,i] = [ebv_agn[ind[0]],ztmp[i],reform(coeff[ind[0],ind[1],*,ind[2]]),chi[ind[0],ind[1],ind[2]],dof[ind[0],ind[1],ind[2]]]
endfor
;toc

;; where negative AGN coefficient was set to zero, also set E(B-V) to zero
param[0,where(param[2,*] eq 0.,/NULL)] = 0.

;; combine and output modeling parameters and object data
fit_str = strjoin(fit_vars+':'+fit_vars,',')
obj_str = strjoin(obj_vars+":"+obj_vars,",")
re = execute('data = soa2aos({'+obj_str+'})')
re = execute('sed_out = {'+fit_str+',data:data}')

return, sed_out


END










;; from when i was but a lowly procedure...

;; save modeling parameters
;sav_vars = ['PARAM','BAND','WAVE',obj_vars]
;sav_str = strjoin(sav_vars,',')

;re = execute('save,'+sav_vars+',/compress,file="fits_"+date_str')
;if (n_elements(suffx) eq 0) then suffx = '.sav' else suffx = '-'+suffx+'.sav'
;re = execute('save,'+sav_str+',/compress,file="fits"+suffx')


