;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;	qsed_zs_multi
;
; PURPOSE:
;	
; CALLING SEQUENCE:
;   qsed_zs_multi, phot, filts
;	
; INPUTS:
;	phot			- Input IDL structure of source photometry and data.
;   filts			- String array of photometric filters, matched to phot.
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   param			- Array of best fit parameters (per object):
;					  [E(B-v), z, C_AGN, C_ELL, C_SFG, C_IRR, chi-square, DoF]
;	band			- String array of photometric filters, matched phot to template.
;	obswv			- Array of central wavelength values, matched phot to template.
;
; OPTIONAL OUTPUTS:
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
PRO qsed_zs_multi, phot, $
				   filts

		
;; directory name
fs = '(I2.2)'
caldat, julday(), mon, d, y, h, m, s
date_str = string(y, format='(I4.2)') + $
          string(mon, format=fs) + $
          string(d, format=fs) + '_' + $
          string(h, format=fs) + '_' + $
          string(m, format=fs) + '_' + $
          string(s, format=fs) + $
          '.sav'
			
;; load object data
phot_tags = tag_names(phot)
;; extract photometry variables
for i = 0,n_elements(phot_tags)-1 do re = execute(phot_tags[i]+' = phot.'+phot_tags[i])

;; load template data
common _galtemp
temp_vars = scope_varname(common='_galtemp')

;; degrees of freedom: color and redshift
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
band = filts[ifilt]
obswv = obswav[ifilt]				;; needed to match properly for plotting!
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
	if (poslen gt 0) then !NULL = min(chi[ipos]/dof[ipos],imin) else stop
	;; find minimum chi-square
	ind = array_indices(chi,ipos[imin])
	;; fill best-fit array
	param[*,i] = [ebv_agn[ind[0]],ztmp[i],reform(coeff[ind[0],ind[1],*,ind[2]]),chi[ind[0],ind[1],ind[2]],dof[ind[0],ind[1],ind[2]]]
endfor
;toc

;; where AGN coefficient set to zero and best fit, ensure E(B-V) is 0
param[0,where(param[2,*] eq 0.,/NULL)] = 0.

;; save modeling parameters
sav_vars = ['PARAM','BAND','OBSWV',phot_tags]
sav_vars = strjoin(sav_vars,',')
re = execute('save,'+sav_vars+',/compress,file="fits_"+date_str')


END





