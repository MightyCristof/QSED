;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   plot_qsed
;   
; PURPOSE:
;   Plot SEDs.
;   
; CALLING SEQUENCE:
;   plot_qsed, obswav, in_flux, in_err, in_bin, in_id, in_fits, [, IND=, SAV= ]
;	
; INPUTS:
;   obswav  		- Vector of central wavelengths.
;	in_flux			- Array of fluxes. First dimension matches obswav.
;   in_err			- Array of flux errors. First dimension matches obswav.
;	in_bin			- Array of byte flags signaling good photometry.
;   in_id           - Vector of object IDs.
;   in_fits         - Array of best-fit modeling parameter output.	
;
; OPTIONAL INPUTS:
;	IND             - Vector containing the indices to plot.
;
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   SAV             - Set this keyword to 'EPS' for encapsulated postscript.
;                     If set as keyword, default is PNG.
;
; COMMENTS:
;
; EXAMPLES:
;	
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2014-Apr-19  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO plot_qsed, obswav, $
			   in_flux, $
			   in_err, $
               in_bin, $
			   in_id, $
               in_fits, $
               IND = ind, $
               SAV = sav       


;; silence math errors
!EXCEPT = 0

;; load template component variables
common _comp
;; determine which template components
temps = ['AGN','ELL','SFG','IRR','DST']		;; all possible templates (SED modeling procedure can handle max=5 templates)
match2,tag_names(comp),temps,icomp,itemp	;; match input components (use MATCH2.PRO to keep named order of TEMPS; MATCH.PRO alphabetizes; important for plotting purposes)
if (total(itemp ne -1) le 0) then stop		;; ensure we contain at least one valid template
temps = temps[where(itemp ne -1)]
ntemps = n_elements(temps)

;; extract indices of sources to plot
if (n_elements(ind) eq 0) then ind = lindgen(n_elements(i_fits[0,*]))
nobj = n_elements(ind)
flux = in_flux[*,ind]
err = in_err[*,ind]
bin = in_bin[*,ind]
id = in_id[ind]
fits = in_fits[*,ind]

;; extract model parameters
ebv = fits[0,*]
z = fits[1,*]
coeff = fits[2:2+ntemps-1,*]
chi = fits[-2:-1,*]

;; calculate rest wavelength and frequency for sources and templates
restwav = obswav#(1+z)^(-1)
objnu = !const.c/(restwav * 1e-6)
tempwav = rebin(comp.wav,n_elements(comp),nobj)
tempnu = !const.c/(tempwav * 1e-6)
	
;; covert from flux density [microjansky] to flux [erg/s/cm2]
err *= 1e-29 * objnu         
flux *= 1e-29 * objnu
;; reconstruct models
agn = 1e-29 * tempnu * (coeff[0,*]##comp.agn) * 10.^(-0.4 * comp.kap # ebv)                         ;; AGN model
for i = 1,ntemps-1 do re = execute(temps[i]+' = 1e-29 * tempnu * (coeff[i,*]##comp.'+temps[i]+')')  ;; galaxy models
re = execute('model = '+strjoin(temps,"+"))                                                         ;; coadded models

;; convert to log scale
err = abs((err)/(flux*alog(10)))
flux = alog10(flux)
for i = 0,ntemps-1 do re = execute(temps[i]+' = alog10('+temps[i]+')')
model = alog10(model)

;; string variables for call to TEXT()
z = strtrim(string(z,format='(d5.3)'),2)
ebv = strtrim(string(ebv,format='(d5.2)'),2)
coeff = reform(strtrim(string(coeff,format='(e10.3)'),2),ntemps,nobj)
chi = strtrim(string(chi[0,*],format='(d0.2)'),2)+'/'+strtrim(string(chi[1,*],format='(i)'),2)

;; plot SEDs
col = ['purple','red','dark green','medium blue']
for i = 0,nobj-1 do begin
    ;; plot good photometry
    ig = where(bin[*,i],/null)
    p = plot(restwav[ig,i],flux[ig,i],XRA=[0.05,30],/XLOG,YRA=[-18.,-5.],/NODATA)                               ;; set plotting window
    for t = 0,ntemps-1 do re = execute('p = plot(tempwav[*,i],'+temps[t]+'[*,i],col=col[t],/ov)')               ;; plot models
    p = plot(tempwav[*,i],model[*,i],/ov)                                                                       ;; plot coadded models
    p = errorplot(restwav[ig,i],flux[ig,i],err[ig,i],'o',/SYM_FILLED,LINESTYLE='',/OV)                          ;; plot data
    p.XTITLE='$Rest wavelength [ \mum ]$' & p.YTITLE='$log( \nu \itF\rm_\nu  /  [ erg s^{-1} cm^{-2} ] )$'
    ;; Model parameters
    !NULL = text(0.18,0.80,'id: '+strtrim(id[i],2),/RELATIVE)
	!NULL = text(0.18,0.76,'$z: $'+z[i],/RELATIVE)
	!NULL = text(0.18,0.72,'E(B-V): '+ebv[i],/RELATIVE)
	!NULL = text(0.18,0.68,'$\chi^2/dof$: '+chi[i],/RELATIVE)
	yp = 0.80
	for t = 0,ntemps-1 do txt = text(0.68,yp-t*0.04,temps[t]+': '+coeff[t,i],col=col[t],/RELATIVE)              ;; template contribution

	if keyword_set(sav) then if (strupcase(sav) eq 'EPS') then p.save,strtrim(id[i],2)+'.eps' else $
															   p.save,strtrim(id[i],2)+'.png'
endfor


END



