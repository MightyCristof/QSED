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
;   TEMP            - String array containing the template components for plotting.
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
               TEMP = temp, $
               IND = ind, $
               SAV = sav, $
               SHOW = show, $
               RESTFRAME = restframe      


;; silence math errors
!EXCEPT = 0

;; load template component variables
common _comp
;; determine which template components
if (n_elements(temp) eq 0) then components = tag_names(comp) else $
                                components = strupcase(temp)
;; all possible templates (SED modeling procedure can handle max=5 templates)
temps = ['AGN','ELL','SFG','IRR','DST']   
;; colors for plotting
col = [[204,121,167],[213,94,0],[0,158,115],[0,114,178],[240,228,66]]
;col = ['purple','purple','red','dark green','dark green','medium blue','brown']
;; match input components (use MATCH2.PRO to keep named order of TEMPS; MATCH.PRO alphabetizes; important for plotting purposes)
match2,components,temps,icomp,itemp
;; ensure we contain at least one valid template and sort
if (total(itemp ne -1) le 0) then stop		           
temps = temps[where(itemp ne -1)]
col = col[*,where(itemp ne -1)]
ntemps = n_elements(temps)

;; extract indices of sources to plot
if (n_elements(ind) eq 0) then ind = lindgen(n_elements(in_fits[0,*]))
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

;; calculate wavelength and frequency for sources and templates
if keyword_set(restframe) then begin
    objwav = obswav#(1+z)^(-1)
    objnu = (!const.c*1e6)/objwav#(1.+z)^(-1)
    tempwav = rebin(comp.wav,n_elements(comp),nobj)
    tempnu = (!const.c*1e6)/tempwav#(1.+z)^(-1)
    xtitle = '$Rest wavelength [ \mum ]$'
endif else begin
    objwav = rebin(obswav,n_elements(obswav),nobj)
    objnu = (!const.c*1e6)/objwav
    tempwav = comp.wav#reform(1+z)
    tempnu = (!const.c*1e6)/tempwav
    xtitle = '$Observed wavelength [ \mum ]$'
endelse

;; covert data from flux density [microjansky] to flux [erg/s/cm2]
err *= 1e-29 * objnu         
flux *= 1e-29 * objnu
;; reconstruct models
;; convert models from flux density [microjansky] to flux [erg/s/cm2]
agn = 1e-29 * tempnu * (coeff[0,*]##comp.(where(strmatch(tag_names(comp),'AGN*')))) * 10.^(-0.4 * comp.kap # ebv)                         ;; AGN model
for i = 1,ntemps-1 do re = execute(temps[i]+' = 1e-29 * tempnu * (coeff[i,*]##comp.'+temps[i]+')')  ;; galaxy models
re = execute('model = '+strjoin(temps,"+"))                                                         ;; coadded models

;; convert to log scale
err = abs((err)/(flux*alog(10)))
flux = alog10(flux)
for i = 0,ntemps-1 do re = execute(temps[i]+' = alog10('+temps[i]+')')
model = alog10(model)
;print, agn[value_locate(tempwav[*,0],6.),0]
;; string variables for call to TEXT()
z = strtrim(string(z,format='(d5.3)'),2)
ebv = strtrim(string(ebv,format='(d5.2)'),2)
coeff = reform(strtrim(string(coeff,format='(e10.3)'),2),ntemps,nobj)
chi = strtrim(string(chi[0,*],format='(d0.2)'),2)+'/'+strtrim(string(chi[1,*],format='(i)'),2)
label = transpose([['ID: '+strtrim(id,2)],['$\itz\rm: $'+z],['$\itE(B-V)\rm$: '+ebv],['$\chi^2 / DoF$: '+chi]])

;; plot SEDs
e = {xr:[0.05,30.],yra:[floor(min(flux[where(finite(flux))]))-1.5,ceil(max(flux[where(finite(flux))]))+2.],xlog:1, $
     xtitle:xtitle, ytitle:'$log( \nu \itF\rm_\nu  /  [ erg s^{-1} cm^{-2} ] )$', $
     nodata:1,buffer:1}
if keyword_set(show) then e.buffer = 0

for i = 0,nobj-1 do begin
    ;; plot good photometry
    ig = where(bin[*,i],/null)
    if keyword_set(sav) then p = plot(objwav[ig,i],flux[ig,i],_extra=e) else $              ;; set plotting window
                             p = plot(objwav[ig,i],flux[ig,i],_extra=e)
    for t = 0,ntemps-1 do re = execute('p = plot(tempwav[*,i],'+temps[t]+'[*,i],col=col[*,t],/ov)')   ;; plot models
    p = plot(tempwav[*,i],model[*,i],/ov)                                                           ;; plot coadded models
    p = errorplot(objwav[ig,i],flux[ig,i],err[ig,i],'o',/SYM_FILLED,LINESTYLE='',/OV)               ;; plot photometry
    ;; Model parameters
	yp = 0.80
	for t = 0,ntemps-1 do begin
	    lab = text(0.18,yp-t*0.04,label[t,i],/RELATIVE)
	    txt = text(0.68,yp-t*0.04,temps[t]+': '+coeff[t,i],col=col[*,t],/RELATIVE)  ;; template contribution
    endfor
    
	if keyword_set(sav) then if (strupcase(sav) eq 'EPS') then p.save,strtrim(id[i],2)+'.eps' else $
															   p.save,strtrim(id[i],2)+'.png'
endfor


END



