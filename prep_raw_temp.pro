;-----------------------------------------------------------------------------------------
; NAME:                                                                   IDL Main Program
;	prep_raw_temp
;
; PURPOSE:
;	Format raw galaxy+AGN templates for SED modeling procedure.
;
; CALLING SEQUENCE:
;   .r prep_raw_temp
;   
; INPUTS:
;	Pointer to directory and template files.   
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;	comp			- IDL structure containing the template components in 
;					  the format required for SED modeling. Here we default to
;					  1 AGN template and 3 galaxy templates (passive, star-forming,
;					  and star-burst).
;   
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   If you prefer to use other templates than those listed here, simply change
;	the pointer to the input data files.
;
;	The SED modeling procedure is equipped to handle up to 5 template components,
;	(only 4 present here) such as a cold dust component (or two AGN components, which
;	I have never tested.) To do so, you need to modify the code and conversions, and
;	add the final component to the output structure (and tags).
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2014-Feb-10  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
print, 'prep_raw_temp'
print, 'Version 1.0'
print, ' '


;; input data files
temp_dir = '~/Research/sed_models/raw_templates/'
a10 = 'assef+10/lrt_templates.dat'
ext = 'assef+10/ext_law_data.dat'
k15 = 'kirkpatrick+15/Comprehensive_library/SFG2.txt'

;; read A10 templates
readcol,temp_dir+a10,wav,agn,agn2,ell,sbc,irr,format='d,d,d,d,d,d',/silent
temps = ['AGN','AGN2','ELL','SBC','IRR']
;; Assef raw data conversion from "lrt_templates.dat"
convrs = [1d-11,1d-14,1d-16,1d-17,1d-15]
;; read extinction law data
readcol,temp_dir+ext,extwav,extkap,format='d,d',/silent
;; convert extinction law wavelengths from Ångstrom to micrometers
extwav *= 1e-4
;; extrapolate extinction coeffienct to A10 wavelengths
kap = spline(extwav,extkap,wav)
;; covert A10 templates to microjansky
erg2mujy = 1e29
for i = 0,n_elements(temps)-1 do re = execute(temps[i]+' *= erg2mujy * convrs[i]')

;; read K15 SFG3 template
readcol,temp_dir+k15,wavlen,lum,dlnu,format='d,d,d'
;; covert K15 templates to microjansky
lum *= 1e7
dlnu *= 1e7
dl = !const.parsec * 100.				;; A10 templates normalized to 10 pc
fnu = lum/(4.*!const.pi*dl^2)
;; match and normalize to A10 SBc
loc = value_locate(wav,wavlen[0])		;; index of wav where wavlen begins
normal = sbc[loc]/fnu[0]				;; flux conversion factor
fnu *= normal							;; fnu normalized to A10 at wav[loc]

;; interpolate fnu to A10 wavelengths from loc to end of wav array
sfg = spline(wavlen,fnu,wav[loc+1:-1])	
;; stitch SBc[0:loc] and SFG[loc:-1]
;; SFG is now A10 until the beginning of the SFG template and spans A10 wavelength
sfg = [sbc[0:loc],sfg]					

;; create template components structure
comp = {wav: 0d, $
        kap: 0d, $
        agn: 0d, $
        agn2: 0d, $
        ell: 0d, $
        sbc: 0d, $
        irr: 0d, $
        sfg: 0d $
        }
comp = replicate(comp,n_elements(wav))

;; fill component structure
tags = tag_names(comp)		  ;; A10 all tempaltes + K15 SFG
for i = 0,n_elements(tags)-1 do re = execute('comp.(i) = '+tags[i])

save, comp, file='comp.sav'


END

