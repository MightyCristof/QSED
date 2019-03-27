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
;; Assef raw data "in units of" from "lrt_templates.dat"
in_units_of = [1e-11,1e-14,1e-16,1e-17,1e-15]
for i = 0,n_elements(temps)-1 do re = execute(temps[i]+' *= in_units_of[i]')
;; read extinction law data
readcol,temp_dir+ext,extwav,extkap,format='d,d',/silent
;; convert extinction law wavelengths from Ångstrom to micrometers
;; 1 A = 10^4 um
extwav /= 1e4
;; extrapolate extinction coeffienct to A10 wavelengths
kap = spline(extwav,extkap,wav)
;; covert A10 templates to microjansky
;; 1 uJy = 1e-29 erg/s/cm2/Hz
cgs2ujy = 1e29
for i = 0,n_elements(temps)-1 do re = execute(temps[i]+' *= cgs2ujy')

;; read K15 SFG3 template
readcol,temp_dir+k15,kirkwav,lum,dlnu,format='d,d,d'
;; covert K15 templates to microjansky
lum *= 1e7                              ;; Watt to cgs
dlnu *= 1e7
dl = (10.*!const.parsec) * 100.			;; A10 templates are normalized to 10 pc; convert
fnu = lum/(4.*!const.pi*dl^2)           ;; to centimeters and use same value for K15
;; match and normalize to A10 SBc
ikirk = 232;0
loc = value_locate(wav,kirkwav[ikirk])		;; index of wav where kirkwav begins
normal = sbc[loc]/fnu[ikirk]				;; flux conversion factor
fnu *= normal							;; fnu normalized to A10 at wav[loc]

;; interpolate fnu to A10 wavelengths from loc to end of wav array
sfg = spline(kirkwav,fnu,wav[loc+1:-1])	
;; stitch SBc[0:loc] and SFG[loc:-1]
;; SFG is now A10 until the beginning of the SFG template and spans A10 wavelength
sfg = [sbc[0:loc],sfg]					

;; create template components structure
comp = {wav: 0d, $
        kap: 0d, $
        agn: 0d, $
;        agn2: 0d, $
        ell: 0d, $
;        sbc: 0d, $
        irr: 0d, $
        sfg: 0d $
        }
comp = replicate(comp,n_elements(wav))

;; fill component structure
tags = tag_names(comp)		  ;; A10 all tempaltes + K15 SFG
for i = 0,n_elements(tags)-1 do re = execute('comp.(i) = '+tags[i])

save, comp, file='comp.sav'


END

