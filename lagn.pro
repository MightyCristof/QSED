;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	lagn
;
; PURPOSE:
;	Calculate AGN luminosity at a desired wavelength using the Asesf+10 AGN tempalate.
; 
; CALLING SEQUENCE:
;   lum = lagn( w0, ebv, z, coeff, [, /LOG ] )
;	
; INPUTS:
;	w0				- Scalar value of desired wavelength in microns.
;	color			- Vector containing the color excess E(B-V) of source AGN.
;	redshift		- Vector containing redshifts of input sources.
;	coefficient		- Vector containing the template contribution/coefficients 
;					  from SED modeling procedure.
;
; OPTIONAL INPUTS:
;   /LOG			- Set keyword to output luminosity in log space.
;
; OUTPUTS:
;	lum				- Luminosity of source in cgs units [erg/s].
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	Function assumes input wavelength is in microns.
;	
;   The call to LUMDIST() is computationally intensive for large input arrays. To 
;	reduce the computation time, a sparse redshift array is used in the call to 
;	LUMDIST(), and the results are then interpolated to the input redshift values.
;	The redshift range covers 0.001 < z < 11.303.
;
; EXAMPLES:
;	IDL> load_comp,'components4.sav',/push
;	IDL> lum6 = lagn(6.,param[0,*],z,param[2,*],/log)
;	IDL> for i = 0,10 do print, param[2,i], lum6[i]
;	       0.0000000       0.0000000
;	       0.0000000       0.0000000
;	   1.4887328e-17       42.994325
;	       0.0000000       0.0000000
;	       0.0000000       0.0000000
;	   9.8133307e-17       41.413075
;	       0.0000000       0.0000000
;	       0.0000000       0.0000000
;	       0.0000000       0.0000000
;	   1.1253022e-17       36.595763
;	       0.0000000       0.0000000
;
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2017-May-12  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION lagn, w0, $
			   color, $
			   redshift, $
			   coefficient, $
			   LOG = log


;; load template component variables
common _comp
;; remove extraneous dimensions
ebv = reform(color)
z = reform(redshift)
coeff = reform(coefficient)

;; flux density Fnu [erg/s/cm2/Hz] at desired wavelength
fnu0 = interpol(comp.agn,comp.wav,w0)
kap0 = interpol(comp.kap,comp.wav,w0)
fnu = coeff * fnu0 * 10d^(-0.4*kap0*ebv) * 1e-29			;; convert microjansky to cgs units
nu = !const.c/(w0*1e-6)									;; convert w0 from micron to m
;; luminosity distance dL [cm2]
testz = 10.^(dindgen(110)/100)-0.999d					;; range of z values to calculate dL
dl = lumdist(testz,h0=70.,omega_m=0.3,lambda0=0.7)		;; luminosity distance in Mpc; dL = (1+z)c ºdz/H(z)
dl = interpol(dl,testz,z)								;; interpolation to input redshift
dl *= 1e6 * !const.parsec * 1e2 						;; luminosity distance converted from Mpc to cm
;; luminosity L [erg/s]
lum = 4d * !const.pi * dl^2 * fnu * nu

;; return in log space
if keyword_set(log) then lum = alog10(lum) > 0.

return, lum


END


