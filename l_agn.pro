;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;   l_agn
;
; PURPOSE:
;   Calculate AGN luminosity at a desired wavelength using the Asesf+10 AGN tempalate.
; 
; CALLING SEQUENCE:
;   lum = l_agn( w0, ebv, z, coeff, [, /LOG ] )
;	
; INPUTS:
;   w0              - Scalar value of desired wavelength in microns.
;   ebv             - Vector containing the color excess E(B-V) of source AGN.
;   z               - Vector containing redshifts of input sources.
;   coeff           - Vector containing the template contribution/coefficients 
;                     from SED modeling procedure.
;
; OPTIONAL INPUTS:
;   /LOG            - Set keyword to output luminosity in log space.
;
; OUTPUTS:
;	lum             - Luminosity of source in cgs units [erg/s].
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   Function assumes input wavelength is in microns.
;	
;   The call to LUMDIST() is computationally intensive for large input arrays. To 
;   reduce the computation time, a sparse redshift array is used in the call to 
;   LUMDIST(), and the results are then interpolated to the input redshift values.
;   The redshift range covers min(z) < z < 29.9.
;
; EXAMPLES:
;	IDL> load_comp,'components4.sav',/push
;	IDL> lum6 = l_agn(6.,param[0,*],z,param[2,*],/log)
;	IDL> for i = 0,10 do print, param[2,i], lum6[i]
;           0.0000000       0.0000000
;           0.0000000       0.0000000
;       1.4887328e-17       42.994325
;           0.0000000       0.0000000
;           0.0000000       0.0000000
;       9.8133307e-17       41.413075
;           0.0000000       0.0000000
;           0.0000000       0.0000000
;           0.0000000       0.0000000
;       1.1253022e-17       36.595763
;           0.0000000       0.0000000
;
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2017-May-12  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION l_agn, w0, $
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

;; determine AGN template
iagn = where(strmatch(tag_names(comp),'AGN'),agnct)    ;; match COMP to find AGN template
if (agnct ne 1) then stop                              ;; ensure there is an AGN template
;; extract AGN template and wavelength from COMP
agn = comp.(iagn)
wav = comp.wav
agn0 = interpol(agn,wav,w0)

;; convert AGN template to Fnu (cgs units), uJy to erg/s/cm2/Hz
fnu0 = coeff * agn0 * 1e-29

;; calculate extinction
kap0 = interpol(comp.kap,wav,w0)
ext = 10.^(-0.4*kap0*ebv)
;; attenuate Fnu [erg/s/cm2/Hz]
fnu0 *= ext

;; frequency at w0 rest-frame, SED fit in observed frame -> nu_obs = nu_emit/(1+z)
nu0 = (!const.c*1e6)/w0/(1.+redshift)

;; calculate flux
nufnu0 = fnu0 * nu0
;print, nufnu0[0]
;; luminosity distance in cgs
testz = 10.^(dindgen(150)/100.)-(1.-min(z)>0.)          ;; range of z values to calculate dL
dl = lumdist(testz,h0=70.,omega_m=0.3,lambda0=0.7)      ;; luminosity distance in Mpc; dL = (1+z)c ºdz/H(z)
dl = interpol(dl,testz,z)                               ;; interpolation to input redshift
dl2 = (dl * 1e6 * !const.parsec * 1e2)^2                ;; luminosity distance converted from Mpc to cm

;; luminosity L in cgs (ergs/s)
lum = 4.*!const.pi*dl2*nufnu0

;; return in log space
if keyword_set(log) then lum = alog10(lum) > 0.

return, lum


END


;wav_cm = wav * 1e-4
;w0_cm = w0 * 1e-4
;;; calculate flux density Fnu [erg/s/cm2/Hz] at desired wavelength
;for i = 0,nsrc-1 do fnu0_cgs[i] = interpol(agn_cgs[*,i],wav_cm,w0_cm)
;fnu0_cgs = interpol(agn_cgs,wav_cm,w0_cm)
;fnu_cgs = coeff * fnu0 * 
;
;fnu0 = interpol(agn,wav,w0)
;kap0 = interpol(comp.kap,comp.wav,w0)
;fnu = coeff * fnu0 * 10d^(-0.4*kap0*ebv) * 1e-29        ;; convert microjansky to erg/s/cm2/Hz
;nu = !const.c/(w0*1e-6)                                 ;; convert w0 from micron to m
;;; luminosity L [erg/s/Hz]
;lum = 4d * !const.pi * dl^2 * fnu * nu
;
;;; return in log space
;if keyword_set(log) then lum = alog10(lum) > 0.
;
;return, lum
;
;
;END


