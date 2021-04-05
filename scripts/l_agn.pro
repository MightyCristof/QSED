;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;   l_agn
;
; PURPOSE:
;   Calculate AGN luminosity at a desired rest wavelength using the Asesf+10 AGN tempalate.
; 
; CALLING SEQUENCE:
;   lum = l_agn( wav0, color, redshift, c_agn, [, /LOG ] )
;	
; INPUTS:
;   wav0            - Scalar value of desired wavelength in microns.
;   ebv             - Vector containing the color excess E(B-V) of source AGN.
;   z               - Vector containing redshifts of input sources.
;   c_agn           - Vector containing the template contribution/coefficients 
;                     from SED modeling procedure.
;
; OPTIONAL INPUTS:
;   /LOG            - Set keyword to output luminosity in log space.
;
; OUTPUTS:
;	nulnu0          - Luminosity of source in cgs units [erg/s].
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
;	IDL> lum6 = l_agn(6.,ebv,z,coeff,/log)
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
FUNCTION l_agn, wav0, $
                color, $
                redshift, $
                c_agn, $
                LOG = log


;; load template component variables
common _comp
;; remove extraneous dimensions
col = reform(color)
red = reform(redshift)
c_a = reform(c_agn)
;; calculate nu*fnu in cgs
kap0 = interpol(comp.kap,comp.wav,wav0)
agn0 = interpol(comp.agn,comp.wav,wav0)*1e-29
fnu0 = c_a*agn0*10.^(-0.4*kap0*col)
nu0 = !const.c/(wav0/1e6)/(1+red)
nufnu0 = nu0*fnu0
;; calculate luminosity
nulnu0 = 4.*!const.pi*dlum(red,/sq)*nufnu0
;; return in log space
if keyword_set(log) then nulnu0 = alog10(nulnu0) > (-9999.)

return, nulnu0


END







