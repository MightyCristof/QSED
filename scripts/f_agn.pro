;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	f_agn
;
; PURPOSE:
;	Calculate the AGN fraction (intrinsic and observed) at a desired wavelength.
;	   
; CALLING SEQUENCE:
;   fa = f_agn( w0, fits, [, MODEL= ] )
;
; INPUTS:
;	w0				- Scalar value of desired wavelength in microns.
;	fits			- Best-fit parameter array from SED modeling (param).
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;	MODEL			- Return template model build (AGN+galaxy.) 
;
; COMMENTS:
;	Function assumes input wavelength is in microns.
;
; EXAMPLES:
;	IDL> load_comp,'components4.sav',/push
;		% Compiled module: LOAD_COMP.
;	IDL> fa = f_agn(1.,param,model=md)
;	IDL> help, fa
;		FA              STRUCT    = -> <Anonymous> Array[5097251]
;	IDL> help, fa, /st
;		** Structure <d6f9ccb8>, 2 tags, length=16, data length=16, refs=1:
;	   		INT             DOUBLE           0.0000000
;	   		OBS             DOUBLE           0.0000000
;	IDL> help, md
;		MD              STRUCT    = -> <Anonymous> Array[5097251]
;	IDL> help, md, /st
;		** Structure <d6f9c348>, 3 tags, length=24, data length=24, refs=1:
;	   		AGN_INT         DOUBLE           0.0000000
;	   		AGN_OBS         DOUBLE           0.0000000
;	   		GAL             DOUBLE           403.98832
;
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2017-Feb-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION f_agn, w0, $
                fits, $
                MODEL = out_model
               

;; load template components variables
common _comp

;; all possible templates (SED modeling procedure can handle max=5 templates)
temps = ['AGN','AGN2','ELL','SBC','SFG','IRR','DST']   
;; match input components
;; NOTE: use MATCH2 to keep named order of temps (MATCH alphabetizes)
;; this is extremely important for plotting purposes!			
match2,tag_names(comp),temps,icomp,itemp
if (total(itemp ne -1) le 0) then stop
temps = temps[where(itemp ne -1)]

;; pull data from fits array
ebv = reform(fits[0,*])
z = reform(fits[1,*])
coeffs = 'c_'+temps
for i = 0,n_elements(coeffs)-1 do re = execute(coeffs[i]+' = reform(fits[i+2,*])')

;; interpolate template components to desired wavelength
wav = comp.wav
nu = (!const.c*1e6)/wav
for i = 0,n_elements(temps)-1 do re = execute(temps[i]+' = interpol(1e-29*nu*comp.'+temps[i]+',wav,w0)/(1+z)')
kap = interpol(comp.kap,wav,w0)

;; compute AGN model contribution
intrinsic_agn = c_agn * agn
observed_agn = intrinsic_agn * 10.^(-0.4 * kap * ebv)
;; coadd galaxy model contribution
re = execute('galaxy = total('+"[["+strjoin(coeffs[1:-1]+" * "+temps[1:-1],"],[")+"]]"+',2)')
;; compute AGN fraction
intrinsic = intrinsic_agn/total([[intrinsic_agn],[galaxy]],2)
observed = observed_agn/total([[observed_agn],[galaxy]],2)

;; preserve model state
out_model = {agn_int:0d,agn_obs:0d,gal:0d}
out_model = replicate(out_model,n_elements(z))
out_model.(0) = intrinsic_agn
out_model.(1) = observed_agn
out_model.(2) = galaxy
model = out_model

;; fa output structure
object = {int:0d,obs:0d}
fa = replicate(object,n_elements(z))
fa.(0) = intrinsic
fa.(1) = observed

return, fa


END





