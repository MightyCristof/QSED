;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   load_comp
;   
; PURPOSE:
;   Load template components into COMMON block "_COMP".
;   
; CALLING SEQUENCE:
;   load_comp, file, [, /PUSH ]
;	
; INPUTS:
;   file			- String containing the name of the template components data.
;	
; OPTIONAL INPUTS:
;   /PUSH			- Push to directory where components file is stored.
;					  Change directory path for personal use.
;	
; OUTPUTS:
;   _COMP			- Loads variables to named common block.
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	Loads COMP, which contains:
;		wavelength (microns)
;		extinction coefficient "kappa"
;		1 AGN templates (AGN: Assef+10)
;		1 Elliptical template (Ell: Assef+10)
;		1 Star-forming template (Sfg: Kirkpatrick+15)
;		1 Irregular or star-burst template (Irr: Assef+10)
;   
; EXAMPLES:
;	IDL> load_comp,'components4.sav',/push
;		% Compiled module: LOAD_COMP.
;	IDL> common _comp
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2015-Jun-23  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO load_comp, file, $
               PUSH = push


common _comp, comp
if keyword_set(push) then restore, '~/Research/sed_models/primed_templates/'+file else $
                          restore, file


END


