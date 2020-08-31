;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   load_gt
;   
; PURPOSE:
;   Load SED template grid into COMMON block "_GALTEMP".
;   
; CALLING SEQUENCE:
;   load_gt, file, [, /PUSH ]
;	
; INPUTS:
;   file			- String containing the name of the template grid data.
;	
; OPTIONAL INPUTS:
;   /PUSH			- Push to directory where template grid file is stored.
;					  Change directory path for personal use.
;	
; OUTPUTS:
;   _GALTEMP		- Loads variables to named common block.
; 
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   Loads temp, wavband, obswav, ztemp, and ebv_agn.
;	temp			- Array of grid templates convolved with instrument bandpasses.
;					  Dimensions of TEMP == [WAVBAND, EBV_AGN, ZTEMP, # of components]
;					  IDL> help, temp
;						TEMP (_GALTEMP) DOUBLE    = Array[13, 120, 8000, 4]
;	wavband			- String array of filter names (e.g., 'SDSS1', 'SDSS2', ...).
;	obswav			- Array of central wavelengths, matched to WAVBAND.
;	ztemp			- Array of redshift values used in constructing TEMP.
;	ebv_agn			- Array of E(B-V) values used in constructing TEMP.
;	
; EXAMPLES:
;	IDL> load_gt,'galtemp_sed4.sav',/push
;		% Compiled module: LOAD_GT.
;	IDL> common _galtemp
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2015-Jun-23  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO load_gt, file, $
			 PUSH = push
			 
			 
common _galtemp, temp, $
                 wavband, $
                 obswav, $
                 ztemp, $
                 ebv_agn

if keyword_set(push) then restore, '~/Research/sed_modeling/galtemps/'+file else $
                          restore, file


END



