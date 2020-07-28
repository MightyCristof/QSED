;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   load_bp
;   
; PURPOSE:
;   Load bandpass wavelength and throughput into COMMON block "_BP".
;   
; CALLING SEQUENCE:
;   load_bp, file, [, /PUSH ]
;	
; INPUTS:
;   file			- String containing the name of the template data.
;	
; OPTIONAL INPUTS:
;   /PUSH			- Push to directory of stored components.
;	
; OUTPUTS:
;   _BP				- Loads variables to named common block.
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	Loads BP, a nested structure which contains the wavelength and throughput for 
;   each instrument with bandpass data.
;   
; EXAMPLES:
;	IDL> load_bp,'bandpass.sav',/push
;		% Compiled module: LOAD_BP.
;	IDL> common _bp
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2015-Jun-23  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO load_bp, file, $
             PUSH = push


common _bp, bp
if keyword_set(push) then restore, '~/Research/sed_modeling/bandpass/'+file else $
						  restore, file


END


