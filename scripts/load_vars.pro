;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   load_vars
;   
; PURPOSE:
;   Load all variables from an IDL .SAV file into a named COMMON block.
;   
; CALLING SEQUENCE:
;   load_vars, file, block
;	
; INPUTS:
;   file			- String containing the name of the IDL .SAV file.
;	block           - String name of common block containing variables.
; OPTIONAL INPUTS:
;	
; OUTPUTS:
;   block			- String name of common block containing variables.
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   Default COMMON block name _VARS if no block name is passed.
;
; EXAMPLES:
;	IDL> load_vars,'fits.sav','_fits'
;       % Compiled module: LOAD_VARS.
;       % Compiled module: RIP_VARS.
;	IDL> common _fits
;   IDL> re = scope_varname(common='_fits')
;       IDL> help, re
;       RE              STRING    = Array[15]
;   IDL> print, strjoin(re,',')
;       BAND,BIN,DEC,E_FLUX,E_MAG,E_Z,E_ZARR,FLUX,MAG,OBJID,PARAM,RA,WAVE,Z,ZARR
;
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   2018-May-16  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO load_vars, file, $
               block


if (n_elements(block) eq 0) then begin
    message, /info, 'NO COMMON BLOCK DESIGNATED'
    retall
endif

;; load first .SAV file as an IDL_Savefile object
object = OBJ_NEW('IDL_Savefile',file)
;; extract variable names
vars = strjoin(object->names(),',')

;; create COMMON block
block = strupcase(block)
re = execute('common '+block+','+vars)
;; load variables
restore, file


END

