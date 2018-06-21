;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	merge_bp
;	
; PURPOSE:
;	Concatenate multiple bandpass structures.
;	
; CALLING SEQUENCE:
;   bp = merge_bp( dirs )
; INPUTS:
;	dirs			- String array of instrument directories containing bandpass
;					  data files.
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   bp				- IDL structure of multi-instrument bandpass data.
;	
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   BP is a nested structure of structures of arrays of structures.
;
; EXAMPLES:
;	IDL> dirs = file_search(/test_dir)
;	IDL> print, dirs
;		flmx galex hsc irac mips pacs sdss spire twom uk wise
;	IDL> bp = merge_bp(dirs)
;		% Compiled module: READCOL.
;		% Compiled module: REMCHAR.
;		% Compiled module: GETTOK.
;		% Compiled module: STRNUMBER.
;	IDL> help, bp
;		** Structure <1a53e08>, 11 tags, length=289216, data length=289216, refs=1:
;   	FLMX            STRUCT    -> <Anonymous> Array[1]
;  		GALEX           STRUCT    -> <Anonymous> Array[1]
;  		HSC             STRUCT    -> <Anonymous> Array[1]
;  		IRAC            STRUCT    -> <Anonymous> Array[1]
;  		MIPS            STRUCT    -> <Anonymous> Array[1]
;  		PACS            STRUCT    -> <Anonymous> Array[1]
;  		SDSS            STRUCT    -> <Anonymous> Array[1]
;  		SPIRE           STRUCT    -> <Anonymous> Array[1]
;  		TWOM            STRUCT    -> <Anonymous> Array[1]
;  		UK              STRUCT    -> <Anonymous> Array[1]
;  		WISE            STRUCT    -> <Anonymous> Array[1]
;   
; PROCEDURES CALLED:
;	CREATE_BP.PRO
;
; REVISION HISTORY:
;   2015-Jun-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION merge_bp, dirs


;; call to CREATE_BP() for each instrument in dirs
for i = 0,n_elements(dirs)-1 do re = execute(dirs[i]+' = create_bp(file_search(dirs[i]+"/*"))')	

bp_out = strjoin(dirs+':'+dirs,',')		;; concatenate output string
re = execute('bp = {'+bp_out+'}')		;; merge multi-instrument bandpass strutures

return, bp


END
