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
;  		SDSS            STRUCT    -> <Anonymous> Array[1]
;  		UK              STRUCT    -> <Anonymous> Array[1]
;  		WISE            STRUCT    -> <Anonymous> Array[1]  
;   IDL> help, bp.sdss
;       ** Structure <1649298>, 5 tags, length=7056, data length=7056, refs=2:
;       FILT1           STRUCT    -> <Anonymous> Array[47]
;       FILT2           STRUCT    -> <Anonymous> Array[89]
;       FILT3           STRUCT    -> <Anonymous> Array[75]
;       FILT4           STRUCT    -> <Anonymous> Array[89]
;       FILT5           STRUCT    -> <Anonymous> Array[141]
;   IDL> help, bp.sdss.filt1,/st
;       ** Structure <21160e8>, 2 tags, length=16, data length=16, refs=6:
;       WAV             DOUBLE          0.29800000
;       THRU            DOUBLE           0.0000000
;
; PROCEDURES CALLED:
;	CREATE_BP.PRO
;
; REVISION HISTORY:
;   2015-Jun-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO merge_bp


;; directory of instruments
pushd,'./bandpass'
inst = file_search('*',/test_dir)
;; call to CREATE_BP() for each instrument in dirs
for i = 0,n_elements(inst)-1 do re = execute(inst[i]+' = create_bp(file_search(inst[i]+"/*"))')	
;; fill structure
bp_str = strjoin(inst+':'+inst,',')		;; concatenate output string
re = execute('bp = {'+bp_str+'}')		;; merge multi-instrument bandpass strutures
save,bp,file='bandpass.sav'
popd


END
