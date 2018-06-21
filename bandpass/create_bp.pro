;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	create_bp
;	
; PURPOSE:
;	Convert a set of transmission curves to an IDL structure.
;	   
; CALLING SEQUENCE:
;   bp = create_bp( files )
;	
; INPUTS:
;	files			- String array containing the instrument bandpass data files.
;	   
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   bp				- IDL structure of single instrument bandpass data.
;	
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	Bandpass files must be in a format recognized by READCOL, with the first two
;	columns being wavelength and throughput. Prefixed files with numbers in ascending 
;	order will ensure that FILT# is sorting by wavelength coverage (short -> long).
;
;	Note: Bandpass files from different instruments may have differing wavelength  
;	units. Please be sure input files are in the units you need!
;
;	Technically, BP is a structure of arrays of structures. See the example below.
;   
; EXAMPLES:
;	Using SDSS:
;	
;	IDL> files = file_search('*.dat')
;	IDL> sdss = create_bp(file)
;	IDL> help, sdss
;		** Structure <31097b8>, 5 tags, length=7056, data length=7056, refs=1:
;	   	FILT1           STRUCT    -> <Anonymous> Array[47]
;	   	FILT2           STRUCT    -> <Anonymous> Array[89]
;	   	FILT3           STRUCT    -> <Anonymous> Array[75]
;	   	FILT4           STRUCT    -> <Anonymous> Array[89]
;	   	FILT5           STRUCT    -> <Anonymous> Array[141]
;
;	Examine SDSS1 (u-band):
;	IDL> help, sdss.filt1, /st
;		** Structure <300b838>, 2 tags, length=16, data length=16, refs=6:
;  		WAV             DOUBLE          0.29800000
;  		THRU            DOUBLE           0.0000000
;	IDL> help, sdss.filt1.wav
;		<Expression>    DOUBLE    = Array[47]
;	IDL> help, sdss.filt1.thru
;		<Expression>    DOUBLE    = Array[47]
;
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2015-Jun-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION create_bp, files


nfilt = n_elements(files)

;; create IDL structure for filter containing wavelength and throughput
filt = {wav:0d, $
		thru:0d}
;; output string array to concatenate instrument bandpass for RETURN
filt_out = !NULL

;; loop over files/filters
for i = 0,nfilt-1 do begin
	readcol,files[i],wav,thru,format='d,d',/silent		;; files should contain wavelength,throughput in columns 1,2
	iw = sort(wav)										;; ensure ascending wavelength
	thru = thru[iw]						
	wav = wav[iw]
	re = execute('filt'+strtrim(i+1,2)+' = replicate(filt,n_elements(wav))')	;; replicate unique filter and fill 
	re = execute('filt'+strtrim(i+1,2)+'.wav = wav')
	re = execute('filt'+strtrim(i+1,2)+'.thru = thru')
	filt_out = [filt_out,'filt'+strtrim(i+1,2)+':filt'+strtrim(i+1,2)]			;; add unique filter name to output string
endfor

filt_out = strjoin(filt_out,',')				;; concatenate output string
re = execute('bp = {'+filt_out+'}')				;; create full instrument bandpass structure

return, bp


END


