;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;	qsed_batch
;
; PURPOSE:
;	Batch large samples of sources for call to qsed_zs_multi, and combine output.
;	
; CALLING SEQUENCE:
;   qsed_batch, dir
;	
; INPUTS:
;	files			- String array of input data files for SED modeling.
; 
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   fits			- Combined output of multi-batch qsed_zs_multi calls.
;	
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   
; EXAMPLES:
;
; PROCEDURES CALLED:
;	LOAD_GT.PRO, QSED_ZS_MULTI.PRO
;
; REVISION HISTORY:
;   2017-Feb-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO qsed_batch, files


;; load template grid variables
load_gt,'galtemp_ebv90.sav'
;load_gt,'galtemp_*.sav',/push

;; create runtime directory
fs = '(I2.2)'
caldat, julday(), mon, d, y, h, m
fit_str = string(y, format='(I4.2)') + $
          string(mon, format=fs) + $
          string(d, format=fs) + '_' + $
          string(h, format=fs) + $
          string(m, format=fs)
fit_dir = 'run_' + fit_str
file_mkdir, fit_dir
pushd, fit_dir

;; directory for all batched output files
file_mkdir, 'batch_output'
pushd, 'batch_output'

;; create date string
date_str = string(y, format='(I4.2)') + $
           string(mon, format=fs) + $
           string(d, format=fs) + '-'

;tic
;; batch sources and run fitting
for i = 0,n_elements(files)-1 do begin
	print, 'Fitting: '+files[i]
	restore,files[i]
	;; create batches for fitting
	nobj = n_elements(obs)						;; number of sources in file
	bsize = 100000.								;; batch size
	nbatches = ceil(nobj/bsize)					;; number of batches
	;; create sky section string
	sky_str = strsplit(files[i],'-_',/extract)
	sky_str = date_str+sky_str[where(strmatch(sky_str,'part*'))]
	;; begin batch fitting
	nleft = nobj								;; number of sources left in batch
	for b = 0,nbatches-1 do begin
		;; create batch file string
		sav_str = sky_str+'-'+string(b,format='(I02)')
		nleft = nleft-bsize > 0								;; track the remaining number of sources
		batch_obs = obs[b*bsize:(b+1)*bsize-1 < (nobj-1)]	;; temporary data structure
		qsed_zs_multi,batch_obs,band,sav_str				;; call SED modeling procedure
	endfor	
endfor
;toc

;; concatenate batched modeling outputs into single variables
fit_file = file_search('fits-*')						;; name of all batched output files
object = OBJ_NEW('IDL_Savefile', fit_file[0])			;; create IDL object
vars = object->names()									;; pull variable names
ivar = where(~strmatch(vars,'BAND') and ~strmatch(vars,'WAVE'),nvars)	;; do not append BAND or WAVE
temp_vars = 'TEMP_'+vars[ivar]							;; create temporary variables
temp_dims = [2,1,1,1,2,2,1,1,1,2,2,1,2,1,1,1]
for i = 0,nvars-1 do re = execute(temp_vars[i]+' = []')

;; loop over batched output files and append variables
for i = 0,n_elements(fit_file)-1 do begin
	restore, fit_file[i]
	for j = 0,nvars-1 do begin
		;; match dimensions of output array
		dim = temp_dims[j]
		case dim of
		    0: re = execute(temp_vars[j]+' = ['+temp_vars[j]+','+vars[ivar[j]]+']')
			1: re = execute(temp_vars[j]+' = ['+temp_vars[j]+','+vars[ivar[j]]+']')
			2: re = execute(temp_vars[j]+' = [['+temp_vars[j]+'],['+vars[ivar[j]]+']]')
		endcase
	endfor
endfor
;; restore variables to original names
for i = 0,nvars-1 do re = execute(vars[ivar[i]]+' = '+temp_vars[i])
;; save concatenated SED modeling variables in top directory
popd
re = execute('save,'+strjoin(vars,",")+',/compress,file="fits.sav"')
popd


END

