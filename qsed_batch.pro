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


;; create runtime directory
fs = '(I2.2)'
caldat, julday(), mon, d, y, h, min
date_str = string(y, format='(I4.2)') + $
          string(mon, format=fs) + $
          string(d, format=fs) + '_' + $
          string(h, format=fs) + $
          string(min, format=fs)
fit_dir = 'run_' + date_str
file_mkdir, fit_dir
pushd, fit_dir

;; load template grid variables
load_gt,'galtemp_sed4.sav',/push

;; directory for all batched output files
file_mkdir, 'fit_output'
pushd, 'fit_output'

;tic
;; batch sources and run fitting
for i = 0,n_elements(files)-1 do begin
	print, 'Fitting: '+files[i]
	restore,files[i]
	;; create batches for fitting
	nobj = n_elements(obs)						;; number of sources in file
	bsize = 100000.								;; batch size
	nbatches = ceil(nobj/bsize)					;; number of batches
	nleft = nobj								;; number of sources left in batch
	for b = 0,nbatches-1 do begin
		nleft = nleft-bsize > 0					;; track the remaining number of sources
		batch_obs = obs[b*bsize:(b+1)*bsize-1 < (nobj-1)]	;; temporary data structure
		qsed_zs_multi,batch_obs,band			;; call SED modeling procedure
	endfor	
endfor
;toc

;; concatenate batched modeling outputs into single variables
fit_file = file_search('fits_*')						;; name of all batched output files
object = OBJ_NEW('IDL_Savefile', fit_file[0])			;; create IDL object
vars = object->names()									;; pull variable names
temp_vars = 'TEMP_'+vars								;; create temporary variables
for i = 0,n_elements(temp_vars)-1 do re = execute(temp_vars[i]+' = []')

;; loop over batched output files and append variables
for i = 0,n_elements(fit_file)-1 do begin
	restore, fit_file[i]
	for j = 0,n_elements(vars)-1 do begin
		re = execute('sz = size('+vars[j]+')')
		if (sz[0] eq 2) then re = execute(temp_vars[j]+' = [['+temp_vars[j]+'],['+vars[j]+']]') else $
		                     re = execute(temp_vars[j]+' = ['+temp_vars[j]+','+vars[j]+']')
	endfor
endfor
;; restore variables to original names
for i = 0,n_elements(vars)-1 do re = execute(vars[i]+' = '+temp_vars[i])
;; save concatenated SED modeling variables in top directory
var_str = strjoin(vars,',')
popd
re = execute('save,'+var_str+',/compress,file="fits.sav"')
popd


END

