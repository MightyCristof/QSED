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
;	dir				- String variable containing directory of formatted flux data.
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
;	
; REVISION HISTORY:
;   2017-Feb-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO qsed_batch, dir


;; directory name
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

flux_file = file_search(dir,+'/*')
load_gt,'../galtemp_qso4.sav'

file_mkdir, 'fit_output'
pushd, 'fit_output'
;; run fitting
tic
for i = 0,n_elements(flux_file)-1 do begin
	print, 'Fitting: '+flux_file[i]
	restore,flux_file[i]
	;;BATCH HERE
	;; create batches for fitting
	nobj = n_elements(obs)
	bsize = 100000.
	nbatches = ceil(nobj/bsize)
	;param_str = make_array(nbatches,value='param_')+string(indgen(nbatches,start=1),format='(i3.3)')
	nleft = nobj
	for b = 0,nbatches-1 do begin
		nleft = nleft-bsize > 0		;; track the number of objects left to fit
		batch_obs = obs[b*bsize:(b+1)*bsize-1 < (nobj-1)]
		qso_zs_multi_github,batch_obs,band
	endfor	
endfor
toc

;; concatenate param outputs into single variable "PARAM"
fit_file = file_search('fits_*')
object = OBJ_NEW('IDL_Savefile', fit_file[0])
vars = object->names()
temp_vars = 'TEMP_'+vars
for i = 0,n_elements(temp_vars)-1 do re = execute(temp_vars[i]+' = []')

for i = 0,n_elements(fit_file)-1 do begin
	restore, fit_file[i]
	for j = 0,n_elements(vars)-1 do begin
		re = execute('sz = size('+vars[j]+')')
		if (sz[0] eq 2) then re = execute(temp_vars[j]+' = [['+temp_vars[j]+'],['+vars[j]+']]') else $
		                     re = execute(temp_vars[j]+' = ['+temp_vars[j]+','+vars[j]+']')
	endfor
endfor
	
for i = 0,n_elements(vars)-1 do re = execute(vars[i]+' = '+temp_vars[i])

var_str = strjoin(vars,',')
popd
re = execute('save,'+var_str+',/compress,file="fits.sav"')
popd


END

