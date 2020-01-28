;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;	qsed_bootstrap
;
; PURPOSE:
;	Batch large samples of sources for call to qsed_zs_multi, and combine output.
;	
; CALLING SEQUENCE:
;   qsed_batch, dir
;	
; INPUTS:
;	files			- String array of input data file for SED modeling.
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
PRO qsed_bootstrap, files, $
                    galtemp, $
                    niter, $
                    TEST = test


;; load template grid variables
load_gt, galtemp
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

;; directory for all batched output file
file_mkdir, 'bootstrap_output'
pushd, 'bootstrap_output'

;; create date string
date_str = string(y, format='(I4.2)') + $
           string(mon, format=fs) + $
           string(d, format=fs) + '-'

;; batch output string
;nfiles = n_elements(files)
;fmt = '(i0'+strtrim(ceil(alog10(nfiles))+1,2)+')'
;if (nfiles eq 1) then batch_str = '' else re = execute('batch_str = string(lindgen(nfiles),format="'+fmt+'")')


;tic
;; batch sources and run fitting
;for f = 0,n_elements(files)-1 do begin
    f = 0
	print, 'Fitting: '+files[f]
	restore,files[f]
	nobj = n_elements(obs)						;; number of sources in file

    if keyword_set(test) then nobj = test
    ebv_sigm = dblarr(8,nobj)
	;; iterate over each object
	for i = 0,nobj-1 do begin
	    ;; pull and replicate individual source
	    this_obs = replicate(obs[i],niter)
	    
	    nbands = n_elements(band)
	    ;; resample photometry 
        for b = 0,nbands-1 do begin
            this_obs.mag[b] += this_obs.e_mag[b]*randomn(seed,niter)
            this_obs.flux[b] = magflux(this_obs.mag[b],this_obs.e_mag[b],band[b],err=this_err)
            this_obs.e_flux[b] = this_err
        endfor
	    
	    sed_out = qsed_fit(this_obs,band)
	    fit_vars = ['PARAM','BAND','WAVE']
	    for v = 0,n_elements(fit_vars)-1 do re = execute(fit_vars[v]+' = SED_OUT.'+fit_vars[v])
	    obj_vars = tag_names(sed_out.obj_data)
	    obj_data = sed_out.obj_data
        ;; best-fit == closest reduced chi-square to 1
        rchi = param[-2,*]/param[-1,*]
        rchi_best = min(abs(rchi-1.),ibest)
        
        ;; construct full SED output array for all objects
        if (i eq 0) then begin
            param_nobj = dblarr(n_elements(param[*,ibest]),nobj)
            band_nobj = band
            wave_nobj = wave
            obj_data_nobj = replicate(obj_data[ibest],nobj)
        endif
        
        ;; best-fit SED for each object        
        param_nobj[*,i] = param[*,ibest]
        obj_data_nobj[i] = obj_data[ibest]

        ebv_sigm[0,i] = param[0,ibest]
        ebv_sigm[1:2,i] = (moment(param[0,*]))[0:1]	    ;; Bootstrap MEAN and STDDEV
        ebv_sigm[3,i] = median(param[0,*])              ;; Bootstrap MEDIAN
        resistant_mean,param[0,*],3.0,mn,sigmn,nrej,goodvec=ig
        ebv_sigm[4:5,i] = (moment(param[0,ig]))[0:1]    ;; Bootstrap RESISTANT MEAN and STDDEV
        ebv_sigm[6,i] = median(param[0,ig])             ;; Bootstrap RESISTANT MEDIAN
        !NULL = min(abs(median(param[0,ig])-param[0,ig]),imin)
        ebv_sigm[7,i] = rchi[imin]
    endfor
    save,ebv_sigm,file='ebv_sigm.sav'
;endfor

;; restore variables to original names
param = param_nobj
for v = 0,n_elements(obj_vars)-1 do re = execute(obj_vars[v]+' = obj_data_nobj.'+obj_vars[v])
;; save concatenated SED modeling variables in top directory
popd
sav_vars = [fit_vars,obj_vars]
sav_str = strjoin(sav_vars,',')
re = execute('save,'+sav_str+',/compress,file="fits.sav"')

popd


END





;p = plot(ebv_sigm[0,*],ebv_sigm[1,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='MEAN')
;p = plot(ebv_sigm[0,*],ebv_sigm[3,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='MEDIAN')
;p = plot(ebv_sigm[0,*],ebv_sigm[4,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='RESISTANT MEAN')
;p = plot(ebv_sigm[0,*],ebv_sigm[6,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='RESISTANT MEDIAN')




