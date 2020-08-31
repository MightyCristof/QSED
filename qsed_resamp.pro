;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;	qsed_resamp
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
PRO qsed_resamp, file, $
                 galtemp, $
                 comp, $
                 niter
                 

;; load template grid variables
load_gt, galtemp
load_comp, comp
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

;; create date string
date_str = string(y, format='(I4.2)') + $
           string(mon, format=fs) + $
           string(d, format=fs) + '-'

;; begin
print, 'Fitting: '+file
restore,file

;; number of sources in file
nobj = n_elements(obs)						
;; arrays for resampling results
ebv_sigm = dblarr(4,nobj)
red_sigm = dblarr(4,nobj)
lir_sigm = dblarr(4,nobj)
fir_sigm = dblarr(4,nobj)

nrej = lonarr(nobj)
;; boolean for bad fits
bad_fit = bytarr(nobj)
perc_agn = dblarr(nobj)

;; iterate over each object
for i = 0,nobj-1 do begin
    ;; pull and replicate individual source
    this_obs = replicate(obs[i],niter)
    
    ;; resample input: photometry
    nbands = n_elements(band)
    for b = 0,nbands-1 do begin
        binct = obs[i].bin[b]
        if (binct gt 0.) then begin
            this_obs.mag[b] += this_obs.e_mag[b]*randomn(seed,niter)
            this_obs.flux[b] = magflux(this_obs.mag[b],this_obs.e_mag[b],band[b],err=this_err)
            this_obs.e_flux[b] = this_err
        endif
    endfor
    
    ;; resample input: redshift
    this_obs.zerr *= randomn(seed,niter)
    this_obs.z += this_obs.zerr
    izlo = where(this_obs.z lt 0.,nzlo)
    if (nzlo gt 0) then this_obs[izlo].z = 0.
    ;izhi = where(this_obs.z gt 0.999,nzhi)
    ;if (nzhi gt 0) then this_obs[izhi].z = 0.999
    
    ;; run SED fitting
    sed_out = qsed_fit(this_obs,band)
    fit_vars = ['PARAM','BAND','WAVE']
    for v = 0,n_elements(fit_vars)-1 do re = execute(fit_vars[v]+' = SED_OUT.'+fit_vars[v])
    obj_vars = tag_names(sed_out.obj_data)
    obj_data = sed_out.obj_data

    ;; construct full SED output array for all objects
    if (i eq 0) then begin
        param_resamp = dblarr(n_elements(param[*,0]),nobj)
        band_resamp = band
        wave_resamp = wave
        obj_data_resamp = replicate(obj_data[0],nobj)
    endif
    
    ;; SED output parameter distributions of realizations
    ebv_dist = reform(param[0,*])
    red_dist = this_obs.z
    c_a = reform(param[2,*])
    lir_dist = l_agn(6.,ebv_dist,red_dist,c_a)
    dl2 = dlum(red_dist,/sq)
    fir_dist = lir_dist/(4.*!const.pi*dl2)

    ;; record the percentage of realizations containing AGN components
    iagn = where(param[2,*] gt 0.,nagn)
    perc_agn[i] = nagn*100./niter
    ;; focus on realizations with AGN component
    if (nagn gt 0.) then begin
    ;; ============================
    ;; AGN presence in realizations
    ;; ============================
        ebv_dist = ebv_dist[iagn]
        red_dist = red_dist[iagn]
        lir_dist = lir_dist[iagn]
        fir_dist = fir_dist[iagn]
        rchi = param[-2,iagn]/param[-1,iagn]
        ;; closest E(B-V) to the mean
        del_ebv = min(abs(median(ebv_dist)-ebv_dist),iloc)
        best_ebv = ebv_dist[iloc]
        ;; find closest realization(s)
        iiebv = ebv_dist eq best_ebv
        ibest = where(iiebv,nbest)
        ;; more than one realization, pick best chi-square
        if (nbest ne 1) then begin
            !NULL = min(rchi[ibest],imin)
            ibest = iagn[ibest[imin]]
            nbest = n_elements(ibest)
        endif
        ;; sanity check
        if (nbest ne 1) then stop
        ;; input resampling results
        if (nagn eq 1) then begin
            ebv_sigm[*,i] = [ebv_dist,-1.,ebv_dist,-1.]
            red_sigm[*,i] = [red_dist,-1.,red_dist,-1.]
            lir_sigm[*,i] = [lir_dist,-1.,lir_dist,-1.]
            fir_sigm[*,i] = [fir_dist,-1.,fir_dist,-1.]
        endif else begin
            ebv_sigm[*,i] = [median(ebv_dist),medabsdev(ebv_dist),mean(ebv_dist),stddev(ebv_dist)]
            red_sigm[*,i] = [median(red_dist),medabsdev(red_dist),mean(red_dist),stddev(red_dist)]
            lir_sigm[*,i] = [median(lir_dist),medabsdev(lir_dist),mean(lir_dist),stddev(lir_dist)]
            fir_sigm[*,i] = [median(fir_dist),medabsdev(fir_dist),mean(fir_dist),stddev(fir_dist)]
        endelse
    endif else begin
    ;; ===================================
    ;; no AGN presence in any realizations
    ;; ===================================
        rchi = param[-2,*]/param[-1,*]
        !NULL = min(rchi,ibest)
        
        ebv_sigm[*,i] = -9999.
        red_sigm[*,i] = [median(red_dist),medabsdev(red_dist),mean(red_dist),stddev(red_dist)]
        lir_sigm[*,i] = -9999.
        fir_sigm[*,i] = -9999.        
    endelse
    ;; best-fit SED for each object        
    param_resamp[*,i] = param[*,ibest]
    obj_data_resamp[i] = obj_data[ibest]
endfor
    
;; save resampled fitting
resamp_vars = ['MAG','FLUX','E_FLUX','Z','ZERR']
for v = 0,n_elements(resamp_vars)-1 do re = execute(resamp_vars[v]+'_RESAMP = obj_data_resamp.'+resamp_vars[v])
sav_vars = [resamp_vars+'_RESAMP',['EBV','RED','LIR','FIR']+'_SIGM','NREJ','BAD_FIT']
sav_str = strjoin(sav_vars,',')
re = execute('save,'+sav_str+',/compress,file="resamp.sav"')
;endfor

;; restore variables to original names and pull original data
param = param_resamp
for v = 0,n_elements(obj_vars)-1 do re = execute(obj_vars[v]+' = obs.'+obj_vars[v])
;; save concatenated SED modeling variables in top directory
sav_vars = [fit_vars,'perc_agn',obj_vars]
sav_str = strjoin(sav_vars,',')
re = execute('save,'+sav_str+',/compress,file="fits.sav"')

popd


END









