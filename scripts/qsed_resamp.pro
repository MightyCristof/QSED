;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;	qsed_resamp
;
; PURPOSE:
;	Batch large samples of sources for call to qsed_zs_multi, and combine output.
;	
; CALLING SEQUENCE:
;   qsed_resamp, phot_file, temp_file, comp_file, niter [, /VERBOSE ] )
;	
; INPUTS:
;	phot_file		- String naming the prepared photometry (e.g., "test_data_READY.sav").
;   temp_file       - String naming the model grid (e.g., "galtemp.sav").
;   comp_file       - String naming the template components (e.g., "comp.sav").
;   niter           - Number of iterations for modeling.
; 
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   fits.sav		- SED modeling output.
;   resamp.sav      - Uncertainy measurements.
;	
; OPTIONAL OUTPUTS:
;   VERBOSE         - Print alert to screen confirming progress.
;   
; COMMENTS:
;   
; EXAMPLES:
;
; PROCEDURES CALLED:
;	LOAD_VARS.PRO, QSED_FIT.PRO
;
; REVISION HISTORY:
;   2017-Feb-17  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO qsed_resamp, phot_file, $
                 temp_file, $
                 comp_file, $
                 niter, $
                 VERBOSE = verbose
                 

;; load template grid variables
load_vars,'models/'+temp_file,'_galtemp'
load_vars,'templates/'+comp_file,'_comp'
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
file_mkdir, 'output/'+fit_dir
pushd, 'output/'+fit_dir

;; create date string
date_str = string(y, format='(I4.2)') + $
           string(mon, format=fs) + $
           string(d, format=fs) + '-'

;; begin
print, '=============================================='
print, 'BEGIN - FITTING: '+phot_file
restore,'../../data/'+phot_file

;; number of sources in file
nobj = n_elements(obs)						
;; arrays for resampling results
sig_ebv = dblarr(4,nobj)
sig_red = dblarr(4,nobj)
sig_lir = dblarr(4,nobj)
sig_fir = dblarr(4,nobj)
;; fraction of realizations which contain AGN
perc_agn = intarr(nobj)

;; run SED fitting
sed_out = qsed_fit(obs,band)
fit_vars = tag_names(sed_out)
for v = 0,n_elements(fit_vars)-1 do re = execute(fit_vars[v]+' = SED_OUT.'+fit_vars[v])
obj_vars = tag_names(data)
for v = 0,n_elements(obj_vars)-1 do re = execute(obj_vars[v]+' = DATA.'+obj_vars[v])
sav_vars = [fit_vars[0:-1],obj_vars]
sav_str = strjoin(sav_vars,',')
re = execute('save,'+sav_str+',/compress,file="fits.sav"')

;; counter for verbose keyword
ncount = ceil(nobj/10.)*10.
;; resample each object and refit for uncertainties
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
    if (nzlo gt 0) then this_obs[izlo].z = abs(this_obs[izlo].z)
    ;izhi = where(this_obs.z gt 2.999,nzhi)
    ;if (nzhi gt 0) then this_obs[izhi].z = 3.
    
    ;; run SED fitting
    sed_out = qsed_fit(this_obs,band)
    fit_vars = tag_names(sed_out)
    for v = 0,n_elements(fit_vars)-1 do re = execute(fit_vars[v]+' = SED_OUT.'+fit_vars[v])
    obj_vars = tag_names(data)

    ;; construct full SED output array for all objects
    if (i eq 0) then begin
        resamp_param = dblarr(n_elements(param[*,0]),nobj)
        resamp_band = band
        resamp_wave = wave
        resamp_data = replicate(data[0],nobj)
    endif
    
    ;; SED output parameter distributions of realizations
    ebv_dist = reform(param[0,*])
    red_dist = this_obs.z
    c_a = reform(param[2,*])
    lir_dist = l_agn(6.,ebv_dist,red_dist,c_a)
    dl2_dist = dlum(red_dist,/sq)
    fir_dist = lir_dist/(4.*!const.pi*dl2_dist)

    ;; record the percentage of realizations containing AGN components
    iagn = where(c_a gt 0.,nagn)
    perc_agn[i] = round(nagn*100./niter)
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
            sig_ebv[*,i] = [ebv_dist,-1.,ebv_dist,-1.]
            sig_red[*,i] = [red_dist,-1.,red_dist,-1.]
            sig_lir[*,i] = [lir_dist,-1.,lir_dist,-1.]
            sig_fir[*,i] = [fir_dist,-1.,fir_dist,-1.]
        endif else begin
            sig_ebv[*,i] = [median(ebv_dist),medabsdev(ebv_dist),mean(ebv_dist),stddev(ebv_dist)]
            sig_red[*,i] = [median(red_dist),medabsdev(red_dist),mean(red_dist),stddev(red_dist)]
            sig_lir[*,i] = [median(lir_dist),medabsdev(lir_dist),mean(lir_dist),stddev(lir_dist)]
            sig_fir[*,i] = [median(fir_dist),medabsdev(fir_dist),mean(fir_dist),stddev(fir_dist)]
        endelse
    endif else begin
    ;; ===================================
    ;; no AGN presence in any realizations
    ;; ===================================
        rchi = param[-2,*]/param[-1,*]
        !NULL = min(rchi,ibest)
        
        sig_ebv[*,i] = -9999.
        sig_red[*,i] = [median(red_dist),medabsdev(red_dist),mean(red_dist),stddev(red_dist)]
        sig_lir[*,i] = -9999.
        sig_fir[*,i] = -9999.        
    endelse
    ;; best-fit SED for each object        
    resamp_param[*,i] = param[*,ibest]
    resamp_data[i] = data[ibest]
    
    ;; print progress to screen
    if keyword_set(verbose) then $
        if ((i ne 0) and (i mod (ncount/10) eq 0)) then print, '      - '+string(100.*i/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FITTING'
print, '=============================================='

;; save resampled fitting
samp_vars = 'RESAMP_'+fit_vars
nsav_vars = [samp_vars,'SIG_'+['EBV','RED','LIR','FIR'],'PERC_AGN']
nsav_str = strjoin(nsav_vars,',')
re = execute('save,'+nsav_str+',/compress,file="resamp.sav"')

popd


END













