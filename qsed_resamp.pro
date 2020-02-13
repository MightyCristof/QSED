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
PRO qsed_resamp, files, $
                 galtemp, $
                 comp, $
                 niter, $
                 IND = ind, $
                 TEST = test, $
                 FLAT = flat


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

;; directory for all batched output file
;; ONLY NECESARRY IF WE ARE ALSO BATCHING LARGE NUMBER OF SOURCES (look at later)
;file_mkdir, 'resamp_output'
;pushd, 'resamp_output'

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

    if keyword_set(ind) then obs = obs[ind]
    if keyword_set(test) then obs = obs[0:test-1]
	nobj = n_elements(obs)						;; number of sources in file
    ;; arrays for resampling results
    ebv_sigm = dblarr(2,nobj)
    red_sigm = dblarr(2,nobj)
    lir_sigm = dblarr(2,nobj)
    flx_sigm = dblarr(2,nobj)
    
    nrej = lonarr(nobj)
    ;; boolean for bad fits
    bad_fit = bytarr(nobj)
    agn_perc = dblarr(nobj)

	;; iterate over each object
	for i = 0,nobj-1 do begin
	    ;; pull and replicate individual source
	    this_obs = replicate(obs[i],niter)
	    
	    nbands = n_elements(band)
	    ;; resample input: photometry
        for b = 0,nbands-1 do begin
            binct = obs[i].bin[b]
            if (binct gt 0.) then begin
                this_obs.mag[b] += this_obs.e_mag[b]*randomn(seed,niter)
                this_obs.flux[b] = magflux(this_obs.mag[b],this_obs.e_mag[b],band[b],err=this_err)
                this_obs.e_flux[b] = this_err
            endif
        endfor
	    ;; resample input: redshift
	    this_obs.z += this_obs.zerr*randomn(seed,niter)
	    izlo = where(this_obs.z lt 0.,nzlo)
	    if (nzlo gt 0) then this_obs[izlo].z = 0.
	    
	    ;; run SED fitting
	    if keyword_set(flat) then sed_out = qsed_fit(this_obs,band,/flat) else $
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
        
        ;; record the number of best-fit models not requiring an AGN component
        !NULL = where(param[2,*] gt 0.,nagn)
        agn_perc[i] = nagn*1./niter
        
        ebv_dist = reform(param[0,*])
        red_dist = this_obs.z
        c_a = reform(param[2,*])
        lir_dist = l_agn(6.,ebv_dist,red_dist,c_a)
        dl2 = dlum(red_dist,/sq)
        flx_dist = lir_dist/(4.*!const.pi*dl2)
        rchi = param[-2,*]/param[-1,*]
        ;; closest E(B-V) to the mean
        !NULL = min(abs(median(ebv_dist)-ebv_dist),iloc)
        best_ebv = ebv_dist[iloc]
            ;; find closest realization(s)
        iiebv = ebv_dist eq best_ebv
        ibest = where(iiebv,nbest)
        ;; more than one realization, pick best chi-square
        if (nbest ne 1) then begin
            if keyword_set(flat) then !NULL = min(abs(1.-rchi[ibest]),imin) else $
                                      !NULL = min(rchi[ibest],imin)
            ibest = ibest[imin]
            nbest = n_elements(ibest)
        endif
        ;; sanity check
        if (nbest ne 1) then stop
        ;; best-fit SED for each object        
        param_resamp[*,i] = param[*,ibest]
        obj_data_resamp[i] = obj_data[ibest]
        ;; input resampling results
        ebv_sigm[*,i] = [median(ebv_dist),medabsdev(ebv_dist)]
        red_sigm[*,i] = [median(red_dist),medabsdev(red_dist)]
        lir_sigm[*,i] = [median(lir_dist),medabsdev(lir_dist)]
        flx_sigm[*,i] = [median(flx_dist),medabsdev(flx_dist)]        
    endfor
    
;; save resampled fitting
resamp_vars = ['MAG','FLUX','E_FLUX','Z']
for v = 0,n_elements(resamp_vars)-1 do re = execute(resamp_vars[v]+'_RESAMP = obj_data_resamp.'+resamp_vars[v])
sav_vars = [resamp_vars+'_RESAMP',['EBV','RED','LIR','FLX']+'_SIGM','NREJ','BAD_FIT']
sav_str = strjoin(sav_vars,',')
re = execute('save,'+sav_str+',/compress,file="resamp_fits.sav"')
;endfor

;; restore variables to original names and pull original data
param = param_resamp
for v = 0,n_elements(obj_vars)-1 do re = execute(obj_vars[v]+' = obs.'+obj_vars[v])
;; save concatenated SED modeling variables in top directory
sav_vars = [fit_vars,'AGN_PERC',obj_vars]
sav_str = strjoin(sav_vars,',')
re = execute('save,'+sav_str+',/compress,file="fits.sav"')

popd


END





;p = plot(ebv_sigm[0,*],ebv_sigm[1,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='MEAN')
;p = plot(ebv_sigm[0,*],ebv_sigm[3,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='MEDIAN')
;p = plot(ebv_sigm[0,*],ebv_sigm[4,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='RESISTANT MEAN')
;p = plot(ebv_sigm[0,*],ebv_sigm[6,*],'o',sym_size=0.5,sym_filled=1,color='blue',xlog=1,ylog=1,xra=[0.01,100],yra=[0.01,100],aspect_ratio=1,ytitle='RESISTANT MEDIAN')





        ;if keyword_set(old_method) then begin
        ;;; remove outliers
        ;resistant_mean,param[0,*],3.,mn,sigmn,nr,goodvec=ig
        ;nrej[i] = nr
        ;if (nr eq niter) then begin
        ;    ;; if resistant mean unsuccessful
        ;    ;; flag source
        ;    bad_fit[i] = 1
        ;    ebv = param[0,*]
        ;    rchi = param[-2,*]/param[-1,*]
        ;    ;; pick best-fit chi-square
        ;    if keyword_set(flat) then min_chi = min(abs(1.-rchi),ibest) else $
        ;                              min_chi = min(rchi,ibest)
        ;endif else begin
        ;    ;; if resistant mean successful
        ;    ebv = param[0,ig]
        ;    red = this_obs[ig].z
        ;    c_a = param[2,ig]
        ;    lir = l_agn(6.,ebv,red,c_a,/log)
        ;    dl2 = dlum(red,/sq)
        ;    flx = lir-alog10(4.*!const.pi*dl2)
        ;    flx[where(lir eq 0.,/null)] = 0.
        ;    rchi = param[-2,ig]/param[-1,ig]
        ;    ;; closest E(B-V) to the mean
        ;    !NULL = min(abs(mn-ebv),iloc)
        ;    best_ebv = ebv[iloc]
        ;    ;; find closest realization(s)
        ;    iiebv = ebv eq best_ebv
        ;    ibest = where(iiebv,nbest)
        ;    ;; more than one realization, pick best chi-square
        ;    if (nbest ne 1) then begin
        ;        if keyword_set(flat) then !NULL = min(abs(1.-rchi[ibest]),imin) else $
        ;                                  !NULL = min(rchi[ibest],imin)
        ;        ibest = ibest[imin]
        ;        nbest = n_elements(ibest)
        ;    endif
        ;    ;; sanity check
        ;    if (nbest ne 1) then stop
        ;endelse
        ;;; best-fit SED for each object        
        ;param_nobj[*,i] = param[*,ibest]
        ;obj_data_nobj[i] = obj_data[ibest]
        ;;; input resampling results
        ;ebv_sigm[*,i] = [median(ebv),moment(ebv)]
        ;red_sigm[*,i] = [median(red),moment(red)]
        ;lir_sigm[*,i] = [median(lir),moment(lir)]
        ;flx_sigm[*,i] = [median(flx),moment(flx)]
        ;endif



