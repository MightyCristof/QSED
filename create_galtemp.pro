;-----------------------------------------------------------------------------------------
; NAME:                                                                      IDL Procedure
;   create_galtemp
;   
; PURPOSE:
;   Create galaxy+AGN template grid over the specified color excess and redshift space.
;
;   Convolves the template components with the bandpass transmission curves and stores 
;   each separately within the TEMP grid. Nuclear obscuration is simulated on the AGN 
;   component only.
;   
; CALLING SEQUENCE:
;   create_galtemp, save_file, ebv_agn, ztemp, instr
;	
; INPUTS:
;   savfile		- String containing the output file name.
;	ebv_agn			- Vector of color excess values for constructing template grid.
;   ztemp			- Vector of redshift values for constructing template grid.
;	instr			- String array containing the desired instruments to be
;                     used in constructing the template grid.
;	
; OPTIONAL INPUTS:
;	
; OUTPUTS:
;   temp			- Template SED grid.
;	wavband			- String array of photometric filters.
;   obswav			- Array of central wavelength values matched to wavband.
;   
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	Instruments recognized:
;       UV: GALEX - 'galex'
;		optical: SDSS - 'sdss', Hyper Suprime-Cam - 'hsc'
;		NIR: UKIDSS - 'uk', 2MASS = 'twom', FLAMEX - 'flmx'
;		MIR: WISE - 'wise', IRAC - 'irac', MIPS - 'mips', PACS - 'pacs'
;
;	Common E(B-V) arrays:
; 		0 <= E(B-V) <= 15           10.^(dindgen(45)/20-1)-0.1d
; 		0 <= E(B-V) < 30	        10.^(dindgen(50)/20-1)-0.1d
;       0 <= E(B-V) <= 50           10.^(dindgen(55)/20-1)-0.1d
; 		0 <= E(B-v) < 90            10.^(dindgen(60)/20-1)-0.1d
;
;	Common redshift arrays:
;		0.000 < z < 7.999           dindgen(8000)/1000.
;       0.000 < z < 1.999           dindgen(2000)/1000.
;       0.000 < z < 1.000           dindgen(1001)/1000.
;
; EXAMPLES:
;	load_comp,'comp_temp4.sav',/push
;   load_bp,'bandpass.sav',/push
;	create_galtemp,'galtemp_sed4.sav',10.^(dindgen(120)/40-1)-0.1d,dindgen(2000)/1000., $
;                  ['sdss','wise','uk']
;	IDL> load_gt,'galtemp_sed4.sav'
;	IDL> common _galtemp
;	IDL> help
;		% At $MAIN$
;		EBV_AGN (_GALTEMP)
;                		DOUBLE    = Array[120]
;		OBSWAV (_GALTEMP)
;                		DOUBLE    = Array[13]
;		TEMP (_GALTEMP) DOUBLE    = Array[13, 120, 8000, 4]
;		WAVBAND (_GALTEMP)
;                		STRING    = Array[13]
;		ZTEMP (_GALTEMP)
;                		DOUBLE    = Array[8000]
;	
; PROCEDURES CALLED:
;	LOAD_COMP.PRO, LOAD_BP.PRO
;	
; REVISION HISTORY:
;   2015-Jul-06  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
PRO create_galtemp, savfile, $
                    ebv_agn, $
                    ztemp, $
                    instr
                                 

;; load template components and bandpass filters
common _comp
common _bp

instr = strupcase(instr)
;; filter identifier + central wavelength array
wavband = []
obswav = []

;; compile tags and wavelengths for template grid
foreach inst,instr do begin
	case inst of
		'GALEX': begin
			wavband = [wavband,'GALEX1','GALEX2']
			obswav = [obswav,0.153862,0.231566]
			end
        'SDSS': begin
			wavband = [wavband,'SDSS1','SDSS2','SDSS3','SDSS4','SDSS5']
			obswav = [obswav,0.354,0.475,0.622,0.763,0.905]
			end
		'HSC': begin
			wavband = [wavband,'HSC1','HSC2','HSC3','HSC4','HSC5']
			obswav = [obswav,0.4754,0.6175,0.7711,0.8898,0.9762]
			end
		'PANS': begin
			wavband = [wavband,'PANS1','PANS2','PANS3','PANS4','PANS5']
			obswav = [obswav,0.481,0.617,0.752,0.866,0.962]
		    end
		'UK': begin
			wavband = [wavband,'UK1','UK2','UK3','UK4']
			obswav = [obswav,1.0305,1.2483,1.6313,2.2010]
			end
		'TWOM': begin
			wavband = [wavband,'TWOM1','TWOM2','TWOM3']
			obswav = [obswav,1.235,1.662,2.159]
			end
		'FLMX': begin
			wavband = [wavband,'FLMX1','FLMX2']
			obswav = [obswav,1.250,2.215]
			end
		'WISE': begin
			wavband = [wavband,'WISE1','WISE2','WISE3','WISE4']
			obswav = [obswav,3.4,4.6,12.,22.]
			end
		'IRAC': begin
			wavband = [wavband,'IRAC1','IRAC2','IRAC3','IRAC4']
			obswav = [obswav,3.6,4.5,5.8,8.]
			end
		'IRAS': begin
		    wavband = [wavband,'IRAS1','IRAS2','IRAS3','IRAS4']
		    obswav = [obswav,12.,25.,60.,100.]
		    end
		'MIPS': begin
			wavband = [wavband,'MIPS1','MIPS2','MIPS3']
			obswav = [obswav,24.,70.,160.]
			end
		'PACS': begin
			wavband = [wavband,'PACS1','PACS2','PACS3']
			obswav = [obswav,70.,100.,160.]
			end
		'SPIRE': begin
			wavband = [wavband,'SPIRE1','SPIRE2','SPIRE3']
			obswav = [obswav,250.,350.,500.]
			end
	endcase
endforeach
		    
;;commonly used array lengths
clen = n_elements(ebv_agn)        	;; color E(B-V)
zlen = n_elements(ztemp)         	;; redshift 
tlen = n_elements(comp)  			;; number of template components (AGN+galaxies)
flen = n_elements(obswav)        	;; number of photometric filters

;; select template components
temps = ['AGN','AGN2','ELL','SBC','SFG','IRR','DST']	;; all possible templates (SED modeling procedure can handle max=5 templates)
match2,tag_names(comp),temps,icomp,itemp	            ;; match input components (use MATCH2.PRO to keep named order of TEMPS; MATCH.PRO alphabetizes; important for plotting purposes)
if (total(itemp ne -1) le 0) then stop		            ;; ensure we contain at least one valid template
temps = temps[where(itemp ne -1)]
                          
;; blank array, named variable
pts = temps+'_PTS'

;; template component shorthand
wav = comp.wav
kap = comp.kap

;; construct templates
for i = 0,n_elements(temps)-1 do begin
	if (strmatch(temps[i],'AGN*')) then begin
	    re = execute(temps[i]+' = rebin(comp.'+temps[i]+',tlen,1,clen)*transpose(10.^(-0.4 * ebv_agn # kap))')
		;agn = rebin(comp.(where(strmatch(tag_names(comp),'AGN*'))),tlen,1,clen)*transpose(10.^(-0.4 * ebv_agn # kap))		
		re = execute(temps[i]+'= rebin('+temps[i]+',tlen,zlen,clen)')
		;agn = rebin(agn,tlen,zlen,clen)											;; add redshift dimension
		agn_pts = dblarr(flen,zlen,clen)											;; blank array for convolved AGN template
		re = execute(pts[i]+' = dblarr(flen,zlen,clen)')
	endif else begin
		re = execute(temps[i]+' = rebin(comp.'+temps[i]+',tlen,zlen)')				;; add redshift dimension
		re = execute(pts[i]+' = dblarr(flen,zlen)')									;; blank arrays for convolved galaxy templates
	endelse
endfor

;; bandpass shorthand
bp_names = tag_names(bp)
for tag = 0,n_elements(bp_names)-1 do re = execute(bp_names[tag]+'=bp.'+bp_names[tag])	;; create bandpass variables
match,instr,bp_names,ibp,iinstr															;; match BP to desired instruments
bp_names = bp_names[iinstr]

;; bandpass frequency bins necessary for convolution
tempwav = wav # (1+ztemp)				;; wavelength x redshift
nu = !const.c/(tempwav*1e-6)			;; frequency x redshift
dnu = nu[0:-2,*]-nu[1:-1,*]				;; frequency bins
dnu = [dnu,dnu[-1,*]]					;; match previous dimensions (shouldn't ever play a part in convolution process)

;; create and normalize bandpass templates
;; determine photometric band/column and convolve galtemps
foreach instr, bp_names do begin
	print, instr
	re = execute('nbands = n_tags('+instr+')')				;; # of filters in instr
	for band = 0,nbands-1 do begin
		thru = dblarr(tlen,zlen)							;; throughput array
		for z = 0,zlen-1 do re = execute('thru[*,z] = interpol('+instr+'.(band).thru,'+instr+'.(band).wave,tempwav[*,z])')	;; step through each redshift and interpolate
		re = execute('mm = minmax('+instr+'.(band).wave)') 	;; determine the upper and lower bounds of throughput from wavelength
		ibounds = where(tempwav lt mm[0] or tempwav gt mm[1], /null)		;; find indices where template is out of bounds
		thru[ibounds] = 0.									;; clear values beyond wavelength boundaries
		thru[where(thru lt 0.,/null)] = 0.					;; ensure no negative throughput values (sanity check)
		normal = total(thru*dnu,1)        					;; normalize throughput
		thru /= rebin(reform(normal,1,zlen),tlen,zlen)		
		thru[where(~finite(thru))] = 0.						;; ensure finite values (sanity check)
		ifilt = where(strmatch(wavband,instr+strtrim(band+1,2)) eq 1,ct)	;; match bandpass filter (BAND) to template output (WAVBAND)
		if (ct eq 0) then stop								;; did you goof?
		for pt = 0,n_elements(pts)-1 do begin				;; fill the convolved template array
			if (strmatch(temps[pt],'AGN*')) then re = execute(pts[pt]+'[ifilt,*,*] = total('+temps[pt]+' * rebin(thru,tlen,zlen,clen) * rebin(dnu,tlen,zlen,clen),1)')
										         re = execute(pts[pt]+'[ifilt,*] = total('+temps[pt]+' * thru * dnu,1)')
		endfor
	endfor
endforeach
print, 'Bandpass templates created and normalized'		

;;create grid array of galaxy templates
temp = dblarr(flen,clen,zlen,n_elements(pts))
for i = 0,n_elements(pts)-1 do begin
	if (strmatch(temps[i],'AGN*')) then for j = 0,flen-1 do re = execute('temp[j,*,*,0] = transpose('+pts[0]+'[j,*,*])') else $
											                re = execute('temp[*,*,*,i] = rebin(reform('+pts[i]+',flen,1,zlen),flen,clen,zlen)')
endfor

save,temp,wavband,obswav,ztemp,ebv_agn,/compress,file=savfile
print, 'OUTPUT: '+savfile


END


