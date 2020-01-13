;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;   MAGFLUX
;   
; PURPOSE:
;   Convert between raw catalog data magnitudes to flux density in microjanskys.
;   
; CALLING SEQUENCE:
;   out_val = magflux( value, error, band, [, ERR=, /FLUX_IN ] )
;	
; INPUTS:
;   value			- Input mag/flux.
;   error			- Input error on mag/flux.
;   band			- Bandpass filter name. 
;   
; OPTIONAL INPUTS:
;   /FLUX_IN		- Set this keyword when input is flux; calculate magnitudes.
;   
; OUTPUTS:
;   out_val			- Output mag/flux. out_val = 0 for value=-9999 or error²0. 
;                     If /FLUX_IN set, out_val = -9999 for value=0 or error<0.
;   
; OPTIONAL OUTPUTS:
;   ERR				- Return output error.
;  
; COMMENTS:
;   All magnitudes are assumed from the raw catalogs. All input and output flux are
;   in microjanskys.
;   
;   Photometric corrections from:
;   
;   SDSS     http://classic.sdss.org/dr7/algorithms/fluxcal.html
;   2MASS    https://iopscience.iop.org/article/10.1086/376474
;            http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
;   UKIDSS   http://www.ukidss.org/technical/photom/hewett-ukidss.pdf
;   WISE     http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEZMA
;   IRAC     http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/17/
;   FLAMEX   ---
;   GALEX    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=GALEX
;            https://www.aanda.org/articles/aa/pdf/2016/12/aa27782-15.pdf  Effective Area == Transmittance
;   HSC      https://www.naoj.org/Observing/Instruments/HSC/sensitivity.html
;   PACS     http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?id=Herschel/
;	SPIRE	 http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=Herschel&gname2=SPIRE
;	
;   Gemini Observatory conversion tool
;   	http://www.gemini.edu/sciops/instruments/midir-resources/imaging-calibrations/fluxmagnitude-conversion
;   Morpheus MIR conversion tool
;   	http://morpheus.phys.lsu.edu/magnitude.html   
;	IPAC
;		http://coolwiki.ipac.caltech.edu/index.php/Central_wavelengths_and_zero_points
;   
; EXAMPLES:
;   Calculate flux density from SDSS r-band magnitudes.
;   
;   IDL> seed = 100
;   IDL> mag = [-2.,3.,8.,13.,18.,23.]
;   IDL> e_mag = randomn(seed,6)/100.
;   IDL> flux = magflux(mag,e_mag,'sdss2',err=e_flux)
;   IDL> print, flux
;      2.2910063e+10   2.2910059e+08       0.0000000       0.0000000       229.10071       0.0000000
;   IDL> print, e_flux
;      4.1734544e+08       1058384.7       0.0000000       0.0000000       1.8852304       0.0000000
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2015-Jun-19  Written by C. M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION magflux, value, $
                  error, $
                  band, $
                  ERR = out_err, $
                  FLUX_IN = flux_in
                       

;; AB zero-pt flux, zero-pt flux error, zero-pt mag
f0 = 3631.
f0e = 0.
m0 = 0.

;; instrument zero points
case strupcase(band) of
	'SDSS1': m0 = 0.04
	'SDSS2': break
	'SDSS3': break
	'SDSS4': break
    'SDSS5': m0 = -0.02 
    'TWOM1': begin
    	;f0 = 1594.
        ;f0e = 27.8
    	m0 = -0.91+0.001
    	end
    'TWOM2': begin
        ;f0 = 1024.
        ;f0e = 20.
        m0 = -1.39-0.019
        end
    'TWOM3': begin
        ;f0 = 666.7
        ;f0e = 12.6
        m0 = -1.85+0.017
        end
    'UK1': m0 = -0.634
    'UK2': m0 = -0.938
    'UK3': m0 = -1.379
    'UK4': m0 = -1.900
    'WISE1': m0 = -2.699
    'WISE2': m0 = -3.339
    'WISE3': m0 = -5.174
    'WISE4': m0 = -6.620
    'IRAC1': begin
    	f0 = 280.9
    	f0e = 4.1
    	end
    'IRAC2': begin
    	f0 = 179.7
    	f0e = 2.6
    	end
    'IRAC3': begin
    	f0 = 115.
    	f0e = 1.7
    	end
    'IRAC4': begin
        f0 = 64.9
        f0e = 0.9
        end
    'FLMX1': begin 
    	f0 = 1594.
        f0e = 27.8
        end
    'FLMX2': begin
        f0 = 666.7
        f0e = 12.6
        end
    'GALEX1': break
    'GALEX2': break
    'HSC1': break
    'HSC2': break
    'HSC3': break
    'HSC4': break
    'HSC5': break
    'MIPS1': m0 = -6.77
    'MIPS2': m0 = -9.18
    'MIPS3': m0 = -10.9
    'PACS1': f0 = 0.7961
    'PACS2': f0 = 0.3868
    'PACS3': f0 = 0.1500
    'SPIRE1': f0 = 0.06
    'SPIRE2': f0 = 0.03
    'SPIRE3': f0 = 0.01
    else: print, 'ASSUMING AB PHOTOMETRY'    
endcase

;; flux2mag
if keyword_set(flux_in) then begin
	out_val = dblarr(n_elements(value))-9999.					;; assume -9999
	out_err = dblarr(n_elements(value))-9999.
	ig = where(error gt 0. and value gt 0.,gct) 				;; good input flux
	if (gct gt 0.) then begin
		out_val[ig] = m0 - 2.5*alog10(value[ig]/10.^6/f0)			;; calculate mag for all good input
		out_err[ig] = abs((2.5*error[ig])/(value[ig]*alog(10.d)))	;; output error
	endif
endif else begin
	out_val = dblarr(n_elements(value))							;; assume 0
	out_err = dblarr(n_elements(value))
	ig = where(error gt 0 and value gt -99.,gct)				;; good input mag
	if (gct gt 0.) then begin
		out_val[ig] = f0 * 10.^((m0-value[ig])/2.5) * 10.^6							;; calculate flux for all good input
		out_err[ig] = sqrt((f0e/f0)^2+((error[ig]*out_val[ig]*alog(10.))/2.5)^2)	;; output error
	endif
endelse

err = out_err													
return, out_val													
	                     
	                     
END



