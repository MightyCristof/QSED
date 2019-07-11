;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;   mw_ext_corr
;   
; PURPOSE:
;   Correct input flux for Galactic extinction using IR dust maps.
;   
; CALLING SEQUENCE:
;   mw_ext_corr, ra, dec, flux_red, band
;
; INPUTS:
;   ra			    - Vector of position right ascension.
;   dec             - Vector of position declination.
;   flux_red        - Vector of input flux to be deredenned.
;   band            - String containing the photometric band for processing.
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   flux_unred      - Vector of output deredenned flux values.
;      
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;   Path to dust map directory must be set prior to running script.
;   
; EXAMPLES:
;
; REVISION HISTORY:
;   2016-Dec-14  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION mw_ext_corr, ra, $
					  dec, $
					  flux_red, $
					  band


case strupcase(band) of
	'SDSS1': wav = 0.354
	'SDSS2': wav = 0.475
	'SDSS3': wav = 0.622
	'SDSS4': wav = 0.763
    'SDSS5': wav = 0.905
    'TWOM1': wav = 1.235
    'TWOM2': wav = 1.662
    'TWOM3': wav = 2.159
    'UK1': wav = 1.0305
    'UK2': wav = 1.2483
    'UK3': wav = 1.6313
    'UK4': wav = 2.2010
    'WISE1': wav = 3.4
    'WISE2': wav = 4.6
    'WISE3': wav = 12.
    'WISE4': wav = 22.
    'IRAC1': wav = 3.6
    'IRAC2': wav = 4.5
    'IRAC3': wav = 5.8
    'IRAC4': wav = 8.
    'FLMX1': wav = 1.250
    'FLMX2': wav = 2.215
    'GALEX1': wav = 0.153862
    'GALEX2': wav = 0.231566
    'HSC1': wav = 0.4754
    'HSC2': wav = 0.6175
    'HSC3': wav = 0.7711
    'HSC4': wav = 0.8898
    'HSC5': wav = 0.9762
    'MIPS1': wav = 24.
    'MIPS2': wav = 70.
    'MIPS3': wav = 160.
    'PACS1': wav = 70.
    'PACS2': wav = 100.
    'PACS3': wav = 160.
    'SPIRE1': wav = 250.
    'SPIRE2': wav = 350.
    'SPIRE3': wav = 500.
    else: begin
        print, 'NO VALID WAVELENGTH'
        return, !NULL
        end
endcase

;;convert to Angstrom
temp_wav = replicate(wav*10000.,n_elements(ra))
;; convert to Galactic coords
euler,ra,dec,gall,galb,select=1
;; obtain E(B-V)_MW from dust map
ebv = dust_getval(gall,galb,/interp,/noloop)
;; correct flux vector for Galactic extinction
ccm_unred,temp_wav,flux_red,ebv,flux_unred

return, flux_unred


end
