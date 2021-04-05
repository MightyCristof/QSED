;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;   zorig
;
; PURPOSE:
;   Given a list of redshifts and string array of all redshift information, 
;   return which catalog the "best" redshift is from.
;   
; CALLING SEQUENCE:
;   zcat = zorig( fullz )
;
; INPUTS:
;   zbin            - String array containing byte flags from multiple catalogs, 
;                     separated by ','.
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   zcat            - String array containing the redshift type/catalog.
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;   Used in conjunction with the QSED modeling. If catalogs change at READ_SED_PHOT,
;   then they need to be updated here as well.
;   
;   It is faster to loop and split than to split and loop using toArray().
;
; EXAMPLES:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   2018-Jun-26  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION zorigin, zbin, $
                  zstr, $
                  ZPREF = zpref


npar = n_params()

zcat = ['ZP','ZPXD','ZS','ZSAGES','ZSUPP']
iz = lonarr(n_elements(zbin))
for i = 0,n_elements(zbin)-1 do iz[i] = max(where(strsplit(zbin[i],',',/extract,/preserve_null) eq 1.))

;; if full redshift string is passed to the function, return preferred redshift
if (npar eq 2) then begin
    zpref = dblarr(n_elements(zbin))
    for i = 0,n_elements(zbin)-1 do zpref[i] = (strsplit(zstr[i],',',/extract,/preserve_null))[iz[i]]
endif

;; where ZP is not valid, replace with negative redshift value
inoz = where(iz eq -1,nnoz)
if (nnoz gt 0) then for i = 0,nnoz-1 do zpref[inoz[i]] = -1.*max(strsplit(zstr[inoz[i]],',',/extract,/preserve_null))

;; return string array of redshift origin ('ZP','PEAKZ','ZS','ZSUPP') for each object
return, zcat[iz]


END