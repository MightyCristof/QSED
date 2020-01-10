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
;   fullz           - String array containing redshifts from multiple catalogs, 
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
FUNCTION zorig, fullz


zcat = ['ZP','PEAKZ','ZS','ZSUPP']
iz = lonarr(n_elements(fullz))
for i = 0,n_elements(fullz)-1 do iz[i] = max(where(strsplit(fullz[i],',',/extract,/preserve_null) gt 0.))

return, zcat[iz]


END