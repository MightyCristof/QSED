;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	ndx2
;	
; PURPOSE:
;	Calculate chi-square statistic in N-dimensional space.
;   
; CALLING SEQUENCE:
;   chi = ndx2( data, uncert, theory, [, DIM= ] )
; INPUTS:
;	data			- Array of measured values.
;	uncert			- Array of measurement errors, the same dimensions as data.
;	theory			- Array of model values, the same dimensions as data.
;
; OPTIONAL INPUTS:
;	DIM				- Scalar indicating the dimension over which to sum in TOTAL().
;					  Default is zero (sum over all elements).
;
; OUTPUTS:
;   chi				- Chi-square statistic.
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
;   2014-Dec-18  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION ndx2, data, $
               uncert, $
               theory, $
               DIM = dim


if (n_elements(dim) eq 0) then dim = 0
return, total(((data - theory)/uncert)^2, dim)


END

