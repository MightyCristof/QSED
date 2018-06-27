;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	nnlls
;	
; PURPOSE:
;	Compute template coefficients using a linear least squares approach.
;	Returns NaN where output is not linearly independent (determinant == 0).
;	   
; CALLING SEQUENCE:
;	call = nnlls( flx, err, tall )
;
; INPUTS:
;	flx				- Array of input photometry.
;	err				- Array of photometric errors.
;	tall			- Array of template grid to compute contribution/coefficients.
;
; OPTIONAL INPUTS:
;   /NNEG			- Keyword to set negative coefficients to zero.
;
; OUTPUTS:
;	call			- Array of template coefficients.
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;	
; EXAMPLES:
;
; PROCEDURES CALLED:
;	LLS.PRO
;
; REVISION HISTORY:
;   2014-Dec-12  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION nnlls, flx, $
                err, $
                tall, $
                NNEG = nneg


sz = size(tall,/dim)							;; dimensions of tall
ncomb = 0										;; number of combinations of templates
for i = sz[3],1,-1 do ncomb += factorial(sz[3])/(factorial(i)*factorial(sz[3]-i))
call = dblarr([sz[1:3],ncomb])					;; array to store all coefficients 

ncomb--											;; use ncomb to index call by combination
;; loop through template combinations
for i = sz[3],2,-1 do begin
	cb = combigen(sz[3],i)						;; set of all combinations
	for s = 0,n_elements(cb[*,0])-1 do begin
		set = cb[s,*]							;; pull the set of templates
		cc = lls(flx,err,tall[*,*,*,set])		;; fit with template set
		call[*,*,set,ncomb] = cc				;; fill coefficient array
		ncomb--
	endfor
endfor
;; call to COMBIGEN() does not handle "n choose 1"
for i = sz[3]-1,0,-1 do begin
	cc = lls(flx,err,tall[*,*,*,i])				;; fit with template set
	call[*,*,i,ncomb] = cc						;; fill coefficient array
	ncomb--
endfor

if keyword_set(nneg) then call[where(finite(call))] >= 0.	;; set negative values to 0

return, call


END


