;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Ferction
;	lls
;	
; PURPOSE:
;	Compute template coefficients using a linear least squares approach.
;	Returns NaN where output is not linearly independent (determinant == 0).
;	   
; CALLING SEQUENCE:
;   cc = lls( fx, er, tt )
;
; INPUTS:
;	fx			- Array of input photometry.
;	er			- Array of photometric errors.
;	tt			- Array of template grid to compute contribution/coefficients.
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;   cc			- Array of template coefficients.
;
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   Photometric band must be the first dimension of all input arrays.
;	
; EXAMPLES:
;
; PROCEDURES CALLED:
;	
; REVISION HISTORY:
;   2014-Jun-19  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION lls, fx, $
              er, $
			  tt


;; size of template grid
sz = size(tt)
;; create coefficient array
cc = dblarr(sz[2:sz[0]])
;; squared ercertainty
er2 = er*er

;; calculate all combinations of templates
;; compute determinant
;; solve for coefficients
if sz[0] eq 2 then begin
	;; case of 1 input template
	t1 = tt[*,*,*,0]
	fxt1 = total(fx*t1/er2,1)
	t1t1 = total(t1*t1/er2,1)
	cc[*,*,0] = fxt1/t1t1
	cc[where(~finite(cc),/null)] = 0.
	return, cc
endif

case sz[sz[0]] of
    4: begin
    	;; case of 4 input templates
    	t1 = tt[*,*,*,0]
		t2 = tt[*,*,*,1]
		t3 = tt[*,*,*,2]
		t4 = tt[*,*,*,3]
    	fxt1 = total(fx*t1/er2,1)
		fxt2 = total(fx*t2/er2,1)
		fxt3 = total(fx*t3/er2,1)
		fxt4 = total(fx*t4/er2,1)
		t1t1 = total(t1*t1/er2,1)
		t1t2 = total(t1*t2/er2,1)
		t1t3 = total(t1*t3/er2,1)
		t1t4 = total(t1*t4/er2,1)
		t2t2 = total(t2*t2/er2,1)
		t2t3 = total(t2*t3/er2,1)
		t2t4 = total(t2*t4/er2,1)
		t3t3 = total(t3*t3/er2,1)
		t3t4 = total(t3*t4/er2,1)
		t4t4 = total(t4*t4/er2,1)
        det = (t1t1*t2t2*t3t3*t4t4 + t1t1*t2t3*t3t4*t2t4 + t1t1*t2t4*t2t3*t3t4 + t1t2*t1t2*t3t4*t3t4 + t1t2*t2t3*t1t3*t4t4 + t1t2*t2t4*t3t3*t1t4 + t1t3*t1t2*t2t3*t4t4 + t1t3*t2t2*t3t4*t1t4 + t1t3*t2t4*t1t3*t2t4 + t1t4*t1t2*t3t3*t2t4 + t1t4*t2t2*t1t3*t3t4 + t1t4*t2t3*t2t3*t1t4 - t1t1*t2t2*t3t4*t3t4 - t1t1*t2t3*t2t3*t4t4 - t1t1*t2t4*t3t3*t2t4 - t1t2*t1t2*t3t3*t4t4 - t1t2*t2t3*t3t4*t1t4 - t1t2*t2t4*t1t3*t3t4 - t1t3*t1t2*t3t4*t2t4 - t1t3*t2t2*t1t3*t4t4 - t1t3*t2t4*t2t3*t1t4 - t1t4*t1t2*t2t3*t3t4 - t1t4*t2t2*t3t3*t1t4 - t1t4*t2t3*t1t3*t2t4)
        cc[*,*,0] = ((t2t2*t3t3*t4t4 + t2t3*t3t4*t2t4 + t2t4*t2t3*t3t4 - t2t2*t3t4*t3t4 - t2t3*t2t3*t4t4 - t2t4*t3t3*t2t4)*fxt1 + (t1t2*t3t4*t3t4 + t1t3*t2t3*t4t4 + t1t4*t3t3*t2t4 - t1t2*t3t3*t4t4 - t1t3*t3t4*t2t4 - t1t4*t2t3*t3t4)*fxt2 + (t1t2*t2t3*t4t4 + t1t3*t2t4*t2t4 + t1t4*t2t2*t3t4 - t1t2*t2t4*t3t4 - t1t3*t2t2*t4t4 - t1t4*t2t3*t2t4)*fxt3 + (t1t2*t2t4*t3t3 + t1t3*t2t2*t3t4 + t1t4*t2t3*t2t3 - t1t2*t2t3*t3t4 - t1t3*t2t4*t2t3 - t1t4*t2t2*t3t3)*fxt4)/det
        cc[*,*,1] = ((t1t2*t3t4*t3t4 + t2t3*t1t3*t4t4 + t2t4*t3t3*t1t4 - t1t2*t3t3*t4t4 - t2t3*t3t4*t1t4 - t2t4*t1t3*t3t4)*fxt1 + (t1t1*t3t3*t4t4 + t1t3*t3t4*t1t4 + t1t4*t1t3*t3t4 - t1t1*t3t4*t3t4 - t1t3*t1t3*t4t4 - t1t4*t3t3*t1t4)*fxt2 + (t1t1*t2t4*t3t4 + t1t3*t1t2*t4t4 + t1t4*t2t3*t1t4 - t1t1*t2t3*t4t4 - t1t3*t2t4*t1t4 - t1t4*t1t2*t3t4)*fxt3 + (t1t1*t2t3*t3t4 + t1t3*t2t4*t1t3 + t1t4*t1t2*t3t3 - t1t1*t2t4*t3t3 - t1t3*t1t2*t3t4 - t1t4*t2t3*t1t3)*fxt4)/det
        cc[*,*,2] = ((t1t2*t2t3*t4t4 + t2t2*t3t4*t1t4 + t2t4*t1t3*t2t4 - t1t2*t3t4*t2t4 - t2t2*t1t3*t4t4 - t2t4*t2t3*t1t4)*fxt1 + (t1t1*t3t4*t2t4 + t1t2*t1t3*t4t4 + t1t4*t2t3*t1t4 - t1t1*t2t3*t4t4 - t1t2*t3t4*t1t4 - t1t4*t1t3*t2t4)*fxt2 + (t1t1*t2t2*t4t4 + t1t2*t2t4*t1t4 + t1t4*t1t2*t2t4 - t1t1*t2t4*t2t4 - t1t2*t1t2*t4t4 - t1t4*t2t2*t1t4)*fxt3 + (t1t1*t2t4*t2t3 + t1t2*t1t2*t3t4 + t1t4*t2t2*t1t3 - t1t1*t2t2*t3t4 - t1t2*t2t4*t1t3 - t1t4*t1t2*t2t3)*fxt4)/det
        cc[*,*,3] = ((t1t2*t3t3*t2t4 + t2t2*t1t3*t3t4 + t2t3*t2t3*t1t4 - t1t2*t2t3*t3t4 - t2t2*t3t3*t1t4 - t2t3*t1t3*t2t4)*fxt1 + (t1t1*t2t3*t3t4 + t1t2*t3t3*t1t4 + t1t3*t1t3*t2t4 - t1t1*t3t3*t2t4 - t1t2*t1t3*t3t4 - t1t3*t2t3*t1t4)*fxt2 + (t1t1*t2t3*t2t4 + t1t2*t1t2*t3t4 + t1t3*t2t2*t1t4 - t1t1*t2t2*t3t4 - t1t2*t2t3*t1t4 - t1t3*t1t2*t2t4)*fxt3 + (t1t1*t2t2*t3t3 + t1t2*t2t3*t1t3 + t1t3*t1t2*t2t3 - t1t1*t2t3*t2t3 - t1t2*t1t2*t3t3 - t1t3*t2t2*t1t3)*fxt4)/det
    end
    3: begin
    	;; case of 3 input templates
        t1 = tt[*,*,*,0]
		t2 = tt[*,*,*,1]
		t3 = tt[*,*,*,2]
    	fxt1 = total(fx*t1/er2,1)
		fxt2 = total(fx*t2/er2,1)
		fxt3 = total(fx*t3/er2,1)
		t1t1 = total(t1*t1/er2,1)
		t1t2 = total(t1*t2/er2,1)
		t1t3 = total(t1*t3/er2,1)
		t2t2 = total(t2*t2/er2,1)
		t2t3 = total(t2*t3/er2,1)
		t3t3 = total(t3*t3/er2,1)
        det = t1t1*(t2t2*t3t3 - t2t3*t2t3) - t1t2*(t1t2*t3t3 - t2t3*t1t3) + t1t3*(t1t2*t2t3 - t2t3*t1t3)
        cc[*,*,0] = ((t2t2*t3t3 - t2t3*t2t3)*fxt1 + (t1t3*t2t3 - t1t2*t3t3)*fxt2 + (t1t2*t2t3 - t1t3*t2t2)*fxt3) / det
        cc[*,*,1] = ((t2t3*t1t3 - t1t2*t3t3)*fxt1 + (t1t1*t3t3 - t1t3*t1t3)*fxt2 + (t1t3*t1t2 - t1t1*t2t3)*fxt3) / det
        cc[*,*,2] = ((t1t2*t2t3 - t2t2*t1t3)*fxt1 + (t1t2*t1t3 - t1t1*t2t3)*fxt2 + (t1t1*t2t2 - t1t2*t1t2)*fxt3) / det
    end
    2: begin
    	;; case of 2 input templates
        t1 = tt[*,*,*,0]
		t2 = tt[*,*,*,1]
    	fxt1 = total(fx*t1/er2,1)
		fxt2 = total(fx*t2/er2,1)
		t1t1 = total(t1*t1/er2,1)
		t1t2 = total(t1*t2/er2,1)
		t2t2 = total(t2*t2/er2,1)
        det = t1t1*t2t2 - t1t2*t1t2
        cc[*,*,0] = (t2t2*fxt1 - t1t2*fxt2) / det
        cc[*,*,1] = (t1t1*fxt2 - t1t2*fxt1) / det
    end
    else: begin
    	;; did you goof again?
        print, 'Insufficient parameters for LINEAR_LS'
        stop
    end
endcase

;; set non-finite values to NaN
cc[where(~finite(cc),/null)] = !VALUES.F_NAN

return, cc


END


