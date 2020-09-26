PRO concat_sed_phot, file, $
                     sav_str


nsrcs = 0ll
ob = []
_sdss = []
_xdqso = []
_zsupp = []
_wise = []
_unwise = []
_ukidss = []
_twom = []
_galex = []
_zrang = []
_zerr = []
_zgood = []
_clean = []
_bands = []
_fourw = []
_nomsk = []
_nodup = []

for i = 0,n_elements(file)-1 do begin
    restore,file[i]
    nsrcs  += n_elements(obs)
    ob      = [ob,obs]    
    _sdss   = [_sdss,iisdss]
    _xdqso  = [_xdqso,iixdqso]
    _zsupp  = [_zsupp,iizsupp]
    _wise   = [_wise,iiwise]
    _unwise = [_unwise,iiunwise]
    _ukidss = [_ukidss,iiukidss]
    _twom   = [_twom,iitwom]
    _galex  = [_galex,iigalex]
    _zrang  = [_zrang,iizrang]
    _zerr   = [_zerr,iizerr]
    _zgood  = [_zgood,iizgood]
    _clean  = [_clean,iiclean]
    _bands  = [_bands,iibands]
    _fourw  = [_fourw,iifourw]
    _nomsk  = [_nomsk,iinomsk]
    _nodup  = [_nodup,iinodup]
endfor

obs = ob
iisdss =   _sdss   
iixdqso =  _xdqso 
iizsupp =  _zsupp 
iiwise =   _wise 
iiunwise = _unwise
iitwom =   _twom
iiukidss = _ukidss
iigalex =  _galex
iizrang =  _zrang
iizerr =   _zerr  
iizgood =  _zgood 
iiclean =  _clean 
iizrang =  _zrang 
iibands =  _bands 
iifourw =  _fourw 
iinomsk =  _nomsk 
iinodup =  _nodup 

save,nsrcs,obs,band, $
     iisdss,iizsupp,iixdqso,iiwise,iiunwise,iiukidss,iitwom,iigalex, $
     iizrang,iizerr,iizgood,iiclean,iibands,iifourw,iinomsk,iinodup, $
     /compress,file=sav_str


END





