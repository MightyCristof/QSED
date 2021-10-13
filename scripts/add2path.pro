;; temporarily add directory to IDL !PATH

cd,current=dir
!PATH = !PATH + ':' + Expand_Path('+'+dir)


END

