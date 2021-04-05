c     This function returns the luminosity distance given a
c     redshift. This function should only be called after initializing
c     KCA or PZA.
c
      real*8 function DL(z)
      implicit real*8 (a-h,o-z)

      real*8 dcom(10000),dlum(10000),dmod(10000),dvol(10000),rstar(10000)
      common /distance/dcom,dlum,dmod,dvol,rstar       

      real*8 zin2(10000)      
      common /distance2/zin2,dz,zinmin,zinmax,nin2
      
      real*8 om,ol,ok,H0,DH
      common /cosmo/om,ol,ok,H0,DH

      real*8 DC,DM

c     Calculate the position of zin that holds the closest, but lower,
c     value of z. Also check that everything is ok.
      i = (z-zinmin)/dz+1
      if(i.gt.nin2) i=nin2

c     Complete the integration up to the redshift needed.
      npt = 1000
      zmin = zin2(i)
      zmax = z
      delz = (zmax-zmin)/float(npt-1) 
      dadd = 0.d0
      do k=1,npt
         ztemp = zmin + delz*float(k-1)
         x   = ztemp
         val = 1.d0/sqrt(om*(1.d0+x)**3+ok*(1.d0+x)**2+ol)
         if (k.ne.1) dadd = dadd + 0.5d0*delz*(val+vold) 
         vold  = val
      enddo
      DC = dcom(i) + DH*dadd
            
c     Calculate DM and DL.
      if(ok.gt.0.d0) then
         DM = DH*sinh(sqrt(ok)*DC/DH)/sqrt(ok)
      else if(ok.lt.0d0) then
         DM = DH*sin(sqrt(-ok)*DC/DH)/sqrt(-ok)
      else
         DM = DC
      endif
      DL = (1.d0+z)*DM

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function returns the comoving volume given a redshift.
c
      real*8 function vc(z)
      implicit real*8 (a-h,o-z)

c     If argument is negative, return -1.
      if(z.lt.0.d0) then
         vc = -1.d0
      else
         vc = (4.d0*3.14159d0/3.d0)*DL(z)**3/(1.d0+z)**3
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Given a certain vector comp, holding the specific luminosities of
c     each component (just like the one returned by the kca function),
c     this subroutine returns the fluxes you would expect at a given
c     redshift. This function if for creating mock galaxies using the
c     templates.
      subroutine get_mags(comp,ebvx,igmx,z,jymodtot)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4)

      real*8 comp(*),jymodtot(*)
      real*8 ebvx,igmx
      real*8 vec(NSMAX)

      real*8 alpha_norm(NSMAX)
      common /alphanorm/alpha_norm

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax 

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

      real*8 jyzero(NCMAX),sat(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,sat,con,lbar

      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      real*8 jymod(NSMAX,NCMAX)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 tigm(NWMAX)


c     Transform the component vectors back to the natural units of the
c     code.
      vecfac = DL(z)**2*1e10*3e-9/(1+z)
      do l = 1,nspec
         vec(l) = comp(l)*alpha_norm(l)/vecfac
      enddo

      do jchan=1,nchan
         do kwave=1,nwave
            wgt(jchan,kwave) = getweight(z,jchan,kwave)
         enddo
         call getrange(jchan)
      enddo

      do kwave=1,nwave
         tigm(kwave) = transmit(bcen(kwave),zuse,guse)
      enddo
      do l=1,nspec
         do j=1,nchan
            jymod(l,j) = 0.d0
            do k=jwmin(j),jwmax(j)
               if(l.ne.1) then
                  dust = 1.d0
               else
                  dust = 10.d0**(-0.4d0*tau(k)*euse)
               endif
               jymod(l,j) = jymod(l,j) + c(j)*spec(l,k)*wgt(j,k)*
     *              dust*tigm(k)
            enddo
         enddo
      enddo

      do j=1,nchan
         jymodtot(j)   = 0.d0
         do l = 1,nspec
            jymodtot(j)   = jymodtot(j) + vec(l)*jymod(l,j)
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Clear the matrix
      subroutine clearmat(a,b,maxdim,nwave)
      real*8 b(maxdim),a(maxdim,maxdim)

      do k1=1,nwave
         b(k1) = 0.0
         do k2=1,nwave
            a(k1,k2) = 0.0
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Symmetrize the matrix
      subroutine symmat(a,b,maxdim,nwave)
      real*8 b(maxdim),a(maxdim,maxdim)
      
      do k1=1,nwave
         do k2=1,k1-1 
            a(k1,k2) = a(k2,k1)
         enddo
      enddo

      return 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Get chi2 of the fit using fluxes.
c
      subroutine get_chi_jy(jy,ejy,jymod,jyuse,nchan,chi2)
      implicit real*8 (a-h,o-z)
      
      real*8 jy(*),jymod(*),ejy(*)
      integer jyuse(*)


      chi2 = 0.d0

      do j=1,nchan
         if(jyuse(j).ge.1) then
            chi2 = chi2 + ((jy(j)-jymod(j))/ejy(j))**2
         endif
      enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     Get chi2 of the fit using mags.
c
      subroutine get_chimag(mag,emag,magmod,maguse,nchan,chi2)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32)

      real*8 mag(*),magmod(*),emag(*)
      integer maguse(*)
      
      real*8 jy(NCMAX),ejy(NCMAX),jymod(NCMAX)

      real*8 jyzero(NCMAX),sat(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,sat,con,lbar


      do jchan=1,nchan
         jy(jchan)   = jyzero(jchan)*10.d0**(-0.4d0*mag(jchan))
         ejy(jchan)  = jy(jchan)*emag(jchan)*2.5d0/dlog(10.d0)
         jymod(jchan)= jyzero(jchan)*10**(-0.4*magmod(jchan))
      enddo


      chi2 = 0.d0

      do j=1,nchan
         if(maguse(j).eq.1) then
            chi2 = chi2 + ((jy(j)-jymod(j))/ejy(j))**2
         endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine noblank(fname,i,j)
      implicit real*8 (a-h,o-z)

      character fname*(*)

      do i=1,100
         if(fname(i:i).ne.' ') goto 100
      enddo
 100  continue
      do j=100,1,-1
         if(fname(j:j).ne.' ') goto 110
      enddo
 110  continue

      return
      end
