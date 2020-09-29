      subroutine sno(deltam,nenu,dm4e,probsnodd,probsnonn,rsno)
      implicit real*4 (a-h,o-z)
      real*4 probsnodd(nenu),probsnonn(nenu)
      real*4 dm4e(nenu),rsno(5)
      real*8 probdd,probnn
      real*8 arg(10),pval(10),perr(10),plo(10),phi(10),gint
      real*8 esno(124),ccd(124),ccn(124),eseld(124),eseln(124),
     %                                   esmud(124),esmun(124)
      real*8 ccsnod(124),ccsnon(124),eselsnod(124),eselsnon(124),
     %                               esmusnod(124),esmusnon(124)
      real*8 cc0,cc1,cc2,aa0,aa1,probddtst,probnntst,aeetst
      common /fluxsno1/esno,ccd,ccn,eseld,eseln,esmud,esmun
      common /fluxsno2/ccsnod,ccsnon,eselsnod,eselsnon,
     %                               esmusnod,esmusnon
      character*10 name(10)
      external minfcn   ! need to pass function name to MINUIT
c............................................................
c............................................................


      
      do i=1,124
         deltam4e=log10(deltam/4.e6/esno(i))
         probdd=divdif(probsnodd,dm4e,nenu,deltam4e,2)
         probnn=divdif(probsnonn,dm4e,nenu,deltam4e,2)
         ccsnod(i)=probdd*ccd(i)
         ccsnon(i)=probnn*ccn(i)
         eselsnod(i)=probdd*eseld(i)
         eselsnon(i)=probnn*eseln(i)
         esmusnod(i)=(1.-probdd)*esmud(i)
         esmusnon(i)=(1.-probnn)*esmun(i)
c         write(89,*)deltam4e,probdd,probnn
      enddo

c set up parameters
      call mnparm(1,'c0', 0.317d0, .001d0,  0.0d0, 0.0d0, ierr)
      call mnparm(2,'c1', 0.0039d0,.0001d0, 0.0d0, 0.0d0, ierr)
      call mnparm(3,'c2',-0.001d0, .001d0,  0.0d0, 0.0d0, ierr)
      call mnparm(4,'a0', 0.046d0, .001d0,  0.0d0, 0.0d0, ierr)
      call mnparm(5,'a1',-0.016d0, .001d0,  0.0d0, 0.0d0, ierr)
c      call mnparm(6,'ph', 1.d0,    .01d0,   0.0d0, 0.0d0, ierr)

c minimize

      call mnexcm(minfcn,'MIGRAD',   arg,0,ierr,0)!perform MIGRAD minimization  

c get parameters and errors
      do i=1,5
         call mnpout(i,name(i),pval(i),perr(i),plo(i),phi(i),ierr)
         rsno(i)=pval(i)
      enddo

      
      return
      end
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine minfcn(npara,grad,fval,xval,iflag,dummy)

c MINUIT objective function
c
c supplied in call to function:
c
c npara - number of parameter modified (or something)
c xval - vector of current parameters
c futil - (optional) name of utility routine [not used here]
c iflag - 1: read input data, etc.
c         2: calculate gradients [not used here]
c         3: final call (print, etc. if needed...)
c
c returned by function:
c
c fval - calculated function value
c grad - (optional) vector of fist derivatives [not used here]

      implicit none
      integer npara
      integer iflag
      real*8 fval,xval(*),grad(*)
      real*8 c0,c1,c2,a0,a1,phib
      real*8 chi2sno,enu,pdd,aee,pnn
      real*8 ccsnodd,ccsnonn,eselsnodd,eselsnonn,esmusnodd,esmusnonn
      real*8 ccsnod(124),ccsnon(124),eselsnod(124),eselsnon(124),
     %                               esmusnod(124),esmusnon(124)
      real*8 esno(124),ccd(124),ccn(124),eseld(124),eseln(124),
     %                                   esmud(124),esmun(124)
      common /fluxsno1/esno,ccd,ccn,eseld,eseln,esmud,esmun
      common /fluxsno2/ccsnod,ccsnon,eselsnod,eselsnon,
     %                               esmusnod,esmusnon
      
      external dummy
      integer i
      
      fval = 0.

      c0=xval(1)
      c1=xval(2)
      c2=xval(3)
      a0=xval(4)
      a1=xval(5)
c      phib=xval(6)

      chi2sno=0.
      do i=1,124
         enu=esno(i)
         pdd=c0+c1*(enu-10.)+c2*(enu-10)**2.
         aee=a0+a1*(enu-10.)
         pnn=pdd*(1.+0.5*aee)/(1.-0.5*aee)
         ccsnodd=pdd*ccd(i)
         ccsnonn=pnn*ccn(i)
         eselsnodd=pdd*eseld(i)
         eselsnonn=pnn*eseln(i)
         esmusnodd=(1.-pdd)*esmud(i)
         esmusnonn=(1.-pnn)*esmun(i)
         chi2sno=chi2sno+(ccsnodd-ccsnod(i))**2.+(ccsnonn-ccsnon(i))**2.
     %        +(eselsnodd-eselsnod(i))**2.+(eselsnonn-eselsnon(i))**2.
     %        +(esmusnodd-esmusnod(i))**2.+(esmusnonn-esmusnon(i))**2.
      enddo

      fval = chi2sno

      return
      end

