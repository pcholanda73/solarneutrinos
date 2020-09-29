      subroutine chi2brxhe(rbrx,chi2brx)
      implicit real*4 (a-h,o-z)
      dimension databrx(6),errsystcorr(6),rbrx(6)
      dimension errstat(6),bkg(6)
      dimension sigma(6,6),ir(6)
c      data databrx/82.025146,56.483402,30.913017,11.866378,
c     %             2.3526304,0.10939147/
      data bkg/34.2,3.1,2.9,1.8,0.,0./
      data databrx/29.,26.,14.,5.,1.,0./
      data errsystcorr/-1.93233769E-02,1.90292820E-02,8.87698382E-02,
     %                0.22328469,0.52661246,1.3736196/ !from brx-bins.f


c.. erro na escala de energia
      do i=1,6
      do j=1,6
         sigma(i,j)=errsystcorr(i)*errsystcorr(j)*rbrx(i)*rbrx(j)
      enddo
      enddo

c.. erro no volume fiducial
      do i=1,6
      do j=1,6
         sigma(i,j)=sigma(i,j)+(0.038)**2.*rbrx(i)*rbrx(j)
      enddo
      enddo

c.. erro estatistico (incluindo background)
      do i=1,6
         sigma(i,i)=sigma(i,i)+databrx(i)+bkg(i)+1.
      enddo

      call rinv(6,sigma,6,ir,ifail)

      chi2brx=0.
      do j1=1,6
      do j2=1,6
	soma = (rbrx(j1)-databrx(j1))*(rbrx(j2)-databrx(j2))
     %         *sigma(j1,j2)
        chi2brx=chi2brx+soma
      enddo
      enddo

c      write(*,*)rbrx
c      write(*,*)databrx

      return
      end
