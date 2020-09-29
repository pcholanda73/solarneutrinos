      subroutine chi2sno(rsno,chi2sno3f,ff)
      implicit real*4 (a-h,o-z)
c      real*8 pval(10)
      dimension corr(6,6),rsno(5)
      data corr/+1.000,-0.723,+0.302,-0.168,+0.028,-0.012,
     %          -0.723,+1.000,-0.299,-0.366,-0.376,+0.129,
     %          +0.302,-0.299,+1.000,-0.206,+0.219,-0.677,
     %          -0.168,-0.366,-0.206,+1.000,+0.008,-0.035,
     %          +0.028,-0.376,+0.219,+0.008,+1.000,-0.297,
     %          -0.012,+0.129,-0.677,-0.035,-0.297,+1.000/
      dimension ppar(6),par(6),sigstat(6),sigsyst(6)
      data par/5.25,0.317,0.0039,-0.001,0.046,-0.016/
      data sigstat/0.16,0.016,0.0066,0.0029,0.031,0.025/
      data sigsyst/0.12,0.009,0.0045,0.0015,0.014,0.011/
c      dimension corr(5,5)
c      data corr/
c     %          +1.000,-0.299,-0.366,-0.376,+0.129,
c     %          -0.299,+1.000,-0.206,+0.219,-0.677,
c     %          -0.366,-0.206,+1.000,+0.008,-0.035,
c     %          -0.376,+0.219,+0.008,+1.000,-0.297,
c     %          +0.129,-0.677,-0.035,-0.297,+1.000/
c      dimension ppar(5),par(5),sigstat(5),sigsyst(5)
c      data par/0.317,0.0039,-0.001,0.046,-0.016/
c      data sigstat/0.016,0.0066,0.0029,0.031,0.025/
c      data sigsyst/0.009,0.0045,0.0015,0.014,0.011/
      dimension sig(6,6),isig(6)
c............................................................
c............................................................
      do i=1,6
         do j=1,6
            sig(i,j)=sqrt(sigstat(i)**2.+sigsyst(i)**2.)*
     %     corr(i,j)*sqrt(sigstat(j)**2.+sigsyst(j)**2.)
         enddo
      enddo
      call rinv(6,sig,6,isig,ifail)
      do k=1,5
         ppar(k+1)=rsno(k)
      enddo

      ppar(1)=ff*5.25
      
      chi2sno3f=0.
      do i=1,6
         do j=1,6
            ppi=ppar(i)
            ppj=ppar(j)
            chi2sno3f=chi2sno3f+(ppi-par(i))*(ppj-par(j))*sig(i,j)
         enddo
      enddo

      return
      end
      
