      implicit real*4 (a-h,o-z)

      open(11,file='chi2all.dat',status='old')

      chi2min=1.e5
      do i=1,10000000
c         read(11,*,end=11)a,b,c,chi2
         read(11,*,end=11)a,b,c,d,e,chi2
         if(chi2.lt.chi2min)then
            chi2min=chi2
            amin=a
            bmin=b
         endif
      enddo

 11   continue
      write(*,*)chi2min,amin,bmin
      end
