      implicit real*4 (a-h,o-z)
c..............................................................
c................... sun profiles .............................
c..............................................................
      dimension rr_b(1268),rho_b(1268)
      dimension r0(501)
      dimension fl1r0(501),fl2r0(501),fl3r0(501),fl4r0(501)
      dimension fl5r0(501),fl6r0(501),fl7r0(501),fl8r0(501)
      common /sun1/rr_b,rho_b,r0,fl1r0,fl2r0,fl3r0,fl4r0
      common /sun2/fl5r0,fl6r0,fl7r0,fl8r0
c..............................................................
c.................. probabilities .............................
c..............................................................
c.. regeneration
      dimension preg(9,46,339),a2e(10000),na2e(2),x2e(2)
      common /preg/a2e,na2e,preg
c.. solar + earth
      parameter (nenu=401)
      dimension dm4e(nenu),enu(nenu)
      dimension pday(8),pnight(8,9)      
c.. different experiments
      dimension probsk(8,nenu,2)
      real*4 probsnodd(nenu),probsnonn(nenu)
      dimension probbe7dd(nenu),probbe7nn(nenu)
      dimension probbxnhep(nenu),probbxnb8(nenu)
      dimension probgach(nenu,8)
            
c..  electron density inside the Sun
      open (11,file='bs2005op.dat',status='old')
      do iji=1,1268
      read(11,*)rr_b(iji),rho_b(iji),rhon
      enddo
      close(11)
c..  production point distribution
      open (12,file='bs05opflux2.dat',status='old')
      do j=1,501
         read(12,*)r0(j),fl1r0(j),fl2r0(j),fl3r0(j)
     1   ,fl4r0(j),fl5r0(j),fl6r0(j),fl7r0(j),fl8r0(j)
      enddo
      close(12)

      open(15,file='probtest.dat',status='unknown')
      
c..............................................................
c............... regeneration probabilities ...................
c..............................................................      
      call pregeneration(a2e,na2e,preg)

c..............................................................
c................... survival probabilities ...................
c..............................................................
      s2_13=0.023
      tg=0.44
c.. probabilities in function of Deltam/4E
      do j=1,nenu
         dm=10.**(-13.+4.*float(j-1)/float(nenu-1))
         dm4e(j)=log10(dm)
         call probsunearth(dm,tg,s2_13,pday,pnight)
         probsk(1,j,1)=pday(5)  !8B
         probsk(1,j,2)=pday(3)  !hep
         do k=1,7
            probsk(k+1,j,1)=pnight(5,k) !8B 
            probsk(k+1,j,2)=pnight(3,k) !hep
         enddo
         probsnodd(j)=pday(5)
         probsnonn(j)=pnight(5,8)
         probbe7dd(j)=pday(4)
         probbe7nn(j)=pnight(4,9)
         probbxnhep(j)=pnight(3,9)
         probbxnb8(j)=pnight(5,9)
         do k=1,8
            probgach(j,k)=0.5*(pday(k)+pnight(k,9))
         enddo
      enddo

      deltam=8.e-5
      do j=1,nenu
         enutst=deltam/10.**dm4e(j)/4.e6
         write(15,*)enutst,(probsk(k,j,1),k=1,7)
      enddo
      
      end


