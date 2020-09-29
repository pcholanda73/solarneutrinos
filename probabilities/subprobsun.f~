      subroutine probsunearth(dm,tg,s2_13,pday,pnight)
      implicit real*4 (a-h,o-z)
c.......................................      
c............. sun profile .............
      dimension rr_b(1268),rho_b(1268)
      dimension r0(501),fmed(8),pmed(8)
      dimension fl1r0(501),fl2r0(501),fl3r0(501),fl4r0(501)
      dimension fl5r0(501),fl6r0(501),fl7r0(501),fl8r0(501)
      common /sun1/rr_b,rho_b,r0,fl1r0,fl2r0,fl3r0,fl4r0
      common /sun2/fl5r0,fl6r0,fl7r0,fl8r0
c.......................................      
      dimension pday(8),p2e(9),pnight(8,9)
      dimension a2e(10000),na2e(2),x2e(2),preg(9,46,339),pregtmp(46,339)
      common /preg/a2e,na2e,preg


      c_2t=(1.-tg)/(1.+tg)
      s_2t=sqrt(1.-c_2t**2.)
      s2_t=tg/(1.+tg)


c........................................................      
c............. survival probability in Sun ..............
c........................................................      
      do i=1,8
         fmed(i)=0.
         pmed(i)=0.
      enddo

c............. averaging on production point ............
      do i=1,501
      rr=r0(i)
      rhosol=divdif(rho_b,rr_b,1268,rr,2)
      rhosol=(1.-s2_13)*10.**(rhosol)

      aa=rhosol-dm*c_2t
      bb=dm*s_2t
      s_2tm = bb/(sqrt(aa**2.+bb**2.))
      c_2tm = -aa/(sqrt(aa**2.+bb**2.))

c      pee=s2_13**2.+(1.-s2_13)**2.*(0.5+0.5*c_2tm*c_2t)
      pee=0.5+0.5*c_2tm*c_2t
      
c............. 8 different production regions ...........
c............. pp, pep, hep, 7Be, 8B, N, O, F ...........
      ff=divdif(fl1r0,r0,501,abs(rr),1)
      fmed(1)=fmed(1)+ff
      pmed(1)=pmed(1)+ff*pee
      ff=divdif(fl2r0,r0,501,abs(rr),1)
      fmed(2)=fmed(2)+ff
      pmed(2)=pmed(2)+ff*pee
      ff=divdif(fl3r0,r0,501,abs(rr),1)
      fmed(3)=fmed(3)+ff
      pmed(3)=pmed(3)+ff*pee
      ff=divdif(fl4r0,r0,501,abs(rr),1)
      fmed(4)=fmed(4)+ff
      pmed(4)=pmed(4)+ff*pee
      ff=divdif(fl5r0,r0,501,abs(rr),1)
      fmed(5)=fmed(5)+ff
      pmed(5)=pmed(5)+ff*pee
      ff=divdif(fl6r0,r0,501,abs(rr),1)
      fmed(6)=fmed(6)+ff
      pmed(6)=pmed(6)+ff*pee
      ff=divdif(fl7r0,r0,501,abs(rr),1)
      fmed(7)=fmed(7)+ff
      pmed(7)=pmed(7)+ff*pee
      ff=divdif(fl8r0,r0,501,abs(rr),1)
      fmed(8)=fmed(8)+ff
      pmed(8)=pmed(8)+ff*pee

      enddo

      do kk=1,8
         pmed(kk)=pmed(kk)/fmed(kk)
         pday(kk)=s2_13**2.+(1.-s2_13)**2.*pmed(kk)
      enddo


c........................................................      
c... P2e transition on Earth (sin^2\theta in vacuum) ....
c........................................................      
      x2e(1)=tg
      x2e(2)=log10(dm)

      
c... 9 different averaged earth profile
c..          1 to 6 + 7  -> 6 zenith bins at SK + averaged night
c..          8           -> averaged night at SNO
c..          9           -> averaged night at GALLEX/GNO and SAGE
      do i=1,9
         do j=1,na2e(1)
         do k=1,na2e(2)
            pregtmp(j,k)=preg(i,j,k)
         enddo
         enddo
         p2e(i)=fint(2,x2e,na2e,a2e,pregtmp)
         if(x2e(2).gt.a2e(na2e(1)+na2e(2)))p2e(i)=s2_t
      enddo
      
      do n=1,8
      do k=1,9
         pnight(n,k)=s2_13**2.+
     %  (1.-s2_13)**2.*(pmed(n)-s2_t+p2e(k)*(1.-2.*pmed(n)))/c_2t  
      enddo
      enddo
      
      return
      end



