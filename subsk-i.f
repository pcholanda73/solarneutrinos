      subroutine ski(nenu,enu,probsk,rsk1b8,rsk1hep)
      implicit real*4 (a-h,o-z)
      dimension rsk1b8(46),rsk1hep(46)
c.
      dimension enu(nenu)
c      dimension E_B8(800),B8(800)
      dimension E_B8(155),B8(155)
      dimension E_hep(1000),hep(1000)
      common /B8/ E_B8,B8
      common /hep/ E_hep,hep
c.
      dimension probsk(8,nenu,2)
      dimension pday(nenu),pnightb(nenu)
      dimension pmanto1(nenu),pmanto2(nenu),pmanto3(nenu)
      dimension pmanto4(nenu),pmanto5(nenu),pnucleo(nenu)
      dimension pdayh(nenu),pnighth(nenu)
      dimension pmanto1h(nenu),pmanto2h(nenu),pmanto3h(nenu)
      dimension pmanto4h(nenu),pmanto5h(nenu),pnucleoh(nenu)
      dimension flboro(0:1000),flhep(0:1000)
      dimension pd(0:1000),pm1(0:1000),pm2(0:1000),pm3(0:1000)
      dimension pm4(0:1000),pm5(0:1000),pnu(0:1000),pn(0:1000)
      dimension pdh(0:1000),pm1h(0:1000),pm2h(0:1000),pm3h(0:1000)
      dimension pm4h(0:1000),pm5h(0:1000),pnuh(0:1000),pnh(0:1000)
c.
      dimension aintday(0:1000),aintnight(0:1000)
      dimension aintmanto1(0:1000),aintmanto2(0:1000)
      dimension aintmanto3(0:1000),aintmanto4(0:1000)
      dimension aintmanto5(0:1000),aintnucleo(0:1000)
      dimension aintdayh(0:1000),aintnighth(0:1000)
      dimension aintmanto1h(0:1000),aintmanto2h(0:1000)
      dimension aintmanto3h(0:1000),aintmanto4h(0:1000)
      dimension aintmanto5h(0:1000),aintnucleoh(0:1000)
c.
c      dimension espec(101),specel(8,101),specmu(8,101)
c      common /cssk/enusk1,specel,specmu
      dimension enusk1(181),cselsk1(8,181),csmusk1(8,181)
      common /cssk1/enusk1,cselsk1,csmusk1
      dimension spectrumel(181),spectrummu(181)
      dimension essm(0:8),specssm(8)
c.
      dimension specday(8),specnight(8)
      dimension specmanto1(8),specmanto2(8),specmanto3(8)
      dimension specmanto4(8),specmanto5(8),specnucleo(8)
      dimension specdayh(8),specnighth(8)
      dimension specmanto1h(8),specmanto2h(8),specmanto3h(8)
      dimension specmanto4h(8),specmanto5h(8),specnucleoh(8)
c.
c      dimension rth(8,46)
c.
      dimension ssmb8(8),ssmhep(8)!,ssmtot(8)
      data ssmb8/1.82866090E-04,3.19422223E-04,3.57678422E-04,
     %   2.23680123E-04,1.43400146E-04,4.43975077E-05,1.24992976E-05,
     %   6.03666990E-07/
      data ssmhep/3.17977083E-07,6.03931426E-07,8.04221202E-07,
     %   6.63202229E-07,6.54652581E-07,4.10399679E-07,3.55849011E-07,
     %   9.48806331E-08/

c      data ssmtot/0.000957152690,0.001653459160,0.001832611040,
c     %               0.001136178270,0.000720176497,0.000219382971,
c     %               4.46319173E-05,2.72680745E-06/ 
c      data specssm/0.00107802858,0.00186217960,0.00206371187,
c     %             0.00127916876,0.00081041234,0.00024652443,
c     %             4.9901722E-05,2.9527821E-06/

c      data specssm/180.059,314.568,352.371,220.522,141.607,44.0513,
c     %             12.6436,0.68882/
c      dimension specssmb8(8),specssmhep(8)
c      data specssmb8/179.739,313.960,351.562,219.855,140.948,
c     %              43.6383,12.2855,0.59334/
c      data specssmhep/0.3200,0.6078,0.8093,0.6674,0.6588,0.4130,
c     %              0.3581,9.5483E-02/
      dimension flux(8)
      common /nuflux/flux
c      dimension fluxo(8)
c      common /ssm/nssm,fluxo
      dimension specssmb8(8),specssmhep(8)!,specssm(46)
c      dimension specssmb82(46),specssmhep2(46)
      data specssmb8/182.9,320.,358.,223.8,143.4,44.4,9.33,0.611/
      data specssmhep/0.309,0.603,0.799,0.653,0.631,0.379,0.211,0.068/

      fb=1.e4*flux(5)
      fh=1.e4*flux(3)
c      fb=5.79
c      fh=7.88e-3
      nssm=0

c..
      do i=1,nenu
         pday(i)=probsk(1,i,1)
         pmanto1(i) =probsk(2,i,1)
         pmanto2(i) =probsk(3,i,1)
         pmanto3(i) =probsk(4,i,1)
         pmanto4(i) =probsk(5,i,1)
         pmanto5(i) =probsk(6,i,1)
         pnucleo(i) =probsk(7,i,1)
         pnightb(i)  =probsk(8,i,1)
         pdayh(i)=probsk(1,i,2)
         pmanto1h(i) =probsk(2,i,2)
         pmanto2h(i) =probsk(3,i,2)
         pmanto3h(i) =probsk(4,i,2)
         pmanto4h(i) =probsk(5,i,2)
         pmanto5h(i) =probsk(6,i,2)
         pnucleoh(i) =probsk(7,i,2)
         pnighth(i)  =probsk(8,i,2)
      enddo

      nint=1000
      do i=0,nint
         energia=2.2+(18.784-2.2)*float(i)/float(nint)
c         flboro(i)=divdif(b8,e_b8,800,energia,2)
         flboro(i)=divdif(b8,e_b8,155,energia,2)
         flhep(i)=divdif(hep,e_hep,1000,energia,2)
         if(energia.gt.16.34)flboro(i)=0.
         if(energia.gt.18.784)flhep(i)=0.
         pd(i)=divdif(pday,enu,nenu,1.e6*energia,2)
         pm1(i)=divdif(pmanto1,enu,nenu,1.e6*energia,2)
         pm2(i)=divdif(pmanto2,enu,nenu,1.e6*energia,2)
         pm3(i)=divdif(pmanto3,enu,nenu,1.e6*energia,2)
         pm4(i)=divdif(pmanto4,enu,nenu,1.e6*energia,2)
         pm5(i)=divdif(pmanto5,enu,nenu,1.e6*energia,2)
         pnu(i)=divdif(pnucleo,enu,nenu,1.e6*energia,2)
         pn(i)=divdif(pnightb,enu,nenu,1.e6*energia,2)
         pdh(i)=divdif(pdayh,enu,nenu,1.e6*energia,2)
         pm1h(i)=divdif(pmanto1h,enu,nenu,1.e6*energia,2)
         pm2h(i)=divdif(pmanto2h,enu,nenu,1.e6*energia,2)
         pm3h(i)=divdif(pmanto3h,enu,nenu,1.e6*energia,2)
         pm4h(i)=divdif(pmanto4h,enu,nenu,1.e6*energia,2)
         pm5h(i)=divdif(pmanto5h,enu,nenu,1.e6*energia,2)
         pnuh(i)=divdif(pnucleoh,enu,nenu,1.e6*energia,2)
         pnh(i)=divdif(pnighth,enu,nenu,1.e6*energia,2)
         if(nssm.eq.1)pd(i)=1.
         if(nssm.eq.1)pm1(i)=1.
         if(nssm.eq.1)pm2(i)=1.
         if(nssm.eq.1)pm3(i)=1.
         if(nssm.eq.1)pm4(i)=1.
         if(nssm.eq.1)pm5(i)=1.
         if(nssm.eq.1)pnu(i)=1.
         if(nssm.eq.1)pn(i)=1.
         if(nssm.eq.1)pdh(i)=1.
         if(nssm.eq.1)pm1h(i)=1.
         if(nssm.eq.1)pm2h(i)=1.
         if(nssm.eq.1)pm3h(i)=1.
         if(nssm.eq.1)pm4h(i)=1.
         if(nssm.eq.1)pm5h(i)=1.
         if(nssm.eq.1)pnuh(i)=1.
         if(nssm.eq.1)pnh(i)=1.
      enddo
c................. loop nos diferentes bins de energia
      do ik=1,8

      do il=1,181
      spectrumel(il)=cselsk1(ik,il)
      spectrummu(il)=csmusk1(ik,il)
      enddo
c..................................................................
      do i=0,nint
      energia=2.2+(18.784-2.2)*float(i)/float(nint)
      rspecel=divdif(spectrumel,enusk1,181,energia,2)
      rspecmu=divdif(spectrummu,enusk1,181,energia,2)
      if(rspecel.lt.0.)rspecel=0.
      if(rspecmu.lt.0.)rspecmu=0.
      aintday(i)=fb*flboro(i)*(rspecel*pd(i)+rspecmu*(1.-pd(i)))
      aintmanto1(i)=fb*flboro(i)*(rspecel*pm1(i)+rspecmu*(1.-pm1(i)))
      aintmanto2(i)=fb*flboro(i)*(rspecel*pm2(i)+rspecmu*(1.-pm2(i)))
      aintmanto3(i)=fb*flboro(i)*(rspecel*pm3(i)+rspecmu*(1.-pm3(i)))
      aintmanto4(i)=fb*flboro(i)*(rspecel*pm4(i)+rspecmu*(1.-pm4(i)))
      aintmanto5(i)=fb*flboro(i)*(rspecel*pm5(i)+rspecmu*(1.-pm5(i)))
      aintnucleo(i)=fb*flboro(i)*(rspecel*pnu(i)+rspecmu*(1.-pnu(i)))
      aintnight(i)=fb*flboro(i)*(rspecel*pn(i)+rspecmu*(1.-pn(i)))
      aintdayh(i)=fh*flhep(i)*(rspecel*pdh(i)+rspecmu*(1.-pdh(i)))
      aintmanto1h(i)=fh*flhep(i)*(rspecel*pm1h(i)+rspecmu*(1.-pm1h(i)))
      aintmanto2h(i)=fh*flhep(i)*(rspecel*pm2h(i)+rspecmu*(1.-pm2h(i)))
      aintmanto3h(i)=fh*flhep(i)*(rspecel*pm3h(i)+rspecmu*(1.-pm3h(i)))
      aintmanto4h(i)=fh*flhep(i)*(rspecel*pm4h(i)+rspecmu*(1.-pm4h(i)))
      aintmanto5h(i)=fh*flhep(i)*(rspecel*pm5h(i)+rspecmu*(1.-pm5h(i)))
      aintnucleoh(i)=fh*flhep(i)*(rspecel*pnuh(i)+rspecmu*(1.-pnuh(i)))
      aintnighth(i)=fh*flhep(i)*(rspecel*pnh(i)+rspecmu*(1.-pnh(i)))
      enddo
c..................................................................

      specday(ik) = simps(aintday,2.2,18.784,nint)
      specmanto1(ik) = simps(aintmanto1,2.2,18.784,nint)
      specmanto2(ik) = simps(aintmanto2,2.2,18.784,nint)
      specmanto3(ik) = simps(aintmanto3,2.2,18.784,nint)
      specmanto4(ik) = simps(aintmanto4,2.2,18.784,nint)
      specmanto5(ik) = simps(aintmanto5,2.2,18.784,nint)
      specnucleo(ik) = simps(aintnucleo,2.2,18.784,nint)
      specnight(ik) = simps(aintnight,2.2,18.784,nint)
      specdayh(ik) = simps(aintdayh,2.2,18.784,nint)
      specmanto1h(ik) = simps(aintmanto1h,2.2,18.784,nint)
      specmanto2h(ik) = simps(aintmanto2h,2.2,18.784,nint)
      specmanto3h(ik) = simps(aintmanto3h,2.2,18.784,nint)
      specmanto4h(ik) = simps(aintmanto4h,2.2,18.784,nint)
      specmanto5h(ik) = simps(aintmanto5h,2.2,18.784,nint)
      specnucleoh(ik) = simps(aintnucleoh,2.2,18.784,nint)
      specnighth(ik) = simps(aintnighth,2.2,18.784,nint)
      enddo
c..................

      rsk1hep(1)=specdayh(1)/ssmhep(1)*specssmhep(1)
      rsk1b8(1)=specday(1)/ssmb8(1)*specssmb8(1)
      rsk1hep(2)=specnighth(1)/ssmhep(1)*specssmhep(1)
      rsk1b8(2)=specnight(1)/ssmb8(1)*specssmb8(1)
      do j=2,7
	ii=(j-2)*7
	rsk1hep(ii+3)=specdayh(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+3)=specday(j)/ssmb8(j)*specssmb8(j)
	rsk1hep(ii+4)=specmanto1h(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+4)=specmanto1(j)/ssmb8(j)*specssmb8(j)
	rsk1hep(ii+5)=specmanto2h(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+5)=specmanto2(j)/ssmb8(j)*specssmb8(j)
	rsk1hep(ii+6)=specmanto3h(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+6)=specmanto3(j)/ssmb8(j)*specssmb8(j)
	rsk1hep(ii+7)=specmanto4h(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+7)=specmanto4(j)/ssmb8(j)*specssmb8(j)
	rsk1hep(ii+8)=specmanto5h(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+8)=specmanto5(j)/ssmb8(j)*specssmb8(j)
	rsk1hep(ii+9)=specnucleoh(j)/ssmhep(j)*specssmhep(j)
	rsk1b8(ii+9)=specnucleo(j)/ssmb8(j)*specssmb8(j)
      enddo
      rsk1hep(45)=specdayh(8)/ssmhep(8)*specssmhep(8)
      rsk1b8(45)=specday(8)/ssmb8(8)*specssmb8(8)
      rsk1hep(46)=specnighth(8)/ssmhep(8)*specssmhep(8)
      rsk1b8(46)=specnight(8)/ssmb8(8)*specssmb8(8)

      if(nssm.eq.1)then
         write(*,*)rsk1b8(1),rsk1hep(1),rsk1b8(1)+rsk1hep(1)  
         write(*,*)rsk1b8(3),rsk1hep(3),rsk1b8(3)+rsk1hep(3)  
         write(*,*)rsk1b8(10),rsk1hep(10),rsk1b8(10)+rsk1hep(10)
         write(*,*)rsk1b8(17),rsk1hep(17),rsk1b8(17)+rsk1hep(17)
         write(*,*)rsk1b8(24),rsk1hep(24),rsk1b8(24)+rsk1hep(24)
         write(*,*)rsk1b8(31),rsk1hep(31),rsk1b8(31)+rsk1hep(31)
         write(*,*)rsk1b8(38),rsk1hep(38),rsk1b8(38)+rsk1hep(38)
         write(*,*)rsk1b8(45),rsk1hep(45),rsk1b8(45)+rsk1hep(45)
         stop
      endif

      return
      end



