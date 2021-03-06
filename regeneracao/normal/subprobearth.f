      subroutine probearth(dm,tg)
      implicit real*4 (a-h,o-z)
      dimension c_t_(100)
      dimension tlsk1_(100)
      dimension tlsg1_(100),tlsg2a_(100),tlsg2b_(100)
      dimension tlsno_(100)
      common /timeexp/c_t_,tlsk1_,tlsg1_,tlsg2a_,tlsg2b_,tlsno_
      dimension P_night(6),pnight(7,5)
      dimension cphi(10000),pphi(10000)
      common /pnight/pnight

      c_2t=(1.-tg)/(1.+tg)
      s_2t=sqrt(1.-c_2t**2.)
      s2_t=tg/(1.+tg)

c..   calculating P(2e) runing through zenith angle, with
c..         840 points for cos\theta=[0:0.840]
c..         675 points for cos\theta=[0.840,0.975]
c..         does not need to go higher due to latitude of detectors
      ncore=5*135
      nmantle=840
      np=nmantle+ncore
      do icos=1,np
         cosphi=float(icos)/1000.
         if(icos.gt.nmantle)then
            cosphi=.84+.135*float(icos-nmantle)/float(ncore)
         endif
         cphi(icos)=cosphi
         pphi(icos)=p2e(tg,dm,cosphi)
      enddo

c..   averaging over zenith angle
c.. SK:
      call  amm(pphi,cphi,ncore,np,c_t_,tlsk1_,P_night,Pn_medio)
      do k=1,6
         pnight(k,1)=P_night(k)
      enddo
      pnight(7,1)=Pn_medio

c.. SNO
      call  amm(pphi,cphi,ncore,np,c_t_,tlsno_,P_night,Pn_medio)
      do k=1,6
         pnight(k,2)=P_night(k)
      enddo
      pnight(7,2)=Pn_medio

c.. SAGE
      call  amm(pphi,cphi,ncore,np,c_t_,tlsg1_,P_night,Pn_medio)
      do k=1,6
         pnight(k,3)=P_night(k)
      enddo
      pnight(7,3)=Pn_medio

c.. GALLEX (summer and winter) 
      call  amm(pphi,cphi,ncore,np,c_t_,tlsg2a_,P_night,Pn_medio)
      do k=1,6
         pnight(k,4)=P_night(k)
      enddo
      pnight(7,4)=Pn_medio

      call  amm(pphi,cphi,ncore,np,c_t_,tlsg2b_,P_night,Pn_medio)
      do k=1,6
         pnight(k,5)=P_night(k)
      enddo
      pnight(7,5)=Pn_medio

      return
      end




c**********************************************************
c**********************************************************
      subroutine amm(pphi_,zenith_,ncore,np,c_t_,timelife_,
     %               P_night,Pn_medio)
      implicit real*4 (a-h,o-z)
      dimension P_night(6)
      dimension zenith_(np),pphi_(np)
      dimension c_t_(100),timelife_(100)
      dimension r(7)
      real*4 pnight1(0:10000),pnight2(0:10000),pnight3(0:10000)
      real*4 pnight4(0:10000),pnight5(0:10000),pnight6(0:10000)
      real*4 pnight7(0:10000)
      real*4 tt1(0:10000),tt2(0:10000),tt3(0:10000)
      real*4 tt4(0:10000),tt5(0:10000),tt6(0:10000)
      real*4 tt7(0:10000)

      cmax1 =.16 
      cmax2 =.33
      cmax3 =.50
      cmax4 =.67
      cmax5 =.84
      cmax6 = 1.
      cmax7 = 1.

      do i=1,6
         P_night(i) = 0.0
      enddo

c..   dividing the night in 6 regions
c..   and each region in <nkb> points in cos(theta)
      delcos = (0.975-0.840)/float(ncore)
      if(delcos.ge.1.e-2)then   !100 points em lt#.dat
         nkb = 100
      else
         nkb = 1./delcos 
      endif

c--------------------------- LOOP --------------------------
      do kb=0,nkb
         c_alpha = float(kb)/float(nkb)
         if((c_alpha.ge.0.0).and.(c_alpha.le.cmax1))then
            nnight=1
            night1=kb
         elseif((c_alpha.gt.cmax1).and.(c_alpha.le.cmax2))then
            nnight=2
            night2=kb
         elseif((c_alpha.gt.cmax2).and.(c_alpha.le.cmax3))then
            nnight=3
            night3=kb
         elseif((c_alpha.gt.cmax3).and.(c_alpha.le.cmax4))then
            nnight=4
            night4=kb
         elseif((c_alpha.gt.cmax4).and.(c_alpha.le.cmax5))then
            nnight=5
            night5=kb
         elseif((c_alpha.gt.cmax5).and.(c_alpha.le.1.0))then
            nnight=6
            night6=kb
         endif


         texposure = divdif(timelife_,c_t_,100,c_alpha,1)
         p2e = divdif(pphi_,zenith_,np,c_alpha,1)
         pnight7(kb)=texposure*p2e
         tt7(kb)=texposure

         if(nnight.eq.1)then
            pnight1(kb)=texposure*p2e
            tt1(kb)=texposure
            pnight2(0)=pnight1(kb)
            tt2(0)=tt1(kb)
         elseif(nnight.eq.2)then
            pnight2(kb-night1)=texposure*p2e
            tt2(kb-night1)=texposure
            pnight3(0)=pnight2(kb-night1)
            tt3(0)=tt2(kb-night1)
         elseif(nnight.eq.3)then
            pnight3(kb-night2)=texposure*p2e
            tt3(kb-night2)=texposure
            pnight4(0)=pnight3(kb-night2)
            tt4(0)=tt3(kb-night2)
         elseif(nnight.eq.4)then
            pnight4(kb-night3)=texposure*p2e
            tt4(kb-night3)=texposure
            pnight5(0)=pnight4(kb-night3)
            tt5(0)=tt4(kb-night3)
         elseif(nnight.eq.5)then
            pnight5(kb-night4)=texposure*p2e
            tt5(kb-night4)=texposure
            pnight6(0)=pnight5(kb-night4)
            tt6(0)=tt5(kb-night4)
         elseif(nnight.eq.6)then
            pnight6(kb-night5)=texposure*p2e
            tt6(kb-night5)=texposure
         endif

      enddo

c---------------------------------------------------------------
      nn1 = night1
      nn2 = night2-night1
      nn3 = night3-night2
      nn4 = night4-night3
      nn5 = night5-night4
      nn6 = night6-night5
      nn7 = nkb

      nt1=2*(nn1/2)
      nt2=2*(nn2/2)
      nt3=2*(nn3/2)
      nt4=2*(nn4/2)
      nt5=2*(nn5/2)
      nt6=2*(nn6/2)
      nt7=2*(nn7/2)
      if(nt1.lt.nn1)then
         nn1=nn1-1
         cmax1 = float(night1-1)/float(nkb)
      endif
      if(nt2.lt.nn2)then
         nn2=nn2-1
         cmax2 = float(night2-1)/float(nkb)
      endif
      if(nt3.lt.nn3)then
         nn3=nn3-1
         cmax3 = float(night3-1)/float(nkb)
      endif
      if(nt4.lt.nn4)then
         nn4=nn4-1
         cmax4 = float(night4-1)/float(nkb)
      endif
      if(nt5.lt.nn5)then
         nn5=nn5-1
         cmax5 = float(night5-1)/float(nkb)
      endif
      if(nt6.lt.nn6)then
         nn6=nn6-1
         cmax6 = float(night6-1)/float(nkb)
      endif
      if(nt7.lt.nn7)then
         nn7=nn7-1
         cmax7 = float(nkb-1)/float(nkb)
      endif

      cmin=0.0
      cmax=cmax1
      rates1 = simps(pnight1,cmin,cmax,nn1)
      texp1 = simps(tt1,cmin,cmax,nn1)
      r(1) = rates1/texp1

      cmin=cmax1
      cmax=cmax2
      rates2 = simps(pnight2,cmin,cmax,nn2)
      texp2 = simps(tt2,cmin,cmax,nn2)
      r(2) = rates2/texp2

      cmin=cmax2
      cmax=cmax3
      rates3 = simps(pnight3,cmin,cmax,nn3)
      texp3 = simps(tt3,cmin,cmax,nn3)
      r(3) = rates3/texp3

      cmin=cmax3
      cmax=cmax4
      rates4 = simps(pnight4,cmin,cmax,nn4)
      texp4 = simps(tt4,cmin,cmax,nn4)
      r(4) = rates4/texp4

      cmin=cmax4
      cmax=cmax5
      rates5 = simps(pnight5,cmin,cmax,nn5)
      texp5 = simps(tt5,cmin,cmax,nn5)
      r(5) = rates5/texp5

      cmin=cmax5
      cmax=cmax6
      rates6 = simps(pnight6,cmin,cmax,nn6)
      texp6 = simps(tt6,cmin,cmax,nn6)
      r(6) = rates6/texp6

      cmin=0.0
      cmax=cmax7
      rates7 = simps(pnight7,cmin,cmax,nn7)
      texp7 = simps(tt7,cmin,cmax,nn7)
      r(7) = rates7/texp7

      do n=1,6
         P_night(n)=r(n)
      enddo
      Pn_medio = r(7)
      
      return
      end


c**********************************************************
c**********************************************************
c.. calculate P2e using analytical results for constant matter layers
c.. PREM model for Earth density
      function p2e(tg,dm,cosphi)
      implicit real *4 (a-h, o-z)
      dimension rr_t(55),rho_t(55)
      dimension rho_(20000),rr_(20000)
      data rr_t/0.0000,0.0313,0.0627,0.0940,0.1254,0.1567,0.1881,0.1915,
     %   0.1916,0.2194,0.2508,0.2822,0.3135,0.3449,0.3762,0.4076,0.4389,
     %   0.4703,0.5017,0.5330,0.5456,0.5457,0.5487,0.5691,0.5957,0.6271,
     %   0.6584,0.6898,0.7212,0.7525,0.7839,0.8152,0.8466,0.8780,0.8938,
     %   0.8939,0.9048,0.9093,0.9250,0.9361,0.9362,0.9407,0.9563,0.9643,
     %   0.9644,0.9720,0.9863,0.9877,0.9950,0.9951,0.9965,0.9966,0.9984,
     %   0.9985,0.9988/
      data rho_t/6.0840,6.0793,6.0654,6.0468,6.0189,5.9817,5.9353,5.9306
     %   ,5.6564,5.6099,5.5541,5.4891,5.4147,5.3357,5.2474,5.1498,5.0429
     %   ,4.9267,4.8012,4.6617,4.6013,2.7492,2.7443,2.7097,2.6702,2.6209
     %   ,2.5715,2.5222,2.4679,2.4185,2.3642,2.3099,2.2507,2.1915,2.1619
     %   ,1.9694,1.9644,1.9447,1.8805,1.8361,1.7472,1.7423,1.7127,1.6979
     %   ,1.6584,1.6584,1.6633,1.6683,1.6683,1.4313,1.4313,1.2833,1.2833
     %   ,0.5034,0.5034/
      complex*8 ucte(2,2),uf(2,2),ufim(2,2),umass(2),umix(2,2)
      complex*8 uevol(2,2)

c-----------------------------------------------------------------
c......... some constants
      pi = acos(-1.0)
      Gf = 1.16632e-23
      R_earth = 6378.14e3

      s2_t = tg/(1.+tg)
      c2_t = 1./(1.+tg)
      s_t = sqrt(s2_t)
      c_t = sqrt(c2_t)
      s_2t = 2.*sqrt(s2_t*c2_t)
      c_2t = c2_t -s2_t

      phi = acos(cosphi)

      rhomax=0.5*sqrt(2.0)*Gf*6.02e23*rho_t(55)/((5.068e4)**3.0)
      rhomantle=0.5*sqrt(2.0)*Gf*6.02e23*rho_t(3)/((5.068e4)**3.0)
      dmtil=sqrt((dm*c_2t-rhomax)**2.+(dm*s_2t)**2.)
      dmtil=max(dmtil,dm)
      alambda=1./(5.068e6*dmtil)
      dr=(alambda/R_earth)/20.
      if(dr.lt.1.e-4)dr=1.e-4

c......... density profile felt by neutrino ...............
      n_passos = cosphi/dr
      if(n_passos.lt.100)n_passos=100
      dr=cosphi/float(n_passos)
      do i=1,n_passos
         rr_(i) = dr*dfloat(i)
         r=sqrt(sin(phi)**2.0+(cos(phi)-2.0*rr_(i))**2.0)
         rho_(i)=divdif(rho_t,rr_t,55,r,1)
         if(rr_(i).gt.rr_t(55))rho_(i) = rho_t(55)
         rho_(i) = 0.5*sqrt(2.0)*Gf*6.02e23*rho_(i)/((5.068e4)**3.0)
      enddo
c..........................................................

       do i=1,2
       do j=1,2
          uf(i,j)=0.
          if(i.eq.j)uf(i,j)=cmplx(1.,0.)
       enddo
       enddo

c..............................
       ddr=dr*2.*R_earth*5.068e6
       do 12 jr = 1, n_passos
          ii=jr

          rho=rho_(ii)
          dm1=rho-sqrt((dm*c_2t-rho)**2.+(dm*s_2t)**2.)
          dm2=rho+sqrt((dm*c_2t-rho)**2.+(dm*s_2t)**2.)
          umass(1)=cmplx(cos(dm1*ddr),-sin(dm1*ddr))
          umass(2)=cmplx(cos(dm2*ddr),-sin(dm2*ddr))
          c_2tm=(dm*c_2t-rho)/sqrt((dm*c_2t-rho)**2.+(dm*s_2t)**2.)
          s_tm=sqrt(0.5*(1.-c_2tm))
          c_tm=sqrt(0.5*(1.+c_2tm))
          umix(1,1)=cmplx(c_tm,0.)
          umix(1,2)=cmplx(s_tm,0.)
          umix(2,1)=cmplx(-s_tm,0.)
          umix(2,2)=cmplx(c_tm,0.)
          do i=1,2
          do k=1,2
          ucte(i,k)=0.
          do j=1,2
          ucte(i,k)=ucte(i,k)+umix(i,j)*umass(j)*umix(k,j)
          enddo
          enddo
          enddo
          do i=1,2
          do j=1,2
             ufim(i,j)=0.
             do k=1,2
                ufim(i,j)=ufim(i,j)+ucte(i,k)*uf(k,j)
             enddo
          enddo
          enddo
          do i=1,2
          do j=1,2
             uf(i,j)=ufim(i,j)
          enddo
          enddo

      p2etmp=abs(s_t*uf(1,1)+c_t*uf(1,2))**2.
12    continue
c.............................
      do i=1,2
      do j=1,2
         uevol(i,j)=uf(i,j)
      enddo
      enddo
      p2e=abs(s_t*uevol(1,1)+c_t*uevol(1,2))**2.

 13   continue
      return
      end
         

