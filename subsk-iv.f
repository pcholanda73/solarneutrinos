      subroutine skiv(nenu,enu,probsk,rsk4b8,rsk4hep)
      implicit real*4 (a-h,o-z)
      parameter (nbins=23)
      dimension rskbssm(nbins),rskhssm(nbins)
      dimension enu(401),phepd(401),phepn(401),pb8d(401),pb8n(401)
      dimension e_b8(155),b8(155)
      dimension e_hep(1000),hep(1000)
      common /B8/ E_B8,B8
      common /hep/ E_hep,hep
      dimension flux(8)
      common /nuflux/flux
      dimension enusk4(196),celsk4(196,nbins),cmusk4(196,nbins)
      common /cssk4/enusk4,celsk4,cmusk4
      dimension probsk(8,nenu,2)
      dimension relb8(196,2*nbins),relhep(196,2*nbins)
      dimension rmub8(196,2*nbins),rmuhep(196,2*nbins)
      dimension rsk4b8(2*nbins),rsk4hep(2*nbins)
      dimension rrelb8(0:194),rrelhep(0:194)
      dimension rrmub8(0:194),rrmuhep(0:194)
      data rskbssm/19122.3555,17722.1621,16253.9531,14742.9980,
     %     13215.9121,11699.7266,10220.4971,8802.34082,7466.53271,
     %     6230.66992,5108.43359,4109.37549,3238.41772,2496.06445,
     %     1878.43237,1377.74512,983.039734,681.042175,457.246460,
     %     296.945953,298.719116,102.002228,38.7173576/
      data rskhssm/34.1817284,33.0443764,31.8033485,30.4627590,
     %     29.0297432,27.5140934,25.9277287,24.2843075,22.5987396,
     %     20.8865662,19.1640129,17.4479160,15.7549229,14.1014280,
     %     12.5030661,10.9744520,9.52873516,8.17738438,6.92985296,
     %     5.79339933,8.64319801,5.50062656,6.24590111/
      dimension ssmb(nbins),ssmh(nbins)
      data ssmb/196.8,182.8,167.8,153.3,137.8,121.9,106.8,
     %     92.1,78.0,65.2,53.4,42.9,33.8,26.0,19.55,14.34,
     %     10.24,7.10,4.80,3.11,3.18,1.117,0.464/
      data ssmh/0.346,0.335,0.323,0.312,0.298,0.282,0.266,
     %     0.250,0.232,0.214,0.197,0.179,0.162,0.144,0.128,
     %     0.112,0.097,0.083,0.070,0.059,0.088,0.056,0.064/
      
c      flb8=5.25e6     !used in SK-IV paper
c      flhep=8.e3      !used in SK-IV paper
      flb8=1.e10*flux(5)
      flhep=1.e10*flux(3)
      

      do j=1,nenu
         pb8d(j)=probsk(1,j,1)
         pb8n(j)=probsk(8,j,1)
         phepd(j)=probsk(1,j,2)
         phepn(j)=probsk(8,j,2)
      enddo
      do ienu=1,196
         eee=enusk4(ienu)
         b8fl=flb8*divdif(b8,e_b8,155,eee,2)
         hepfl=flhep*divdif(hep,e_hep,1000,eee,2)
         if(eee.gt.e_b8(155))b8fl=0.
         if(eee.gt.e_hep(1000))hepfl=0.
         if(eee.lt.e_b8(1))b8fl=0.
         if(eee.lt.e_hep(1))hepfl=0.
         dm=log10(deltam/4.e6/eee)
         probb8d=divdif(pb8d,enu,401,1.e6*eee,2)
         probb8n=divdif(pb8n,enu,401,1.e6*eee,2)
         probhepd=divdif(phepd,enu,401,1.e6*eee,2)
         probhepn=divdif(phepn,enu,401,1.e6*eee,2)
c.................ssm (to obtain normalization) .............
c         probb8d=1.
c         probhepd=1.
c         probb8n=1.
c         probhepn=1.
c.............................................................
         do ispec=1,nbins
         relb8(ienu,ispec)=celsk4(ienu,ispec)*b8fl*probb8d
         relhep(ienu,ispec)=celsk4(ienu,ispec)*hepfl*probhepd
         rmub8(ienu,ispec)=cmusk4(ienu,ispec)*b8fl*(1.-probb8d)
         rmuhep(ienu,ispec)=cmusk4(ienu,ispec)*hepfl*(1.-probhepd)
         relb8(ienu,ispec+nbins)=celsk4(ienu,ispec)*b8fl*probb8n
         relhep(ienu,ispec+nbins)=celsk4(ienu,ispec)*hepfl*probhepn
         rmub8(ienu,ispec+nbins)=cmusk4(ienu,ispec)*b8fl*(1.-probb8n)
         rmuhep(ienu,ispec+nbins)=cmusk4(ienu,ispec)*hepfl*(1.-probhepn)
         enddo
      enddo

c...................................................................
c.............  runing on electron kinetic energy bins .............
c...................................................................
      rtot=0.
      do ispec=1,nbins 
         tini = (4.0-0.511)+0.5*float(ispec-1)
         tend = (4.0-0.511)+0.5*float(ispec)
         if(ispec.eq.21)then
            tini=14.-0.511
            tend=15.-0.511
         elseif(ispec.eq.22)then
            tini=15.-0.511
            tend=16.-0.511
         elseif(ispec.eq.23)then
            tini=16.-0.511
            tend=20.-0.511
         endif
c...................................................................
c.............  integrating on neutrino energy .....................
c...................................................................
         do ienu=0,194 !eliminando o primeiro ponto e rearranjando
            rrelb8(ienu)=relb8(ienu+2,ispec)
            rrelhep(ienu)=relhep(ienu+2,ispec)
            rrmub8(ienu)=rmub8(ienu+2,ispec)
            rrmuhep(ienu)=rmuhep(ienu+2,ispec)
         enddo
         rrrelb8=simps(rrelb8,0.6,20.,194)
         rrrelhep=simps(rrelhep,0.6,20.,194)
         rrrmub8=simps(rrmub8,0.6,20.,194)
         rrrmuhep=simps(rrmuhep,0.6,20.,194)
         rsk4b8(ispec)=(rrrelb8+rrrmub8)*ssmb(ispec)/rskbssm(ispec)
         rsk4hep(ispec)=(rrrelhep+rrrmuhep)*ssmh(ispec)/rskhssm(ispec)

         do ienu=0,194 !eliminando o primeiro ponto e rearranjando
            rrelb8(ienu)=relb8(ienu+2,ispec+nbins)
            rrelhep(ienu)=relhep(ienu+2,ispec+nbins)
            rrmub8(ienu)=rmub8(ienu+2,ispec+nbins)
            rrmuhep(ienu)=rmuhep(ienu+2,ispec+nbins)
         enddo
         rrrelb8=simps(rrelb8,0.6,20.,194)
         rrrelhep=simps(rrelhep,0.6,20.,194)
         rrrmub8=simps(rrmub8,0.6,20.,194)
         rrrmuhep=simps(rrmuhep,0.6,20.,194)
         rsk4b8(ispec+nbins)=(rrrelb8+rrrmub8)
     %                                       *ssmb(ispec)/rskbssm(ispec)
         rsk4hep(ispec+nbins)=(rrrelhep+rrrmuhep)
     %                                       *ssmh(ispec)/rskhssm(ispec)
c...................................................................
c...................................................................

c         write(*,100)tini,rsk4b8(ispec),rsk4hep(ispec),
c     %               rsk4b8(ispec+nbins),rsk4hep(ispec+nbins)
c 100  format(f8.4,1x,4(f9.4,1x))
      enddo
c      stop
c...................................................................
c...................................................................
      
      return
      end
      
