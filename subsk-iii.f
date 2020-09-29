      subroutine skiii(nenu,enu,probsk,rsk3b8,rsk3hep)
      implicit real*4 (a-h,o-z)
      parameter (nbins=21)
      dimension rskbssm(nbins),rskhssm(nbins)
      dimension enu(nenu),phepd(nenu),phepn(nenu),pb8d(nenu),pb8n(nenu)
      dimension e_b8(155),b8(155)
      dimension e_hep(1000),hep(1000)
      common /B8/ E_B8,B8
      common /hep/ E_hep,hep
      dimension flux(8)
      common /nuflux/flux
      dimension enusk3(196),celsk3(196,nbins),cmusk3(196,nbins)
      common /cssk3/enusk3,celsk3,cmusk3
      dimension probsk(8,nenu,2)
      dimension relb8(196,2*nbins),relhep(196,2*nbins)
      dimension rmub8(196,2*nbins),rmuhep(196,2*nbins)
      dimension rsk3b8(2*nbins),rsk3hep(2*nbins)
      dimension rrelb8(0:194),rrelhep(0:194)
      dimension rrmub8(0:194),rrmuhep(0:194)
      data rskbssm/16246.9912,14736.4590,13210.1875,11694.9678,
     %     10216.8193,8799.71973,7464.80420,6229.68018,5107.99414,
     %     4109.22852,3238.43579,2496.07349,1878.34875,1377.52539,
     %     982.678894,680.578430,456.730621,296.430481,297.854309,
     %     101.465950,38.3304329/
      data rskhssm/31.7931805,30.4524746,29.0198669,27.5050259,
     %     25.9198418,24.2778511,22.5937614,20.8831959,19.1622353,
     %     17.4475384,15.7558861,14.1034470,12.5059481,10.9778719,
     %     9.53240585,8.18098927,6.93311787,5.79609537,8.64636230,
     %     5.50050831,6.23191881/
      dimension ssmb(nbins),ssmh(nbins)
      data ssmb/189.7,172.2,155.2,134.3,117.1,101.2,85.8,71.7,58.5,
     % 47.1,37.0,28.5,21.45,15.76,11.21,7.79,5.22,3.39,3.49,1.227,0.513/
      data ssmh/.334,.321,.310,.289,.271,.257,.240,.223,.205,.186,
     %          .169,.151,.134,.118,.102,.088,.074,.062,.092,.059,.068/

c      flb8=5.25e6
c      flhep=8.e3
      flb8=1.e10*flux(5)
      flhep=1.e10*flux(3)
      
      do j=1,nenu
         pb8d(j)=probsk(1,j,1)
         pb8n(j)=probsk(8,j,1)
         phepd(j)=probsk(1,j,2)
         phepn(j)=probsk(8,j,2)
      enddo
      do ienu=1,196
         eee=enusk3(ienu)
         b8fl=flb8*divdif(b8,e_b8,155,eee,2)
         hepfl=flhep*divdif(hep,e_hep,1000,eee,2)
         if(eee.gt.e_b8(155))b8fl=0.
         if(eee.gt.e_hep(1000))hepfl=0.
         if(eee.lt.e_b8(1))b8fl=0.
         if(eee.lt.e_hep(1))hepfl=0.
         dm=log10(deltam/4.e6/eee)
         probb8d=divdif(pb8d,enu,nenu,1.e6*eee,2)
         probb8n=divdif(pb8n,enu,nenu,1.e6*eee,2)
         probhepd=divdif(phepd,enu,nenu,1.e6*eee,2)
         probhepn=divdif(phepn,enu,nenu,1.e6*eee,2)
c.................ssm (to obtain normalization) .............
c         probb8d=1.
c         probhepd=1.
c         probb8n=1.
c         probhepn=1.
c.............................................................
         do ispec=1,nbins
         relb8(ienu,ispec)=celsk3(ienu,ispec)*b8fl*probb8d
         relhep(ienu,ispec)=celsk3(ienu,ispec)*hepfl*probhepd
         rmub8(ienu,ispec)=cmusk3(ienu,ispec)*b8fl*(1.-probb8d)
         rmuhep(ienu,ispec)=cmusk3(ienu,ispec)*hepfl*(1.-probhepd)
         relb8(ienu,ispec+nbins)=celsk3(ienu,ispec)*b8fl*probb8n
         relhep(ienu,ispec+nbins)=celsk3(ienu,ispec)*hepfl*probhepn
         rmub8(ienu,ispec+nbins)=cmusk3(ienu,ispec)*b8fl*(1.-probb8n)
         rmuhep(ienu,ispec+nbins)=cmusk3(ienu,ispec)*hepfl*(1.-probhepn)
         enddo
      enddo

c...................................................................
c.............  runing on electron kinetic energy bins .............
c...................................................................
      rtot=0.
      do ispec=1,nbins 
         tini = (5.0-0.511)+0.5*float(ispec-1)
         tend = (5.0-0.511)+0.5*float(ispec)
         if(ispec.eq.19)then
            tini=14.-0.511
            tend=15.-0.511
         elseif(ispec.eq.20)then
            tini=15.-0.511
            tend=16.-0.511
         elseif(ispec.eq.21)then
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
         rsk3b8(ispec)=(rrrelb8+rrrmub8)*ssmb(ispec)/rskbssm(ispec)
         rsk3hep(ispec)=(rrrelhep+rrrmuhep)*ssmh(ispec)/rskhssm(ispec)

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
         rsk3b8(ispec+nbins)=(rrrelb8+rrrmub8)
     %                                      *ssmb(ispec)/rskbssm(ispec)
         rsk3hep(ispec+nbins)=(rrrelhep+rrrmuhep)
     %                                       *ssmh(ispec)/rskhssm(ispec)
c...................................................................
c...................................................................

c         write(*,100)tini,rsk3b8(ispec),rsk3hep(ispec),
c     %               rsk3b8(ispec+nbins),rsk3hep(ispec+nbins)
c 100  format(f8.4,1x,4(f9.4,1x))
      enddo
c      write(*,*)
c      stop
c...................................................................
c...................................................................
      
      return
      end
      
