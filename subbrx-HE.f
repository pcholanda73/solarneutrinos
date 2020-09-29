      subroutine brx_he(nenu,enu,probbxnhep,probbxnb8,rbrx)
      implicit real*4 (a-h,o-z)
      parameter (nbins=6)
      parameter (n_prob=1000)   ! deve ser > nenu
      dimension ebrx(196),cbrxel(196,nbins),cbrxmu(196,nbins)
      dimension rbins(nbins),rbinspl(nbins),rbinsmi(nbins)
c      dimension E_B8(800),B8(800)
      dimension E_B8(155),B8(155)
      dimension E_hep(1000),hep(1000)
      dimension enu(n_prob)
      dimension probbxnhep(n_prob),probbxnb8(n_prob)
      dimension relb8(196,nbins),relhep(196,nbins)
      dimension rmub8(196,nbins),rmuhep(196,nbins)
      dimension rrelb8(0:194),rrelhep(0:194)
      dimension rrmub8(0:194),rrmuhep(0:194)
c      dimension fluxo(8)
      dimension rbrx(6)
      common /csbrx/ebrx,cbrxel,cbrxmu
c      common /B8brx/ E_B8brx,B8brx
      common /B8/ E_B8,B8
      common /hep/ E_hep,hep
c      common /ssm/nssm,fluxo
      dimension flux(8)
      common /nuflux/flux

      flb8=1.e4*flux(5)
      flhep=1.e4*flux(3)
c      flb8=5.88
c      flhep=7.91e-3

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      do ienu=1,196
         eee=ebrx(ienu)
         b8fl=flb8*divdif(B8,E_B8,155,eee,2)
         hepfl=flhep*divdif(hep,E_hep,1000,eee,2)
         if(eee.gt.E_B8(155))b8fl=0.
         if(eee.gt.E_hep(1000))hepfl=0.
         if(eee.lt.E_B8(1))b8fl=0.
         if(eee.lt.E_hep(1))hepfl=0.
         probb8=divdif(probbxnb8,enu,nenu,1.e6*eee,2)
         probhep=divdif(probbxnhep,enu,nenu,1.e6*eee,2)
         do ispec=1,nbins
            relb8(ienu,ispec)=cbrxel(ienu,ispec)*b8fl*probb8
            relhep(ienu,ispec)=cbrxel(ienu,ispec)*hepfl*probhep
            rmub8(ienu,ispec)=cbrxmu(ienu,ispec)*b8fl*(1.-probb8)
            rmuhep(ienu,ispec)=cbrxmu(ienu,ispec)*hepfl*(1.-probhep)
         enddo
      enddo

c..................... normalizacao obtida com ssm (P=1)
c..................... sem oscilacao espera-se 75 eventos * (2.4/5.88)
c.....................
      anorm=75.*5.88/2.4 /0.19002149 !ultimo numero obtido com anorm=1
c.....................
c.....................
      rtot=0.
      do ispec=1,nbins
         tini=3.+2.*float(ispec-1)
         tend=3.+2.*float(ispec)
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
         r=anorm*(rrrelb8+rrrelhep+rrrmub8+rrrmuhep)
         rbrx(ispec)=r
      enddo

c      write(*,*)rbrx
c      stop

      return
      end
