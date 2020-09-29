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
c.................. neutrino flux, different sources ..........
c..............................................................
      dimension E_pp(84),pp(84)
      real*4 n13(502)
      dimension E_cno(502),o15(502),f17(502)
      dimension e_b8(155),b8(155)
      dimension e_hep(1000),hep(1000)
      common /pp/ E_pp,pp
      common /CNO/ E_cno,n13,o15,f17
      common /B8/ E_B8,B8
      common /hep/ E_hep,hep
      dimension flux(8)
      data flux/5.98,1.44e-2,7.98e-7,4.93e-1,
     %              5.46e-4,2.78e-2,2.05e-2,5.28d-4/
c      data flux/5.991,1.421E-02,7.930E-07,4.844E-01,
c     %      5.691E-04,3.066E-02,2.331E-02,5.836E-04/ 
      common /nuflux/flux
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
c..............................................................
c................ rates and chi2 SK-I .........................
c..............................................................
      dimension probsk(8,nenu,2)
      dimension enusk1(181),cselsk1(8,181),csmusk1(8,181)
      common /cssk1/enusk1,cselsk1,csmusk1
      dimension rsk1b8(46),rsk1hep(46)
c..............................................................
c................ rates and chi2 SK-III .......................
c..............................................................
      parameter (nbsk4=23) 
      dimension enusk4(196),celsk4(196,nbsk4),cmusk4(196,nbsk4)
      common /cssk4/enusk4,celsk4,cmusk4
      dimension rsk4b8(2*nbsk4),rsk4hep(2*nbsk4)
c..............................................................
c................ rates and chi2 SK-IV ........................
c..............................................................
      parameter (nbsk3=21) 
      dimension enusk3(196),celsk3(196,nbsk3),cmusk3(196,nbsk3)
      common /cssk3/enusk3,celsk3,cmusk3
      dimension rsk3b8(2*nbsk3),rsk3hep(2*nbsk3)
c..............................................................
c.................. SNO-all phases  ...........................
c..............................................................
      real*4 probsnodd(nenu),probsnonn(nenu)
      real*4 rsno(5)
      real*8 esno(124),ccd(124),ccn(124),eseld(124),eseln(124),
     %                                   esmud(124),esmun(124)
      common /fluxsno1/esno,ccd,ccn,eseld,eseln,esmud,esmun
c..............................................................
c.................. Borexino, 7Be line ........................
c..............................................................
      dimension probbe7dd(nenu),probbe7nn(nenu)
c..............................................................
c.................. Borexino, high-energy .....................
c..............................................................
      dimension probbxnhep(nenu),probbxnb8(nenu)
      dimension ebrx(196),cbrxel(196,6),cbrxmu(196,6),rbrx(6)
      common /csbrx/ebrx,cbrxel,cbrxmu
c..............................................................
c.................. Gallex/GNO, Sage, Homestake ...............
c..............................................................
      dimension probgach(nenu,8)
      dimension e_gal(58),cross_gal(58)
      dimension e_chl(19),cross_chl(19)
      dimension rr(8,3)
      common /bahcall_1/ E_gal,cross_gal
      common /bahcall_2/ E_chl,cross_chl
c..............................................................
c.................. KamLAND ...................................
c..............................................................
c      dimension chi2kamland(131,62)
c      dimension akl(10000),nakl(2),xkl(2)
      dimension chi2kamland(41,201,81)
      dimension akl(10000),nakl(3),xkl(3)
c..............................................................
c................. Projections ................................
c..............................................................
      parameter (ntg=50,ndm=90,ns13=10)
      dimension chi2tgdm(0:ntg,0:ndm,2),chi2t12t13(0:ntg,0:ns13,2)
      dimension chi2dmt13(0:ndm,0:ns13,2)


c..............................................................
c................... cross sections ...........................
c..............................................................
c. sk1 ..............................................     
      open(11,file='skspecel.zenithbins.dat',status='old')
      open(12,file='skspecmu.zenithbins.dat',status='old')
      do i=1,181
         read(11,*)enusk1(i),a,b,(cselsk1(j,i),j=1,8)
         read(12,*)enusk1tmp,a,b,(csmusk1(j,i),j=1,8)
      enddo
      close(11)
      close(12)
c. sk3 ..............................................     
      open(11,file='sk-iii.specel.dat',status='old')
      open(12,file='sk-iii.specmu.dat',status='old')
      do i=1,196
         read(11,*)enusk3(i),(celsk3(i,k),k=1,nbsk3)
         read(12,*)enutmpppp,(cmusk3(i,k),k=1,nbsk3)
      enddo
      close(11)
      close(12)
c. sk4 ..............................................     
      open(11,file='sk-iv.specel.dat',status='old')
      open(12,file='sk-iv.specmu.dat',status='old')
      do i=1,196
         read(11,*)enusk4(i),(celsk4(i,k),k=1,nbsk4)
         read(12,*)enutmpppp,(cmusk4(i,k),k=1,nbsk4)
      enddo
      close(11)
      close(12)
c. sno ..............................................     
      open(12,file='sno3phase.in',status='old')
      do i=1,124
         read(12,*)esno(i),ccd(i),ccn(i)
     %               ,eseld(i),eseln(i),esmud(i),esmun(i)
      enddo
      close(12)
      open(86,file='minorm.out')
      call mninit(5,86,7)
c. borexino high energy...................................     
      open(11,file='borexino-el.dat',status='old')
      open(12,file='borexino-mu.dat',status='old')
      do i=1,196
         read(11,*)ebrx(i),(cbrxel(i,k),k=1,6)
         read(12,*)ebrx(i),(cbrxmu(i,k),k=1,6)
      enddo
      close(11)
      close(12)
c. gallex/GNO and sage ...................................     
      open(11,file='galliumcross.dat',status='old')
      do i=1,58
         read(11,*)E_gal(i),cross_gal(i)
      enddo
      close(11)
c. homestake .............................................     
      open(12,file='chlorinecross.dat',status='old')
      do i=1,19
         read(12,*)E_chl(i),cross_chl(i)
      enddo
      close(12)
      
c..............................................................
c................... solar neutrino inputs ....................
c..............................................................
c.. neutrino spectrum
      open(11,file='ppenergy.dat',status='old')
      do i=1,84
         read(11,*)E_pp(i),pp(i)
      enddo
      close(11)
      open(13,file='cno.dat',status='old')
      do i=1,502
         read(13,*)E_cno(i),n13(i),o15(i),f17(i)
      enddo
      open(12,file='b8spec-2006.dat',status='old')
      do i=1,155
         read(12,*)e_b8(i),b8(i)
         b8(i)=1.e-3*b8(i)
      enddo
      close(12)	
      open(14,file='hepspectrum.dat',status='old')
      do i=1,1000
         read(14,*)e_hep(i),hep(i)
      enddo
      close(14)      
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
      
c..............................................................
c................... KamLAND final chi2 .......................
c..............................................................
c      open(41,file='chi2_KL.dat',status='old')
c      ntgkl=131
c      ndmkl=61       
c      nakl(1)=ntgkl
c      nakl(2)=ndmkl
c      chi2min=1.e5
c      do itg=1,ntgkl
c         do idm=1,ndmkl
c            read(41,*)tgkl,dmkl,chi2kamland(itg,idm)
c            tgkl=10.**tgkl
c            dmkl=10.**dmkl
c            akl(itg)=tgkl
c            akl(ntgkl+idm)=dmkl
c         enddo
c         read(41,*)
c      enddo
c      close(41)
      open(41,file='../KamLAND/chi2pois.dat',status='old')
      ntgkl=201
      ndmkl=81
      ns13kl=41
      nakl(1)=ns13kl
      nakl(2)=ntgkl
      nakl(3)=ndmkl
      do is13=1,ns13kl
      do itg=1,ntgkl
      do idm=1,ndmkl
         read(41,*)s13kl,tgkl,dmkl,chi2kamland(is13,itg,idm)
         akl(is13)=s13kl
         akl(ns13kl+itg)=tgkl
         akl(ns13kl+ntgkl+idm)=dmkl
      enddo
      read(41,*)
      enddo
      enddo
      close(41)
c      chi2min=1.e5
c      s2_13=0.
c      do  itg=0,100
c         tg=10.**(-1.+2.*float(itg)/100.)
c         do idm=0,40
c            deltam=(6.+4.*float(idm)/40.)*1.e-5
c            xkl(1)=s2_13
c            xkl(2)=tg
c            xkl(3)=deltam
c            chi2kl=fint(3,xkl,nakl,akl,chi2kamland)
c            if(chi2kl.lt.chi2min)chi2min=chi2kl
c            write(99,*)tg,deltam,chi2kl
c         enddo
c         write(99,*)
c      enddo
c      write(*,*)chi2min
c      stop

c..............................................................
c................... survival probabilities ...................
c..............................................................      
c.. Earth regeneration
      call pregeneration(a2e,na2e,preg)
c.. Survival probability
      open(21,file='chi2all.tmp',status='unknown')
      open(22,file='chi2exp.tmp',status='unknown')


c..............................................................
c................... runing neutrino parameters ...............
c..............................................................
      chi2mingl1=1.e5
      chi2mingl2=1.e5
      s213min=0.
      s213max=0.1
      tgmin=0.2
      tgmax=0.7
      dmmin=2.e-5
      dmmax=2.e-4
      do is13=0,ns13
         do itg=0,ntg
            chi2t12t13(itg,is13,1)=1.e5
            chi2t12t13(itg,is13,2)=1.e5
         enddo
      enddo
      do is13=0,ns13
         do idm=0,ndm
            chi2dmt13(idm,is13,1)=1.e5
            chi2dmt13(idm,is13,2)=1.e5
         enddo
      enddo
      do itg=0,ntg
         do idm=0,ndm
            chi2tgdm(itg,idm,1)=1.e5
            chi2tgdm(itg,idm,2)=1.e5
         enddo
      enddo

c.. runing 13 mixing angle
      do is13=0,ns13
      s2_13=s213min+(s213max-s213min)*float(is13)/float(ns13) !0.023
      write(*,*)is13   
c.. parameters to boron flux normalization
      fini0=1.
      fini=1.
      delf=1.e-4
      delff=delf
c.. runing 12 mixing angle
      do  itg=0,ntg
         tg=tgmin+(tgmax-tgmin)*float(itg)/float(ntg)
         write(*,*)'    ',tg,fini0
c.. probabilities in function of Deltam/4E
         do j=1,nenu
            dm=10.**(-13.+4.*float(j-1)/float(nenu-1))
            dm4e(j)=log10(dm)
            call probsunearth(dm,tg,s2_13,pday,pnight)
c            write(99,*)dm,(pday(k),k=1,8)
c            write(98,*)dm,(pnight(5,k),k=1,9)
            probsk(1,j,1)=pday(5) !8B
            probsk(1,j,2)=pday(3) !hep
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
c.. runing deltam
        do idm=0,ndm
            deltam=dmmin+(dmmax-dmmin)*float(idm)/float(ndm)
            do j=1,nenu
               enu(j)=deltam/10.**(dm4e(j))/4.
            enddo
c.. calculating detection rates
            call ski(nenu,enu,probsk,rsk1b8,rsk1hep)
            call skiii(nenu,enu,probsk,rsk3b8,rsk3hep)
            call skiv(nenu,enu,probsk,rsk4b8,rsk4hep)
            call sno(deltam,nenu,dm4e,probsnodd,probsnonn,rsno)
            dmbe7=log10(deltam/4.e6/0.8631)
            pbe7d=divdif(probbe7dd,dm4e,nenu,dmbe7,2)
            pbe7n=divdif(probbe7nn,dm4e,nenu,dmbe7,2)
            pprob=0.5*(pbe7d+pbe7n)
            rbe7=(pprob + (1.-pprob)*0.22)*74. !prevision in 0805.3843
            call brx_he(nenu,enu,probbxnhep,probbxnb8,rbrx)
            call rates(nenu,enu,probgach,rr)
c.. runing boron neutrino normalization
            f1=fini
            if(idm.eq.0)f1=fini0
            f2=f1+delf
c.. setting ff direction of change
            do ipar=1,2
               if(ipar.eq.1)ff=f1
               if(ipar.eq.2)ff=f2
c.. chi2 calculation
               call chi2ski(rsk1b8,rsk1hep,chi2sk1,ff)
               call chi2skiii(rsk3b8,rsk3hep,chi2sk3,ff)
               call chi2skiv(rsk4b8,rsk4hep,chi2sk4,ff)
               call chi2sno(rsno,chi2sno3f,ff)
               call chi2brxhe(rbrx,chi2brx)
               call chi2lowenu(rr,rbe7,chi2le,chi2gach,chi2be7,ff)
               chi2=chi2sk1+chi2sk3+chi2sk4+chi2sno3f+chi2brx+chi2le
               if(ipar.eq.1)chi2a=chi2
               if(ipar.eq.2)chi2b=chi2
            enddo
            if(chi2a.lt.chi2b)delff=-delf
c.. runing phib until minimum found
            chi2min=1.e5
            do ipar=0,10000
               ff=f1+delff*float(ipar)
               call chi2ski(rsk1b8,rsk1hep,chi2sk1,ff)
               call chi2skiii(rsk3b8,rsk3hep,chi2sk3,ff)
               call chi2skiv(rsk4b8,rsk4hep,chi2sk4,ff)
               call chi2sno(rsno,chi2sno3f,ff)
               call chi2brxhe(rbrx,chi2brx)
               call chi2lowenu(rr,rbe7,chi2le,chi2gach,chi2be7,ff)
               chi2=chi2sk1+chi2sk3+chi2sk4+chi2sno3f+chi2brx+chi2le
               if(chi2.le.chi2min)then
                  chi2min=chi2
                  chi2lemin=chi2le
                  chi2brxmin=chi2brx
                  chi2snomin=chi2sno3f
                  chi2skmin=chi2sk1+chi2sk3+chi2sk4
                  ffmin=ff
               else
                  goto 11
               endif
            enddo
 11         continue
            fini=ffmin
            if(idm.eq.0)fini0=ffmin
c.. KamLAND chi2
c            xkl(1)=tg
c            xkl(2)=deltam
c            chi2kl=fint(2,xkl,nakl,akl,chi2kamland)
            xkl(1)=s2_13
            xkl(2)=tg
            xkl(3)=deltam
            chi2kl=fint(3,xkl,nakl,akl,chi2kamland)
c..
            write(21,*)s2_13,tg,deltam,chi2min,ffmin,chi2min+chi2kl
            write(22,*)s2_13,tg,deltam,chi2skmin,chi2snomin,
     %                           chi2brxmin,chi2lemin
c            write(*,*)s2_13,tg,deltam,chi2min,ffmin,chi2min+chi2kl
c            write(*,*)tg,deltam,chi2skmin,chi2snomin,
c     %                           chi2brxmin,chi2lemin
c            stop
            if(chi2min.lt.chi2mingl1)chi2mingl1=chi2min
            if(chi2min+chi2kl.lt.chi2mingl2)
     %                             chi2mingl2=chi2min+chi2kl
            if(chi2min.lt.chi2tgdm(itg,idm,1))
     %                             chi2tgdm(itg,idm,1)=chi2min
            if(chi2min.lt.chi2t12t13(itg,is13,1))
     %                             chi2t12t13(itg,is13,1)=chi2min
            if(chi2min.lt.chi2dmt13(idm,is13,1))
     %                             chi2dmt13(idm,is13,1)=chi2min
            if(chi2min+chi2kl.lt.chi2tgdm(itg,idm,2))
     %                             chi2tgdm(itg,idm,2)=chi2min+chi2kl
            if(chi2min+chi2kl.lt.chi2t12t13(itg,is13,2))
     %                             chi2t12t13(itg,is13,2)=chi2min+chi2kl
            if(chi2min+chi2kl.lt.chi2dmt13(idm,is13,2))
     %                             chi2dmt13(idm,is13,2)=chi2min+chi2kl

         enddo
         write(21,*)
         write(22,*)
      enddo
      enddo

      open(31,file='chi2tgdm.dat',status='unknown')
      open(32,file='chi2t12t13.dat',status='unknown')
      open(33,file='chi2dmt13.dat',status='unknown')

      do is13=0,ns13
         s2_13=s213min+(s213max-s213min)*float(is13)/float(ns13)
         do  itg=0,ntg
            tg=tgmin+(tgmax-tgmin)*float(itg)/float(ntg)
            write(32,*)tg,s2_13,chi2t12t13(itg,is13,1)-chi2mingl1,
     %                          chi2t12t13(itg,is13,2)-chi2mingl2
         enddo
         write(32,*)
      enddo
      do is13=0,ns13
         s2_13=s213min+(s213max-s213min)*float(is13)/float(ns13)
         do idm=0,ndm
            deltam=dmmin+(dmmax-dmmin)*float(idm)/float(ndm)
            write(33,*)deltam,s2_13,chi2dmt13(idm,is13,1)-chi2mingl1,
     %                              chi2dmt13(idm,is13,2)-chi2mingl2
         enddo
         write(33,*)
      enddo
      do  itg=0,ntg
         tg=tgmin+(tgmax-tgmin)*float(itg)/float(ntg)
         do idm=0,ndm
            deltam=dmmin+(dmmax-dmmin)*float(idm)/float(ndm)
            write(31,*)tg,deltam,chi2tgdm(itg,idm,1)-chi2mingl1,
     %                           chi2tgdm(itg,idm,2)-chi2mingl2
         enddo
         write(31,*)
      enddo

      
      end


