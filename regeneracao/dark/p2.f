c.. all in one, except for probabilities
        implicit real*4 (a-h,o-z)
c..
        dimension c_t_(100)
        dimension tlsk1_(100)!,tlsk2a_(100),tlsk2b_(100)
c        dimension tlsk8a_(100),tlsk8b_(100),tlsk8c_(100),tlsk8d_(100)
        dimension tlsg1_(100),tlgn2a_(100),tlgn2b_(100)
c        dimension tlsg8a_(100),tlsg8b_(100),tlsg8c_(100),tlsg8d_(100)
        dimension tlsno_(100)
        common /timeexp/c_t_,tlsk1_,tlsg1_,tlgn2a_,tlgn2b_,tlsno_
c        external amm
        dimension pnight(7,5)
        common /pnight/pnight
        dimension pregsk1(8,139),pregsk2(8,139),pregsk3(8,139)
        dimension pregsk4(8,139),pregsk5(8,139),pregsk6(8,139)
        dimension pregsk7(8,139),pregsno(8,139),pregsgg(8,139)
        dimension preggn1(8,139),preggn2(8,139)
        dimension a2e(10000),na2e(2),x2e(2)
        common /preg/a2e,na2e,pregsk1,pregsk2,pregsk3,pregsk4,pregsk5,
     %               pregsk6,pregsk7,pregsno,pregsgg,preggn1,preggn2


c......... lendo a exposicao zenital de SK, sage e SNO2
      open(21,file='ltsk1.dat',status='old')
c      open(22,file='ltsk2.dat',status='old')
c      open(23,file='ltsk8.dat',status='old')
      open(31,file='ltsage1.dat',status='old')
      open(32,file='ltsage2.dat',status='old')
c      open(33,file='ltsage8.dat',status='old')
      open(34,file='livetimesno2.dat',status='old')
        do j=1,100
        read(21,*)c_t_(j),tlsk1_(j)
c        read(22,*)c_t_(j),tlsk2a_(j),tlsk2b_(j)
c        read(23,*)c_t_(j),tlsk8a_(j),tlsk8b_(j),tlsk8c_(j),tlsk8d_(j)
        read(31,*)c_t_(j),tlsg1_(j)
        read(32,*)c_t_(j),tlgn2a_(j),tlgn2b_(j)
c        read(33,*)c_t_(j),tlsg8a_(j),tlsg8b_(j),tlsg8c_(j),tlsg8d_(j)
        read(34,*)c_t_(j),tlsno_(j)
        enddo
      close(21)
c      close(22)
c      close(23)
      close(31)
c      close(32)
c      close(33)
      close(34)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        open(11,file='pregsk.2dat',status='unknown')
        open(12,file='pregsno.2dat',status='unknown')
        open(13,file='pregsg.2dat',status='unknown')
        open(14,file='preggno.2dat',status='unknown')
        ntg=9
        ndm=200       
        do itg=5,ntg
           tg=0.1+0.1*float(itg)
           write(*,*)itg,tg
c           a2e(itg+1)=tg
c           jjdm=0
           do jdm=ndm,-200,-1
              if((jdm.gt.100).and.(2*(jdm/2).lt.jdm))goto 13
              if((jdm.gt.150).and.(4*(jdm/4).lt.jdm))goto 13
              dm=-10.**(-13.+2.*float(jdm)/200.)
              write(*,*)'    ',jdm,dm
c              jjdm=jjdm+1
c              a2e(ntg+1+jjdm)=log10(dm)

              call probearth(dm,tg)
              write(11,101)tg,dm,(pnight(k,1),k=1,7)
              write(12,101)tg,dm,pnight(7,2)
              write(13,101)tg,dm,pnight(7,3)
              write(14,101)tg,dm,pnight(7,4),pnight(7,5)
 13           continue
           enddo
           write(11,*)
           write(12,*)
           write(13,*)
           write(14,*)
        enddo
 101    format(10(e14.7,1x))
c        na2e(1)=ntg+1
c        na2e(2)=jjdm
c        x2e(1)=0.45
c        x2e(2)=-11.1
c        pp2e=fint(2,x2e,na2e,a2e,pregsk1)
c        write(*,*)x2e(1),10.**x2e(2),pp2e
c        stop

        end


