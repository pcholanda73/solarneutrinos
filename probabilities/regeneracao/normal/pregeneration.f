c........................................................................
c..   calculate the regeneratioin probabilities for specific experiments
c........................................................................
      implicit real*4 (a-h,o-z)
      dimension c_t_(100)
      dimension tlsk1_(100)
      dimension tlsg1_(100),tlgn2a_(100),tlgn2b_(100)
      dimension tlsno_(100)
      common /timeexp/c_t_,tlsk1_,tlsg1_,tlgn2a_,tlgn2b_,tlsno_
      dimension pnight(7,5)
      common /pnight/pnight


c........................................................................
c..   reading zenithal exposition for SK, SAGE/GALLEX and SNO
c..   some calculated elsewhere, some given by collaborations
c........................................................................

      open(21,file='ltsk1.dat',status='old')
      open(31,file='ltsage1.dat',status='old')
      open(32,file='ltgallex2.dat',status='old')
      open(34,file='livetimesno2.dat',status='old')
      do j=1,100
         read(21,*)c_t_(j),tlsk1_(j)
         read(31,*)c_t_(j),tlsg1_(j)
         read(32,*)c_t_(j),tlgn2a_(j),tlgn2b_(j)
         read(34,*)c_t_(j),tlsno_(j)
      enddo
      close(21)
      close(31)
      close(32)
      close(34)

c........................................................................

      open(11,file='pregsk.tmp',status='unknown')
      open(12,file='pregsno.tmp',status='unknown')
      open(13,file='pregsg.tmp',status='unknown')
      open(14,file='preggno.tmp',status='unknown')
      ntg=45
      ndm=400
      do itg=0,ntg
         tg=0.1+0.02*float(itg)
         write(*,*)itg,tg
         do jdm=0,ndm
            if((jdm.gt.300).and.(2*(jdm/2).lt.jdm))goto 13
            if((jdm.gt.350).and.(4*(jdm/4).lt.jdm))goto 13
            dm=10.**(-15.+4.*float(jdm)/400.)
            call probearth(dm,tg)
            write(11,101)tg,dm,(pnight(k,1),k=1,7)
            write(12,101)tg,dm,pnight(7,2)
            write(13,101)tg,dm,pnight(7,3)
            write(14,101)tg,dm,pnight(7,4),pnight(7,5)
 13         continue
         enddo
         write(11,*)
         write(12,*)
         write(13,*)
         write(14,*)
      enddo
 101  format(10(e14.7,1x))

      end


