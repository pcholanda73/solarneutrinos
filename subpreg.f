      subroutine pregeneration(a2e,na2e,preg)
      implicit real*4 (a-h,o-z)
      dimension a2e(10000),na2e(2),x2e(2)
c      dimension pregsk1(10,339),pregsk2(10,339),pregsk3(10,339)
c      dimension pregsk4(10,339),pregsk5(10,339),pregsk6(10,339)
c      dimension pregsk7(10,339),pregsno(10,339),pregsgg(10,339)
c.. preg(i,j,k)
c..      i:  1 to 6 + 7  -> 6 zenith bins at SK + averaged night
c..          8           -> averaged night at SNO
c..          9           -> averaged night at GALLEX/GNO and SAGE
c..      j: 10 points in tg    
c..      k: 339 points in deltam/4E, (not equally spaced)
      dimension preg(9,46,339) 

      open(11,file='regeneracao/normal/pregsk.dat',status='old')
      open(12,file='regeneracao/normal/pregsno.dat',status='old')
      open(13,file='regeneracao/normal/pregsg.dat',status='old')
      ntg=46
      ndm=339
      na2e(1)=ntg
      na2e(2)=ndm

      do itg=1,ntg
         do jdm=1,ndm
            read(11,*)tg,dm,(preg(k,itg,jdm),k=1,7)
            read(12,*)tg,dm,preg(8,itg,jdm)
            read(13,*)tg,dm,preg(9,itg,jdm)
            a2e(ntg+jdm)=log10(dm)
         enddo
         a2e(itg)=tg
         read(11,*)
         read(12,*)
         read(13,*)
      enddo
      
      return
      end
      
