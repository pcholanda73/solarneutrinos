      subroutine chi2skiii(rskb8,rskhep,chi2sk3,ff)
      implicit real*4 (a-h,o-z)
      parameter (nbins=21) 
      dimension rskb8(2*nbins),rskhep(2*nbins),rsk(2*nbins)
      dimension sigma(2*nbins,2*nbins),ir(2*nbins)
c      dimension datask3(nbins),sigstd(nbins),sigstn(nbins)
c      dimension sigst1(nbins),sigst2(nbins)
c      dimension datask3d(nbins),datask3n(nbins),
      dimension datask3(2*nbins)
      dimension sigstd1(nbins),sigstd2(nbins)
      dimension sigstn1(nbins),sigstn2(nbins)
c      dimension ssmb8(nbins),ssmhep(nbins)
c      data ssmb8/189.7,172.2,155.2,134.3,117.1,101.2,85.8,71.7,58.5,
c     % 47.1,37.0,28.5,21.45,15.76,11.21,7.79,5.22,3.39,3.49,1.227,0.513/
c      data ssmhep/.334,.321,.310,.289,.271,.257,.240,.223,.205,.186,
c     %          .169,.151,.134,.118,.102,.088,.074,.062,.092,.059,.068/
      
c........  systematic energy-correlated errors 
c........  (energy scale, energy resolution and b8 spectrum)
c........  graphic reduction
      dimension sigsyscorr1(nbins),sigsyscorr2(nbins),sigsyscorr3(nbins)
      data sigsyscorr1/0.11,0.18,0.32,0.47,0.61,0.83,1.05,1.26,1.55
     %    ,1.84,2.20,2.64,3.07,3.50,4.01,4.66,5.23,5.96,6.97,8.70,11.37/
      data sigsyscorr2/-0.11,-0.11,-0.11,-0.18,-0.25,-0.18,-0.11,-0.04,
     %    0.04,.11,.32,0.54,0.76,1.12,1.55,2.06,2.71,3.43,.87,.33,11.95/
      data sigsyscorr3/0.18,0.25,0.25,0.25,0.32,0.40,0.40,0.47,0.47
     %    ,0.54,0.54,0.69,0.69,0.76,0.83,0.83,0.90,0.90,0.97,0.12,0.41/
      dimension sigsyscorr(nbins)
      dimension sigsysuncorr(nbins)
      data sigsysuncorr/4.3,3.0,2.6,1.5,1.3,.8,.8,.8,.8,.8,
     %                    .8,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8/
c........ 
      data datask3/93.4,73.7,55.3,50.8,55.6,39.6,37.2,28.4,19.8,17.7,
     %     15.0,14.7,9.36,5.24,4.08,2.67,1.59,1.13,2.00,0.381,0.244,
     %             72.6,59.9,70.4,58.7,52.1,41.1,35.7,32.6,24.9,20.3,
     %     13.6,12.9,9.44,6.04,5.69,3.38,2.25,1.48,2.31,1.208,0.000/
      data sigstd1/15.7,9.8,7.0,3.8,3.6,3.1,2.7,2.3,1.9,1.7,1.5,1.4,
     %             1.17,0.90,0.79,0.61,0.47,0.39,0.51,0.289,0.238/
 
      data sigstd2/14.9,9.3,6.5,3.7,3.5,3.0,2.6,2.2,1.8,1.6,1.4,1.3,
     %             1.03,0.76,0.66,0.49,0.35,0.27,0.40,0.158,0.117/
      data sigstn1/13.7,8.4,7.1,3.8,3.5,3.1,2.6,2.4,2.1,1.8,1.4,1.3,
     %             1.11,0.94,0.85,0.65,0.55,0.47,0.53,0.385,0.123/   
      data sigstn2/13.0,7.9,6.7,3.7,3.3,2.9,2.5,2.2,1.9,1.7,1.3,1.2,
     %             0.98,0.81,0.73,0.53,0.43,0.35,0.42,0.275,0.401/

      do i=1,2*nbins
         do j=1,2*nbins
            sigma(i,j)=0.
         enddo
      enddo
      do i=1,nbins
c     do j=1,nbins
         j=i
         sigma(i,j)=1.e-4*(sigsyscorr1(i)*sigsyscorr1(j)+
     %                     sigsyscorr2(i)*sigsyscorr2(j)+
     %                     sigsyscorr3(i)*sigsyscorr3(j))
     %                                      *datask3(i)*datask3(j)
c         sigma(i,j+nbins)=1.e-4*(sigsyscorr1(i)*sigsyscorr1(j)+
c     %                     sigsyscorr2(i)*sigsyscorr2(j)+
c     %                     sigsyscorr3(i)*sigsyscorr3(j))
c     %                                      *datask3(i)*datask3(j+nbins)
c         sigma(i+nbins,j)=sigma(i,j+nbins)
         sigma(i+nbins,j+nbins)=1.e-4*(sigsyscorr1(i)*sigsyscorr1(j)+
     %                     sigsyscorr2(i)*sigsyscorr2(j)+
     %                     sigsyscorr3(i)*sigsyscorr3(j))
     %                      *datask3(i+nbins)*datask3(j+nbins)
c      enddo
         sigma(i,i)=sigma(i,i)+
     %               (0.5*(sigstd1(i)+sigstd2(i)))**2.+
     %               1.e-4*sigsysuncorr(i)**2.*datask3(i)**2.
c         sigma(i,i+nbins)=sigma(i,i+nbins)+
c     %               (0.5*(sigstd1(i)+sigstd2(i)))**2.+
c     %               1.e-4*sigsysuncorr(i)**2.*datask3(i)*datask3(i+nbins)
c         sigma(i+nbins,i)=sigma(i,i+nbins)
         sigma(i+nbins,i+nbins)=sigma(i+nbins,i+nbins)+
     %                      (0.5*(sigstn1(i)+sigstn2(i)))**2.+
     %                   1.e-4*sigsysuncorr(i)**2.*datask3(i+nbins)**2.
      enddo

c      do i=1,2*nbins
c         write(*,*)(sigma(i,j),j=1,2*nbins)
c         write(*,*)
c      enddo
c      write(*,*)
      call rinv(2*nbins,sigma,2*nbins,ir,ifail)
c      do i=1,2*nbins
c         write(*,*)(sigma(i,j),j=1,2*nbins)
c         write(*,*)
c      enddo
c      stop
      
c      chi2min=1.e5
c      do ifac=0,800
c      fator=0.6+0.8*float(ifac)/800.
      fator=ff
      
      do i=1,nbins
         rsk(i)=fator*rskb8(i)+rskhep(i)
         rsk(i+nbins)=fator*rskb8(i+nbins)+rskhep(i+nbins)
      enddo
      
      chi2sk3=0.
      do i=1,2*nbins
      do j=1,2*nbins
         chi2sk3=chi2sk3+
     %           (rsk(i)-datask3(i))*(rsk(j)-datask3(j))*sigma(i,j)
      enddo
      enddo

c      if(chi2.lt.chi2min)then
c         chi2min=chi2
c         fmin=fator
c      endif
c      enddo
c      chi2=chi2min
c      ff=fmin
      
      return
      end
