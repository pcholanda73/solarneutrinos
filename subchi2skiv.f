      subroutine chi2skiv(rskb8,rskhep,chi2sk4,ff)
      implicit real*4 (a-h,o-z)
      parameter (nbins=23) 
      dimension rskb8(2*nbins),rskhep(2*nbins),rsk(2*nbins)
      dimension sigma(2*nbins,2*nbins),ir(2*nbins)
c      dimension datask4(nbins),sigstd(nbins),sigstn(nbins)
c      dimension sigst1(nbins),sigst2(nbins)
c      dimension datask4d(nbins),datask4n(nbins),
      dimension datask4(2*nbins)
      dimension sigstd1(nbins),sigstd2(nbins)
      dimension sigstn1(nbins),sigstn2(nbins)
c      dimension ssmb8(nbins),ssmhep(nbins)
c      data ssmb8/196.8,182.8,167.8,153.3,137.8,121.9,106.8,
c     %     92.1,78.0,65.2,53.4,42.9,33.8,26.0,19.55,14.34,
c     %     10.24,7.10,4.80,3.11,3.18,1.117,0.464/
c      data ssmhep/0.346,0.335,0.323,0.312,0.298,0.282,0.266,
c     %     0.250,0.232,0.214,0.197,0.179,0.162,0.144,0.128,
c     %     0.112,0.097,0.083,0.070,0.059,0.088,0.056,0.064/
c     
c........  systematic energy-correlated errors 
c........  (energy scale, energy resolution and b8 spectrum)
c........  graphic reduction
      dimension sigsyscorr1(nbins),sigsyscorr2(nbins),sigsyscorr3(nbins)
      data sigsyscorr1/-0.27,-0.20,-0.16,0.18,0.29,0.44,0.55,0.78,
     %     0.96,1.23,1.49,1.83,2.16,2.54,3.02,3.51,4.03,4.63,5.27,
     %     5.98,7.03,8.86,8.86/
      data sigsyscorr2/-0.12,-0.08,-0.08,-0.01,-0.01,0.03,0.03,0.03,
     %     0.03,0.03,0.03,0.03,0.07,0.10,0.14,0.18,0.37,0.48,0.63,0.78,
     %     1.15,1.71,1.71/
      data sigsyscorr3/0.07,0.07,0.14,0.18,0.18,0.22,0.25,0.25,0.29,
     %     0.33,0.37,0.44,0.52,0.55,0.59,0.59,0.67,0.74,0.81,0.81,
     %     0.96,1.08,1.08/
      dimension sigsyscorr(nbins)
      dimension sigsysuncorr(nbins)
      data sigsysuncorr/5.0,2.6,2.4,0.9,0.9,2.0,2.0,1.9,0.9,0.9,0.9,0.9,
     %                  0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9/
c........ 
c      data datask4/92.2,76.7,82.1,69.3,59.6,54.2,47.8,40.6,35.7,
c     %     29.1,24.0,18.5,14.5,10.7,8.43,6.60,4.40,3.04,
c     %     2.14,1.47,1.59,0.469,0.186/
c      data datask4d/96.0,64.6,79.4,65.9,58.3,51.0,45.7,41.8,35.0,
c     %     28.6,24.1,17.9,14.5,10.2,7.73,6.60,
c     %     3.83,3.04,2.41,1.48,1.54,0.486,0.150/
c      data datask4n/81.5,85.2,84.6,72.5,60.5,56.9,49.9,39.5,36.1,
c     %     28.9,23.7,19.2,14.4,11.1,9.23,6.72,
c     %     4.89,3.06,1.93,1.47,1.63,0.493,0.203/      
       data datask4/96.0,64.6,79.4,65.9,58.3,51.0,45.7,41.8,35.0,
     %     28.6,24.1,17.9,14.5,10.2,7.73,6.60,
     %     3.83,3.04,2.41,1.48,1.54,0.486,0.150,
     %              81.5,85.2,84.6,72.5,60.5,56.9,49.9,39.5,36.1,
     %     28.9,23.7,19.2,14.4,11.1,9.23,6.72,
     %     4.89,3.06,1.93,1.47,1.63,0.493,0.203/      
c     data sigst1/10.8,5.2,3.4,2.1,1.6,1.4,1.3,1.2,1.0,0.9,0.8,
c     %  0.7,0.6,0.5,0.43,0.37,0.30,0.25,0.20,0.17,0.17,0.102,0.072/
c      data sigst2/10.6,5.1,3.3,2.1,1.6,1.4,1.3,1.1,1.0,0.9,0.8,
c     %  0.7,0.6,0.5,0.41,0.35,0.28,0.23,0.18,0.15,0.15,0.082,0.051/
      data sigstd1/16.8,7.9,5.1,3.1,2.3,2.1,1.9,1.7,1.5,1.3,1.2,
     %  1.0,0.9,0.7,0.61,0.54,0.41,0.35,0.31,0.25,0.25,0.151,0.108/
      data sigstd2/16.3,7.6,5.0,3.0,2.2,2.0,1.8,1.7,1.5,1.3,1.1,
     %  0.9,0.8,0.7,0.56,0.49,0.37,0.31,0.27,0.21,0.21,0.112,0.065/
      data sigstn1/14.0,6.9,4.6,3.0,2.2,2.0,1.9,1.6,1.5,1.3,1.1,
     %  1.0,0.8,0.7,0.64,0.53,0.44,0.36,0.29,0.25,0.25,0.161,0.113/
      data sigstn2/13.6,6.7,4.5,2.9,2.2,2.0,1.8,1.6,1.4,1.2,1.1,
     %  0.9,0.8,0.7,0.60,0.49,0.40,0.32,0.25,0.21,0.22,0.121,0.071/

      
c      do i=1,23
c         sigsyscorr(i)=sqrt(sigsyscorr1(i)**2.+sigsyscorr2(i)**2.+
c     %                      sigsyscorr3(i)**2.)
c      enddo

      do i=1,46
         do j=1,46
            sigma(i,j)=0.
         enddo
      enddo
      do i=1,23
         do j=1,23
         sigma(i,j)=1.e-4*(sigsyscorr1(i)*sigsyscorr1(j)+
     %                     sigsyscorr2(i)*sigsyscorr2(j)+
     %                     sigsyscorr3(i)*sigsyscorr3(j))
     %                                      *datask4(i)*datask4(j)
         sigma(i+23,j+23)=1.e-4*(sigsyscorr1(i)*sigsyscorr1(j)+
     %                     sigsyscorr2(i)*sigsyscorr2(j)+
     %                     sigsyscorr3(i)*sigsyscorr3(j))
     %                                      *datask4(i+23)*datask4(j+23)
         enddo
         sigma(i,i+23)=1.e-4*(sigsyscorr1(i)*sigsyscorr1(i)+
     %                     sigsyscorr2(i)*sigsyscorr2(i)+
     %                     sigsyscorr3(i)*sigsyscorr3(i))
     %                                      *datask4(i)*datask4(i+23)
         sigma(i+23,i)=sigma(i,i+23)

         sigma(i,i)=sigma(i,i)+
     %               (0.5*(sigstd1(i)+sigstd2(i)))**2.+
     %               1.e-4*sigsysuncorr(i)**2.*datask4(i)**2.
         sigma(i+23,i+23)=sigma(i+23,i+23)+
     %                      (0.5*(sigstn1(i)+sigstn2(i)))**2.+
     %                      1.e-4*sigsysuncorr(i)**2.*datask4(i+23)**2.
      enddo

      call rinv(2*nbins,sigma,2*nbins,ir,ifail)
      fator=ff
      
      do i=1,nbins
c         rsk(i)=fator*(rskb8(i)*ssmb8(i)+rskhep(i)*ssmhep(i))
c         rsk(i+nbins)=fator*(rskb8(i+nbins)*ssmb8(i)+
c     %                     rskhep(i+nbins)*ssmhep(i))
         rsk(i)=fator*rskb8(i)+rskhep(i)
         rsk(i+nbins)=fator*rskb8(i+nbins)+rskhep(i+nbins)
      enddo
      
      chi2sk4=0.
      do i=1,2*nbins
      do j=1,2*nbins
         chi2sk4=chi2sk4+
     %           (rsk(i)-datask4(i))*(rsk(j)-datask4(j))*sigma(i,j)
      enddo
      enddo

      
      return
      end
