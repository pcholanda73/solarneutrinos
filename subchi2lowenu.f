      subroutine chi2lowenu(rr,rbe7,chi2le,chi2gach,chi2be7,ff)
      implicit real*4 (a-h,o-z)
c.. incertezas astrofisicas hep-ph/9912231
      dimension corr_astr(8,8),sigma_astr(8)
      data corr_astr/
     %         1.000,0.969,0.111,-.163,0.111,0.459,0.416,0.580,
     %         0.969,1.000,0.123,-.062,0.216,0.588,0.546,0.682,
     %         0.111,0.123,1.000,0.038,0.057,0.097,0.090,0.109,
     %         -.163,-.062,0.038,1.000,0.889,0.507,0.497,0.509,
     %         0.111,0.216,0.057,0.889,1.000,0.721,0.705,0.719,
     %         0.459,0.588,0.097,0.507,0.721,1.000,0.997,0.857,
     %         0.416,0.546,0.090,0.497,0.705,0.997,1.000,0.826,
     %         0.580,0.682,0.109,0.509,0.719,0.857,0.826,1.000/
      data sigma_astr/
c     %   0.0105,0.0167,0.1549,0.1052,0.1631,0.2947,0.3075,0.5212/
     %   0.0105,0.0167,0.1549,0.1052,0.0000,0.2947,0.3075,0.5212/

c.....................................................................
c.....................................................................
c.. gallex/GNO, sage e homestake
      dimension rr(8,3),Del_lnC(3,8),del_lnX0(11)
      dimension rexgach(3)
cccccc
cccccc   Gallex + GNO, Sage and Homestake experimental data
      data rexgach/67.6,65.4,2.56/  
cccccc   Gallex + GNO, Sage and Homestake statystical and systematical errors
      data siggno_stt/4./
      data siggno_sys/3.2/
      data sigsage_sys/2.7/
      data sigsage_stt/3.1/
      data sigcl_sys/.16/
      data sigcl_stt/.16/
cccccc
      data Del_lnC/ .023, .023,  0., 
     %               .17,  .17,  .02, 
     %               .32,  .32,  .037,
     %               .07,  .07,  .02, 
     %               .32,  .32,  .032,
     %               .06,  .06,  .02, 
     %               .12,  .12,  .02, 
     %               .12,  .12,  .02/ 
      data Del_lnX0/.017,.06,.094,.143,.040,.004,.061,.004,.02,.02,.02/
c.. chi2
      dimension rall(8,4),rthall(4),rexall(4)
      dimension s_st(4,4),s_ap(4,4),s_sys(4,4)
      dimension ir(4),sigma(4,4)
      dimension siggach(3,3),irgach(3)
c
      do j1=1,4
      do j2=1,4
         s_st(j1,j2)=0. !sigma_st**2.
         sigma(j1,j2)=0.
      enddo
      enddo
      do j=1,4
      do k=1,8
      rall(k,j)=0.
      enddo
      enddo
      do i=1,4
      do j=1,4
         s_sys(i,j)=0.
      enddo
      enddo



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  experimental data and statistical errors squared  ccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.. gallex+gno
      s_st(1,1)=siggno_stt
      rexall(1)=rexgach(1)
c.. sage and homestake
      s_st(2,2)=sigsage_stt
      s_st(3,3)=sigcl_stt
      rexall(2)=rexgach(2)
      rexall(3)=rexgach(3)
c.. borexino
      s_st(4,4)=1.5**2.
      rexall(4)=46.0 !from 1104.1816

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc    theoretical prevision and systematic errors       ccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..  gallex+gno, sage and homestake
      do k1=1,8
         rall(k1,1)=rr(k1,1) !*rssm(k1,1)
         rall(k1,2)=rr(k1,2) !*rssm(k1,1)
         rall(k1,3)=rr(k1,3) !*rssm(k1,2)
         if(k1.eq.5)then
            do k2=1,3
               rall(k1,k2)=ff*rall(k1,k2)
            enddo
         endif
      enddo
      s_sys(1,1)=siggno_sys**2.
      s_sys(2,2)=sigsage_sys**2.
      s_sys(3,3)=sigcl_sys**2.
      do j=1,3
         do i1=1,8
            do i2=1,8
               n1=1
               n2=1
               if((i1.eq.3).or.(i1.eq.5))n1=-1
               if((i2.eq.3).or.(i2.eq.5))n2=-1
               if(n1*n2.eq.-1)icor=0
               if(n1*n2.eq.1)icor=1
               s_sys(j,j)=s_sys(j,j)+rall(i1,j)*Del_lnC(j,i1)*
     %              rall(i2,j)*Del_lnC(j,i2)*icor
            enddo
         enddo
      enddo
c.. borexino
      rall(4,4)=rbe7
      s_sys(4,4)=1.55**2.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc astrophysical theoretical flux correlation cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do j1=1,4
         do j2=1,4
            s_ap(j1,j2)=0.
            do i1=1,8
            do i2=1,8
               s_ap(j1,j2) = s_ap(j1,j2) + rall(i1,j1)*rall(i2,j2)*
     %             sigma_astr(i1)*sigma_astr(i2)*corr_astr(i1,i2)
c     1         1.e-2*sigma_astr(i1,i2)
            enddo
            enddo
         enddo
      enddo

      do j1=1,4
         rthall(j1)=0.
         do k=1,8
            rthall(j1)=rthall(j1)+rall(k,j1)
         enddo
      enddo

      do j1=1,4
      do j2=1,4
      sigma(j1,j2)=s_st(j1,j2)+s_sys(j1,j2)+s_ap(j1,j2)
      enddo
      enddo

      do j1=1,3
      do j2=1,3
         siggach(j1,j2)=sigma(j1,j2)
      enddo
      enddo

      sigbe7=sigma(4,4)

      call rinv(4,sigma,4,ir,ifail)
      call rinv(3,siggach,3,irgach,ifail)

      chi2le=0.
      do j1=1,4
      do j2=1,4
	soma = (rthall(j1)-rexall(j1))*(rthall(j2)-rexall(j2))
     %         *sigma(j1,j2)
        chi2le=chi2le+soma
      enddo
      enddo

      chi2gach=0.
      do j1=1,3
      do j2=1,3
	soma = (rthall(j1)-rexall(j1))*(rthall(j2)-rexall(j2))
     %         *siggach(j1,j2)
        chi2gach=chi2gach+soma
      enddo
      enddo

      chi2be7=(rthall(4)-rexall(4))**2./sigbe7

      return
      end






