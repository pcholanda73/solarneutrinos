      subroutine chi2ski(rsk1b8,rsk1hep,chi2sk1,ff)
      implicit real*4 (a-h,o-z)
      dimension rsk1b8(46),rsk1hep(46)
c.....................................................................
      dimension datask1(46)
      dimension statsk(46),sysunc(46)
      dimension syscor1(46),syscor2(46),syscor3(46)
      data datask1(1)/72.1/
      data datask1(2)/77.1/
      data (datask1(i),i=3,9)/127.,124.,106.,132.,146.,140.,119./
      data (datask1(i),i=10,16)/149.,166.,158.,137.,150.,141.,137./
      data (datask1(i),i=17,23)/87.8,90.7,92.1,90.5,99.8,90.3,88.5/
      data (datask1(i),i=24,30)/57.1,56.5,63.3,56.8,59.6,60.1,60.9/
      data (datask1(i),i=31,37)/18.7,20.0,13.8,15.3,19.5,17.0,20.4/
      data (datask1(i),i=38,44)/4.28,4.78,6.97,5.82,5.58,3.70,3.93/
      data datask1(45)/0.240/
      data datask1(46)/0.423/
c.......
      data statsk(1)/9.5/
      data statsk(2)/9.1/
      data (statsk(i),i=3,9)/6.,15.,14.,13.,13.,14.,15./
      data (statsk(i),i=10,16)/4.,11.,10.,8.,8.,9.,10./
      data (statsk(i),i=17,23)/2.6,7.0,6.6,5.7,5.8,6.2,6.8/
      data (statsk(i),i=24,30)/1.9,5.0,4.8,4.0,4.1,4.5,5.0/
      data (statsk(i),i=31,37)/1.0,2.6,2.2,1.9,2.1,2.2,2.5/
      data (statsk(i),i=38,44)/.46,1.27,1.45,1.1,1.07,1.0,1.1/
      data statsk(45)/.120/
      data statsk(46)/.132/
c.......
      data sysunc/
     % 3.2,3.2,
     % 1.7,1.7,1.7,1.7,1.7,1.7,1.7,
     % 1.5,1.5,1.5,1.5,1.5,1.5,1.5,
     % 1.4,1.4,1.4,1.4,1.4,1.4,1.4,
     % 1.4,1.4,1.4,1.4,1.4,1.4,1.4,
     % 1.4,1.4,1.4,1.4,1.4,1.4,1.4,
     % 1.4,1.4,1.4,1.4,1.4,1.4,1.4,
     % 1.4,1.4/
c.......
      data syscor1/
     % -0.05,-0.05,
     %  0.18,0.18,0.18,0.18,0.18,0.18,0.18,
     %  0.67,0.67,0.67,0.67,0.67,0.67,0.67,
     %  1.42,1.42,1.42,1.42,1.42,1.42,1.42,
     %  2.55,2.55,2.55,2.55,2.55,2.55,2.55,
     %  4.19,4.19,4.19,4.19,4.19,4.19,4.19,
     %  6.47,6.47,6.47,6.47,6.47,6.47,6.47,
     %  10.05,10.05/
      data syscor2/
     % -0.20,-0.20,
     % -0.20,-0.20,-0.20,-0.20,-0.20,-0.20,-0.20,
     % -0.20,-0.20,-0.20,-0.20,-0.20,-0.20,-0.20,
     % -0.13,-0.13,-0.13,-0.13,-0.13,-0.13,-0.13,
     %  0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20,
     %  1.23, 1.23, 1.23, 1.23, 1.23, 1.23, 1.23,
     %  3.43, 3.43, 3.43, 3.43, 3.43, 3.43, 3.43,
     %  8.15, 8.15/
      data syscor3/
     % -0.05,-0.05,
     %  0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
     %  0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38,
     %  0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
     %  1.76, 1.76, 1.76, 1.76, 1.76, 1.76, 1.76,
     %  3.04, 3.04, 3.04, 3.04, 3.04, 3.04, 3.04,
     %  4.78, 4.78, 4.78, 4.78, 4.78, 4.78, 4.78,
     %  7.15, 7.15/
c.....................................................................
c.....................................................................
c.. all
      dimension rsk1(46),rexall(46)
c      dimension s_st(46,46),s_ap(46,46),s_sys(46,46)
      dimension s_sys(46,46)
      dimension ir(46),sigma(46,46)
      dimension sigsk(46,46),irsk(46)

      fator=ff

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  experimental data and statistical errors squared  ccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      do j1=1,46
c         s_st(j1,j1)=statsk(j1)**2.
c      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc    theoretical prevision and systematic errors       ccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do j1=1,46
         rsk1(j1)=fator*rsk1b8(j1)+rsk1hep(j1)
      enddo
      
      do j1=1,46
      do j2=1,46
      n=0
      if(((j1.ge.1).and.(j1.le.2)).and.((j2.ge.1).and.(j2.le.2)))n=1
      if(((j1.ge.3).and.(j1.le.9)).and.((j2.ge.3).and.(j2.le.9)))n=1
      if(((j1.ge.10).and.(j1.le.16)).and.((j2.ge.10).and.(j2.le.16)))n=1
      if(((j1.ge.17).and.(j1.le.23)).and.((j2.ge.17).and.(j2.le.23)))n=1
      if(((j1.ge.24).and.(j1.le.30)).and.((j2.ge.24).and.(j2.le.30)))n=1
      if(((j1.ge.31).and.(j1.le.37)).and.((j2.ge.31).and.(j2.le.37)))n=1
      if(((j1.ge.38).and.(j1.le.44)).and.((j2.ge.38).and.(j2.le.44)))n=1
      if(((j1.ge.45).and.(j1.le.46)).and.((j2.ge.45).and.(j2.le.46)))n=1
         s_sys(j1,j2)=float(n)*sysunc(j1)**2.0 +
     %         (syscor1(j1)*syscor1(j2)+syscor2(j1)*syscor2(j2)+
     %         syscor3(j1)*syscor3(j2))*rsk1(j1)*rsk1(j2)/1.e4
      enddo
      enddo

      do j1=1,46
      do j2=1,46
         sigsk(j1,j2)=s_sys(j1,j2)
      enddo
      sigsk(j1,j1)=sigsk(j1,j1)+statsk(j1)**2.
      enddo

c................................

      call rinv(46,sigsk,46,irsk,ifail)
      
      chi2sk1=0.
      do j1=1,46
      do j2=1,46
	soma = (rsk1(j1)-datask1(j1))
     %         *(rsk1(j2)-datask1(j2))*sigsk(j1,j2)
        chi2sk1=chi2sk1+soma
      enddo
      enddo

      return
      end






