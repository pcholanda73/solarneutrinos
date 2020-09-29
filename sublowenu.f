	subroutine rates(nenu,enu,probgach,rr)
	implicit real*4 (a-h,o-z)
c.. prob. survival
	dimension enu(nenu),probgach(nenu,8)
	dimension prob_ppsg(nenu),prob_pepsg(nenu),prob_hepsg(nenu)
        dimension prob_besg(nenu),prob_bsg(nenu),prob_nsg(nenu)
	dimension prob_osg(nenu),prob_fsg(nenu)
c.. rates
	dimension RR(8,3)
c	dimension rssm(9,2)
cc.. bahcall
c	data rssm/ 71.3874,  2.9300,  0.0566, 34.7825, 14.0562,  
c     %      1.8500,  2.6512,  0.0668,127.7807,0.0000,0.2375,0.0350,
c     %      1.1626,6.7412,0.0504,0.1623,0.0041,8.3931/ !2005op
	dimension flux(8)
	common /nuflux/flux
c.. cs
	dimension E_gal(58),cross_gal(58)
	dimension E_chl(19),cross_chl(19)
	common /bahcall_1/ E_gal,cross_gal
	common /bahcall_2/ E_chl,cross_chl
c.. fluxos
	dimension E_pp(84),pp(84)
	dimension E_B8(155),B8(155)
	dimension E_cno(502),n13(502),o15(502),f17(502)
	dimension E_hep(1000),hep(1000)
	common /pp/ E_pp,pp
	common /B8/ E_B8,B8
	common /CNO/ E_cno,n13,o15,f17
	common /hep/ E_hep,hep
c...........................................
	parameter (nint = 100)
	dimension fintegra(0:nint)

	nssm=0

	do i=1,8
	do j=1,3
	   rr(i,j)=0.
	enddo
	enddo

c------------------------------
	do i=1,nenu
	   prob_ppsg(i)= probgach(i,1)  
	   prob_pepsg(i)=probgach(i,2)	
	   prob_hepsg(i)=probgach(i,3)	
	   prob_besg(i)= probgach(i,4)	
	   prob_bsg(i)=  probgach(i,5)	
	   prob_nsg(i)=  probgach(i,6)	
	   prob_osg(i)=  probgach(i,7)	
	   prob_fsg(i)=  probgach(i,8)	
	enddo


c***********************************************************************
c******************  Gallex, Sage e Homestake **************************
c	n_exper = 1:  Gallex, GNO and Sage
c	n_exper = 2:  Homestake
c***********************************************************************

	do n_source=1,8
	do n_exper=1,2


c====================  fontes mono-energeticas  ========================
	if(n_source.eq.2) goto 32
	if(n_source.eq.4) goto 34

c====================  fontes espectrais ===============================
	nn = n_source

c..................... minimo dado pelo experimento considerado
	if(n_exper.eq.1)then
	   emin = .233
	elseif(n_exper.eq.2)then
	   emin = 0.814
	endif
c..................... maximo dado pela fonte de neutrinos
	if(n_source.eq.1)then
	   emax = .42341
	elseif(n_source.eq.3)then
	   emax = 18.784
	elseif(n_source.eq.5)then
	   emax = 16.36
	elseif(n_source.eq.6)then
	   emax = 1.202
	elseif(n_source.eq.7)then
	   emax = 1.732
	elseif(n_source.eq.8)then
	   emax = 1.7389
	endif

	if(emax.lt.emin)then
	   RR(n_source,n_exper) = 0.0
	   goto 40
	endif

c..................................
	do ii=0,nint
	energy = emin+(emax-emin)*float(ii)/float(nint)
c........
	if(n_exper.eq.1)then
	   ccross = divdif(cross_gal,e_gal,58,energy,3)
	elseif(n_exper.eq.2)then
	   ccross = divdif(cross_chl,e_chl,19,energy,3)
	endif
	if(ccross.lt.0.)ccross=0.
c........
	if(n_source.eq.1)then
	   ffluxo = divdif(pp,e_pp,84,energy,3)
	   pprob=divdif(prob_ppsg,enu,nenu,1.e6*energy,3)
	elseif(n_source.eq.3)then
	   ffluxo = divdif(hep,e_hep,1000,energy,3)
	   pprob=divdif(prob_hepsg,enu,nenu,1.e6*energy,3)
	elseif(n_source.eq.5)then
	   ffluxo = divdif(b8,e_b8,155,energy,3)	
	   pprob=divdif(prob_bsg,enu,nenu,1.e6*energy,3)
	elseif(n_source.eq.6)then
	   ffluxo = divdif(n13,e_cno,502,energy,3)	
	   pprob=divdif(prob_nsg,enu,nenu,1.e6*energy,3)
	elseif(n_source.eq.7)then
	   ffluxo = divdif(o15,e_cno,502,energy,3)	
	   pprob=divdif(prob_osg,enu,nenu,1.e6*energy,3)
	elseif(n_source.eq.8)then
	   ffluxo = divdif(f17,e_cno,502,energy,3)	
	   pprob=divdif(prob_fsg,enu,nenu,1.e6*energy,3)
	endif
c........
	if(nssm.eq.1)pprob=1.
c........
	fintegra(ii) = ffluxo*ccross*pprob
	enddo
c..................................

	RR(n_source,n_exper)=
     %            flux(n_source)*simps(fintegra,emin,emax,nint)

	goto 40

c======================== pep ======================================
32	continue
	energy = 1.445

	if(nssm.eq.1)then
	   pprob=1.
	else
	   pprob=divdif(prob_pepsg,enu,nenu,1.e6*energy,3)
	endif

	if(n_exper.eq.1)ccross = divdif(cross_gal,e_gal,58,energy,3)
	if(n_exper.eq.2)ccross = divdif(cross_chl,e_chl,19,energy,3)

	RR(2,n_exper) = flux(2) * pprob * ccross
	goto 40

c======================== Be =======================================
34	continue
	energy1 = 0.8631
	energy2 = 0.3855

	if(nssm.eq.1)then
	   pprob1 = 1.
	   pprob2 = 1.
	else
	   pprob1=divdif(prob_besg,enu,nenu,1.e6*energy1,3)
	   pprob2=divdif(prob_besg,enu,nenu,1.e6*energy2,3)
	endif

	if(n_exper.eq.1)then
	   ccross1=divdif(cross_gal,e_gal,58,energy1,3)
	   ccross2=divdif(cross_gal,e_gal,58,energy2,3)
	   rex = 0.897*pprob1*ccross1+0.103*pprob2*ccross2
	elseif(n_exper.eq.2)then
	   ccross1 = 2.4
	   ccross2 = 0.0
	   rex = pprob1*ccross1
	endif
	RR(4,n_exper) = flux(n_source) * rex
	goto 40

	
c***********************************************************************
 40	continue

	enddo
	enddo

	do i=1,8
	   RR(i,3)=RR(i,2) !sage = gallex+gno
	   RR(i,2)=RR(i,1) !homestake
	enddo

	return
	end
	




