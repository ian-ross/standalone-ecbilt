












*********************************************************************
              SUBROUTINE CCSTAT(fracgr,darea)
*********************************************************************
        implicit none


*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
        INCLUDE 'comunit.h'
*********************************************************************

         real*8 fracgr(nlat,nlon), darea(nlat), tempor1
         integer indxv
C...1) calculation of initial carbon cycle parameters

	call CCPARAM
* calculation of equilibrium storages

* leaves biomass
* b1t is leaves phytomass for trees, b1g - for grass (kg C/m2)
* t1t is residence time of carbon in trees, t1g - in grass (years)

	b1t(lat,lon)=k1t*t1t*nppt
	b1g(lat,lon)=k1g*t1g*nppg

*   stems and roots biomass

	b2t(lat,lon)=(1-k1t)*t2t*nppt
	b2g(lat,lon)=(1-k1g)*t2g*nppg

*   litter

	b3t(lat,lon)=(k0t*b1t(lat,lon)/t1t+k2t/t2t*b2t(lat,lon))*t3t
	b3g(lat,lon)=(k0g*b1g(lat,lon)/t1g+k2g/t2g*b2g(lat,lon))*t3g

* mortmass and soil organic matter

	b4t(lat,lon)=(k3t/t3t*b3t(lat,lon))*t4t
	b4g(lat,lon)=(k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon))*t4g

* initialization of fraction dynamic variables
        
c
c until year 1992: distribution min & max is prescribed:
c after year 1992: vegetation is allowed to grow in the grid cell where no deforestation
c takes place in 1992. Vegetation can decline everywhere
c
        indxv=0

        if ((iscendef.eq.1).AND.((irunlabelf+iyear).ge.ivegstrt)) then
         if (irunlabelf+iyear.eq.ivegstrt) st_const(lat,lon)=st(lat,lon)
         tempor1=-99.99
         if((irunlabelf+iyear).le.1992) then
          indxv=(irunlabelf+iyear)
          if(indxv.eq.0) indxv=1

          if(farea(lon,lat,indxv).gt.0.D0) then
            tempor1=st_const(lat,lon)-farea(lon,lat,indxv) !orig
            if((lon.eq.47).and.(lat.eq.22))write(iuo+99, *) 
     &        "Veget(47,22,",VegetTime(indxv),")=", farea(47,22,indxv)
          endif
c        write(299,*) 'tempor1',tempor1,
c     &  'farea',farea(25,1,indxv),'indxv',indxv
         endif
         st(lat,lon)=min(forshare_st,tempor1)
        else
         st(lat,lon)=forshare_st
        endif
        if(st(lat,lon).lt.0.) st(lat,lon)=0.  

        sd(lat,lon)=desshare_st

        snlt(lat,lon)=nlshare_st
        if(sd(lat,lon).lt.0.) sd(lat,lon)=0.
        sg(lat,lon)=1.-st(lat,lon)-sd(lat,lon)
        if(sg(lat,lon).lt.0.) sg(lat,lon)=0.

	call CLIMPAR(fracgr,darea)
	
	return
	end
	
*********************************************************************
              SUBROUTINE CCDYN(fracgr,darea)
*********************************************************************
        implicit none
*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
        INCLUDE 'comunit.h'
*********************************************************************
* temporal var*/
        integer indxv
	real*8 tempor1,tempor2,db2,fd,dst,dd,nld,dstime
        real*8 dsd,temp_sg,temp_st
        real*8 fracgr(nlat,nlon), darea(nlat)
        
* calculation of current carbon cycle parameters

	call CCPARAM
	
* calculation of fraction dynamic variables

	fd=forshare_st-st(lat,lon)
	dd=desshare_st-sd(lat,lon)
	nld=nlshare_st-snlt(lat,lon)
        temp_st=st(lat,lon)
        temp_sg=sg(lat,lon)
        

* calculation of forest dynamics; exponential filtre
        dst=fd*(1.d0-exp(-1./t2t))
        st(lat,lon)=temp_st+dst
        if (st(lat,lon).lt.0.) st(lat,lon)=0.

* desert dynamics; exponential filtre
        dsd=dd*(1.-exp(-1./t2g))
* calculation of characteristic time of desert propagation
        tempor1=sd(lat,lon)+dsd+st(lat,lon)
        if (tempor1.gt.0.9) then
             dstime=(t2g*(1.-tempor1)+t2t*(tempor1-0.9))*10.
             dsd=dd*(1.-exp(-1./dstime))
        endif
        sd(lat,lon)=sd(lat,lon)+dsd
        if (sd(lat,lon).lt.0.) sd(lat,lon)=0.
*
c
c
        indxv=irunlabelf+iyear-i0dfor
        tempor1=0.
c
c ----------------------
c Constant vegetation:
c ----------------------
        if ((iscendef.eq.-1).AND.(indxv.ge.0)) then
c Defines ref. tree and desert fractions:
         if (indxv.eq.0) then
          st_const(lat,lon)=st(lat,lon)
          sd_const(lat,lon)=sd(lat,lon)
	 else
	  st(lat,lon)=st_const(lat,lon)
	  sd(lat,lon)=sd_const(lat,lon)
	 endif
	endif
c
c Available fraction for vegetation
        tempor2=1.-sd(lat,lon)
        if(tempor2.lt.0.) tempor2=0.
c
c ----------------------
c Deforestation scenario:
c ----------------------
c Forests/trees are replaced with grassland (cropland) in accordance with R&F 1700-1992 scenario.
c
c Once the end of the deforestation scenario file is reached the land use is kept at its state
c  as in the last year.
c
c Three versions: version A1 & A2 (AM & MFL, july 2008);
c                 version B (VB & AM, sept. 2008)
c                 version C (ED & AM for MILMO);
c
c SCENARIO version A ================= goes from here ====>>
c
c        if ((iscendef.eq.1).AND.(indxv.gt.0)) then
c         if(indxv.gt.ndfor) indxv=ndfor
c Version A1
c Scenario conc/efor
c actual fraction of trees is the minimum of
c     [potential fraction , fraction not occupied by cultures]
c         if(farea(lat,lon,indxv).lt.9.0E+19)
c     &       tempor1=tempor2-farea(lat,lon,indxv)
c         if(tempor1.lt.0.) tempor1=0.
c         st(lat,lon)=min(st(lat,lon),tempor1)
c Version A2
c Scenario Conc/Efor
c the tree fraction is reduced by the crop fraction
c         if(farea(lat,lon,indxv).lt.9.0E+19)
c     &    st(lat,lon)=st(lat,lon)-farea(lat,lon,indxv)
c         if(st(lat,lon).lt.0.) st(lat,lon)=0.
c        endif
c
c SCENARIO version A ======================== to there ====<<
c
c If replacing one version by the other be careful to comment and decomment
c  all instructions comprised between the tags "from here" and "to there".
c
c SCENARIO version B ================= goes from here ====>>
c
c Scenario as implemented by Victor Brovkin,
c  complies with the rules edicted in EMIC intercomparison project.
c Deforestation is computed with respect to a reference vegetation distribution;
c the reference distribution corresponds to that as in year ivgstrt.
c It goes like this:
c  - crop = 0 => sd=sd_ref; st=st_ref & sg=1-st-sd=sg_ref;
c  - crop > 0 => sd=sd_ref; st=max(0,st_ref-crop) ; sg=1-sd-st;
c

        if ((iscendef.eq.1).AND.(indxv.ge.0)) then
c Defines ref. tree and desert fractions:
         if((lon.eq.47).and.(lat.eq.22)) write(iuo+99, *) "Veget(47,22,",VegetTime(indxv),")=", farea(47,22,indxv)
         if (indxv.eq.0) then
          st_const(lat,lon)=st(lat,lon)
          sd_const(lat,lon)=sd(lat,lon)
         else
          if(indxv.gt.ndfor) indxv=ndfor
          sd(lat,lon)=sd_const(lat,lon)
          tempor2=1.-sd_const(lat,lon)
          st(lat,lon)=st_const(lat,lon)
          if(farea(lon,lat,indxv).gt.0.)
     &     st(lat,lon)=st_const(lat,lon)-farea(lon,lat,indxv)
c          write(299,*) 'farea',farea(25,1,indxv)
          if(st(lat,lon).lt.0.) st(lat,lon)=0.
         endif
        endif
c
c SCENARIO version B ======================== to there ====<<
c
c If replacing one version by the other be careful to comment and decomment
c  all instructions comprised between the tags "from here" and "to there".
c
c SCENARIO version C ================= goes from here ====>>
c
c Scenario as formerly written (MILMO)
c The tree fraction is the mimimum among:
c      - potential tree fraction
c      - reference tree-fraction (st_const) less the crop fraction (farea)
c This is not exactly the scenario preconised by V.B. for EMICs 
c        (since veg. cover does not necessarily correspond to the reference cover: st smaller, sd evolves...)
c It may produce deforestation fluxes even if farea=0 (if st < st_const)
c
c this part NEEDs TO BE VERIFIED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c        if ((iscendef.eq.1).AND.(indxv.ge.0)) then
c         if (indxv.eq.0) then
c Defines ref. tree fraction:
c          st_const(lat,lon)=st(lat,lon)
c         else
c          if(indxv.le.ndfor) then
c           if(farea(lat,lon,indxv).lt.9.0E+19)
c     &       tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
c          else
c           indxv=ndfor
c           if((farea(lat,lon,indxv).lt.9.0E+19).and.(farea(lat,lon,indxv).gt.0.D0))
c     &       tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
c          endif
c          if(tempor1.lt.0.) tempor1=0.
c          st(lat,lon)=min(st(lat,lon),tempor1)
c         endif
c
c SCENARIO version C ======================== to there ====<<
c
c
c do not forget to update dst (needed for corrections below)
        dst=st(lat,lon)-temp_st
c
        sg(lat,lon)=tempor2-st(lat,lon)
        if (sg(lat,lon).lt.0) sg(lat,lon)=0.
c
        snlt(lat,lon)=nlshare_st-nld*exp(-1./t2t)
c
* calculation of dynamics of storages

* calculation of changes of storages due to conservation law 

* correction for trees

        tempor1=b4t(lat,lon)
        tempor2=b3t(lat,lon)

        if(st(lat,lon).gt.0) then
cv         if(dst.gt.0) then  
            b4t(lat,lon)=(b4t(lat,lon)*temp_st
     >         +b4g(lat,lon)*dst)/st(lat,lon)
            b3t(lat,lon)=(b3t(lat,lon)*temp_st
     >         +b3g(lat,lon)*dst)/st(lat,lon)
cv         endif         
          b2t(lat,lon)=b2t(lat,lon)*temp_st/st(lat,lon)
          b1t(lat,lon)=b1t(lat,lon)*temp_st/st(lat,lon)
        endif
        
* correction for grass

cv       if(sg(lat,lon).gt.0) then	      
cv               if(dst.gt.0) then  
cv                b4g(lat,lon)=b4g(lat,lon)*(temp_sg-dst)
cv    >           /sg(lat,lon)
cv                b3g(lat,lon)=b3g(lat,lon)*(temp_sg-dst)
cv    >           /sg(lat,lon)
cv               else 
cv          b4g(lat,lon)=(b4g(lat,lon)*temp_sg-tempor1*dst)
cv    >     /sg(lat,lon)
cv          b3g(lat,lon)=(b3g(lat,lon)*temp_sg-tempor2*dst)
cv    >     /sg(lat,lon)
cv               endif
cv        b2g(lat,lon)=b2g(lat,lon)*temp_sg/sg(lat,lon)
cv        b1g(lat,lon)=b1g(lat,lon)*temp_sg/sg(lat,lon)
cv       endif
              
* slow soil organic matter

	b4t(lat,lon)=b4t(lat,lon)
     &                 +k3t/t3t*b3t(lat,lon)-b4t(lat,lon)/t4t
        b4g(lat,lon)=b4g(lat,lon)+k4g/t2g*b2g(lat,lon)
     &                 +k3g/t3g*b3g(lat,lon)-b4g(lat,lon)/t4g

*   fast soil organic matter

	b3t(lat,lon)=b3t(lat,lon)+b1t(lat,lon)/t1t*k0t
     &                 +k2t/t2t*b2t(lat,lon)-b3t(lat,lon)/t3t
	b3g(lat,lon)=b3g(lat,lon)+b1g(lat,lon)/t1g*k0g
     &                 +k2g/t2g*b2g(lat,lon)-b3g(lat,lon)/t3g

* leaves biomass

	b1t(lat,lon)=b1t(lat,lon)+k1t*nppt-b1t(lat,lon)/t1t
	b1g(lat,lon)=k1g*nppg*t1g

*   stems and roots biomass

	b2t(lat,lon)=b2t(lat,lon)+(1-k1t)*nppt-b2t(lat,lon)/t2t
	b2g(lat,lon)=b2g(lat,lon)+(1-k1g)*nppg-b2g(lat,lon)/t2g


	call CLIMPAR(fracgr,darea)
	
	return
	end

*********************************************************************
	    SUBROUTINE CLIMPAR(fracgr,darea)
*********************************************************************
        implicit none
*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
        include 'veget.h'
C
C Local variables
       REAL tempor1
       real*8 fracgr(nlat,nlon), darea(nlat)
       integer KTVM
*********************************************************************
         KTVM=2
* calculation of annual averaged LAI - lai

	 laig=b1g(lat,lon)*deng
         lait=b1t(lat,lon)*(dentn*snlt(lat,lon)+
     *   dentd*(1-snlt(lat,lon)))
     
  	 BLAI(lat,lon,1)=lait
  	 BLAI(lat,lon,2)=laig

* calculation of annual carbon uptake

	if (KTVM.eq.2) then 
	   tempor1=b1(lat,lon)+b2(lat,lon)+
     >          b3(lat,lon)+b4(lat,lon)
        else 
           tempor1=0.D0
        endif
        
        b1(lat,lon)=b1t(lat,lon)*st(lat,lon)+b1g(lat,lon)*sg(lat,lon)
        b2(lat,lon)=b2t(lat,lon)*st(lat,lon)+b2g(lat,lon)*sg(lat,lon)
        b3(lat,lon)=b3t(lat,lon)*st(lat,lon)+b3g(lat,lon)*sg(lat,lon)
        b4(lat,lon)=b4t(lat,lon)*st(lat,lon)+b4g(lat,lon)*sg(lat,lon)
        b12(lat,lon)=b1(lat,lon)+b2(lat,lon)
        b34(lat,lon)=b3(lat,lon)+b4(lat,lon)
	anup(lat,lon)=(b1(lat,lon)+b2(lat,lon)
     >       +b3(lat,lon)+b4(lat,lon)-tempor1)

        stock(lat,lon)=b1(lat,lon)+b2(lat,lon)+b3(lat,lon)
     &  +b4(lat,lon)
        stockloch(lat,lon)=fracgr(lat,lon)*darea(lat)*stock(lat,lon)*
     &   1E-12
        anuploch(lat,lon)=fracgr(lat,lon)*darea(lat)*anup(lat,lon)*1E-12
	anuploch(lat,lon)=fco2veg*anuploch(lat,lon)

c       anup_moy(lat,lon)=anup(lat,lon)*darea(lat)
c    &     *fracgr(lat,lon)

C...      NET PRIMARY PRODUCTION

c       pnpp(lat,lon)=npp*(1-sd(lat,lon))
        pnpp(lat,lon)=nppt*st(lat,lon)+nppg*sg(lat,lon)

	return
	end

*********************************************************************
	    SUBROUTINE CCPARAM
*********************************************************************

        implicit none
c
*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
        include 'veget.h'
        real*8  npp1,npp2,avefor,differ,pcr
        real*8  db1,db2,db3
*********************************************************************
* calculation of current cycle parameters

* potential trees share

       avefor=ave_pr*ave_pr*ave_pr*ave_pr
       differ=gdd0-gdd0_min
       db1=-bet*differ
       db2=gamm*differ
       db3=differ*differ

       if(differ.lt.0) then 
          forshare_st=0
       else 
          forshare_st=(1-exp(db1))*avefor/(avefor+a*db3*exp(db2))
       endif
       if (forshare_st.gt.fmax) forshare_st=fmax

* potential desert share - desshare_st

       desshare_st=0
	 
* northern deserts 

       if(gdd0.lt.100) desshare_st=1
         
       if(gdd0.ge.100.and.gdd0.lt.gdd0_min)
     >      desshare_st=(gdd0_min-gdd0)/(gdd0_min-100.)

* southern deserts

 	 if (gdd0.ge.gdd0_max) then
 	 
            pcr=acr*exp(gamm2/2.*differ)

            if (ave_pr05.le.pcr) then 
                desshare_st=1
                forshare_st=0
            else
                db2=(ave_pr05-pcr)/exp(gamm2*differ)
                desshare_st=1.03/(1+ades*db2*db2)-0.03
                if (desshare_st.lt.0) desshare_st=0
            endif
            
         endif

* calculation of NPP, Lieth''s formula

	db1=-v1*ave_pr
c	if (gdd0.ge.gdd0_max) db1=-v1*ave_pr05
	db2=-v2*ave_t
	npp1=(1.-exp(db1))
	npp2=1./(1.+v3*exp(db2))
	if(npp1.lt.npp2) then
                npp=nppmax*npp1
	    else
                npp=nppmax*npp2
	endif

* CO2 enrichment factor
c
c-AM (may 2007)
c Modified in order to account for different enhancement factors for tree and grass
c  -> betat & betag
c also division by log(2) is now included in beta_x (see veget.f)
c       npp=npp*(1.0+((0.25/LOG(2.))*LOG(co2ghg/280.)))
        nppt=npp*(1.0+(betat*LOG(co2ghg/280.)))
        nppg=npp*(1.0+(betag*LOG(co2ghg/280.)))

* allocation factors and residence time of leaves biomass

	k1t=c1t+c2t/(1+c3t*nppt)
	k1g=c1g+c2g/(1+c3g*nppg)

	t1t=d1t+d2t/(1+d3t*nppt)
	t1g=d1g+d2g/(1+d3g*nppg)

*   residence time of stems and roots biomass

	t2t=e1t+e2t/(1+e3t*nppt)
	t2g=e1g+e2g/(1+e3g*nppg)

*   residence time of fast carbon pool

        t3t=16.*exp(-ps5*(ave_t-soilt))  
        t3g=40.*exp(-ps5*(ave_t-soilt))  

* residence time of slow soil organic matter

        t4t=900.*exp(-ps5*(ave_t-soilt))
        t4g=t4t

*calculation of potential nedleleaves trees ratio

	nlshare_st=(t1t-t1td)/(t1tn-t1td)
	if (nlshare_st.gt.1) nlshare_st=1
	if (nlshare_st.lt.0) nlshare_st=0

	return
	end

*********************************************************************
	    SUBROUTINE INITCPAR
*********************************************************************

        implicit none
*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
       include 'veget.h'
*********************************************************************

* initialisation of variables
         ades=0.0011
         acr=28
         a=7000.
         bet=0.002
         gamm=0.00017
c        gamm2=0.00025
         gdd0_min=1000.
         gdd0_max=1800.
Cjmc0    gdd0_min=1000.
Cjmc0    gdd0_max=1700.
         fmax=0.9
	 nppmax=1.3
	 v1=0.000664
	 v2=0.119
	 v3=3.73
	 c1t=0.046
	 c2t=0.58
	 c3t=1.6
	 c1g=0.069
	 c2g=0.38
	 c3g=1.6
	 d1t=0.22
	 d2t=7.19
	 d3t=5.5
	 d1g=0.6
	 d2g=0.41
	 d3g=6.0
	 e1t=17.9
	 e2t=167.3
	 e3t=15.
	 e1g=0.67
	 e2g=50.5
	 e3g=100.
	 f1t=0.43
	 f2t=24.3
	 f3t=13.
	 f1g=0.34
	 f2g=17.8
	 f3g=50.
         k2t=1.
         k3t=0.025
         k0t=0.6
         k0g=0.2
         k2g=0.55
         k4g=0.025
         k3g=0.025
         t3g=1.
    	 t1tn=4
	 t1td=1
	 deng=20
	 dentd=20
    	 dentn=6
	 ps5=0.04
	 soilt=5
         acwd=100
         acwt=100
         acwg=100
         acwn=100
         zrd=1.
	 zrt=1.
	 zrg=0.6
	 zrn=0.6
         rsd=0.
	 rst=300.
	 rsg=130.
	 rsn=160.

	 return
	 end
*********************************************************************
              SUBROUTINE CCSTATR(fracgr,darea)
*********************************************************************
        implicit none
*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
*********************************************************************

         real*8 fracgr(nlat,nlon), darea(nlat), tempor1

        call CCPARAM


        stR(lat,lon)=forshare_st
        sdR(lat,lon)=desshare_st
        snltR(lat,lon)=nlshare_st
        if(sdR(lat,lon).lt.0.) sdR(lat,lon)=0.
        sgR(lat,lon)=1.-stR(lat,lon)-sdR(lat,lon)
        if(sgR(lat,lon).lt.0.) sgR(lat,lon)=0.

        call CLIMPAR(fracgr,darea)

        return
        end

*********************************************************************
              SUBROUTINE CCDYNR(fracgr,darea)
*********************************************************************
        implicit none
*       INCLUDE 'declar.inc'
*       INCLUDE 'params.inc'
C       INCLUDE 'bio.inc'
C       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
*********************************************************************
* temporal var*/
        integer indxv
        real*8 tempor1,tempor2,db2,fd,dst,dd,nld,dstime
        real*8 dsd,temp_sg,temp_st
        real*8 fracgr(nlat,nlon), darea(nlat)

* calculation of current carbon cycle parameters

        call CCPARAM

* calculation of fraction dynamic variables

        fd=forshare_st-stR(lat,lon)
        dd=desshare_st-sdR(lat,lon)
        nld=nlshare_st-snltR(lat,lon)
        temp_st=stR(lat,lon)
        temp_sg=sgR(lat,lon)


* calculation of forest dynamics; exponential filtre
c-AM    dst=forshare_st-fd*exp(-1./t2t)-st(lat,lon)
        dst=fd*(1.d0-exp(-1./t2t))
c
c
         stR(lat,lon)=stR(lat,lon)+dst
c
        if (stR(lat,lon).lt.0.) stR(lat,lon)=0.
c
        snltR(lat,lon)=nlshare_st-nld*exp(-1./t2t)

* desert dynamics; exponential filtre
        dsd=desshare_st-dd*exp(-1./t2g)-sdR(lat,lon)
        tempor1=sdR(lat,lon)+dsd+stR(lat,lon)

* calculation of characteristic time of desert propagation
        if (tempor1.gt.0.9) then
             dstime=t2g*(1-tempor1)*10.+t2t*(tempor1-0.9)*10.
             dsd=desshare_st-dd*exp(-1./dstime)-sdR(lat,lon)
        endif

        sdR(lat,lon)=sdR(lat,lon)+dsd
        if (sdR(lat,lon).lt.0) sdR(lat,lon)=0.
c
*
        sgR(lat,lon)=1.-stR(lat,lon)-sdR(lat,lon)
        if (sgR(lat,lon).lt.0) sgR(lat,lon)=0.



        return
        end

