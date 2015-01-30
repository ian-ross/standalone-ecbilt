!*********************************************************************
              SUBROUTINE CCSTAT(fracgr,darea)
!*********************************************************************
        implicit none


!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
        INCLUDE 'comunit.h'
!*********************************************************************

         real*8 fracgr(nlat,nlon), darea(nlat), tempor1
         integer indxv
!...1) calculation of initial carbon cycle parameters

	call CCPARAM
! calculation of equilibrium storages

! leaves biomass
! b1t is leaves phytomass for trees, b1g - for grass (kg C/m2)
! t1t is residence time of carbon in trees, t1g - in grass (years)

	b1t(lat,lon)=k1t*t1t*nppt
	b1g(lat,lon)=k1g*t1g*nppg

!   stems and roots biomass

	b2t(lat,lon)=(1-k1t)*t2t*nppt
	b2g(lat,lon)=(1-k1g)*t2g*nppg

!   litter

	b3t(lat,lon)=(k0t*b1t(lat,lon)/t1t+k2t/t2t*b2t(lat,lon))*t3t
	b3g(lat,lon)=(k0g*b1g(lat,lon)/t1g+k2g/t2g*b2g(lat,lon))*t3g

! mortmass and soil organic matter

	b4t(lat,lon)=(k3t/t3t*b3t(lat,lon))*t4t
	b4g(lat,lon)=(k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon))*t4g

! initialization of fraction dynamic variables

!
! until year 1992: distribution min & max is prescribed:
! after year 1992: vegetation is allowed to grow in the grid cell where no deforestation
! takes place in 1992. Vegetation can decline everywhere
!
        indxv=0

        if ((iscendef.eq.1).AND.((irunlabelf+iyear).ge.ivegstrt)) then
         if (irunlabelf+iyear.eq.ivegstrt) st_const(lat,lon)=st(lat,lon)
         tempor1=-99.99
         if((irunlabelf+iyear).le.1992) then
          indxv=(irunlabelf+iyear)
          if(indxv.eq.0) indxv=1

          if(farea(lon,lat,indxv).gt.0.D0) then
            tempor1=st_const(lat,lon)-farea(lon,lat,indxv) !orig
            if((lon.eq.47).and.(lat.eq.22))write(iuo+99, *) &
     &        "Veget(47,22,",VegetTime(indxv),")=", farea(47,22,indxv)
          endif
!        write(299,*) 'tempor1',tempor1,
!     &  'farea',farea(25,1,indxv),'indxv',indxv
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

!*********************************************************************
              SUBROUTINE CCDYN(fracgr,darea)
!*********************************************************************
        implicit none
!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
        INCLUDE 'comunit.h'
!*********************************************************************
! temporal var*/
        integer indxv
	real*8 tempor1,tempor2,db2,fd,dst,dd,nld,dstime
        real*8 dsd,temp_sg,temp_st
        real*8 fracgr(nlat,nlon), darea(nlat)

! calculation of current carbon cycle parameters

	call CCPARAM

! calculation of fraction dynamic variables

	fd=forshare_st-st(lat,lon)
	dd=desshare_st-sd(lat,lon)
	nld=nlshare_st-snlt(lat,lon)
        temp_st=st(lat,lon)
        temp_sg=sg(lat,lon)


! calculation of forest dynamics; exponential filtre
        dst=fd*(1.d0-exp(-1./t2t))
        st(lat,lon)=temp_st+dst
        if (st(lat,lon).lt.0.) st(lat,lon)=0.

! desert dynamics; exponential filtre
        dsd=dd*(1.-exp(-1./t2g))
! calculation of characteristic time of desert propagation
        tempor1=sd(lat,lon)+dsd+st(lat,lon)
        if (tempor1.gt.0.9) then
             dstime=(t2g*(1.-tempor1)+t2t*(tempor1-0.9))*10.
             dsd=dd*(1.-exp(-1./dstime))
        endif
        sd(lat,lon)=sd(lat,lon)+dsd
        if (sd(lat,lon).lt.0.) sd(lat,lon)=0.

        indxv=irunlabelf+iyear-i0dfor
        tempor1=0.
!
! ----------------------
! Constant vegetation:
! ----------------------
        if ((iscendef.eq.-1).AND.(indxv.ge.0)) then
! Defines ref. tree and desert fractions:
         if (indxv.eq.0) then
          st_const(lat,lon)=st(lat,lon)
          sd_const(lat,lon)=sd(lat,lon)
	 else
	  st(lat,lon)=st_const(lat,lon)
	  sd(lat,lon)=sd_const(lat,lon)
	 endif
	endif
!
! Available fraction for vegetation
        tempor2=1.-sd(lat,lon)
        if(tempor2.lt.0.) tempor2=0.
!
! ----------------------
! Deforestation scenario:
! ----------------------
! Forests/trees are replaced with grassland (cropland) in accordance with R&F 1700-1992 scenario.
!
! Once the end of the deforestation scenario file is reached the land use is kept at its state
!  as in the last year.
!
! Three versions: version A1 & A2 (AM & MFL, july 2008);
!                 version B (VB & AM, sept. 2008)
!                 version C (ED & AM for MILMO);
!
! SCENARIO version A ================= goes from here ====>>
!
!        if ((iscendef.eq.1).AND.(indxv.gt.0)) then
!         if(indxv.gt.ndfor) indxv=ndfor
! Version A1
! Scenario conc/efor
! actual fraction of trees is the minimum of
!     [potential fraction , fraction not occupied by cultures]
!         if(farea(lat,lon,indxv).lt.9.0E+19)
!     &       tempor1=tempor2-farea(lat,lon,indxv)
!         if(tempor1.lt.0.) tempor1=0.
!         st(lat,lon)=min(st(lat,lon),tempor1)
! Version A2
! Scenario Conc/Efor
! the tree fraction is reduced by the crop fraction
!         if(farea(lat,lon,indxv).lt.9.0E+19)
!     &    st(lat,lon)=st(lat,lon)-farea(lat,lon,indxv)
!         if(st(lat,lon).lt.0.) st(lat,lon)=0.
!        endif
!
! SCENARIO version A ======================== to there ====<<
!
! If replacing one version by the other be careful to comment and decomment
!  all instructions comprised between the tags "from here" and "to there".
!
! SCENARIO version B ================= goes from here ====>>
!
! Scenario as implemented by Victor Brovkin,
!  complies with the rules edicted in EMIC intercomparison project.
! Deforestation is computed with respect to a reference vegetation distribution;
! the reference distribution corresponds to that as in year ivgstrt.
! It goes like this:
!  - crop = 0 => sd=sd_ref; st=st_ref & sg=1-st-sd=sg_ref;
!  - crop > 0 => sd=sd_ref; st=max(0,st_ref-crop) ; sg=1-sd-st;
!

        if ((iscendef.eq.1).AND.(indxv.ge.0)) then
! Defines ref. tree and desert fractions:
         if((lon.eq.47).and.(lat.eq.22)) write(iuo+99, *) "Veget(47,22,",VegetTime(indxv),")=", farea(47,22,indxv)
         if (indxv.eq.0) then
          st_const(lat,lon)=st(lat,lon)
          sd_const(lat,lon)=sd(lat,lon)
         else
          if(indxv.gt.ndfor) indxv=ndfor
          sd(lat,lon)=sd_const(lat,lon)
          tempor2=1.-sd_const(lat,lon)
          st(lat,lon)=st_const(lat,lon)
          if(farea(lon,lat,indxv).gt.0.) &
     &     st(lat,lon)=st_const(lat,lon)-farea(lon,lat,indxv)
!          write(299,*) 'farea',farea(25,1,indxv)
          if(st(lat,lon).lt.0.) st(lat,lon)=0.
         endif
        endif
!
! SCENARIO version B ======================== to there ====<<
!
! If replacing one version by the other be careful to comment and decomment
!  all instructions comprised between the tags "from here" and "to there".
!
! SCENARIO version C ================= goes from here ====>>
!
! Scenario as formerly written (MILMO)
! The tree fraction is the mimimum among:
!      - potential tree fraction
!      - reference tree-fraction (st_const) less the crop fraction (farea)
! This is not exactly the scenario preconised by V.B. for EMICs
!        (since veg. cover does not necessarily correspond to the reference cover: st smaller, sd evolves...)
! It may produce deforestation fluxes even if farea=0 (if st < st_const)
!
! this part NEEDs TO BE VERIFIED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        if ((iscendef.eq.1).AND.(indxv.ge.0)) then
!         if (indxv.eq.0) then
! Defines ref. tree fraction:
!          st_const(lat,lon)=st(lat,lon)
!         else
!          if(indxv.le.ndfor) then
!           if(farea(lat,lon,indxv).lt.9.0E+19)
!     &       tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
!          else
!           indxv=ndfor
!           if((farea(lat,lon,indxv).lt.9.0E+19).and.(farea(lat,lon,indxv).gt.0.D0))
!     &       tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
!          endif
!          if(tempor1.lt.0.) tempor1=0.
!          st(lat,lon)=min(st(lat,lon),tempor1)
!         endif
!
! SCENARIO version C ======================== to there ====<<
!
!
! do not forget to update dst (needed for corrections below)
        dst=st(lat,lon)-temp_st
!
        sg(lat,lon)=tempor2-st(lat,lon)
        if (sg(lat,lon).lt.0) sg(lat,lon)=0.
!
        snlt(lat,lon)=nlshare_st-nld*exp(-1./t2t)
!
! calculation of dynamics of storages

! calculation of changes of storages due to conservation law

! correction for trees

        tempor1=b4t(lat,lon)
        tempor2=b3t(lat,lon)

        if(st(lat,lon).gt.0) then
!v         if(dst.gt.0) then
            b4t(lat,lon)=(b4t(lat,lon)*temp_st+b4g(lat,lon)*dst)/st(lat,lon)
            b3t(lat,lon)=(b3t(lat,lon)*temp_st+b3g(lat,lon)*dst)/st(lat,lon)
!v         endif
          b2t(lat,lon)=b2t(lat,lon)*temp_st/st(lat,lon)
          b1t(lat,lon)=b1t(lat,lon)*temp_st/st(lat,lon)
        endif

! correction for grass

!v       if(sg(lat,lon).gt.0) then
!v               if(dst.gt.0) then
!v                b4g(lat,lon)=b4g(lat,lon)*(temp_sg-dst)
!v    >           /sg(lat,lon)
!v                b3g(lat,lon)=b3g(lat,lon)*(temp_sg-dst)
!v    >           /sg(lat,lon)
!v               else
!v          b4g(lat,lon)=(b4g(lat,lon)*temp_sg-tempor1*dst)
!v    >     /sg(lat,lon)
!v          b3g(lat,lon)=(b3g(lat,lon)*temp_sg-tempor2*dst)
!v    >     /sg(lat,lon)
!v               endif
!v        b2g(lat,lon)=b2g(lat,lon)*temp_sg/sg(lat,lon)
!v        b1g(lat,lon)=b1g(lat,lon)*temp_sg/sg(lat,lon)
!v       endif

! slow soil organic matter

	b4t(lat,lon)=b4t(lat,lon) &
     &                 +k3t/t3t*b3t(lat,lon)-b4t(lat,lon)/t4t
        b4g(lat,lon)=b4g(lat,lon)+k4g/t2g*b2g(lat,lon) &
     &                 +k3g/t3g*b3g(lat,lon)-b4g(lat,lon)/t4g

!   fast soil organic matter

	b3t(lat,lon)=b3t(lat,lon)+b1t(lat,lon)/t1t*k0t &
     &                 +k2t/t2t*b2t(lat,lon)-b3t(lat,lon)/t3t
	b3g(lat,lon)=b3g(lat,lon)+b1g(lat,lon)/t1g*k0g &
     &                 +k2g/t2g*b2g(lat,lon)-b3g(lat,lon)/t3g

! leaves biomass

	b1t(lat,lon)=b1t(lat,lon)+k1t*nppt-b1t(lat,lon)/t1t
	b1g(lat,lon)=k1g*nppg*t1g

!   stems and roots biomass

	b2t(lat,lon)=b2t(lat,lon)+(1-k1t)*nppt-b2t(lat,lon)/t2t
	b2g(lat,lon)=b2g(lat,lon)+(1-k1g)*nppg-b2g(lat,lon)/t2g


	call CLIMPAR(fracgr,darea)

	return
	end

!********************************************************************
	    SUBROUTINE CLIMPAR(fracgr,darea)
!********************************************************************
        implicit none
!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
        include 'veget.h'
!
! Local variables
       REAL tempor1
       real*8 fracgr(nlat,nlon), darea(nlat)
       integer KTVM
!********************************************************************
         KTVM=2
! calculation of annual averaged LAI - lai

	 laig=b1g(lat,lon)*deng
         lait=b1t(lat,lon)*(dentn*snlt(lat,lon)+dentd*(1-snlt(lat,lon)))

  	 BLAI(lat,lon,1)=lait
  	 BLAI(lat,lon,2)=laig

! calculation of annual carbon uptake

	if (KTVM.eq.2) then
	   tempor1=b1(lat,lon)+b2(lat,lon)+b3(lat,lon)+b4(lat,lon)
        else
           tempor1=0.D0
        endif

        b1(lat,lon)=b1t(lat,lon)*st(lat,lon)+b1g(lat,lon)*sg(lat,lon)
        b2(lat,lon)=b2t(lat,lon)*st(lat,lon)+b2g(lat,lon)*sg(lat,lon)
        b3(lat,lon)=b3t(lat,lon)*st(lat,lon)+b3g(lat,lon)*sg(lat,lon)
        b4(lat,lon)=b4t(lat,lon)*st(lat,lon)+b4g(lat,lon)*sg(lat,lon)
        b12(lat,lon)=b1(lat,lon)+b2(lat,lon)
        b34(lat,lon)=b3(lat,lon)+b4(lat,lon)
	anup(lat,lon)=(b1(lat,lon)+b2(lat,lon)+b3(lat,lon)+b4(lat,lon)-tempor1)

        stock(lat,lon)=b1(lat,lon)+b2(lat,lon)+b3(lat,lon)+b4(lat,lon)
        stockloch(lat,lon)=fracgr(lat,lon)*darea(lat)*stock(lat,lon)*1E-12
        anuploch(lat,lon)=fracgr(lat,lon)*darea(lat)*anup(lat,lon)*1E-12
	anuploch(lat,lon)=fco2veg*anuploch(lat,lon)

!       anup_moy(lat,lon)=anup(lat,lon)*darea(lat)
!    &     *fracgr(lat,lon)

!...      NET PRIMARY PRODUCTION

!       pnpp(lat,lon)=npp*(1-sd(lat,lon))
        pnpp(lat,lon)=nppt*st(lat,lon)+nppg*sg(lat,lon)

	return
	end

!********************************************************************
	    SUBROUTINE CCPARAM
!********************************************************************

        implicit none
!
!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
        include 'veget.h'
        real*8  npp1,npp2,avefor,differ,pcr
        real*8  db1,db2,db3
!********************************************************************
! calculation of current cycle parameters

! potential trees share

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

! potential desert share - desshare_st

       desshare_st=0

! northern deserts

       if(gdd0.lt.100) desshare_st=1

       if(gdd0.ge.100.and.gdd0.lt.gdd0_min) &
            & desshare_st=(gdd0_min-gdd0)/(gdd0_min-100.)

! southern deserts

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

! calculation of NPP, Lieth''s formula

	db1=-v1*ave_pr
!	if (gdd0.ge.gdd0_max) db1=-v1*ave_pr05
	db2=-v2*ave_t
	npp1=(1.-exp(db1))
	npp2=1./(1.+v3*exp(db2))
	if(npp1.lt.npp2) then
                npp=nppmax*npp1
	    else
                npp=nppmax*npp2
	endif

! CO2 enrichment factor
!
!-AM (may 2007)
! Modified in order to account for different enhancement factors for tree and grass
!  -> betat & betag
! also division by log(2) is now included in beta_x (see veget.f)
!       npp=npp*(1.0+((0.25/LOG(2.))*LOG(co2ghg/280.)))
        nppt=npp*(1.0+(betat*LOG(co2ghg/280.)))
        nppg=npp*(1.0+(betag*LOG(co2ghg/280.)))

! allocation factors and residence time of leaves biomass

	k1t=c1t+c2t/(1+c3t*nppt)
	k1g=c1g+c2g/(1+c3g*nppg)

	t1t=d1t+d2t/(1+d3t*nppt)
	t1g=d1g+d2g/(1+d3g*nppg)

!   residence time of stems and roots biomass

	t2t=e1t+e2t/(1+e3t*nppt)
	t2g=e1g+e2g/(1+e3g*nppg)

!   residence time of fast carbon pool

        t3t=16.*exp(-ps5*(ave_t-soilt))
        t3g=40.*exp(-ps5*(ave_t-soilt))

! residence time of slow soil organic matter

        t4t=900.*exp(-ps5*(ave_t-soilt))
        t4g=t4t

!calculation of potential nedleleaves trees ratio

	nlshare_st=(t1t-t1td)/(t1tn-t1td)
	if (nlshare_st.gt.1) nlshare_st=1
	if (nlshare_st.lt.0) nlshare_st=0

	return
	end

!********************************************************************
	    SUBROUTINE INITCPAR
!********************************************************************

        implicit none
!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
       include 'veget.h'
!********************************************************************

! initialisation of variables
         ades=0.0011
         acr=28
         a=7000.
         bet=0.002
         gamm=0.00017
!        gamm2=0.00025
         gdd0_min=1000.
         gdd0_max=1800.
!jmc0    gdd0_min=1000.
!jmc0    gdd0_max=1700.
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
!********************************************************************
              SUBROUTINE CCSTATR(fracgr,darea)
!********************************************************************
        implicit none
!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
!********************************************************************

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

!********************************************************************
              SUBROUTINE CCDYNR(fracgr,darea)
!********************************************************************
        implicit none
!       INCLUDE 'declar.inc'
!       INCLUDE 'params.inc'
!       INCLUDE 'bio.inc'
!       INCLUDE 'buffer.inc'
        INCLUDE 'veget.h'
        INCLUDE 'comrunlabel.h'
        INCLUDE 'comemic.h'
!********************************************************************
! temporal var*/
        integer indxv
        real*8 tempor1,tempor2,db2,fd,dst,dd,nld,dstime
        real*8 dsd,temp_sg,temp_st
        real*8 fracgr(nlat,nlon), darea(nlat)

! calculation of current carbon cycle parameters

        call CCPARAM

! calculation of fraction dynamic variables

        fd=forshare_st-stR(lat,lon)
        dd=desshare_st-sdR(lat,lon)
        nld=nlshare_st-snltR(lat,lon)
        temp_st=stR(lat,lon)
        temp_sg=sgR(lat,lon)


! calculation of forest dynamics; exponential filtre
!-AM    dst=forshare_st-fd*exp(-1./t2t)-st(lat,lon)
        dst=fd*(1.d0-exp(-1./t2t))
!
!
         stR(lat,lon)=stR(lat,lon)+dst
!
        if (stR(lat,lon).lt.0.) stR(lat,lon)=0.
!
        snltR(lat,lon)=nlshare_st-nld*exp(-1./t2t)

! desert dynamics; exponential filtre
        dsd=desshare_st-dd*exp(-1./t2g)-sdR(lat,lon)
        tempor1=sdR(lat,lon)+dsd+stR(lat,lon)

! calculation of characteristic time of desert propagation
        if (tempor1.gt.0.9) then
             dstime=t2g*(1-tempor1)*10.+t2t*(tempor1-0.9)*10.
             dsd=desshare_st-dd*exp(-1./dstime)-sdR(lat,lon)
        endif

        sdR(lat,lon)=sdR(lat,lon)+dsd
        if (sdR(lat,lon).lt.0) sdR(lat,lon)=0.

        sgR(lat,lon)=1.-stR(lat,lon)-sdR(lat,lon)
        if (sgR(lat,lon).lt.0) sgR(lat,lon)=0.



        return
        end
