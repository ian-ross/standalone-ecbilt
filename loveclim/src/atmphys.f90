!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iatmphys
!-----------------------------------------------------------------------
! *** initializes variables used in subroutines of atmphys.f
!-----------------------------------------------------------------------
      USE NETCDF    ! netcdf module to read forcing data

      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'netcdf.inc'

      integer ios

      integer i,j,k,l,ireg,im,nn,is,j1,i1,ii,jj,ism
      real*8  beta,draganr,draglar,dum(2),asum,spv
      integer jyear,kyear,ilat,jmonth,m,status
      real*8 ksw
      character*6 numyear
      character*3 numday
      integer tmp_imonth
      real*8 globalmean

      !write(numyear,'(i6.6)') irunlabel+int((irunlabeld)/360)
      !write(numday,'(i3.3)') mod(irunlabeld,360)
      write(numyear,'(i6.6)') irunlabel
      write(numday,'(i3.3)') irunlabeld
!
! *** initial atmospheric temperatures and moisture
!
      if (irunlabel .eq. 0) then
        tempm(0)=tzero-35d0
        tempm(1)=tzero-35d0
        tempm(2)=tzero-8d0
        do j=1,nlon
          do i=1,nlat
            temp0g(i,j)=tempm(0)
            temp2g(i,j)=tempm(1)
            temp4g(i,j)=tempm(2)
            tempsg(i,j)=290d0
            rmoisg(i,j)=0d0
            relhum(i,j)=0d0
            q10(i,j)=0d0
            qsurf(i,j)=0d0
            do nn=1,ntyps
              q10n(i,j,nn)=0.0d0
              qsurfn(i,j,nn)=0.0d0
              tempsgn(i,j,nn)=290d0
              pgroundn(i,j,nn)=p0
            enddo
            pground(i,j)=p0
            geopg(i,j,2)=0d0
            tcc(i,j)=ccisccp(i,j,1)
          enddo
        enddo


        do j=1,nlon
          do i=1,nlat
            tsurfn (i,j,noc)= tempsg(i,j)
            tsurfn (i,j,nse)= tzero
            tsurfn (i,j,nld)= tempsg(i,j)
            tsurf (i,j)= tempsg(i,j)
          enddo
        enddo
        call atmphyszero
      else
	ios=0
        open(iuo+95,file='startdata/inatphy'//numyear//'_'//numday//'.dat', &
     &        form='unformatted')
        read(iuo+95) tsurfn,tempm,temp0g
        read(iuo+95) rmoisg,torain,tosnow
      endif

      do j=1,nlon
        do i=1,nlat
          rmountn(i,j,nld)=rmount(i,j)
          rmountn(i,j,nse)=0d0
          rmountn(i,j,noc)=0d0
          if (rmountn(i,j,nld).lt.0d0) rmountn(i,j,nld)=0d0
          if (fractn(i,j,nld).lt.epss) then
            rmountn(i,j,nld)=0d0
            qmount(i,j)=0d0
          else
            qmount(i,j)=rmountn(i,j,nld)
          endif
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
          asum=0d0
          do i1=-1,1
            do j1=-1,1
              ii=i+i1
              jj=j+j1
              if (ii.lt.1) then
                ii=1
                jj=jj+nlon/2
              endif
              if (ii.gt.nlat) then
                ii=nlat
                jj=jj+nlon/2
              endif
              if (jj.lt.1) jj=jj+nlon
              if (jj.gt.nlon) jj=jj-nlon
              asum=asum+rmountn(ii,jj,nld)
            enddo
          enddo
          qmount(i,j)=asum/9d0
        enddo
      enddo


!***  longwave radiation parameterisation, based on a linearisation
!***  of KRCM with respect to reference T and q profiles and other
!***  variables (greenhousegases)

      read(iuo+16) irn,ipl,pisccp,pncep,z500ncep
      read(iuo+16) tncep,qancep,ghgipcc,ccisccp
      read(iuo+16) lwrref

      read(iuo+17) lwrt,lwrts,lwrqts,lwrqa,lwrghg

!**   amplifcation of freshwater feedback
      lwrqa(:,:,:,:)=lwrqa(:,:,:,:)*AMPWIR
!**   amplifcation of feedback at the equator
      lwrqa(:,9,:,:)=lwrqa(:,9,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,10,:,:)=lwrqa(:,10,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,11,:,:)=lwrqa(:,11,:,:)*AMPEQIR
      lwrqa(:,12,:,:)=lwrqa(:,12,:,:)*AMPEQIR
      lwrqa(:,13,:,:)=lwrqa(:,13,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,14,:,:)=lwrqa(:,14,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,23,:,:)=lwrqa(:,23,:,:)*AMPANIR
      lwrqa(:,1,:,:)=lwrqa(:,1,:,:)*AMPANIR2
      lwrqa(:,2,:,:)=lwrqa(:,2,:,:)*AMPANIR2
      lwrqa(:,3,:,:)=lwrqa(:,3,:,:)*AMPANIR2

!     write(iuo+99,*) "Modif IR scheme"
!     write(iuo+99,*) AMPWIR,AMPEQIR,expIR
!     write(iuo+99,*) HPROFW,HPROFTROP,1.0+(HPROFTROP-1.0)/2.0,HPROFEQ

!***  update moisture profile used in the linearization of the
!***  radiative scheme in the tropics:easy surrogate to a change
!***  mean IR flux in the model without affecting sensitivity
      do ism=1,12
        qancep(1,ism)=qancep(1,ism)*HPROFAN2
        qancep(2,ism)=qancep(2,ism)*HPROFAN2
        qancep(3,ism)=qancep(3,ism)*HPROFAN2
        qancep(4,ism)=qancep(4,ism)*HPROFW
        qancep(5,ism)=qancep(5,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(6,ism)=qancep(6,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(7,ism)=qancep(7,ism)*HPROFTROP
        qancep(8,ism)=qancep(8,ism)*HPROFTROP
        qancep(9,ism)=qancep(9,ism)*HPROFTROP
        qancep(10,ism)=qancep(10,ism)*HPROFTROP

        qancep(11,ism)=qancep(11,ism)*HPROFEQ
        qancep(12,ism)=qancep(12,ism)*HPROFEQ

        qancep(13,ism)=qancep(13,ism)*HPROFTROP
        qancep(14,ism)=qancep(14,ism)*HPROFTROP
        qancep(15,ism)=qancep(15,ism)*HPROFTROP
        qancep(16,ism)=qancep(16,ism)*HPROFTROP
        qancep(17,ism)=qancep(17,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(18,ism)=qancep(18,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(19,ism)=qancep(19,ism)*HPROFW
        qancep(20,ism)=qancep(20,ism)*HPROFW
        qancep(21,ism)=qancep(21,ism)*HPROFW
        qancep(22,ism)=qancep(22,ism)*HPROFW
        qancep(23,ism)=qancep(23,ism)*HPROFAN
        qancep(24,ism)=qancep(24,ism)*HPROFW
        qancep(25,ism)=qancep(25,ism)*HPROFW
        qancep(26,ism)=qancep(26,ism)*HPROFW
        qancep(27,ism)=qancep(27,ism)*HPROFW
      enddo


!read sulfates optical depths
      status = nf_open("inputdata/SUL.nc", nf_nowrite, ireg)
      status = nf_inq_varid(ireg, "Sul", j)
      status = nf_get_vara_real(ireg, j, (/1,1/), (/64,32/), sulopt)

      call ghgupdate

!* set reference value for CO2 concentration

      do ireg=1,27
         do im=1,12
            beta=(tncep(11,ireg,im)-tncep(10,ireg,im))/alog(400./300.)
            tncep12(1,ireg,im)=tncep(10,ireg,im)+beta*alog(350./300.)
            beta=(tncep(14,ireg,im)-tncep(13,ireg,im))/alog(700./600.)
            tncep12(2,ireg,im)=tncep(13,ireg,im)+beta*alog(650./600.)
         enddo
      enddo

      do k=1,17
         rlogtl(k)=log(pncep(k)/65000.)
      enddo

      do ireg=1,27
         rlogts(ireg)=log(pisccp(ireg)/65000.)
      enddo

! *** shortwave radiation parameters

      read(iuo+18) costref,salbref
      read(iuo+18) swrref
      read(iuo+19) swrcost,swrsalb

      call detqmtabel

310   format(9f7.2)
330   format(2f5.2)
340   format(f5.2)

! *** computation of the effective turning angle

      draglar=3.141592654/180.0*dragla
      draganr=3.141592654/180.0*dragan
      do i=1,nlat
        if (phi(i).lt.(-1.0*draglar)) then
           dragane(i)=-1.0*draganr
        else
           if (phi(i).gt.draglar) then
               dragane(i)=draganr
           else
               dragane(i)=draganr*phi(i)/draglar
           endif
        endif
      enddo

! *** evaporation factor

      evfac=1d0
      ksw=0.042

      solarcl = solarm + tsi
      solarcl(1:10)  = solarcl(1:10)  + solarvol(imonth,1)
      solarcl(11:16) = solarcl(11:16) + solarvol(imonth,2)
      solarcl(17:22) = solarcl(17:22) + solarvol(imonth,3)
      solarcl(23:32) = solarcl(23:32) + solarvol(imonth,4)
      IF (o3 /= 25.0) THEN
         solarcl(1:16)  = solarcl(1:16)  + (4 * 0.247 * ksw * (o3 - 25.0))
         solarcl(17:32) = solarcl(17:32) + (4 * 0.324 * ksw * (o3 - 25.0))
      endif

      RETURN
    END SUBROUTINE iatmphys

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ghgupdate
!-----------------------------------------------------------------------
! *** Setup from GHG concentrations
!-----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'comatm.h'
      INCLUDE 'comphys.h'
      INCLUDE 'comemic.h'

      INTEGER s, r, k, l, m, h
      REAL*8 logco2, sqrch4, sqrn2o
      REAL*8 alpho3lw(2)

      PGACO2 = ghg(1)
!*** Update LW reference radiation fluxes using new GHG concentrations
      logco2 = LOG(ghg(1) / ghgipcc(1))
      IF (ghg(1) == 0) STOP
      sqrch4 = SQRT(ghg(2)) - SQRT(ghgipcc(2))
      sqrn2o = SQRT(ghg(3)) - SQRT(ghgipcc(3))
      alpho3lw(1) = 153.6
      alpho3lw(2) = 201.2
      DO h = 1, 2
         DO l = 0, 1
            DO s = 1, 4
               DO r = 1, 27
                  DO k = 1, 7
                     lwrflux(k,r,s,l,h) = lwrref(k,r,s,l) + &
                          & lwrghg(k,1,r,s,l) * logco2 + &
                          & lwrghg(k,2,r,s,l) * sqrch4 + &
                          & lwrghg(k,3,r,s,l) * sqrn2o
                     DO m = 4, 19
                        lwrflux(k,r,s,l,h) = lwrflux(k,r,s,l,h) + &
                             & lwrghg(k,m,r,s,l) * (ghg(m) - ghgipcc(m))
                     END DO
                     lwrflux(k,r,s,l,h) = lwrflux(k,r,s,l,h) + &
                          & lwrghg(k,4,r,s,l) * alpho3lw(h) * (o3 - 25.0)
                  END DO
               END DO
            END DO
         END DO
      END DO
    END SUBROUTINE ghgupdate



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmphyszero
!-----------------------------------------------------------------------
! *** initializes data arrays to zero for physics routines
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comrunlabel.h'

      integer i,j


      do j=1,nlon
        do i=1,nlat
          cormois(i,j)=0.d0
          torain(i,j)=0.d0
          tosnow(i,j)=0.d0
          dyrain(i,j)=0.d0
          corain(i,j)=0.d0
          dysnow(i,j)=0.d0
          cosnow(i,j)=0.d0
          thforg0(i,j)=0.d0
          thforg1(i,j)=0.d0
          thforg2(i,j)=0.d0
          vhforg0(i,j)=0.d0
          vhforg1(i,j)=0.d0
          vhforg2(i,j)=0.d0
        enddo
      enddo


      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sensrad
!-----------------------------------------------------------------------
! *** computes atmospheric forcing due to sensible heat and radiation
! *** the forcing terms are computed in dimensional units (K/s)
! *** input
! ***       hflux : sensible heat flux between atmosphere and earth
! ***       hesws : short wave solar radiation in layer 1 or 2
! ***       ulrad : upward long wave radiation in layer 1 or 2
! *** output
! ***       thforg: temperature forcing in layer 1 or 2 in K/s
! ***       vhforg: diabatic forcing in layer 1 or 2 in K/s
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  halpha,halpha1,halpha2,sum1,sum2,sum0,halpha0

      halpha=grav/cpair


      halpha0 =halpha/dp0
      halpha1 =halpha/dp1
      halpha2 =halpha/dp2



!
! *** summation of forcing terms
!
      do j=1,nlon
        do i=1,nlat


          sum0 = (hesw0(i,j) - ulrad0(i,j) + ulrad1(i,j))*halpha0
          sum1 = (hesw1(i,j) - ulrad1(i,j) + ulrad2(i,j))*halpha1
          sum2 = (hflux(i,j) + hesw2(i,j) - ulrad2(i,j) + &
               & ulrads(i,j) - dlrads(i,j))*halpha2



          thforg0(i,j) = thforg0(i,j) + sum0
          thforg1(i,j) = thforg1(i,j) + sum1
          thforg2(i,j) = thforg2(i,j) + sum2


          vhforg0(i,j) = vhforg0(i,j) + sum0
          vhforg1(i,j) = vhforg1(i,j) + sum1
          vhforg2(i,j) = vhforg2(i,j) + sum2


        enddo
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine solar(istep)
!-----------------------------------------------------------------------
! Calculates incoming solar radiation as a function of latitude
! for each day of the year, given the orbital parameters (see PMIP)
! One year has 360 days.  Reference: A. Berger, JAS, 35, 2362-2367,1978
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'

      integer i,j,l,NVE


      real*8 beta,alam,alam0,ala,ala0
      real*8 fac1,fac2,fac3,ro,roref
      real*8 deltal, sindl, cosdl, tandl
      real*8 rkosz, rkosha1, ha1
      real*8 deg2rad, day2rad
      real*8 solard, solarcf(nlat)
      real*8 ksw
      real*8 alpho3sw(2)
      integer istep

      deg2rad=pi/180.d0
      day2rad=pi/180.d0


! Present-day orbital parameters: eccentricity ecc, obliquity obl and
! angle om between Vernal Equinox and Perihelion (angles all given
! in degrees and converted to radians). Solarc is the solar constant.
! NVE is day of the Vernal Equinox, set at 21 MARCH
! Implementatie van Nanne

      ! ecc=0.016724
      ! obl=23.446*deg2rad
      ! omweb=(102.04+180.00)*deg2rad
      ecc = eccf
      obl = oblf * deg2rad
      omweb = (omwebf + 180.00) * deg2rad
      NVE=30+30+21
!
! In old SW-routine of ECBilt-model values were as follows:
!
!      ecc=0.0
!      solarc=1353.
!      NVE=90
!
! At 6000 years BP (Before Present) values were as follows:
!
!      ecc=0.018682
!      obl=24.105*deg2rad
!      omweb=(0.87+180.00)*deg2rad
!
! :0:
!     ecc=0.018994
!     obl=22.949*deg2rad
!     omweb=(114.42+180.00)*deg2rad


! First compute alam0 (the starting point). Then convert days to
! the true longitude ala, using Berger's formulae.
! Longitude ala loops from 0 (Vernal Equinox) to 359, ro is earth-sun
! distance relative to the major axis, del is declination.
!
      ala0=0.
      beta=(1.-ecc**2.)**0.5
      fac1=(0.5*ecc+0.125*ecc**3.)*(1.+beta)*sin(ala0-omweb)
      fac2=0.25*ecc**2.*(0.5+beta)*sin(2.*(ala0-omweb))
      fac3=0.125*ecc**3.*(1./3.+beta)*sin(3.*(ala0-omweb))
      alam0=ala0-2.*(fac1-fac2+fac3)


      l=(imonth-1)*30+iday-NVE
      if (l.lt.0) l=l+360
      alam=alam0+l*1.0*pi/180.


      fac1=(2.*ecc-0.25*ecc**3.)*sin(alam-omweb)
      fac2=1.25*ecc**2.*sin(2.*(alam-omweb))
      fac3=(13./12.)*ecc**3.*sin(3.*(alam-omweb))
      ala=alam+fac1+fac2+fac3
      ro=(1.-ecc**2.)/(1.+ecc*cos(ala-omweb))
      deltal=asin(sin(obl)*sin(ala))


      sindl=sin(deltal)
      cosdl=cos(deltal)
      tandl=tan(deltal)


! factor voor variable afstand Aarde-Zon (Berger, p.2367; Velds, p. 99)
      solard=1./ro**2.
!     solardref=1./roref**2.
      solardref=solard
      ksw=0.042
      alpho3sw(1)=0.247
      alpho3sw(2)=0.324
      solarc = solarm + tsi
      solarcl(:) = solarc
      solarcl(1:10)  = solarcl(1:10)  + solarvol(imonth,1)
      solarcl(11:16) = solarcl(11:16) + solarvol(imonth,2)
      solarcl(17:22) = solarcl(17:22) + solarvol(imonth,3)
      solarcl(23:32) = solarcl(23:32) + solarvol(imonth,4)
      IF (o3 /= 25.0) THEN
        solarcl(1:16)  = solarcl(1:16)  + (4 * alpho3sw(1) * ksw * (o3 - 25.0))
        solarcl(17:32) = solarcl(17:32) + (4 * alpho3sw(2) * ksw * (o3 - 25.0))
      END IF
      !write(*,*) day, iatm, nstpyear, initialization
      if((mod(nint(day*real(iatm)),30*iatm).eq.0).or. &
           & (initialization.eqv..true.)) then
       write(iuo+99,12) 'vol forcing ',iyear,imonth, &
            & solarcl(15), solarvol(imonth,:)
      endif
11      format(A12,2i6,2f12.3)
12      format(A12,3i6,5f12.3)
      do i=1,nlat
! zonneconstante is 1370 in sw parameterisatie
       solarcf(i)=solarcl(i)/1370.d0
! beide effecten samen
       solarf(i)=solarcf(i)*solard
      enddo


      do i=1,nlat
         rkosha1=-tanfi(i)*tandl
         rkosha1=sign(min(abs(rkosha1),1.d0),rkosha1)
         ha1=acos(rkosha1)
         rkosz=sinfi(i)*sindl*(ha1 - tan(ha1))/pi
         if (ha1 .ne. 0.d0) then
            kosz(i) = rkosz*pi/ha1
         else
            kosz(i) = rkosz
         endif
         dayfr(i) = ha1/pi
         q0(i)=rkosz*solarcl(i)*solard
      enddo


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fluxes(nn)
!-----------------------------------------------------------------------
! *** computes energy fluxes above ocean surface
! *** short wave radiation, long wave radiation, sensible heat flux
! *** latent heat flux, evaporation
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comsurf.h'
      include 'comrunlabel.h'
      include 'comdiag.h'

      integer nn


      if (irad.eq.1) call swaverad2(nn)
      call swaverad(nn)
      call lwaverad(nn)
      if (irad.eq.1) call lwaverad2(nn)
      call dragcoef(nn)
      call surfmois(nn)
      call sensibheat(nn)
      call latentheat(nn)
      if (nn.eq.noc) call momentflux(nn)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine swaverad(nn)
!-----------------------------------------------------------------------
! *** computes short wave radiation
! *** linearization of RCM with ISCCP D2 1990 clouds
!-----------------------------------------------------------------------



      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'comdiag.h'


      integer i,j,k,l,ireg
      integer m, d, r, nn , nol


      real*8 f0,f1,ftot(8),fn(8,0:1)
      real*8 drs, drs2, drs3
      real*8 dcost, df,sk,sr,x,y,dfs,smsc
      real*8 fswdtoa,fswutoa,fswdsfc(nlat,nlon),fswusfc
      real*8 fswutoa2(nlat,nlon),fswutoa0(nlat,nlon)
      real*8 fswdtoa0(nlat,nlon)
      real*8 globalmean
      real*8 fswutoaG0,fswutoaGA,fswdtoa2,fswdtoaG0,fswdtoaGA
      real*8 fswutoa_diff,fswutoaG,df_test,fswdtoa_diff,fswdtoaG


      integer nreg(2)
      real*8 zac(2),asup
!     real*8 zac(2),asup,bup
      common /rad_sul0 /fswutoaG,df_test,fswdtoaG
      common /rad_sul1 /fswutoa0,fswutoaGA,fswdtoa0,fswdtoaGA

! *** aerosol scattering included as a correction on the upward
! *** clear sky fluxes
! *** sk,sr: empirical coefficients Dorland et al, J. Geophys. Res.,102,
! *** 28079-28100, 1997.
! *** smsc: mass scattering coefficient [m2/g]
! *** dso4: change in sulfate aerosol column integrated concentration since
! *** pre-industrial times [g/m2]


      sk=0.058d0*1370d0
      sr=0.05d0
      smsc=8.0
!     write(iuo+99,*) 'bup=',bup

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2

      do j=1,nlon
       do i=1,nlat
         tas1(i,j) = sulopt(j,i)
         nreg(nol)=irn(i,j,nol)
         zac(nol)=dble(costref(nreg(nol),imonth))
         if (zac(nol).GT.0.) then
           asup = (bup*tas1(i,j)*(1-albesn(i,j,nn))**2)/zac(nol)
         else
           asup =0.
         endif

         alb2esn(i,j,nn) = albesn(i,j,nn)+ asup
!        else
!          alb2esn(i,j,nn) = albesn(i,j,nn)
!        endif
          if (alb2esn(i,j,nn).ge.1.) then
            alb2esn(i,j,nn)=1.
          endif
          df=dayfr(i)*solarf(i)
          ireg=irn(i,j,nol)
          dcost=kosz(i)-costref(ireg,imonth)
          do l=1,8
            do k=0,1
              fn(l,k) = swrref(l,ireg,imonth,k) &
     &               +  swrcost(l,ireg,imonth,k)*dcost
            enddo
          enddo

          x=sqrt(kosz(i))
          if(kosz(i).lt.0.) x=0d0
          y=sqrt(1-alb2esn(i,j,nn))
          dfs=sk*(4d0*x*y*(y-x)-sr)*dso4(i,j)*smsc
!         WRITE(*,*)dso4(i,j)
          if (dfs.gt.0d0.and.kosz(i).lt.0.05) dfs=0d0
          drs=alb2esn(i,j,nn)-salbref(ireg,imonth)
          drs2=drs*drs
          drs3=drs2*drs


          do l=1,4
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs+dfs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
!           write(*,*)f0,f1
!           if (l.eq.1)ftot_test=f0
          enddo
          do l=5,8
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo


! alternative calculation of upward flux at ground:
! in parameterisation no cross terms are accounted for, which are important for
! upward shortwave radiation at surface and therefore also for net flux
! heswsn(i,j)

          ftot(4)=-alb2esn(i,j,nn)*ftot(8)

          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*df
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*df
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*df
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*df

          if (irad.eq.1) then
! for diagnostic purposes:
! (1) downward shortwave radiation at TOA
             if (nn.eq.1)fswdtoa0(i,j)=0.
             fswdtoa2=-ftot(5)*df
             fswdtoa0(i,j)=fswdtoa0(i,j)+(fractn(i,j,nn)*fswdtoa2)
! (2) upward shortwave radiation at TOA
             if (nn.eq.1)fswutoa0(i,j)=0.
             fswutoa2(i,j)=ftot(1)*df
             fswutoa0(i,j)=fswutoa0(i,j)+(fractn(i,j,nn)*fswutoa2(i,j))
! (3) downward shortwave radiation at SURFACE
!            fswdsfc(i,j)=-ftot(8)*df
! (4) upward shortwave radiation at SURFACE
!            fswusfc=(heswsn(i,j)+ftot(8)*df)

         endif
        enddo
      enddo

      if (irad.eq.1) then
       if (nn.eq.3) then
        fswutoaG0=globalmean(fswutoa0)
        fswdtoaG0=globalmean(fswdtoa0)
        fswutoa_diff=fswutoaG0-fswutoaG
        fswdtoa_diff=fswdtoaG0-fswdtoaG
        !if (iyear.eq.0) fswutoaGA=0.
        !if (iyear.eq.0) fswdtoaGA=0.
        if (initialization.eqv..true.) fswutoaGA=0.
        if (initialization.eqv..true.) fswdtoaGA=0.
        fswutoaGA=fswutoaGA+(fswutoa_diff/(360.*6.))
        fswdtoaGA=fswdtoaGA+(fswdtoa_diff/(360.*6.))
       endif
      else
       fswutoaGA=0.
       fswdtoaGA=0.
      endif



      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lwaverad(nn)
!-----------------------------------------------------------------------
! *** computes long wave radiation according to the parameterization of
! *** Chao Chou and Neelin and substantially adapted and extended
! *** for global scale and more
! *** specific ECBILT application by the one and only Michiel Schaeffer
! ***
! *** parameters: nlat   = number of gridpoints in the meridional
! ***                      direction (32)
! ***             nlon   = number of gridpoints in the zonal
! ***                      direction (64)
! ***
! *** input : dtemp(19,nlat,nlon): temperature anomalies [K] wrt ncep
! ***                              climatology tncep in common lwrscheme
! ***         dqa(nlat,nlon) : anomalies of total prec. water cont. below
! ***                          500 hPa wrt ncep climatology
! ***         tcc(nlat,nlon)  : total cloud cover
! ***         ghg(19) : concentrations of well mixed ghg's (see comphys.h)
! ***
! *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             toa
! ***          ulrad2(nlat,nlon): net longwave radiation [Wm-2] at
! ***                             (500 hPa)
! ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
! ***                             the surface
! ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             the surface
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'comdiag.h'


      integer i,j,l,k,m,is,ism,nol,nn,ireg,h
      real*8  lwr(7,0:1),dumts
      real*8  dqa,dqreg(27)
      real*8  ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm
      real*8  globalmean
      real*8  ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  ulrad2nz(nlat,nlon),ulradsnz(nlat,nlon)
      real*8  dlradsnz(nlat,nlon), ulrad0nUz(nlat,nlon)
      real*8  ulrad1nUz(nlat,nlon)
      common / radO3 / ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      common / radO32 / ulrad0nUz,ulrad1nUz



      is=imonth/3+1
      if (is.gt.4) is=1
      ism=(is-1)*3+1

      do i=1,27
!dqa    dqreg(i)=qancep(i,ism)**0.3333
        dqreg(i)=qancep(i,ism)
      enddo

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
        do i=1,nlat
          ireg=irn(i,j,nol)

!-hemispheric dependence of tropospheric ozone forcing
          if (i.le.16) then
           h=1
          else
           h=2
          endif

!dqa      dqa=lwrmois(i,j)-dqreg(ireg)
!dqa      q**1/3-qm**1/3=qm**(1/3-n)*(q**n-qm**n)
          dqa=dqreg(ireg)**(0.3333-EXPIR)* &
               & (lwrmois(i,j)**EXPIR-dqreg(ireg)**EXPIR)
!dqa      write(iuo+99,'(i3,3F12.5)') i,lwrmois(i,j)-
!dqa &      dqreg(ireg),dqa,lwrmois(i,j)**0.333-dqreg(ireg)**0.33333
          do l=0,1
            do k=1,7
              lwr(k,l)=lwrflux(k,ireg,is,l,h)+lwrqa(k,ireg,is,l)*dqa
!                   write(*,*) "lwr=",lwrflux(k,ireg,is,l,h),"+",lwrqa(k,ireg,is,l),"*",dqa
              do m=1,ipl(ireg)-1
                lwr(k,l)=lwr(k,l)+lwrt(k,m,ireg,is,l)*dtemp(m,i,j,nol)
              enddo
              lwr(k,l)=lwr(k,l)+lwrt(k,18,ireg,is,l)*dtemp(18,i,j,nol)
            enddo

            dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism)
            do m=1,4
              do k=1,3
                lwr(k,l)=lwr(k,l)+ &
                     & (lwrts(k,m,ireg,is,l)+lwrqts(k,m,ireg,is,l)*dqa)*dumts
              enddo
              lwr(7,l)=lwr(7,l)+ &
                   & (lwrts(7,m,ireg,is,l)+lwrqts(7,m,ireg,is,l)*dqa)*dumts
              dumts=dumts*(tsurfn(i,j,nn)-tncep(19,ireg,ism))
            enddo


          enddo

          ulrad0n(i,j,nn)=(lwr(1,0)*(1-tcc(i,j))+lwr(1,1)*tcc(i,j))
          ulrad1n(i,j,nn)=(lwr(2,0)+lwr(5,0))*(1-tcc(i,j)) + &
               & (lwr(2,1)+lwr(5,1))*tcc(i,j)
          ulrad2n(i,j,nn)=(lwr(3,0)+lwr(6,0))*(1-tcc(i,j)) + &
               & (lwr(3,1)+lwr(6,1))*tcc(i,j)

	  ulradsn(i,j,nn)=emisn(nn)*sboltz*tsurfn(i,j,nn)**4
          dlradsn(i,j,nn)=-lwr(7,0)*(1-tcc(i,j))-lwr(7,1)*tcc(i,j)


         if (irad.eq.1) then
          if(nn.eq.1) then
            ulrad0nUz(i,j)=0.
            ulrad1nUz(i,j)=0.
           endif

          ulrad0nUz(i,j)=ulrad0nUz(i,j)+(ulrad0n(i,j,nn)*fractn(i,j,nn))
          ulrad1nUz(i,j)=ulrad1nUz(i,j)+(ulrad1n(i,j,nn)*fractn(i,j,nn))

         endif
       enddo
      enddo

      if (irad.eq.1) then
       ulrad0nU=globalmean(ulrad0nUz)
       ulrad1nU=globalmean(ulrad1nUz)
      endif

!     ulrad0nm=globalmean(ulrad0n)
!     ulrad1nm=globalmean(ulrad1n)
!     ulrad2nm=globalmean(ulrad2n)
!     ulradsnm=globalmean(ulradsn)
!     dlradsnm=globalmean(dlradsn)

!     write(iuo+37,*)ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm

! *** that's all folks


      return
      end


!2345678901234567890123456789012345678901234567890123456789012345678901
      subroutine dragcoef(nn)
!----------------------------------------------------------------------
! *** drag coefficient
! *** depending on stability
!----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comsurf.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,nn


      real*8 cdrags,cdragl,tdif,cred,cdum

      cdrags=cdrag
      cdragl=cdrag
      cred=0.2
      cdum=(1-cred)*0.5
      do j=1,nlon
        do i=1,nlat
          cdragw(i,j)=cwdrag

          tdif=tsurfn(i,j,nn)-tempsgn(i,j,nn)
          cdragvn(i,j,nn)=cdrags*max(cred,cred+min(cdum+cdum*tdif,1-cred))

        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sensibheat(nn)
!-----------------------------------------------------------------------
! *** sensible heatflux between surface and atmosphere
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn


! *** sensible heatflux (watt/m**2)
      if (nn.eq.2) then
       do j=1,nlon
        do i=1,nlat
          hficof(i,j)=alphad*cdragvn(i,j,2)*uv10(i,j)
        enddo
       enddo
       do j=1,nlon
        do i=1,nlat
          hfluxn(i,j,2)=hficof(i,j)*(tsurfn(i,j,2)-tempsgn(i,j,2))
        enddo
       enddo
      else
       do j=1,nlon
        do i=1,nlat
           hfluxn(i,j,nn)=alphad*cdragvn(i,j,nn)*uv10(i,j)* &
     &                      (tsurfn(i,j,nn)-tempsgn(i,j,nn))
        enddo
       enddo
      endif


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine latentheat(nn)
!-----------------------------------------------------------------------
! *** latent heatflux between surface and atmosphere
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn,n,NSTAT
      real*8  qsat,db,emois,esubf,evapf,esnow,sfrac,edum,psilai

      real*8    bmoisg(nlat,nlon),resist(3),lai(2),k0(2)
      real*8    bmoism(nlat,nlon),rs
      real*8    rainf(nlat,nlon),snowf(nlat,nlon)
      real*8    fswdsfc(nlat,nlon)
      real*8    dc,eflux_t,eflux_g,eflux_bare
      real*8    fswdsfcM

      real*8 st,sg,sd,snlt,anup,blai,pnpp,b12,b34,b1,b2,b3,b4, &
     &       anup_moy,stock,st_moy
      real*8 stR,sgR,sdR,snltR

      common /lbmbmois/ bmoisg,bmoism,rainf,snowf
      common /pr_evap /fswdsfc
      COMMON /BIOTA/ &
     &   ST(nlat,nlon), SG(nlat,nlon), SD(nlat,nlon), SNLT(nlat,nlon), &
     &   BLAI(nlat,nlon,2), PNPP(nlat,nlon), &
     &   B12(nlat,nlon),   B34(nlat,nlon), &
     &   B1(nlat,nlon), B2(nlat,nlon), B3(nlat,nlon), B4(nlat,nlon), &
     &   ANUP_MOY(nlat,nlon),ANUP(nlat,nlon), STOCK(nlat,nlon), &
     &   st_moy(nlat,nlon), &
     &   NSTAT(nlat,nlon),STR(nlat,nlon), SGR(nlat,nlon),SDR(nlat,nlon), &
     &   SNLTR(nlat,nlon)


! *** latent heatflux due to evaporation from surface (watt/m**2)
! *** and evaporation rate (m/s)
! *** limit evaporation to available snow or bottom moisture
! *** evaporation factor =1 over snow, over wet land maximal 1


      if (nn.ne.nld) then
        do j=1,nlon
          do i=1,nlat
            evfacan(i,j,nn)=evfac
            efluxn(i,j,nn)=alphav*cdragvn(i,j,nn)*uv10(i,j)* &
     &                 (qsurfn(i,j,nn)-q10n(i,j,nn))
            efluxn(i,j,nn)=evfacan(i,j,nn)*max(0.d0,efluxn(i,j,nn))
            evapn(i,j,nn)=efluxn(i,j,nn)/(rowat*rlatvap)
          enddo
        enddo
      else
        do j=1,nlon
          do i=1,nlat



! *** evaporation factor =1 over snow, over wet land maximal 1
! *** care has to be taken in case of a snowcover in order to
! *** conserve heat: the evaporation is constraint to the amount
! *** of snow and moisture available; it can happen that in
! *** one timestep, the remaining snow is sublimated and part
! *** of the latent heat flux is also used to evaporate bottom
! *** moisture: Eflux = (1-sfrac)*Esub + sfrac*Evap


            if (fractn(i,j,nld).gt.epss) then


              if (adsnow(i,j).gt.0.) then
                evfacan(i,j,nld)=evfac
                edum=cdragvn(i,j,nld)*uv10(i,j)*(qsurfn(i,j,nld)-q10n(i,j,nld))
                edum=evfacan(i,j,nld)*max(edum,0d0)
                esubf=alphas*edum
                evapf=alphav*edum
                esnow=min(rowat*adsnow(i,j)*rlatsub/dtime,esubf)
                if (esnow.lt.esubf) then
                  sfrac=(esubf-esnow)/esubf
                  emois=min(rowat*abmoisg(i,j)*rlatvap/dtime,sfrac*evapf)
                  efluxn(i,j,nld)=esnow+emois
                  evapn(i,j,nld)=esnow/(rowat*rlatsub)+emois/(rowat*rlatvap)
                else
                  efluxn(i,j,nld)=esubf
                  evapn(i,j,nld)=esubf/(rowat*rlatsub)
                endif
              else
                evfacan(i,j,nld)=evfac*min(1d0,abmoisg(i,j)/max(1.0E-10,abmoism(i,j)))
                dc=qsat(pgroundn(i,j,nld),tempsgn(i,j,nld))-q10n(i,j,nld)
                psilai=0.2+(0.08*min(max(tempsgn(i,j,nld)-tzero,0.),10.))
                lai(1)=6*psilai
                lai(2)=2*psilai
                k0(1)=30.0E-5
                k0(2)=25.0E-5
                fswdsfcM=max(fswdsfc(i,j),1.0d0)
                do n=1,2
                rs=((fswdsfcM+125)/fswdsfcM)*(23E-3+(1.5*dc))*(1/(lai(n)*k0(n)))
                rs=max(0d0,rs)
                resist(n)=rs+(1/(cdragvn(i,j,nld)*uv10(i,j)))
                resist(n)=1/resist(n)
                enddo
                resist(3)=1/(cdragvn(i,j,nld)*uv10(i,j))
                resist(3)=1/resist(3)
!	write(*,*)cdragvn(i,j,nld)*uv10(i,j), rs(i,j)

!               efluxn(i,j,nld)=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld))*
!    &              ((st(i,j)*resist(1))+(sg(i,j)*resist(2))+
!    &              (sd(i,j)*resist(3)))
                eflux_bare=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld)) &
     &           *cdragvn(i,j,nld)*uv10(i,j)*(10./30.)
                eflux_g=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld)) &
     &           *resist(2)*(15./30.)
                eflux_t=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld)) &
     &           *resist(1)*(30./30.)
                efluxn(i,j,nld)=eflux_bare+(sg(i,j)*eflux_g) &
     &          + (st(i,j)*eflux_t)


                efluxn(i,j,nld)=evfacan(i,j,nld)* &
     &                          max(efluxn(i,j,nld),0d0)
                efluxn(i,j,nld)=min(abmoisg(i,j)*rowat*rlatvap/dtime &
     &                       ,efluxn(i,j,nld))
                evapn(i,j,nld)=efluxn(i,j,nld)/(rowat*rlatvap)
              endif


            endif
          enddo
        enddo
      endif

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine momentflux(nn)
!-----------------------------------------------------------------------
! *** computation of windstress
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn
      real*8  uv,costt,sintt,facstr

      facstr=roair * uv10rws

      do j=1,nlon
        do i=1,nlat
          costt=cos(dragane(i))
          sintt=sin(dragane(i))
          winstu(i,j)=cdragw(i,j)*facstr*uvw10(i,j)* &
               & (utot(i,j,3)*costt-vtot(i,j,3)*sintt)
          winstv(i,j)=cdragw(i,j)*facstr*uvw10(i,j)* &
               & (utot(i,j,3)*sintt+vtot(i,j,3)*costt)
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine totwind10
!-----------------------------------------------------------------------
! *** computation of strength of 10-meter winds (uv10r* 800 hPa wind)
! *** input  u800, v800 , udivg, vdivg
! *** output uv10 strength of 10 m wind at gaussian grid with minimum
! ***        uvw10 strength of 10 m wind at gaussian grid for windstress
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,k
      real*8  uv

! *** bug fix 27 march 97: uv was declared integer

      do j=1,nlon
        do i=1,nlat
          uv=sqrt((utot(i,j,3))**2 + (vtot(i,j,3))**2)
          uv10(i,j)=uv10rfx*uv
          uvw10(i,j)=uv10rws*uv
! *** minimum value of uv10 set to uv10m m/s
          if (uv10(i,j).lt.uv10m) uv10(i,j)=uv10m
        enddo
      enddo


      return
      end
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine tracer
!-----------------------------------------------------------------------
! *** advection of tracer field
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  hdivmg(nlat,nlon)
      real*8  co2sp(nsh2)



      call rggtosp(co2,co2sp)
      call sptogg(co2sp,co2,pp)



! *** horizontal divergence of tracer


      call trafluxdiv(hdivmg,co2sp,co2)


!
! *** time stepping forward time stepping
!


      do j=1,nlon
        do i=1,nlat
          co2(i,j)=co2(i,j)-dtime*(hdivmg(i,j))
          if (co2(i,j).lt.0d0) co2(i,j)=0d0
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisfields
!-----------------------------------------------------------------------
! *** calculates relative humidity of the moised layer
! *** and specific humidity above the surface and at the surface
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn
      real*8  qsat,pmount,tmount,qmax,dqmdt

      do j=1,nlon
        do i=1,nlat


          call ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)


          if (qmax.gt.0d0) then


            relhum(i,j)=min(1d0,rmoisg(i,j)/qmax)


          else

            relhum(i,j)=0d0


          endif


          q10(i,j)= 0.d0
          do nn=1,ntyps
            q10n(i,j,nn)=relhum(i,j) * qsat(pgroundn(i,j,nn),tempsgn(i,j,nn))
          enddo

! *** lwrmois is used in the lwr parameterization

!dqa      lwrmois(i,j)=rmoisg(i,j)**0.3333
          lwrmois(i,j)=rmoisg(i,j)


        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine surfmois(nn)
!-----------------------------------------------------------------------
! *** calculates specific humidity at the surface
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'



      integer i,j,nn
      real*8  qsat


      do j=1,nlon
        do i=1,nlat


          qsurfn(i,j,nn)=qsat(pgroundn(i,j,nn),tsurfn(i,j,nn))

        enddo
      enddo


      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine cloud

!----------------------------------------------------------------------
! ***calculates the total cloud cover tcc using Axel's new scheme
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! *** calculates total cloud cover tcc, based on a global climatology
! *** from isccpd2. Michiel Schaeffer may 1998.
! *** diagnostic cloud scheme based on relative humidity and omega
! *** and stability. Jules Beersma
! *** depending on iradcloud, cloud climatology of diagnostic clouds are
! *** used in the calculation of radiative fluxes.
!-----------------------------------------------------------------------

      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j
      real*8 rhc, rhfac
      real*8 cc


      do j=1,nlon
        do i=1,nlat

           rhc = relhcrit
           rhfac=relhfac
! enhance clouds in case of vertical motion
           if (omegg(i,j,2) .lt.  0.0 )  rhfac=0.95d0
           if (omegg(i,j,2) .lt. -0.04)  rhfac=0.90d0
! enhance clouds in areas of subsidence inversions
           if ((fractoc(i,j) .gt. epss) .and. (omegg(i,j,2) .gt.  0.03)) &
     &         rhfac=0.7d0

           cc=(relhum(i,j)/rhfac - rhc)/(1.0d0 - rhc)
           cc=max(cc,0.0d0)
           cc=min(cc,1.0d0)


           tccd(i,j)=cc
           if (iradcloud.eq.1) then
             tcc(i,j)=tccd(i,j)
           else
             tcc(i,j)=ccisccp(i,j,imonth)
           endif
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine detqmtabel
!-----------------------------------------------------------------------
!***  calculate tabel of maximum content of water in the layer below
!***  500 hPa using the Clausius-Clapeyron relation and assuming a
!***  temperature profile linear in log(p): T(p)= Tr + alpha*log(p/pr)
!***  where alpha is given by (T350-T650)/log(350/650)
!***  for given groundtemperature and 650 and 350 hPa temperatures
!***  the maximum water content in [m] is given by qmtabel
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k
      real*8  tmount,t500,b,qmax,expint,z1,z2,bz1,bz2,hulpx
      real*8  t350,t650,rlogp500,alpha
      real*8  detqmax,detqmaxexact
      real*4  hulp(0:iqmtab,0:jqmtab,0:kqmtab)


      rlogp500=log(500d0/650d0)
      b=cc2*cc3-cc2*tzero


!      call system('rm outputdata/atmos/qmtabel.dat')


!      open(88,file='qmtabel.dat')
!      open(89,file='qmtabel.test')
!     *  form='unformatted',access='direct',recl=51*21*21)


      do i=0,iqmtab
        tqmi(i)=tqmimin + i*dtqmi
      enddo
      do j=0,jqmtab
        tqmj(j)=tqmjmin + j*dtqmj
      enddo
      do k=0,kqmtab
        tqmk(k)=tqmkmin + k*dtqmk
      enddo


      hulpx=cc1*exp(cc2)/(rowat*grav)


      do i=0,iqmtab
        t650=tqmi(i)


        do j=0,jqmtab
          tmount=tqmj(j)+t650


          do k=0,kqmtab
            t350=t650-tqmk(k)


            alpha=(t350-t650)*rlogtl12


            t500=t650+alpha*rlogp500


            z1=1/(tmount-cc3)
            z2=1/(t500-cc3)


            bz1=b*z1
            bz2=b*z2


            qmax=(exp(bz1)+expint(1,-bz1)*bz1)/z1 - &
                 & (exp(bz2)+expint(1,-bz2)*bz2)/z2


            qmax=qmax*hulpx/alpha


            if (qmax.lt.0d0) qmax=0d0
            qmtabel(i,j,k)=qmax
!            write(88,111) tqmi(i),tqmj(j),tqmk(k),qmtabel(i,j,k)
          enddo
        enddo
      enddo


!      write(89) (((qmtabel(i,j,k),i=0,iqmtab),j=0,jqmtab),k=0,kqmtab)



!      do i=0,iqmtab
!        temp4g(1,1)=tqmi(i)
!        do j=0,jqmtab
!          tmount=tqmj(j)+temp4g(1,1)
!          do k=0,kqmtab
!            temp2g(1,1)=temp4g(1,1)-tqmk(k)
!            hulp(i,j,k)=detqmaxexact(tmount,1,1)
!            write(89,111) tqmi(i),tqmj(j),tqmk(k),hulp(i,j,k)
!          enddo
!        enddo
!      enddo



 111  format(4F14.7)
!      close(88)
!      close(89)
!      stop


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function detqmaxexact(tmount,i,j)
!-----------------------------------------------------------------------
! *** determines the maximum water content in latlon point
! *** i,j for given ground- and 650 and 350 hPa temperature
! *** by linear interpolation in qmtabel
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comunit.h'


      integer i,j
      real*8  tmount,alpha,t500,z1,z2,bz1,bz2,b,hulpx
      real*8  qmax,detqmaxexact,expint


      b=cc2*cc3-cc2*tzero


      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12


      t500=temp4g(i,j)+alpha*log(500d0/650d0)


      z1=1/(tmount-cc3)
      z2=1/(t500-cc3)


      bz1=b*z1
      bz2=b*z2


      hulpx=cc1*exp(cc2)/(rowat*grav*alpha)


      qmax=hulpx*(exp(bz1)+expint(1,-bz1)*bz1)/z1 - &
           & hulpx*(exp(bz2)+expint(1,-bz2)*bz2)/z2


      if (qmax.lt.0d0) qmax=0d0


      if (qmax.gt.0.2) then
        write(iuo+29,*) 'in latlon ',i,j,' qmax ',qmax
        call error(121)
      endif


      detqmaxexact=qmax


      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function detqmax(tmount,i,j,dqmdt)
!-----------------------------------------------------------------------
! *** determines the maximum water content in latlon point
! *** i,j for given ground- and 650 and 350 hPa temperature
! *** by linear interpolation in qmtabel
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comunit.h'


      integer i,j,ii,jj,kk
      real*8  tmount,ti,tj,tk
      real*8  qmax,detqmax,dqmdi,dqmdj,dqmdk
      real*8  dqmdt,hmount,z500,t500,alpha,dtgdt

      ti=temp4g(i,j)
      tj=tmount-temp4g(i,j)
      tk=temp4g(i,j)-temp2g(i,j)


      if (ti.lt.tqmi(0)) then
        ti=tqmi(0)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif
      if (ti.gt.tqmi(iqmtab)) then
        ti=tqmi(iqmtab)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif

      if (tj.lt.tqmj(0)) then
        tj=tqmj(0)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif
      if (tj.gt.tqmj(jqmtab)) then
        tj=tqmj(jqmtab)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif

      if (tk.lt.tqmk(0)) then
        tk=tqmk(0)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif
      if (tk.gt.tqmk(kqmtab)) then
        tk=tqmk(kqmtab)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif

      ii=min(iqmtab-1,int((ti-tqmimin)*rdtqmi))
      jj=min(jqmtab-1,int((tj-tqmjmin)*rdtqmj))
      kk=min(kqmtab-1,int((tk-tqmkmin)*rdtqmk))
!       if( ii.lt.0) then
!         PRINT *,ii, jj, kk
!         PRINT *,"min(",iqmtab,"-1,int((",ti,"-",tqmimin,")*",rdtqmi,"))"
!         write(*,*) temp4g
!       endif

      dqmdi=(qmtabel(ii+1,jj,kk)-qmtabel(ii,jj,kk))*rdtqmi
      dqmdj=(qmtabel(ii,jj+1,kk)-qmtabel(ii,jj,kk))*rdtqmj
      dqmdk=(qmtabel(ii,jj,kk+1)-qmtabel(ii,jj,kk))*rdtqmk

      qmax = qmtabel(ii,jj,kk) + (ti-tqmi(ii))*dqmdi + &
           & (tj-tqmj(jj))*dqmdj + (tk-tqmk(kk))*dqmdk
      if (qmax.lt.0d0) qmax=0d0

      if (qmax.gt.0.2) then
        write(iuo+29,*) 'in latlon ',i,j,' qmax ',qmax
        call error(121)
      endif


      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12
      t500=temp4g(i,j)+alpha*alogpl2tl2
      z500=gpm500*grav
      hmount=qmount(i,j)*hmoisr*grav


      dtgdt=(rgas*t500*alogtl1pl2 + (hmount-geopg(i,j,2)-z500))/ &
           & (rgas*tmount*alogtl12)


      dqmdt=dqmdi + dqmdj * (dtgdt - 1d0) + dqmdk


      detqmax=0.9*qmax
!     detqmax=0.85*qmax


      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function expint(n,x)
      implicit none
      integer n,maxit
      real*8 expint,x,eps,fpmin,euler
      parameter (maxit=100,eps=1.e-10,fpmin=1.e-30,euler=.5772156649)
      integer i,ii,nm1
      real*8 a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
!        pause 'bad arguments in expint'
        call error(20)
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0.)then
        expint=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./fpmin
        d=1./b
        h=d
        do 11 i=1,maxit
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.eps)then
            expint=h*exp(-x)
            return
          endif
11      continue
!        pause 'continued fraction failed in expint'
        call error(20)
      else
        if(nm1.ne.0)then
          expint=1./nm1
        else
          expint=-log(x)-euler
        endif
        fact=1.
        do 13 i=1,maxit
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-euler
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*eps) return
13      continue
!        pause 'series failed in expint'
        call error(20)
      endif
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisture
!-----------------------------------------------------------------------
! *** acvection and sources and sinks of atmospheric moisture
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  hdmoisg(nlat,nlon),d1(nlat,nlon),d3(nlat,nlon)
      real*8  d2(nlat,nlon),qstar,hdivmg(nlat,nlon),globalmean
      real*8  levtempgp,factemv,factems,omegg500,t500,qsat,gm1,gm2

      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair


      call rggtosp(rmoisg,rmoiss)
      call sptogg(rmoiss,rmoisg,pp)

      do j=1,nlon
        do i=1,nlat
          cormois(i,j)=0d0
          if (rmoisg(i,j).lt.0d0) then
            cormois(i,j)=cormois(i,j)-rmoisg(i,j)
            rmoisg(i,j)= 0d0
          endif
        enddo
      enddo


! *** horizontal divergence of moisture


      call trafluxdiv(hdivmg,rmoiss,rmoisg)

! *** vertical advection of moisture


      do j=1,nlon
        do i=1,nlat
          omegg500=(omegg(i,j,1)+omegg(i,j,2))/2.d0
!         d1(i,j)=omegg500
          vemoisg(i,j)=0d0
!         d2(i,j)=0d0
!          t500=levtempgp(plevel(2),i,j)
!          qstar=relhum(i,j)*qsat(plevel(2),t500)
!           d3(i,j)=qstar
!           d2(i,j)=t500
          if (omegg500.lt.0.d0) then
            t500=levtempgp(plevel(2),i,j)
!           d2(i,j)=t500
            qstar=relhum(i,j)*qsat(plevel(2),t500)
!           d3(i,j)=qstar
            vemoisg(i,j)=-omegg500*qstar/(grav*rowat)
            vemoisg(i,j)=min(vemoisg(i,j),rmoisg(i,j)/dtime)
            if (vemoisg(i,j).LT.0.d0) then
              vemoisg(i,j) = rmoisg(i,j)/dtime
            endif
            if (tsurf(i,j).ge.tzero) then
              dyrain(i,j) = dyrain(i,j) + vemoisg(i,j)
              thforg1(i,j)=thforg1(i,j) + factemv*vemoisg(i,j)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factemv*vemoisg(i,j)/dp1
            else
              dysnow(i,j) = dysnow(i,j) + vemoisg(i,j)
              thforg1(i,j)=thforg1(i,j) + factems*vemoisg(i,j)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factems*vemoisg(i,j)/dp1
            endif
          endif
        enddo
      enddo


! *** horizontal diffusion of moisture


      call hdiff(hdmoisg)


!
! *** time stepping forward time stepping
!


      do j=1,nlon
        do i=1,nlat
          rmoisg(i,j)=rmoisg(i,j)+dtime*(-ihavm*hdivmg(i,j) &
               & -ivavm*vemoisg(i,j) + hdmoisg(i,j) + imsink*evap(i,j))
          if (rmoisg(i,j).lt.0d0) then
            cormois(i,j)=cormois(i,j)-rmoisg(i,j)
            rmoisg(i,j)= 0d0
          endif
        enddo
      enddo

!      if (iyear.eq.7.and.imonth.eq.1) then
!        write(350) ((real(rmoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(hdivmg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(vemoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(hdmoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(cormois(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(d1(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(d2(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(d3(i,j)),j=1,nlon),i=1,nlat)
!      endif
      call moisbalance()
!      if (iyear.eq.7.and.imonth.eq.1) then
!        write(350) ((real(rmoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(cormois(i,j)),j=1,nlon),i=1,nlat)
!      endif


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine trafluxdiv(tfdiv,ctrasp,ctra)
!-----------------------------------------------------------------------
! *** computes horizontal divergence of tracer flux
!-----------------------------------------------------------------------

      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k
      real*8  ctrasp(nsh2),vv(nsh2),ww(nsh2)
      real*8  dcdl(nlat,nlon),dcdm(nlat,nlon)
      real*8  tfdiv(nlat,nlon),ctra(nlat,nlon)


! *** 800 hPa winds are reduced with umoisr in the advection of the
! *** tracer field

! *** spatial derivatives of tracer


      call ddl (ctrasp,vv)
      call sptogg (vv,dcdl,pp)
      call sptogg (ctrasp,dcdm,pd)


! *** advection of tracer by total wind + convergence of tracer


      do j=1,nlon
        do i=1,nlat
          tfdiv(i,j)=dcdl(i,j)*(u800(i,j) + udivg(i,j,3))/ &
               & (radius*cosfi(i)) + dcdm(i,j)*(v800(i,j) + vdivg(i,j,3))/ &
               & (radius/cosfi(i)) + ctra(i,j)*divg(i,j,3)
          tfdiv(i,j)=tfdiv(i,j)*umoisr


        enddo
      enddo


      call rggtosp(tfdiv,vv)
      vv(1)=0d0
      call sptogg (vv,tfdiv,pp)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine hdivspec(hduvg,ug,vg)
!-----------------------------------------------------------------------
! *** computes horizontal divergence
!-----------------------------------------------------------------------

      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  hduvg(nlat,nlon),ug(nlat,nlon),vg(nlat,nlon)
      real*8  dugdl(nlat,nlon),dvgdm(nlat,nlon)
      real*8  usp(nsh2),vv(nsh2),vsp(nsh2)
      real*8  dx,dy(nlat)


      do j=1,nlon
        do i=1,nlat
          vg(i,j)=vg(i,j)*cosfi(i)
        enddo
      enddo


      call rggtosp(ug,usp)
      call rggtosp(vg,vsp)
      call ddl (usp,vv)
      call sptogg (vv,dugdl,pp)
      call sptogg (vsp,dvgdm,pd)


      do j=1,nlon
        do i=1,nlat
          hduvg(i,j)= (dugdl(i,j)/cosfi(i)+dvgdm(i,j))/radius
        enddo
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine hdiff(hdmg)
!-----------------------------------------------------------------------
! *** horizontal diffusion of moisture
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer idifq,k
      real*8  hdmoiss(nsh2),hdmg(nlat,nlon)
      real*8  difq,rll

      difq=max(0.d0,1.d0/(tdifq*24d0*3600d0))


      call lap(rmoiss,hdmoiss)


      hdmoiss(1)=0d0


      do k=2,nsh2
        hdmoiss(k)=difq*hdmoiss(k)
      enddo


      call sptogg(hdmoiss,hdmg,pp)


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine convec
!-----------------------------------------------------------------------
! *** moist convective adjustment
!-----------------------------------------------------------------------
      implicit none



      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comrunlabel.h'


      integer ncmax,iconvn,i,j
      real*8  qsatcr,tsatcr,pref,t500,qsat500,pot2g,pot4g,dcmoisg
      real*8  fachulp,facteta,factemv,factems,pmount,tmount
      real*8  qmax,qsat,hulp,redrain
      real*8  temp2go,temp4go,levtempgp
      real*8  detqmax,drainm,crainm,dqmdt


      fachulp=0.622d0*(rlatvap**2)/(cpair*rgas)
      facteta=0.6d0*rgas*(2.d0**rkappa)/grav
      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair
      ncmax=0
      crainm=0.5d0*rainmax
      drainm=rainmax-crainm



      do j=1,nlon
        do i=1,nlat
          iconvn=0
  10      continue


! ***     calculate pressure and temperature at the ground
! ***     and the maximum water content


            call ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)


! ***     relhmax defines the relative humidity at which oversaturation
! ***     occurs


            qmax=relhmax*qmax


            pot2g=temp2g(i,j)/potfac1
            pot4g=temp4g(i,j)/potfac2
            teta(i,j)=0.5d0*(pot2g-pot4g)
! ***       dry adiabatic lapse rate
            tetacr(i,j)=0d0
            dcmoisg=0d0


            if (rmoisg(i,j).gt.qmax) then


! ***     calculate rain reduction factor to account for increased
! ***     moisture capacity due to latent heat release


              redrain=1d0+dqmdt*relhmax*rowat*rlatvap*grav/(cpair*dp2)


              dcmoisg=(rmoisg(i,j)-qmax)/(redrain*dtime)


! ***         if the air is supersaturated initially, this is due to
! ***         large scale convergence of moisture and large scale
! ***         temperature changes. Excessive moisture
! ***         is then removed as dynamic rain

              if (iconvn.eq.0) then
                if (tsurf(i,j).ge.tzero) then
                  dyrain(i,j)=dyrain(i,j) + dcmoisg
                  if (dyrain(i,j).gt.drainm) then
                    dcmoisg=drainm-dyrain(i,j)+dcmoisg
                    dyrain(i,j)=drainm
                  endif
                  thforg2(i,j)=thforg2(i,j) + factemv*dcmoisg/dp2
!                   write(*,*) "thforg2", thforg2(i,j),"=", factemv,"*",dcmoisg,"/",dp2
                  vhforg2(i,j)=vhforg2(i,j) + factemv*dcmoisg/dp2
                else
                  dysnow(i,j)=dysnow(i,j) + dcmoisg
                  if (dysnow(i,j).gt.drainm) then
                    dcmoisg=drainm-dysnow(i,j)+dcmoisg
                    dysnow(i,j)=drainm
                  endif
                  thforg2(i,j)=thforg2(i,j) + factems*dcmoisg/dp2
                  vhforg2(i,j)=vhforg2(i,j) + factems*dcmoisg/dp2
                endif
              else
                if (tsurf(i,j).ge.tzero) then
                  corain(i,j)=corain(i,j)+ dcmoisg
                  if (corain(i,j).gt.crainm) then
                    dcmoisg=crainm-corain(i,j)+dcmoisg
                    corain(i,j)=crainm
                  endif
                  temp4g(i,j) =temp4g(i,j) + factemv*dcmoisg*dtime/dp2
                  vhforg2(i,j)=vhforg2(i,j)+ factemv*dcmoisg/dp2
                else
                  cosnow(i,j)=cosnow(i,j)+ dcmoisg
                  if (cosnow(i,j).gt.crainm) then
                    dcmoisg=crainm-cosnow(i,j)+dcmoisg
                    cosnow(i,j)=crainm
                  endif
                  temp4g(i,j) =temp4g(i,j) + factems*dcmoisg*dtime/dp2
                  vhforg2(i,j)=vhforg2(i,j)+ factems*dcmoisg/dp2
                endif
              endif


! ***         moisture changes due to dynamic and convective rain


              rmoisg(i,j)=rmoisg(i,j)-ivavm*dcmoisg*dtime

              if (rmoisg(i,j).lt.0d0) then
!                write(99,*) 'moisture less than zero in convec !!!'
!                write(99,*) i,j,rmoisg(i,j),qmax,redrain
                rmoisg(i,j)= 0d0
              endif



! ***         calculate moist adiabatic lapse rate


              pot2g=temp2g(i,j)/potfac1
              pot4g=temp4g(i,j)/potfac2
              teta(i,j)=0.5d0*(pot2g-pot4g)


!              t500 = levtempgp(plevel(2),i,j)


              t500 = 0.5d0 * (temp2g(i,j) + temp4g(i,j))


              qsat500=qsat(plevel(2),t500)


              hulp=1d0 + fachulp*qsat500/(t500**2)
              gams(i,j)=gamd*(1+(rlatvap*qsat500)/(rgas*t500))/hulp
              tetacr(i,j)=0.5d0*t500*facteta*(gamd-gams(i,j))

            endif
            if (teta(i,j) .lt. (tetacr(i,j)-0.1)) then


              pot2g=(dp1*temp2g(i,j)+dp2*temp4g(i,j)+ &
     &                                 dp2*potfac2*2.*tetacr(i,j))/ &
     &              (potfac1*dp1+potfac2*dp2)
              pot4g=pot2g - 2.d0*tetacr(i,j)
              temp2go=temp2g(i,j)
              temp4go=temp4g(i,j)
              temp2g(i,j)=pot2g*potfac1
              temp4g(i,j)=pot4g*potfac2
!               write(*,*) "ici", temp4g(1,1),"=", pot4g,"*",potfac2
              vhforg1(i,j)=vhforg1(i,j) + (temp2g(i,j)-temp2go)/dtime
              vhforg2(i,j)=vhforg2(i,j) + (temp4g(i,j)-temp4go)/dtime
              iconvn=iconvn+1
              if (iconvn.lt.10) then
                if (dyrain(i,j).eq.drainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' dyrain'
                  call error(122)
                  goto 20
                endif
                if (corain(i,j).eq.crainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' corain'
                  call error(122)
                  goto 20
                endif
                if (dysnow(i,j).eq.drainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' dysnow'
                  call error(122)
                  goto 20
                endif
                if (cosnow(i,j).eq.crainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' cosnow'
                  call error(122)
                  goto 20
                endif
                goto 10
              else
                write(iuo+29,*) 'warning in lat-lon point: ',i,j
                write(iuo+29,*) temp2g(i,j),temp4g(i,j)
                call error(120)
                goto 20
              endif
            else
              goto 20
            endif
  20      continue
          if (iconvn.gt.ncmax) ncmax=iconvn
          torain(i,j)= dyrain(i,j) + corain(i,j)
          tosnow(i,j)= dysnow(i,j) + cosnow(i,j)
        enddo
      enddo
!      write(99,*) 'maximum iterations of convection : ',ncmax


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisbalance
!-----------------------------------------------------------------------
! *** fix moisture balance in the atmosphere: due to advection and
! *** spectral truncation, the moisture at a specific point can become
! *** negative. All moisture additions to reset the moisture at zero
! *** have been accumulated in cormois. In this routine these additions
! *** are subtracted from rmoisg. This is done for each latitude
! *** separately to prevent an artificial meridional transport of
! *** moisture. If at a given latitude not enough moisture is available
! *** to accomodate the subtraction, the moisture balance is fixed at
! *** the neighbouring equatorward latitude. The weights pw(nlat,1) are
! *** used to correct for changes in the gridsize with latitude.
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comrunlabel.h'


      integer i,j,nn
      real*8  gmc,gmm,gfac


      do i=1,nlat/2
        gmc=0d0
        gmm=0d0
        do j=1,nlon
          gmc=gmc+cormois(i,j)
          gmm=gmm+rmoisg(i,j)
        enddo
        if (gmm.gt.0d0.and.gmc.gt.0d0) then
          if (gmm.gt.gmc) then
            gfac=gmc/gmm
            do j=1,nlon
              rmoisg(i,j)=rmoisg(i,j)-gfac*rmoisg(i,j)
            enddo
          else
            gfac=gmm/gmc
            gmc=(gmc-gmm)*pw(i,1)/(pw(i+1,1)*dble(nlon))
            do j=1,nlon
              rmoisg(i,j)=0d0
              cormois(i,j)=gfac*cormois(i,j)
              cormois(i+1,j)=cormois(i+1,j)+gmc
            enddo
          endif
        endif
      enddo

      do i=nlat,1+nlat/2,-1
        gmc=0d0
        gmm=0d0
        do j=1,nlon
          gmc=gmc+cormois(i,j)
          gmm=gmm+rmoisg(i,j)
        enddo
        if (gmm.gt.0d0.and.gmc.gt.0d0) then
          if (gmm.gt.gmc) then
            gfac=1d0-gmc/gmm
            do j=1,nlon
              rmoisg(i,j)=rmoisg(i,j)*gfac
            enddo
          else
            gfac=gmm/gmc
            gmc=(gmc-gmm)*pw(i,1)/(pw(i-1,1)*dble(nlon))
            do j=1,nlon
              rmoisg(i,j)=0d0
              cormois(i,j)=gfac*cormois(i,j)
              cormois(i-1,j)=cormois(i-1,j)+gmc
            enddo
          endif
        endif
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function qsat(press,temp)
!-----------------------------------------------------------------------
! *** saturation mixing ratio
! *** input press in [Pa], temp in K
! *** output qsat: saturation mixing ratio
!-----------------------------------------------------------------------
      implicit none
      include 'comatm.h'
      include 'comphys.h'


      real*8  press,temp,qsat


      qsat=cc1*exp(cc2*(temp-tzero)/(temp-cc3))/press

      end



!23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine levtemp(tlev,plev)
!-----------------------------------------------------------------------
! *** computation temperatures at level p [Pa] assuming a constant
! *** temperature lapse rate : dt/dlnp = constant
! *** input: plev
! *** ouput: tlev
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  tlev(nlat,nlon)
      real*8  r,plev


      r=log(plev/65000.d0)*rlogtl12
      do j=1,nlon
        do i=1,nlat
          tlev(i,j)=temp4g(i,j)+r*(temp2g(i,j)-temp4g(i,j))
        enddo
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function levtempgp(plev,i,j)
!-----------------------------------------------------------------------
! *** computation temperatures at level p [Pa] assuming a constant
! *** temperature lapse rate : dt/dlnp = constant
! *** input: plev [Pa],i,j
! *** output: levtempgp [K]
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'


      integer i,j
      real*8  levtempgp
      real*8  r,plev


      r=log(plev/65000.d0)*rlogtl12
      levtempgp=temp4g(i,j)+r*(temp2g(i,j)-temp4g(i,j))

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fortemp
!-----------------------------------------------------------------------
! *** computes  new atmospheric temperatures due to heatforcing
! *** apply diffusion in stratosphere: timestep smaller than:
! *** 2*diffusiontimescale/462(=21*22) fe 10 days, 1 hour timestep
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comrunlabel.h'


      integer i,j,k,it,ipd
      real*8  temp0sp(nsh2),tdifc,globalmean,tstep,tdifday

      tdifday=100d0
      tdifc=1.0d0/(tdifday*24.*3600.)

      tstep=2d0*tdifday*24.*3600./462.


      ipd=1+int(dtime)/int(tstep)


      tstep=dtime/ipd

      call ggtosp(temp0g,temp0sp)

      do it=1,ipd

        do k=1,nsh2
          temp0sp(k)=temp0sp(k) + tstep*rinhel(k,0)*temp0sp(k)*tdifc
        enddo
      enddo

      call sptogg(temp0sp,temp0g,pp)


      do j=1,nlon
        do i=1,nlat
          temp0g(i,j)=temp0g(i,j) + dtime*thforg0(i,j)
          temp2g(i,j)=temp2g(i,j) + dtime*thforg1(i,j)
          temp4g(i,j)=temp4g(i,j) + dtime*thforg2(i,j)
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine tempprofile
!-----------------------------------------------------------------------
! *** computation of vertical temperature profile
! *** based on reference profiles from NCEP reanalysis data
! *** also used in the LWR parameterisation
! *** assuming temperature anomalies wrt these temperature profiles
! *** vary linearly with log of the pressure
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comrunlabel.h'

      integer i,j,k,l,ireg(2),is,ism,nn,k1,k2,k2_tmp
      real*8 ro,ro1,ro2,z,z0,dt350,dt650(2),beta(2),tsref,dt100,z1,z2
      real*8 dtemp_tmp, pground_tmp,tsref_tmp,z_tmp,ro_tmp,z0_tmp
      real*8 beta_tmp

! *** Example reference profile for
! *** zonal band between 15s and 15n
! *** Month  4 : tncep(19,27,12)
! ***  nr.   pfl     tfl
! ***        (mb)    (K)
! ***  1   10.00  235.113
! ***  2   20.00  223.086
! ***  3   30.00  216.125
! ***  4   50.00  206.119
! ***  5   70.00  198.420
! ***  6  100.00  196.311
! ***  7  150.00  208.144
! ***  8  200.00  221.461
! ***  9  250.00  232.611
! *** 10  300.00  242.385
! *** 11  400.00  257.836
! *** 12  500.00  268.271
! *** 13  600.00  276.160
! *** 14  700.00  282.859
! *** 15  850.00  290.298
! *** 16  925.00  294.405
! *** 17 1000.00  299.345
! *** 18 1011.99  300.156
! *** Ps 1013.00  301.265 Ts


      is=imonth/3+1
      if (is.gt.4) is=1
      ism=(is-1)*3+1
      do j=1,nlon
        do i=1,nlat
          ireg(1)=irn(i,j,1)
          ireg(2)=irn(i,j,2)

! *** logarithmic interpolation of temperature anomalies; 200 hPa or
! *** higher the anomalies approach T100 temperature within 3
! *** pressure levels


          do nn=1,2
            dt100=temp0g(i,j)-tncep(6,ireg(nn),imonth)
            dt350=temp2g(i,j)-tncep12(1,ireg(nn),imonth)
            dt650(nn)=temp4g(i,j)-tncep12(2,ireg(nn),imonth)
            beta(nn)=(dt350-dt650(nn))*rlogtl12
            do k=1,6
              dtemp(k,i,j,nn)=dt100
            enddo
            do k=7,9
              dtemp(k,i,j,nn)=((10-k)*dt100+(k-6)*(dt650(nn) &
                   & + beta(nn)*rlogtl(k)))*0.25
            enddo
            do k=10,17
              dtemp(k,i,j,nn)=dt650(nn) + beta(nn)*rlogtl(k)
            enddo
          enddo
! *** from a mean height of 500 hPa from NCEP reanalysis data, the
! *** surface pressure is found using hydrostatic equilibrium
! *** and ideal gas law.


          z=z500ncep(ireg(1),imonth)
          k1=12
          DO WHILE (z.GT.rmountn(i,j,noc).AND.k1.LT.17)
            z0=z
            ro1=pncep(k1)/(rgas*(dtemp(k1,i,j,1)+tncep(k1,ireg(1),imonth)))
            ro2=pncep(k1+1)/(rgas*(dtemp(k1+1,i,j,1)+ &
                 & tncep(k1+1,ireg(1),imonth)))
            ro=(ro1+ro2)*0.5
            z=z-(pncep(k1+1)-pncep(k1))/(grav*ro)
            k1=k1+1
          END DO
          z1=z
          pgroundn(i,j,noc)=ro*grav*(z0-rmountn(i,j,noc))+pncep(k1-1)
          pgroundn(i,j,nse)=pgroundn(i,j,noc)


          z=z500ncep(ireg(2),imonth)
          k2=12

          DO WHILE ((z.GT.rmountn(i,j,nld)).AND.(k2.LT.17))
            z0=z
            ro1=pncep(k2)/(rgas*(dtemp(k2,i,j,2)+tncep(k2,ireg(2),imonth)))
            ro2=pncep(k2+1)/(rgas*(dtemp(k2+1,i,j,2)+ &
                 & tncep(k2+1,ireg(2),imonth)))
            ro=(ro1+ro2)*0.5
            z=z-(pncep(k2+1)-pncep(k2))/(grav*ro)
            k2=k2+1
          END DO

          z2=z
          pgroundn(i,j,nld)=ro*grav*(z0-rmountn(i,j,nld))+pncep(k2-1)
          pground(i,j)=0.0
          do nn=1,ntyps
            pground(i,j)=pground(i,j)+fractn(i,j,nn)*pgroundn(i,j,nn)
          enddo



! *** Temperature of air near surface is then found by
! *** interpolating for this pressure level the temperature profile calculated
! *** above.

          dtemp(18,i,j,1)=dt650(1)+beta(1)*log(pgroundn(i,j,noc)/65000.)
          beta(1)=(tncep(k1,ireg(1),imonth)-tncep((k1-1),ireg(1),imonth))/ &
               & log(pncep(k1)/pncep(k1-1))
          tsref=tncep((k1-1),ireg(1),imonth)+ &
               & beta(1)*log(pgroundn(i,j,noc)/pncep(k1-1))
          tempsgn(i,j,noc)=tsref+dtemp(18,i,j,1)
          tempsgn(i,j,nse)=tempsgn(i,j,noc)
          dtemp(18,i,j,2)=dt650(2)+beta(2)*log(pgroundn(i,j,nld)/65000.)
          beta(2)=(tncep(k2,ireg(2),imonth)-tncep((k2-1),ireg(2),imonth))/ &
               & log(pncep(k2)/pncep(k2-1))
          tsref=tncep((k2-1),ireg(2),imonth)+ &
               & beta(2)*log(pgroundn(i,j,nld)/pncep(k2-1))
          tempsgn(i,j,nld)=tsref+dtemp(18,i,j,2)
          tempsg(i,j)=0.0
          do nn=1,ntyps
            tempsg(i,j)=tempsg(i,j)+fractn(i,j,nn)*tempsgn(i,j,nn)
          enddo

! *** Temperature of pressure levels between diagnosed surface pressure and
! *** reference surface pressure are set equal to surface air temp.


          IF(z1.LE.rmountn(i,j,noc))THEN
            DO l=ipl(ireg(1)),k1,-1
               dtemp(l,i,j,1)=tempsgn(i,j,noc)-tncep(l,ireg(1),imonth)
            ENDDO
          ENDIF

          IF(z2.LE.rmountn(i,j,nld))THEN
            DO l=ipl(ireg(2)),k2,-1
               dtemp(l,i,j,2)=tempsgn(i,j,nld)-tncep(l,ireg(2),imonth)
            ENDDO
          ENDIF

! *** for LWR parameterisation temperature anomalies wrt seasonal mean are
! *** required:
          do nn=1,2
            DO k=1,18
              dtemp(k,i,j,nn)=tncep(k,ireg(nn),imonth) + &
                   & dtemp(k,i,j,nn)-tncep(k,ireg(nn),ism)
            ENDDO
          enddo


         enddo
       enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)
!-----------------------------------------------------------------------
! *** computation of ground pressure and temperature in order to
! *** to calculate the maximum precipitable water content in latlon i,j
! *** qmount contains the topography for this purpose
! *** assuming temperature varies linearly with log of the pressure
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comunit.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  t500,levtempgp,hmount,hred,z500,dqmdt
      real*8  alpha,pfac,hfac,pmount,tmount,qmax,detqmax


      z500=gpm500*grav
      hfac=2/rgas
      hred=hmoisr*grav
      pfac=log(plevel(2)/tlevel(2))


! *** calculate temperature at t500 assuming the temperature varies
! *** linearly with log(p) : T = Tref + alpha * log (p/pref)


      alpha=(temp2g(i,j) - temp4g(i,j))*rlogtl12
      t500 =temp4g(i,j) + alpha*pfac


! *** calculate reduced ground height in decameters
! *** reduction occurs in order to tune the amount of moisture which
! *** is allowed to pass a topographic barier


      hmount=qmount(i,j)*hred
      if (hmount.lt.0d0) hmount=0d0


! *** calculate the groundpressure assuming that the mean geopotential
! *** height at 500 hPa is gpm500 decameter
! *** calculate 10 mtr temperature in K


      tmount=t500**2 - hfac*alpha*(hmount-geopg(i,j,2)-z500)
      if (tmount.lt.0) then
        write(iuo+29,*) 'in latlon ',i,j
        write(iuo+29,*) tmount,hmount,t500,geopg(i,j,2)
        call error(18)
      else
        tmount=sqrt(tmount)
      endif


!      pmount=plevel(2)*exp((tmount-t500)/alpha)

      qmax=detqmax(tmount,i,j,dqmdt)


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function globalmean(gfield)
!-----------------------------------------------------------------------
! *** computes global mean of gfield
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'


      integer i,j
      real*8  gfield(nlat,nlon),sum(nlat),globsum,globfac,globalmean


      globfac=1d0/dsqrt(dble(nlon))


      do i=1,nlat
        sum(i)=0d0
      enddo


      do j=1,nlon
        do i=1,nlat
          sum(i)=sum(i)+gfield(i,j)
        enddo
      enddo


      globsum=0d0


      do i=1,nlat
        globsum=globsum+sum(i)*pw(i,1)
      enddo


      globalmean=globsum*globfac

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine meantemp
!-----------------------------------------------------------------------
! *** computes mean atmospheric temperatures
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      real*8  globalmean


! *** mean temperatures

      tempm(0)=globalmean(temp0g)
      tempm(1)=globalmean(temp2g)
      tempm(2)=globalmean(temp4g)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine dyntemp
!-----------------------------------------------------------------------
! *** computes temperature distribution in K from geopotential
! *** the mean level is given by tempm
! *** input:  geopg,tempm
! *** output: temp2g,temp4g
!-----------------------------------------------------------------------
      implicit none



      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,l
      real*8  tempfac(ntl)
      real*8  geogt(nlat,nlon,ntl),tempgt(nlat,nlon,ntl)


      tempfac(1)=350.d0/(rgas*300.d0)
      tempfac(2)=650.d0/(rgas*300.d0)


      do j=1,nlon
        do i=1,nlat
          geogt(i,j,1)=geopg(i,j,1)-geopg(i,j,2)
          geogt(i,j,2)=geopg(i,j,2)-geopg(i,j,3)
        enddo
      enddo

      do l=1,ntl


! ***  calculate temperatures and add mean temperature level


        do j=1,nlon
          do i=1,nlat
            tempgt(i,j,l)=tempfac(l)*geogt(i,j,l)+tempm(l)
          enddo
        enddo


      enddo

      do j=1,nlon
        do i=1,nlat
          temp2g(i,j)=tempgt(i,j,1)
          temp4g(i,j)=tempgt(i,j,2)
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine vortfor
!-----------------------------------------------------------------------
! *** computes vorticity forcing
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)
      real*8  forcgg(nlat,nlon),forcgg2(nlat,nlon)
      real*8  zetas(nsh2,3),zetag(nlat,nlon,3)
      real*8  dimfac


      call rggtosp(vhforg1,dfor1)
      call rggtosp(vhforg2,dfor2)


! ***  compute the relative vorticity zeta from steamfunction psi
! ***  psi is dimensionless, zeta with dimension


      do l=1,nvl
        do k=1,nsh2
          zetas(k,l)=rinhel(k,0)*psi(k,l)*om
        enddo
        call sptogg(zetas(1,l),zetag(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=0d0
          enddo
        enddo
      enddo


! *** potential vorticity (pv) forcing (1/(second**2)) due to
! *** the diabatic heating


      if (ipvf1.ne.0) call pvf1(vforg)


! *** pv forcing due to advection of the planetary vorticy f
! *** by the divergent wind

      if (ipvf2.ne.0) call pvf2(vforg)


! *** pv forcing due to d*zeta. d is the divergence.
! *** zeta is the relative vorticity


      if (ipvf3.ne.0) call pvf3(vforg,zetag)


! *** pv forcing due to advection of the relative vorticy
! *** zeta by the divergent wind


      if (ipvf4.ne.0) call pvf4(vforg,zetas)


! *** pv forcing due to advection of temperature by
! *** the divergent wind


      if (ipvf5.ne.0) call pvf5(vforg)


! *** the total pv forcing in nondimensional units


      dimfac=1./(om**2)


      call ggtosp(vforg(1,1,1),vfor1)
      call ggtosp(vforg(1,1,2),vfor2)
      call ggtosp(vforg(1,1,3),vfor3)


! *** adding the artificial forcing (forcgg1)

      vfor1(1)=0.d0
      vfor2(1)=0.d0
      vfor3(1)=0.d0


      do k=2,nsh2
        for(k,1)=dimfac*vfor1(k)
        for(k,2)=dimfac*vfor2(k)
        for(k,3)=dimfac*vfor3(k)
      enddo


! *** transfer to grid point


      call sptogg(vfor1,vforg1,pp)
      call sptogg(vfor2,vforg2,pp)
      call sptogg(vfor3,vforg3,pp)


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf1(vforg)
!-----------------------------------------------------------------------
! *** potential vorticity (pv) forcing (1/(second**2)) due to
! *** the adiabatic heating
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  pvf1s(nsh2,3),vforg(nlat,nlon,nvl)
      real*8  vhfor1(nsh2),vhfor2(nsh2)
      real*8  vhforg1x(nlat,nlon),vhforg2x(nlat,nlon)
      real*8  vorfac1,vorfac2,drdef1,drdef2


      vorfac1=+rgas*dp/3.5d+4
      vorfac2=+rgas*dp/6.5d+4
      drdef1=1./( om*fzero*(radius*rrdef1)**2 )
      drdef2=1./( om*fzero*(radius*rrdef2)**2 )


      do j=1,nlon
        do i=1,nlat
          vhforg1x(i,j)=vhforg1(i,j)*sinfi(i)/fzero
          vhforg2x(i,j)=vhforg2(i,j)*sinfi(i)/fzero
        enddo
      enddo


      call rggtosp(vhforg1x,vhfor1)
      call rggtosp(vhforg2x,vhfor2)


      pvf1s(1,1)=0d0
      pvf1s(1,2)=0d0
      pvf1s(1,3)=0d0


      do k=2,nsh2
        pvf1s(k,1)=-drdef1*vorfac1*vhfor1(k)
        pvf1s(k,2)=drdef1*vorfac1*vhfor1(k)-drdef2*vorfac2*vhfor2(k)
        pvf1s(k,3)=drdef2*vorfac2*vhfor2(k)
      enddo


      do l=1,nvl
        call sptogg(pvf1s(1,l),vforg(1,1,l),pp)
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf2(vforg)
!-----------------------------------------------------------------------
! *** pv forcing due to advection of the planetary vorticy f
! *** by the divergent wind
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)

      do l=1,nvl
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)-vdivg(i,j,l)*om*cosfi(i)/radius
          enddo
        enddo
      enddo


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf3(vforg,zetag)
!-----------------------------------------------------------------------
! *** pv forcing due to d*zeta. d is the divergence.
! *** zeta is the relative vorticity
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),zetag(nlat,nlon,nvl)


      do l=1,nvl
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)-divg(i,j,l)*zetag(i,j,l)
          enddo
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf4(vforg,zetas)
!-----------------------------------------------------------------------
! *** pv forcing due to advection of the relative vorticy
! *** zeta by the divergent wind
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),zetas(nsh2,nvl)
      real*8  x(nsh2),xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)


      do l=1,nvl
        do k=1,nsh2
          x(k)=zetas(k,l)
        enddo
        call ddl(x,xhelp)
        call sptogg(xhelp,dxdl,pp)
        call sptogg(x,dxdm,pd)
        do i=1,nlat
          do j=1,nlon
             vforg(i,j,l)=vforg(i,j,l) &
     &                  -udivg(i,j,l)*dxdl(i,j)/(radius*cosfi(i)) &
     &                  -vdivg(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf5(vforg)
!-----------------------------------------------------------------------
! *** computes vorticity forcing due to advection of temperature by
! *** divergent wind
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),pvf7g(nlat,nlon,nvl)
      real*8  ud(nlat,nlon,2),vd(nlat,nlon,2)
      real*8  x(nsh2,2),y(nlat,nlon,3)
      real*8  xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)
      real*8  rr1,rr2,sinfact


      rr1=1.d0/(rrdef1*radius)**2
      rr2=1.d0/(rrdef2*radius)**2


      do j=1,nlon
        do i=1,nlat
          ud(i,j,1)=(udivg(i,j,1)+udivg(i,j,2))/2.d0
          ud(i,j,2)=(udivg(i,j,2)+udivg(i,j,3))/2.d0
          vd(i,j,1)=(vdivg(i,j,1)+vdivg(i,j,2))/2.d0
          vd(i,j,2)=(vdivg(i,j,2)+vdivg(i,j,3))/2.d0
        enddo
      enddo

      do k=1,nsh2
        x(k,1)=(psi(k,1)-psi(k,2))*om*radius**2
        x(k,2)=(psi(k,2)-psi(k,3))*om*radius**2
      enddo


      do l=1,2
        call ddl(x(1,l),xhelp)
        call sptogg(xhelp,dxdl,pp)
        call sptogg(x(1,l),dxdm,pd)

        do j=1,nlon
          do i=1,nlat
             y(i,j,l)=ud(i,j,l)*dxdl(i,j)/(radius*cosfi(i)) &
     &               +vd(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo


      do j=1,nlon
        do i=1,nlat
          pvf7g(i,j,1)=rr1*y(i,j,1)
          pvf7g(i,j,2)=-rr1*y(i,j,1)+rr2*y(i,j,2)
          pvf7g(i,j,3)=-rr2*y(i,j,2)
        enddo
      enddo


      do i=1,nlat
        sinfact=(sinfi(i)/fzero)**2
        do j=1,nlon
          vforg(i,j,1)=vforg(i,j,1)+pvf7g(i,j,1)*sinfact
          vforg(i,j,2)=vforg(i,j,2)+pvf7g(i,j,2)*sinfact
          vforg(i,j,3)=vforg(i,j,3)+pvf7g(i,j,3)*sinfact
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lwaverad2(nn)
!-----------------------------------------------------------------------
! *** computes long wave radiation according to the parameterization of
! *** Chao Chou and Neelin and substantially adapted and extended
! *** for global scale and more
! *** specific ECBILT application by the one and only Michiel Schaeffer
! ***
! *** parameters: nlat   = number of gridpoints in the meridional
! ***                      direction (32)
! ***             nlon   = number of gridpoints in the zonal
! ***                      direction (64)
! ***
! *** input : dtemp(19,nlat,nlon): temperature anomalies [K] wrt ncep
! ***                              climatology tncep in common lwrscheme
! ***         dqa(nlat,nlon) : anomalies of total prec. water cont. below
! ***                          500 hPa wrt ncep climatology
! ***         tcc(nlat,nlon)  : total cloud cover
! ***         ghg(19) : concentrations of well mixed ghg's (see comphys.h)
! ***
! *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             toa
! ***          ulrad2(nlat,nlon): net longwave radiation [Wm-2] at
! ***                             (500 hPa)
! ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
! ***                             the surface
! ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             the surface
!-----------------------------------------------------------------------
           implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'


      integer i,j,l,k,m,is,ism,nol,nn,ireg,h,r,s,igas
      real*8  lwrz(7,0:1),dumts
      real*8  dqa,dqreg(27)
      real*8  ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm
      real*8  ulrad0nmm,ulrad1nmm
      real*8  ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  ulrad1nzz(nlat,nlon,3)
      real*8  ulrad2nz(nlat,nlon),ulradsnz(nlat,nlon)
      real*8  dlradsnz(nlat,nlon)
      real*8  ulrad0nT,ulrad1nT,ulrad2nT,ulradsnT,dlradsnT
      real*8  globalmean
      real*8  logco2T,sqrch4T,sqrn2oT,ghgz(20)
      real*8  alpho3lw(2)
      real*4  lwrfluxz(7,27,4,0:1,2)
      real*8  moc,tmc,tmc0,tsurfmean,cland,thex
      common / radO3 / ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      common /rad031/ulrad0nz,ulrad1nz,ulrad0nT,ulrad1nT
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex



        ghgz(1)=280.
        do igas=2,20
         ghgz(igas)=ghg(igas)
        enddo
        logco2T=log(ghgz(1)/ghgipcc(1))
        sqrch4T=sqrt(ghgz(2))-sqrt(ghgipcc(2))
        sqrn2oT=sqrt(ghgz(3))-sqrt(ghgipcc(3))
        alpho3lw(1)=153.6
        alpho3lw(2)=201.2
        do h=1,2
        do l=0,1
         do s=1,4
          do r=1,27
           do k=1,7
            lwrfluxz(k,r,s,l,h)=lwrref(k,r,s,l)+lwrghg(k,1,r,s,l)*logco2T+ &
                 & lwrghg(k,2,r,s,l)*sqrch4T+lwrghg(k,3,r,s,l)*sqrn2oT
            do m=4,19
             lwrfluxz(k,r,s,l,h)=lwrfluxz(k,r,s,l,h)+ &
                  & lwrghg(k,m,r,s,l)*(ghgz(m)-ghgipcc(m))
            enddo
              lwrfluxz(k,r,s,l,h)=lwrfluxz(k,r,s,l,h)+ &
                   & lwrghg(k,4,r,s,l)*alpho3lw(h)*(ghgz(20)-25.)
           enddo
          enddo
         enddo
        enddo
        enddo

       is=imonth/3+1
       if (is.gt.4) is=1
       ism=(is-1)*3+1

      do i=1,27
!dqa    dqreg(i)=qancep(i,ism)**0.3333
        dqreg(i)=qancep(i,ism)
      enddo

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
        do i=1,nlat
          ireg=irn(i,j,nol)

!-Hemispheric dependence of tropospheric ozone forcing
          if (i.le.16) then
           h=1
          else
           h=2
          endif

!dqa      dqa=lwrmois(i,j)-dqreg(ireg)
!dqa      q**1/3-qm**1/3=qm**(1/3-n)*(q**n-qm**n)
          dqa=dqreg(ireg)**(0.3333-EXPIR)* &
               & (lwrmois(i,j)**EXPIR-dqreg(ireg)**EXPIR)

          do l=0,1
            do k=1,7
              lwrz(k,l)=lwrfluxz(k,ireg,is,l,h)+lwrqa(k,ireg,is,l)*dqa
              do m=1,ipl(ireg)-1
                lwrz(k,l)=lwrz(k,l)+lwrt(k,m,ireg,is,l)*dtemp(m,i,j,nol)
              enddo
              lwrz(k,l)=lwrz(k,l)+lwrt(k,18,ireg,is,l)*dtemp(18,i,j,nol)
            enddo

            dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism)
            do m=1,4
              do k=1,3
                lwrz(k,l)=lwrz(k,l)+ &
                     & (lwrts(k,m,ireg,is,l)+lwrqts(k,m,ireg,is,l)*dqa)*dumts
              enddo
!             lwrz(7,l)=lwrz(7,l)+
!    *        (lwrts(7,m,ireg,is,l)+lwrqts(7,m,ireg,is,l)*dqa)
!    *        *dumts
              dumts=dumts*(tsurfn(i,j,nn)-tncep(19,ireg,ism))
            enddo


          enddo

          if (nn.eq.1) then
            ulrad1nz(i,j)=0.
            if(initialization.eqv..true.) then
            !if (iyear.eq.0) then
             ulrad1nT=0.
            endif
          endif


          ulrad1nzz(i,j,nn)=(lwrz(2,0)+lwrz(5,0))*(1-tcc(i,j)) + &
               & (lwrz(2,1)+lwrz(5,1))*tcc(i,j)
         ulrad1nz(i,j)=ulrad1nz(i,j)+(ulrad1nzz(i,j,nn)*fractn(i,j,nn))

        enddo
      enddo


      if (nn.eq.3) then
      ulrad1nm=globalmean(ulrad1nz)

      ulrad1nmm=ulrad1nm-ulrad1nU

      ulrad1nT=ulrad1nT+(ulrad1nmm/(360.*6.))

      endif

! *** that's all folks
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine swaverad2(nn)
!-----------------------------------------------------------------------
! *** computes short wave radiation
! *** linearization of RCM with ISCCP D2 1990 clouds
!-----------------------------------------------------------------------



      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'

      integer i,j,k,l,ireg
      integer m, d, r, nn , nol


      real*8 f0,f1,ftot(8),fn(8,0:1)
      real*8 drs, drs2, drs3
      real*8 dcost, df,sk,sr,x,y,dfs,smsc,df2
      real*8 fswutoa(nlat,nlon),fswdsfc(nlat,nlon),fswusfc
      real*8 fswutoa2(nlat,nlon),fswdtoa(nlat,nlon)
      real*8 fswutoaG,fswdtoa2,fswdtoaG


      integer nreg(2)
!     real*8 zac(2),asup,bup
      real*8 zac(2),asup
      real*8  globalmean
      real*8 fswutoaGA,fswutoaG0
      real*8 fswutoa_diff,df_test
      common /rad_sul2 /fswutoa,fswdtoa
      common /rad_sul0 /fswutoaG,df_test,fswdtoaG
      common /pr_evap /fswdsfc
! *** aerosol scattering included as a correction on the upward
! *** clear sky fluxes
! *** sk,sr: empirical coefficients Dorland et al, J. Geophys. Res.,102,
! *** 28079-28100, 1997.
! *** smsc: mass scattering coefficient [m2/g]
! *** dso4: change in sulfate aerosol column integrated concentration since
! *** pre-industrial times [g/m2]


      sk=0.058d0*1370d0
      sr=0.05d0
      smsc=8.0

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
       do i=1,nlat
           alb2esn(i,j,nn)=albesn(i,j,nn)
           alb2esn(i,j,3) = albesnR(i,j)
          if (alb2esn(i,j,nn).ge.1.) then
            alb2esn(i,j,nn)=1.
          endif
          df=dayfr(i)*solarf(i)
          df2=dayfr(i)*solarm*solardref/1370.d0
          ireg=irn(i,j,nol)
          dcost=kosz(i)-costref(ireg,imonth)
          do l=1,8
            do k=0,1
              fn(l,k) = swrref(l,ireg,imonth,k) &
     &               +  swrcost(l,ireg,imonth,k)*dcost
            enddo
          enddo

          x=sqrt(kosz(i))
          y=sqrt(1-alb2esn(i,j,nn))
          dfs=sk*(4d0*x*y*(y-x)-sr)*dso4(i,j)*smsc
          if (dfs.gt.0d0.and.kosz(i).lt.0.05) dfs=0d0
          drs=alb2esn(i,j,nn)-salbref(ireg,imonth)
          drs2=drs*drs
          drs3=drs2*drs


          do l=1,4
           f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs+dfs
           f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo
          do l=5,8
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo


! alternative calculation of upward flux at ground:
! in parameterisation no cross terms are accounted for, which are important for
! upward shortwave radiation at surface and therefore also for net flux
! heswsn(i,j)

          ftot(4)=-alb2esn(i,j,nn)*ftot(8)

          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*df
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*df
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*df
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*df

! for diagnostic purposes:
! (1) downward shortwave radiation at TOA
             if (nn.eq.1)fswdtoa(i,j)=0.
             fswdtoa2=-ftot(5)*df2
             fswdtoa(i,j)=fswdtoa(i,j)+(fractn(i,j,nn)*fswdtoa2)
! (2) upward shortwave radiation at TOA
            if (nn.eq.1)fswutoa(i,j)=0.
            fswutoa2(i,j)=ftot(1)*df2
            fswutoa(i,j)=fswutoa(i,j)+(fractn(i,j,nn)*fswutoa2(i,j))
! (3) downward shortwave radiation at SURFACE
            fswdsfc(i,j)=-ftot(8)*df
! (4) upward shortwave radiation at SURFACE
!            fswusfc=(heswsn(i,j)+ftot(8)*df)

        enddo
      enddo
      if (nn.eq.3) then
        fswutoaG=globalmean(fswutoa)
        fswdtoaG=globalmean(fswdtoa)
      endif

      return
      end
