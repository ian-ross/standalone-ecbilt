!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmout(istep)
      implicit none

      include 'comatm.h'
      include 'comemic.h'

      integer istep

      if ( mod(istep,iatm) .eq. 0) then
        call selectout(istep)
      endif

      call outamean(istep)
      call outiocht(istep)

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine selectout(istep)
!-----------------------------------------------------------------------
! *** this routine selects which kind of outputs should be written to output
! *** itel counts the number of days in the output interval
! *** meantype = 1: output monthly mean fields.
! *** meantype = 2: output seasonal mean fields.
! *** meantot = 1: computes whole period monthly or seasonal mean fields.
! *** meanyl  = 1: computes yearly monthly or seasonal mean fields.
! *** ioutdaily = 1 output instantanous fields.
! *** instcount: counter for output instantaneous fields.
! *** ixout: frequency for output instantanous fields in days.
! *** written by xueli wang.
!-----------------------------------------------------------------------
      USE Atmosphere_Output
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comdiag.h'
      include 'comoutlocal.h'

      integer i,j,k,l,istep
      integer idmean,idstd
      real*8  facstr,costt,sintt,pfac,psifac,qpfac

      ivlevel(1) = 200
      ivlevel(2) = 500
      ivlevel(3) = 800
      itlevel(0) = 100
      itlevel(1) = 350
      itlevel(2) = 650
      itlevel(3) = 1000

!** some computation and unit transformation

      call sptogg(psi(1,1),psig(1,1,1),pp)
      call sptogg(psi(1,2),psig(1,1,2),pp)
      call sptogg(psi(1,3),psig(1,1,3),pp)
      call sptogg(qprime(1,1),qgpv(1,1,1),pp)
      call sptogg(qprime(1,2),qgpv(1,1,2),pp)
      call sptogg(qprime(1,3),qgpv(1,1,3),pp)

!  *** compute the precipitation, evaporation and runoffs in cm/year
      facstr = roair*uv10rws
      pfac=100.*3600.*24.*360.
      psifac=radius*radius*om
      qpfac=om
      do i=1,nlat
        costt=cos(dragane(i))
        sintt=sin(dragane(i))
        do j=1,nlon
	  do k=1,3
            qgpv(i,j,k) = qgpv(i,j,k)*qpfac
            psig(i,j,k) = psig(i,j,k)*psifac
	  enddo
          dyrain1(i,j) = (dyrain(i,j)+dysnow(i,j))*pfac
          corain1(i,j) = (corain(i,j)+cosnow(i,j))*pfac
          torain1(i,j) = (torain(i,j)+tosnow(i,j))*pfac
	  snow1(i,j)   = tosnow(i,j)*pfac
          evap1(i,j)   = evap(i,j)*pfac
          runofl1(i,j) = arunofl(i,j)*pfac
          runofo1(i,j) = arunofo(i,j)*pfac
          eminp1(i,j)  = evap1(i,j)-torain1(i,j)
          hesw(i,j)    = hesw0(i,j)+hesw1(i,j)+hesw2(i,j)+hesws(i,j)
          nlrads(i,j)  = ulrads(i,j)-dlrads(i,j)
          if (q0(i).gt.0d0) then
            albep(i,j)   = 1d0 - hesw(i,j)/q0(i)
          else
            albep(i,j)   = 1d0
          endif
          winstu1(i,j)= cdragw(i,j)*facstr*uvw10(i,j)* &
               & (utot(i,j,3)*costt-vtot(i,j,3)*sintt)
          winstv1(i,j)= cdragw(i,j)*facstr*uvw10(i,j)* &
               & (utot(i,j,3)*sintt+vtot(i,j,3)*costt)
        enddo
      enddo
!  *** change in compute temperature
      do i=1,nlat
        do j=1,nlon
          tsurf1(i,j) = tsurf(i,j)+newtotvar(Surface_Temperature,6)
          temp4g1(i,j)=temp4g(i,j)+newtotvar(Surface_Temperature,6)
          temp2g1(i,j)=temp2g(i,j)+newtotvar(Surface_Temperature,6)
          tempsg1(i,j)=tempsg(i,j)+newtotvar(Surface_Temperature,6)
          temp0g1(i,j)=temp0g(i,j)+newtotvar(Surface_Temperature,6)
        enddo
      enddo

      if(meantype.eq.2) then
        if(istep.gt.(11*30*iatm)) then
          iseason = iseason + 1
          if (iseason.gt.90) iseason = 1
          itel = itel + 1
        endif
      else
        itel = itel + 1
      endif

      instcount = instcount + 1

!  *** write the instant data
!      if (ioutdaily .eq. 1.and.instcount.eq.ixout) then
!        call outputinst
!      endif

      if (meantype .eq. 1) then
        minterv=30
        if (meantot .eq. 1) then
          call outputmtl
        endif
        if (meanyl .eq. 1) then
          call outputmyl
        endif
      endif

!      if (meantype .eq. 2 .and. istep .gt. 11*30*iatm) then
!        minterv=90
!        if (meantot .eq. 1) then
!          call outputmtl
!        endif
!        if (meanyl .eq. 1) then
!          call outputmyl
!        endif
!      endif

      if (ioutyearly .eq. 1) then
        call outputyrl
      endif

      if (itel.eq.minterv) itel = 0
      if (instcount.eq.ixout) instcount = 0

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
!      SUBROUTINE outputinst
!-----------------------------------------------------------------------
! *** this routine output instantanous fields.
! *** written by xueli wang, september 1995.
! *** adapted to netcdf format by camiel severijns, march 2000
!-----------------------------------------------------------------------
!      USE Atmosphere_Output
!      IMPLICIT NONE

!      INCLUDE 'comatm.h'
!      INCLUDE 'comdyn.h'
!      INCLUDE 'comphys.h'
!      INCLUDE 'comsurf.h'
!      INCLUDE 'comemic.h'
!      INCLUDE 'comcoup.h'
!      INCLUDE 'comdiag.h'
!      INCLUDE 'comoutlocal.h'

!      INTEGER, PARAMETER :: IOFlag = 1

!      REAL, DIMENSION(1:nlat,1:nlon,1:4) :: tmp

!      CALL open(Instantaneous_Data)

!      IF (output(newts(IOFlag)))
!     &     CALL write(Surface_Temperature,tsurf1(1:nlat,1:nlon))

!      IF (output(newt(1))) THEN
!         tmp(1:nlat,1:nlon,1) = temp0g1(:,:)
!         tmp(1:nlat,1:nlon,2) = temp2g1(:,:)
!         tmp(1:nlat,1:nlon,3) = temp4g1(:,:)
!         tmp(1:nlat,1:nlon,4) = tempsg1(:,:)
!         CALL write(Temperature,tmp(1:nlat,1:nlon,1:4))
!      END IF

!      IF (output(newtstrat(IOFlag)))
!     &     CALL write(Stratospheric_Temperature,temp0g1(1:nlat,1:nlon))

!      IF (output(newt2m(IOFlag)))
!     &     CALL write(Two_Meter_Temperature,tempsg1(1:nlat,1:nlon))

!      IF (output(newu(IOFlag))) THEN
!         tmp(1:nlat,1:nlon,1) = u800(:,:)
!         tmp(1:nlat,1:nlon,2) = u500(:,:)
!         tmp(1:nlat,1:nlon,3) = u200(:,:)
!         CALL write(Wind_U,tmp(1:nlat,1:nlon,1:3))
!      END IF

!      IF (output(newv(IOFlag))) THEN
!         tmp(1:nlat,1:nlon,1) = v800(:,:)
!         tmp(1:nlat,1:nlon,2) = v500(:,:)
!         tmp(1:nlat,1:nlon,3) = v200(:,:)
!         CALL write(Wind_V,tmp(1:nlat,1:nlon,1:3))
!      END IF

!      IF (output(newsp(IOFlag)))
!     &     CALL write(Surface_Pressure,pground(1:nlat,1:nlon))

!      IF (output(newomega(IOFlag)))
!     &     CALL write(Vertical_Pressure_Wind,omegg(1:nlat,1:nlon,1:nvl))

!      IF (output(newustress(IOFlag)))
!     &     CALL write(U_Stress,winstu1(1:nlat,1:nlon))
!      IF (output(newvstress(IOFlag)))
!     &     CALL write(V_Stress,winstv1(1:nlat,1:nlon))

!      IF (output(newuv10(IOFlag)))
!     &     CALL write(Wind_at_10_Meter,uv10(1:nlat,1:nlon))

!      IF (output(newdivu(IOFlag)))
!     &     CALL write(Ageostrophic_Wind_U,udivg(1:nlat,1:nlon,1:nvl))
!      IF (output(newdivv(IOFlag)))
!     &     CALL write(Ageostrophic_Wind_V,vdivg(1:nlat,1:nlon,1:nvl))

!      IF (output(newpsi(IOFlag)))
!     &     CALL write(Stream_Function,psig(1:nlat,1:nlon,1:nvl))

!      IF (output(newchi(IOFlag)))
!     &     CALL write(Velocity_Potential,chig(1:nlat,1:nlon,1:nvl))

!      IF (output(newqgpv(IOFlag)))
!     &     CALL write(QG_Potential_Vorticity,qgpv(1:nlat,1:nlon,1:nvl))

!      IF (output(newgh(IOFlag)))
!     &     CALL write(Geopotential_Height,geopg(1:nlat,1:nlon,1:nvl))

!      IF (output(newhforg(IOFlag))) THEN
!         tmp(1:nlat,1:nlon,1) = vhforg1(:,:)
!         tmp(1:nlat,1:nlon,2) = vhforg2(:,:)
!         CALL write(Heating_Force,tmp(1:nlat,1:nlon,1:2))
!      END IF

!      IF (output(newvforg(IOFlag))) THEN
!         tmp(1:nlat,1:nlon,1) = vforg1(:,:)
!         tmp(1:nlat,1:nlon,2) = vforg2(:,:)
!         tmp(1:nlat,1:nlon,3) = vforg3(:,:)
!         CALL write(Potential_Vorticity_Forcing,tmp(1:nlat,1:nlon,1:3))
!      END IF

!      IF (output(newdyrain(IOFlag)))
!     &     CALL write(Large_Scale_Precipitation,dyrain1(1:nlat,1:nlon))

!      IF (output(newcorain(IOFlag)))
!     &     CALL write(Convective_Precipitation,corain1(1:nlat,1:nlon))

!      IF (output(newtorain(IOFlag)))
!     &     CALL write(Total_Precipitation,torain1(1:nlat,1:nlon))

!      IF (output(newsnow(IOFlag)))
!     &     CALL write(Total_Snow_Fall,snow1(1:nlat,1:nlon))

!      IF (output(newevap(IOFlag)))
!     &     CALL write(Surface_Evaporation,evap1(1:nlat,1:nlon))

!      IF (output(neweminp(IOFlag)))
!     &     CALL write(Evap_Minus_Precip,eminp1(1:nlat,1:nlon))

!      IF (output(newhflux(IOFlag)))
!     &     CALL write(Surface_Sensible_Heat_Flux,hflux(1:nlat,1:nlon))

!      IF (output(neweflux(IOFlag)))
!     &     CALL write(Surface_Latent_Heat_Flux,eflux(1:nlat,1:nlon))

!      IF (output(newtsr(IOFlag)))
!     &     CALL write(Top_Solar_Radiation,hesw(1:nlat,1:nlon))

!      IF (output(newssr(IOFlag)))
!     &     CALL write(Surface_Solar_Radiation,hesws(1:nlat,1:nlon))

!      IF (output(newstr(IOFlag)))
!     &     CALL write(Surface_Thermal_Radiation,nlrads(1:nlat,1:nlon))

!      IF (output(newttr(IOFlag)))
!     &     CALL write(Top_Thermal_Radiation,ulrad0(1:nlat,1:nlon))

!      IF (output(newalbs(IOFlag)))
!     &     CALL write(Surface_Albedo,alb2es(1:nlat,1:nlon))

!      IF (output(newalbp(IOFlag)))
!     &     CALL write(Planetary_Albedo,albep(1:nlat,1:nlon))

!      IF (output(newbmoisg(IOFlag)))
!     &     CALL write(Bottom_Moisture,abmoisg(1:nlat,1:nlon))

!      IF (output(newrunoffl(IOFlag)))
!     &     CALL write(Land_Surface_Runoff,runofl1(1:nlat,1:nlon))

!      IF (output(newrunoffo(IOFlag)))
!     &     CALL write(Ocean_Surface_Runoff,runofo1(1:nlat,1:nlon))

!      IF (output(newsdl(IOFlag)))
!     &     CALL write(Land_Snow_Depth,adsnow(1:nlat,1:nlon))

!      IF (output(newhic(IOFlag)))
!     &     CALL write(Sea_Ice_Thickness,ahic(1:nlat,1:nlon))

!      IF (output(newrmoisg(IOFlag)))
!     &     CALL write(Specific_Humidity,rmoisg(1:nlat,1:nlon))

!      IF (output(newrelhum(IOFlag)))
!     &     CALL write(Relative_Humidity,relhum(1:nlat,1:nlon))

!      IF (output(newcdragw(IOFlag)))
!     &     CALL write(Drag_Coefficient_W,cdragw(1:nlat,1:nlon))

!      IF (output(newcdragv(IOFlag)))
!     &     CALL write(Drag_Coefficient_V,cdragv(1:nlat,1:nlon))

!      IF (output(newtcc(IOFlag)))
!     &     CALL write(Total_Cloud_Cover,tccd(1:nlat,1:nlon))

!      IF (output(newdumt1(IOFlag)))
!     &     CALL write(User_Assigned_T1,dumt1(1:nlat,1:nlon,1:nvl))

!      IF (output(newdumt2(IOFlag)))
!     &     CALL write(User_Assigned_T2,dumt2(1:nlat,1:nlon,1:nvl))

!      IF (output(newdumu1(IOFlag)))
!     &     CALL write(User_Assigned_U1,dumu1(1:nlat,1:nlon,1:nvl))

!      IF (output(newdumu2(IOFlag)))
!     &     CALL write(User_Assigned_U2,dumu2(1:nlat,1:nlon,1:nvl))

!      CALL close

!      RETURN
!      END SUBROUTINE outputinst


!23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE outputmtl
!-----------------------------------------------------------------------
! *** this routine calls another routine which computes the  whole period
! *** monthly or seasonal mean fields.
! *** written by xueli wang and nanne weber, april 1995.
!-----------------------------------------------------------------------
      USE Atmosphere_Output
      IMPLICIT NONE

      INCLUDE 'comatm.h'
      INCLUDE 'comdyn.h'
      INCLUDE 'comphys.h'
      INCLUDE 'comsurf.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comcoup.h'
      INCLUDE 'comdiag.h'
      INCLUDE 'comoutlocal.h'

      INTERFACE
         SUBROUTINE totstat(xx,sumxx,sumxxsq,xxm,xxdev,compute)
         IMPLICIT NONE
         REAL*8, DIMENSION(:,:),   INTENT(in)    :: xx
         REAL*8, DIMENSION(:,:,:), INTENT(inout) :: sumxx, sumxxsq
         REAL*8, DIMENSION(:,:),   INTENT(out)   :: xxm, xxdev
         LOGICAL,                  INTENT(in)    :: compute
         END SUBROUTINE totstat
      END INTERFACE

      INTEGER, PARAMETER :: IOFlag = 3

      LOGICAL :: need_to_write
      REAL*8, DIMENSION(1:nlat,1:nlon,1:4) :: mean, stddev
      INTEGER :: realyear,realmonth

      realyear=iyear
      realmonth=imonth
      if ( irunlabeld/=360 ) realyear=iyear+1
      if ( irunlabeld/=360 ) realmonth=(irunlabeld/30)+1+(ndays/30)

      IF (itel == 0) RETURN

      IF (meantype == 1) THEN
         need_to_write = ( realyear == nyears .AND. imonth .GE. realmonth .AND. iday == 30 )
         IF (.NOT. need_to_write) need_to_write = ( iyear == nyears .AND. imonth .LT. realmonth .AND. iday == 30 )
         IF (need_to_write) CALL open(Total_Monthly_Means)
      ELSE
         need_to_write = ( realyear == nyears .AND. imonth .GE. realmonth .AND. iseason == 90 )
         IF (.NOT. need_to_write) need_to_write = ( iyear == nyears .AND. imonth .LT. realmonth .AND. iseason == 90 )
         IF (need_to_write) CALL open(Total_Seasonal_Means)
      END IF

      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
         CALL totstat(tsurf1(1:nlat,1:nlon), &
     &        sxtsurf(1:nlat,1:nlon,1:12), sytsurf(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),    stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         CALL totstat(temp0g1(1:nlat,1:nlon), &
     &        sxtstrat(1:nlat,1:nlon,1:12),sytstrat(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(temp2g1(1:nlat,1:nlon), &
     &        sxtemp2g(1:nlat,1:nlon,1:12),sytemp2g(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(temp4g1(1:nlat,1:nlon), &
     &        sxtemp4g(1:nlat,1:nlon,1:12),sytemp4g(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         CALL totstat(tempsg1(1:nlat,1:nlon), &
     &        sxtempsg(1:nlat,1:nlon,1:12),sytempsg(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,4),   stddev  (1:nlat,1:nlon,4), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
         END IF
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         CALL totstat(temp0g1(1:nlat,1:nlon), &
     &        sxtstrat(1:nlat,1:nlon,1:12), sytstrat(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),    stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         CALL totstat(tempsg1(1:nlat,1:nlon), &
     &        sxt2m(1:nlat,1:nlon,1:12),syt2m (1:nlat,1:nlon,1:12), &
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         CALL totstat(u200(1:nlat,1:nlon), &
     &        sxu200(1:nlat,1:nlon,1:12),syu200(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(u500(1:nlat,1:nlon), &
     &        sxu500(1:nlat,1:nlon,1:12),syu500(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(u800(1:nlat,1:nlon), &
     &        sxu800(1:nlat,1:nlon,1:12),syu800(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         CALL totstat(v200(1:nlat,1:nlon), &
     &        sxv200(1:nlat,1:nlon,1:12),syv200(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(v500(1:nlat,1:nlon), &
     &        sxv500(1:nlat,1:nlon,1:12),syv500(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(v800(1:nlat,1:nlon), &
     &        sxv800(1:nlat,1:nlon,1:12),syv800(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         CALL totstat(pground(1:nlat,1:nlon), &
     &        sxpground(1:nlat,1:nlon,1:12),sypground(1:nlat,1:nlon,1:12), &
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         CALL totstat(omegg(1:nlat,1:nlon,1), &
     &        sxomeg1(1:nlat,1:nlon,1:12),syomeg1(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(omegg(1:nlat,1:nlon,2), &
     &        sxomeg2(1:nlat,1:nlon,1:12),syomeg2(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,2),   stddev (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(omegg(1:nlat,1:nlon,3), &
     &        sxomeg3(1:nlat,1:nlon,1:12),syomeg3(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,3),   stddev (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         CALL totstat(winstu1(1:nlat,1:nlon), &
     &        sxwinstu1(1:nlat,1:nlon,1:12),sywinstu1(1:nlat,1:nlon,1:12), &
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         CALL totstat(winstv1(1:nlat,1:nlon), &
     &        sxwinstv1(1:nlat,1:nlon,1:12),sywinstv1(1:nlat,1:nlon,1:12), &
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         CALL totstat(uv10(1:nlat,1:nlon), &
     &        sxuv10(1:nlat,1:nlon,1:12),syuv10(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         CALL totstat(udivg(1:nlat,1:nlon,1), &
     &        sxudivg1(1:nlat,1:nlon,1:12),syudivg1(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(udivg(1:nlat,1:nlon,2), &
     &        sxudivg2(1:nlat,1:nlon,1:12),syudivg2(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(udivg(1:nlat,1:nlon,3), &
     &        sxudivg3(1:nlat,1:nlon,1:12),syudivg3(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         CALL totstat(vdivg(1:nlat,1:nlon,1), &
     &        sxvdivg1(1:nlat,1:nlon,1:12),syvdivg1(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(vdivg(1:nlat,1:nlon,2), &
     &        sxvdivg2(1:nlat,1:nlon,1:12),syvdivg2(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(vdivg(1:nlat,1:nlon,3), &
     &        sxvdivg3(1:nlat,1:nlon,1:12),syvdivg3(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         CALL totstat(psig(1:nlat,1:nlon,1), &
     &        sxgrpsi1(1:nlat,1:nlon,1:12),sygrpsi1(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(psig(1:nlat,1:nlon,2), &
     &        sxgrpsi2(1:nlat,1:nlon,1:12),sygrpsi2(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(psig(1:nlat,1:nlon,3), &
     &        sxgrpsi3(1:nlat,1:nlon,1:12),sygrpsi3(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         CALL totstat(chig(1:nlat,1:nlon,1), &
     &        sxchi1(1:nlat,1:nlon,1:12),sychi1(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(chig(1:nlat,1:nlon,1), &
     &        sxchi2(1:nlat,1:nlon,1:12),sychi2(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(chig(1:nlat,1:nlon,1), &
     &        sxchi3(1:nlat,1:nlon,1:12),sychi3(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         CALL totstat(qgpv(1:nlat,1:nlon,1), &
     &        sxqgpv1(1:nlat,1:nlon,1:12),syqgpv1(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(qgpv(1:nlat,1:nlon,2), &
     &        sxqgpv2(1:nlat,1:nlon,1:12),syqgpv2(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,2),   stddev (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(qgpv(1:nlat,1:nlon,3), &
     &        sxqgpv3(1:nlat,1:nlon,1:12),syqgpv3(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,3),   stddev (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         CALL totstat(geopg(1:nlat,1:nlon,1), &
     &        sxgh1(1:nlat,1:nlon,1:12),sygh1 (1:nlat,1:nlon,1:12), &
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(geopg(1:nlat,1:nlon,2), &
     &        sxgh2(1:nlat,1:nlon,1:12),sygh2 (1:nlat,1:nlon,1:12), &
     &        mean (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(geopg(1:nlat,1:nlon,3), &
     &        sxgh3(1:nlat,1:nlon,1:12),sygh3 (1:nlat,1:nlon,1:12), &
     &        mean (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         CALL totstat(vhforg1(1:nlat,1:nlon), &
     &        sxvhforg1(1:nlat,1:nlon,1:12),syvhforg1(1:nlat,1:nlon,1:12), &
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(vhforg2(1:nlat,1:nlon), &
     &        sxvhforg2(1:nlat,1:nlon,1:12),syvhforg2(1:nlat,1:nlon,1:12), &
     &        mean     (1:nlat,1:nlon,2),   stddev   (1:nlat,1:nlon,2), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
         END IF
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         CALL totstat(vforg1(1:nlat,1:nlon), &
     &        sxvforg1(1:nlat,1:nlon,1:12),syvforg1(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(vforg2(1:nlat,1:nlon), &
     &        sxvforg2(1:nlat,1:nlon,1:12),syvforg2(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(vforg3(1:nlat,1:nlon), &
     &        sxvforg3(1:nlat,1:nlon,1:12),syvforg3(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         CALL totstat(dyrain1(1:nlat,1:nlon), &
     &        sxdyrain(1:nlat,1:nlon,1:12),sydyrain(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         CALL totstat(corain1(1:nlat,1:nlon), &
     &        sxcorain(1:nlat,1:nlon,1:12),sycorain(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         CALL totstat(torain1(1:nlat,1:nlon), &
     &        sxtorain(1:nlat,1:nlon,1:12),sytorain(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         CALL totstat(snow1(1:nlat,1:nlon), &
     &        sxsnow(1:nlat,1:nlon,1:12),sysnow(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         CALL totstat(evap1(1:nlat,1:nlon), &
     &        sxevap(1:nlat,1:nlon,1:12),syevap(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         CALL totstat(eminp1(1:nlat,1:nlon), &
     &        sxeminp(1:nlat,1:nlon,1:12),syeminp(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         CALL totstat(hflux(1:nlat,1:nlon), &
     &        sxhflux(1:nlat,1:nlon,1:12),syhflux(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
          CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         CALL totstat(eflux(1:nlat,1:nlon), &
     &        sxeflux(1:nlat,1:nlon,1:12),syeflux(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         CALL totstat(hesws(1:nlat,1:nlon), &
     &        sxhesws(1:nlat,1:nlon,1:12),syhesws(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         CALL totstat(hesw(1:nlat,1:nlon), &
     &        sxhesw(1:nlat,1:nlon,1:12),syhesw(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         CALL totstat(alb2es(1:nlat,1:nlon), &
     &        sxalbes(1:nlat,1:nlon,1:12),syalbes(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         CALL totstat(albep(1:nlat,1:nlon), &
     &        sxalbep(1:nlat,1:nlon,1:12),syalbep(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         CALL totstat(nlrads(1:nlat,1:nlon), &
     &        sxnlrads(1:nlat,1:nlon,1:12),synlrads(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         CALL totstat(ulrad0(1:nlat,1:nlon), &
     &        sxulrad1(1:nlat,1:nlon,1:12),syulrad1(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
         CALL totstat(abmoisg(1:nlat,1:nlon), &
     &        sxbmoisg(1:nlat,1:nlon,1:12),sybmoisg(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
         CALL totstat(adsnow(1:nlat,1:nlon), &
     &        sxdsnow(1:nlat,1:nlon,1:12),sydsnow(1:nlat,1:nlon,1:12), &
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         CALL totstat(ahic(1:nlat,1:nlon), &
     &        sxhic(1:nlat,1:nlon,1:12),syhic (1:nlat,1:nlon,1:12), &
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         CALL totstat(runofl1(1:nlat,1:nlon), &
     &        sxrunofl(1:nlat,1:nlon,1:12),syrunofl(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         CALL totstat(runofo1(1:nlat,1:nlon), &
     &        sxrunofo(1:nlat,1:nlon,1:12),syrunofo(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
             CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
         CALL totstat(rmoisg(1:nlat,1:nlon), &
     &        sxrmoisgw3(1:nlat,1:nlon,1:12),syrmoisgw3(1:nlat,1:nlon,1:12), &
     &        mean      (1:nlat,1:nlon,1),   stddev    (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         CALL totstat(relhum(1:nlat,1:nlon), &
     &        sxrelhum(1:nlat,1:nlon,1:12),syrelhum(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         CALL totstat(cdragw(1:nlat,1:nlon), &
     &        sxcdragw(1:nlat,1:nlon,1:12),sycdragw(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         CALL totstat(cdragv(1:nlat,1:nlon), &
     &        sxcdragv(1:nlat,1:nlon,1:12),sycdragv(1:nlat,1:nlon,1:12), &
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         CALL totstat(tccd(1:nlat,1:nlon), &
     &        sxtcc(1:nlat,1:nlon,1:12),sytcc (1:nlat,1:nlon,1:12), &
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         CALL totstat(dumt1(1:nlat,1:nlon,1), &
     &        sxdt11(1:nlat,1:nlon,1:12),sydt11(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(dumt1(1:nlat,1:nlon,2), &
     &        sxdt12(1:nlat,1:nlon,1:12),sydt12(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(dumt1(1:nlat,1:nlon,3), &
     &        sxdt13(1:nlat,1:nlon,1:12),sydt13(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         CALL totstat(dumt2(1:nlat,1:nlon,1), &
     &        sxdt21(1:nlat,1:nlon,1:12),sydt21(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(dumt2(1:nlat,1:nlon,2), &
     &        sxdt22(1:nlat,1:nlon,1:12),sydt22(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(dumt2(1:nlat,1:nlon,3), &
     &        sxdt23(1:nlat,1:nlon,1:12),sydt23(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         CALL totstat(dumu1(1:nlat,1:nlon,1), &
     &        sxdu11(1:nlat,1:nlon,1:12),sydu11(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(dumu1(1:nlat,1:nlon,2), &
     &        sxdu12(1:nlat,1:nlon,1:12),sydu12(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(dumu1(1:nlat,1:nlon,3), &
     &        sxdu13(1:nlat,1:nlon,1:12),sydu13(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         CALL totstat(dumu2(1:nlat,1:nlon,1), &
     &        sxdu21(1:nlat,1:nlon,1:12),sydu21(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL totstat(dumu2(1:nlat,1:nlon,2), &
     &        sxdu22(1:nlat,1:nlon,1:12),sydu22(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL totstat(dumu2(1:nlat,1:nlon,3), &
     &        sxdu23(1:nlat,1:nlon,1:12),sydu23(1:nlat,1:nlon,1:12), &
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (need_to_write) CALL close

      RETURN
      END

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE totstat(x,sumx1,sumy1,xmean,xstd,compute)
! *** --------------------------------------------------------------------
! *** This routine computes whole period monthly or seasonal mean and standard
! *** deviation around it.
! *** all arrays are assumed to have the same shape.
! *** ------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:),   INTENT(in)    :: x
      REAL*8, DIMENSION(:,:,:), INTENT(inout) :: sumx1, sumy1
      REAL*8, DIMENSION(:,:),   INTENT(out)   :: xmean, xstd
      LOGICAL,                  INTENT(in)    :: compute

      INCLUDE 'comatm.h'
      INCLUDE 'comcoup.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comdiag.h'

      INTEGER :: k,ncase,nk
      INTEGER :: realyear,realmonth

      realyear=iyear
      realmonth=1
      if ( irunlabeld/=360 ) realyear=iyear+1
      if ( irunlabeld/=360 ) realmonth=(irunlabeld/30)+1
      IF (meantype == 1) THEN
         k     = imonth
         ncase = nyears*30
         nk    = 12
         IF (realyear == 1 .AND. imonth == realmonth .AND. iday == 1) THEN
            sumx1(:,:,1:nk)   = 0.0
            sumy1(:,:,1:nk)   = 0.0
         END IF
      ELSE                      ! IF (meantype == 2) THEN
         k     = Mod(imonth/3,4) + 1
         ncase = (nyears-1)*90
         nk    = 4
         IF (realyear == 1 .AND. imonth == 12 .AND. iday == 1) then
            sumx1(:,:,1:nk)   = 0.0
            sumy1(:,:,1:nk)   = 0.0
         END IF
      END IF

      sumx1(:,:,k) = sumx1(:,:,k) + x(:,:)
      sumy1(:,:,k) = sumy1(:,:,k) + x(:,:)**2

      IF (compute) THEN
         !IF ( iyear == nyears .AND. imonth .GE. realmonth .AND. iday == 30 ) ncase=ncase+30
         xmean(:,:) = sumx1(:,:,k)/ncase
         xstd (:,:) = sumy1(:,:,k) - sumx1(:,:,k)**2/ncase
         WHERE (xstd(:,:) <= 0.0)
            xstd(:,:) = 0.0
         ELSE WHERE
            xstd(:,:) = Sqrt(xstd(:,:)/(ncase - 1.0))
         END WHERE
      END IF

      RETURN
      END SUBROUTINE totstat

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE outputmyl
!-----------------------------------------------------------------------
! *** this routine calls another routine which computes the yearly monthly
! *** mean fields.
! *** written by xueli wang, september 1995.
!-----------------------------------------------------------------------
      USE Atmosphere_Output
      IMPLICIT NONE

      INCLUDE 'comatm.h'
      INCLUDE 'comdyn.h'
      INCLUDE 'comphys.h'
      INCLUDE 'comsurf.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comcoup.h'
      INCLUDE 'comdiag.h'
      INCLUDE 'comoutlocal.h'

      INTERFACE
         SUBROUTINE stat(xx,sumxx,sumxxsq,xxm,xxdev,compute)
         IMPLICIT NONE
         REAL*8, DIMENSION(:,:), INTENT(in)    :: xx
         REAL*8, DIMENSION(:,:), INTENT(inout) :: sumxx, sumxxsq
         REAL*8, DIMENSION(:,:), INTENT(out)   :: xxm, xxdev
         LOGICAL,                INTENT(in)    :: compute
         END SUBROUTINE stat
      END INTERFACE

      INTEGER, PARAMETER :: IOFlag = 2

      LOGICAL :: need_to_write
      REAL*8, DIMENSION(1:nlat,1:nlon,1:4) :: mean, stddev

      IF (itel == 0) RETURN

      need_to_write = ( itel == minterv )

      IF (need_to_write) THEN
         IF (meantype == 1) THEN
            CALL open(Monthly_Means)
         ELSE
            CALL open(Seasonal_Means)
         END IF
      END IF

      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
         CALL stat(tsurf1(1:nlat,1:nlon), &
     &        s1tsurf(1:nlat,1:nlon),  s2tsurf(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         CALL stat(temp0g1(1:nlat,1:nlon), &
     &        s1tstrat(1:nlat,1:nlon),  s2tstrat(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(temp2g1(1:nlat,1:nlon), &
     &        s1temp2g(1:nlat,1:nlon),  s2temp2g(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(temp4g1(1:nlat,1:nlon), &
     &        s1temp4g(1:nlat,1:nlon),  s2temp4g(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         CALL stat(tempsg1(1:nlat,1:nlon), &
     &        s1tempsg(1:nlat,1:nlon),  s2tempsg(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,4),stddev  (1:nlat,1:nlon,4), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
         END IF
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         CALL stat(temp0g1(1:nlat,1:nlon), &
     &        s1tstrat(1:nlat,1:nlon),  s2tstrat(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         CALL stat(tempsg1(1:nlat,1:nlon), &
     &        s1t2m(1:nlat,1:nlon),  s2t2m (1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         CALL stat(u200(1:nlat,1:nlon), &
     &        s1u200(1:nlat,1:nlon),  s2u200(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(u500(1:nlat,1:nlon), &
     &        s1u500(1:nlat,1:nlon),  s2u500(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(u800(1:nlat,1:nlon), &
     &        s1u800(1:nlat,1:nlon),  s2u800(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         CALL stat(v200(1:nlat,1:nlon), &
     &        s1v200(1:nlat,1:nlon),  s2v200(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(v500(1:nlat,1:nlon), &
     &        s1v500(1:nlat,1:nlon),  s2v500(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(v800(1:nlat,1:nlon), &
     &        s1v800(1:nlat,1:nlon),  s2v800(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         CALL stat(pground(1:nlat,1:nlon), &
     &        s1pground(1:nlat,1:nlon),  s2pground(1:nlat,1:nlon), &
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         CALL stat(omegg(1:nlat,1:nlon,1), &
     &        s1omeg(1:nlat,1:nlon,1),s2omeg(1:nlat,1:nlon,1), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(omegg(1:nlat,1:nlon,2), &
     &        s1omeg(1:nlat,1:nlon,2),s2omeg(1:nlat,1:nlon,2), &
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(omegg(1:nlat,1:nlon,3), &
     &        s1omeg(1:nlat,1:nlon,3),s2omeg(1:nlat,1:nlon,3), &
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         CALL stat(winstu1(1:nlat,1:nlon), &
     &        s1winstu1(1:nlat,1:nlon),  s2winstu1(1:nlat,1:nlon), &
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         CALL stat(winstv1(1:nlat,1:nlon), &
     &        s1winstv1(1:nlat,1:nlon),  s2winstv1(1:nlat,1:nlon), &
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         CALL stat(uv10(1:nlat,1:nlon), &
     &        s1uv10(1:nlat,1:nlon),  s2uv10(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         CALL stat(udivg(1:nlat,1:nlon,1), &
     &        s1udivg(1:nlat,1:nlon,1),s2udivg(1:nlat,1:nlon,1), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(udivg(1:nlat,1:nlon,2), &
     &        s1udivg(1:nlat,1:nlon,2),s2udivg(1:nlat,1:nlon,2), &
     &        mean   (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(udivg(1:nlat,1:nlon,3), &
     &        s1udivg(1:nlat,1:nlon,3),s2udivg(1:nlat,1:nlon,3), &
     &        mean   (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         CALL stat(vdivg(1:nlat,1:nlon,1), &
     &        s1vdivg(1:nlat,1:nlon,1),s2vdivg(1:nlat,1:nlon,1), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(vdivg(1:nlat,1:nlon,2), &
     &        s1vdivg(1:nlat,1:nlon,2),s2vdivg(1:nlat,1:nlon,2), &
     &        mean   (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(vdivg(1:nlat,1:nlon,3), &
     &        s1vdivg(1:nlat,1:nlon,3),s2vdivg(1:nlat,1:nlon,3), &
     &        mean   (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         CALL stat(psig(1:nlat,1:nlon,1), &
     &        s1psi(1:nlat,1:nlon,1),s2psi (1:nlat,1:nlon,1), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(psig(1:nlat,1:nlon,2), &
     &        s1psi(1:nlat,1:nlon,2),s2psi (1:nlat,1:nlon,2), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(psig(1:nlat,1:nlon,3), &
     &        s1psi(1:nlat,1:nlon,3),s2psi (1:nlat,1:nlon,3), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         CALL stat(chig(1:nlat,1:nlon,1), &
     &        s1chi(1:nlat,1:nlon,1),s2chi (1:nlat,1:nlon,1), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(chig(1:nlat,1:nlon,2), &
     &        s1chi(1:nlat,1:nlon,2),s2chi (1:nlat,1:nlon,2), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(chig(1:nlat,1:nlon,3), &
     &        s1chi(1:nlat,1:nlon,3),s2chi (1:nlat,1:nlon,3), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         CALL stat(qgpv(1:nlat,1:nlon,1), &
     &        s1qgpv(1:nlat,1:nlon,1),s2qgpv(1:nlat,1:nlon,1), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(qgpv(1:nlat,1:nlon,2), &
     &        s1qgpv(1:nlat,1:nlon,2),s2qgpv(1:nlat,1:nlon,2), &
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(qgpv(1:nlat,1:nlon,3), &
     &        s1qgpv(1:nlat,1:nlon,3),s2qgpv(1:nlat,1:nlon,3), &
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         CALL stat(geopg(1:nlat,1:nlon,1), &
     &        s1gh(1:nlat,1:nlon,1),s2gh  (1:nlat,1:nlon,1), &
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(geopg(1:nlat,1:nlon,2), &
     &        s1gh(1:nlat,1:nlon,2),s2gh  (1:nlat,1:nlon,2), &
     &        mean(1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(geopg(1:nlat,1:nlon,3), &
     &        s1gh(1:nlat,1:nlon,3),s2gh  (1:nlat,1:nlon,3), &
     &        mean(1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         CALL stat(vhforg1(1:nlat,1:nlon), &
     &        s1vhforg1(1:nlat,1:nlon),  s2vhforg1(1:nlat,1:nlon), &
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(vhforg2(1:nlat,1:nlon), &
     &        s1vhforg2(1:nlat,1:nlon),  s2vhforg2(1:nlat,1:nlon), &
     &        mean     (1:nlat,1:nlon,2),stddev   (1:nlat,1:nlon,2), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
         END IF
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         CALL stat(vforg1(1:nlat,1:nlon), &
     &        s1vforg1(1:nlat,1:nlon),  s2vforg1(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(vforg2(1:nlat,1:nlon), &
     &        s1vforg2(1:nlat,1:nlon),  s2vforg2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(vforg3(1:nlat,1:nlon), &
     &        s1vforg3(1:nlat,1:nlon),  s2vforg3(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         CALL stat(dyrain1(1:nlat,1:nlon), &
     &        s1dyrain(1:nlat,1:nlon),  s2dyrain(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         CALL stat(corain1(1:nlat,1:nlon), &
     &        s1corain(1:nlat,1:nlon),  s2corain(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         CALL stat(torain1(1:nlat,1:nlon), &
     &        s1torain(1:nlat,1:nlon),  s2torain(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         CALL stat(snow1(1:nlat,1:nlon), &
     &        s1snow(1:nlat,1:nlon),  s2snow(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
             CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         CALL stat(evap1(1:nlat,1:nlon), &
     &        s1evap(1:nlat,1:nlon),  s2evap(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         CALL stat(eminp1(1:nlat,1:nlon), &
     &        s1eminp(1:nlat,1:nlon),  s2eminp(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         CALL stat(hflux(1:nlat,1:nlon), &
     &        s1hflux(1:nlat,1:nlon),  s2hflux(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
          CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         CALL stat(eflux(1:nlat,1:nlon), &
     &        s1eflux(1:nlat,1:nlon),  s2eflux(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         CALL stat(hesws(1:nlat,1:nlon), &
     &        s1hesws(1:nlat,1:nlon),  s2hesws(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         CALL stat(hesw(1:nlat,1:nlon), &
     &        s1hesw(1:nlat,1:nlon),  s2hesw(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         CALL stat(ulrad0(1:nlat,1:nlon), &
     &        s1ulrad1(1:nlat,1:nlon),  s2ulrad1(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         CALL stat(nlrads(1:nlat,1:nlon), &
     &        s1nlrads(1:nlat,1:nlon),  s2nlrads(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         CALL stat(alb2es(1:nlat,1:nlon), &
     &        s1albes(1:nlat,1:nlon),  s2albes(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         CALL stat(albep(1:nlat,1:nlon), &
     &        s1albep(1:nlat,1:nlon),  s2albep(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
         CALL stat(abmoisg(1:nlat,1:nlon), &
     &        s1bmoisg(1:nlat,1:nlon),  s2bmoisg(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
         CALL stat(adsnow(1:nlat,1:nlon), &
     &        s1dsnow(1:nlat,1:nlon),  s2dsnow(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         CALL stat(ahic(1:nlat,1:nlon), &
     &        s1hic(1:nlat,1:nlon),  s2hic (1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         CALL stat(runofl1(1:nlat,1:nlon), &
     &        s1runofl(1:nlat,1:nlon),  s2runofl(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         CALL stat(runofo1(1:nlat,1:nlon), &
     &        s1runofo(1:nlat,1:nlon),  s2runofo(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
         CALL stat(rmoisg(1:nlat,1:nlon), &
     &        s1rmoisgw3(1:nlat,1:nlon),  s2rmoisgw3(1:nlat,1:nlon), &
     &        mean      (1:nlat,1:nlon,1),stddev    (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         CALL stat(relhum(1:nlat,1:nlon), &
     &        s1relhum(1:nlat,1:nlon),  s2relhum(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         CALL stat(cdragw(1:nlat,1:nlon), &
     &        s1cdragw(1:nlat,1:nlon),  s2cdragw(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         CALL stat(cdragv(1:nlat,1:nlon), &
     &        s1cdragv(1:nlat,1:nlon),  s2cdragv(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         CALL stat(tccd(1:nlat,1:nlon), &
     &        s1tcc(1:nlat,1:nlon),  s2tcc (1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         CALL stat(dumt1(1:nlat,1:nlon,1), &
     &        s1dt1(1:nlat,1:nlon,1),s2dt1 (1:nlat,1:nlon,1), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(dumt1(1:nlat,1:nlon,2), &
     &        s1dt1(1:nlat,1:nlon,2),s2dt1 (1:nlat,1:nlon,2), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(dumt1(1:nlat,1:nlon,3), &
     &        s1dt1(1:nlat,1:nlon,3),s2dt1 (1:nlat,1:nlon,3), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         CALL stat(dumt2(1:nlat,1:nlon,1), &
     &        s1dt2(1:nlat,1:nlon,1),s2dt2 (1:nlat,1:nlon,1), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(dumt2(1:nlat,1:nlon,2), &
     &        s1dt2(1:nlat,1:nlon,2),s2dt2 (1:nlat,1:nlon,2), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(dumt2(1:nlat,1:nlon,3), &
     &        s1dt2(1:nlat,1:nlon,3),s2dt2 (1:nlat,1:nlon,3), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         CALL stat(dumu1(1:nlat,1:nlon,1), &
     &        s1du1(1:nlat,1:nlon,1),s2du1 (1:nlat,1:nlon,1), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(dumu1(1:nlat,1:nlon,2), &
     &        s1du1(1:nlat,1:nlon,2),s2du1 (1:nlat,1:nlon,2), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(dumu1(1:nlat,1:nlon,3), &
     &        s1du1(1:nlat,1:nlon,3),s2du1 (1:nlat,1:nlon,3), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         CALL stat(dumu2(1:nlat,1:nlon,1), &
     &        s1du2(1:nlat,1:nlon,1),s2du2 (1:nlat,1:nlon,1), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write)
         CALL stat(dumu2(1:nlat,1:nlon,2), &
     &        s1du2(1:nlat,1:nlon,2),s2du2 (1:nlat,1:nlon,2), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write)
         CALL stat(dumu2(1:nlat,1:nlon,3), &
     &        s1du2(1:nlat,1:nlon,3),s2du2 (1:nlat,1:nlon,3), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (need_to_write) CALL close

      RETURN
      END SUBROUTINE outputmyl

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE stat(xx,sumxx,sumxxsq,xxm,xxdev,compute)
!-----------------------------------------------------------------------
! *** this function sums xx and its square and computes the means and
! *** standard deviation at the end of a summation interval. The return
! *** value of this function indicates whether the output arrays xxm and
! *** xxdev contain valid data.
! *** all arrays are assumed to have the same shape.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:), INTENT(in)    :: xx
      REAL*8, DIMENSION(:,:), INTENT(inout) :: sumxx, sumxxsq
      REAL*8, DIMENSION(:,:), INTENT(out)   :: xxm, xxdev
      LOGICAL,                INTENT(in)    :: compute

      INCLUDE 'comatm.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comdiag.h'

      IF (itel /= 0) THEN
         IF (itel == 1) THEN
            sumxx  (:,:) = xx(:,:)
            sumxxsq(:,:) = xx(:,:)**2
         ELSE
            sumxx  (:,:) = sumxx  (:,:) + xx(:,:)
            sumxxsq(:,:) = sumxxsq(:,:) + xx(:,:)**2
         END IF

         IF (compute) THEN
            xxm  (:,:) = sumxx  (:,:)/minterv
            xxdev(:,:) = sumxxsq(:,:) - sumxx(:,:)**2/minterv
            WHERE (xxdev(:,:) < 0.0)
               xxdev(:,:) = 0.0
            ELSE WHERE
               xxdev (:,:) = Sqrt(xxdev(:,:)/(minterv - 1.0))
            END WHERE
         END IF
      END IF

      RETURN
      END SUBROUTINE stat


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iniocbas
!-----------------------------------------------------------------------

! *** initialises the oceanic bassin : based on iniocout
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comdiag.h'
      include 'comcoup.h'
      include 'comunit.h'

      integer     i,j,k,ias,ibas
      character*1 ch(nlon),space
      real*8 sum(nlat,nbasa)

!     (1)='ATL N  '
!     (2)='PAC N  '
!     (3)='ARCTIC '
!     (4)='INDIAN '
!     (5)='ATL S  '
!     (6)='PAC S  '
!     (7)='ANTAR  '


! *** computation of basin masks

! *** asci number of small letter a

      ias=ichar('a') - 1
      do j=nlat,1,-1
        read (iuo+14,100) k,space,(ch(i),i=1,nlon)
        do i=1,nlon
          iocbasa(j,i)=ichar(ch(i)) - ias
          if (iocbasa(j,i).lt.1.or.iocbasa(j,i).gt.nbasa) then
            if (fractoc(j,i).gt.0.01) then
              write(iuo+29,*) 'no ocean basin defined in',i,j
            endif
            iocbasa(j,i)=0
          endif
        enddo
      enddo

! *** computation of area of ocean basins
      do ibas=1,nbasa
        arocbasa(ibas)=0.
        do i=1,nlat
          sum(i,ibas)=0d0
        enddo
      enddo

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0) then
          sum(i,iocbasa(i,j))=sum(i,iocbasa(i,j))+fractoc(i,j)
         endif
         sum(i,nbasa)=sum(i,nbasa)+fractoc(i,j)
        enddo
      enddo

      do i=1,nlat
       do ibas=1,nbasa
        arocbasa(ibas)=arocbasa(ibas)+sum(i,ibas)*darea(i)
       enddo
      enddo

 100  format(i4,65A1)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine outamean(istep)
!-----------------------------------------------------------------------

! *** computes mean values of ocean basins: m/yr
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comdiag.h'
      include 'comemic.h'
      include 'comcoup.h'
      include 'comunit.h'

      integer i,j,k,l,istep,ibas
      real*8  precmean(nbasa),evapmean(nbasa),runomean(nbasa)
      real*8 sum1(nlat,nbasa),sum2(nlat,nbasa),sum3(nlat,nbasa)
      real*8 ttest,unit

      unit=86400.

      if (istep.eq.1) then
        do ibas=1,nbasa
          precan(ibas)=0.0
          evapan(ibas)=0.0
          runoan(ibas)=0.0
        enddo
      endif

      do l=1,nbasa
         precmean(l)=0.
         evapmean(l)=0.
         runomean(l)=0.
      enddo

! *** calculate mean values of all basins
      do ibas=1,nbasa
        do i=1,nlat
          sum1(i,ibas)=0d0
          sum2(i,ibas)=0d0
          sum3(i,ibas)=0d0
        enddo
      enddo

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0) then
          sum1(i,iocbasa(i,j))=sum1(i,iocbasa(i,j))+fractoc(i,j) &
     &         *(torain(i,j)+tosnow(i,j))
          sum2(i,iocbasa(i,j))=sum2(i,iocbasa(i,j))+fractoc(i,j) &
     &         *evap(i,j)
          sum3(i,iocbasa(i,j))=sum3(i,iocbasa(i,j))+fractoc(i,j) &
     &         *arunofo(i,j)
         endif

         sum1(i,nbasa)=sum1(i,nbasa)+fractoc(i,j)*(torain(i,j)+tosnow(i,j))
         sum2(i,nbasa)=sum2(i,nbasa)+fractoc(i,j)*evap(i,j)
         sum3(i,nbasa)=sum3(i,nbasa)+fractoc(i,j)*arunofo(i,j)
        enddo
      enddo

      do i=1,nlat
       do ibas=1,nbasa
        precmean(ibas)=precmean(ibas)+sum1(i,ibas)*darea(i)
        evapmean(ibas)=evapmean(ibas)+sum2(i,ibas)*darea(i)
        runomean(ibas)=runomean(ibas)+sum3(i,ibas)*darea(i)
       enddo
      enddo
      do ibas=1,nbasa
        precan(ibas)=precan(ibas) + precmean(ibas)
        evapan(ibas)=evapan(ibas) + evapmean(ibas)
        runoan(ibas)=runoan(ibas) + runomean(ibas)
      enddo

      ttest=float(istep)/float(nstpyear)

      if (ttest.eq.int(ttest)) then
        do ibas=1,nbasa
          precan(ibas)=precan(ibas)*unit/(arocbasa(ibas)*float(iatm))
          evapan(ibas)=evapan(ibas)*unit/(arocbasa(ibas)*float(iatm))
          runoan(ibas)=runoan(ibas)*unit/(arocbasa(ibas)*float(iatm))
        enddo
        write(iuo+26,*) 'precip'
        write(iuo+26,*)  precan
        write(iuo+26,*) 'evap'
        write(iuo+26,*)  evapan
        write(iuo+26,*) 'runoff'
        write(iuo+26,*)  runoan
        write(iuo+26,*)
        call flush(iuo+26)
        do ibas=1,nbasa
          precan(ibas)=0.0
          evapan(ibas)=0.0
          runoan(ibas)=0.0
        enddo
      endif
      end

!23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine outiocht(istep)
!-----------------------------------------------------------------------

! *** computes annual mean implied ocean heat transport: PW (10^15 W)
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comdiag.h'
      include 'comemic.h'
      include 'comcoup.h'
      include 'comunit.h'

      integer i,j,k,l,istep,ibas
      real*8  hfmean(nlat,nbasa)
      real*8  hatl(nlat),hpac(nlat),hind(nlat),hall(nlat)
      real*8  ttest,unit,dum

      unit=1e-15
      if (istep.eq.1) then
        do ibas=1,nbasa
          do i=1,nlat
            hfmeanan(i,ibas)=0d0
          enddo
        enddo
      endif
      do ibas=1,nbasa
        do i=1,nlat
          hfmean(i,ibas)=0d0
        enddo
      enddo

! *** calculate zonal mean net surface heatflux of all basins

      do j=1,nlon
        do i=1,nlat
          dum=fractoc(i,j)*(heswsn(i,j,noc)+dlradsn(i,j,noc)- &
     &        ulradsn(i,j,noc)-efluxn(i,j,noc)-hfluxn(i,j,noc))
          if (iocbasa(i,j).gt.0) then
            hfmean(i,iocbasa(i,j))=hfmean(i,iocbasa(i,j))+dum
          endif
          hfmean(i,nbasa)=hfmean(i,nbasa)+dum
        enddo
      enddo

      do ibas=1,nbasa
        do i=1,nlat
          hfmeanan(i,ibas)=hfmeanan(i,ibas) + hfmean(i,ibas)*darea(i)
        enddo
      enddo

      if(mod(nint(day*real(iatm))+1,nstpyear).eq.0) then
!if (mod(istep,nstpyear).eq.0) then
        do ibas=1,nbasa
          do i=1,nlat
            hfmeanan(i,ibas)=hfmeanan(i,ibas)/real(nstpyear)
          enddo
        enddo
        do i=1,nlat
          hatl(i)=0d0
          hpac(i)=0d0
          hind(i)=0d0
          hall(i)=0d0
        enddo
        hatl(nlat)=-hfmeanan(nlat,1)
        hpac(nlat)=-hfmeanan(nlat,2)
        hind(nlat)=-hfmeanan(nlat,4)
        hall(nlat)=-hfmeanan(nlat,nbasa)
        do i=nlat-1,1,-1
          hatl(i)=hatl(i+1)-hfmeanan(i,1)-hfmeanan(i,5)
          hpac(i)=hpac(i+1)-hfmeanan(i,2)-hfmeanan(i,6)
          hind(i)=hind(i+1)-hfmeanan(i,4)
          hall(i)=hall(i+1)-hfmeanan(i,nbasa)
        enddo
        write(iuo+27) (real(hatl(i)*unit),i=1,nlat)
        write(iuo+27) (real(hpac(i)*unit),i=1,nlat)
        write(iuo+27) (real(hind(i)*unit),i=1,nlat)
        write(iuo+27) (real(hall(i)*unit),i=1,nlat)
        call flush(iuo+27)
        do ibas=1,nbasa
          do i=1,nlat
            hfmeanan(i,ibas)=0d0
          enddo
        enddo
      endif
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE outputyrl
!-----------------------------------------------------------------------
! *** this routine calls another routine which computes the yearly mean
! *** written by camiel severijns, march 2001.
!-----------------------------------------------------------------------
      USE Atmosphere_Output
      IMPLICIT NONE

      INCLUDE 'comatm.h'
      INCLUDE 'comdyn.h'
      INCLUDE 'comphys.h'
      INCLUDE 'comsurf.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comcoup.h'
      INCLUDE 'comdiag.h'
      INCLUDE 'comoutlocal.h'

      REAL*8, DIMENSION(1:nlat,1:nlon), SAVE :: &
     &     tsurfs, tsurfs2, tstrats, tstrats2, &
     &     temp2gs, temp2gs2, temp4gs, temp4gs2, &
     &     tempsgs, tempsgs2, tstrat2s, tstrat2s2, &
     &     t2ms, t2ms2, u200s, u200s2, u500s, u500s2, &
     &     u800s, u800s2, v200s, v200s2, v500s, v500s2, &
     &     v800s, v800s2, pgrounds, pgrounds2, omeg1s, &
     &     omeg1s2, omeg2s, omeg2s2, omeg3s, omeg3s2, &
     &     winstu1s, winstu1s2, winstv1s, winstv1s2, &
     &     uv10s, uv10s2, udivg1s, udivg1s2, udivg2s, &
     &     udivg2s2, udivg3s, udivg3s2, vdivg1s, vdivg1s2, &
     &     vdivg2s, vdivg2s2, vdivg3s, vdivg3s2, grpsi1s, &
     &     grpsi1s2, grpsi2s, grpsi2s2, grpsi3s, grpsi3s2, &
     &     chi1s, chi1s2, chi2s, chi2s2, chi3s, chi3s2, &
     &     qgpv1s, qgpv1s2, qgpv2s, qgpv2s2, qgpv3s, &
     &     qgpv3s2, gh1s, gh1s2, gh2s, gh2s2, gh3s, gh3s2, &
     &     vhforg1s, vhforg1s2, vhforg2s, vhforg2s2
      REAL*8, DIMENSION(1:nlat,1:nlon), SAVE :: &
     &     vforg1s, vforg1s2, vforg2s, vforg2s2, vforg3s, &
     &     vforg3s2, dyrains, dyrains2, corains, corains2, &
     &     torains, torains2, snows, snows2, evaps, evaps2, &
     &     eminps, eminps2, hfluxs, hfluxs2, efluxs, &
     &     efluxs2, heswss, heswss2, hesw_s, hesw_s2, &
     &     albess, albess2, albeps, albeps2, nlradss, &
     &     nlradss2, ulrad1s, ulrad1s2, bmoisgs, bmoisgs2, &
     &     dsnows, dsnows2, hics, hics2, runofls, runofls2, &
     &     runofos, runofos2, rmoisgw3s, rmoisgw3s2, &
     &     relhums, relhums2, cdragws, cdragws2, cdragvs, &
     &     cdragvs2, tccs, tccs2, dt11s, dt11s2, dt12s, &
     &     dt12s2, dt13s, dt13s2, dt21s, dt21s2, dt22s, &
     &     dt22s2, dt23s, dt23s2, du11s, du11s2, du12s, &
     &     du12s2, du13s, du13s2, du21s, du21s2, du22s, &
     &     du22s2, du23s, du23s2

      INTERFACE
         SUBROUTINE yrlstat(xx,sumxx,sumxxsq,xxm,xxdev,compute,samples)
         IMPLICIT NONE
         REAL*8, DIMENSION(:,:), INTENT(in)    :: xx
         REAL*8, DIMENSION(:,:), INTENT(inout) :: sumxx, sumxxsq
         REAL*8, DIMENSION(:,:), INTENT(out)   :: xxm, xxdev
         LOGICAL,                INTENT(in)    :: compute
         INTEGER,                INTENT(in)    :: samples
         END SUBROUTINE yrlstat
      END INTERFACE

      INTEGER, PARAMETER :: IOFlag = 4

      LOGICAL :: need_to_write
      REAL*8, DIMENSION(1:nlat,1:nlon,1:4) :: mean, stddev
      INTEGER, SAVE :: samples = 0
      INTEGER :: realmonth

      realmonth=12
      if ( irunlabeld/=360 ) realmonth=(irunlabeld/30)

      samples = samples + 1 ! Increment sampel counter

      need_to_write = ( imonth == realmonth .AND. iday == 30 )

      IF (need_to_write) CALL open(Yearly_Means)

      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
         CALL yrlstat(tsurf1(1:nlat,1:nlon), &
     &        tsurfs(1:nlat,1:nlon),   tsurfs2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1), stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         CALL yrlstat(temp0g1(1:nlat,1:nlon), &
     &        tstrats(1:nlat,1:nlon),  tstrats2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(temp2g1(1:nlat,1:nlon), &
     &        temp2gs(1:nlat,1:nlon),  temp2gs2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(temp4g1(1:nlat,1:nlon), &
     &        temp4gs(1:nlat,1:nlon),  temp4gs2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         CALL yrlstat(tempsg1(1:nlat,1:nlon), &
     &        tempsgs(1:nlat,1:nlon),  tempsgs2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,4),stddev  (1:nlat,1:nlon,4), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
         END IF
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         CALL yrlstat(temp0g1(1:nlat,1:nlon), &
     &        tstrat2s(1:nlat,1:nlon),  tstrat2s2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
           CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         CALL yrlstat(tempsg1(1:nlat,1:nlon), &
     &        t2ms(1:nlat,1:nlon),  t2ms2 (1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         CALL yrlstat(u200(1:nlat,1:nlon), &
     &        u200s(1:nlat,1:nlon),  u200s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(u500(1:nlat,1:nlon), &
     &        u500s(1:nlat,1:nlon),  u500s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(u800(1:nlat,1:nlon), &
     &        u800s(1:nlat,1:nlon),  u800s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         CALL yrlstat(v200(1:nlat,1:nlon), &
     &        v200s(1:nlat,1:nlon),  v200s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(v500(1:nlat,1:nlon), &
     &        v500s(1:nlat,1:nlon),  v500s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(v800(1:nlat,1:nlon), &
     &        v800s(1:nlat,1:nlon),  v800s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         CALL yrlstat(pground(1:nlat,1:nlon), &
     &        pgrounds(1:nlat,1:nlon),  pgrounds2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         CALL yrlstat(omegg(1:nlat,1:nlon,1), &
     &        omeg1s(1:nlat,1:nlon),  omeg1s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(omegg(1:nlat,1:nlon,2), &
     &        omeg2s(1:nlat,1:nlon),  omeg2s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(omegg(1:nlat,1:nlon,3), &
     &        omeg3s(1:nlat,1:nlon),  omeg3s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         CALL yrlstat(winstu1(1:nlat,1:nlon), &
     &        winstu1s(1:nlat,1:nlon),  winstu1s2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         CALL yrlstat(winstv1(1:nlat,1:nlon), &
     &        winstv1s(1:nlat,1:nlon),  winstv1s2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         CALL yrlstat(uv10(1:nlat,1:nlon), &
     &        uv10s(1:nlat,1:nlon),  uv10s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         CALL yrlstat(udivg(1:nlat,1:nlon,1), &
     &        udivg1s(1:nlat,1:nlon),  udivg1s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(udivg(1:nlat,1:nlon,2), &
     &        udivg2s(1:nlat,1:nlon),  udivg2s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(udivg(1:nlat,1:nlon,3), &
     &        udivg3s(1:nlat,1:nlon),  udivg3s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         CALL yrlstat(vdivg(1:nlat,1:nlon,1), &
     &        vdivg1s(1:nlat,1:nlon),  vdivg1s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(vdivg(1:nlat,1:nlon,2), &
     &        vdivg2s(1:nlat,1:nlon),  vdivg2s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(vdivg(1:nlat,1:nlon,3), &
     &        vdivg3s(1:nlat,1:nlon),  vdivg3s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         CALL yrlstat(psig(1:nlat,1:nlon,1), &
     &        grpsi1s(1:nlat,1:nlon),  grpsi1s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(psig(1:nlat,1:nlon,2), &
     &        grpsi2s(1:nlat,1:nlon),  grpsi2s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(psig(1:nlat,1:nlon,3), &
     &        grpsi3s(1:nlat,1:nlon),  grpsi3s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         CALL yrlstat(chig(1:nlat,1:nlon,1), &
     &        chi1s(1:nlat,1:nlon),  chi1s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(chig(1:nlat,1:nlon,1), &
     &        chi2s(1:nlat,1:nlon),  chi2s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(chig(1:nlat,1:nlon,1), &
     &        chi3s(1:nlat,1:nlon),  chi3s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         CALL yrlstat(qgpv(1:nlat,1:nlon,1), &
     &        qgpv1s(1:nlat,1:nlon),  qgpv1s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(qgpv(1:nlat,1:nlon,2), &
     &        qgpv2s(1:nlat,1:nlon),  qgpv2s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(qgpv(1:nlat,1:nlon,3), &
     &        qgpv3s(1:nlat,1:nlon),  qgpv3s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         CALL yrlstat(geopg(1:nlat,1:nlon,1), &
     &        gh1s(1:nlat,1:nlon),  gh1s2 (1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(geopg(1:nlat,1:nlon,2), &
     &        gh2s(1:nlat,1:nlon),  gh2s2 (1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(geopg(1:nlat,1:nlon,3), &
     &        gh3s(1:nlat,1:nlon),  gh3s2 (1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         CALL yrlstat(vhforg1(1:nlat,1:nlon), &
     &        vhforg1s(1:nlat,1:nlon),  vhforg1s2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(vhforg2(1:nlat,1:nlon), &
     &        vhforg2s(1:nlat,1:nlon),  vhforg2s2(1:nlat,1:nlon), &
     &        mean    (1:nlat,1:nlon,2),stddev   (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
         END IF
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         CALL yrlstat(vforg1(1:nlat,1:nlon), &
     &        vforg1s(1:nlat,1:nlon),  vforg1s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(vforg2(1:nlat,1:nlon), &
     &        vforg2s(1:nlat,1:nlon),  vforg2s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(vforg3(1:nlat,1:nlon), &
     &        vforg3s(1:nlat,1:nlon),  vforg3s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         CALL yrlstat(dyrain1(1:nlat,1:nlon), &
     &        dyrains(1:nlat,1:nlon),  dyrains2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         CALL yrlstat(corain1(1:nlat,1:nlon), &
     &        corains(1:nlat,1:nlon),  corains2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         CALL yrlstat(torain1(1:nlat,1:nlon), &
     &        torains(1:nlat,1:nlon),   torains2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1), stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         CALL yrlstat(snow1(1:nlat,1:nlon), &
     &        snows(1:nlat,1:nlon),  snows2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         CALL yrlstat(evap1(1:nlat,1:nlon), &
     &        evaps(1:nlat,1:nlon),  evaps2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         CALL yrlstat(eminp1(1:nlat,1:nlon), &
     &        eminps(1:nlat,1:nlon),  eminps2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         CALL yrlstat(hflux(1:nlat,1:nlon), &
     &        hfluxs(1:nlat,1:nlon),  hfluxs2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         CALL yrlstat(eflux(1:nlat,1:nlon), &
     &        efluxs(1:nlat,1:nlon),  efluxs2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         CALL yrlstat(hesws(1:nlat,1:nlon), &
     &        heswss(1:nlat,1:nlon),  heswss2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         CALL yrlstat(hesw(1:nlat,1:nlon), &
     &        hesw_s(1:nlat,1:nlon),  hesw_s2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         CALL yrlstat(alb2es(1:nlat,1:nlon), &
     &        albess(1:nlat,1:nlon),  albess2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         CALL yrlstat(albep(1:nlat,1:nlon), &
     &        albeps(1:nlat,1:nlon),  albeps2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         CALL yrlstat(nlrads(1:nlat,1:nlon), &
     &        nlradss(1:nlat,1:nlon),  nlradss2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         CALL yrlstat(ulrad0(1:nlat,1:nlon), &
     &        ulrad1s(1:nlat,1:nlon),  ulrad1s2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
         CALL yrlstat(abmoisg(1:nlat,1:nlon), &
     &        bmoisgs(1:nlat,1:nlon),  bmoisgs2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
         CALL yrlstat(adsnow(1:nlat,1:nlon), &
     &        dsnows(1:nlat,1:nlon),  dsnows2(1:nlat,1:nlon), &
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         CALL yrlstat(ahic(1:nlat,1:nlon), &
     &        hics(1:nlat,1:nlon),  hics2 (1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         CALL yrlstat(runofl1(1:nlat,1:nlon), &
     &        runofls(1:nlat,1:nlon),  runofls2(1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         CALL yrlstat(runofo1(1:nlat,1:nlon), &
     &        runofos(1:nlat,1:nlon),  runofos2(1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
         CALL yrlstat(rmoisg(1:nlat,1:nlon), &
     &        rmoisgw3s(1:nlat,1:nlon),  rmoisgw3s2(1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev    (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         CALL yrlstat(relhum(1:nlat,1:nlon), &
     &        relhums(1:nlat,1:nlon),  relhums2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         CALL yrlstat(cdragw(1:nlat,1:nlon), &
     &        cdragws(1:nlat,1:nlon),  cdragws2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         CALL yrlstat(cdragv(1:nlat,1:nlon), &
     &        cdragvs(1:nlat,1:nlon),  cdragvs2(1:nlat,1:nlon), &
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         CALL yrlstat(tccd(1:nlat,1:nlon), &
     &        tccs(1:nlat,1:nlon),  tccs2 (1:nlat,1:nlon), &
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         CALL yrlstat(dumt1(1:nlat,1:nlon,1), &
     &        dt11s(1:nlat,1:nlon),  dt11s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(dumt1(1:nlat,1:nlon,2), &
     &        dt12s(1:nlat,1:nlon),  dt12s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(dumt1(1:nlat,1:nlon,3), &
     &        dt13s(1:nlat,1:nlon),  dt13s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         CALL yrlstat(dumt2(1:nlat,1:nlon,1), &
     &        dt21s(1:nlat,1:nlon),  dt21s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(dumt2(1:nlat,1:nlon,2), &
     &        dt22s(1:nlat,1:nlon),  dt22s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(dumt2(1:nlat,1:nlon,3), &
     &        dt23s(1:nlat,1:nlon),  dt23s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         CALL yrlstat(dumu1(1:nlat,1:nlon,1), &
     &        du11s(1:nlat,1:nlon),  du11s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(dumu1(1:nlat,1:nlon,2), &
     &        du12s(1:nlat,1:nlon),  du12s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(dumu1(1:nlat,1:nlon,3), &
     &        du13s(1:nlat,1:nlon),  du13s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         CALL yrlstat(dumu2(1:nlat,1:nlon,1), &
     &        du21s(1:nlat,1:nlon),  du21s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1), &
     &        need_to_write,samples)
         CALL yrlstat(dumu2(1:nlat,1:nlon,2), &
     &        du22s(1:nlat,1:nlon),  du22s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2), &
     &        need_to_write,samples)
         CALL yrlstat(dumu2(1:nlat,1:nlon,3), &
     &        du23s(1:nlat,1:nlon),  du23s2(1:nlat,1:nlon), &
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3), &
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (need_to_write) THEN
         CALL close
         samples = 0 ! Restart counting samples
      END IF

      RETURN
      END SUBROUTINE outputyrl

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE yrlstat(x,sumx1,sumy1,xmean,xstd,compute,samples)
! *** --------------------------------------------------------------------
! *** This routine computes yearly mean and standard  deviation around it.
! *** all arrays are assumed to have the same shape.
! *** ------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:), INTENT(in)    :: x
      REAL*8, DIMENSION(:,:), INTENT(inout) :: sumx1, sumy1
      REAL*8, DIMENSION(:,:), INTENT(out)   :: xmean, xstd
      LOGICAL,                INTENT(in)    :: compute
      INTEGER,                INTENT(in)    :: samples

      INCLUDE 'comatm.h'
      INCLUDE 'comcoup.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comdiag.h'

      IF (imonth == 1 .AND. iday == 1) THEN
         sumx1(:,:) = 0.0
         sumy1(:,:) = 0.0
      END IF

      sumx1(:,:) = sumx1(:,:) + x(:,:)
      sumy1(:,:) = sumy1(:,:) + x(:,:)**2

      IF (compute) THEN
         xmean(:,:) = sumx1(:,:)/samples
         xstd (:,:) = sumy1(:,:) - sumx1(:,:)**2/samples
         WHERE (xstd(:,:) <= 0.0)
            xstd(:,:) = 0.0
         ELSE WHERE
            xstd(:,:) = Sqrt(xstd(:,:)/(samples - 1.0))
         END WHERE
      END IF

      RETURN
      END SUBROUTINE yrlstat
