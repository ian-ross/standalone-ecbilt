!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inioceanfixed
!-----------------------------------------------------------------------
! *** initialises and sets parameters of the fixed ocean model
!-----------------------------------------------------------------------
      implicit none

      include 'comunit.h'
      include 'comland.h'
      include 'comemic.h'
      include 'comocfix.h'

      integer      i,j,k,id,ia1,ia9,ica,icd,iap
      real*4       sstday4(nlat,nlon)
      character*1  ch(nlon),space

!
! *** read sst on daily basis
!

      do id=1,360
        read(iuo+41) (( sstday4(i,j),j=1,nlon),i=1,nlat)
        do i=1,nlat
          do j=1,nlon
            sstday(i,j,id)=sstday4(i,j)
            if (fractl(i, j) .lt. epsl) then
              if (sstday(i, j, id) .ge. 310.or. &
                   & sstday(i, j, id) .lt. tzero - 1.) then
                 write(iuo+29,*) 'sst out of range ',i,j,sstday(i,j,id)
                 if (sstday(i, j, id) >= 310) sstday(i, j, id) = 310
                 if (sstday(i, j, id) < tzero - 1.0) &
                      & sstday(i, j, id) = tzero - 1.0
              endif
            endif
          enddo
        enddo
      enddo

! *** read climatological seaice cover

      ia1=ichar('1')
      ia9=ichar('9')
      ica=ichar('A')
      icd=ichar('D')
      iap=ichar('.')

      do i=nlat,1,-1
         read (iuo+43,100) k,space,(ch(j),j=1,nlon)
         do j=1,nlon
            lsicebirth(i,j)=ichar(ch(j))
            if (lsicebirth(i,j) .eq. iap .and. fractl(i, j) .le. epsl) then
               write(iuo+29,*) 'birth icemask latlon ',i,j,ch(j), fractl(i,j)
               call error(17)
            endif
            if (lsicebirth(i,j) .ne. iap .and. fractl(i, j) .gt. epsl) then
               write(iuo+29,*) 'birth icemask latlon ',i,j,ch(j), fractl(i,j)
               call error(17)
            endif
            if (lsicebirth(i, j) .ge. ia1 .and. lsicebirth(i, j) .le. ia9) then
               lsicebirth(i,j)=lsicebirth(i,j)-ia1+1
            else
               if (lsicebirth(i,j) .ge. ica .and. lsicebirth(i,j) .le. icd) then
                  lsicebirth(i,j)=lsicebirth(i,j)-ica+10
               else
                  lsicebirth(i,j)=0
               endif
            endif
         enddo
      enddo

      do i=nlat,1,-1
        read (iuo+43,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          lsicedeath(i,j)=ichar(ch(j))
          if (lsicedeath(i,j).eq.iap.and. fractl(i, j) .le. epsl) then
            write(iuo+29,*) 'death icemask latlon ',i,j
            call error(17)
          endif
          if (lsicedeath(i,j).ne.iap.and. fractl(i, j) .gt. epsl) then
            write(iuo+29,*) 'death icemask latlon ',i,j
            call error(17)
          endif
          if (lsicedeath(i,j).ge.ia1.and.lsicedeath(i,j).le.ia9) then
            lsicedeath(i,j)=lsicedeath(i,j)-ia1+1
          else
            if (lsicedeath(i,j).ge.ica.and.lsicedeath(i,j).le.icd) then
              lsicedeath(i,j)=lsicedeath(i,j)-ica+10
            else
              lsicedeath(i,j)=0
            endif
          endif
          if (lsicedeath(i,j).gt.0) then
            if (lsicebirth(i,j).eq.0.or. &
                 & lsicebirth(i,j).eq.lsicedeath(i,j)) then
              write(iuo+29,*) 'icemask latlon ',i,j
              call error(17)
            endif
          endif
          if (lsicebirth(i,j).gt.0) then
            if (lsicedeath(i,j).eq.0.or. &
                 & lsicebirth(i,j).eq.lsicedeath(i,j)) then
              write(iuo+29,*) 'icemask latlon ',i,j
              call error(17)
            endif
          endif
        enddo
      enddo

100   format(i4,65A1)

      call initseaalb

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine oceanfixed(index)
!-----------------------------------------------------------------------
! *** prescribe sst's and seaice cover over ocean and lakes
! *** calculate seaice temperature assuming all surface fluxes balance
!-----------------------------------------------------------------------
      implicit none

      include 'comland.h'
      include 'comemic.h'
      include 'comocfix.h'
      include 'comcoup.h'
      include 'comsurf.h'
      include 'comice.h'

      integer index,i,j
      real*8 facwin, facsum
      real*8 seaalb(nlat, nlon)
      real*8 zalb,zalbp

! *** prescription of seaice cover

      facwin = dabs(180.d0 - day) / 180.d0
      facsum = 1d0 - facwin

      do i = 1, nlat
        do j = 1, nlon
          lseaice(i, j) = 0
          if (fractl(i, j) .le. epsl) then
            if (lsicebirth(i, j) .gt. 0) then
              if (lsicebirth(i, j) .lt. lsicedeath(i, j)) then
                if (imonth .ge. lsicebirth(i, j).and. &
                     & imonth .lt. lsicedeath(i, j)) lseaice(i, j) = 1
              endif
              if (lsicebirth(i, j) .gt. lsicedeath(i, j)) then
                if (imonth .ge. lsicebirth(i, j).or. &
                     & imonth .lt. lsicedeath(i, j)) lseaice(i, j) = 1
              endif
            endif
          endif
        enddo
      enddo

! *** calculate seaice temperatures

      call seaicetemp

! *** prescription of sst

      call detseaalb(seaalb)

      do i = 1, nlat
         do j = 1, nlon
            fractn(i, j, noc) = (1.0d0 - lseaice(i, j)) * fracto(i, j)
            fractn(i, j, nse) = lseaice(i, j) * fracto(i, j)
            tsurfn(i, j, nse) = min(tzero, tijs(i, j))
            tsurfn(i, j, noc) = max(tzero - 1.8d0, sstday(i, j, index))
            call shine(tzero - 0.15, tzero - 0.25, tsurfn(i, j, nse), &
     &           hic(i, j), hsn(i, j), zalb, zalbp)
            albesn(i, j, nse) = (1.0 - couptcc(i, j)) * zalbp + &
     &           couptcc(i, j) * zalb
            albesn(i, j, noc) = albocef * seaalb(i, j)
         end do
      end do

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** calculates seaice temperatures assuming that all surface fluxes
! *** balance
!-----------------------------------------------------------------------
      SUBROUTINE seaicetemp
      IMPLICIT NONE
      INCLUDE 'comunit.h'
      INCLUDE 'comatm.h'
      INCLUDE 'comdyn.h'
      INCLUDE 'comphys.h'
      INCLUDE 'comsurf.h'
      INCLUDE 'comice.h'
      INCLUDE 'comocfix.h'
      INCLUDE 'comemic.h'

      INTEGER i, j, itetel, mxtetel, il, jl
      REAL*8 stice, stice1, stice2, fluxsumice, zbrent
! *** tol is the wanted accuracy of the seaice temperature in degrees
      REAL*8, PARAMETER :: tol = 0.1
      EXTERNAL fluxsumice
      COMMON /landpoint/ il, jl

      mxtetel = 0
      DO j = 1, nlon
         DO i = 1, nlat
            IF (lseaice(i,j) == 1) THEN
               il = i
               jl = j
               stice = tzero
               stice1 = stice - 1.0
               stice2 = stice
               CALL zbrac(fluxsumice, stice1, stice2, itetel)
               IF (itetel == 100) CALL error(7)
               tijs(i,j) = zbrent(fluxsumice, stice1, stice2, tol, itetel)
               IF (itetel == 100) CALL error(8)
               IF (itetel > mxtetel) mxtetel = itetel

! *** in case of temperatures above zero, set
! *** surface temperature to meltpoint
               IF (tijs(i,j) > 330 .OR. tijs(i,j) < 200) THEN
                  WRITE(iuo+29,*) 'ice temperature out of range ', &
                       & iyear,imonth,iday,i,j,tijs(i,j)
                  WRITE(iuo+29,*) lsicebirth(i,j),lsicedeath(i,j),lseaice(i,j)
                  WRITE(iuo+29,*) uv10(i,j),tempsg(i,j),q10(i,j)
                  WRITE(iuo+29,*) dlrads(i,j),hesws(i,j)
               END IF

               IF (tijs(i,j) < 200.0) tijs(i,j) = 200.0
               IF (tijs(i,j) > tzero) tijs(i,j) = tzero
            END IF
         END DO
      END DO

      CALL flush(iuo+29)

      RETURN
      END SUBROUTINE seaicetemp

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine detseaalb(seaalb)
!-----------------------------------------------------------------------
! *** calculates albedos as a function of the time of the year, linearly
! *** interpolating between seasonal mean values
! *** albsea is the albedo of the open sea
!-----------------------------------------------------------------------

      implicit none

      include 'comcouphelp.h'
      include 'comemic.h'

      integer i,j,id1,is1,is2
      real*8  albseaz(nlat),albsea(nlat,4),seaalb(nlat,nlon)
      real*8  sfrac

      common /albedoclio/albsea

! *** interpolate between seasonal means

      id1=(imonth-1)*30+iday-14
      if (id1.lt.1) id1=id1+360

      is1=(id1+89)/90
      is2=is1+1
      if (is2.eq.5) is2=1

      sfrac=(id1-((is1-1)*90.+1.))/90.

      do j=1,nlat
        albseaz(j)=albsea(j,is1)+(albsea(j,is2)-albsea(j,is1))*sfrac
      enddo

      do j=1,nlon
        do i=1,nlat
          seaalb(i,j)=albseaz(i)
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
    SUBROUTINE initseaalb
      IMPLICIT NONE
      INCLUDE 'comcouphelp.h'
      INCLUDE 'comunit.h'

      INTEGER i, is
      REAL*8  albsea(nlat, 4)
      COMMON /albedoclio/ albsea

! *** read climatological zonal mean albedos for each season

      OPEN(iuo+49, file='inputdata/sea_albedo.dat')

      READ (iuo+49,*)
      DO i = 1, nlat
         READ (iuo+49,45) (albsea(i, is), is = 1, 4)
      END DO
45    FORMAT(4(2x,f7.4))
      CLOSE(iuo+49)
    END SUBROUTINE initseaalb

!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** computes sum of fluxes between the seaice and the atmosphere
!-----------------------------------------------------------------------
      FUNCTION fluxsumice(stice)
      IMPLICIT NONE
      INCLUDE 'comatm.h'
      INCLUDE 'comdyn.h'
      INCLUDE 'comphys.h'
      INCLUDE 'comsurf.h'

      INTEGER il, jl
      REAL*8 fluxsumice, stice, qsatss, qsat, efluxgp, hfluxgp
      COMMON /landpoint/ il, jl

! *** sensible heatflux
      hfluxgp = alphad * cdragv(il,jl) * uv10(il,jl) * (stice - tempsg(il,jl))

! *** latent heat flux
      qsatss = qsat(1.0D+5, stice)
      IF (stice > tzero) THEN
         efluxgp = alphav * cdragv(il,jl) * uv10(il,jl) * (qsatss - q10(il,jl))
      ELSE
         efluxgp = alphas * cdragv(il,jl) * uv10(il,jl) * (qsatss - q10(il,jl))
      END IF
      IF (efluxgp < 0.0) efluxgp = 0.0

! *** sum of all fluxes
      fluxsumice = hesws(il,jl) - hfluxgp + dlrads(il,jl) - &
           & sboltz * stice**4 - efluxgp

      RETURN
      END FUNCTION fluxsumice

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine shine(tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp)
!-----------------------------------------------------------------------
! *** This subroutine computes albedo of snow-sea ice following SHINE &
! *** HENDERSSON-SELLERS [1985]
!-----------------------------------------------------------------------
!
!  tfsn   : melting point temperature of snow (273.15 in CLIO)
!  tfsg   : melting point temperature of ice (273.05 in CLIO)
!  ts     : surface temperature
!  hgbq   : ice thickness
!  hnbq   : snow thickness
!  zalb   : ice/snow albedo for overcast sky
!  zalbp  : ice/snow albedo for clear sky
      IMPLICIT NONE
      include 'comatm.h'
      include 'comphys.h'
      real*8 tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp,al
!     real*8 albin,albis,albice,alphd,alphdi,alphs,cgren



!driess    albin = 0.45
!driess    albis = 0.45
!driess    albice = 0.45
!     albice = 0.53
!     alphd  = 0.80
!     alphdi = 0.72
!     alphs  = 0.65
!driess    alphd  = 0.72
!driess    alphdi = 0.64
!driess    alphs  = 0.55
!     albin = 0.43
!     albis = 0.43
!     albice = 0.43
!     alphd  = 0.70
!     alphdi = 0.62
!     alphs  = 0.53
!     cgren = 0.04

!  albin: Albedo of melting ice in the arctic.
!  albis: Albedo of melting ice in the antarctic (SHINE
!         & HENDERSSON-SELLERS, 1985).
!  albice: Albedo of melting ice.
!  alphd : Albedo of snow (thickness > 0.05m)
!  alphdi: Albedo of thick bare ice
!  alphs : Albedo of melting snow
!  cgren: Correction of the snow or ice albedo to take into account
!         effects of cloudiness (GRENFELL & PEROVICH, 1984)

      if (hnbq.gt.0.0) then

! ***  Case of ice covered by snow.

        if (ts.lt.tfsn) then

! ***    Freezing snow.

          IF (hnbq > 0.05) THEN
            zalbp = alphd
          ELSE
            IF (hgbq > 1.5) THEN
               zalbp = alphdi + (hnbq * (alphd - alphdi) / 0.05)
            ELSE
               IF (hgbq > 1.0 .AND. hgbq <= 1.5) THEN
                  al = 0.472 + 2.0 * (alphdi - 0.472) * (hgbq - 1.0)
               ELSE IF (hgbq > 0.05 .AND. hgbq <= 1.0) THEN
                  al = 0.2467 + (0.7049 * hgbq) - (0.8608 * (hgbq * hgbq)) + &
                       & (0.3812 * (hgbq * hgbq * hgbq))
               ELSE
                  al = 0.1 + 3.6 * hgbq
               END IF
               IF (hgbq <= 1.5) zalbp = al + (hnbq * (alphd - al) / 0.05)
            END IF
          ENDIF
        else
!
! ***    Melting snow.
!
          if (hnbq.ge.0.1) then
            zalbp = alphs
          else
            zalbp = albice+((alphs-albice)/0.1)*hnbq
          endif
        endif
      else
!
! *** Case of ice free of snow.
!
        if (ts.lt.tfsg) then
!
! ***    Freezing ice.
!
          if (hgbq.gt.1.5) then
            zalbp = alphdi
          else if (hgbq.gt.1..and.hgbq.le.1.5) then
             zalbp = 0.472+2.*(alphdi-0.472)*(hgbq-1.)
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then
             zalbp = 0.2467+ &
                  & (0.7049*hgbq)-(0.8608*(hgbq*hgbq))+ &
                  & (0.3812*(hgbq*hgbq*hgbq))
          else
            zalbp = 0.1+3.6*hgbq
          endif
        else
!
! *** Melting ice.
!
          if (hgbq.gt.1.5) then
            zalbp = albice
          else if (hgbq.gt.1..and.hgbq.le.1.5)  then
                 zalbp = 0.472+(2.*(albice-0.472)*(hgbq-1.))
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then
                 zalbp = 0.2467+0.7049*hgbq &
     &                  -(0.8608*(hgbq*hgbq)) &
     &                  +(0.3812*(hgbq*hgbq*hgbq))
          else
            zalbp = 0.1+3.6*hgbq
          endif
        endif
      endif
      zalb=zalbp+cgren

      return
      end
