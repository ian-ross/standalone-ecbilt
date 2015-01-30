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
                write(100,*) 'sst out of range ',i,j,sstday(i,j,id)
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
               write(29,*) 'birth icemask latlon ',i,j
               call error(17)
            endif
            if (lsicebirth(i,j) .ne. iap .and. fractl(i, j) .gt. epsl) then
               write(29,*) 'birth icemask latlon ',i,j
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
            write(29,*) 'death icemask latlon ',i,j
            call error(17)
          endif
          if (lsicedeath(i,j).ne.iap.and. fractl(i, j) .gt. epsl) then
            write(29,*) 'death icemask latlon ',i,j
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
              write(29,*) 'icemask latlon ',i,j
              call error(17)
            endif
          endif
          if (lsicebirth(i,j).gt.0) then
            if (lsicedeath(i,j).eq.0.or. &
                 & lsicebirth(i,j).eq.lsicedeath(i,j)) then
              write(29,*) 'icemask latlon ',i,j
              call error(17)
            endif
          endif
        enddo
      enddo

100   format(i4,65A1)
310   format(i4,i2,90i1)

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
      real*8 dum, facwin, facsum
      real*8 seaalb(nlat, nlon)
      real*8 zalb,zalbp

!      write(100,*) 'index in oceanfixed ',index

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
      subroutine seaicetemp
!-----------------------------------------------------------------------
! *** calculates seaice temperatures assuming that all surface fluxes
! *** balance
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comice.h'
      include 'comocfix.h'

      integer  i,j,itetel,mxtetel,il,jl
      real*8   stice,stice1,stice2,tol,fluxsumice,zbrent

      external fluxsumice

! *** tol is the wanted accuracy of the seaice temperature in degrees

      parameter (tol=0.1)

      common /landpoint/il,jl

      mxtetel = 0

      do j=1,nlon
        do i=1,nlat

          if (lseaice(i,j).eq.1) then

            il=i
            jl=j
            stice=tzero
            stice1=stice - 1.
            stice2=stice

            call zbrac(fluxsumice,stice1,stice2,itetel)
            if (itetel.eq.100) call error(7)
            tijs(i,j)=zbrent(fluxsumice,stice1,stice2,tol,itetel)
            if (itetel.eq.100) call error(8)
            if (itetel.gt.mxtetel) mxtetel=itetel

! *** in case of temperatures above zero, set
! *** surface temperature to meltpoint

            if (tijs(i,j).gt.330.or.tijs(i,j).lt.200) then
              write(100,*) 'ice temperature out of range ',i,j
              write(100,*) lsicebirth(i,j),lsicedeath(i,j),lseaice(i,j)
              write(100,*) uv10(i,j),tempsg(i,j),q10(i,j)
              write(100,*) dlrads(i,j),hesws(i,j)
            endif

            if (tijs(i,j).gt.tzero) then
	      tijs(i,j)=tzero
            endif

          endif
        enddo
      enddo


!      write(100,*) 'in oceanfixed over ice ',mxtetel
      call flush(100)

      return
      end

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
      subroutine initseaalb

      implicit none

      include 'comcouphelp.h'
      include 'comunit.h'

      integer i,is
      real*8  albsea(nlat,4)
      common /albedoclio/albsea

! *** read climatological zonal mean albedos for each season

      open (iuo+49,file='inputdata/albedo.dat')

      read(iuo+49,*)
      do i=1,nlat
	read(iuo+49,45)(albsea(i,is),is=1,4)
      enddo
45    format(4(2x,f7.4))
      close(iuo+49)
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      function fluxsumice(stice)
!-----------------------------------------------------------------------
! *** computes sum of fluxes between the seaice and the atmosphere
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'

      integer il,jl
      real*8  fluxsumice,stice,qsatss,qsat,efluxgp,hfluxgp

      common /landpoint/il,jl

! *** sensible heatflux

      hfluxgp=alphad*cdragv(il,jl)*uv10(il,jl)*(stice-tempsg(il,jl))

! *** latent heat flux

      qsatss=qsat(1.d+5,stice)

      if (stice.gt.tzero) then
        efluxgp=alphav*cdragv(il,jl)*uv10(il,jl)*(qsatss-q10(il,jl))
      else
        efluxgp=alphas*cdragv(il,jl)*uv10(il,jl)*(qsatss-q10(il,jl))
      endif

      if (efluxgp.lt.0) efluxgp=0.

! *** sum of all fluxes

      fluxsumice=hesws(il,jl) - hfluxgp + dlrads(il,jl) - &
     &            sboltz*stice**4 - efluxgp


      return
      end

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

      include 'comatm.h'
      include 'comphys.h'
      integer ih
      real*8 tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp
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

          if (hnbq.gt.0.05) then
            zalbp = alphd
          else
            if (hgbq.gt.1.5) then
              zalbp = alphdi+(hnbq*(alphd-alphdi)/0.05)
            else if (hgbq.gt.1.0.and.hgbq.le.1.5) then
                   al = 0.472+2.0*(alphdi-0.472)*(hgbq-1.0)
            else if (hgbq.gt.0.05.and.hgbq.le.1.0) then
                   al = 0.2467+(0.7049*hgbq)-(0.8608*(hgbq*hgbq))+ &
     &                 (0.3812*(hgbq*hgbq*hgbq))
            else
              al = 0.1+3.6*hgbq
            endif
            if (hgbq.le.1.5) zalbp=al+(hnbq*(alphd-al)/0.05)
          endif
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
     &                   (0.7049*hgbq)-(0.8608*(hgbq*hgbq))+ &
     &                   (0.3812*(hgbq*hgbq*hgbq))
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
