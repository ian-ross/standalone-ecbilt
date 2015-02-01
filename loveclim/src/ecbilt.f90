!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ecbilt(istep)


!-----------------------------------------------------------------------
! *** this routine performs one timestep of ECBilt which is a
! *** QG3L atmospheric model with simple physical parameterizations
! ***
! *** 6 april 1999 KNMI, De Bilt
! ***
! *** joint project Hugues Goosse
! ***               Rein Haarsma
! ***               Theo Opsteegh
! ***               Thijs van Reenen
! ***               Michiel Schaeffer
! ***               Frank Selten
! ***               Xueli Wang
! ***               Nanne Weber
! ***
! *** Modified code : P. Mathiot (01/2012) Merge topo with atmdyn0.f
! ***                                      Allow interannual topo forcing
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'

      integer istep

! *** atmospheric physics (in file atmphys.f)

      if (mod(nint(day*real(iatm)),nstpyear).eq.0) call topo

      call atmout(istep)
      call checks(istep)

      if (iaphys.eq.1) then

        call atmphyszero
        call sensrad
        call moisture
        call convec
        call fortemp
        call meantemp
      endif

! *** atmospheric dynamics (in file atmdyn.f)

      if (iadyn.eq.1) then
        call forward
      endif

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine update(istep)
!--------------------------------------------------------------------------------
! ***
! *** This routine updates the model date, incoming solar radiation and
! *** the atmospheric state and does outputting and checking
! ***
!-------------------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'

      integer istep

      call mdldate(istep)
      call solar
      call atmstate
      call vortfor


      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine checks(istep)
!-----------------------------------------------------------------------
! *** this routine performs checks on the model formulation
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'

      integer istep

      if ( mod(istep,iatm) .eq. 0) then
        call testecbilt
      endif

!      call moischeck(istep)
!      call heatcheck(istep)

      call flush(99)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine testecbilt
!-----------------------------------------------------------------------
! *** testing if model variables are inside prescribed ranges
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'

      integer i,j
      real*8  tsurfmean,globalmean
      character*12 chts
      real*8  moc,tmc,tmc0,cland,thex
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex

      do j=1,nlon
        do i=1,nlat

          if (uv10(i,j).gt.60) then
            write(iuo+99,*) 'uv10 out of range'
            write(iuo+99,*) i,j,uv10(i,j)
          endif
          if (tsurf(i,j).gt.400.or.tsurf(i,j).lt.150) then
            write(iuo+99,*) 'tsurf out of range in test'
            write(iuo+99,*) i,j,tsurf(i,j)
          endif
          if (eflux(i,j).gt.2000.or.hflux(i,j).gt.2000) then
            write(iuo+99,*) 'surface flux out of range in test'
            write(iuo+99,*) i,j,eflux(i,j),hflux(i,j),tsurf(i,j)
          endif

        enddo
      enddo

      tsurfmean=globalmean(tsurf)-tzero

      write(chts,900) tsurfmean
 900  format(E12.5)
      if (chts(1:3).eq.'nan') call error(99)

      write(iuo+20,110) iyear,int((day+0.5*dt)/(iatm*dt))+1,tsurfmean
      call flush(iuo+20)

      if (tsurfmean.gt.40..or.tsurfmean.lt.-10.) then
        write(iuo+29,*) 'mean surface temperature ',tsurfmean
        call error(3)
      endif

  110 format(i8,i8,f7.2)

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmstate
!-----------------------------------------------------------------------
! *** calculates atmospheric fields from potential vorticity
! *** field, mean atmospheric temperatures and the moisture field
!-----------------------------------------------------------------------
      implicit none

      call qtopsi
      call psitogeo
      call dyntemp
      call tempprofile
      call geowin
      call omega3
      call diver
      call divwin
      call totwind
      call totwind10
      call moisfields
      call cloud

      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wrenddyn
!-----------------------------------------------------------------------
! *** output atmosphere for start new run
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comunit.h'

      write(iuo+95) qprime,for

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wrendphy
!-----------------------------------------------------------------------
! *** output atmosphere for start new run
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comemic.h'

      write(iuo+95) tsurfn,tempm,temp0g
      write(iuo+95) rmoisg,torain,tosnow

      return
      end
