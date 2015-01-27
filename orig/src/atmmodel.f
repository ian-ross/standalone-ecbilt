c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine initatmmodel
c-----------------------------------------------------------------------
c *** initialisation of the atmospheric dynamics and physics
c-----------------------------------------------------------------------
      implicit none

      call iatmdyn
      call iatmphys
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmstate
c-----------------------------------------------------------------------
c *** calculates atmospheric fields from potential vorticity
c *** field, mean atmospheric temperatures and the moisture field
c-----------------------------------------------------------------------
      implicit none

      call qtopsi
      call psitogeo
      call dyntemp
      call ptground
      call geowin
      call omega3
      call diver
      call divwin
      call totwind
      call totwind10
      call moisfields
      call dragcoef

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmmodel
c-----------------------------------------------------------------------
c *** atmosphere model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j

c *** atmospheric physics (in file atmphys.f)


      if (iaphys.eq.1.and.irunatm.eq.1) then

        call atmphyszero
        call sensrad
c        call tracer
        call moisture
        call convec
        call fortemp
        call meantemp
      endif

c *** atmospheric dynamics (in file atmdyn.f)

      if (iadyn.eq.1.and.irunatm.eq.1) then
        call forward
      endif

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine surfacetemp
c-----------------------------------------------------------------------
c *** this routine assigns values to tsurf over oceans, lakes,
c *** seaice and land
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comice.h'
      include 'comcoup.h'
      include 'comland.h'
      include 'comglobal.h'

      integer i,j

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then
            if (lseaice(i,j).eq.0) then
              tsurf(i,j)=stmix(i,j)
            else
              tsurf(i,j)=tijs(i,j)
            endif
          else
            if (irunatm.eq.1) tsurf(i,j)=tland(i,j)
          endif
        enddo
      enddo

      return
      end

c2345678901234567890123456789012345678901234567890123456789012345678901
      subroutine dragcoef
c----------------------------------------------------------------------
c *** drag coefficient
c *** depending on surface roughness and stability (Richardson number)
c----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      integer i,j

      real*8 cdrags,cdragl

      cdrags=cdrag
      cdragl=cdrag
      do j=1,nlon
        do i=1,nlat
          cdragw(i,j) = cwdrag
          if (lsmask(i,j).eq.1) then
            cdragv(i,j)=cdrags
          else
            cdragv(i,j)=cdragl
          endif
          if (tsurf(i,j).lt.tempsg(i,j)) then
            cdragv(i,j)=0.2*cdragv(i,j)
          endif
        enddo
      enddo

      return
      end




c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wrendatmos
c-----------------------------------------------------------------------
c *** output atmosphere for start new run
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'

      integer i,j

      write(95) qprime,for
      write(95) tsurf,tempm
      write(95) rmoisg,torain

      return
      end
