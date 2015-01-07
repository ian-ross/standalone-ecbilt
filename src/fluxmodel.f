c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fluxmodel
c-----------------------------------------------------------------------
c *** computes radiational fluxes, fluxes between the surface 
c *** and the atmosphere and potential vorticity forcing due to
c *** diabatical heating and estimated ageostrophic terms
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      call radiation
      call senhflux
      call lathflux
      call vortfor

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine senhflux
c-----------------------------------------------------------------------
c *** sensible heatflux between surface and atmosphere
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j

c *** sensible heatflux (watt/m**2)

      do j=1,nlon
        do i=1,nlat

c frank flux ook over het ijs en over land

          hflux(i,j)=alphad*cdragv(i,j)*uv10(i,j)*
     &               (tsurf(i,j)-tempsg(i,j))
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lathflux
c-----------------------------------------------------------------------
c *** latent heatflux between surface and atmosphere
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comland.h'

      integer i,j
      real*8  qsat
 
c *** latent heatflux due to evaporation from surface (watt/m**2)
c *** and evaporation rate (m/s)
       
      do j=1,nlon
        do i=1,nlat

          if (lseaice(i,j).eq.0.and.lsmask(i,j).eq.1) then
            eflux(i,j)=alphav*cdragv(i,j)*uv10(i,j)*
     &                 (qsurf(i,j)-q10(i,j))
            if (eflux(i,j).lt.0) eflux(i,j)=0.
            evap(i,j)=eflux(i,j)/(rowat*rlatvap)
          else
            if (tsurf(i,j).ge.tzero) then
              eflux(i,j)=alphav*cdragv(i,j)*uv10(i,j)*
     &                   (qsurf(i,j)-q10(i,j))
              if (eflux(i,j).lt.0) eflux(i,j)=0.
              evap(i,j)=eflux(i,j)/(rowat*rlatvap)
            else
              eflux(i,j)=alphas*cdragv(i,j)*uv10(i,j)*
     &                   (qsurf(i,j)-q10(i,j))
              if (eflux(i,j).lt.0) eflux(i,j)=0.
              evap(i,j)=eflux(i,j)/(rowat*rlatsub)
            endif
          endif

c *** evaporation factor =1 over snow, over wet land maximal 1

          evfaca(i,j)=evfac
          if (lsmask(i,j).eq.0) then

	    if (dsnow(i,j).gt.0.) then 
              evfaca(i,j)=evfac
            else
              evfaca(i,j)=evfac*min(1d0,bmoisg(i,j)/bmoism)
            endif

            eflux(i,j)=evfaca(i,j)*eflux(i,j)
            evap (i,j)=evfaca(i,j)*evap (i,j)
          endif
        enddo
      enddo        
 
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ocfluxes
c-----------------------------------------------------------------------
c *** computes fluxes between the ocean and the atmosphere
c *** for forcing the ocean model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      call osenhflux
      call olathflux
      call olradflux

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine osenhflux
c-----------------------------------------------------------------------
c *** sensible heatflux between ocean/seaice and atmosphere
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j

c *** sensible heatflux (watt/m**2)

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then
            hflux(i,j)=alphad*cdragv(i,j)*uv10(i,j)*
     &                 (tsurf(i,j)-tempsg(i,j))
          endif
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine olathflux
c-----------------------------------------------------------------------
c *** latent heatflux between ocean/seaice and atmosphere
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
 
      real*8  qsat,qsatss
 
c *** latent heatflux due to evaporation from ocean (watt/m**2)
c *** evaporation rate over ocean/seaice is prescribed from inato
       
      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then
            qsatss=qsat(1.d+5,tsurf(i,j))

            if (lseaice(i,j).eq.0) then
              eflux(i,j)=alphav*cdragv(i,j)*uv10(i,j)*
     &                   (qsatss-q10(i,j))
              if (eflux(i,j).lt.0) eflux(i,j)=0.
            else
              if (tsurf(i,j).ge.tzero) then
                eflux(i,j)=alphav*cdragv(i,j)*uv10(i,j)*
     &                     (qsatss-q10(i,j))
                if (eflux(i,j).lt.0) eflux(i,j)=0.
              else
                eflux(i,j)=alphas*cdragv(i,j)*uv10(i,j)*
     &                     (qsatss-q10(i,j))
                if (eflux(i,j).lt.0) eflux(i,j)=0.
              endif
            endif

          endif
        enddo
      enddo        
 
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine olradflux
c-----------------------------------------------------------------------
c *** upward longwave radiation over ocean and seaice
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j

c *** upward longwave radiation (watt/m**2)

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then
            ulrads(i,j)=sboltz*tsurf(i,j)**4
          endif
        enddo
      enddo

      return
      end



