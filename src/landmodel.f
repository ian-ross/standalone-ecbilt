c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iniland
c-----------------------------------------------------------------------
c *** initialises and sets parameters of the land model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comland.h'
      include 'comglobal.h'

      integer     i,j

c *** initialisation evaporation, bottom moisture and snow cover

      if (irunlabel.eq.0) then
        do j=1,nlon
          do i=1,nlat
            evfaca(i,j)=1d0
            bmoisg(i,j)=bmoism
            landsnow(i,j)=0
            dsnow(i,j)=0d0
            runofo(i,j)=0d0
          enddo
        enddo
      else
        read(91) bmoisg,runofo,dsnow,landsnow
      endif
       
      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then
          
c *** mountain heights set to zero above water surfaces

            rmount(i,j)=0d0
            qmount(i,j)=0d0
            bmoisg(i,j)=undef
          else
            tland(i,j)=tsurf(i,j)
          endif
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inirunoff
c-----------------------------------------------------------------------
c *** initialises the runof basins
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comland.h'
      include 'comglobal.h'

      integer     i,j,k,ias,ibas,iac
      character*1 ch(nlon),space

c *** computation of runoff masks
 
c *** asci number of small letter a

      ias=ichar('a') - 1
      iac=ichar('A') - 1
      do i=nlat,1,-1
        read (10,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          if (lsmask(i,j).eq.0) then      
            ilabas(i,j)=ichar(ch(j)) - ias
            if (ilabas(i,j).lt.1.or.ilabas(i,j).gt.24) then
              write(29,*) 'in lat-lon point ',i,j
              write(29,*) ch(j),ilabas(i,j)
              call error(16)
            endif
            iocbas(i,j)=0
          endif
        enddo
      enddo

      do i=nlat,1,-1
        read (10,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          if (lsmask(i,j).eq.1) then      
            iocbas(i,j)=ichar(ch(j)) - ias
            if (iocbas(i,j).lt.1.or.iocbas(i,j).gt.24) then
              iocbas(i,j)=0
            endif
            ilabas(i,j)=0
          endif
        enddo
      enddo

      do i=nlat,1,-1
        do j=1,nlon
          if (lsmask(i,j).eq.0) then      
            ch(j)=char(ilabas(i,j)+ias) 
          else
            if (iocbas(i,j).eq.0) iocbas(i,j)=ichar('0')-iac
            ch(j)=char(iocbas(i,j)+iac)
          endif
        enddo
        write(10,100) i-1,space,(ch(j),j=1,nlon)
      enddo

      rewind(10)

c *** computation of area of ocean runoff basins

      do ibas=1,nbasins
        arocbas(ibas)=0.
      enddo
      do i=1,nlat
        do j=1,nlon
          if (iocbas(i,j).gt.0) then
            arocbas(iocbas(i,j))=arocbas(iocbas(i,j)) + darea*cosfi(i)
          endif
        enddo
      enddo


 100  format(i4,65A1)
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine landtemp
c-----------------------------------------------------------------------
c *** computes surface land temperature
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comland.h'

      integer  i,j,itetel,mitetel,il,jl
      real*8   tsland,tsland1,tsland2,tol,fluxsumland,zbrent

      external fluxsumland

c *** tol is the wanted accuracy of the land temperature in degrees

      parameter (tol=0.1)

      common /landpoint/il,jl

      mitetel = 0

      do j=1,nlon
        do i=1,nlat

          if (lsmask(i,j).eq.0) then    

            il=i
            jl=j
            tsland=tland(il,jl)
            tsland1=tsland
            tsland2=tsland + 1.
            
            call zbrac(fluxsumland,tsland1,tsland2,itetel)
            if (itetel.eq.100) call error(7)
            tland(i,j)=zbrent(fluxsumland,tsland1,tsland2,tol,itetel)
            if (tland(i,j).lt.200.or.tland(i,j).gt.350) then
              write(100,*) 'tland out of range'       
              write(100,*) i,j,tland(i,j)
            endif
            if (itetel.eq.100) call error(8)
            if (itetel.gt.mitetel) mitetel=itetel

c *** in case of temperatures above zero in case of snowcover, set
c *** surface temperature to meltpoint

            if (dsnow(i,j).gt.0d0.and.tland(i,j).gt.tzero) then
	      tland(i,j)=tzero
            endif
            
          endif
        enddo
      enddo


c      write(100,*) mitetel
c      call flush(100)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function fluxsumland(tsland) 
c-----------------------------------------------------------------------
c *** computes sum of fluxes between the land and the atmosphere
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer il,jl
      real*8  fluxsumland,tsland,qsatss,qsat,efluxgp,hfluxgp

      common /landpoint/il,jl

c *** sensible heatflux

      hfluxgp=alphad*cdragv(il,jl)*uv10(il,jl)*(tsland-tempsg(il,jl))

c *** latent heat flux

      qsatss=qsat(1.d+5,tsland)

      if (tsland.gt.tzero) then
        efluxgp=alphav*cdragv(il,jl)*uv10(il,jl)*(qsatss-q10(il,jl))
      else
        efluxgp=alphas*cdragv(il,jl)*uv10(il,jl)*(qsatss-q10(il,jl))
      endif      
      efluxgp=evfaca(il,jl)*efluxgp  
       
      if (efluxgp.lt.0) efluxgp=0.

c *** sum of all fluxes

      fluxsumland=hesws(il,jl) - hfluxgp + dlrads(il,jl) -
     &            sboltz*tsland**4 - efluxgp


      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine landmois
c-----------------------------------------------------------------------
c *** landmodel
c *** computes: bottom moisture, snow coverage and runoff
c-----------------------------------------------------------------------

      call landprecip
      call runoff

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine landprecip
c-----------------------------------------------------------------------
c *** computes snow coverage and bottom moisture
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comland.h'


      integer i,j
      real*8  betam,dsnomel,dsnovap
      real*8  qsat

      betam=1./(rlatfus*rowat)

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.0) then

c *** precipitation

            if (tland(i,j).lt.tzero) then
              dsnow (i,j)=dsnow (i,j) + dtime*torain(i,j)
            else
              bmoisg(i,j)=bmoisg(i,j) + dtime*torain(i,j)
            endif

c ***       sublimation or evaporation

            if (dsnow(i,j).le.0.) then
              bmoisg(i,j)=bmoisg(i,j)-dtime*evap(i,j)
            else
	      dsnovap=dtime*evap(i,j)
 
              if (dsnovap.gt.dsnow(i,j)) then
                bmoisg(i,j)=bmoisg(i,j)-(dsnovap-dsnow(i,j))
                dsnow(i,j)=0.
              else
                dsnow(i,j)=dsnow(i,j)-dsnovap
              endif
            endif

c *** melting

            if (tland(i,j).eq.tzero.and.dsnow(i,j).gt.0.) then
              dsnomel=dtime*betam*(hesws(i,j)+dlrads(i,j)-
     *                ulrads(i,j)-eflux(i,j)-hflux(i,j))
	      if (dsnomel.gt.0.) then
                if (dsnomel.gt.dsnow(i,j)) dsnomel=dsnow(i,j)
                dsnow(i,j)=dsnow(i,j)-dsnomel
	      
c *** melt water to bottom moisture

                bmoisg(i,j)=bmoisg(i,j)+dsnomel
              
	      endif
	    endif

c *** if snowdepth above a thresshold, remove excessive snow 
c *** artificially through the bottom moisture. 
c *** neglect the heat involved

            if (dsnow(i,j).gt.dsnm) then
              bmoisg(i,j)=bmoisg(i,j) + dsnow(i,j)-dsnm
              dsnow(i,j) = dsnm
            endif

            if (bmoisg(i,j).lt.0d0) bmoisg(i,j)=0d0

          endif
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
c *** landsnow is used in the definition of albedo at the surface
          landsnow(i,j)=0
          if (lsmask(i,j).eq.0) then
          if (tland(i,j).lt.tzero.and.dsnow(i,j).gt.0.) landsnow(i,j)=1
          endif
        enddo
      enddo

      return
      end
        
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine runoff
c-----------------------------------------------------------------------
c *** computes runoff from rivers and distributes it over the ocean
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comland.h'

      real*8    runo(nbasins)
      integer   i,j,ibas



c ***  total runoff for each land basin

      do ibas=1,nbasins
        runo(ibas)=0d0
      enddo

      do i=1,nlat
        do j=1,nlon

          runofl(i,j)=0.

          if (lsmask(i,j).eq.0) then       

            if (bmoisg(i,j).gt.bmoism) then
              runofl(i,j)=(bmoisg(i,j)-bmoism)/dtime
              bmoisg(i,j)=bmoism
            endif

            runo(ilabas(i,j))=runo(ilabas(i,j)) + 
     *                            runofl(i,j)*darea*cosfi(i)
          endif
        enddo
      enddo

c ***  distribution of land runoff over the ocean

      do i=1,nlat
        do j=1,nlon
          if (iocbas(i,j).gt.0) then
            runofo(i,j)=runo(iocbas(i,j))/arocbas(iocbas(i,j))
          endif
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wrendland
c-----------------------------------------------------------------------
c *** output atmosphere for start new run
c-----------------------------------------------------------------------
      implicit none
      
      include 'comatm.h'
      include 'comphys.h'
      include 'comland.h'

      write(96) bmoisg,runofo,dsnow,landsnow

      return
      end


      

