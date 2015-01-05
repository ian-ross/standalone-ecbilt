c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_ecbilt(ist,jst)

#define LGM 0
#define LOCH 0

c-----------------------------------------------------------------------
c *** this routine performs one timestep of ECBilt which is a
c *** QG3L atmospheric model with simple physical parameterizations 
c *** 
c *** 6 april 1999 KNMI, De Bilt
c ***
c *** joint project Hugues Goosse
c ***               Rein Haarsma
c ***               Theo Opsteegh
c ***               Thijs van Reenen
c ***               Michiel Schaeffer
c ***               Frank Selten 
c ***               Xueli Wang
c ***               Nanne Weber
c ***
c *** Modified code : P. Mathiot (01/2012) Merge ec_topo with atmdyn0.f
c ***                                      Allow interannual topo forcing
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'

      integer ist,jst,istep
 
c *** atmospheric physics (in file atmphys.f)
 
      istep=(ist-1)*iatm+jst

c *** modif P.M.
c *** Update topo when ism model is used (each year and for initialization)
c *** Initialisation = .FALSE. every time step here, is not usefull 
c      if(flgism.AND.((mod(nint(day*real(iatm)),nstpyear).eq.0).or.(initialization.eq..true.))) call ec_topo(.TRUE.)
c *** Update top when ism model is not used
c *** UPDATE topo each year (flag flgism is include in ec_topo
      if (mod(nint(day*real(iatm)),nstpyear).eq.0) call ec_topo
c *** End modif P.M.
c
c      if(flgism.AND.((mod(nint(day*real(iatm)),nstpyear).eq.0).or.(initialization.eq..true.))) call ec_topo
c      !if (flgism.AND.(mod(istep,nstpyear).eq.1)) call ec_topo

      call ec_atmout(istep)
      call ec_checks(istep)

      if (iaphys.eq.1) then

        call ec_atmphyszero
        call ec_sensrad
c        call ec_tracer
        call ec_moisture 
        call ec_convec
        call ec_fortemp
        call ec_meantemp
      endif

c *** atmospheric dynamics (in file atmdyn.f) 
 
      if (iadyn.eq.1) then
        if (iartif.eq.1) call ec_forcdaily
        call ec_forward
      endif
      
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_update(ist,jst)
C--------------------------------------------------------------------------------
C ***
C *** This routine updates the model date, incoming solar radiation and
c *** the atmospheric state and does outputting and checking
C ***
C-------------------------------------------------------------------------------
      implicit none
      
      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'
      
      integer ist,jst,istep
      
      istep=(ist-1)*iatm+jst

      call ec_mdldate(istep)
      call ec_ghgupdate(istep) 
      call ec_solar(istep)
      call ec_atmstate
      call ec_vortfor

#if LGM == 1
c *** UPDATE irn (used by fluxes and vertical profile)
      if (mod(nint(day*real(iatm)),nstpyear).eq.0) call ec_irn
#endif

      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_checks(istep)
c-----------------------------------------------------------------------
c *** this routine performs checks on the model formulation
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'

      integer i,j,istep

      if ( mod(istep,iatm) .eq. 0) then
        call ec_testecbilt(istep)
      endif

c      call ec_moischeck(istep)
c      call ec_heatcheck(istep)

      call flush(99)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_testecbilt(istep)
c-----------------------------------------------------------------------
c *** testing if model variables are inside prescribed ranges
c-----------------------------------------------------------------------
      implicit none
  
      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'

      integer i,j,istep
      real*8  tsurfmean,ec_globalmean,dum1,dum2
      character*12 chts
      real*8  moc,tmc,tmc0,cland,thex
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex

c      dum1=ec_globalmean(ulrads)
c      dum2=ec_globalmean(dlrads)
c      write(*,*) istep,dum1,dum2

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

      tsurfmean=ec_globalmean(tsurf)-tzero

      write(chts,900) tsurfmean
 900  format(E12.5)
      if (chts(1:3).eq.'nan') call ec_error(99)

      write(iuo+20,110) iyear,int((day+0.5*dt)/(iatm*dt))+1,tsurfmean
      call flush(iuo+20)

      if (tsurfmean.gt.40..or.tsurfmean.lt.-10.) then
        write(iuo+29,*) 'mean surface temperature ',tsurfmean
        call ec_error(3)
      endif

  110 format(i8,i8,f7.2)
  100 format(f7.2)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_atmstate
c-----------------------------------------------------------------------
c *** calculates atmospheric fields from potential vorticity
c *** field, mean atmospheric temperatures and the moisture field
c-----------------------------------------------------------------------
      implicit none

      call ec_qtopsi 
      call ec_psitogeo
      call ec_dyntemp
      call ec_tempprofile
      call ec_geowin  
      call ec_omega3
      call ec_diver
      call ec_divwin
      call ec_totwind
      call ec_totwind10
      call ec_moisfields
      call ec_cloud

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_wrenddyn
c-----------------------------------------------------------------------
c *** output atmosphere for start new run
c-----------------------------------------------------------------------
      implicit none
      
      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comunit.h'

      write(iuo+95) qprime,for

      return
      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_wrendphy
c-----------------------------------------------------------------------
c *** output atmosphere for start new run
c-----------------------------------------------------------------------
      implicit none
      
      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comemic.h'

      write(iuo+95) tsurfn,tempm,temp0g
      write(iuo+95) rmoisg,torain,tosnow
#if LOCH == 1
      write(iuo+95) PGACO2,PCO2ref
#endif

      return
      end

 
