      subroutine ec_inierror

      implicit none

      integer nwarns

      common /cerror/nwarns

      nwarns = 0

      return
      end

      subroutine ec_error(ierr)

      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comsurf.h'
      include 'comcoup.h'
      include 'comunit.h'

      integer ierr,nwarns,i,j
      common /cerror/nwarns

      if (ierr.lt.100) then

      if (ierr.eq.1) then
        write(iuo+29,*) 'error in routine surftl of landmodel'
        write(iuo+29,*) 'too many iterations (> 20)'
      endif

      if (ierr.eq.2) then
        write(iuo+29,*) 'error in routine surftl of landmodel'
        write(iuo+29,*) 'divergence of solution'
      endif

      if (ierr.eq.3) then
        write(iuo+29,*) 'error in routine test of atmdiag0.f'
        write(iuo+29,*) 'surface temperature out of range'
      endif

      if (ierr.eq.4) then
        write(iuo+29,*) 'error in routine iniglobal of initial0.f'
        write(iuo+29,*) 'error in reading gausspoints'
      endif

      if (ierr.eq.5) then
        write(iuo+29,*) 'error in routine rooster of oceandyn0.f'
        write(iuo+29,*) 'error in reading landseamask'
      endif

      if (ierr.eq.6) then
        write(iuo+29,*) 'error in routine seaice of icemodel0.f'
        write(iuo+29,*) 'tijs undefined'
      endif

      if (ierr.eq.7) then
        write(iuo+29,*) 'error in routine surftl of landmodel0.f'
        write(iuo+29,*) 'too many iterations (> 100) in zbrac'
      endif

      if (ierr.eq.8) then
        write(iuo+29,*) 'error in routine surftl of landmodel0.f'
        write(iuo+29,*) 'too many iterations (> 100) in zbrent'
      endif

      if (ierr.eq.9) then
        write(iuo+29,*) 'error in routine surfacetemp of icemodel0.f'
        write(iuo+29,*) 'too many iterations (> 100)'
      endif

      if (ierr.eq.10) then
        write(iuo+29,*) 'error in routine roostl of lakemodel0.f'
        write(iuo+29,*) 'error in reading lakemask'
      endif

      if (ierr.eq.11) then
        write(iuo+29,*) 'error in routine surfacetemp of icemodel0.f'
        write(iuo+29,*) 'too many iterations (> 100) in zbrac'
      endif

      if (ierr.eq.12) then
        write(iuo+29,*) 'error in routine surfacetemp of icemodel0.f'
        write(iuo+29,*) 'too many iterations (> 100) in zbrent'
      endif

      if (ierr.eq.13) then
        write(iuo+29,*) 'error in routine zbrent of root0.f'
        write(iuo+29,*) 'root not in specified interval'
      endif

      if (ierr.eq.15) then
        write(iuo+29,*) 'lake salinity out of range (<15)'
      endif

      if (ierr.eq.16) then
        write(iuo+29,*) 'error in routine iniland of landmodel0.f'
        write(iuo+29,*) 'land point not in a specified landbasin'
      endif

      if (ierr.eq.17) then
        write(iuo+29,*) 'error in reading ice mask in oceanfixed0.f'
      endif

      if (ierr.eq.18) then
        write(iuo+29,*) 'unrealistic ground pressure'
      endif

      if (ierr.eq.19) then
        write(iuo+29,*) 'too many ocean-lake neighbours in inioclacp'
      endif

      if (ierr.eq.20) then
        write(iuo+29,*) 'failure in expint called by detqmax'
      endif

      if (ierr.eq.99) then
        write(iuo+29,*) 'not a number detected in surfacetemperature'
      endif

      else

        if (ierr.eq.117) then
          write(iuo+29,*) ' longwave par. out of range too low temp'
        endif

        if (ierr.eq.118) then
          write(iuo+29,*) ' longwave par. out of range too high temp'
        endif

        if (ierr.eq.120) then
          write(iuo+29,*) 'convec: 10 iterations, still unstable'
        endif

        if (ierr.eq.121) then
          write(iuo+29,*) 'qmax out of range'
        endif

        if (ierr.eq.122) then
          write(iuo+29,*) 'rain larger than maximum set by rainmax'
        endif

        if (ierr.eq.123) then
          write(iuo+29,*) 'failure in netcdf export'
        endif

        if (ierr.eq.999) then
          write(iuo+29,*) 'integration succesfully completed'
        endif

        call flush(iuo+29)
        nwarns=nwarns+1
        if (nwarns.gt.100000) goto 10

        return

      endif

 10   continue

      write(iuo+29,*) 'program aborted due to soft error'
      write(iuo+29,*) 'see state.grads for state of the system'

      open(iuo+46,file='state.grads',form='unformatted')

      write(iuo+46) ((real(u200(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(v200(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(u500(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(v500(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(u800(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(v800(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(udivg(i,j,1)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(vdivg(i,j,1)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(udivg(i,j,2)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(vdivg(i,j,2)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(udivg(i,j,3)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(vdivg(i,j,3)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(temp0g(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(temp2g(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(temp4g(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(tempsg(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(tsurf(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(torain(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(tosnow(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(rmoisg(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(cormois(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(adsnow(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(ahic(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(eflux(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(hflux(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(hesws(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(dlrads(i,j)),j=1,nlon),i=1,nlat)
      write(iuo+46) ((real(ulrads(i,j)),j=1,nlon),i=1,nlat)
      close(iuo+46)
      stop

      return
      end
