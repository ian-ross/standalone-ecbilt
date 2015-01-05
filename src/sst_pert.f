c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_perttsurfn
c----------------------------------------------------
c *** calculates noise
c *** and adds it to tsurfn
c----------------------------------------------------

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'netcdf.inc'
      include 'st_noise.com'

C     Parameters and variables, the same as for nudging
      integer dimnudge,Nens
      character(len=256) Clioarea
      character(len=256) Cliogrid
      character(len=256) eoffile
      character(len=256) resfile
      character(len=256) alphafile
      character(len=256) fcostEPFfile
      character(len=256) stoutfile
      character(len=256) nudgefile
      real*8 relaxcoefNoise,relaxcoef,sigmaEPF, nudgemax
      integer nudgeperiod, nudgeFileStartY, nudgeFileStartD
      integer xpStartY, xpStartD, wrtstnoise
      integer inudge


C     Parameters for netcdf reading
      integer fileID, tsID, recID, Nerror, status, tdim, ntime
      parameter (ndims=3) ! number of dimensions of the netcdf file
      integer start(ndims), count(ndims)
      integer NCStartTime
      real*8  stnoiseSST(nlat,nlon)

C     Parameters for perturbations
      real*8 rnum
      integer idum, ryear, dumy

      common / tsurfn_nudgingParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,relaxcoefNoise,sigmaEPF,
     &                           fcostEPFfile,stoutfile,wrtstnoise,
     &                           nudgefile,relaxcoef,nudgemax,
     &                           nudgeperiod,nudgefileStartY,
     &                           nudgeFileStartD,xpStartY,xpStartD,Nens

      common /tsurfn_stnoise/ stnoiseSST

      write(iuo+66,*) " "
      write(iuo+66,*) 'ensemble',iens,numens
      if (iens.eq.1) then
c--- Reading of the parameters, the same as for the nudging
	call read_tsurfnParam('nudging.param')

c--- STOCHASTIC NOISE >>>
	if((eoffile.ne.'none').and.(resfile.ne.'none').and.(alphafile.ne.'none').and.(relaxcoefNoise.ne.0)) then
c--- Calculation of the proper time step for the time series st_noise, in months
	  NCStartTime=(((irunlabel+iyear)*360+(imonth-1)*30+iday-1)-
     &                 (nudgefileStartY*360+nudgeFileStartD-1))/30+1

	  call def_stnoiseSST(NCStartTime)
	else
	  do i=1,nlat
	    do j=1,nlon
              stnoiseSST(i,j)=0.
	    enddo
	  enddo
	endif
!	write(*,*)'stnoiseSST(10,10) =',stnoiseSST(10,10)
c--- <<< STOCHASTIC NOISE

c--- Adding perturbation to tsurfn according to a land/ocean/sea ice fraction of a grid cell
        do j=1,nlon
          do i=1,nlat
	    do nn=1,ntyps
              tsurfn(i,j,nn)=tsurfn(i,j,nn)+stnoiseSST(i,j)*fractn(i,j,nn)
	    enddo
          enddo
        enddo
      endif

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine read_tsurfnParam(nudgingParamFile)
c----------------------------------------------------
c *** Reads the parameters, the same as for the nudging
c----------------------------------------------------
      character(len=*) nudgingParamFile

C     Parameters and variables, the same as for the nudging
      integer dimnudge,Nens
      character(len=256) Clioarea
      character(len=256) Cliogrid
      character(len=256) eoffile
      character(len=256) resfile
      character(len=256) alphafile
      character(len=256) fcostEPFfile
      character(len=256) stoutfile
      character(len=256) nudgefile
      real*8 relaxcoefNoise,relaxcoef,sigmaEPF, nudgemax
      integer nudgeperiod, nudgeFileStartY, nudgeFileStartD
      integer xpStartY, xpStartD, wrtstnoise

      common / tsurfn_nudgingParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,relaxcoefNoise,sigmaEPF,
     &                           fcostEPFfile,stoutfile,wrtstnoise,
     &                           nudgefile,relaxcoef,nudgemax,
     &                           nudgeperiod,nudgefileStartY,
     &                           nudgeFileStartD,xpStartY,xpStartD,Nens

      open(iuo+60,file=nudgingParamFile,status='old',form='formatted')

      read(iuo+60,*)
      read(iuo+60,*) dimnudge
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') Clioarea
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') Cliogrid
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') eoffile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') resfile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') alphafile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) relaxcoefNoise
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) sigmaEPF
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') fcostEPFfile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') stoutfile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) wrtstnoise
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') nudgefile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) relaxcoef
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgemax
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgeperiod
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgeFileStartY
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgeFileStartD
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) xpStartY
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) xpStartD
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) Nens

      close(iuo+60)

!       if((eoffile.eq.'none').or.(resfile.eq.'none').or.(alphafile.eq.'none').or.(relaxcoefNoise.eq.0)) then
! 	write(*,*) "No sst stochastic noise"
!       else
! 	write(*,*) " "
! 	write(*,*) "Adding noise to the surface temperature (sst)"
! 	write(*,*) "Relaxation coefficient = ",relaxcoefNoise
!         write(*,*) "Surface stochastic noise files are "
!         write(*,*) trim(eoffile)
!         write(*,*) trim(resfile)
!         write(*,*) trim(alphafile)
!       endif

      end

c------------------------------
c ***Subroutine def_stnoise***
c------------------------------

      subroutine def_stnoiseSST(NCStartTime)

      implicit double precision (a-h,o-z)
      include 'para.com'
      include 'st_noise.com'
      include 'comcouphelp.h'
      include 'comemic.h'


      integer dimnudge,Nens
      character(len=256) Clioarea
      character(len=256) Cliogrid
      character(len=256) eoffile
      character(len=256) resfile
      character(len=256) alphafile
      character(len=256) fcostEPFfile
      character(len=256) stoutfile
      character(len=256) nudgefile
      real*8 relaxcoefNoise,relaxcoef,sigmaEPF,nudgemax
      integer nudgeperiod, nudgeFileStartY, nudgeFileStartD
      integer xpStartY, xpStartD, wrtstnoise

      integer NCStartTime, nmodes
      real*8  alpha_stnoise(nmodesMax_st), lambda_stnoise(nmodesMax_st)
      real*8  stnoiseSST(nlat,nlon)

      common / tsurfn_nudgingParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,relaxcoefNoise,sigmaEPF,
     &                           fcostEPFfile,stoutfile,wrtstnoise,
     &                           nudgefile,relaxcoef,nudgemax,
     &                           nudgeperiod,nudgefileStartY,
     &                           nudgeFileStartD,xpStartY,xpStartD,Nens

      common /tsurfn_stnoise/ stnoiseSST

      real*8  scal_relaxcoef
      integer fileObsID,tempObsID,Nerror
      integer start(4), count(4)
      integer i,j,k,ii,jj,idum,ntmax
      integer*4 timeArray(3)

      real*4 temp_resNC(nlon,nlat,1,1)
      real*4 temp_eofNC(nlon,nlat,1,1)

      call itime(timeArray)     ! Get the current time
      idum = timeArray(1)+timeArray(2)+timeArray(3)+Nens
!      idum = Nens

      open(10,file=alphafile,form='formatted')
      read(10,*) nmodes
      read(10,*) ntmax
      do i=1,nmodes
        read(10,*) alpha_stnoise(i)
      enddo
      close(10)

!       write(*,*)'NCStartTime, nmodes, ntmax =',NCStartTime, nmodes,ntmax
!       write(*,*)'alpha(1), alpha(nmodes) =',alpha_stnoise(1), alpha_stnoise(nmodes)

      if (nmodes.gt.nmodesMax_st)then
       write(*,*) "nmodes is greater than its maximum allowed value",nmodesMax_st
       write(*,*) "In this case, nmodes is equal to the maximum"
       nmodes=nmodesMax_st
      endif

      do i=1,nmodes
        lambda_stnoise(i)=gasdev2DSST(idum)
      enddo

      write(*,*)'lambda_stnoise(1) =',lambda_stnoise(1)

c--- Reading the residual file
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=NCStartTime
      if(NCStartTime.gt.ntmax) then
	start(4)=mod(NCStartTime,ntmax)
! 	write(*,*) "THIS SIMULATION MIGHT BE WRONG, since"
! 	write(*,*) "current time step exceeds the time dimension of the st.noise."
! 	write(*,*) "Therefore st.noise will be periodic with a period",ntmax
!         write(*,*) 'NCStartTime for noise is',start(4),'instead of',NCStartTime
!         write(*,*) " "
!         write(*,*) " "
      endif
      if(NCStartTime.le.0) then
	start(4)=mod(NCStartTime,ntmax)+ntmax
! 	write(*,*) "THIS SIMULATION MIGHT BE WRONG, since"
! 	write(*,*) "current time step exceeds the time dimension of the st.noise."
! 	write(*,*) "Therefore st.noise will be periodic with a period",ntmax
!         write(*,*) 'NCStartTime for noise is',start(4),'instead of',NCStartTime
!         write(*,*) " "
!         write(*,*) " "
      endif


c--- Remarks stnoiseTs(nlat,nlon) but temp_resNC(nlon,nlat,1,1) and temp_eofNC(nlon,nlat,1,1)
      count(1)=nlon
      count(2)=nlat
      count(3)=1
      count(4)=1

      fileObsID = NCOPN( resfile, NCnowrit , NCtest )
      tempObsID = NCVID( fileObsID, 'ts', Nerror )

      CALL NCVGT(fileObsID,tempObsID,start,count,temp_resNC,NCtest)
      CALL NCCLOS(fileObsID, NCtest)

      do i=1,nlat
        do j=1,nlon
!            stnoiseSST(i,j)=temp_resNC(j,i,1,1)
           stnoiseSST(i,j)=0.0
        enddo
      enddo

c--- Reading the EOF file
      fileObsID = NCOPN( eoffile, NCnowrit , NCtest )
      tempObsID = NCVID( fileObsID, 'ts', Nerror )
      do ii=1,nmodes
        start(4)=ii

        CALL NCVGT(fileObsID,tempObsID,start,count,temp_eofNC,NCtest)

        do i=1,nlat
          do j=1,nlon
              stnoiseSST(i,j)=stnoiseSST(i,j)
     *+lambda_stnoise(ii)*temp_eofNC(j,i,1,1)
          enddo
        enddo
      enddo
      CALL NCCLOS(fileObsID, NCtest)

      if (wrtstnoise.ge.360) then
	scal_relaxcoef=1.0
      else
	scal_relaxcoef=360.0/wrtstnoise
      endif

      do i=1,nlat
        do j=1,nlon
c--- Daily noise is monthly mean + random
!            stnoiseSST(i,j)=relaxcoefNoise/scal_relaxcoef*stnoiseSST(i,j)+gasdev2DSST(idum)
           stnoiseSST(i,j)=relaxcoefNoise/scal_relaxcoef*stnoiseSST(i,j)
        enddo
      enddo


      end

      double precision function usran2DSST(ir)
c
c   this subroutine generates random values between 0.0 and 1.0 using
c   an integer seed
c   it is based on the imsl routine ggubs.
c
c   double precision version
c
      implicit double precision (a-h,o-z)
      parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
      ir=abs(mod(da*ir,db)+0.5d0)
      usran2DSST=dfloat(ir)/dc
      return
      end

      double precision function gasdev2DSST(idum)
c
c   function gasdev2D
c   (Numerical Recipes, W.T.Vetterling, et.al., p.203)
c
c   description
c   ===========
c   returns a normally distributed deviate with zero mean and unit variance
c   using usran2DSST(idum) in file ranims.f
c   (formerly used ran1(idum) in file ran1ol.f)
c
c 09/01/88:converted to double precision
c 07/06/89:added save for gset (most computers dont need these save statements)
      implicit double precision (a-h,o-z)
      save iset, gset
      data iset/0/
      if (iset.eq.0) then
c***we don't have an extra deviate handy,
c***so pick two uniform random numbers in the square extending from -1 to +1
c***in each direction.
1       v1=2.d0*usran2DSST(idum)-1.d0
        v2=2.d0*usran2DSST(idum)-1.d0
        r=v1**2+v2**2
        if(r.ge.1.d0)go to 1
        fac=sqrt(-2.d0*log(r)/r)
        gset=v1*fac
        gasdev2DSST=v2*fac
        iset=1
      else
c***we do have an extra deviate handy,
c***so return it, and unset the flag.
        gasdev2DSST=gset
        iset=0
      endif
      return
      end
