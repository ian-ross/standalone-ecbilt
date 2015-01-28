












c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine nudging(ist,jst)
c----------------------------------------------------
c *** calculates the value of the nudging 
c *** added to sumohfx in the routine ec_co2oc 
c----------------------------------------------------

      include 'comcouphelp.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'netcdf.inc'
      include 'st_noise.com'

      integer ist,jst

C     Parameters and variables for the nudging      
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

      real*8  tsurf4nudge(nlat,nlon)
      integer  count4nudge
      real*8  ts_obs(nlat,nlon), diff4nudge(nlat,nlon), nudgingTs(nlat,nlon)
      real*8 diff4nudgeMax

C     Parameters for netcdf reading             
      integer fileID, tsID, recID, Nerror, status, tdim, ntime
      parameter (ndims=3) ! number of dimensions of the netcdf file
      integer start(ndims), count(ndims)
      integer NCStartTime
      real*4 ts_obsNC(nlon, nlat, 1)
      real*4 undefVal,undefValMin,undefValMax

      real*8  alpha_stnoise(nmodesMax_st), lambda_stnoise(nmodesMax_st) 
      real*8  stnoiseTs(nlat,nlon)

      common / ec_nudgingParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,relaxcoefNoise,sigmaEPF,
     &                           fcostEPFfile,stoutfile,wrtstnoise, 
     &                           nudgefile,relaxcoef,nudgemax,
     &                           nudgeperiod,nudgefileStartY,
     &                           nudgeFileStartD,xpStartY,xpStartD,Nens,
     &                           undefVal,undefValMin,undefValMax

      common /ec_nudging/ tsurf4nudge,count4nudge,nudgingTs

      common /ec_defstnoise/ alpha_stnoise, lambda_stnoise,
     &                       nmodes, ntmax
      common /ec_stnoise/ stnoiseTs

      real*8 tzero

      inudge=6
      tzero=273.15



      if((ist.eq.1).and.(jst.eq.1)) then

C***    Reading of the parameters for the nudging

        call read_nudgingParam('nudging.param')

C***    Reading of the special values of the netcdf file
C       containing the observations

        if((relaxcoef.ne.0).and.(nudgefile.ne.'none')) then

          fileID = NCOPN( nudgefile, NCnowrit , NCtest )  
	    
	  tsID = NCVID( fileID, 'ts', Nerror )

          call GetUndef4Var2D(fileID,tsID,undefVal,undefValMin,undefValMax)

          tsurf4nudge(:,:)=0.
          count4nudge=0

        endif

      endif

c--- STOCHASTIC NOISE >>>
      if((eoffile.ne.'none').and.(resfile.ne.'none').and.(alphafile.ne.'none').and.(relaxcoefNoise.ne.0)) then
        if(jst.eq.inudge) then
c--- Calculation of the proper time step for the time series st_noise, in months
	  NCStartTime=(((irunlabel+iyear)*360+(imonth-1)*30+iday-1)-
     &                 (nudgefileStartY*360+nudgeFileStartD-1))/30+1

	  call def_stnoise(ist,NCStartTime)
	endif
      else
        do i=1,nlat
          do j=1,nlon
              stnoiseTs(i,j)=0.
          enddo
        enddo
      endif
c--- <<< STOCHASTIC NOISE 

      do i=1,nlat
        do j=1,nlon
          nudgingTs(i,j)=0.
        enddo
      enddo

      if((relaxcoef.ne.0).and.(nudgefile.ne.'none')) then

C***    Calcul of the daily mean of the surface temperature used for the nudging
        
        do ilat=1,nlat
          do ilon=1,nlon
            tsurf4nudge(ilat,ilon)=tsurf4nudge(ilat,ilon)+tsurf(ilat,ilon)
          enddo
        enddo

        count4nudge=count4nudge+1

        if(jst.eq.inudge) then

C ***     Calcul of the nudging value that will be added to sumohfx
	
	  NCStartTime=(((irunlabel+iyear)*360+(imonth-1)*30+iday-1)-
     &                (nudgeFileStartY*360+nudgeFileStartD-1))/nudgeperiod+1

	  diff4nudgeMax=abs(nudgemax/relaxcoef)

	        
c ***     Opening and reading of the netcdf file containing the observations  
	  
	  fileID = NCOPN( nudgefile, NCnowrit , NCtest )  
	    
	  tsID = NCVID( fileID, 'ts', Nerror )
	  
	  start(1)=1
	  start(2)=1
	  start(3)=NCStartTime
	  
	  count(1)=64
	  count(2)=32
	  count(3)=1	  

	  status=nf_inq_dimid(fileID,'time',tdim)
	  status=nf_inq_dimlen(fileID,tdim,ntime)
	  if (start(3).gt.ntime) then
	    do i=1, nlat
	      do j=1, nlon
	        ts_obs(i,j)=undefVal 
		stnoiseTs(i,j)=0.
	      enddo
	    enddo  
	  else
	    CALL NCVGT( fileID, tsID, start, count, ts_obsNC, NCtest )
 	    do i=1, nlat
	      do j=1, nlon
	        ts_obs(i,j)=ts_obsNC(j,i,1)
	      enddo
	    enddo 
	    CALL NCCLOS( fileID, NCtest)
	  endif

  
	  do i=1, nlat
	    do j=1, nlon
	    
	      tsurf4nudge(i,j)=tsurf4nudge(i,j)/count4nudge

c Nudging and noise are applied only over the same domain, which is defined in the nudging_file!!!
	      if(((ts_obs(i,j).lt.(undefValMax)).and.(ts_obs(i,j).gt.(undefValMin))).or.
     &            (isnan(ts_obs(i,j)))) then
	        
		diff4nudge(i,j)=0.
		stnoiseTs(i,j)=0.
	
	      else
	      
	        diff4nudge(i,j)=tsurf4nudge(i,j)-tzero-ts_obs(i,j)
              endif

	      if(abs(diff4nudge(i,j)).gt.diff4nudgeMax) then
                if (diff4nudge(i,j).gt.0.) then
                  nudgingTs(i,j)=relaxcoef*diff4nudgeMax
                else if (diff4nudge(i,j).lt.0.) then 
                  nudgingTs(i,j)=-relaxcoef*diff4nudgeMax
                endif
	      
              else
		nudgingTs(i,j)=relaxcoef*diff4nudge(i,j)  
	      endif  
		 
	    enddo
	  enddo 

          do i=1,nlat
            do j=1,nlon
	      tsurf4nudge(i,j)=0.
            enddo
          enddo

	  count4nudge=0

          
        endif

      endif

      end 

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine read_nudgingParam(nudgingParamFile)
c----------------------------------------------------
c *** Reads the parameters for the nudging
c----------------------------------------------------
      character(len=13) nudgingParamFile

C     Parameters and variables for the nudging      
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

      common / ec_nudgingParam / dimnudge, Clioarea, Cliogrid,
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

      if (dimnudge.ne.0) write(*,*) " "
      if (dimnudge.ne.2) then 
        nudgefile='none'
        relaxcoef=0.
	eoffile='none'
        relaxcoefNoise=0.
      endif
      if((nudgefile.eq.'none').or.(relaxcoef.eq.0)) then
!         if (dimnudge.ne.0) write(*,*) "No surface nudging"
      else
	write(*,*) "Nudging of the surface temperature"
! 	write(*,*) "Relaxation coefficient = ",relaxcoef
! 	write(*,*) "Max value of the nudging (in W/m2)= ",nudgemax
! 	write(*,*) "Observation file = ",nudgefile
      endif      
      if((eoffile.eq.'none').or.(resfile.eq.'none').or.(alphafile.eq.'none').or.(relaxcoefNoise.eq.0)) then
!         if (dimnudge.ne.0) write(*,*) "No surface stochastic noise"
      else
	write(*,*) "Adding noise to the surface temperature"
! 	write(*,*) "Relaxation coefficient = ",relaxcoefNoise
!         write(*,*) "Surface stochastic noise files are "
!         write(*,*) trim(eoffile)
!         write(*,*) trim(resfile)
!         write(*,*) trim(alphafile)
      endif
!        write(*,*) "Writing surface st.noise and nudging after",wrtstnoise,"days to a file",trim(stoutfile)

      end

c------------------------------
c ***Subroutine GetUndef4Var2D***
c------------------------------

      subroutine GetUndef4Var2D(inputID, varID, undef, undefMin, undefMax)
 
      include 'netcdf.inc'

      integer inputID, varID
      real*4 undef, undefMin, undefMax
      integer status
      real*4 undefMinMax(2)
      write(*,*) " "
      write(*,*) "Reading of the special value of the netcdf containing the obs. for the nudging"
      
      call NCAGT(inputID, varID, 'missing_value', undef, status)

      if(status.ne.0) then
        call NCAGT(inputID, varID, '_FillValue', undef, status)
        if(status.ne.0) then
          
          undef=-99.99

        endif
      endif
          
      if(undef.lt.0) then
        undefMin=undef*1.01
        undefMax=undef*0.99
      else
        undefMin=undef*0.99
        undefMax=undef*1.01
      endif

      write(*,*) "undefVal = ",undef
      write(*,*) "undefValMin = ",undefMin
      write(*,*) "undefValMax = ",undefMax
      write(*,*) " "

      end subroutine GetUndef4Var2D


c------------------------------
c ***Subroutine def_stnoise***
c------------------------------

      subroutine def_stnoise(ist,NCStartTime)

      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'netcdf.inc'
      include 'st_noise.com'

      include 'comcouphelp.h'
      include 'comemic.h'

      integer ist

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
      real*8  stnoiseTs(nlat,nlon)

      common /ec_defstnoise/ alpha_stnoise, lambda_stnoise,
     &                       nmodes, ntmax

      common / ec_nudgingParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,relaxcoefNoise,sigmaEPF,
     &                           fcostEPFfile,stoutfile,wrtstnoise, 
     &                           nudgefile,relaxcoef,nudgemax,
     &                           nudgeperiod,nudgefileStartY,
     &                           nudgeFileStartD,xpStartY,xpStartD,Nens

      common /ec_stnoise/ stnoiseTs
    
      integer fileObsID,tempObsID,Nerror
      integer start(4), count(4)
      integer i,j,k,ii,jj,idum,ntmax 
      integer*4 timeArray(3)
      real*8 diff4nudgeMax

      real*4 temp_resNC(nlon,nlat,1,1)
      real*4 temp_eofNC(nlon,nlat,1,1)


      call itime(timeArray)     ! Get the current time
      idum = timeArray(1)+timeArray(2)+timeArray(3)+Nens
      if(ist.eq.1) then
        open(10,file=alphafile,form='formatted')
        read(10,*) nmodes
        read(10,*) ntmax
        do i=1,nmodes
          read(10,*) alpha_stnoise(i)
        enddo
        close(10)
!         write(*,*) " "
!         write(*,*)'nmodes, ntmax =',nmodes,ntmax
!         write(*,*)'alpha(1), alpha(nmodes) =',alpha_stnoise(1), alpha_stnoise(nmodes)

        if (nmodes.gt.nmodesMax_st)then
         write(*,*) "nmodes is greater than its maximum allowed value",nmodesMax_st
         write(*,*) "In this case, nmodes is equal to the maximum"
         nmodes=nmodesMax_st
        endif

        do i=1,nmodes
          lambda_stnoise(i)=gasdev2D(idum)
        enddo

      else
	! The same lambda_stnoise over a month
	! So that daily noise is the monthly mean + random
! 	if (iday.eq.1) then
	  do i=1,nmodes
	    lambda_stnoise(i)=alpha_stnoise(i)*lambda_stnoise(i)
     *+sqrt(1-alpha_stnoise(i)*alpha_stnoise(i))*gasdev2D(idum)
	  enddo
! 	endif
      endif

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
!         write(*,*)'NCStartTime for noise is',start(4),'instead of',NCStartTime
      endif
      if(NCStartTime.le.0) then
	start(4)=mod(NCStartTime,ntmax)+ntmax
! 	write(*,*) "THIS SIMULATION MIGHT BE WRONG, since"
! 	write(*,*) "current time step exceeds the time dimension of the st.noise."
! 	write(*,*) "Therefore st.noise will be periodic with a period",ntmax
!         write(*,*)'NCStartTime for noise is',start(4),'instead of',NCStartTime
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
           stnoiseTs(i,j)=temp_resNC(j,i,1,1)
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
              stnoiseTs(i,j)=stnoiseTs(i,j)
     *+lambda_stnoise(ii)*temp_eofNC(j,i,1,1)
          enddo
        enddo
      enddo 
      CALL NCCLOS(fileObsID, NCtest)

      do i=1,nlat
        do j=1,nlon
! c--- The same noise over a month
!            stnoiseTs(i,j)=relaxcoefNoise*stnoiseTs(i,j)/30
! c--- Daily noise is monthly mean + random
!            stnoiseTs(i,j)=relaxcoefNoise*stnoiseTs(i,j)+gasdev2D(idum)
           stnoiseTs(i,j)=relaxcoefNoise*stnoiseTs(i,j)
        enddo
      enddo 


      end

      double precision function usran2D(ir)
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
      usran2D=dfloat(ir)/dc
      return
      end

      double precision function gasdev2D(idum)
c
c   function gasdev2D
c   (Numerical Recipes, W.T.Vetterling, et.al., p.203)
c
c   description
c   ===========
c   returns a normally distributed deviate with zero mean and unit variance
c   using usran2D(idum) in file ranims.f
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
1       v1=2.d0*usran2D(idum)-1.d0
        v2=2.d0*usran2D(idum)-1.d0
        r=v1**2+v2**2
        if(r.ge.1.d0)go to 1
        fac=sqrt(-2.d0*log(r)/r)
        gset=v1*fac
        gasdev2D=v2*fac
        iset=1
      else
c***we do have an extra deviate handy,
c***so return it, and unset the flag.
        gasdev2D=gset
        iset=0
      endif
      return
      end

