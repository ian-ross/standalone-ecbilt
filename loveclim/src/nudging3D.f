












c234567890123456789012345678901234567I89012345678901234567890123456789012
      subroutine nudging3D(ist)
c-------------------------------------------------------
c *** calculates the value of the nudging 
c *** added to scal (temp and salt)  in the routine clio 
c-------------------------------------------------------   
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'netcdf.inc'
      include 'st_noise.com'

      include 'comcouphelp.h'
      include 'comemic.h'

      include 'dynami.com'

      integer ist

C     Parameters and variables for the nudging      
      integer dimnudge,Nens
      character(len=256) Clioarea
      character(len=256) Cliogrid
      character(len=256) eoffile
      character(len=256) resfile
      character(len=256) alphafile
      character(len=256) stoutfile
      integer wrtstnoise
      character(len=256) obsfile
      character(len=256) coeffile
      real*8  nudgemaxSalt, nudgemaxTemp
      real*8  relaxcoefSalt(kmax) 
      real*8  relaxcoefTemp(kmax)
      integer nudgeperiod, obsfileStartY, obsfileStartD
      integer xpStartY, xpStartD
      integer inudge

      real*8  temp_obs(imax-2,jmax,kmax), salt_obs(imax-2,jmax,kmax) 
      real*8  diff4nudge_temp(imax,jmax,kmax), nudging3DTemp(imax,jmax,kmax)
      real*8  diff4nudgeMaxTemp
      real*8  diff4nudge_salt(imax,jmax,kmax), nudging3DSalt(imax,jmax,kmax)
      real*8  diff4nudgeMaxSalt
      real*8  stnoise3DTemp(imax,jmax,kmax), stnoise3DSalt(imax,jmax,kmax)
      real*8  alpha_stnoise(nmodesMax_st), lambda_stnoise(nmodesMax_st)
      real*8  stnoise3DTempAn(imax-2,jmax,kmax,12),  stnoise3DSaltAn(imax-2,jmax,kmax,12)
      real*8  nudging3DTempAn(imax-2,jmax,kmax,12), nudging3DSaltAn(imax-2,jmax,kmax,12)
      real*8  numbnois(12), numbnudg(12)


C     Parameters for netcdf reading             
      integer fileObsID, tempObsID,saltObsID, recID, Nerror
      parameter (ndims=4) ! number of dimensions of the netcdf file
      integer start(4), count(4)
      integer NCStartTime
      double precision temp_obsNC(imax-2, jmax, kmax, 1)
      double precision salt_obsNC(imax-2, jmax, kmax, 1)
      double precision undefValObsTemp,undefValObsTempMin,undefValObsTempMax
      double precision undefValObsSalt,undefValObsSaltMin,undefValObsSaltMax
      double precision undefValCoefTemp,undefValCoefTempMin,undefValCoefTempMax
      double precision undefValCoefSalt,undefValCoefSaltMin,undefValCoefSaltMax
      integer nmodes,ntmax,idum
      integer*4 timeArray(3) 

      common / ec_nudging3DParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,Nens,stoutfile,
     &                           wrtstnoise, 
     &                           obsfile,coeffile,
     &                           nudgemaxSalt,nudgemaxTemp,
     &                           relaxcoefSalt,relaxcoefTemp,
     &                           nudgeperiod,obsfileStartY,
     &                           obsfileStartD,xpStartY,xpStartD
     &                           undefValObsTemp,undefValObsTempMin,
     &                           undefValObsTempMax,
     &                           undefValObsSalt,undefValObsSaltMin,
     &                           undefValObsSaltMax,
     &                           undefValCoefTemp,undefValCoefTempMin,
     &                           undefValCoefTempMax,
     &                           undefValCoefSalt,undefValCoefSaltMin,
     &                           undefValCoefSaltMax

      common /ec_nudging3D/ nudging3DTemp, nudging3DSalt,
     &                      nudging3DTempAn, nudging3DSaltAn
      common /ec_defstnoise3D/ alpha_stnoise, lambda_stnoise,
     &                       nmodes, ntmax
      common /ec_stnoise3D/ stnoise3DTemp, stnoise3DSalt,
     &                    stnoise3DTempAn, stnoise3DSaltAn,
     &                    numbnois,numbnudg

      real*8 tzero
      tzero=273.15

!      write(*,*) 'time step',ist

      if(ist.eq.1) then

C***    Reading of the parameters for the nudging

        call read_nudging3DParam('nudging3D.param')
      endif

c--- STOCHASTIC NOISE
      if((eoffile.ne.'none').and.(resfile.ne.'none').and.(alphafile.ne.'none')) then
c--- Calculation of the proper time step for the time series st_noise
	 NCStartTime=(((irunlabel+iyear)*360+(imonth-1)*30+iday-1)-
     &                (obsfileStartY*360+obsfileStartD-1))/30+1

         call def_stnoise3D(ist,NCStartTime)
      else
        do i=1,imax-2
          do j=1,jmax
            do k=1,kmax
              stnoise3DTemp(i,j,k)=0.
              stnoise3DSalt(i,j,k)=0.
            enddo
          enddo
        enddo
      endif


      if(ist.eq.1) then
 
        if(obsfile.ne.'none') then
 
          fileObsID = NCOPN( obsfile, NCnowrit , NCtest )  
 	    
 	  tempObsID = NCVID( fileObsID, 'temp', Nerror )
 
          saltObsID = NCVID( fileObsID, 'salt', Nerror )
 
          call GetUndef4Var3D(fileObsID,tempObsID,undefValObsTemp,undefValObsTempMin,undefValObsTempMax)
 
          call GetUndef4Var3D(fileObsID,saltObsID,undefValObsSalt,undefValObsSaltMin,undefValObsSaltMax)
 
        endif
 
      endif
 
      nudging3DTemp=0.
      nudging3DSalt=0.
 
      if(obsfile.ne.'none') then
 
C ***     Calcul of the nudging values that will be added
C ***     to scal (ns=1 and ns=2)
 
 	 NCStartTime=(((irunlabel+iyear)*360+(imonth-1)*30+iday-1)-
     &                (obsfileStartY*360+obsfileStartD-1))/nudgeperiod+1
       
c ***     Opening and reading of the netcdf file containing the observations  
 	  
 	 fileObsID = NCOPN( obsfile, NCnowrit , NCtest )  
 	    
 	 tempObsID = NCVID( fileObsID, 'temp', Nerror )
           
         saltObsID = NCVID( fileObsID, 'salt', Nerror )
 
 	 start(1)=1
 	 start(2)=1
         start(3)=1
 	 start(4)=NCStartTime
 
 	 count(1)=imax-2
 	 count(2)=jmax
 	 count(3)=kmax
         count(4)=1
 
 	 CALL NCVGT( fileObsID, tempObsID, start, count, temp_obsNC, NCtest)
         CALL NCVGT( fileObsID, saltObsID, start, count, salt_obsNC, NCtest)

 	 do i=1, imax-2
 	   do j=1, jmax
             do k=1, kmax
 	    
 	       temp_obs(i,j,k)=temp_obsNC(i,j,k,1)
               salt_obs(i,j,k)=salt_obsNC(i,j,k,1)
 
 	     enddo
 	   enddo
 	 enddo    
 
 	 CALL NCCLOS( fileObsID, NCtest)
 
 	 do i=1, imax-2
 	   do j=1, jmax
             do k=1,kmax
 	    
 	      if(((temp_obs(i,j,k).lt.(undefValObsTempMax)).and.
     &            (temp_obs(i,j,k).gt.(undefValObsTempMin))).or.
     &           (isnan(temp_obs(i,j,k)))) then
 	        
 		diff4nudge_temp(i,j,k)=0.
 	
 	      else
 	      
 	        diff4nudge_temp(i,j,k)=temp_obs(i,j,k)-(scal(i,j,k,1)-tzero)
              endif
 
              if(((salt_obs(i,j,k).lt.(undefValObsTempMax)).and.
     &            (salt_obs(i,j,k).gt.(undefValObsTempMin))).or.
     &           (isnan(salt_obs(i,j,k)))) then
         
                 diff4nudge_salt(i,j,k)=0.
         
              else
         
                diff4nudge_salt(i,j,k)=salt_obs(i,j,k)-scal(i,j,k,2)
              endif
 
 	      if(abs(diff4nudge_temp(i,j,k)).gt.(nudgeMaxTemp/relaxcoefTemp(k))) then
                 if (diff4nudge_temp(i,j,k).gt.0.) then
                   nudging3DTemp(i,j,k)=nudgeMaxTemp
                 else if (diff4nudge_temp(i,j,k).lt.0.) then 
                   nudging3DTemp(i,j,k)=-nudgeMaxTemp
                 endif
 	      
              else
                 nudging3DTemp(i,j,k)=relaxcoefTemp(k)*diff4nudge_temp(i,j,k)  
 	      endif  
 
 
              if(abs(diff4nudge_salt(i,j,k)).gt.(nudgeMaxSalt/relaxcoefSalt(k))) then
                 if (diff4nudge_salt(i,j,k).gt.0.) then
                   nudging3DSalt(i,j,k)=nudgeMaxSalt
                 else if (diff4nudge_salt(i,j,k).lt.0.) then
                   nudging3DSalt(i,j,k)=-nudgeMaxSalt
                 endif
         
              else
                 nudging3DSalt(i,j,k)=relaxcoefSalt(k)*diff4nudge_salt(i,j,k)
              endif
 
             enddo
 	   enddo
	 enddo  

         if(ist.eq.1) then
           do ii=1,12
	     numbnudg(ii)=0.
             do i=1,imax-2
               do j=1,jmax
                 do k=1,kmax
                   nudging3DTempAn(i,j,k,ii)=0.
                   nudging3DSaltAn(i,j,k,ii)=0.
                 enddo
               enddo
             enddo
           enddo
	 endif
	 numbnudg(imonth)=numbnudg(imonth)+1.
         do i=1,imax-2
           do j=1,jmax
             do k=1,kmax
               nudging3DTempAn(i,j,k,imonth)=nudging3DTempAn(i,j,k,imonth)+nudging3DTemp(i,j,k)
               nudging3DSaltAn(i,j,k,imonth)=nudging3DSaltAn(i,j,k,imonth)+nudging3DSalt(i,j,k)
             enddo
           enddo
         enddo

      else
         do i=1,imax-2
           do j=1,jmax
             do k=1,kmax
               nudging3DTemp(i,j,k)=0.
               nudging3DSalt(i,j,k)=0.
             enddo
           enddo
         enddo
      endif
      
      if (ist.eq.wrtstnoise) then 
	 !Writing 3D nudging, noise
         call write_3Dstnoise
      endif 

      end subroutine nudging3D

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine read_nudging3DParam(nudgingParamFile)

      include 'para.com'

c----------------------------------------------------
c *** Reads the parameters for the nudging
c----------------------------------------------------
      character(len=15) nudgingParamFile

C     Parameters and variables for the nudging      
      integer nmodes
      integer dimnudge,Nens
      character(len=256) Clioarea
      character(len=256) Cliogrid
      character(len=256) eoffile
      character(len=256) resfile
      character(len=256) alphafile
      character(len=256) stoutfile
      integer wrtstnoise
      character(len=256) obsfile
      character(len=256) coeffile
      real*8  nudgemaxSalt, nudgemaxTemp
      real*8  relaxcoefSalt(kmax) 
      real*8  relaxcoefTemp(kmax)
      integer nudgeperiod, obsfileStartY, obsfileStartD
      integer xpStartY, xpStartD
      integer i 

      double precision undefValObsTemp,undefValObsTempMin,undefValObsTempMax
      double precision undefValObsSalt,undefValObsSaltMin,undefValObsSaltMax
      double precision undefValCoefTemp,undefValCoefTempMin,undefValCoefTempMax
      double precision undefValCoefSalt,undefValCoefSaltMin,undefValCoefSaltMax

      common / ec_nudging3DParam / dimnudge, Clioarea, Cliogrid,
     &                           eoffile,resfile,alphafile,Nens,stoutfile,
     &                           wrtstnoise, 
     &                           obsfile,coeffile,
     &                           nudgemaxSalt,nudgemaxTemp,
     &                           relaxcoefSalt,relaxcoefTemp,
     &                           nudgeperiod,obsfileStartY,
     &                           undefValObsTempMax,
     &                           undefValObsSalt,undefValObsSaltMin,
     &                           undefValObsSaltMax,
     &                           undefValCoefTemp,undefValCoefTempMin,
     &                           undefValCoefTempMax,
     &                           undefValCoefSalt,undefValCoefSaltMin,
     &                           undefValCoefSaltMax


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
      read(iuo+60,'(A256)') stoutfile
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) wrtstnoise
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,'(A256)') obsfile
      read(iuo+60,*)
      read(iuo+60,*)
      do i=1,kmax
        read(iuo+60,*) relaxcoefTemp(i)
      enddo
      read(iuo+60,*)
      read(iuo+60,*)
      do i=1,kmax
        read(iuo+60,*) relaxcoefSalt(i)
      enddo
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgemaxTemp
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgemaxSalt
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) nudgeperiod
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) obsfileStartY
      read(iuo+60,*)
      read(iuo+60,*)
      read(iuo+60,*) obsfileStartD
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
      if (dimnudge.ne.3) then 
        obsfile='none'
	eoffile='none'
      endif
      if(obsfile.eq.'none') then
!         if (dimnudge.ne.0) write(*,*) "No 3D nudging"
      else
        write(*,*) "Nudging of 3D temperature and salinity"
        write(*,*) "Max value of the nudging on temp. = ",nudgemaxTemp
        write(*,*) "Max value of the nudging on salt = ",nudgemaxSalt
        write(*,*) "Observation file = ",obsfile
      endif

      if((eoffile.eq.'none').or.(resfile.eq.'none').or.(alphafile.eq.'none')) then
!         if (dimnudge.ne.0) write(*,*) "No 3D stochastic noise"
      else
        write(*,*) "Stochastic noise applied to 3D temperature and salinity"
        write(*,*) trim(eoffile)
        write(*,*) trim(resfile)
        write(*,*) trim(alphafile)
      endif
      if (dimnudge.ne.0) write(*,*) " "
!       write(*,*) "Writing 3D st.noise and nudging after",wrtstnoise,"days to a file",trim(stoutfile)


      end

c------------------------------
c ***Subroutine GetUndef4Var3D***
c------------------------------

      subroutine GetUndef4Var3D(inputID, varID, undef, undefMin, undefMax)
 
      include 'netcdf.inc'

      integer inputID, varID
      double precision undef, undefMin, undefMax
      integer status
      double precision undefMinMax(2)

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

      end subroutine GetUndef4Var3D

c------------------------------
c ***Subroutine def_stnoise3D***
c------------------------------

      subroutine def_stnoise3D(ist,NCStartTime)

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
      integer NCStartTime, nmodes
      real*8  alpha_stnoise(nmodesMax_st), lambda_stnoise(nmodesMax_st) 
      real*8  stnoise3DTemp(imax,jmax,kmax), stnoise3DSalt(imax,jmax,kmax)
      real*8  stnoise3DTempAn(imax-2,jmax,kmax,12), stnoise3DSaltAn(imax-2,jmax,kmax,12)
      real*8  numbnois(12),numbnudg(12)

      common /ec_defstnoise3D/ alpha_stnoise, lambda_stnoise,
     &                       nmodes, ntmax

      common /ec_nudging3DParam/ dimnudge, Clioarea, Cliogrid,
     &                    eoffile,resfile,alphafile,Nens
     
      common /ec_stnoise3D/ stnoise3DTemp, stnoise3DSalt,
     &                    stnoise3DTempAn, stnoise3DSaltAn,
     &                    numbnois,numbnudg
    
      integer fileObsID,tempObsID,saltObsID,Nerror
      integer start(4), count(4)
      integer i,j,k,ii,jj,idum,ntmax 

      real*4 temp_resNC(imax-2,jmax,kmax,1)
      real*4 salt_resNC(imax-2,jmax,kmax,1)
      real*4 temp_eofNC(imax-2,jmax,kmax,1)
      real*4 salt_eofNC(imax-2,jmax,kmax,1)

      integer*4 timeArray(3) 

      call itime(timeArray)     ! Get the current time
      idum = timeArray(1)+timeArray(2)+timeArray(3)+Nens
      if(ist.eq.1) then
        do jj=1,12
	  numbnois(jj)=0.
          do i=1,imax-2
	     do j=1,jmax
	       do k=1,kmax
		 stnoise3DTempAn(i,j,k,jj)=0.
		 stnoise3DSaltAn(i,j,k,jj)=0.
	       enddo
	     enddo
	  enddo
        enddo

        open(10,file=alphafile,form='formatted')
        read(10,*) nmodes
        read(10,*) ntmax
        do i=1,nmodes
          read(10,*) alpha_stnoise(i)
        enddo
        close(10)
        write(*,*) " "
        write(*,*)'nmodes, ntmax =',nmodes,ntmax
        write(*,*)'alpha(1), alpha(nmodes) =',alpha_stnoise(1), alpha_stnoise(nmodes)

        if (nmodes.gt.nmodesMax_st)then
         write(*,*) "nmodes is greater than its maximum allowed value",nmodesMax_st
         write(*,*) "In this case, nmodes is equal to the maximum"
         nmodes=nmodesMax_st
        endif

        do i=1,nmodes
          lambda_stnoise(i)=gasdev3D(idum)
        enddo

      else
	! The same lambda_stnoise over a month
	! So that daily noise is the monthly mean + random
	if (iday.eq.1) then
	  do i=1,nmodes
	    lambda_stnoise(i)=alpha_stnoise(i)*lambda_stnoise(i)
     *+sqrt(1-alpha_stnoise(i)*alpha_stnoise(i))*gasdev3D(idum)
	  enddo
	endif
      endif

c--- Reading the residual file
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=NCStartTime
      if(NCStartTime.gt.ntmax) then
	write(*,*) "THIS SIMULATION MIGHT BE WRONG, since"
	write(*,*) "current time step exceeds the time dimension of the st.noise."
	write(*,*) "Therefore st.noise will be periodic with a period",ntmax
	start(4)=NCStartTime-ntmax
        write(*,*)'NCStartTime for noise is',start(4),'instead of',NCStartTime
      endif
      if(NCStartTime.lt.0) then
	write(*,*) "THIS SIMULATION MIGHT BE WRONG, since"
	write(*,*) "current time step exceeds the time dimension of the st.noise."
	write(*,*) "Therefore st.noise will be periodic with a period",ntmax
	start(4)=NCStartTime+ntmax
        write(*,*)'NCStartTime for noise is',start(4),'instead of',NCStartTime
      endif


      count(1)=imax-2
      count(2)=jmax
      count(3)=kmax
      count(4)=1
 
      fileObsID = NCOPN( resfile, NCnowrit , NCtest )  
      tempObsID = NCVID( fileObsID, 'temp', Nerror )
      saltObsID = NCVID( fileObsID, 'salt', Nerror )

      CALL NCVGT(fileObsID,tempObsID,start,count,temp_resNC,NCtest)
      CALL NCVGT(fileObsID,saltObsID,start,count,salt_resNC,NCtest)
      CALL NCCLOS(fileObsID, NCtest)

      do i=1,imax-2
        do j=1,jmax
          do jj=1,kmax
           stnoise3DTemp(i,j,jj)=temp_resNC(i,j,jj,1)
           stnoise3DSalt(i,j,jj)=salt_resNC(i,j,jj,1)
          enddo
        enddo
      enddo


c--- Reading the EOF file
      fileObsID = NCOPN( eoffile, NCnowrit , NCtest )
      tempObsID = NCVID( fileObsID, 'temp', Nerror )
      saltObsID = NCVID( fileObsID, 'salt', Nerror )
      do ii=1,nmodes
        start(4)=ii

        CALL NCVGT(fileObsID,tempObsID,start,count,temp_eofNC,NCtest)
        CALL NCVGT(fileObsID,saltObsID,start,count,salt_eofNC,NCtest)

        do i=1,imax-2
          do j=1,jmax
            do jj=1,kmax
              stnoise3DTemp(i,j,jj)=stnoise3DTemp(i,j,jj)
     *+lambda_stnoise(ii)*temp_eofNC(i,j,jj,1)
              stnoise3DSalt(i,j,jj)=stnoise3DSalt(i,j,jj)
     *+lambda_stnoise(ii)*salt_eofNC(i,j,jj,1)
            enddo
          enddo
        enddo
      enddo 
      CALL NCCLOS(fileObsID, NCtest)

c--- Daily noise is monthly mean + random
      do i=1,imax-2
        do j=1,jmax
          do k=1,kmax
! c--- The same noise over a month
! 	    stnoise3DTemp(i,j,k)=stnoise3DTemp(i,j,k)/30
! 	    stnoise3DSalt(i,j,k)=stnoise3DSalt(i,j,k)/30
	    stnoise3DTemp(i,j,k)=stnoise3DTemp(i,j,k)+gasdev3D(idum)
	    stnoise3DSalt(i,j,k)=stnoise3DSalt(i,j,k)+gasdev3D(idum)
          enddo
        enddo
      enddo 
      numbnois(imonth)=numbnois(imonth)+1.
      do i=1,imax-2
	do j=1,jmax
	   do k=1,kmax
	     stnoise3DTempAn(i,j,k,imonth)=stnoise3DTempAn(i,j,k,imonth)+stnoise3DTemp(i,j,k)
	     stnoise3DSaltAn(i,j,k,imonth)=stnoise3DSaltAn(i,j,k,imonth)+stnoise3DSalt(i,j,k)
	   enddo
	enddo
      enddo

      end

c------------------------------
c ***Subroutine write_3Dstnoise***
c------------------------------

      subroutine write_3Dstnoise

      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'netcdf.inc'

      include 'comcouphelp.h'
      include 'comemic.h'

      integer dimnudge,Nens
      character(len=256) Clioarea
      character(len=256) Cliogrid
      character(len=256) eoffile
      character(len=256) resfile
      character(len=256) alphafile
      character(len=256) stoutfile
      real*8  stnoise3DTemp(imax,jmax,kmax), stnoise3DSalt(imax,jmax,kmax)
      real*8  stnoise3DTempAn(imax-2,jmax,kmax,12), stnoise3DSaltAn(imax-2,jmax,kmax,12)
      real*8  nudging3DTemp(imax,jmax,kmax), nudging3DSalt(imax,jmax,kmax)
      real*8  nudging3DTempAn(imax-2,jmax,kmax,12), nudging3DSaltAn(imax-2,jmax,kmax,12)
      real*8  numbnois(12),numbnudg(12)

      common /ec_nudging3DParam/ dimnudge, Clioarea, Cliogrid,
     &                    eoffile,resfile,alphafile,Nens,stoutfile
      common /ec_stnoise3D/ stnoise3DTemp, stnoise3DSalt,
     &                    stnoise3DTempAn, stnoise3DSaltAn,
     &                    numbnois,numbnudg
      common /ec_nudging3D/ nudging3DTemp, nudging3DSalt,
     &                    nudging3DTempAn, nudging3DSaltAn
      

      integer fileNoiID,tempNoiID,saltNoiID,Nerror,tempNudID,saltNudID
      integer start(4),count(4),i,j,k,ii
      integer dimidx,dimidy,dimidz,dimidt,dd(4)
      integer varidx,varidy,varidz,varidt
      real*8  xvals(imax-2),yvals(jmax),zvals(kmax),tvals(12)

      do ii=1,12
        do i=1,imax-2
          do j=1,jmax
            do k=1,kmax
              if ((numbnudg(ii)).ne.0) then
		nudging3DTempAn(i,j,k,ii)=nudging3DTempAn(i,j,k,ii)/numbnudg(ii)
		nudging3DSaltAn(i,j,k,ii)=nudging3DSaltAn(i,j,k,ii)/numbnudg(ii)
	      endif
              if ((numbnois(ii)).ne.0) then
	        stnoise3DTempAn(i,j,k,ii)=stnoise3DTempAn(i,j,k,ii)/numbnois(ii)
	        stnoise3DSaltAn(i,j,k,ii)=stnoise3DSaltAn(i,j,k,ii)/numbnois(ii)
	      endif
            enddo
          enddo
        enddo
      enddo

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1

      count(1)=imax-2
      count(2)=jmax
      count(3)=kmax
      count(4)=12
      
c---  Get the grid
      Nerror=NF_OPEN(Cliogrid,NF_NOWRITE,fileNoiID)
c---  Get variable IDs
      Nerror=NF_INQ_VARID (fileNoiID,'lon',varidx)
      Nerror=NF_INQ_VARID (fileNoiID,'lat',varidy)
      Nerror=NF_INQ_VARID (fileNoiID,'tdepth',varidz)
c --- Get contents of dimension variables
      Nerror=NF_GET_VARA_DOUBLE (fileNoiID,varidx,start(1),count(1),xvals)
      Nerror=NF_GET_VARA_DOUBLE (fileNoiID,varidy,start(2),count(2),yvals)
      Nerror=NF_GET_VARA_DOUBLE (fileNoiID,varidz,start(3),count(3),zvals)
      Nerror=NF_CLOSE(fileNoiID)
      do i=1,12
        tvals(i)=1.0*i
      enddo

c---  Create a file 
      Nerror=NF_CREATE(stoutfile,NF_CLOBBER,fileNoiID)

c---  Define the dimensions
      Nerror=NF_DEF_DIM(fileNoiID,'lon',imax-2,dimidx)
      Nerror=NF_DEF_DIM(fileNoiID,'lat',jmax,dimidy)
      Nerror=NF_DEF_DIM(fileNoiID,'tdepth',kmax,dimidz)
      Nerror=NF_DEF_DIM(fileNoiID,'time',12,dimidt)

c---  Define the variables 
      dd(1)=dimidx
      Nerror=NF_DEF_VAR(fileNoiID,'lon',nf_double,1,dd,varidx)
      dd(1)=dimidy
      Nerror=NF_DEF_VAR(fileNoiID,'lat',nf_double,1,dd,varidy)
      dd(1)=dimidz
      Nerror=NF_DEF_VAR(fileNoiID,'tdepth',nf_double,1,dd,varidz)
      dd(1)=dimidt
      Nerror=NF_DEF_VAR(fileNoiID,'time',nf_double,1,dd,varidt)

      dd(1)=dimidx
      dd(2)=dimidy
      dd(3)=dimidz
      dd(4)=dimidt
      Nerror=NF_DEF_VAR(fileNoiID,'temp_noise',nf_double,4,dd,tempNoiID)
      Nerror=NF_DEF_VAR(fileNoiID,'salt_noise',nf_double,4,dd,saltNoiID)
      Nerror=NF_DEF_VAR(fileNoiID,'temp_nudg',nf_double,4,dd,tempNudID)
      Nerror=NF_DEF_VAR(fileNoiID,'salt_nudg',nf_double,4,dd,saltNudID)

c---  Change mode of netCDF operation
      Nerror=NF_ENDDEF(fileNoiID)

c---  Output the values of the variables (include dimension variables)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,varidx,start(1),count(1),xvals)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,varidy,start(2),count(2),yvals)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,varidz,start(3),count(3),zvals)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,varidt,start(4),count(4),tvals)

c---  Main variables values
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,tempNoiID,start,count,stnoise3DTempAn)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,saltNoiID,start,count,stnoise3DSaltAn)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,tempNudID,start,count,nudging3DTempAn)
      Nerror=NF_PUT_VARA_DOUBLE(fileNoiID,saltNudID,start,count,nudging3DSaltAn)

c---  Close the file
      Nerror=NF_CLOSE(fileNoiID)

      end

      double precision function usran3D(ir)
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
      usran3D=dfloat(ir)/dc
      return
      end

      double precision function gasdev3D(idum)
c
c   function gasdev3D
c   (Numerical Recipes, W.T.Vetterling, et.al., p.203)
c
c   description
c   ===========
c   returns a normally distributed deviate with zero mean and unit variance
c   using usran3D(idum) in file ranims.f
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
1       v1=2.d0*usran3D(idum)-1.d0
        v2=2.d0*usran3D(idum)-1.d0
        r=v1**2+v2**2
        if(r.ge.1.d0)go to 1
        fac=sqrt(-2.d0*log(r)/r)
        gset=v1*fac
        gasdev3D=v2*fac
        iset=1
      else
c***we do have an extra deviate handy,
c***so return it, and unset the flag.
        gasdev3D=gset
        iset=0
      endif
      return
      end
