












c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_initlbm
c-----------------------------------------------------------------------
c *** initialises and sets parameters of the land model
c-----------------------------------------------------------------------
      implicit none


      include 'comland.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comsurf.h'
      
      integer     i,j,ija,is
      real*8      ds,db,dumwei
      real*8  phi(nlat)
      real*8  cosfi(nlat),sinfi(nlat),tanfi(nlat)
      real*8  dt,dtime,dtt,rdtime
      character*6 numyear
      character*3 numday

      common /ec_ctstep/ dt,dtime,dtt,rdtime

      NAMELIST /landpar/ bmoismfix,dsnm,lhcap,iscenland,islndstrt

      include 'openlainpfiles.h'


!       write(numyear,'(i6.6)') irunlabel+int((irunlabeld)/360)
!       write(numday,'(i3.3)') mod(irunlabeld,360)  
      write(numyear,'(i6.6)') irunlabel
      write(numday,'(i3.3)') irunlabeld 

c *** land fraction

      do j=1,nlon
        do i=1,nlat
          fractl(i,j)=1.-fracto(i,j)
	enddo
      enddo


c *** set land parameters


      pi=4d0*datan(1d0)
      radius=6.37e+6 
      rowat=1000.
      rlatvap=2.5e+06
      rlatsub=2.8e+06
      rlatfus=0.3e+06
      tzero  = 273.15
      
      bmoismfix = 0.15
      dsnm      = 1000.
      lhcap     = 2e6
      iscenland = 0
      islndstrt = 1990
      
      read(iuo+47, NML = landpar)
      
      write(iuo+30, 910) 'bmoism    =', bmoismfix
      write(iuo+30, 910) 'dsnm      =', dsnm
      write(iuo+30, 910) 'lhcap     =', lhcap
      write(iuo+30, 900) 'iscenland =', iscenland
      write(iuo+30, 900) 'islndstrt =', islndstrt
      
      call flush(iuo+30)
      
      rlhcap=1d0/lhcap

c *** grid of the surface is gaussian, read from dareafac
c *** which is initialised in initemic from file darea.dat

      do i=1,nlat
        dareas(i)= dareafac(i)
      enddo
      
      tareas=0d0
      do i=1,nlat
	tareas=tareas+nlon*dareas(i)
      enddo
      

c *** initialisation of evaporation, bottom moisture and snow cover

      if (irunlabel.eq.0) then
        do j=1,nlon
          do i=1,nlat
            bmoisg(i,j)=bmoismfix
            dsnow(i,j)=0d0
            runofo(i,j)=0d0
            tland(i,j)=tzero+10.
          enddo
        enddo
        heatsnown=0.0
        heatsnows=0.0
      else
        open(iuo+95,file='startdata/inland'//numyear//'_'//numday//'.dat',
     *    form='unformatted')
        read(iuo+95) bmoisg,runofo,dsnow,tland
 	close(iuo+95)
        do j=1,nlon
          do i=1,nlat
            dsnow(i,j)=min(dsnow(i,j),dsnm)
          enddo
        enddo
        heatsnown=0.0
        heatsnows=0.0

      endif
       
      do j=1,nlon
        do i=1,nlat
          if (fractl(i,j).lt.epsl) then
            bmoisg(i,j)=0d0
          endif
        enddo
      enddo
      
c *** read snow albedos
      
      read(iuo+4,*)
      do i=1,nlat
        read(iuo+4,*) albsnow(i)
      enddo
       
      call ec_landcoverupdate(1,1)

! Because of this flag, there are no initial conditions for the surface vegetation at all when VECODE is put to .FALSE.
! By commenting the if/endif lines in the initialization part of the landmodel0.f, the model can run with a non-dynamical constant vegetation cover.

!      if (flgveg) then

C Note: the value of patmCO2 has no impact during this call.
       call veget(1,1,dtime,1d-10,280,fractn(1,1,nld),
     &         dareas,tempsgn(1,1,nld))

!      endif



      call ec_landalbedoR(1)
      call ec_landalbedo(1)
      call ec_inirunoff

      
c *** land time step
      
      dtland=3600.*24./(iatm*ilan)
      rdtland=1d0/dtland

900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_inirunoff
c-----------------------------------------------------------------------
c *** initialises the runof basins
c-----------------------------------------------------------------------
      implicit none

      include 'comland.h'
      include 'comunit.h'
      
      integer     i,j,k,ias,ibas,iac
      character*1 ch(nlon),space

c *** computation of runoff masks
 
c *** asci number of small letter a
      do i=1,nlat
       do j=1,nlon
         iocbas(i,j)=0
         ilabas(i,j)=0
       enddo
      enddo

      nbasins=1
      ias=ichar('a') - 1
      iac=ichar('A') - 1
      do i=nlat,1,-1
        read (iuo+10,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          if (fractl(i,j).gt.epsl) then      
            ilabas(i,j)=ichar(ch(j)) - ias
            if (ilabas(i,j).lt.1.or.ilabas(i,j).gt.26) then
              write(iuo+29,*) 'in lat-lon point ',i,j
              write(iuo+29,*) ch(j),ilabas(i,j),fractl(i,j)
              call ec_error(16)
            endif
	    if (nbasins.lt.ilabas(i,j)) nbasins=ilabas(i,j)
          endif
        enddo
      enddo
      
      if (nbasins.gt.mbasins) then
        write(iuo+29,*) 
     & 'Error inirunoff: number of land basins greater than mbasins'
        STOP
      endif  

      do i=nlat,1,-1
        read (iuo+10,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          if ((1d0-fractl(i,j)).gt.epsl) then      
            iocbas(i,j)=ichar(ch(j)) - ias
            if (iocbas(i,j).lt.1.or.iocbas(i,j).gt.26) then
              iocbas(i,j)=0
            endif
          endif
        enddo
      enddo

      do i=nlat,1,-1
        do j=1,nlon
          if (fractl(i,j).gt.epsl) then
            ch(j)=char(ilabas(i,j)+ias) 
          else
            if (iocbas(i,j).eq.0) iocbas(i,j)=ichar('0')-iac
            ch(j)=char(iocbas(i,j)+iac)
          endif
        enddo
        write(iuo+10,100) i-1,space,(ch(j),j=1,nlon)
      enddo

      rewind(iuo+10)

c *** computation of area of ocean runoff basins

      do ibas=1,nbasins
        arocbas(ibas)=0.
      enddo
      do i=1,nlat
        do j=1,nlon
          if (iocbas(i,j).gt.0) then
            arocbas(iocbas(i,j))=arocbas(iocbas(i,j)) + dareas(i)
     &                               *(1.-fractl(i,j))
          endif
          if (iocbas(i,j).gt.nbasins) then
	    write(iuo+29,*) 'Error inirunof: iocbas out of range'
	    STOP
          endif
        enddo
      enddo
      do ibas=1,nbasins
        if (arocbas(ibas).eq.0.) then
	  write(iuo+29,*) 'Error inirunof: ocean basin empty ',ibas
	  STOP
        endif	
      enddo
      

 100  format(i4,65A1)
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_lbm(ist,jst,kst)
c-----------------------------------------------------------------------
c *** landmodel
c *** computes: bottom moisture, snow coverage, runoff, landtemperature
c-----------------------------------------------------------------------
      include 'comland.h'
      include 'comemic.h'

      integer ist,jst,kst,istep
      
      istep=(ist-1)*iatm+jst

      if (iscenland.eq.1) call ec_landcoverupdate(istep,kst)      
      call ec_landtemp
      call ec_landprecip
      call ec_runoff(istep)
      call ec_landalbedoR(istep)
      call ec_landalbedo(istep)
      call ec_landcheck(istep)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_landtemp
c-----------------------------------------------------------------------
c *** computes surface land temperature
c-----------------------------------------------------------------------
      implicit none

      include 'comland.h'

      integer  i,j
      real*8 betam

      betam=1./(rlatfus*rowat)
      do j=1,nlon
        do i=1,nlat
	
          meltheat(i,j)=0d0           
          if (fractl(i,j).gt.epsl) then    

            tland(i,j)=tland(i,j)+dtland*rlhcap*
     &         (nethfxland(i,j)+landheat(i,j))
         
c *** in case of temperatures above zero in case of snowcover, set
c *** surface temperature to meltpoint

            if (dsnow(i,j).gt.0d0.and.tland(i,j).gt.tzero) then
	      meltheat(i,j)=
     &          min(dsnow(i,j)/betam,lhcap*(tland(i,j)-tzero))
	      tland(i,j)=tland(i,j)-meltheat(i,j)/lhcap
	      meltheat(i,j)=meltheat(i,j)*rdtland
            endif
          endif
        enddo
      enddo

      return
      end



c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_landprecip
c-----------------------------------------------------------------------
c *** computes snow coverage and bottom moisture
c-----------------------------------------------------------------------
      implicit none

      include 'comland.h'
      include 'comunit.h'
      include 'comemic.h'

      integer i,j
      real*8  betam,dsnomel,dsnovap
      real*8  qsat

      betam=1./(rlatfus*rowat)
      heatsnows=0.0
      heatsnown=0.0
      
      do j=1,nlon
        do i=1,nlat
          landheat(i,j)=0d0
          if (fractl(i,j).gt.epsl) then

c ***       sublimation or evaporation

            if (dsnow(i,j).le.0.) then
              bmoisg(i,j)=bmoisg(i,j)-dtland*evapl(i,j)
            else
	      dsnovap=dtland*evapl(i,j)
 
              if (dsnovap.gt.dsnow(i,j)) then
                bmoisg(i,j)=bmoisg(i,j)-(dsnovap-dsnow(i,j))
                dsnow(i,j)=0.
              else
                dsnow(i,j)=dsnow(i,j)-dsnovap
              endif
            endif


c *** precipitation

            dsnow(i,j) =dsnow(i,j)  + dtland*snowf(i,j)
            bmoisg(i,j)=bmoisg(i,j) + dtland*rainf(i,j)

c *** melting: melt water to bottom moisture
c *** other imbalances are stored in landheat

            if (meltheat(i,j).gt.0d0) then
              dsnomel=dtland*betam*meltheat(i,j)
              if (dsnomel.gt.dsnow(i,j)) then
                landheat(i,j)=landheat(i,j)+
     *              (dsnomel-dsnow(i,j))/(dtland*betam)
                  dsnomel=dsnow(i,j)
              endif
              dsnow(i,j)=dsnow(i,j)-dsnomel
              bmoisg(i,j)=bmoisg(i,j)+dsnomel

c-Accumulation of ice sheet melt: dsnomeln et dsnomels en m**3
            endif

c *** if snowdepth above a thresshold, remove excessive snow
c *** artificially through the bottom moisture.
c *** account for the heat involved in landheat or in the ocean...

              if (dsnow(i,j).gt.dsnm) then
C               landheat(i,j)=landheat(i,j)-
C    *          (dsnow(i,j)-dsnm)/(dtland*betam)
                if (i.gt.15) then 
                  heatsnown=heatsnown+
     *             (dsnow(i,j)-dsnm)*dareas(i)*fractl(i,j)/betam
                else
                  heatsnows=heatsnows+
     *             (dsnow(i,j)-dsnm)*dareas(i)*fractl(i,j)/betam
                endif
                bmoisg(i,j)=bmoisg(i,j) + dsnow(i,j)-dsnm
                dsnow(i,j) = dsnm
              endif

            if (bmoisg(i,j).lt.0d0) then
	      if (bmoisg(i,j).lt.-1e-10) then
                write(iuo+99,*) 'bottom moisture less than zero '
                write(iuo+99,*) i,j,bmoisg(i,j)
                write(iuo+99,*) dtland*betam*meltheat(i,j)
                write(iuo+99,*) dtland*snowf(i,j)
                write(iuo+99,*) dtland*rainf(i,j)
                write(iuo+99,*) dtland*evapl(i,j)
                write(iuo+99,*) dsnow(i,j),dsnovap
              endif
              bmoisg(i,j)=0d0
            endif

          endif
        enddo
      enddo
C      write(iuo+28,*) 'heatsnow 1',heatsnown,heatsnows,
C    &                 heatsnows*betam
      
      return
      end
        
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_runoff(istep)
c-----------------------------------------------------------------------
c *** computes runoff from rivers and distributes it over the ocean
c-----------------------------------------------------------------------
      implicit none

      include 'comland.h'
      include 'comunit.h'
      include 'comemic.h'

      real*8    runo(nbasins)
      integer   i,j,ibas,istep



c ***  total runoff for each land basin

      do ibas=1,nbasins
        runo(ibas)=0d0
      enddo

      do i=1,nlat
        do j=1,nlon

          runofl(i,j)=0.

          if (fractl(i,j).gt.epsl) then

            if (bmoisg(i,j).gt.bmoism(i,j)) then
              runofl(i,j)=(bmoisg(i,j)-bmoism(i,j))/dtland
              bmoisg(i,j)=bmoism(i,j)
            endif
              runo(ilabas(i,j))=runo(ilabas(i,j)) + dareas(i)*
     *             runofl(i,j)*fractl(i,j)
          endif
        enddo
      enddo

c ***  distribution of land runoff over the ocean
C correction test Greenland, Antartcica
C     runo(5)=runo(5)*0.80
C     runo(23)=runo(23)*0.80

      do i=1,nlat
        do j=1,nlon
          runofo(i,j)=0.
          if (iocbas(i,j).gt.0) then
          runofo(i,j)=runo(iocbas(i,j))/arocbas(iocbas(i,j))
          endif
        enddo
      enddo
C outputs for Greenland and Antarctica
c     write(iuo+28,*) runo(5),runo(23)
      if (mod(istep,nstpyear).eq.1) then
       runo_yn=0.
       runo_ys=0.    
      endif
       runo_yn=runo_yn+(runo(5)*dtland)
       runo_ys=runo_ys+(runo(23)*dtland)

c      if (mod(istep,nstpyear).eq.0.) then
c        write(iuo+28,'(A25,I3,2E15.5)')'runoff land (m**3/y)',
c    &   iyear,runo_yn,runo_ys
c      endif

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_landcoverupdate(istep,kst)
c-----------------------------------------------------------------------
c *** updates land surface albedo and forest fraction
c *** once every 5 years
c *** starts when simulation year is equivalent to 1970, remains
c *** unchanged if simulation year exceeds 2100
c-----------------------------------------------------------------------
      implicit none

      include 'comland.h'
      include 'comemic.h'
      include 'comunit.h'

      integer i,j,is,scenyr,yr,istep,iy,kst
      real*8  d
      
      
c *** update once every 5 years
      
      if (mod(istep,5*nstpyear).eq.1.and.kst.eq.1) then
      
c *** scenario starts in 1970, after 2100 no updates
        if (iyear.eq.0) then
          scenyr=islndstrt
	else
	  scenyr=islndstrt+iyear-1
	endif
        if (scenyr.ge.1970.and.scenyr.le.2100) then
	
c *** read seasonal albedo data
	  read(iuo+31,90)yr
	  read(iuo+31,*)
	  if(yr.ne.scenyr) then
	    rewind(iuo+31)
	    read(iuo+31,90)yr
	    read(iuo+31,*)
	  endif
	  if(yr.ne.scenyr) then
	    call ec_forwardfile(iuo+31,scenyr,yr)
	  endif
          do i=1,nlat
            do j=1,nlon
              read(iuo+31,100)(albland(i,j,is),is=1,4)
            enddo
          enddo

          alblandR(:,:,:)=albland(:,:,:)
	  
c *** read yearly forest fraction data
	  read(iuo+32,90)yr
	  if(yr.ne.scenyr) then
	    rewind(iuo+32)
	    read(iuo+32,90)yr
	  endif
	  if(yr.ne.scenyr) then
	    call ec_forwardfile(iuo+32,scenyr,yr)
	  endif
          do i=1,nlat
            do j=1,nlon
              read(iuo+32,110)forestfr(i,j)
            enddo
          enddo
        endif


        d=0d0
        do j=2,25
          d=d+albland(27,j,1)
        enddo
        d=d/24.
        write(iuo+99,120) 'landcover update year: ',iyear,scenyr,yr,d
        do i=1,nlat
          do j=1,nlon
	    if (fractl(i,j).gt.epsl) then
	      do is=1,4
                if (albland(i,j,is).lt.0.01.or.
     *            albland(i,j,is).gt.0.99) then
                  write(iuo+99,130) 
     *            'Albedo of land out of range ',i,j,is,albland(i,j,is)
                endif
	      enddo
	    endif
          enddo
        enddo
		  
      endif
      call flush(iuo+99)
C *** FORMATS:

90    FORMAT(I5)
100   FORMAT(4(2X,F8.4))
110   FORMAT(F8.4)
120   FORMAT(A23,3I8,F12.5)
130   FORMAT(A28,3I4,F12.5)
      
      
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_forwardfile(filenr,year,yr)
c-----------------------------------------------------------------------
c *** forwards scenario ASCII file to find scenario year in 
c *** file with number filenr
c-----------------------------------------------------------------------
      implicit none
      
      include 'comland.h'
      include 'comunit.h'
      
      integer year,filenr,yr,i,j
      
      do i=1,nlat
        do j=1,nlon
          read(filenr,*)
        enddo
      enddo
      read(filenr,90)yr
      if (filenr.eq.iuo+31) read(filenr,*)
      dowhile(yr.ne.year.and.yr.lt.2100)
        do i=1,nlat
          do j=1,nlon
            read(filenr,*)
          enddo
        enddo
        read(filenr,90)yr
        if (filenr.eq.iuo+31) read(filenr,*)
      enddo
      
90    FORMAT(I5)
      
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_landalbedo(istep)
c-----------------------------------------------------------------------
c *** calculates albedos as a function of the time of the year, linearly
c *** interpolating between seasonal mean values
c *** albes  is the albedo of the earth surface, its value depends on 
c ***      whether it is a land or sea point, and on whether the grid 
c ***      point is snow or ice covered or not
c-----------------------------------------------------------------------
  
      implicit none

      include 'comland.h'
      include 'comemic.h'
      include 'comunit.h'

      integer i,j,id1,is1,is2,istep
      real*8  sfrac,albforfac,asnow,anosnow,snm,rsnm,snd

      integer imonth2,ktype,spv
      real*8  fracm,snowfrac,wcontent,forestfrac,hsnow
      real*8  albvec,albvecnosnow,albvecsnow,sndism
C     real*8  albsnow2

c *** reduction factor of snow albedo over forests

      albforfac=0.5
      
c *** if snowdepth exceeds snm, the albedo of snow
c *** is taken; for smaller snowdepths the albedo is
c *** a linear interpolation between the landalbedo
c *** and the albedo of snow.

      snm=0.05d0
      rsnm=1d0/snm
      
c *** interpolate between seasonal means: first shift day
c *** one to the middle of the first season (15 januari) and
c *** then interpolate

      id1=(imonth-1)*30+iday-14
      if (id1.lt.1) id1=id1+360

      is1=(id1+89)/90
      is2=is1+1
      if (is2.eq.5) is2=1

      sfrac=(id1-((is1-1)*90.+1.))/90.
      
      do j=1,nlon
        do i=1,nlat
	  anosnow=albland(i,j,is1)+
     *          (albland(i,j,is2)-albland(i,j,is1))*sfrac
          asnow=albsnow(i)*(1.-forestfr(i,j))+
     *                    forestfr(i,j)*albsnow(i)*albforfac
          snd=min(dsnow(i,j),snm)
	  
          alblbm(i,j)=rsnm*((snm-snd)*anosnow+snd*asnow)

        enddo
      enddo
            
      return
      end


      subroutine ec_landcheck(istep)
c-----------------------------------------------------------------------
c *** this routine checks the conservation of the moisture and
c *** heat budget over land
c-----------------------------------------------------------------------
      implicit none

      include 'comland.h'
      include 'comemic.h'
      include 'comunit.h'

      integer i,j,istep,nn
      real*8 iflux,oflux,labufch,labufcht,moisbdif
      real*8 totrunl,totruno,difrun,dum1,dum2,dum3
      real*8 ifluxt,ofluxt,abufchn,abufchnt,facmois,totrain,totevap
      real*8 difatmois,obufchnt,ocbud,totcor,totmois,globalmean
      
      real*8 bmoisgo(nlat,nlon),dsnowo(nlat,nlon),tlando(nlat,nlon)
      
      common /rlhch/bmoisgo,dsnowo,tlando

      if (istep.gt.2) then

        do j=1,nlon
          do i=1,nlat
            if (fractl(i,j).gt.epsl) then
              if (evapl(i,j)*dtland-bmoisgo(i,j)-dsnowo(i,j).gt.1E-15) 
     &          then
                write(iuo+99,*) 'error in evaporation over land '
                write(iuo+99,*) 'latlon ',i,j,bmoisgo(i,j)+
     &                dsnowo(i,j),evapl(i,j)*dtland
              endif
            endif
          enddo
        enddo


c *** local control for surface moisture budget

c       if (.not.flgism) then
        labufcht=0d0
        do j=1,nlon
          do i=1,nlat
            if (fractl(i,j).gt.epsl) then
              iflux=(rainf(i,j)+snowf(i,j)-evapl(i,j))
              oflux=runofl(i,j)
              labufch=(bmoisg(i,j)+dsnow(i,j))-(bmoisgo(i,j)+dsnowo(i,j))
              labufcht=labufcht + labufch*dareas(i)*fractl(i,j)
              moisbdif=labufch - dtland*(iflux-oflux)
              if (abs(labufch).gt.1d-6) then
                if (abs(moisbdif/labufch).gt.1d-5) then
                  write(iuo+99,*) 'error in landmoisture budget',
     &                        fractl(i,j)
                  write(iuo+99,*) 'latlon',i,j,moisbdif,labufch,
     &                     dtland*(iflux-oflux)
          	  write(iuo+99,*) rainf(i,j),evapl(i,j),runofl(i,j)
		  write(iuo+99,*) moisbdif/dtland,bmoisg(i,j)
		  write(iuo+99,*) 
                endif
              else
                if (abs(moisbdif).gt.1d-5) then
                  write(iuo+99,*) 'error in landmoisture budget'
                  write(iuo+99,*) 'latlon',i,j,moisbdif,labufch
                endif
              endif
            endif
          enddo
        enddo
c       endif
       
c *** global control for runoff budget

c- Not ok in the ice sheet-coupled version because
        totrunl=0d0
        totruno=0d0
        do j=1,nlon
          do i=1,nlat
            if (fractl(i,j).gt.epsl) then
              totrunl=totrunl + runofl(i,j)*dareas(i)*fractl(i,j)
	    endif
	    if ((1d0-fractl(i,j)).gt.epsl) then
              totruno=totruno + runofo(i,j)*dareas(i)*(1-fractl(i,j))
            endif
          enddo
        enddo

        totrunl=totrunl/tareas
        totruno=totruno/tareas
        difrun=totrunl-totruno

        if (abs(totruno).gt.0d0) then
          if (abs(difrun/totruno).gt.1d-5) then
            write(iuo+99,*) 'error in runoff budget'
            write(iuo+99,*) difrun,totrunl,totruno
          endif
        endif
       
      endif
      
c *** water conservation in the 0 coupled version (see forism.f)
       

c *** local control for surface heat budget

        labufcht=0d0
        do j=1,nlon
          do i=1,nlat
            if (fractl(i,j).gt.epsl) then
              iflux=(rainf(i,j)+snowf(i,j)-evapl(i,j))
              oflux=runofl(i,j)
              labufch=(tland(i,j)-tlando(i,j))*lhcap
              labufcht=labufcht + labufch*dareas(i)*fractl(i,j)
              moisbdif=labufch - dtland*(iflux-oflux)
              if (abs(labufch).gt.1d-6) then
                if (abs(moisbdif/labufch).gt.1d-5) then
c                  write(iuo+99,*) 'error in landmoisture budget',
c     &                        fractl(i,j)
c                  write(iuo+99,*) 'latlon',i,j,moisbdif,labufch,
c     &                     dtland*(iflux-oflux)
c          	  write(iuo+99,*) rainf(i,j),evapl(i,j),runofl(i,j)
c		  write(iuo+99,*) moisbdif/dtland,bmoisg(i,j)
c		  write(iuo+99,*) 
                endif
              else
                if (abs(moisbdif).gt.1d-5) then
c                  write(iuo+99,*) 'error in landmoisture budget'
c                  write(iuo+99,*) 'latlon',i,j,moisbdif,labufch
                endif
              endif
            endif
          enddo
        enddo

c *** storing of moisture variables 

       do j=1,nlon
         do i=1,nlat
	   bmoisgo(i,j)=bmoisg(i,j)
	   dsnowo(i,j)=dsnow(i,j)
	   tlando(i,j)=tland(i,j)
	 enddo
       enddo
      
       return
       end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_co2la
c-----------------------------------------------------------------------
c *** communicate between coupler and land
c-----------------------------------------------------------------------

      implicit none
      
      include 'comland.h'
      include 'comcoup.h'
      include 'comsurf.h'
      include 'comemic.h'
      
      integer i,j
      real*8 rsnrai
      real*8 temp_pr,temp_ev,tot_pr,pr_ism,tot_prnet
      real*8 betam

      betam=1./(rlatfus*rowat)

      evapoc(:,:)=evapn(:,:,noc)

      do j=1,nlon
        do i=1,nlat
	  nethfxland(i,j)=(heswsn(i,j,nld)+dlradsn(i,j,nld)
     &                    -ulradsn(i,j,nld)-efluxn(i,j,nld)
     &                    -hfluxn(i,j,nld))
	  evapl(i,j)=evapn(i,j,nld)
	  rainf(i,j)=couprf(i,j)
	  snowf(i,j)=coupsf(i,j)
          
c-Only the precip not taken into account by 0 are considered here

	enddo
      enddo
      
      end
	  

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_la2co
c-----------------------------------------------------------------------
c *** communicate land data to coupler
c-----------------------------------------------------------------------

      implicit none
      
      include 'comland.h'
      include 'comcoup.h'
      include 'comsurf.h'
      
      integer i,j
      
      do j=1,nlon
        do i=1,nlat
	  albesn(i,j,nld)=alblbm(i,j)
          albesnR(i,j)=alblbmR(i,j)
                          
	  tsurfn(i,j,nld)=tland(i,j)
	  abmoisg(i,j)=bmoisg(i,j)
	  abmoism(i,j)=bmoism(i,j)
	  adsnow(i,j)=dsnow(i,j)
	enddo
      enddo
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_lae2co
c-----------------------------------------------------------------------
c *** communicate land data to coupler
c-----------------------------------------------------------------------

      implicit none
      
      include 'comland.h'
      include 'comcoup.h'
      include 'comsurf.h'
      
      integer i,j
      
      do j=1,nlon
        do i=1,nlat
	  coupruno(i,j)=runofo(i,j)
	  couprunl(i,j)=runofl(i,j)
	enddo
      enddo
      couphsnn=heatsnown
      couphsns=heatsnows
      return
      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_wrendland
c-----------------------------------------------------------------------
c *** output land for start new run
c-----------------------------------------------------------------------
      implicit none
      
      include 'comland.h'
      include 'comunit.h'

      write(iuo+95) bmoisg,runofo,dsnow,tland

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_landalbedoR(istep)
c-----------------------------------------------------------------------
c *** calculates albedos as a function of the time of the year, linearly
c *** interpolating between seasonal mean values
c *** albes  is the albedo of the earth surface, its value depends on
c ***      whether it is a land or sea point, and on whether the grid
c ***      point is snow or ice covered or not
c-----------------------------------------------------------------------

      implicit none

      include 'comland.h'
      include 'comemic.h'
      include 'comunit.h'

      integer i,j,id1,is1,is2,istep
      real*8  sfrac,albforfac,asnow,anosnow,snm,rsnm,snd

      integer imonth2,ktype,spv
      real*8  fracm,snowfrac,wcontent,forestfrac,hsnow
      real*8  albvec,albvecnosnow,albvecsnow,sndism

c *** reduction factor of snow albedo over forests

      albforfac=0.5

      
*** if snowdepth exceeds snm, the albedo of snow
c *** is taken; for smaller snowdepths the albedo is
c *** a linear interpolation between the landalbedo
c *** and the albedo of snow.

      snm=0.05d0
      rsnm=1d0/snm

c *** interpolate between seasonal means: first shift day
c *** one to the middle of the first season (15 januari) and
c *** then interpolate

      id1=(imonth-1)*30+iday-14
      if (id1.lt.1) id1=id1+360

      is1=(id1+89)/90
      is2=is1+1
      if (is2.eq.5) is2=1

      sfrac=(id1-((is1-1)*90.+1.))/90.

      do j=1,nlon
        do i=1,nlat
          anosnow=alblandR(i,j,is1)+
     *          (alblandR(i,j,is2)-alblandR(i,j,is1))*sfrac
          asnow=albsnow(i)*(1.-forestfrR(i,j))+
     *                    forestfrR(i,j)*albsnow(i)*albforfac
          snd=min(dsnow(i,j),snm)

          alblbmR(i,j)=rsnm*((snm-snd)*anosnow+snd*asnow)

        enddo
      enddo

      return
      end
