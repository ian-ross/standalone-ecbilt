












c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_co2oc(ist)
c----------------------------------------------------------------------
c *** communicate data from coupler to ocean
c----------------------------------------------------------------------


      include 'type.com'

      include 'comcouphelp.h'
      include 'comcoup.h'

      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'const.com'

      dimension zlatt(imax,jmax),zlont(imax,jmax)
      real*8 frwism(imax,jmax),deltamass,massGold,massAold
      real*8  day
      integer iyear,imonth,iday,iseason,nbclins,nbtrops,ntotday
      integer ntstep,nstpyear,nocstpyear
      integer nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,iobclin,
     *            nwrskip,nwrskip_days,end_year,end_day
      logical flgveg,flgicb,flgism,flgisma,flgismg

      COMMON/AG2CLIO2/frwism,deltamass,massGold,massAold 
      common /ec_coupl/ flgveg,flgicb,flgisma,flgismg
      common /ec_timectl/ day,ntstep,nstpyear,nocstpyear,
     *                 iyear,imonth,iday,iseason,ntotday,nbclins,nbtrops
      common /ec_timepar/ nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,iobclin,
     *                 nwrskip,nwrskip_days,end_year,end_day
      common / coord /  zlont,zlatt 


      real*8 zfld(imax,jmax)
      integer ix,jy,ist


c *** zonal wind stress 
      call at2oc(winstua,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
           tairox(ix,jy) = zfld(ix,jy)
           tenagx(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** meridional wind stress 
      call at2oc(winstva,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
           tairoy(ix,jy) = zfld(ix,jy)
           tenagy(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** short wave radiation over ocean 
      call at2oc(sumoswr,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
          fsolcn(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** net fresh water flux over ocean/sea-ice
      call at2oc(sumofwf,zfld)    
           areat=0.
           x1=280*pi/180
           x2=360*pi/180
           y1=20*pi/180
           y2=50*pi/180
      do ix = 1, imax
        do jy = 1, jmax
          fwat(ix,jy) = -zfld(ix,jy)*86400.*1000.
c ecrespin:mettre la valeur du flux eau douce dans Sv  
        if (ihyster.ne.0) then
          sv=xfreshsv(irunlabel+iyear-ihyster+1)
!  write(99,*) 'sv xfreshsv',sv, irunlabel+iyear-ihyster+1
        else
          sv=0.0
        endif
          if ((zlont(ix,jy).le.x2).and.(zlont(ix,jy).ge.x1).and.
     &    (zlatt(ix,jy).le.y2).and.(zlatt(ix,jy).ge.y1)) then
           areat=areat+(area(ix,jy)*tms(ix,jy,ks2))
          endif
        enddo
      enddo
c     write(*,*)'total',areat
      do ix = 1, imax
        do jy = 1, jmax
          if ((zlont(ix,jy).le.x2).and.(zlont(ix,jy).ge.x1).and.
     &    (zlatt(ix,jy).le.y2).and.(zlatt(ix,jy).ge.y1)) then
           fwat(ix,jy)=fwat(ix,jy)+((sv*1.E+6/areat)*86400.*1000.)
          endif
        enddo
      enddo

      call at2oc(sumohfx,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
          ratbqo(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** net snowflux over sea-ice 
      call at2oc(sumisno,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
          hnplbq(ix,jy) = max (zero,zfld(ix,jy)*86400.*1000.)
        enddo
      enddo
c *** short wave radiation over sea-ice 
      call at2oc(sumiswr,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
          fsolg(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** net heat flux over sea-ice without solar flux
      call at2oc(sumihfx,zfld)    
      do ix = 1, imax
        do jy = 1, jmax
          ratbqg(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** coefficients sensible-heat flux for implicit ice-temperature calculation
      call at2oc(sumicof,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          vabq(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** 10 meter height Wind Magnitude (m/s) for iceberg model  
      call at2oc(sumuv10,zfld)
      do ix = 1, imax
        do jy = 1, jmax
         temp_vent(ix,jy) = zfld(ix,jy)
        enddo
      enddo

c *** Surface pressure
      call at2oc(sumpress,zfld)
      do ix = 1, imax
        do jy = 1, jmax
         temp_press(ix,jy) = zfld(ix,jy)
        enddo
      enddo
c *** heat due to iceberg melting
      ficebergn=sumohsn
      ficebergs=sumohss
C     write(99,*) 'iceberg',ficebergn,ficebergs

C Rotation of the wind stress

      do ix = 1, imax
        do jy = 1, jmax
         tairox(ix,jy)  = tenagx(ix,jy)*xang2(ix,jy)-
     &                  tenagy(ix,jy)*xang1(ix,jy)
         tairoy(ix,jy)  = tenagy(ix,jy)*xang2(ix,jy)+
     &                  tenagx(ix,jy)*xang1(ix,jy)
        enddo
      enddo
      do ix = 1, imax
        do jy = 1, jmax
         tenagx(ix,jy)  = tairox(ix,jy)*1.4
         tenagy(ix,jy)  = tairoy(ix,jy)*1.4
        enddo
      enddo

      do ix = 1, imax
        do jy = 1, jmax
          fevabq(ix,jy) = 0.0
          fcscn(ix,jy) = 0.0
          flecn(ix,jy) = 0.0
        enddo
      enddo
      
c      zsumz1=0.0
c      do j=1,jmax
c        do i=is1(j),is2(j)
c          zsumz1=zsumz1+area(i,j)*tms(i,j,ks2)
c        enddo
c      enddo
c      zsumz=0.0
c      do j=1,jmax
c        do i=is1(j),is2(j)
c          zsumz=zsumz+area(i,j)*tms(i,j,ks2)*fwat(i,j)
c        enddo
c      enddo
c      write(97,*) zsumz/zsumz1,zsumz1

      return
 
      end

      subroutine at2oc(fin,fout)
      
      implicit none
      
      include 'comcouphelp.h'
      include 'comcoup.h'
      
      integer ji,jk,i,j
      real*8  fin(nlat,nlon),fout(ijocn),finf(ijatm),zsum
      
      ji=0
      do i=1,nlat
        do j=1,nlon
	  ji=ji+1
	  finf(ji)=fin(i,j)
	enddo
      enddo
      
      do ji = 1, ijocn
      zsum = 0.
        do jk = 1, komax
          zsum = zsum + wa2o(ji,jk) * finf(inda2o(ji,jk))
        enddo 
        fout(ji) = zsum
      enddo

      end

