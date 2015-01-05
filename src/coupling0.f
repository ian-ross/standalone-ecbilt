c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_initcoup
c-----------------------------------------------------------------------
c *** initialises the coupling of ocean with atmosphere model
c-----------------------------------------------------------------------
      implicit none

#define ISOTOPE 0

      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comsurf.h'
      include 'comcoup.h'
      include 'comunit.h'

      integer i,j,nn,ijo,ija,k
      real*8  areafac
      character*6 numyear
      character*3 numday


!       write(numyear,'(i6.6)') irunlabel+int((irunlabeld)/360)
!       write(numday,'(i3.3)') mod(irunlabeld,360)
      write(numyear,'(i6.6)') irunlabel
      write(numday,'(i3.3)') irunlabeld

      if (irunlabel.eq.0) then
        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = 0.
            clhesw0(i,j) = 0.
            clhesw1(i,j) = 0.
            clhesw2(i,j) = 0.
            clulrad0(i,j)= 0.
            clulrad1(i,j)= 0.
            clulrad2(i,j)= 0.
            clulrads(i,j)= 0.
            cldlrads(i,j)= 0.
            clhflux(i,j) = 0.
            cleflux(i,j) = 0.
            clevap(i,j)  = 0.
            cohesws(i,j) = 0.
            cohesw0(i,j) = 0.
            cohesw1(i,j) = 0.
            cohesw2(i,j) = 0.
            coulrad0(i,j)= 0.
            coulrad1(i,j)= 0.
            coulrad2(i,j)= 0.
            coulrads(i,j)= 0.
            codlrads(i,j)= 0.
            cohflux(i,j) = 0.
            coeflux(i,j) = 0.
            coevap(i,j)  = 0.
            sumohfx(i,j) = 0.
            sumoswr(i,j) = 0.
            sumihfx(i,j) = 0.
            sumiswr(i,j) = 0.
            sumofwf(i,j) = 0.
#if ISOTOPE == 1
            sumoiof(i,j) = 0.
            sumoidf(i,j) = 0.
#endif
            sumisno(i,j) = 0.
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0.
            winstua(i,j) = 0.
            winstva(i,j) = 0.
            sumtx(i,j)   = 0.
            sumty(i,j)   = 0.
          enddo
        enddo
        sumohsn=0.
        sumohss=0.

      else

         open(iuo+95,file='startdata/incoup' // numyear // '_' //
     &               numday // '.dat', form='unformatted')
         read(iuo+95) sumohfx,sumoswr,sumihfx,sumiswr,sumofwf,sumisno,
     &              winstua,winstva,sumtx,sumty
         read(iuo+95) cohesws,cohesw0,cohesw1,cohesw2,coulrad0,coulrad1,
     &         coulrad2,coulrads,codlrads,cohflux,coeflux,coevap
        close(iuo+95)
        sumohsn=0.0
        sumohss=0.0

        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = 0.
            clhesw0(i,j) = 0.
            clhesw1(i,j) = 0.
            clhesw2(i,j) = 0.
            clulrad0(i,j)= 0.
            clulrad1(i,j)= 0.
            clulrad2(i,j)= 0.
            clulrads(i,j)= 0.
            cldlrads(i,j)= 0.
            clhflux(i,j) = 0.
            cleflux(i,j) = 0.
            clevap(i,j)  = 0.
#if ISOTOPE == 1
            sumoiof(i,j) = 0.
            sumoidf(i,j) = 0.
#endif
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0.
          enddo
        enddo
      endif

c *** ISOTOPE added without modifying the restart file.

c *** calculate total area for atmospheric gaussian grid darea
c *** dareafac is read in initemic.f from file darea.dat

      areafac=4*pi*radius**2/dsqrt(dble(nlon))
      tarea=0d0
      tareas=0d0

      do i=1,nlat
        darea(i)=dareafac(i)
        dareas(i)=darea(i)
        tarea=tarea + nlon*darea(i)
        tareas=tareas + nlon*dareas(i)
      enddo

c *** read adresses of neighbouring points and weights to
c *** interpolate between CLIO and ECBilt grid
c *** unit 45 connected to mozaic.w
c *** interpolate from ocean to atmosphere:
c *** do ji = 1, ijatm
c ***   zsum = 0.
c ***   do jk = 1, kamax
c ***     zsum = zsum + wo2a(ji,jk) * ocean(indo2a(ji,jk))
c ***   enddo
c ***   atmos(ji) = zsum
c *** enddo

c *** interpolate from atmosphere to ocean:
c *** do ji = 1, ijocn
c ***   zsum = 0.
c ***   do jk = 1, komax
c ***     zsum = zsum + wa2o(ji,jk) * atmos(indo2a(ji,jk))
c ***   enddo
c ***   ocean(ji) = zsum
c *** enddo

c *** kamax=14 komax=17 (see comcoup.h)

      read(iuo+45)
      read(iuo+45) ((indo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(iuo+45)
      read(iuo+45) ((wo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(iuo+45)
      read(iuo+45) ((inda2o(ijo,k),k=1,komax),ijo=1,ijocn)
      read(iuo+45)
      read(iuo+45) ((wa2o(ijo,k),k=1,komax),ijo=1,ijocn)


      do k=1,kamax
        do ija=1,ijatm
          if (indo2a(ija,k).eq.0) then
            indo2a(ija,k)=1
            wo2a(ija,k)=0d0
          endif
        enddo
      enddo

      do k=1,komax
        do ijo=1,ijocn
          if (inda2o(ijo,k).eq.0) then
            inda2o(ijo,k)=1
            wa2o(ijo,k)=0d0
          endif
        enddo
      enddo

      call ec_iniocbas

c *** initialisation of fluxes

      do nn=1,ntyps
        call ec_fluxes(nn)
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_co2at
c-----------------------------------------------------------------------
c     *** communicate coupler data to the atmosphere
!
! ######################################################################
!
! THIS IS WHERE THE FIXED SURFACE BOUNDARY CONDITIONS SHOULD BE IMPOSED
! ON THE ATMOSPHERE.
!
! ######################################################################
!
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comcoup.h'

      integer i,j,nn

      do j=1,nlon
        do i=1,nlat
          tsurf(i,j)=0.0
          albes(i,j) =0.0
          alb2es(i,j) =0.0
          qsurf(i,j)=0.0
          cdragv(i,j)=0.0
          do nn=1,ntyps
            tsurf(i,j) =tsurf(i,j)+fractn(i,j,nn)*tsurfn(i,j,nn)
            albes(i,j) =albes(i,j)+fractn(i,j,nn)*albesn(i,j,nn)
            alb2es(i,j)=alb2es(i,j)+fractn(i,j,nn)*alb2esn(i,j,nn)
            qsurf(i,j) =qsurf(i,j)+fractn(i,j,nn)*qsurfn(i,j,nn)
            cdragv(i,j)=cdragv(i,j)+fractn(i,j,nn)*cdragvn(i,j,nn)
          enddo
          hesws(i,j) =cohesws(i,j)  + clhesws(i,j)
          hesw0(i,j) =cohesw0(i,j)  + clhesw0(i,j)
          hesw1(i,j) =cohesw1(i,j)  + clhesw1(i,j)
          hesw2(i,j) =cohesw2(i,j)  + clhesw2(i,j)
          ulrad0(i,j)=coulrad0(i,j) + clulrad0(i,j)
          ulrad1(i,j)=coulrad1(i,j) + clulrad1(i,j)
          ulrad2(i,j)=coulrad2(i,j) + clulrad2(i,j)
          ulrads(i,j)=coulrads(i,j) + clulrads(i,j)
          dlrads(i,j)=codlrads(i,j) + cldlrads(i,j)
          hflux(i,j) =cohflux(i,j)  + clhflux(i,j)
          eflux(i,j) =coeflux(i,j)  + cleflux(i,j)
          evap(i,j)  =coevap(i,j)   + clevap(i,j)
c *** these two variables are communicated to the atmosphere only
c *** for the output routines
          arunofo(i,j)=sumro(i,j)
          arunofl(i,j)=sumrl(i,j)
        enddo
      enddo

      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_sumfluxland(kst)
c-----------------------------------------------------------------------
c *** accumulate surface fluxes between atmosphere and land
c-----------------------------------------------------------------------

      implicit none

      include 'comcouphelp.h'
      include 'comcoup.h'
      include 'comsurf.h'
      include 'comemic.h'

      integer kst,i,j
      real*8  rlan

      if (kst.eq.1) then
        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = 0.
            clhesw0(i,j) = 0.
            clhesw1(i,j) = 0.
            clhesw2(i,j) = 0.
            clulrad0(i,j)= 0.
            clulrad1(i,j)= 0.
            clulrad2(i,j)= 0.
            clulrads(i,j)= 0.
            cldlrads(i,j)= 0.
            clhflux(i,j) = 0.
            cleflux(i,j) = 0.
            clevap(i,j)  = 0.
            sumrl(i,j)   = 0.
            sumro(i,j)   = 0.
          enddo
        enddo
        sumhsn=0.
        sumhss=0.
      endif
      do j=1,nlon
        do i=1,nlat
          clhesws(i,j) = clhesws(i,j) + fractn(i,j,nld)*heswsn(i,j,nld)
          clhesw0(i,j) = clhesw0(i,j) + fractn(i,j,nld)*hesw0n(i,j,nld)
          clhesw1(i,j) = clhesw1(i,j) + fractn(i,j,nld)*hesw1n(i,j,nld)
          clhesw2(i,j) = clhesw2(i,j) + fractn(i,j,nld)*hesw2n(i,j,nld)
          clulrad0(i,j)= clulrad0(i,j)+ fractn(i,j,nld)*ulrad0n(i,j,nld)
          clulrad1(i,j)= clulrad1(i,j)+ fractn(i,j,nld)*ulrad1n(i,j,nld)
          clulrad2(i,j)= clulrad2(i,j)+ fractn(i,j,nld)*ulrad2n(i,j,nld)
          clulrads(i,j)= clulrads(i,j)+ fractn(i,j,nld)*ulradsn(i,j,nld)
          cldlrads(i,j)= cldlrads(i,j)+ fractn(i,j,nld)*dlradsn(i,j,nld)
          clhflux(i,j) = clhflux(i,j) + fractn(i,j,nld)*hfluxn(i,j,nld)
          cleflux(i,j) = cleflux(i,j) + fractn(i,j,nld)*efluxn(i,j,nld)
          clevap(i,j)  = clevap(i,j)  + fractn(i,j,nld)*evapn(i,j,nld)
          sumrl(i,j)   = sumrl(i,j)   + couprunl(i,j)
          sumro(i,j)   = sumro(i,j)   + coupruno(i,j)
        enddo
      enddo
      sumhsn=sumhsn+couphsnn
      sumhss=sumhss+couphsns

      if (kst.eq.ilan) then
        rlan=1d0/float(ilan)
        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = clhesws(i,j) *rlan
            clhesw0(i,j) = clhesw0(i,j) *rlan
            clhesw1(i,j) = clhesw1(i,j) *rlan
            clhesw2(i,j) = clhesw2(i,j) *rlan
            clulrad0(i,j)= clulrad0(i,j)*rlan
            clulrad1(i,j)= clulrad1(i,j)*rlan
            clulrad2(i,j)= clulrad2(i,j)*rlan
            clulrads(i,j)= clulrads(i,j)*rlan
            cldlrads(i,j)= cldlrads(i,j)*rlan
            clhflux(i,j) = clhflux(i,j) *rlan
            cleflux(i,j) = cleflux(i,j) *rlan
            clevap(i,j)  = clevap(i,j)  *rlan
            sumrl(i,j)   = sumrl(i,j)   *rlan
            sumro(i,j)   = sumro(i,j)   *rlan
          enddo
        enddo
        sumhsn=sumhsn *rlan
        sumhss=sumhss *rlan
C       sumhsn=sumhsn
C       sumhss=sumhss
      endif

      end


c23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine preccor
c-----------------------------------------------------------------------

c *** computes mean values of ocean basins
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdiag.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comcoup.h'

      integer i,j
      real*8 amcor(nbasa)
      real*8 sum1rf,sum2rf,sum1sf,sum2sf

C     (1)='ATL N  '
      amcor(1)=corAN
C     (2)='PAC N  '
      amcor(2)=corPN
C     (3)='ARCTIC '
      amcor(3)=corAC
C     (4)='INDIAN '
      amcor(4)=corID
C     (5)='ATL S  '
      amcor(5)=corAS
C     (6)='PAC S  '
      amcor(6)=corPS
C     (7)='ANTAR  '
      amcor(7)=corAA

c *** apply correction to rainfall (couprf)
c *** sum volume of water that is subtracted at higher latitudes over
c *** gridpoints that cover ocean

      sum1rf=0.0
      do j=1,nlon
       do i=1,nlat
C        if (iocbasa(i,j).gt.0.and.fracto(i,j).gt.epss) then
         if (iocbasa(i,j).gt.0.and.fracto(i,j).gt.0.1) then
           if (amcor(iocbasa(i,j)).lt.0 ) then
             sum1rf=sum1rf+couprf(i,j)*darea(i)*amcor(iocbasa(i,j))
             couprf(i,j)=couprf(i,j)*(1.+amcor(iocbasa(i,j)))
           endif
         endif
       enddo
      enddo

c *** sum area over which this water is to be spread, over full ocean
c *** grid points only
      sum2rf=0.0

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
             sum2rf=sum2rf+darea(i)*amcor(iocbasa(i,j))
           endif
         endif
       enddo
      enddo
c *** calculate water flux per meter squared

      sum1rf=sum1rf/sum2rf

c *** reduce rainfall accordingly

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
             couprf(i,j)=couprf(i,j)-sum1rf*amcor(iocbasa(i,j))
           endif
         endif
       enddo
      enddo

c *** apply correction to snowfall (coupsf)

      sum1sf=0.0
      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).gt.0.1) then
           if (amcor(iocbasa(i,j)).lt.0 ) then
             sum1sf=sum1sf+coupsf(i,j)*darea(i)*amcor(iocbasa(i,j))
             coupsf(i,j)=coupsf(i,j)*(1.+amcor(iocbasa(i,j)))
           endif
         endif
       enddo
      enddo

      sum2sf=0.0

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
             sum2sf=sum2sf+darea(i)*amcor(iocbasa(i,j))
           endif
         endif
       enddo
      enddo

      sum1sf=sum1sf/sum2sf

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
             coupsf(i,j)=coupsf(i,j)-sum1sf*amcor(iocbasa(i,j))
           endif
         endif
       enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_wrendcoup
c-----------------------------------------------------------------------
c *** output fluxes for start new run
c-----------------------------------------------------------------------
      implicit none

      include 'comcouphelp.h'
      include 'comcoup.h'
      include 'comunit.h'

      integer i,j
      real*8  dum(nlat,nlon,12)

      do i=1,nlat
        do j=1,nlon
           dum(i,j,1)   =cohesws(i,j)  + clhesws(i,j)
           dum(i,j,2)   =cohesw0(i,j)  + clhesw0(i,j)
           dum(i,j,3)   =cohesw1(i,j)  + clhesw1(i,j)
           dum(i,j,4)   =cohesw2(i,j)  + clhesw2(i,j)
           dum(i,j,5)   =coulrad0(i,j) + clulrad0(i,j)
           dum(i,j,6)   =coulrad1(i,j) + clulrad1(i,j)
           dum(i,j,7)   =coulrad2(i,j) + clulrad2(i,j)
           dum(i,j,8)   =coulrads(i,j) + clulrads(i,j)
           dum(i,j,9)   =codlrads(i,j) + cldlrads(i,j)
           dum(i,j,10)  =cohflux(i,j)  + clhflux(i,j)
           dum(i,j,11)  =coeflux(i,j)  + cleflux(i,j)
           dum(i,j,12)  =coevap(i,j)   + clevap(i,j)
        enddo
      enddo

      write(iuo+95) sumohfx,sumoswr,sumihfx,sumiswr,sumofwf,sumisno,
#if ISOTOPE == 1
     &winstua,winstva,sumtx,sumty,sumoiof,sumoidf
c      write(99,*)'sumoiof',sumoiof
c      write(99,*)'sumoidf',sumoidf
c      write(99,*)'sumofwf',sumofwf
c      write(99,*)'sumofwf88',sumofwf(8,8)
c      write(99,*)'sumoiof88',sumoiof(8,8)
c      write(99,*)'sumoidf88',sumoidf(8,8)
#else
     &winstua,winstva,sumtx,sumty
#endif
      write(iuo+95) dum
      return
      end
