












c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_initcoup
c-----------------------------------------------------------------------
c *** initialises the coupling of ocean with atmosphere model
c-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comsurf.h'
      include 'comcoup.h'
      include 'comunit.h'
     
      integer i,j,g05dyf,nn,ijo,ija,k
      real*8  areafac,dumwei,ds,db
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
      
	open(iuo+95,file='startdata/incoup'//numyear//'_'//numday//'.dat',
     *              form='unformatted')
	read(iuo+95) sumohfx,sumoswr,sumihfx,sumiswr,sumofwf,sumisno,
     &winstua,winstva,sumtx,sumty
	read(iuo+95) cohesws,cohesw0,cohesw1,cohesw2,coulrad0,coulrad1,
     *         coulrad2,coulrads,codlrads,cohflux,coeflux,coevap
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
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0.  
          enddo
	enddo
      endif

c *** 0 added without modifying the restart file.

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
     
      call ec_la2co
      call ec_lae2co
      call ec_oc2co(1)
      call ec_at2co

c *** initialisation of fluxes
      
      do nn=1,ntyps
        call ec_fluxes(nn)
      enddo
      

900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_at2co
c-----------------------------------------------------------------------
c *** communicate coupler data to the atmosphere
c-----------------------------------------------------------------------
      implicit none
      
      include 'comatm.h'
      include 'comphys.h'
      include 'comcoup.h'
      
      integer i,j
      
      do j=1,nlon
        do i=1,nlat
	  couprf(i,j)=torain(i,j)
	  coupsf(i,j)=tosnow(i,j)
	  couptcc(i,j)=tcc(i,j)
	  
 	enddo
      enddo
      
c *** calculate flux correction to snow and rain over oceans
c *** to deal with too much precipitation in the arctic in ecbilt

      call preccor
      
      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_co2at
c-----------------------------------------------------------------------
c *** communicate coupler data to the atmosphere
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
      subroutine ec_sumfluxocean(ist,jst)
c-----------------------------------------------------------------------
c *** accumulate surface fluxes for the ocean-atmosphere interface
c-----------------------------------------------------------------------

      implicit none
      
      include 'comatm.h'
      include 'comphys.h'
      include 'comcoup.h'
      include 'comsurf.h'
      include 'comemic.h'
      
      integer ist,jst,i,j
      real*8  ratm,facstr,rndws,dwsm1,facsnow
      
      facsnow=rowat*rlatfus
      
      do j=1,nlon
	do i=1,nlat
          cohesws(i,j) = fractn(i,j,nse)*heswsn(i,j,nse)+
     &                                  fractn(i,j,noc)*heswsn(i,j,noc)
          cohesw0(i,j) = fractn(i,j,nse)*hesw0n(i,j,nse)+
     &                                  fractn(i,j,noc)*hesw0n(i,j,noc)
          cohesw1(i,j) = fractn(i,j,nse)*hesw1n(i,j,nse)+
     &                                  fractn(i,j,noc)*hesw1n(i,j,noc)
          cohesw2(i,j) = fractn(i,j,nse)*hesw2n(i,j,nse)+
     &                                  fractn(i,j,noc)*hesw2n(i,j,noc)
          coulrad0(i,j)= fractn(i,j,nse)*ulrad0n(i,j,nse)+
     &                                  fractn(i,j,noc)*ulrad0n(i,j,noc)
          coulrad1(i,j)= fractn(i,j,nse)*ulrad1n(i,j,nse)+
     &                                  fractn(i,j,noc)*ulrad1n(i,j,noc)
          coulrad2(i,j)= fractn(i,j,nse)*ulrad2n(i,j,nse)+
     &                                  fractn(i,j,noc)*ulrad2n(i,j,noc)
          coulrads(i,j)= fractn(i,j,nse)*ulradsn(i,j,nse)+
     &                                  fractn(i,j,noc)*ulradsn(i,j,noc)
          codlrads(i,j)= fractn(i,j,nse)*dlradsn(i,j,nse)+
     &                                  fractn(i,j,noc)*dlradsn(i,j,noc)
          cohflux(i,j) = fractn(i,j,nse)*hfluxn(i,j,nse)+
     &                                  fractn(i,j,noc)*hfluxn(i,j,noc)
	  coeflux(i,j) = fractn(i,j,nse)*efluxn(i,j,nse)+
     &                                  fractn(i,j,noc)*efluxn(i,j,noc)
	  coevap(i,j)  = fractn(i,j,nse)*evapn(i,j,nse)+
     &                                  fractn(i,j,noc)*evapn(i,j,noc)
        enddo
      enddo

      if ((jst.eq.1.and.ist.eq.1).or.iobclint.eq.0) then
        do j=1,nlon
	  do i=1,nlat
            sumohfx(i,j) = 0. 
            sumoswr(i,j) = 0. 
            sumihfx(i,j) = 0. 
            sumiswr(i,j) = 0. 
            sumofwf(i,j) = 0. 
            sumisno(i,j) = 0. 
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0. 
	  enddo
	enddo
        sumohsn=0.0
        sumohss=0.0
	iobclint=0
      endif
      
      if ((jst.eq.1.and.ist.eq.1).or.iobtropt.eq.0) then
        do j=1,nlon
	  do i=1,nlat
            winstua(i,j) = 0.
            winstva(i,j) = 0.
	  enddo
	enddo
	iobtropt=0
      endif

      iobtropt=iobtropt+1
      iobclint=iobclint+1
      
      do j=1,nlon
	do i=1,nlat
	  if (fracto(i,j).gt.epss) then
            sumohfx(i,j) = sumohfx(i,j)+
     &      		   dlradsn(i,j,noc)-ulradsn(i,j,noc)
     &      		  -efluxn(i,j,noc)-hfluxn(i,j,noc)-
     &      		  facsnow*coupsf(i,j)
            sumoswr(i,j) = sumoswr(i,j)+ heswsn(i,j,noc)
            sumihfx(i,j) = sumihfx(i,j)+
     &      		   dlradsn(i,j,nse)
     &      		  -efluxn(i,j,nse)-hfluxn(i,j,nse)
            sumiswr(i,j) = sumiswr(i,j)+ heswsn(i,j,nse)
            sumofwf(i,j) = sumofwf(i,j)+
     &      		   (evapn(i,j,noc)*fractn(i,j,noc)+
     &      		   evapn(i,j,nse)*fractn(i,j,nse))/fracto(i,j)-
     &      		   couprf(i,j)-coupsf(i,j)-sumro(i,j)
            sumisno(i,j) = sumisno(i,j)+ coupsf(i,j)-
     &      		   evapn(i,j,nse)*fractn(i,j,nse)/fracto(i,j)
            sumicof(i,j) = sumicof(i,j)+hficof(i,j)
          endif
c-For the iceberg coupling
          sumuv10(i,j) = sumuv10(i,j)+uv10(i,j)
          sumpress(i,j)= sumpress(i,j)+pground(i,j) 
c-icb
          winstua(i,j) = winstua(i,j) + winstu(i,j)
          winstva(i,j) = winstva(i,j) + winstv(i,j)
	enddo
      enddo
      sumohsn      = sumohsn+sumhsn
      sumohss      = sumohss+sumhss
      
      if (iobclint.eq.nbclins) then
        ratm=1d0/float(nbclins)
	do j=1,nlon
	  do i=1,nlat
            sumohfx(i,j) = sumohfx(i,j)*ratm 
            sumoswr(i,j) = sumoswr(i,j)*ratm
            sumihfx(i,j) = sumihfx(i,j)*ratm 
            sumiswr(i,j) = sumiswr(i,j)*ratm
            sumofwf(i,j) = sumofwf(i,j)*ratm 
            sumisno(i,j) = sumisno(i,j)*ratm 
            sumicof(i,j) = sumicof(i,j)*ratm
            sumuv10(i,j) = sumuv10(i,j)*ratm
            sumpress(i,j)=sumpress(i,j)*ratm   
	  enddo
        enddo
C       sumohsn=sumohsn*ratm
C       sumohss=sumohss*ratm
        sumohsn=sumohsn
        sumohss=sumohss
	iobclint=0
	if (iclimflux.eq.1) call ec_climflux(ist,jst)
      endif
      
      if (iobtropt.eq.nbtrops) then
        rndws=1d0/dble(ndayws)
        dwsm1=ndayws-1d0
        ratm=1d0/float(nbtrops)
	do j=1,nlon
	  do i=1,nlat
	    winstua(i,j) = winstua(i,j)*ratm
	    winstva(i,j) = winstva(i,j)*ratm
            sumtx(i,j)=(dwsm1*sumtx(i,j)+winstua(i,j))*rndws
            sumty(i,j)=(dwsm1*sumty(i,j)+winstva(i,j))*rndws
	  enddo
        enddo
	iobtropt=0
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
      subroutine ec_climflux(ist,jst)
c-----------------------------------------------------------------------
c *** this routine calculates climatogical SST's and heatfluxes
c *** from year 1 onwards and outputs these to a file
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comcoup.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comsurf.h'

      integer i,j,k,l,index,istep,ist,jst
      
      real*4 hulph(nlat,nlon),hulpt(nlat,nlon)
      real*8 hefxcl(nlat,nlon,360),tmixcl(nlat,nlon,360),ryear,riatm
      real*8 hefx(nlat,nlon),hefxo(nlat,nlon)
      
      common /mixlayer/ hefxcl,tmixcl,hefx,hefxo
      
      
      if ( ist.eq.iobclin) then
        do k=1,360
          do j=1,nlon
            do i=1,nlat
              hefxcl(i,j,k)=0d0
              tmixcl(i,j,k)=0d0
            enddo
          enddo
	enddo
      endif
      
      index=(imonth-1)*30+iday
            
      do i=1,nlat
        do j=1,nlon
          hefxcl(i,j,index)=hefxcl(i,j,index)+sumohfx(i,j)+sumoswr(i,j)
          tmixcl(i,j,index)=tmixcl(i,j,index)+tsurfn(i,j,noc)
        enddo
      enddo

      if (index.eq.360.and.iyear.eq.nyears) then
	ryear=1/real(nyears)
        do k=1,360
	  do i=1,nlat
            do j=1,nlon
              if (fracto(i,j).gt.epss) then
                hefxcl(i,j,k)=hefxcl(i,j,k)*ryear
                tmixcl(i,j,k)=tmixcl(i,j,k)*ryear
              else
                hefxcl(i,j,k)=undef
                tmixcl(i,j,k)=undef
              endif
            enddo
	  enddo
	enddo
	index=0
	do k=1,360/iobclin
	  do l=1,iobclin
	    index=index+1
            write(51) ((real(hefxcl(i,j,index)),j=1,nlon),i=1,nlat)
            write(52) ((real(tmixcl(i,j,index)),j=1,nlon),i=1,nlat)
	  enddo
	enddo
      endif

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
     &winstua,winstva,sumtx,sumty
      write(iuo+95) dum
      return
      end

