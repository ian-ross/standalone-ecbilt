












      subroutine forcat(iyear,xjour)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c-This routine calculates the heat fluxes from atmospheric data
c read by forcng.
c Derniere modification 04/06/95
c Differences entre forcat.f et formel.f tcn(i,j,1) - scal(i,j,ks2,1)
c et le runoff, oceacn.com
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
C     include 'oceacn.com'
c
      dimension ws(96),zmue(96),zalcnp(96)
c     write(iuo+66,*) 'debut de forcat.f'
      data iflgaa /0/
      data iflgab /0/
      data iflago /0/
      data nintsr /24/
c
      zeps  = 1d-20
      zeps0 = 1.d-13
      zeps1 = 1.d-06
      dpi   = 2.*pi
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1) READ ATMOSPHERIC FORCING.                                        |
c-----------------------------------------------------------------------
c
          call forcng(iyear,xjour)
c
	    do jj=jcl1,jcl2
	       tenagx(ims1-1,jj) = tenagx(ims2,jj)    
	       tenagx(ims2+1,jj) = tenagx(ims1,jj)
	       tenagy(ims1-1,jj) = tenagy(ims2,jj)
	       tenagy(ims2+1,jj) = tenagy(ims1,jj)
            enddo
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2) COMPUTATION OF SNOW PRECIPITATION                                |
c-----------------------------------------------------------------------
c
c--2.1. FRACTION OF SNOW PRECIPITATION (LEDLEY, 1985).
c-----------------------------------------------------
c
	  do j=js1,js2
	    do i=is1(j),is2(j)
	      zind1       = max(zero,sign(one,253.0-tabq(i,j)))
	      zind2       = max(zero,sign(one,272.0-tabq(i,j)))
	      zind3       = max(zero,sign(one,281.0-tabq(i,j)))
	      hnplbq(i,j) = (zind1+(1.0-zind1)*
     &                      (zind2*(0.5+(272.0-tabq(i,j))/38.0)+
     &		            (1.0-zind2)
     &                       *zind3*((281.0-tabq(i,j))/18.0)))*
     &                      fwat(i,j)
 	      fwat(i,j)=fwat(i,j)+runoff(i,j)*86400*1000
	    enddo
	  enddo
c
c--2.1. FRACTION OF SNOW PRECIPITATION (ROSS & WALSH, 1987).
c-----------------------------------------------------------
c
c	  do j=js1,js2
c	    do i=is1(j),is2(j)
c	      if (tabq(i,j).lt.268.15) then
c		hnplbq(i,j) = fwat(i,j)
c	      else 
c		if (tabq(i,j).lt.276.15) then
c		  hnplbq(i,j) = ((276.15-tabq(i,j))/8.0)*fwat(i,j)
c	        else 
c		  hnplbq(i,j) = 0.0
c		endif
c	      endif
c	    enddo
c	  enddo
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3) COMPUTATION OF SOLAR IRRADIANCE.                                 |
c-----------------------------------------------------------------------
c
          indaet = 1
	  ijour  = int(xjour)
          dday = xjour*2.*pi/yeaday
          dcor=1.000d0+0.0013d0*sin(dday)+0.0342d0*cos(dday)
C         write(96,*) dcor,dday,xjour,pi,yeaday
C         write(96,*) dcor

C	  write(iuo+66,*) 'dcor: ',dcor,' jour : ', ijour
	  dec  = pdecli(indaet,ijour)
	  dec  = dec*radian
	  sdec = sin(dec)
	  cdec = cos(dec)
          do j=js1,js2
	    do i=is1(j),is2(j)
	    slat   = covrai(i,j)
	    zps    = slat*sdec
	    zpc    = cos(asin(slat))*cdec
	    zljour = acos(-sign(one,zps)
     &               *min(one,sign(one,zps)*(zps/zpc)))
	    dws    = (2.0*zljour)/real(nintsr)
	    zlmidi = asin(( zps +zpc ))/radian
	    zalcnq = 0.0
	    do k=1,nintsr
	      ws(k)     = zljour-(real(k)-0.5)*dws
	      zmue(k)   = max(zero,zps+zpc*cos(ws(k)))
	      zalcnp(k) = 0.05/(1.1*zmue(k)**1.4+0.15)
	      zalcnq    = zalcnq+zalcnp(k)*dws
	    enddo
	    zalcnq = zalcnq/max(2.0*zljour,zeps0)
	    zmudum = 0.4
c
c
              ih=max0(0,isign(1,j-jeq))
c
c--ALBEDO.
c
      	      call shine(ih,zmudum,tfsn,tfsg,ts(i,j),hgbq(i,j)
     &          	,hnbq(i,j),zalb,zalcn,zalbp,zaldum)
c
c--SHORTWAVE RADIATION.
c
C             evg=611.0*10.0**(9.5*(tdew(i,j)-273.16)
C    &            /(tdew(i,j)-7.660))
C             evo=611.0*10.0**(7.5*(tdew(i,j)-273.16)
C    &            /(tdew(i,j)-35.86))
C             es1=611.0*10.0**(7.5*(tabq(i,j)-273.16)
C    &            /(tabq(i,j)-7.660))
C             es2=611.0*10.0**(7.5*(tabq(i,j)-273.16)
C    &           /(tabq(i,j)-35.86))
C             evg=tdew(i,j)*es1
C             evo=tdew(i,j)*es2
              es =611.0*exp(min(sign(17.269*one,
     &            tabq(i,j)-too),sign(21.875*one,
     &            tabq(i,j)-too))*abs(tabq(i,j)-too)/
     &            (tabq(i,j)-35.86+max(zero,sign(28.2*one,
     &            too-tabq(i,j)))))
              evg=tdew(i,j)*es
              evo=tdew(i,j)*es
              e  =tdew(i,j)*es
              qabq(i,j) = (0.622*e)/(psbq(i,j)-(1.0-0.622)*e)
c
c--ZILLMAN.
c
 	      albg =  (1.0-cloud(i,j))*zalbp+cloud(i,j)*zalb
              frsdtg = 0.0
              frsdto = 0.0
              do k=1,nintsr
                albo   = (1.0-cloud(i,j))*zalcnp(k)+cloud(i,j)*zalcn
C               frsdtg = frsdtg+dws*(1.0-albg)*
C    &                   (1368.0*zmue(k)*zmue(k))/
C    &	                 ((zmue(k)+2.7)*evg*1.0e-05+1.085*zmue(k)+0.10)
 	        frsdto = frsdto+dws*(1.0-albo)*
     &                   (1368.0*zmue(k)*zmue(k))/
     &	                 ((zmue(k)+2.7)*evo*1.0e-05+1.085*zmue(k)+0.10)
 	      enddo
c
c--Computation of the solar heat flux absorbed at
c  the snow/ice and ocean surfaces.
c
C             fsolg(i,j) =(1.0-0.6*cloud(i,j)
C    &                     *cloud(i,j)*cloud(i,j))*frsdtg/dpi
C             fsolcn(i,j)=(1.0-0.6*cloud(i,j)
C    &                    *cloud(i,j)*cloud(i,j))*frsdto/dpi
C             fsolg(i,j) = 0.9*min(one,(1-.62*cloud(i,j)+.0019*zlmidi))*
C    &                     frsdtg/dpi
              fsolcn(i,j)= 0.9*min(one,(1-.62*cloud(i,j)+.0019*zlmidi))*
     &                     frsdto/dpi
 
c--Computation of the effective sea-water albedo.
c
C             albege(i,j) = 1.0-fsolg(i,j)/
C    &                      max(fsolg(i,j)/(1.0-albg),zeps0)
              albecn(i,j) = 1.0-fsolcn(i,j)/
     &                      max(fsolcn(i,j)/(1.0-((1.0-cloud(i,j))
     &                      *zalcnq+cloud(i,j)*zalcn)),zeps0)
c
c--SHINE AND CRANE.
c
              frsdrg = 0.0
              frsdfg = 0.0
              frsdro = 0.0
              frsdfo = 0.0
              do k=1,nintsr
 	        frsdrg = frsdrg+dws*
     &	                 (1368.0*zmue(k)*zmue(k)*(1.0-zalbp))/
     &	                 (1.2*zmue(k)+(1.0+zmue(k))*evg*1.0e-05+0.0455)
                frsdfg  = frsdfg+dws*
     &	                 ((53.5+1274.5*zmue(k))*sqrt(zmue(k))
     &                   *(1.0-0.996*zalb))/
     &                   (1.0+0.139*(1.0-0.9435*zalb)*tauc(i,j))
C               frsdro = frsdro+dws*
C    &	                 (1368.0*zmue(k)*zmue(k)*(1.0-zalcnp(k)))/
C    &	                 (1.2*zmue(k)+(1.0+zmue(k))*evo*1.0e-05+0.0455)
C               frsdfo = frsdfo+dws*
C    &	                 ((53.5+1274.5*zmue(k))*sqrt(zmue(k))
C    &                   *(1.0-0.996*zalcn))/
C    &                   (1.0+0.139*(1.0-0.9435*zalcn)*tauc(i,j))
              enddo
c
c--Computation of the solar heat flux absorbed at
c  the snow/ice and ocean surfaces.
C
              zinda       = 1.0-max(zero,sign(one,-(-0.5-albq(i,j))))
              fsolg(i,j)  = ((1.0-cloud(i,j))*frsdrg
     &                      +cloud(i,j)*frsdfg)/dpi
              fsolcn(i,j) = (((1.0-cloud(i,j))*frsdro
     &                      +cloud(i,j)*frsdfo)/dpi)*zinda
     &                      +(1.-zinda)*fsolcn(i,j)
 
c--Computation of the effective snow/ice albedo.
c
              albege(i,j) = (((frsdrg/(1.0-zalbp))*zalbp+
     &          	    (frsdfg/(1.0-zalb))*zalb)/
     &                      max(frsdrg/(1.0-zalbp)
     &                      +frsdfg/(1.0-zalb),zeps0))
     &                      *(1-max(zero,sign(one,-hgbq(i,j))))

c--Taking into account the ellipsity of the earth orbit
              fsolg(i,j)  =fsolg(i,j)*dcor
              fsolcn(i,j) =fsolcn(i,j)*dcor
	      ii1=85
	      ii2=55
	      jj1=2
	      jj2=2

c
c--Computation of the effective sea-water albedo.
c
C             albecn(i,j) = ((frsdro/(1.0-zalcnq))*zalcnq+
C    &    		    (frsdfo/(1.0-zalcn))*zalcn)/
C    &                      max(frsdro/(1.0-zalcnq)
C    &                      +frsdfo/(1.0-zalcn),zeps0)
c
	    enddo
	  enddo
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  4) PARTIAL COMPUTATION OF HEAT, WATER AND MOMENTUM FLUXES           |
c-----------------------------------------------------------------------
c
          do 170 j=js1,js2
            do 160 i=is1(j),is2(j)
c
c rhoa: Air density, obtained from the equation of
c	state for dry air.
c
              rhoa = psbq(i,j)/(287.04*tabq(i,j))
c
c--4.1. Calculate turbulent heat fluxes over water.
c--------------------------------------------------
c
cNCAR: bulk sensible and latent heat fluxes.
c
c  es: Saturation vapor pressure.
c  e : Vapor pressure.
c
 	      es         =  611.0*10.0**(7.5*(scal(i,j,ks2,1)
     &                      -273.16)/(scal(i,j,ks2,1)-35.86))
C	      es         =  611.0*10.0**(7.5*(tcn(i,j,1)
C    &                      -273.16)/(tcn(i,j,1)-35.86))
C             es2        =  611.0*10.0**(7.5*(tabq(i,j)
C    &                      -273.16)/(tabq(i,j)-35.86))
C	      e          =  611.0*10.0**(7.5*(tdew(i,j)
C    &                      -273.16)/(tdew(i,j)-35.86))
C             e          = tdew(i,j)*es2
	      zqsat      = (0.622*es)/(psbq(i,j)-(1.0-0.622)*es)
c             qabq(i,j)  = (0.622*e)/(psbq(i,j)-(1.0-0.622)*e)
c

c--computation of drag coefficients following Bunker
c
c	      tvirt      = tabq(i,j)*(1+.608*qabq(i,j))
c	      tvirtsu    = scal(i,j,ks2,1)*(1+.608*zqsat)
c	      ce         = ceb(vabq(i,j),tabq(i,j)-scal(i,j,ks2,1))
c	      ceorg      = ceb(vabq(i,j),tvirt-tvirtsu)
c	      ce         = .92*ceorg
c	      ch         = .87*ceorg

C	      ce         = ceb(vabq(i,j),tabq(i,j)-tcn(i,j,1))
C	      ce         = 1.75e-03
C             ch         = ce

c--drag coefficients from Large and Pond
      
c  Stability parameters
              dteta  = scal(i,j,ks2,1)-tabq(i,j)
C             dteta  = tcn(i,j,1)-tabq(i,j)
              deltaq = qabq(i,j)-zqsat
              ztvmoy = tabq(i,j)*(1.+2.2E-3*tabq(i,j)*qabq(i,j))
              Obouks = -70.0*10.*(dteta+3.2E-3*ztvmoy*ztvmoy*deltaq)/
     &                 (vabq(i,j)*vabq(i,j)*ztvmoy)
	      Obouks =max(zero,Obouks)
              psims  = -7.0*Obouks
              psihs  =  psims
              psils  =  psims
              Obouki =  -100.0*10.0*(dteta
     &                    +2.2E-3*ztvmoy*ztvmoy*deltaq)/
     &                 (vabq(i,j)*vabq(i,j)*ztvmoy)
	      Obouki = min(zero,Obouki)
              xins   = (1-16.*Obouki)**.25
              psimi  = 2*log((1+xins)/2)+
     &                 log((1+xins**2)/2)-2*atan(xins)+pi/2
              psihi  = 2*log((1+xins**2)/2)
              psili  = psihi
              stab   = max(zero,sign(one,dteta))
              psim   = stab*psimi+(1.0-stab)*psims
              psih   = stab*psihi+(1.0-stab)*psihs
              psil   = stab*psili+(1.0-stab)*psils

c  computation of intermediate values
              zzero  = .032*1.5E-3*vabq(i,j)*vabq(i,j)/gpes
              cmn    = (vkarmn/log(10./zzero))**2
              chn    = 0.0327*vkarmn/log(10./zzero)
              cln    = 0.0346*vkarmn/log(10./zzero)
              cmcmn  = 1/(1-sqrt(cmn)*psim/vkarmn)

c  ch and ce
              ch     = chn*cmcmn/(1-chn*psih/(vkarmn*sqrt(cmn)))
              ce     = cln*cmcmn/(1-cln*psil/(vkarmn*sqrt(cmn)))
c
cEnd compuation of ch and ce
c
	      drghce     = rhoa*ce*vabq(i,j)
	      drghch     = rhoa*ch*vabq(i,j)
C	      fcscn(i,j) = drghch*1004.0*(tcn(i,j,1)-tabq(i,j))
 	      fcscn(i,j) = drghch*1004.0*(scal(i,j,ks2,1)-tabq(i,j))
	      flecn(i,j) = drghce*2.5e+06*(zqsat-qabq(i,j))
c
c
c--4.3. CALCULATE EVAPORATION OVER WATER.
c----------------------------------------
c
              fevabq(i,j) = flecn(i,j)/cevap
c
160         continue
170       continue
c
c
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine forcat -
      end
c
