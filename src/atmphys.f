c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iatmphys
c-----------------------------------------------------------------------
c *** initializes variables used in subroutines of atmphys.f
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j,k,l
      real*8  albtune

c
c *** tuneable albedo at top of the atmosphere
c
      albtune=1.05


c
c *** initial atmospheric temperatures and moisture
c
      if (irunlabel.eq.0) then
        tempm(1)=tzero-35d0
        tempm(2)=tzero-8d0
        do j=1,nlon
          do i=1,nlat
            temp2g(i,j)=tempm(1)
            temp4g(i,j)=tempm(2)
            tsurf (i,j)=290d0
            tempsg(i,j)=290d0
            rmoisg(i,j)=0d0
            relhum(i,j)=0d0
            q10(i,j)=0d0
            qsurf(i,j)=0d0
            pground(i,j)=p0
            dp2(i,j)=dp1
            geopg(i,j,2)=0d0
          enddo
        enddo
        call ptground
        do j=1,nlon
          do i=1,nlat
            tsurf (i,j)= tempsg(i,j)
          enddo
        enddo
        call atmphyszero
      else
        read(90) tsurf,tempm
        read(90) rmoisg,torain
      endif

c
c *** longwave radiation parameters (held and suarez)
c
      do i=1,21
        read(9,310) pteta(i),pa(1,i),pb(1,i),pc(1,i),pa(2,i),pb(2,i),
     *              pc(2,i),pa(3,i),pb(3,i)
        pteta(i)=pteta(i) + tzero
        pc(3,i)=0.d0
      enddo

c *** calculation of interpolation factors for the Held and Suarez
c *** parameterization

      do l=1,3
        do i=2,21
          pafac(l,i)=(pa(l,i)-pa(l,i-1))/(pteta(i)-pteta(i-1))
          pbfac(l,i)=(pb(l,i)-pb(l,i-1))/(pteta(i)-pteta(i-1))
          pcfac(l,i)=(pc(l,i)-pc(l,i-1))/(pteta(i)-pteta(i-1))
        enddo
      enddo

c
c *** shortwave radiation parameters
c
      do i=1,nlat
        read(4,330) albeaw(i),albeas(i)
        albeaw(i)=albtune*albeaw(i)
        albeas(i)=albtune*albeas(i)
      enddo
      do i=1,nlat
        read(4,340) albland(i)
      enddo
      do i=1,nlat
       read(4,340) albsnow(i)
      enddo
      do i=1,nlat
        albice(i)=albsnow(i)
      enddo
      do i=1,nlat
        read(4,330) albseaw(i),albseas(i)
      enddo
      do i=1,nlat
       read(4,330) abstow(i),abstos(i)
      enddo

c
c *** evaporation factor
c
      evfac=1d0

      call detqmtabel

310   format(9f7.2)
330   format(2f5.2)
340   format(f5.2)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmphyszero
c-----------------------------------------------------------------------
c *** initializes data arrays to zero for physics routines
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      integer i,j

      do j=1,nlon
        do i=1,nlat
          torain(i,j)=0.d0
          dyrain(i,j)=0.d0
          corain(i,j)=0.d0
          thforg1(i,j)=0.d0
          thforg2(i,j)=0.d0
          vhforg1(i,j)=0.d0
          vhforg2(i,j)=0.d0
        enddo
      enddo

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine radiation
c-----------------------------------------------------------------------
c *** computes albedos and short and long wave radiational fluxes
c-----------------------------------------------------------------------
      implicit none

      call solar
      call albedo
      call swaverad
      call lwaverad

      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sensrad
c-----------------------------------------------------------------------
c *** computes atmospheric forcing due to sensible heat and radiation
c *** the forcing terms are computed in dimensional units (K/s)
c *** input
c ***       hflux : sensible heat flux between atmosphere and earth
c ***       hesws : short wave solar radiation in layer 1 or 2
c ***       ulrad : upward long wave radiation in layer 1 or 2
c ***       ulrads: upward long wave radiation in lower layer (2)
c ***       dlrads: downward long wave radiation in lower layer (2)
c *** output
c ***       thforg: temperature forcing in layer 1 or 2 in K/s
c ***       vhforg: diabatic forcing in layer 1 or 2 in K/s
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j
      real*8  halpha,halpha1,halpha2(nlat,nlon),sum1,sum2

      halpha=grav/cpair

      halpha1 =halpha/dp1

      do j=1,nlon
        do i=1,nlat
          halpha2(i,j)=halpha/dp2(i,j)
        enddo
      enddo

c
c *** summation of forcing terms
c
      do j=1,nlon
        do i=1,nlat

          sum1 = (hesw1(i,j) - ulrad1(i,j) + ulrad2(i,j))*halpha1
          sum2 = (hflux(i,j) + hesw2(i,j) - ulrad2(i,j) +
     *           ulrads(i,j) - dlrads(i,j))*halpha1

          thforg1(i,j) = thforg1(i,j) + sum1
          thforg2(i,j) = thforg2(i,j) + sum2

          vhforg1(i,j) = vhforg1(i,j) + sum1
          vhforg2(i,j) = vhforg2(i,j) + sum2

        enddo
      enddo

      return
      end



c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine solar
c-----------------------------------------------------------------------
c *** calculates incoming solar radiation as a function
c *** of the day of the year
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i
      real*8  deltal
      real*8  rkosz(nlat),rkosha1(nlat),ha1(nlat),ha2(nlat),qi(nlat)

      do i=1,nlat
        deltal=23.44d0*(pi/180.d0)*cos((180.d0-day)*pi/180.d0)
        rkosz(i)=sinfi(i)*sin(deltal)-cosfi(i)*cos(deltal)
        if (rkosz(i).lt.0.) then
          rkosha1(i)=-tanfi(i)*tan(deltal)
          if ((rkosha1(i).le.1).and.(rkosha1(i).ge.-1)) then
            ha1(i)=acos(rkosha1(i))
            ha2(i)=-ha1(i)
            rkosz(i)=sinfi(i)*sin(deltal)*(ha1(i)-ha2(i))/(2d0*pi)+
     *               cosfi(i)*cos(deltal)*(sin(ha1(i))-
     *               sin(ha2(i)))/(2d0*pi)
          else
            ha1(i)=pi
            rkosz(i)=sinfi(i)*sin(deltal)
          endif
        else
          ha1(i)=pi
          rkosz(i)=sinfi(i)*sin(deltal)
        endif
        qi(i)=solarc*rkosz(i)
        if  (qi(i).le.0.) qi(i)=0.d0
        q0(i)=qi(i)
      enddo

      return
      end



c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine albedo
c-----------------------------------------------------------------------
c *** calculates albedos as a function of the time of the year, linearly
c ***      interpolating between winter and summer values
c *** albea is the albedo at the top of the atmosphere
c *** albsea is the albedo of the open sea
c *** albes is the albedo of the earth surface, its value depends on whether
c ***      it is a land or sea point, and on whether the grid point is ice
c ***      or snow covered or not
c *** abso1,2 are the absorption coefficients for the upper and lower
c ***      atmoshperic layer
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j,nli,nsi
      real*8  abstot,alblat,delalbcls(nlat),delalbclw(nlat)
      real*8  delalbcl,albealand,albeasea

c *** frank: increase albedo at the top of the atmosphere over sea
c *** to account for higher cloudiness over sea: reduce albedo over
c *** land accordingly to keep the zonal mean albedo at the original
c *** value: include a seasonal cycle : northern hemisphere only

c      delalbcls(1)=0.1
c      delalbcls(2)=0.1
c      delalbcls(3)=0.2
c      delalbcls(4)=0.3
c      delalbcls(5)=0.5
c      delalbcls(6)=0.5
c      delalbcls(7)=0.5
c      delalbcls(8)=0.1
c      delalbcls(9)=0.1
c      delalbcls(10)=0.1
c      delalbcls(11)=0.1
c      delalbcls(12)=0.1
c      delalbcls(13)=0.1
c      delalbcls(14)=0.0
c      delalbcls(15)=0.0
c      delalbcls(16)=0.0
c      delalbcls(17)=0.0
c      delalbcls(18)=0.0
c      delalbcls(19)=0.1
c      delalbcls(20)=0.1
c      delalbcls(21)=0.1
c      delalbcls(22)=0.1
c      delalbcls(23)=0.1
c      delalbcls(24)=0.2
c      delalbcls(25)=0.6
c      delalbcls(26)=0.5
c      delalbcls(27)=0.5
c      delalbcls(28)=0.4
c      delalbcls(29)=0.3
c      delalbcls(30)=0.3
c      delalbcls(31)=0.1
c      delalbcls(32)=0.1


c      delalbclw(1)=0.1
c      delalbclw(2)=0.1
c      delalbclw(3)=0.3
c      delalbclw(4)=0.3
c      delalbclw(5)=0.3
c      delalbclw(6)=0.3
c      delalbclw(7)=0.3
c      delalbclw(8)=0.3
c      delalbclw(9)=0.2
c      delalbclw(10)=0.1
c      delalbclw(11)=0.1
c      delalbclw(12)=0.1
c      delalbclw(13)=0.1
c      delalbclw(14)=0.0
c      delalbclw(15)=0.0
c      delalbclw(16)=0.0
c      delalbclw(17)=0.0
c      delalbclw(18)=0.0
c      delalbclw(19)=0.0
c      delalbclw(20)=0.1
c      delalbclw(21)=0.1
c      delalbclw(22)=0.1
c      delalbclw(23)=0.1
c      delalbclw(24)=0.1
c      delalbclw(25)=0.1
c      delalbclw(26)=0.1
c      delalbclw(27)=0.5
c      delalbclw(28)=0.5
c      delalbclw(29)=0.3
c      delalbclw(30)=0.2
c      delalbclw(31)=0.1
c      delalbclw(32)=0.1

      do i=1,nlat

        alblat=abs(180.d0-day)*albeaw(i)/180.d0 +
     *           albeas(i) - abs(180.d0-day)*albeas(i)/180.d0

        do j=1,nlon
          albea(i,j)=alblat
        enddo


c        if (i.gt.16) then
c          delalbcl=abs(180.d0-day)*delalbclw(i)/180.d0 +
c     *           delalbcls(i) - abs(180.d0-day)*delalbcls(i)/180.d0
c          nli=0
c          nsi=0
c          do j=1,nlon
c            if (lsmask(i,j).eq.1) then
c              nli=nli+1
c            else
c              nsi=nsi+1
c            endif
c          enddo
c          if (nli.gt.0) then
c            albealand=alblat*(1d0-nsi*(1d0+delalbcl))/dble(nli)
c            albeasea=alblat*(1d0+delalbcl)
c          else
c            albealand=alblat
c            albeasea=alblat
c          endif
c        else
c          albealand=alblat
c          albeasea = alblat
c        endif

c        do j=1,nlon
c          if (lsmask(i,j).eq.1) then
c            albea(i,j)=albeasea
c          else
c            albea(i,j)=albealand
c          endif
c        enddo

      enddo

      do i=1,nlat
        albsea(i)=abs(180.d0-day)*albseaw(i)/180.d0 +
     *            albseas(i) - abs(180.d0-day)*albseas(i)/180.d0
      enddo

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then
            if (lseaice(i,j).eq.1) then
               albes(i,j)=albice(i)
            else
               albes(i,j)=albsea(i)
            endif
          else
            if (landsnow(i,j).eq.1) then
               albes(i,j)=albsnow(i)
            else
               albes(i,j)=albland(i)
            endif
          endif
        enddo
      enddo

      do i=1,nlat
        abstot=abs(180.d0-day)*abstow(i)/180.d0 +
     *         abstos(i) - abs(180.d0-day)*abstos(i)/180.d0
        do j=1,nlon
          abso1(i,j)=0.29d0*abstot
          abso2(i,j)=0.71d0*abstot
        enddo
      enddo

      return
      end



c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine swaverad
c-----------------------------------------------------------------------
c *** computes short wave radiation
c *** very simple, similar to approach of sergio franchito
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
      real*8  help

      do j=1,nlon
        do i=1,nlat
          hesw1(i,j)=(1.d0-albea(i,j))*abso1(i,j)*q0(i)
          hesw2(i,j)=(1.d0-albea(i,j))*abso2(i,j)*q0(i)
          help=(1.d0-albea(i,j))*(1.d0-abso1(i,j)-abso2(i,j))
          hesws(i,j)=help*(1-albes(i,j))*q0(i)
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lwaverad
c-----------------------------------------------------------------------
c *** computes long wave radiation according to the parameterization of
c *** held and suarez (jas 1978) and Idso and Jackson (1969) for the
c *** downward longwave radiation at the surface
c ***
c *** parameters: nlat   = number of gridpoints in the meridional
c ***                      direction (32)
c ***             nlon   = number of gridpoints in the zonal
c ***                      direction (64)
c ***             p0     = surface pressure (= 100000 Pa)
c ***             tzero  = 273.15 K
c ***             sboltz = constant of stefan boltzman
c ***
c *** input : temp4g(nlat,nlon): temperature [K] at 650 hPa
c ***         temp2g(nlat,nlon): temperature [K] at 350 hPa
c ***         tempsg(nlat,nlon): temperature [K] at 10 mtr
c ***         tsurf(nlat,nlon) : temperature [K] at the ground
c ***         pteta(21)        : contains 30 values of the mean
c ***                            potential temperatures running from
c ***                            -30 C to 70 C in steps of 5 degrees
c ***         pa(21,3)         : regression coefficients of Held and
c ***         pb(21,3)           Suarez for each of the potential
c ***         pc(21,3)           temperatures in pteta(21). index 3 is
c ***                            not used but instead use if made of the
c ***                            parameterization of Idso and Jackson
c ***                            for the downward flux at the surface
c ***         clfrac(nlat)     : climatological zonal mean cloudfraction
c ***                            used in the parameterization of Idso
c ***                            and Jackson for the downward flux at
c ***                            the surface
c ***
c *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] in
c ***                             the top layer (above 500 hPa)
c ***          ulrad2(nlat,nlon): upward longwave radiation [Wm-2] in
c ***                             the lower layer (below 500 hPa)
c ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
c ***                             the surface
c ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
c ***                             the surface
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,l,k,ix0
      real*8  facpot1,facpot2,dif,help,p750,p250
      real*8  pot1(nlat,nlon),pot2(nlat,nlon),avpot(nlat,nlon)
      real*8  depot(nlat,nlon),detes(nlat,nlon),ra(3),rb(3),rc(3)

      p250=25000d0
      p750=75000d0

c *** computation of potential temperatures at 250 and 750 mb
c *** given temperatures at 650 hPa and 350 hPa
c *** pot1 is potential temperature at 250 hPa
c *** pot2 is potential temperature at 750 hPa

      facpot1=1.d0/((p250/p0)**rkappa)
      facpot2=1.d0/((p750/p0)**rkappa)

      call levtemp(pot1,p250)
      call levtemp(pot2,p750)

      do j=1,nlon
        do i=1,nlat
          pot1(i,j)=pot1(i,j)*facpot1
          pot2(i,j)=pot2(i,j)*facpot2
        enddo
      enddo
c
c *** computation of mean and deviation of potential temperature, and
c *** difference between surface temperature and temperature of the
c *** ground.
c
      do j=1,nlon
        do i=1,nlat
          avpot(i,j)=0.5d0*(pot1(i,j)+pot2(i,j))
          depot(i,j)=0.5d0*(pot1(i,j)-pot2(i,j))
          detes(i,j)=(tsurf(i,j)-tempsg(i,j))
        enddo
      enddo
c
c *** computation of coefficients ai,bi and ci of regression formulae
c *** of held and suarez using linear interpolation in between the
c *** potential temperature intervals
c
      do j=1,nlon
        do i=1,nlat

          ix0=int((avpot(i,j)-tzero+40.d0)/5.d0)

          if ((ix0.ge.2).and.(ix0.le.21)) then
             dif=avpot(i,j)-pteta(ix0-1)
             do l=1,3
                ra(l)=pa(l,ix0-1) + pafac(l,ix0)*dif
                rb(l)=pb(l,ix0-1) + pbfac(l,ix0)*dif
                rc(l)=pc(l,ix0-1) + pcfac(l,ix0)*dif
             enddo
          else if (ix0.le.1) then
c ***   the potential temperature is out of range (too low)
c ***   the longwave radiation is calculated for the lowest potential
c ***   temperature ( -30 C)
c             call error(117)
             do l=1,3
                ra(l)=pa(l,1)
                rb(l)=pb(l,1)
                rc(l)=pc(l,1)
             enddo
          else if (ix0.ge.22) then
c ***   the potential temperature is out of range (too high)
c ***   the longwave radiation is calculated for the highest potential
c ***   temperature ( 70 C)
c             call error(118)
             do l=1,3
                ra(l)=pa(l,21)
                rb(l)=pb(l,21)
                rc(l)=pc(l,21)
             enddo
          endif

c ***     parameterization of Held and Suarez for the upward fluxes

          ulrad1(i,j)=ra(1) + rb(1)*depot(i,j) + rc(1)*detes(i,j)
          ulrad2(i,j)=ra(2) + rb(2)*depot(i,j) + rc(2)*detes(i,j)

c ***     parameterization of Idso and Jackson for the downward flux

c          help=1.-0.261d0*exp(-7.77d-04*(tzero-tempsg(i,j))**2)
c          dlrads(i,j)=sboltz*(tempsg(i,j)**4)*help*
c     &                (1.d0+0.275d0*clfrac(i))

c rein
c *** alternative paramterisation of dlrads (Ulden en Holtslag)
c
          dlrads(i,j)=9.35d-6*sboltz*(tempsg(i,j)**6) + 60.*clfrac(i)

c ***     black body radiation at the surface (does not as yet depend
c ***     on surface conditions)

          ulrads(i,j)=sboltz*tsurf(i,j)**4

        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine tracer
c-----------------------------------------------------------------------
c *** advection of tracer field
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j
      real*8  hdivmg(nlat,nlon)
      real*8  co2sp(nsh2)


      call rggtosp(co2,co2sp)
      call sptogg(co2sp,co2,pp)


c *** horizontal divergence of tracer

      call trafluxdiv(hdivmg,co2sp,co2)

c
c *** time stepping forward time stepping
c

      do j=1,nlon
        do i=1,nlat
          co2(i,j)=co2(i,j)-dtime*(hdivmg(i,j))
          if (co2(i,j).lt.0d0) co2(i,j)=0d0
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisfields
c-----------------------------------------------------------------------
c *** calculates relative humidity of the moised layer
c *** and specific humidity above the surface and at the surface
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      integer i,j
      real*8  qsat,pmount,tmount,qmax,dqmdt

      do j=1,nlon
        do i=1,nlat

          call ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)

c          dp2(i,j)=pmount-plevel(2)

          if (qmax.gt.0d0) then

            relhum(i,j)=rmoisg(i,j)/qmax

          else

            relhum(i,j)=1d0

          endif

          q10(i,j)=relhum(i,j) * qsat(pground(i,j),tempsg(i,j))

          qsurf(i,j) =qsat(pground(i,j),tsurf(i,j))

        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine detqmtabel
c-----------------------------------------------------------------------
c***  calculate tabel of maximum content of water in the layer below
c***  500 hPa using the Clausius-Clapeyron relation and assuming a
c***  temperature profile linear in log(p): T(p)= Tr + alpha*log(p/pr)
c***  where alpha is given by (T350-T650)/log(350/650)
c***  for given groundtemperature and 650 and 350 hPa temperatures
c***  the maximum water content in [m] is given by qmtabel
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      integer i,j,k
      real*8  tmount,t500,b,qmax,expint,z1,z2,bz1,bz2,hulpx
      real*8  t350,t650,rlogp500,alpha
      real*8  detqmax,detqmaxexact
      real*4  hulp(0:iqmtab,0:jqmtab,0:kqmtab)

      rlogp500=log(500d0/650d0)
      b=cc2*cc3-cc2*tzero

c      call system('rm -f outputdata/atmos/qmtabel.dat')

c      open(88,file='qmtabel.dat')
c      open(89,file='qmtabel.test')
c     *  form='unformatted',access='direct',recl=51*21*21)

      do i=0,iqmtab
        tqmi(i)=tqmimin + i*dtqmi
      enddo
      do j=0,jqmtab
        tqmj(j)=tqmjmin + j*dtqmj
      enddo
      do k=0,kqmtab
        tqmk(k)=tqmkmin + k*dtqmk
      enddo

      hulpx=cc1*exp(cc2)/(rowat*grav)

      do i=0,iqmtab
        t650=tqmi(i)

        do j=0,jqmtab
          tmount=tqmj(j)+t650

          do k=0,kqmtab
            t350=t650-tqmk(k)

            alpha=(t350-t650)*rlogtl12

            t500=t650+alpha*rlogp500

            z1=1/(tmount-cc3)
            z2=1/(t500-cc3)

            bz1=b*z1
            bz2=b*z2

            qmax=(exp(bz1)+expint(1,-bz1)*bz1)/z1 -
     #           (exp(bz2)+expint(1,-bz2)*bz2)/z2

            qmax=qmax*hulpx/alpha

            if (qmax.lt.0d0) qmax=0d0
            qmtabel(i,j,k)=qmax
c            write(88,111) tqmi(i),tqmj(j),tqmk(k),qmtabel(i,j,k)
          enddo
        enddo
      enddo

c      write(89) (((qmtabel(i,j,k),i=0,iqmtab),j=0,jqmtab),k=0,kqmtab)


c      do i=0,iqmtab
c        temp4g(1,1)=tqmi(i)
c        do j=0,jqmtab
c          tmount=tqmj(j)+temp4g(1,1)
c          do k=0,kqmtab
c            temp2g(1,1)=temp4g(1,1)-tqmk(k)
c            hulp(i,j,k)=detqmaxexact(tmount,1,1)
c            write(89,111) tqmi(i),tqmj(j),tqmk(k),hulp(i,j,k)
c          enddo
c        enddo
c      enddo


 111  format(4F14.7)
c      close(88)
c      close(89)
c      stop

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function detqmaxexact(tmount,i,j)
c-----------------------------------------------------------------------
c *** determines the maximum water content in latlon point
c *** i,j for given ground- and 650 and 350 hPa temperature
c *** by linear interpolation in qmtabel
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      integer i,j
      real*8  tmount,alpha,t500,z1,z2,bz1,bz2,b,hulpx
      real*8  qmax,detqmaxexact,expint

      b=cc2*cc3-cc2*tzero

      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12

      t500=temp4g(i,j)+alpha*log(500d0/650d0)

      z1=1/(tmount-cc3)
      z2=1/(t500-cc3)

      bz1=b*z1
      bz2=b*z2

      hulpx=cc1*exp(cc2)/(rowat*grav*alpha)

      qmax=hulpx*(exp(bz1)+expint(1,-bz1)*bz1)/z1 -
     #     hulpx*(exp(bz2)+expint(1,-bz2)*bz2)/z2

      if (qmax.lt.0d0) qmax=0d0

      if (qmax.gt.0.2) then
        write(29,*) 'in latlon ',i,j,' qmax ',qmax
        call error(121)
      endif

      detqmaxexact=qmax

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function detqmax(tmount,i,j,dqmdt)
c-----------------------------------------------------------------------
c *** determines the maximum water content in latlon point
c *** i,j for given ground- and 650 and 350 hPa temperature
c *** by linear interpolation in qmtabel
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'

      integer i,j,ii,jj,kk
      real*8  tmount,ti,tj,tk
      real*8  qmax,detqmax,dqmdi,dqmdj,dqmdk
      real*8  dqmdt,hmount,z500,t500,alpha,dtgdt

      ti=temp4g(i,j)
      tj=tmount-temp4g(i,j)
      tk=temp4g(i,j)-temp2g(i,j)

      if (ti.lt.tqmi(0)) then
        ti=tqmi(0)
c        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
c        call error(121)
      endif
      if (ti.gt.tqmi(iqmtab)) then
        ti=tqmi(iqmtab)
c        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
c        call error(121)
      endif

      if (tj.lt.tqmj(0)) then
        tj=tqmj(0)
c        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
c        call error(121)
      endif
      if (tj.gt.tqmj(jqmtab)) then
        tj=tqmj(jqmtab)
c        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
c        call error(121)
      endif

      if (tk.lt.tqmk(0)) then
        tk=tqmk(0)
c        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
c        call error(121)
      endif
      if (tk.gt.tqmk(kqmtab)) then
        tk=tqmk(kqmtab)
c        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
c        call error(121)
      endif

      ii=min(iqmtab-1,int((ti-tqmimin)*rdtqmi))
      jj=min(jqmtab-1,int((tj-tqmjmin)*rdtqmj))
      kk=min(kqmtab-1,int((tk-tqmkmin)*rdtqmk))

      dqmdi=(qmtabel(ii+1,jj,kk)-qmtabel(ii,jj,kk))*rdtqmi
      dqmdj=(qmtabel(ii,jj+1,kk)-qmtabel(ii,jj,kk))*rdtqmj
      dqmdk=(qmtabel(ii,jj,kk+1)-qmtabel(ii,jj,kk))*rdtqmk

      qmax = qmtabel(ii,jj,kk) + (ti-tqmi(ii))*dqmdi +
     *     (tj-tqmj(jj))*dqmdj + (tk-tqmk(kk))*dqmdk
      if (qmax.lt.0d0) qmax=0d0

      if (qmax.gt.0.2) then
        write(29,*) 'in latlon ',i,j,' qmax ',qmax
        call error(121)
      endif

      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12
      t500=temp4g(i,j)+alpha*alogpl2tl2
      z500=gpm500*grav
      hmount=qmount(i,j)*hmoisr*grav

      dtgdt=(rgas*t500*alogtl1pl2 + (hmount-geopg(i,j,2)-z500))/
     *      (rgas*tmount*alogtl12)

      dqmdt=dqmdi + dqmdj * (dtgdt - 1d0) + dqmdk

      detqmax=qmax

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function expint(n,x)
      implicit none
      integer n,maxit
      real*8 expint,x,eps,fpmin,euler
      parameter (maxit=100,eps=1.e-10,fpmin=1.e-30,euler=.5772156649)
      integer i,ii,nm1
      real*8 a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
         write (*,*) 'bad arguments in expint'
         stop
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0.)then
        expint=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./fpmin
        d=1./b
        h=d
        do 11 i=1,maxit
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.eps)then
            expint=h*exp(-x)
            return
          endif
11      continue
c        pause 'continued fraction failed in expint'
        call error(20)
      else
        if(nm1.ne.0)then
          expint=1./nm1
        else
          expint=-log(x)-euler
        endif
        fact=1.
        do 13 i=1,maxit
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-euler
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*eps) return
13      continue
c        pause 'series failed in expint'
        call error(20)
      endif
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisture
c-----------------------------------------------------------------------
c *** acvection and sources and sinks of atmospheric moisture
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j
      real*8  hdmoisg(nlat,nlon)
      real*8  fomoisg(nlat,nlon),qstar,hdivmg(nlat,nlon)
      real*8  levtempgp,factemv,factems,omegg500,t500,qsat

      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair

      call rggtosp(rmoisg,rmoiss)
      call sptogg(rmoiss,rmoisg,pp)

      do j=1,nlon
        do i=1,nlat
          if (rmoisg(i,j).lt.0.) then
            rmoisg(i,j)= 0d0
c            write(100,*) '1. in lat-lon ',j,i,' rmoisg lt 0'
          endif
        enddo
      enddo

c *** horizontal divergence of moisture

      call trafluxdiv(hdivmg,rmoiss,rmoisg)
c      call trafluxdivupw(hdivmg,rmoisg)

c *** vertical advection of moisture

      do j=1,nlon
        do i=1,nlat
          omegg500=(omegg(i,j,1)+omegg(i,j,2))/2.d0
          vemoisg(i,j)=0d0
          if (omegg500.lt.0.d0) then
            t500=levtempgp(plevel(2),i,j)
            qstar=relhum(i,j)*qsat(plevel(2),t500)
            vemoisg(i,j)=-omegg500*qstar/(grav*rowat)
            dyrain(i,j) = dyrain(i,j) + vemoisg(i,j)

            if (tsurf(i,j).ge.tzero) then
              thforg1(i,j)=thforg1(i,j) + factemv*vemoisg(i,j)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factemv*vemoisg(i,j)/dp1
            else
              thforg1(i,j)=thforg1(i,j) + factems*vemoisg(i,j)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factems*vemoisg(i,j)/dp1
            endif
          endif
        enddo
      enddo

c *** horizontal diffusion of moisture

      call hdiff(hdmoisg)

c
c *** time stepping forward time stepping
c

      do j=1,nlon
        do i=1,nlat
          rmoisg(i,j)=rmoisg(i,j)+dtime*(-ihavm*hdivmg(i,j)
     *      -ivavm*vemoisg(i,j) + hdmoisg(i,j) + imsink*evap(i,j))
          if (rmoisg(i,j).lt.0.) then
            rmoisg(i,j)= 0d0
c            write(100,*) '2. in lat-lon ',j,i,' rmoisg lt 0'
          endif
        enddo
      enddo

c      call trasmooth(rmoisg)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine trafluxdiv(tfdiv,ctrasp,ctra)
c-----------------------------------------------------------------------
c *** computes horizontal divergence of tracer flux
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k
      real*8  ctrasp(nsh2),vv(nsh2),ww(nsh2)
      real*8  dcdl(nlat,nlon),dcdm(nlat,nlon)
      real*8  tfdiv(nlat,nlon),ctra(nlat,nlon)

c *** 800 hPa winds are reduced with umoisr in the advection of the
c *** tracer field

c *** spatial derivatives of tracer

      call ddl (ctrasp,vv)
      call sptogg (vv,dcdl,pp)
      call sptogg (ctrasp,dcdm,pd)

c *** advection of tracer by total wind + convergence of tracer

      do j=1,nlon
        do i=1,nlat
          tfdiv(i,j)=dcdl(i,j)*(u800(i,j) + udivg(i,j,3))/
     *               (radius*cosfi(i)) +
     *               dcdm(i,j)*(v800(i,j) + vdivg(i,j,3))/
     *               (radius/cosfi(i)) +
     *               ctra(i,j)*divg(i,j,3)

          tfdiv(i,j)=tfdiv(i,j)*umoisr

        enddo
      enddo

      call rggtosp(tfdiv,vv)
      vv(1)=0d0
      call sptogg (vv,tfdiv,pp)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine trafluxdivupw(tfdiv,ctra)
c-----------------------------------------------------------------------
c *** computes horizontal divergence of tracer flux
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k
      real*8  dcdl(nlat,nlon),dcdf(nlat,nlon)
      real*8  tfdiv(nlat,nlon),ctra(nlat,nlon)

c *** 800 hPa winds are reduced with umoisr in the advection of the
c *** tracer field

c *** spatial derivatives of tracer

      call tragradupw(ctra,dcdl,dcdf)

c *** advection of tracer by total wind + convergence of tracer

      do j=1,nlon
        do i=1,nlat
          tfdiv(i,j)=dcdl(i,j)*(u800(i,j) + udivg(i,j,3))/
     *               (radius*cosfi(i)) +
     *               dcdf(i,j)*(v800(i,j) + vdivg(i,j,3))/
     *               (radius) +
     *               ctra(i,j)*divg(i,j,3)

          tfdiv(i,j)=tfdiv(i,j)*umoisr

        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine tragradupw(ctra,dcdl,dcdf)
c-----------------------------------------------------------------------
c *** computes horizontal gradient of tracer field
c *** using an upwind differencing scheme
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,jdum,jdlab
      real*8  ctra(nlat,nlon),dcdl(nlat,nlon),dcdf(nlat,nlon)
      real*8  rdyu(nlat),rdyd(nlat),ctranp,ctrasp,dlabj

      do i=2,nlat-1
       rdyu(i)=phi(i)-phi(i-1)
       rdyd(i)=phi(i+1)-phi(i)
      enddo

      rdyu(1)=(phi(1)+ pi/2.)
      rdyd(1)=phi(2)-phi(1)
      rdyu(nlat)=phi(nlat)-phi(nlat-1)
      rdyd(nlat)=(pi/2. - phi(nlat))

      do i=2,nlat-1

        do j=2,nlon-1
          if (utot(i,j,3).ge.0.) then
            dcdl(i,j)=(ctra(i,j)-ctra(i,j-1))/dlab
          else
            dcdl(i,j)=(ctra(i,j+1)-ctra(i,j))/dlab
          endif
        enddo

        if (utot(i,1,3).ge.0.) then
          dcdl(i,1)=(ctra(i,1)-ctra(i,nlon))/dlab
        else
          dcdl(i,1)=(ctra(i,2)-ctra(i,1))/dlab
        endif

        if (utot(i,nlon,3).ge.0.) then
          dcdl(i,nlon)=(ctra(i,nlon)-ctra(i,nlon-1))/dlab
        else
          dcdl(i,nlon)=(ctra(i,1)-ctra(i,nlon))/dlab
        endif

      enddo

c *** determine zonal gradients at the most southern and northern
c *** latitudes from differences jdlab gridpoints away in order
c *** to avoid violation of the cfl cricterium

      jdlab=2
      dlabj=jdlab*dlab

      do i=1,nlat,nlat-1

        do j=jdlab+1,nlon-jdlab
          if (utot(i,j,3).ge.0.) then
            dcdl(i,j)=(ctra(i,j)-ctra(i,j-jdlab))/dlabj
          else
            dcdl(i,j)=(ctra(i,j+jdlab)-ctra(i,j))/dlabj
          endif
        enddo

        do j=1,jdlab
          if (utot(i,1,3).ge.0.) then
            dcdl(i,j)=(ctra(i,j)-ctra(i,j-jdlab+nlon))/dlabj
          else
            dcdl(i,j)=(ctra(i,j+jdlab)-ctra(i,j))/dlabj
          endif
        enddo

        do j=nlon-jdlab+1,nlon
          if (utot(i,nlon,3).ge.0.) then
            dcdl(i,j)=(ctra(i,j)-ctra(i,j-jdlab))/dlabj
          else
            dcdl(i,j)=(ctra(i,j+jdlab-nlon)-ctra(i,j))/dlabj
          endif
        enddo
      enddo

      do j=1,nlon
        do i=2,nlat-1
          if (vtot(i,j,3).ge.0.) then
            dcdf(i,j)=(ctra(i,j)-ctra(i-1,j))/rdyu(i)
          else
            dcdf(i,j)=(ctra(i+1,j)-ctra(i,j))/rdyd(i)
          endif
        enddo
      enddo

      ctranp=0d0
      ctrasp=0d0

      do j=1,nlon
        ctrasp=ctrasp+ctra(1,j)
        ctranp=ctranp+ctra(nlat,j)
      enddo

      ctrasp=ctrasp/dble(nlon)
      ctranp=ctranp/dble(nlon)

      do j=1,nlon

        if (vtot(1,j,3).ge.0.) then
c          jdum=mod(j+nlon/2-1,nlon)+1
          dcdf(1,j)=(ctra(1,j)-ctrasp)/rdyu(1)
        else
          dcdf(1,j)=(ctra(2,j)-ctra(1,j))/rdyd(1)
        endif

        if (vtot(nlat,j,3).ge.0.) then
          dcdf(nlat,j)=(ctra(nlat,j)-ctra(nlat-1,j))/rdyu(nlat)
        else
c          jdum=mod(j+nlon/2-1,nlon)+1
          dcdf(nlat,j)=(ctranp-ctra(nlat,j))/rdyd(nlat)
        endif

      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine trasmooth(co2)
c-----------------------------------------------------------------------
c *** this routine smooths the tracer field co2 using a jmean points
c *** running mean in the zonal direction at the latitude circles
c *** closest to the poles
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'

      integer i,j,j1,jdel,jmean,ilat(4)
      real*8  clat(nlon),co2(nlat,nlon)

      jdel=5
      jmean=2*jdel+1

      ilat(1)=1
      ilat(2)=2
      ilat(3)=nlat-1
      ilat(4)=nlat

      do i=1,4
        do j=1,nlon
          clat(j)=co2(ilat(i),j)
        enddo

        do j=1,jdel
          co2(ilat(i),j)=0d0
          do j1=j-jdel+nlon,nlon
            co2(ilat(i),j)=co2(ilat(i),j)+clat(j1)
          enddo
          do j1=1,j+jdel
            co2(ilat(i),j)=co2(ilat(i),j)+clat(j1)
          enddo
          co2(ilat(i),j)=co2(ilat(i),j)/dble(jmean)
        enddo

        do j=jdel+1,nlon-jdel
          co2(ilat(i),j)=0d0
          do j1=j-jdel,j+jdel
            co2(ilat(i),j)=co2(ilat(i),j)+clat(j1)
          enddo
          co2(ilat(i),j)=co2(ilat(i),j)/dble(jmean)
        enddo

        do j=nlon-jdel+1,nlon
          co2(ilat(i),j)=0d0
          do j1=j-jdel,nlon
            co2(ilat(i),j)=co2(ilat(i),j)+clat(j1)
          enddo
          do j1=1,j+jdel-nlon
            co2(ilat(i),j)=co2(ilat(i),j)+clat(j1)
          enddo
          co2(ilat(i),j)=co2(ilat(i),j)/dble(jmean)
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine hdivspec(hduvg,ug,vg)
c-----------------------------------------------------------------------
c *** computes horizontal divergence
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
      real*8  hduvg(nlat,nlon),ug(nlat,nlon),vg(nlat,nlon)
      real*8  dugdl(nlat,nlon),dvgdm(nlat,nlon)
      real*8  usp(nsh2),vv(nsh2),vsp(nsh2)
      real*8  dx,dy(nlat)

      do j=1,nlon
        do i=1,nlat
          vg(i,j)=vg(i,j)*cosfi(i)
        enddo
      enddo

      call rggtosp(ug,usp)
      call rggtosp(vg,vsp)
      call ddl (usp,vv)
      call sptogg (vv,dugdl,pp)
      call sptogg (vsp,dvgdm,pd)

      do j=1,nlon
        do i=1,nlat
          hduvg(i,j)= (dugdl(i,j)/cosfi(i)+dvgdm(i,j))/radius
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine hdiff(hdmg)
c-----------------------------------------------------------------------
c *** horizontal diffusion of moisture
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer idifq,k
      real*8  hdmoiss(nsh2),hdmg(nlat,nlon)
      real*8  difq,rll

      rll=dble(ll(nsh))
      difq=max(0.d0,1.d0/(rll*(rll+1)*tdifq*24d0*3600d0))

      call lap(rmoiss,hdmoiss)

      hdmoiss(1)=0d0

      do k=2,nsh2
        hdmoiss(k)=difq*hdmoiss(k)
      enddo

      call sptogg(hdmoiss,hdmg,pp)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine convec
c-----------------------------------------------------------------------
c *** moist convective adjustment
c-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer ncmax,iconvn,i,j
      real*8  qsatcr,tsatcr,pref,t500,qsat500,pot2g,pot4g,dcmoisg
      real*8  fachulp,facteta,factemv,factems,pmount,tmount
      real*8  qmax,qsat,hulp,theat,hulp1,hulp2,redrain
      real*8  temp2go,temp4go,levtempgp
      real*8  detqmax,drainm,crainm,dqmdt

      fachulp=0.622d0*(rlatvap**2)/(cpair*rgas)
      facteta=0.6d0*rgas*(2.d0**rkappa)/grav
      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair
      ncmax=0
      crainm=0.5d0*rainmax
      drainm=rainmax-crainm


      do j=1,nlon
        do i=1,nlat
          iconvn=0
  10      continue

c ***     calculate pressure and temperature at the ground
c ***     and the maximum water content

            call ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)

c            dp2(i,j)=pmount-plevel(2)

c ***     relhmax defines the relative humidity at which oversaturation
c ***     occurs

            qmax=relhmax*qmax

            pot2g=temp2g(i,j)/potfac1
            pot4g=temp4g(i,j)/potfac2
            teta(i,j)=0.5d0*(pot2g-pot4g)
c ***       dry adiabatic lapse rate
            tetacr(i,j)=0d0
            dcmoisg=0d0

            if (rmoisg(i,j).gt.qmax) then

c ***     calculate rain reduction factor to account for increased
c ***     moisture capacity due to latent heat release

              redrain=1d0+dqmdt*relhmax*rowat*rlatvap*grav/(cpair*dp1)

              dcmoisg=(rmoisg(i,j)-qmax)/(redrain*dtime)

c ***         if the air is supersaturated initially, this is due to
c ***         large scale convergence of moisture and large scale
c ***         temperature changes. Excessive moisture
c ***         is then removed as dynamic rain

              if (iconvn.eq.0) then
                dyrain(i,j)=dyrain(i,j) + dcmoisg
                if (dyrain(i,j).gt.drainm) then
                  dcmoisg=drainm-dyrain(i,j)+dcmoisg
                  dyrain(i,j)=drainm
                endif
                if (tsurf(i,j).ge.tzero) then
                  thforg2(i,j)=thforg2(i,j) + factemv*dcmoisg/dp1
                  vhforg2(i,j)=vhforg2(i,j) + factemv*dcmoisg/dp1
                else
                  thforg2(i,j)=thforg2(i,j) + factems*dcmoisg/dp1
                  vhforg2(i,j)=vhforg2(i,j) + factems*dcmoisg/dp1
                endif
              else
                corain(i,j)=corain(i,j)+ dcmoisg
                if (corain(i,j).gt.crainm) then
                  dcmoisg=crainm-corain(i,j)+dcmoisg
                  corain(i,j)=crainm
                endif
                if (tsurf(i,j).ge.tzero) then
                  temp4g(i,j) =temp4g(i,j) + factemv*dcmoisg*dtime
     *                                       /dp1
                  vhforg2(i,j)=vhforg2(i,j)+ factemv*dcmoisg/dp1
                else
                  temp4g(i,j) =temp4g(i,j) + factems*dcmoisg*dtime
     *                                       /dp1
                  vhforg2(i,j)=vhforg2(i,j)+ factems*dcmoisg/dp1
                endif
              endif

c ***         moisture changes due to dynamic and convective rain

              rmoisg(i,j)=rmoisg(i,j)-ivavm*dcmoisg*dtime


c ***         calculate moist adiabatic lapse rate

              pot2g=temp2g(i,j)/potfac1
              pot4g=temp4g(i,j)/potfac2
              teta(i,j)=0.5d0*(pot2g-pot4g)

c              t500 = levtempgp(plevel(2),i,j)

              t500 = 0.5d0 * (temp2g(i,j) + temp4g(i,j))

              qsat500=qsat(plevel(2),t500)

              hulp=1d0 + fachulp*qsat500/(t500**2)
              gams(i,j)=gamd*(1+(rlatvap*qsat500)/(rgas*t500))/hulp
              tetacr(i,j)=0.5d0*t500*facteta*(gamd-gams(i,j))

            endif
            if (teta(i,j) .lt. (tetacr(i,j)-0.1)) then


c ***  frank: verdisconteer nog de varierende dikte dp2 van de onderste
c ***         laag (onderste laag is niet even dik als bovenste !!!)

c              hulp1=potfac1*dp1
c              hulp2=potfac2*dp2(i,j)

              hulp1=potfac1
              hulp2=potfac2

              theat = 0.5 * ( pot2g*hulp1 + pot4g*hulp2)

              pot2g=2.d0*(theat + hulp2*tetacr(i,j))/(hulp1+hulp2)
              pot4g=pot2g - 2.d0*tetacr(i,j)
              temp2go=temp2g(i,j)
              temp4go=temp4g(i,j)
              temp2g(i,j)=pot2g*potfac1
              temp4g(i,j)=pot4g*potfac2
              vhforg1(i,j)=vhforg1(i,j) + (temp2g(i,j)-temp2go)/dtime
              vhforg2(i,j)=vhforg2(i,j) + (temp4g(i,j)-temp4go)/dtime
              iconvn=iconvn+1
              if (iconvn.lt.10) then
                if (dyrain(i,j).eq.drainm) then
                  write(29,*) 'in latlon ',i,j, ' dyrain'
                  call error(122)
                  goto 20
                endif
                if (corain(i,j).eq.crainm) then
                  write(29,*) 'in latlon ',i,j, ' corain'
                  call error(122)
                  goto 20
                endif
                goto 10
              else
                write(29,*) 'warning in lat-lon point: ',i,j
                write(29,*) temp2g(i,j),temp4g(i,j)
                call error(120)
                goto 20
              endif
            else
              goto 20
            endif
  20      continue
          if (iconvn.gt.ncmax) ncmax=iconvn
          torain(i,j)= dyrain(i,j) + corain(i,j)
        enddo
      enddo
c      write(100,*) 'maximum iterations of convection : ',ncmax

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      function qsat(press,temp)
c-----------------------------------------------------------------------
c *** saturation mixing ratio
c *** input press in [Pa], temp in K
c *** output qsat: saturation mixing ratio
c-----------------------------------------------------------------------
      implicit none
      include 'comatm.h'
      include 'comphys.h'

      real*8  press,temp,qsat

      qsat=cc1*exp(cc2*(temp-tzero)/(temp-cc3))
     &     /press

      end


c23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine levtemp(tlev,plev)
c-----------------------------------------------------------------------
c *** computation temperatures at level p [Pa] assuming a constant
c *** temperature lapse rate : dt/dlnp = constant
c *** input: plev
c *** ouput: tlev
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
      real*8  tlev(nlat,nlon)
      real*8  r,plev

      r=log(plev/65000.d0)*rlogtl12
      do j=1,nlon
        do i=1,nlat
          tlev(i,j)=temp4g(i,j)+r*(temp2g(i,j)-temp4g(i,j))
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function levtempgp(plev,i,j)
c-----------------------------------------------------------------------
c *** computation temperatures at level p [Pa] assuming a constant
c *** temperature lapse rate : dt/dlnp = constant
c *** input: plev [Pa],i,j
c *** output: levtempgp [K]
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
      real*8  levtempgp
      real*8  r,plev

      r=log(plev/65000.d0)*rlogtl12
      levtempgp=temp4g(i,j)+r*(temp2g(i,j)-temp4g(i,j))

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fortemp
c-----------------------------------------------------------------------
c *** computes  new atmospheric temperatures due to heatforcing
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      integer i,j

      do j=1,nlon
        do i=1,nlat
          temp2g(i,j)=temp2g(i,j) + dtime*thforg1(i,j)
          temp4g(i,j)=temp4g(i,j) + dtime*thforg2(i,j)
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine ptground
c-----------------------------------------------------------------------
c *** computation of ground pressure and temperature
c *** assuming temperature varies linearly with log of the pressure
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
      real*8  t500(nlat,nlon),hulpx,levtempgp,hmount,hred
      real*8  alpha(nlat,nlon),pfac,hfac,z500

      z500=gpm500*grav
      hfac=2/rgas
      hred=1d0*grav
      pfac=log(plevel(2)/tlevel(2))

c *** calculate temperature at t500 assuming the temperature varies
c *** linearly with log(p) : T = Tref + alpha * log (p/pref)

      do j=1,nlon
        do i=1,nlat
          alpha(i,j)=(temp2g(i,j) - temp4g(i,j))*rlogtl12
          t500(i,j) = temp4g(i,j) + alpha(i,j)*pfac
        enddo
      enddo


      do j=1,nlon
        do i=1,nlat

c *** calculate ground height in decameters

          hmount=rmount(i,j)*hred

c *** calculate 10 mtr temperature in K assuming that the mean
c *** geopotential height at 500 hPa is gpm500 meter

          tempsg(i,j)=t500(i,j)*t500(i,j) -
     *                    hfac*alpha(i,j)*(hmount-geopg(i,j,2)-z500)
          tempsg(i,j)=max(0d0,tempsg(i,j))
          tempsg(i,j)=sqrt(tempsg(i,j))
        enddo
      enddo


      do j=1,nlon
        do i=1,nlat

          if (tempsg(i,j).eq.0d0) then
            write(29,*) 'in latlon ',i,j
            call error(18)
          endif

        enddo
      enddo

      do j=1,nlon
        do i=1,nlat

c *** calculate the groundpressure

          pground(i,j)=plevel(2)*exp((tempsg(i,j)-t500(i,j))/alpha(i,j))
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)
c-----------------------------------------------------------------------
c *** computation of ground pressure and temperature in order to
c *** to calculate the maximum precipitable water content in latlon i,j
c *** qmount contains the topography for this purpose
c *** assuming temperature varies linearly with log of the pressure
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j
      real*8  t500,levtempgp,hmount,hred,z500,dqmdt
      real*8  alpha,pfac,hfac,pmount,tmount,qmax,detqmax

      z500=gpm500*grav
      hfac=2/rgas
      hred=hmoisr*grav
      pfac=log(plevel(2)/tlevel(2))

c *** calculate temperature at t500 assuming the temperature varies
c *** linearly with log(p) : T = Tref + alpha * log (p/pref)

      alpha=(temp2g(i,j) - temp4g(i,j))*rlogtl12
      t500 =temp4g(i,j) + alpha*pfac

c *** calculate reduced ground height in decameters
c *** reduction occurs in order to tune the amount of moisture which
c *** is allowed to pass a topographic barier

      hmount=qmount(i,j)*hred
      if (hmount.lt.0d0) hmount=0d0

c *** calculate the groundpressure assuming that the mean geopotential
c *** height at 500 hPa is gpm500 decameter
c *** calculate 10 mtr temperature in K

      tmount=t500**2 - hfac*alpha*(hmount-geopg(i,j,2)-z500)
      if (tmount.lt.0) then
        write(29,*) 'in latlon ',i,j
        call error(18)
      else
        tmount=sqrt(tmount)
      endif

c      pmount=plevel(2)*exp((tmount-t500)/alpha)

      qmax=detqmax(tmount,i,j,dqmdt)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function globalmean(gfield)
c-----------------------------------------------------------------------
c *** computes global mean of gfield
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j
      real*8  gfield(nlat,nlon),sum(nlat),globsum,globfac,globalmean

      globfac=1d0/dsqrt(dble(nlon))

      do i=1,nlat
        sum(i)=0d0
      enddo

      do j=1,nlon
        do i=1,nlat
          sum(i)=sum(i)+gfield(i,j)
        enddo
      enddo

      globsum=0d0

      do i=1,nlat
        globsum=globsum+sum(i)*pw(i,1)
      enddo

      globalmean=globsum*globfac

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine meantemp
c-----------------------------------------------------------------------
c *** computes mean atmospheric temperatures
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'

      real*8  globalmean

c *** mean temperatures

      temp4gm=globalmean(temp4g)
      temp2gm=globalmean(temp2g)

      tempm(1)=temp2gm
      tempm(2)=temp4gm

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine dyntemp
c-----------------------------------------------------------------------
c *** computes temperature in K from streamfunction by solving the
c *** linear balance equation:
c *** del phi = (1 - mu**2 ) d psi/dmu + mu del psi
c *** the mean level is given by tempm
c *** input:  psi,tempm
c *** output: temp2g,temp4g,temp2gm,temp4gm,geopg(nlat,nlon,nvl)
c-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,l
      real*8  tempfac(ntl)
      real*8  geogt(nlat,nlon,ntl),tempgmt(ntl),tempgt(nlat,nlon,ntl)

      tempfac(1)=350.d0/(rgas*300.d0)
      tempfac(2)=650.d0/(rgas*300.d0)

      do j=1,nlon
        do i=1,nlat
          geogt(i,j,1)=geopg(i,j,1)-geopg(i,j,2)
          geogt(i,j,2)=geopg(i,j,2)-geopg(i,j,3)
        enddo
      enddo

      do l=1,ntl

c ***  calculate temperatures and add mean temperature level

        do j=1,nlon
          do i=1,nlat
            tempgt(i,j,l)=tempfac(l)*geogt(i,j,l)+tempm(l)
          enddo
        enddo

      enddo

      do j=1,nlon
        do i=1,nlat
          temp2g(i,j)=tempgt(i,j,1)
          temp4g(i,j)=tempgt(i,j,2)
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine vortfor
c-----------------------------------------------------------------------
c *** computes vorticity forcing
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)
      real*8  forcgg(nlat,nlon),forcgg2(nlat,nlon)
      real*8  zetas(nsh2,3),zetag(nlat,nlon,3)
      real*8  dimfac

      call rggtosp(vhforg1,dfor1)
      call rggtosp(vhforg2,dfor2)

c ***  compute the relative vorticity zeta from steamfunction psi
c ***  psi is dimensionless, zeta with dimension

      do l=1,nvl
        do k=1,nsh2
          zetas(k,l)=rinhel(k,0)*psi(k,l)*om
        enddo
        call sptogg(zetas(1,l),zetag(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=0d0
          enddo
        enddo
      enddo

c *** potential vorticity (pv) forcing (1/(second**2)) due to
c *** the diabatic heating

      if (ipvf1.ne.0) call pvf1(vforg)

c *** pv forcing due to advection of the planetary vorticy f
c *** by the divergent wind

      if (ipvf2.ne.0) call pvf2(vforg)

c *** pv forcing due to d*zeta. d is the divergence.
c *** zeta is the relative vorticity

      if (ipvf3.ne.0) call pvf3(vforg,zetag)

c *** pv forcing due to advection of the relative vorticy
c *** zeta by the divergent wind

      if (ipvf4.ne.0) call pvf4(vforg,zetas)

c *** pv forcing due to

      if (ipvf5.ne.0) call pvf5(vforg,zetag)

c *** pv forcing due to

      if (ipvf6.ne.0) call pvf6(vforg)

c *** pv forcing due to

      if (ipvf7.ne.0) call pvf7(vforg)

c *** the total pv forcing in nondimensional units

      dimfac=1./(om**2)

      call ggtosp(vforg(1,1,1),vfor1)
      call ggtosp(vforg(1,1,2),vfor2)
      call ggtosp(vforg(1,1,3),vfor3)

c *** adding the artificial forcing (forcgg1)

      vfor1(1)=0.d0
      vfor2(1)=0.d0
      vfor3(1)=0.d0

      do k=2,nsh2
        for(k,1)=dimfac*vfor1(k)
        for(k,2)=dimfac*vfor2(k)
        for(k,3)=dimfac*vfor3(k)
      enddo

c *** transfer to grid point

      call sptogg(vfor1,vforg1,pp)
      call sptogg(vfor2,vforg2,pp)
      call sptogg(vfor3,vforg3,pp)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf1(vforg)
c-----------------------------------------------------------------------
c *** potential vorticity (pv) forcing (1/(second**2)) due to
c *** the adiabatic heating
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  pvf1s(nsh2,3),vforg(nlat,nlon,nvl)
      real*8  vhfor1(nsh2),vhfor2(nsh2)
      real*8  vhforg1x(nlat,nlon),vhforg2x(nlat,nlon)
      real*8  vorfac1,vorfac2,drdef1,drdef2

      vorfac1=+rgas*dp/3.5d+4
      vorfac2=+rgas*dp/6.5d+4
      drdef1=1./( om*fzero*(radius*rrdef1)**2 )
      drdef2=1./( om*fzero*(radius*rrdef2)**2 )

      do j=1,nlon
        do i=1,nlat
          vhforg1x(i,j)=vhforg1(i,j)*sinfi(i)/fzero
          vhforg2x(i,j)=vhforg2(i,j)*sinfi(i)/fzero
        enddo
      enddo

      call rggtosp(vhforg1x,vhfor1)
      call rggtosp(vhforg2x,vhfor2)

      pvf1s(1,1)=0d0
      pvf1s(1,2)=0d0
      pvf1s(1,3)=0d0

      do k=2,nsh2
        pvf1s(k,1)=-drdef1*vorfac1*vhfor1(k)
        pvf1s(k,2)=drdef1*vorfac1*vhfor1(k)-drdef2*vorfac2*vhfor2(k)
        pvf1s(k,3)=drdef2*vorfac2*vhfor2(k)
      enddo

      do l=1,nvl
        call sptogg(pvf1s(1,l),vforg(1,1,l),pp)
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf2(vforg)
c-----------------------------------------------------------------------
c *** pv forcing due to advection of the planetary vorticy f
c *** by the divergent wind
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)

      do l=1,nvl
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)-vdivg(i,j,l)*om*cosfi(i)/radius
          enddo
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf3(vforg,zetag)
c-----------------------------------------------------------------------
c *** pv forcing due to d*zeta. d is the divergence.
c *** zeta is the relative vorticity
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),zetag(nlat,nlon,nvl)

      do l=1,nvl
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)-divg(i,j,l)*zetag(i,j,l)
          enddo
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf4(vforg,zetas)
c-----------------------------------------------------------------------
c *** pv forcing due to advection of the relative vorticy
c *** zeta by the divergent wind
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),zetas(nsh2,nvl)
      real*8  x(nsh2),xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)

      do l=1,nvl
        do k=1,nsh2
          x(k)=zetas(k,l)
        enddo
        call ddl(x,xhelp)
        call sptogg(xhelp,dxdl,pp)
        call sptogg(x,dxdm,pd)
        do i=1,nlat
          do j=1,nlon
             vforg(i,j,l)=vforg(i,j,l)
     &                  -udivg(i,j,l)*dxdl(i,j)/(radius*cosfi(i))
     &                  -vdivg(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf5(vforg,zetag)
c-----------------------------------------------------------------------
c *** computes vorticity forcing due to
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)
      real*8  omeg200(nlat,nlon),omeg500(nlat,nlon),omeg800(nlat,nlon)
      real*8  zetag(nlat,nlon,nvl),z350(nlat,nlon),z650(nlat,nlon)

      do j=1,nlon
        do i=1,nlat
          omeg200(i,j)=omegg(i,j,1)/2.d0
          omeg500(i,j)=(omegg(i,j,1)+omegg(i,j,2))/2.d0
          omeg800(i,j)=(omegg(i,j,2)+omegg(i,j,3))/2.d0
          z350(i,j)=(zetag(i,j,1)+zetag(i,j,2))/2.d0
          z650(i,j)=(zetag(i,j,2)+zetag(i,j,3))/2.d0

          vforg(i,j,1)=vforg(i,j,1)-
     &                 omeg200(i,j)*(z350(i,j)-0.d0)/35000.
          vforg(i,j,2)=vforg(i,j,2)-
     &                 omeg500(i,j)*(z650(i,j)-z350(i,j))/30000.d0
          vforg(i,j,3)=vforg(i,j,3)-
     &                 omeg800(i,j)*(0.d0-z650(i,j))/35000.d0
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf6(vforg)
c-----------------------------------------------------------------------
c *** computes vorticity forcing due to
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)
      real*8  o(nsh2,3),ohelp(nsh2),dodl(nlat,nlon),dodm(nlat,nlon)
      real*8  u350(nlat,nlon),u650(nlat,nlon),v350(nlat,nlon),v650(nlat,nlon)
      real*8  up(nlat,nlon,3),vp(nlat,nlon,3)

      do j=1,nlon
        do i=1,nlat
          u350(i,j)=(u200(i,j)+u500(i,j))/2.d0
          u650(i,j)=(u500(i,j)+u800(i,j))/2.d0
          v350(i,j)=(v200(i,j)+v500(i,j))/2.d0
          v650(i,j)=(v500(i,j)+v800(i,j))/2.d0
          up(i,j,1)=(u350(i,j)-0.d0)/35000.d0
          up(i,j,2)=(u650(i,j)-u350(i,j))/30000.d0
          up(i,j,3)=(0.d0-u650(i,j))/35000.d0
          vp(i,j,1)=(v350(i,j)-0.d0)/35000.d0
          vp(i,j,2)=(v650(i,j)-v350(i,j))/30000.d0
          vp(i,j,3)=(0.d0-v650(i,j))/35000.d0
        enddo
      enddo

      do k=1,nsh2
        o(k,1)=omegs(k,1)/2.d0
        o(k,2)=(omegs(k,1)+omegs(k,2))/2.d0
        o(k,3)=(omegs(k,2)+omegs(k,3))/2.d0
      enddo

      do l=1,nvl

        call ddl(o(1,l),ohelp)
        call sptogg(ohelp,dodl,pp)
        call sptogg(o(1,l),dodm,pd)

        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)+
     &                  up(i,j,l)*dodm(i,j)*cosfi(i)/radius
     &                 -vp(i,j,l)*dodl(i,j)/(radius*cosfi(i))
          enddo
        enddo
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf7(vforg)
c-----------------------------------------------------------------------
c *** computes vorticity forcing due to
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),pvf7g(nlat,nlon,nvl)
      real*8  ud(nlat,nlon,2),vd(nlat,nlon,2)
      real*8  x(nsh2,2),y(nlat,nlon,3)
      real*8  xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)
      real*8  rr1,rr2,sinfact

      rr1=1.d0/(rrdef1*radius)**2
      rr2=1.d0/(rrdef2*radius)**2

      do j=1,nlon
        do i=1,nlat
          ud(i,j,1)=(udivg(i,j,1)+udivg(i,j,2))/2.d0
          ud(i,j,2)=(udivg(i,j,2)+udivg(i,j,3))/2.d0
          vd(i,j,1)=(vdivg(i,j,1)+vdivg(i,j,2))/2.d0
          vd(i,j,2)=(vdivg(i,j,2)+vdivg(i,j,3))/2.d0
        enddo
      enddo

      do k=1,nsh2
        x(k,1)=(psi(k,1)-psi(k,2))*om*radius**2
        x(k,2)=(psi(k,2)-psi(k,3))*om*radius**2
      enddo

      do l=1,2
        call ddl(x(1,l),xhelp)
        call sptogg(xhelp,dxdl,pp)
        call sptogg(x(1,l),dxdm,pd)

        do j=1,nlon
          do i=1,nlat
             y(i,j,l)=ud(i,j,l)*dxdl(i,j)/(radius*cosfi(i))
     &               +vd(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
          pvf7g(i,j,1)=rr1*y(i,j,1)
          pvf7g(i,j,2)=-rr1*y(i,j,1)+rr2*y(i,j,2)
          pvf7g(i,j,3)=-rr2*y(i,j,2)
        enddo
      enddo

      do i=1,nlat
        sinfact=(sinfi(i)/fzero)**2
        do j=1,nlon
          vforg(i,j,1)=vforg(i,j,1)+pvf7g(i,j,1)*sinfact
          vforg(i,j,2)=vforg(i,j,2)+pvf7g(i,j,2)*sinfact
          vforg(i,j,3)=vforg(i,j,3)+pvf7g(i,j,3)*sinfact
        enddo
      enddo

      return
      end
