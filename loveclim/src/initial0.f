












c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_initecbilt
c-----------------------------------------------------------------------
c *** initialisation of ECBilt
c-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comphys.h'
      

c *** open files
 
      include 'openatinpfiles.h'

      call ec_iniparameterat
      

c******change H. Renssen, 25-02-2003
      if (iscencel.eq.1) call celest
      if (iscencel.eq.2) call bretagnon
c******end of change   

      call ec_inierror
      call ec_inimdldate
      call ec_iatmpar
      call ec_iatmdyn
      call ec_iatmphys
      call ec_inioutparat
      call ec_atmstate
      return
      end
      

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_iatmpar
      implicit none
      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      
      real*8  ari(nlat/2)
      real*8  rj,rw,dumwei

      integer ilat,k,ird,i,j

c *** constants
c *** rlatvap: latent heat of condensation in J/kg
c *** rlatsub: latent heat of sublimation in J/kg
c *** rlatfus: latent heat of fusion in J/kg
 
      pi=4d0*datan(1d0)
      fzero=sin(0.25*pi)
c     let op om is 2 keer de hoeksnelheid van de aarde !!!!
      om=4d0*pi/(24d0*3600d0)
      grav=9.8
      rgas=287.
      radius=6.37e+6 
      rowat=1000.
      roair=1.25
      rlatvap=2.5e+06
      rlatsub=2.8e+06
      rlatfus=0.3e+06
      sboltz=5.67e-08
      cwater=4180.
      cpair=1004.
      cvair=717.
      gamma=cpair/cvair
      rkappa=(gamma-1.)/gamma

      plevel(1)=2.0d+4
      plevel(2)=5.0d+4
      plevel(3)=8.0d+4

      tlevel(1)=3.5d+4
      tlevel(2)=6.5d+4

      dp=plevel(2)-plevel(1)

      dp0=2.0d+4
      dp1=3.0d+4
      dp2=5.0d+4

      p0=1d5

      rlogtl12=1d0/log(tlevel(1)/tlevel(2))
      alogtl12=log(tlevel(1)/tlevel(2))
      alogtl1pl2=log(tlevel(1)/plevel(2))
      alogpl2tl2=log(plevel(2)/tlevel(2))

      potfac1=(tlevel(1)/p0)**rkappa
      potfac2=(tlevel(2)/p0)**rkappa

      gamd=grav/cpair
      tzero=273.15d0
      alphad=roair*cpair
      alphas=roair*rlatsub
      alphav=roair*rlatvap
      
c *** constants in clausius clapeyron relation

      cc1=0.662*611.2
      cc2=17.67
      cc3=29.66

c *** definitions for tabel of qmax values used in iatmphys
c *** i corresponds to temperature at 650 hPa
c *** j corresponds to temperature difference between ground and 650 hPa
c *** k corresponds to temperature difference between 650 and 350 hPa

      tqmimin=200d0
      dtqmi=2d0
      rdtqmi=1d0/dtqmi
      tqmjmin=-10d0
      dtqmj=2d0
      rdtqmj=1d0/dtqmj
      tqmkmin=5d0
      dtqmk=2d0
      rdtqmk=1d0/dtqmk
      
c *** time step of the atmospheric model: 
c *** dt    : fraction of one day
c *** dtime : in seconds
c *** dtt   : dimensionless
 
      dt     = 1d0/dble(iatm)
      dtime  = dt*(24d0*3600d0)
      rdtime = 1d0/dtime
      dtt    = dt*pi*4d0
       
c *** gauss points and weights
      ilat=nlat/2
      rewind(iuo+7)
  10  continue
        read(iuo+7,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
          do i=1,ird
            read(iuo+7,220) ari(i),dumwei
          enddo
          goto 20
        else
          goto 10
        endif 
  15    continue
        call ec_error(4)
  20  continue
      
      do i=1,ilat
        phi(i)=-ari(ilat+1-i)
        phi(ilat+i)=ari(i)
      enddo

      do i=1,nlat
        phi(i)=asin(phi(i))
      enddo

      do i=1,nlat
        cosfi(i)=cos(phi(i))
        sinfi(i)=sin(phi(i))
        tanfi(i)=tan(phi(i))
      enddo

c *** land/sea/sea-ice fraction

      
      do j=1,nlon
        do i=1,nlat
	  fractoc(i,j)=fracto(i,j)
          fractn(i,j,noc)=fractoc(i,j)
          fractn(i,j,nld)=1-fractoc(i,j)
          fractn(i,j,nse)=0.0
	enddo
      enddo 
      
c      do j=nlat,1,-1
c         write (6, '(i3,x,64(i1))') j, (int(fracto(j,i)), i = 1, nlon)
c      enddo

      epss=1d-10


220   format(f18.10,f17.10)

      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_iniparameterat
c-----------------------------------------------------------------------
c *** initialisation of parameters from the namelist 
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comatfor.h'
      include 'comunit.h'
      include 'comrunlabel.h'

      NAMELIST /runatctl/ iadyn,iaphys,ipert,initfield,initdate
      NAMELIST /dispar/   tdis,addisl,addish,trel,tdif,idif,initfield,initdate
      NAMELIST /dfmpar/  rrdef1,rrdef2,h0
      NAMELIST /moipar/  ihavm,ivavm,imsink,tdifq,gpm500,relhmax,
     *                   hmoisr,umoisr,rainmax
      NAMELIST /cloudpar/ relhcrit,relhfac
      NAMELIST /forpar/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      NAMELIST /radpar/    solarc,iradcloud,iscenghg,iscenghg2s,
     *                     iscentsi,iscenvol,
     *                     iscensul,isceno3,
     *                     iscencel,iens,numens,emisoc,emisse,emisld,
     *                     albin,albis,albice,alphd,alphdi,alphs,cgren,
     *                     albocef,facttsi,bup,AMPWIR,EXPIR,HPROFTROP,
     *                     HPROFEQ,HPROFW,AMPEQIR,HPROFAN,AMPANIR,
     *                     eccf,oblf,omwebf,AMPANIR2,HPROFAN2
      NAMELIST /satfor/   isatfor,nbsatfor,nafyear,iclimflux
      NAMELIST /fluxpar/ cdrag,cwdrag,dragan,dragla,uv10rfx,uv10m,
     *                   uv10rws,ndayws
      NAMELIST /fluxcorw/ corAN,corPN,corAC,corID,corAS,corPS,
     *                   corAA


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C intgrtn parameter:                                                   C
C iadyn:      with (1) or without (0) atmospheric dynamics.            C
C iaphys:     with (1) or without (0) atmospheric physics.             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C disspt parameter:                                                    C
C tdis:       Ekman dissipation time scale atlower level               C
C addisl:     parameter for computing dissipation time scale for land. C
C addish:     parameter for computing dissipation time scale over      C
C             mountains                                                C
C trel:       temperature relaxation time scale.                       C
C tdif:       damping time scale of hyperviscosity of the smallest     C
C             waves at all levels                                      C
C idif:       determines scale-selectivity of hyperviscosity: power of C
C             Laplace operator on PV                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C dfmtion parameter:                                                   C
C rrdef1:     rossby radius of deformation layer 1 (200 - 500 hPa).    C
C rrdef2:     rossby radius of deformation layer 2 (500 - 800 hPa).    C
C h0:         mountain scale height.                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C moisture Parameter:                                                  C
C ihavm:      with (1) or without (0) horizontal divergence of moistureC
C ivavm:      with (1) or without (0) vertical divergence of moisture. C
C imsink:     with (1) or without (0) the source and sink of moisture. C
C tdifq:      diffusion time scale for moisture in days                C
C gpm500  :   mean 500 hPa geopotential height [m] used in the         C
C             calculation of the groundpressure and temperature        C
C relhmax:    maximum relative humidity before saturation.             C
C hmoisr:     reduction factor of mountain heights in order to tune    C
C             the amount of water that is allowed to pass a            C
C             topographic barier                                       C
C umoisr:     reduction factor of 800 hPa winds used in advecting the  C
C             moisture                                                 C
C rainmax:    maximum rate of dynamic and convective rain in m/s.      C
C             if the maximum is reached by both type of rain           C
C             oversaturation occurs.                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Surface flux Parameters:                                             C
C cdrag:      coefficient in sensible and latent air-surface heat flux C
C cwdrag:     coefficient in winds stress                              C
C uv10rfx:    reduction factor of 800 hPa winds to obtain              C
C             10 mtr wind used in the definition of the surface        C
C             heat fluxes                                              C
C             used in subroutine wind10 of atmdyn.f                    C
C uv10m:      minimum value of 10 mtr wind used in all surface fluxes  C
C             used in subroutine wind10 of atmdyn.f                    C
C uv10rws:    reduction factor of 800 hPa winds to obtain              C
C             10 mtr wind used in the definition of the wind stress    C
C             to drive the ocean currents                              C
C             used in subroutine wind10 of atmdyn.f                    C
C dsnm:       maximum land snow depth                                  C
C bmoism:     maximum bottom moisture.                                 C
C ndayws:     averaging period of wind stresses in days                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Flux correction parameters                                           C
C corAN       flux correction in the North Atlantic                    C
C corPN       flux correction in the North Pacific                     C
C corAC       flux correction in the Arctic                            C
C corID       flux correction in the Indian                            C
C corAS       flux correction in the South Atlantic                    C
C corPS       flux correction in the South Pacific                     C
C corAA       flux correction in the Southern Ocean                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C forcing parameter:                                                   C
C iartif:     with (1) or without (0) artificial forcing               C 
C ipvf1 :     with (1) or without (0) diabatic heating                 C
C ipvf2 :     with (1) or without (0) advection of f by divergent wind C
C ipvf3 :     with (1) or without (0) stretching term                  C
C ipvf4 :     with (1) or without (0) advection of zeta by divergent   C
C             wind                                                     C
C ipvf5 :     with (1) or without (0) advection of temperature by      C
C                                  divergent wind                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C tcont parameter:                                                     C
C solarc:     solar constant.                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C satfor parameter:                                                    C
C isatfor:    if (1) saves nafyear of atmospheric data to disk to be   C 
C             used to drive the ocean in uncoupled mode                C
C nbsatfor:   first year in the integrations to start saving nafyears  C
C             of data                                                  C
C nafyear:    number of years of data to save                          C
C iclimflux:  if (1) daily climatological sst and heatflux are output  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      iadyn    = 1
      iaphys   = 1
      ipert    = 0
      initfield= 1
      initdate = 0

      tdis   = 3.0
      addisl = 0.5
      addish = 0.5
      trel   = 50.0
      tdif   = 1.0
      idif   = 4

      rrdef1 = .110
      rrdef2 = .070
      h0     = 3.

      ihavm  = 1
      ivavm  = 1
      imsink = 1
      tdifq  = 5d0
      gpm500 = 5500
      relhmax= 1.0
      hmoisr = 1d0
      umoisr = 0.8
      rainmax= 1e-5
      
      cdrag  = 1.4e-03
      cwdrag = 2.0e-03
      dragan = 10.0
      dragla = 10.0
      uv10rfx= 0.8
      uv10m  = 4.
      uv10rws= 0.8
      ndayws = 30

      corAN  = -0.075d0
      corPN  =  1.0d0
      corAC  = -0.20d0
      corID  =  0.0d0
      corAS  = -0.075d0
      corPS  =  0.d0
      corAA  =  0.d0

      relhcrit = 0.5d0
      relhfac  = 1.0d0


      iartif = 0
      ipvf1  = 1
      ipvf2  = 1
      ipvf3  = 1
      ipvf4  = 1
      ipvf5  = 1

C     solarc = 1353
      solarc = 1365
      solarm=solarc
      iradcloud= 0
      iscenghg = 0
      iscentsi = 0
      facttsi =1.0
      iscenvol = 0
      iscensul = 0
      isceno3 = 0
      iscencel = 0
      eccf   = 0.016724
      oblf   = 23.446
      omwebf = 102.04 
      iens=0
      bup=0.13
      numens=1
      emisoc=1.0
      emisse=1.0
      emisld=1.0
      albin = 0.43
      albis = 0.43
      albice = 0.43
      alphd  = 0.70
      alphdi = 0.62
      alphs  = 0.53
      cgren = 0.04
      albocef =1.0
      AMPWIR    = 1.0
      AMPEQIR   = 1.0
      EXPIR     = 0.3333
      HPROFTROP = 1.0
      HPROFEQ   = 1.0
      HPROFW    = 1.0
      HPROFAN   = 1.0
      AMPANIR   = 1.0
      HPROFAN2  = 1.0
      AMPANIR2  = 1.0


      isatfor = 0
      nbsatfor= 0
      nafyear = 0
      iclimflux = 0


      read(iuo+15, NML = runatctl)
      read(iuo+15, NML = dispar)
      read(iuo+15, NML = dfmpar)
      read(iuo+15, NML = moipar)
      read(iuo+15, NML = cloudpar)
      read(iuo+15, NML = forpar)
      read(iuo+15, NML = radpar)
      read(iuo+15, NML = satfor)
      read(iuo+15, NML = fluxpar)
      read(iuo+15, NML = fluxcorw)

      write(iuo+30, 900) 'iaphys   =', iaphys
      write(iuo+30, 900) 'iadyn    =', iadyn
      write(iuo+30, 900) 'ipert    =', ipert
      write(iuo+30, 900) 'initfield=', initfield
      write(iuo+30, 900) 'initdate =', initdate

      irunlabelf=irunlabel-initdate 
      write(iuo+30, 900) 'irunlabelf=', irunlabelf

      write(iuo+30, 910) 'tdis     =', tdis
      write(iuo+30, 910) 'addisl   =', addisl
      write(iuo+30, 910) 'addish   =', addish
      write(iuo+30, 910) 'trel     =', trel
      write(iuo+30, 910) 'tdif     =', tdif
      write(iuo+30, 900) 'idif     =', idif

      write(iuo+30, 910) 'h0       =', h0
      write(iuo+30, 910) 'rrdef1   =', rrdef1
      write(iuo+30, 910) 'rrdef2   =', rrdef2

      write(iuo+30, 910) 'relhcrit =', relhcrit
      write(iuo+30, 910) 'relhfac  =', relhfac

      write(iuo+30, 910) 'cdrag    =', cdrag
      write(iuo+30, 910) 'cwdrag   =', cwdrag
      write(iuo+30, 910) 'dragan   =', dragan
      write(iuo+30, 910) 'dragla   =', dragla
      write(iuo+30, 910) 'uv10rfx  =', uv10rfx
      write(iuo+30, 910) 'uv10m    =', uv10m
      write(iuo+30, 910) 'uv10rws  =', uv10rws
      write(iuo+30, 900) 'ndayws   =', ndayws
      write(iuo+30, 910) 'corAN    =', corAN
      write(iuo+30, 910) 'corPN    =', corPN
      write(iuo+30, 910) 'corAC    =', corAC
      write(iuo+30, 910) 'corID    =', corID
      write(iuo+30, 910) 'corAS    =', corAS
      write(iuo+30, 910) 'corPS    =', corPS
      write(iuo+30, 910) 'corAA    =', corAA


      write(iuo+30, 900) 'ihavm    =', ihavm
      write(iuo+30, 900) 'ivavm    =', ivavm
      write(iuo+30, 900) 'imsink   =', imsink
      write(iuo+30, 910) 'tdifq    =', tdifq
      write(iuo+30, 910) 'gpm500   =', gpm500
      write(iuo+30, 910) 'relhmax  =', relhmax
      write(iuo+30, 910) 'hmoisr   =', hmoisr
      write(iuo+30, 910) 'umoisr   =', umoisr
      write(iuo+30, 910) 'rainmax  =', rainmax


      write(iuo+30, 900) 'iartif   =', iartif
      write(iuo+30, 900) 'ipvf1    =', ipvf1
      write(iuo+30, 900) 'ipvf2    =', ipvf2
      write(iuo+30, 900) 'ipvf3    =', ipvf3
      write(iuo+30, 900) 'ipvf4    =', ipvf4
      write(iuo+30, 900) 'ipvf5    =', ipvf5

      write(iuo+30, 910) 'solarc   =', solarc
      write(iuo+30, 900) 'iradcloud=', iradcloud
      write(iuo+30, 900) 'iscenghg =', iscenghg
      write(iuo+30, 900) 'iscentsi =', iscentsi
      write(iuo+30, 910) 'facttsi  =', facttsi
      write(iuo+30, 900) 'iscenvol =', iscenvol
      write(iuo+30, 900) 'iscensul =', iscensul
      write(iuo+30, 910) 'bup      =', bup

      write(iuo+30, 900) 'isatfor  =', isatfor
      write(iuo+30, 900) 'nbsatfor =', nbsatfor
      write(iuo+30, 900) 'nafyear  =', nafyear
      write(iuo+30, 900) 'iclimflux=', iclimflux

      call flush(iuo+30)
      emisn(noc)=emisoc
      emisn(nse)=emisse
      emisn(nld)=emisld


900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end



c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_inioutparat
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comdiag.h'
      include 'comunit.h'

      character*60 part1
      integer      totvar(80,20)
      integer      i,l


      NAMELIST /outatctl/ixout,ifrendat
      NAMELIST /flgout/irad
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ixout:        output frequency for instantanous field in days        C
c meantype      output monthly mean fields.			       C
c meantype      output seasonal mean fields.			       C
c meantot       computes whole period monthly or seasonal mean fields. C
c meanyl        computes yearly monthly or seansonal mean fields.      C
c ioutdaily     output instantanous fields.			       C
c ioutyearly    output yearly mean fields.			       C
C ts:           output surface temperature.                            C
C t:            output temperature.                                    C
C u             output wind component U.                               C
C v:            output wind component V.                               C
C om:           output wind component omega.                           C
C psi:          output stream function.                                C
C shf:          output surface sensible heat flux.                     C
C lhf:          output surface latent heat flux.                       C
C lsp:          output large scale precipitation.                      C
C cp:           output convective precipitation.                       C
C q:            output specific humidity.                              C
C r:            output relative humidity.                              C
C ageu:         output ageostrophic wind component U.                  C
C agev:         output ageostrophic wind component V.                  C
C ssr:          output surface solar radiation (downward).             C        
C tsr:          output top solar radiation (downward).                 C
C ttr:          output top thermal radiation (upward).		       C
C str:          output surface thermal radiation (upward).             C
C pp:		output total percipitation.			       C 
C evap:		output evaporation.				       C
C eminp:	output evaporation minus precipitation.                C 
c albs:         output surface albedo.                                 C
c albp:         output planetary albedo.                               C
c ustress:      output u wind stress.                                  C
c vstress:      output v wind stress.                                  C
c runoffo:      output runoff over ocean.                              C
c runoffl:      output runoff over land.                               C
c sdl:          output snow depth over land.                           C
c chi:          output velocity potential.                             C
c sp:           output surface pressure                                C
c uv10:         output wind magnitude at 10 meter height.              C
c cdragw:       output drag coefficient of wind.                       C
c cdragv:       output drag coefficient of sensible heat flux.         C
c richar:       output richardson number.                              C
c qgpv:         output quasi geostrophic potential vorticity.          C
c gh:           output geopotential height.                            C
c tcc:          output total cloud cover.                              c
c dumt1:        free output variable at 350 hPa, 650 hPa and 1000 hPa  c
c dumt2:        free output variable at 350 hPa, 650 hPa and 1000 hPa  c
c dumu1:        free output variable at 200 hPa, 500 hPa and 800 hPa   c
c dumu2:        free output variable at 200 hPa, 500 hPa and 800 hPa   c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C** Here following are default values for above mentioned parameters. 
C** These parameters can be updated in the namelist.

      itel=0
      instcount=0

      ixout  = 30
      ioutdaily = 0
      ioutyearly = 0
      meantype = 0
      meantot  = 0
      meanyl   = 0
      irad     = 1
      thirdd(:)= 0

      open(iuo+65,file='outp_atmos.param')
      read(iuo+65,'(/,/,/,/,/,/,/,/)')
      read(iuo+65,*) numtotvar,fill_value,missing_value
      read(iuo+65,*)

      do i=1,numtotvar

         read(iuo+65,"(A)") part1
         nametotvar(i,1)=trim(part1)
         read(iuo+65,*) (nametotvar(i,l),l=2,5)
         read(iuo+65,*) (newtotvar(i,l),l=1,7)
         do l=1,6
            newtotvar(i,l)=newtotvar(i,l+1)
         enddo

         IF ( newtotvar(i,4)==1 ) ioutyearly = 1
         IF ((newtotvar(i,3)==1).OR.(newtotvar(i,2)==1) ) meantype = 1
         IF ( newtotvar(i,2)==1 ) meanyl     = 1
         IF ( newtotvar(i,3)==1 ) meantot    = 1
         IF ( newtotvar(i,1)==1 ) ioutdaily  = 1

         IF ((newtotvar(i,2)==1).OR.(newtotvar(i,3)==1)
     &                          .OR.(newtotvar(i,4)==1) ) THEN
	         SELECT CASE ( nametotvar(i,5) )
	         CASE ("T2")
	            thirdd(1)=1
	         CASE ("T3")
	            thirdd(2)=1
	         CASE ("T4")
	            thirdd(3)=1
	         CASE ("U3")
	            thirdd(4)=1
	         CASE ("N")
	            l=0
	         CASE DEFAULT
	            call ec_error(123)
	         END SELECT
	      END IF

      enddo

      close(iuo+65)

      read(iuo+15, NML = outatctl)
c      read(iuo+15, NML = wratpar)
      read(iuo+15, NML = flgout)

      return
      end
