c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine initecbilt
c-----------------------------------------------------------------------
c *** initialisation of ECBILT
c-----------------------------------------------------------------------
      implicit none

      include 'comglobal.h'

c *** caution !
c *** do not change the order of calling initialisation routines

      call iniparameterat
      call iniparameteroc
      call iniglobal
      call inierror
      call inimdldate
      call initatmmodel
      call inioceanfixed
      call inirunoff
      call iniland
      call inioutparat
C xueli
c      call inioutparoc

      if (irunatm.eq.1) then
        if (irunlabel.gt.0) call atmstate
        call radiation
        call senhflux
        call lathflux
      endif

      if (isatfor.eq.1) call outato(0)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iniparameterat
c-----------------------------------------------------------------------
c *** initialisation of parameters from the namelist
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'
      include 'comcoup.h'
      include 'comatfor.h'
      include 'comocfix.h'

      NAMELIST /runatctl/  nyears,irunlabel,iatm,nwrskip,iadyn,iaphys
      NAMELIST /dispar/   tdis,addisl,addish,trel,tdif,idif
      NAMELIST /dfmpar/  rrdef1,rrdef2,h0
      NAMELIST /moipar/  ihavm,ivavm,imsink,tdifq,gpm500,relhmax,
     *                   hmoisr,umoisr,rainmax
      NAMELIST /fluxpar/ cdrag,uv10rfx,uv10m,uv10rws,dsnm,bmoism,ndayws
      NAMELIST /forpar/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5,ipvf6,
     *                    ipvf7
      NAMELIST /radpar/    solarc
      NAMELIST /satfor/   isatfor,nbsatfor,nafyear,irunatm


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C intgrtn parameter:                                                   C
C nyears:     integration period in years                              C
C irunlabel:  number of startfiles inat###.dat and inoc###.dat         C
C             inat001.dat is default                                   C
C iatm:       number of atmospheric timesteps in one day               C
C nwrskip:    number of years between writing the model state to disk  C
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
C forcing parameter:                                                   C
C iartif:     with (1) or without (0) artificial forcing               C
C ipvf1 :     with (1) or without (0) diabatic heating                 C
C ipvf2 :     with (1) or without (0) advection of f by divergent wind C
C ipvf3 :     with (1) or without (0) stretching term                  C
C ipvf4 :     with (1) or without (0) advection of zeta by divergent   C
C             wind                                                     C
C ipvf5 :     with (1) or without (0) vertical advection of zeta       C
C ipvf6 :     with (1) or without (0) solenoidal term                  C
C ipvf7 :     with (1) or without (0) advection of temperature by      C
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
C irunatm:    if (1) the atmospheric model is run, if (0) the ocean is C
C             forced with atmospheric data read from disk              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      nyears   = 1
      irunlabel= 1
      iatm     = 6
      nwrskip  = 10
      iadyn    = 1
      iaphys   = 1


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
      uv10rfx= 0.8
      uv10m  = 4.
      uv10rws= 0.8
      bmoism = 0.15
      dsnm   = 2000.
      ndayws = 30

      iartif = 0
      ipvf1  = 1
      ipvf2  = 1
      ipvf3  = 1
      ipvf4  = 1
      ipvf5  = 1
      ipvf6  = 1
      ipvf7  = 0

      solarc = 1353

      isatfor = 0
      nbsatfor= 0
      nafyear = 0
      irunatm = 1


      read(15,NML = runatctl)
      read(15,NML = dispar)
      read(15,NML = dfmpar)
      read(15,NML = moipar)
      read(15,NML = fluxpar)
      read(15,NML = forpar)
      read(15,NML = radpar)
      read(15,NML = satfor)

      write(fini,1) irunlabel
      write(fend,1) irunlabel+nyears
  1   format(i4.4)

      include 'openatstartfiles.h'
      include 'openatoutfiles.h'

      write(30, *) 'parameters of atmosphere only run'
      write(30, *)
      write(30, 900) 'nyears   =', nyears
      write(30, 900) 'irunlabel=', irunlabel
      write(30, 900) 'iatm     =', iatm
      write(30, 900) 'nwrskip  =', nwrskip
      write(30, 900) 'iaphys   =', iaphys
      write(30, 900) 'iadyn    =', iadyn

      write(30, 910) 'tdis     =', tdis
      write(30, 910) 'addisl   =', addisl
      write(30, 910) 'addish   =', addish
      write(30, 910) 'trel     =', trel
      write(30, 910) 'tdif     =', tdif
      write(30, 900) 'idif     =', idif

      write(30, 910) 'h0       =', h0
      write(30, 910) 'rrdef1   =', rrdef1
      write(30, 910) 'rrdef2   =', rrdef2

      write(30, 900) 'ihavm    =', ihavm
      write(30, 900) 'ivavm    =', ivavm
      write(30, 900) 'imsink   =', imsink
      write(30, 910) 'tdifq    =', tdifq
      write(30, 910) 'gpm500   =', gpm500
      write(30, 910) 'relhmax  =', relhmax
      write(30, 910) 'hmoisr   =', hmoisr
      write(30, 910) 'umoisr   =', umoisr
      write(30, 910) 'rainmax  =', rainmax

      write(30, 910) 'cdrag    =', cdrag
      write(30, 910) 'uv10rfx  =', uv10rfx
      write(30, 910) 'uv10m    =', uv10m
      write(30, 910) 'uv10rws  =', uv10rws
      write(30, 910) 'bmoism   =', bmoism
      write(30, 910) 'dsnm     =', dsnm
      write(30, 900) 'ndayws   =', ndayws

      write(30, 900) 'iartif   =', iartif
      write(30, 900) 'ipvf1    =', ipvf1
      write(30, 900) 'ipvf2    =', ipvf2
      write(30, 900) 'ipvf3    =', ipvf3
      write(30, 900) 'ipvf4    =', ipvf4
      write(30, 900) 'ipvf5    =', ipvf5
      write(30, 900) 'ipvf6    =', ipvf6
      write(30, 900) 'ipvf7    =', ipvf7

      write(30, 910) 'solarc   =', solarc

      write(30, 900) 'isatfor  =', isatfor
      write(30, 900) 'nbsatfor =', nbsatfor
      write(30, 900) 'nafyear  =', nafyear
      write(30, 900) 'irunatm  =', irunatm

      call flush(30)

900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iniparameteroc
c-----------------------------------------------------------------------
c *** initialisation of parameters from the namelist of the ocean part
c *** of the model.
c-----------------------------------------------------------------------
      implicit none

      include 'comglobal.h'
      include 'comcouphelp.h'
      include 'comocean.h'
      include 'comcoup.h'
      include 'comice.h'

      NAMELIST /runoctl/  idtbclin,idtbtrop
      NAMELIST /couppar/  ndreloct,ndrelocs,icoupleh,icouplew,icouples,
     *                    ltflux,mico
      NAMELIST /icepar/   rks,rki,abottom
      NAMELIST /ocpar/    rkap,rkapv,rkaph,accstr
      NAMELIST /levpar/   h


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C runoctl parameter:                                                   C
C idtbclin:   time step of baroclinic part of ocean in days    	       C
C idtbtrop:   time step of barotropic part of ocean in days            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C couppar parameter:                                                   C
C icoupleh:   with or without heat flux exchange between atm and ocean.C
C icouplew:   with or without salinity exchange between atm and ocean. C
C icouples:   with or without windstress forcing from atm and ocean.   C
C ltflux:     true:  heat flux prescribed.                             C
C             false: relaxation.                                       C
C mico:       true:  salt flux prescribed.                             C
C             false: relaxation.                                       C
C ndreloct:   relaxaton constant for temperature in days               C
C ndrelocs:   relaxaton constant for salinity in days                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ocpar parameter:                                                     C
C rks:       thermal conductivity of snow.                             C
C rki:       thermal conductivity of ice.                              C
C abottom:   thermal exchange coefficient between ocean and atmosphere.C
C rkap:      stommel friction.                                         C
C rkapv:     vertical diffusion coefficient.                           C
C rkaph:     horizontal diffusion coefficient.                         C
C accstr:    strength of acc current.                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C levpar parameter:                                                    C
C h:          array containing depth's of ocean levels                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      idtbclin = 1
      idtbtrop = 1

      ndreloct = 20
      ndrelocs = 40
      icoupleh = 1
      icouplew = 1
      icouples = 1
      ltflux   = .false.
      mico     = .false.

      rks      = 0.3097
      rki      = 2.034
      abottom  = 20d0

      rkap     = 8d-6
      rkapv    = 3d-5
      rkaph    = 1d3
      accstr   = 100d6

      h(1)=30d0
      h(2)=50d0
      h(3)=80d0
      h(4)=140d0
      h(5)=250d0
      h(6)=350d0
      h(7)=400d0
      h(8)=450d0
      h(9)=500d0
      h(10)=550d0
      h(11)=600d0
      h(12)=600d0

      read(15,NML = runoctl)
      read(15,NML = couppar)
      read(15,NML = icepar)
      read(15,NML = ocpar)
      read(15,NML = levpar)

      write(30, 900) 'idtbclin =', idtbclin
      write(30, 900) 'idtbtrop =', idtbtrop

      write(30, 910) 'ndreloct =', ndreloct
      write(30, 910) 'ndrelocs =', ndrelocs
      write(30, 900) 'icoupleh =', icoupleh
      write(30, 900) 'icouplew =', icouplew
      write(30, 900) 'icouples =', icouples
      write(30, 930) 'ltflux   =', ltflux
      write(30, 930) 'mico     =', mico

      write(30, 910) 'rks      =',      rks
      write(30, 910) 'rki      =',      rki
      write(30, 910) 'abottom  =', abottom

      write(30, 910) 'rkap     =',   rkap
      write(30, 910) 'rkapv    =',   rkapv
      write(30, 910) 'rkaph    =',   rkaph
      write(30, 910) 'accstr   =',  accstr

      write(30, 910) 'h(1)     =',   h(1)
      write(30, 910) 'h(2)     =',   h(2)
      write(30, 910) 'h(3)     =',   h(3)
      write(30, 910) 'h(4)     =',   h(4)
      write(30, 910) 'h(5)     =',   h(5)
      write(30, 910) 'h(6)     =',   h(6)
      write(30, 910) 'h(7)     =',   h(7)
      write(30, 910) 'h(8)     =',   h(8)
      write(30, 910) 'h(9)     =',   h(9)
      write(30, 910) 'h(10)    =',   h(10)
      write(30, 910) 'h(11)    =',   h(11)
      write(30, 910) 'h(12)    =',   h(12)
      close(30)

900   format(a12,1x,i6)
910   format(a12,1x,e12.5)
920   format(a12,1x,a9)
930   format(a12,1x,L9)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inioutparat
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'
      include 'comdiag.h'


      integer ts(20),t(20),u(20),v(20),omega(20),psi(20),hforc(20),
     &        vforc(20),shf(20),lhf(20),lsp(20),cp(20),ageu(20),agev(20),
     &        ssr(20),tsr(20),str(20),ttr(20),bm(20),r(20),q(20),pp(20),
     &        evap(20),eminp(20),albs(20),ustress(20),vstress(20),sdl(20),
     &        runoffo(20),runoffl(20),uv10(20),chi(20),sp(20),albp(20),
     &        cdragw(20),cdragv(20),richar(20),qgpv(20),gh(20),tcc(20),
     &        dumt1(20),dumt2(20),dumu1(20),dumu2(20)

      integer i

      NAMELIST /outatctl/ixout,ioutdaily,meantype,
     &                    meanyl,meantot,ifrendat
      NAMELIST /wratpar/ ts,t,u,v,omega,psi,hforc,vforc,shf,lhf,pp,
     &                    lsp,cp,q,r,ageu,agev,ssr,tsr,ttr,str,bm,
     &                    evap,eminp,albs,ustress,vstress,sdl,
     &                    runoffo,runoffl,uv10,chi,sp,albp,
     &                    cdragw,cdragv,richar,qgpv,gh,tcc,dumu1,dumu2,
     &                    dumt1,dumt2


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ixout:        output frequency for instantanous field in days        C
c meantype      output monthly mean fields.			       C
c meantype      output seasonal mean fields.			       C
c meantot       computes whole period monthly or seasonal mean fields. C
c meanyl        computes yearly monthly or seansonal mean fields.      C
c ioutdaily     output instantanous fields.			       C
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

      ixout  = 30
      ioutdaily = 0
      meantype = 1
      meantot  = 0
      meanyl   = 1

      do i = 1, 20
        ts(i)      = 0
        t(i)       = 0
        u(i)       = 0
        v(i)       = 0
        omega(i)   = 0
        psi(i)     = 0
        shf(i)     = 0
        lhf(i)     = 0
        lsp(i)     = 0
        hforc(i)   = 0
        vforc(i)   = 0
        cp(i)      = 0
        q(i)       = 0
        r(i)       = 0
        ageu(i)    = 0
        agev(i)    = 0
        ssr(i)     = 0
        tsr(i)     = 0
        ttr(i)     = 0
        str(i)     = 0
        bm(i)      = 0
        pp(i)      = 0
        evap(i)    = 0
        eminp(i)   = 0
        albs(i)    = 0
        ustress(i) = 0
        vstress(i) = 0
        sdl(i)     = 0
        uv10(i)    = 0
        runoffo(i) = 0
        runoffl(i) = 0
        chi(i)     = 0
        sp(i)      = 0
        albp(i)    = 0
        cdragw(i)  = 0
        cdragv(i)  = 0
        richar(i)  = 0
        qgpv(i)    = 0
        gh (i)     = 0
        tcc (i)    = 0
        dumt1(i)   = 0
        dumt2(i)   = 0
        dumu1(i)   = 0
        dumu2(i)   = 0

      enddo

      read(15,NML = outatctl)
      read(15,NML = wratpar)

      do i = 1, 20
        newts(i)     = ts(i)
        newdyrain(i) = lsp(i)
        newcorain(i) = cp(i)
        neweflux(i)  = lhf(i)
        newhflux(i)  = shf(i)
        newrelhum(i) = r(i)
        newrmoisg(i) = q(i)
        newbmoisg(i) = bm(i)
        newtorain(i) = pp(i)
        newevap(i)   = evap(i)
        neweminp(i)  = eminp(i)
        newt(i)      = t(i)
        newu(i)      = u(i)
        newv(i)      = v(i)
        newomega(i)  = omega(i)
        newdivu(i)   = AGEU(i)
        newdivv(i)   = AGEV(i)
        newpsi(i)    = psi(i)
        newvhforg(i) = hforc(i)
        newvforg(i)  = vforc(i)
        newssr(i)    = ssr(i)
        newtsr(i)    = tsr(i)
        newttr(i)    = ttr(i)
        newstr(i)    = str(i)
        newalbs(i)   = albs(i)
        newustress(i)= ustress(i)
        newvstress(i)= vstress(i)
        newsdl(i)    = sdl(i)
        newuv10(i)   = uv10(i)
        newrunoffo(i)= runoffo(i)
        newrunoffl(i)= runoffl(i)
        newsp(i)     = sp(i)
        newchi(i)    = chi(i)
        newalbp(i)   = albp(i)
        newcdragw(i) = cdragw(i)
        newcdragv(i) = cdragv(i)
        newrichar(i) = richar(i)
        newqgpv(i)   = qgpv(i)
        newgh(i)     = gh(i)
        newtcc(i)    = tcc(i)
        newdumt1(i)  = dumt1(i)
        newdumt2(i)  = dumt2(i)
        newdumu1(i)  = dumu1(i)
        newdumu2(i)  = dumu2(i)
      enddo

      return
      end
C123456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inioutparoc
c-----------------------------------------------------------------------
      implicit none

      include 'comdiago.h'

      integer to(20),sa(20),uo(20),vo(20),wo(20),ro(20),hadvt(20),
     &        hadvs(20),vadvt(20),vadvs(20),hdift(20),hdifs(20),
     &        vdift(20),vdifs(20),convs(20),convt(20),relaxt(20),relaxs(20)

      integer heex(20),saex(20),brin(20),dice(20),
     &        sds(20), tice(20),heico(20)

      integer i

      NAMELIST /outoctl/isnapshot,ifreq
      NAMELIST /wrocpar/ to,sa,uo,vo,wo,ro,hadvt,hadvs,vadvt,
     &         vadvs,hdift,hdifs,vdift,vdifs,convt,convs,relaxt,relaxs,
     &         heex,saex,brin,dice,sds,tice,heico

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C isnapshot:    output (1) or not (0) instantanous field               C
c ifreq:	output frenqency for snapshot fields.      	       C
c to:           output ocean temperature.                              C
c sa:           output ocean salinity.                                 C
c uo:           output ocean u velocity.                               C
c vo:           output ocean v velocity.                               c
c wo:           output ocean w velocity.                               c
c ro:           output ocean density.                                  c
c hadvt:        output ocean horizontal advection of temperature.      c
c hadvs:        output ocean horizontal advection of salinity.         c
c vadvt:        output ocean vertical advection of temperature.        c
c vadvs:        output ocean vertical advection of salinity.           c
c hdift:        output ocean horizontal diffusion of temperature.      c
c hdifs:        output ocean horizontal diffusion of salinity.         c
c vdift:        output ocean vertical diffusion of temperature.        c
c vdifs:        output ocean vertical diffusion of salinity.           c
c convt:        output ocean convective adjustment of temperature.     c
c convs:        output ocean convective adjustment of salinity.        c
c relaxt:       output ocean temperature forceing.                     c
c relaxs:       output ocean salinity forceing.                        c
c brin:         output fresh water release or update due to ice change.c
c dice:         output depth of ice.                                   c
c sds:          output snow depth over seaice                          c
c tice:         output ice temperature.                                c
c heico:        output heat flux between ice and ocean.                c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C** Here following are default values for above mentioned parameters.
C** These parameters can be updated in the namelist.


      isnapshot = 0
      ifreq    = 30

      do i = 1, 20
        to(i)      = 0
        sa(i)      = 0
        uo(i)      = 0
        vo(i)      = 0
        wo(i)      = 0
        ro(i)      = 0
        hadvt(i)   = 0
        hadvs(i)   = 0
        vadvt(i)   = 0
        vadvs(i)   = 0
        hdift(i)   = 0
        hdifs(i)   = 0
        vdift(i)   = 0
        vdifs(i)   = 0
        hdift(i)   = 0
        convt(i)   = 0
        convs(i)   = 0
        relaxs(i)  = 0
        relaxt(i)  = 0
        heex(i)    = 0
        saex(i)    = 0
        brin(i)    = 0
        dice(i)    = 0
        sds(i)     = 0
        tice(i)    = 0
        heico(i)   = 0
      enddo

      read(15,NML = outoctl)
      read(15,NML = wrocpar)

      do i = 1, 20
        nto(i)      = to(i)
        nsa(i)      = sa(i)
        nuo(i)      = uo(i)
        nvo(i)      = vo(i)
        nwo(i)      = wo(i)
        nro(i)      = ro(i)
        nhadvt(i)   = hadvt(i)
        nhadvs(i)   = hadvs(i)
        nvadvt(i)   = vadvt(i)
        nvadvs(i)   = vadvs(i)
        nhdift(i)   = hdift(i)
        nhdifs(i)   = hdifs(i)
        nvdift(i)   = vdift(i)
        nvdifs(i)   = vdifs(i)
        nhdift(i)   = hdift(i)
        nconvt(i)   = convt(i)
        nconvs(i)   = convs(i)
        nrelaxt(i)  = relaxt(i)
        nrelaxs(i)  = relaxs(i)
        nheex(i)    = heex(i)
        nsaex(i)    = saex(i)
        nbrin(i)    = brin(i)
        ndice(i)    = dice(i)
        nsds(i)     = sds(i)
        ntice(i)    = tice(i)
        nheico(i)   = heico(i)
      enddo

c *** initialise counters

      iss   = 0
      icont = 0

      return
      end
