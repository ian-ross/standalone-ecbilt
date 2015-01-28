












      subroutine defcst(nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C Initialisation of the various physical variables and parameters
c  nn99, : 1 => open "mouchard" and nn99 <- 2
c   0 => if kfond=-1,-3 , open "mouchard" and nn99 <- 2
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 23/08/02
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'reper.com'
      include 'varno.com'
      include 'comunit.h'
 
c--local variables :
      dimension coef2(kmax), coef3(kmax)
      dimension coef4(kmax), coef5(kmax), coef6(kmax), coef7(kmax)
      dimension coef8(kmax), coef9(kmax)
      dimension coef10(kmax),coef11(kmax)
C     dimension zrr(5),zd1(5),zd2(5)
Cic0  dimension acrit(2), hgcrit(2)
 
      character*70 line
      real  tmpzfluxm,tmpzfluxms
      integer ios
c--instructions "data" :
      include 'datadc.com'
      real*8 moc,tmc,tmc0,tsurfmean,cland,thex
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Reading of the parameter of the run ; arrays of associated parameters |
c-----------------------------------------------------------------------
 
       open(unit=10,file='run.param',status='old')
c- header 1 :
        read(10,*)
        read(10,*)
c- Checking at one particular point :
        icheck = 0
        jcheck = 0
        kcheck = 0
        read(10,'(A)') line
        if (line(14:14).ne.' ') read(line(15:),*) icheck, jcheck, kcheck
Cray    read(10,'(15x,3i5)') icheck, jcheck, kcheck
c- header 2 :
      do 100 n=1,4
        read(10,*)
 100  continue
c- option of the code :
c---
c kfond > 0 => Nb_niv.=kfond (flat); < -2 => DownSloping , -1,-3 => mouchard
c lstab > 1 => AjConv by pairs, lstab times ; -1 total mixing; 0 : nothing
c        -3 =>  permutation & mixing of a fraction ajcmix
      read(10,*)
      read(10,*) kfond, lstab, nsewfr, bering, ajcmix, xslop
c- initial value (by level) for each scalar :
      do 105 n=1,nsmax
        read(10,*)
        read(10,*) (scal0(k,n),k=kmax,1,-1)
 105  continue
c- conversion Celcius -> Kelvin :
      do 110 k=1,kmax
        scal0(k,1) = scal0(k,1) + 273.15d0
 110  continue
c- restoring :
      nitrap = 0
      ahrap = 0.
      read(10,*)
      read(10,*) unstyr, (rapp0(k),k=0,(nsmax+2)), ahrap, nitrap
      read(10,*)
      read(10,*) (rapp1(k),k=kmax,1,-1)
c- time steps :
      read(10,*)
      read(10,*) dts(kmax), dtu, dtb
      read(10,*)
      read(10,*) (coef2(k),k=kmax,1,-1)
c- Parameters associated with the coupling with an atmospheric model
c icoupl=0 : not coupled
c icoupl=1 : following of a preceding run
c icoupl=2 : set time to zero
c icoupl=3 : Beginning of a coupled run with starting state
c Coupling type
c icoutp=1 : LMD 5bis
c icoutp=2 : LMD Z
c icoutp=3 : ECBILT
c icoutp=4 : (to be continued)
 
      read(10,*)
      read(10,*) icoupl, icoutp, itau_slow
c-----
c- nb-iterations :
      read(10,*)
      read(10,*) nitrun, nsplit, nsplaj
      read(10,*)
      read(10,*) nclin, nclmoy
c- option start/read/write :
      read(10,*)
      read(10,*) kstart, kinput, koutpu
c- frequ. writing :
      read(10,*)
      read(10,*) nsav, ninfo, ntmoy,ntout
c- option for the sea-ice ...
c       nwtest=0, WRITE OUTPUTS OF YEARS:
c                              nw(..), 2*nw(..), ...
c       nwtest=1, WRITE OUTPUTS OF YEARS:
c                              1, nw(..)+1, 2*nw(..)+1, ...
      read(10,*)
      read(10,*) nwjl, nwm, nwa, nwtoa, nwtom
      read(10,*)
      read(10,*) nwtest, nwtal
 
c- Bottom drag coefficient
      read(10,*)
c     read(10,*) cdbot
      read(10,*) cdbot, txifcb, txifcu
 
c- diffusivity & visc. H. :
      read(10,*)
      read(10,*) ahs(kmax), ahu, ahe
      read(10,*)
      read(10,*) (coef3(k),k=kmax,1,-1)
c- diffusivity & visc. V. :
      read(10,*)
      read(10,*) avnub(kmax), avnu0(kmax), rifumx
      read(10,*)
         read(10,*) avkb(kmax),  avk0(kmax),  rifsmx
      read(10,*)
      read(10,*) (coef4(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) (coef5(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) (coef6(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) (coef7(k),k=kmax,1,-1)
c- upwind rate :
      read(10,*)
      read(10,*) alphxu, alphxv, alphyu, alphyv
      read(10,*)
      read(10,*) alphah(1), alphah(2),
     &           (alphgr(ns),ns=1,nsmax), (algrmn(ns),ns=1,nsmax)
      read(10,*)
      read(10,*) (alphmi(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) (alphaz(k),k=kmax,1,-1)
c- reading of the parameters for the formulation of Avs(N2)
      read(10,*)
      read(10,*) kavsmx, qavs, avsn2, ccfmn, ccfmx
c- rate Implic. / Explic.
      read(10,*)
      read(10,*) txiadu, txiads, txidfu, txidfs, txeflx
c- modif forcing :
      read(10,*)
      read(10,'(A)') line
Cray  read(10,'(I4,6X,A40)') mdforc,filcor
c- which forcing files read :
c kforc : (unite-> constant, dizaine-> saisonnier, centaine ...)
      read(10,*)
      read(10,*) kforc, (yforc(ns),ns=0,nsmax)
 
c- value of precipitations and value at surface for the scalars
      read(10,*)
      read(10,*) (scpme(ns),ns=1,nsmax)
      read(10,*)
      read(10,*) (scssv(ns),ns=1,nsmax)
 
c- constants and parameters for the turbulence model.
      read(10,*)
      read(10,*) q2tmin,zlotur,vlmin,varfor,kajul
      vkappa  = 0.4d0
      ghmax  = 0.0233d0
      ghmin  = -0.28d0
      sqrghm = sqrt(-ghmin)

c- volume correction + bottom heat flux + hysteresis or not
      read(10,*)
      read(10,*) vcor,bheat,ihyster
 
c- Parameters for isoslope.f, isodiffu.f, isoadvec.f
      read(10,*)
      read(10,*) ai(kmax),slopemax(kmax)
      read(10,*)
      read(10,*) (coef8(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) (coef9(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) aitd(kmax),slopmgm(kmax),afilt,ahh,avv
      read(10,*)
      read(10,*) (coef10(k),k=kmax,1,-1)
      read(10,*)
      read(10,*) (coef11(k),k=kmax,1,-1)
 
      close(10)

      xfreshsv(:)=0.0
      if (ihyster.ne.0) then
        open(10,file='inputdata/fresh_for.dat')
         do i=1,5000
           read(10,*) jj,xfreshsv(i)
         enddo
      endif

      close(10)
   
 
      icheck = min(imax,icheck)
      jcheck = min(jmax,jcheck)
      kcheck = min(kmax,kcheck)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--1.2 computation of the other parameters :
c--------------
 
      do 120 k=1,kmax
        dts(k) = coef2(k) * dts(kmax)
        ahs(k) = coef3(k) * ahs(kmax)
        avnub(k) = coef4(k) * avnub(kmax)
        avnu0(k) = coef5(k) * avnu0(kmax)
           avkb(k)  = coef6(k) * avkb(kmax)
        avk0(k)  = coef7(k) * avk0(kmax)
        ai(k)=coef8(k)*ai(kmax)
        slopemax(k)=coef9(k)*slopemax(kmax)
        aitd(k)=coef10(k)*aitd(kmax)
        slopmgm(k)=coef11(k)*slopmgm(kmax)
 120  continue

c- unchanged values unchanged at the bottom :
        avnub(1) = coef4(1)
        avnu0(1) = coef5(1)
           avkb(1)  = coef6(1)
        avk0(1)  = coef7(1)
c-
      if(nsplit.le.nsplaj.or.nsplaj.lt.0) then
        write(iuo+66,*) 'STOP in "defcst", nsplit & nsplaj =',
     &             nsplit, nsplaj
        stop
      endif
      if(nsewfr.eq.0) nsewfr = 1
      nn = 500 * nsav
      if (nn.ge.1.and.nn.le.numit) then
        write(iuo+66,*) 'STOP because Pb for parameter "nsav" :'
     &           //' too frequent output !'
        stop
      endif
      if ( koutpu.ge.2 ) then
c- test for the computation of the averages of Ajust.Conv. :
        if (ninfo.eq.0) then
          write(iuo+66,*) 'Defcst : modif of the parameter "ninfo" :'
          write(iuo+66,*) 'Old , New :', ninfo, nsav
          ninfo = nsav
        elseif ( mod(nsav,ninfo).ne.0 ) then
          nninfo = abs(ninfo)
          do 130 nn=nninfo,nsav
            if (mod(nsav,nn).eq.0) goto 135
 130      continue
          nn = nsav
 135      continue
          nn = sign(nn,ninfo)
          write(iuo+66,*) 'Defcst : modif of the parameter "ninfo" :'
          write(iuo+66,*) 'Old , New :', ninfo, nn
          ninfo = nn
        endif
      endif
 
c--Bering Strait opening:
      ibera = 0
      iberp = 0
      jbera = jmax
      jberp = jsepar
      if ( ltest.eq.3 .and. bering.ne.0. ) iberp = 1
 
c--Modificarion of the forcaging if needed :
      read(line,*) mdforc
      if (mdforc.ne.0) read(line,'(10x,A40)') filcor
 
c--Opening of the file. "mouchard" (=unit 99) and modification of "nn99" for the outputs
      if (nn99.eq.0 .and. mod(kfond,2).eq.-1 ) nn99 = 1
      if(nn99.eq.1) then
        nn99 = 2
        open(99,file='mouchard',status='unknown')
      endif
 
c--Initialisation of the arrsy "nvrl" (depending on  "koutpu"),
c   nvrl(nv)=1 => write the variable "nv" in the ouput/restart file.
      nvrl(nvret ) = 1
      nvrl(nvrub ) = 1
      nvrl(nvrvb ) = 1
      nvrl(nvru  ) = 1
      nvrl(nvrv  ) = 1
      nvrl(nvrtke) = 1
      do 170 ns=1,nsmax
        nvrl(ns  ) = 1
        if (koutpu.ge.2) nvrl(ns+nvrfw) = 1
 170  continue
      if (koutpu.ge.2) then
        nvrl(nvrfw ) = 1
        nvrl(nvrajc) = 1
C       nvrl(nvrb ) = 1
C       nvrl(nvrn2 ) = 1
        nvrl(nvras ) = 1
        nvrl(nvrau ) = 1
      endif
      nvrl(nvrum)  = 1
      nvrl(nvrvm)  = 1
      nvrl(nvrhg)  = 1
      nvrl(nvrfq)  = 1
      nvrl(nvrqs)  = 1
      nvrl(nvrhn)  = 1
      nvrl(nvral)  = 1
      nvrl(nvrts)  = 1
      nvrl(nvrug)  = 1
      nvrl(nvrvg)  = 1
      nvrl(nvrtbq) = 1
      nvrl(nvrxzo) = 1
      nvrl(nvrtgx) = 1
      nvrl(nvrtgy) = 1
      nvrl(nvrmom) = 1
      do 180 nv=1,nvmax
        if (ltyp(nv).eq.99) nvrl(nv) = 0
 180  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Initiailisation of   the numerical/physical constants  |
c-----------------------------------------------------------------------
 
      zero = 0.d0
      one = 1.d0
      epsil = 0.1d-9
      cstmin = -1.d20
      cstmax =  1.d20
      untour = 360.d0
 
      pi = 4.0 * atan(one)
      radian = pi / 180.0
      omega = 2.0 * pi / 86164.0
      rterre = 6371000.0
      unsrt = 1.0 / rterre
      gpes = 9.80d0
      rho0 = 1030.0
      svrdrp = 1.0d-6
      cpo = 4002.
c--number of day in a year :
Cic0  yeaday = 365.25d0
      yeaday = 365.0
      yeaday = 360.0
 
c--used in redforc, informe : change unit(model)  -> standard  unit 
      unitfx(0) = 1.
      do 290 ns=1,nsmax
        scalwr(ns) = 0.
        unitfx(ns) = 1.
 290  continue
c- unitfx(0) = Year(s) ;
c- unitfx(1) = rho.Cp : Fx. en W/m2 ; unitfx(2) = Year / 34.7 g/l : Fx. en m/y
      unitfx(0) = yeaday * 86400.0
      unitfx(1) = -4.1d+06
      unitfx(2) =  9.1d+05
 
c- initialisation of the bottom stress  <- moved to "defgrid"
 
      if (nn99.eq.2) write(99,*) 'year : yeaday =', yeaday
 
Cic0  return
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Values linked with the sea-ice thermodynamics   |
c-----------------------------------------------------------------------
 
c--3.1 Physical constants.
c--------------------------
 
      tfsn=273.15d0
      tfsg=273.05d0
      xkn=0.22d0
      xkg=2.034396d0
      rcpn=6.9069e+05
      rcpg=1.8837e+06
      xlg=300.33e+06
      xln=110.121e+06
      xsn=2.8e+06
      rhog=900.0
      rhon=330.0
c
C     emig=0.97d0
C     emig=1.0d0
C     sglace=4.
      sglace=6.
c
      rhoesn=1/rhon
c
c--3.2. Physical constants for the heat fluxes.
c-----------------------------------------------
c
      stefan=5.67e-08
c
      too=273.16d0
      vkarmn=0.4d0
      cevap=2.5e+06
      zemise=0.97d0
      zemise=1.00d0
c
c--3.3.Physical parameters.
c--------------------------
c
      open(53,file='thermo.param', status='old')
      read(53,*)
      read(53,*)
      read(53,*)
      read(53,*)
      read(53,*)
      read(53,*) hmelt
      read(53,*)
      read(53,*) acrit(1)
      read(53,*)
      read(53,*) acrit(2)
      read(53,*)
      read(53,*) hgcrit(1)
      read(53,*)
      read(53,*) hgcrit(2)
      read(53,*)
      read(53,*) emig
c
c--3.4. Numerical parameters.
c----------------------------
c
      read(53,*)
      read(53,*) hgmin
      read(53,*)
      read(53,*) hndif
      read(53,*)
      read(53,*) hgdif
      read(53,*)
      read(53,*) hglim
      read(53,*)
      read(53,*) amax
      read(53,*)
      read(53,*) swiqst
c
      uscomi=1.0/(1.0-amax)
c
c--3.5. Numerical caracteristics.
c---------------------------------
c
      read(53,*)
      read(53,*)
      read(53,*) beta
      read(53,*)
      read(53,*) ddtb
      read(53,*)
      read(53,*) nbits
      read(53,*)
      read(53,*) parlat
      read(53,*)
      read(53,*) hakspl
      read(53,*)
      read(53,*) hibspl
      read(53,*)
      read(53,*) exld
      read(53,*)
      read(53,*) hakdif
      read(53,*)
      read(53,*) hth
      read(53,*)
      read(53,*) hnzst
      read(53,*)
      read(53,*) parsub
      read(53,*)
      read(53,*) alphs
c
      xkn=hakdif*xkn
      xkg=hakdif*xkg
      if ((hndif.gt.100.0).or.(hgdif.gt.100.0)) then
        cnscg = 0.0
      else
        cnscg = rcpn/rcpg
      endif
c
      close(53)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Values for ice-dynamic model.                                   |
c-----------------------------------------------------------------------
c
c--3.1. begin reading file of parameters.
c-----------------------------------------
c
      open(54,file='dynami.param',status='old')
c
      read(54,*)
      read(54,*)
      read(54,*)
      read(54,*)
      read(54,*)
c
c--3.2. numerical characteristics.
c---------------------------------
c
      read(54,*)
      read(54,*) idyn
      read(54,*)
      read(54,*) zepsd1
      read(54,*)
      read(54,*) zepsd2
      read(54,*)
      read(54,*) nlminn
      read(54,*)
      read(54,*) nlmaxs
      nlmaxn=jmax-1
      nlmins=3
      read(54,*)
      read(54,*) usdt
      read(54,*)
      read(54,*) alpha
      read(54,*)
      read(54,*) bound
      read(54,*)
      read(54,*) dm
      read(54,*)
      read(54,*) nbitdf
      read(54,*)
      read(54,*) nbiter
      read(54,*)
      read(54,*) nbitdr
      read(54,*)
      read(54,*) om
      read(54,*)
      read(54,*) resl
c
c--3.3. physical parameters.
c---------------------------
c
      read(54,*)
      read(54,*) cw
      read(54,*)
      read(54,*) angvg
      read(54,*)
      read(54,*) pstar
      read(54,*)
      read(54,*) c
      read(54,*)
      read(54,*) zetamn
      read(54,*)
      read(54,*) creepl
      read(54,*)
      read(54,*) ecc
      read(54,*)
      read(54,*) uvdif
      read(54,*)
      read(54,*) ren
      read(54,*)
      read(54,*) gridsz
      read(54,*)
      read(54,*) iameth
c
      close(54)
c
      usecc2 = 1.0/(ecc*ecc)
      rhoco  = rho0*cw
      rhoco2 = rhoco*rhoco
      angvg  = angvg*radian
      sangvg = sin(angvg)
      cangvg = cos(angvg)
      pstarh = pstar/2.0
      sber   = 1.-max(zero,sign(one,0.1-bering))
c
c- 3.4 Mean value of w
      
      zfluxm=0.0
      zfluxms=0.0
      if (vcor.gt.zero) then
        open(36,file='correcw.dat',form='formatted')
        ios=0
        do while (ios==0)
c S'il n'y a qu'une valeur sur la ligne => tmc & ios>0 => va ligne 974
        read(36,*,end=974,iostat=ios) tmpzfluxm,tmpzfluxms
        zfluxm=tmpzfluxm
        zfluxms=tmpzfluxms
       enddo
974    if (ios.NE.0)tmc0=tmpzfluxm  
       close(36)
       write(iuo+66,*) 'correction vol applied',vcor,zfluxm,zfluxms
      endif
     
c--3.4. Geostrophic velocities
c------------------------------
c
C     open (89,file='/u14/grpastr/hgs/data/geost/ugeost.mig',
C    &      form='unformatted',status='old')
C     read(89) uost
C     close(89)
C     open (89,file='/u14/grpastr/hgs/data/geost/vgeost.mig',
C    &      form='unformatted',status='old')
C     read(89) vost
C     close(89)
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- end of the routine defcst -
      end
