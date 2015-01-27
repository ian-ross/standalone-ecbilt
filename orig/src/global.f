c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iniglobal
c-----------------------------------------------------------------------
c *** initialisation of climate model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comatfor.h'

      real*8  ari(nlat/2)
      real*8  clfrcs(nlat,2),rj,rw,dumphi,dumwei
      integer indu(0:nlon,0:nlat),indt(0:nlon-1,0:nlat-1)
      integer ilat,k,irn,i,j


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

      dp1=5.0d+4

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

      undef = 9.99E10

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

c *** nstpyear is number of atmospheric timesteps per year
c *** nocstpyear is number of ocean timesteps per year
c *** ntstep is total number of timesteps
c *** nbclins is number of atmospheric time steps per baroclinic ocean
c *** timestep
c *** nbtrops is number of atmospheric time steps per barotropic ocean
c *** timestep


      nstpyear   = iatm*360
      nocstpyear = 360/idtbclin
      ntstep     = nstpyear*nyears
      nbclins    = iatm*idtbclin
      nbtrops    = iatm*idtbtrop

c *** time step of the atmospheric model:
c *** dt    : fraction of one day
c *** dtime : in seconds
c *** dtt   : dimensionless

      dt     = 1d0/dble(iatm)
      dtime  = dt*(24d0*3600d0)
      dtt    = dt*pi*4d0

c *** gauss points and weights

      ilat=nlat/2
  10  continue
        read(7,220,end=15) rj,rw
        irn=int(rj)
        if (irn.eq.ilat) then
          do i=1,irn
            read(7,220) ari(i),dumwei
          enddo
          goto 20
        else
          goto 10
        endif
  15    continue
        call error(4)
  20  continue

      do i=1,ilat
        phi(i)=-ari(ilat+1-i)
        phi(ilat+i)=ari(i)
      enddo

      do i=1,nlat
        phi(i)=asin(phi(i))
      enddo

      dphi=(phi(nlat)-phi(ilat+1))/(ilat-1)
      dlab=pi/dble(nlon/2)
      darea=dphi*dlab
      tarea=0.

      do i=1,nlat
        cosfi(i)=cos(phi(i))
        sinfi(i)=sin(phi(i))
        tanfi(i)=tan(phi(i))
        cosfid(i)=cosfi(i)*darea
        sinfid(i)=sinfi(i)*darea
        do j=1,nlon
          tarea=tarea + cosfid(i)
        enddo
      enddo
c
c *** landsea mask
c
      read (40,120)
      do j=nlat,0,-1
        read(40,120) (indu(i,j),i=0,nlon)
      enddo
      do j=nlat-1,0,-1
        read (40,310) k,(indt(i,j),i=0,nlon-1)
      enddo
      do i=1,nlat
        do j=1,nlon
          lsmask(i,j)=indt(j-1,i-1)
        enddo
      enddo
      rewind(40)

c
c *** climatological cloud fractions
c
      do i=1,nlat
        read(8,100) clfrcs(i,1),clfrcs(i,2)
      enddo
      do i=1,nlat
        clfrac(i)=0.5*(clfrcs(i,1)+clfrcs(i,2))
      enddo

      do j=1,nlon
        do i=1,nlat
          otempsg(i,j)=0.d0
          odlrads(i,j)=0.d0
          ouv10(i,j)=0.d0
          owinstu(i,j)=0.d0
          owinstv(i,j)=0.d0
          oq10(i,j)=0.d0
          otorain(i,j)=0.d0
          orunofo(i,j)=0.d0
          oevap(i,j)=0.d0
        enddo
      enddo

310   format(i4,i2,90i1)
120   format(65i1)
220   format(f18.10,f17.10)
100   format(2f5.2)


      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inimdldate
C--------------------------------------------------------------------------------
C ***
C *** This routine initialises the day, month, year of the model run
C ***
C-------------------------------------------------------------------------------
      implicit none

      include 'comglobal.h'

      day    = 0
      iyear  = 1
      imonth = 1
      iday   = 0
      itel   = 0
      iseason= 0
      instcount=0

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine mdldate(istep)
C--------------------------------------------------------------------------------
C ***
C *** This routine calculates the day, month, year of the model run from
C *** integration steps.
C *** Written by xueli wang April 1995.
C ***
C-------------------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comglobal.h'

      integer iy,im,idate,istep


      day = (mod(istep-1,nstpyear)) * dt

      if (mod(istep-1,iatm).eq.0) then
        iday = iday + 1
        if (iday .gt. 30) then
          iday = 1
          imonth = imonth + 1
          if (imonth.gt.12) then
            imonth = 1
            iyear = iyear + 1
          endif
        endif
      endif


      iy = iyear * 10000
      im = imonth * 100

      idate = 20000000 + iy + im +iday

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine checks(istep)
c-----------------------------------------------------------------------
c *** this routine performs checks on the model formulation
c-----------------------------------------------------------------------
      implicit none

      include 'comglobal.h'

      integer i,j,istep

      if ( mod(istep,iatm) .eq. 0) then
        call testocean(istep)
        call test(istep)
      endif

c      call checkbm
c      call moiscntrl(istep)

      call flush(100)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine checkbm
c-----------------------------------------------------------------------
c *** this routine performs checks on the model formulation
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'
      include 'comphys.h'
      include 'comland.h'

      integer i,j

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.0) then
            if (bmoisg(i,j).eq.0d0.and.evap(i,j).gt.0d0) then
              write(100,*) 'error in land moisture '
              write(100,*) 'lonlat ',i,j,bmoisg(i,j),evfaca(i,j)
            endif
          endif
        enddo
      enddo

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine test(istep)
c-----------------------------------------------------------------------
c *** testing if model variables are inside prescribed ranges
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'

      integer i,j,istep
      real*8  tsurfmean,globalmean
      character*12 chts


      do j=1,nlon
        do i=1,nlat

          if (uv10(i,j).gt.60) then
            write(100,*) 'uv10 out of range'
            write(100,*) i,j,uv10(i,j)
          endif
          if (tsurf(i,j).gt.400.or.tsurf(i,j).lt.150) then
            write(100,*) 'tsurf out of range in test'
            write(100,*) i,j,tsurf(i,j)
          endif
          if (eflux(i,j).gt.2000.or.hflux(i,j).gt.2000) then
            write(100,*) 'surface flux out of range in test'
            write(100,*) i,j,eflux(i,j),hflux(i,j)
          endif

        enddo
      enddo

      tsurfmean=globalmean(tsurf)-tzero

      write(chts,900) tsurfmean
 900  format(E12.5)
      if (chts(1:3).eq.'nan') call error(99)

      write(20,110) iyear,int((day+0.5*dt)/(iatm*dt))+1,tsurfmean
      call flush(20)

      if (tsurfmean.gt.40..or.tsurfmean.lt.-10.) then
        write(29,*) 'mean surface temperature ',tsurfmean
        call error(3)
      endif

  110 format(i8,i8,f7.2)
  100 format(f7.2)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine testocean(istep)
c-----------------------------------------------------------------------
c *** testing if model variables are inside prescribed ranges
c-----------------------------------------------------------------------
      implicit none

      include 'comocean.h'
      include 'comglobal.h'
      include 'comlake.h'

      integer      i,j,k,istep
      character*12 chts


      do j=-1,nm+1
        do i=-1,nl+1
          if (iwater(i,j).eq.1) then
            do k=1,nzl
              if (te(i,j,k).lt.-4.or.te(i,j,k).gt.45) then
                write(100,*) 'T ocean out of range ',istep
                write(100,10) i,j,k,te(i,j,k)
              endif
            enddo
            if (isea(i,j).eq.1) then
              do k=nzl+1,nz
                if (te(i,j,k).lt.-4.or.te(i,j,k).gt.45) then
                  write(100,*) 'T ocean out of range ',istep
                  write(100,10) i,j,k,te(i,j,k)
                endif
              enddo
            endif
          endif
        enddo
      enddo

 10   format(3i4,E14.7)
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moiscntrl(istep)
c-----------------------------------------------------------------------
c *** this routine controls on the conservation of the moisture budget
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'
      include 'comphys.h'
      include 'comland.h'
      include 'comcontrl.h'

      integer i,j,istep
      real*8 iflux,oflux,labufch,labufcht,moisbdif
      real*8 totrunl,totruno,difrun
      real*8 ifluxt,ofluxt,abufchn,abufchnt,facmois
      real*8 difatmois,obufchnt,ocbud

      if (istep.gt.2) then

c *** local control for surface moisture budget

      labufcht=0.
      do j=1,nlon
        do i=1,nlat
           if (lsmask(i,j).eq.0) then
             iflux=toraino(i,j)-evapo(i,j)
             oflux=runofl(i,j)
             labufch=(bmoisg(i,j)+dsnow(i,j))-(bmoisgo(i,j)+dsnowo(i,j))
             labufcht=labufcht + labufch*cosfi(i)
             moisbdif=labufch - dtime*(iflux-oflux)
               if (abs(moisbdif).gt.1d-5) then
                 write(100,*) 'error in landmoisture budget'
                 write(100,*) 'lonlat',i,j,moisbdif
               endif
           endif
         enddo
       enddo

c *** global control for runoff budget

       totrunl=0.
       totruno=0.
       do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.0) then
            totrunl=totrunl + runofl(i,j)*darea*cosfi(i)
          else
            totruno=totruno + runofo(i,j)*darea*cosfi(i)
          endif
        enddo
       enddo

       difrun=totrunl-totruno

       if (abs(difrun).gt.1d-5) then
         write(100,*) 'error in runoff budget'
         write(100,*) difrun
       endif

c ***  atmospheric moisture budget

       ifluxt=0.
       ofluxt=0.
       abufchnt=0.
       do j=1,nlon
         do i=1,nlat
           ifluxt=ifluxt + evapo(i,j)*cosfi(i)
           ofluxt=ofluxt + toraino(i,j)*cosfi(i)
           abufchnt=abufchnt + (rmoisg(i,j)-rmoisgo(i,j))*cosfi(i)
         enddo
       enddo

       difatmois=abufchnt - dtime*(ifluxt-ofluxt)

       if (abs(difatmois).gt.1d-5) then
         write(100,*) 'error in atmospheric budget'
         write(100,*) difatmois
       endif

c ***  budget of total hydrological cycle
c ***  net moisture flux to the oceans is due to
c ***  change in total moisture content of the atmosphere
c ***  and land surface

       obufchnt=0.
       do j=1,nlon
         do i=1,nlat
          if (lsmask(i,j).eq.1) then
            obufchnt=obufchnt +
     *      dtime*(runofo(i,j)+torain(i,j)-evap(i,j))*cosfi(i)
          endif
         enddo
       enddo

       ocbud=abufchnt + labufcht + obufchnt


       write(100,*) 'moisture budgets'
       write(100,*) 'atmos',abufchnt,'land', labufcht,'ocean',obufchnt
       write(100,*) 'difference',ocbud

       endif


c *** storing of moisture variables


       do j=1,nlon
         do i=1,nlat
           toraino(i,j)=torain(i,j)
           evapo(i,j)=evap(i,j)
           runoflo(i,j)=runofl(i,j)
           runofoo(i,j)=runofo(i,j)
           bmoisgo(i,j)=bmoisg(i,j)
           dsnowo(i,j)=dsnow(i,j)
           rmoisgo(i,j)=rmoisg(i,j)
         enddo
       enddo


       return
       end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine outato(istep)
c-----------------------------------------------------------------------
c *** saving atmospheric data for forcing the ocean only model
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'
      include 'comatfor.h'
      include 'comcoup.h'

c *** outato stores the daily mean atmospheric fields
c *** for driving the ocean
c *** ntboa: start of the period
c *** nteoa: end of the period

      integer     i,j,istep,ntboa,nteoa
      real*8      facstr,riatm
      character*2 fnoc

      facstr = cdrag*roair*uv10rws

      ntboa=nstpyear*(nbsatfor-1)
      nteoa=nstpyear*(nbsatfor+nafyear-1)


      if (istep.gt.ntboa.and.istep.le.nteoa) then

        if (mod(istep,nstpyear).eq.1) then
          if (istep.gt.ntboa+1) then
            close (85)
          else
            nocfile=0
          endif
          nocfile=nocfile+1
          write(fnoc,10) nocfile
   10     format(i2.2)
          if (nocfile.le.nafyear) then
            open(85,file='outputdata/ocean/atofor'//fnoc//'.dat',
     *      form='unformatted')
          endif
        endif

        do j=1,nlon
          do i=1,nlat
            otempsg(i,j)=otempsg(i,j)+ tempsg(i,j)
            odlrads(i,j)=odlrads(i,j)+ dlrads(i,j)
            ouv10(i,j)  =ouv10(i,j)  + uv10(i,j)
            owinstu(i,j)=owinstu(i,j)+ facstr*uvw10(i,j)*utot(i,j,3)
            owinstv(i,j)=owinstv(i,j)+ facstr*uvw10(i,j)*vtot(i,j,3)
            oq10(i,j)   =oq10(i,j)   + q10(i,j)
            otorain(i,j)=otorain(i,j)+ torain(i,j)
            orunofo(i,j)=orunofo(i,j)+ runofo(i,j)
            oevap(i,j)  =oevap(i,j)  + evap(i,j)
          enddo
        enddo

        if (mod(istep,iatm).eq.0) then

          riatm=1d0/iatm
          do j=1,nlon
            do i=1,nlat
              otempsg(i,j)=otempsg(i,j)*riatm
              odlrads(i,j)=odlrads(i,j)*riatm
              ouv10(i,j)  =ouv10(i,j)  *riatm
              owinstu(i,j)=owinstu(i,j)*riatm
              owinstv(i,j)=owinstv(i,j)*riatm
              oq10(i,j)   =oq10(i,j)   *riatm
              otorain(i,j)=otorain(i,j)*riatm
              orunofo(i,j)=orunofo(i,j)*riatm
              oevap(i,j)  =oevap(i,j)  *riatm
            enddo
          enddo

          write(85) ((otempsg(i,j),j=1,nlon),i=1,nlat)
          write(85) ((odlrads(i,j),j=1,nlon),i=1,nlat)
          write(85) ((ouv10(i,j)  ,j=1,nlon),i=1,nlat)
          write(85) ((owinstu(i,j),j=1,nlon),i=1,nlat)
          write(85) ((owinstv(i,j),j=1,nlon),i=1,nlat)
          write(85) ((oq10(i,j)   ,j=1,nlon),i=1,nlat)
          write(85) ((otorain(i,j),j=1,nlon),i=1,nlat)
          write(85) ((orunofo(i,j),j=1,nlon),i=1,nlat)
          write(85) ((oevap(i,j)  ,j=1,nlon),i=1,nlat)

          do j=1,nlon
            do i=1,nlat
              otempsg(i,j)=0d0
              odlrads(i,j)=0d0
              ouv10(i,j)  =0d0
              owinstu(i,j)=0d0
              owinstv(i,j)=0d0
              oq10(i,j)   =0d0
              otorain(i,j)=0d0
              orunofo(i,j)=0d0
              oevap(i,j)  =0d0
            enddo
          enddo

        endif
      endif

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inato(istep)
c-----------------------------------------------------------------------
c *** read output of a previous coupled run for forcing ocean model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comland.h'
      include 'comatfor.h'
      include 'comcoup.h'

      integer     i,ifail,j,g05dyf,istep
      character*2 fnoc

      if (mod(istep,nstpyear).eq.0) then
        if (istep.gt.0) close(85)
        nocfile=g05dyf(1,nafyear)
        write(fnoc,10) nocfile
   10   format(i2.2)
        open(85,file='outputdata/ocean/atofor'//fnoc//'.dat',
     *          form='unformatted')
      endif

      if (mod(istep,iatm).eq.0) then

        read(85) ((otempsg(i,j),j=1,nlon),i=1,nlat)
        read(85) ((odlrads(i,j),j=1,nlon),i=1,nlat)
        read(85) ((ouv10(i,j)  ,j=1,nlon),i=1,nlat)
        read(85) ((owinstu(i,j),j=1,nlon),i=1,nlat)
        read(85) ((owinstv(i,j),j=1,nlon),i=1,nlat)
        read(85) ((oq10(i,j)   ,j=1,nlon),i=1,nlat)
        read(85) ((otorain(i,j),j=1,nlon),i=1,nlat)
        read(85) ((orunofo(i,j),j=1,nlon),i=1,nlat)
        read(85) ((oevap(i,j)  ,j=1,nlon),i=1,nlat)

        do j=1,nlon
          do i=1,nlat
            tempsg(i,j) =otempsg(i,j)
            dlrads(i,j) =odlrads(i,j)
            uv10(i,j)   =ouv10(i,j)
            winstua(i,j)=owinstu(i,j)
            winstva(i,j)=owinstv(i,j)
            q10(i,j)    =oq10(i,j)
            torain(i,j) =otorain(i,j)
            runofo(i,j) =orunofo(i,j)
            evap(i,j)   =oevap(i,j)
          enddo
        enddo

        call solar
        call albedo
        call swaverad

        call ocfluxes

      endif


      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine outcheck
c-----------------------------------------------------------------------
c *** output check file
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'
      include 'comcoup.h'
      include 'comphys.h'

      integer  i,j,istep
      real*8   sstm,sssm,hefxm,safxm,tcosf

      sstm=0d0
      sssm=0d0
      hefxm=0d0
      safxm=0d0
      tcosf=0d0

      do j=1,nlon
        do i=1,nlat
          if (lsmask(i,j).eq.1) then

            sstm = sstm + stmix(i,j) * cosfid(i)
            sssm = sssm + samix(i,j) * cosfid(i)
            safxm = safxm + safxo(i,j) * cosfid(i)
            hefxm = hefxm + hefxo(i,j) * cosfid(i)
            tcosf = tcosf + cosfid(i)

          endif
        enddo
      enddo

      sstm =sstm/tcosf
      sssm =sssm/tcosf
      hefxm=hefxm/tcosf
      safxm=safxm/tcosf

      open(28,file='checkrun'//fini)

      write(28,900) 'global mean sea surface temperature : ',sstm
      write(28,900) 'global mean sea surface salinity    : ',sssm
      write(28,900) 'global mean air-sea heat exchange   : ',hefxm
      write(28,900) 'global mean air-sea water exchange  : ',safxm
      write(28,*)
      write(28,*) 'meridional profile at 180 degrees'
      write(28,*)
      write(28,910) '    SST    ','    SSS    ',' Heat flux '
      write(28,*)
      do i=1,nlat
        write(28,920) stmix(i,35),samix(i,35),hefxo(i,35)
      enddo

      close(28)

 900  format(A38,E14.8)
 910  format(3A11)
 920  format(3F10.4)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine writestate(istep)
c-----------------------------------------------------------------------
c *** this routine writes the current state of ecbilt to datafiles
c-----------------------------------------------------------------------
      implicit none

      include 'comglobal.h'

      integer istep
      character*4 chf

      if (mod(istep,nstpyear).eq.0) then
        if (mod(iyear,nwrskip).eq.0.or.iyear.eq.nyears) then
          write(chf,1) iyear+irunlabel
    1     format(i4.4)

          if (irunatm.eq.1) then
            open(95,file='startdata/inatmos'//chf//'.dat'
     *           ,form='unformatted')
            open(96,file='startdata/inland'//chf//'.dat'
     *           ,form='unformatted')
            call wrendatmos
            call wrendland
            close(95)
            close(96)
          endif
        endif
      endif

      return
      end
