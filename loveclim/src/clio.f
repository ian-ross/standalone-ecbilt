












      subroutine init_CLIO
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 23/08/02
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'comunit.h'
c     include 'iceberg.com'
 
c- nn99=2 => writing in the file "mouchard", unit=99
      common / mchd99 / nn99
      common /clio_control_int/ ktvar,ntrmax,mixage
      common /clio_control_fp/ dtsd2,yrsec,daysec,unsplt
 
c-driess
      logical flgveg,flgicb,flgisma,flgismg

      common /ec_coupl/ flgveg,flgicb,flgisma,flgismg 

c--local variables :
      dimension irn(imax,8), jrn(jmax,8)
      character*4 fchnum
 
      open (iuo+66,file='clio3.out')
      write(iuo+66,*) 'clio3 : Adv. Alterne X,Y (nsewfr)'
Cic0  write(iuo+66,*) 'om : Adv. Alterne X,Y (nsewfr)'
 
      call foroutp(irn,jrn)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1) Preparation of the run.                                          |
c-----------------------------------------------------------------------
 
      jflag = 0
      nn99 = 0
      call defcst(nn99)
      call defgrid(nn99)
      call redforc(nn99)
      call initseaalb
c--Zonaly uniform, time dependant forcing :
      ktvar = abs(kforc) / 100
 
c--kstart = 0 start from routine start   (ocean at rest, no ice)
c--       = 1                    redrun* (follow up of a run, same conditions)
c--       > 1 transition :
c--       = 2                    staoc*  (ocean not at rest, sea ice prescribed)
c--kinput = 0/2 (nn2t=0) restart from the binairy file (*b)
c--       = 1/3 (nn2t=1)                  NetCDF       (*c)
c--koutpu = 0/2 (nn3t=0) output on a binary file (*b)
c--       = 1/3 (nn3t=1)             NetCDF  (*c)
 
      nn2t = mod(kinput,2)
      nn3t = mod(koutpu,2)
      ntrmax = int(dts(ks2)/ddtb)
      if ((ntrmax .ne. 1) .and. (icoupl.ge.1)) then
          write(iuo+66,*) 'problem with ntrmax in  coupled mode ',ntrmax
          stop
      endif
      if (kstart.eq.0) then
c- case kstart = 0 :
        call start
      elseif (kstart.eq.2) then
        if (nn2t.eq.0) then
          call staocb(kinput,'resto.om')
        else
C         call staocc(kinput,'resto.ncdf')
        endif
      else
        if (nn2t.eq.0) then
          call redrunb(kinput,nn99,'rest.om')
        else
Cncd      call redrunc(kinput,nn99,'rest.ncdf')
          write(iuo+66,*) 'om : Version without NetCDF '//
     &               '=> reading rest.ncdf impossible !'
          stop
        endif
      endif
 
      if (icoupl.ge.2)  then
        numit = 0
        tpstot = 0.0
      endif
      nstart = numit + 1
      nlast  = numit + nitrun
      yrsec  = yeaday*86400.
      daysec  = 86400.
      dtsd2 = 0.5*dts(ks2)
      iyear = 1+int(tpstot/yrsec)
      xday = mod(tpstot,yrsec)/daysec
      write(iuo+66,'(A,I8,A,F7.2,A,I10)') ' start of the run ; year=', iyear,
     &    ' +', xday, ' days ; iteration=', nstart
      if (ninfo.eq.0) ninfo = nlast + 1
      if (nsav.eq.0)  nsav  = nlast + 1
      mixage = max(min(lstab+lstab,2),-lstab)
      unsplt = 1.0 / DFLOAT(nsplit - nsplaj)
 
Csai  call splsclr1
      call etat(0, nn99)
      call informe(nn99)
      if (nstart.eq.1) call inforun(nn99)
C     call conti3d
 
c- Initialisation of the coupling
      if (icoupl.ge.1) then
Ccp1     call mytid_ocn
Ccp1     call ocn_rfm_cpl
         numcpl=0
Ccp2     call inicmo2(nitrun,freqcpl,int(dts))
c         call inicmo3(nitrun,freqcpl,int(dts))
c         write(iuo+66,*) ' ocn recv initi forc from cpl, itr :', nstart
      end if

c- Initialisation of the netcdf outputs
      call clio_inioutparat

c     if (flgicb) call deficeb
      end
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      subroutine CLIO(istep)
            
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'comunit.h'
c     include 'iceberg.com'
      
      character*4 fchnum
      common /clio_control_int/ ktvar,ntrmax,mixage
      common /clio_control_fp/ dtsd2,yrsec,daysec,unsplt
      
c-driess
      logical  flgveg,flgicb,flgisma,flgismg
      common / ec_coupl/ flgveg,flgicb,flgisma,flgismg

      
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2) Daily do loop.                                                   |
c-----------------------------------------------------------------------
 
        numit = numit + 1
c modif for  the counting of the interations in couple mode
        nsew = mod(numit,nsewfr)
c--The date corresponds to the middle of the time step.
c--iyear is the year, the first year is the year 1.
c--The integer part of xjour is the day of the year(1-yeaday) and the
c--non-integer one is the fraction of this day.
        tpstot = tpstot + dtsd2
        iyear  = 1+int(tpstot/yrsec)
        xjour  = mod(tpstot,yrsec)/daysec +1.0
        xjour1 = mod(tpstot+dts(ks2),yrsec)/daysec +1.0
        tpstot = tpstot + dtsd2
 
c-- 2.1. start of an iteration :
c-------------------------------
 
Cic0    if (ktvar.ne.0) call tvforc(nn99)
Csai    call splsclr2(jflag,iyear,xjour)
        if (icoupl.eq.0) then
           call forcat(iyear,xjour)
        elseif ((icoupl.eq.3) .and. (numit.eq.nstart)) then
Ccp2       call rfmcpl2(numcpl)
           call forcat(iyear,xjour)
C       elseif (icoupl.ge.1) then
        else
Ccp2       call rfmcpl2(numcpl)
c           call rfmcpl3(numcpl)
        endif
        call icdyna
        call icdadv(xjour)
c        call cfc_flx(nn99)
        do ither=1,ntrmax
           call thersf(ither,ntrmax)
Ccp1       if (icoupl.ge.1.and.ither.lt.ntrmax) then
Ccp1          call ocn_s2_cpl
Ccp1          call ocn_rfm_cpl
Ccp1       endif
        enddo
C       call ocesla
Ctk0    call vdiffu(nn99)
        call diftur
        call engtur
        call conti3d
        call flucor
        call isoslope
        if (aitd(kmax).gt.zero) call isoadvec

        call scale(nsew,nn99)
        call scali

        call etat(mixage, nn99)
        do 250 nuclin=1,nclin
          call uve
          call barot
          ccsplt = 0.0
          do 200 numspl=1,nsplit
            if (numspl.gt.nsplaj) ccsplt = unsplt
            if (ahe.eq.zero) then
              call uvb0et(ccsplt)
            else
              call uvbfet(ccsplt)
            endif
 200      continue
C         call uvi   ! --> ds uvm
          call uvm
 250    continue
C       write(99,*) 'clio',toticesm
 
c--partie ICEBERG de l'ITERATION
c--------------------
         

c       if (flgicb) then 
c        call iceberg(istep,iyear,nn99)
c       endif

c--end of an iteration 


c-- 2.2. Outputs.
c-----------------
        if (ntout.eq.1) then
        call outnet(iyear,xjour,xjour1)
        else
        call outave(iyear,xjour,xjour1)
        endif
        if (mod(numit,ninfo).eq.0 .or. ntmoy.ge.1) call inforun(nn99)
 
        if (mod(numit,nsav).eq.0) then
c--definition of the root f or the name of the restart file:
          indicf = (nlast - numit) / nsav
          ncfch = max(1,indicf)
          ncfch = 4 - int( log10( DFLOAT(ncfch) ) )
          write(fchnum,'(I4)') indicf
 
c--vertical velocity "w" n agreement with (u,v) :
          if (koutpu.ge.2) call conti3d
 
c--save of teh results
          nn3t = mod(koutpu,2) 
          if (nn3t.eq.1) then
Cncd        call savrunc(koutpu,nn99,'res'//fchnum(ncfch:)//'.ncdf')
            write(iuo+66,*) 'om : Version without NetCDF '//
     &       '=> writing of  res'//fchnum(ncfch:)//'.ncdf impossible !'
            call savrunb(koutpu,nn99,'res'//fchnum(ncfch:)//'.om')
          else
            call savrunb(koutpu,nn99,'res'//fchnum(ncfch:)//'.om')
          endif
          write(iuo+66,*) 'deriv:',deriv
          write(iuo+66,*) 'day:',xjour,' ------ ','year:',iyear
c--buffeur associed with evolu emptied (No=90) :
          if (numit+ninfo.le.nlast) call flush(90)
        endif
 
c-receive  and/or  send forcing in coupled mode
        if (icoupl.ge.1) then
Ccp1       call ocn_s2_cpl
Ccp2       call s2cpl2(numcpl)
c           call s2cpl3(numcpl)
c           write(iuo+66,*) 'sst and ice has been sent to cpl', numit
Ccp1       call ocn_rfm_cpl
Ccp1       write(iuo+66,*) 'fluxes have been received by sioclm', numit
        endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3) End of the run.                                                  |
c-----------------------------------------------------------------------
 
c      xday = xjour - 1. + dtsd2/daysec
c      write(iuo+66,'(A,I8,A,F7.2,A,I10)') '  end of the run ; year=', iyear,
c     &    ' +', xday, ' days ; iteration=', numit
c      if(nn99.eq.2) close(99)
 
c-deconnection of sioclm from PVM in coupled mode
      if (icoupl.ge.1) then
Ccp1     call ocn_stop_ocn
Ccp2     call quitcmo2
c         call quitcmo3
      endif
c     if (flgicb) call iceb_out
 
c      stop
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- end of the routine clio -
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine clio_inioutparat
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comclio.h'
      include 'comunit.h'

      character*60 part1
      integer      i,l

C** Here following are default values for above mentioned parameters. 
C** These parameters can be updated in the namelist.

      flagmon=.FALSE.
      flagyear=.FALSE.
      numtotvaro=0
      newtotvaro(:,:)=0.
      nametotvaro(:,:)=""
      dimtotvaro(:,:)=""

      do i=1,26
        globalatto(i,1)=globalatt(i,1)
        globalatto(i,2)=globalatt(i,2)
      enddo

      open(unit=30,file='outp_ocean.param',status='old',err=128)
      write(iuo+66,*) 'MESSAGE : file outp_ocean.param found,',
     &    '  checking for monthly/annual mean dumps to netcdf file:'
c
c -- switches for monthly/annual mean dumps
c
      read(30,'(/,/,/,/,/,/,/,/,/)')
      read(30,*) numtotvaro,maxmrecs,maxarecs,missing_valo
      read(30,*)

      do i=1,numtotvaro
        read(30,"(A)") part1
        nametotvaro(i,1)=trim(part1)
        read(30,*) (dimtotvaro(i,l),l=1,4)
        read(30,*) (nametotvaro(i,l),l=2,4)
        read(30,*) (newtotvaro(i,l),l=1,5)
        if(newtotvaro(i,2).eq.1) flagmon=.true.
        if(newtotvaro(i,3).eq.1) flagyear=.true.
      enddo

      if (flagmon) then
        write(iuo+66,*) '   monthly mean netcdf dumps will be written'//
     &               ' using the native model grid.'
      endif
      if (flagyear) then
        write(iuo+66,*) '   annual mean netcdf dumps will be written'//
     &               ' using the native model grid.'
      endif
      goto 129
  128 write(iuo+66,'(a,/)') 'WARNING : file outp_ocean.param not found',
     &    '  outave will perform dumps to cresum etc. as usual.'
  129 close(30)

      return
      end
