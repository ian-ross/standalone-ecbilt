c23456789012345678901234567890123456789012345678901234567890123456789012
      program ECbilt
c-----------------------------------------------------------------------
c *** coupled ocean-atmosphere-seaice model.
c *** 10 augustus 1995 KNMI, De Bilt
c ***
c *** joint project Rein Haarsma
c ***               Geert Lenderink
c ***               Theo Opsteegh
c ***               Liu Qing
c ***               Frank Selten
c ***               Xueli Wang
c ***               Nanne Weber
c-----------------------------------------------------------------------

      implicit none

      include 'comglobal.h'

      integer icount,istep,i,j

c *** open files

      include 'openatinpfiles.h'
      include 'openocfixfiles.h'

c *** initialisation of parameters and initial state

      call initecbilt


c *** forward time integration

      do icount=1,ntstep
        write (*,*) 'icount=', icount

        istep=icount
        call mdldate(istep)

        call oceanfixed

        if (irunatm.eq.1) then
          call landmois
          call atmmodel
          call landtemp
          call surfacetemp
          call atmstate
          call fluxmodel
          if (isatfor.eq.1) call outato(istep)
        else

        endif

        call checks(istep)

        if (irunatm.eq.1) call atmout(istep)

        call writestate(istep)

      enddo

      call writestate(ntstep)
      call error(999)

      end
