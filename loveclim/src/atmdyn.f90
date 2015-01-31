!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iatmdyn
!-----------------------------------------------------------------------
! *** initialise parameters and operators and read initial state
!-----------------------------------------------------------------------

      USE NETCDF    ! netcdf module to read forcing data


      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comcoup.h'
      include 'comunit.h'
      INCLUDE 'comrunlabel.h'  ! usefull to detect the right position in berg.nc

      integer i,j,k1,k2,k,l,m,n,ifail,ii,jj,i1,j1,nn
      real*8  pigr4,dis,dif,rll,ininag(nlat,nlon),asum
      real*8  r1,a,b,c,d,e,sqn,rsqn
      real*8  rnorm,rh0,dd
      real*8  agg(nlat,nlon), agg1(nlat,nlon), agg2(nlat,nlon)
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zagg1, zagg2
      real*8  areafac
      real*8  spv
      character*6 numyear
      character*3 numday

      ! netcdf file
      !LGM variable
      INTEGER :: idd_time, idf_berg, idv_time, idv_h, idv_sfric, istatus
      INTEGER :: ntime_berg, itime_berg
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: nvtime_berg

      read (iuo+1) nshm, ll

! *** legendre associated functions and derivatives

      read (iuo+1) pp
      read (iuo+1) pd
      read (iuo+1) pw

! *** compensation for normalization in nag fft routines

      sqn=sqrt(dble(nlon))
      rsqn=1d0/sqn
      do k=1,nsh
        do i=1,nlat
          pp(i,k)=pp(i,k)*sqn
          pd(i,k)=pd(i,k)*sqn
          pw(i,k)=pw(i,k)*rsqn
        enddo
      enddo

!       write(numyear,'(i6.6)') irunlabel+int((irunlabeld)/360)
!       write(numday,'(i3.3)') mod(irunlabeld,360)
      write(numyear,'(i6.6)') irunlabel
      write(numday,'(i3.3)') irunlabeld

! *** real parameter
! *** zonal derivative operator
! *** Laplace/Helmoltz
! *** Fft coefficients
! *** Orography and dissipation term
! *** Heigh orography in m
! *** Surface dependent friction
      CALL topo

! *** input initial qprime and for

      if (irunlabel .eq. 0) then
        do k=1,nsh2
          do l=1,nvl
            qprime(k,l)=0d0
            psi(k,l)=0d0
            for(k,l)=0d0
          enddo
          do l=1,ntl
            psit(k,l)=0d0
          enddo
        enddo
        do j=1,nlon
          do i=1,nlat
            u800(i,j)=0d0
            u500(i,j)=0d0
            u200(i,j)=0d0
            uv10(i,j)=0d0
            uvw10(i,j)=0d0
            do l=1,nvl
              utot(i,j,l)=0d0
              udivg(i,j,l)=0d0
              divg(i,j,l)=0d0
            enddo
          enddo
        enddo
      else

      open(iuo+95,file='startdata/inatdyn'//numyear//'_'//numday//'.dat', &
           & form='unformatted')
        read(iuo+95) qprime,for
        close(iuo+95)
      endif

  910 format(10e12.5)

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ddt
!----------------------------------------------------------------------
! *** computation of time derivative of the potential vorticity fields
! *** input qprime, psi, psit
! *** output dqprdt
! *** NOTE psit is destroyed
!----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer k,l
      real*8  dum1,dum2

! *** advection of potential vorticity at upper level

      call jacob (psi(1,1),qprime(1,1),dqprdt(1,1))

! *** advection of potential vorticity at middle level

      call jacob (psi(1,2),qprime(1,2),dqprdt(1,2))

! *** advection of potential vorticity and dissipation at lower level

      call jacobd (psi(1,3),qprime(1,3),dqprdt(1,3))

! *** relaxation of temperature and forcing

      do k=1,nsh2
        dum1=relt1*psit(k,1)
        dum2=relt2*psit(k,2)
        dqprdt(k,1)=dqprdt(k,1)+dum1              +for(k,1)
        dqprdt(k,2)=dqprdt(k,2)-dum1+dum2         +for(k,2)
        dqprdt(k,3)=dqprdt(k,3)          -dum2    +for(k,3)
      enddo

! *** explicit horizontal diffusion

      do l=1,3
        do k=1,nsh2
          dqprdt(k,l)=dqprdt(k,l)+diss(k,1)*qprime(k,l)
        enddo
      enddo
      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine jacob (psiloc,pvor,sjacob)
!----------------------------------------------------------------------
! *** advection of potential vorticity
! *** input psiloc, pvor
! *** output sjacob
!----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon), &
           & dvordm(nlat,nlon), gjacob(nlat,nlon), dpsidls(nsh2)

! *** space derivatives of potential vorticity

      call ddl (pvor,vv)
      call sptogg (vv,dvordl,pp)
      call sptogg (pvor,dvordm,pd)

! *** space derivatives of streamfunction

      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dvordl(i,j)-dpsidl(i,j)*dvordm(i,j)
        enddo
      enddo

      call ggtosp (gjacob,sjacob)

! *** planetary vorticity advection

      do k=1,nsh2
        sjacob(k)=sjacob(k)-dpsidls(k)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine jacobd (psiloc,pvor,sjacob)
!----------------------------------------------------------------------
! *** advection of potential vorticity and dissipation on gaussian grid
! *** input psiloc, pvor
! *** output sjacob
!----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon), &
           & dvordm(nlat,nlon), gjacob(nlat,nlon), vv(nsh2), &
           & azeta(nlat,nlon),dpsidls(nsh2)

! *** space derivatives of potential vorticity

      call ddl (pvor,vv)
      call sptogg (vv,dvordl,pp)
      call sptogg (pvor,dvordm,pd)

! *** space derivatives of streamfunction

      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)

! *** jacobian term + orographic forcing

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*(dvordl(i,j)+sinfi(i)*dorodl(i,j))- &
               & dpsidl(i,j)*(dvordm(i,j)+sinfi(i)*dorodm(i,j))
        enddo
      enddo

! *** dissipation


      if (lgdiss) then

! ***   spatially varying dissipation

        do k=1,nsh2
          vv(k)=diss(k,2)*psiloc(k)
        enddo

        call sptogg (vv,azeta,pp)

        do j=1,nlon
          do i=1,nlat
            gekdis(i,j)=-dpsidm(i,j)*ddisdy(i,j)-dpsidl(i,j)*ddisdx(i,j)+ &
                 & rdiss(i,j)*azeta(i,j)
            gjacob(i,j)=gjacob(i,j) + gekdis(i,j)
          enddo
        enddo

        call ggtosp (gjacob,sjacob)

      else

! ***   uniform dissipation

        call ggtosp (gjacob,sjacob)

        do k=1,nsh2
          sjacob(k)=sjacob(k)+diss(k,2)*psi(k,3)
        enddo

      endif

! *** planetary vorticity advection

      do k=1,nsh2
        sjacob(k)=sjacob(k)-dpsidls(k)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine jacobr (psiloc,pvor,sjacob)
!-----------------------------------------------------------------------
! *** computation of jacobian without planetary vorticity
! *** input psiloc, pvor
! *** output sjacob
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon), &
           & dvordm(nlat,nlon), gjacob(nlat,nlon), dpsidls(nsh2)

! *** space derivatives of potential vorticity

      call ddl (pvor,vv)
      call sptogg (vv,dvordl,pp)
      call sptogg (pvor,dvordm,pd)

! *** space derivatives of streamfunction

      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dvordl(i,j)-dpsidl(i,j)*dvordm(i,j)
        enddo
      enddo

      call ggtosp (gjacob,sjacob)


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine omoro (psiloc,sjacob)
!-----------------------------------------------------------------------
! *** computation of jacobian without planetary vorticity
! *** input psiloc, pvor
! *** output sjacob
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k
      real*8  psiloc(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), &
           & gjacob(nlat,nlon), dpsidls(nsh2)

! *** space derivatives of streamfunction

      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dorodl(i,j)-dpsidl(i,j)*dorodm(i,j)
          gjacob(i,j)=sinfi(i)*gjacob(i,j)
        enddo
      enddo

      call ggtosp (gjacob,sjacob)


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ddl (as,dadl)
!-----------------------------------------------------------------------
! *** zonal derivative in spectral space
! *** input spectral field as
! *** output spectral field dadl which is as differentiated wrt lambda
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer k
      real*8 as(nsh,2), dadl(nsh,2)

      do k=1,nsh
        dadl(k,1)=-rm(k)*as(k,2)
        dadl(k,2)= rm(k)*as(k,1)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sptogg (as,agg,pploc)

!-----------------------------------------------------------------------
! *** conversion from spectral coefficients to gaussian grid
! *** input  spectral field as, legendre polynomials pploc (pp or pd)
! ***        where pp are legendre polynomials and pd derivatives with
! ***        respect to sin(fi)
! *** output gaussian grid agg
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,ifail,j,k,k1,k2,m,mi,mr,nlon1
      real*8  as(nsh,2), agg(nlat,nlon), pploc(nlat,nsh)

! *** inverse legendre transform

      do j=1,nlon
        do i=1,nlat
          agg(i,j)=0.0d0
        enddo
      enddo

      nlon1=nlon+1
      k2=nshm(0)

      do k=1,k2
        do i=1,nlat
          agg(i,1)=agg(i,1)+as(k,1)*pploc(i,k)
        enddo
      enddo

      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            agg(i,mr)=agg(i,mr)+as(k,1)*pploc(i,k)
          enddo
          do i=1,nlat
            agg(i,mi)=agg(i,mi)-as(k,2)*pploc(i,k)
          enddo
        enddo
      enddo

! *** inverse fourier transform

      ifail=0
      call c06fqf (nlat,nlon,agg,'r',trigi,wgg,ifail)

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ggtosp (agg,as)
!-----------------------------------------------------------------------
! *** conversion from gaussian grid (agg) to spectral coefficients (as)
! *** input array agg is destroyed
! *** output as contains spectral coefficients
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer ir,ifail,j,k,k1,k2,m,mi,mr,nlon1,i
      real*8  as(nsh,2), agg(nlat,nlon)
!
! *** fourier transform
!
      ifail=0
      call c06fpf (nlat,nlon,agg,'r',trigd,wgg,ifail)
!
! *** legendre transform
!
      do ir=1,2
        do k=1,nsh
          as(k,ir)=0.0d0
        enddo
      enddo

      nlon1=nlon+1

      k2=nshm(0)

      do k=1,k2
        do i=1,nlat
          as(k,1)=as(k,1)+agg(i,1)*pw(i,k)
        enddo
      enddo

      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            as(k,1)=as(k,1)+agg(i,mr)*pw(i,k)
            as(k,2)=as(k,2)+agg(i,mi)*pw(i,k)
          enddo
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine rggtosp (agg,as)
!-----------------------------------------------------------------------
! *** conversion from gaussian grid (agg) to spectral coefficients (as)
! *** input array agg is conserved
! *** output as contains spectral coefficients
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,ifail,ir,j,k,k1,k2,m,mi,mr,nlon1
      real*8 as(nsh,2), agg(nlat,nlon)
      real*8 store(nlat,nlon)

      do j=1,nlon
        do i=1,nlat
          store(i,j)=agg(i,j)
        enddo
      enddo

! *** fourier transform

      ifail=0
      call c06fpf (nlat,nlon,store,'r',trigd,wgg,ifail)

! *** legendre transform

      do ir=1,2
        do k=1,nsh
          as(k,ir)=0.0d0
        enddo
      enddo

      nlon1=nlon+1

      k2=nshm(0)

      do k=1,k2
        do i=1,nlat
          as(k,1)=as(k,1)+store(i,1)*pw(i,k)
        enddo
      enddo

      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            as(k,1)=as(k,1)+store(i,mr)*pw(i,k)
            as(k,2)=as(k,2)+store(i,mi)*pw(i,k)
          enddo
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine qtopsi
!-----------------------------------------------------------------------
! *** computation of streamfunction from potential vorticity
! *** input  qprime which is potential vorticity field
! *** output psi, the streamfunction and psit, the layer thicknesses
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer k
      real*8  r3

      do k=1,nsh2
        ws(k)=qprime(k,1)+qprime(k,3)
        psi(k,1)=rinhel(k,1)*(ws(k)+qprime(k,2))
        psi(k,2)=ws(k)-2.d0*qprime(k,2)
        psi(k,3)=qprime(k,1)-qprime(k,3)
      enddo

      do k=1,nsh2
        psit(k,1)=rinhel(k,2)*psi(k,2)+rinhel(k,3)*psi(k,3)
        psit(k,2)=rinhel(k,4)*psi(k,2)+rinhel(k,5)*psi(k,3)
      enddo

      r3=1./3
      do k=1,nsh2
        psi(k,2)=r3*(psi(k,1)-psit(k,1)+psit(k,2))
        psi(k,1)=psi(k,2)+psit(k,1)
        psi(k,3)=psi(k,2)-psit(k,2)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine psitoq
!-----------------------------------------------------------------------
! *** computation of potential vorticity from stream function
! *** input psi streamfunction
! *** output qprime, the potential vorticity and psit, the layer thick.
!-----------------------------------------------------------------------

      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      integer k

      do k=1,nsh2
        psit(k,1)=psi(k,1)-psi(k,2)
        psit(k,2)=psi(k,2)-psi(k,3)
        qprime(k,1)=rinhel(k,0)*psi(k,1)-rl1*psit(k,1)
        qprime(k,2)=rinhel(k,0)*psi(k,2)+rl1*psit(k,1)-rl2*psit(k,2)
        qprime(k,3)=rinhel(k,0)*psi(k,3)+rl2*psit(k,2)
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine psiq(sfin,qout)
!-----------------------------------------------------------------------
! ***  computation of potential vorticity qout from stream function sfin
!-----------------------------------------------------------------------

      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      integer k
      real*8  sfin(nsh2,nvl),qout(nsh2,nvl),tus(nsh2)

      do k=1,nsh2
        tus(k)=rl1*sfin(k,1)-rl1*sfin(k,2)
      enddo

      do k=1,nsh2
        qout(k,1)=rinhel(k,0)*sfin(k,1)-tus(k)
        qout(k,2)=rinhel(k,0)*sfin(k,2)+tus(k)
      enddo

      do k=1,nsh2
        tus(k)=rl2*sfin(k,2)-rl2*sfin(k,3)
      enddo

      do k=1,nsh2
        qout(k,2)=qout(k,2)-tus(k)
        qout(k,3)=rinhel(k,0)*sfin(k,3)+tus(k)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine qpsi(qin,sfout)
!-----------------------------------------------------------------------
! *** computation of streamfunction bb from potential vorticity qin
!-----------------------------------------------------------------------

      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      real*8  qin(nsh2,nvl),sfout(nsh2,nvl), tus(nsh2,ntl), r3
      integer k

      do k=1,nsh2
        ws(k)=qin(k,1)+qin(k,3)
        sfout(k,1)=rinhel(k,1)*(ws(k)+qin(k,2))
        sfout(k,2)=ws(k)-2.*qin(k,2)
        sfout(k,3)=qin(k,1)-qin(k,3)
      enddo

      do k=1,nsh2
        tus(k,1)=rinhel(k,2)*sfout(k,2)+rinhel(k,3)*sfout(k,3)
        tus(k,2)=rinhel(k,4)*sfout(k,2)+rinhel(k,5)*sfout(k,3)
      enddo

      r3=1./3
      do k=1,nsh2
        sfout(k,2)=r3*(sfout(k,1)-tus(k,1)+tus(k,2))
        sfout(k,1)=sfout(k,2)+tus(k,1)
        sfout(k,3)=sfout(k,2)-tus(k,2)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine qpsit(qin,tus)
!-----------------------------------------------------------------------
! *** computation of streamfunction bb from potential vorticity qin
!-----------------------------------------------------------------------

      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      real*8  qin(nsh2,nvl),tus(nsh2,ntl), r3,sfout(nsh2,nvl)
      integer k

      do k=1,nsh2
        ws(k)=qin(k,1)+qin(k,3)
        sfout(k,1)=rinhel(k,1)*(ws(k)+qin(k,2))
        sfout(k,2)=ws(k)-2.*qin(k,2)
        sfout(k,3)=qin(k,1)-qin(k,3)
      enddo

      do k=1,nsh2
        tus(k,1)=rinhel(k,2)*sfout(k,2)+rinhel(k,3)*sfout(k,3)
        tus(k,2)=rinhel(k,4)*sfout(k,2)+rinhel(k,5)*sfout(k,3)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fmtofs (y,z)
!-----------------------------------------------------------------------
! *** transforms franco's format to the french format for global fields
! *** input  y spectral coefficients in franco's format
! *** output z spectral coefficients in french format
!-----------------------------------------------------------------------

      implicit none
      include'comatm.h'

      integer   m,n,k,indx,l
      real*8    y(nsh2,nvl),z(nsh2,nvl)

      do l=1,nvl
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if (m.eq.0) then
              indx=n**2
            else
              indx=n**2+2*m-1
            end if
            z(indx,l)=y(k,l)
            if (m.ne.0) z(indx+1,l)=y(k+nsh,l)
          enddo
        enddo
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fstofm (y,z,ntr)
!-----------------------------------------------------------------------
! *** transforms the french format to franco's format for global fields
! *** input  y spectral coef. in french format, ntr is truncation limit
! *** output z spectral coefficients in franco's format
!-----------------------------------------------------------------------

      implicit none
      include 'comatm.h'

      integer   m,n,k,indx,i,l,ntr
      real*8    y(nsh2,nvl),z(nsh2,nvl)

      do l=1,nvl
        do i=1,nsh2
          z(i,l)=0d0
        enddo
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if ((m.le.ntr).and.(n.le.ntr)) then
              if (m.eq.0) then
                indx=n**2
              else
                indx=n**2+2*m-1
              end if
              z(k,l)=y(indx,l)
              if (m.ne.0) z(k+nsh,l)=y(indx+1,l)
            endif
          enddo
        enddo
      enddo
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine forward
!-----------------------------------------------------------------------
! *** performs a fourth order runge kutta time step at truncation nm
! *** with time step dt
! *** dqdt calculates the time derivative
! *** input  qprime at current time
! *** output qprime at current time plus dt
!-----------------------------------------------------------------------
      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      integer  k,l,nvar
      real*8   dt2,dt6
      real*8   y(nsh2,nvl),dydt(nsh2,nvl),yt(nsh2,nvl)
      real*8   dyt(nsh2,nvl),dym(nsh2,nvl)

      nvar=(nm+2)*nm
      dt2=dtt*0.5d0
      dt6=dtt/6d0
      call fmtofs(qprime,y)
      call dqdt(y,dydt)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dt2*dydt(k,l)
        enddo
      enddo
      call dqdt(yt,dyt)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dt2*dyt(k,l)
        enddo
      enddo
      call dqdt(yt,dym)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dtt*dym(k,l)
          dym(k,l)=dyt(k,l)+dym(k,l)
        enddo
      enddo
      call dqdt(yt,dyt)
      do l=1,nvl
        do k=1,nvar
          y(k,l)=y(k,l)+dt6*(dydt(k,l)+dyt(k,l)+2.*dym(k,l))
        enddo
      enddo
      call fstofm(y,qprime,nm)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine dqdt(y,dydt)
!-----------------------------------------------------------------------
! *** computation of time derivative of the potential vorticity field
! *** input  y potential vorticity in french format
! *** output dydt time derivative of y in french format
! *** values of qprime, psi and psit are changed
!-----------------------------------------------------------------------

      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      real*8  y(nsh2,nvl),dydt(nsh2,nvl)

      call fstofm(y,qprime,nm)
      call qtopsi
      call ddt
      call fmtofs(dqprdt,dydt)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine psitogeo
!-----------------------------------------------------------------------
! *** computes geopotential in [m2/s2] from the streamfunction
! *** by solving the linear balance equation:
! *** del phi = (1 - mu**2 ) d psi/dmu + mu del psi
! *** the global mean value is not determined and set to zero
! *** input:  psi
! *** output: grpsi1,grpsi2,grpsi3,geopg(nlat,nlon,nvl)
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer i,j,l,k
      real*8  dmu(nlat),cdim,sdim,tempfac(ntl)
      real*8  delpsis(nsh2),delpsig(nlat,nlon)
      real*8  dmupsig(nlat,nlon),delgeog(nlat,nlon)
      real*8  delgeos(nsh2),geos(nsh2)


      do i=1,nlat
        dmu(i)=1-sinfi(i)**2
      enddo

      cdim=(om**2)*(radius**2)

      do l=1,nvl

! *** solve linear balance equation

        call lap(psi(1,l),delpsis)
        call sptogg(delpsis,delpsig,pp)
        call sptogg(psi(1,l),dmupsig,pd)

        do j=1,nlon
          do i=1,nlat
            delgeog(i,j)=dmu(i)*dmupsig(i,j)+sinfi(i)*delpsig(i,j)
          enddo
        enddo
        call ggtosp(delgeog,delgeos)
        call lapinv(delgeos,geos)
        geos(1)=0.d0
        call sptogg(geos,geopg(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            geopg(i,j,l)=cdim*geopg(i,j,l)
          enddo
        enddo

      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine omega3
!-----------------------------------------------------------------------
! *** computes the vertical velocity at the two temperature levels
! *** and the surface in pa/sec
! *** input psi,psit
! *** output omegs is vertical velocity at three levels in spectral form
! ***        omegg at gaussian grid
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comemic.h'

      integer i,j,k,l
      real*8  adoro(nsh2),adpsit(nsh2,ntl),dpsitdt(nsh2,ntl)
      real*8  facom(nvl),facoc
      real*8  facd1,facd2
      real*8  facekm
      real*8  omegsd(nsh2)

      facom(1)=(dp*om)/(rrdef1**2*fzero)
      facom(2)=(dp*om)/(rrdef2**2*fzero)
      facom(3)=(dp*om)

      call jacobr(psi(1,1),psit(1,1),adpsit(1,1))
      call jacobr(psi(1,2),psit(1,2),adpsit(1,2))
      call omoro(psi(1,3),adoro)

      call ddt

      call qpsit(dqprdt,dpsitdt)

      facoc=trel*4.d0*pi


! *** ekman dissipation induced omega at ground level
! *** rein

      do j=1,nlon
        do i=1,nlat
!          if (i.gt.14.and.i.lt.19) then
!            gekdis(i,j)=gekdis(i,j)/abs(sinfi(14))
!          else
!            gekdis(i,j)=gekdis(i,j)/abs(sinfi(i))
!          endif
          gekdis(i,j)=gekdis(i,j)/fzero
        enddo
      enddo

      call ggtosp(gekdis,omegsd)

! *** first contributions that are odd in the equator

      do k=1,nsh2

! ***   level1  (350 hpa)

        omegs(k,1)= dpsitdt(k,1) - adpsit(k,1) + psit(k,1)/facoc

! ***   level2  (650 hpa)

        omegs(k,2)= dpsitdt(k,2) - adpsit(k,2) + psit(k,2)/facoc

! ***   level 3 (surface)

!        omegs(k,3)= omegsd(k) + adoro(k)
        omegs(k,3)= 0d0

      enddo

      do l=1,nvl
        do k=1,nsh2
          omegs(k,l)=omegs(k,l)*facom(l)
        enddo
      enddo

      do l=1,nvl
        call sptogg(omegs(1,l),omegg(1,1,l),pp)
      enddo

      do l=1,nvl
        do j=1,nlon
          do i=1,nlat/2
            omegg(i,j,l)=-omegg(i,j,l)
          enddo
        enddo
      enddo

      do l=1,nvl
        call ggtosp(omegg(1,1,l),omegs(1,l))
      enddo


! *** second contributions that are even in the equator
! *** diabatic forcing and lift due to orography

      facd1=(rgas*(dp**2))/((rrdef1*radius*om*fzero)**2*tlevel(1))
      facd2=(rgas*(dp**2))/((rrdef2*radius*om*fzero)**2*tlevel(2))

      do k=1,nsh2
        omegs(k,1)=omegs(k,1) - dfor1(k)*facd1
        omegs(k,2)=omegs(k,2) - dfor2(k)*facd2
!        omegs(k,3)=omegs(k,3) + adoro(k)*facom(3)
      enddo

      do l=1,nvl
        omegs(1,l)=0.
        call sptogg(omegs(1,l),omegg(1,1,l),pp)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine diver
!-----------------------------------------------------------------------
! *** computes divergence from omega using conservation of mass
! *** input omegs is vertical velocity in pa/s
! *** output divs is divergence in spectral form
! ***        divg at gaussian grid
!-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer k,l

      do k=1,nsh2
        divs(k,1)=-omegs(k,1)/dp
        divs(k,2)=(omegs(k,1)-omegs(k,2))/dp
        divs(k,3)=(omegs(k,2)-omegs(k,3))/dp
      enddo

      do l=1,nvl
        call sptogg(divs(1,l),divg(1,1,l),pp)
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine divwin
!-----------------------------------------------------------------------
! *** computes divergent wind from the divergence
! *** input  divs is divergence at the three pressure levels
! *** output udivg zonal divergent wind at gaussian grid
! ***        vdivg meridional divergent wind at gaussian grid
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k,l
      real*8  r1,r2
      real*8  x(nsh2),xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)

      r1=radius
      r2=radius**2

      do l=1,nvl

         do k=1,nsh2
           x(k)=r2*divs(k,l)*rinhel(k,1)
           chi(k,l)=x(k)
	 enddo

         call sptogg(x,chig(1,1,l),pp)

         call ddl(x,xhelp)
         call sptogg(xhelp,dxdl,pp)
         call sptogg(x,dxdm,pd)

         do j=1,nlon
           do i=1,nlat
             udivg(i,j,l)=dxdl(i,j)/(r1*cosfi(i))
             vdivg(i,j,l)=dxdm(i,j)*cosfi(i)/r1
           enddo
         enddo

      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine geowin
!-----------------------------------------------------------------------
! *** computation of geostrophic winds at all levels
! *** input psi streamfunction in spectral form
! *** output u200,v200 wind components at 200 hPa
! ***        u500,v500 wind components at 500 hPa
! ***        u800,v800 wind components at 800 hPa
!
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k,l
      real*8  facwin2
      real*8  dpsdl(nlat,nlon),dpsdm(nlat,nlon),psik(nsh2),vv(nsh2)

! *** space derivatives of streamfunction

      facwin2=radius*om

      do l=1,nvl

        do k=1,nsh2
          psik(k)=psi(k,l)
        enddo

        call ddl (psik,vv)
        call sptogg (vv,dpsdl,pp)
        call sptogg (psik,dpsdm,pd)

        if (l.eq.1) then
          do j=1,nlon
            do i=1,nlat
              u200(i,j)=-facwin2*dpsdm(i,j)*cosfi(i)
              v200(i,j)=+facwin2*dpsdl(i,j)/cosfi(i)
            enddo
          enddo
        endif

        if (l.eq.2) then
          do j=1,nlon
            do i=1,nlat
              u500(i,j)=-facwin2*dpsdm(i,j)*cosfi(i)
              v500(i,j)=+facwin2*dpsdl(i,j)/cosfi(i)
            enddo
          enddo
        endif

        if (l.eq.3) then
          do j=1,nlon
            do i=1,nlat
              u800(i,j)=-facwin2*dpsdm(i,j)*cosfi(i)
              v800(i,j)=+facwin2*dpsdl(i,j)/cosfi(i)
            enddo
          enddo
        endif

      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine totwind
!-----------------------------------------------------------------------
! *** computation of total wind at all levels
!
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer i,j,k,l

      do l=1,nvl

        if (l.eq.1) then
          do j=1,nlon
            do i=1,nlat
              utot(i,j,l)=u200(i,j) + udivg(i,j,l)
              vtot(i,j,l)=v200(i,j) + vdivg(i,j,l)
            enddo
          enddo
        endif

        if (l.eq.2) then
          do j=1,nlon
            do i=1,nlat
              utot(i,j,l)=u500(i,j) + udivg(i,j,l)
              vtot(i,j,l)=v500(i,j) + vdivg(i,j,l)
            enddo
          enddo
        endif

        if (l.eq.3) then
          do j=1,nlon
            do i=1,nlat
              utot(i,j,l)=u800(i,j) + udivg(i,j,l)
              vtot(i,j,l)=v800(i,j) + vdivg(i,j,l)
            enddo
          enddo
        endif

      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lap(xs,xsl)
!-----------------------------------------------------------------------
! *** computation of laplace operator in spectral domain
! *** input  xs  field in spectral form
! *** output xsl laplace of xs in spectral form
!-----------------------------------------------------------------------
      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      integer k
      real*8  xs(nsh2),xsl(nsh2)

      do k=1,nsh2
        xsl(k)=xs(k)*rinhel(k,0)
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lapinv(xsl,xs)
!-----------------------------------------------------------------------
! *** computation of laplace operator in spectral domain
! *** input  xsl field in spectral form
! *** output xs  inverse laplace of xs in spectral form
!-----------------------------------------------------------------------
      implicit none
      include 'comatm.h'
      include 'comdyn.h'

      integer k
      real*8  xs(nsh2),xsl(nsh2)

      do k=1,nsh2
        xs(k)=xsl(k)*rinhel(k,1)
      enddo

      return
      end

!123456789012345678901234567890123456789012345678901234567890123456789012
      subroutine topo
!-----------------------------------------------------------------------
! *** Purpose : replace topography by 0 topography
! *** Original code : Cedric Van Meerbeek
! *** Modified code : P.Mathiot (01/2012) Allow interannual modification
!-----------------------------------------------------------------------
      USE NETCDF

      implicit none

      include'comatm.h'
      include'comsurf.h'
      include'comphys.h'
      include'comdyn.h'
      include'comemic.h'
      include'comcoup.h'
      include'comunit.h'
      INCLUDE 'comrunlabel.h'

      real*8 asum,spv
      integer i,j,i1,j1,ii,jj

      integer k1,k2,k,l,m,n,ifail,nn
      real*8  pigr4,dis,dif,rll,ininag(nlat,nlon)
      real*8  r1,a,b,c,d,e,sqn,rsqn
      real*8  rnorm,rh0,dd
      real*8  agg(nlat,nlon), agg1(nlat,nlon), agg2(nlat,nlon)
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zagg1, zagg2
      real*8  fmu(nlat,2)
      real*8  wsx(nsh2),areafac

      INTEGER :: idd_time, idf_berg, idv_time, idv_h, idv_sfric, istatus
      INTEGER :: ntime_berg, itime_berg
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: nvtime_berg

!-pour changer topo (dyn+thermo), mettre le nouveau champ dans
!-rmount_ism(verifier les flgism):

! *** real parameters
      pigr4=4.d0*pi
      rl1=1.0d0/rrdef1**2
      rl2=1.0d0/rrdef2**2
      relt1=max(0.0d0,rl1/(trel*pigr4))
      relt2=max(0.0d0,rl2/(trel*pigr4))
      dis=max(0.0d0,1.0d0/(tdis*pigr4))
      rll=dble(ll(nsh))
      dif=max(0.0d0,1.0d0/(tdif*pigr4*(rll*(rll+1))**idif))

! *** zonal derivative operator
      k2=0
      do m=0,nm
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          rm(k)=dble(m)
        enddo
      enddo

! *** laplace/helmholtz direct and inverse operators
      do j=0,5
        rinhel(1,j)=0.0d0
      enddo

      diss(1,1)=0.0d0
      diss(1,2)=0.0d0

      do k=2,nsh
        r1=dble(ll(k)*(ll(k)+1))
        a=-r1-3.0d0*rl1
        b=-r1-3.0d0*rl2
        c=-r1-rl1
        d=-r1-rl2
        e=a*d+b*c
        rinhel(k,0)=-r1
        rinhel(k,1)=-1.0d0/r1
        rinhel(k,2)= d/e
        rinhel(k,3)= b/e
        rinhel(k,4)=-c/e
        rinhel(k,5)= a/e
        diss(k,2)=dis*r1
        diss(k,1)=-dif*r1**idif
      enddo

      do j=0,5
        do k=1,nsh
          rinhel(k+nsh,j)=rinhel(k,j)
        enddo
      enddo

      do j=1,2
        do k=1,nsh
          diss(k+nsh,j)=diss(k,j)
        enddo
      enddo

! *** computation of the weight needed in OASIS
!     include 'weight.h'

! *** initialization of coefficients for fft
      do j=1,nlon
        do i=1,nlat
          ininag(i,j)=1.0d0
        enddo
      enddo

      ifail=0
      call c06fpf (nlat,nlon,ininag,'i',trigd,wgg,ifail)

      ifail=0
      call c06fqf (nlat,nlon,ininag,'i',trigi,wgg,ifail)

! *** orography and dissipation terms
! *** fmu(i,1): sin(phi(i))
! *** fmu(i,2): 1-sin**2(phi(i))
      rnorm=1.0d0/sqrt(3.0d0*nlon)
      do i=1,nlat
        fmu(i,1)=rnorm*pp(i,2)
        fmu(i,2)=1.d0-fmu(i,1)**2
      enddo

! *** height of orography in meters

      ! LOAD netcdf forcing file
      istatus=NF90_OPEN("inputdata/berg.nc", NF90_NOWRITE, idf_berg)
      istatus=NF90_INQ_DIMID(idf_berg, 'time', idd_time)
      istatus=NF90_INQUIRE_DIMENSION(idf_berg, idd_time, len = ntime_berg)
      istatus=nf90_inq_varid(idf_berg, 'time', idv_time)
      istatus=nf90_inq_varid(idf_berg, "h", idv_h)
      istatus=nf90_inq_varid(idf_berg, "sfric", idv_sfric)

      ! ALLOCATE variable
      ALLOCATE(nvtime_berg(ntime_berg), zagg1(nlon, nlat))

      ! READ time variable
      istatus=nf90_get_var(idf_berg, idv_time, nvtime_berg)

      ! SELECT the right time
      itime_berg=0
      DO i=1,ntime_berg
         IF (INT(ABS((nvtime_berg(i)-(irunlabelf+iyear)))) == &
              & INT(MINVAL(ABS(nvtime_berg(:)-(irunlabelf+iyear))))) &
              & itime_berg = i
      END DO
      IF (itime_berg==0) THEN
         PRINT *, "topo = error in detection time slide : ", itime_berg
         STOP
      END IF

      ! READ variable agg1
      istatus=nf90_get_var(idf_berg, idv_h, zagg1, &
           & start = (/1,1,itime_berg/), count = (/nlon,nlat,1/))
      agg1=TRANSPOSE(zagg1)

      ! UPDATE rmount
      rh0=max(0.0d0,0.001d0/h0)
        do j=1,nlon
          do i=1,nlat
!          agg(i,j)=fmu(i,1)*agg1(i,j)*rh0
          agg(i,j) = agg1(i,j)*rh0
          rmount(i,j)=agg1(i,j)
          if (rmount(i,j).lt.0d0) rmount(i,j)=0d0
         enddo
       enddo

      WRITE(iuo+99,*) "===================== topo ==========================="
      WRITE(iuo+99,*) "Read berg in file ",TRIM('inputdata/berg.nc')
      WRITE(iuo+99,*) "For year (AD) = ",irunlabelf+iyear
      WRITE(iuo+99,*) "selected time is ", nvtime_berg(itime_berg), &
           & " index is ",itime_berg,"/",ntime_berg
      WRITE(iuo+99,*) "======================================================"
      CALL FLUSH(iuo+99)

      DEALLOCATE(zagg1, nvtime_berg)

! *** surface dependent friction
      lgdiss=((addisl.gt.0.0).or.(addish.gt.0.0))

      call ggtosp (agg,orog)
      call ddl (orog,ws)
      call sptogg (ws,dorodl,pp)
      call sptogg (orog,dorodm,pd)

      if (lgdiss) then
        ! READ friction variable
        ALLOCATE(zagg2(nlon, nlat))
        istatus=nf90_get_var(idf_berg, idv_sfric, zagg2, &
             & start = (/1,1,itime_berg/), count = (/nlon,nlat,1/))
        agg2=TRANSPOSE(zagg2)
        DEALLOCATE(zagg2)

        do j=1,nlon
          do i=1,nlat
            agg(i,j)=1.0d0+addisl*agg2(i,j)+ &
     &                addish*(1.0d0-exp(-0.001d0*agg1(i,j)))
          enddo
        enddo

        call ggtosp (agg,ws)
        call ddl (ws,wsx)

        call sptogg (ws,rdiss,pp)
        call sptogg (wsx,ddisdx,pp)
        call sptogg (ws,ddisdy,pd)

        dd=0.5d0*diss(2,2)
        do j=1,nlon
          do i=1,nlat
            ddisdx(i,j)=dd*ddisdx(i,j)/fmu(i,2)
            ddisdy(i,j)=dd*ddisdy(i,j)*fmu(i,2)
          enddo
        enddo

      endif

       ! PRINT *, "topo rmount = ",SUM(rmount)
       ! PRINT *, "topo agg    = ",SUM(agg)
       ! Close the file
       istatus=nf90_close(idf_berg)
      end
