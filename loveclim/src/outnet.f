












      subroutine outnet(ja,xjour,xjour1)


      USE Ocean_Output

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  This routine computes the average of some variables and write it
c  on the ouput files.
c
c  This is a major revision of outave that implements netcdf-output of
c  monthly and/or annual means following as closely as possible the
c  "NetCDF Climate and Forecast (CF) Metadata Conventions" (refer to
c  http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html for more
c  information). The update retains full compatibility with the former
c  output to cresum.dat etc. Outside outave.f only minor changes have
c  been applied to defgrid.f, where true gridpoint longitudes and
c  latitudes of the combined WW- and AA-grids are computed, and to
c  bloc0.com for some additional arrays that needed to be defined.
c
c  This work has been funded by the Deutsche ForschungsGemeinschaft
c  (DFG) within the research framework "Passagen" (please visit
c  http://www.passagen.uni-kiel.de). Extensions of the code have
c  been enclosed by comment lines

cDFG start ...

cDFG end
c
c  Even though this revision of the code is a bit lengthy due to the
c  administrative nature of netcdf, it is fairly well commented and
c  should be easily understood and flexible enough for modifications
c  according to other users' needs.
c
c  Oct-2002, Peter Herrmann (ph@phsck.de)
c
c---
c Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
c---
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'thermo.com'
      include 'densit.com'
      include 'isoslope.com'
      include 'comclio.h'

c- nn99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
      integer  nn99
      real*8   zfluxmt,zfluxmts

      common / mchd99 / nn99
      common / zfluc / zfluxmt,zfluxmts

c--Output variables and arrays.
c
      integer njcum(0:12),zlim1(1:4),zlim2(1:4),noumax
      parameter (noumax=47)

c Choice noumax :18,34,43,47, noumax2 depending of the choice
      real*8 cmoyan(imax,jmax,noumax),cmoymo(imax,jmax,noumax),
     &       cmoymap(imax,jmax,noumax,12)
cDFG start real*8 define for lahey compatibility
cDFG original define is plain real
      real cmulti(noumax),cadd(noumax)
cDFG end
      integer ncor(noumax),nmap(12),nstreamout

c streamfunction part
c options: nstreamupdate=0(no output), 1(monthly) or 2(annually)
cDFG start
c In order to interfere as little as possible with streamfunc.f,
c we here let nstreamout unchanged to yield monthly means computed
c by that original routine. Those transports are exported from
c streamfunc.f through common streamflu. We will then add them up
c here for annual means... Please leave nstreamout untouched!
cDFG end
      parameter (nstreamout=1)

      real*8 umm(imax,jmax,kmax),vmm(imax,jmax,kmax),
     &       tmm(imax,jmax,kmax),smm(imax,jmax,kmax),
     &       uma(imax,jmax,kmax),vma(imax,jmax,kmax),
     &       tma(imax,jmax,kmax),sma(imax,jmax,kmax)

      common /totalanu/ umm,vmm,tmm,smm,uma,vma,
     &                  tma,sma

cDFG start additonal fields for monthly and annual means of w
      real*8  wma(imax,jmax,kmax+1),wmm(imax,jmax,kmax+1)
      common /totalanu/ wma,wmm
c
      integer klevm(4,2),klevv,njm,noumef,jmoisi,jmoifl

      common/moyout/ cmoyan,cmoymo,cmoymap
      common/parout/ klevm,klevv,njm
      common/parout2/ noumef
      common/parout3/ ncor,nmap,jmois,nptj
      common/parout4/ cmulti,cadd,jmoisi,jmoifl
      common/parout5/ znumas(imax,jmax,4),znumav(imax,jmax,5)

cDFG start
c
c declarations for netcdf output
c
c ncidm, mcida: integer IDs for netcdf files
c lcdf[mon,ann]: logical switch for [monthly,annual] mean netcdf dumps
c l[v,m,a]cdf: logical switches for variables to dump to netcdf file
c
c copied from streamfunc.f: jmtt, kmtt, uuu, vvv
c
      integer start(4),count(4),jmtt,kmtt
      parameter (jmtt=57, kmtt=kmax+1)
      integer iMtimrec,iAtimrec
      common /icdfout/iMtimrec,iAtimrec
c
      real*8 f2d(imax,jmax),etamo(imax,jmax),ubmo(imax,jmax),
     &       vbmo(imax,jmax),etaan(imax,jmax),uban(imax,jmax),
     &       vban(imax,jmax)

      real*4 uuu(jmtt,0:kmtt,0:nbsmax),vvv(jmtt,0:nbsmax+4)
      common /streamflu/ uuu,vvv
cDFG start define meridional arrays of minimum bathymetry
c    for masking of streamfunction in netcdf output
c
      integer kfloor(jmtt,0:nbsmax)
      common/streambath/ kfloor
cDFG end
      real*8 moc_m(0:nbsmax,jmtt,0:kmtt),mht_m(0:nbsmax,jmtt),
     &       mst_m(0:nbsmax,jmtt),
     &       moc_a(0:nbsmax,jmtt,0:kmtt),mht_a(0:nbsmax,jmtt),
     &       mst_a(0:nbsmax,jmtt)
cDFG end
c
      integer ijour,ijour1,jmois,ijouri,jjoufl
      integer i,j,k,kk
      real*8  xjour,xjour1,undim6

      data njcum/0,30,60,90,120,150,180,210,240,270,300,330,360/
      data zlim1/0,500,1500,3500/
      data zlim2/500,1500,3500,5600/

      ijour=int(xjour)
      ijour1=int(xjour1)
C     jours=xjour+dts(ks2)/86400
      undim6=1.D-6
C     write(iuo+66,*) 'debut outave.f',ja,jmois,ijour

      if(numit.eq.nstart) then

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1) INITIALIZATIONS.                                                 |
c-----------------------------------------------------------------------

        zfluxmt=0.0
        zfluxmts=0.0

c--1.1 OPEN FILES FOR INITIALIZATION AND FOR OUTPUTS.
c----------------------------------------------------
c
        open(42,file='correcw.dat',form='formatted')
c
c--1.3 Computation of the month.
c---------------------------------
c
        jmois=1
        do 10 i=1,ijour
           if (ijour.gt.njcum(jmois)) jmois=jmois+1
 10     continue
        jmoisi=jmois
        ijouri=ijour-njcum(jmois-1)-1
        jmoifl=0
        jjoufl=0
        njm=0
c
c--1.4. INITIALIZATION OF CUMULATIVE ARRAYS
c-------------------------------------------
c

c streamfunction part
        if (nstreamout.ne.0) call streamfunc(-1,0.d0,nstreamout)

        cmoyan(:,:,:)=0.0
        cmoymo(:,:,:)=0.0
        cmoymap(:,:,:,:)=0.0

        umm(:,:,:)= 0.0
        vmm(:,:,:)= 0.0
        tmm(:,:,:)= 0.0
        smm(:,:,:)= 0.0

        uma(:,:,:)= 0.0
        vma(:,:,:)= 0.0
        tma(:,:,:)= 0.0
        sma(:,:,:)= 0.0

c
cDFG start initialisation of arrays for netcdf-output
c
        f2d(:,:)=0.0
        etamo(:,:)=0.0
        ubmo(:,:)=0.0
        vbmo(:,:)=0.0
        etaan(:,:)=0.0
        uban(:,:)=0.0
        vban(:,:)=0.0

        wmm(:,:,:)=0.0
        wma(:,:,:)=0.0

        moc_m(:,:,:)=0.0
        moc_a(:,:,:)=0.0

        mht_m(:,:)=0.0
        mst_m(:,:)=0.0
        mht_a(:,:)=0.0
        mst_a(:,:)=0.0
c
c---------------------------------
c
        noumef = 47
c
c Computation of the levels for averages
        do k=ks2,ks1,-1
          if ((-60).ge.zw(k)) then
             klevv=k
             goto 141
          endif
        enddo
141     continue
        klevm(:,:)=0
        if (nn99.eq.2) then
          write(99,*)
          write(99,'(A,I3)')
     &     ' Vertical averages for outputs(outave) : klevv=', klevv
        endif
        do kk=1,4
          do k=ks2,ks1,-1
            if ((-zlim1(kk)).ge.zw(k)) then
              klevm(kk,1)=k
              goto 142
            endif
          enddo
142       continue
          do k=ks2,ks1,-1
            if ((-zlim2(kk)).ge.zw(k)) then
              klevm(kk,2)=k+1
              goto 143
            endif
          enddo
143       if (klevm(kk,2).eq.0) klevm(kk,2)=1

          if (nn99.eq.2) then
            write(99,*) 'Layer number',kk
            write(99,*) 'lim1 and lim2 :',zlim1(kk),zlim2(kk)
            write(99,*) 'lev1 and lev2 :',klevm(kk,1),klevm(kk,2)
          endif
        enddo
c  Number of value for each average
        do j=2,jmax-1
          do i=2,imax-1
            znumav(i,j,5)=0.0
            do k=klevv,ku2
              znumav(i,j,5)=znumav(i,j,5)+tmu(i,j,k)
            enddo
            znumav(i,j,5)=max(znumav(i,j,5),undim6)
          enddo
        enddo

        do j=ju1,ju2
          do i=iu1(j),iu2(j)
            do k=1,4
              znumav(i,j,k)=0.0
              do kk=klevm(k,2),klevm(k,1)
                znumav(i,j,k)=znumav(i,j,k)+tmu(i,j,kk)
              enddo
              znumav(i,j,k)=max(znumav(i,j,k),undim6)
            enddo
          enddo
        enddo

        do j=js1,js2
          do i=is1(j),is2(j)
            do k=1,4
              znumas(i,j,k)=0.0
              do kk=klevm(k,2),klevm(k,1)
                znumas(i,j,k)=znumas(i,j,k)+tms(i,j,kk)
              enddo
              znumas(i,j,k)=max(znumas(i,j,k),undim6)
            enddo
          enddo
        enddo

c
        spval=spvr
c
c end of preparation if numit=nstart (INITIALIZATION)
c
      endif
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2. Computation of averages                                          |
c-----------------------------------------------------------------------
c
c computation of annula mean surface velocity
      zfluxmt=zfluxmt+zflux0
      zfluxmts=zfluxmts+zflux0s
c
      njm=njm+1
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.2. CUMULATION OF GLOBAL MONTHLY OUTPUTS.
c--------------------------------------------
c streamfunction part
c cumulation for monthly and/or yearly outputs !
      if (nstreamout.ne.0) call streamfunc(0,0.d0,nstreamout)
c
cDFG start accumulation of ssh, ubar, vbar
c
        do j=1,jmax
          do i=1,imax
            etamo(i,j) = etamo(i,j) + eta(i,j)
             ubmo(i,j) =  ubmo(i,j) +  ub(i,j)
             vbmo(i,j) =  vbmo(i,j) +  vb(i,j)
          end do
        end do
cDFG end
c
cDFG start
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              umm(i,j,k) = umm(i,j,k)+u(i,j,k)
              vmm(i,j,k) = vmm(i,j,k)+v(i,j,k)
              tmm(i,j,k) = tmm(i,j,k)+scal(i,j,k,1)
              smm(i,j,k) = smm(i,j,k)+scal(i,j,k,2)
              uma(i,j,k) = uma(i,j,k)+u(i,j,k)
              vma(i,j,k) = vma(i,j,k)+v(i,j,k)
              tma(i,j,k) = tma(i,j,k)+scal(i,j,k,1)
              sma(i,j,k) = sma(i,j,k)+scal(i,j,k,2)
            enddo
          enddo
        enddo
cDFG end
c
cDFG start
        do k=1,kmax+1
          do j=1,jmax
            do i=1,imax
              wmm(i,j,k) = wmm(i,j,k)+w(i,j,k)
              wma(i,j,k) = wma(i,j,k)+w(i,j,k)
            enddo
          enddo
        enddo
cDFG end

        if (flagmon.or.flagyear) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          do j=2,jmax-1
            do i=2,imax-1
                zindh       = max(zero,sign(one,hgbq(i,j)*
     &                        (1.0-albq(i,j))-0.10))
                zinda       = max(zero,sign(one,
     &                         (1.0-albq(i,j))-0.10))
                zindb       = zindh*zinda
                cmoymo(i,j,1)= cmoymo(i,j,1)+hnbq(i,j)
                cmoymo(i,j,2)= cmoymo(i,j,2)+hgbq(i,j)
                cmoymo(i,j,3)= cmoymo(i,j,3)+hgbqp(i,j)
                cmoymo(i,j,4)= cmoymo(i,j,4)+albq(i,j)
                cmoymo(i,j,5)= cmoymo(i,j,5)+ts(i,j)
                cmoymo(i,j,6)= cmoymo(i,j,6)+fbbq(i,j)
                cmoymo(i,j,7)= cmoymo(i,j,7)+
     &                ug(i,j)*tmu(i,j,ks2)*zindb
                cmoymo(i,j,8)= cmoymo(i,j,8)+
     &                vg(i,j)*tmu(i,j,ks2)*zindb
                cmoymo(i,j,9) = cmoymo(i,j,9)+scal(i,j,ks2,1)
                cmoymo(i,j,10)= cmoymo(i,j,10)+scal(i,j,ks2,2)
                cmoymo(i,j,15)= cmoymo(i,j,15) - phiss(i,j,1)*
     &                         (rho0*cpo)/(dts(ks2)*unsdz(ks2))
     &                         +fcm1(i,j)
                cmoymo(i,j,16)= cmoymo(i,j,16)+86400./dts(ks2)*rho0*
     &                  ( phiss(i,j,2)/(unsdz(ks2)*34.7)
     &                    -phiss(i,j,0)
C    &                   +rappes(i,j,0)*dz(ks2)
C    &             *(scal(i,j,ks2,2)-scalr(i,j,ks2,2))/scal(i,j,ks2,2)
     &                  )
              cmoymo(i,j,17) = cmoymo(i,j,17) + ai(ks2) * sqrt(
     &          c4x(i,j,ks2)*c4x(i,j,ks2)+c4y(i,j,ks2)*c4y(i,j,ks2) )
                cmoymo(i,j,18)=cmoymo(i,j,18)+ eta(i,j)
Cage            cmoymo(i,j,18)=cmoymo(i,j,18)+ageg(i,j)
            enddo
          enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          do j=2,jmax-1
            do i=2,imax-1
              zztmp1=0.0
              zztmp2=0.0
              do kk=klevv,ku2
                zztmp1=zztmp1+tmu(i,j,kk)*u(i,j,kk)
                zztmp2=zztmp2+tmu(i,j,kk)*v(i,j,kk)
              enddo
              cmoymo(i,j,11)= cmoymo(i,j,11) +
     &                      zztmp1/znumav(i,j,5)
              cmoymo(i,j,12)= cmoymo(i,j,12) +
     &                      zztmp2/znumav(i,j,5)
            enddo
          enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          do j=js1,js2
            do  i=is1(j),is2(j)
              kmin3 = ks2
              ztest3 = 0.0
              kmin4 = ks2
              ztest4 = 0.0

              tloc=scal(i,j,ks2,1)-273.15
              sloc=scal(i,j,ks2,2)
              ccb1 = cstrho(4)*sloc
     &         + (cstrho(2)-cstrho(3)*tloc)*tloc
              ccb2 = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
              sigsurf=1000.0*(1.0/(cstrho(0)+
     &                 ccb2/(ccb1+cfb1z4(ks2+1)) ) - 1.0 )
C             write(155,*) sigsurf*tms(i,j,ks2)

              do 245 k=ks2,kfs(i,j),-1
                ztest3= ztest3+max(zero,sign(one,
     &                    bvf(i,j,k) ) )
                ztest3=min(ztest3,one)
                kmin3 = int(1.0-ztest3)*k+int(ztest3)*kmin3

                tloc=scal(i,j,k,1)-273.15
                sloc=scal(i,j,k,2)
                ccb1 = cstrho(4)*sloc
     &         +  (cstrho(2)-cstrho(3)*tloc)*tloc
                ccb2 = cstrho(5)
     &               + (cstrho(6) - cstrho(7)*tloc)*tloc
     &               - (cstrho(8) + cstrho(9)*tloc)*sloc
                sigz=1000.0*(1.0/(cstrho(0)+ ccb2
     &              /(ccb1+cfb1z4(ks2+1)) ) -1.0 )
                ztest4= ztest4+max(zero,sign(one,
     &                     -0.02+abs(sigz-sigsurf) ) )
                ztest4=min(ztest4,one)
                kmin4 = int(1.0-ztest4)*k+int(ztest4)*kmin4
  245         continue
              cmoymo(i,j,14)= cmoymo(i,j,14)-zw(kmin3-1)
              cmoymo(i,j,13)= cmoymo(i,j,13)-zw(kmin4)
            enddo
          enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        endif
c-------
        if (flagmon.or.flagyear) then
          do k=1,4
            kji=21+(k-1)*4
            do kk=klevm(k,2),klevm(k,1)
              do j=ju1,ju2
                do i=iu1(j),iu2(j)
                  cmoymo(i,j,kji)=cmoymo(i,j,kji)+
     &                tmu(i,j,kk)*u(i,j,kk)/znumav(i,j,k)
                  cmoymo(i,j,kji+1)=cmoymo(i,j,kji+1)+
     &                tmu(i,j,kk)*v(i,j,kk)/znumav(i,j,k)
                enddo
              enddo
            enddo
          enddo
          do k=1,4
            kji=19+(k-1)*4
            do kk=klevm(k,2),klevm(k,1)
              do j=js1,js2
                do i=is1(j),is2(j)
                  cmoymo(i,j,kji)=cmoymo(i,j,kji)+
     &               tms(i,j,kk)*scal(i,j,kk,1)/znumas(i,j,k)
                  cmoymo(i,j,kji+1)=cmoymo(i,j,kji+1)+
     &               tms(i,j,kk)*scal(i,j,kk,2)/znumas(i,j,k)
                enddo
              enddo
            enddo
          enddo
        endif
c-------
        if (flagmon.or.flagyear) then
          do j=js1,js2
            do  i=is1(j),is2(j)
              cmoymo(i,j,35)= cmoymo(i,j,35)+fcm1(i,j)
              cmoymo(i,j,36)= cmoymo(i,j,36)+tenagx(i,j)
              cmoymo(i,j,37)= cmoymo(i,j,37)+tenagy(i,j)
              cmoymo(i,j,38)= cmoymo(i,j,38)+fsolg(i,j)*
     &                   (1.0-albq(i,j))+albq(i,j)*fsolcn(i,j)
              cmoymo(i,j,39)= cmoymo(i,j,39)+(1.0-albq(i,j))
     &               *firg(i,j)+albq(i,j)*ratbqo(i,j)
              cmoymo(i,j,40)= cmoymo(i,j,40)+(1.0-albq(i,j))
     &               *fcsg(i,j)+albq(i,j)*fcscn(i,j)
              cmoymo(i,j,41)= cmoymo(i,j,41)+(1.0-albq(i,j))
     &               *fleg(i,j)+albq(i,j)*flecn(i,j)
              cmoymo(i,j,42)=cmoymo(i,j,42)+phiss(i,j,2)*rho0*
     &                  86400./(unsdz(ks2)*34.7*dts(ks2))
              cmoymo(i,j,43)= cmoymo(i,j,43)+albege(i,j)*
     &               (1-albq(i,j))+ albecn(i,j)*albq(i,j)
            enddo
          enddo
        endif
c-------
        if (flagmon.or.flagyear) then
          do j=js1,js2
            do  i=is1(j),is2(j)
              cmoymo(i,j,44)= cmoymo(i,j,44)+tabq(i,j)
              cmoymo(i,j,45)= cmoymo(i,j,45)+qabq(i,j)
              cmoymo(i,j,46)= cmoymo(i,j,46)+vabq(i,j)
              cmoymo(i,j,47)= cmoymo(i,j,47)+hnplbq(i,j)
            enddo
          enddo
        endif




c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- COMPUTES MONTHLY MEAN
c
      if (ijour.eq.njcum(jmois).and.ijour.ne.ijour1) then
c
        jjoufl = 1
        usnjm = 1.0/dble(njm)
        njm=0
c
        do 300 n=1,noumef
          do 300 j=1,jmax
            do 300 i=1,imax
              cmoymo(i,j,n)=cmoymo(i,j,n)*usnjm
  300   continue
c
        do jj=1,jmax
          do ij=1,imax
            zindp=max(0.0,sign(1.0D0,cmoymo(ij,jj,2)-0.2))
            cmoymo(ij,jj,6)=cmoymo(ij,jj,6)*zindp
          enddo
        enddo
c
        do 320 n=1,noumef
          do 320 j=1,jmax
            do 320 i=1,imax
              zindo       = tms(i,j,ks2)
              cmoymo(i,j,n)=zindo*cmoymo(i,j,n)+(1.0-zindo)*spvr
  320   continue
c
        do 325 kk=1,4
          do 325 j=1,jmax
            do 325 i=1,imax
              zindo1      = tms(i,j,klevm(kk,1))
              zindo2      = tmu(i,j,klevm(kk,1))
              n=19+(kk-1)*4
              cmoymo(i,j,n)=zindo1*cmoymo(i,j,n)+(1.0-zindo1)*spvr
              cmoymo(i,j,n+1)=zindo1*cmoymo(i,j,n+1)+(1.0-zindo1)*spvr
              cmoymo(i,j,n+2)=zindo2*cmoymo(i,j,n+2)+(1.0-zindo2)*spvr
              cmoymo(i,j,n+3)=zindo2*cmoymo(i,j,n+3)+(1.0-zindo2)*spvr
  325   continue
c
cDFG start computing of rest of 2d monthly means not yet contained
c    in cmoymo
c
        do j=1,jmax
          do i=1,imax
            etamo(i,j) = etamo(i,j)*usnjm
             ubmo(i,j) =  ubmo(i,j)*usnjm*0.01
             vbmo(i,j) =  vbmo(i,j)*usnjm*0.01
          end do
        end do
cDFG end
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c WRITES MONTHLY MEAN
cDFG start
c
        CALL open(Monthly_Means,refexp)
        if(flagmon) spval=missing_value
c
c Update time record. Data are dumped at the end of the month, the
c timestamp represents the entire averaging period. See also the
c timerecord definition with the 360 day calendar and the 30 day
c averaging period attributes.
c
        iMtimrec=iMtimrec+1
        if(flagmon) call incretime(ja,ijour,iMtimrec)

cDFG end
        if (mod(ja,nwm).eq.nwtest) then

c streamfunction part
          if (nstreamout.eq.1) call streamfunc(1,usnjm,nstreamout)
cDFG start streamfunction part

          do nb=0,nbsmax
            do j=1,jmtt
              do k=0,kmtt
                moc_m(nb,j,k) = uuu(j,k,nb)
                moc_a(nb,j,k) = moc_a(nb,j,k)+  moc_m(nb,j,k)
                uuu(j,k,nb) = 0.0
              end do
              do k=0,kfloor(j,nb)
                moc_m(nb,j,k) = spval
                moc_a(nb,j,k) = spval
              end do
              mht_m(nb,j) = vvv(j,nb)
              mst_m(nb,j) = vvv(j,nb+4)
              mht_a(nb,j) = mht_a(nb,j)+ mht_m(nb,j)
              mst_a(nb,j) = mst_a(nb,j)+ mst_m(nb,j)
              vvv(j,nb) = 0.0
              vvv(j,nb+4) = 0.0
            end do
          end do
c
c -- apply mask
          do nb=1,nbsmax
            do j=1,16
              do k=0,kmtt
                moc_m(nb,j,k) = spval
              end do
              mht_m(nb,j) = spval
              mst_m(nb,j) = spval
            end do
          end do

c
c dump of moc fields
c
          start(:)=1
          start(4)=iMtimrec
          count(1)=nbsmax+1
          count(2)=jmtt
          count(3)=kmtt
          count(4)=1
          IF (output(newtotvaro(meridional_overturning_streamfunction,2)))
     &       CALL writet(meridional_overturning_streamfunction,start,count,moc_m)

          start(3)=iMtimrec
          start(4)=0
          count(3)=1
          count(4)=0
          IF (output(newtotvaro(meridional_heat_transport,2)))
     &       CALL writed(meridional_heat_transport,start,count,mht_m)

          IF (output(newtotvaro(meridional_salt_transport,2)))
     &       CALL writed(meridional_salt_transport,start,count,mst_m)

c
c -- clean up
          do j=1,jmtt
            do nb=0,nbsmax
              do k=0,kmtt
                moc_m(nb,j,k) = 0.0
              end do
              mht_m(nb,j) = 0.0
              mst_m(nb,j) = 0.0
            end do
          end do

            do k=1,kmax
              do j=1,jmax
                do i=1,imax
                  umm(i,j,k)= tmu(i,j,k)*umm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spval
                  vmm(i,j,k)= tmu(i,j,k)*vmm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spval
                  tmm(i,j,k)= tms(i,j,k)*(tmm(i,j,k)*usnjm-273.15)
     &                     +(1-tms(i,j,k))*spval
                  smm(i,j,k)= tms(i,j,k)*smm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  wmm(i,j,k)= tms(i,j,k)*wmm(i,j,k)*usnjm*86400.
     &                   +(1-tms(i,j,k))*spval
                enddo
              enddo
            enddo
            k=kmax+1
            do j=1,jmax
              do i=1,imax
                wmm(i,j,k)= tms(i,j,k-1)*wmm(i,j,k)*usnjm*86400.
     &                 +(1-tms(i,j,k-1))*spval
              end do
            end do
c
c dump of 3d monthly means
c
            IF (output(newtotvaro(Potential_temperature,2)))
     &        CALL write_3d(Potential_temperature,tmm,kmax,iMtimrec)
            IF (output(newtotvaro(Salinity,2)))
     &        CALL write_3d(Salinity,smm,kmax,iMtimrec)

            IF (output(newtotvaro(Zonal_velocity,2)))
     &        CALL write_3d(Zonal_velocity,umm,kmax,iMtimrec)
            IF (output(newtotvaro(Meridional_velocity,2)))
     &        CALL write_3d(Meridional_velocity,vmm,kmax,iMtimrec)
            IF (output(newtotvaro(Vertical_velocity,2)))
     &        CALL write_3d(Vertical_velocity,wmm,kmax+1,iMtimrec)
c

            do k=1,kmax
              do j=1,jmax
                do i=1,imax
cDFG end
                  umm(i,j,k)= 0.0
                  vmm(i,j,k)= 0.0
                  tmm(i,j,k)= 0.0
                  smm(i,j,k)= 0.0
                enddo
              enddo
            enddo
            do k=1,kmax+1
              do j=1,jmax
                do i=1,imax
                  wmm(i,j,k)= 0.0
                end do
              end do
            end do
c

cDFG start output of 2d monthly means
c
c We need to care about the size of the floats defined for
c netcdf-output. Whereas everything in the netcdf-file is
c defined as double precision, cmoymo is single (even that
c may sometimes yield problems of accuracy...). Nevertheless,
c we someimes need to convert units anyway, so we simply
c copy all the data to the double precision work arry f2d
c first...
c
c ice model data:
c
c - snow cover above ice
c
            if (output(newtotvaro(Snow_thickness,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 1)*newtotvaro(Snow_thickness,4)+
     &                                  newtotvaro(Snow_thickness,5)
                end do
              end do
              call write_2d(Snow_thickness,999,f2d,.true.,.false.,iMtimrec)
            endif
c
c - snow precipitation
c
            if (output(newtotvaro(Snow_precipitation,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,47)*newtotvaro(Snow_precipitation,4)+
     &                                  newtotvaro(Snow_precipitation,5)
                end do
              end do
              call write_2d(Snow_precipitation,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - ice thickness
c
            if (output(newtotvaro(Ice_thickness,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 2)*newtotvaro(Ice_thickness,4)+
     &                                  newtotvaro(Ice_thickness,5)
                end do
              end do
              call write_2d(Ice_thickness,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - ice production
c
            if (output(newtotvaro(Ice_production,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 3)*newtotvaro(Ice_production,4)+
     &                                  newtotvaro(Ice_production,5)
                end do
              end do
              call write_2d(Ice_production,999,f2d,.true. ,.false.,iMtimrec)
            endif

c - lead fraction
c
            if (output(newtotvaro(Lead_fraction,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,4)*newtotvaro(Lead_fraction,4)+newtotvaro(Lead_fraction,5)
                end do
              end do
              call write_2d(Lead_fraction,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - ice temperature
c
            if (output(newtotvaro(Ice_temperature,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,5)*newtotvaro(Ice_temperature,4)+newtotvaro(Ice_temperature,5)
                end do
              end do
              call write_2d(Ice_temperature,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - oceanic heat flux at ice base
c
            if (output(newtotvaro(Heat_flux_ice_base,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 6)*newtotvaro(Heat_flux_ice_base,4)+
     &                                  newtotvaro(Heat_flux_ice_base,5)
                end do
              end do
              call write_2d(Heat_flux_ice_base,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - zonal ice velocity
c
            if (output(newtotvaro(Zonal_ice_velocity,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 7)*newtotvaro(Zonal_ice_velocity,4)+
     &                                  newtotvaro(Zonal_ice_velocity,5)
                end do
              end do
              call write_2d(Zonal_ice_velocity,999,f2d,.false.,.true. ,iMtimrec)
            endif
c
c -meridional ice velocity
c
            if (output(newtotvaro(Meridional_ice_velocity,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 8)*newtotvaro(Meridional_ice_velocity,4)+
     &                                  newtotvaro(Meridional_ice_velocity,5)
                end do
              end do
              call write_2d(Meridional_ice_velocity,999,f2d,.false.,.true. ,iMtimrec)
            endif
c
c 2-dimensional ocean (or ice/ocean) data:
c
c - sea surface height
c
            if (output(newtotvaro(Sea_surface_height,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=etamo(i,j)*newtotvaro(Sea_surface_height,4)+
     &                                  newtotvaro(Sea_surface_height,5)
                end do
              end do
              call write_2d(Sea_surface_height,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - ubar
c
            if (output(newtotvaro(Zonal_barotropic_momentum,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=ubmo(i,j)*newtotvaro(Zonal_barotropic_momentum,4)+
     &                                  newtotvaro(Zonal_barotropic_momentum,5)
                end do
              end do
              call write_2d(Zonal_barotropic_momentum,999,f2d,.false.,.true. ,iMtimrec)
            endif
c
c - vbar
c
            if (output(newtotvaro(Meridional_barotropic_momentum,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=vbmo(i,j)*newtotvaro(Meridional_barotropic_momentum,4)+
     &                                  newtotvaro(Meridional_barotropic_momentum,5)
                end do
              end do
              call write_2d(Meridional_barotropic_momentum,999,f2d,.false.,.true. ,iMtimrec)
            endif
c
c - sea surface temperature
c
            if (output(newtotvaro(Surface_temperature,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,9)*newtotvaro(Surface_temperature,4)+newtotvaro(Surface_temperature,5)
                end do
              end do
              call write_2d(Surface_temperature,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - sea surface salinity
c
            if (output(newtotvaro(Surface_salinity,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,10)*newtotvaro(Surface_salinity,4)+
     &                                  newtotvaro(Surface_salinity,5)
                end do
              end do
              call write_2d(Surface_salinity,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - depth of mixed layer
c
            if (output(newtotvaro(Depth_surface_mixed_layer,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,13)*newtotvaro(Depth_surface_mixed_layer,4)+
     &                                  newtotvaro(Depth_surface_mixed_layer,5)
                end do
              end do
              call write_2d(Depth_surface_mixed_layer,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - depth of convection
c
            if (output(newtotvaro(Depth_convection,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,14)*newtotvaro(Depth_convection,4)+
     &                                  newtotvaro(Depth_convection,5)
                end do
              end do
              call write_2d(Depth_convection,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - surface heat flux
c
            if (output(newtotvaro(Heat_flux_ice_base,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,15)*newtotvaro(Heat_flux_ice_base,4)+
     &                                  newtotvaro(Heat_flux_ice_base,5)
                end do
              end do
              call write_2d(Heat_flux_ice_base,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - surface freshwater flux
c
            if (output(newtotvaro(Surface_freshwater_flux,2))) then
c
c cmoymo(i,j,16) is in cm/day = 3.6m/year
c
              fac=1.0/3.6
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,16)*fac*newtotvaro(Surface_freshwater_flux,4)+
     &                                  newtotvaro(Surface_freshwater_flux,5)
                end do
              end do
              call write_2d(Surface_freshwater_flux,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - G-M slope parameter
c
            if (output(newtotvaro(G_M_slope,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,17)*newtotvaro(G_M_slope,4)+
     &                                  newtotvaro(G_M_slope,5)
                end do
              end do
              call write_2d(G_M_slope,999,f2d,.true. ,.false.,iMtimrec)
            endif
c
c - zonal windstress component
c   CLIO3 uses dynes/cm**2 for the stresses, but we want N/m**2 !
c
            if (output(newtotvaro(Zonal_wind_stress,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,36)*0.1*newtotvaro(Zonal_wind_stress,4)+
     &                                  newtotvaro(Zonal_wind_stress,5)
                end do
              end do
              call write_2d(Zonal_wind_stress,999,f2d,.false.,.true. ,iMtimrec)
            endif
c
c - meridional windstress component
c
            if (output(newtotvaro(Meridional_wind_stress,2))) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,37)*0.1*newtotvaro(Meridional_wind_stress,4)+
     &                                  newtotvaro(Meridional_wind_stress,5)
                end do
              end do
              call write_2d(Meridional_wind_stress,999,f2d,.false.,.true. ,iMtimrec)
            endif
c
c sync file (such that user may access it with analysis tools while
c model is still running)
c
            CALL close
c
c close netcdf-file at end of run or when record limit is reached
c
            if (iMtimrec.eq.maxmrecs) iMtimrec=0

        endif
c
c
c       UPDATE jmois.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.4. CUMULATION FOR GLOBAL YEARLY OUTPUTS AND ALL PERIOD MEAN
c---------------------------------------------------------------
c
        do 340 n=1,noumef
          do 340 j=1,jmax
            do 340 i=1,imax
              cmoymap(i,j,n,jmois)=cmoymap(i,j,n,jmois)+
     &                                 cmoymo(i,j,n)
  340   continue
        nmap(jmois)=nmap(jmois)+1

        if (mod(ja,nwa).eq.nwtest) then
c
          do 350 n=1,noumef
            do 350 j=1,jmax
              do 350 i=1,imax
                cmoyan(i,j,n)=cmoyan(i,j,n)+cmoymo(i,j,n)
  350     continue
c
cDFG start accumulation of 2d momentum and SSH
c
          do j=1,jmax
            do i=1,imax
              etaan(i,j) = etaan(i,j) + etamo(i,j)
               uban(i,j) =  uban(i,j) +  ubmo(i,j)
               vban(i,j) =  vban(i,j) +  vbmo(i,j)
            end do
          end do
cDFG end
        endif
c
c       UPDATE jmois.
        jmois = jmois+1
        if (jmois.gt.12) jmois = 1
cDFG start
c
c now that we have updated all annual sums, we can clear all
c monthly means
c
cDFG end
        do 343 n=1,noumef
          do 343 j=1,jmax
            do 343 i=1,imax
              cmoymo(i,j,n) =0.0
  343   continue
c-----
cDFG start cleanup
c
        do j=1,jmax
          do i=1,imax
            etamo(i,j) = 0.0
             ubmo(i,j) = 0.0
             vbmo(i,j) = 0.0
          end do
        end do
cDFG end

c-----
      endif




c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- WRITE GLOBAL YEARLY OUTPUTS.
c
      if (mod(ja,nwa).eq.nwtest.and.ijour.eq.int(yeaday)
     &               .and.ijour.ne.ijour1) then
c
        usanms = 1.0/real((12-(jmoisi-1)*(1-jmoifl)))
        jmoifl = 1
        do 370 n=1,noumef
          do 370 j=1,jmax
            do 370 i=1,imax
              zindo        = tms(i,j,ks2)
              cmoyan(i,j,n)=zindo*cmoyan(i,j,n)*usanms+(1.0-zindo)*spvr
  370   continue
c
cDFG start netcdf output of 2d annual fields
c
        do j=1,jmax
          do i=1,imax
            etaan(i,j) = etaan(i,j)*usanms
             uban(i,j) =  uban(i,j)*usanms
             vban(i,j) =  vban(i,j)*usanms
          end do
        end do

c
c Do we need to create a new NetCDF-file?
c generate a unique filename
c
        CALL open(Yearly_Means,refexp)
        if(flagyear) spval=missing_value
c
c Update time record. Data are dumped at the end of the year, the
c timestamp represents the entire averaging period. See also the
c timerecord definition with the 360 day calendar and the
c 360 day averaging period attributes.
c
        iAtimrec=iAtimrec+1
        if(flagyear) call incretime(ja,ijour,iAtimrec)
c
c ice model data:
c
c - snow cover above ice
c
          if (output(newtotvaro(Snow_thickness,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 1)*newtotvaro(Snow_thickness,4)+newtotvaro(Snow_thickness,5)
              end do
            end do
            call write_2d(Snow_thickness,999,f2d,.true.,.false.,iAtimrec)
          endif
c
c - snow precipitation
c
          if (output(newtotvaro(Snow_precipitation,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,47)*newtotvaro(Snow_precipitation,4)+newtotvaro(Snow_precipitation,5)
              end do
            end do
            call write_2d(Snow_precipitation,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - ice thickness
c
          if (output(newtotvaro(Ice_thickness,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 2)*newtotvaro(Ice_thickness,4)+newtotvaro(Ice_thickness,5)
              end do
            end do
            call write_2d(Ice_thickness,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - ice production
c
          if (output(newtotvaro(Ice_production,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 3)*newtotvaro(Ice_production,4)+newtotvaro(Ice_production,5)
              end do
            end do
            call write_2d(Ice_production,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - lead fraction
c
          if (output(newtotvaro(Lead_fraction,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,4)*newtotvaro(Lead_fraction,4)+newtotvaro(Lead_fraction,5)
              end do
            end do
            call write_2d(Lead_fraction,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - ice temperature
c
          if (output(newtotvaro(Ice_temperature,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,5)*newtotvaro(Ice_temperature,4)+newtotvaro(Ice_temperature,5)
              end do
            end do
            call write_2d(Ice_temperature,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - oceanic heat flux at ice base
c
          if (output(newtotvaro(Heat_flux_ice_base,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 6)*newtotvaro(Heat_flux_ice_base,4)+newtotvaro(Heat_flux_ice_base,5)
              end do
            end do
            call write_2d(Heat_flux_ice_base,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - zonal ice velocity
c
          if (output(newtotvaro(Zonal_ice_velocity,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 7)*newtotvaro(Zonal_ice_velocity,4)+newtotvaro(Zonal_ice_velocity,5)
              end do
            end do
            call write_2d(Zonal_ice_velocity,999,f2d,.false.,.true. ,iAtimrec)
          endif
c
c -meridional ice velocity
c
          if (output(newtotvaro(Meridional_ice_velocity,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 8)*newtotvaro(Meridional_ice_velocity,4)+
     &                                  newtotvaro(Meridional_ice_velocity,5)
              end do
            end do
            call write_2d(Meridional_ice_velocity,999,f2d,.false.,.true. ,iAtimrec)
          endif
c
c 2-dimensional ocean (or ice/ocean) data:
c
c - sea surface height
c
          if (output(newtotvaro(Sea_surface_height,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=etaan(i,j)*newtotvaro(Sea_surface_height,4)+
     &                                  newtotvaro(Sea_surface_height,5)
              end do
            end do
            call write_2d(Sea_surface_height,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - ubar
c
          if (output(newtotvaro(Zonal_barotropic_momentum,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=uban(i,j)*newtotvaro(Zonal_barotropic_momentum,4)+
     &                                  newtotvaro(Zonal_barotropic_momentum,5)
              end do
            end do
            call write_2d(Zonal_barotropic_momentum,999,f2d,.false.,.true. ,iAtimrec)
          endif
c
c - vbar
c
          if (output(newtotvaro(Meridional_barotropic_momentum,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=vban(i,j)*newtotvaro(Meridional_barotropic_momentum,4)+
     &                                  newtotvaro(Meridional_barotropic_momentum,5)
              end do
            end do
            call write_2d(Meridional_barotropic_momentum,999,f2d,.false.,.true. ,iAtimrec)
          endif
c
c - sea surface temperature
c
          if (output(newtotvaro(Surface_temperature,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 9)*newtotvaro(Surface_temperature,4)+newtotvaro(Surface_temperature,5)
              end do
            end do
            call write_2d(Surface_temperature,999,f2d,.true.,.false.,iAtimrec)
          endif
c
c - sea surface salinity
c
          if (output(newtotvaro(Surface_salinity,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,10)*newtotvaro(Surface_salinity,4)+
     &                                  newtotvaro(Surface_salinity,5)
              end do
            end do
            call write_2d(Surface_salinity,999,f2d,.true.,.false.,iAtimrec)
          endif
c
c - depth of mixed layer
c
          if (output(newtotvaro(Depth_surface_mixed_layer,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,13)*newtotvaro(Depth_surface_mixed_layer,4)+
     &                                  newtotvaro(Depth_surface_mixed_layer,5)
              end do
            end do
            call write_2d(Depth_surface_mixed_layer,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - depth of convection
c
          if (output(newtotvaro(Depth_convection,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,14)*newtotvaro(Depth_convection,4)+
     &                                  newtotvaro(Depth_convection,5)
              end do
            end do
            call write_2d(Depth_convection,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - surface heat flux
c
          if (output(newtotvaro(Surface_heat_flux,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,15)*newtotvaro(Surface_heat_flux,4)+
     &                                  newtotvaro(Surface_heat_flux,5)
              end do
            end do
            call write_2d(Surface_heat_flux,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - surface freshwater flux
c
          if (output(newtotvaro(Surface_freshwater_flux,3))) then
            fac=1.0/3.6
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,16)*fac*newtotvaro(Surface_freshwater_flux,4)+
     &                                  newtotvaro(Surface_freshwater_flux,5)
              end do
            end do
            call write_2d(Surface_freshwater_flux,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - G-M slope parameter
c
          if (output(newtotvaro(G_M_slope,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,17)*newtotvaro(G_M_slope,4)+
     &                                  newtotvaro(G_M_slope,5)
              end do
            end do
            call write_2d(G_M_slope,999,f2d,.true. ,.false.,iAtimrec)
          endif
c
c - zonal windstress component
c
          if (output(newtotvaro(Zonal_wind_stress,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,36)*newtotvaro(Zonal_wind_stress,4)+
     &                                  newtotvaro(Zonal_wind_stress,5)
              end do
            end do
            call write_2d(Zonal_wind_stress,999,f2d,.false.,.true. ,iAtimrec)
          endif
c
c - meridional windstress component
c
          if (output(newtotvaro(Meridional_wind_stress,3))) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,37)*newtotvaro(Meridional_wind_stress,4)+
     &                                  newtotvaro(Meridional_wind_stress,5)
              end do
            end do
            call write_2d(Meridional_wind_stress,999,f2d,.false.,.true. ,iAtimrec)
          endif
c
c cleanup of annual files that have been written to disk meanwhile
c
        do j=1,jmax
          do i=1,imax
            etaan(i,j) = 0.0
             uban(i,j) = 0.0
             vban(i,j) = 0.0
          end do
        end do
cDFG end netcdf output
c
        do 390 n=1,noumef
          do 390 j=1,jmax
            do 390 i=1,imax
              cmoyan(i,j,n)=0.0
  390   continue

        zfluxm=zfluxmt/yeaday
        zfluxms=zfluxmts/yeaday
        write(42,*) zfluxm,zfluxms
        zfluxmt=0.0
        zfluxmts=0.0

c streamfunction part
        if (nstreamout.eq.2) call streamfunc(1,usanms*usnjm,nstreamout)
cDFG start streamfunction part
c
c -- time average
        usnja=1./12.
        do j=1,jmtt
          do nb=0,nbsmax
            do k=0,kmtt
              moc_a(nb,j,k) =  moc_a(nb,j,k)*usnja
            end do
            do k=0,kfloor(j,nb)
              moc_a(nb,j,k) = spval
            end do
            mht_a(nb,j) = mht_a(nb,j)*usnja
            mst_a(nb,j) = mst_a(nb,j)*usnja
          end do
        end do
c
c -- apply mask
        do nb=1,nbsmax
          do j=1,16
            do k=0,kmtt
              moc_a(nb,j,k) = spval
            end do
            mht_a(nb,j) = spval
            mst_a(nb,j) = spval
          end do
        end do
c
c --write to netcdf file

        start(:)=1
        start(4)=iAtimrec
        count(1)=nbsmax+1
        count(2)=jmtt
        count(3)=kmtt
        count(4)=1
        IF (output(newtotvaro(meridional_overturning_streamfunction,3)))
     &     CALL writet(meridional_overturning_streamfunction,start,count,moc_a)

        start(3)=iAtimrec
        start(4)=0
        count(3)=1
        count(4)=0
        IF (output(newtotvaro(meridional_heat_transport,3)))
     &     CALL writed(meridional_heat_transport,start,count,mht_a)

        IF (output(newtotvaro(meridional_salt_transport,3)))
     &     CALL writed(meridional_salt_transport,start,count,mst_a)

c
c -- wipe arrays
        do j=1,jmtt
          do nb=0,nbsmax
            do k=0,kmtt
               moc_a(nb,j,k) = 0.0
            end do
            mht_a(nb,j) = 0.0
            mst_a(nb,j) = 0.0
          end do
        end do
cDFG end

cDFG start computing annual mean 3d fields
          usnja=1/yeaday
          do k=1,kmax
            do j=1,jmax
              do i=1,imax
                uma(i,j,k)= tmu(i,j,k)*uma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spval
                vma(i,j,k)= tmu(i,j,k)*vma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spval
                tma(i,j,k)= tms(i,j,k)*(tma(i,j,k)*usnja-273.15)
     &                     +(1-tms(i,j,k))*spval
                sma(i,j,k)= tms(i,j,k)*sma(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                wma(i,j,k)= tms(i,j,k)*wma(i,j,k)*usnja*86400.
     &                       +(1-tms(i,j,k))*spval
              end do
            end do
          end do
          k=kmax+1
          do j=1,jmax
            do i=1,imax
              wma(i,j,k)= tms(i,j,k-1)*wma(i,j,k)*usnja*86400.
     &                   +(1-tms(i,j,k-1))*spval
            end do
          end do
c
c netcdf output of 3d annual means
c

          IF (output(newtotvaro(Potential_temperature,3)))
     &      CALL write_3d(Potential_temperature,tma,kmax,iAtimrec)
          IF (output(newtotvaro(Salinity,3)))
     &      CALL write_3d(Salinity,sma,kmax,iAtimrec)

          IF (output(newtotvaro(Zonal_velocity,3)))
     &      CALL write_3d(Zonal_velocity,uma,kmax,iAtimrec)
          IF (output(newtotvaro(Meridional_velocity,3)))
     &      CALL write_3d(Meridional_velocity,vma,kmax,iAtimrec)
          IF (output(newtotvaro(Vertical_velocity,3)))
     &      CALL write_3d(Vertical_velocity,wma,kmax+1,iAtimrec)
c
          do k=1,kmax
            do j=1,jmax
              do i=1,imax
cDFG end
                uma(i,j,k)= 0.0
                vma(i,j,k)= 0.0
                tma(i,j,k)= 0.0
                sma(i,j,k)= 0.0
              enddo
            enddo
          enddo
cDFG start
          do k=1,kmax+1
            do j=1,jmax
              do i=1,imax
                wma(i,j,k)= 0.0
              end do
            end do
          end do
c
        call close
c
c close netcdf-file at end of run or when record limit is reached
c
        if (iAtimrec.eq.maxarecs) iAtimrec=0
cDFG end
c-----
      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- WRITE GLOBAL ALL PERIOD OUTPUTS.
      if (numit.eq.nlast.and.nwtal.eq.1) then

        do 410 k=1,12
          do 400 n=1,noumef
            do 400 j=1,jmax
              do 400 i=1,imax
                zindo        = tms(i,j,ks2)
                cmoymap(i,j,n,k)=zindo*cmoymap(i,j,n,k)
     &                          /max(dfloat(nmap(k)),undim6)
     &                          +(1.0-zindo)*spvr
  400     continue
  410   continue
      endif
c
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine outnet -
      end
c
cDFG end
