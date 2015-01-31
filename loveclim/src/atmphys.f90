!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine iatmphys
!-----------------------------------------------------------------------
! *** initializes variables used in subroutines of atmphys.f
!-----------------------------------------------------------------------
      USE NETCDF    ! netcdf module to read forcing data

      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'netcdf.inc'

      real*8 fvolcan
      integer ios

      integer i,j,k,l,ireg,im,nn,is,j1,i1,ii,jj,ism
      real*8  beta,draganr,draglar,dum(2),asum,spv
      integer jyear,kyear,ilat,jmonth,m,indxvol,indxtsi,status
      real*8 tsi,ksw,valVolc1,valVolc2,valVolc3,valVolc4
      character*6 numyear
      character*3 numday
      integer tmp_imonth
      real*8 globalmean

      !write(numyear,'(i6.6)') irunlabel+int((irunlabeld)/360)
      !write(numday,'(i3.3)') mod(irunlabeld,360)
      write(numyear,'(i6.6)') irunlabel
      write(numday,'(i3.3)') irunlabeld
!
! *** initial atmospheric temperatures and moisture
!
      if (initfield.eq.0) then
        tempm(0)=tzero-35d0
        tempm(1)=tzero-35d0
        tempm(2)=tzero-8d0
        do j=1,nlon
          do i=1,nlat
            temp0g(i,j)=tempm(0)
            temp2g(i,j)=tempm(1)
            temp4g(i,j)=tempm(2)
            tempsg(i,j)=290d0
            rmoisg(i,j)=0d0
            relhum(i,j)=0d0
            q10(i,j)=0d0
            qsurf(i,j)=0d0
            do nn=1,ntyps
              q10n(i,j,nn)=0.0d0
              qsurfn(i,j,nn)=0.0d0
              tempsgn(i,j,nn)=290d0
              pgroundn(i,j,nn)=p0
            enddo
            pground(i,j)=p0
            geopg(i,j,2)=0d0
            tcc(i,j)=ccisccp(i,j,1)
          enddo
        enddo


        do j=1,nlon
          do i=1,nlat
            tsurfn (i,j,noc)= tempsg(i,j)
            tsurfn (i,j,nse)= tzero
            tsurfn (i,j,nld)= tempsg(i,j)
            tsurf (i,j)= tempsg(i,j)
          enddo
        enddo
        call atmphyszero
      else
	ios=0
        open(iuo+95,file='startdata/inatphy'//numyear//'_'//numday//'.dat', &
     &        form='unformatted')
        read(iuo+95) tsurfn,tempm,temp0g
        read(iuo+95) rmoisg,torain,tosnow
   99   if((ios.ne.0).and.((iscenghg.ne.3).and.(iscenghg.ne.0))) &
     &              print *,'CO2 ref from value in emic.f:',PCO2ref
        close(iuo+95)
      endif

      do j=1,nlon
        do i=1,nlat
          rmountn(i,j,nld)=rmount(i,j)
          rmountn(i,j,nse)=0d0
          rmountn(i,j,noc)=0d0
          if (rmountn(i,j,nld).lt.0d0) rmountn(i,j,nld)=0d0
          if (fractn(i,j,nld).lt.epss) then
            rmountn(i,j,nld)=0d0
            qmount(i,j)=0d0
          else
            qmount(i,j)=rmountn(i,j,nld)
          endif
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
          asum=0d0
          do i1=-1,1
            do j1=-1,1
              ii=i+i1
              jj=j+j1
              if (ii.lt.1) then
                ii=1
                jj=jj+nlon/2
              endif
              if (ii.gt.nlat) then
                ii=nlat
                jj=jj+nlon/2
              endif
              if (jj.lt.1) jj=jj+nlon
              if (jj.gt.nlon) jj=jj-nlon
              asum=asum+rmountn(ii,jj,nld)
            enddo
          enddo
          qmount(i,j)=asum/9d0
        enddo
      enddo


!***  longwave radiation parameterisation, based on a linearisation
!***  of KRCM with respect to reference T and q profiles and other
!***  variables (greenhousegases)

      read(iuo+16) irn,ipl,pisccp,pncep,z500ncep
      read(iuo+16) tncep,qancep,ghgipcc,ccisccp
      read(iuo+16) lwrref

! *** UPDATE land surface each year

      read(iuo+17) lwrt,lwrts,lwrqts,lwrqa,lwrghg

!**   amplifcation of freshwater feedback
      lwrqa(:,:,:,:)=lwrqa(:,:,:,:)*AMPWIR
!**   amplifcation of feedback at the equator
      lwrqa(:,9,:,:)=lwrqa(:,9,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,10,:,:)=lwrqa(:,10,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,11,:,:)=lwrqa(:,11,:,:)*AMPEQIR
      lwrqa(:,12,:,:)=lwrqa(:,12,:,:)*AMPEQIR
      lwrqa(:,13,:,:)=lwrqa(:,13,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,14,:,:)=lwrqa(:,14,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,23,:,:)=lwrqa(:,23,:,:)*AMPANIR
      lwrqa(:,1,:,:)=lwrqa(:,1,:,:)*AMPANIR2
      lwrqa(:,2,:,:)=lwrqa(:,2,:,:)*AMPANIR2
      lwrqa(:,3,:,:)=lwrqa(:,3,:,:)*AMPANIR2

!     write(iuo+99,*) "Modif IR scheme"
!     write(iuo+99,*) AMPWIR,AMPEQIR,expIR
!     write(iuo+99,*) HPROFW,HPROFTROP,1.0+(HPROFTROP-1.0)/2.0,HPROFEQ

!***  update moisture profile used in the linearization of the
!***  radiative scheme in the tropics:easy surrogate to a change
!***  mean IR flux in the model without affecting sensitivity
      do ism=1,12
        qancep(1,ism)=qancep(1,ism)*HPROFAN2
        qancep(2,ism)=qancep(2,ism)*HPROFAN2
        qancep(3,ism)=qancep(3,ism)*HPROFAN2
        qancep(4,ism)=qancep(4,ism)*HPROFW
        qancep(5,ism)=qancep(5,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(6,ism)=qancep(6,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(7,ism)=qancep(7,ism)*HPROFTROP
        qancep(8,ism)=qancep(8,ism)*HPROFTROP
        qancep(9,ism)=qancep(9,ism)*HPROFTROP
        qancep(10,ism)=qancep(10,ism)*HPROFTROP

        qancep(11,ism)=qancep(11,ism)*HPROFEQ
        qancep(12,ism)=qancep(12,ism)*HPROFEQ

        qancep(13,ism)=qancep(13,ism)*HPROFTROP
        qancep(14,ism)=qancep(14,ism)*HPROFTROP
        qancep(15,ism)=qancep(15,ism)*HPROFTROP
        qancep(16,ism)=qancep(16,ism)*HPROFTROP
        qancep(17,ism)=qancep(17,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(18,ism)=qancep(18,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(19,ism)=qancep(19,ism)*HPROFW
        qancep(20,ism)=qancep(20,ism)*HPROFW
        qancep(21,ism)=qancep(21,ism)*HPROFW
        qancep(22,ism)=qancep(22,ism)*HPROFW
        qancep(23,ism)=qancep(23,ism)*HPROFAN
        qancep(24,ism)=qancep(24,ism)*HPROFW
        qancep(25,ism)=qancep(25,ism)*HPROFW
        qancep(26,ism)=qancep(26,ism)*HPROFW
        qancep(27,ism)=qancep(27,ism)*HPROFW
      enddo


!read GHG concentrations
      ghgscen(:,:)=0.0
      i=1
      k=1
      do
        read(iuo+33,*,iostat=k) (ghgscen(j,i),j=1,20)
        if(k.lt.0) exit
        i=i+1
      enddo
      y1scenghg=ghgscen(1,1)
      nyscenmaxghg=i-1
      if (nyscenmaxghg.eq.0) nyscenmaxghg=1

      if (iscenghg.eq.2) then
!2 times CO2 concentrations
!-> other GHGs at their "standard" concentrations (1st line of file)
        i=iscenghg2s
        k=1
        do
          read(iuo+39,*,iostat=k) j, ghgscen(2,i)
          if(k.lt.0) exit
          i=i+1
        enddo
        do k=i,nyscenmaxghg
          ghgscen(2,k)=ghgscen(2,i-1)
        enddo
        i=iscenghg2s+1
        do k=i,nyscenmaxghg
          do j=3,20
            ghgscen(j,k)=ghgscen(j,iscenghg2s)
          enddo
        enddo
      endif
      write(iuo+99,*) 'scen GHG start=',y1scenghg, "AD nbline=",nyscenmaxghg

!read O3 concentrations
      o3scen(1,:)=0 !year
      o3scen(2,:)=25.0 !value
      if (isceno3.eq.1) then
        i=1
        k=1
        do
          read(iuo+37,*,iostat=k) (o3scen(j,i),j=1,2)
          if(k.lt.0) exit
          i=i+1
        enddo
        y1sceno3=o3scen(1,1)
        nyscenmaxo3=i-1
        if (nyscenmaxo3.eq.0) nyscenmaxo3=1
        write(iuo+99,*) 'scen O3 start=',y1sceno3, "AD nbline=",nyscenmaxo3
      endif

!read sulfates optical depths
      sulopt(:,:,:)=0.0
      if (iscensul.eq.1) then
        status=nf_open("inputdata/SUL.nc", nf_nowrite, ireg) !ouvre le fichier
        status=nf_inq_dimid(ireg, 'time', i) !recupère l'id du temps
        status=nf_inq_dimlen(ireg, i, nyscenmaxsul) !recupère à partir de l'id, le nombre de pas de temps
        status=nf_inq_varid(ireg, 'time', i) !recupère l'id du temps
        status=nf_get_vara_double(ireg, i, (/1/), (/nyscenmaxsul/), suloptTime)
        issulstrt=int(suloptTime(1)) !la 1ere valeur du temps donne la date de debut du forcage
        status=nf_inq_varid(ireg, "Sul", j) !récupère l'id de la variable Sul
        status=nf_get_vara_real(ireg, j, (/1,1,1/), (/64,32,nyscenmaxsul/), sulopt) !charge les valeurs dans la variable sulopt
        oldmonth=-1
        write(iuo+99,*) 'scen Sul start=',issulstrt, "AD nbline(month)=",nyscenmaxsul
      endif

!***  read 0-2000 TSI anomalies

      if (iscentsi.eq.1) then
        tsiscen(:,:)=0.0 !utilisé lorsque l'on demande un forçage précédent le debut du fichier
        i=1
        k=1
        do
         read(iuo+34,*,iostat=k) (tsiscen(j,i),j=1,2)
         if(k.lt.0) exit
         i=i+1
        enddo
        y1scentsi=tsiscen(1,1)
        nyscenmaxtsi=i-1
	if (nyscenmaxtsi.eq.0) nyscenmaxtsi=1
        write(iuo+99,*) 'scen TSI start=',y1scentsi, "AD nbline=",nyscenmaxtsi
      else
        tsiscen(:,:)=0.0
      endif

!***  read 0-2000 anomalies associated with volcanos
      solarvol(:,:,:)=0.
      if (iscenvol.eq.1) then
!       changement fait par Yoann Sallaz-Damaz pour permettre l'utilisation
!       d'un fichier de forçage volcanique évolutif au cours d'une assimilation
        i=1
        k=1
        do
          read(iuo+35,'(I5,I3,1X,4(F8.3,1X))',iostat=k) jyear,jmonth,valVolc1,valVolc2,valVolc3,valVolc4
          solarvol(i,jmonth,1)=valVolc1
          solarvol(i,jmonth,2)=valVolc2
          solarvol(i,jmonth,3)=valVolc3
          solarvol(i,jmonth,4)=valVolc4
          if(i.eq.1) then
            y1scenvol=jyear
            m1scenvol=jmonth
          endif
          if(k.lt.0) exit
          if(jmonth.eq.12) i=i+1
        enddo
        nyscenmaxvol=i-1
	if (nyscenmaxvol.eq.0) nyscenmaxvol=1
        write(iuo+99,*) 'scen VOLC starty=',y1scenvol,"AD startm=",m1scenvol, "nbline=",nyscenmaxvol
      endif
100   FORMAT(5X,3(X,F7.2),9(2X,F6.2),X,F7.2,6(2X,F6.2))


      i=1
      call ghgupdate(i)

!* set reference value for CO2 concentration

      do ireg=1,27
        do im=1,12
          beta=(tncep(11,ireg,im)-tncep(10,ireg,im))/alog(400./300.)
          tncep12(1,ireg,im)=tncep(10,ireg,im)+beta*alog(350./300.)
          beta=(tncep(14,ireg,im)-tncep(13,ireg,im))/alog(700./600.)
          tncep12(2,ireg,im)=tncep(13,ireg,im)+beta*alog(650./600.)
        enddo
      enddo

      do k=1,17
         rlogtl(k)=log(pncep(k)/65000.)
      enddo

      do ireg=1,27
        rlogts(ireg)=log(pisccp(ireg)/65000.)
      enddo


!
! *** shortwave radiation parameters
!


      read(iuo+18) costref,salbref
      read(iuo+18) swrref


      read(iuo+19) swrcost,swrsalb



      call detqmtabel


310   format(9f7.2)
330   format(2f5.2)
340   format(f5.2)


! *** computation of the effective turning angle


      draglar=3.141592654/180.0*dragla
      draganr=3.141592654/180.0*dragan
      do i=1,nlat
        if (phi(i).lt.(-1.0*draglar)) then
           dragane(i)=-1.0*draganr
        else
           if (phi(i).gt.draglar) then
               dragane(i)=draganr
           else
               dragane(i)=draganr*phi(i)/draglar
           endif
        endif
      enddo

!
! *** evaporation factor
!
      evfac=1d0

      ksw=0.042

      do i=1,32
        solarcl(i)=solarm
      enddo
      indxtsi=1
      indxvol=1
      if (iscentsi.eq.1) then
        indxtsi=irunlabelf+iyear-y1scentsi+1
        if ((initialization.eqv..true.).and.(irunlabeld.eq.360)) indxtsi=indxtsi+1
        if (indxtsi.gt.nyscenmaxtsi) indxtsi=nyscenmaxtsi
        if (indxtsi.lt.1) indxtsi=1
        do i=1,32
          solarcl(i)=solarcl(i)+tsiscen(2,indxtsi)*facttsi
        enddo
        write(iuo+99,*) 'TSI y=',tsiscen(1,indxtsi), "AD values=",tsiscen(1:2,indxtsi)
      endif

      if (iscenvol.eq.1) then
	tmp_imonth=imonth
        !indxvol=irunlabelf
        indxvol=irunlabelf+iyear-y1scenvol+1
        if (initialization.eqv..true.) then
	  if (irunlabeld.eq.360) then
	    indxvol=indxvol+1
	    imonth=1
	  else
	    if ((mod(irunlabeld,30)).eq.0) imonth=imonth+1
	  endif
	endif

        if (indxvol.gt.nyscenmaxvol) indxvol=nyscenmaxvol
        if (indxvol.lt.1) indxvol=1
        do i=1,10
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,1)
        enddo
        do i=11,16
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,2)
        enddo
        do i=17,22
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,3)
        enddo
        do i=23,32
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,4)
        enddo
	imonth=tmp_imonth
      endif

      if (isceno3.eq.1) then
        do i=1,16
        solarcl(i)=solarcl(i)+(4*0.247*ksw*(o3-25.))
        enddo
        do i=17,32
        solarcl(i)=solarcl(i)+(4*0.324*ksw*(o3-25.))
        enddo
      endif

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ghgupdate(istep)
!-----------------------------------------------------------------------
! *** updates ghg concentrations: indxghg 1 corresponds to y1scenghg AD
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'


      integer i,istep,indxghg,s,r,k,l,m,indxo3,h
      real*8  logco2,sqrch4,sqrn2o
      real*8 alpho3lw(2)

      if((mod(nint(day*real(iatm)),nstpyear).eq.0).or.(initialization.eqv..true.)) then
      !if(mod(istep,nstpyear).eq.1) then

      write(iuo+99,*) 'iscenghg',iscenghg

      indxghg=0
      if ((iscenghg.eq.1).or.(iscenghg.eq.2)) then
         indxghg=irunlabelf+iyear-y1scenghg+1
         if ((initialization.eqv..true.).and.(irunlabeld.eq.360)) indxghg=indxghg+1
         if (indxghg.gt.nyscenmaxghg) indxghg=nyscenmaxghg
      else
         indxghg=1
      endif
      if (indxghg.lt.1) then
        write(*,*) "GHG error :"
        write(*,*) "your GHG.dat file doesn't include the simulation period"
        stop
      endif
      write(iuo+99,*) 'GHG y=',ghgscen(1,indxghg), "AD values=",ghgscen(2:4,indxghg)

      ghg(1:19)=ghgscen(2:20,indxghg)

        PGACO2=ghg(1)

      indxo3=1
      if (isceno3.eq.1) then
         indxo3=irunlabelf+iyear-y1sceno3+1
         if ((initialization.eqv..true.).and.(irunlabeld.eq.360)) indxo3=indxo3+1
         if (indxo3.gt.nyscenmaxo3) indxo3=nyscenmaxo3
      endif

      if (indxo3.lt.1) then
        indxo3=1
      endif

      o3=o3scen(2,indxo3)
!      write(iuo+99,*) 'O3 y=',o3scen(1,indxo3),"AD value=", o3
!10     format(A12,1i6,4f12.3)
!      write(iuo+99,*) 'Aure - indxo3', indxo3
      write(iuo+99,*) 'O3 value', o3-25.
       call flush(iuo+99)
!*** Update LW reference radiation fluxes using new GHG concentrations
        logco2=log(ghg(1)/ghgipcc(1))
        !write(*,*) "logco2=log(",ghg(1),"/",ghgipcc(1),")"
        if(ghg(1).eq.0) stop
        sqrch4=sqrt(ghg(2))-sqrt(ghgipcc(2))
        sqrn2o=sqrt(ghg(3))-sqrt(ghgipcc(3))
        alpho3lw(1)=153.6
        alpho3lw(2)=201.2
        do h=1,2
        do l=0,1
         do s=1,4
          do r=1,27
           do k=1,7
            lwrflux(k,r,s,l,h)=lwrref(k,r,s,l)+ lwrghg(k,1,r,s,l)*logco2+ &
                 & lwrghg(k,2,r,s,l)*sqrch4+ lwrghg(k,3,r,s,l)*sqrn2o
!                   write(*,*) "lwrfluxA+=",lwrref(k,r,s,l),"+",lwrghg(k,1,r,s,l),"*",logco2,"+",lwrghg(k,2,r,s,l),"*",sqrch4,"+",
!      *                 lwrghg(k,3,r,s,l),"*",sqrn2o
            do m=4,19
             lwrflux(k,r,s,l,h)=lwrflux(k,r,s,l,h)+ &
                  & lwrghg(k,m,r,s,l)*(ghg(m)-ghgipcc(m))
!                   write(*,*) "lwrfluxB+=",lwrflux(k,r,s,l,h),"+",lwrghg(k,m,r,s,l),"*(",ghg(m),"-",ghgipcc(m),")"
            enddo
              lwrflux(k,r,s,l,h)=lwrflux(k,r,s,l,h)+ &
                   & lwrghg(k,4,r,s,l)*alpho3lw(h)*(o3-25.)
!                   write(*,*) "lwrfluxC+=",lwrghg(k,4,r,s,l),"*",alpho3lw(h),"*(",o3,"-",25.,")"
!            write(iuo+99,*) 'Aure2 - o3 moins 25', o3-25.
           enddo
          enddo
         enddo
        enddo
        enddo
      endif
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmphyszero
!-----------------------------------------------------------------------
! *** initializes data arrays to zero for physics routines
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comrunlabel.h'

      integer i,j


      do j=1,nlon
        do i=1,nlat
          cormois(i,j)=0.d0
          torain(i,j)=0.d0
          tosnow(i,j)=0.d0
          dyrain(i,j)=0.d0
          corain(i,j)=0.d0
          dysnow(i,j)=0.d0
          cosnow(i,j)=0.d0
          thforg0(i,j)=0.d0
          thforg1(i,j)=0.d0
          thforg2(i,j)=0.d0
          vhforg0(i,j)=0.d0
          vhforg1(i,j)=0.d0
          vhforg2(i,j)=0.d0
        enddo
      enddo


      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sensrad
!-----------------------------------------------------------------------
! *** computes atmospheric forcing due to sensible heat and radiation
! *** the forcing terms are computed in dimensional units (K/s)
! *** input
! ***       hflux : sensible heat flux between atmosphere and earth
! ***       hesws : short wave solar radiation in layer 1 or 2
! ***       ulrad : upward long wave radiation in layer 1 or 2
! *** output
! ***       thforg: temperature forcing in layer 1 or 2 in K/s
! ***       vhforg: diabatic forcing in layer 1 or 2 in K/s
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  halpha,halpha1,halpha2,sum1,sum2,sum0,halpha0

      halpha=grav/cpair


      halpha0 =halpha/dp0
      halpha1 =halpha/dp1
      halpha2 =halpha/dp2



!
! *** summation of forcing terms
!
      do j=1,nlon
        do i=1,nlat


          sum0 = (hesw0(i,j) - ulrad0(i,j) + ulrad1(i,j))*halpha0
          sum1 = (hesw1(i,j) - ulrad1(i,j) + ulrad2(i,j))*halpha1
          sum2 = (hflux(i,j) + hesw2(i,j) - ulrad2(i,j) + &
               & ulrads(i,j) - dlrads(i,j))*halpha2



          thforg0(i,j) = thforg0(i,j) + sum0
          thforg1(i,j) = thforg1(i,j) + sum1
          thforg2(i,j) = thforg2(i,j) + sum2


          vhforg0(i,j) = vhforg0(i,j) + sum0
          vhforg1(i,j) = vhforg1(i,j) + sum1
          vhforg2(i,j) = vhforg2(i,j) + sum2


        enddo
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine solar(istep)
!-----------------------------------------------------------------------
! Calculates incoming solar radiation as a function of latitude
! for each day of the year, given the orbital parameters (see PMIP)
! One year has 360 days.  Reference: A. Berger, JAS, 35, 2362-2367,1978
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'

      integer i,j,l,NVE


      real*8 beta,alam,alam0,ala,ala0
      real*8 fac1,fac2,fac3,ro,roref
      real*8 deltal, sindl, cosdl, tandl
      real*8 rkosz, rkosha1, ha1
      real*8 deg2rad, day2rad
      real*8 solard, solarcf(nlat)
      real*8 tsi,ksw
      real*8 alpho3sw(2)
      integer indxvol,istep,indxtsi

      deg2rad=pi/180.d0
      day2rad=pi/180.d0


! Present-day orbital parameters: eccentricity ecc, obliquity obl and
! angle om between Vernal Equinox and Perihelion (angles all given
! in degrees and converted to radians). Solarc is the solar constant.
! NVE is day of the Vernal Equinox, set at 21 MARCH
! Implementatie van Nanne

      if (iscencel .eq. 1) then
         if ((mod(nint(day * real(iatm)), nstpyear) .eq. 0) &
              & .or. (initialization .eqv. .true.)) &
              & call celest
         ecc = ecc2
         obl = asin(so)
         omweb = (perh + 180.00) * deg2rad
      elseif (iscencel .eq. 2) then
         if ((mod(nint(day * real(iatm)), nstpyear) .eq. 0) &
              & .or. (initialization .eqv. .true.)) &
              & call bretagnon
         ecc = ecc2
         obl = asin(so)
         omweb = (perh + 180.00) * deg2rad
      else
         ! iscencel == 0
!        ecc=0.016724
!        obl=23.446*deg2rad
!        omweb=(102.04+180.00)*deg2rad
         ecc = eccf
         obl = oblf * deg2rad
         omweb = (omwebf + 180.00) * deg2rad
!        write(iuo+99,*) 'Orbital parameters',eccf,oblf,omwebf
!        solarc = 1365.
      endif
      NVE=30+30+21
!
! In old SW-routine of ECBilt-model values were as follows:
!
!      ecc=0.0
!      solarc=1353.
!      NVE=90
!
! At 6000 years BP (Before Present) values were as follows:
!
!      ecc=0.018682
!      obl=24.105*deg2rad
!      omweb=(0.87+180.00)*deg2rad
!
! :0:
!     ecc=0.018994
!     obl=22.949*deg2rad
!     omweb=(114.42+180.00)*deg2rad


! First compute alam0 (the starting point). Then convert days to
! the true longitude ala, using Berger's formulae.
! Longitude ala loops from 0 (Vernal Equinox) to 359, ro is earth-sun
! distance relative to the major axis, del is declination.
!
      ala0=0.
      beta=(1.-ecc**2.)**0.5
      fac1=(0.5*ecc+0.125*ecc**3.)*(1.+beta)*sin(ala0-omweb)
      fac2=0.25*ecc**2.*(0.5+beta)*sin(2.*(ala0-omweb))
      fac3=0.125*ecc**3.*(1./3.+beta)*sin(3.*(ala0-omweb))
      alam0=ala0-2.*(fac1-fac2+fac3)


      l=(imonth-1)*30+iday-NVE
      if (l.lt.0) l=l+360
      alam=alam0+l*1.0*pi/180.


      fac1=(2.*ecc-0.25*ecc**3.)*sin(alam-omweb)
      fac2=1.25*ecc**2.*sin(2.*(alam-omweb))
      fac3=(13./12.)*ecc**3.*sin(3.*(alam-omweb))
      ala=alam+fac1+fac2+fac3
      ro=(1.-ecc**2.)/(1.+ecc*cos(ala-omweb))
      deltal=asin(sin(obl)*sin(ala))


      sindl=sin(deltal)
      cosdl=cos(deltal)
      tandl=tan(deltal)


! factor voor variable afstand Aarde-Zon (Berger, p.2367; Velds, p. 99)
      solard=1./ro**2.
!     solardref=1./roref**2.
      solardref=solard
      ksw=0.042
      alpho3sw(1)=0.247
      alpho3sw(2)=0.324
      solarcl(:)=solarm
      indxtsi=1
      indxvol=1
      if (iscentsi.eq.1) then
        indxtsi=irunlabelf+iyear-y1scentsi+1
        if ((initialization.eqv..true.).and.(irunlabeld.eq.360)) indxtsi=indxtsi+1
        if (indxtsi.gt.nyscenmaxtsi) indxtsi=nyscenmaxtsi
        if (indxtsi.lt.1) indxtsi=1
        solarc=solarm+tsiscen(2,indxtsi)*facttsi
        solarcl(:)=solarc
      endif

        if((mod(nint(day*real(iatm)),nstpyear).eq.0).or.(initialization.eqv..true.)) then
          write(iuo+99,*) 'tsi year=', tsiscen(1,indxtsi), "val=", tsiscen(2,indxtsi)
          write(iuo+99,11) ' sol forcing ',iyear,indxtsi,solarcl(15), tsiscen(2,indxtsi)*facttsi
        endif


      if (iscenvol.eq.1) then
!         write(*,*) "yo", nyears, ndays,  iyear, imonth, irunlabelf, irunlabel, irunlabeld
        indxvol=irunlabelf+iyear-y1scenvol+1
        if ((initialization.eqv..true.).and.(irunlabeld.eq.360)) indxvol=indxvol+1

        if (indxvol.gt.nyscenmaxvol) indxvol=nyscenmaxvol
        if (indxvol.lt.1) indxvol=1
        do i=1,10
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,1)
        enddo
        do i=11,16
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,2)
        enddo
        do i=17,22
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,3)
        enddo
        do i=23,32
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,4)
        enddo
      endif
      if (isceno3.eq.1) then
        do i=1,16
        solarcl(i)=solarcl(i)+(4*alpho3sw(1)*ksw*(o3-25.))
        enddo
        do i=17,32
        solarcl(i)=solarcl(i)+(4*alpho3sw(2)*ksw*(o3-25.))
        enddo
      endif
      !write(*,*) day, iatm, nstpyear, initialization
      if((mod(nint(day*real(iatm)),30*iatm).eq.0).or.(initialization.eqv..true.)) then
       write(iuo+99,12) 'vol forcing ',iyear,imonth,indxvol,solarcl(15), solarvol(indxvol,imonth,:)
      endif
11      format(A12,2i6,2f12.3)
12      format(A12,3i6,5f12.3)
      do i=1,nlat
! zonneconstante is 1370 in sw parameterisatie
       solarcf(i)=solarcl(i)/1370.d0
! beide effecten samen
       solarf(i)=solarcf(i)*solard
      enddo


      do i=1,nlat
         rkosha1=-tanfi(i)*tandl
         rkosha1=sign(min(abs(rkosha1),1.d0),rkosha1)
         ha1=acos(rkosha1)
         rkosz=sinfi(i)*sindl*(ha1 - tan(ha1))/pi
         if (ha1 .ne. 0.d0) then
            kosz(i) = rkosz*pi/ha1
         else
            kosz(i) = rkosz
         endif
         dayfr(i) = ha1/pi
         q0(i)=rkosz*solarcl(i)*solard
      enddo


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fluxes(nn)
!-----------------------------------------------------------------------
! *** computes energy fluxes above ocean surface
! *** short wave radiation, long wave radiation, sensible heat flux
! *** latent heat flux, evaporation
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comsurf.h'
      include 'comrunlabel.h'
      include 'comdiag.h'

      integer nn


      if (irad.eq.1) call swaverad2(nn)
      call swaverad(nn)
      call lwaverad(nn)
      if (irad.eq.1) call lwaverad2(nn)
      call dragcoef(nn)
      call surfmois(nn)
      call sensibheat(nn)
      call latentheat(nn)
      if (nn.eq.noc) call momentflux(nn)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine swaverad(nn)
!-----------------------------------------------------------------------
! *** computes short wave radiation
! *** linearization of RCM with ISCCP D2 1990 clouds
!-----------------------------------------------------------------------



      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'comdiag.h'


      integer i,j,k,l,ireg
      integer m, d, r, nn , nol


      real*8 f0,f1,ftot(8),fn(8,0:1)
      real*8 drs, drs2, drs3
      real*8 dcost, df,sk,sr,x,y,dfs,smsc
      real*8 fswdtoa,fswutoa,fswdsfc(nlat,nlon),fswusfc
      real*8 fswutoa2(nlat,nlon),fswutoa0(nlat,nlon)
      real*8 fswdtoa0(nlat,nlon)
      real*8 globalmean
      real*8 fswutoaG0,fswutoaGA,fswdtoa2,fswdtoaG0,fswdtoaGA
      real*8 fswutoa_diff,fswutoaG,df_test,fswdtoa_diff,fswdtoaG


      integer nreg(2),indxsul
      real*8 zac(2),asup
!     real*8 zac(2),asup,bup
      common /rad_sul0 /fswutoaG,df_test,fswdtoaG
      common /rad_sul1 /fswutoa0,fswutoaGA,fswdtoa0,fswdtoaGA

! *** aerosol scattering included as a correction on the upward
! *** clear sky fluxes
! *** sk,sr: empirical coefficients Dorland et al, J. Geophys. Res.,102,
! *** 28079-28100, 1997.
! *** smsc: mass scattering coefficient [m2/g]
! *** dso4: change in sulfate aerosol column integrated concentration since
! *** pre-industrial times [g/m2]


      sk=0.058d0*1370d0
      sr=0.05d0
      smsc=8.0
!     write(iuo+99,*) 'bup=',bup

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2

      indxsul=irunlabelf+iyear-issulstrt+1
      if ((initialization.eqv..true.).and.(irunlabeld.eq.360)) then
        indxsul=indxsul+1
        indxsul=12*(indxsul-1)+1
      else
        indxsul=12*(indxsul-1)+imonth
      endif
      if (indxsul.gt.nyscenmaxsul) indxsul=nyscenmaxsul


      do j=1,nlon
       do i=1,nlat
         if(indxsul.lt.1) then
          tas1(i,j) =0.
         else
           tas1(i,j) = sulopt(j,i,indxsul)
         endif

         nreg(nol)=irn(i,j,nol)
         zac(nol)=dble(costref(nreg(nol),imonth))
         if (zac(nol).GT.0.) then
           asup = (bup*tas1(i,j)*(1-albesn(i,j,nn))**2)/zac(nol)
         else
           asup =0.
         endif

         alb2esn(i,j,nn) = albesn(i,j,nn)+ asup
!        else
!          alb2esn(i,j,nn) = albesn(i,j,nn)
!        endif
          if (alb2esn(i,j,nn).ge.1.) then
!           write(*,*)alb2esn(i,j,nn),i,j,iyear,imonth,iday,zac(nol),
!    &       asup,albesn(i,j,nn),iscensul,tas1(i,j)
            alb2esn(i,j,nn)=1.
          endif
          df=dayfr(i)*solarf(i)
          ireg=irn(i,j,nol)
          dcost=kosz(i)-costref(ireg,imonth)
          do l=1,8
            do k=0,1
              fn(l,k) = swrref(l,ireg,imonth,k) &
     &               +  swrcost(l,ireg,imonth,k)*dcost
            enddo
          enddo

          x=sqrt(kosz(i))
          if(kosz(i).lt.0.) x=0d0
          y=sqrt(1-alb2esn(i,j,nn))
          dfs=sk*(4d0*x*y*(y-x)-sr)*dso4(i,j)*smsc
!         WRITE(*,*)dso4(i,j)
          if (dfs.gt.0d0.and.kosz(i).lt.0.05) dfs=0d0
          drs=alb2esn(i,j,nn)-salbref(ireg,imonth)
          drs2=drs*drs
          drs3=drs2*drs


          do l=1,4
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs+dfs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
!           write(*,*)f0,f1
!           if (l.eq.1)ftot_test=f0
          enddo
          do l=5,8
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo


! alternative calculation of upward flux at ground:
! in parameterisation no cross terms are accounted for, which are important for
! upward shortwave radiation at surface and therefore also for net flux
! heswsn(i,j)

          ftot(4)=-alb2esn(i,j,nn)*ftot(8)

          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*df
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*df
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*df
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*df

          if (irad.eq.1) then
! for diagnostic purposes:
! (1) downward shortwave radiation at TOA
             if (nn.eq.1)fswdtoa0(i,j)=0.
             fswdtoa2=-ftot(5)*df
             fswdtoa0(i,j)=fswdtoa0(i,j)+(fractn(i,j,nn)*fswdtoa2)
! (2) upward shortwave radiation at TOA
             if (nn.eq.1)fswutoa0(i,j)=0.
             fswutoa2(i,j)=ftot(1)*df
             fswutoa0(i,j)=fswutoa0(i,j)+(fractn(i,j,nn)*fswutoa2(i,j))
! (3) downward shortwave radiation at SURFACE
!            fswdsfc(i,j)=-ftot(8)*df
! (4) upward shortwave radiation at SURFACE
!            fswusfc=(heswsn(i,j)+ftot(8)*df)

         endif
        enddo
      enddo
      if((nn.eq.1).and.(oldmonth.ne.imonth)) then
        if(indxsul.ge.1) write(iuo+99, *) "Sulfate(25,2,",suloptTime(indxsul),")=", tas1(2,25)
        oldmonth=imonth
      endif

      if (irad.eq.1) then
       if (nn.eq.3) then
        fswutoaG0=globalmean(fswutoa0)
        fswdtoaG0=globalmean(fswdtoa0)
        fswutoa_diff=fswutoaG0-fswutoaG
        fswdtoa_diff=fswdtoaG0-fswdtoaG
        !if (iyear.eq.0) fswutoaGA=0.
        !if (iyear.eq.0) fswdtoaGA=0.
        if (initialization.eqv..true.) fswutoaGA=0.
        if (initialization.eqv..true.) fswdtoaGA=0.
        fswutoaGA=fswutoaGA+(fswutoa_diff/(360.*6.))
        fswdtoaGA=fswdtoaGA+(fswdtoa_diff/(360.*6.))
       endif
      else
       fswutoaGA=0.
       fswdtoaGA=0.
      endif



      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lwaverad(nn)
!-----------------------------------------------------------------------
! *** computes long wave radiation according to the parameterization of
! *** Chao Chou and Neelin and substantially adapted and extended
! *** for global scale and more
! *** specific ECBILT application by the one and only Michiel Schaeffer
! ***
! *** parameters: nlat   = number of gridpoints in the meridional
! ***                      direction (32)
! ***             nlon   = number of gridpoints in the zonal
! ***                      direction (64)
! ***
! *** input : dtemp(19,nlat,nlon): temperature anomalies [K] wrt ncep
! ***                              climatology tncep in common lwrscheme
! ***         dqa(nlat,nlon) : anomalies of total prec. water cont. below
! ***                          500 hPa wrt ncep climatology
! ***         tcc(nlat,nlon)  : total cloud cover
! ***         ghg(19) : concentrations of well mixed ghg's (see comphys.h)
! ***
! *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             toa
! ***          ulrad2(nlat,nlon): net longwave radiation [Wm-2] at
! ***                             (500 hPa)
! ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
! ***                             the surface
! ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             the surface
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'comdiag.h'


      integer i,j,l,k,m,is,ism,nol,nn,ireg,h
      real*8  lwr(7,0:1),dumts
      real*8  dqa,dqreg(27)
      real*8  ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm
      real*8  globalmean
      real*8  ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  ulrad2nz(nlat,nlon),ulradsnz(nlat,nlon)
      real*8  dlradsnz(nlat,nlon), ulrad0nUz(nlat,nlon)
      real*8  ulrad1nUz(nlat,nlon)
      common / radO3 / ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      common / radO32 / ulrad0nUz,ulrad1nUz



      is=imonth/3+1
      if (is.gt.4) is=1
      ism=(is-1)*3+1

      do i=1,27
!dqa    dqreg(i)=qancep(i,ism)**0.3333
        dqreg(i)=qancep(i,ism)
      enddo

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
        do i=1,nlat
          ireg=irn(i,j,nol)

!-hemispheric dependence of tropospheric ozone forcing
          if (i.le.16) then
           h=1
          else
           h=2
          endif

!dqa      dqa=lwrmois(i,j)-dqreg(ireg)
!dqa      q**1/3-qm**1/3=qm**(1/3-n)*(q**n-qm**n)
          dqa=dqreg(ireg)**(0.3333-EXPIR)* &
               & (lwrmois(i,j)**EXPIR-dqreg(ireg)**EXPIR)
!dqa      write(iuo+99,'(i3,3F12.5)') i,lwrmois(i,j)-
!dqa &      dqreg(ireg),dqa,lwrmois(i,j)**0.333-dqreg(ireg)**0.33333
          do l=0,1
            do k=1,7
              lwr(k,l)=lwrflux(k,ireg,is,l,h)+lwrqa(k,ireg,is,l)*dqa
!                   write(*,*) "lwr=",lwrflux(k,ireg,is,l,h),"+",lwrqa(k,ireg,is,l),"*",dqa
              do m=1,ipl(ireg)-1
                lwr(k,l)=lwr(k,l)+lwrt(k,m,ireg,is,l)*dtemp(m,i,j,nol)
              enddo
              lwr(k,l)=lwr(k,l)+lwrt(k,18,ireg,is,l)*dtemp(18,i,j,nol)
            enddo

            dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism)
            do m=1,4
              do k=1,3
                lwr(k,l)=lwr(k,l)+ &
                     & (lwrts(k,m,ireg,is,l)+lwrqts(k,m,ireg,is,l)*dqa)*dumts
              enddo
              lwr(7,l)=lwr(7,l)+ &
                   & (lwrts(7,m,ireg,is,l)+lwrqts(7,m,ireg,is,l)*dqa)*dumts
              dumts=dumts*(tsurfn(i,j,nn)-tncep(19,ireg,ism))
            enddo


          enddo

          ulrad0n(i,j,nn)=(lwr(1,0)*(1-tcc(i,j))+lwr(1,1)*tcc(i,j))
          ulrad1n(i,j,nn)=(lwr(2,0)+lwr(5,0))*(1-tcc(i,j)) + &
               & (lwr(2,1)+lwr(5,1))*tcc(i,j)
          ulrad2n(i,j,nn)=(lwr(3,0)+lwr(6,0))*(1-tcc(i,j)) + &
               & (lwr(3,1)+lwr(6,1))*tcc(i,j)

	  ulradsn(i,j,nn)=emisn(nn)*sboltz*tsurfn(i,j,nn)**4
          dlradsn(i,j,nn)=-lwr(7,0)*(1-tcc(i,j))-lwr(7,1)*tcc(i,j)


         if (irad.eq.1) then
          if(nn.eq.1) then
            ulrad0nUz(i,j)=0.
            ulrad1nUz(i,j)=0.
           endif

          ulrad0nUz(i,j)=ulrad0nUz(i,j)+(ulrad0n(i,j,nn)*fractn(i,j,nn))
          ulrad1nUz(i,j)=ulrad1nUz(i,j)+(ulrad1n(i,j,nn)*fractn(i,j,nn))

         endif
       enddo
      enddo

      if (irad.eq.1) then
       ulrad0nU=globalmean(ulrad0nUz)
       ulrad1nU=globalmean(ulrad1nUz)
      endif

!     ulrad0nm=globalmean(ulrad0n)
!     ulrad1nm=globalmean(ulrad1n)
!     ulrad2nm=globalmean(ulrad2n)
!     ulradsnm=globalmean(ulradsn)
!     dlradsnm=globalmean(dlradsn)

!     write(iuo+37,*)ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm

! *** that's all folks


      return
      end


!2345678901234567890123456789012345678901234567890123456789012345678901
      subroutine dragcoef(nn)
!----------------------------------------------------------------------
! *** drag coefficient
! *** depending on stability
!----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comsurf.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,nn


      real*8 cdrags,cdragl,tdif,cred,cdum

      cdrags=cdrag
      cdragl=cdrag
      cred=0.2
      cdum=(1-cred)*0.5
      do j=1,nlon
        do i=1,nlat
          cdragw(i,j)=cwdrag

          tdif=tsurfn(i,j,nn)-tempsgn(i,j,nn)
          cdragvn(i,j,nn)=cdrags*max(cred,cred+min(cdum+cdum*tdif,1-cred))

        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sensibheat(nn)
!-----------------------------------------------------------------------
! *** sensible heatflux between surface and atmosphere
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn


! *** sensible heatflux (watt/m**2)
      if (nn.eq.2) then
       do j=1,nlon
        do i=1,nlat
          hficof(i,j)=alphad*cdragvn(i,j,2)*uv10(i,j)
        enddo
       enddo
       do j=1,nlon
        do i=1,nlat
          hfluxn(i,j,2)=hficof(i,j)*(tsurfn(i,j,2)-tempsgn(i,j,2))
        enddo
       enddo
      else
       do j=1,nlon
        do i=1,nlat
           hfluxn(i,j,nn)=alphad*cdragvn(i,j,nn)*uv10(i,j)* &
     &                      (tsurfn(i,j,nn)-tempsgn(i,j,nn))
        enddo
       enddo
      endif


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine latentheat(nn)
!-----------------------------------------------------------------------
! *** latent heatflux between surface and atmosphere
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn,n,NSTAT
      real*8  qsat,db,emois,esubf,evapf,esnow,sfrac,edum,psilai

      real*8    bmoisg(nlat,nlon),resist(3),lai(2),k0(2)
      real*8    bmoism(nlat,nlon),rs
      real*8    rainf(nlat,nlon),snowf(nlat,nlon)
      real*8    fswdsfc(nlat,nlon)
      real*8    dc,eflux_t,eflux_g,eflux_bare
      real*8    fswdsfcM

      real*8 st,sg,sd,snlt,anup,blai,pnpp,b12,b34,b1,b2,b3,b4, &
     &       anup_moy,stock,st_moy
      real*8 stR,sgR,sdR,snltR

      common /lbmbmois/ bmoisg,bmoism,rainf,snowf
      common /pr_evap /fswdsfc
      COMMON /BIOTA/ &
     &   ST(nlat,nlon), SG(nlat,nlon), SD(nlat,nlon), SNLT(nlat,nlon), &
     &   BLAI(nlat,nlon,2), PNPP(nlat,nlon), &
     &   B12(nlat,nlon),   B34(nlat,nlon), &
     &   B1(nlat,nlon), B2(nlat,nlon), B3(nlat,nlon), B4(nlat,nlon), &
     &   ANUP_MOY(nlat,nlon),ANUP(nlat,nlon), STOCK(nlat,nlon), &
     &   st_moy(nlat,nlon), &
     &   NSTAT(nlat,nlon),STR(nlat,nlon), SGR(nlat,nlon),SDR(nlat,nlon), &
     &   SNLTR(nlat,nlon)


! *** latent heatflux due to evaporation from surface (watt/m**2)
! *** and evaporation rate (m/s)
! *** limit evaporation to available snow or bottom moisture
! *** evaporation factor =1 over snow, over wet land maximal 1


      if (nn.ne.nld) then
        do j=1,nlon
          do i=1,nlat
            evfacan(i,j,nn)=evfac
            efluxn(i,j,nn)=alphav*cdragvn(i,j,nn)*uv10(i,j)* &
     &                 (qsurfn(i,j,nn)-q10n(i,j,nn))
            efluxn(i,j,nn)=evfacan(i,j,nn)*max(0.d0,efluxn(i,j,nn))
            evapn(i,j,nn)=efluxn(i,j,nn)/(rowat*rlatvap)
          enddo
        enddo
      else
        do j=1,nlon
          do i=1,nlat



! *** evaporation factor =1 over snow, over wet land maximal 1
! *** care has to be taken in case of a snowcover in order to
! *** conserve heat: the evaporation is constraint to the amount
! *** of snow and moisture available; it can happen that in
! *** one timestep, the remaining snow is sublimated and part
! *** of the latent heat flux is also used to evaporate bottom
! *** moisture: Eflux = (1-sfrac)*Esub + sfrac*Evap


            if (fractn(i,j,nld).gt.epss) then


              if (adsnow(i,j).gt.0.) then
                evfacan(i,j,nld)=evfac
                edum=cdragvn(i,j,nld)*uv10(i,j)*(qsurfn(i,j,nld)-q10n(i,j,nld))
                edum=evfacan(i,j,nld)*max(edum,0d0)
                esubf=alphas*edum
                evapf=alphav*edum
                esnow=min(rowat*adsnow(i,j)*rlatsub/dtime,esubf)
                if (esnow.lt.esubf) then
                  sfrac=(esubf-esnow)/esubf
                  emois=min(rowat*abmoisg(i,j)*rlatvap/dtime,sfrac*evapf)
                  efluxn(i,j,nld)=esnow+emois
                  evapn(i,j,nld)=esnow/(rowat*rlatsub)+emois/(rowat*rlatvap)
                else
                  efluxn(i,j,nld)=esubf
                  evapn(i,j,nld)=esubf/(rowat*rlatsub)
                endif
              else
                evfacan(i,j,nld)=evfac*min(1d0,abmoisg(i,j)/max(1.0E-10,abmoism(i,j)))
                dc=qsat(pgroundn(i,j,nld),tempsgn(i,j,nld))-q10n(i,j,nld)
                psilai=0.2+(0.08*min(max(tempsgn(i,j,nld)-tzero,0.),10.))
                lai(1)=6*psilai
                lai(2)=2*psilai
                k0(1)=30.0E-5
                k0(2)=25.0E-5
                fswdsfcM=max(fswdsfc(i,j),1.0d0)
                do n=1,2
                rs=((fswdsfcM+125)/fswdsfcM)*(23E-3+(1.5*dc))*(1/(lai(n)*k0(n)))
                rs=max(0d0,rs)
                resist(n)=rs+(1/(cdragvn(i,j,nld)*uv10(i,j)))
                resist(n)=1/resist(n)
                enddo
                resist(3)=1/(cdragvn(i,j,nld)*uv10(i,j))
                resist(3)=1/resist(3)
!	write(*,*)cdragvn(i,j,nld)*uv10(i,j), rs(i,j)

!               efluxn(i,j,nld)=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld))*
!    &              ((st(i,j)*resist(1))+(sg(i,j)*resist(2))+
!    &              (sd(i,j)*resist(3)))
                eflux_bare=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld)) &
     &           *cdragvn(i,j,nld)*uv10(i,j)*(10./30.)
                eflux_g=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld)) &
     &           *resist(2)*(15./30.)
                eflux_t=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld)) &
     &           *resist(1)*(30./30.)
                efluxn(i,j,nld)=eflux_bare+(sg(i,j)*eflux_g) &
     &          + (st(i,j)*eflux_t)


                efluxn(i,j,nld)=evfacan(i,j,nld)* &
     &                          max(efluxn(i,j,nld),0d0)
                efluxn(i,j,nld)=min(abmoisg(i,j)*rowat*rlatvap/dtime &
     &                       ,efluxn(i,j,nld))
                evapn(i,j,nld)=efluxn(i,j,nld)/(rowat*rlatvap)
              endif


            endif
          enddo
        enddo
      endif

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine momentflux(nn)
!-----------------------------------------------------------------------
! *** computation of windstress
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn
      real*8  uv,costt,sintt,facstr

      facstr=roair * uv10rws

      do j=1,nlon
        do i=1,nlat
          costt=cos(dragane(i))
          sintt=sin(dragane(i))
          winstu(i,j)=cdragw(i,j)*facstr*uvw10(i,j)* &
               & (utot(i,j,3)*costt-vtot(i,j,3)*sintt)
          winstv(i,j)=cdragw(i,j)*facstr*uvw10(i,j)* &
               & (utot(i,j,3)*sintt+vtot(i,j,3)*costt)
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine totwind10
!-----------------------------------------------------------------------
! *** computation of strength of 10-meter winds (uv10r* 800 hPa wind)
! *** input  u800, v800 , udivg, vdivg
! *** output uv10 strength of 10 m wind at gaussian grid with minimum
! ***        uvw10 strength of 10 m wind at gaussian grid for windstress
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,k
      real*8  uv

! *** bug fix 27 march 97: uv was declared integer

      do j=1,nlon
        do i=1,nlat
          uv=sqrt((utot(i,j,3))**2 + (vtot(i,j,3))**2)
          uv10(i,j)=uv10rfx*uv
          uvw10(i,j)=uv10rws*uv
! *** minimum value of uv10 set to uv10m m/s
          if (uv10(i,j).lt.uv10m) uv10(i,j)=uv10m
        enddo
      enddo


      return
      end
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine tracer
!-----------------------------------------------------------------------
! *** advection of tracer field
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  hdivmg(nlat,nlon)
      real*8  co2sp(nsh2)



      call rggtosp(co2,co2sp)
      call sptogg(co2sp,co2,pp)



! *** horizontal divergence of tracer


      call trafluxdiv(hdivmg,co2sp,co2)


!
! *** time stepping forward time stepping
!


      do j=1,nlon
        do i=1,nlat
          co2(i,j)=co2(i,j)-dtime*(hdivmg(i,j))
          if (co2(i,j).lt.0d0) co2(i,j)=0d0
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisfields
!-----------------------------------------------------------------------
! *** calculates relative humidity of the moised layer
! *** and specific humidity above the surface and at the surface
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j,nn
      real*8  qsat,pmount,tmount,qmax,dqmdt

      do j=1,nlon
        do i=1,nlat


          call ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)


          if (qmax.gt.0d0) then


            relhum(i,j)=min(1d0,rmoisg(i,j)/qmax)


          else

            relhum(i,j)=0d0


          endif


          q10(i,j)= 0.d0
          do nn=1,ntyps
            q10n(i,j,nn)=relhum(i,j) * qsat(pgroundn(i,j,nn),tempsgn(i,j,nn))
          enddo

! *** lwrmois is used in the lwr parameterization

!dqa      lwrmois(i,j)=rmoisg(i,j)**0.3333
          lwrmois(i,j)=rmoisg(i,j)


        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine surfmois(nn)
!-----------------------------------------------------------------------
! *** calculates specific humidity at the surface
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comrunlabel.h'



      integer i,j,nn
      real*8  qsat


      do j=1,nlon
        do i=1,nlat


          qsurfn(i,j,nn)=qsat(pgroundn(i,j,nn),tsurfn(i,j,nn))

        enddo
      enddo


      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine cloud

!----------------------------------------------------------------------
! ***calculates the total cloud cover tcc using Axel's new scheme
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! *** calculates total cloud cover tcc, based on a global climatology
! *** from isccpd2. Michiel Schaeffer may 1998.
! *** diagnostic cloud scheme based on relative humidity and omega
! *** and stability. Jules Beersma
! *** depending on iradcloud, cloud climatology of diagnostic clouds are
! *** used in the calculation of radiative fluxes.
!-----------------------------------------------------------------------

      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comrunlabel.h'


      integer i,j
      real*8 rhc, rhfac
      real*8 cc


      do j=1,nlon
        do i=1,nlat

           rhc = relhcrit
           rhfac=relhfac
! enhance clouds in case of vertical motion
           if (omegg(i,j,2) .lt.  0.0 )  rhfac=0.95d0
           if (omegg(i,j,2) .lt. -0.04)  rhfac=0.90d0
! enhance clouds in areas of subsidence inversions
           if ((fractoc(i,j) .gt. epss) .and. (omegg(i,j,2) .gt.  0.03)) &
     &         rhfac=0.7d0

           cc=(relhum(i,j)/rhfac - rhc)/(1.0d0 - rhc)
           cc=max(cc,0.0d0)
           cc=min(cc,1.0d0)


           tccd(i,j)=cc
           if (iradcloud.eq.1) then
             tcc(i,j)=tccd(i,j)
           else
             tcc(i,j)=ccisccp(i,j,imonth)
           endif
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine detqmtabel
!-----------------------------------------------------------------------
!***  calculate tabel of maximum content of water in the layer below
!***  500 hPa using the Clausius-Clapeyron relation and assuming a
!***  temperature profile linear in log(p): T(p)= Tr + alpha*log(p/pr)
!***  where alpha is given by (T350-T650)/log(350/650)
!***  for given groundtemperature and 650 and 350 hPa temperatures
!***  the maximum water content in [m] is given by qmtabel
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k
      real*8  tmount,t500,b,qmax,expint,z1,z2,bz1,bz2,hulpx
      real*8  t350,t650,rlogp500,alpha
      real*8  detqmax,detqmaxexact
      real*4  hulp(0:iqmtab,0:jqmtab,0:kqmtab)


      rlogp500=log(500d0/650d0)
      b=cc2*cc3-cc2*tzero


!      call system('rm outputdata/atmos/qmtabel.dat')


!      open(88,file='qmtabel.dat')
!      open(89,file='qmtabel.test')
!     *  form='unformatted',access='direct',recl=51*21*21)


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


            qmax=(exp(bz1)+expint(1,-bz1)*bz1)/z1 - &
                 & (exp(bz2)+expint(1,-bz2)*bz2)/z2


            qmax=qmax*hulpx/alpha


            if (qmax.lt.0d0) qmax=0d0
            qmtabel(i,j,k)=qmax
!            write(88,111) tqmi(i),tqmj(j),tqmk(k),qmtabel(i,j,k)
          enddo
        enddo
      enddo


!      write(89) (((qmtabel(i,j,k),i=0,iqmtab),j=0,jqmtab),k=0,kqmtab)



!      do i=0,iqmtab
!        temp4g(1,1)=tqmi(i)
!        do j=0,jqmtab
!          tmount=tqmj(j)+temp4g(1,1)
!          do k=0,kqmtab
!            temp2g(1,1)=temp4g(1,1)-tqmk(k)
!            hulp(i,j,k)=detqmaxexact(tmount,1,1)
!            write(89,111) tqmi(i),tqmj(j),tqmk(k),hulp(i,j,k)
!          enddo
!        enddo
!      enddo



 111  format(4F14.7)
!      close(88)
!      close(89)
!      stop


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function detqmaxexact(tmount,i,j)
!-----------------------------------------------------------------------
! *** determines the maximum water content in latlon point
! *** i,j for given ground- and 650 and 350 hPa temperature
! *** by linear interpolation in qmtabel
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comunit.h'


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


      qmax=hulpx*(exp(bz1)+expint(1,-bz1)*bz1)/z1 - &
           & hulpx*(exp(bz2)+expint(1,-bz2)*bz2)/z2


      if (qmax.lt.0d0) qmax=0d0


      if (qmax.gt.0.2) then
        write(iuo+29,*) 'in latlon ',i,j,' qmax ',qmax
        call error(121)
      endif


      detqmaxexact=qmax


      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function detqmax(tmount,i,j,dqmdt)
!-----------------------------------------------------------------------
! *** determines the maximum water content in latlon point
! *** i,j for given ground- and 650 and 350 hPa temperature
! *** by linear interpolation in qmtabel
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comunit.h'


      integer i,j,ii,jj,kk
      real*8  tmount,ti,tj,tk
      real*8  qmax,detqmax,dqmdi,dqmdj,dqmdk
      real*8  dqmdt,hmount,z500,t500,alpha,dtgdt

      ti=temp4g(i,j)
      tj=tmount-temp4g(i,j)
      tk=temp4g(i,j)-temp2g(i,j)


      if (ti.lt.tqmi(0)) then
        ti=tqmi(0)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif
      if (ti.gt.tqmi(iqmtab)) then
        ti=tqmi(iqmtab)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif

      if (tj.lt.tqmj(0)) then
        tj=tqmj(0)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif
      if (tj.gt.tqmj(jqmtab)) then
        tj=tqmj(jqmtab)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif

      if (tk.lt.tqmk(0)) then
        tk=tqmk(0)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif
      if (tk.gt.tqmk(kqmtab)) then
        tk=tqmk(kqmtab)
!        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!        call error(121)
      endif

      ii=min(iqmtab-1,int((ti-tqmimin)*rdtqmi))
      jj=min(jqmtab-1,int((tj-tqmjmin)*rdtqmj))
      kk=min(kqmtab-1,int((tk-tqmkmin)*rdtqmk))
!       if( ii.lt.0) then
!         PRINT *,ii, jj, kk
!         PRINT *,"min(",iqmtab,"-1,int((",ti,"-",tqmimin,")*",rdtqmi,"))"
!         write(*,*) temp4g
!       endif

      dqmdi=(qmtabel(ii+1,jj,kk)-qmtabel(ii,jj,kk))*rdtqmi
      dqmdj=(qmtabel(ii,jj+1,kk)-qmtabel(ii,jj,kk))*rdtqmj
      dqmdk=(qmtabel(ii,jj,kk+1)-qmtabel(ii,jj,kk))*rdtqmk

      qmax = qmtabel(ii,jj,kk) + (ti-tqmi(ii))*dqmdi + &
           & (tj-tqmj(jj))*dqmdj + (tk-tqmk(kk))*dqmdk
      if (qmax.lt.0d0) qmax=0d0

      if (qmax.gt.0.2) then
        write(iuo+29,*) 'in latlon ',i,j,' qmax ',qmax
        call error(121)
      endif


      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12
      t500=temp4g(i,j)+alpha*alogpl2tl2
      z500=gpm500*grav
      hmount=qmount(i,j)*hmoisr*grav


      dtgdt=(rgas*t500*alogtl1pl2 + (hmount-geopg(i,j,2)-z500))/ &
           & (rgas*tmount*alogtl12)


      dqmdt=dqmdi + dqmdj * (dtgdt - 1d0) + dqmdk


      detqmax=0.9*qmax
!     detqmax=0.85*qmax


      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function expint(n,x)
      implicit none
      integer n,maxit
      real*8 expint,x,eps,fpmin,euler
      parameter (maxit=100,eps=1.e-10,fpmin=1.e-30,euler=.5772156649)
      integer i,ii,nm1
      real*8 a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
!        pause 'bad arguments in expint'
        call error(20)
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
!        pause 'continued fraction failed in expint'
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
!        pause 'series failed in expint'
        call error(20)
      endif
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisture
!-----------------------------------------------------------------------
! *** acvection and sources and sinks of atmospheric moisture
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  hdmoisg(nlat,nlon),d1(nlat,nlon),d3(nlat,nlon)
      real*8  d2(nlat,nlon),qstar,hdivmg(nlat,nlon),globalmean
      real*8  levtempgp,factemv,factems,omegg500,t500,qsat,gm1,gm2

      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair


      call rggtosp(rmoisg,rmoiss)
      call sptogg(rmoiss,rmoisg,pp)

      do j=1,nlon
        do i=1,nlat
          cormois(i,j)=0d0
          if (rmoisg(i,j).lt.0d0) then
            cormois(i,j)=cormois(i,j)-rmoisg(i,j)
            rmoisg(i,j)= 0d0
          endif
        enddo
      enddo


! *** horizontal divergence of moisture


      call trafluxdiv(hdivmg,rmoiss,rmoisg)

! *** vertical advection of moisture


      do j=1,nlon
        do i=1,nlat
          omegg500=(omegg(i,j,1)+omegg(i,j,2))/2.d0
!         d1(i,j)=omegg500
          vemoisg(i,j)=0d0
!         d2(i,j)=0d0
!          t500=levtempgp(plevel(2),i,j)
!          qstar=relhum(i,j)*qsat(plevel(2),t500)
!           d3(i,j)=qstar
!           d2(i,j)=t500
          if (omegg500.lt.0.d0) then
            t500=levtempgp(plevel(2),i,j)
!           d2(i,j)=t500
            qstar=relhum(i,j)*qsat(plevel(2),t500)
!           d3(i,j)=qstar
            vemoisg(i,j)=-omegg500*qstar/(grav*rowat)
            vemoisg(i,j)=min(vemoisg(i,j),rmoisg(i,j)/dtime)
            if (vemoisg(i,j).LT.0.d0) then
              vemoisg(i,j) = rmoisg(i,j)/dtime
            endif
            if (tsurf(i,j).ge.tzero) then
              dyrain(i,j) = dyrain(i,j) + vemoisg(i,j)
              thforg1(i,j)=thforg1(i,j) + factemv*vemoisg(i,j)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factemv*vemoisg(i,j)/dp1
            else
              dysnow(i,j) = dysnow(i,j) + vemoisg(i,j)
              thforg1(i,j)=thforg1(i,j) + factems*vemoisg(i,j)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factems*vemoisg(i,j)/dp1
            endif
          endif
        enddo
      enddo


! *** horizontal diffusion of moisture


      call hdiff(hdmoisg)


!
! *** time stepping forward time stepping
!


      do j=1,nlon
        do i=1,nlat
          rmoisg(i,j)=rmoisg(i,j)+dtime*(-ihavm*hdivmg(i,j) &
               & -ivavm*vemoisg(i,j) + hdmoisg(i,j) + imsink*evap(i,j))
          if (rmoisg(i,j).lt.0d0) then
            cormois(i,j)=cormois(i,j)-rmoisg(i,j)
            rmoisg(i,j)= 0d0
          endif
        enddo
      enddo

!      if (iyear.eq.7.and.imonth.eq.1) then
!        write(350) ((real(rmoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(hdivmg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(vemoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(hdmoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(cormois(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(d1(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(d2(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(d3(i,j)),j=1,nlon),i=1,nlat)
!      endif
      call moisbalance()
!      if (iyear.eq.7.and.imonth.eq.1) then
!        write(350) ((real(rmoisg(i,j)),j=1,nlon),i=1,nlat)
!        write(350) ((real(cormois(i,j)),j=1,nlon),i=1,nlat)
!      endif


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine trafluxdiv(tfdiv,ctrasp,ctra)
!-----------------------------------------------------------------------
! *** computes horizontal divergence of tracer flux
!-----------------------------------------------------------------------

      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k
      real*8  ctrasp(nsh2),vv(nsh2),ww(nsh2)
      real*8  dcdl(nlat,nlon),dcdm(nlat,nlon)
      real*8  tfdiv(nlat,nlon),ctra(nlat,nlon)


! *** 800 hPa winds are reduced with umoisr in the advection of the
! *** tracer field

! *** spatial derivatives of tracer


      call ddl (ctrasp,vv)
      call sptogg (vv,dcdl,pp)
      call sptogg (ctrasp,dcdm,pd)


! *** advection of tracer by total wind + convergence of tracer


      do j=1,nlon
        do i=1,nlat
          tfdiv(i,j)=dcdl(i,j)*(u800(i,j) + udivg(i,j,3))/ &
               & (radius*cosfi(i)) + dcdm(i,j)*(v800(i,j) + vdivg(i,j,3))/ &
               & (radius/cosfi(i)) + ctra(i,j)*divg(i,j,3)
          tfdiv(i,j)=tfdiv(i,j)*umoisr


        enddo
      enddo


      call rggtosp(tfdiv,vv)
      vv(1)=0d0
      call sptogg (vv,tfdiv,pp)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine hdivspec(hduvg,ug,vg)
!-----------------------------------------------------------------------
! *** computes horizontal divergence
!-----------------------------------------------------------------------

      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine hdiff(hdmg)
!-----------------------------------------------------------------------
! *** horizontal diffusion of moisture
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer idifq,k
      real*8  hdmoiss(nsh2),hdmg(nlat,nlon)
      real*8  difq,rll

      difq=max(0.d0,1.d0/(tdifq*24d0*3600d0))


      call lap(rmoiss,hdmoiss)


      hdmoiss(1)=0d0


      do k=2,nsh2
        hdmoiss(k)=difq*hdmoiss(k)
      enddo


      call sptogg(hdmoiss,hdmg,pp)


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine convec
!-----------------------------------------------------------------------
! *** moist convective adjustment
!-----------------------------------------------------------------------
      implicit none



      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'
      include 'comrunlabel.h'


      integer ncmax,iconvn,i,j
      real*8  qsatcr,tsatcr,pref,t500,qsat500,pot2g,pot4g,dcmoisg
      real*8  fachulp,facteta,factemv,factems,pmount,tmount
      real*8  qmax,qsat,hulp,redrain
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


! ***     calculate pressure and temperature at the ground
! ***     and the maximum water content


            call ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)


! ***     relhmax defines the relative humidity at which oversaturation
! ***     occurs


            qmax=relhmax*qmax


            pot2g=temp2g(i,j)/potfac1
            pot4g=temp4g(i,j)/potfac2
            teta(i,j)=0.5d0*(pot2g-pot4g)
! ***       dry adiabatic lapse rate
            tetacr(i,j)=0d0
            dcmoisg=0d0


            if (rmoisg(i,j).gt.qmax) then


! ***     calculate rain reduction factor to account for increased
! ***     moisture capacity due to latent heat release


              redrain=1d0+dqmdt*relhmax*rowat*rlatvap*grav/(cpair*dp2)


              dcmoisg=(rmoisg(i,j)-qmax)/(redrain*dtime)


! ***         if the air is supersaturated initially, this is due to
! ***         large scale convergence of moisture and large scale
! ***         temperature changes. Excessive moisture
! ***         is then removed as dynamic rain

              if (iconvn.eq.0) then
                if (tsurf(i,j).ge.tzero) then
                  dyrain(i,j)=dyrain(i,j) + dcmoisg
                  if (dyrain(i,j).gt.drainm) then
                    dcmoisg=drainm-dyrain(i,j)+dcmoisg
                    dyrain(i,j)=drainm
                  endif
                  thforg2(i,j)=thforg2(i,j) + factemv*dcmoisg/dp2
!                   write(*,*) "thforg2", thforg2(i,j),"=", factemv,"*",dcmoisg,"/",dp2
                  vhforg2(i,j)=vhforg2(i,j) + factemv*dcmoisg/dp2
                else
                  dysnow(i,j)=dysnow(i,j) + dcmoisg
                  if (dysnow(i,j).gt.drainm) then
                    dcmoisg=drainm-dysnow(i,j)+dcmoisg
                    dysnow(i,j)=drainm
                  endif
                  thforg2(i,j)=thforg2(i,j) + factems*dcmoisg/dp2
                  vhforg2(i,j)=vhforg2(i,j) + factems*dcmoisg/dp2
                endif
              else
                if (tsurf(i,j).ge.tzero) then
                  corain(i,j)=corain(i,j)+ dcmoisg
                  if (corain(i,j).gt.crainm) then
                    dcmoisg=crainm-corain(i,j)+dcmoisg
                    corain(i,j)=crainm
                  endif
                  temp4g(i,j) =temp4g(i,j) + factemv*dcmoisg*dtime/dp2
                  vhforg2(i,j)=vhforg2(i,j)+ factemv*dcmoisg/dp2
                else
                  cosnow(i,j)=cosnow(i,j)+ dcmoisg
                  if (cosnow(i,j).gt.crainm) then
                    dcmoisg=crainm-cosnow(i,j)+dcmoisg
                    cosnow(i,j)=crainm
                  endif
                  temp4g(i,j) =temp4g(i,j) + factems*dcmoisg*dtime/dp2
                  vhforg2(i,j)=vhforg2(i,j)+ factems*dcmoisg/dp2
                endif
              endif


! ***         moisture changes due to dynamic and convective rain


              rmoisg(i,j)=rmoisg(i,j)-ivavm*dcmoisg*dtime

              if (rmoisg(i,j).lt.0d0) then
!                write(99,*) 'moisture less than zero in convec !!!'
!                write(99,*) i,j,rmoisg(i,j),qmax,redrain
                rmoisg(i,j)= 0d0
              endif



! ***         calculate moist adiabatic lapse rate


              pot2g=temp2g(i,j)/potfac1
              pot4g=temp4g(i,j)/potfac2
              teta(i,j)=0.5d0*(pot2g-pot4g)


!              t500 = levtempgp(plevel(2),i,j)


              t500 = 0.5d0 * (temp2g(i,j) + temp4g(i,j))


              qsat500=qsat(plevel(2),t500)


              hulp=1d0 + fachulp*qsat500/(t500**2)
              gams(i,j)=gamd*(1+(rlatvap*qsat500)/(rgas*t500))/hulp
              tetacr(i,j)=0.5d0*t500*facteta*(gamd-gams(i,j))

            endif
            if (teta(i,j) .lt. (tetacr(i,j)-0.1)) then


              pot2g=(dp1*temp2g(i,j)+dp2*temp4g(i,j)+ &
     &                                 dp2*potfac2*2.*tetacr(i,j))/ &
     &              (potfac1*dp1+potfac2*dp2)
              pot4g=pot2g - 2.d0*tetacr(i,j)
              temp2go=temp2g(i,j)
              temp4go=temp4g(i,j)
              temp2g(i,j)=pot2g*potfac1
              temp4g(i,j)=pot4g*potfac2
!               write(*,*) "ici", temp4g(1,1),"=", pot4g,"*",potfac2
              vhforg1(i,j)=vhforg1(i,j) + (temp2g(i,j)-temp2go)/dtime
              vhforg2(i,j)=vhforg2(i,j) + (temp4g(i,j)-temp4go)/dtime
              iconvn=iconvn+1
              if (iconvn.lt.10) then
                if (dyrain(i,j).eq.drainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' dyrain'
                  call error(122)
                  goto 20
                endif
                if (corain(i,j).eq.crainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' corain'
                  call error(122)
                  goto 20
                endif
                if (dysnow(i,j).eq.drainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' dysnow'
                  call error(122)
                  goto 20
                endif
                if (cosnow(i,j).eq.crainm) then
                  write(iuo+29,*) 'in latlon ',i,j, ' cosnow'
                  call error(122)
                  goto 20
                endif
                goto 10
              else
                write(iuo+29,*) 'warning in lat-lon point: ',i,j
                write(iuo+29,*) temp2g(i,j),temp4g(i,j)
                call error(120)
                goto 20
              endif
            else
              goto 20
            endif
  20      continue
          if (iconvn.gt.ncmax) ncmax=iconvn
          torain(i,j)= dyrain(i,j) + corain(i,j)
          tosnow(i,j)= dysnow(i,j) + cosnow(i,j)
        enddo
      enddo
!      write(99,*) 'maximum iterations of convection : ',ncmax


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine moisbalance
!-----------------------------------------------------------------------
! *** fix moisture balance in the atmosphere: due to advection and
! *** spectral truncation, the moisture at a specific point can become
! *** negative. All moisture additions to reset the moisture at zero
! *** have been accumulated in cormois. In this routine these additions
! *** are subtracted from rmoisg. This is done for each latitude
! *** separately to prevent an artificial meridional transport of
! *** moisture. If at a given latitude not enough moisture is available
! *** to accomodate the subtraction, the moisture balance is fixed at
! *** the neighbouring equatorward latitude. The weights pw(nlat,1) are
! *** used to correct for changes in the gridsize with latitude.
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comrunlabel.h'


      integer i,j,nn
      real*8  gmc,gmm,gfac


      do i=1,nlat/2
        gmc=0d0
        gmm=0d0
        do j=1,nlon
          gmc=gmc+cormois(i,j)
          gmm=gmm+rmoisg(i,j)
        enddo
        if (gmm.gt.0d0.and.gmc.gt.0d0) then
          if (gmm.gt.gmc) then
            gfac=gmc/gmm
            do j=1,nlon
              rmoisg(i,j)=rmoisg(i,j)-gfac*rmoisg(i,j)
            enddo
          else
            gfac=gmm/gmc
            gmc=(gmc-gmm)*pw(i,1)/(pw(i+1,1)*dble(nlon))
            do j=1,nlon
              rmoisg(i,j)=0d0
              cormois(i,j)=gfac*cormois(i,j)
              cormois(i+1,j)=cormois(i+1,j)+gmc
            enddo
          endif
        endif
      enddo

      do i=nlat,1+nlat/2,-1
        gmc=0d0
        gmm=0d0
        do j=1,nlon
          gmc=gmc+cormois(i,j)
          gmm=gmm+rmoisg(i,j)
        enddo
        if (gmm.gt.0d0.and.gmc.gt.0d0) then
          if (gmm.gt.gmc) then
            gfac=1d0-gmc/gmm
            do j=1,nlon
              rmoisg(i,j)=rmoisg(i,j)*gfac
            enddo
          else
            gfac=gmm/gmc
            gmc=(gmc-gmm)*pw(i,1)/(pw(i-1,1)*dble(nlon))
            do j=1,nlon
              rmoisg(i,j)=0d0
              cormois(i,j)=gfac*cormois(i,j)
              cormois(i-1,j)=cormois(i-1,j)+gmc
            enddo
          endif
        endif
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function qsat(press,temp)
!-----------------------------------------------------------------------
! *** saturation mixing ratio
! *** input press in [Pa], temp in K
! *** output qsat: saturation mixing ratio
!-----------------------------------------------------------------------
      implicit none
      include 'comatm.h'
      include 'comphys.h'


      real*8  press,temp,qsat


      qsat=cc1*exp(cc2*(temp-tzero)/(temp-cc3))/press

      end



!23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine levtemp(tlev,plev)
!-----------------------------------------------------------------------
! *** computation temperatures at level p [Pa] assuming a constant
! *** temperature lapse rate : dt/dlnp = constant
! *** input: plev
! *** ouput: tlev
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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


!23456789012345678901234567890123456789012345678901234567890123456789012
      function levtempgp(plev,i,j)
!-----------------------------------------------------------------------
! *** computation temperatures at level p [Pa] assuming a constant
! *** temperature lapse rate : dt/dlnp = constant
! *** input: plev [Pa],i,j
! *** output: levtempgp [K]
!-----------------------------------------------------------------------
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


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fortemp
!-----------------------------------------------------------------------
! *** computes  new atmospheric temperatures due to heatforcing
! *** apply diffusion in stratosphere: timestep smaller than:
! *** 2*diffusiontimescale/462(=21*22) fe 10 days, 1 hour timestep
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comphys.h'
      include 'comdyn.h'
      include 'comrunlabel.h'


      integer i,j,k,it,ipd
      real*8  temp0sp(nsh2),tdifc,globalmean,tstep,tdifday

      tdifday=100d0
      tdifc=1.0d0/(tdifday*24.*3600.)

      tstep=2d0*tdifday*24.*3600./462.


      ipd=1+int(dtime)/int(tstep)


      tstep=dtime/ipd

      call ggtosp(temp0g,temp0sp)

      do it=1,ipd

        do k=1,nsh2
          temp0sp(k)=temp0sp(k) + tstep*rinhel(k,0)*temp0sp(k)*tdifc
        enddo
      enddo

      call sptogg(temp0sp,temp0g,pp)


      do j=1,nlon
        do i=1,nlat
          temp0g(i,j)=temp0g(i,j) + dtime*thforg0(i,j)
          temp2g(i,j)=temp2g(i,j) + dtime*thforg1(i,j)
          temp4g(i,j)=temp4g(i,j) + dtime*thforg2(i,j)
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine tempprofile
!-----------------------------------------------------------------------
! *** computation of vertical temperature profile
! *** based on reference profiles from NCEP reanalysis data
! *** also used in the LWR parameterisation
! *** assuming temperature anomalies wrt these temperature profiles
! *** vary linearly with log of the pressure
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comrunlabel.h'

      integer i,j,k,l,ireg(2),is,ism,nn,k1,k2,k2_tmp
      real*8 ro,ro1,ro2,z,z0,dt350,dt650(2),beta(2),tsref,dt100,z1,z2
      real*8 dtemp_tmp, pground_tmp,tsref_tmp,z_tmp,ro_tmp,z0_tmp
      real*8 beta_tmp

! *** Example reference profile for
! *** zonal band between 15s and 15n
! *** Month  4 : tncep(19,27,12)
! ***  nr.   pfl     tfl
! ***        (mb)    (K)
! ***  1   10.00  235.113
! ***  2   20.00  223.086
! ***  3   30.00  216.125
! ***  4   50.00  206.119
! ***  5   70.00  198.420
! ***  6  100.00  196.311
! ***  7  150.00  208.144
! ***  8  200.00  221.461
! ***  9  250.00  232.611
! *** 10  300.00  242.385
! *** 11  400.00  257.836
! *** 12  500.00  268.271
! *** 13  600.00  276.160
! *** 14  700.00  282.859
! *** 15  850.00  290.298
! *** 16  925.00  294.405
! *** 17 1000.00  299.345
! *** 18 1011.99  300.156
! *** Ps 1013.00  301.265 Ts


      is=imonth/3+1
      if (is.gt.4) is=1
      ism=(is-1)*3+1
      do j=1,nlon
        do i=1,nlat
          ireg(1)=irn(i,j,1)
          ireg(2)=irn(i,j,2)

! *** logarithmic interpolation of temperature anomalies; 200 hPa or
! *** higher the anomalies approach T100 temperature within 3
! *** pressure levels


          do nn=1,2
            dt100=temp0g(i,j)-tncep(6,ireg(nn),imonth)
            dt350=temp2g(i,j)-tncep12(1,ireg(nn),imonth)
            dt650(nn)=temp4g(i,j)-tncep12(2,ireg(nn),imonth)
            beta(nn)=(dt350-dt650(nn))*rlogtl12
            do k=1,6
              dtemp(k,i,j,nn)=dt100
            enddo
            do k=7,9
              dtemp(k,i,j,nn)=((10-k)*dt100+(k-6)*(dt650(nn) &
                   & + beta(nn)*rlogtl(k)))*0.25
            enddo
            do k=10,17
              dtemp(k,i,j,nn)=dt650(nn) + beta(nn)*rlogtl(k)
            enddo
          enddo
! *** from a mean height of 500 hPa from NCEP reanalysis data, the
! *** surface pressure is found using hydrostatic equilibrium
! *** and ideal gas law.


          z=z500ncep(ireg(1),imonth)
          k1=12
          DO WHILE (z.GT.rmountn(i,j,noc).AND.k1.LT.17)
            z0=z
            ro1=pncep(k1)/(rgas*(dtemp(k1,i,j,1)+tncep(k1,ireg(1),imonth)))
            ro2=pncep(k1+1)/(rgas*(dtemp(k1+1,i,j,1)+ &
                 & tncep(k1+1,ireg(1),imonth)))
            ro=(ro1+ro2)*0.5
            z=z-(pncep(k1+1)-pncep(k1))/(grav*ro)
            k1=k1+1
          END DO
          z1=z
          pgroundn(i,j,noc)=ro*grav*(z0-rmountn(i,j,noc))+pncep(k1-1)
          pgroundn(i,j,nse)=pgroundn(i,j,noc)


          z=z500ncep(ireg(2),imonth)
          k2=12

          DO WHILE ((z.GT.rmountn(i,j,nld)).AND.(k2.LT.17))
            z0=z
            ro1=pncep(k2)/(rgas*(dtemp(k2,i,j,2)+tncep(k2,ireg(2),imonth)))
            ro2=pncep(k2+1)/(rgas*(dtemp(k2+1,i,j,2)+ &
                 & tncep(k2+1,ireg(2),imonth)))
            ro=(ro1+ro2)*0.5
            z=z-(pncep(k2+1)-pncep(k2))/(grav*ro)
            k2=k2+1
          END DO

          z2=z
          pgroundn(i,j,nld)=ro*grav*(z0-rmountn(i,j,nld))+pncep(k2-1)
          pground(i,j)=0.0
          do nn=1,ntyps
            pground(i,j)=pground(i,j)+fractn(i,j,nn)*pgroundn(i,j,nn)
          enddo



! *** Temperature of air near surface is then found by
! *** interpolating for this pressure level the temperature profile calculated
! *** above.

          dtemp(18,i,j,1)=dt650(1)+beta(1)*log(pgroundn(i,j,noc)/65000.)
          beta(1)=(tncep(k1,ireg(1),imonth)-tncep((k1-1),ireg(1),imonth))/ &
               & log(pncep(k1)/pncep(k1-1))
          tsref=tncep((k1-1),ireg(1),imonth)+ &
               & beta(1)*log(pgroundn(i,j,noc)/pncep(k1-1))
          tempsgn(i,j,noc)=tsref+dtemp(18,i,j,1)
          tempsgn(i,j,nse)=tempsgn(i,j,noc)
          dtemp(18,i,j,2)=dt650(2)+beta(2)*log(pgroundn(i,j,nld)/65000.)
          beta(2)=(tncep(k2,ireg(2),imonth)-tncep((k2-1),ireg(2),imonth))/ &
               & log(pncep(k2)/pncep(k2-1))
          tsref=tncep((k2-1),ireg(2),imonth)+ &
               & beta(2)*log(pgroundn(i,j,nld)/pncep(k2-1))
          tempsgn(i,j,nld)=tsref+dtemp(18,i,j,2)
          tempsg(i,j)=0.0
          do nn=1,ntyps
            tempsg(i,j)=tempsg(i,j)+fractn(i,j,nn)*tempsgn(i,j,nn)
          enddo

! *** Temperature of pressure levels between diagnosed surface pressure and
! *** reference surface pressure are set equal to surface air temp.


          IF(z1.LE.rmountn(i,j,noc))THEN
            DO l=ipl(ireg(1)),k1,-1
               dtemp(l,i,j,1)=tempsgn(i,j,noc)-tncep(l,ireg(1),imonth)
            ENDDO
          ENDIF

          IF(z2.LE.rmountn(i,j,nld))THEN
            DO l=ipl(ireg(2)),k2,-1
               dtemp(l,i,j,2)=tempsgn(i,j,nld)-tncep(l,ireg(2),imonth)
            ENDDO
          ENDIF

! *** for LWR parameterisation temperature anomalies wrt seasonal mean are
! *** required:
          do nn=1,2
            DO k=1,18
              dtemp(k,i,j,nn)=tncep(k,ireg(nn),imonth) + &
                   & dtemp(k,i,j,nn)-tncep(k,ireg(nn),ism)
            ENDDO
          enddo


         enddo
       enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)
!-----------------------------------------------------------------------
! *** computation of ground pressure and temperature in order to
! *** to calculate the maximum precipitable water content in latlon i,j
! *** qmount contains the topography for this purpose
! *** assuming temperature varies linearly with log of the pressure
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comunit.h'
      include 'comrunlabel.h'


      integer i,j
      real*8  t500,levtempgp,hmount,hred,z500,dqmdt
      real*8  alpha,pfac,hfac,pmount,tmount,qmax,detqmax


      z500=gpm500*grav
      hfac=2/rgas
      hred=hmoisr*grav
      pfac=log(plevel(2)/tlevel(2))


! *** calculate temperature at t500 assuming the temperature varies
! *** linearly with log(p) : T = Tref + alpha * log (p/pref)


      alpha=(temp2g(i,j) - temp4g(i,j))*rlogtl12
      t500 =temp4g(i,j) + alpha*pfac


! *** calculate reduced ground height in decameters
! *** reduction occurs in order to tune the amount of moisture which
! *** is allowed to pass a topographic barier


      hmount=qmount(i,j)*hred
      if (hmount.lt.0d0) hmount=0d0


! *** calculate the groundpressure assuming that the mean geopotential
! *** height at 500 hPa is gpm500 decameter
! *** calculate 10 mtr temperature in K


      tmount=t500**2 - hfac*alpha*(hmount-geopg(i,j,2)-z500)
      if (tmount.lt.0) then
        write(iuo+29,*) 'in latlon ',i,j
        write(iuo+29,*) tmount,hmount,t500,geopg(i,j,2)
        call error(18)
      else
        tmount=sqrt(tmount)
      endif


!      pmount=plevel(2)*exp((tmount-t500)/alpha)

      qmax=detqmax(tmount,i,j,dqmdt)


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      function globalmean(gfield)
!-----------------------------------------------------------------------
! *** computes global mean of gfield
!-----------------------------------------------------------------------
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


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine meantemp
!-----------------------------------------------------------------------
! *** computes mean atmospheric temperatures
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      real*8  globalmean


! *** mean temperatures

      tempm(0)=globalmean(temp0g)
      tempm(1)=globalmean(temp2g)
      tempm(2)=globalmean(temp4g)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine dyntemp
!-----------------------------------------------------------------------
! *** computes temperature distribution in K from geopotential
! *** the mean level is given by tempm
! *** input:  geopg,tempm
! *** output: temp2g,temp4g
!-----------------------------------------------------------------------
      implicit none



      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,l
      real*8  tempfac(ntl)
      real*8  geogt(nlat,nlon,ntl),tempgt(nlat,nlon,ntl)


      tempfac(1)=350.d0/(rgas*300.d0)
      tempfac(2)=650.d0/(rgas*300.d0)


      do j=1,nlon
        do i=1,nlat
          geogt(i,j,1)=geopg(i,j,1)-geopg(i,j,2)
          geogt(i,j,2)=geopg(i,j,2)-geopg(i,j,3)
        enddo
      enddo

      do l=1,ntl


! ***  calculate temperatures and add mean temperature level


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


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine vortfor
!-----------------------------------------------------------------------
! *** computes vorticity forcing
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)
      real*8  forcgg(nlat,nlon),forcgg2(nlat,nlon)
      real*8  zetas(nsh2,3),zetag(nlat,nlon,3)
      real*8  dimfac


      call rggtosp(vhforg1,dfor1)
      call rggtosp(vhforg2,dfor2)


! ***  compute the relative vorticity zeta from steamfunction psi
! ***  psi is dimensionless, zeta with dimension


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


! *** potential vorticity (pv) forcing (1/(second**2)) due to
! *** the diabatic heating


      if (ipvf1.ne.0) call pvf1(vforg)


! *** pv forcing due to advection of the planetary vorticy f
! *** by the divergent wind

      if (ipvf2.ne.0) call pvf2(vforg)


! *** pv forcing due to d*zeta. d is the divergence.
! *** zeta is the relative vorticity


      if (ipvf3.ne.0) call pvf3(vforg,zetag)


! *** pv forcing due to advection of the relative vorticy
! *** zeta by the divergent wind


      if (ipvf4.ne.0) call pvf4(vforg,zetas)


! *** pv forcing due to advection of temperature by
! *** the divergent wind


      if (ipvf5.ne.0) call pvf5(vforg)


! *** the total pv forcing in nondimensional units


      dimfac=1./(om**2)


      call ggtosp(vforg(1,1,1),vfor1)
      call ggtosp(vforg(1,1,2),vfor2)
      call ggtosp(vforg(1,1,3),vfor3)


! *** adding the artificial forcing (forcgg1)

      vfor1(1)=0.d0
      vfor2(1)=0.d0
      vfor3(1)=0.d0


      do k=2,nsh2
        for(k,1)=dimfac*vfor1(k)
        for(k,2)=dimfac*vfor2(k)
        for(k,3)=dimfac*vfor3(k)
      enddo


! *** transfer to grid point


      call sptogg(vfor1,vforg1,pp)
      call sptogg(vfor2,vforg2,pp)
      call sptogg(vfor3,vforg3,pp)


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf1(vforg)
!-----------------------------------------------------------------------
! *** potential vorticity (pv) forcing (1/(second**2)) due to
! *** the adiabatic heating
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf2(vforg)
!-----------------------------------------------------------------------
! *** pv forcing due to advection of the planetary vorticy f
! *** by the divergent wind
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf3(vforg,zetag)
!-----------------------------------------------------------------------
! *** pv forcing due to d*zeta. d is the divergence.
! *** zeta is the relative vorticity
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf4(vforg,zetas)
!-----------------------------------------------------------------------
! *** pv forcing due to advection of the relative vorticy
! *** zeta by the divergent wind
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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
             vforg(i,j,l)=vforg(i,j,l) &
     &                  -udivg(i,j,l)*dxdl(i,j)/(radius*cosfi(i)) &
     &                  -vdivg(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine pvf5(vforg)
!-----------------------------------------------------------------------
! *** computes vorticity forcing due to advection of temperature by
! *** divergent wind
!-----------------------------------------------------------------------
      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comrunlabel.h'


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
             y(i,j,l)=ud(i,j,l)*dxdl(i,j)/(radius*cosfi(i)) &
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

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine celest
! ---
! --- Name and arguments:
! --- ^^^^^^^^^^^^^^^^^^^
! ---   celest
! ---
! --- Parents:
! --- ^^^^^^^^
! ---   insol2
! --- Children:
! --- ^^^^^^^^^
! ---   none
! --- Input:
! --- ^^^^^^
! ---   - time elapsed since 1950 (>0 after 1950) (annee, units : kyr)
! ---
! --- Output:
! --- ^^^^^^
! ---   - longitude of the perihelion measured from the moving vernal
! ---     point (perh, units : deg)
! ---   - eccentricity (ecc)
! ---   - sine of the obliquity (so)
! ---
! --- Read:
! --- ^^^^^
! ---   File: mbcs2_cor
! ---   ~~~~~
! ---     - amplitude for dvlpt of (e-pi) system (cfr. theoretical part)
! ---       (ae)
! ---     - mean rate for dvlpt of (e-pi) system (cfr. theoretical part)
! ---       (y -> be, units : rad.yr-1)
! ---     - phase for dvlpt of (e-pi) system (cfr. theoretical part)
! ---       (z -> ce, units : rad)
! ---
! ---     - amplitude for dvlpt of obliquity (aob, units : arcsecond)
! ---     - mean rate for dvlpt of obliquity
! ---       (y -> bob, units : rad.yr-1)
! ---     - phase for dvlpt of obliquity
! ---       (z -> cob, units : rad)
! ---
! ---     - amplitude for dvlpt of general precession in longitude
! ---       (aop, units : arcsecond)
! ---     - mean rate for dvlpt of general precession in longitude
! ---       (y -> bop, units : rad.yr-1)
! ---     - phase for dvlpt of general precession in longitude
! ---       (z -> cop, units : rad)
! ---
! --- Write:
! --- ^^^^^^
! ---   File: mout.d (iwr1)
! ---   ~~~~~
! ---     - number of terms for dvlpt of (e-pi) system
! ---       (cfr. theoretical part) (nef)
! ---     - number of terms for dvlpt of obliquity (nob)
! ---     - number of terms for dvlpt of general precession in
! ---       longitude (nop)
! ---     - climatic precession (pre)
! ---     - absolute value of year (or number of page-AB)
! ---       (ipage, units : kyr)
! ---     - year (ikyr, units : kyr)
! ---     - eccentricity (ecc)
! ---     - longitude of the perihelion measured from the moving vernal
! ---       point (perh, units : deg)
! ---     - obliquity (xob, units : deg)
! ---  ?  - caloric equator (xeq, units : deg)
! ---
! ---   File: bcmlog (iwr6)
! ---   ~~~~~
! ---     - year (ikyr, units : kyr)
! ---     - eccentricity (ecc)
! ---     - longitude of the perihelion measured from the moving vernal
! ---       point (perh, units : deg)
! ---     - obliquity (xob, units : deg)
! ---  ?  - caloric equator (xeq, units : deg)
! ---
! --- Computation:
! --- ^^^^^^^^^^^^
! ---   - pi (pi, units : rad)
! ---   - conversion factor deg --> rad (pir, units : rad.deg-1)
! ---   - conversion factor arcsecond --> rad (pirr, units : rad.arcsecond-1)
! ---   - independent term for dvlpt of obliquity (xod, units : deg)
! ---   - independent term for dvlpt of general precession in longitude
! ---     (xop, units : deg)
! ---   - linear term (precessional constant) for dvlpt of general
! ---     precession in longitude (prm, units : arcsecond.yr-1)
! ---   - eccentricity * sine (longitude of perihelion in a fixed
! ---     reference frame) (xes)
! ---   - eccentricity * cosine (longitude of perihelion in a fixed
! ---     reference frame) (xec)
! ---   - argument of the trigonometrical function in the expansions of
! ---     (e-pi) system or general precession in longitude or obliquity
! ---     (arg, units : rad)
! ---   - longitude of perihelion in a fixed reference frame
! ---     (rp, units : rad; drp, units : deg)
! ---   - general precession in longitude (prg, units : arcsecond;
! ---     dprg, units : deg)
! ---
! --- Bibliography:
! --- ^^^^^^^^^^^^^
! ---   - Berger, A., Long-term variations of daily insolation and
! ---     Quaternary climatic changes
! ---     J. Atmos. Sc., 35 (12), 2362-2367, 1978.
! ---
! --- Last changes made (by,when):
! --- ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! ---   DUTRIEUX Alexis and LOUTRE Marie-France, 05/11/1992
! ---
! --- ===================================================================
!
       implicit none
!      implicit double precision (a-h,o-z)
!
!      include 'div08.inc'
!      include 'ray15.inc'
!      include 'tst01.inc'
      integer  iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
      common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
!      double precision perh,ecc2,so
!      common /so10/ perh,ecc2,so
      double precision errtst(-16:+16)
      common/tst/errtst

      include 'comatm.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
!
      real*8    ae(19),be(19),ce(19),aob(104),bob(104),cob(104), &
           & aop(177),bop(177),cop(177)
      real*8 pi1,pir,pirr,annee,y,z,xod,xop,prm,t,xes,xec,arg, &
           & rp,drp,prg,dprg,pre,xob,xeq
      integer nef,nob,nop,ipage,ikyr,i
! --- -------------------------------------------------------------------
!

!*****change by H. Renssen, 25-02-03
!     annee=(19550.0d0-irunlabelf-iyear)/-1000.0d0
      annee=(+irunlabelf+iyear-1950)/1000.0d0
!*****
! --- 0. Computation of errtst (= needed to resolve error test)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=-16,16,1
         errtst(i)=10.d0**(i)
      enddo
! ----

      pi1 = dacos(-1.0d0)
      pir = pi1 / 180.0d0
      pirr = pir / 3600.0d0
!
! --- 1.earth orbital elements :
! --- ^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!     open(unit=47,status='old',
!    +file='/home/astr/hgs/ecclio_vecode/ecbilt/inputdata/mbcs2_cor')
!     rewind  47
      rewind (iuo+38)
!
! --- Eccentricity:
! --- ~~~~~~~~~~~~~
!
      nef = 19
      do 10 i=1,nef
!       read(47,9000) ae(i),y,z
        read(iuo+38,9000) ae(i),y,z
        be(i) = y * pirr
        ce(i) = z * pir
 10   continue
!
! --- Obliquity:
! --- ~~~~~~~~~~
!
      xod = 23.320556d0
      nob = 104
      do 20 i=1,nob
        read(iuo+38,9005) aob(i),y,z
!       read(47,9005) aob(i),y,z
        bob(i) = y * pirr
        cob(i) = z * pir
 20   continue
!
! --- General precession in longitude:
! --- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      xop = 3.392506d0
      prm = 50.439273d0
      nop = 177
      do 30 i=1,177
!       read(47,9005) aop(i),y,z
        read(iuo+38,9005) aop(i),y,z
        bop(i) = y * pirr
        cop(i) = z * pir
 30   continue
!
!     close(unit=47)
!
!
! --- 2.numerical value for ecc pre xob:
! --- ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      t = annee * 1000.0d0
      xes = 0.0d0
      xec = 0.0d0
      do 40 i=1,nef
        arg = be(i) * t + ce(i)
        xes = xes + ae(i) * dsin(arg)
        xec = xec + ae(i) * dcos(arg)
 40   continue
      ecc2 = dsqrt(xes * xes + xec * xec)
      if ((dabs(xec).lt.errtst(-8)).and.(dabs(xes).lt.errtst(-8))) then
        rp = 0.0d0
      else
        rp = datan2(xes,xec)
        if (rp.lt.0.0d0) then
          rp = rp + 2.0d0 * pi1
        endif
      endif
      drp = rp / pir
!
      prg = prm * t
      do 50 i=1,nop
        arg = bop(i) * t + cop(i)
        prg = prg + aop(i) * dsin(arg)
 50   continue
      dprg = prg / 3600.0d0 + xop
      perh = drp + dprg
      perh = dmod(perh,360.0d0)
      if (perh.lt.0.0d0) then
        perh = perh + 360.0d0
      endif
!
      pre = ecc2 * dsin(perh * pir)
!
      xob = xod
      do 60 i=1,nob
        arg = bob(i) * t + cob(i)
        xob = xob + aob(i) / 3600.0d0 * dcos(arg)
 60   continue
!
      so = dsin(xob * pir)
      xeq = (datan(4.0d0 * pre / (pi1 * so))) / pir
!
      ipage = dabs(t / 1000.)
!
      ikyr = t / 1000.
!
!     *** change by Hans Renssen, June 13 2001
      write(iuo+99,*) 'celest: annee, irunlabelf, iyear, ecc2(eccf), so(oblf), perh(omwebf)'
      write(iuo+99,*) annee,irunlabelf,iyear,ecc2,so,perh
        call flush(iuo+99)
!     *** end change HR
!
 9000 format(13x,f11.8,f20.7,f20.6)
 9005 format(7x,f13.7,2x,f10.6,2x,f10.4)
 9010 format(//,1x,'long term daily insolation',/,1x, &
           & 'number of terms in:',1x,'eccentricity',i5,2x, &
           & 'obliquity',i5,2x,'general precession',i5,//)
 9015 format('precession climatique= ',f8.5,5x,'page = ',i4)
 9020 format(1x,'date = ',i6,3x,'eccen = ',f9.6,3x,'long per = ', &
           & f7.2,3x,'obliq = ',f7.3,3x,'cal eq = ',f5.2,/)
!
      return
      end
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lwaverad2(nn)
!-----------------------------------------------------------------------
! *** computes long wave radiation according to the parameterization of
! *** Chao Chou and Neelin and substantially adapted and extended
! *** for global scale and more
! *** specific ECBILT application by the one and only Michiel Schaeffer
! ***
! *** parameters: nlat   = number of gridpoints in the meridional
! ***                      direction (32)
! ***             nlon   = number of gridpoints in the zonal
! ***                      direction (64)
! ***
! *** input : dtemp(19,nlat,nlon): temperature anomalies [K] wrt ncep
! ***                              climatology tncep in common lwrscheme
! ***         dqa(nlat,nlon) : anomalies of total prec. water cont. below
! ***                          500 hPa wrt ncep climatology
! ***         tcc(nlat,nlon)  : total cloud cover
! ***         ghg(19) : concentrations of well mixed ghg's (see comphys.h)
! ***
! *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             toa
! ***          ulrad2(nlat,nlon): net longwave radiation [Wm-2] at
! ***                             (500 hPa)
! ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
! ***                             the surface
! ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             the surface
!-----------------------------------------------------------------------
           implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'


      integer i,j,l,k,m,is,ism,nol,nn,ireg,h,r,s,igas
      real*8  lwrz(7,0:1),dumts
      real*8  dqa,dqreg(27)
      real*8  ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm
      real*8  ulrad0nmm,ulrad1nmm
      real*8  ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  ulrad1nzz(nlat,nlon,3)
      real*8  ulrad2nz(nlat,nlon),ulradsnz(nlat,nlon)
      real*8  dlradsnz(nlat,nlon)
      real*8  ulrad0nT,ulrad1nT,ulrad2nT,ulradsnT,dlradsnT
      real*8  globalmean
      real*8  logco2T,sqrch4T,sqrn2oT,ghgz(20)
      real*8  alpho3lw(2)
      real*4  lwrfluxz(7,27,4,0:1,2)
      real*8  moc,tmc,tmc0,tsurfmean,cland,thex
      common / radO3 / ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      common /rad031/ulrad0nz,ulrad1nz,ulrad0nT,ulrad1nT
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex



        ghgz(1)=280.
        do igas=2,20
         ghgz(igas)=ghgscen(igas,1)
        enddo
        logco2T=log(ghgz(1)/ghgipcc(1))
        sqrch4T=sqrt(ghgz(2))-sqrt(ghgipcc(2))
        sqrn2oT=sqrt(ghgz(3))-sqrt(ghgipcc(3))
        alpho3lw(1)=153.6
        alpho3lw(2)=201.2
        do h=1,2
        do l=0,1
         do s=1,4
          do r=1,27
           do k=1,7
            lwrfluxz(k,r,s,l,h)=lwrref(k,r,s,l)+lwrghg(k,1,r,s,l)*logco2T+ &
                 & lwrghg(k,2,r,s,l)*sqrch4T+lwrghg(k,3,r,s,l)*sqrn2oT
            do m=4,19
             lwrfluxz(k,r,s,l,h)=lwrfluxz(k,r,s,l,h)+ &
                  & lwrghg(k,m,r,s,l)*(ghgz(m)-ghgipcc(m))
            enddo
              lwrfluxz(k,r,s,l,h)=lwrfluxz(k,r,s,l,h)+ &
                   & lwrghg(k,4,r,s,l)*alpho3lw(h)*(ghgz(20)-25.)
           enddo
          enddo
         enddo
        enddo
        enddo

       is=imonth/3+1
       if (is.gt.4) is=1
       ism=(is-1)*3+1

      do i=1,27
!dqa    dqreg(i)=qancep(i,ism)**0.3333
        dqreg(i)=qancep(i,ism)
      enddo

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
        do i=1,nlat
          ireg=irn(i,j,nol)

!-Hemispheric dependence of tropospheric ozone forcing
          if (i.le.16) then
           h=1
          else
           h=2
          endif

!dqa      dqa=lwrmois(i,j)-dqreg(ireg)
!dqa      q**1/3-qm**1/3=qm**(1/3-n)*(q**n-qm**n)
          dqa=dqreg(ireg)**(0.3333-EXPIR)* &
               & (lwrmois(i,j)**EXPIR-dqreg(ireg)**EXPIR)

          do l=0,1
            do k=1,7
              lwrz(k,l)=lwrfluxz(k,ireg,is,l,h)+lwrqa(k,ireg,is,l)*dqa
              do m=1,ipl(ireg)-1
                lwrz(k,l)=lwrz(k,l)+lwrt(k,m,ireg,is,l)*dtemp(m,i,j,nol)
              enddo
              lwrz(k,l)=lwrz(k,l)+lwrt(k,18,ireg,is,l)*dtemp(18,i,j,nol)
            enddo

            dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism)
            do m=1,4
              do k=1,3
                lwrz(k,l)=lwrz(k,l)+ &
                     & (lwrts(k,m,ireg,is,l)+lwrqts(k,m,ireg,is,l)*dqa)*dumts
              enddo
!             lwrz(7,l)=lwrz(7,l)+
!    *        (lwrts(7,m,ireg,is,l)+lwrqts(7,m,ireg,is,l)*dqa)
!    *        *dumts
              dumts=dumts*(tsurfn(i,j,nn)-tncep(19,ireg,ism))
            enddo


          enddo

          if (nn.eq.1) then
            ulrad1nz(i,j)=0.
            if(initialization.eqv..true.) then
            !if (iyear.eq.0) then
             ulrad1nT=0.
            endif
          endif


          ulrad1nzz(i,j,nn)=(lwrz(2,0)+lwrz(5,0))*(1-tcc(i,j)) + &
               & (lwrz(2,1)+lwrz(5,1))*tcc(i,j)
         ulrad1nz(i,j)=ulrad1nz(i,j)+(ulrad1nzz(i,j,nn)*fractn(i,j,nn))

        enddo
      enddo


      if (nn.eq.3) then
      ulrad1nm=globalmean(ulrad1nz)

      ulrad1nmm=ulrad1nm-ulrad1nU

      ulrad1nT=ulrad1nT+(ulrad1nmm/(360.*6.))

      endif

! *** that's all folks
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine swaverad2(nn)
!-----------------------------------------------------------------------
! *** computes short wave radiation
! *** linearization of RCM with ISCCP D2 1990 clouds
!-----------------------------------------------------------------------



      implicit none


      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comemic.h'
      include 'comsurf.h'
      include 'comunit.h'

      integer i,j,k,l,ireg
      integer m, d, r, nn , nol


      real*8 f0,f1,ftot(8),fn(8,0:1)
      real*8 drs, drs2, drs3
      real*8 dcost, df,sk,sr,x,y,dfs,smsc,df2
      real*8 fswutoa(nlat,nlon),fswdsfc(nlat,nlon),fswusfc
      real*8 fswutoa2(nlat,nlon),fswdtoa(nlat,nlon)
      real*8 fswutoaG,fswdtoa2,fswdtoaG


      integer nreg(2),indxsul
!     real*8 zac(2),asup,bup
      real*8 zac(2),asup
      real*8  globalmean
      real*8 fswutoaGA,fswutoaG0
      real*8 fswutoa_diff,df_test
      common /rad_sul2 /fswutoa,fswdtoa
      common /rad_sul0 /fswutoaG,df_test,fswdtoaG
      common /pr_evap /fswdsfc
! *** aerosol scattering included as a correction on the upward
! *** clear sky fluxes
! *** sk,sr: empirical coefficients Dorland et al, J. Geophys. Res.,102,
! *** 28079-28100, 1997.
! *** smsc: mass scattering coefficient [m2/g]
! *** dso4: change in sulfate aerosol column integrated concentration since
! *** pre-industrial times [g/m2]


      sk=0.058d0*1370d0
      sr=0.05d0
      smsc=8.0

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
       do i=1,nlat
           alb2esn(i,j,nn)=albesn(i,j,nn)
           alb2esn(i,j,3) = albesnR(i,j)
          if (alb2esn(i,j,nn).ge.1.) then
            alb2esn(i,j,nn)=1.
          endif
          df=dayfr(i)*solarf(i)
          df2=dayfr(i)*solarm*solardref/1370.d0
          ireg=irn(i,j,nol)
          dcost=kosz(i)-costref(ireg,imonth)
          do l=1,8
            do k=0,1
              fn(l,k) = swrref(l,ireg,imonth,k) &
     &               +  swrcost(l,ireg,imonth,k)*dcost
            enddo
          enddo

          x=sqrt(kosz(i))
          y=sqrt(1-alb2esn(i,j,nn))
          dfs=sk*(4d0*x*y*(y-x)-sr)*dso4(i,j)*smsc
          if (dfs.gt.0d0.and.kosz(i).lt.0.05) dfs=0d0
          drs=alb2esn(i,j,nn)-salbref(ireg,imonth)
          drs2=drs*drs
          drs3=drs2*drs


          do l=1,4
           f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs+dfs
           f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo
          do l=5,8
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs &
     &                     +swrsalb(l,ireg,imonth,2)*drs2 &
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo


! alternative calculation of upward flux at ground:
! in parameterisation no cross terms are accounted for, which are important for
! upward shortwave radiation at surface and therefore also for net flux
! heswsn(i,j)

          ftot(4)=-alb2esn(i,j,nn)*ftot(8)

          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*df
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*df
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*df
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*df

! for diagnostic purposes:
! (1) downward shortwave radiation at TOA
             if (nn.eq.1)fswdtoa(i,j)=0.
             fswdtoa2=-ftot(5)*df2
             fswdtoa(i,j)=fswdtoa(i,j)+(fractn(i,j,nn)*fswdtoa2)
! (2) upward shortwave radiation at TOA
            if (nn.eq.1)fswutoa(i,j)=0.
            fswutoa2(i,j)=ftot(1)*df2
            fswutoa(i,j)=fswutoa(i,j)+(fractn(i,j,nn)*fswutoa2(i,j))
! (3) downward shortwave radiation at SURFACE
            fswdsfc(i,j)=-ftot(8)*df
! (4) upward shortwave radiation at SURFACE
!            fswusfc=(heswsn(i,j)+ftot(8)*df)

        enddo
      enddo
      if (nn.eq.3) then
        fswutoaG=globalmean(fswutoa)
        fswdtoaG=globalmean(fswdtoa)
      endif


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine bretagnon
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension dkhqp(4),psiom(2)
      dimension epi(2)
!
      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'

2000  format(/)

      pi314=3.1415926535897932d0
      datzer=2451545.d0-500*365.25d0
      pas=0.25d0
      datzer=2451545.d0-2000*365.25d0
      pas=1.d0
        i=(+irunlabelf+iyear)
        dj=datzer+365.25d0*pas*i
        date=(dj-2451545.d0)/365.25d0+2000
        call subkhqp(dj,dkhqp)
        call subpsiom(dj,psiom)
        call obliquity(dj,dkhqp,psiom,obliq)
        call expi(dj,dkhqp,epi)
!       write(iuo+66,4000) dj,date,epi,obliq
4000    format(f10.2,f10.2,f11.8,f10.6,f11.8)
        ecc2=epi(1)
        perh=epi(2)*180.0D0/pi314
        obl=obliq*180.0D0/pi314
        so=dsin(obliq)

      return
      end

      SUBROUTINE obliquity(dj,dkhqp,psiom,obliq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension dkhqp(4),psiom(2)

!alcul de l'inclinaison et du noeud de l'ecliptique puis du vecteur unite
!  normal a l'ecliptique vrai de la date

      cosisur2=sqrt(1.d0-dkhqp(3)**2-dkhqp(4)**2)
      xi= 2*dkhqp(4)*cosisur2
      yi=-2*dkhqp(3)*cosisur2
      zi=sqrt(1.d0-2*dkhqp(3)**2-2*dkhqp(4)**2)
      di=2*dasin(sqrt(dkhqp(3)**2+dkhqp(4)**2))
      do=datan2(dkhqp(4),dkhqp(3))

!alcul du vecteur unite normal a l'equateur vrai de la date

      xj=sin(psiom(1))*sin(psiom(2))
      yj=-cos(psiom(1))*sin(psiom(2))
      zj=cos(psiom(2))

!alcul du produit scalaire et de l'obliquite

      pscal=xi*xj+yi*yj+zi*zj
      obliq=dacos(pscal)

      return
      end

      SUBROUTINE EXPI(dj,dkhqp,epi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension dkhqp(4),epi(2)
      data dj2000/2451545.d0/

      pi314=3.1415926535897932d0

      t=(dj-dj2000)/365250.d0

!alcul de e, pi, i, Omega rapportes a l'ecliptique J2000

      e=sqrt(dkhqp(1)**2+dkhqp(2)**2)
      pi0=atan2(dkhqp(2),dkhqp(1))
      di=2*dasin(sqrt(dkhqp(3)**2+dkhqp(4)**2))
      omega0=atan2(dkhqp(4),dkhqp(3))
!      write(iuo+66,3000) e,pi0,di,omega0
3000  format('e =',f11.8,'  pi0 =',f10.6,'  di =',f11.8,'  Om0 =',f10.6)

!alcul de pi_A, Pi_A, P_A

      dq=-0.00113472340d0*t+0.00001236850d0*t**2+0.00000126542d0*t**3
      dp=0.00010179800d0*t+0.00004702790d0*t**2-0.00000054230d0*t**3
      pa=0.2438016495d0*t+539.3230d-6*t**2+373.3d-9*t**3-1138.3d-9*t**4 &
           & -8.6d-9*t**5
      ppia=2*dasin(sqrt(dq**2+dp**2))
      gpia=atan2(dp,dq)

4001  format('dq,dp,pa,ppia,gpia',5d15.8)

      z=cos(di)*cos(ppia)+sin(di)*sin(ppia)*cos(gpia-omega0)
      y1=sin(di)*sin(gpia-omega0)
      x1=-cos(di)*sin(ppia)+sin(di)*cos(ppia)*cos(gpia-omega0)
      y2=sin(ppia)*sin(gpia-omega0)
      x2=sin(di)*cos(ppia)-cos(di)*sin(ppia)*cos(gpia-omega0)
!      write(iuo+66,4002) z,y1,x1,y2,x2
4002  format('z,y1,x1,y2,x2',5d15.8)

!      did=dacos(z)
      omegad=pa+gpia-atan2(y1,x1)
      deltaomega=atan2(y2,x2)
      pid=pi0+omegad-omega0+deltaomega
      pid=mod(pid,2*pi314)
!      write(iuo+66,3001) pid,did,omegad
3001  format('pid =',f10.6,' did =',f11.8,'  Omd =',f10.6)

      epi(1)=e
      epi(2)=pid

      return
      end

      SUBROUTINE SUBKHQP(dj,dkhqp)
      IMPLICIT NONE
      DOUBLE PRECISION t, dj, dj2000
      DOUBLE PRECISION dkhqp(4)
      DOUBLE PRECISION dk0(3,102),dk1(3,10),dk2(3),dk3(3), &
           & dh0(3,102),dh1(3,10),dh2(3),dh3(3),&
           & dq0(3,16),dq1(3,2),dq2(3),dq3(3), &
           & dp0(3,16),dp1(3,2),dp2(3),dp3(3)
      INTEGER dimkh(4),dimqp(4)
      INTEGER j
      data dimkh/102,10,1,1/
      data dimqp/16,2,1,1/
      data dj2000/2451545.d0/
      data dk0/ &
     &   0.00374081650d0, 3.14159265359d0,       0.00000000000d0, &
     &   0.00001988852d0, 4.23374621009d0,    1577.34354244780d0, &
     &   0.00001859231d0, 0.55463591479d0,    5223.69391980220d0, &
     &   0.00001497439d0, 3.72409379834d0,     529.69096509460d0, &
     &   0.00000823038d0, 0.33112005725d0,    2352.86615377180d0, &
     &   0.00000483421d0, 3.17751155482d0,   10213.28554621100d0, &
     &   0.00000483174d0, 5.65660621901d0,    5507.55323866740d0, &
     &   0.00000441134d0, 1.21138303352d0,     398.14900340820d0, &
     &   0.00000354179d0, 0.13054037265d0,    4694.00295470760d0, &
     &   0.00000278456d0, 5.27885763953d0,    1059.38193018920d0, &
     &   0.00000294442d0, 3.92692187453d0,     775.52261132400d0, &
     &   0.00000229665d0, 0.79626727662d0,    9437.76293488700d0, &
     &   0.00000211663d0, 1.69085049749d0,   10977.07880469900d0, &
     &   0.00000178132d0, 4.06312103648d0,   17789.84561978500d0, &
     &   0.00000095668d0, 6.10260944466d0,     796.29800681640d0, &
     &   0.00000128829d0, 2.21903872695d0,   13367.97263110660d0, &
     &   0.00000096199d0, 4.01090243080d0,     213.29909543800d0, &
     &   0.00000079084d0, 0.00566665201d0,    5856.47765911540d0, &
     &   0.00000077808d0, 3.64174224007d0,   17298.18232732620d0, &
     &   0.00000074517d0, 1.76395145858d0,    6283.07584999140d0, &
     &   0.00000071832d0, 0.64118659399d0,    3154.68708489560d0, &
     &   0.00000066753d0, 4.20078421387d0,    5753.38488489680d0, &
     &   0.00000054305d0, 1.27898233220d0,   10447.38783960440d0, &
     &   0.00000047317d0, 5.98541885945d0,    4164.31198961300d0, &
     &   0.00000053850d0, 3.75922237603d0,    2544.31441988340d0, &
     &   0.00000041341d0, 5.15630887259d0,       7.11354700080d0, &
     &   0.00000048964d0, 5.06438171801d0,   21228.39202354580d0, &
     &   0.00000038553d0, 1.62992893340d0,     801.82093112380d0, &
     &   0.00000033520d0, 5.14670482356d0,    2146.16541647520d0, &
     &   0.00000044839d0, 5.35602472529d0,    4705.73230754360d0, &
     &   0.00000030006d0, 5.82045023860d0,    1589.07289528380d0, &
     &   0.00000039353d0, 4.40875564745d0,   26087.90314157420d0, &
     &   0.00000025552d0, 4.70498282243d0,    1194.44701022460d0, &
     &   0.00000034288d0, 2.87629134458d0,   16730.46368959580d0, &
     &   0.00000032230d0, 5.19482255479d0,   23543.23050468179d0, &
     &   0.00000031620d0, 0.20377436978d0,   25158.60171976540d0, &
     &   0.00000029988d0, 2.06647143238d0,    7084.89678111520d0, &
     &   0.00000022914d0, 0.26684134231d0,   11506.76976979360d0, &
     &   0.00000028371d0, 3.55121155276d0,   17260.15465469040d0, &
     &   0.00000024525d0, 5.33410545785d0,    6812.76681508600d0, &
     &   0.00000022400d0, 6.18159425612d0,     426.59819087600d0, &
     &   0.00000018808d0, 1.65360362287d0,    3738.76143010800d0, &
     &   0.00000017863d0, 0.26479139813d0,    1748.01641306700d0, &
     &   0.00000016727d0, 5.32163043011d0,      26.29831979980d0, &
     &   0.00000020783d0, 1.62629258894d0,   29088.81141598500d0, &
     &   0.00000019032d0, 3.07717685025d0,    3930.20969621960d0, &
     &   0.00000018419d0, 5.55814968279d0,    5486.77784317500d0, &
     &   0.00000018697d0, 0.89438266755d0,   13521.75144159140d0, &
     &   0.00000017980d0, 3.02828450256d0,    3340.61242669980d0, &
     &   0.00000017383d0, 3.49102971235d0,   11015.10647733480d0, &
     &   0.00000012241d0, 0.68301850852d0,    5088.62883976680d0, &
     &   0.00000015652d0, 0.75677278647d0,    5643.17856367740d0, &
     &   0.00000015213d0, 0.56985565245d0,    8635.94200376320d0, &
     &   0.00000014789d0, 2.84922701747d0,   18073.70493865020d0, &
     &   0.00000014393d0, 4.27640159260d0,   22003.91463486980d0, &
     &   0.00000013642d0, 6.00388395797d0,   12036.46073488820d0, &
     &   0.00000013834d0, 3.04875242870d0,   33019.02111220460d0, &
     &   0.00000010346d0, 2.96548563796d0,    3128.38876509580d0, &
     &   0.00000010485d0, 2.44285226960d0,   16200.77272450120d0, &
     &   0.00000011665d0, 5.69855413828d0,   25934.12433108940d0, &
     &   0.00000010374d0, 0.42204428363d0,    5230.80746680300d0, &
     &   0.00000010593d0, 3.78906764969d0,    5216.58037280140d0, &
     &   0.00000011214d0, 4.91512650529d0,   14945.31617355440d0, &
     &   0.00000008850d0, 0.85735777525d0,    9917.69687450980d0, &
     &   0.00000009429d0, 1.66843094483d0,    1349.86740965880d0, &
     &   0.00000007963d0, 1.17522727831d0,   52175.80628314840d0, &
     &   0.00000009297d0, 4.47115496366d0,   36949.23080842420d0, &
     &   0.00000008844d0, 0.83625571936d0,   29864.33402730900d0, &
     &   0.00000008447d0, 1.06163760227d0,    8429.24126646660d0, &
     &   0.00000006441d0, 0.25331823447d0,    4136.91043351620d0, &
     &   0.00000006339d0, 2.49567231737d0,    8031.09226305840d0, &
     &   0.00000006299d0, 2.08920832990d0,    4690.47983635860d0, &
     &   0.00000008297d0, 3.51264577993d0,   18422.62935909819d0, &
     &   0.00000007403d0, 4.71075536490d0,   23013.53953958720d0, &
     &   0.00000007598d0, 0.05572271846d0,   18875.52586977400d0, &
     &   0.00000007446d0, 3.30242876494d0,    1592.59601363280d0, &
     &   0.00000005472d0, 5.55412164007d0,    3634.62102451840d0, &
     &   0.00000006081d0, 3.07880786015d0,    4732.03062734340d0, &
     &   0.00000005289d0, 3.07331745335d0,     951.71840625060d0, &
     &   0.00000004781d0, 2.95678402973d0,     955.59974160860d0, &
     &   0.00000006498d0, 2.25653898617d0,   33794.54372352860d0, &
     &   0.00000004739d0, 4.39334857882d0,    2942.46342329160d0, &
     &   0.00000004611d0, 3.67855085148d0,   24072.92146977640d0, &
     &   0.00000004465d0, 0.55652742216d0,    7860.41939243920d0, &
     &   0.00000006294d0, 5.89350110582d0,   40879.44050464380d0, &
     &   0.00000006063d0, 4.04895816524d0,   22483.84857449259d0, &
     &   0.00000005750d0, 0.04561508156d0,   29296.61538957860d0, &
     &   0.00000004666d0, 4.42336639384d0,    7058.59846131540d0, &
     &   0.00000005286d0, 1.47927631912d0,   22805.73556599360d0, &
     &   0.00000004248d0, 1.22058351301d0,   18319.53658487960d0, &
     &   0.00000005216d0, 0.88886576578d0,   11926.25441366880d0, &
     &   0.00000004823d0, 5.88544271733d0,     155.42039943420d0, &
     &   0.00000004571d0, 1.35743621290d0,   14143.49524243060d0, &
     &   0.00000004245d0, 0.43954994331d0,     536.80451209540d0, &
     &   0.00000004254d0, 0.54107269293d0,    7238.67559160000d0, &
     &   0.00000004696d0, 3.67645105718d0,   37724.75341974820d0, &
     &   0.00000004020d0, 1.70026746262d0,   12566.15169998280d0, &
     &   0.00000004304d0, 2.82777079450d0,   11371.70468975820d0, &
     &   0.00000004059d0, 3.85797352990d0,     522.57741809380d0, &
     &   0.00000004420d0, 5.12393027283d0,     103.09277421860d0, &
     &   0.00000004286d0, 1.03260631515d0,   44809.65020086340d0, &
     &   0.00000004053d0, 1.14006806984d0,    4590.91018048900d0/
      data dk1/ &
     &   0.00082267418d0, 3.14159265359d0,       0.00000000000d0, &
     &   0.00000026605d0, 2.07229275260d0,     775.52261132400d0, &
     &   0.00000021081d0, 3.57698980914d0,    1059.38193018920d0, &
     &   0.00000019125d0, 1.01508722580d0,    4694.00295470760d0, &
     &   0.00000008380d0, 3.44586991840d0,       7.11354700080d0, &
     &   0.00000007457d0, 4.57985332526d0,     796.29800681640d0, &
     &   0.00000006463d0, 2.49473878135d0,    3154.68708489560d0, &
     &   0.00000005146d0, 0.58714547534d0,    4164.31198961300d0, &
     &   0.00000004067d0, 3.20010724578d0,    1194.44701022460d0, &
     &   0.00000004067d0, 3.49983542863d0,    4705.73230754360d0/
      data dk2/ &
     &   0.00002766670d0, 0.00000000000d0,       0.00000000000d0/
      data dk3/ &
     &   0.00000117290d0, 0.00000000000d0,       0.00000000000d0/

      data dh0/ &
     &   0.01628447663d0, 0.00000000000d0,       0.00000000000d0, &
     &   0.00001986929d0, 5.80464886318d0,    1577.34354244780d0, &
     &   0.00001864029d0, 2.12650300196d0,    5223.69391980220d0, &
     &   0.00001510978d0, 2.16070229051d0,     529.69096509460d0, &
     &   0.00000819544d0, 5.04224333254d0,    2352.86615377180d0, &
     &   0.00000483355d0, 0.94384676328d0,    5507.55323866740d0, &
     &   0.00000480730d0, 1.60400966048d0,   10213.28554621100d0, &
     &   0.00000448935d0, 5.94987750309d0,     398.14900340820d0, &
     &   0.00000354605d0, 1.70188277221d0,    4694.00295470760d0, &
     &   0.00000278987d0, 3.70742582004d0,    1059.38193018920d0, &
     &   0.00000294243d0, 2.35624478692d0,     775.52261132400d0, &
     &   0.00000230089d0, 2.36629249490d0,    9437.76293488700d0, &
     &   0.00000210039d0, 3.28388935733d0,   10977.07880469900d0, &
     &   0.00000178133d0, 2.49083562024d0,   17789.84561978500d0, &
     &   0.00000096628d0, 4.53886845461d0,     796.29800681640d0, &
     &   0.00000129280d0, 3.78871222562d0,   13367.97263110660d0, &
     &   0.00000096975d0, 2.43961766700d0,     213.29909543800d0, &
     &   0.00000078833d0, 1.57611478459d0,    5856.47765911540d0, &
     &   0.00000078224d0, 5.21110078908d0,   17298.18232732620d0, &
     &   0.00000073030d0, 0.17047151152d0,    6283.07584999140d0, &
     &   0.00000071889d0, 2.21128742509d0,    3154.68708489560d0, &
     &   0.00000054085d0, 2.86035616077d0,   10447.38783960440d0, &
     &   0.00000047355d0, 1.27342614649d0,    4164.31198961300d0, &
     &   0.00000042340d0, 3.63204564682d0,       7.11354700080d0, &
     &   0.00000055026d0, 5.31305130919d0,    2544.31441988340d0, &
     &   0.00000049329d0, 0.35027085630d0,   21228.39202354580d0, &
     &   0.00000038591d0, 3.20013618718d0,     801.82093112380d0, &
     &   0.00000044725d0, 3.78457637829d0,    4705.73230754360d0, &
     &   0.00000033434d0, 0.44799033642d0,    2146.16541647520d0, &
     &   0.00000030028d0, 4.24943152250d0,    1589.07289528380d0, &
     &   0.00000034257d0, 2.08062585709d0,   17260.15465469040d0, &
     &   0.00000038922d0, 2.82554883886d0,   26087.90314157420d0, &
     &   0.00000025568d0, 3.12836763274d0,    1194.44701022460d0, &
     &   0.00000032236d0, 3.62097537831d0,   23543.23050468179d0, &
     &   0.00000031929d0, 1.77259235162d0,   25158.60171976540d0, &
     &   0.00000030052d0, 3.63609619943d0,    7084.89678111520d0, &
     &   0.00000026998d0, 4.37231025874d0,   16730.46368959580d0, &
     &   0.00000024484d0, 3.76449325575d0,    6812.76681508600d0, &
     &   0.00000022368d0, 4.61049763041d0,     426.59819087600d0, &
     &   0.00000019819d0, 4.80088598059d0,    5753.38488489680d0, &
     &   0.00000018985d0, 0.09504486018d0,    3738.76143010800d0, &
     &   0.00000017742d0, 1.85524859155d0,    1748.01641306700d0, &
     &   0.00000016759d0, 0.60835730976d0,      26.29831979980d0, &
     &   0.00000021039d0, 3.19487995539d0,   29088.81141598500d0, &
     &   0.00000014503d0, 0.19585627636d0,   11506.76976979360d0, &
     &   0.00000018583d0, 1.51359089364d0,    3340.61242669980d0, &
     &   0.00000018586d0, 2.46923185056d0,   13521.75144159140d0, &
     &   0.00000018063d0, 0.89567248704d0,    5486.77784317500d0, &
     &   0.00000017442d0, 5.06023707430d0,   11015.10647733480d0, &
     &   0.00000016550d0, 1.33631601302d0,   18073.70493865020d0, &
     &   0.00000015631d0, 2.32467237202d0,    5643.17856367740d0, &
     &   0.00000015893d0, 2.75526767021d0,   22003.91463486980d0, &
     &   0.00000012020d0, 2.29584383924d0,    5088.62883976680d0, &
     &   0.00000015103d0, 5.28027396060d0,    8635.94200376320d0, &
     &   0.00000014044d0, 4.61713404227d0,   33019.02111220460d0, &
     &   0.00000013794d0, 4.46145453857d0,   12036.46073488820d0, &
     &   0.00000012899d0, 4.17819697737d0,   25934.12433108940d0, &
     &   0.00000010345d0, 1.39496302569d0,    3128.38876509580d0, &
     &   0.00000010621d0, 5.36114548056d0,    5216.58037280140d0, &
     &   0.00000010349d0, 1.99887324181d0,    5230.80746680300d0, &
     &   0.00000011266d0, 0.20076961236d0,   14945.31617355440d0, &
     &   0.00000008825d0, 2.43493655720d0,    9917.69687450980d0, &
     &   0.00000009711d0, 3.99457031445d0,    9917.69687450980d0, &
     &   0.00000008726d0, 3.23652502460d0,   23013.53953958720d0, &
     &   0.00000009835d0, 5.60203721863d0,   29864.33402730900d0, &
     &   0.00000009345d0, 3.26194362349d0,    1349.86740965880d0, &
     &   0.00000009466d0, 6.03935520837d0,   36949.23080842420d0, &
     &   0.00000007903d0, 5.87711198557d0,   52175.80628314840d0, &
     &   0.00000006484d0, 4.96979335761d0,    4136.91043351620d0, &
     &   0.00000006188d0, 3.69907141335d0,    4690.47983635860d0, &
     &   0.00000007859d0, 2.77416225550d0,    8429.24126646660d0, &
     &   0.00000008300d0, 1.94153363328d0,   18422.62935909819d0, &
     &   0.00000006081d0, 4.14906498513d0,    8031.09226305840d0, &
     &   0.00000007643d0, 1.62419959918d0,   18875.52586977400d0, &
     &   0.00000007426d0, 1.72027148458d0,    1592.59601363280d0, &
     &   0.00000007282d0, 0.74305710251d0,   33794.54372352860d0, &
     &   0.00000005475d0, 0.84206361823d0,    3634.62102451840d0, &
     &   0.00000006094d0, 4.64871534287d0,    4732.03062734340d0, &
     &   0.00000006571d0, 6.23951366290d0,   14143.49524243060d0, &
     &   0.00000005234d0, 4.66867318754d0,     951.71840625060d0, &
     &   0.00000004844d0, 4.51875245164d0,     955.59974160860d0, &
     &   0.00000006430d0, 1.17835886689d0,   40879.44050464380d0, &
     &   0.00000004612d0, 2.10710399676d0,   24072.92146977640d0, &
     &   0.00000004800d0, 2.44768544974d0,    7860.41939243920d0, &
     &   0.00000005753d0, 4.75282240478d0,   29296.61538957860d0, &
     &   0.00000004660d0, 2.85251167320d0,    7058.59846131540d0, &
     &   0.00000004244d0, 2.98876821470d0,    2942.46342329160d0, &
     &   0.00000005325d0, 3.04742256903d0,   22805.73556599360d0, &
     &   0.00000005308d0, 2.16743188943d0,   37724.75341974820d0, &
     &   0.00000004250d0, 5.93279864167d0,   18319.53658487960d0, &
     &   0.00000004766d0, 1.19968271915d0,     155.42039943420d0, &
     &   0.00000003742d0, 2.03947338751d0,    6681.22485339960d0, &
     &   0.00000004966d0, 2.44601751167d0,   11926.25441366880d0, &
     &   0.00000004293d0, 2.10316987651d0,    7238.67559160000d0, &
     &   0.00000004285d0, 5.15927478446d0,     536.80451209540d0, &
     &   0.00000003938d0, 0.10955967703d0,   12566.15169998280d0, &
     &   0.00000004712d0, 0.43341694722d0,     103.09277421860d0, &
     &   0.00000004053d0, 2.25469129195d0,     522.57741809380d0, &
     &   0.00000004394d0, 2.60051647343d0,   44809.65020086340d0, &
     &   0.00000004066d0, 2.71215397963d0,    4590.91018048900d0, &
     &   0.00000003834d0, 3.59188920421d0,   41654.96311596780d0, &
     &   0.00000003768d0, 4.47048984887d0,   26735.94526221320d0/
      data dh1/ &
     &   0.00062029655d0, 3.14159265359d0,       0.00000000000d0, &
     &   0.00000026631d0, 0.50260243452d0,     775.52261132400d0, &
     &   0.00000021214d0, 2.00566682331d0,    1059.38193018920d0, &
     &   0.00000019195d0, 2.58846455447d0,    4694.00295470760d0, &
     &   0.00000008642d0, 1.96182223710d0,       7.11354700080d0, &
     &   0.00000007680d0, 3.04144704600d0,     796.29800681640d0, &
     &   0.00000006469d0, 4.06548157609d0,    3154.68708489560d0, &
     &   0.00000005154d0, 2.15924623854d0,    4164.31198961300d0, &
     &   0.00000004073d0, 1.62252345380d0,    1194.44701022460d0, &
     &   0.00000004073d0, 1.93158131461d0,    4705.73230754360d0/
      data dh2/ &
     &   0.00003387470d0, 3.14159265359d0,       0.00000000000d0/
      data dh3/ &
     &   0.00000085620d0, 0.00000000000d0,       0.00000000000d0/

      data dq0/ &
     &   0.00000046990d0, 1.03836320801d0,     775.52261132400d0, &
     &   0.00000037030d0, 2.58501310328d0,    1059.38193018920d0, &
     &   0.00000023828d0, 2.93938256767d0,    3930.20969621960d0, &
     &   0.00000014714d0, 4.46415660357d0,    7860.41939243920d0, &
     &   0.00000009899d0, 2.46044981348d0,    4705.73230754360d0, &
     &   0.00000009873d0, 1.94962975227d0,     529.69096509460d0, &
     &   0.00000008903d0, 3.52488566870d0,    3154.68708489560d0, &
     &   0.00000007864d0, 2.90239730110d0,     426.59819087600d0, &
     &   0.00000009005d0, 5.16294437877d0,   12566.15169998280d0, &
     &   0.00000007970d0, 5.89298627379d0,   11790.62908865880d0, &
     &   0.00000006750d0, 3.88287479975d0,    8635.94200376320d0, &
     &   0.00000005196d0, 4.14011267868d0,    1577.34354244780d0, &
     &   0.00000004252d0, 2.99484359241d0,    1589.07289528380d0, &
     &   0.00000004395d0, 4.81752364842d0,     801.82093112380d0, &
     &   0.00000004735d0, 1.03635433956d0,   15720.83878487840d0, &
     &   0.00000004753d0, 0.95915716757d0,     337.81426319640d0/
      data dq1/ &
     &   0.00113472340d0, 3.14159265359d0,       0.00000000000d0, &
     &   0.00000004145d0, 2.55097496137d0,     775.52261132400d0/
      data dq2/ &
     &   0.00001236850d0, 0.00000000000d0,       0.00000000000d0/
      data dq3/ &
     &   0.00000126542d0, 0.00000000000d0,       0.00000000000d0/

      data dp0/ &
     &   0.00000048408d0, 5.76054381234d0,     775.52261132400d0, &
     &   0.00000036656d0, 1.01572916759d0,    1059.38193018920d0, &
     &   0.00000010147d0, 0.89874990533d0,    4705.73230754360d0, &
     &   0.00000009225d0, 5.08523515352d0,    3154.68708489560d0, &
     &   0.00000007856d0, 0.62393917012d0,    3930.20969621960d0, &
     &   0.00000007976d0, 1.32938474979d0,     426.59819087600d0, &
     &   0.00000009098d0, 3.59816659677d0,   12566.15169998280d0, &
     &   0.00000006891d0, 2.32046988968d0,    8635.94200376320d0, &
     &   0.00000004566d0, 0.57013141989d0,    7860.41939243920d0, &
     &   0.00000004236d0, 1.42387417254d0,    1589.07289528380d0, &
     &   0.00000004541d0, 0.09108817855d0,     801.82093112380d0, &
     &   0.00000003670d0, 0.53038616499d0,       7.11354700080d0, &
     &   0.00000004753d0, 2.52994390566d0,     337.81426319640d0, &
     &   0.00000003523d0, 0.30506773876d0,   20426.57109242200d0, &
     &   0.00000003152d0, 0.22265746056d0,    7084.89678111520d0, &
     &   0.00000003219d0, 2.50965126428d0,   11506.76976979360d0/
      data dp1/ &
     &   0.00010179800d0, 0.00000000000d0,       0.00000000000d0, &
     &   0.00000004328d0, 1.01402968363d0,     775.52261132400d0/
      data dp2/ &
     &   0.00004702790d0, 0.00000000000d0,       0.00000000000d0/
      data dp3/ &
     &   0.00000054230d0, 3.14159265359d0,       0.00000000000d0/

      t=(dj-dj2000)/365250.d0

      dkhqp(1)=0.d0
      do j=1,dimkh(1)
        dkhqp(1)=dkhqp(1)+dk0(1,j)*cos(dk0(2,j)+dk0(3,j)*t)
      enddo
      do j=1,dimkh(2)
        dkhqp(1)=dkhqp(1)+t*dk1(1,j)*cos(dk1(2,j)+dk1(3,j)*t)
      enddo
      dkhqp(1)=dkhqp(1)+t**2*dk2(1)*cos(dk2(2)+dk2(3)*t)
      dkhqp(1)=dkhqp(1)+t**3*dk3(1)*cos(dk3(2)+dk3(3)*t)

      dkhqp(2)=0.d0
      do j=1,dimkh(1)
        dkhqp(2)=dkhqp(2)+dh0(1,j)*cos(dh0(2,j)+dh0(3,j)*t)
      enddo
      do j=1,dimkh(2)
        dkhqp(2)=dkhqp(2)+t*dh1(1,j)*cos(dh1(2,j)+dh1(3,j)*t)
      enddo
      dkhqp(2)=dkhqp(2)+t**2*dh2(1)*cos(dh2(2)+dh2(3)*t)
      dkhqp(2)=dkhqp(2)+t**3*dh3(1)*cos(dh3(2)+dh3(3)*t)

      dkhqp(3)=0.d0
      do j=1,dimqp(1)
        dkhqp(3)=dkhqp(3)+dq0(1,j)*cos(dq0(2,j)+dq0(3,j)*t)
      enddo
      do j=1,dimqp(2)
        dkhqp(3)=dkhqp(3)+t*dq1(1,j)*cos(dq1(2,j)+dq1(3,j)*t)
      enddo
      dkhqp(3)=dkhqp(3)+t**2*dq2(1)*cos(dq2(2)+dq2(3)*t)
      dkhqp(3)=dkhqp(3)+t**3*dq3(1)*cos(dq3(2)+dq3(3)*t)

      dkhqp(4)=0.d0
      do j=1,dimqp(1)
        dkhqp(4)=dkhqp(4)+dp0(1,j)*cos(dp0(2,j)+dp0(3,j)*t)
      enddo
      do j=1,dimqp(2)
        dkhqp(4)=dkhqp(4)+t*dp1(1,j)*cos(dp1(2,j)+dp1(3,j)*t)
      enddo
      dkhqp(4)=dkhqp(4)+t**2*dp2(1)*cos(dp2(2)+dp2(3)*t)
      dkhqp(4)=dkhqp(4)+t**3*dp3(1)*cos(dp3(2)+dp3(3)*t)

      return
      end

      SUBROUTINE SUBPSIOM(dj,psiom)
      IMPLICIT NONE
      DOUBLE PRECISION dj, dj2000, pi, t, xxx
      DOUBLE PRECISION psiom(2)
      DOUBLE PRECISION psi0(3,33),psi1(3,6),psi2(3,7),psi3(3,2),psi4(3)
      DOUBLE PRECISION om0(3,33),om1(3,6),om2(3,7),om3(3,2),om4(3)
      INTEGER j, dimpo(5)
      data dimpo/33,6,7,2,1/
      data dj2000/2451545.d0/
      data psi0/ &
     &   17206641.474d0, 0.611445269d0,     -337.5704708d0, &
     &    1318626.533d0, 1.936645871d0,    12566.6393056d0, &
     &     227626.194d0, 6.048628418d0,   167994.1822020d0, &
     &     207461.312d0, 5.935320674d0,     -675.1409416d0, &
     &     128455.042d0, 1.461501441d0,     6283.3196528d0, &
     &      71119.766d0, 3.925112589d0,    83286.9142477d0, &
     &      51707.381d0, 1.893812078d0,    18849.9589584d0, &
     &      38729.439d0, 3.866282942d0,   168331.7526728d0, &
     &      30142.671d0, 2.119466369d0,   251281.0964497d0, &
     &      15699.729d0, 3.327895623d0,    72140.6286486d0, &
     &      12814.834d0, 2.897015352d0,    12904.2097764d0, &
     &      12352.571d0, 0.552584902d0,    84707.2679542d0, &
     &       6347.520d0, 5.682920310d0,   155427.5428964d0, &
     &       6310.847d0, 6.109231943d0,    82949.3437769d0, &
     &       5962.992d0, 1.522300656d0,   240134.8108506d0, &
     &       5799.404d0, 4.542669292d0,   -83624.4847185d0, &
     &       5160.805d0, 6.220315307d0,   251618.6669205d0, &
     &       4771.281d0, 4.115622859d0,   -11146.2855991d0, &
     &       4591.630d0, 2.298450889d0,     1757.9241773d0, &
     &       3856.313d0, 3.876270881d0,   323421.7250983d0, &
     &       3104.231d0, 4.473252858d0,   334568.0106974d0, &
     &       2925.209d0, 6.279521207d0,   166573.8284954d0, &
     &       2862.235d0, 1.149500341d0,    95853.5535533d0, &
     &       2589.629d0, 4.823896749d0,   168669.3231436d0, &
     &       2177.827d0, 2.429254786d0,   -13241.7802472d0, &
     &       2044.959d0, 4.653771385d0,    85044.8384250d0, &
     &       1708.476d0, 0.528529857d0,    -6620.8901236d0, &
     &       1584.065d0, 1.852104693d0,    25133.2786112d0, &
     &       1516.871d0, 5.512100210d0,    71803.0581778d0, &
     &       1406.399d0, 0.561668193d0,     5945.7491820d0, &
     &       1287.304d0, 5.139807391d0,   -72478.1991195d0, &
     &       1102.087d0, 3.258399637d0,     2095.4946481d0, &
     &       1020.267d0, 5.622978530d0,   240472.3813214d0/
      data psi1/ &
     & -50384564880.869d0, 0.000000000d0,        0.0000000d0, &
     &        85078.637d0, 0.769402051d0,     -337.5704708d0, &
     &        49893.061d0, 6.227018463d0,     6283.3196528d0, &
     &        15568.378d0, 0.224782731d0,    18849.9589584d0, &
     &        14103.294d0, 0.925955348d0,  2301216.7536515d0, &
     &         6043.891d0, 5.059057568d0,    12566.6393056d0/
      data psi2/ &
     &    107194829.468d0, 0.000000000d0,        0.0000000d0, &
     &        35665.735d0, 2.215176555d0,     -337.5704708d0, &
     &         6230.483d0, 4.368719914d0,     6283.3196528d0, &
     &         2402.183d0, 4.852828549d0,    18849.9589584d0, &
     &         1535.294d0, 4.467680770d0,   167994.1822020d0, &
     &         1240.113d0, 1.235424136d0,     -675.1409416d0, &
     &         1083.476d0, 5.500284858d0,    83286.9142477d0/
      data psi3/ &
     &      1143656.495d0, 0.000000000d0,        0.0000000d0, &
     &         1344.571d0, 1.236660169d0,     -337.5704708d0/
      data psi4/ &
     &     -1328317.884d0, 0.000000000d0,        0.0000000d0/

      data om0/ &
     & -84381409000.000d0, 0.000000000d0,        0.0000000d0, &
     &      9205151.971d0, 5.323867953d0,     -337.5704708d0, &
     &       573040.588d0, 0.366174057d0,    12566.6393056d0, &
     &        97840.400d0, 4.477657022d0,   167994.1822020d0, &
     &        89755.029d0, 4.364557132d0,     -675.1409416d0, &
     &        22439.282d0, 0.323015641d0,    18849.9589584d0, &
     &        20072.758d0, 2.294821149d0,   168331.7526728d0, &
     &        16647.291d0, 3.373699722d0,     6283.3196528d0, &
     &        12902.410d0, 0.548538206d0,   251281.0964497d0, &
     &         6901.459d0, 1.324873162d0,    12904.2097764d0, &
     &         5330.577d0, 5.265015033d0,    84707.2679542d0, &
     &         3322.456d0, 4.537664710d0,    82949.3437769d0, &
     &         3142.820d0, 2.970988600d0,   -83624.4847185d0, &
     &         2636.682d0, 4.649134152d0,   251618.6669205d0, &
     &         2554.007d0, 6.234595942d0,   240134.8108506d0, &
     &         2422.688d0, 0.728105741d0,     1757.9241773d0, &
     &         1645.201d0, 2.305355427d0,   323421.7250983d0, &
     &         1323.691d0, 2.902507060d0,   334568.0106974d0, &
     &         1233.743d0, 5.862182822d0,    95853.5535533d0, &
     &         1075.803d0, 3.082469592d0,    85044.8384250d0, &
     &          857.264d0, 5.286182958d0,     5945.7491820d0, &
     &          800.493d0, 3.940488735d0,    71803.0581778d0, &
     &          695.709d0, 3.568204712d0,   -72478.1991195d0, &
     &          685.864d0, 0.280494104d0,    25133.2786112d0, &
     &          674.703d0, 2.408976642d0,    83286.9142477d0, &
     &          522.227d0, 4.051664597d0,   240472.3813214d0, &
     &          497.627d0, 4.184531777d0,       -7.8449519d0, &
     &          414.461d0, 5.650818169d0,    -6620.8901236d0, &
     &          335.193d0, 0.122363251d0,   323759.2955691d0, &
     &          327.112d0, 3.152993847d0,   155089.9724256d0, &
     &          326.954d0, 4.658649906d0,   406708.6393460d0, &
     &          325.028d0, 1.292657215d0,   174277.5018547d0, &
     &          305.745d0, 4.522535824d0,   161710.8625492d0/
      data om1/ &
     &       265011.262d0, 0.000000000d0,        0.0000000d0, &
     &        10168.539d0, 5.770898235d0,     -337.5704708d0, &
     &         6753.746d0, 4.934337891d0,    18849.9589584d0, &
     &         5609.966d0, 5.638344460d0,  2301216.7536515d0, &
     &         3011.032d0, 3.380408067d0,    12566.6393056d0, &
     &         1656.785d0, 6.112431360d0,     6283.3196528d0/
      data om2/ &
     &     -5127689.571d0, 0.000000000d0,        0.0000000d0, &
     &        25453.279d0, 0.611923477d0,     -337.5704708d0, &
     &         1039.165d0, 3.272875387d0,    18849.9589584d0, &
     &          723.538d0, 0.169638986d0,     6283.3196528d0, &
     &          642.263d0, 2.906202613d0,   167994.1822020d0, &
     &          554.189d0, 5.935823188d0,     -675.1409416d0, &
     &          200.721d0, 0.723759607d0,   168331.7526728d0/
      data om3/ &
     &      7727159.659d0, 0.000000000d0,        0.0000000d0, &
     &          300.348d0, 0.778609438d0,     -337.5704708d0/
      data om4/ &
     &         4917.017d0, 0.000000000d0,        0.0000000d0/

      pi=3.1415926535897932d0
      xxx=6.48d11/pi
      t=(dj-dj2000)/365250.d0

      psiom(1)=0.d0
      do j=1,dimpo(1)
        psiom(1)=psiom(1)+psi0(1,j)*cos(psi0(2,j)+psi0(3,j)*t)
      enddo
      do j=1,dimpo(2)
        psiom(1)=psiom(1)+t*psi1(1,j)*cos(psi1(2,j)+psi1(3,j)*t)
      enddo
      do j=1,dimpo(3)
        psiom(1)=psiom(1)+t**2*psi2(1,j)*cos(psi2(2,j)+psi2(3,j)*t)
      enddo
      do j=1,dimpo(4)
        psiom(1)=psiom(1)+t**3*psi3(1,j)*cos(psi3(2,j)+psi3(3,j)*t)
      enddo
      psiom(1)=psiom(1)+t**4*psi4(1)*cos(psi4(2)+psi4(3)*t)

      psiom(2)=0.d0
      do j=1,dimpo(1)
        psiom(2)=psiom(2)+om0(1,j)*cos(om0(2,j)+om0(3,j)*t)
      enddo
      do j=1,dimpo(2)
        psiom(2)=psiom(2)+t*om1(1,j)*cos(om1(2,j)+om1(3,j)*t)
      enddo
      do j=1,dimpo(3)
        psiom(2)=psiom(2)+t**2*om2(1,j)*cos(om2(2,j)+om2(3,j)*t)
      enddo
      do j=1,dimpo(4)
        psiom(2)=psiom(2)+t**3*om3(1,j)*cos(om3(2,j)+om3(3,j)*t)
      enddo
      psiom(2)=psiom(2)+t**4*om4(1)*cos(om4(2)+om4(3)*t)

      psiom(1)=psiom(1)/xxx
      psiom(2)=psiom(2)/xxx

      return
      end

      SUBROUTINE update_irn
! *** =========================================================================
! *** Purpose : update irn variable (surface type) in case of icemask forcing
! *** Original code : C. Van Meerbeek (XXXX)
! *** Modified code : P. Mathiot (01/2012) Allow interannual simulation
! *** =========================================================================

      USE NETCDF    ! netcdf module to read forcing data

      IMPLICIT NONE

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'icemask.h'

!0 variable
!      INTEGER, DIMENSION(nlat, nlon):: icemask
!      REAL(KIND=4) :: totalprecipcases
      INTEGER :: idd_time, idf_imsk, idv_time, idv_imsk, istatus        ! id netcdf file 0 icemask
      INTEGER :: ntime_lgm, itime_lgm                                   ! time dimension of icemask file
      INTEGER, DIMENSION(:)  , ALLOCATABLE :: nvtime_lgm                ! time of each icemask data

!      common /lev_tarasov_forcing/icemask,totalprecipcases

! temporary variable
      INTEGER :: i                                                      ! Loop index
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: zicemask                  ! temporary icemask data
      integer zipl(27)
      real*4  zpisccp(27),zpncep(17),zz500ncep(27,12)

! netcdf version + interannual version : Pierre Mathiot (12/11/11) UCL
!*** 0 : adapt LWR scheme ==>
!          use k=27 Greenland values for icecaps in Northern Hemisphere
!dmr 0--- With use of the T21 ICE-5G icemask ...

! LOAD netcdf icemask file
      istatus=NF90_OPEN("inputdata/icemask.nc", NF90_NOWRITE, idf_imsk) !ouvre le fichier
      istatus=NF90_INQ_DIMID(idf_imsk, 'time', idd_time) !recupère l'id du temps
      istatus=NF90_INQUIRE_DIMENSION(idf_imsk, idd_time, len = ntime_lgm) !recupère à partir de l'id, le nombre de pas de temps
      istatus=nf90_inq_varid(idf_imsk, 'time', idv_time) !recupère l'id du temps
      istatus=nf90_inq_varid(idf_imsk, "imask", idv_imsk) !récupère l'id de la variable Sul

! ALLOCATE variable
      ALLOCATE(nvtime_lgm(ntime_lgm), zicemask(nlon, nlat))

! READ time variable
      istatus=nf90_get_var(idf_imsk, idv_time, nvtime_lgm)

! SELECT right time
      itime_lgm=0
      DO i=1,ntime_lgm
         IF (INT(ABS((nvtime_lgm(i)-(irunlabelf+iyear)))) == INT(MINVAL(ABS(nvtime_lgm(:)-(irunlabelf+iyear))))) itime_lgm=i
      END DO
      IF (itime_lgm==0) THEN
         PRINT *, "ic_irn : error in detection time slide = ", itime_lgm
         STOP
      END IF

! READ icemask variable
      istatus=nf90_get_var(idf_imsk, idv_imsk, zicemask, start = (/1,1,itime_lgm/), count = (/nlon,nlat,1/)) !charge les valeurs dans la variable sulopt
      icemask=TRANSPOSE(zicemask)
! 6 is not to take antarctica, but the World Everywhere otherwise
      WHERE (fracto(6:nlat,:) .GT. 0.99)
          icemask(6:nlat,:) = 0
      END WHERE

! READ irn before applied icemask file
      open(111,file='inputdata/lwrref.dat',form='unformatted')
      read(111) irn,zipl,zpisccp,zpncep,zz500ncep
      close(111)

! UPDATE irn
      WHERE (INT(icemask(6:nlat,:)) == 1)
          irn(6:nlat,:,2) = 27
      END WHERE

! CLOSE netcdf file
      istatus=nf90_close(idf_imsk) ! Close icemask file

      WRITE(iuo+99,*) "================================================================="
      WRITE(iuo+99,*) "========================== irn ==============================="
      WRITE(iuo+99,*) "Read icemask in file ",TRIM('inputdata/icemask.nc')
      WRITE(iuo+99,*) "For year (AD) = ",irunlabelf+iyear
      WRITE(iuo+99,*) "selected time is ", nvtime_lgm(itime_lgm)," index is ",itime_lgm,"/",ntime_lgm
      WRITE(iuo+99,*) "================================================================="
      CALL FLUSH(iuo+99)

! DEALLOCATE variable
      DEALLOCATE(zicemask, nvtime_lgm)

      END
