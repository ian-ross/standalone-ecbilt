












c- issu de la routine "TVM" ; modif (JMC) : 08/10/00; inclu
c- dans emic (driess) :30/08/2002
*********************************************************************
*      Terrestrial Vegetation annual Model
*
*  Purpose: Computation of forest/grass/desert ratio;
*           npp and living and soil biomass
*
*  By: V.Brovkin
*  Modifyed by: A.Ganopolski
*  Last modification: 17.11.97
**********************************************************************
c- nyvegt = frequence (in Years) of "call veget"
c      if = 0 --> return flgveg=False and do not "call veget" anymore.
c- nwrveg = frequence (in Years) of ASCII output (0 => no output)
c---------
c-First call of "veget" (iyear=0) : read "veget.par"
c kveget < 0  => compute initial vegetation using T,Pre,GDD0 from "veget.init"
c kveget > 0  => read vegetation on "veget.rest"
c   after 1rst call : kveget =abs(kveget)  and other param. are kept constant
c-In Running Loop (iyear > 0) :
c kveget = 0  => no effect on Atmospheric Model (Albedo not changed)
c         but write T_moy,Pre_moy & GDD0 on file "veget.init" & "veget.outp"
c kveget = 1  => synchronous run of vegetation model
c kveget > 1  => do "kveget" iterations of vegetation to get faster equilibrium
c----------
c   prcmin= minimun daily precip (in m) available for vegetation in warm area.
c bmtdry = Soil-Water threshold(m) to compute "tpsdry"=Nb_days/yr of dry soil C.
c tmxdry : above this limit(Nb_days/yr), soil_dryness start to affect vegetation
c dtrdry = time interval(days) for linear transition [tmxdry,tmxdry+dtrdry] ;
c  If tpsdry(=Nb_day_dry) > tmxdry+dtrdry => Precip(Veget) is reduced by rpfdry
c rpfdry = reduction factor applied to Precip(Veget) if dry condition satisfied

*********************************************************************
c input data: annual mean temperature in degr. Celc. - ave_t
c             annual mean precipitation, mm/year - ave_pr
c             growing degree days above 0, degr. Celc - gdd0
c             CO2 concentration in the atmosphere, ppm - co2
c
c output data: st(i,j) - forest's share in cell i,j
c              sg(i,j) - grass share in cell i,j
c              sd(i,j) - desert's share in cell i,j
c              snlt(i,j) -needle leaved trees ratio
c              st(i,j)+sg(i,j)+sd(i,j)=1, 0<= snlt(i,j)<= 1
c              plai(i,j) - annual average of leaf area index, m2/m2
c              anup(i,j) - annual uptake of carbon in cell, Gt -> kg/m2 ???
c
c COMMON block described in buffer.inc
*********************************************************************
      subroutine veget(ist,jst,dtime,epss,patmCO2,fracgr,darea,tempgr) 
**************************************************************************
      implicit none


*       include 'declar.inc'
*       include 'params.inc'
C       include 'bio.inc'
C       include 'buffer.inc'
      include 'veget.h'
      include 'comphys.h'
      include 'comrunlabel.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comdiag.h'
      include 'netcdf.inc'
      real*8 veg_albd
*********************************************************************
      real*8    albsnow(nlat),albland(nlat,nlon,4)
      real*8    bmoisg(nlat,nlon),forestfr(nlat,nlon)
      real*8    bmoism(nlat,nlon)
C     real*8    rs(nlat,nlon)
      real*8    alblbm(nlat,nlon),rainf(nlat,nlon),snowf(nlat,nlon)
      real*8    alblbmR(nlat,nlon),alblandR(nlat,nlon,4)
      real*8    alblandismR(nlat,nlon,4,3),forestfrR(nlat,nlon)
      common /ec_lbmbmois/ bmoisg,bmoism,rainf,snowf
      common /ec_lbmcalbedo/albsnow,albland,forestfr,alblbm
      common /ec_lbmcalbedo2/alblandR,forestfrR,alblbmR,alblandismR

c--dummy variables :
c- input (except 1rst call) :
      integer ist,jst, istep
      real*8 dtime, epss
      real*8 fracgr(nlat,nlon), darea(nlat), tempgr(nlat,nlon)
c- output :
      character*6 endyear
      character*3 endday
      logical existe
      
c--local variables saved (in common) from one call to the other :
      logical flgwrt
      integer kveget, nyvegt, nwrveg, ncumvg, iyr0vg, ns0inv
      integer ieq
      real*8 tmxdry, dtrdry, rpfdry, unsdry
      real*8 temveg(nlat,nlon), gd0veg(nlat,nlon),
     &     prcveg(nlat,nlon,2), tpsdry(nlat,nlon)
      real*8 prcmin, bmtdry, prcday(nlat,nlon)
      character*15 titveg, titv00
      common / com0veg / kveget, nyvegt, nwrveg, ncumvg,
     &                   iyr0vg, ns0inv(4,0:1)
      common / com1veg / tmxdry,dtrdry,rpfdry,unsdry,
     &                   temveg,gd0veg,prcveg,tpsdry
      common / com2veg / prcmin, bmtdry, prcday
      common / comcveg / titveg

        
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c--local variables :
      integer i,j,k,ns,ii,kk2,nlatd2, iyrloc, ireg
      real*8 zero, one, xx, xxalb, zlai, ddc, bmtd00
      real*8 dzero, ttemp, prcrit, zmoy1, zmoy2, zmoy3, zmoy4
      real*8 albet(4), albeg(4), albed(4), albegc(4), albedc(4)
      real*8 xxalbR
      integer istart,igroen 
      integer ios,ios2
      double precision patmCO2
      character*60 part1
      CHARACTER(len=256) :: output_filename
      real*8  ari(nlat/2)
      real*8  rj,rw,dumwei
      integer ilat,ird, status
 
c- Albedo of tree, grass, dessert for Winter,Spring,Summer and Fall:
c     data albet / 0.13 , 0.13 , 0.13 , 0.13 /
c     data albeg / 0.20 , 0.20 , 0.20 , 0.20 /
c     data albed / 0.33 , 0.33 , 0.33 , 0.33 /
c- Albedo of grass over cold area (Steppe) = albeg + albegc
c     data albegc / -.06 , -.04 , -.02 , -.04 /
c- Albedo of bright sand desert (Sahara) = albed + albedc
c     data albedc / 0.07 , 0.07 , 0.07 , 0.07 /
c--For output :
      data flgwrt / .TRUE. /
      logical veg_act
    
      zero = 0.
      one  = 1.
      dzero = 0.d0

c-----
      if(initialization.eqv..true.) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) Read "veget.par" and Define parameter for this routine.         |
c-----------------------------------------------------------------------

c--Initialization of Cumul.array

c--ouverture de "veget.uptake" pour le calcul du flux de carbon
c--entre latmosphere et la biomasse
 
c      open(iveg+3,file='veget.uptake',status='unknown')

c- Read Vegetation parameter :
      open(iveg+1,file='veget.par',status='old')
      read(iveg+1,*)
      read(iveg+1,*)
      read(iveg+1,'(A)') titveg
      read(iveg+1,*)
      read(iveg+1,*) kveget, nyvegt, nwrveg
      read(iveg+1,*)
      read(iveg+1,*) prcmin, bmtdry, tmxdry, dtrdry, rpfdry
      read(iveg+1,*)
      read(iveg+1,*) ieqveg,ieqvegwr
      read(iveg+1,*)
      read(iveg+1,*) iscendef
      read(iveg+1,*)
      read(iveg+1,*) (albet(i),i=1,4)
      read(iveg+1,*)
      read(iveg+1,*) (albeg(i),i=1,4)
      read(iveg+1,*)
      read(iveg+1,*) (albed(i),i=1,4)
      read(iveg+1,*)
      read(iveg+1,*) (albegc(i),i=1,4)
      read(iveg+1,*)
      read(iveg+1,*) (albedc(i),i=1,4)
      read(iveg+1,*)
      read(iveg+1,*) gamm2
      read(iveg+1,*)
      read(iveg+1,*) betag,betat
      read(iveg+1,*)
      read(iveg+1,*) fco2veg
      close(iveg+1)

C *** gauss points and weights
      rewind(iuo+7)
      ilat=nlat/2
   10 continue
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
   15 continue
C     call ec_error(4)
   20 continue

      do i=1,ilat
        phi(i)=-ari(ilat+1-i)
        phi(ilat+i)=ari(i)
      enddo

      do i=1,nlat
        phi(i)=asin(phi(i))
      enddo

 220  format(f18.10,f17.10)

      open(iveg+65,file='outp_veget.param')
      read(iveg+65,'(/,/,/,/,/,/,/,/)')
      read(iveg+65,*) numvegvar,veg_fill_value,veg_missing_value
      read(iveg+65,*)

      do i=1,numvegvar

         read(iveg+65,"(A)") part1
         namevegvar(i,1)=trim(part1)
         read(iveg+65,*) (namevegvar(i,k),k=2,5)
         read(iveg+65,*) (newvegvar(i,k),k=1,7)
         do k=1,6
            newvegvar(i,k)=newvegvar(i,k+1)
         enddo

C         IF ( newvegvar(i,4)==1 ) ioutyearly = 1
C         IF ((newvegvar(i,3)==1).OR.(newvegvar(i,2)==1) ) meantype = 1
C         IF ( newvegvar(i,3)==1 ) meanyl     = 1
C         IF ( newvegvar(i,2)==1 ) meantot    = 1
C         IF ( newvegvar(i,1)==1 ) ioutdaily  = 1

      enddo

      close(iveg+65)

c-AM 
C beta_i divided once  for all by log(2)
      betat=betat/log(2.0)
      betag=betag/log(2.0)

c-AM (2008)
c***  read deforestation
c-YSD (2011)
c***  Netcdfisation
      i0dfor=0
      ndfor=0
      if (iscendef.eq.1) then
        status=nf_open("inputdata/VEGET.nc", nf_nowrite, ireg) !ouvre le fichier
        status=nf_inq_dimid(ireg, 'time', i) !recupère l'id du temps
        status=nf_inq_dimlen(ireg, i, ndfor) !recupère à partir de l'id, le nombre de pas de temps
        if(ndfor.gt.mdfor) STOP 'Please adjust  mdfor  in veget.h'
        status=nf_inq_varid(ireg, 'time', i) !recupère l'id du temps
        status=nf_get_vara_real(ireg, i, (/1/), (/ndfor/), VegetTime)
        ivegstrt=int(VegetTime(1)) !la 1ere valeur du temps donne la date de debut du forcage
        status=nf_inq_varid(ireg, "vegfrac", j) !récupère l'id de la variable Sul
        status=nf_get_vara_real(ireg, j, (/1,1,1/), (/64,32,ndfor/), farea) !charge les valeurs dans la variable sulopt
        ivegstrt=max(i0dfor,ivegstrt)
        write(iuo+99,*) 'scen Veget start=',ivegstrt, "AD nbline=",ndfor
      endif

c
c-AM (2008)
c  in case of constant vegetation distribution
      if (iscendef.eq.-1) i0dfor=ivegstrt

c- check param :
      unsdry = 0.0
      if (dtrdry.gt.1.0e-6) unsdry = 1.0 / dtrdry

      if (nyvegt.eq.0) then
        nyvegt = nyears + 1
        kveget = abs(kveget)
        flgveg = .FALSE.
        return
      endif
      nwrveg = (nwrveg/nyvegt) * nyvegt
      if (nwrveg.eq.0) then
        nwrveg = nyears + 1
        flgwrt = .FALSE.
      endif
      iyr0vg = 0

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Before the 1rst time step of the run : Setup Initial Vegetation.|
c-----------------------------------------------------------------------

c--reverse season index in South.H. = ns0inv(ns,0)
      do ns=1,4
        ns0inv(ns,0) = 1 + mod(ns+1,4)
        ns0inv(ns,1) = ns
      enddo

c- Initialisation of Vegetation Parameters :
      if (kveget.ne.0) call initcpar

      if (kveget.gt.0) then
c- Read Vegetation from file "veget.rest".
      open(iveg+2,file='veget.rest',status='old', form='unformatted')
      read(iveg+2) iyr0vg,bmtd00,titv00
      read(iveg+2) st
      read(iveg+2) sg
      read(iveg+2) sd
      read(iveg+2) snlt
      read(iveg+2) temveg
      read(iveg+2) b1t
      read(iveg+2) b1g
      read(iveg+2) b2t
      read(iveg+2) b2g
      read(iveg+2) b3t
      read(iveg+2) b3g
      read(iveg+2) b4t
      read(iveg+2) b4g
      read(iveg+2) b1
      read(iveg+2) b2
      read(iveg+2) b3
      read(iveg+2) b4
      read(iveg+2) anup
c
c iscendef = -1 : veget does not evolve (maintained at its initial distribution)
c          = +1 : land-use scenario 
c          =  0 : no constraint
c
      if ((iscendef.eq.1).or.(iscendef.eq.-1)) then
c Reference vegetation (note no need to save sd_const <- AM)
c ivegstrt define the reference year for anomalies of cropland
       read(iveg+2,end=99,iostat=ios) st_const
   99  if ((ios.ne.0).or.(irunlabelf.eq.ivegstrt)) st_const(:,:)=st(:,:)
       sd_const(:,:)=sd(:,:)
       read(iveg+2,end=98,iostat=ios2)str
       read(iveg+2,end=98,iostat=ios2)sgr
       read(iveg+2,end=98,iostat=ios2)sdr
      endif 
   98 if ((ios2.ne.0).or.(iscendef.ne.1)) then 
       sdr(:,:)=sd(:,:)
       sgr(:,:)=sg(:,:)
       str(:,:)=st(:,:)
      endif 
      close(iveg+2) 
      write(iuo+66,'(A)') 'Initial Vegetation <- Read file "veget.rest".'
      write(iuo+66,'(A,I6,A,F6.3,2A)') ' year:', iyr0vg,
     &        ' Dry.Soil Lim=', bmtd00, ' ; title: ', titv00
      snltr(:,:)=snlt(:,:)

      do i=1,nlat
        do j=1,nlon
          stock(i,j)=b1(i,j)+b2(i,j)+b3(i,j)+b4(i,j)
          stockloch(i,j)=fracgr(i,j)*darea(i)*stock(i,j)*1E-12
          anuploch(i,j)=fracgr(i,j)*darea(i)*anup(i,j)*1E-12
          anuploch(i,j)=fco2veg*anuploch(i,j)
        enddo
      enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-----
      elseif (kveget.lt.0) then
c- Read from file=veget.init : An_Mean_Temp, gdd0 & An_Mean_Precip.
      open(iveg+2,file='veget.init',status='old', form='unformatted')
      read(iveg+2) temveg
      read(iveg+2) gd0veg
      read(iveg+2) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
      read(iveg+2) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
      read(iveg+2) tpsdry
      read(iveg+2) iyrloc, bmtd00, titv00
      close(iveg+2)
      write(iuo+66,'(A)') 'Initial Veg. computed from file "veget.init" :'
      write(iuo+66,'(A,I6,A,F6.3,2A)') ' year:', iyrloc,
     &        ' Dry.Soil Lim=', bmtd00, ' ; title: ', titv00
c-----
      endif

      else
c------------------------------------------------
c- during time integration : <--> kveget > or = 0


*************************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) current time step : cum. atmospheric fields used for vegetation :
c------------------------------------------------

      istep=(ist-1)*iatm+jst

      ncumvg = ncumvg + 1
      do j=1,nlon
       do i=1,nlat
         ttemp = tempgr(i,j)-273.15d0
         temveg(i,j)=temveg(i,j)+ttemp
c- torain is in m/s ; dtime = length of 1 iter, in second
         prcday(i,j)=prcday(i,j)+dtime*(torain(i,j)+tosnow(i,j))
         gd0veg(i,j)=gd0veg(i,j)+max(ttemp,dzero)
       enddo
      enddo
 
c- last iter. of the day :
      if (mod(istep,iatm).eq.0) then
        prcrit = abs(prcmin)
        do j=1,nlon
         do i=1,nlat
c- total precipitation :
           prcveg(i,j,1) = prcveg(i,j,1)+prcday(i,j)
c- precipitation above the daily threshold "prcmin" :
           if (prcday(i,j).ge.prcrit)
     &       prcveg(i,j,2)=prcveg(i,j,2)+prcday(i,j)
           prcday(i,j) = 0.
c- cumul Nb_days of dry soil conditions :
           if (bmoisg(i,j).le.bmtdry)
     &        tpsdry(i,j) = tpsdry(i,j) + 1.
         enddo
        enddo
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!       if (mod(istep,nstpyear*nyvegt).eq.0) then
! SDubinkina, needed for restart every N month, N<12
      if (mod(irunlabeld*iatm+istep,nstpyear*nyvegt).eq.0) then
c-----------------------------------------------------------
c--Last time step of the year : compute Anual_Mean Climato :
        zmoy1 = 0.
        if (ncumvg.ge.1) zmoy1 = 1.d0 / DBLE(ncumvg)
        zmoy2 = zmoy1 * 360.
c- conversion : hauteur cumulee (m) -> precip. (mm/yr)
        zmoy3 = zmoy2 * DBLE(iatm*1000)
c- conversion Nb_day/yr
        zmoy4 = zmoy2 * DBLE(iatm)
c- precip for veget: Sum precip > prcmin (if prcmin > 0) + Reduction fct(tpsdry)
        kk2 = 2
        if (prcmin.lt.dzero) kk2 = 1
        write(iuo+66,'(A,I2,1P4E14.6)') ' veget:k2,zmoy1,2,3,4=', kk2,
     &                               zmoy1, zmoy2, zmoy3, zmoy4
        do 270 j=1,nlon
         do 270 i=1,nlat
           temveg(i,j) = temveg(i,j) * zmoy1
           gd0veg(i,j) = gd0veg(i,j) * zmoy2
           prcveg(i,j,1) = prcveg(i,j,1) * zmoy3
           prcveg(i,j,2) = prcveg(i,j,2) * zmoy3
           tpsdry(i,j) = tpsdry(i,j) * zmoy4
 270    continue
c- Reduce Precip by 0 -> rpfdry if Nb_Day_Dry(=tpsdry) > tmxdry -> tmxdry+dtrdry
        do 280 j=1,nlon
         do 280 i=1,nlat
           prcveg(i,j,2) = prcveg(i,j,kk2) *
     &      (one-rpfdry*min(one,max(zero,(tpsdry(i,j)-tmxdry)*unsdry)))
 280    continue
 
      else
c-----------------------------------------------------------
c-Not [last time_step of the year] => return directly (no call veget.routines).
!         return
! SDubinkina, needed for restart every N month, N<12
c-Not [last time_step of the year] => go to the saving files part
c-Since it could be the last time_step of the calculation
	     goto 25
      endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

! SDubinkina, needed for restart every N month, N<12
      veg_act=.false.
      if ((iyear.eq.0).and.(mod(irunlabeld*iatm+istep,nstpyear*nyvegt).eq.0)) veg_act=.true.
      do 281 ieq=1,ieqveg

!       if (kveget.lt.0 .or. iyear*kveget.ne.0)  then
! SDubinkina, needed for restart every N month, N<12
      if (kveget.lt.0 .or. iyear*kveget.ne.0 .or. veg_act) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Call the Vegetation Routines :
c-----------------------------------------------------------------------
 
c- co2 enrichment
c-AM      co2ghg = ghg(1)
          co2ghg = patmCO2
 
c...    SPATIAL LOOP
 
      do j=1,nlon
       do i=1,nlat
        if (fracgr(i,j).gt.epss) then
c- Transfert var. for Veget. model (& modif Precip_2) :
          lat = i
          lon = j
          ave_t    = temveg(i,j)
          gdd0     = gd0veg(i,j)
          ave_pr   = prcveg(i,j,1)
          ave_pr05 = prcveg(i,j,2)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c... calculation of dynamics of carbon pools and forest ratio
 
          if (kveget.lt.0) then
            if (irad.eq.1) call ccstatR(fracgr,darea)
            call ccstat(fracgr,darea)
          else
            do k=1,kveget
              if (irad.eq.1) call ccdynR(fracgr,darea)
              call ccdyn(fracgr,darea)
            enddo
          endif
        endif
       enddo
      enddo

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      if (kveget.ne.0) then
      
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Compute albedo of land (except over Greenland & Antartica )
c      Compute bmoism en fonction de veget
c-----------------------------------------------------------------------
 
c- Antarctica Ice i: 1--> 5 [90.S -> 61.875 S]
c- Greenland Ice : continent & i >= 26 + |j-57| (inclus)
      nlatd2 = nlat / 2
      istart=6

C      bmoism(i,j)=bmoisg(i,j)
C      rs(i,j)=0.

      do i=1,nlat
         do j=1,nlon
C           bmoism(i,j)=bmoisg(i,j)
            bmoism(i,j)=0.15
C           rs(i,j)=0.
         enddo
      enddo


      do 450 j=1,nlon
       do 450 i=istart,nlat
      igroen=26+abs(j-57)

        if ( fracgr(i,j).gt.epss .and. i.lt.igroen ) then
          ii =i / nlatd2
          if (ii.eq.2) ii=1
c- transition -> Tundra & Steppe (xx=1) : T_ann entre 4 et 0 deg.C
          xx = (4. - temveg(i,j)) * 0.25
          xx = min(1.,max(0.,xx))
c--Bright sand desert (Sahara) : + albedc :
          ddc = 0.
Csah1     if (  (i.ge.19.and.i.le.21) .and. (j.le.10.or.j.ge.62)
Csah1&     .or. (i.eq.22.and.j.le.7) ) ddc = 1.
          if ( (i.ge.19.and.i.le.23).and.(j.le.12.or.j.ge.62) ) ddc=1.
          if ( i.eq.23 .and. j.ge.5 .and. j.le.12 ) ddc = 0.5
c--test Albedo funct(LAI) :
C_2       zlai = blai(i,j,1)*st(i,j)+blai(i,j,2)*sg(i,j)
C_2       xxalb = veg_albd(i,zlai,snlt(i,j),temveg(i,j),ddc)
          do 430 ns=1,4
            xxalb = st(i,j)*albet(ns)
     &            + sd(i,j)*( albed(ns) + ddc*albedc(ns) )
     &            + sg(i,j)*( albeg(ns) + xx*albegc(ns) )
           xxalbR=str(i,j)*albet(ns)
     &            + sdr(i,j)*( albed(ns) + ddc*albedc(ns) )
     &            + sgr(i,j)*( albeg(ns) + xx*albegc(ns) )
c!Attention pour calcul FR, creer bmoismR, rsR et adapter ds lbm
           bmoism(i,j)=((st(i,j)*(1-snlt(i,j)))*0.25)
     &       + (sd(i,j)*0.1)+(sg(i,j)*0.15)+(st(i,j)*snlt(i,j)*0.25)
clbm2      bmoism(i,j)=((st(i,j)*(1-snlt(i,j)))*acwt*zrt/1020.)
clbm2&            + (sd(i,j)*acwd*zrd/1020.)   
clbm2&            + (sg(i,j)*acwg*zrg/1020.)   
clbm2&            + (st(i,j)*snlt(i,j)*acwn*zrn/1020.)   
clbm1      rs(i,j)=((st(i,j)*(1-snlt(i,j)))*rst)+(sd(i,j)*rsd)
clbm1&            +(sg(i,j)*rsg)+(st(i,j)*snlt(i,j)*rsn)

c--Test Eq.Rain.Forest Low Albedo :
C_1       if (i.ge.15.and.i.le.18) xxalb = .09
          if (i.eq.16.or.i.eq.17) xxalb = .09
          if (i.eq.15.or.i.eq.18) xxalb = 0.5*xxalb + .045
          if (i.eq.16.or.i.eq.17) xxalbR = .09
          if (i.eq.15.or.i.eq.18) xxalbR = 0.5*xxalbR + .045
C_2       do 430 ns=1,4
            albland(i,j,ns0inv(ns,ii)) = xxalb
            alblandR(i,j,ns0inv(ns,ii)) = xxalbR

 430      continue
          forestfr(i,j) = st(i,j)
          forestfrR(i,j) = stR(i,j)
        endif
 450  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Ouput on ASCII file.
c-----------------------------------------------------------------------
 
      if (ieqveg.ne.1) then
      if((mod(ieq,ieqvegwr).eq.0).or.((initialization.eqv..true.).and.(ieq.eq.1))) then
         flgwrt=.TRUE.
      else
         flgwrt=.FALSE.
      endif
      endif
c- Compute & write Global & Zonal mean Var. ; Write 2.D vegetation Var.
      if (flgwrt) call veget_wr(
     &             kveget, nyvegt, nwrveg,iyr0vg, iyear, nyears,
     &             temveg,gd0veg,prcveg,tpsdry,
     &             prcmin,bmtdry, epss,fracgr,darea, 
     &             titveg)
 
 281    continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Prepare the next run (write restart) or the next year (reset).  |
c-----------------------------------------------------------------------
! SDubinkina, needed for restart every N month, N<12
      if (istep.eq.(ntotday*iatm)) then
c--Last call : Write restart files veget.init & veget.rest
        open(iveg+2,file='veget.init',status='unknown',form='unformatted')
        write(iveg+2) temveg
        write(iveg+2) gd0veg
        write(iveg+2) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
        write(iveg+2) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
        write(iveg+2) tpsdry
        write(iveg+2) iyear+iyr0vg, bmtdry, titveg
        close(iveg+2)
 
        if (kveget.ne.0) then
        open(iveg+2,file='veget.rest',status='unknown',form='unformatted')
        write(iveg+2) iyear+iyr0vg, bmtdry, titveg
        write(iveg+2) st
        write(iveg+2) sg
        write(iveg+2) sd
        write(iveg+2) snlt
        write(iveg+2) temveg  
        write(iveg+2) b1t
        write(iveg+2) b1g
        write(iveg+2) b2t
        write(iveg+2) b2g
        write(iveg+2) b3t
        write(iveg+2) b3g
        write(iveg+2) b4t
        write(iveg+2) b4g 
        write(iveg+2) b1
        write(iveg+2) b2
        write(iveg+2) b3
        write(iveg+2) b4
        write(iveg+2) anup 
        if (iscendef.eq.1) then
          write(iveg+2) st_const
          write(iveg+2) str
          write(iveg+2) sgr
          write(iveg+2) sdr
        else
          write(iveg+2) st
        endif
        
        close(iveg+2)
        write(iuo+66,'(2A,I8,A,F8.3)') ' Last call: write "veget.init"',
     &      ' & "veget.rest" ; yr=', iyear+iyr0vg, ' ; CO_2=', co2ghg
      else
        write(iuo+66,'(2A,I8,A,F8.3)') ' Last call: write "veget.init"',
     &                     ' ; yr=', iyear+iyr0vg, ' ; CO_2=', co2ghg
      endif
c-----
      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Before starting a new year, reset Cum.Var (Climato) at Zero :
      ncumvg = 0
      do j=1,nlon
       do i=1,nlat
         temveg(i,j) = 0.0
         gd0veg(i,j) = 0.0
         prcveg(i,j,1) = 0.0
         prcveg(i,j,2) = 0.0
         tpsdry(i,j) = 0.0
       enddo
      enddo
      kveget = abs(kveget)
! SDubinkina, needed for restart every N month, N<12
c--read Cum.Var (Climato) from a file :
      if(initialization.eqv..true.) then
        open(iveg+2,file='veget_accu.rest',status='old',form='unformatted',iostat=ios)
	if (ios.eq.0) then
	  read(iveg+2) ncumvg
          read(iveg+2) temveg
          read(iveg+2) prcday
          read(iveg+2) gd0veg
          read(iveg+2) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
          read(iveg+2) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
          read(iveg+2) tpsdry
	endif
        close(iveg+2)
!           write(*,*) " "
!           write(*,*) " "
!           write(*,*) " "
!           write(*,*) "Vegetation check veget, 1st"
!           write(*,*) "iostat=",ios
!           write(*,*) initialization
!           write(*,*)ist,jst,irunlabeld
!           write(*,*) kveget
!           write(*,*) iyear,nyvegt,nyears
!           write(*,*) istep,nstpyear*nyvegt
!           write(*,*) "ncumvg= ",ncumvg
!           write(*,*) "prcday(10,12)= ",prcday(10,12)
!           write(*,*) "tpsdry(10,12)= ",tpsdry(10,12)
!           write(*,*) " "
!           write(*,*) " "
!           write(*,*) " "
      endif
! SDubinkina, needed for restart every N month, N<12
c-Not [last time_step of the year] and [last time_step of the year]
   25 continue
      if ((istep.eq.(ntotday*iatm)).and.(initialization.eqv..false.)) then
c--Last call : Write restart file veget_accu.rest
!           write(*,*) " "
!           write(*,*) " "
!           write(*,*) " "
!           write(*,*) "Vegetation check veget, 2st"
!           write(*,*) initialization
!           write(*,*)ist,jst,irunlabeld
!           write(*,*) kveget
!           write(*,*) iyear,nyvegt,nyears
!           write(*,*) istep,nstpyear*nyvegt
!           write(*,*) "ncumvg= ",ncumvg
!           write(*,*) "prcday(10,12)= ",prcday(10,12)
!           write(*,*) "tpsdry(10,12)= ",tpsdry(10,12)
!           write(*,*) " "
!           write(*,*) " "
!           write(*,*) " "
        open(iveg+2,file='veget_accu.rest',status='unknown',form='unformatted')
        write(iveg+2) ncumvg
        write(iveg+2) temveg
	write(iveg+2) prcday
        write(iveg+2) gd0veg
        write(iveg+2) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
        write(iveg+2) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
        write(iveg+2) tpsdry
        write(iveg+2) iyear+iyr0vg, bmtdry, titveg
        close(iveg+2)

      else
c-Not [last time_step of the run] => return without saving veget_accu.rest
	return
      endif


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine veget -
      end
