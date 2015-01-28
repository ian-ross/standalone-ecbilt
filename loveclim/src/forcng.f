












      subroutine forcng(ja,xjour)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--This routine reads the and atmospheric forcing (NCAR)
c  needed by the model. It computes surface fluxes and inter-
c  polates every quantity in time (daily) by means of a cubic spline.
c  Derniere modification 30/05/95
c  difference entre version glace et version CLIO:
c  include 'oceacn.com',tcn-scal
c-----
c  modif : 04/10/99

 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'trace.com'
 
c- nn99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
      common / mchd99 / nn99
 
      parameter(nsamt=12)
 
c iwind.eq.1, use already computed geostrophic winds.
c iwind.ne.1, compute wind stress from derived geostrophic winds.
c uwind     Wind velocity in the x direction
c vwind     Wind velocity in the y direction
c
        parameter(iwind=1)
c
        real*4 rtt(imax,jmax,2,2),rtr(imax,jmax,2,2),
     &         rpr(imax,jmax,2,2),rdw(imax,jmax,2,2,2),
     &         rrn(imax,jmax,2,2),rcd(imax,jmax,2,2),
     &         rhr(imax,jmax,2,2,2),rvv(imax,jmax,2,2),
     &         rru(imax,jmax,2,2),rdv(imax,jmax,2,2)
c
        logical flgveg,flgicb,flgisma,flgismg
        dimension rhcg2(nsamt,2),theta(nsamt,2)
        dimension psl(imax,jmax),ugw(imax,jmax),vgw(imax,jmax)
C       dimension cmartc(nsamt,2),cmarsh(nsamt,imax,jmax)
        dimension bunker(19),budyko(19)
        dimension trans(imax,jmax),tenhrx(imax,jmax),tenhry(imax,jmax)
        dimension vabqec(imax,jmax), zmodfo(imax,jmax)
c-uwind et vwind en common pour les icebergs
c       dimension uwind(imax,jmax),vwind(imax,jmax)
c
        common /ec_coupl/ flgveg,flgicb,flgisma,flgismg
        common/splno3/ rtt,rtr,rpr,rdw,rrn,rcd,rhr,rvv,rru,rdv
        common/splno4/ amtt(imax,jmax),amtr(imax,jmax),ampr(imax,jmax),
     &                 amwx(imax,jmax),amwy(imax,jmax),amrn(imax,jmax),
     &                 amcd(imax,jmax),amhrx(imax,jmax),amhry(imax,jmax)
     &                ,amvv(imax,jmax),amru(imax,jmax),amdv(imax,jmax)
        common/splno5/ jnddt,mit,mft
        common/splno6/ cmarsh,bunker,budyko
        common/tranz/ trans,zmodfo
c
        data rhcg2  /0.00085,0.00085,0.000850,0.000933,0.001016,0.00110,
     &               0.00110,0.00110,0.001016,0.000933,0.000850,0.00085,
     &               0.00110,0.00110,0.001016,0.000933,0.000850,0.00085,
     &               0.00085,0.00085,0.000850,0.000933,0.001016,0.00110/
        data theta  /33.0,33.0,33.0,29.6,26.3,23.0,
     &               23.0,23.0,26.3,29.6,33.0,33.0,
     &               23.0,23.0,26.3,29.6,33.0,33.0,
     &               33.0,33.0,33.0,29.6,26.3,23.0/
c
c  MARSHUNOVA's cp coefficient in the central Arctic:
c               source: Shine and Crane, 1984.
c
C       data cmartc /0.310,0.285,0.260,0.235,0.210,0.185,
C    &                     0.160,0.185,0.210,0.235,0.260,0.285,
C    &                     0.160,0.185,0.210,0.235,0.260,0.285,
C    &                     0.310,0.285,0.260,0.235,0.210,0.185/
c
c               source: Doronin, 1969.
c
C       data cmartc /0.300,0.300,0.300,0.280,0.270,0.240,
C    &                     0.220,0.230,0.270,0.290,0.300,0.300,
C    &                     0.220,0.230,0.270,0.290,0.300,0.300,
C    &                     0.300,0.300,0.300,0.280,0.270,0.240/
c
c  BUNKER's coefficient (cloudiness effect on LW radiation):
c
C       data bunker /1.00,1.00,0.95,0.90,0.86,0.81,0.75,0.70,0.62,0.60,
C    &                    0.62,0.70,0.75,0.81,0.86,0.90,0.95,1.00,1.00/
c
c  BUDYKO's coefficient (cloudiness effect on LW radiation):
c
C       data budyko /0.82,0.82,0.82,0.82,0.80,0.78,0.76,0.74,0.72,0.70,
C    &               0.68,0.65,0.63,0.61,0.59,0.57,0.55,0.52,0.50/
C       data budyko /1.00,0.98,0.95,0.92,0.89,0.86,0.83,0.80,0.78,0.75,
C    &               0.72,0.69,0.67,0.64,0.61,0.58,0.56,0.53,0.50/
        data budyko /1.00,0.95,0.89,0.83,0.78,0.72,0.67,0.61,0.56,0.50,
     &               0.57,0.62,0.68,0.75,0.90,1.00,1.00,1.00,1.00/
 
c
c  Begin.
c
      dpi   = 2.0*pi
      dtt=yeaday/float(nsamt)
      dtts6=dtt/6.0
C     jour  = int(xjour)
C     datet = (dble((ja-1)*yeaday+jour)-0.5)/dtt
      datet = ((ja-1)*yeaday+xjour-1.0)/dtt
      ttbt  = datet-int(datet)
      ttat  = 1.0-ttbt
      jndt  = int(datet)
      jm    = mod(jndt,nsamt)+1
      jmm1  = jm-1+nsamt*(1/jm)
      jmp1  = jm+1-nsamt*(jm/nsamt)
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c       READ NCAR DATA.                                                |
c-----------------------------------------------------------------------
c
      if(numit.eq.nstart) then
        if (flgicb) then
          open(92,file='wind.int.dat',
     &        form='unformatted')
          open(93,file='wind.int.d.dat',
     &        form='unformatted')
        endif
c
C       open(71,file='/u14/grpastr/hgs/data/ncar/ncar.int.dat',
C    &        form='unformatted')
C       open(72,file='/u14/grpastr/hgs/data/ncar/ncar.int.d.dat',
C    &        form='unformatted')
        open(71,file='pres.int.dat',
     &  form='unformatted')
        open(72,file='pres.int.d.dat',
     &  form='unformatted')
        open(73,file='prcp.int.dat',
     &  form='unformatted')
        open(74,file='prcp.int.d.dat',
     &  form='unformatted')
C       open(73,file='/u14/grpastr/hgs/data/precip/pleg.int.dat',
C    &  form='unformatted')
C       open(74,file='/u14/grpastr/hgs/data/precip/pleg.int.d.dat',
C    &  form='unformatted')
        open(75,file='cldsa.int.dat',
     &  form='unformatted')
        open(76,file='cldsa.int.d.dat',
     &  form='unformatted')
C       open(77,file='tenhr.int.dat',
C    &  form='unformatted')
C       open(78,file='tenhr.int.d.dat',
C    &  form='unformatted')
        open(77,file='ten5.int.dat',
     &  form='unformatted')
        open(78,file='ten5.int.d.dat',
     &  form='unformatted')
C       open(79,file='runoff.anu2.dat',
C    &  form='unformatted')
        open(80,file='vvent.int.dat',
     &  form='unformatted')
        open(81,file='vvent.int.d.dat',
     &  form='unformatted')
        open(82,file='temp.int.dat',
     &  form='unformatted')
        open(83,file='temp.int.d.dat',
     &  form='unformatted')
C       open(84,file='/u14/grpastr/hgs/data/ncar/tdew.int.dat',
C    &  form='unformatted')
C       open(85,file='/u14/grpastr/hgs/data/ncar/tdew.int.d.dat',
C    &  form='unformatted')
        open(84,file='humid.int.dat',
     &  form='unformatted')
        open(85,file='humid.int.d.dat',
     &  form='unformatted')
        open(86,file='runo.int.dat',
     &  form='unformatted')
        open(87,file='runo.int.d.dat',
     &  form='unformatted')
        open(88,file='ust.int.dat',
     &  form='unformatted')
        open(89,file='ust.int.d.dat',
     &  form='unformatted')
 
        rewind 71
        rewind 72
        rewind 73
        rewind 74
        rewind 75
        rewind 76
        rewind 77
        rewind 78
        rewind 79
        rewind 80
        rewind 81
        rewind 82
        rewind 83
        rewind 84
        rewind 85
        rewind 86
        rewind 87
        rewind 88
        rewind 89
        rewind 92
        rewind 93
c
        jjt=mod(jndt,nsamt)
        do 2 j=1,jmax
        read(71) (rpr(i,j,1,1),i=1,imax)
C       read(71) (rtt(i,j,1,1),i=1,imax)
C       read(71) (rtr(i,j,1,1),i=1,imax)
C       read(71) (rdw(i,j,1,1,1),i=1,imax)
C       read(71) (rdw(i,j,2,1,1),i=1,imax)
        read(72) (rpr(i,j,2,1),i=1,imax)
C       read(72) (rtt(i,j,2,1),i=1,imax)
C       read(72) (rtr(i,j,2,1),i=1,imax)
C       read(72) (rdw(i,j,1,2,1),i=1,imax)
C       read(72) (rdw(i,j,2,2,1),i=1,imax)
        if (flgicb) then
          read(92) (rdw(i,j,1,1,1),i=1,imax)
          read(92) (rdw(i,j,2,1,1),i=1,imax)
          read(93) (rdw(i,j,1,2,1),i=1,imax)
          read(93) (rdw(i,j,2,2,1),i=1,imax)
        endif
        read(73) (rrn(i,j,1,1),i=1,imax)
        read(74) (rrn(i,j,2,1),i=1,imax)
        read(75) (rcd(i,j,1,1),i=1,imax)
        read(76) (rcd(i,j,2,1),i=1,imax)
        read(77) (rhr(i,j,1,1,1),i=1,imax)
        read(77) (rhr(i,j,2,1,1),i=1,imax)
        read(78) (rhr(i,j,1,2,1),i=1,imax)
        read(78) (rhr(i,j,2,2,1),i=1,imax)
        read(80) (rvv(i,j,1,1),i=1,imax)
        read(81) (rvv(i,j,1,2),i=1,imax)
        read(82) (rtt(i,j,1,1),i=1,imax)
        read(83) (rtt(i,j,2,1),i=1,imax)
        read(84) (rtr(i,j,1,1),i=1,imax)
        read(85) (rtr(i,j,2,1),i=1,imax)
        read(86) (rru(i,j,1,1),i=1,imax)
        read(87) (rru(i,j,2,1),i=1,imax)
        read(88) (rdv(i,j,1,1),i=1,imax)
        read(89) (rdv(i,j,2,1),i=1,imax)
        do 1 i=1,imax
        amtt(i,j)=rtt(i,j,1,1)
        amtr(i,j)=rtr(i,j,1,1)
        ampr(i,j)=rpr(i,j,1,1)
C       amwx(i,j)=rdw(i,j,1,1,1)
C       amwy(i,j)=rdw(i,j,2,1,1)
        if (flgicb) then
          amwx(i,j)=rdw(i,j,1,1,1)
          amwy(i,j)=rdw(i,j,2,1,1)
        endif
        amrn(i,j)=rrn(i,j,1,1)
        amcd(i,j)=rcd(i,j,1,1)
        amhrx(i,j)=rhr(i,j,1,1,1)
        amhry(i,j)=rhr(i,j,2,1,1)
        amvv(i,j)=rvv(i,j,1,1)
        amru(i,j)=rru(i,j,1,1)
        amdv(i,j)=rdv(i,j,1,1)
 1      continue
 2        continue
        do 4 l=1,jjt
        do 3 j=1,jmax
        read(71) (rpr(i,j,1,1),i=1,imax)
C       read(71) (rtt(i,j,1,1),i=1,imax)
C       read(71) (rtr(i,j,1,1),i=1,imax)
C       read(71) (rdw(i,j,1,1,1),i=1,imax)
C       read(71) (rdw(i,j,2,1,1),i=1,imax)
        if (flgicb) then
          read(92) (rdw(i,j,1,1,1),i=1,imax)
          read(92) (rdw(i,j,2,1,1),i=1,imax)
          read(93) (rdw(i,j,1,2,1),i=1,imax)
          read(93) (rdw(i,j,2,2,1),i=1,imax)
        endif
        read(72) (rpr(i,j,2,1),i=1,imax)
C       read(72) (rtt(i,j,2,1),i=1,imax)
C       read(72) (rtr(i,j,2,1),i=1,imax)
C       read(72) (rdw(i,j,1,2,1),i=1,imax)
C       read(72) (rdw(i,j,2,2,1),i=1,imax)
        read(73) (rrn(i,j,1,1),i=1,imax)
        read(74) (rrn(i,j,2,1),i=1,imax)
        read(75) (rcd(i,j,1,1),i=1,imax)
        read(76) (rcd(i,j,2,1),i=1,imax)
        read(77) (rhr(i,j,1,1,1),i=1,imax)
        read(77) (rhr(i,j,2,1,1),i=1,imax)
        read(78) (rhr(i,j,1,2,1),i=1,imax)
        read(78) (rhr(i,j,2,2,1),i=1,imax)
        read(80) (rvv(i,j,1,1),i=1,imax)
        read(81) (rvv(i,j,2,1),i=1,imax)
        read(82) (rtt(i,j,1,1),i=1,imax)
        read(83) (rtt(i,j,2,1),i=1,imax)
        read(84) (rtr(i,j,1,1),i=1,imax)
        read(85) (rtr(i,j,2,1),i=1,imax)
        read(86) (rru(i,j,1,1),i=1,imax)
        read(87) (rru(i,j,2,1),i=1,imax)
        read(88) (rdv(i,j,1,1),i=1,imax)
        read(89) (rdv(i,j,2,1),i=1,imax)
 3      continue
 4      continue
C       write(95,*) ampr(25,32)
C       write(95,*) amtt(20,40),amtr(20,40),ampr(20,40),amwx(20,40),
C    &              amwy(20,40),amrn(20,40),amcd(20,40),amhrx(20,40),
C    &              amhry(20,40),runoff(20,40)
c modif precip
        nzmdif = 0
        open(91,file='modif.pre',status='old',err=5)
        read(91,*)
        read(91,*)
        read(91,*)
        read(91,*) nzmdif
        do n=1,nzmdif
          read(91,*)
          read(91,*) i1,i2
          read(91,*) j1,j2
          read(91,*) zmopr
          zmopr=zmopr/(100.*365.0*86400.0)
          do j=j1-1,j2+1
           do i=i1-1,i2+1
            zmodfo(i,j)=zmopr/2.0
           enddo
          enddo
          do j=j1,j2
           do i=i1,i2
            zmodfo(i,j)=zmopr
           enddo
          enddo
        enddo
        close(unit=91)
 5      continue
c
c  Data is given for 66 years, and for the SH (1) and NH (2).
c  Note that year 1989 data is my own extrapolation to ensure last
c    few months of 1988 are forced properly.  STOP PRESS: Now have 0
c    data going until 1995 based on Elkins et al., 1993 Nature Vol 364
c    pp 780-783.
 
c=======================================================================
c    Read in atmospheric 0 histories
c=======================================================================
Ccfc   if (nn99.eq.2) write(99,977)
Ccfc   open(91,form='formatted',file='CFC_atmos_data')
Ccfc   do 7 nnyear=1930,1995
Ccfc       read(91,979) mmyr, cfc11(nnyear-1929,2),
Ccfc &             cfc12(nnyear-1929,2), xcfcn,
Ccfc &             cfc11(nnyear-1929,1), cfc12(nnyear-1929,1), xcfcs
Ccfc       if (nn99.eq.2) write(99,978) mmyr, cfc11(nnyear-1929,2),
Ccfc & cfc12(nnyear-1929,2),cfc11(nnyear-1929,1),cfc12(nnyear-1929,1)
7      continue
Ccfc   close(unit=91)
Ccfc   if (nn99.eq.2) write(99,*)
979    format(i4, 2f8.2, f7.3, 2f8.2, f7.3)
978    format(i5, 4f11.2)
977    format('Year:   CFC-11(N): CFC-12(N): CFC-11(S): CFC-12(S):')
 
c
c   Coefficients for geostrophic derivation of wind stress.
c   (Overland & Colony, 1994; Hibler, 1972).
c
        do k=1,nsamt
          theta(k,1) = radian*theta(k,1)
          theta(k,2) = radian*theta(k,2)
        enddo
c
c  Coefficient for cloudiness effect on LW radiation.
c
C       do j=1,jmax
C         do i=1,imax
C           alat    = asin(covrai(i,j))/radian
C           clat    = (90.0-alat)/10.0
C           indx    = 1+int(clat)
C           ihm     = 1+max(0,sign(1,indx-10))
C           indxp1  = indx+1
C           ihmp1   = 1+max(0,sign(1,indxp1-10))
C           dl      = clat-int(clat)
C           dr      = 1.0-dl
C           do k=1,nsamt
C             cmarsh(k,i,j) = dr*bunker(indx)*cmartc(k,ihm)
C    &                       +dl*bunker(indxp1)*cmartc(k,ihmp1)
C             cmarsh(k,i,j) = dr*cmartc(k,ihm)+dl*cmartc(k,ihmp1)
C           enddo
C         enddo
C       enddo
 
c      coef for the transition between world ocean and the poles
       do 11 j=1,jmax
         do 11 i=1,imax
C        do 11 i=is1(j),is2(j)
           xlati=asin(covrai(i,j))/radian
           if (xlati.lt.(-65.0)) then
             trans(i,j)= 1.0
           else
             if (xlati.lt.(-50.0)) then
               trans(i,j)=1.0 -(xlati+65.0)/15.0
             else
               if (xlati.lt.55.0) then
                 trans(i,j)=0.0
               else
                 if (xlati.lt.70.0) then
                   trans(i,j)=0.0+(xlati-55)/15.0
                 else
                   if (xlati.lt.90.0) then
                     trans(i,j)=1.0
                   else
                     write(iuo+66,*) 'arret dans latitude',i,j,xlati
                   endif
                 endif
               endif
             endif
           endif
11     continue
c
        mit=2
        mft=1
        jnddt=jndt
c
      endif
c
      do 15 l=jnddt,jndt
      mit=3-mit
      mft=3-mft
      if (mod(jndt+1,nsamt).eq.0) then
        rewind 71
        rewind 72
        rewind 73
        rewind 74
        rewind 75
        rewind 76
        rewind 77
        rewind 78
        rewind 80
        rewind 81
        rewind 82
        rewind 83
        rewind 84
        rewind 85
        rewind 86
        rewind 87
        rewind 88
        rewind 89
        rewind 92
        rewind 93 
      endif
      do 25 j=1,jmax
        read(71) (rpr(i,j,1,mft),i=1,imax)
C       read(71) (rtt(i,j,1,mft),i=1,imax)
C       read(71) (rtr(i,j,1,mft),i=1,imax)
C       read(71) (rdw(i,j,1,1,mft),i=1,imax)
C       read(71) (rdw(i,j,2,1,mft),i=1,imax)
        if (flgicb) then
          read(92) (rdw(i,j,1,1,mft),i=1,imax)
          read(92) (rdw(i,j,2,1,mft),i=1,imax)
          read(93) (rdw(i,j,1,2,mft),i=1,imax)
          read(93) (rdw(i,j,2,2,mft),i=1,imax)
        endif
        read(72) (rpr(i,j,2,mft),i=1,imax)
C       read(72) (rtt(i,j,2,mft),i=1,imax)
C       read(72) (rtr(i,j,2,mft),i=1,imax)
C       read(72) (rdw(i,j,1,2,mft),i=1,imax)
C       read(72) (rdw(i,j,2,2,mft),i=1,imax)
        read(73) (rrn(i,j,1,mft),i=1,imax)
        read(74) (rrn(i,j,2,mft),i=1,imax)
        read(75) (rcd(i,j,1,mft),i=1,imax)
        read(76) (rcd(i,j,2,mft),i=1,imax)
        read(77) (rhr(i,j,1,1,mft),i=1,imax)
        read(77) (rhr(i,j,2,1,mft),i=1,imax)
        read(78) (rhr(i,j,1,2,mft),i=1,imax)
        read(78) (rhr(i,j,2,2,mft),i=1,imax)
        read(80) (rvv(i,j,1,mft),i=1,imax)
        read(81) (rvv(i,j,2,mft),i=1,imax)
        read(82) (rtt(i,j,1,mft),i=1,imax)
        read(83) (rtt(i,j,2,mft),i=1,imax)
        read(84) (rtr(i,j,1,mft),i=1,imax)
        read(85) (rtr(i,j,2,mft),i=1,imax)
        read(86) (rru(i,j,1,mft),i=1,imax)
        read(87) (rru(i,j,2,mft),i=1,imax)
        read(88) (rdv(i,j,1,mft),i=1,imax)
        read(89) (rdv(i,j,2,mft),i=1,imax)
 25   continue
      jnddt=jndt+1
 15   continue
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c       INTERPOLATE DATA.                                              |
c-----------------------------------------------------------------------
c
      do 35 j=1,jmax
        ihm = max(0,sign(1,j-jeq))
        hem  = real(ihm)
      do 30 i=1,imax
C       alat = asin(covrai(i,j))/radian
C       indx = 1+int((90.0-abs(alat))/5.0)
        alat = asin(covrai(i,j))/radian
        clat = (95.0-alat)/10.0
        indx = 1+int(clat)
 
C     fwat(i,j)=max(zero,((rrn(i,j,1,mft)-rrn(i,j,1,mit))/dtt-
C    &         ((3.0*ttat*ttat-1.0)*rrn(i,j,2,mit)-
C    &          (3.0*ttbt*ttbt-1.0)*rrn(i,j,2,mft))*dtts6+
C    &               amrn(i,j))/dtt)
      fwat(i,j)=max(zero,(rrn(i,j,1,mft)-rrn(i,j,1,mit))/dtt-
     &         ((3.0*ttat*ttat-1.0)*rrn(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rrn(i,j,2,mft))*dtts6+
     &         amrn(i,j))
      cloud(i,j)=max(zero,(rcd(i,j,1,mft)-rcd(i,j,1,mit))/dtt-
     &         ((3.0*ttat*ttat-1.0)*rcd(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rcd(i,j,2,mft))*dtts6+
     &         amcd(i,j))
      tabq(i,j)=(rtt(i,j,1,mft)-rtt(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rtt(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rtt(i,j,2,mft))*dtts6+
     &          amtt(i,j)+273.15
C     tdew(i,j)=(rtr(i,j,1,mft)-rtr(i,j,1,mit))/dtt-
C    &          ((3.0*ttat*ttat-1.0)*rtr(i,j,2,mit)-
C    &          (3.0*ttbt*ttbt-1.0)*rtr(i,j,2,mft))*dtts6+
C    &                amtr(i,j)+273.15
      tdew(i,j)=(rtr(i,j,1,mft)-rtr(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rtr(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rtr(i,j,2,mft))*dtts6+
     &          amtr(i,j)
      psl(i,j)=((rpr(i,j,1,mft)-rpr(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rpr(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rpr(i,j,2,mft))*dtts6+
     &          ampr(i,j))*100.0
C     psbq(i,j)=hem*101400.0+(1.0-hem)*98800.0
      psbq(i,j)=psl(i,j)
      vabqec(i,j)=(rvv(i,j,1,mft)-rvv(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rvv(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rvv(i,j,2,mft))*dtts6+
     &          amvv(i,j)
      runoff(i,j)=(rru(i,j,1,mft)-rru(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rru(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rru(i,j,2,mft))*dtts6+
     &          amru(i,j)
      sdvt(i,j)=((rdv(i,j,1,mft)-rdv(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rdv(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rdv(i,j,2,mft))*dtts6+
     &          amdv(i,j))
c modification du forcing P-E
C     if (j.eq.6) then
C     write(112,*) i,runoff(i,6),zmodfo(i,6)
C     endif
      runoff(i,j)=runoff(i,j)+zmodfo(i,j)
c
c
c  Longwave radiation:
c----------------------
c
c  vapor pressures (in millibars).
c
C     evg=611.0*10.0**(9.5*(tdew(i,j)-273.16)/(tdew(i,j)-7.660))*0.01
C     evo=611.0*10.0**(7.5*(tdew(i,j)-273.16)/(tdew(i,j)-35.86))*0.01
C     es1=611.0*10.0**(7.5*(tabq(i,j)-273.16)/(tabq(i,j)-7.660))*0.01
C     es2=611.0*10.0**(7.5*(tabq(i,j)-273.16)/(tabq(i,j)-35.86))*0.01
C     evg=tdew(i,j)*es1
C     evo=tdew(i,j)*es2
      es =611.0*exp(min(sign(17.269*one,tabq(i,j)-too),sign(21.875*one,
     &    tabq(i,j)-too))*abs(tabq(i,j)-too)/(tabq(i,j)-35.86+
     &    max(zero,sign(28.2*one,too-tabq(i,j)))))
      evg=tdew(i,j)*es*0.01
      evo=tdew(i,j)*es*0.01
c
c  1. IDSO and JACKSON (all latitudes). (from Parkinson and Washington).
c
c     ratbqg(i,j)=(1.0+cloud(i,j)*0.275)*
c    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
c    &            (1.0-0.261*exp(-7.77e-04*(273.0-tabq(i,j))
c    &             *(273.0-tabq(i,j))))
c     ratbqo(i,j)=(1.0+cloud(i,j)*0.275)*
c    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
c    &            (1.0-0.261*exp(-7.77e-04*(273.0-tabq(i,j))
c    &             *(273.0-tabq(i,j))))
c
c  1.b IDSO, 1981 (all latitudes).
c
c     ratbqg(i,j)=(1.0+cloud(i,j)*0.275)*
c    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
c    &            (0.70+5.95e-05*evg*exp(1500.0/tabq(i,j)))
c     ratbqo(i,j)=(1.0+cloud(i,j)*0.275)*
c    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
c    &            (0.70+5.95e-05*evo*exp(1500.0/tabq(i,j)))
c
c  2. MAYKUT AND CHURCH (Arctic).
c
c     ratbqg(i,j)=0.7855*(1.0+0.2232*cloud(i,j)**2.75)*
c    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))
c     ratbqo(i,j)=0.7855*(1.0+0.2232*cloud(i,j)**2.75)*
c    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))
c
c  3. MARSHUNOVA (Arctic).
c
C     cmar        = 0.5*((cmarsh(jmm1,i,j)+cmarsh(jm,i,j))*ttat+
C    &                   (cmarsh(jm,i,j)+cmarsh(jmp1,i,j))*ttbt)
C     ratbqg(i,j) = stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
C    &              (0.67+0.05*sqrt(evg))*(1.0+cloud(i,j)*cmar)
C     ratbqo(i,j) = stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
C    &              (0.67+0.05*sqrt(evo))*(1.0+cloud(i,j)*cmar)
c
c  3.b MARSHUNOVA with cloud correction depending on latitude.
c
c     ratbqg(i,j)=stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
c    &            (1.0-(0.33-0.05*sqrt(evg))*
c    &                 (1.0-budyko(indx)*cloud(i,j)*cloud(i,j)))
c     ratbqo(i,j)=stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
c    &            (1.0-(0.33-0.05*sqrt(evo))*
c    &                 (1.0-budyko(indx)*cloud(i,j)*cloud(i,j)))
c
c  4. BERLIAND (all latitudes).
c
C     ratbqg(i,j)=stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
C    &            (1.0-(0.39-0.05*sqrt(evg))*(1.0-budyko(indx)*
C    &              cloud(i,j)*cloud(i,j)))
C     ratbqo(i,j)=stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
C    &            (1.0-(0.39-0.05*sqrt(evo))*(1.0-budyko(indx)
C    &             *cloud(i,j)*cloud(i,j)))
      ratbqg(i,j)=-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
     &            (0.39-0.05*sqrt(evg))*(1.0-budyko(indx)*
     &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
     &             tabq(i,j)*tabq(i,j)*(ts(i,j)-tabq(i,j))+
     &             stefan*ts(i,j)*ts(i,j)*ts(i,j)*ts(i,j)
C     ratbqo(i,j)=-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
C    &            (0.39-0.05*sqrt(evo))*(1.0-budyko(indx)*
C    &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
C    &             tabq(i,j)*tabq(i,j)*(tcn(i,j,1)-tabq(i,j))+
C    &             stefan*tcn(i,j,1)*tcn(i,j,1)*tcn(i,j,1)*tcn(i,j,1)
C     ratbqo(i,j)=-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
C    &            (0.39-0.05*sqrt(evo))*(1.0-budyko(indx)*
C    &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
C    &             tabq(i,j)*tabq(i,j)*(scal(i,j,ks2,1)-tabq(i,j))+
C    &             stefan*scal(i,j,ks2,1)*scal(i,j,ks2,1)*
C    &             scal(i,j,ks2,1)*scal(i,j,ks2,1)
      ratbqo(i,j)=zemise*(-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)
     &             *tabq(i,j)*(0.39-0.05*sqrt(evo))*(1.0-budyko(indx)*
     &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
     &             tabq(i,j)*tabq(i,j)*(scal(i,j,ks2,1)-tabq(i,j)) )
C  Comparison with the param of Andreas
C     if (i.eq.50.and.j.lt.12) then
C      toto=(0.601+5.95e-5*evo*exp(1500./tabq(i,j)))*
C    &        stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)
C      titi=stefan*scal(i,j,ks2,1)*scal(i,j,ks2,1)*
C    &             scal(i,j,ks2,1)*scal(i,j,ks2,1)
C      fnet=(toto-titi)*(1-.11*8.0*cloud(i,j))
C      write(iuo+66,*) 'IR',j,ratbqo(i,j),fnet
C     endif
      uwind(i,j)=(rdw(i,j,1,1,mft)-rdw(i,j,1,1,mit))/dtt-
     &            ((3.0*ttat*ttat-1.0)*rdw(i,j,1,2,mit)-
     &            (3.0*ttbt*ttbt-1.0)*rdw(i,j,1,2,mft))*dtts6+
     &            amwx(i,j)
      write(*,*)'forcng',uwind(i,j)
      vwind(i,j)=(rdw(i,j,2,1,mft)-rdw(i,j,2,1,mit))/dtt-
     &            ((3.0*ttat*ttat-1.0)*rdw(i,j,2,2,mit)-
     &            (3.0*ttbt*ttbt-1.0)*rdw(i,j,2,2,mft))*dtts6+
     &            amwy(i,j)
 
      tenhrx(i,j)=(rhr(i,j,1,1,mft)-rhr(i,j,1,1,mit))/dtt-
     &             ((3.0*ttat*ttat-1.0)*rhr(i,j,1,2,mit)-
     &             (3.0*ttbt*ttbt-1.0)*rhr(i,j,1,2,mft))*dtts6+
     &             amhrx(i,j)
      tenhry(i,j)=(rhr(i,j,2,1,mft)-rhr(i,j,2,1,mit))/dtt-
     &             ((3.0*ttat*ttat-1.0)*rhr(i,j,2,2,mit)-
     &             (3.0*ttbt*ttbt-1.0)*rhr(i,j,2,2,mft))*dtts6+
     &             amhry(i,j)
 
 30   continue
 35   continue
c
c Set up 0-11 and 0-12 atmospheric concentrations as a function
c   of hemisphere and time (linearly interpolate yearly data).
c 1st some year variables and time weights for atmospheric cfc conc'ns:
 
Ccfc    cfcyear = 1929.0 + float(ja)
Ccfc    intyear = ja
 
Ccfc    if(cfcyear.eq.1930) then
Ccfc       atmcfc11(1) = 0.0
Ccfc       atmcfc11(2) = 0.0
Ccfc       atmcfc12(1) = 0.0
Ccfc       atmcfc12(2) = 0.0
Ccfc       goto 41
Ccfc    endif
 
Ccfc    if(cfcyear.ge.1931 .and. cfcyear.le.1978) then
c During and before 1978 data is for 31st December:
Ccfc       cfcweighta = (xjour)/yeaday
Ccfc       cfcweightb = 1.0 - cfcweighta
Ccfc       goto 51
Ccfc    endif
 
Ccfc    if(cfcyear.eq.1979 .and. xjour.le.180.0) then
c  In first part of 1979, interpolation is with only 6 months gap!
Ccfc       cfcweighta =  (xjour)/182.5
Ccfc       cfcweightb = 1.0 - cfcweighta
Ccfc       goto 51
Ccfc    endif
 
Ccfc    if(cfcyear.ge.1979) then
c  Later than July 1st 1979, data is 12 monthly centred on July 1:
Ccfc        if(xjour.le.183.0) then
Ccfc            cfcweighta = (xjour + 181.0)/yeaday
Ccfc            cfcweightb = 1.0 - cfcweighta
Ccfc        endif
Ccfc        if(xjour.gt.183.0) then
c  For months after July, interpolation is between year *ahead* and tyear:
Ccfc            intyear = intyear + 1
Ccfc            cfcweightb = (548.0 - xjour)/yeaday
Ccfc            cfcweighta = 1.0 - cfcweightb
Ccfc        endif
Ccfc    endif
 
51      continue
 
c Atmospheric 0 computed as a function of time and hemisphere:
 
Ccfc       ihemcfc = 1
Ccfc       atmcfc11(ihemcfc) =
Ccfc $          cfc11(intyear,ihemcfc)*cfcweighta +
Ccfc $          cfc11(intyear-1,ihemcfc)*cfcweightb
Ccfc       atmcfc12(ihemcfc) =
Ccfc $          cfc12(intyear,ihemcfc)*cfcweighta +
Ccfc $          cfc12(intyear-1,ihemcfc)*cfcweightb
 
Ccfc       ihemcfc = 2
Ccfc       atmcfc11(ihemcfc) =
Ccfc $          cfc11(intyear,ihemcfc)*cfcweighta +
Ccfc $          cfc11(intyear-1,ihemcfc)*cfcweightb
Ccfc       atmcfc12(ihemcfc) =
Ccfc $          cfc12(intyear,ihemcfc)*cfcweighta +
Ccfc $          cfc12(intyear-1,ihemcfc)*cfcweightb
 
c Test output some 0 weightings to check time interpolation:
c       if(j.eq.20 .and. prntsi .and. eots) then
c          write(stdout,888) intyear, ihemcfc,
c     $ cfc12(intyear,ihemcfc), cfc12(intyear-1,ihemcfc),
c     $ cfcweighta, cfcweightb,tmonth, atmcfc12
c888       format('Year+1=',i3,' hem=',i2,'cfc12(t1)=',f6.2,
c     $ 'cfc12(t0)=',f6.2,'weight to t1=',f4.2,'weight to t0=',
c     $ f4.2,'since tmonth=',f4.1,'==> CFC-12=',f6.2)
c       endif
 
41    continue
 
c  Set some 0 variables for printing:
C       if(j.eq.15 .and. prntsi .and. eots) then
C          cfc11atm = atmcfc11
C          cfc12atm = atmcfc12
C          write(stdout,889) atmcfc11, atmcfc12
C889        format('Calculated values for Southern Hemisphere CFCs:',
C    $ ' CFC-11 =',f8.3,'  CFC-12=',f8.3)
C       endif
C      if ((xjour/50.0).eq.int(xjour/50.)) then
C      write(117,*) ja,xjour,atmcfc11
C      write(118,*) ja,xjour,atmcfc12
C      endif
 
      if (iwind.eq.1) go to 400
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3) Wind velocity and wind stress (OVERLAND and COLONY, 1994).       |
c-----------------------------------------------------------------------
c
c--3.1. Calculate geostrophic wind at the corners of the grid squares.
c---------------------------------------------------------------------
c
      rhlim=2.0*omega*sin(15.0*radian)
      do 330 j=2,jmax
          jm1 = j-1
        do 320 i=1,imax
          im1      =  (i-1)+(imax-2)*(1/i)
          rhoa     =  psbq(i,j)/(287.04*tabq(i,j))
C         rhoafn   =  rhoa*zfn(i,j)+1d-18
          rhoafn   =  rhoa*sign(one,zfn(i,j))*max(abs(zfn(i,j)),rhlim)
c         write(89,*) i,j,rhoa,zfn(i,j)
          ugw(i,j) = -(-(alambd(i,j,1,1,2,1)+
     &                       alambd(i,j,1,2,2,1))*psl(i,jm1)-
     &                 (alambd(i,j,1,1,1,1)
     &                      +alambd(i,j,1,2,1,1))*psl(im1,jm1)+
     &                 (alambd(i,j,1,1,2,2)
     &                      -alambd(i,j,1,2,2,2))*psl(i,j)+
     &                 (alambd(i,j,1,1,1,2)
     &                      -alambd(i,j,1,2,1,2))*psl(im1,j))/
     &                rhoafn
 
          vgw(i,j) =  ((alambd(i,j,2,2,2,1)
     &                     -alambd(i,j,2,1,2,1))*psl(i,jm1)+
     &                 (alambd(i,j,2,2,2,2)
     &                     -alambd(i,j,2,1,2,2))*psl(i,j)-
     &                 (alambd(i,j,2,2,1,1)
     &                     +alambd(i,j,2,1,1,1))*psl(im1,jm1)-
     &                 (alambd(i,j,2,2,1,2)
     &                     +alambd(i,j,2,1,1,2))*psl(im1,j))/
     &                rhoafn
320     continue
330   continue
c  South pole: never used.
c
      do 360 i=1,imax
         ugw(i,1) = 0.0
         vgw(i,1) = 0.0
360   continue
c
c  Equator.
c
      njh = jeq-1
      do 370 i=1,imax
        ugw(i,njh+1) =0.5*(ugw(i,njh)+ugw(i,njh+2))
        vgw(i,njh+1) =0.5*(vgw(i,njh)+vgw(i,njh+2))
370   continue
c
c  Does not use coastal points for determining wind
c  at oceanic grid.
c
      do j=js1,js2
        jp1 = j+1
        do i=is1(j),is2(j)
          ip1         = (i+1)
          usp         = 1.0/
     &    max(1.0e-12*one,tmu(i,j,ks2)+tmu(ip1,j,ks2)
     &                   +tmu(i,jp1,ks2)+tmu(ip1,jp1,ks2))
          uwind(i,j)  = (ugw(i,j)*tmu(i,j,ks2)+
     &                   ugw(ip1,j)*tmu(ip1,j,ks2)+
     &                   ugw(i,jp1)*tmu(i,jp1,ks2)+
     &                   ugw(ip1,jp1)*tmu(ip1,jp1,ks2))*usp
          vwind(i,j)  = (vgw(i,j)*tmu(i,j,ks2)+
     &                   vgw(ip1,j)*tmu(ip1,j,ks2)+
     &                   vgw(i,jp1)*tmu(i,jp1,ks2)+
     &                   vgw(ip1,jp1)*tmu(ip1,jp1,ks2))*usp
        enddo
      enddo
c
400   continue
c
c--3.2. Calculate wind and wind stress at the center of the grid squares.
c------------------------------------------------------------------------
c
c  gam:  enhancing factor of OVERLAND and COLONY, 1994.
c
      gam = 1.3
      do 390 j=js1,js2
C       ihm   = max(0,sign(1,j-jeq))
C       hem   = 2.0*real(ihm)-1.0
C       thet  = 0.5*((theta(jmm1,2-ihm)+theta(jm,2-ihm))*ttat+
C    &               (theta(jm,2-ihm)+theta(jmp1,2-ihm))*ttbt)
C       rhc2  = 0.5*((rhcg2(jmm1,2-ihm)+rhcg2(jm,2-ihm))*ttat+
C    &               (rhcg2(jm,2-ihm)+rhcg2(jmp1,2-ihm))*ttbt)
C       sint  = sin(thet)
C       cost  = cos(thet)
C       costa = cos(15.0*radian)
C       sinta = sin(15.0*radian)
        do 380 i=is1(j),is2(j)
C         vabq(i,j)   = max(0.1*one,sqrt(uwind(i,j)*uwind(i,j)
C    &                       +vwind(i,j)*vwind(i,j)))
C         uwind(i,j)  = gam*uwind(i,j)
C         vwind(i,j)  = gam*vwind(i,j)
C         vabq(i,j)   = gam*vabq(i,j)
C         tenagx(i,j) = rhc2*vabq(i,j)*
C    &                  (uwind(i,j)*cost-vwind(i,j)*sint*hem)
C         tenagy(i,j) = rhc2*vabq(i,j)*
C    &                  (vwind(i,j)*cost+uwind(i,j)*sint*hem)
c
C         rhoa     =  psbq(i,j)/(287.04*tabq(i,j))
C         tairox(i,j) = cdb(vabq(i,j),tabq(i,j)
C    &                  -scal(i,j,ks2,1))*rhoa*vabq(i,j)*
C    &                  (uwind(i,j)*costa-vwind(i,j)*sinta*hem)
C         tairoy(i,j) = cdb(vabq(i,j),tabq(i,j)
C    &                  -scal(i,j,ks2,1))*rhoa*vabq(i,j)*
C    &                  (vwind(i,j)*costa+uwind(i,j)*sinta*hem)
C         tairox(i,j) = 1.3e-03*rhoa*vabq(i,j)*
C    &                  (uwind(i,j)*costa-vwind(i,j)*sinta*hem)
C         tairoy(i,j) = 1.3e-03*rhoa*vabq(i,j)*
C    &                  (vwind(i,j)*costa+uwind(i,j)*sinta*hem)
c
C         tairox(i,j) = trans(i,j)*tairox(i,j)
C    &                 +(1-trans(i,j))*tenhrx(i,j)
C         tairoy(i,j) = trans(i,j)*tairoy(i,j)
C    &                 +(1-trans(i,j))*tenhry(i,j)
C vent de Esbensen et Kushnir
 
          vabq(i,j)=vabqec(i,j)
c
c tension Trenberth
          zmod=1.0
          angt=0.*radian
          costt=cos(angt)
          sintt=sin(angt)
          tairox(i,j) = zmod*(tenhrx(i,j)*costt-tenhry(i,j)*sintt*hem)
          tairoy(i,j) = zmod*(tenhrx(i,j)*sintt*hem+tenhry(i,j)*costt)
          angt=0.*radian
          costb=cos(angt)
          sintb=sin(angt)
          tenagx(i,j) = zmod*(tenhrx(i,j)*costb-tenhry(i,j)*sintb*hem)
          tenagy(i,j) = zmod*(tenhrx(i,j)*sintb*hem+tenhry(i,j)*costb)
 
380     continue
390   continue
c
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine forcng -
      end
