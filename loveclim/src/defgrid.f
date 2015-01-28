












      subroutine defgrid(nflag)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c initatilisation of the masque and of the array linked with the grid.
c If nflag=2 : write testing varaibles on the file "mouchard" ;
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 23/08/02
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'reper.com'
      include 'dynami.com'
C     include 'met.inc.com'
 
c ntypcn    Water type
c rrro      --
c dd1o       | Parameter for the abosorption of solar radiation in the ocean
c dd2o      __
c dx1       Relative length of a grid square (direction x). Gen. in radian
c dx2       Relative length of a grid square (direction y). Gen. in radian
c h1        Metric coefficient in the direction x (defined at the center
c           of the grid )
c h2        Metric coefficient in the direction y (defined at the center
c           of the grid )
c d2d1      Derivatite of h2 in the x direction (defined at the center)
c d1d2      Derivatite of h1 in the y direction (defined at the center)
c h1p       Idem h1 for the bottom left corner of the grid
c h2p       Idem h2 for the bottom left corner of the grid
c d2d1p     Idem d2d1 for the bottom left corner of the grid
c d1d2p     Idem d1d2 for the bottom left corner of the grid
c h1pp      Idem h1 for the middle up of the grid
c h2pp      Idem h2 for the middle up of the grid
c
c--local variables :
      dimension zbath(kmax), dzbath(kmax)
cDFG hs now defined in bloc.com
c      dimension hs(imax,jmax), unszw(kmax)
      dimension unszw(kmax)
      dimension csxsa(jmax), csxua(jmax)
C     dimension snxua(jmax)
      dimension zrr(5),zd1(5),zd2(5)
      dimension tauco(20),zlatt(imax,jmax),zlont(imax,jmax)
      dimension ipt0v(10), jpt0v(10)
      dimension ntypcn(imax,jmax),rrro(imax,jmax)
      dimension dd1o(imax,jmax),dd2o(imax,jmax)
 
      dimension dx1(imax,jmax),dx2(imax,jmax),
     &          h1(imax,jmax),h2(imax,jmax),
     &          h1p(imax,jmax),h2p(imax,jmax),
     &          d1d2p(imax,jmax),d2d1p(imax,jmax),
     &          d1d2(imax,jmax),d2d1(imax,jmax)
     &         ,h1pp(imax,jmax),h2pp(imax,jmax)

cph
cph jsep2   aux. array to mask unneeded regions in rotated grid 
cph
      integer jsep2(83:imax),isep2(29:jmax)
cph
      character*1  cc1(0:9)
      character*8  fmt1
      character*30 fmt, fmtrep

      common / coord /  zlont,zlatt

      fotr(rf,d1f,d2f,h) = rf*exp(-h/d1f)+(1.0-rf)*exp(-h/d2f)
      data tauco  /6.6,6.6,7.0,7.2,7.1,6.8,6.5,6.6,7.1,7.6,
     &             6.6,6.1,5.6,5.5,5.8,5.8,5.6,5.6,5.6,5.6/
C    &             6.6,6.1,5.6,5.5,5.8,6.1,6.2,6.0,5.6,5.6/
cph
      data jsep2 / 49, 50, 50, 50, 50, 49, 47, 46, 46, 47, 47, 47,
     &             56, 56, 56, 56, 63, 64, 65, 65, 65, 65, 65, 65,
     &             65, 64, 64, 61, 61, 60, 59, 46, 46, 46, 46, 46,
     &             46, 46, 46, 29 /
      data isep2 /116,109,108,107,107,107,107,108,109,109,110,113,
     &            118,120,121,121,110,109,112,112,110,110,110,110,
     &            111,113,113,113,113,112,111,109,109,109,107,107,
     &            105 /
c      data jsep2 / 49, 50, 50, 50, 50, 49, 47, 46, 46, 47, 47, 47,
c     &             47, 56, 56, 55, 55, 64, 65, 65, 65, 65, 64, 64,
c     &             63, 63, 60, 60, 59, 58, 45, 45, 45, 45, 45, 45,
c     &             45, 45, 45, 29 /
cph 
      kmaxp1 = kmax + 1
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation of local and global variables (common) .         |
c-----------------------------------------------------------------------
 
      do 1 n=1,9
        cc1(n) = '-'
 1    continue
      cc1(0) = '|'
      cc1(5) = '5'
 
      do 11 j=1,jmax
       is1(j) = imax
       is2(j) = 1
       iu1(j) = imax
       iu2(j) = 1
       iuf1(j) = imax
       iuf2(j) = 1
       isf1(j) = imax
       isf2(j) = 1
 11   continue
 
      do 12 k=1,kmax
        z(k)  = 0.0
        zw(k) = 0.0
        dz(k)  = 0.0
        dzw(k) = 0.0
        unsdz(k)  = 0.0
        unsdzw(k) = 0.0
 12   continue
      z(kmax+1) = 0.0
      zw(kmax+1) = 0.0
      dzw(kmax+1) = 0.0
      unsdzw(kmax+1) = 0.0
 
      do 20 j=1,jmax
       do 20 i=1,imax
        hs(i,j) = 0.0
        huy(i,j) = 0.0
        hux(i,j) = 0.0
        hu(i,j) = 0.0
        unshu(i,j) = 0.0
C       kfs(i,j) = kmax
C       kfu(i,j) = kmax
        kniv(i,j,-1) = kmaxp1
        kniv(i,j, 0) = kmaxp1
        kniv(i,j, 1) = kmaxp1
        phifu(i,j)  = 0.0
        phifv(i,j)  = 0.0
        xang1(i,j)=0.0
        xang2(i,j)=1.0
cph
cph preset t- and u-grid longitudes and latitudes, and the
cph respective gridcell corners
cph
        xslon(i,j  )=90.0
        yslat(i,j  )=90.0
        xsedg(i,j,1)=90.0
        xsedg(i,j,2)=90.0
        xsedg(i,j,3)=90.0
        xsedg(i,j,4)=90.0
        ysedg(i,j,1)=90.0
        ysedg(i,j,2)=90.0
        ysedg(i,j,3)=90.0
        ysedg(i,j,4)=90.0
        xulon(i,j  )=90.0
        yulat(i,j  )=90.0
        xuedg(i,j,1)=90.0
        xuedg(i,j,2)=90.0
        xuedg(i,j,3)=90.0
        xuedg(i,j,4)=90.0
        yuedg(i,j,1)=90.0
        yuedg(i,j,2)=90.0
        yuedg(i,j,3)=90.0
        yuedg(i,j,4)=90.0
        angle(i,j)=0.0
cph
 20   continue
      do j=1,jmax+1
        do i=1,imax+1
          xslonp(i,j)=90.0
          yslatp(i,j)=90.0
          xslonp(i,j)=90.0
          yslatp(i,j)=90.0
        end do
      end do

      do 21 ns=0,nsmax
       do 21 j=1,jmax
        do 21 i=1,imax
         fss(i,j,ns)  = 0.0
         phifs(i,j,ns)  = 0.0
         phiss(i,j,ns)  = 0.0
         rappes(i,j,ns) = 0.0
         phimnx(i,j,0,ns) = cstmin
         phimnx(i,j,1,ns) = cstmax
 21   continue
 
      do 22 nn=1,6
       do 22 j=1,jmax
        do 22 i=1,imax
         phihhh(i,j,nn) = 0.0
 22   continue
 
      do 30 k=1,kmax
       do 30 j=1,jmax
        do 30 i=1,imax
         fub(i,j,k) = 0.0
         fvb(i,j,k) = 0.0
         q(i,j,k) = 0.0
         bvf(i,j,k) = 0.0
         avsdz(i,j,k) = 0.0
         avudz(i,j,k) = 0.0
         fqajc(i,j,k) = 0.0
         tms(i,j,k) = 0.0
         tmu(i,j,k) = 0.0
         rappel(i,j,k) = 0.0
 30   continue
 
      do 31 k=1,kmax+1
       do 31 j=1,jmax
        do 31 i=1,imax
         w(i,j,k) = 0.0
         phizzz(i,j,k,1) = 0.0
         phizzz(i,j,k,2) = 0.0
 31   continue
 
c--Initialisation of the scalairs,  forcings, and of the
c  turbuent kinetic energy
      do 32 ns=1,nsmax
       do 32 k=1,kmax
        nrap(k,ns) = 0
        do 32 j=1,jmax
         do 32 i=1,imax
          scalr(i,j,k,ns)=0.0
          scal(i,j,k,ns) = scal0(k,ns)
          phivs(i,j,k,ns) = 0.0
 32   continue
 
      do 33 k=1,kmax
       do 33 j=1,jmax
        do 33 i=1,imax
         q2turb(i,j,k)=q2tmin
 33   continue
 
c--Initialisation of the scalars associated with the freshwater flux
      do 34 ns=1,nsmax
       do 34 j=1,jmax
        do 34 i=1,imax
         scs(i,j,ns) = scpme(ns)
 34   continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Initialosarion of the constants ( resolution, grid ...)          |
c-----------------------------------------------------------------------
 
c--size (maximum)  of the basin (pts type scalar) :
      ims1 = 2
      ims2 = imax - 1
      js1 = 2
      js2 = jmax - 1
      ks1 = 1
      ks2 = kmax
c--others :
      ku2 = ks2
      ku1 = ks1
 
c--constantes linked with the resolution / axes /  Geography :
      call geogra(js1, js2, jeq, jdl1, jdl2, ijsdl, ijudl,
     &            iberp, ibera)
 
      if (ltest.lt.1) then
c--basin type : rectangular
        jcl1 = js2
        jcl2 = js1
      else
c--basin with cyclic boundaries, definition of the linked areas :
        jcl1 = js1
        jcl2 = jeq
      endif
 
c--initialisation of the constants
      if (ltest.lt.0) then
        dx = 1000.
        dy = 1000.
      else
        dx = dlong * radian * rterre
        dy = dlat  * radian * rterre
      endif
      unsdx = 1.0 / dx
      unsdy = 1.0 / dy
      uns2dx = 0.5 / dx
      uns2dy = 0.5 / dy
 
c- Warning : Bering Stait closed => (iberp,a)=0,0. But not when living "defgrid"
      jberpm = jberp - 1
      iberpm = iberp - 1
      jberam = jbera - 1
      iberam = ibera - 1
 
c--initialisation of island with no surface :
c- Spitzberg :
      ipt0v(1) = 109
      jpt0v(1) =  55
Cic0  npt0v =  1
      npt0v =  0
 
c--point with ocean velocity but no sea-ice velocity 
      npo1i0 = 1
      npo1i0 = 2
      ipo1i0(1) = 102
      jpo1i0(1) = 56
      ipo1i0(2) = 111
      jpo1i0(2) = 56
 
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) For each grid, initialisation of metric coefficients .   |
c-----------------------------------------------------------------------
 
c-------------------
c Conventions :    |
c-------------------
c  index 1 \ 2 <-> location in the direction x \ y
c    3eme indice <-> location on the grid :
c  0=Centre, 1=Side W x(i-1/2), 2=side S y(j-1/2), 3=corner SW x(i-1/2),y(j-1/2)
c-------------------
c Metric coeficients : cmx \ cmy = metric. coef. . direction x (h1) \ y (h2)
c Derivatives of the metriq. coef. : cmxdy \ cmydx = d(cmx)/dy  \ d(cmy)/dx located.3
c
c Other coeficients built from metric. coef. :
c  cmxy(i,j,0\3) = cmx(i,j,0\3) * cmy(i,j,0\3) = element de surface
c  cmxy(i,j,1\2) = cmx(i,j,1\2) / cmy(i,j,1\2)
c  smx = 1 / cmx ; smy = 1 / cmy ; smxy = 1 / cmxy
c-------------------
c loacl conventions :
c  cs <- cosine ; sc <- 1 / cosine ; sn <- sine [ou 1st letter doubled]
c  x\y_s\u_w\a : related to a
c   longitude (x) \ latitude (y) <-- differs from the convention metric coef. 
c   for a scalar point (=Centre) (s) \ velocity (corner) (u) [default = s]
c   on the grid WW (w) \ AA (a)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--3.1 1st grid classical spheric , grid World (WW) :
c-------------------------------------------------------
 
c--initialisation of the geometric arrays cos,sin, etc ...
      do 320 j=1,jsepar
        ysdeg = ylat1 + dlat * DFLOAT(j-1)
cDFG
        yudeg = ysdeg - 0.5 * dlat
cDFG
        yurad = (ysdeg - 0.5 * dlat) * radian
c-
        if (ltest.lt.0) then
c--Cartesian Grid :
          yu0rad = 45.0 * radian
          ccsysw = 1.0
          sscysw = 1.0
          ccsyuw = 1.0
          sscyuw = 1.0
          ccmxdy = 0.0
          ccfcor = omega * sin(yu0rad)
          ccfcu8 = -0.25 * omega * cos(yu0rad)
        else
c--Spherical grid:
          ccsysw = cos(ysdeg * radian)
          sscysw = 1.0 / ccsysw
          ccsyuw = cos(yurad)
          sscyuw = 1.0 / ccsyuw
          ssnyuw = sin(yurad)
          ssnysw = sin(ysdeg*radian)
          ccmxdy = -ssnyuw * unsrt
          ccfcor = omega * ssnyuw
          ccfcu8 = -0.25 * omega * ccsyuw
c--end of the difference between  cartesian/spheric .
        endif
 
c--metric : spherical grid long x lat :
        do 310 i=1,imax
          zlatt(i,j)=yurad
cph
cph simple lon-lat-grid with longitude interval dlong and
cph latitudinal width dlat. corners are indexed as follows:
cph  1) SW  2) SE  3) NE  4) NW  of respective gridpoint
cph momentum-gridpoint i,j is at SW corner of corresponding
cph tracer point.
cph
          yslat(i,j)=ysdeg
          ysedg(i,j,1)=yudeg
          ysedg(i,j,2)=yudeg
          ysedg(i,j,3)=yudeg+dlat
          ysedg(i,j,4)=yudeg+dlat
          yulat(i,j)=yudeg
          yuedg(i,j,1)=ysdeg-dlat
          yuedg(i,j,2)=ysdeg-dlat
          yuedg(i,j,3)=ysdeg
          yuedg(i,j,4)=ysdeg
          xsdeg = xlon1 + dlong * DFLOAT(i-1)
          xurad = (xsdeg - 0.5 * dlong) * radian
          xudeg = xsdeg - 0.5 * dlong  
          xslon(i,j)=xsdeg
          xsedg(i,j,1)=xudeg
          xsedg(i,j,2)=xudeg+dlong
          xsedg(i,j,3)=xudeg+dlong
          xsedg(i,j,4)=xudeg
          xulon(i,j)=xudeg
          if ((xsdeg-dlong).le.0.0d0) then
            xuedg(i,j,1)=xsdeg-dlong+360.0d0
          else
            xuedg(i,j,1)=xsdeg-dlong
          endif
          xuedg(i,j,2)=xsdeg
          xuedg(i,j,3)=xsdeg
          if ((xsdeg-dlong).le.0.0d0) then
            xuedg(i,j,4)=xsdeg-dlong+360.0d0
          else
            xuedg(i,j,4)=xsdeg-dlong
          endif

          zlont(i,j)=xurad
c- coef :
          cmx(i,j,0) = ccsysw
          cmx(i,j,1) = ccsysw
          cmx(i,j,2) = ccsyuw
          cmx(i,j,3) = ccsyuw
          cmy(i,j,0) = 1.
          cmy(i,j,1) = 1.
          cmy(i,j,2) = 1.
          cmy(i,j,3) = 1.
c- inverse :
          smx(i,j,0) = sscysw
          smx(i,j,1) = sscysw
          smx(i,j,2) = sscyuw
          smx(i,j,3) = sscyuw
          smy(i,j,0) = 1.
          smy(i,j,1) = 1.
          smy(i,j,2) = 1.
          smy(i,j,3) = 1.
c- product, ratio of 2 coefs :
          cmxy(i,j,0) = ccsysw
          cmxy(i,j,1) = ccsysw
          cmxy(i,j,2) = ccsyuw
          cmxy(i,j,3) = ccsyuw
          smxy(i,j,0) = sscysw
          smxy(i,j,1) = sscysw
          smxy(i,j,2) = sscyuw
          smxy(i,j,3) = sscyuw
c- derivatives :
          cmxdy(i,j) = ccmxdy
          cmydx(i,j) = 0.
c- metric coefficients for sea ice dynamic:
          dx1(i,j)   = dlong*radian
          yu2rad     = (ysdeg + 0.5*dlat)*radian
          dx2(i,j)   = sin(yu2rad) -sin(yurad)
          h1(i,j)    = rterre*ccsysw
          h2(i,j)    = rterre/ccsysw
          d1d2(i,j)  = -rterre*sin(ysdeg*radian)*sscysw
          d2d1(i,j)  = 0.0
 310    continue

cph - rotation of grid (should yield zero for simple lon-lat grid)
        do i=1,imax-1
          angle(i,j)=atan2((yulat(i+1,j)-yulat(i,j)),
     &                     (xulon(i+1,j)-xulon(i,j)))
        end do
        angle(imax,j)=angle(imax-1,j)

c--Coriolis factor  :
        do 315 i=1,imax
          fs2cor(i,j) = ccfcor
Cfcc      fcucor(i,j) = ccfcu8
Cfcc      fcvcor(i,j) = 0.0
c- Sine (Geog. Lat.) :
          covrai(i,j) = ssnysw
 315    continue
 320  continue
 
      if (ltest.ge.3) then
 
c--3.3 North Atlantic + Arctic (AA), sherical grid N-S / E-W inversed
c-----------------------------------------------------------------------------

C       xudeg0 = xaj1  - 1.5 * dxaj
        xu0m90 = xaj1  - 1.5 * dxaj - 90.
        do 330 j=jeq,jmax
C         xudeg = xudeg0 + dxaj * DFLOAT(j)
          xudeg = xu0m90 + dxaj * DFLOAT(j)
c- Caution : change of  sign ! csxua = -cos(Long_AA)
C         csxua(j) = -cos(xudeg * radian)
          csxua(j) = sin(xudeg * radian)
Cfcc      snxua(j) = cos(xudeg * radian)
          csxsa(j) = sin((xudeg+0.5*dxaj)*radian)
 330    continue
        yyalim = 90.0 - abs(dyai)
        do 350 i=1,imax
          ysdeg = yai1 + dyai * DFLOAT(i-1)
          if (abs(ysdeg).ge.yyalim) then
            ccsysa = 1.
            yurad =  0.
            ysrad =  0.
            ccsyua = 1.
            ssnyua = 0.
          else
            ccsysa = cos(ysdeg * radian)
            yurad  = (ysdeg - 0.5 * dyai) * radian
            ysrad  = ysdeg * radian
            ccsyua = cos(yurad)
Cfcc        ssnyua = sin(yurad)
          endif
          sscysa = 1.0 / ccsysa
          sscyua = 1.0 / ccsyua
          ccmydx =  sin(yurad) * unsrt
c--Metrique : 2nd spherical grid lat x long :
          do 340 j=jsep(i),jmax
c- coef :
            cmx(i,j,0) = 1.
            cmx(i,j,1) = 1.
            cmx(i,j,2) = 1.
            cmx(i,j,3) = 1.
            cmy(i,j,0) = ccsysa
            cmy(i,j,1) = ccsyua
            cmy(i,j,2) = ccsysa
            cmy(i,j,3) = ccsyua
c- inverse :
            smx(i,j,0) = 1.
            smx(i,j,1) = 1.
            smx(i,j,2) = 1.
            smx(i,j,3) = 1.
            smy(i,j,0) = sscysa
            smy(i,j,1) = sscyua
            smy(i,j,2) = sscysa
            smy(i,j,3) = sscyua
 
c- product, ratio de 2 coefs :
            cmxy(i,j,0) = ccsysa
            cmxy(i,j,1) = sscyua
            cmxy(i,j,2) = sscysa
            cmxy(i,j,3) = ccsyua
            smxy(i,j,0) = sscysa
            smxy(i,j,1) = ccsyua
            smxy(i,j,2) = ccsysa
            smxy(i,j,3) = sscyua
c- derivatives :
            cmxdy(i,j) = 0.
            cmydx(i,j) = ccmydx
c- metric coefficients for sea ice dynamic:
            dx2(i,j)   = dlat*radian
            yu2rad     = (ysdeg + 0.5*dyai)*radian
            dx1(i,j)   = -sin(yu2rad) + sin(yurad)
            h1(i,j)    = rterre/ccsysa
            h2(i,j)    = rterre*ccsysa
            d1d2(i,j)  = 0.0
            d2d1(i,j)  = +rterre*sin(ysdeg*radian)*sscysa
            xsrad=(xaj1+dxaj*DFLOAT(j-1))*radian
            xurad=(xaj1-1.5*dxaj+dxaj*DFLOAT(j))*radian
            yslat(i,j)=-asin(cos(ysrad)*cos(xsrad))/radian
            yulat(i,j)=-asin(cos(yurad)*cos(xurad))/radian
            zlatt(i,j)=-asin(cos(yurad)*cos(xurad))
            szphi=(cos(yurad)*cos(xurad))
            zlam=atan2(cos(ysrad)*sin(xsrad),sin(ysrad))
            xslon(i,j)=(zlam+pi)/radian+69.0
            zlam=atan2(cos(yurad)*sin(xurad),sin(yurad))
            xulon(i,j)=(zlam+pi)/radian+69.0 
            zlont(i,j)=zlam+pi+69.0*pi/180.0
CCice        xang1(i,j)=(-sin(zphi)*cos(zlam))/max(zepsd1,cos(yurad))
            xang1(i,j)=(szphi*cos(zlam))/max(zepsd1,cos(yurad))
            xang2(i,j)=(sin(zlam))/max(zepsd1,cos(yurad))
 340      continue
c--coef " between the two grids 2 grids :
          if (jsep(i).eq.jeq) then
            cmy(i,jeq,2) = 0.5 * ( 1. + cmy(i,jeq,0) )
            cmy(i,jeq,3) = 0.5 * ( 1. + cmy(i,jeq,1) )
            smy(i,jeq,2) = 1. / cmy(i,jeq,2)
            smy(i,jeq,3) = 1. / cmy(i,jeq,3)
            cmxy(i,jeq,2) = smy(i,jeq,2)
            cmxy(i,jeq,3) = cmy(i,jeq,3)
            smxy(i,jeq,2) = cmy(i,jeq,2)
            smxy(i,jeq,3) = smy(i,jeq,3)
            cmydx(i,jeq) = 0.5 * ccmydx
          endif
c--Coriolis factor :
          do 345 j=jsep(i),jmax
            fs2cor(i,j) = omega * ccsyua * csxua(j)
Cfcc        fcucor(i,j) = omega * snxua(j) * -0.25
Cfcc        fcvcor(i,j) = omega * ssnyua * csxua(j) * 0.25
c- Sine (Geog. Lat) :
            covrai(i,j) = ccsysa * csxsa(j)
 345      continue
 350    continue
cph
cph set momentum gridpoints as corners of tracer cells
cph and vice-versa in rotated grid
cph
        do i=1,imax
          do j=jsep(i),jmax
            ipm=min(i+1,imax)
            jpm=min(j+1,jmax)
            xsedg(i,j,1)=xulon(i  ,j  )
            xsedg(i,j,2)=xulon(ipm,j  )
            xsedg(i,j,3)=xulon(ipm,jpm)
            xsedg(i,j,4)=xulon(i  ,jpm)
            ysedg(i,j,1)=yulat(i  ,j  )
            ysedg(i,j,2)=yulat(ipm,j  )
            ysedg(i,j,3)=yulat(ipm,jpm)
            ysedg(i,j,4)=yulat(i  ,jpm)
            ipm=max(i-1,1)
            jpm=max(j-1,1)
            xuedg(i,j,1)=xslon(ipm,jpm)
            xuedg(i,j,2)=xslon(i  ,jpm)
            xuedg(i,j,3)=xslon(i  ,j  )
            xuedg(i,j,4)=xslon(ipm,j  )
            yuedg(i,j,1)=yslat(ipm,jpm)
            yuedg(i,j,2)=yslat(i  ,jpm)
            yuedg(i,j,3)=yslat(i  ,j  )
            yuedg(i,j,4)=yslat(ipm,j  )
          end do
        end do
cph
cph avoid unsteady, huge coordinates in the outer space
cph between the two grids
cph
        do i=1,82
          yfac=0.25/float(jmax-jsep(i))
          jpm=jsep(i)-1
          do j=jsep(i),jmax
            xslon(i,j)=xslon(i,jpm)
            yslat(i,j)=yslat(i,jpm)+yfac*(j-jpm)
            xsedg(i,j,1)=xsedg(i,jpm,1)
            xsedg(i,j,2)=xsedg(i,jpm,2)
            xsedg(i,j,3)=xsedg(i,jpm,3)
            xsedg(i,j,4)=xsedg(i,jpm,4)
            ysedg(i,j,1)=yslat(i,j)-0.5*yfac
            ysedg(i,j,2)=yslat(i,j)-0.5*yfac
            ysedg(i,j,3)=yslat(i,j)+0.5*yfac
            ysedg(i,j,4)=yslat(i,j)+0.5*yfac
            xulon(i,j)=xulon(i,jpm)
            yulat(i,j)=yulat(i,jpm)+yfac*(float(j-jpm)-0.5)
            xuedg(i,j,1)=xuedg(i,jpm,1)
            xuedg(i,j,2)=xuedg(i,jpm,2)
            xuedg(i,j,3)=xuedg(i,jpm,3)
            xuedg(i,j,4)=xuedg(i,jpm,4)
            yuedg(i,j,1)=yulat(i,j)-0.5*yfac
            yuedg(i,j,2)=yulat(i,j)-0.5*yfac
            yuedg(i,j,3)=yulat(i,j)+0.5*yfac
            yuedg(i,j,4)=yulat(i,j)+0.5*yfac
          end do
        end do
cph
cph same procedure in the arctic ocean region
cph
        do i=83,imax
          if (jsep2(i).ne.jmax) then
            yfac=0.25/float(jmax-jsep2(i))
            jpm=jsep2(i)-1
            do j=jsep2(i),jmax
              xslon(i,j)=xslon(i,jpm)
              yslat(i,j)=yslat(i,jpm)+yfac*(j-jpm)
              xsedg(i,j,1)=xsedg(i,jpm,1)
              xsedg(i,j,2)=xsedg(i,jpm,2)
              xsedg(i,j,3)=xsedg(i,jpm,3)
              xsedg(i,j,4)=xsedg(i,jpm,4)
              ysedg(i,j,1)=yslat(i,j)-0.5*yfac
              ysedg(i,j,2)=yslat(i,j)-0.5*yfac
              ysedg(i,j,3)=yslat(i,j)+0.5*yfac
              ysedg(i,j,4)=yslat(i,j)+0.5*yfac
              xulon(i,j)=xulon(i,jpm)
              yulat(i,j)=yulat(i,jpm)+yfac*(float(j-jpm)-0.5)
              xuedg(i,j,1)=xuedg(i,jpm,1)
              xuedg(i,j,2)=xuedg(i,jpm,2)
              xuedg(i,j,3)=xuedg(i,jpm,3)
              xuedg(i,j,4)=xuedg(i,jpm,4)
              yuedg(i,j,1)=yulat(i,j)-0.5*yfac
              yuedg(i,j,2)=yulat(i,j)-0.5*yfac
              yuedg(i,j,3)=yulat(i,j)+0.5*yfac
              yuedg(i,j,4)=yulat(i,j)+0.5*yfac
            end do
          endif
        end do

cph
cph angle of rotation of grid
cph
        do i=1,imax-1
          do j=jsep(i),jmax
            angle(i,j)=atan2((yulat(i+1,j)-yulat(i,j)),
     &                       (xulon(i+1,j)-xulon(i,j)))
          end do
        end do
        i=imax
        do j=jsep(i),jmax
          angle(i,j)=angle(i-1,j)
        end do

c-- end of the initialisation of the AA grid .
c-----
      endif
cph
cph second set of edges for the 2nd form of the 3-argument shade in 
cph ferret (the four point algorithm)
cph
      do j=1,jmax
        do i=1,imax
          xslonp(i,j)=xsedg(i,j,1)
          yslatp(i,j)=ysedg(i,j,1)
          xulonp(i,j)=xuedg(i,j,1)
          yulatp(i,j)=yuedg(i,j,1)
        end do
        xslonp(imax+1,j)=xsedg(imax,j,2)
        yslatp(imax+1,j)=ysedg(imax,j,2)
        xulonp(imax+1,j)=xuedg(imax,j,2)
        yulatp(imax+1,j)=yuedg(imax,j,2)
      end do
      do i=1,imax
        xslonp(i,jmax+1)=xsedg(i,jmax,4)
        yslatp(i,jmax+1)=ysedg(i,jmax,4)
        xulonp(i,jmax+1)=xuedg(i,jmax,4)
        yulatp(i,jmax+1)=yuedg(i,jmax,4)
      end do
      xslonp(imax+1,jmax+1)=xsedg(imax,jmax,3)
      yslatp(imax+1,jmax+1)=ysedg(imax,jmax,3)
      xulonp(imax+1,jmax+1)=xuedg(imax,jmax,3)
      yulatp(imax+1,jmax+1)=yuedg(imax,jmax,3)

      if (ltest.eq.3 .and. iberp.ne.0) then
C       smy(iberp-1,jberp,2) = 0.
C       smy(iberp , jberp,2) = 0.
c--Bering : local modification (AA) of the Metric Coefs :
        cmx(ibera, jbera,2) = cmx(iberpm,jberp,2)
        smx(ibera, jbera,2) = smx(iberpm,jberp,2)
        cmy(ibera, jbera,2) = cmy(iberpm,jberp,2)
        smy(ibera, jbera,2) = smy(iberpm,jberp,2)
        cmxy(ibera, jbera,2) = cmxy(iberpm,jberp,2)
        smxy(ibera, jbera,2) = smxy(iberpm,jberp,2)
c-
        cmx(iberam,jbera,2) = cmx(iberp, jberp,2)
        smx(iberam,jbera,2) = smx(iberp, jberp,2)
        smx(iberam,jbera,2) = smx(iberp, jberp,2)
        cmy(iberam,jbera,2) = cmy(iberp, jberp,2)
        smy(iberam,jbera,2) = smy(iberp, jberp,2)
        cmxy(iberam,jbera,2) = cmxy(iberp, jberp,2)
        smxy(iberam,jbera,2) = smxy(iberp, jberp,2)
c--Only for the diagnostics :
        cmxy(ibera,jbera,3) = cmxy(iberp,jberp,3)
        smxy(ibera,jbera,3) = smxy(iberp,jberp,3)
        cmxy(iberp, jberp,0) = cmxy(iberam,jberam,0)
        cmxy(iberpm,jberp,0) = cmxy(ibera, jberam,0)
        cmxy(ibera, jbera,0) = cmxy(iberpm,jberpm,0)
        cmxy(iberam,jbera,0) = cmxy(iberp, jberpm,0)
        smxy(iberp, jberp,0) = smxy(iberam,jberam,0)
        smxy(iberpm,jberp,0) = smxy(ibera, jberam,0)
        smxy(ibera, jbera,0) = smxy(iberpm,jberpm,0)
        smxy(iberam,jbera,0) = smxy(iberp, jberpm,0)
      endif
 
c--Modification of the metric  coeff. close to Gibraltar.
C     include 'defgrid.gibr.inc'
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Computation of metric coefficients needed in ice dynamics       |
c-----------------------------------------------------------------------
c
c
c--4.1. WEIGHTS FOR COMPUTING INTERPOLATED VALUES AT THE CORNERS
c       OF THE GRID SQUARES.
c----------------------------------------------------------------
c
      do 420 j=2,jmax
        jm1 = j-1
        do 410 i=1,imax
          im1           = (i-1)+(imax-2)*(1/i)
          usden         = 1.0/((dx1(i,j)+dx1(im1,j))
     &                    *(dx2(i,j)+dx2(i,jm1)))
          wght(i,j,1,1) = usden*dx1(i,j)*dx2(i,j)
          wght(i,j,1,2) = usden*dx1(i,j)*dx2(i,jm1)
          wght(i,j,2,2) = usden*dx1(im1,j)*dx2(i,jm1)
          wght(i,j,2,1) = usden*dx1(im1,j)*dx2(i,j)
 410    continue
 420  continue
c
c--4.2. metric coefficients at the grid squares corners,
c       and mid-sides points, and their derivatives (h1p,  h2p,
c       h1pp, h2pp, d1d2p,d2d1p), corilis coefficient (zfn)
c       and coefficients for the strain rate tensor (akappa).
c----------------------------------------------------------------
c
      do 430 i=1,imax
        h1p(i,1)   = 0.0
        h2p(i,1)   = 1.0e+20
        d1d2p(i,1) = 1.0e+20
        d2d1p(i,1) = 0.0
c
        zfn(i,1)   = fs2cor(i,1)*2
 430  continue
c
      do 450 j=1,jmax-1
        do 440 i=1,imax
          im1          = (i-1)+(imax-2)*(1/i)
          ip1          = (i+1)-(imax-2)*(i/imax)
          h1p(i,j+1)   = h1(i,j+1)*wght(i,j+1,2,2)+
     &                   h1(im1,j+1)*wght(i,j+1,1,2)+
     &                   h1(i,j)*wght(i,j+1,2,1)+
     &                   h1(im1,j)*wght(i,j+1,1,1)
          h2p(i,j+1)   = h2(i,j+1)*wght(i,j+1,2,2)+
     &                   h2(im1,j+1)*wght(i,j+1,1,2)+
     &                   h2(i,j)*wght(i,j+1,2,1)+
     &                   h2(im1,j)*wght(i,j+1,1,1)
          d1d2p(i,j+1) = 2.0*
     &                   (dx1(i,j)*(-h1(im1,j)+h1(im1,j+1))+
     &                    dx1(im1,j)*(-h1(i,j)+h1(i,j+1)))/
     &                   ((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,j+1)))
          d2d1p(i,j+1) = 2.0*
     &                   (dx2(i,j+1)*(h2(i,j)-h2(im1,j))+
     &                    dx2(i,j)*(h2(i,j+1)-h2(im1,j+1)))/
     &                   ((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,j+1)))
c
          h1pp(i,j)    = (dx2(i,j)*h1(i,j+1)+dx2(i,j+1)*h1(i,j))/
     &                   (dx2(i,j)+dx2(i,j+1))
          h2pp(i,j)    = (dx1(i,j)*h2(ip1,j)+dx1(ip1,j)*h2(i,j))/
     &                   (dx1(i,j)+dx1(ip1,j))
c
          zfn(i,j+1)   = 2*fs2cor(i,j+1)
 440    continue
 450  continue
      do 460 i=1,imax
        ip1          = (i+1)-(imax-2)*(i/imax)
        h1pp(i,jmax) = 0.0
        h2pp(i,jmax) = (dx1(i,j)*h2(ip1,j)+dx1(ip1,j)*h2(i,j))/
     &                   (dx1(i,j)+dx1(ip1,j))
 460  continue
c
      do 475 j=1,jmax
        do 470 i=1,imax
          ip1             = (i+1)-(imax-2)*(i/imax)
          akappa(i,j,1,1) = 1.0/(2.0*h1(i,j)*dx1(i,j))
          akappa(i,j,1,2) = d1d2(i,j)/(4.0*h1(i,j)*h2(i,j))
          akappa(i,j,2,1) = d2d1(i,j)/(4.0*h1(i,j)*h2(i,j))
          akappa(i,j,2,2) =  1.0/(2.0*h2(i,j)*dx2(i,j))
 470    continue
 475  continue
c
c--4.3. COEFFICIENTS FOR DIVERGENCE OF THE STRESS TENSOR.
c---------------------------------------------------------
c
      do 485 j=2,jmax
        jm1 = j-1
        do 480 i=1,imax
          im1=(i-1)+imax*(1/i)
          usden               =
     &           1.0/(h1p(i,j)*h2p(i,j)*(dx1(i,j)+dx1(im1,j))
     &           *(dx2(i,j)+dx2(i,jm1)))
          alambd(i,j,2,2,2,1) = usden*2.0*dx2(i,j)*h2(i,jm1)
          alambd(i,j,2,2,2,2) = usden*2.0*dx2(i,jm1)*h2(i,j)
          alambd(i,j,2,2,1,1) = usden*2.0*dx2(i,j)*h2(im1,jm1)
          alambd(i,j,2,2,1,2) = usden*2.0*dx2(i,jm1)*h2(im1,j)
          alambd(i,j,1,1,2,1) = usden*2.0*dx1(im1,j)*h1(i,jm1)
          alambd(i,j,1,1,1,1) = usden*2.0*dx1(i,j)*h1(im1,jm1)
          alambd(i,j,1,1,2,2) = usden*2.0*dx1(im1,j)*h1(i,j)
          alambd(i,j,1,1,1,2) = usden*2.0*dx1(i,j)*h1(im1,j)
          alambd(i,j,1,2,1,1) = usden*d1d2p(i,j)*dx2(i,j)*dx1(i,j)
          alambd(i,j,1,2,2,1) = usden*d1d2p(i,j)*dx2(i,j)*dx1(im1,j)
          alambd(i,j,1,2,1,2) = usden*d1d2p(i,j)*dx2(i,jm1)*dx1(i,j)
          alambd(i,j,1,2,2,2) = usden*d1d2p(i,j)*dx2(i,jm1)*dx1(im1,j)
          alambd(i,j,2,1,1,1) = usden*d2d1p(i,j)*dx2(i,j)*dx1(i,j)
          alambd(i,j,2,1,2,1) = usden*d2d1p(i,j)*dx2(i,j)*dx1(im1,j)
          alambd(i,j,2,1,1,2) = usden*d2d1p(i,j)*dx2(i,jm1)*dx1(i,j)
          alambd(i,j,2,1,2,2) = usden*d2d1p(i,j)*dx2(i,jm1)*dx1(im1,j)
 480    continue
 485  continue
c
c--4.4 METRICS FOR ADVECTION AND SCALAR DIFFUSION COEFFICIENTS.
c--------------------------------------------------------------
c
c  dxs1, dxs2: LENGHT OF THE GRID SQUARES SIDES.
c  dxc1, dxc2: LENGHT OF THE GRID SQUARES CENTRES.
c  area: SURFACE OF A GRID SQUARE.
c
      do 495 j=1,jmax
        do 490 i=1,imax
          dxs1(i,j)  = h1pp(i,j)*dx1(i,j)
          dxs2(i,j)  = h2pp(i,j)*dx2(i,j)
          dxc1(i,j)  = h1(i,j)*dx1(i,j)
          dxc2(i,j)  = h2(i,j)*dx2(i,j)
          area(i,j)  = dxc1(i,j)*dxc2(i,j)
 490    continue
 495  continue
c
      do 497 j=2,jmax
       do 497 i=1,imax
         im1=(i-1)+(imax-2)*(1/i)
         bkappa(i,j,1,1) = dxs2(i,j)/area(i,j)
         bkappa(i,j,1,2) = dxs2(im1,j)/area(i,j)
         bkappa(i,j,2,2) = dxs1(i,j)/area(i,j)
         bkappa(i,j,2,1) = dxs1(i,j-1)/area(i,j)
 497  continue
c
      do 498 i=1,imax
         j=1
         bkappa(i,j,1,1) = dxs2(i,j)/area(i,j)
         bkappa(i,j,1,2) = dxs2(im1,j)/area(i,j)
         bkappa(i,j,2,2) = dxs1(i,j)/area(i,j)
         bkappa(i,j,2,1) = dxs1(i,1)/area(i,j)
 498  continue
c--modification of the metric coefficients for Bering (Arctic)
      if (ltest.eq.3 .and. iberp.ne.0) then
        dxs1(ibera-1,jbera-1) = dxs1(iberp,jberp-1)
        dxs1(ibera,jbera-1)   = dxs1(iberp-1,jberp-1)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Definition of the domain : reading of the bathymetry               |
c-----------------------------------------------------------------------
 
c--reading of the bathymetry file (thickness of each level (m) )
c        and of the map tha contains the number of levels)
      open(unit=20,file='bath.om',status='old',form='UNFORMATTED')
      read (20) zbath
      read (20) dzbath
      read (20) kbath
      close(20)
 
c--initialisation of the depth scale :
      zw(ks2+1) = 0.0
      do 510 k=ks2,ks1,-1
        z(k)  = -zbath(ks1+ks2-k)
        dz(k) = dzbath(ks1+ks2-k)
        zw(k)    = zw(k+1) - dz(k)
        unsdz(k) = 1.0 / dz(k)
        if (zw(k).ne.zero) unszw(k) = 1.0 / zw(k)
 510  continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 520 k=ks1+1,ks2
        dzw(k) = z(k)  - z(k-1)
        unsdzw(k) = 1.0 / dzw(k)
 520  continue

      if(kfond.gt.0) then
c--flat bottom:
      do 530 j=1,jmax
       do 530 i=1,imax
        kbath(i,j) = min(kfond,(kbath(i,j)*kmax))
 530  continue
      endif
 
c--closing the basin North and South (after save) :
      do 540 i=1,imax
        kbath1(i) = kbath(i,1)
        kbath(i,1) = 0
        kbath2(i) = kbath(i,jmax)
        kbath(i,jmax) = 0
 540  continue
 
c--closing the basin East and Ouest :
      do 550 j=1,jmax
        kbath(1,j) = 0
        kbath(imax,j) = 0
 550  continue
 
c--closing between the 2 grids  (except at the Equator) :
      if (ltest.ge.2) then
        do 560 i=1,imax
          if (jsep(i).ne.jeq) kbath(i,jsep(i)-1) = 0
 560    continue
      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Initialisation of the mask and of the indexes of beginning/end of domain.
c-----------------------------------------------------------------------
 
c--6.1 Centre of the grid  (pt. type scalar) - land/sea mask  pts Scal. :
c---------------------------------------------------------------------------
 
c- initialisation of the depth and computation of the index of beginning 
c  and end of basin  :
      do 615 j=js1,js2
       is1(j) = imax
       is2(j) = 1
       do 615 i=ims1,ims2
         km = kbath(i,j)
         if(km.gt.0) then
           kkm = ks1 + ks2 - km
           kniv(i,j,0) = kkm
           hs(i,j) = -zw(kkm)
           if(is1(j).eq.imax) is1(j) = i
           is2(j) = i
           do 610 k=kkm,kmax
             tms(i,j,k) = 1.0
 610       continue
         endif
 615  continue
 
c--6.2 definition of the "openings" for the E & S boundaries :
c---------------------------------------------------------------
 
      if (ltest.ge.1) then
c--half-cycli correspondance for kbath :
      do 620 j=jcl1,jcl2
        kbath(ims1-1,j) = kbath(ims2,j)
 620  continue
      endif
 
      if (ltest.eq.3 .and. iberp.ne.0) then
c--demi correspondance  for kbath, Bering Strait :
        kbath(iberp, jberp) = kbath(iberam,jberam)
        kbath(iberpm,jberp) = kbath(ibera, jberam)
      endif
 
c--6.3 Corners of the grid (pt. type velocity) - land/sea mask for velocity :
c--------------------------------------------------------------------------
 
c- Computation of the depths and computation of the index of beginning end of basin :
      do 635 j=js1,js2
       iu1(j) = imax
       iu2(j) = 1
       do 635 i=ims1,ims2
c--for each point (i,j) :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         km3 = min(kbath(i-1,j-1),kbath(i,j),kbath(i-1,j),kbath(i,j-1))
c--
c--Suppresion  of the island with no surface  :
         nn0vit = 1
         do 633 n=1,npt0v
           nn0vit = min( nn0vit , abs(i-ipt0v(n)) + abs(j-jpt0v(n)) )
 633     continue
c--
         if(km3.gt.0 .and. nn0vit.eq.1) then
           kkm3 = ks1 + ks2 - km3
           kniv(i,j,-1) = kkm3
           hu(i,j) = -zw(kkm3)
           unshu(i,j) = -unszw(kkm3)
           if(iu1(j).eq.imax) iu1(j) = i
           iu2(j) = i
           do 630 k=kkm3,kmax
             tmu(i,j,k) = 1.0
 630       continue
         endif
c--end of the work for each point (i,j).
 635  continue
 
c--6.4 definition of the  "openings" for  the W & N  boundaries (correspondance) :
c---------------------------------------------------------------------
 
      if (ltest.ge.1) then
      do 640 j=jcl1,jcl2
c--half cyclic correspondance  kbath :
       kbath(ims2+1,j) = kbath(ims1,j)
c--cyclic correspondance for kniv, hu, unshu, tms et tmu :
       kniv(ims1-1,j,0) = kniv(ims2,j,0)
       kniv(ims2+1,j,0) = kniv(ims1,j,0)
       kniv(ims1-1,j,-1) = kniv(ims2,j,-1)
       kniv(ims2+1,j,-1) = kniv(ims1,j,-1)
       hu(ims1-1,j) = hu(ims2,j)
       hu(ims2+1,j) = hu(ims1,j)
       unshu(ims1-1,j) = unshu(ims2,j)
       unshu(ims2+1,j) = unshu(ims1,j)
       do 640 k=ks1,ks2
        tms(ims1-1,j,k) = tms(ims2,j,k)
        tms(ims2+1,j,k) = tms(ims1,j,k)
        tmu(ims1-1,j,k) = tmu(ims2,j,k)
        tmu(ims2+1,j,k) = tmu(ims1,j,k)
 640  continue
      endif
c--
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Bering Strait, end of the closing :
      if (ltest.eq.3 .and. iberp.eq.0) then
c- if in the middle of 4 "land point"=> no connexion :
        do 647 j=2,jmax
          if (iberp.eq.0) then
            if (is1(j-1).gt.ims1 .and. is1(j).gt.ims1) then
              iberp = ims1
              jberp = j
            elseif (is2(j-1).lt.ims2 .and. is2(j).lt.ims2 ) then
              iberp = ims2
              jberp = j
            endif
          endif
 647    continue
        if (iberp.eq.0) then
          write(iuo+66,*) 'STOP in the routine "defgrid" :'
          write(iuo+66,*) 'Problem in Bering closing !'
          stop
        endif
        ibera = iberp
        jbera = jberp
        jberpm = jberp - 1
        iberpm = iberp - 1
        jberam = jbera - 1
        iberam = ibera - 1
      endif
      if (ltest.eq.3 .and. iberp.ne.ibera) then
c--half correspondance  for kbath, Bering Strait :
        kbath(ibera, jbera) = kbath(iberpm,jberpm)
        kbath(iberam,jbera) = kbath(iberp, jberpm)
c--correspondance  Bering for kniv, hu, unshu, tms et tmu :
        kniv(iberp, jberp,0) = kniv(iberam,jberam,0)
        kniv(iberpm,jberp,0) = kniv(ibera, jberam,0)
        kniv(ibera, jbera,0) = kniv(iberpm,jberpm,0)
        kniv(iberam,jbera,0) = kniv(iberp, jberpm,0)
        kniv(ibera,jbera,-1) = kniv(iberp,jberp,-1)
        hu(ibera,jbera) = hu(iberp,jberp)
        unshu(ibera,jbera) = unshu(iberp,jberp)
        do 650 k=ks1,ks2
          tms(iberp, jberp,k) = tms(iberam,jberam,k)
          tms(iberpm,jberp,k) = tms(ibera ,jberam,k)
          tms(ibera, jbera,k) = tms(iberpm,jberpm,k)
          tms(iberam,jbera,k) = tms(iberp ,jberpm,k)
          tmu(ibera,jbera,k) = tmu(iberp,jberp,k)
 650    continue
c--Correction of  "iu1(jberp),iu2(jberp)" (necessairy ?) :
C       j = jberp
C       do 651 i=iberp+1,imax
C       if (iu1(j).eq.iberp.and.tmu(i,j,ks2).eq.one) iu1(j) = i
C651    continue
C       do 652 i=iberp-1,1,-1
C       if (iu2(j).eq.iberp.and.tmu(i,j,ks2).eq.one) iu2(j) = i
C652    continue
c--Computation of  the coeff. "alpha" (-> bering = coeff for geostrophic control )
        bering = bering * gpes * hu(iberp,jberp)
     &                  * 0.5 / fs2cor(iberp,jberp)
        bering = max(zero, bering)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c--6.6 Corners of the grid (pt. type velocity) - mask sea+coast : kniv(1)
c--------------------------------------------------------------------------
 
c- computation of the maximum depth :
      do 660 j=js1,jmax
       do 660 i=ims1,imax
c--for each point (i,j) :
         km1 = max(kbath(i-1,j-1),kbath(i,j),kbath(i-1,j),kbath(i,j-1))
         kniv(i,j,1) = kmaxp1 - km1
 660  continue
c--
      if (ltest.ge.1) then
        do 662 j=jcl1,jcl2
c--cyclic boundary for kniv(1) :
          kniv(ims1-1,j,1) = kniv(ims2,j,1)
          kniv(ims2+1,j,1) = kniv(ims1,j,1)
 662    continue
      endif
      if (ltest.eq.3 .and. iberp.ne.ibera) then
c--correspondance  Bering for kniv(1) :
        kniv(ibera,jbera,1) = kniv(iberp,jberp,1)
        kniv(ibera+1,jbera,1) = kniv(iberp-1,jberp,1)
        kniv(ibera-1,jbera,1) = kniv(iberp+1,jberp,1)
        do 666 ii=-1,1
          kniv(iberp+ii,jberp+1,1) = kmaxp1
 666    continue
      endif
 
c--6.7 definition of the indexes for the computation of the fluxes :
c-------------------------------------------------------------------
 
c--initialisation of  isf1 and isf2 :
      do 670 j=js1,1+js2
       isf1(j) = min(is1(j-1),is1(j))
       isf2(j) = max(is2(j-1),is2(j))
 670  continue
 
c--initialisation of iuf1 and iuf2 , huy and hux :
      do 675 j=1,jmax - 1
        iuf1(j) = min(iu1(j),iu1(j+1))
        iuf2(j) = max(iu2(j),iu2(j+1))
        do 675 i=1,imax - 1
          huy(i,j) = min(hu(i,j),hu(i,j+1))
          hux(i,j) = min(hu(i,j),hu(i+1,j))
 675  continue
 
c--initialisation of imu1 et imu2 , ju1 et ju2 ,
c- correction (empty lign) of is2, iu2, isf2 & iuf2 :
      imu1 = imax
      imu2 = 1
      ju1 = jmax
      ju2 = 1
      do 677 j=1,jmax
        imu1 = min(imu1, iu1(j))
        imu2 = max(imu2, iu2(j))
        if (iu1(j) .gt.iu2(j) ) then
          iu2(j) = ims2
        else
          if (ju1.eq.jmax) ju1 = j
          ju2 = j
        endif
        if (is1(j) .gt.is2(j) ) is2(j) = ims2
        if (isf1(j).gt.isf2(j)) isf2(j) = ims2
        if (iuf1(j).gt.iuf2(j)) iuf2(j) = ims2
 677  continue
 
c--6.9 initialisation of the indexes indicating the bottom : kfs & kfu
c-----------------------------------------------------------------------
 
c--the level corresponding to the bottom : kfs (pts Sclal), kfu (pts Vites.)
      do 690 j=1,jmax
       do 690 i=1,imax
        kfs(i,j) = min( kmax, kniv(i,j,0) )
        kfu(i,j) = min( kmax, kniv(i,j,-1))
 690  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) Other Indices - Surface for Integration over the whole domain:        |
c-----------------------------------------------------------------------
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--7.1 Initialisation of the arrays indicating the corners :
 
      do 710 k=ks1,ks2
       n1coin(k) = 0
       n2coin(k) = 0
       n3coin(k) = 0
       n4coin(k) = 0
       do 710 j=js1,js2
        ij = (j-1) * imax
        do 710 i=is1(j),is2(j)
          xx4tms =  tms(i-1,j,k) + tms(i+1,j,k)
     &            + tms(i,j-1,k) + tms(i,j+1,k)
          if (xx4tms.eq.2.0d0 .and. tms(i,j,k).eq.one ) then
            if (tmu(i,j,k).eq.one) then
              n = 1 + n1coin(k)
              if (n.gt.ncomax) goto 1710
              i1coin(n,k) = i + ij
              n1coin(k) = n
            elseif (tmu(i+1,j,k).eq.one) then
              n = 1 + n2coin(k)
              if (n.gt.ncomax) goto 1710
              i2coin(n,k) = i + ij
              n2coin(k) = n
            elseif (tmu(i,j+1,k).eq.one) then
              n = 1 + n3coin(k)
              if (n.gt.ncomax) goto 1710
              i3coin(n,k) = i + ij
              n3coin(k) = n
            elseif (tmu(i+1,j+1,k).eq.one) then
              n = 1 + n4coin(k)
              if (n.gt.ncomax) goto 1710
              i4coin(n,k) = i + ij
              n4coin(k) = n
            endif
          endif
 710  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--7.3 Initialisation of the arrays indcating the slopes:
 
c- direction X :
      nl = 0
      do 730 k=ks1+1,ks2
       do 730 j=js1,js2
        do 730 i=is1(j),is2(j)+1
         if (tms(i-1,j,k).ne.one .or. tms(i,j,k).ne.one) goto 730
         if (tms(i-1,j,k-1).eq.one .and. tms(i,j,k-1).eq.zero) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730
           ijslp(nl) = i - 1 + (j-1) * imax
           kslp(nl) = k
           lslp(nl) = 1
         endif
         if (tms(i-1,j,k-1).eq.zero .and. tms(i,j,k-1).eq.one) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730
           ijslp(nl) = i + (j-1) * imax
           kslp(nl) = k
           lslp(nl) = -1
         endif
 730  continue
      nxslp = nl
 
c- direction Y :
      do 740 k=ks1+1,ks2
       do 740 j=js1,js2+1
        do 740 i=isf1(j),isf2(j)
         if (tms(i,j-1,k).ne.one .or. tms(i,j,k).ne.one) goto 740
         if (tms(i,j-1,k-1).eq.one .and. tms(i,j,k-1).eq.zero) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730
           ijslp(nl) = i + (j-2) * imax
           kslp(nl) = k
           lslp(nl) = imax
         endif
         if (tms(i,j-1,k-1).eq.zero .and. tms(i,j,k-1).eq.one) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730
           ijslp(nl) = i + (j-1) * imax
           kslp(nl) = k
           lslp(nl) = -imax
         endif
 740  continue
      nxyslp = nl
 
      nxslp0 = nxslp
      nyslp0 = nl - nxslp
      if (kfond.gt.-2) then
        nxslp  = 0
        nxyslp = 0
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--7.5 Initialisation of the Surfaces for the  Integration over the whole domain :
 
c--initialisation :
      do 750 llm=0,1
       do 750 k=1,kmax
        do 750 j=1,jmax
         do 750 i=1,imax
           ctmi(i,j,k,llm) = 0.
 750  continue
 
c--Point type Scalar :
      do 760 k=1,kmax
       do 760 j=js1,js2
        do 760 i=is1(j),is2(j)
          ctmi(i,j,k,0) = tms(i,j,k) * cmxy(i,j,0)
 760  continue
 
c--Point type velocity :
      do 770 k=1,kmax
       do 770 j=ju1,ju2
        do 770 i=iu1(j),iu2(j)
          ctmi(i,j,k,1) = tmu(i,j,k) * cmxy(i,j,3)
 770  continue
C     if (ltest.eq.3 .and. iberp.ne.ibera) then
c--Bering : ajoute 1 pt vitesse : <- si necessaire ?
C       i = iberp
C       j = jberp
C       do 775 k=1,kmax
C         ctmi(i,j,k,1) = tmu(i,j,k) * cmxy(i,j,3)
C775    continue
C     endif
 
c--computation of the ocean "volume" (in fact : volume /(dx*dy)) in metres :
      volz = 0.0
      do 780 j=js1,js2
        do 780 i=is1(j),is2(j)
          volz = volz + cmxy(i,j,0)*hs(i,j)
 780  continue
      unsvol = 1.0 / volz
c--end of the computation of the  volume.
      ccdxdy = dx * dy / 1.0d+12
      do 790 j=1,jmax
         do 790 i=1,imax
C           aire(i,j) = ctmi(i,j,ks2,0)*dx*dy/1.0d+12
            aire(i,j) = ctmi(i,j,ks2,0) * ccdxdy
790   continue

C surface for mean w 
        zurfow=0.0
        do j=js1,js2
         do i=is1(j),is2(j)
           zurfow=zurfow+aire(i,j)*tms(i,j,ks2)
         enddo
        enddo

 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  8 ) Checking - writing of parametres and indices of the domain .       |
c-----------------------------------------------------------------------
 
      if(nflag.ne.2) goto 899
c--checking
      write(99,*) 'dtb, dtu, dts(k) :'
      write(99,*)  dtb, dtu, dts
      write(99,*) 'dx, dy, pi :'
      write(99,*)  dx, dy, pi
      write(99,*) 'icheck, jcheck, kcheck :'
      write(99,*)  icheck, jcheck, kcheck
 
      if (kcheck.eq.0 .and. icheck.ge.1 .and. icheck.le.imax) then
c-----
c- Metric Coef. : writing for  i=icheck :
      do 820 jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(99,'(3(A,I3))') 'cmx(0), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(99,*) (cmx(icheck,j,0),j=jj1,jj2)
 820  continue
      do 821 jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(99,'(3(A,I3))') 'cmy(0), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(99,*) (cmy(icheck,j,0),j=jj1,jj2)
 821  continue
      do 826 jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(99,'(3(A,I3))') 'cmx(3), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(99,*) (cmx(icheck,j,3),j=jj1,jj2)
 826  continue
      do 827 jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(99,'(3(A,I3))') 'cmy(3), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(99,*) (cmy(icheck,j,3),j=jj1,jj2)
 827  continue
c-----
      endif
 
      if (kcheck.eq.0 .and. jcheck.ge.1 .and. jcheck.le.jmax) then
c-----
c- Metric coef. : wrting for j=jcheck :
      do 830 ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(99,'(3(A,I3))') 'cmx(0), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(99,*) (cmx(i,jcheck,0),i=ii1,ii2)
 830  continue
      do 831 ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(99,'(3(A,I3))') 'cmy(0), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(99,*) (cmy(i,jcheck,0),i=ii1,ii2)
 831  continue
 
      do 836 ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(99,'(3(A,I3))') 'cmx(3), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(99,*) (cmx(i,jcheck,3),i=ii1,ii2)
 836  continue
      do 837 ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(99,'(3(A,I3))') 'cmy(3), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(99,*) (cmy(i,jcheck,3),i=ii1,ii2)
 837  continue
c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      write(99,*) 'z :'
      write(99,*)  z
      write(99,*) 'dz :'
      write(99,*)  dz
      write(99,*) 'zw :'
      write(99,*)  zw
      write(99,*) 'dzw :'
      write(99,*)  dzw
 
      write(99,*) 'ahu, ahe :'
      write(99,*)  ahu, ahe
      write(99,*) 'ahs :'
      write(99,*)  ahs
      write(99,*) 'ai  :'
      write(99,*)  ai
      write(99,*) 'slopemax  :'
      write(99,*)  slopemax
      write(99,*) 'aitd  :'
      write(99,*)  aitd
      write(99,*) 'slopmgm  :'
      write(99,*)  slopmgm
      write(99,*) 'afilt,ahh,avv  :',afilt,ahh,avv
      write(99,*) 'avnub :'
      write(99,*)  avnub
      write(99,*) 'avnu0 :'
      write(99,*)  avnu0
      write(99,*) 'avkb :'
         write(99,*)  avkb
      write(99,*) 'avk0 :'
      write(99,*)  avk0
      write(99,*) 'rifumx, rifsmx :'
      write(99,*)  rifumx, rifsmx
      write(99,*) 'alphah :'
      write(99,*)  alphah
      write(99,*) 'alphgr :'
      write(99,*)  alphgr
      write(99,*) 'algrmn :'
      write(99,*)  algrmn
      write(99,*) 'alphmi :'
      write(99,*)  alphmi
      write(99,*) 'alphaz :'
      write(99,*)  alphaz
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      write(99,'(A)') 'Indices linked with the domain :'
      write(99,*) 'jsep(i) :'
      write(99,'(20I4)') (jsep(i),i=1,imax)
      write(99,*) 'jnorth(i) :'
      write(99,'(20I4)') (jnorth(i),i=1,imax)
 
      write(99,*) 'dxwi, dywj, xwi1, ywj1 :'
      write(99,*)  dxwi, dywj, xwi1, ywj1
      write(99,*) 'dxaj, dyai, xaj1, yai1, xwpoln :'
      write(99,*)  dxaj, dyai, xaj1, yai1, xwpoln
      write(99,*) 'jsepar, jeq, ijsdl, ijudl :'
      write(99,*)  jsepar,jeq,ijsdl,ijudl
      write(99,*) 'jcl1, jcl2, jdl1, jdl2 :'
      write(99,*)  jcl1, jcl2, jdl1, jdl2
c-----
      write(99,'(A,2I4)') 'Detroits : pour fichier "evolu" et "*.x" :'
     &                  //' nvhsf,ndhsf =', nvhsf, ndhsf
      do 855 nv=1,max(nvhsf,ndhsf)
        write(99,'(A,4I4)') tithsf(nv), ishsf(nv), iehsf(nv),
     &                      jshsf(nv), jehsf(nv)
 855  continue
      write(99,*)
c-----
      write(99,*) 'bering :', bering
      write(99,*) 'iberp, jberp, ibera, jbera :'
      write(99,*)  iberp, jberp, ibera, jbera
      write(99,*) 'cmx,cmy(iberp-1,jberp,2) :'
      write(99,*)  cmx(iberp-1,jberp,2), cmy(iberp-1,jberp,2)
      write(99,*) 'cmx,cmy(iberp,jberp,2) :'
      write(99,*)  cmx(iberp,jberp,2), cmy(iberp,jberp,2)
      write(99,*) 'cmx,cmy(ibera-1,jbera,2) :'
      write(99,*)  cmx(ibera-1,jbera,2), cmy(ibera-1,jbera,2)
      write(99,*) 'cmx,cmy(ibera,jbera,2) :'
      write(99,*)  cmx(ibera,jbera,2), cmy(ibera,jbera,2)
 
      do 891 n=1,npt0v
        write(99,*) 'Ajout ile surf=0 en (i,j)=', ipt0v(n), jpt0v(n)
 891  continue
      write(99,*) 'n1coin(k) :'
      write(99,'(20I4)') (n1coin(k),k=1,kmax)
      write(99,*) 'n2coin(k) :'
      write(99,'(20I4)') (n2coin(k),k=1,kmax)
      write(99,*) 'n3coin(k) :'
      write(99,'(20I4)') (n3coin(k),k=1,kmax)
      write(99,*) 'n4coin(k) :'
      write(99,'(20I4)') (n4coin(k),k=1,kmax)
C     do 892 k=1,kmax
C       write(99,*) 'Coins 1, k=', k
C       write(99,'(20I5)') (i1coin(n,k),n=1,n1coin(k))
C       write(99,*) 'Coins 2, k=', k
C       write(99,'(20I5)') (i2coin(n,k),n=1,n2coin(k))
C       write(99,*) 'Coins 3, k=', k
C       write(99,'(20I5)') (i3coin(n,k),n=1,n3coin(k))
C       write(99,*) 'Coins 4, k=', k
C       write(99,'(20I5)') (i4coin(n,k),n=1,n4coin(k))
C892  continue
      write(99,*) 'nxyslp, nXslope, nYslope, kfond :'
      write(99,*) nxyslp, nxslp0, nyslp0, kfond
 
      write(99,*) 'ju1, ju2, imu1, imu2 :'
      write(99,*)  ju1, ju2, imu1, imu2
      write(99,'(A)') 'Indices Debut(=1)/Fin(=2) :'
      write(99,'(A)') '  j | is1,is2   iu1,iu2  iuf1,iuf2'
     &          //' isf1,isf2     (iszon,iezon)(0/1/2/3) :'
      do 893 j=jmax,1,-1
        write(99,'(I3,1X,A1,4(2I4,2X),3(2I4,1X),2I4)') j, '|'
     &   , is1(j),  is2(j),  iu1(j),  iu2(j)
     &   , iuf1(j), iuf2(j), isf1(j), isf2(j)
     &   , (iszon(j,n),iezon(j,n),n=0,nbsmax)
 893  continue
 
c-- writing of the file  kbath (=sum[tm(k)]) :
c      <-- transfered in the routine  "local"
 
      write(99,*) 'array kfs(i,j) :'
      if (kmax.le.15) then
        nfrc = 125
        fmt1 = 'Z1'
      else
        nfrc = 41
        fmt1 = 'I3'
      endif
      write(fmtrep,'(A,I3,A)') '(',iabs(nfrc),'A1)'
      do 894 ii1=1,imax,nfrc
       ii2=min(ii1+nfrc-1,imax)
       write(fmt,'(A,I3,2A)') '(',(ii2-ii1+1),fmt1,',A1,I3)'
       write(99,*) 'portion de i=',ii1,' a i=',ii2
       write(99,fmtrep) (cc1(mod(i,10)),i=ii1,ii2)
       do 894 j=jmax,1,-1
       write(99,fmt) (kfs(i,j),i=ii1,ii2),'|',j
 894  continue
      write(99,*) 'volz, volume-real :'
      write(99,*)  volz, volz*dx*dy
 
 899  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  9 ) Needed in ice dynamic (comprise tms).                           |
c-----------------------------------------------------------------------
c
c--9.1. Ice thickness and advection masks.
c-------------------------------------------
c
      do 905 j=1,jmax
        do 905 i=1,imax
          zindfa(i,j) = tms(i,j,ks2)
 905  continue
 
c  Diffusion coefficients.
c
       ah = (uvdif/ren)*gridsz
       do 920 j=1,jmax-1
        jp1 = (j+1)
        do 910 i=1,imax
          ip1       = (i+1)-(imax-2)*(i/imax)
          dfhu(i,j) = tms(i,j,ks2)*tms(ip1,j,ks2)*ah
          dfhv(i,j) = tms(i,j,ks2)*tms(i,jp1,ks2)*ah
 910    continue
 920   continue
c     write (89,*) ah
c     write (89,*) dfhv
c
c--9.2. Determine cloud optical depths as a function of
c       latitude (CHOU ET AL., 1981).
c------------------------------------------------------
c
      do 930 j=js1,js2
         do 930 i=is1(j),is2(j)
         alat    = asin(covrai(i,j))/radian
         clat    = (95.0-alat)/10.0
         indx    = 1+int(clat)
         indxp1  = indx+1
         dl      = clat-int(clat)
         dr      = 1.0-dl
         tauc(i,j) = dr*tauco(indx)+dl*tauco(indxp1)
 930    continue
Cic0  return
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 10 ) Solar absorption at the various levels .                        |
c-----------------------------------------------------------------------
 
      open(11,file='typeaux.dat',status='old')
      do 1010 j=1,jmax
C        read(11,'(122(I3))') (ntypcn(i,j),i=1,imax)
         read(11,*) (ntypcn(i,j),i=1,imax)
 1010  continue
 
      data zrr /0.58,0.62,0.67,0.77,0.78/
      data zd1 /0.35,0.60,1.00,1.50,1.40/
      data zd2 /23.0,20.0,17.0,14.0,7.90/
      do 1020 j=1,jmax
        do 1020 i=1,imax
          rrro(i,j) = zrr(ntypcn(i,j))
          dd1o(i,j) = zd1(ntypcn(i,j))
          dd2o(i,j) = zd2(ntypcn(i,j))
 1020 continue
 
      do 1050 j=1,jmax
        do 1040 i=1,imax
          do k=1,kmax
            reslum(i,j,k) = 0.0
          enddo
          znivin = 0.0
          klim=max(11,kfs(i,j))
          do 1030 k=klim,kmax
            znivsp =fotr(rrro(i,j),dd1o(i,j),dd2o(i,j),-zw(k+1))
            reslum(i,j,k) = (znivsp - znivin)*tms(i,j,k)
     &                      *dts(k)/dts(ks2)*unsdz(k)/(rho0*cpo)
            znivin = znivsp
 1030     continue
 1040   continue
 1050 continue

cDFG start
c
c set horizontal tracer and momentum grid masks
c
      do j=1,jmax
        hs(1,j)=hs(imax-1,j)
        hs(imax,j)=hs(2,j)
        hu(1,j)=hu(imax-1,j)
        hu(imax,j)=hu(2,j)
        do i=1,imax
          msks(i,j)=1
          if (hs(i,j).eq.0) msks(i,j)=0
          msku(i,j)=1
          if (hu(i,j).eq.0) msku(i,j)=0
          angle(i,j)=angle(i,j)*float(msku(i,j))
        end do
      end do
cDFG end

C output of latitude and longitude
      open (61,file='lat.dat')
      do j=1,jmax
         write(61,'(122( F10.5))' ) (zlatt(i,j),i=1,imax)
      enddo
      write(61,*)
      close (61)
      open (61,file='long.dat')
      do j=1,jmax
         write(61,'(122( F10.5))' ) (zlont(i,j),i=1,imax)
      enddo
      write(61,*)
      close (61)

c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 11 ) Ice shelves and icebergs .                                         |
c-----------------------------------------------------------------------
c melting of iceberg in the North
c     
      iicebern1=98
      iicebern2=115
      jicebern1=46
      jicebern2=55
      areicebn=0.0
      do i=iicebern1,iicebern2
       do j=jicebern1,jicebern2
        areicebn=areicebn+area(i,j)*tms(i,j,ks2)
       enddo
      enddo

      iicebers1=2
      iicebers2=121
      jicebers1=1
      jicebers2=9
      areicebs=0.0
      do i=iicebers1,iicebers2
       do j=jicebers1,jicebers2
        areicebs=areicebs+area(i,j)*tms(i,j,ks2)
       enddo
      enddo

      write (99,*) 'surface iceberg melting N,S',areicebn,areicebs

C Definition of the ice-shelves
C     effecta=10E3
      effecta=5E3
      do i=1,imax
       do j=1,jmax
         tmics(i,j)=0.0
       enddo
      enddo
C 1. Rhonne Ice-shelf
      do i=93,103
        tmics(i,2)=tms(i,2,ks2)*dxs1(i,2)*effecta
      enddo
      do i=93,96
        tmics(i,3)=tms(i,3,ks2)*dxs1(i,2)*effecta
      enddo
C 2. East-Weddel 1
      do i=104,106
        tmics(i,3)=tms(i,3,ks2)*dxs1(i,3)*effecta
      enddo
      do i=107,109
        tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
      enddo
C 3. East-Weddel 2
      do i=110,114
        tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
      enddo
      do i=115,117
        tmics(i,5)=tms(i,5,ks2)*dxs1(i,5)*effecta
      enddo
C 4. Larsen
      do j=4,5
        tmics(93,j)=tms(93,j,ks2)*dxs2(93,j)*effecta
      enddo
C 5. Aimery
      do i=14,17
        tmics(i,5)=tms(i,5,ks2)*dxs1(i,5)*effecta
      enddo
C 6. Ross
      do i=48,61
        tmics(i,2)=tms(i,2,ks2)*dxs1(i,2)*effecta
      enddo
C 7. Getz
      do i=73,77
        tmics(i,3)=tms(i,3,ks2)*dxs1(i,3)*effecta
      enddo
C 8. George VI
      do i=84,87
        tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
      enddo

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 13 ) Pathological cases .                              |
c-----------------------------------------------------------------------
 
 1710 continue
      write(iuo+66,'(A,2I5)') 'Stop in routine "defgrid" :'
     &    //' ncomax too small ! (k,ncomax)= ', k, ncomax
      stop
 
 1730 continue
      write(iuo+66,'(A,4I5)') 'Stop in routine "defgrid" :'
     &    //' nlpmax too small ! (i,j,k,nlpmax)= ', i,j,k, ncomax
      stop
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- end of the routine defgrid -
      end
