












      subroutine streaminit
c
c A)
c i1(2)zon arrays determining basins in model grid coordinates
c South of Bering Strait: 3 zones 
c nb=1: Indian, nb=2: Pacific, nb=3: Atlantic, nb=0: Global Ocean
c North of Bering Strait included in Atlantic part
c B)
c reading data file: jmaxl.dat ==> declared in 'reper.com'
c connects the model j-coordinates to true lattitudes and is used in computing 
c contributions to streamfunctions for North Atlantic and Arctic Ocean: 
c v.dx_(model j) and/or u.dy_(model j) ==> v_eff.dy_(true latitude)
c this part is independant of land-sea masks !
c
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'reper.com'


      integer i1zon(jmax,0:nbsmax),i2zon(jmax,0:nbsmax)
      common/streamzon/ i1zon,i2zon


      character*30 fmtl
      character*70 line


c A) --------------------- A)
c code taken from geogra.f for defining borders of basins i1zon/i2zon
c may be different from arrays iszon and iezon of geogra.f !
      dlat  = 3.0
      dlong = 3.0
      xlon1 =  25.5
      ylat1 = -79.5
      do j=1,jmax
        i1zon(j,0) = 2
        i2zon(j,0) = imax - 1
        i2zon(j,nbsmax) = imax - 1
        do nb=1,nbsmax-1
          i2zon(j,nb) = 1
        enddo
      enddo
      ybering = 67.0
      jbering = nint( (ybering - ylat1) / dlat ) + 1
      do j=js1,jbering
        yy = ylat1 + dlat * DFLOAT(j-1)
c- border Indian/Pacific: xx = f(yy) ; then i2zon(-,1)
        if (yy.le.-31.0) then
          xx = 143.5
        elseif (yy.le.-6.0) then
          xx = 112.5 - yy
        elseif (yy.le.9.0) then
          xx = 103.5 - 0.25 * yy
        else
          xx = 105.0
        endif
        i2zon(j,1) = nint( (xx-xlon1) / dlong - 0.5 ) + 1
c- border Pacific/Atlantic: i2zon(-,2)
        if (yy.le.-60.0) then
          xx = 297.0
        elseif (yy.le.-54.0) then
          xx = 237.0 - yy
        elseif (yy.le.-42.0) then
          xx = 291.0
        elseif (yy.le.-36.0) then
          xx = 333.0 + yy
        elseif (yy.le.3.0) then
          xx = 297.0
        elseif (yy.le.22.5) then
          xx = 303.0 - yy - yy
        else
          xx = 291.0
          xx = 258.0
        endif
        i2zon(j,2) = nint( (xx-xlon1) / dlong - 0.5 ) + 1
      enddo
      do j=1,jmax
        i1zon(j,1)  = 2
        do nb=2,nbsmax
          i1zon(j,nb) = i2zon(j,nb-1) + 1
        enddo
      enddo
c B) --------------------- B)
c read file : jmaxl.dat
      do j=1,jmax
       do i=1,imax
        jmaxl(i,j)=j
       enddo
      enddo
      jmvlat = 0



      open(52, file='jmaxl.dat', status='old', err=480)
      read(52,'(2A)',err=475) fmtl, line
      read(line,*) spvl, iml, jml, kml, nfrcl
      read(52,*)
      read(52,*)
      read(52,'(A)',err=475) ttvlat
c      write(99,*)'File jmaxl.dat'
c      write(99,*) spvl, iml, jml, kml, nfrcl
      jj1 = max(1,-jml)
      jj2 = max(jml,1)
      jj3 = sign(1,jml)
      do j=jj1,jj2,jj3
        read(52,fmtl,err=475) (jmaxl(i,j),i=1,iml)
c        write(99,fmtl,err=475) (jmaxl(i,j),i=1,iml)
      enddo
      close(52)
      jmvlat = abs(jml)
      write(99,*)'File jmaxl.dat; jmvlat =', jmvlat
      return
 475  continue
      STOP ' File jmaxl.dat not properly read '
 480  continue
      STOP ' File jmaxl.dat not properly opened ' 
      return     
      end
c------------------------------------------------------------
      subroutine streaminitout(nstreamout)
c
c (A) write GrADS str*.ctl files:
c stream.ctl; heatsaltflux.ctl
c
c (B) opening outputfiles:
c unit = 38: stream.dat for streamfunction meridional velocity v         
c unit = 39:  heatsaltflux.dat for meridional heat flux and salt flux
c
      include 'para.com'
      include 'reper.com'

      integer end_year,end_day

      common /ec_timepar/ nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,iobclin,
     &                     nwrskip,nwrskip_days,end_year,end_day

c 
      parameter (jmtt=57,kmtt=kmax+1)
      real*4 uuu(jmtt,0:kmtt,0:nbsmax)
      real*4 vvv(jmtt,0:nbsmax+4)
      common/streamflu/ uuu,vvv



c A) ------------------------ c write str*.ctl for GrADS



c) for meridional overturning + (stream.ctl; stream.dat)
      open(38,file='stream.ctl')
      write(38,fmt="('dset   ^stream.dat')")    
      write(38,fmt="('options big_endian')")
      write(38,fmt="('undef ',1p,e12.4)") -1.0e20
      write(38,fmt="('title Clio streamfunctions for merid. overturning')")
      write(38,fmt="('xdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(38,fmt="('ydef ',i3,' linear ',2f7.2)") 57,-81.0,3.0
      write(38,fmt="('zdef ',i3,' levels')") kmtt+1
      write(38,fmt="(' 5500     5126.18 4385.22 3661.11 2963.25 2307.36')") 
      write(38,fmt="(' 1717.90  1225.11  850.19  588.88  415.07  299.29')") 
      write(38,fmt="('  219.86   163.28  121.52   89.75   64.96   45.20')")
      write(38,fmt="('   29.17    15.98    5.00    0.00')")
      if (nstreamout.eq.1) then
       write(38,fmt="('tdef ',i5,' linear 1jan0001  30dy')") nyears*12
      else
       write(38,fmt="('tdef ',i5,' linear 1jan0001  1YR')") nyears
      endif
      write(38,fmt="('vars 4')") 
      write(38,fmt="('psig       22  99 Streamfunction Global   Sv')")
      write(38,fmt="('psii       22  99 Streamfunction Indian   Sv')")
      write(38,fmt="('psip       22  99 Streamfunction Pacific  Sv')")
      write(38,fmt="('psia       22  99 Streamfunction Atlantic Sv')")
      write(38,fmt="('endvars')")       
      close(38)
c) for heat flux and salt flux  + (heatsaltflux.ctl; heatsaltflux.dat)
      open(38,file='heatsaltflux.ctl')
      write(38,fmt="('dset   ^heatsaltflux.dat')")    
      write(38,fmt="('options big_endian')")
      write(38,fmt="('undef ',1p,e12.4)") -1.0e20
      write(38,fmt="('title Clio streamfunctions for heat flux')")
      write(38,fmt="('xdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(38,fmt="('ydef ',i3,' linear ',2f7.2)") 57,-81.0,3.0
      write(38,fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      if (nstreamout.eq.1) then
       write(38,fmt="('tdef ',i4,' linear 1jan0001  30dy')") nyears*12
      else
       write(38,fmt="('tdef ',i4,' linear 1jan0001  1YR')") nyears
      endif
      write(38,fmt="('vars 8')") 
      write(38,fmt="('hg       1  99 Heat Flux Global   PW')")
      write(38,fmt="('hi       1  99 Heat Flux Indian   PW')")
      write(38,fmt="('hp       1  99 Heat Flux Pacific  PW')")
      write(38,fmt="('ha       1  99 Heat Flux Atlantic PW')")
      write(38,fmt="('sg       1  99 Salt Flux Global   psuSv')")
      write(38,fmt="('si       1  99 Salt Flux Indian   psuSv')")
      write(38,fmt="('sp       1  99 Salt Flux Pacific  psuSv')")
      write(38,fmt="('sa       1  99 Salt Flux Atlantic psuSv')")
      write(38,fmt="('endvars')")       
      close(38)


c B) ------------------------ c opening outputfiles 



      open(38,file='stream.dat',form='unformatted',access='direct',
     &         recl=Size(uuu)*Kind(uuu(1,0,0)))
      open(39,file='heatsaltflux.dat',form='unformatted',access='direct',
     &         recl=Size(vvv)*Kind(vvv(1,0)))


      return     
      end
c------------------------------------------------------------
