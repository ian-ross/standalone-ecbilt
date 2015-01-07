c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine atmout(istep)
      implicit none

      include 'comglobal.h'

      integer istep

      if ( mod(istep,iatm) .eq. 0) then
c        if (iyear.eq.nyears) call outliu
        call selectout(istep)
      endif

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine selectout(istep)
c-----------------------------------------------------------------------
c *** this routine selects which kind of outputs should be written to output
c *** itel counts the number of days in the output interval
c *** meantype = 1: output monthly mean fields.
c *** meantype = 2: output seasonal mean fields.
c *** meantot = 1: computes whole period monthly or seasonal mean fields.
c *** meanyl  = 1: computes yearly monthly or seasonal mean fields.
c *** ioutdaily = 1 output instantanous fields.
c *** instcount: counter for output intantaneous fields.
c *** ixout: frequency for output instantanous fields in days.
c *** written by xueli wang.
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'
      include 'comdiag.h'
      include 'comoutlocal.h'

      integer i,j,k,l,istep
      integer idmean,idstd
      real*8  facstr

      ivlevel(1) = 200
      ivlevel(2) = 500
      ivlevel(3) = 800
      itlevel(1) = 350
      itlevel(2) = 650
      itlevel(3) = 1000

c** some computation and unit transformation

      call sptogg(psi(1,1),psig(1,1,1),pp)
      call sptogg(psi(1,2),psig(1,1,2),pp)
      call sptogg(psi(1,3),psig(1,1,3),pp)
      call sptogg(qprime(1,1),qgpv(1,1,1),pp)
      call sptogg(qprime(1,2),qgpv(1,1,2),pp)
      call sptogg(qprime(1,3),qgpv(1,1,3),pp)

c  *** compute the precipitation, evaporation and runoffs in cm/year
      facstr = roair*uv10rws
      do i=1,nlat
        do j=1,nlon
          dyrain1(i,j) = dyrain(i,j)*100.*3600.*24.*360.
          corain1(i,j) = corain(i,j)*100.*3600.*24.*360.
          torain1(i,j) = torain(i,j)*100.*3600.*24.*360.
          evap1(i,j)   = evap(i,j)*100.*3600.*24.*360.
          runofl1(i,j) = runofl(i,j)*100.*3600.*24.*360.
          runofo1(i,j) = runofo(i,j)*100.*3600.*24.*360.
          eminp1(i,j)  = evap1(i,j)-torain1(i,j)
          hesw(i,j)    = hesw1(i,j)+hesw2(i,j)+hesws(i,j)
          nlrads(i,j)  = ulrads(i,j)-dlrads(i,j)
          if (q0(i).gt.0d0) then
            albep(i,j)   = 1d0 - hesw(i,j)/q0(i)
          else
            albep(i,j)   = 1d0
          endif
          winstu1(i,j)=cdragw(i,j)*facstr*uvw10(i,j)*utot(i,j,3)
          winstv1(i,j)=cdragw(i,j)*facstr*uvw10(i,j)*vtot(i,j,3)
        enddo
      enddo
c  *** compute temperature in C
      do i=1,nlat
        do j=1,nlon
          tsurf1(i,j)=tsurf(i,j)-tzero
          temp4g1(i,j)=temp4g(i,j)-tzero
          temp2g1(i,j)=temp2g(i,j)-tzero
          tempsg1(i,j)=tempsg(i,j)-tzero
        enddo
      enddo

      if(meantype.eq.2) then
        if(istep.gt.(11*30*iatm)) then
          iseason = iseason + 1
          if (iseason.gt.90) iseason = 1
          itel = itel + 1
        endif
      else
        itel = itel + 1
      endif

      instcount = instcount + 1

c  *** write the instant data
      if (ioutdaily .eq. 1.and.instcount.eq.ixout) then
        call outputinst
      endif

      if (meantype .eq. 1) then
        minterv=30
        if (meantot .eq. 1) then
          call outputmtl
        endif
        if (meanyl .eq. 1) then
          call outputmyl
        endif
      endif

      if (meantype .eq. 2 .and. istep .gt. 11*30*iatm) then
        minterv=90
        if (meantot .eq. 1) then
          call outputmtl
        endif
        if (meanyl .eq. 1) then
          call outputmyl
        endif
      endif

      if (itel.eq.minterv) itel = 0
      if (instcount.eq.ixout) instcount = 0

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine outputinst
c-----------------------------------------------------------------------
c *** this routine output instantanous fields.
c *** written by xueli wang, september 1995.
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'
      include 'comcoup.h'
      include 'comdiag.h'
      include 'comoutlocal.h'


      integer i,j,k,l

c** write surface temperature, tempertature at 650, 350 mb and 10 meter
c** temperature.

      call wrout(139,0,tsurf1,0,newts(1))

      call wrout(130,1000,tempsg1,0,newt(1))
      call wrout(130,itlevel(2),temp4g1,0,newt(1))
      call wrout(130,itlevel(1),temp2g1,0,newt(1))

c** write wind components u, v

      call wrout(131,ivlevel(3),u800,0,newu(1))
      call wrout(131,ivlevel(2),u500,0,newu(1))
      call wrout(131,ivlevel(1),u200,0,newu(1))

      call wrout(132,ivlevel(3),v800,0,newv(1))
      call wrout(132,ivlevel(2),v500,0,newv(1))
      call wrout(132,ivlevel(1),v200,0,newv(1))

c** write surface pressure

      call wrout(134,0,pground,0,newsp(1))

c** write vertical velocity

      do l=nvl,1,-1
        call wrout(135,itlevel(l),omegg(1,1,l),0,newomega(1))
      enddo

c** write wind stress components u, v.

      call wrout(180,0,winstu1,0,newustress(1))
      call wrout(181,0,winstv1,0,newvstress(1))

c** write 10 meter wind magnitude.

      call wrout(159,0,uv10,0,newuv10(1))

c**  write ageostrophic wind u and v

      do l=nvl,1,-1
        call wrout(303,ivlevel(l),udivg(1,1,l),0,newdivu(1))
      enddo
      do l=nvl,1,-1
        call wrout(304,ivlevel(l),vdivg(1,1,l),0,newdivv(1))
      enddo

c** write stream function for 800 mb, 500 mb, 200 mb.

      do l=nvl,1,-1
        call wrout(148,ivlevel(l),psig(1,1,l),0,newpsi(1))
      enddo

c**  write velocity potential.

      do l=nvl,1,-1
        call wrout(149,ivlevel(l),chig(1,1,l),0,newchi(1))
      enddo

c**  write quasi geostrophic potential vorticity.

      do l=nvl,1,-1
        call wrout(310,ivlevel(l),qgpv(1,1,l),0,newqgpv(1))
      enddo

c**  write geopotential height.

      do l=nvl,1,-1
        call wrout(156,ivlevel(l),geopg(1,1,l),0,newgh(1))
      enddo

c** write heating force for 650, 350 mb.

      call wrout(301,itlevel(2),vhforg2,0,newvhforg(1))
      call wrout(301,itlevel(1),vhforg1,0,newvhforg(1))

c** write potential vorticity forceing for 800, 500, 200 mb.

      call wrout(302,ivlevel(3),vforg3,0,newvforg(1))
      call wrout(302,ivlevel(2),vforg2,0,newvforg(1))
      call wrout(302,ivlevel(1),vforg1,0,newvforg(1))

c** write dynamic rain in cm/year

      call wrout(142,0,dyrain1,0,newdyrain(1))

c** write convective rain cm/year

      call wrout(143,0,corain1,0,newcorain(1))

c** write total rain cm/year

      call wrout(260,0,torain1,0,newtorain(1))

c** write evaporation

      call wrout(182,0,evap1,0,newevap(1))

c** write evaporation minus percipitation.

      call wrout(309,0,eminp1,0,neweminp(1))

c** write sensible heat flux from the surface.

      call wrout(146,0,hflux,0,newhflux(1))

c** write latent heat flux due to surface evaporation and condensation.

      call wrout(147,0,eflux,0,neweflux(1))

C** Write top solar radiaton and surface solar radiation.

      call wrout(178,0,hesw,0,newtsr(1))
      call wrout(176,0,hesws,0,newssr(1))


C** Write surface thermal radiation and top thermal radiation.

      call wrout(177,0,nlrads,0,newstr(1))
      call wrout(179,0,ulrad1,0,newttr(1))

      call wrout(175,0,albes,0,newalbs(1))
      call wrout(174,0,albep,0,newalbp(1))

C** Write bottom moisture.

      call wrout(140,0,bmoisg,0,newbmoisg(1))

C** Write runoff over land and over ocean.

      call wrout(161,0,runofl1,0,newrunoffl(1))
      call wrout(160,0,runofo1,0,newrunoffo(1))

C** Write snow depth.

      call wrout(141,0,dsnow,0,newsdl(1))

C** Write specific humidity.

      call wrout(133,0,rmoisg(1,1),0,newrmoisg(1))

C** Write relative humidity.

      call wrout(157,0,relhum,0,newrelhum(1))

C** Write drag coefficient and richardson number

      call wrout(305,0,cdragw,0,newcdragw(1))
      call wrout(306,0,cdragv,0,newcdragv(1))
      call wrout(307,0,richar,0,newrichar(1))

C** write cloud cover

      call wrout(164,0,tcc,0,newtcc(1))

C** write dummy variable dumt1

      do l=nvl,1,-1
        call wrout(996,itlevel(l),dumt1(1,1,l),0,newdumt1(1))
      enddo

C** write dummy variable dumt2

      do l=nvl,1,-1
        call wrout(997,itlevel(l),dumt2(1,1,l),0,newdumt2(1))
      enddo

C** write dummy variable dumu1

      do l=nvl,1,-1
        call wrout(998,ivlevel(l),dumu1(1,1,l),0,newdumu1(1))
      enddo

C** write dummy variable dumu2

      do l=nvl,1,-1
        call wrout(999,ivlevel(l),dumu2(1,1,l),0,newdumu2(1))
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wrout(icode,ilevel,xxout,ifield,inx1)
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'
      include 'comdiag.h'

      integer  icode,ilevel,ifield,i,j,kunit,inx1
      real*8   xxout(nlat,nlon)
      real*4   yyout(nlat,nlon)

      if (inx1 .eq. 0 ) goto 100
      do i = 1, nlat
        do j = 1, nlon
          yyout(i,j) = xxout(i,j)
        enddo
      enddo

      write(21) icode,ilevel,iyear,imonth,iday,nlon,nlat,ifield
      write(21) ((yyout(i,j),j=1,nlon),i=1,nlat)

900   format(10e12.6)
910   format(8i10)

100   continue
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine outputmtl
c-----------------------------------------------------------------------
c *** this routine calls another routine which computes the  whole period
c *** monthly or seasonal mean fields.
c *** written by xueli wang and nanne weber, april 1995.
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'
      include 'comcoup.h'
      include 'comdiag.h'
      include 'comoutlocal.h'


      integer i,j,k,l

      k = imonth

c** write surface temperature, tempertature at 650, 350 mb and 10 meter
c** temperature.

      call tomean(139,0,tsurf1,sxtsurf,sytsurf,newts(2))
      call tomean(130,1000,tempsg1,sxtempsg,sytempsg,newt(2))
      call tomean(130,itlevel(2),temp4g1,sxtemp4g,sytemp4g,newt(2))
      call tomean(130,itlevel(1),temp2g1,sxtemp2g,sytemp2g,newt(2))

c** write wind components u, v

      call tomean(131,ivlevel(3),u800,sxu800,syu800,newu(2))
      call tomean(131,ivlevel(2),u500,sxu500,syu500,newu(2))
      call tomean(131,ivlevel(1),u200,sxu200,syu200,newu(2))

      call tomean(132,ivlevel(3),v800,sxv800,syv800,newv(2))
      call tomean(132,ivlevel(2),v500,sxv500,syv500,newv(2))
      call tomean(132,ivlevel(1),v200,sxv200,syv200,newv(2))

c**  write surface pressure

      call tomean(134,0,pground,sxpground,sypground,newsp(2))

c**  write vertical velocity

      call tomean(135,1000,omegg(1,1,3),sxomeg3,syomeg3,newomega(2))
      call tomean(135,650,omegg(1,1,2),sxomeg2,syomeg2,newomega(2))
      call tomean(135,350,omegg(1,1,1),sxomeg1,syomeg1,newomega(2))

c** write wind stress components.

      call tomean(180,0,winstu1,sxwinstu1,sywinstu1,newustress(2))
      call tomean(181,0,winstv1,sxwinstv1,sywinstv1,newvstress(2))

c** write magnitude of 10 meter height wind.

      call tomean(159,0,uv10,sxuv10,syuv10,newuv10(2))

c**  write ageostrophic wind u and v

      call tomean(303,ivlevel(3),udivg(1,1,3),sxudivg3,syudivg3,newdivu(2))
      call tomean(303,ivlevel(2),udivg(1,1,2),sxudivg2,syudivg2,newdivu(2))
      call tomean(303,ivlevel(1),udivg(1,1,1),sxudivg1,syudivg1,newdivu(2))
      call tomean(304,ivlevel(3),vdivg(1,1,3),sxvdivg3,syvdivg3,newdivv(2))
      call tomean(304,ivlevel(2),vdivg(1,1,2),sxvdivg2,syvdivg2,newdivv(2))
      call tomean(304,ivlevel(1),vdivg(1,1,1),sxvdivg1,syvdivg1,newdivv(2))

c** write stream function for 800 mb, 500 mb, 200 mb.

      call tomean(148,ivlevel(3),psig(1,1,3),sxgrpsi3,sygrpsi3,newpsi(2))
      call tomean(148,ivlevel(2),psig(1,1,2),sxgrpsi2,sygrpsi2,newpsi(2))
      call tomean(148,ivlevel(1),psig(1,1,1),sxgrpsi1,sygrpsi1,newpsi(2))

c**  write velocity potential.

      call tomean(149,ivlevel(3),chig(1,1,3),sxchi3,sychi3,newchi(2))
      call tomean(149,ivlevel(2),chig(1,1,2),sxchi2,sychi2,newchi(2))
      call tomean(149,ivlevel(1),chig(1,1,1),sxchi1,sychi1,newchi(2))

c**  write quasi geostrophic potential vorticity.

      call tomean(310,ivlevel(3),qgpv(1,1,3),sxqgpv3,syqgpv3,newqgpv(2))
      call tomean(310,ivlevel(2),qgpv(1,1,2),sxqgpv2,syqgpv2,newqgpv(2))
      call tomean(310,ivlevel(1),qgpv(1,1,1),sxqgpv1,syqgpv1,newqgpv(2))

c**  write geopotential height.

      call tomean(156,ivlevel(3),geopg(1,1,3),sxgh3,sygh3,newgh(2))
      call tomean(156,ivlevel(2),geopg(1,1,2),sxgh2,sygh2,newgh(2))
      call tomean(156,ivlevel(1),geopg(1,1,1),sxgh1,sygh1,newgh(2))

c** write heating force for 650, 350 mb.

      call tomean(301,itlevel(2),vhforg2,sxvhforg2,syvhforg2,newvhforg(2))
      call tomean(301,itlevel(1),vhforg1,sxvhforg1,syvhforg1,newvhforg(2))

c** write potential vorticity forceing for 800, 500, 200 mb.

      call tomean(302,ivlevel(3),vforg3,sxvforg3,syvforg3,newvforg(2))
      call tomean(302,ivlevel(2),vforg2,sxvforg2,syvforg2,newvforg(2))
      call tomean(302,ivlevel(1),vforg1,sxvforg1,syvforg1,newvforg(2))

c** write dynamic rain in cm/year

      call tomean(142,0,dyrain1,sxdyrain,sydyrain,newdyrain(2))

c** write convective rain cm/year

      call tomean(143,0,corain1,sxcorain,sycorain,newcorain(2))

c** write total rain cm/year

      call tomean(260,0,torain1,sxtorain,sytorain,newtorain(2))

c** write evaporation

      call tomean(182,0,evap1,sxevap,syevap,newevap(2))

c** write evaporation minus percipitation.

      call tomean(309,0,eminp1,sxeminp,syeminp,neweminp(2))

c** write sensible heat flux from the surface.

      call tomean(146,0,hflux,sxhflux,syhflux,newhflux(2))

c** write latent heat flux due to surface evaporation and condensation.

      call tomean(147,0,eflux,sxeflux,syeflux,neweflux(2))

C** Write surface solar radiation and top solar radiaton heating.

      call tomean(176,0,hesws,sxhesws,syhesws,newssr(2))
      call tomean(178,0,hesw,sxhesw,syhesw,newtsr(2))
      call tomean(175,0,albes,sxalbes,syalbes,newalbs(2))
      call tomean(174,0,albep,sxalbep,syalbep,newalbp(2))

C** Write surface thermal radiation and top thermal radiation.

      call tomean(177,0,nlrads,sxnlrads,synlrads,newstr(2))
      call tomean(179,0,ulrad1,sxulrad1,syulrad1,newttr(2))


C** Write bottom moisture.

      call tomean(140,0,bmoisg,sxbmoisg,sybmoisg,newbmoisg(2))

C** Write snow depth over land.

      call tomean(141,0,dsnow,sxdsnow,sydsnow,newsdl(2))

C** Write runoff over land and over sea.

      call tomean(161,0,runofl1,sxrunofl,syrunofl,newrunoffl(2))
      call tomean(160,0,runofo1,sxrunofo,syrunofo,newrunoffo(2))

C** Write specific humidity.

      call tomean(133,0,rmoisg(1,1),sxrmoisgw3,syrmoisgw3,newrmoisg(2))

C** Write relative humidity.

      call tomean(157,0,relhum,sxrelhum,syrelhum,newrelhum(2))

C** Write drag coefficient and richardson number

      call tomean(305,0,cdragw,sxcdragw,sycdragw,newcdragw(2))
      call tomean(306,0,cdragv,sxcdragv,sycdragv,newcdragv(2))
      call tomean(307,0,richar,sxrichar,syrichar,newrichar(2))

C** write cloud cover

      call tomean(164,0,tcc,sxtcc,sytcc,newtcc(2))

C** write dummy variable dumt1

      call tomean(996,itlevel(3),dumt1(1,1,3),sxdt13,sydt13,newdumt1(2))
      call tomean(996,itlevel(2),dumt1(1,1,2),sxdt12,sydt12,newdumt1(2))
      call tomean(996,itlevel(1),dumt1(1,1,1),sxdt11,sydt11,newdumt1(2))

C** write dummy variable dumt2

      call tomean(997,itlevel(3),dumt2(1,1,3),sxdt23,sydt23,newdumt2(2))
      call tomean(997,itlevel(2),dumt2(1,1,2),sxdt22,sydt22,newdumt2(2))
      call tomean(997,itlevel(1),dumt2(1,1,1),sxdt21,sydt21,newdumt2(2))

C** write dummy variable dumu1

      call tomean(998,ivlevel(3),dumu1(1,1,3),sxdu13,sydu13,newdumu1(2))
      call tomean(998,ivlevel(2),dumu1(1,1,2),sxdu12,sydu12,newdumu1(2))
      call tomean(998,ivlevel(1),dumu1(1,1,1),sxdu11,sydu11,newdumu1(2))

C** write dummy variable dumu2

      call tomean(999,ivlevel(3),dumu2(1,1,3),sxdu23,sydu23,newdumu2(2))
      call tomean(999,ivlevel(2),dumu2(1,1,2),sxdu22,sydu22,newdumu2(2))
      call tomean(999,ivlevel(1),dumu2(1,1,1),sxdu21,sydu21,newdumu2(2))

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine tomean(icode,ilevel,x,sumx1,sumy1,inx1)
c *** --------------------------------------------------------------------
c *** This routine computes whole period monthly or seasonal mean and standard
c *** deviation around it.
c *** ------------------------------------------------------------------------


      implicit none

      include 'comatm.h'
      include 'comcoup.h'
      include 'comglobal.h'

      integer i,j,ncase,k,m,inx1,nk
      integer  icode,ilevel,idmean,idstd,kunit
      real*8  sumx1(nlat,nlon,12),sumy1(nlat,nlon,12)
      real*8  x(nlat,nlon),xmean(nlat,nlon,12)
      real*8  xstd(nlat,nlon,12)
      real*4  xstd4(nlat,nlon,12),xmean4(nlat,nlon,12)

      if (inx1.eq.0) goto 100
      if (itel.eq.0) goto 100

      if(meantype.eq.1) then
        k = imonth
        ncase = nyears * 30
        nk = 12
        if(iyear.eq.1.and.imonth.eq.1.and.iday.eq.1) then
          do m = 1,nk
            do j = 1, nlon
              do i = 1, nlat
                sumx1(i,j,m)   = 0.0
                sumy1(i,j,m)   = 0.0
              enddo
            enddo
          enddo
        endif
      endif

      if(meantype.eq.2) then
        if(imonth.eq.12.or.imonth.eq.1.or.imonth.eq.2) k = 1
        if(imonth.ge.3.and.imonth.le.5)  k = 2
        if(imonth.ge.6.and.imonth.le.8)  k = 3
        if(imonth.ge.9.and.imonth.le.11) k = 4
        ncase = (nyears-1) * 90
        nk = 4
        if(iyear.eq.1.and.imonth.eq.12.and.iday.eq.1) then
          do m = 1,nk
            do j = 1, nlon
              do i = 1, nlat
                sumx1(i,j,m)   = 0.0
                sumy1(i,j,m)   = 0.0
              enddo
            enddo
          enddo
        endif
      endif

        do i = 1, nlat
          do j = 1, nlon
            sumx1(i,j,k) = sumx1(i,j,k) + x(i,j)
            sumy1(i,j,k) = sumy1(i,j,k) + x(i,j) ** 2
          enddo
        enddo

      do m = 1, nk
        do i = 1, nlat
          do j = 1, nlon
            xmean(i,j,m) = 0.0
            xstd (i,j,m) = 0.0
          enddo
        enddo
      enddo


      idmean = 1
      idstd  = 2

      if (meantype.eq.1) then
        if(iyear.eq.nyears.and.iday.eq.30) then
          do j = 1, nlon
            do i = 1, nlat
              xmean(i,j,k) = sumx1(i,j,k)/ncase
              xstd(i,j,k)  = sumy1(i,j,k) - sumx1(i,j,k)**2/ncase
              if (xstd(i,j,k) .le. 0.) xstd(i,j,k)=0.
              xstd(i,j,k) = sqrt (xstd(i,j,k)/(ncase - 1.0))
              xmean4(i,j,k) = xmean(i,j,k)
              xstd4(i,j,k) = xstd(i,j,k)
            enddo
          enddo

          kunit=22

          write(kunit) icode,ilevel,iyear,imonth,iday,nlon,nlat,idmean
          write(kunit) ((xmean4(i,j,k),j=1,nlon),i=1,nlat)
          if (inx1 .eq. 2) then
            write(kunit) icode,ilevel,iyear,imonth,iday,nlon,nlat,idstd
            write(kunit) ((xstd4(i,j,k),j=1,nlon),i=1,nlat)
          endif
        endif
      endif

      if (meantype.eq.2) then
        if(iyear.eq.nyears.and.iseason.eq.90) then
          do j = 1, nlon
            do i = 1, nlat
              xmean(i,j,k) = sumx1(i,j,k)/ncase
              xstd(i,j,k)  = sumy1(i,j,k) - sumx1(i,j,k)**2/ncase
              if (xstd(i,j,k) .le. 0.) xstd(i,j,k)=0.
              xstd(i,j,k) = sqrt (xstd(i,j,k)/(ncase - 1.0))
              xmean4(i,j,k) = xmean(i,j,k)
              xstd4(i,j,k) = xstd(i,j,k)
            enddo
          enddo

          kunit = 24

          write(kunit) icode,ilevel,iyear,imonth,iday,nlon,nlat,idmean
          write(kunit) ((xmean4(i,j,k),j=1,nlon),i=1,nlat)
          if (inx1 .eq. 2) then
            write(kunit) icode,ilevel,iyear,imonth,iday,nlon,nlat,idstd
            write(kunit) ((xstd4(i,j,k),j=1,nlon),i=1,nlat)
          endif
        endif
      endif


900   format(10e12.6)
910   format(9i10)

100   continue
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine outputmyl
c-----------------------------------------------------------------------
c *** this routine calls another routine which computes the yearly monthly
c *** mean fields.
c *** written by xueli wang, september 1995.
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comland.h'
      include 'comcoup.h'
      include 'comdiag.h'
      include 'comoutlocal.h'

      integer i,j,k,l

c** write surface temperature, tempertature at 650, 350 mb and 10 meter
c** temperature.

      call meanout(139,0,tsurf1,s1tsurf,s2tsurf,newts(3))
      call meanout(130,1000,tempsg1,s1tempsg,s2tempsg,newt(3))
      call meanout(130,itlevel(2),temp4g1,s1temp4g,s2temp4g,newt(3))
      call meanout(130,itlevel(1),temp2g1,s1temp2g,s2temp2g,newt(3))


c** write wind components u, v

      call meanout(131,ivlevel(3),u800,s1u800,s2u800,newu(3))
      call meanout(131,ivlevel(2),u500,s1u500,s2u500,newu(3))
      call meanout(131,ivlevel(1),u200,s1u200,s2u200,newu(3))

      call meanout(132,ivlevel(3),v800,s1v800,s2v800,newv(3))
      call meanout(132,ivlevel(2),v500,s1v500,s2v500,newv(3))
      call meanout(132,ivlevel(1),v200,s1v200,s2v200,newv(3))

c** write surface pressure

      call meanout(134,0,pground,s1pground,s2pground,newsp(3))

c** write vertical velocity

      call meanout(135,1000,omegg(1,1,3),s1omeg(1,1,3),
     *               s2omeg(1,1,3),newomega(3))
      call meanout(135,650,omegg(1,1,2),s1omeg(1,1,2),
     *               s2omeg(1,1,2),newomega(3))
      call meanout(135,350,omegg(1,1,1),s1omeg(1,1,1),
     *               s2omeg(1,1,1),newomega(3))

c** write u, v wind stress.

      call meanout(180,0,winstu1,s1winstu1,s2winstu1,newustress(3))
      call meanout(181,0,winstv1,s1winstv1,s2winstv1,newvstress(3))

c** write magnitude of 10 meter height wind.

      call meanout(159,0,uv10,s1uv10,s2uv10,newuv10(3))

c**  write ageostrophic wind u and v

      do l=nvl,1,-1
         call meanout(303,ivlevel(l),udivg(1,1,l),s1udivg(1,1,l),
     *                s2udivg(1,1,l),newdivu(3))
      enddo
      do l=nvl,1,-1
         call meanout(304,ivlevel(l),vdivg(1,1,l),s1vdivg(1,1,l),
     *                s2vdivg(1,1,l),newdivv(3))
      enddo

c** write stream function for 800 mb, 500 mb, 200 mb.

      do l=nvl,1,-1
         call meanout(148,ivlevel(l),psig(1,1,l),s1psi(1,1,l),
     *                s2psi(1,1,l),newpsi(3))
      enddo

c**  write velocity potential.

      do l=nvl,1,-1
         call meanout(149,ivlevel(l),chig(1,1,l),s1chi(1,1,l),
     *                s2chi(1,1,l),newchi(3))
      enddo

c**  write quasi geostrophic potential vorticity.

      do l=nvl,1,-1
         call meanout(310,ivlevel(l),qgpv(1,1,l),s1qgpv(1,1,l),
     *                s2qgpv(1,1,l),newqgpv(3))
      enddo

c**  write geopotential height.

      do l=nvl,1,-1
         call meanout(156,ivlevel(l),geopg(1,1,l),s1gh(1,1,l),
     *                s2gh(1,1,l),newgh(3))
      enddo

c** write heating force for 650, 350 mb.

      call meanout(301,itlevel(2),vhforg2,s1vhforg2,s2vhforg2,newvhforg(3))
      call meanout(301,itlevel(1),vhforg1,s1vhforg1,s2vhforg1,newvhforg(3))

c** write potential vorticity forcing for 800, 500, 200 mb.

      call meanout(302,ivlevel(3),vforg3,s1vforg3,s2vforg3,newvforg(3))
      call meanout(302,ivlevel(2),vforg2,s1vforg2,s2vforg2,newvforg(3))
      call meanout(302,ivlevel(1),vforg1,s1vforg1,s2vforg1,newvforg(3))

c** write dynamic rain in cm/year

      call meanout(142,0,dyrain1,s1dyrain,s2dyrain,newdyrain(3))

c** write convective rain cm/year

      call meanout(143,0,corain1,s1corain,s2corain,newcorain(3))

c** write total rain cm/year

      call meanout(260,0,torain1,s1torain,s2torain,newtorain(3))

c** write evaporation

      call meanout(182,0,evap1,s1evap,s2evap,newevap(3))

c** write evaporation minus percipitation.

      call meanout(309,0,eminp1,s1eminp,s2eminp,neweminp(3))

c** write sensible heat flux from the surface.

      call meanout(146,0,hflux,s1hflux,s2hflux,newhflux(3))

c** write latent heat flux due to surface evaporation and condensation.

      call meanout(147,0,eflux,s1eflux,s2eflux,neweflux(3))

C** Write surface solar radiaton and to solar radiation.

       call meanout(176,0,hesws,s1hesws,s2hesws,newssr(3))
       call meanout(178,0,hesw,s1hesw,s2hesw,newtsr(3))

C** Write top thermal radiation and surface thermal radiation.
       call meanout(179,ivlevel(1),ulrad1,s1ulrad1,s2ulrad1,newttr(3))
       call meanout(177,0,nlrads,s1nlrads,s2nlrads,newstr(3))

       call meanout(175,0,albes,s1albes,s2albes,newalbs(3))
       call meanout(174,0,albep,s1albep,s2albep,newalbp(3))

C** Write bottom moisture.

       call meanout(140,0,bmoisg,s1bmoisg,s2bmoisg,newbmoisg(3))

C** Write snow depth over land.

       call meanout(141,0,dsnow,s1dsnow,s2dsnow,newsdl(3))

C** Write runoff over land and over ocean.

       call meanout(161,0,runofl1,s1runofl,s2runofl,newrunoffl(3))
       call meanout(160,0,runofo1,s1runofo,s2runofo,newrunoffo(3))

C** Write specific humidity.

       call meanout(133,0,rmoisg(1,1),s1rmoisgw3,s2rmoisgw3,newrmoisg(3))

C** Write relative humidity.

       call meanout(157,0,relhum,s1relhum,s2relhum,newrelhum(3))

C** Write drag coefficient and richardson number

       call meanout(305,0,cdragw,s1cdragw,s2cdragw,newcdragw(3))
       call meanout(306,0,cdragv,s1cdragv,s2cdragv,newcdragv(3))
       call meanout(307,0,richar,s1richar,s2richar,newrichar(3))

C** write cloud cover

       call meanout(164,0,tcc,s1tcc,s2tcc,newtcc(3))


C** write dummy variable dumt1

      do l=nvl,1,-1
        call meanout(996,itlevel(l),dumt1(1,1,l),
     *     s1dt1(1,1,l),s2dt1(1,1,l),newdumt1(2))
      enddo

C** write dummy variable dumt2

      do l=nvl,1,-1
        call meanout(997,itlevel(l),dumt1(1,1,l),
     *     s1dt1(1,1,l),s2dt1(1,1,l),newdumt1(2))
      enddo

C** write dummy variable dumu1

      do l=nvl,1,-1
        call meanout(998,ivlevel(l),dumu1(1,1,l),
     *     s1du1(1,1,l),s2du1(1,1,l),newdumu1(2))
      enddo

C** write dummy variable dumu2

      do l=nvl,1,-1
        call meanout(999,ivlevel(l),dumu2(1,1,l),
     *     s1du2(1,1,l),s2du2(1,1,l),newdumu2(2))
      enddo

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine meanout(icode,ilevel,xx,sumxx,sumxxsq,inx1)
c-----------------------------------------------------------------------
c *** this routine computes mean xxm over the output interval minterv
c *** (in days),and standard deviation xxdev for the variable xx.
c *** itel counts the number of days in the output interval
c *** written by xueli wang and nanne weber, april 1995.
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comglobal.h'

      integer i,j,inx1
      integer icode,ilevel,idmean,idstd,kunit

      real*8  sumxx(nlat,nlon),sumxxsq(nlat,nlon)
      real*8  xx(nlat,nlon),xxm(nlat,nlon)
      real*8  xxdev(nlat,nlon)
      real*4  xxdev4(nlat,nlon),xxm4(nlat,nlon)

      if (inx1 .eq. 0 ) goto 100
      if (itel .eq. 0 ) goto 100

      if(itel .eq. 1)  then
        do j = 1, nlon
          do i = 1, nlat
            sumxx(i,j)   = 0.0
            sumxxsq(i,j)   = 0.0
          enddo
        enddo
      endif

      do j = 1, nlon
        do i = 1, nlat
          sumxx(i,j)   = sumxx(i,j) + xx(i,j)
          sumxxsq(i,j) = sumxxsq(i,j) + xx(i,j) ** 2
        enddo
      enddo

      do j = 1, nlon
        do i = 1, nlat
          xxm(i,j)   = 0.0
          xxdev(i,j) = 0.0
        enddo
      enddo

      if (itel .eq. minterv) then
        do j = 1, nlon
          do i = 1, nlat
            xxm(i,j)   = sumxx(i,j)/minterv
            xxdev(i,j) = sumxxsq(i,j) - sumxx(i,j)**2/minterv
            if (xxdev(i,j) .le. 0.) xxdev(i,j)=0.
            xxdev(i,j) = sqrt (xxdev(i,j)/(minterv - 1.0))
            xxm4(i,j) = xxm(i,j)
            xxdev4(i,j) = xxdev(i,j)
          enddo
        enddo

        idmean=1
        idstd =2

        if (meantype .eq. 1) then
          kunit=23
        endif
        if (meantype .eq. 2) then
          kunit=25
        endif

        write(kunit) icode,ilevel,iyear,imonth,iday,nlon,nlat,idmean
        write(kunit) ((xxm4(i,j),j=1,nlon),i=1,nlat)
        if (inx1 .eq. 2) then
          write(kunit) icode,ilevel,iyear,imonth,iday,nlon,nlat,idstd
          write(kunit) ((xxdev4(i,j),j=1,nlon),i=1,nlat)
        endif
      endif

900   format(10e12.6)
910   format(9i10)


100   continue
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine outliu

      implicit none

      include 'comatm.h'
      include 'comdyn.h'

      integer l,k

      do l=1,nvl
        write(26) (psi(k,l),k=1,nsh2)
      enddo
      return
      end
