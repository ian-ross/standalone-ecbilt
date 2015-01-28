













      subroutine cfc_flx(nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  This routine compute the air-sea exchange of 0 11 & 12.
c  (Rq : initialement incorpore a "icdyna.Fom" et separe le 08/11/98)
c  modif : 22/04/99
c  modif MFL : 10/6/2011
c---
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'comrunlabel.h'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'trace.com'
c
c c Also define Coefficients for Schmidt numbers (Wanninkhof, JGR, 1992)
      parameter (a11=4039.8, b11=264.7, c11=8.2552, d11=0.10359)
      parameter (a12=3713.2, b12=243.3, c12=7.5879, d12=0.095215)
 
c----------------------------------------------------------------------
c   Set up 0-fluxes from the atmosphere into the ocean, taking into
c  account the solubility of CFCs in seawater (Warner and Weiss, 1985).
c----------------------------------------------------------------------
 
      do 20 j=js1,js2
         do 10 i=is1(j),is2(j)
c Set integer variable to indicate whether SH or NH:
        if(j.lt.jeq) then
            ihemcfc = 1
        else
            ihemcfc = 2
        endif
 
c Set up 0-11 and 0-12 solubilities (tsol,ssol used for T,S
c       with T in degrees K, and S in parts per thousand):
 
C         mon = nint(tmonth)
C         tsol = sstobs(i,j,mon) + 273.0
C         ssol = (salobs(i,j,mon) * 1000.) + 35.0
 
          tsol = scal(i,j,ks2,1)
          ssol = scal(i,j,ks2,2)
 
          cfcsol11 = exp(a1cfc11 + a2cfc11*(100./tsol) +
     $          a3cfc11*log(tsol/100.) + a4cfc11*((tsol/100.)**2) +
     $          ssol*(b1cfc11 + b2cfc11*(tsol/100.) +
     $          b3cfc11*((tsol/100.)**2)) )
          cfcsol12 = exp(a1cfc12 + a2cfc12*(100./tsol) +
     $          a3cfc12*log(tsol/100.) + a4cfc12*((tsol/100.)**2) +
     $          ssol*(b1cfc12 + b2cfc12*(tsol/100.) +
     $          b3cfc12*((tsol/100.)**2)) )
 
c  Test output some 0 solubilities (near the dateline):
C         if(first .and. .not. mxpas2) then
C            if(i.eq.48 .and. j.eq.2) write(stdout,9088)
C9088         format(/,'Some examples of CFC solubilities (as a test):'
C    $,//,'   Temperature   Salinity   CFC-11*100  CFC-12*1000')
C            if((kmt(i,j).ge.1) .and. (i.eq.48))
C    $   write(stdout,9089) t(i,1,jc,nm,1), ssol,
C    $                          cfcsol11*100., cfcsol12*1000.
C         endif
C9089    format(2f11.3, 2f13.3)
 
C      if (i.eq.100.and.j.eq.50) then
C        write(120,*) tsol,ssol,cfcsol11,cfcsol12
C      endif
 
c Use a variety of techniques to get the ultimate "piston velocities":
 
 
c       cfcrestore = restore with e-folding time of 20 days (for example)
c               (i.e., surface level is restored towards atmospheric
c               concentration * solubility, factoring in a time scale
c               of 1/(20 days)
 
c  Value for the restoring
c      scalr(i,j,ks2,3) =cfcsol11*atmcfc11(ihemcfc)
c      scalr(i,j,ks2,4) =cfcsol12*atmcfc12(ihemcfc)
c      phiss(i,j,3) =0.0
c      phiss(i,j,4) =0.0
 
c       cfcwinds = Use complete Wanninkhof formulation with schmidt number
c                  and wind speed dependencies.  The wind speed dependent
c                  gas exchange is derived from E&K wind speed data (note
c                  interpolation over land and ice may give weird
c                  GFDL_prep_data values).
 
c       cfcliss = Use complete Liss and Merlivat (1986) formulation
c                  with schmidt number and wind speed dependencies.
 
c  Set up gammacfc to follow Wanninkhof (1992, JGR) formula:
        xsst = scal(i,j,ks2,1)-273.15
 
 
c  The Schmidt numbers (sc11, sc12) are slightly different (7% or so)
c    depending on which 0 is being considered:
        sc11 = a11 - b11*xsst + c11*(xsst**2) - d11*(xsst**3)
        sc12 = a12 - b12*xsst + c12*(xsst**2) - d12*(xsst**3)
 
C ifndef cfcliss
c  Wanninkhof's formula:
        xk11 = 0.31 * (vabq(i,j)**2) * sqrt(660/sc11)
        xk12 = 0.31 * (vabq(i,j)**2) * sqrt(660/sc12)
C endif
 
C ifdef cfcliss
c  Use Liss and Merlivat (1986) formulation instead ifdef cfcliss
C       if(vabq(i,j).le.3.6) xk11 = 0.17 * vabq(i,j)
C       if(vabq(i,j).gt.3.6 .and. vabq(i,j).le.13.)
C    $                  xk11 = 2.85 * vabq(i,j) - 9.65
C       if(vabq(i,j).gt.13.) xk11 = 5.9 * vabq(i,j) - 49.3
C       xk12 = xk11
C       if(vabq(i,j).le.3.6) then
C          xk11 = xk11 * ((sc11/600)**(-2./3.))
C          xk12 = xk12 * ((sc12/600)**(-2./3.))
C       else
C          xk11 = xk11 * ((sc11/600)**(-1./2.))
C          xk12 = xk12 * ((sc12/600)**(-1./2.))
C       endif
C endif
 
c Then the k (cm/hour) into time-scale (days) for level 1:
c Computation of the flux
C       t11(i,j) = (zw(1)/xk11) / 24.0
C       t12(i,j) = (zw(1)/xk12) / 24.0
        flcfc11  = xk11/(3600.*100.)*
     $    (cfcsol11*atmcfc11(ihemcfc) - scal(i,j,ks2,3))
        flcfc12  = xk12/(3600.*100.)*
     $    (cfcsol12*atmcfc12(ihemcfc) - scal(i,j,ks2,4))
        phiss(i,j,3) = - flcfc11 * dts(ks2) * unsdz(ks2)
        phiss(i,j,4) = - flcfc12 * dts(ks2) * unsdz(ks2)
c------------------------------------------------------------------------
c NB: Now limit the air-sea 0 flux in direct proportion to the
c     percentage coverage of sea-ice:
c------------------------------------------------------------------------
        phiss(i,j,3) = phiss(i,j,3) * albq(i,j)
        phiss(i,j,4) = phiss(i,j,4) * albq(i,j)
 
c  Then just add the 0 fluxes to the surface level 0-11 and 0-12....
 10      continue
 20   continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine cfc_flx -
      end
