












      subroutine veget_wr(kveget,nyvegt, nwrveg,iyr0vg, iyear, nyears,
     &             temveg,gd0veg,prcveg,tpsdry,
     &             prcmin,bmtdry, epss,fracgr,darea, 
     &             titveg)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c Write vegetation output on ASCI formated files "veget.zav" & "veget.outp" :
c       zonal+global mean: every "nyvegt" (=every call veget)
c       2.D maps: every "nwrveg" year(s).
c----
c  iyear = 0 : (=1rst call) open & write Header (2 lines) of output files.
c            and if kveget < 0 write Equilibr. Vegetation on output files.
c  WARNING : put Special Values in array:
c                          st,sg,sd,temveg,gd0veg,prcveg,tpsdry,snlt,blai
c  use file_unit iveg+7 & iveg+8 (remains opened during the run).
c---------
c  modif : 01/10/00
 
      USE Vegetation_Output
      implicit none


      include 'veget.h'
      real*8 forestfr(nlat,nlon),albsnow(nlat),albland(nlat,nlon,4)
      real*8 alblbm(nlat,nlon)
      real*8 moc,tmc,tmc0,tsurfmean,cland,thex 
      real*8 alblandR(nlat,nlon,4),forestfrR(nlat,nlon),alblbmR(nlat,nlon)
      real*8 alblandismR(nlat,nlon,4,3) 
      common /ec_lbmcalbedo/albsnow,albland,forestfr,alblbm
      common /ec_lbmcalbedo2/alblandR,forestfrR,alblbmR,alblandismR
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex

 
c--dummy variables :
c- input :
      integer kveget, nyvegt, nwrveg, iyr0vg, iyear, nyears
      real*8 temveg(nlat,nlon), gd0veg(nlat,nlon),
     &     prcveg(nlat,nlon,2), tpsdry(nlat,nlon)
      real*8 epss, prcmin, bmtdry
      real*8 fracgr(nlat,nlon), darea(nlat)
      character*15 titveg
c- output :
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--local variables saved (in common) from one call to the other :
      integer nnctr
      common / cm0wvg / nnctr
 
      integer nvgmax, kwradd
      parameter ( nvgmax = 19 , kwradd = 0 )
      real*8 vegsum, vegmap
      common / cm1wvg / vegsum, vegmap(nlat,nlon,nvgmax),vegsum2
 
      logical flgveg,flgicb,flgisma,flgismg
      common /ec_coupl/ flgveg,flgicb,flgisma,flgismg
c--data + local variables :
      integer i,j,k,n, kk,ns, nbwr1,nbwr2
      real*8 zero, one, var(nlon)
      real*8 dzero, zavsum, totsum,vegsum2
      real*8 vegzav(0:nlat,nvgmax), cfmzav(nvgmax), vegspv, ttscun, tts1
      character*6 ccny
      character*15 titzav(nvgmax)
      character*30 titmap(nvgmax)
      character*32 fmtveg, fmtzav
      real*8 soiltype2(nlat,nlon) 
 
c--For output :
      data vegspv / 999.0 /
      data fmtveg / '(65E12.5)                      ' /
      data fmtzav / '(65E12.5)                      ' /
      data cfmzav / 4*100. , 2*1. , 25. , 1. , 3*0.001, 7*1E-12,  100/
      data titzav / 'tree o/o', 'grass   ', 'desert  ', 't.needle',
     &              't_lai   ', 'g_lai   ', 'albe o/o',
     &              'temp oC ', 'gdd0 e-3', 'prc m/y ', 'prc_veg ', 
     &              'uptake GtC/y ', 'stock GtC' , 'b1 GtC' ,
     &              ' b2 GtC' , 'b3 GtC' ,' b4 GtC' ,'NPP GtC/y ',
     &               'icesheet o/o'  /
      data titmap / 
     & 'Fraction of Tree   (o/o) ; ' , 'Fraction of Grass  (o/o) ; ' ,
     & 'Fraction of Desert (o/o) ; ' , 'Fract. of Needle L.(o/o) ; ' ,
     & ' L.A.I. for Tree    (1)  ; ' , ' L.A.I. for Grass   (1)  ; ' ,
     & 'An. Albedo, snow=0 (o/o) ; ' ,
     & ' Annual Temperature (oC) ; ' , ' Annual G.D.D.0  (/1000) ; ' ,
     & ' Annual Precipit.  (m/y) ; ' , 'Reduced Prc. (Min&Dry_S) ; ' ,
     & ' Annual C Uptake (kgC/y/m2) ; ','Carbon Stock (kgC/m2) ; ' ,
     & ' B1 (kgC/m2) ; ' , ' B2 (kgC/m2) ; ' , 'B3 (kgC/m2) ; ',
     & 'B4 (kgC/m2) ; ' ,
     & ' NPP (kgC/y/m2) ','icesheet fraction o/o' /
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 1112 format(3(F7.2,1X,F7.3,1X),I3,A)
 
      zero = 0.
      one  = 1.
      dzero = 0.d0
c-----
 
      if (iyear.eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) 1rst call : Open & Write Header of veget. output files :
c-----------------------------------------------------------------------
 
        nnctr = 0
        do 100 n=1,len(titveg)
          if (titveg(n:n).ne.' ') nnctr = n
 100    continue
        vegsum = 0.
        vegsum2 = 0.
        vegmap(:,:,:)= 0.
        st_moy(:,:)= 0. 
 
c- time scale unit of "veget.zav" Header = kyr
        ttscun = one/1000.
        tts1 = (nyvegt+iyr0vg)*ttscun
        nbwr1 = nyears / nyvegt
        nbwr2 = nyears / nwrveg
        if (kveget.eq.0) then
c--No vegetation : write Imposed Tree Fraction & Albedo (<- from EcBilt) :
          titmap(1) = 'EcBilt Origin. Tree Frac.; '
          nbwr2 = (5-nvgmax)*nbwr2
        elseif (kveget.lt.0) then
          tts1 = iyr0vg*ttscun
          nbwr1 = nbwr1 + 1
          nbwr2 = -nvgmax*(1+nbwr2)
        else
          nbwr2 = -nvgmax*nbwr2
        endif
 
        open(iveg+7,file='veget.zav',status='unknown')
        write(iveg+7,1000) fmtzav, vegspv, 1+nlat, nvgmax, nbwr1, 65
        write(iveg+7,1112) -87.19-5.625, 5.625, 8., 0.,
     &                    tts1, nyvegt*ttscun, -5
 
        open(iveg+8,file='veget.outp',status='unknown')
        write(iveg+8,1000) fmtveg, vegspv, nlon, nlat, nbwr2, 65
        write(iveg+8,1111) 0., 5.625, -87.19, 5.625, 1., 1., 0
 
        if (kveget.ge.0) return
c-------
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      if (kveget.eq.0) then
c--No vegetation : fill in "st" with EcBilt "Original" Tree Fraction  :
        do 150 j=1,nlon
         do 150 i=1,nlat
          if (fracgr(i,j).gt.epss) st(i,j) = forestfr(i,j)
 150    continue
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Compute & Write Global + Zonal mean values.
c-----------------------------------------------------------------------
 
c-- compute Global & Zonal mean (veget+albedo+Veg_input_var) :
 
c-- Initialize :
      totsum = 0.
        vegzav(:,:)=0.
 
c-- Compute sum :
      do 250 i=1,nlat
        zavsum = 0.
        do 230 j=1,nlon
          if (fracgr(i,j).gt.epss) then
            zavsum = zavsum + fracgr(i,j)
             soiltype2(i,j)=0.
            vegzav(i,1) = vegzav(i,1) + fracgr(i,j)*st(i,j)
            vegzav(i,2) = vegzav(i,2) + fracgr(i,j)*sg(i,j)
            vegzav(i,3) = vegzav(i,3) + fracgr(i,j)*sd(i,j)
            vegzav(i,4) = vegzav(i,4) + fracgr(i,j)*st(i,j)*snlt(i,j)
            vegzav(i,5) = vegzav(i,5) + fracgr(i,j)*st(i,j)*blai(i,j,1)
            vegzav(i,6) = vegzav(i,6) + fracgr(i,j)*sg(i,j)*blai(i,j,2)
            vegzav(i,7) = vegzav(i,7) + fracgr(i,j)*
     &  (albland(i,j,1)+albland(i,j,2)+albland(i,j,3)+albland(i,j,4))
            vegzav(i,8) = vegzav(i,8) + fracgr(i,j)*temveg(i,j)
            vegzav(i,9) = vegzav(i,9) + fracgr(i,j)*gd0veg(i,j)
            vegzav(i,10)= vegzav(i,10)+ fracgr(i,j)*prcveg(i,j,1)
            vegzav(i,11)= vegzav(i,11)+ fracgr(i,j)*prcveg(i,j,2)
            vegzav(i,12)= vegzav(i,12)+ fracgr(i,j)*(anup(i,j)*darea(i))
            vegzav(i,13)= vegzav(i,13)+ fracgr(i,j)*(stock(i,j)*darea(i))
            vegzav(i,14)= vegzav(i,14)+ fracgr(i,j)*(b1(i,j)*darea(i))
            vegzav(i,15)= vegzav(i,15)+ fracgr(i,j)*(b2(i,j)*darea(i))
            vegzav(i,16)= vegzav(i,16)+ fracgr(i,j)*(b3(i,j)*darea(i))
            vegzav(i,17)= vegzav(i,17)+ fracgr(i,j)*(b4(i,j)*darea(i))
            vegzav(i,18)= vegzav(i,18)+ fracgr(i,j)*(pnpp(i,j)*darea(i))  
            vegzav(i,19)= vegzav(i,19)+ fracgr(i,j)*soiltype2(i,j)
          endif
 230    continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- Compute Zonal Mean Value :
        totsum = totsum + zavsum*darea(i)
        do 240 k=nvgmax,1,-1
          if ((k.ge.12).AND.(k.le.18)) then
            vegzav(i,k) = vegzav(i,k)*cfmzav(k)
            vegzav(0,k) = vegzav(0,k) + vegzav(i,k)
            goto 240
          endif
          if (k.ge.4 .and. k.le.6) then
            kk = (k-2)/2
            if (vegzav(i,kk).gt.epss) then
              vegzav(0,k) = vegzav(0,k) + vegzav(i,k)*darea(i)
              vegzav(i,k) = cfmzav(k)*vegzav(i,k)/vegzav(i,kk)
            else
              vegzav(i,k) = vegspv
            endif
          elseif (zavsum.gt.epss) then
            vegzav(0,k) = vegzav(0,k) + vegzav(i,k)*darea(i)
            vegzav(i,k) = cfmzav(k)*vegzav(i,k)/zavsum
          else
            vegzav(i,k) = vegspv
          endif
 240    continue
 250  continue

c     WRITE(*,*) 'NPP tot:',vegzav(0,18)
c     WRITE(*,'(E16.6)') (vegzav(i,18),i=1,nlat)
 
c-- Compute Global Mean Value :
      do 260 k=nvgmax,1,-1
        if ((k.ge.12).AND.(k.le.18)) goto 260 
        if (k.ge.4 .and. k.le.6) then
          kk = (k-2)/2
          if (vegzav(0,kk).gt.epss) then
            vegzav(0,k) = cfmzav(k)*vegzav(0,k)/vegzav(0,kk)
          else
            vegzav(0,k) = vegspv
          endif
        elseif (totsum.gt.epss) then
          vegzav(0,k) = cfmzav(k)*vegzav(0,k)/totsum
        else
          vegzav(0,k) = vegspv
        endif
 260  continue
 
c-- Write Global & Zonal_mean values :
      write(iveg+7,'(3A,I6)') 'Veget. Zonal_mean Output ; ',titveg(:nnctr)
     &                   , ' ; year=', iyear+iyr0vg
      write(iveg+7,'(99A)') (titzav(k),k=1,nvgmax)
      do k=1,nvgmax
        write(iveg+7,fmtzav) (vegzav(i,k),i=0,nlat)
      enddo
      cland=vegzav(0,13)
      write(iveg+7,*)
 
c-- compute & write Global + Zonal mean : End. ----------
c-----------------------------------------------------------------------
      if (iyear+nyvegt.gt.nyears) then
c       close(iveg+7)
      elseif (mod(iyear,nwrveg).eq.0) then
        call flush(iveg+7)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Sum 2.D vegetation map ; Prepare output .
c-----------------------------------------------------------------------
 
      vegsum = vegsum + 1.
      do 330 j=1,nlon
       do 330 i=1,nlat
         if (fracgr(i,j).gt.epss) then
           vegmap(i,j,1) = vegmap(i,j,1) + st(i,j)
           vegmap(i,j,2) = vegmap(i,j,2) + sg(i,j)
           vegmap(i,j,3) = vegmap(i,j,3) + sd(i,j)
           vegmap(i,j,4) = vegmap(i,j,4) + snlt(i,j)
           vegmap(i,j,5) = vegmap(i,j,5) + blai(i,j,1)
           vegmap(i,j,6) = vegmap(i,j,6) + blai(i,j,2)
           vegmap(i,j,7) = vegmap(i,j,7) +
     &  (albland(i,j,1)+albland(i,j,2)+albland(i,j,3)+albland(i,j,4))
           vegmap(i,j,8) = vegmap(i,j,8) + temveg(i,j)
           vegmap(i,j,9) = vegmap(i,j,9) + gd0veg(i,j)
           vegmap(i,j,10)= vegmap(i,j,10)+ prcveg(i,j,1)
           vegmap(i,j,11)= vegmap(i,j,11)+ prcveg(i,j,2)
           vegmap(i,j,12)= vegmap(i,j,12)+ (anup(i,j)*1E+12)
           vegmap(i,j,13)= vegmap(i,j,13)+ (stock(i,j)*1E+12) 
           vegmap(i,j,14)= vegmap(i,j,14)+ (b1(i,j)*1E+12) 
           vegmap(i,j,15)= vegmap(i,j,15)+ (b2(i,j)*1E+12)
           vegmap(i,j,16)= vegmap(i,j,16)+ (b3(i,j)*1E+12)
           vegmap(i,j,17)= vegmap(i,j,17)+ (b4(i,j)*1E+12)
           vegmap(i,j,18)= vegmap(i,j,18)+ (pnpp(i,j)*1E+12)
           vegmap(i,j,19)= vegmap(i,j,19)+ soiltype2(i,j)
         else
c          vegmap(i,j,19)= vegmap(i,j,19)+ soiltype2(i,j) 
           gd0veg(i,j) = -900.
         endif
 330  continue
 
      if (mod(iyear,nwrveg).eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c ==> decide to write 2.D map on file ; Prepare output :
 
c- x scaling Factor (if chg units) ; put SPecial_Value over non-land grid point
      vegsum = 1. / vegsum
      st_moy(:,:)=st_moy(:,:)+vegmap(:,:,1)*vegsum 
      vegsum2=vegsum2+1
      if (iyear+nwrveg.gt.nyears) then
       vegsum2= 1. / vegsum2
       st_moy(:,:)=st_moy(:,:)*vegsum2
      endif   
      do 370 k=1,nvgmax
       do 370 j=1,nlon
        do 370 i=1,nlat
          if (fracgr(i,j).gt.epss) then
            vegmap(i,j,k) = cfmzav(k)*vegmap(i,j,k)*vegsum
          else
            st_moy(i,j) = 0.
            vegmap(i,j,k) = vegspv
          endif
 370  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) write 2.D vegetation map on ascii file "veget.outp" :
c-----------------------------------------------------------------------
 
      CALL open(Yearly_Means)
      do k=1,numvegvar
         IF (output(newvegvar(k,4)))
     &   CALL write(k,vegmap(1:nlat,1:nlon,k))
      enddo

      do 450 k=1,nvgmax
        if (kveget.eq.0 .and. k.ne.1 .and. k.le.6) goto 450
        if (kveget.eq.0 .and. k.eq.19) goto 450 
c--Write 2 titles :
        if (k.eq.11) then
          write(iveg+8,'(3A,F6.3,A)') titmap(k), titveg(:nnctr),
     &                            ' (>', prcmin*1000., ' mm/day)'
        else
          write(iveg+8,'(2A)') titmap(k), titveg(:nnctr)
        endif
        write(iveg+8,*) 'year=', iyear+iyr0vg
c--Write vegetation distribution :
        do i=1,nlat
          write(iveg+8,fmtveg) (vegmap(i,j,k),j=1,nlon)
        enddo
        write(iveg+8,*)
 450  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- End of writing :
      if (iyear+nwrveg.gt.nyears) then
        close(iveg+8)
      else
        call flush(iveg+8)
      endif
      CALL close()
 
c- Reset to zero (for next Sum):
      vegsum = 0.
      do 470 k=1,nvgmax
       do 470 j=1,nlon
        do 470 i=1,nlat
          vegmap(i,j,k) = 0.
 470  continue
 
      if (kwradd.ne.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) write additional 2.D map on ascii file "veget.+[no_yr]" :
c-----------------------------------------------------------------------
      if (mod(iyear,kwradd).eq.0) then
 
c--open & write 2.l Header :
        write(ccny,'(I6)') iyear+iyr0vg
        do k=1,len(ccny)
          if (ccny(k:k).eq.' ') ccny(k:k) = '0'
        enddo
 
        open(85, file='veget.+'//ccny, status='unknown')
        write(85,1000) fmtveg, vegspv, nlon, nlat, 5, 65
        write(85,1111) 0., 5.625, -87.19, 5.625, 1., 1., 0
 
c- write Number of days/yr with Dry Soil Condition, i.e. Soil.W < bmtdry
        write(85,'(3A,F5.3,A)') 'Nb_Day/yr with Dry soil cond.; ',
     &                  titveg(:nnctr), ' (Soil.W <', bmtdry,'m)'
        write(85,*) 'year=', iyear+iyr0vg
        do 550 i=1,nlat
          do j=1,nlon
            var(j) = vegspv
            if (fracgr(i,j).gt.epss) var(j) = tpsdry(i,j)
          enddo
          write(85,fmtveg) (var(j),j=1,nlon)
 550    continue
        write(85,*)
 
c--
c- write 4 seasonal albedo map :
        do 570 ns=1,4
          write(85,'(2A)') 'Season. Albedo (o/o) (no snow) ; ',
     &                     titveg(:nnctr)
          write(85,*) 'year=', iyear+iyr0vg, '  saison=', ns
          do 560 i=1,nlat
            do j=1,nlon
              var(j) = vegspv
              if (fracgr(i,j).gt.epss) var(j) = 100.0*albland(i,j,ns)
            enddo
            write(85,fmtveg) (var(j),j=1,nlon)
 560      continue
          write(85,*)
 570    continue
 
        close(85)
 
      endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- write additional 2.D map : End. -----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- write 2.D map : End. ----------------
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine veget_wr -
      end
