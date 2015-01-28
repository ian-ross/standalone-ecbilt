












      subroutine outave(ja,xjour,xjour1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  This routine computes the average of some variables and write it
c  on the ouput files.
c  ATTENTION cette routine devrait etre valable meme si le pas de temps est
c  different de 1 jour. Pour plus de generalite il faudrait introduire la
c  variable xjours indiquant le jour au pas de temps suivant pour savoir si
c  on restera dans le meme mois ou la meme annee. De plus il faudra calculer
c  njucum par un cumul lors de chaque passage.
c  modif : 24/12/99
 
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


c- nn99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
      common / mchd99 / nn99
      common / zfluc / zfluxmt,zfluxmts


c--Output variables and arrays.
c
      integer njcum(0:12)
      integer iwl(30),jwl(30)
      parameter (noumax=47,noumax2=30)


c Choice noumax :18,34,43,47, noumax2 depending of the choice
      real*4 cmoyan(imax,jmax,noumax),cmoymo(imax,jmax,noumax),
     &       cmoymap(imax,jmax,noumax2,12)
      real cmulti(noumax),cadd(noumax)


      real*4 sort1(kmax),sort2(kmax),sort3(kmax)
      integer ntypou(noumax),nmm(noumax),nmal(noumax),nma(noumax),
     &        ncor(noumax2),nmap(12)
      character*60 titn(noumax)
      real*4 umm,vmm,tmm,smm,uma,vma,tma,sma,wma,wmm,rmm,rma
      real*4 rtemp,T,S


c streamfunction part
c options: nstreamupdate=0(no output), 1(monthly) or 2(annually)
      parameter (nstreamout=1)



      common/totalanu/ umm(imax,jmax,kmax),vmm(imax,jmax,kmax),
     &           tmm(imax,jmax,kmax),smm(imax,jmax,kmax),
     &           uma(imax,jmax,kmax),vma(imax,jmax,kmax),
     &           tma(imax,jmax,kmax),sma(imax,jmax,kmax),
     &           wma(imax,jmax,kmax),wmm(imax,jmax,kmax),
     &           rma(imax,jmax,kmax),rmm(imax,jmax,kmax)
 


c 
      common/moyout/ cmoyan,cmoymo,cmoymap
      common/parout/ titn,npack,klevm(4,2),klevv,njm
      common/parout2/ noumef,nmm,nmal,nma,nmmef,nmaef,nmalef
      common/parout3/ ncor,nmap,ntypou,jmois,nptj,iwl,jwl
      common/parout4/ cmulti,cadd,jmoisi,jmoifl
      common/parout5/ znumas(imax,jmax,4),znumav(imax,jmax,5)
 
c
Ccp0  data njcum/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data njcum/0,30,60,90,120,150,180,210,240,270,300,330,360/
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
c--cresujl.dat : CONTAINS ICE LOCAL DAILY OUTPUTS.
c--cresujg.dat : CONTAINS ICE GLOBAL DAILY OUTPUTS.
c--cresum.dat  : CONTAINS MONTHLY OUTPUTS.
c--cresua.dat  : CONTAINS YEARLY OUTPUTS.
c
      open(31,file='cresujl.dat',form='unformatted')
      nmtot=365
      write(31) nmtot
      open(32,file='cresujg.dat',form='unformatted')
      open(33,file='cresum.dat',form='unformatted')
      open(34,file='cresua.dat',form='unformatted')
      open(35,file='cresal.dat',form='unformatted')
      open(36,file='crestm.dat',form='unformatted')
      open(37,file='cresta.dat',form='unformatted')
      open(42,file='correcw.dat',form='formatted')


c
c--1.2. SELECTION OF POINTS FOR LOCAL DAILY OUTPUTS.
c---------------------------------------------------
c
c--cpointj.dat: CONTAINS POINTS FOR LOCAL DAILY OUTPUTS.
c
      nptj=0
      open(30,file='cpointj.dat',status='old',err=6)
      read(30,*) nptj
      do 5 npw=1,nptj
        read(30,*) iwl(npw),jwl(npw)
 5    continue
      close(30)
 6    continue
c
c--1.3 Computation of the month.
c---------------------------------
c
      jmois=1
      do 10 ii=1,ijour
         if (ijour.gt.njcum(jmois)) jmois=jmois+1
 10   continue
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


      do 140 n=1,noumax
         do 130 j=1,jmax
           do 120 i=1,imax
            cmoyan(i,j,n)=0.0
            cmoymo(i,j,n)=0.0
120         continue
130      continue
140   continue
      do  n=1,noumax2
         do  j=1,jmax
           do  i=1,imax
            do jj=1,12
              cmoymap(i,j,n,jj)=0.0
            enddo
           enddo
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
              umm(i,j,k)= 0.0
              vmm(i,j,k)= 0.0
              wmm(i,j,k)= 0.0
              tmm(i,j,k)= 0.0
              smm(i,j,k)= 0.0
              rmm(i,j,k)= 0.0
              uma(i,j,k)= 0.0
              vma(i,j,k)= 0.0
              wma(i,j,k)= 0.0
              tma(i,j,k)= 0.0
              sma(i,j,k)= 0.0
              rma(i,j,k)= 0.0
            enddo
         enddo
      enddo


c
c--1.5. Reading output.param
c---------------------------
       open(30,file='output.param')
c
       read(30,*)
       read(30,*)
       read(30,*)
       read(30,*) 
       read(30,*)
C      read(30,*) noumef
       read(30,*) npack
c
c Correspondance between npack and he number of fields
       noumef = 0
       if (npack.eq.1) noumef=18
       if (npack.eq.2) noumef=34
       if (npack.eq.3) noumef=43
       if (npack.eq.4) noumef=47
c
       if (noumef.gt.noumax) then
         write(iuo+66,*) 'STOP : Too much fields for monthly mean'
         write(iuo+66,*) 'noumef=',noumef,'  noumax=',noumax
         STOP
       endif
       if (noumef.lt.noumax) then
         write(iuo+66,*) 'WARNING : You can reduce the memory requirements'
         write(iuo+66,*) '  by reducing noumax in outave : noumax=',noumax
         write(iuo+66,*) '  and it can de reduced to noumef=', noumef
       endif
c
c Computation of the levels for averages
       do kz=ks2,ks1,-1
         if ((-60).ge.zw(kz)) then
            klevv=kz
            goto 141
         endif
       enddo
141    continue
       do kk= 1,4
         klevm(kk,1)=0
         klevm(kk,2)=0
       enddo
       if (nn99.eq.2) then
         write(99,*)
         write(99,'(A,I3)')
     &     ' Vertical averages for outputs(outave) : klevv=', klevv
       endif
       read(30,*)
       do kk=1,4
         read(30,*) zlim1,zlim2
          do kz=ks2,ks1,-1
            if ((-zlim1).ge.zw(kz)) then
              klevm(kk,1)=kz
              goto 142
            endif
          enddo
142       continue
          do kz=ks2,ks1,-1
            if ((-zlim2).ge.zw(kz)) then
              klevm(kk,2)=kz+1
              goto 143
            endif
          enddo
143       if (klevm(kk,2).eq.0) klevm(kk,2)=1
 
          if (nn99.eq.2) then
            write(99,*) 'Layer number',kk
            write(99,*) 'lim1 and lim2 :',zlim1,zlim2
            write(99,*) 'lev1 and lev2 :',klevm(kk,1),klevm(kk,2)
          endif
       enddo
c  Number of value for each average
      do j=2,jmax-1
       do i=2,imax-1
         znumav(i,j,5)=0.0
         do kk=klevv,ku2
           znumav(i,j,5)=znumav(i,j,5)+tmu(i,j,kk)
         enddo
         znumav(i,j,5)=max(znumav(i,j,5),undim6)
       enddo
      enddo
 
      do j=ju1,ju2
       do i=iu1(j),iu2(j)
         do kz=1,4
          znumav(i,j,kz)=0.0
          do kk=klevm(kz,2),klevm(kz,1)
            znumav(i,j,kz)=znumav(i,j,kz)+tmu(i,j,kk)
          enddo
          znumav(i,j,kz)=max(znumav(i,j,kz),undim6)
         enddo
       enddo
      enddo
 
      do j=js1,js2
       do i=is1(j),is2(j)
         do kz=1,4
          znumas(i,j,kz)=0.0
          do kk=klevm(kz,2),klevm(kz,1)
            znumas(i,j,kz)=znumas(i,j,kz)+tms(i,j,kk)
          enddo
          znumas(i,j,kz)=max(znumas(i,j,kz),undim6)
         enddo
       enddo
      enddo
 
       read(30,*)
       read(30,*)
       read(30,*)
       do 150 npw=1,noumef
         read(30,'(A)') titn(npw)
         read(30,*) i,ntypou(npw),nmm(npw),nmal(npw),nma(npw),
     &               cmulti(npw),cadd(npw)
150    continue
       close(30)
       nmmef=0
       nmalef=0
       nmaef=0
       do 160 npw=1,noumef
          if (nmm(npw).eq.1) then
            nmmef=nmmef+1
          endif
          if (nma(npw).eq.1) then
            nmaef=nmaef+1
          endif
          if (nmal(npw).eq.1) then
            nmalef=nmalef+1
            ncor(nmalef)=npw
          endif
160    continue
C      write(iuo+66,*) 'nmmef etc'
C      write(iuo+66,*) nmmef,nmaef,nmalef,ncor
       if (nmalef.gt.noumax2) then
         write(iuo+66,*) 'STOP : Too much fields for all period mean'
         STOP
       endif
       if (nmalef.lt.noumax2) then
         write(iuo+66,*) 'WARNING : You can reduce the memory requirements'
         write(iuo+66,*) '  by reducing noumax2 in outave : noumax2=',noumax2
         write(iuo+66,*) '  and it can de reduced to nmalef=', nmalef
       endif
 
      endif
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2. Computation of averages                                          |
c-----------------------------------------------------------------------
c
c--2.1. WRITE LOCAL DAILY OUTPUTS.
c---------------------------------
c
      if (mod(ja,nwjl).eq.nwtest) then
c
            nom2=4
            nom3=3
            write(31) nom2,nom3
            write(31) ja,ijour
            do 220 np=1,nptj
              write(31) 'Epaisseur de neige   (71)'
              write(31) hnbq(iwl(np),jwl(np))
              write(31) 'Epaisseur de glace   (72)'
              write(31) hgbq(iwl(np),jwl(np))
              write(31) 'Epaisseur cree       (73)'
              write(31) hgbqp(iwl(np),jwl(np))
              write(31) 'Aire des leads       (74)'
              write(31) albq(iwl(np),jwl(np))
C             write(31) 'Temp de surface      (75)'
C             write(31) ts(iwl(np),jwl(np))
C             write(31) 'Temp de neige        (76)'
C             write(31) tbq(iwl(np),jwl(np),1)
C             write(31) 'Temp de glace 1      (77)'
C             write(31) tbq(iwl(np),jwl(np),2)
C             write(31) 'Temp de glace 2      (78)'
C             write(31) tbq(iwl(np),jwl(np),3)
C             write(31) 'var volume surf      (79)'
C             write(31) dvosbq(iwl(np),jwl(np))
C             write(31) 'var volume bott      (80)'
C             write(31) dvobbq(iwl(np),jwl(np))
C             write(31) 'var volume lead      (81)'
C             write(31) dvolbq(iwl(np),jwl(np))
C             write(31) 'var volume neige     (82)'
C             write(31) dvonbq(iwl(np),jwl(np))
C             write(31) 'flux base glace      (83)'
C             write(31) fbbq(iwl(np),jwl(np))
C             write(31) 'temp air             (84)'
C             write(31) tabq(iwl(np),jwl(np))
C             write(31) 'humid air            (85)'
C             write(31) qabq(iwl(np),jwl(np))
C             write(31) 'vit vent             (86)'
C             write(31) vabq(iwl(np),jwl(np))
C             write(31) 'flux sol             (87)'
C             write(31) fsolg(iwl(np),jwl(np))
C             write(31) 'flux IR              (88)'
C             write(31) firg(iwl(np),jwl(np))
C             write(31) 'flux sensible        (89)'
C             write(31) fcsg(iwl(np),jwl(np))
C             write(31) 'flux latent          (90)'
C             write(31) fleg(iwl(np),jwl(np))
 
              do kk=1,kmax
                 zms=tms(iwl(np),jwl(np),kk)
                 sort1(kk)=zms*(scal(iwl(np),jwl(np),kk,1)
     &                    -273.15) + (1.0-zms)*(-100.0)
                 sort2(kk)=zms*scal(iwl(np),jwl(np),kk,2)
     &                     + (1.0-zms)*(-100.0)
                 sort3(kk)=zms*q2turb(iwl(np),jwl(np),kk)
     &                     + (1.0-zms)*(-100.0)
Ctk0             sort3(kk)=zms*bvf(iwl(np),jwl(np),kk)
Ctk0 &                     + (1.0-zms)*(-100.0)
              enddo
              write(31) 'temperature ocean    (91)'
C             write(31) (scal(iwl(np),jwl(np),kv,1),kv=1,kmax)
              write(31) (sort1(kv),kv=1,kmax)
              write(31) 'Salinite ocean       (92)'
C             write(31) (scal(iwl(np),jwl(np),kv,2),kv=1,kmax)
              write(31) (sort2(kv),kv=1,kmax)
              write(31) 'Q2turb               (93)'
C             write(31) 'n2                   (93)'
C             write(31) (q2turb(iwl(np),jwl(np),kv),kv=2,kmax+1)
              write(31) (sort3(kv),kv=1,kmax)
 
220         continue
c
      endif
 
c
c computation of annula mean surface velocity
      zfluxmt=zfluxmt+zflux0
      zfluxmts=zfluxmts+zflux0s
c
      if (nwtal.eq.1.or.mod(ja,nwm).eq.nwtest) then
        njm=njm+1
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.2. CUMULATION OF GLOBAL MONTHLY OUTPUTS.
c--------------------------------------------
c streamfunction part
c cumulation for monthly and/or yearly outputs !
       if (nstreamout.ne.0) call streamfunc(0,0.d0,nstreamout)


       do k=1,kmax
          do j=2,jmax-1
             do i=2,imax-1
               T=scal(i,j,k,1)-273.15
               S=scal(i,j,k,2)
               rtemp=999.842594+(T*6.793952e-02)- 
     &         (9.095290e-3*T**2)+(1.001685e-4*T**3)- 
     &         (1.120083e-6*T**4)+(6.536332e-9*T**5)+ 
     &         S*(0.824493-(4.0899e-3*T)+(7.6438e-5* 
     &         T**2)-(8.2467e-7*T**3)+(5.3875e-9* 
     &         T**4))+(S*1.5)*(0-5.72466e-3+(1.0227e-4* 
     &         T)-(1.6546e-6*T**2))+4.8314e-4* 
     &         S**2

               umm(i,j,k) = umm(i,j,k)+u(i,j,k)
               vmm(i,j,k) = vmm(i,j,k)+v(i,j,k)
               wmm(i,j,k) = wmm(i,j,k)+w(i,j,k)
               tmm(i,j,k) = tmm(i,j,k)+scal(i,j,k,1)
               smm(i,j,k) = smm(i,j,k)+scal(i,j,k,2)
               rmm(i,j,k) = rmm(i,j,k)+rtemp
               uma(i,j,k) = uma(i,j,k)+u(i,j,k)
               vma(i,j,k) = vma(i,j,k)+v(i,j,k)
               wma(i,j,k) = wma(i,j,k)+w(i,j,k)
               tma(i,j,k) = tma(i,j,k)+scal(i,j,k,1)
               sma(i,j,k) = sma(i,j,k)+scal(i,j,k,2)
               rma(i,j,k) = rma(i,j,k)+rtemp
             enddo
          enddo
       enddo


        if (npack.ge.1) then


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
     &             )
              cmoymo(i,j,17) = cmoymo(i,j,17) + ai(ks2) * sqrt(
     &          c4x(i,j,ks2)*c4x(i,j,ks2)+c4y(i,j,ks2)*c4y(i,j,ks2) )
                cmoymo(i,j,18)=cmoymo(i,j,18)+eta(i,j)
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
C               write(155,*) sigsurf*tms(i,j,ks2)
 
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
 245          continue
              cmoymo(i,j,14)= cmoymo(i,j,14)-zw(kmin3-1)
              cmoymo(i,j,13)= cmoymo(i,j,13)-zw(kmin4)
            enddo
          enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        endif
c-------
        if (npack.ge.2) then
            do kz=1,4
             kji=21+(kz-1)*4
             do kk=klevm(kz,2),klevm(kz,1)
              do j=ju1,ju2
               do i=iu1(j),iu2(j)
                  cmoymo(i,j,kji)=cmoymo(i,j,kji)+
     &                tmu(i,j,kk)*u(i,j,kk)/znumav(i,j,kz)
                  cmoymo(i,j,kji+1)=cmoymo(i,j,kji+1)+
     &                tmu(i,j,kk)*v(i,j,kk)/znumav(i,j,kz)
               enddo
              enddo
             enddo
            enddo
            do kz=1,4
            kji=19+(kz-1)*4
             do kk=klevm(kz,2),klevm(kz,1)
              do j=js1,js2
               do i=is1(j),is2(j)
                   cmoymo(i,j,kji)=cmoymo(i,j,kji)+
     &               tms(i,j,kk)*scal(i,j,kk,1)/znumas(i,j,kz)
                   cmoymo(i,j,kji+1)=cmoymo(i,j,kji+1)+
     &               tms(i,j,kk)*scal(i,j,kk,2)/znumas(i,j,kz)
               enddo
              enddo
             enddo
            enddo
        endif
c-------
        if (npack.ge.3) then
            do j=js1,js2
              do  i=is1(j),is2(j)
                cmoymo(i,j,35)= cmoymo(i,j,35)+fcm1(i,j)
                cmoymo(i,j,36)= cmoymo(i,j,36)+tenagx(i,j)
                cmoymo(i,j,37)= cmoymo(i,j,37)+tenagy(i,j)
                cmoymo(i,j,38)= cmoymo(i,j,38)+fsolg(i,j)*
     &                     (1.0-albq(i,j))+albq(i,j)*fsolcn(i,j)
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
        if (npack.ge.4) then
          do j=js1,js2
            do  i=is1(j),is2(j)
                cmoymo(i,j,44)= cmoymo(i,j,44)+tabq(i,j)
                cmoymo(i,j,45)= cmoymo(i,j,45)+qabq(i,j)
                cmoymo(i,j,46)= cmoymo(i,j,46)+vabq(i,j)
                cmoymo(i,j,47)= cmoymo(i,j,47)+hnplbq(i,j)
            enddo
          enddo
        endif
c-------
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- COMPUTES MONTHLY MEAN
c
C         if (mod(ja,nwm).eq.nwtest.and.ijour.eq.njcum(jmois)) then
C         if (ijour.eq.njcum(jmois)) then



      if (ijour.eq.njcum(jmois).and.ijour.ne.ijour1) then
c
C           njm = (njcum(jmois)-njcum(jmois-1))-ijouri*(1-jjoufl)
            jjoufl = 1
            usnjm = 1.0/dble(njm)
            njm=0
c
            do 300 n=1,noumef
             do 300 j=1,jmax
              do 300 i=1,imax
C              cmoymo(i,j,n)=cmoymo(i,j,n)*usnjm
               cmoymo(i,j,n)=(cmoymo(i,j,n)*usnjm)*cmulti(n)+cadd(n)
300         continue
c
            do jj=1,jmax
               do ij=1,imax
                  zindp=max(0.0,sign(1.0,cmoymo(ij,jj,2)-0.2))
                  cmoymo(ij,jj,6)=cmoymo(ij,jj,6)*zindp
               enddo
            enddo
c
            do 320 n=1,noumef
             do 320 j=1,jmax
              do 320 i=1,imax
                zindo       = tms(i,j,ks2)
                cmoymo(i,j,n)=zindo*cmoymo(i,j,n)+(1.0-zindo)*spvr
320         continue
        if (npack.ge.2) then
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
325         continue
        endif
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c WRITES MONTHLY MEAN
 
        if (mod(ja,nwm).eq.nwtest) then
              if (nn99.eq.2) write(99,'(2A,3I6)') 'Write monthly mean',
     &        ' on file cresum.out ; year,month,day=', ja, jmois, ijour
              write(33) ja,jmois
              write(33) nmmef
              do n=1,noumef
                if (nmm(n).eq.1) then
                  write(33) n,ntypou(n),titn(n)
                  do j=1,jmax
                      write(33) (cmoymo(i,j,n),i=1,imax)
                  enddo
                endif
              enddo


c streamfunction part
          if (nstreamout.eq.1) call streamfunc(1,usnjm,nstreamout)


          if (nwtom.eq.1) then
            do j=2,jmax-1
             do i=2,imax-1
              do k=1,kmax
                umm(i,j,k)= tmu(i,j,k)*umm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spvr
                vmm(i,j,k)= tmu(i,j,k)*vmm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spvr
                wmm(i,j,k)= tmu(i,j,k)*wmm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spvr
                tmm(i,j,k)= tms(i,j,k)*tmm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spvr
                smm(i,j,k)= tms(i,j,k)*smm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spvr
                rmm(i,j,k)= tms(i,j,k)*rmm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spvr
              enddo
             enddo
            enddo
            ntotoc=999
            write(36) ja,jmois
            write(36) ntotoc
            write(36) umm
            write(36) vmm
            write(36) wmm  
            write(36) tmm
            write(36) smm
	    write(36) rmm  
            do j=2,jmax-1
             do i=2,imax-1
              do k=1,kmax
                umm(i,j,k)= 0.0
                vmm(i,j,k)= 0.0
                wmm(i,j,k)= 0.0  
                tmm(i,j,k)= 0.0
                smm(i,j,k)= 0.0
                rmm(i,j,k)= 0.0
              enddo
             enddo
            enddo
          endif


        endif
 
c
c                       UPDATE jmois.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.4. CUMULATION FOR GLOBAL YEARLY OUTPUTS AND ALL PERIOD MEAN
c---------------------------------------------------------------
c
              do 340 n=1,nmalef
               do 340 j=1,jmax
                do 340 i=1,imax
                  cmoymap(i,j,n,jmois)=cmoymap(i,j,n,jmois)+
     &                                 cmoymo(i,j,ncor(n))
340           continue
              nmap(jmois)=nmap(jmois)+1
 
            if (mod(ja,nwa).eq.nwtest) then
c
              do 350 n=1,noumef
               do 350 j=1,jmax
                do 350 i=1,imax
                  cmoyan(i,j,n)=cmoyan(i,j,n)+cmoymo(i,j,n)
350           continue
c
            endif
c
c                       UPDATE jmois.
            jmois = jmois+1
            if (jmois.gt.12) jmois = 1



            do 343 n=1,noumef
             do 343 j=1,jmax
              do 343 i=1,imax
                cmoymo(i,j,n) =0.0
343           continue


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
370       continue
c
          if (nn99.eq.2) write(99,'(2A,3I6)') 'Write annual mean ',
     &     ' on file cresua.out ; year,month,day=', ja, jmois, ijour
          write(iuo+66,*) 'moyenne annuelle',ja
          write(34) ja,jmois
          write(34) nmaef
          do n=1,noumef
            if (nma(n).eq.1) then
              write(34) n,ntypou(n),titn(n)
              do j=1,jmax
                  write(34) (cmoyan(i,j,n),i=1,imax)
              enddo
            endif
          enddo
c
          do 390 n=1,noumef
           do 390 j=1,jmax
             do 390 i=1,imax
                cmoyan(i,j,n)=0.0
390       continue


          zfluxm=zfluxmt/yeaday
          zfluxms=zfluxmts/yeaday
          write(42,*) zfluxm,zfluxms
          zfluxmt=0.0
          zfluxmts=0.0



c streamfunction part
          if (nstreamout.eq.2) call streamfunc(1,usanms*usnjm,nstreamout)


          if (nwtoa.eq.1) then
            usnja=1/yeaday
            do j=2,jmax-1
             do i=2,imax-1
              do k=1,kmax
                uma(i,j,k)= tmu(i,j,k)*uma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spvr
                vma(i,j,k)= tmu(i,j,k)*vma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spvr
                wma(i,j,k)= tmu(i,j,k)*wma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spvr
                tma(i,j,k)= tms(i,j,k)*tma(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spvr
                sma(i,j,k)= tms(i,j,k)*sma(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spvr
                rma(i,j,k)= tms(i,j,k)*rma(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spvr
              enddo
             enddo
            enddo
            ntotoc=999
            write(37) ja,jmois
            write(37) ntotoc
            write(37) uma
            write(37) vma
            write(37) wma
            write(37) tma
            write(37) sma
            write(37) rma
          endif
            do j=2,jmax-1
             do i=2,imax-1
              do k=1,kmax
                uma(i,j,k)= 0.0
                vma(i,j,k)= 0.0
                wma(i,j,k)= 0.0 
                tma(i,j,k)= 0.0
                sma(i,j,k)= 0.0
                rma(i,j,k)= 0.0  
              enddo
             enddo
            enddo



c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- WRITE GLOBAL ALL PERIOD OUTPUTS.
      if (numit.eq.nlast.and.nwtal.eq.1) then
 
          do 410 k=1,12
           do 400 n=1,nmalef
            do 400 j=1,jmax
             do 400 i=1,imax
               zindo        = tms(i,j,ks2)
               cmoymap(i,j,n,k)=zindo*cmoymap(i,j,n,k)
     &                          /max(dfloat(nmap(k)),undim6)
     &                          +(1.0-zindo)*spvr
400        continue
           if (nn99.eq.2) write(99,'(2A,3I6)') 'Write monthly mean',
     &       ' on all period on cresal.out ; yy,mm,dd=',ja,jmois,ijour
           write(35) ja,k
           write(35) nmalef
           do n=1,nmalef
             write(35) ncor(n),ntypou(ncor(n)),titn(ncor(n))
             do j=1,jmax
                  write(35) (cmoymap(i,j,n,k),i=1,imax)
             enddo
           enddo
 410      continue
      endif
c
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine outave -
      end
