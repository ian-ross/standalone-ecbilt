












      subroutine icdyna


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  This routine calls ice dynamic and advcect sea ice properties.
c---
c Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
c---
c  modif : 28/09/99
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
 
      integer     iyear,imonth,iday,iseason,ntotday,nbclins,nbtrops
      common /ec_timectl/ day,ntstep,nstpyear,nocstpyear,
     *                 iyear,imonth,iday,iseason,ntotday,nbclins,nbtrops 
      dimension trest(5)
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--1. Dynamics of sea ice-initialisation                               |
c-----------------------------------------------------------------------
c
      do 10 j=js1,js2
        do 10 i=is1(j),is2(j)
          tfu(i,j)   = abs(273.15-
     &                         0.0575*scal(i,j,ks2,2)+
     &                         1.710523e-3*sqrt(scal(i,j,ks2,2))**3-
     &                         2.154996e-4*scal(i,j,ks2,2)**2)
        
          phiss(i,j,0) = 0.0
          phiss(i,j,1) = 0.0
          phiss(i,j,2) = 0.0
          hgbqp(i,j)   = 0.0
 10    continue
       do k=ks1,ks2
        do j=js1,js2
          do i=is1(j),is2(j)
            phivs(i,j,k,1) = 0.0
          enddo
        enddo
       enddo
 
 
      if (idyn.eq.1) then
c
c     initialise array alcd for tracking rate of lead closure/opening 
c     due to convervenge/divergence.
c     initialise array alcr for extra amount of lead opening/closing 
c     due to shearing deformation. 
c
      do j=1,jmax
        do i=1,imax
          alcd(i,j) = zero
          alcr(i,j) = zero
        enddo
      enddo
c
c--1.1. Mean ice and snow thicknesses.
c--------------------------------------
c
        do 20 j=1,jmax
          do 20 i=1,imax
            hnm(i,j) = (1.0-albq(i,j))*hnbq(i,j)
            hgm(i,j) = (1.0-albq(i,j))*hgbq(i,j)
C           uo(i,j)  = uost(i,j)*tmu(i,j,ku2)
C           vo(i,j)  = vost(i,j)*tmu(i,j,ku2)
            uo(i,j)  = u(i,j,ku2)*tmu(i,j,ku2)
            vo(i,j)  = v(i,j,ku2)*tmu(i,j,ku2)
C           uo(i,j)  = umoy(i,j)*tmu(i,j,ku2)
C           vo(i,j)  = vmoy(i,j)*tmu(i,j,ku2)
20        continue
c
c--1.2. Call to dynamics routine.
c---------------------------------
c
c-no-ice velocity point
        do n=1,npo1i0
         trest(n)=tmu(ipo1i0(n),jpo1i0(n),ks2)
         tmu(ipo1i0(n),jpo1i0(n),ks2)= 0.0
        enddo
c
c  Northern hemisphere.
c
c Cdzh dynamics according to Zhang and Hibler, 1997.
c Cdzr dynamics according to Zhang and Rothrock, 2000.
c
        call dynami_zh(+1)
Cdzr    call dynami_zr(+1)
c
c  Southern hemisphere.
c
        call dynami_zh(-1)
Cdzr    call dynami_zr(-1)
c
c-no-ice velocity point
        do n=1,npo1i0
         tmu(ipo1i0(n),jpo1i0(n),ks2)= trest(n)
        enddo
c
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2. Computation of flux to the ocean.                                |
c-----------------------------------------------------------------------
c
        do 200 j=ju1,ju2
           ih=sign(1,j-jeq)
           sang = real(ih)*sangvg
           do 200 i=iu1(j),iu2(j)
Ccp1        if (icoupl .ne. 0) then
Ccp1         tairx = albq(i,j)*tairox(i,j)
Ccp1 &              +albq(i-1,j)*tairox(i-1,j)
Ccp1 &              +albq(i-1,j-1)*tairox(i-1,j-1)
Ccp1 &              +albq(i,j-1)*tairox(i,j-1)
C
Ccp1         tairy = albq(i,j)*tairoy(i,j)
Ccp1 &              +albq(i-1,j)*tairoy(i-1,j)
Ccp1 &              +albq(i-1,j-1)*tairoy(i-1,j-1)
Ccp1 &              +albq(i,j-1)*tairoy(i,j-1)
Ccp1        else
             tairx = (albq(i,j)+albq(i-1,j)
     &               +albq(i-1,j-1)+albq(i,j-1))*tairox(i,j)
C
             tairy = (albq(i,j)+albq(i-1,j)
     &               +albq(i-1,j-1)+albq(i,j-1))*tairoy(i,j)
Ccp1        endif
             zmod  = sqrt((ug(i,j)-uo(i,j))**2+(vg(i,j)-vo(i,j))**2)
             tglx  = (4-albq(i,j)-albq(i-1,j)-albq(i-1,j-1)-albq(i,j-1))
C    &               *rhoco*zmod*(ug(i,j)-uo(i,j))
     &               *rhoco*zmod*(cangvg*(ug(i,j)-uo(i,j))-
     &                            sang*(vg(i,j)-vo(i,j)))
             tgly  = (4-albq(i,j)-albq(i-1,j)-albq(i-1,j-1)-albq(i,j-1))
C    &               *rhoco*zmod*(vg(i,j)-vo(i,j))
     &               *rhoco*zmod*(cangvg*(vg(i,j)-vo(i,j))+
     &                            sang*(ug(i,j)-uo(i,j)))
             phisu(i,j) = -(tairx+1.0*tglx)/(4*rho0)
             phisv(i,j) = -(tairy+1.0*tgly)/(4*rho0)
200     continue
c
        do 210 j=js1,js2
           do 210 i=is1(j),is2(j)
              t11=rhoco*((ug(i-1,j-1)-uo(i-1,j-1))**2
     &                  +(vg(i-1,j-1)-vo(i-1,j-1))**2)
              t12=rhoco*((ug(i-1,j)-uo(i-1,j))**2
     &                  +(vg(i-1,j)-vo(i-1,j))**2)
              t21=rhoco*((ug(i,j-1)-uo(i,j-1))**2
     &                  +(vg(i,j-1)-vo(i,j-1))**2)
              t22=rhoco*((ug(i,j)-uo(i,j))**2
     &                  +(vg(i,j)-vo(i,j))**2)
              tmoy=0.25*(t11+t12+t21+t22)
              tot=tmoy*(1-albq(i,j))
     &            +albq(i,j)*sqrt(tairox(i,j)**2+tairoy(i,j)**2)
              zustm=tot/rho0
              ust2s(i,j)=zustm*(one+sdvt(i,j))
210      continue
      else
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--3. If no ice dynamics                                               |
c-----------------------------------------------------------------------
c
        do 1140 j=ju1,ju2
          do 1130 i=iu1(j),iu2(j)
Ccp1        if (icoupl .ne. 0) then
Ccp1           phisu(i,j) = -(tairox(i,j) + tairox(i-1,j)
Ccp1 &                      + tairox(i-1,j-1) + tairox(i,j-1))/(rho0*4)
Ccp1           phisv(i,j) = -(tairoy(i,j) + tairoy(i-1,j)
Ccp1 &                      + tairoy(i-1,j-1) + tairoy(i,j-1))/(rho0*4)
Ccp1        else
              phisu(i,j) = -(tairox(i,j))/(rho0)
              phisv(i,j) = -(tairoy(i,j))/(rho0)
Ccp1        endif
 1130     continue
 1140   continue
c
        do 1120 j=js1,js2
           do 1110 i=is1(j),is2(j)
              tot=sqrt(tairox(i,j)**2+tairoy(i,j)**2)
              zustm=tot/rho0
              ust2s(i,j)=zustm*(one+sdvt(i,j))
 1110      continue
 1120   continue
c
      endif
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--4. Set umoy, vmoy at zero                                           |
c-----------------------------------------------------------------------
 
      do j=ju1,ju2
        do i=iu1(j),iu2(j)
C         umoy(i,j)=0.0
C         vmoy(i,j)=0.0
         enddo
      enddo
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine icdyna -
      end
