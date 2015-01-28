












c
      subroutine advx(dt,ut,crh,sm,s0,sx,sxx,sy,syy,sxy)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
c  PRATHER'S (1986) Advection scheme (periodic conditions).
c
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'dynami.com'
c
      dimension ut(imax,jmax)
      dimension sm(imax,jmax)
      dimension s0(imax,jmax),sx(imax,jmax),sxx(imax,jmax),
     &          sy(imax,jmax),syy(imax,jmax),sxy(imax,jmax)
      dimension f0(imax,jmax),fx(imax,jmax),fxx(imax,jmax),
     &          fy(imax,jmax),fyy(imax,jmax),fxy(imax,jmax),
     &          fm(imax,jmax)
c
c  bet: Coefficient with values 0/1 to avoid over-writing
c	of arrays.
c
      dimension alg(imax,jmax),alg1(imax,jmax),alg1q(imax,jmax),
     &		bet(imax,jmax)
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1) Limitation of moments.                                           |
c-----------------------------------------------------------------------
c
      do 5 j=1,jmax
         do 5 i=1,imax
         bet(i,j) = 0.
 5    continue
      do 15 j=1,jmax
        do 10 i=1,imax
          slpmax   = max(zero,s0(i,j))
          s1max    = 1.5*slpmax
          s1new    = min(s1max,max(-s1max,sx(i,j)))
          s2new    = min((2.0*slpmax-0.3334*abs(s1new)),
     &               max(abs(s1new)-slpmax,sxx(i,j)))
	  s0(i,j)  = slpmax  
          sx(i,j)  = s1new
          sxx(i,j) = s2new
          sxy(i,j) = min(slpmax,max(-slpmax,sxy(i,j)))
10      continue
15    continue
c
c--1.1. CASE OF EMPTY BOXES. 
c---------------------------
c
      do 29 j=1,jmax
	do 27 i=1,imax
          slpmax   = max(zero,s0(i,j))
          zin0     = 1.0-max(zero,sign(one,-slpmax))
          sx(i,j)  = sx(i,j)*zin0
          sxx(i,j) = sxx(i,j)*zin0
          sy(i,j)  = sy(i,j)*zin0
          syy(i,j) = syy(i,j)*zin0
          sxy(i,j) = sxy(i,j)*zin0
27      continue
29    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2) Apply mask.                                                      |
c-----------------------------------------------------------------------
c
      do 40 j=1,jmax
        do 30 i=1,imax
  	  sx(i,j)  = sx(i,j)*zindfa(i,j)
	  sxx(i,j) = sxx(i,j)*zindfa(i,j)
	  sy(i,j)  = sy(i,j)*zindfa(i,j)
	  syy(i,j) = syy(i,j)*zindfa(i,j)
	  sxy(i,j) = sxy(i,j)*zindfa(i,j)
30      continue
40    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3) Initialize volumes of boxes.                                     |
c-----------------------------------------------------------------------
c
      do j=1,jmax
	do i=1,imax
  	  sm(i,j) = crh*area(i,j)+(1.0-crh)*sm(i,j)
	enddo
      enddo
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  4) Calculate fluxes and moments between boxes i<-->i+1              |
c-----------------------------------------------------------------------
c
      do 80 j=js1,js2
        do 70 i=is1(j)-1,is2(j)
c
c--4.1. FLUX FROM i TO i+1 WHEN u.GT.0.0. 
c----------------------------------------
c     
          bet(i,j) = max(zero,sign(one,ut(i,j)))
          alf      = (max(zero,ut(i,j))*dt*dxs2(i,j))/sm(i,j)
          alfq     = alf*alf
          alf1     = 1.0-alf
          alf1q    = alf1*alf1
          fm(i,j)  = alf*sm(i,j)
          f0(i,j)  = alf*(s0(i,j)+alf1*(sx(i,j)+(alf1-alf)*sxx(i,j)))
          fx(i,j)  = alfq*(sx(i,j)+3.0*alf1*sxx(i,j))
          fxx(i,j) = alf*alfq*sxx(i,j)
          fy(i,j)  = alf*(sy(i,j)+alf1*sxy(i,j))
          fxy(i,j) = alfq*sxy(i,j)
          fyy(i,j) = alf*syy(i,j)
c
c  4.1.1. Readjust moments remaining in the box.
c
          sm(i,j)  = sm(i,j)-fm(i,j)
          s0(i,j)  = s0(i,j)-f0(i,j)
          sx(i,j)  = alf1q*(sx(i,j)-3.0*alf*sxx(i,j))
          sxx(i,j) = alf1*alf1q*sxx(i,j)
          sy(i,j)  = sy(i,j)-fy(i,j)
          syy(i,j) = syy(i,j)-fyy(i,j)
          sxy(i,j) = alf1q*sxy(i,j)
70	continue
80    continue
c
      do 100 j=js1,js2
        do 90 i=is1(j)-1,is2(j)
c
c--4.2. FLUX FROM i+1 TO i WHEN u.LT.0.0.
c----------------------------------------
c
	  ii         = i+1
          alf        = (max(zero,-ut(i,j))*dt*dxs2(i,j))/sm(ii,j)
          alg(i,j)   = alf
          alfq       = alf*alf
          alf1       = 1.0-alf
	  alg1(i,j)  = alf1
          alf1q      = alf1*alf1
	  alg1q(i,j) = alf1q
	  fm(i,j)    = fm(i,j)+alf*sm(ii,j)
          f0(i,j)    = f0(i,j)+
     &	               alf*(s0(ii,j)-alf1*(sx(ii,j)
     &                 -(alf1-alf)*sxx(ii,j)))
          fx(i,j)    = fx(i,j)+alfq*(sx(ii,j)-3.0*alf1*sxx(ii,j))
          fxx(i,j)   = fxx(i,j)+alf*alfq*sxx(ii,j)
          fy(i,j)    = fy(i,j)+alf*(sy(ii,j)-alf1*sxy(ii,j))
          fxy(i,j)   = fxy(i,j)+alfq*sxy(ii,j)
          fyy(i,j)   = fyy(i,j)+alf*syy(ii,j)
90      continue
100   continue 
      do 120 j=js1,js2
        do 110 i=is1(j),is2(j)
c
c  4.2.1. Readjust moments remaining in the box. 
c
	  ii       = i-1
          sm(i,j)  = bet(ii,j)*sm(i,j)
     &               +(1.0-bet(ii,j))*(sm(i,j)-fm(ii,j))
          s0(i,j)  = bet(ii,j)*s0(i,j)
     &               +(1.0-bet(ii,j))*(s0(i,j)-f0(ii,j))
          sx(i,j)  = alg1q(ii,j)*(sx(i,j)+3.0*alg(ii,j)*sxx(i,j))
          sxx(i,j) = alg1(ii,j)*alg1q(ii,j)*sxx(i,j)
          sy(i,j)  = bet(ii,j)*sy(i,j)
     &               +(1.0-bet(ii,j))*(sy(i,j)-fy(ii,j))
          syy(i,j) = bet(ii,j)*syy(i,j)
     &               +(1.0-bet(ii,j))*(syy(i,j)-fyy(ii,j))
          sxy(i,j) = alg1q(ii,j)*sxy(i,j)
110     continue
120   continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  5) Put the temporary moments into appropriate neighboring boxes.    |
c-----------------------------------------------------------------------
c
      do 140 j=js1,js2
        do 130 i=is1(j),is2(j)
c
c--5.1. FLUX FROM i to i+1 IF u.GT.0.0.
c--------------------------------------
c
	  ii       = i-1
          sm(i,j)  = bet(ii,j)*(sm(i,j)+fm(ii,j))
     &               +(1.0-bet(ii,j))*sm(i,j)
          alf      = bet(ii,j)*fm(ii,j)/sm(i,j)
          alf1     = 1.0-alf
          temp     = alf*s0(i,j)-alf1*f0(ii,j)
          s0(i,j)  = bet(ii,j)*(s0(i,j)+f0(ii,j))
     &               +(1.0-bet(ii,j))*s0(i,j)
          sx(i,j)  = bet(ii,j)*(alf*fx(ii,j)
     &               +alf1*sx(i,j)+3.0*temp)+
     &               (1.0-bet(ii,j))*sx(i,j)
          sxx(i,j) = bet(ii,j)*(alf*alf*fxx(ii,j)+alf1*alf1*sxx(i,j)+
     &               5.0*(alf*alf1*(sx(i,j)
     &                    -fx(ii,j))-(alf1-alf)*temp))+
     &		     (1.0-bet(ii,j))*sxx(i,j)
          sxy(i,j) = bet(ii,j)*(alf*fxy(ii,j)+alf1*sxy(i,j)+
     & 	             3.0*(-alf1*fy(ii,j)+alf*sy(i,j)))+
     &		     (1.0-bet(ii,j))*sxy(i,j)
          sy(i,j)  = bet(ii,j)*(sy(i,j)+fy(ii,j))
     &               +(1.0-bet(ii,j))*sy(i,j)
          syy(i,j) = bet(ii,j)*(syy(i,j)+fyy(ii,j))
     &               +(1.0-bet(ii,j))*syy(i,j)
130     continue
140   continue
c
      do 160 j=js1,js2
	do 150 i=is1(j),is2(j)
c
c--5.2. FLUX FROM i+1 to i IF u.LT.0.0.
c--------------------------------------
c
          sm(i,j)  = bet(i,j)*sm(i,j)
     &               +(1.0-bet(i,j))*(sm(i,j)+fm(i,j))
          alf      = (1.0-bet(i,j))*fm(i,j)/sm(i,j)
          alf1     = 1.0-alf
          temp     = -alf*s0(i,j)+alf1*f0(i,j)
          s0(i,j)  = bet(i,j)*s0(i,j)+
     &               (1.0-bet(i,j))*(s0(i,j)+f0(i,j))
          sx(i,j)  = bet(i,j)*sx(i,j)+
     &	             (1.0-bet(i,j))*(alf*fx(i,j)
     &                +alf1*sx(i,j)+3.0*temp)
          sxx(i,j) = bet(i,j)*sxx(i,j)+
     &  	     (1.0-bet(i,j))*(alf*alf*fxx(i,j)
     &               +alf1*alf1*sxx(i,j)+
     &               5.0*(alf*alf1*(-sx(i,j)+fx(i,j))+(alf1-alf)*temp))
          sxy(i,j) = bet(i,j)*sxy(i,j)+
     &	             (1.0-bet(i,j))*(alf*fxy(i,j)+alf1*sxy(i,j)+
     &               3.0*(alf1*fy(i,j)-alf*sy(i,j)))
          sy(i,j)  = bet(i,j)*sy(i,j)+(1.0-bet(i,j))*(sy(i,j)+fy(i,j))
          syy(i,j) = bet(i,j)*syy(i,j)
     &               +(1.0-bet(i,j))*(syy(i,j)+fyy(i,j))
150     continue
160   continue
c
c--Raccord cyclique
c
        call raccord(sm(1,1),zero,1,8)
        call raccord(s0(1,1),zero,1,8)
        call raccord(sx(1,1),zero,1,8)
        call raccord(sxx(1,1),zero,1,8)
        call raccord(sy(1,1),zero,1,8)
        call raccord(syy(1,1),zero,1,8)
        call raccord(sxy(1,1),zero,1,8)
      return
      end
C
