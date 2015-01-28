












      subroutine icdadv(xjour)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  This routine advcect sea ice properties.
c  
c  modif : 29/11/99

 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'moment.com'
      include 'comunit.h'
c


      dimension ut(imax,jmax),vt(imax,jmax)
      dimension sm(imax,jmax)
      dimension s1(imax,jmax)
      dimension amskx(imax,jmax),amsky(imax,jmax),
     &          difhx(imax,jmax),difhy(imax,jmax),
     &          fld0(imax,jmax),fld1(imax,jmax)
      dimension s0g(imax,jmax),s0n(imax,jmax),
     &          s0a(imax,jmax),s0c0(imax,jmax),
     &          s0c1(imax,jmax),s0c2(imax,jmax),s0st(imax,jmax)
Cage &         ,s0an(imax,jmax),s0ag(imax,jmax)

c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c--1. Initialisation                                                   |
c-----------------------------------------------------------------------
c
      jour=int(xjour)
      zeps0 = 1.0d-16
      zeps1 = 1.0d-20
      onelo =  one
      zerolo = zero
      do j=1,jmax
       do i=1,imax
        sm(i,j)=area(i,j)
       enddo
      enddo

      if (idyn.eq.1) then
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c--2.Advection of sea ice properties.                                  |
c-----------------------------------------------------------------------
c
c--2.1. Computation of velocities for advection.
c-----------------------------------------------
c
c vbord factor between 1 and 2 to take into account slip
c or no-slip boundary conditions.
c
        vbord = 1.0+(1.0-bound)
        do 80 j=1,jmax-1
          do 70 i=1,imax
            ip1     =  i+1-(imax-2)*(i/imax)
C           ut(i,j) =  0.5*(ug(ip1,j)+ug(ip1,j+1))
C           vt(i,j) =  0.5*(vg(i,j+1)+vg(ip1,j+1))
            ut(i,j) = (ug(ip1,j)+ug(ip1,j+1))/
     &          (max(tmu(ip1,j,ku2)+tmu(ip1,j+1,ku2),vbord))
            vt(i,j) =  (vg(i,j+1)+vg(ip1,j+1))/
     &          (max(tmu(i,j+1,ku2)+tmu(ip1,j+1,ku2),vbord))
70        continue
80	continue
	    
c
c  2.1.1. Boundary conditions at the top of the grids.
c
        do 90 i=1,imax
          ut(i,jmax) = 0.0
          vt(i,jmax) = 0.0
90      continue
c
c  2.1.2. CFL test for stability.	
c
 	cfl    = 0.0
 	do 110 j=1,jmax
 	  do 100 i=1,imax
	    ip1  = (i+1)-(imax-2)*(i/imax)
 	    cfl  = max(cfl,(abs(ut(i,j))*dts(ks2))/
     &             (0.5*(dxc1(i,j)+dxc1(ip1,j))))
100	  continue
110	continue
 	do 130 j=1,jmax-1
	  jp1 = (j+1)
 	  do 120 i=1,imax
 	    cfl = max(cfl,(abs(vt(i,j))*dts(ks2))/
     &            (0.5*(dxc2(i,j)+dxc2(i,jp1))))
120	  continue
130	continue
 	if (cfl.gt.0.5) then 
	  write(iuo+66,'(a31,i3,a14,f10.6)') 
     &	  'violation of cfl criterion the ',jour,'th day, cfl = ',cfl
	endif
c
c--2.2. Transported properties.
c-------------------------------
c       
        do 150 j=1,jmax
          do 140 i=1,imax
c
c  Snow volume.	
c
            s0n(i,j)  = hnm(i,j)*area(i,j)
c
c  Ice volume.
c
            s0g(i,j)  = hgm(i,j)*area(i,j)
c
c  Surface covered by ice.
c
            s0a(i,j)  = (1.0-albq(i,j))*area(i,j)
c
c  Heat content of the snow layer.
c
            s0c0(i,j) = (tbq(i,j,1)/tfsn)*s0n(i,j)
c
c  Heat content of the first ice layer.
c
            s0c1(i,j) = (tbq(i,j,2)/tfsg)*s0g(i,j)
c
c  Heat content of the second ice layer.
c
            s0c2(i,j) = (tbq(i,j,3)/tfsg)*s0g(i,j)
c
c  Heat reservoir for brine pockets.
c
            s0st(i,j) = (qstobq(i,j)/xlg)*s0a(i,j)
c
c  Age of snow
c
Cage        s0an(i,j) = agen(i,j)*s0n(i,j)
c
c  Age of ice
c
Cage        s0ag(i,j) = ageg(i,j)*s0g(i,j)
c
140       continue
150	continue
c
c--2.3. Calls to advection and diffusion routines.
c-------------------------------------------------
c
c  Advection. 
c	
c	If ice drift field is too fast, use an
c	appropriate time step for advection.
c
        nitad = 1+int(max(zero,sign(one,cfl-0.5)))
	usnit = 1.0/real(nitad)
	dta   = dts(ks2)
	if (iameth.eq.1) then
	  do k=1,nitad
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0g(i,j)
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0g,s1)
	    do j=1,jmax
	      do i=1,imax
	       s1(i,j) = 0.5*(s0g(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0g,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0g(i,j) = s1(i,j)
	     enddo
	    enddo
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0n(i,j)
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0n,s1)
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = 0.5*(s0n(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0n,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0n(i,j) = s1(i,j)
	      enddo
	    enddo
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0a(i,j)
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0a,s1)
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = 0.5*(s0a(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0a,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0a(i,j) = s1(i,j)
	      enddo
	    enddo
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0c0(i,j)
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0c0,s1)
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = 0.5*(s0c0(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0c0,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0c0(i,j) = s1(i,j)
	      enddo
	    enddo
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0c1(i,j)
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0c1,s1)
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = 0.5*(s0c1(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0c1,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0c1(i,j) = s1(i,j)
	      enddo
	    enddo
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0c2(i,j)
      	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0c2,s1)
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = 0.5*(s0c2(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0c2,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0c2(i,j) = s1(i,j)
	      enddo
	    enddo
c
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = s0st(i,j)
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0st,s1)
	    do j=1,jmax
	      do i=1,imax
	        s1(i,j) = 0.5*(s0st(i,j)+s1(i,j))
	      enddo
	    enddo
	    call adv(dta*usnit,ut,vt,s0st,s1)
	    do j=1,jmax
	      do i=1,imax
	        s0st(i,j) = s1(i,j)
	      enddo
	    enddo
c
Cage        do j=1,jmax
Cage          do i=1,imax
Cage            s1(i,j) = s0an(i,j)
Cage          enddo
Cage        enddo
Cage        call adv(dta*usnit,ut,vt,s0an,s1)
Cage        do j=1,jmax
Cage          do i=1,imax
Cage            s1(i,j) = 0.5*(s0an(i,j)+s1(i,j))
Cage          enddo
Cage        enddo
Cage        call adv(dta*usnit,ut,vt,s0an,s1)
Cage        do j=1,jmax
Cage          do i=1,imax
Cage            s0an(i,j) = s1(i,j)
Cage          enddo
Cage        enddo
c
Cage         do j=1,jmax
Cage          do i=1,imax
Cage            s1(i,j) = s0ag(i,j)
Cage          enddo
Cage        enddo
Cage        call adv(dta*usnit,ut,vt,s0ag,s1)
Cage        do j=1,jmax
Cage          do i=1,imax
Cage            s1(i,j) = 0.5*(s0ag(i,j)+s1(i,j))
Cage          enddo
Cage        enddo
Cage        call adv(dta*usnit,ut,vt,s0ag,s1)
Cage        do j=1,jmax
Cage          do i=1,imax
Cage            s0ag(i,j) = s1(i,j)
Cage          enddo
Cage        enddo
c
	  enddo
	else
	  if (mod(jour,2).eq.0) then
	    do k=1,nitad
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
c    
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
c
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
c
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
Cage          call advx(dta*usnit,ut,onelo,
Cage &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
Cage          call advy(dta*usnit,vt,zerolo,
Cage &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
Cage          call advx(dta*usnit,ut,onelo,
Cage &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
Cage          call advy(dta*usnit,vt,zerolo,
Cage &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
	    enddo
	  else
	    do k=1,nitad
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
c    
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
c
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
c
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
Cage          call advy(dta*usnit,vt,onelo,
Cage &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
Cage          call advx(dta*usnit,ut,zerolo,
Cage &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
Cage          call advy(dta*usnit,vt,onelo,
Cage &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
Cage          call advx(dta*usnit,ut,zerolo,
Cage &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
	    enddo
	  endif
	endif
c
c     extra amount of lead opening/closing due to shearing
c     deformation (hibler, 1980, 1984; flato and hibler, 1991,1995; 
c     harder and lemke, 1994) and ridging (schulkes, 1995).
c     also stern et al. (1995).
c
      do j=1,jmax
        do i=1,imax
          s0a(i,j) = s0a(i,j)-alcr(i,j)*dta*area(i,j)
        enddo
      enddo
c
c  Diffusion.
c  
      do j=js1-1,js2
	jp1 = j+1   
	do i=is1(j)-1,is2(j)
	  ip1        = (i+1)-(imax-2)*(i/imax) 
	  amskx(i,j) = (1.0-max(zero,sign(one,-s0a(i,j))))*
     &	  	       (1.0-max(zero,sign(one,-s0a(ip1,j))))
	  amsky(i,j) = (1.0-max(zero,sign(one,-s0a(i,j))))*
     &	  	       (1.0-max(zero,sign(one,-s0a(i,jp1))))
c
	enddo
      enddo
      do j=1,jmax
         do i=1,imax
            difhx(i,j) = 0.
            difhy(i,j) = 0.
         enddo
      enddo
      do j=js1-1,js2-1
	jp1 = (j+1)
	do i=is1(j)-1,is2(j)
	  ip1          = (i+1)
	  difhx(i,j)   = amskx(i,j)*dfhu(i,j)
     &                   /(0.5*(dxc1(i,j)+dxc1(ip1,j)))
	  difhy(i,j)   = amsky(i,j)*dfhv(i,j)
     &                   /(0.5*(dxc2(i,j)+dxc2(i,jp1)))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0g(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0g(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0n(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0n(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0a(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0a(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0c0(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0c0(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0c1(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0c1(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0c2(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0c2(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
      do j=1,jmax
	do i=1,imax
	  fld0(i,j) = s0st(i,j)/area(i,j)
	enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
	do i=1,imax
	  s0st(i,j) = max(zero,fld1(i,j)*area(i,j))
	enddo
      enddo
c
Cage  do j=1,jmax
Cage    do i=1,imax
Cage      fld0(i,j) = s0an(i,j)/area(i,j)
Cage    enddo
Cage  enddo
Cage  call diffus(dta,difhx,difhy,fld0,fld1)
Cage  do j=1,jmax
Cage    do i=1,imax
Cage      s0an(i,j) = fld1(i,j)*area(i,j)
Cage    enddo
Cage  enddo
c
Cage  do j=1,jmax
Cage    do i=1,imax
Cage      fld0(i,j) = s0ag(i,j)/area(i,j)
Cage    enddo
Cage  enddo
Cage  call diffus(dta,difhx,difhy,fld0,fld1)
Cage  do j=1,jmax
Cage    do i=1,imax
Cage      s0ag(i,j) = fld1(i,j)*area(i,j)
Cage    enddo
Cage  enddo

c
c--2.5. Up-dating and limitation of sea ice properties
c       after transport.
c------------------------------------------------------
c
        do 190 j=js1,js2
	  zindhe = real(max(0,isign(1,j-jeq)))
          do 180 i=is1(j),is2(j)
c
c  2.5.1. Recover mean values over the grid squares.
c
	    s0n(i,j)    = max(zero,s0n(i,j)/area(i,j))
	    s0g(i,j)    = max(zero,s0g(i,j)/area(i,j))
	    s0a(i,j)    = max(zero,s0a(i,j)/area(i,j))
	    s0c0(i,j)   = max(zero,s0c0(i,j)/area(i,j))
	    s0c1(i,j)   = max(zero,s0c1(i,j)/area(i,j))
	    s0c2(i,j)   = max(zero,s0c2(i,j)/area(i,j))
	    s0st(i,j)   = max(zero,s0st(i,j)/area(i,j))
Cage        s0an(i,j)   = max(zero,s0an(i,j)/area(i,j))
Cage        s0ag(i,j)   = max(zero,s0ag(i,j)/area(i,j))
c
c  2.5.2. Recover in situ values.
c
            zindb       = max(zero,sign(one,s0a(i,j)-1.0e-06))
            acrith      = 1.0-(zindhe*acrit(1)+(1.0-zindhe)*acrit(2))
            s0a(i,j)    = zindb*min(s0a(i,j),acrith)
            hnbq(i,j)   = zindb*(s0n(i,j)/max(s0a(i,j),zeps0))
            hgbq(i,j)   = zindb*(s0g(i,j)/max(s0a(i,j),zeps0))
            zindn       = max(zero,sign(one,hnbq(i,j)-1.0e-06))
            zindg       = max(zero,sign(one,hgbq(i,j)-1.0e-03))
            zindb       = max(zindn,zindg)
            s0a(i,j)    = zindb*s0a(i,j)
            albq_prev   = albq(i,j)
            albq(i,j)   = 1.0-s0a(i,j)
            alcd(i,j)   = (albq(i,j)-albq_prev)/dta
     &                   -alcr(i,j)
            hnbq(i,j)   = zindn*hnbq(i,j)
            hgbq(i,j)   = zindg*hgbq(i,j)
            zusvon      = 1.0/max(hnbq(i,j)*s0a(i,j),zeps0)
            zusvog      = 1.0/max(hgbq(i,j)*s0a(i,j),zeps0)
	    zignm       = max(zero,sign(one,hndif-hnbq(i,j)))
            tbq(i,j,1)  = zindn*(zignm*tbq(i,j,1)+(1.0-zignm)*
     &                    min(max(173.15*one,
     &                    tfsn*zusvon*s0c0(i,j)),tfu(i,j)))+
     &	                  (1.0-zindn)*tfu(i,j)
            tbq(i,j,2)  = zindg*min(max(173.15*one,tfsg*zusvog*
     &                    s0c1(i,j)),tfu(i,j))+
     &	                  (1.0-zindg)*tfu(i,j)
            tbq(i,j,3)  = zindg*min(max(173.15*one,tfsg*zusvog*
     &                    s0c2(i,j)),tfu(i,j))+
     &	                  (1.0-zindg)*tfu(i,j)
            qstobq(i,j) = zindb*xlg*s0st(i,j)/max(s0a(i,j),zeps0)
Cage        agen(i,j)   = zindn*s0an(i,j)*zusvon
Cage        ageg(i,j)   = zindg*s0ag(i,j)*zusvog
180       continue
190	continue
c
c  2.5.3. Raccord cyclique (not neccesary here)
c
	do 200 jj=jcl1,jcl2
           hnbq(ims1-1,jj) = hnbq(ims2,jj)
           hnbq(ims2+1,jj) = hnbq(ims1,jj)
           hgbq(ims1-1,jj) = hgbq(ims2,jj)
           hgbq(ims2+1,jj) = hgbq(ims1,jj)
           albq(ims1-1,jj) = albq(ims2,jj)
           albq(ims2+1,jj) = albq(ims1,jj)
           alcd(ims1-1,jj) = alcd(ims2,jj)
           alcd(ims2+1,jj) = alcd(ims1,jj)
           tbq(ims1-1,jj,1) = tbq(ims2,jj,1)
           tbq(ims2+1,jj,1) = tbq(ims1,jj,1)
	   tbq(ims1-1,jj,2) = tbq(ims2,jj,2)
	   tbq(ims2+1,jj,2) = tbq(ims1,jj,2)
	   tbq(ims1-1,jj,3) = tbq(ims2,jj,3)
	   tbq(ims2+1,jj,3) = tbq(ims1,jj,3)
           qstobq(ims1-1,jj) = qstobq(ims2,jj)
           qstobq(ims2+1,jj) = qstobq(ims1,jj)
Cage       agen(ims1-1,jj) = agen(ims2,jj)
Cage       agen(ims2+1,jj) = agen(ims1,jj)
Cage       ageg(ims1-1,jj) = ageg(ims2,jj)
Cage       ageg(ims2+1,jj) = ageg(ims1,jj)
 200    continue
c
      endif

      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine icdadv -
      end
