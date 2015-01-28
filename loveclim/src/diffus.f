












      subroutine diffus(dt,difhx,difhy,fld0,fld1)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----
c
	include 'type.com'
	include 'para.com'
	include 'const.com'
	include 'bloc.com'
        include 'dynami.com'
c
	dimension difhx(imax,jmax),difhy(imax,jmax),
     &	          fld0(imax,jmax),fld1(imax,jmax)
c
	dimension fld2(imax,jmax)
        dimension grxh(imax,jmax),gryh(imax,jmax),
     &            div0(imax,jmax),div(imax,jmax),fact(imax,jmax)
c
c  Fully implicit (alfa = 1) or CRANK-NICHOLSON (alfa = 0.5),
c  or fully  explicit (alfa = 0.0).
c
	alfa = 0.5
c
c  Relaxation constant.
c
	om = 0.5
c
        do 10 j=1,jmax
	  do 10 i=1,imax
	    fld1(i,j) = fld0(i,j)
	    grxh(i,j) = 0.0
	    gryh(i,j) = 0.0
10      continue
c
	do 90 k=1,100
	  zindk = real(1/k)
	  do 20 j=js1-1,js2
	    jp1 = (j+1)
	    do 15 i=is1(j)-1,is2(j)
	      ip1       = (i+1)
	      grxh(i,j) = difhx(i,j)*(fld1(ip1,j)-fld1(i,j))
	      gryh(i,j) = difhy(i,j)*(fld1(i,jp1)-fld1(i,j))
15          continue
20        continue
	  do 30 j=js1,js2
	    jm1=j-1
	    do 25 i=is1(j),is2(j)
	      im1=i-1
	      fact(i,j) = bkappa(i,j,1,1)+bkappa(i,j,1,2)
     &		         +bkappa(i,j,2,2)+bkappa(i,j,2,1)	
c	      fact(i,j) = 2.0*
c     &                    (difhx(i,j)*(akappa(i,j,1,1)+akappa(i,j,2,1))+
c     &                     difhx(im1,j)*(akappa(i,j,1,1)-akappa(i,j,2,1))+
c     &                     difhy(i,jm1)*(-akappa(i,j,1,2)+akappa(i,j,2,2))+
c     &                     difhy(i,j)*(akappa(i,j,2,2)+akappa(i,j,1,2)))
25            continue
30        continue
c
          do 40 j=js1,js2
	    jm1 = j-1
	    do 35 i=is1(j),is2(j)
	      im1       = i-1
	      div(i,j)  = bkappa(i,j,1,1)*grxh(i,j)
     &         	         -bkappa(i,j,1,2)*grxh(im1,j)
     &                   +bkappa(i,j,2,2)*gryh(i,j)
     &                   -bkappa(i,j,2,1)*gryh(i,jm1)
c	      div(i,j)  = 2.0*
c     &                    (akappa(i,j,1,1)*(grxh(i,j)-grxh(im1,j))+
c     &                     akappa(i,j,1,2)*(gryh(i,jm1)+gryh(i,j))+
c     &                     akappa(i,j,2,2)*(gryh(i,j)-gryh(i,jm1))+
c     &                     akappa(i,j,2,1)*(grxh(i,j)+grxh(im1,j)))
	      div0(i,j) = zindk*div(i,j)+(1.0-zindk)*div0(i,j)
35          continue
40        continue
c
          do 50 j=js1,js2
	    do 50 i=is1(j),is2(j)
	      fldnw     = (fld0(i,j)+dt*
     &	                   (alfa*(div(i,j)+fact(i,j)*fld1(i,j))+
     &		            (1.0-alfa)*div0(i,j)))/
     &	      	          (1.0+alfa*dt*fact(i,j))
	      fld2(i,j) = fld1(i,j)+om*(fldnw-fld1(i,j))
 50       continue
c
	  amx =  0.0
	  do 60 j=js1,js2
	    do 60 i=is1(j),is2(j)
	      amx = max(amx,abs(fld2(i,j)-fld1(i,j)))
 60       continue

          do 70 j=js1,js2
	    do 70 i=is1(j),is2(j)
	      fld1(i,j) = fld2(i,j)
70        continue
	  do 80 jj=jcl1,jcl2
	     fld1(ims1-1,jj) = fld1(ims2,jj)
	     fld1(ims2+1,jj) = fld1(ims1,jj)
80        continue

	  if (amx.lt.2.0e-04) go to 100
90      continue
c
100     continue	
C       write(53,*) amx,k
	return
	end
