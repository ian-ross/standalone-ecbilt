












      subroutine adv(dt,ut,vt,s0,s1)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
c    Avection methode Hibler.
C    Attention : pas conservatif au front de glace.
C    Pas toujours tres stable (instable si pad diffusion).
C
	include 'type.com'
	include 'para.com'
	include 'const.com'
	include 'bloc.com'
	include 'dynami.com'
c
	dimension ut(imax,jmax),vt(imax,jmax),s0(imax,jmax)
	dimension s1(imax,jmax)
	dimension div(imax,jmax),huc(imax,jmax),hvc(imax,jmax)
c
	do j=js1-1,js2
	   do i=is1(j)-1,is2(j)
	      ip1=i+1
              huc(i,j)=(s1(i,j)/area(i,j)*dxc1(ip1,j)
     &		      +s1(ip1,j)/area(ip1,j)*dxc1(i,j))
     &                /(dxc1(ip1,j)+dxc1(i,j))
              huc(i,j)=.5*(s1(i,j)/area(i,j)+s1(ip1,j)/area(ip1,j))
	   enddo
	enddo
c
	do j=js1-1,js2
	   jp1=j+1
 	   do i=is1(j)-1,is2(j)
	      hvc(i,j)=(s1(i,j)/area(i,j)*dxc2(i,jp1)
     & 	              +s1(i,jp1)/area(i,jp1)*dxc2(i,j))
     &                 /(dxc2(i,j)+dxc2(i,jp1))
 	      hvc(i,j)=.5*(s1(i,j)/area(i,j)+s1(i,jp1)/area(i,jp1))
	   enddo
	enddo
c
       do j=js1,js2
 	  jm1 = j-1
 	  do i=is1(j),is2(j)
 	     im1      = i-1
 	     div(i,j) = 2.0*
     &      (akappa(i,j,1,1)*(ut(i,j)*huc(i,j)-ut(im1,j)*huc(im1,j))-
     &        akappa(i,j,1,2)*(vt(i,jm1)*hvc(i,jm1)+vt(i,j)*hvc(i,j))-
     &        akappa(i,j,2,2)*(-vt(i,jm1)*hvc(i,jm1)+vt(i,j)*hvc(i,j))+
     &        akappa(i,j,2,1)*(ut(i,j)*huc(i,j)+ut(im1,j)*huc(im1,j)))
C	     div(i,j)  = bkappa(i,j,1,1)*ut(i,j)*huc(i,j)
C     &                  -bkappa(i,j,1,2)*ut(im1,j)*huc(im1,j)
C     &                  +bkappa(i,j,2,2)*vt(i,j)*hvc(i,j)
C     &                  -bkappa(i,j,2,1)*vt(i,jm1)*hvc(i,jm1)
	     s1(i,j)  = max(zero,s0(i,j)-area(i,j)*dt*div(i,j))
         enddo
       enddo
       do 200 jj=jcl1,jcl2
	   s1(ims1-1,jj) = s1(ims2,jj)
           s1(ims2+1,jj) = s1(ims1,jj)
200    continue
c
       return
       end
