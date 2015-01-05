      program write
      implicit none
      
      integer i,j,k,l,m,ireg,im

c SW coefficienten die gedeclareerd en ingelezen moeten worden:
c index 1: maand
c index 2: regio; 1-11 land, 12-22 sea, 23-27 mountains
c index 3: level; 1 toa, 2 200hPa, 3 500hPa, 4 sfc
c index 4: up/down; 1 down, 2 up
c index 5: (alleen voor dfuocdrs) coefficienten van 3e orde polynoom fit

      real*8 cost_ref(12,27), salb_ref(12,27)
      real*8 fclrref(12,27,4,2)
      real*8 dfclrdcost(12,27,4,2)
      real*8 dfclrdrs1(12,27,4,2)
      real*8 fuocref(12,27,4,2)
      real*8 dfuocdcost(12,27,4,2)
      real*8 dfuocdrs(12,27,4,2,3)

C index 1: flux levels: 1 toa up, 2 200 up 3 500 up, 4 sfc up 5 toa down etc.
C index 2: regions
C index 3: month's
C index 4: 0 clear sky, 1 cloudy sky, except for surface albedo: 1, 2,3 
c          correspond to coefficients of 3rd order polynomial fit
      
      real*4 costref(27,12), salbref(27,12)
      real*4 swrref(8,27,12,0:1)
      real*4 swrcost(8,27,12,0:1)
      real*4 swrsalb(8,27,12,0:3)
      
      
      character ls
      integer  irind(27)
      DATA irind /22,11,21,10,20,9,19,8,18,7,17,6,16,5,15,4,14,3,13,2,
     *            12,1,27,26,25,24,23/
     
      open(10,file='cost_ref',form='formatted')
      read(10,*)
      do i=1, 27
         read(10,*) (cost_ref(m,i), m=1, 12)
      enddo
      close(10)

      open(10,file='salb_ref',form='formatted')
      read(10,*)
      do i=1, 27
         read(10,*) (salb_ref(m,i), m=1, 12)
      enddo
      close(10)

      open(10,file='fclrref',form='formatted')
      read(10,*)    fclrref
      close(10)

      open(10,file='dfclrdcost',form='formatted')
      read(10,*)    dfclrdcost
      close(10)

      open(10,file='dfclrdrs1',form='formatted')
      read(10,*)    dfclrdrs1
      close(10)

      open(10,file='fuocref',form='formatted')
      read(10,*)    fuocref
      close(10)

      open(10,file='dfuocdcost',form='formatted')
      read(10,*)    dfuocdcost
      close(10)

      open(10,file='dfuocdrs',form='formatted')
      read(10,*)    dfuocdrs
      close(10)



      do ireg=1,27
        do im=1,12 
        
          costref(ireg,im)=cost_ref(im,irind(ireg))
          salbref(ireg,im)=salb_ref(im,irind(ireg))
          do k=1,4
            swrref(k,ireg,im,0)=fclrref(im,irind(ireg),k,2)
            swrref(k,ireg,im,1)=fuocref(im,irind(ireg),k,2)
            swrcost(k,ireg,im,0)=dfclrdcost(im,irind(ireg),k,2)
            swrcost(k,ireg,im,1)=dfuocdcost(im,irind(ireg),k,2)
            swrsalb(k,ireg,im,0)=dfclrdrs1(im,irind(ireg),k,2)
            do l=1,3
              swrsalb(k,ireg,im,l)=dfuocdrs(im,irind(ireg),k,2,l)
            enddo
          enddo
          do k=5,8
            swrref(k,ireg,im,0)=fclrref(im,irind(ireg),k-4,1)
            swrref(k,ireg,im,1)=fuocref(im,irind(ireg),k-4,1)
            swrcost(k,ireg,im,0)=dfclrdcost(im,irind(ireg),k-4,1)
            swrcost(k,ireg,im,1)=dfuocdcost(im,irind(ireg),k-4,1)
            swrsalb(k,ireg,im,0)=dfclrdrs1(im,irind(ireg),k-4,1)
            do l=1,3
              swrsalb(k,ireg,im,l)=dfuocdrs(im,irind(ireg),k-4,1,l)
            enddo
          enddo
        enddo
      enddo
      
      open(2,file='swrref.dat',form='unformatted')
      
      write(2) costref,salbref
      write(2) swrref
      close(2)
      
      open(2,file='swrcoef.dat',form='unformatted')
      
      write(2) swrcost,swrsalb
      close(2)

      end
      
