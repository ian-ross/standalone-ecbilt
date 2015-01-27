c23456789012345678901234567890123456789012345678901234567890123456789012      
      program plotmaps
c-----------------------------------------------------------------------
c ***
c *** this program produces a postscript file of the world map with
c *** the masks of ecbilt superposed
c ***
c-----------------------------------------------------------------------

      implicit none

      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      real*4    xsc,ysc

      common /scal / xsc , ysc

      xsc=6d0/90d0
      ysc=8d0/90d0

      call rooster
      open(1,file='lakes.ps')
      call header
      call readgaust
      call mklsmaskt
      call mkgridt
      call mkwrld
      call mktextt
      call mklakent
      call trailer
      close(1)
      open(1,file='rivers.ps')
      call header
      call readgaust
      call mklsmaskt
      call mkgridt
      call mkwrld
      call mktextt
      call mkrunoft
      call trailer
      close(1)
      open(1,file='oceans.ps')
      call header
      call readgaust
      call mklsmaskt
      call mkgridt
      call mkwrld
      call mktextt
      call mkocbast
      call trailer
      close(1)
      open(1,file='worldmapv.ps')
      call header
      call readgausu
      call mklsmasku
      call mkgridu
      call mkwrld
      call mktextu
      call trailer
      close(1)
      open(1,file='icebirth.ps')
      call header
      call readgaust
      call mklsmaskt
      call mkgridt
      call mkwrld
      call mktextt
      call mkicebt
      call trailer
      close(1)
      open(1,file='icedeath.ps')
      call header
      call readgaust
      call mklsmaskt
      call mkgridt
      call mkwrld
      call mktextt
      call mkicedt
      call trailer
      close(1)
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine readgaust
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      integer   ilat,i,j
      real*4    g,ag(nlat),rlat(0:nlat+1),rlon(1:nlon+1),dlon,dfac

      common /gridt/ rlon,rlat

c *** gauss points

      dfac=90d0/asin(1d0)
       
      open(2,file='../inputdata/atmos/gauss.asc',
     *      status='old')
      ilat=nlat/2      
  10  continue
        read(2,220,end=20) g
        if (int(g).eq.ilat) then
          do i=1,int(g)
            read(2,220) ag(i)
          enddo
          goto 20
        else
          goto 10
        endif 
  20  continue
      close(2)
      
      do i=1,ilat
        rlat(i)=-ag(ilat+1-i)
        rlat(ilat+i)=ag(i)
      enddo

      do i=1,nlat
        rlat(i)=dfac*asin(rlat(i))
c        write(100,*) i,rlat(i)
      enddo
      rlat(0)=-180-rlat(1)
      rlat(nlat+1)=-rlat(0)

      dlon=360d0/nlon

      do i=1,nlon+1
        rlon(i)=(i-1)*dlon-180d0
      enddo

220   format (f18.10,f17.10)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine readgausu
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      integer   ilat,i,j
      real*4    g,ag(nlat),rlat(0:nlat+2),rlon(1:nlon+1),dlon,dfac

      common /gridu/ rlon,rlat

c *** gauss points

      dfac=90d0/asin(1d0)
       
      open(2,file='../inputdata/atmos/gauss.asc',
     *         status='old')
      ilat=nlat/2      
  10  continue
        read(2,220,end=20) g
        if (int(g).eq.ilat) then
          do i=1,int(g)
            read(2,220) ag(i)
          enddo
          goto 20
        else
          goto 10
        endif 
  20  continue
      close(2)
      
      do i=1,ilat
        rlat(i)=-ag(ilat+1-i)
        rlat(ilat+i)=ag(i)
      enddo

      do i=1,nlat
        rlat(i)=dfac*asin(rlat(i))
c        write(100,*) i,rlat(i)
      enddo

      do i=nlat,2,-1
        rlat(i)=rlat(i)-0.5*abs(rlat(i)-rlat(i-1))
      enddo

      rlat(nlat+1)= 90d0
      rlat(nlat+2)= 180d0-rlat(nlat)
      rlat(1) =    -90d0
      rlat(0) = -180 - rlat(2)

      dlon=360d0/nlon

      do i=1,nlon+1
        rlon(i)=(i-1.5)*dlon-180d0
      enddo

220   format (f18.10,f17.10)
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mklsmaskt
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      integer   ilat,i,j,k,indu(0:nlon,0:nlat),indt(0:nlon,0:nlat)
      integer   lsmask(1:nlon+1,1:nlat),inx
      real*4    rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,hdx,x1,x2,y1,y2
      character*1 cn

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask

c *** landsea mask

      open(2,file='../inputdata/ocean/mask.dat',
     *     status='old')

      read (2,120) 
      do j=nlat,0,-1
        read(2,120) (indu(i,j),i=0,nlon)
      enddo
      do j=nlat-1,0,-1
        read (2,310) k,(indt(i,j),i=0,nlon-1)
      enddo

      do i=1,nlat
        do j=1,nlon
          lsmask(j,i)=indt(j-1,i-1)
        enddo
        lsmask(nlon+1,i)=lsmask(1,i)
      enddo
      close(2)

      open(2,file='../inputdata/ocean/lakemask.dat',
     *     status='old')

      do j=nlat-1,0,-1
        read (2,310) k,(indt(i,j),i=0,nlon-1)
      enddo

      do i=1,nlat
        do j=1,nlon
          if (indt(j-1,i-1).gt.1) then
            if (lsmask(j,i).eq.0) then
              lsmask(j,i)=indt(j-1,i-1)
            else
              write(*,*) 'error in reading lakemask'
              write(*,*) i-1,j-1,' already sea'
              stop
            endif
          endif
        enddo
        lsmask(nlon+1,i)=lsmask(1,i)
      enddo

      do j=nlat-1,0,-1
        write (2,310) j,(lsmask(i+1,j+1),i=0,nlon-1)
      enddo
      close(2)

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.1 scalefont'
      write(1,200) 'setfont'

      hdx=(rlon(2)-rlon(1))/2.
      do i=1,nlon+1
        x1=rlon(i)-hdx
        x2=rlon(i)+hdx
        do j=1,nlat
          y1=rlat(j)-(rlat(j)-rlat(j-1))/2.
          y2=rlat(j)+(rlat(j+1)-rlat(j))/2.
          write (1,200) 'newpath'
          write (1,100) rxc(x1), ryc(y1),' M'
          write (1,100) rxc(x2), ryc(y1),' L'
          write (1,100) rxc(x2), ryc(y2),' L'
          write (1,100) rxc(x1), ryc(y2),' L'
          write (1,200) 'closepath'
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (lsmask(inx,j).eq.0) then
            write (1,200) '1. setgray'
          else
            if (lsmask(inx,j).ge.2) then
              write (1,200) '.95 setgray'
            else
              write (1,200) '0.8 setgray'
            endif
          endif
          write (1,200) 'fill'
          write (1,200) 'stroke'
          if (lsmask(inx,j).ge.2) then
            write (cn,300) lsmask(inx,j)
            write (1,100) rxc(x1), ryc(y2),' M'
            write (1,200) '('//cn//') show'
          endif
        enddo
      enddo
      
 100  format (2F7.2,A)
 200  format (A)
 300  format (i1)
220   format (f18.10,f17.10)
120   format(65i1)
310   format(i4,i2,90i1)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mklsmasku
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      integer   ilat,i,j,k,indu(0:nlon,0:nlat),indt(0:nlon,0:nlat)
      integer   lsmask(1:nlon+1,0:nlat+1),inx
      real*4    rlat(0:nlat+2),rlon(1:nlon+1),rxc,ryc,hdx,x1,x2,y1,y2

      common /gridu/ rlon,rlat

c *** landsea mask

      open(2,file='../inputdata/ocean/mask.dat',
     *     status='old')

      read (2,120) 
      do j=nlat,0,-1
        read(2,120) (indu(i,j),i=0,nlon)
      enddo

      do i=1,nlat+1
        do j=1,nlon
          lsmask(j,i)=indu(j-1,i-1)
        enddo
        lsmask(nlon+1,i)=lsmask(1,i)
      enddo

      close(2)


      hdx=(rlon(2)-rlon(1))/2.
      do i=1,nlon+1
        x1=rlon(i)-hdx
        x2=rlon(i)+hdx
        do j=1,nlat+1
          y1=rlat(j)-(rlat(j)-rlat(j-1))/2.
          y2=rlat(j)+(rlat(j+1)-rlat(j))/2.
          write (1,200) 'newpath'
          write (1,100) rxc(x1), ryc(y1),' M'
          write (1,100) rxc(x2), ryc(y1),' L'
          write (1,100) rxc(x2), ryc(y2),' L'
          write (1,100) rxc(x1), ryc(y2),' L'
          write (1,200) 'closepath'
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (lsmask(inx,j).eq.0) then
            write (1,200) '.8 setgray'
          else
            write (1,200) '1. setgray'
          endif
          write (1,200) 'fill'
          write (1,200) 'stroke'
        enddo
      enddo
      write (1,200) '0. setgray'
      
 100  format (2F7.2,A)
 200  format (A)
220   format (f18.10,f17.10)
120   format(65i1)
310   format(i4,i2,90i1)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkgridt
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      integer   i,j
      real*4    rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,hdx,hdy

      common /gridt/ rlon,rlat

      
      hdx=(rlon(2)-rlon(1))/2.
      hdy=(rlat(1)-rlat(0))/2.

      write (1,200) '0 setgray'
      write (1,200) '0.01 setlinewidth'

      do i=1,nlon+1
        write(1,200) 'newpath'
        write(1,100) rxc(rlon(i)),  ryc(rlat(1)-hdy),' M'
        write(1,100) rxc(rlon(i)),  ryc(rlat(nlat)+hdy),' L'
        write(1,200) 'stroke'
      enddo

      do i=1,nlat
        write(1,200) 'newpath'
        write(1,100) rxc(rlon(1)-hdx),   ryc(rlat(i)),' M'
        write(1,100) rxc(rlon(nlon+1)+hdx),ryc(rlat(i)),' L'
        write(1,200) 'stroke'
      enddo
      

 100  format (2F7.2,A2)
 200  format (A)
220   format (f18.10,f17.10)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkgridu
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)
      integer   i,j
      real*4    rlat(0:nlat+2),rlon(1:nlon+1),rxc,ryc,hdx,hdy

      common /gridu/ rlon,rlat

      
      hdx=(rlon(2)-rlon(1))/2.
      hdy=(rlat(1)-rlat(0))/2.

      write (1,200) '0 setgray'
      write (1,200) '0.01 setlinewidth'

      do i=1,nlon+1
        write(1,200) 'newpath'
        write(1,100) rxc(rlon(i)),  ryc(rlat(1)-hdy),' M'
        write(1,100) rxc(rlon(i)),  ryc(rlat(nlat+1)+hdy),' L'
        write(1,200) 'stroke'
      enddo

      do i=1,nlat+1
        write(1,200) 'newpath'
        write(1,100) rxc(rlon(1)-hdx),   ryc(rlat(i)),' M'
        write(1,100) rxc(rlon(nlon+1)+hdx),ryc(rlat(i)),' L'
        write(1,200) 'stroke'
      enddo
      

 100  format (2F7.2,A2)
 200  format (A)
220   format (f18.10,f17.10)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mktextt
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j
      real*4      rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,dx,dy
      integer     lsmask(1:nlon+1,1:nlat),inx
      real*4      x1,y1
      character*1 cn
      character*3 fn
      character*4 fm

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'
      
      dx=-10.
      dy=-1.
      do i=1,nlat
        write(fn,300) i
        write(1,100) rxc(rlon(1)+dx),ryc(rlat(i)+dy),' M'
        write(1,200) '('//fn//') show'
      enddo

      dx=5.
      dy=-1.
      do i=1,nlat
        write(fm,400) nint(rlat(i))
        write(1,100) rxc(rlon(nlon+1)+dx),ryc(rlat(i)+dy),' M'
        write(1,200) '('//fm//') show'
      enddo

      dx=-3.5
      dy=-2.

      do i=1,nlon/2+1
        j=i+nlon/2
        write(fn,300) i
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(0)+dy),' M'
        write(1,200) '('//fn//') show'
      enddo

      do i=nlon/2+1,nlon
        j=i-nlon/2
        write(fn,300) i
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(0)+dy),' M'
        write(1,200) '('//fn//') show'
      enddo

      dx=-4.
      dy=7.

      do i=1,nlon/2+1,2
        j=i+nlon/2
        write(fm,400) nint(rlon(j))
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(nlat)+dy),' M'
        write(1,200) '('//fm//') show'
      enddo

      do i=nlon/2+1,nlon+1,2
        j=i-nlon/2
        write(fm,400) nint(rlon(j))
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(nlat)+dy),' M'
        write(1,200) '('//fm//') show'
      enddo


 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mklakent
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j
      real*4      rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,dx,dy
      integer     lsmask(1:nlon+1,1:nlat),inx
      real*4      x1,y1
      character*1 cn
      character*3 fn
      character*4 fm

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask


      write(1,200) '/Times-Roman findfont'
      write(1,200) '.8 scalefont'
      write(1,200) 'setfont'

      write(1,100) rxc(-50.),ryc(100.),' M'
      write(1,200) '(Land-sea mask and lakes) show'

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'

      dx=(rlon(2)-rlon(1))/6.
      do i=1,nlon+1
        x1=rlon(i)-dx
        do j=1,nlat
          y1=rlat(j)-(rlat(j)-rlat(j-1))/6.
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (lsmask(inx,j).ge.2) then
            write(1,200) '.95 setgray'
            write(1,200) 'newpath'
            write(1,600) rxc(rlon(i)),ryc(rlat(j)),rxc(2.),
     *                   0.,360.,' arc fill'
            write(1,200) '0. setgray'
            write(1,100) rxc(rlon(i)), ryc(rlat(j)),' M'
            write(cn,500) lsmask(inx,j)
            write(1,200) '('//cn//') stringwidth pop'
            write(1,200) '2 div neg dup'
            write(1,800) ryc(-1.),' lt'
            write(1,700) '{dup rmoveto} {',ryc(-1.),' rmoveto} ifelse'
            write(1,200) '('//cn//') show'
          endif
        enddo
      enddo

 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
 500  format (i1)
 600  format (5F7.2,A)
 700  format (A,F7.2,A)
 800  format (F7.2,A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkrunoft
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j,k
      real*4      rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,dx,dy
      integer     lsmask(1:nlon+1,1:nlat),inx
      real*4      x1,y1
      character*1 cn,ch(1:nlon+1,1:nlat)
      character*3 fn
      character*4 fm

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask

      open(2,file='../inputdata/land/labas.dat',status='old')

      do i=nlat,1,-1
        read(2,*)
      enddo

      do i=nlat,1,-1
        read(2,*)
      enddo

      do i=nlat,1,-1
        read(2,50) k,cn,(ch(j,i),j=1,nlon)
        ch(nlon+1,i)=ch(1,i)
      enddo

      close(2)

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.8 scalefont'
      write(1,200) 'setfont'

      write(1,100) rxc(-50.),ryc(100.),' M'
      write(1,200) '(Land-sea mask and runoff basins) show'

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'

      dx=(rlon(2)-rlon(1))/6.
      do i=1,nlon+1
        x1=rlon(i)-dx
        do j=1,nlat
          y1=rlat(j)-(rlat(j)-rlat(j-1))/6.
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (ch(inx,j).ne.'0') then
            write(1,200) '1. setgray'
            write(1,200) 'newpath'
            write(1,600) rxc(rlon(i)),ryc(rlat(j)),rxc(2.),
     *                   0.,360.,' arc fill'
            write(1,200) '0. setgray'
            write (1,100) rxc(rlon(i)), ryc(rlat(j)),' M'
            cn=ch(inx,j)
            write(1,200) '('//cn//') stringwidth pop'
            write(1,200) '2 div neg dup'
            write(1,800) ryc(-1.),' lt'
            write(1,700) '{dup rmoveto} {',ryc(-1.),' rmoveto} ifelse'
            write (1,200) '('//cn//') show'
          endif
        enddo
      enddo

  50  format(i4,65A1)
 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
 500  format (i1)
 600  format (5F7.2,A)
 700  format (A,F7.2,A)
 800  format (F7.2,A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkocbast
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j,k
      real*4      rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,dx,dy
      integer     lsmask(1:nlon+1,1:nlat),inx
      real*4      x1,y1
      character*1 cn,ch(1:nlon+1,1:nlat)
      character*3 fn
      character*4 fm

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask

      open(2,file='../inputdata/ocean/ocbas.dat',status='old')

      do i=nlat,0,-1
        read(2,*)
      enddo

      do i=nlat,1,-1
        read(2,50) k,cn,(ch(j,i),j=1,nlon)
        ch(nlon+1,i)=ch(1,i)
      enddo

      close(2)

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.8 scalefont'
      write(1,200) 'setfont'

      write(1,100) rxc(-50.),ryc(100.),' M'
      write(1,200) '(Land-sea mask and ocean basins) show'

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'

      dx=(rlon(2)-rlon(1))/6.
      do i=1,nlon+1
        x1=rlon(i)-dx
        do j=1,nlat
          y1=rlat(j)-(rlat(j)-rlat(j-1))/6.
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (ch(inx,j).ne.'0') then
            write(1,200) '1. setgray'
            write(1,200) 'newpath'
            write(1,600) rxc(rlon(i)),ryc(rlat(j)),rxc(2.),
     *                   0.,360.,' arc fill'
            write(1,200) '0. setgray'
            write (1,100) rxc(rlon(i)), ryc(rlat(j)),' M'
            cn=ch(inx,j)
            write(1,200) '('//cn//') stringwidth pop'
            write(1,200) '2 div neg dup'
            write(1,800) ryc(-1.),' lt'
            write(1,700) '{dup rmoveto} {',ryc(-1.),' rmoveto} ifelse'
            write (1,200) '('//cn//') show'
          endif
        enddo
      enddo

  50  format(i4,65A1)
 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
 500  format (i1)
 600  format (5F7.2,A)
 700  format (A,F7.2,A)
 800  format (F7.2,A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkicebt
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j,k,ia0,ia9,ica,icd,iap
      real*4      rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,dx,dy
      integer     lsmask(1:nlon+1,1:nlat),inx,lsb(nlat,nlon+1)
      real*4      x1,y1
      character*1 cn,ch(nlat,nlon+1),space
      character*3 fn
      character*4 fm

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask


      write(1,200) '/Times-Roman findfont'
      write(1,200) '.8 scalefont'
      write(1,200) 'setfont'

      write(1,100) rxc(-50.),ryc(100.),' M'
      write(1,200) '(Prescribed seaice birth) show'

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'


      ia0=ichar('0') 
      ia9=ichar('9') 
      ica=ichar('A') 
      icd=ichar('D') 
      iap=ichar('.') 

      open(2,file='../inputdata/ocfix/icemask.dat')

      do i=nlat,1,-1
        read (2,10) k,space,(ch(i,j),j=1,nlon)
        do j=1,nlon
          lsb(i,j)=ichar(ch(i,j))
          if (lsb(i,j).ge.ia0.and.lsb(i,j).le.ia9) then
            lsb(i,j)=lsb(i,j)-ia0+1
          else
            if (lsb(i,j).ge.ica.and.lsb(i,j).le.icd) then
              lsb(i,j)=lsb(i,j)-ica+10
            else
              lsb(i,j)=0
            endif
          endif
        enddo
        lsb(i,nlon+1)=lsb(i,1)
        ch(i,nlon+1)=ch(i,1)
      enddo

      close(2)

      dx=(rlon(2)-rlon(1))/6.
      do i=1,nlon+1
        x1=rlon(i)-dx
        do j=1,nlat
          y1=rlat(j)-(rlat(j)-rlat(j-1))/6.
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (lsb(j,inx).ge.1) then
            write(1,200) '.95 setgray'
            write(1,200) 'newpath'
            write(1,600) rxc(rlon(i)),ryc(rlat(j)),rxc(2.),
     *                   0.,360.,' arc fill'
            write(1,200) '0. setgray'
            write(1,100) rxc(rlon(i)), ryc(rlat(j)),' M'
            cn=ch(j,inx)
            write(1,200) '('//cn//') stringwidth pop'
            write(1,200) '2 div neg dup'
            write(1,800) ryc(-1.),' lt'
            write(1,700) '{dup rmoveto} {',ryc(-1.),' rmoveto} ifelse'
            write(1,200) '('//cn//') show'
          endif
        enddo
      enddo

 10   format(i4,65A1)
 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
 500  format (i1)
 600  format (5F7.2,A)
 700  format (A,F7.2,A)
 800  format (F7.2,A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkicedt
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j,k,ia0,ia9,ica,icd,iap
      real*4      rlat(0:nlat+1),rlon(1:nlon+1),rxc,ryc,dx,dy
      integer     lsmask(1:nlon+1,1:nlat),inx,lsb(nlat,nlon+1)
      real*4      x1,y1
      character*1 cn,ch(nlat,nlon+1),space
      character*3 fn
      character*4 fm

      common /gridt/ rlon,rlat
      common /lmaskt/ lsmask


      write(1,200) '/Times-Roman findfont'
      write(1,200) '.8 scalefont'
      write(1,200) 'setfont'

      write(1,100) rxc(-50.),ryc(100.),' M'
      write(1,200) '(Prescribed seaice death) show'

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'


      ia0=ichar('0') 
      ia9=ichar('9') 
      ica=ichar('A') 
      icd=ichar('D') 
      iap=ichar('.') 

      open(2,file='../inputdata/ocfix/icemask.dat')

      do i=nlat,1,-1
        read (2,10) k,space,(ch(i,j),j=1,nlon)
      enddo
      do i=nlat,1,-1
        read (2,10) k,space,(ch(i,j),j=1,nlon)
        do j=1,nlon
          lsb(i,j)=ichar(ch(i,j))
          if (lsb(i,j).ge.ia0.and.lsb(i,j).le.ia9) then
            lsb(i,j)=lsb(i,j)-ia0+1
          else
            if (lsb(i,j).ge.ica.and.lsb(i,j).le.icd) then
              lsb(i,j)=lsb(i,j)-ica+10
            else
              lsb(i,j)=0
            endif
          endif
        enddo
        lsb(i,nlon+1)=lsb(i,1)
        ch(i,nlon+1)=ch(i,1)
      enddo

      close(2)

      dx=(rlon(2)-rlon(1))/6.
      do i=1,nlon+1
        x1=rlon(i)-dx
        do j=1,nlat
          y1=rlat(j)-(rlat(j)-rlat(j-1))/6.
          inx=i+nlon/2
          if (inx.ge.nlon+1) inx=inx-nlon
          if (lsb(j,inx).ge.1) then
            write(1,200) '.95 setgray'
            write(1,200) 'newpath'
            write(1,600) rxc(rlon(i)),ryc(rlat(j)),rxc(2.5),
     *                   0.,360.,' arc fill'
            write(1,200) '0. setgray'
            write(1,100) rxc(rlon(i)), ryc(rlat(j)),' M'
            cn=ch(j,inx)
            write(1,200) '('//cn//') stringwidth pop'
            write(1,200) '2 div neg dup'
            write(1,800) ryc(-1.),' lt'
            write(1,700) '{dup rmoveto} {',ryc(-1.),' rmoveto} ifelse'
            write(1,200) '('//cn//') show'
          endif
        enddo
      enddo

 10   format(i4,65A1)
 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
 500  format (i1)
 600  format (5F7.2,A)
 700  format (A,F7.2,A)
 800  format (F7.2,A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mktextu
      
      implicit none
      
      integer   nlat,nlon
      parameter (nlon=64,nlat=32)

      integer     i,j
      real*4      rlat(0:nlat+2),rlon(1:nlon+1),rxc,ryc,dx,dy
      character*3 fn
      character*4 fm

      common /gridu/ rlon,rlat

      write(1,200) '/Times-Roman findfont'
      write(1,200) '.8 scalefont'
      write(1,200) 'setfont'

      write(1,100) rxc(-50.),ryc(100.),' M'
      write(1,200) '(Land-sea mask at v-points) show'


      write(1,200) '/Times-Roman findfont'
      write(1,200) '.3 scalefont'
      write(1,200) 'setfont'
      
      dx=-10.
      dy=-1.
      do i=1,nlat+1
        write(fn,300) i
        write(1,100) rxc(rlon(1)+dx),ryc(rlat(i)+dy),' M'
        write(1,200) '('//fn//') show'
      enddo

      dx=5.
      dy=-1.
      do i=1,nlat+1
        write(fm,400) nint(rlat(i))
        write(1,100) rxc(rlon(nlon+1)+dx),ryc(rlat(i)+dy),' M'
        write(1,200) '('//fm//') show'
      enddo

      dx=-3.5
      dy=-2.

      do i=1,nlon/2+1
        j=i+nlon/2
        write(fn,300) i
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(0)+dy),' M'
        write(1,200) '('//fn//') show'
      enddo

      do i=nlon/2+1,nlon
        j=i-nlon/2
        write(fn,300) i
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(0)+dy),' M'
        write(1,200) '('//fn//') show'
      enddo

      dx=-4.
      dy=5.

      do i=1,nlon/2+1,2
        j=i+nlon/2
        write(fm,400) nint(rlon(j))
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(nlat+1)+dy),' M'
        write(1,200) '('//fm//') show'
      enddo

      do i=nlon/2+1,nlon+1,2
        j=i-nlon/2
        write(fm,400) nint(rlon(j))
        write(1,100) rxc(rlon(j)+dx),ryc(rlat(nlat+1)+dy),' M'
        write(1,200) '('//fm//') show'
      enddo

 100  format (2F7.2,A)
 200  format (A)
 300  format (i3)
 400  format (i4)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine mkwrld
      
      implicit none
      
      integer          i,j,nlen

      parameter (nlen=28000)
      logical   far,eor,nh,newp
      integer*2 x(nlen),m
      integer   kk,ll,nyp
      real*4    lon,lat,dlon,dlat,v(2),olon,reflon,rxc,ryc,dmax
      data dmax/32768.0/
      data m/-32768/
     

      dlon=-1.0/dmax
      dlat=-1.0/dmax
      reflon=0.
      nyp=0
      open(unit=2,file='wereld.asc')
      read(2,*) x
      close(2)
      kk=x(2)

      write (1,200) '0.02 setlinewidth'
      write (1,200) '0 setgray'
      write (1,200) 'newpath'
      newp=.true.

      nh=x(4).lt.0.
      lon=real(x(3))*dlon
      lon=(lon+1.)/2.
      lon=lon+reflon
      if (lon.ge.1.) then 
        lon=lon-1.
      endif
      olon=2.*lon-1.
      do i=3,kk,2
        lon=real(x(i))*dlon
        lon=(lon+1.)/2.
        lon=lon+reflon
        if (lon.ge.1.) then 
          lon=lon-1.
        endif
        lon=2.*lon-1.
        if (abs(olon-lon).gt.1.0) then
          write (1,200) 'stroke'
          write (1,200) 'newpath'
          newp=.true.
        endif
        olon=lon
        lat=real(x(i+1))*dlat
        if (nh.and.(lat.lt.0).and.(nyp.eq.1)) then
          nh=.false.
          write (1,200) 'stroke'
        end if
        if (.not.nh.and.(lat.gt.0).and.(nyp.eq.1)) then
          nh=.true.
          write (1,200) 'newpath'
          newp=.true.
        end if
        far=abs(x(i)-x(i-2)).gt.1000
        eor=x(i).eq.m
        if (far.or.eor)  then
          write (1,200) 'stroke'
          write (1,200) 'newpath'
          newp=.true.
        endif
        if (.not.eor) then
          lon=lon*180.
          lat=lat*90.
          if (newp) then
            newp=.false.
            write(1,100) rxc(lon),ryc(lat),' M'
          else
            write(1,100) rxc(lon),ryc(lat),' L'
          endif
        endif
      enddo


 100  format (2F7.2,A2)
 200  format (A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      function rxc(x)
      implicit none

      real*4   rxc,xsc,ysc,x

      common /scal/ xsc, ysc
      rxc=xsc*x
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      function ryc(y)
      implicit none

      real*4   ryc,xsc,ysc,y

      common /scal/ xsc, ysc
      ryc=ysc*y
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine header
      write (1,100) '%!PS-Adobe-1'
      write (1,100) '28.35 28.35 scale'
      write (1,100) '10.5 14.9 translate'
      write (1,100) '90 rotate'
      write (1,100) '/M {moveto} def'
      write (1,100) '/L {lineto} def'
      write (1,100) '/R {rlineto} def'
      write (1,100) '/S {stroke} def'
      write (1,100) '/W {setdash} def'
  100 format(A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012      
      subroutine trailer
      write (1,100) 'stroke'
      write (1,100) 'showpage'
      write (1,100) '%%Trailer'
      write (1,100) '%%Pages: 1'
  100 format(A)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine rooster
c-----------------------------------------------------------------------
c *** subroutine reads land-sea mask
c-----------------------------------------------------------------------
      implicit none

      integer   i,j,l,itel,ihulp,jhulp,nl2,nltot,nlhalf
      integer   nl,nm
      parameter (nl=63,nm=31)
      integer   indu(0:nl+1,0:nm+1),isea(-1:nl+1,-1:nm+1)

c *** reading land sea mask

      open(2,file='../inputdata/ocean/mask.dat',
     *     status='old')
      read (2,120) 
      do j=nm+1,0,-1
        read(2,120) (indu(i,j),i=0,nl+1)
      enddo

      do i=0,nl+1
        do j=0,nm+1
          if (indu(i,j).eq.0) then
            indu(i,j)=1
          else
            indu(i,j)=0
          endif
        enddo
      enddo

      do j=0,nm+1
        if (indu(0,j).ne.indu(nl+1,j)) then
          write(*,*) 'error in reading mask'
          write(*,*) j,indu(0,j),indu(nl+1,j)
          stop
        endif
      enddo
       
      do j=0,nm       
        do i=0,nl
          isea(i,j)=0
          do ihulp=0,1
            do jhulp=0,1
              if (indu(i+ihulp,j+jhulp).eq.1) isea(i,j)=1
            enddo
          enddo
        enddo
        isea(nl+1,j)=isea(0,j)
        isea(-1,j)=isea(nl,j)
      enddo

c *** meridional cyclic conditions for isea 

      nltot=nl+1
      nlhalf=nltot/2

      do i=-1,nlhalf
        isea(i,-1)=isea(i+nlhalf,0)
        isea(i+nlhalf,-1)=isea(i,0)
        isea(i,nm+1)=isea(i+nlhalf,nm)
        isea(i+nlhalf,nm+1)=isea(i,nm)
      enddo

      do j=nm,0,-1
        write (2,310) j,(isea(i,j),i=0,nl)
      enddo

      close(2)

120   format(65i1)
310   format(i4,i2,90i1)


      return 
      end



          
