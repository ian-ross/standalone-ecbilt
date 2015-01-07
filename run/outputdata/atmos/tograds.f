c23456789012345678901234567890123456789012345678901234567890123456789012
      program tograds6

c-----------------------------------------------------------------------
c**
c** this program seperates the big dataset from climate model run into 
c** 2 small dataset for grads plotting program. T level variables into
c** one datafile for grads. U level variables into another datafile for
c** grads. 
c** to run the program : 
c** compile this program : f77 togradsv0.6.f
c** excute:                a.out inputdataset
c** 
c** This version differs from the version 0.5 in a way that it includes 
c** two dumy variables for  T leve variable and two dumy variables for
c** U level variables.
c**
c** written by xueli wang june, 1995. last revision sep. 1996.
c-----------------------------------------------------------------------
      implicit none
  
      integer nx,ny

      parameter(nx=32,ny=64)

      integer        i,j,k,fl,n1,n3,nday,itel,idateold,icode0
      integer        icode,ilevel,iyear,imonth,iday,ix,iy,ifield,idate
      integer        idate0
      integer        indy,indx
      real*4         x(nx,ny),y(nx,ny,10)
      integer        irec(3),level1(15),level3(15)
      character*8    cvar1(50),cvar3(50)
      character*256  fn,fn1,fn2
    

      character*20   title


      n1 = 1
      n3 = 1

      nday = 0
      indx = 0
      indy = 0

      level1(1)  = 1000
      level1(2)  = 650
      level1(3)  = 350

      level3(1) = 800
      level3(2) = 500
      level3(3) = 200

      call getarg(1,fn)
      if (fn .eq. '') then
        write(*,*)'enter input data file name'
        read(5,930) fn
        write(*,930) fn
      endif
      fl=index(fn,' ')-1
      fn1=fn(1:fl)//'grdstl.dat'
      fn2=fn(1:fl)//'grdsul.dat'
      
      open(999, file=fn, form='unformatted')

      open(110,file=fn1, form='unformatted', 
     &                access='direct',recl=ny*nx)

      open(112,file=fn2, form='unformatted',
     &                access='direct',recl=ny*nx)

 
      open(10,file=fn(1:fl)//'grdstl.ctl', form='formatted') 
      open(11,file=fn(1:fl)//'grdsul.ctl', form='formatted')


      itel     = 1
      nday     = 1
      idateold = 1111
      icode0   = 888
      do i = 1, 3
        irec(i) = 1
      enddo
 

      read(999,end=100) icode,ilevel,iyear,imonth,iday,ix,iy,ifield
      read(999,end=100) ((x(i,j),j=1,ny),i=1,nx)

      idate = 20000000+iyear*10000+imonth*100+iday

      if (itel.eq.1) then
        idate0 = idate
      endif

      call wridatax(itel,icode,x,irec,indx)

      call count1(itel,icode,ifield,cvar1,n1)
      call count3(itel,icode,ifield,cvar3,n3)
      indx      = indx + 1
      icode0    = icode
      idateold  = idate

  5   continue
      read(999,end=100) icode,ilevel,iyear,imonth,iday,ix,iy,ifield
      read(999,end=100) ((x(i,j),j=1,ny),i=1,nx)
c      write(*,*) icode,ilevel,iyear,imonth,iday,ix,iy,ifield
      itel = itel + 1

      idate = 20000000+iyear*10000+imonth*100+iday
      if (idate .ne. idateold) nday = nday + 1

      if (icode.eq.icode0) then
        if (ifield.eq.1.or.ifield.eq.0) then    
          call wridatax(itel,icode,x,irec,indx)
          indx = indx + 1 
          if (idate.eq.idate0.and.(indx.eq.1))
     &       then
            call count1(itel,icode,ifield,cvar1,n1)
            call count3(itel,icode,ifield,cvar3,n3)
          endif
        endif
        if (ifield.eq.2) then    
          indy = indy + 1
          do i = 1, nx
            do j= 1, ny
               y(i,j,indy) = x(i,j)
            enddo
          enddo
          if (idate.eq.idate0.and.(indy.eq.1))
     &        then
            call count1(itel,icode,ifield,cvar1,n1)
            call count3(itel,icode,ifield,cvar3,n3)
          endif
        endif
      endif
      if (icode.ne.icode0) then
       if(indy.gt.0) then
         call wridatay2(icode0,y,irec,indy)
       endif 
        indy = 0
        indx = 0

        if (ifield.eq.1.or.ifield.eq.0) then    
          call wridatax(itel,icode,x,irec,indx)
          indx = indx + 1

          if (idate.eq.idate0.and.(indx.eq.1))
     &       then
            call count1(itel,icode,ifield,cvar1,n1)
            call count3(itel,icode,ifield,cvar3,n3)
          endif
        endif
      endif
      idateold = idate
      icode0   = icode
      goto 5

100   continue
      if (ifield.eq.2) then
        call wridatay2(icode,y,irec,indy)
      endif 
      call wrctl(10,'1 level variables plus t and omega',1,level1,
     &           nday,cvar1,n1-1,fn1)
      call wrctl(11,'3 level variables',3,level3,nday,cvar3,
     &             n3-1,fn2)

    
  900   format(10e12.6)
  910   format(8i10)
  930   format(a4,4x,i2,a7,2x,a13)
  940   format(8I10)

  950   format(5i10)
      
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wrctl(iunit,title,nlevel,level,nday,cvar,nvar,ctlname)

      implicit none
  
      integer nx,ny

      parameter(nx=32,ny=64)

      integer       i,iunit,nlevel,nday,nvar,fl
      character*20  title
      character*256 ctlname
      integer       level(3)
      character*8   cvar(50)

      fl=index(ctlname,' ')-1
      write(iunit,900)'dset ^', ctlname(1:fl)
      write(iunit, '("undef 9.99e+10")')
      write(iunit,'("title ", a15)') title
      write(iunit,'("*")')
      write(iunit,'("xdef",i4, " linear",2f10.3)') ny, 0.0, 5.625
      write(iunit,'("*")')
      write(iunit,'("ydef",i4, " linear",2f10.3)') nx, -85.76, 5.54
      write(iunit,'("*")')
      if (nlevel.eq.1) then
          write(iunit,'("zdef",i3, " levels",50i5)') 3,1000,650,350
        else
          write(iunit,'("zdef",i3, " levels",50i5)') nlevel,
     &                       (level(i),i=1,nlevel)
      endif
        
      write(iunit,'("*")')
      write(iunit,930)'tdef',nday, ' linear', '1jan0001 1dy'
      write(iunit,'("*")')
      write(iunit,'("vars",i3)') nvar
     
      do i= 1, nvar
        if(cvar(i).eq.'t'    .or.cvar(i).eq.'td'    .or.
     &     cvar(i).eq.'omega'.or.cvar(i).eq.'omegad'.or.
     &     cvar(i).eq.'hforc'.or.cvar(i).eq.'hforcd'.or.
     &     cvar(i).eq.'dumt1'.or.cvar(i).eq.'dumt1d'.or.
     &     cvar(i).eq.'dumt2'.or.cvar(i).eq.'dumt2d') then
           write(iunit,'(a8, "  " i3,"  9   99")')cvar(i),3
        else
           write(iunit,'(a8, "  " i3,"  9   99")')cvar(i),nlevel
        endif
      enddo
      write(iunit,'("endvars")')

900   format(a6,a)
930   format(a4,4x,i6,a7,2x,a13)
   
      return
      end
       
c23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine count1(itel,icode,ifield,cvar1,n1) 

c------------------------------------------------------------------------------
c**
c** this routine counts variable's number for one level variables.
c** 
c**
c------------------------------------------------------------------------------
      implicit none

  
      integer nx,ny

      parameter(nx=32,ny=64)
      integer      itel,icode,ifield,n1
      real*4       x(nx,ny)
      character*8  cvar1(50)


      if(itel.eq.1) n1 = 1

        if(icode.eq.134.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "sp  "
          n1 = n1 + 1
        endif
        if(icode.eq.134.and.(ifield.eq.2)) then
          cvar1(n1) = "spd  "
          n1 = n1 + 1
        endif
        if(icode.eq.139.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "ts  "
          n1 = n1 + 1
        endif
        if(icode.eq.139.and.(ifield.eq.2)) then
          cvar1(n1) = "tsd  "
          n1 = n1 + 1
        endif
        if(icode.eq.140.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "bmu   "
          n1 = n1 + 1
        endif
        if(icode.eq.140.and.(ifield.eq.2)) then
          cvar1(n1) = "bmud   "
          n1 = n1 + 1
        endif
        if(icode.eq.308.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "bml   "
          n1 = n1 + 1
        endif
        if(icode.eq.308.and.(ifield.eq.2)) then
          cvar1(n1) = "bmld   "
          n1 = n1 + 1
        endif
        if(icode.eq.141.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "sdl   "
          n1 = n1 + 1
        endif
        if(icode.eq.141.and.(ifield.eq.2)) then
          cvar1(n1) = "sdld   "
          n1 = n1 + 1
        endif
        if(icode.eq.142.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) ="lsp "
          n1 = n1 + 1
        endif
        if(icode.eq.142.and.(ifield.eq.2)) then
          cvar1(n1) ="lspd "
          n1 = n1 + 1
        endif
        if(icode.eq.143.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "cp  "
          n1 = n1 + 1
        endif
        if(icode.eq.143.and.(ifield.eq.2)) then
          cvar1(n1) = "cpd  "
          n1 = n1 + 1
        endif
        if(icode.eq.146.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "shf "
          n1 = n1 + 1
        endif
        if(icode.eq.146.and.(ifield.eq.2)) then
          cvar1(n1) = "shfd "
          n1 = n1 + 1
        endif
        if(icode.eq.147.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "lhf "
          n1 = n1 + 1
        endif
        if(icode.eq.147.and.(ifield.eq.2)) then
          cvar1(n1) = "lhfd "
          n1 = n1 + 1
        endif
        if(icode.eq.133.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "q   "
          n1 = n1 + 1
        endif
        if(icode.eq.133.and.(ifield.eq.2)) then
          cvar1(n1) = "qd   "
          n1 = n1 + 1
        endif
        if(icode.eq.157.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "r"
          n1 = n1 + 1
        endif
        if(icode.eq.157.and.(ifield.eq.2)) then
          cvar1(n1) = "rd"
          n1 = n1 + 1
        endif
        if(icode.eq.159.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "uv10"
          n1 = n1 + 1
        endif
        if(icode.eq.159.and.(ifield.eq.2)) then
          cvar1(n1) = "uv10d"
          n1 = n1 + 1
        endif
        if(icode.eq.160.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "runoffo"
          n1 = n1 + 1
        endif
        if(icode.eq.160.and.(ifield.eq.2)) then
          cvar1(n1) = "runoffod"
          n1 = n1 + 1
        endif
        if(icode.eq.161.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "runoffl"
          n1 = n1 + 1
        endif
        if(icode.eq.161.and.(ifield.eq.2)) then
          cvar1(n1) = "runoffld"
          n1 = n1 + 1
        endif
        if(icode.eq.164.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "tcc"
          n1 = n1 + 1
        endif
        if(icode.eq.164.and.(ifield.eq.2)) then
          cvar1(n1) = "tccd"
          n1 = n1 + 1
        endif

        if(icode.eq.182.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "evap"
          n1 = n1 + 1
        endif
        if(icode.eq.182.and.(ifield.eq.2)) then
          cvar1(n1) = "evapd"
          n1 = n1 + 1
        endif
        if(icode.eq.174.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "palb  "
            n1 = n1 + 1
        endif
        if(icode.eq.174.and.(ifield.eq.2)) then
            cvar1(n1) = "palbd  "
            n1 = n1 + 1
        endif
        if(icode.eq.175.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "alb  "
            n1 = n1 + 1
        endif
        if(icode.eq.175.and.(ifield.eq.2)) then
            cvar1(n1) = "albd  "
            n1 = n1 + 1
        endif
        if(icode.eq.176.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "ssr  "
            n1 = n1 + 1
        endif
        if(icode.eq.176.and.(ifield.eq.2)) then
            cvar1(n1) = "ssrd"
            n1 = n1 + 1
        endif
        if(icode.eq.177.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "str  "
            n1 = n1 + 1
        endif
        if(icode.eq.177.and.(ifield.eq.2)) then
            cvar1(n1) = "strd  "
            n1 = n1 + 1
        endif
        if(icode.eq.178.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "tsr "
            n1 = n1 + 1
        endif
        if(icode.eq.178.and.(ifield.eq.2)) then
            cvar1(n1) = "tsrd  "
            n1 = n1 + 1
        endif
        if(icode.eq.179.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "ttr  "
            n1 = n1 + 1
        endif
        if(icode.eq.179.and.(ifield.eq.2)) then
            cvar1(n1) = "ttrd  "
            n1 = n1 + 1
        endif
        if(icode.eq.180.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "ustress"
          n1 = n1 + 1
        endif
        if(icode.eq.180.and.(ifield.eq.2)) then
          cvar1(n1) = "ustressd"
          n1 = n1 + 1
        endif
        if(icode.eq.181.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "vstress"
          n1 = n1 + 1
        endif
        if(icode.eq.181.and.(ifield.eq.2)) then
          cvar1(n1) = "vstressd"
          n1 = n1 + 1
        endif
        if(icode.eq.260.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "pp"
          n1 = n1 + 1
        endif
        if(icode.eq.260.and.(ifield.eq.2)) then
          cvar1(n1) = "ppd "
          n1 = n1 + 1
        endif
         if(icode.eq.309.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar1(n1) = "eminp   "
          n1 = n1 + 1
        endif
        if(icode.eq.309.and.(ifield.eq.2)) then
          cvar1(n1) = "eminpd"
          n1 = n1 + 1
        endif
        if(icode.eq.130.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "t "
            n1 = n1 + 1
        endif
        if(icode.eq.130.and.(ifield.eq.2)) then
            cvar1(n1) = "td  "
            n1 = n1 + 1
        endif
        if(icode.eq.135.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "omega"
            n1 = n1 + 1
        endif
        if(icode.eq.135.and.(ifield.eq.2)) then
            cvar1(n1) = "omegad  "
            n1 = n1 + 1
        endif

        if(icode.eq.301.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "hforc "
            n1 = n1 + 1
        endif
        if(icode.eq.301.and.(ifield.eq.2)) then
            cvar1(n1) = "hforcd  "
            n1 = n1 + 1
        endif

        if(icode.eq.307.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "richarson"
            n1 = n1 + 1
        endif
        if(icode.eq.307.and.(ifield.eq.2)) then
            cvar1(n1) = "richarsond"
            n1 = n1 + 1
        endif

        if(icode.eq.305.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "cdragw"
            n1 = n1 + 1
        endif
        if(icode.eq.305.and.(ifield.eq.2)) then
            cvar1(n1) = "cdragwd"
            n1 = n1 + 1
        endif
        if(icode.eq.306.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "cdragv"
            n1 = n1 + 1
        endif
        if(icode.eq.306.and.(ifield.eq.2)) then
            cvar1(n1) = "cdragvd"
            n1 = n1 + 1
        endif
        if(icode.eq.996.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "dumt1"
            n1 = n1 + 1
        endif
        if(icode.eq.996.and.(ifield.eq.2)) then
            cvar1(n1) = "dumt1d"
            n1 = n1 + 1
        endif
        if(icode.eq.997.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar1(n1) = "dumt2"
            n1 = n1 + 1
        endif
        if(icode.eq.997.and.(ifield.eq.2)) then
            cvar1(n1) = "dumt2d"
            n1 = n1 + 1
        endif

       return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine count3(itel,icode,ifield,cvar3,n3)

c------------------------------------------------------------------------------
c**
c** this routine counts variable's number for three level variables.
c** 
c------------------------------------------------------------------------------
      implicit none

  
      integer nx,ny

      parameter(nx=32,ny=64)
      integer     itel,icode,ifield,n3
      real*4      x(nx,ny)
      character*8 cvar3(50)


      if(itel.eq.1) n3 = 1

        if(icode.eq.131.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "u   "
            n3 = n3 + 1
        endif
        if(icode.eq.131.and.(ifield.eq.2)) then
            cvar3(n3) = "ud   "
            n3 = n3 + 1
        endif
        if(icode.eq.132.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "v   "
            n3 = n3 + 1
        endif
        if(icode.eq.132.and.(ifield.eq.2)) then
            cvar3(n3) = "vd   "
            n3 = n3 + 1
        endif
        if(icode.eq.148.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "psi "
            n3 = n3 + 1
        endif
        if(icode.eq.148.and.(ifield.eq.2)) then
            cvar3(n3) = "psid "
            n3 = n3 + 1
        endif
        if(icode.eq.149.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "chi "
            n3 = n3 + 1
        endif
        if(icode.eq.149.and.(ifield.eq.2)) then
            cvar3(n3) = "chid "
            n3 = n3 + 1
        endif
        if(icode.eq.156.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "gh"
            n3 = n3 + 1
        endif
        if(icode.eq.156.and.(ifield.eq.2)) then
            cvar3(n3) = "ghd "
            n3 = n3 + 1
        endif
        if(icode.eq.302.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "vforc "
            n3 = n3 + 1
        endif
        if(icode.eq.302.and.(ifield.eq.2)) then
            cvar3(n3) = "vforcd"
            n3 = n3 + 1
        endif
        if(icode.eq.303.and.(ifield.eq.1.or.ifield.eq.0)) then
            cvar3(n3) = "ageu"
            n3 = n3 + 1
        endif
        if(icode.eq.303.and.(ifield.eq.2)) then
            cvar3(n3) = "ageud"
            n3 = n3 + 1
        endif
        if(icode.eq.304.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar3(n3) = "agev"
          n3 = n3 + 1
        endif
        if(icode.eq.304.and.(ifield.eq.2)) then
          cvar3(n3) = "agevd"
          n3 = n3 + 1
        endif
        if(icode.eq.310.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar3(n3) = "pv"
          n3 = n3 + 1
        endif
        if(icode.eq.310.and.(ifield.eq.2)) then
          cvar3(n3) = "pvd"
          n3 = n3 + 1
        endif
        if(icode.eq.998.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar3(n3) = "dumu1"
          n3 = n3 + 1
        endif
        if(icode.eq.998.and.(ifield.eq.2)) then
          cvar3(n3) = "dumu1d"
          n3 = n3 + 1
        endif
        if(icode.eq.999.and.(ifield.eq.1.or.ifield.eq.0)) then
          cvar3(n3) = "dumu2"
          n3 = n3 + 1
        endif
        if(icode.eq.999.and.(ifield.eq.2)) then
          cvar3(n3) = "dumu2d"
          n3 = n3 + 1
        endif

      return  
      end
c23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine wridatax(itel,icode,x,irec,indx)

      implicit none
  
      integer nx,ny

      parameter(nx=32,ny=64)
      
      integer i,j,k,indx
      integer icode,icode0,ifield,irec(3),itel
      real*4  x(nx,ny),y(nx,ny,10)
      real*4  dumy(nx,ny)

      if(itel.eq.1) then
        do i = 1, 3
          irec(i) = 1
        enddo
      endif
      do i = 1, nx
        do j= 1, ny
          dumy(i,j) = 0.0
        enddo
      enddo
  
      if (icode.eq.133.or.icode.eq.139.or.icode.eq.142.or.
     &  icode.eq.143.or.icode.eq.146.or.icode.eq.147.or.
     &  icode.eq.157.or.icode.eq.182.or.icode.eq.178.or.
     &  icode.eq.175.or.icode.eq.176.or.icode.eq.177.or.
     &  icode.eq.179.or.icode.eq.140.or.icode.eq.134.or.
     &  icode.eq.260.or.icode.eq.309.or.icode.eq.159.or.
     &  icode.eq.141.or.icode.eq.160.or.icode.eq.161.or.
     &  icode.eq.180.or.icode.eq.181.or.icode.eq.164.or.
     &  icode.eq.130.or.icode.eq.135.or.icode.eq.174.or.
     &  icode.eq.305.or.icode.eq.306.or.icode.eq.307.or.
     &  icode.eq.308.or.icode.eq.996.or.icode.eq.997) then
        write(110,rec=irec(1)) ((x(i,j),j=1,ny),i=1,nx)
        irec(1) = irec(1) + 1
      endif

      if (icode.eq.301.and.indx.eq.0) then
        write(110,rec=irec(1)) ((dumy(i,j),j=1,ny),i=1,nx)
        irec(1) = irec(1) + 1
        write(110,rec=irec(1)) ((x(i,j),j=1,ny),i=1,nx)
        irec(1) = irec(1) + 1
      endif
      if (icode.eq.301.and.indx.ne.0) then
        write(110,rec=irec(1)) ((x(i,j),j=1,ny),i=1,nx)
        irec(1) = irec(1) + 1
      endif
 
      if (icode.eq.131.or.icode.eq.132.or.icode.eq.148.or.
     &  icode.eq.302.or.icode.eq.303.or.icode.eq.304.or.
     &  icode.eq.149.or.icode.eq.156.or.icode.eq.310.or.
     &  icode.eq.998.or.icode.eq.999) then
        write(112,rec=irec(3)) ((x(i,j),j=1,ny),i=1,nx)
        irec(3) = irec(3) + 1
      endif 
      
      return
      end
c123456789012345678901234567890123456789012345678901234567890123456789012
      subroutine wridatay2(icode,y,irec,indy)
       
      implicit none
    
  
      integer nx,ny

      parameter(nx=32,ny=64)
      integer i,j,k
      integer icode, irec(3),indy
      real*4  y(nx,ny,10)
      real*4  dumy(nx,ny)

      do i = 1, nx
        do j= 1, ny
          dumy(i,j) = 0.0
        enddo
      enddo
        if (icode.eq.133.or.icode.eq.139.or.icode.eq.142.or.
     &      icode.eq.143.or.icode.eq.146.or.icode.eq.147.or.
     &      icode.eq.157.or.icode.eq.182.or.icode.eq.178.or.
     &      icode.eq.175.or.icode.eq.176.or.icode.eq.177.or.
     &      icode.eq.179.or.icode.eq.140.or.icode.eq.134.or.
     &      icode.eq.260.or.icode.eq.309.or.icode.eq.159.or.
     &      icode.eq.141.or.icode.eq.160.or.icode.eq.161.or.
     &      icode.eq.180.or.icode.eq.181.or.icode.eq.164.or.
     &      icode.eq.305.or.icode.eq.306.or.icode.eq.307.or.
     &      icode.eq.174.or.icode.eq.308) then
            write(110,rec=irec(1)) ((y(i,j,1),j=1,ny),i=1,nx)
            irec(1) = irec(1) + 1
        endif

        if (icode.eq.301) then
          write(110,rec=irec(1)) ((dumy(i,j),j=1,ny),i=1,nx)
          irec(1) = irec(1) + 1
          do k = 1,indy
            write(110,rec=irec(1)) ((y(i,j,k),j=1,ny),i=1,nx)
            irec(1) = irec(1) + 1
          enddo
        endif
        if (icode.eq.130.or.icode.eq.135.or.icode.eq.996.or.
     &      icode.eq.997) then
          do k = 1,indy
            write(110,rec=irec(1)) ((y(i,j,k),j=1,ny),i=1,nx)
            irec(1) = irec(1) + 1
          enddo
        endif 
 
        if (icode.eq.131.or.icode.eq.132.or.icode.eq.148.or.
     &      icode.eq.302.or.icode.eq.303.or.icode.eq.304.or.
     &      icode.eq.149.or.icode.eq.156.or.icode.eq.310.or.
     &      icode.eq.998.or.icode.eq.999) then
          do k = 1,indy
            write(112,rec=irec(3)) ((y(i,j,k),j=1,ny),i=1,nx)
            irec(3) = irec(3) + 1
          enddo
        endif 
  
      return
      end
