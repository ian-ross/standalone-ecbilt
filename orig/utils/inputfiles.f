      program inputfiles
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comocfix.h'


      real*4       sstday4(nlat,nlon)
      real*8       agg1(nlat,nlon), agg2(nlat,nlon)
      integer      id,i,j

      open(unit=1,file='inputdata/atmos/coef.dat',
     &     status='unknown',form='unformatted')
      open(unit=3,file='inputdata/atmos/berg.dat',
     &     status='unknown',form='unformatted')
      open(unit=4,file='inputdata/ocfix/sst_daily.dat',
     &              status='unknown',form='unformatted')

      open(unit=11,file='inputdata/atmos/coef.asc',
     &      status='unknown',form='formatted')
      open(unit=13,file='inputdata/atmos/berg.asc',
     &      status='unknown',form='formatted')
      open(unit=14,file='inputdata/ocfix/sst_daily.asc',
     &      status='unknown',form='formatted')


      read (11,*) nshm, ll
      read (11,*) pp
      read (11,*) pd
      read (11,*) pw
      read (13,*) agg1
      read (13,*) agg2

      do id=1,360
        read(14,*) (( sstday4(i,j),j=1,nlon),i=1,nlat)
        write(4) (( sstday4(i,j),j=1,nlon),i=1,nlat)
      enddo

      write (1) nshm, ll
      write (1) pp
      write (1) pd
      write (1) pw
      write (3) agg1
      write (3) agg2

      end
