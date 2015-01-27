      program convert_sst_ascii_to_dat
      implicit none

      include 'comatm.h'

      real*4  sstday4(nlat, nlon)
      integer id, i, j

      open (unit=4,  file='sst_daily.dat', form='unformatted')
      open (unit=14, file='sst_daily.asc', form='formatted')
      do id = 1, 360
        read(14, *) ((sstday4(i, j), j = 1, nlon), i = 1, nlat)
        write(4)    ((sstday4(i, j), j = 1, nlon), i = 1, nlat)
      end do

      end
