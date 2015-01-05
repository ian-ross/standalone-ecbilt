         PROGRAM lirelw

         integer irn(32,64,2)
         open(1,file='lwrref.dat',form='unformatted')
         open(2,file='irnout.dat',status='unknown')

         read(1)irn
         do k=1,2
         write(2,*)k,"--------------------------------"
         do i=1,32
           write(2,'(64I3)')(irn(i,j,k),j=1,64)
         enddo 
         enddo 
         end
