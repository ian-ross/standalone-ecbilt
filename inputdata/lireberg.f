        PROGRAM lireberg

        integer lrec,irec
        real*8  agg2(32,64), agg1(32,64)
        real*4  agg(32,64)
        lrec=32*64
c       open(10,file='berg.dat',form='unformatted')
        open(10,file='toporefism.dat',status='old')
        open(11,file='berg2.dat',form='unformatted',access='direct',
     *       recl=lrec)

        do i=1,32
         read(10,'(64F9.3)')(agg2(i,j),j=1,64)
        enddo

c       read(10)agg2
c       read(10)agg1
        do i=1,32
         do j=1,64
           agg(i,j)=agg2(i,j)
         enddo
        enddo

        irec=1
        write(11,rec=irec) ((agg(i,j),j=1,64),i=1,32)
        do i=1,32
         do j=1,64
           write(*,*)agg(i,j)
         enddo
        enddo
        close(10)  
        close(11)
        end
