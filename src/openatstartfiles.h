c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of files containing initial state of ecbilt: 
c-----------------------------------------------------------------------
      if (irunlabel.gt.0) then
      open(90,file='startdata/inatmos'//fini//'.dat',form='unformatted')
      open(91,file='startdata/inland'//fini//'.dat',form='unformatted')

      endif
