c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of files containing initial state of ecbilt: 
c-----------------------------------------------------------------------
      if (irunlabel.gt.0) then
      open(92,file='startdata/inice'//fini//'.dat',form='unformatted')
      open(93,file='startdata/inocean'//fini//'.dat',form='unformatted')
      endif
