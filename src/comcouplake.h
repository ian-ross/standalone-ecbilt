c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comcouplake.h                                                    
c *** Contents: Common declarations for coupling lakes to ocean ECbilt
c-----------------------------------------------------------------------


      logical inter(0:nll,0:nml)
      integer iwater(-1:nll+1,-1:nml+1)
      real*8  relalt(-1:nll+1,-1:nml+1,nzl)
      real*8  relals(-1:nll+1,-1:nml+1,nzl)
      real*8  reloct(-1:nll+1,-1:nml+1,nzl)
      real*8  relocs(-1:nll+1,-1:nml+1,nzl)


      common /lcouplaoc/ inter
      common /icouplaoc/ iwater
      common /rcouplaoc/ relalt,relals,reloct,relocs

	
