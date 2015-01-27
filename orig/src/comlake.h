c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comlake.h                                                    
c *** Contents: Common declarations lakes of ECbilt
c-----------------------------------------------------------------------

      integer nml,nll,nzl      

      parameter (nml=31,nll=63,nzl=6)

      real*8  rkapvl
      real*8  hmsl,hml(0:nzl)
      real*8  rkaphla(9),rkaphin(9)

      integer ilaken(-1:nll+1,-1:nml+1),ilake(-1:nll+1,-1:nml+1) 
      integer iwater(-1:nll+1,-1:nml+1) 

      logical lake(0:nll,0:nml),water(0:nll,0:nml)

      common /ilakem/ lake,water
      common /imaskla/ ilaken,ilake,iwater
      common /hmaskla/ hmsl,hml
      common /paral/ rkapvl,rkaphla,rkaphin
