c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comocouthelp.h                                                    
c *** Contents: Common declarations for ocean output routines
c-----------------------------------------------------------------------

      real*8  fi(0:nmh),fiv(0:nmh+1)
      real*8  cofi(0:nmh),cofiv(0:nmh+1),dx(0:nmh),dy

      common /posit / fi,fiv,dx,dy,cofi,cofiv
