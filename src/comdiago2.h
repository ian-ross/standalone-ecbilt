c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comdiago2.h                                                   
c *** Contents: Common declarations for diagnostics of ocean
c-----------------------------------------------------------------------
      real*8 s1te(0:nl,0:nm,nz),s2te(0:nl,0:nm,nz),
     &       s1sa(0:nl,0:nm,nz),s2sa(0:nl,0:nm,nz),
     &       s1ro(0:nl,0:nm,nz),s2ro(0:nl,0:nm,nz),
     &       s1hfluxt(0:nl,0:nm,nz),s2hfluxt(0:nl,0:nm,nz),
     &       s1hfluxs(0:nl,0:nm,nz),s2hfluxs(0:nl,0:nm,nz),
     &       s1vfluxt(0:nl,0:nm,nz),s2vfluxt(0:nl,0:nm,nz),
     &       s1vfluxs(0:nl,0:nm,nz),s2vfluxs(0:nl,0:nm,nz),
     &       s1hdift(0:nl,0:nm,nz),s2hdift(0:nl,0:nm,nz),
     &       s1hdifs(0:nl,0:nm,nz),s2hdifs(0:nl,0:nm,nz),
     &       s1vdift(0:nl,0:nm,nz),s2vdift(0:nl,0:nm,nz),
     &       s1vdifs(0:nl,0:nm,nz),s2vdifs(0:nl,0:nm,nz),
     &       s1dcont(0:nl,0:nm,nz),s2dcont(0:nl,0:nm,nz),
     &       s1dcons(0:nl,0:nm,nz),s2dcons(0:nl,0:nm,nz),
     &       s1u(0:nl,0:nm,nz),s2u(0:nl,0:nm,nz),
     &       s1v(0:nl,0:nm,nz),s2v(0:nl,0:nm,nz),
     &       s1w(0:nl,0:nm,nz),s2w(0:nl,0:nm,nz),
     &       s1relaxt(0:nl,0:nm,nz),s2relaxt(0:nl,0:nm,nz),
     &       s1relaxs(0:nl,0:nm,nz),s2relaxs(0:nl,0:nm,nz)

      common /meanocx2/s1te,s2te, s1sa,s2sa,s1ro,s2ro,
     &       s1hfluxt,s2hfluxt, s1hfluxs,s2hfluxs, s1vfluxt,s2vfluxt,
     &       s1vfluxs,s2vfluxs, s1hdift,s2hdift, s1hdifs,s2hdifs,
     &       s1vdift,s2vdift, s1vdifs,s2vdifs, s1dcont,s2dcont,
     &       s1dcons,s2dcons, s1u,s2u, s1v,s2v, s1w,s2w,
     &       s1relaxt,s2relaxt,s1relaxs,s2relaxs



