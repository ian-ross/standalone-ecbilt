c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** this routine from numerical recipes determines an interval with
c *** bounds x1 and x2 that contains the root of the function func
c-----------------------------------------------------------------------
      SUBROUTINE zbrac(func, x1, x2, iter)
      IMPLICIT NONE

      INTEGER j, iter
      REAL*8 f1, f2, x1, x2, func
      EXTERNAL func
      INTEGER, PARAMETER :: ntry = 100
      REAL*8, PARAMETER :: factor = 1.6

      f1 = func(x1)
      f2 = func(x2)
      DO j = 1, ntry
         IF (f1 <= 0.0 .AND. f2 >= 0.0 .OR.
     $       f1 >= 0.0 .AND. f2 <= 0.0) THEN
            iter = j
            RETURN
         END IF
         IF (ABS(f1) < ABS(f2)) THEN
            x1 = x1 + factor * (x1 - x2)
            f1 = func(x1)
         ELSE
            x2 = x2 + factor * (x2 - x1)
            f2 = func(x2)
         END IF
      END DO
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012
      function zbrent(func,x1,x2,tol,iter)
c-----------------------------------------------------------------------
c *** this routine from numerical recipes determines the root of  the
c *** function func which is contained in the interval with bounds
c *** x1 and x2
c-----------------------------------------------------------------------

      implicit none
      include 'comunit.h'

      integer   iter,itmax
      real*8    eps
      parameter (itmax=100,eps=1e-8)
      real*8    func,x1,x2,tol,a,b,fa,fb,fc,c,d,e,tol1,xm,zbrent
      real*8    p,q,r,s
      external  func

      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      d=0
      e=0
      if ((fa.gt.0.and.fb.gt.0.).or.(fa.lt.0.and.fb.lt.0.)) then
         call error(13)
      endif
      c=b
      fc=fb
      do iter=1,itmax
        if ((fb.gt.0.and.fc.gt.0.).or.(fb.lt.0.and.fc.lt.0.)) then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if (abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*eps*abs(b)+0.5*tol
        xm=0.5*(c-b)
        if (abs(xm).le.tol1.or.fb.eq.0.) then
          zbrent=b
          if (iter.gt.8) write(iuo+29,*) 'zbrent ',iter,b
          return
        endif
        if (abs(e).ge.tol1.and.abs(fa).gt.abs(fb)) then
          s=fb/fa
          if (a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if (p.gt.0.) q=-q
          p=abs(p)
          if (2.*p.lt.min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if (abs(d).gt.tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
        if (iter.gt.8) write(iuo+29,*) 'zbrent ',iter,b
      enddo
      zbrent=b
      return
      end
