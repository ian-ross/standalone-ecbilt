c23456789012345678901234567890123456789012345678901234567890123456789012
      function newton(func,x,iter)
      implicit none

      integer   j,itmax,iter
      real*8    tol,dx
      parameter (tol=0.1,itmax=100,dx=0.01)
      real*8    f,df,x,func,dfdx,delx,newton
      external  func

      iter=0
10    continue
        iter=iter+1
        f=func(x)
        df=func(x+dx)
        dfdx=(df-f)/dx
        delx=f/dfdx

15      if (abs(delx).gt.10.) then
          write(100,*) 'converg. problem in surfacetemp of icemodel'
          write(100,*) iter,f,dfdx,delx
          write(100,*) 'icemodel',iter,x
          delx=0.5*delx
          write(100,*) iter,f,dfdx,delx
          goto 15
        endif
        x=x-delx
        if (iter.lt.2) then
          goto 10
        else
          if (iter.gt.5) then
            write(100,*) 'icemodel',iter,x
          endif
          if (iter.gt.100) then
            call error(9)
          else
	    if (abs(delx).lt.tol) then
	      goto 20
            else
              goto 10
            endif
          endif
        endif
20    continue
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine zbrac(func,x1,x2,iter)
c-----------------------------------------------------------------------
c *** this routine from numerical recipes determines an interval with
c *** bounds x1 and x2 that contains the root of the function func
c-----------------------------------------------------------------------

      implicit none

      integer   j,ntry,iter
      real*8    factor
      parameter (factor=1.6,ntry=100)
      real*8    f1,f2,x1,x2,func
      external  func

      f1=func(x1)
      f2=func(x2)
      do j=1,ntry
        if ((f1.le.0..and.f2.ge.0.).or.(f1.ge.0..and.f2.le.0.)) then
          iter=j
          return
        endif
        if (abs(f1).lt.abs(f2)) then
          x1=x1+factor*(x1-x2)
          f1=func(x1)
        else
          x2=x2+factor*(x2-x1)
          f2=func(x2)
        endif
      enddo
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function zbrent(func,x1,x2,tol,iter)
c-----------------------------------------------------------------------
c *** this routine from numerical recipes determines the root of  the
c *** function func which is contained in the interval with bounds 
c *** x1 and x2
c-----------------------------------------------------------------------

      implicit none

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
          if (iter.gt.8) write(100,*) 'zbrent ',iter,b
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
        if (iter.gt.8) write(100,*) 'zbrent ',iter,b
      enddo
      zbrent=b
      return
      end
