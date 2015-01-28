












      subroutine inter_cfc(xjour)
c
      include 'comcouphelp.h'
c      include 'comcfchelp.h'
      include 'comemic.h'
      include 'comrunlabel.h'
      include 'const.com'
      include 'trace.com'
c
c Set up CFC-11 and CFC-12 atmospheric concentrations as a function
c of hemisphere and time (linearly interpolate yearly data).
c
      nn99=2
      if(nn99.eq.2) write(99,*)'pendant inter',iyear,xjour
      cfcyear=irunlabel+iyear
      x=cfcyear+xjour/yeaday
c
      if(x.lt.1931.5)then
c     CFC=0 before 1931.5
        atmcfc11(1) = 0.0
        atmcfc11(2) = 0.0
        atmcfc12(1) = 0.0
        atmcfc12(2) = 0.0
      else if(x.le.1997.5)then
c     linear interpolation between two dates
        x1=int(x-0.5)+0.5
        i1=int(x1-1930.5)
        do ihemcfc=1,2
          atmcfc11(ihemcfc) = cfc11(i1,ihemcfc)+
     &           (cfc11(i1+1,ihemcfc)-cfc11(i1,ihemcfc))*(x-x1)
          atmcfc12(ihemcfc) = cfc12(i1,ihemcfc)+
     &           (cfc12(i1+1,ihemcfc)-cfc12(i1,ihemcfc))*(x-x1)
        enddo
      else
c     CFC are kept constant to the value at 1997 after that date. 
        i1=68
        do ihemcfc=1,2
          atmcfc11(ihemcfc) = cfc11(i1,ihemcfc)
          atmcfc12(ihemcfc) = cfc12(i1,ihemcfc)
        enddo
      endif
c
c     write
c
      if (nn99.eq.2) then
        write(99,977)
        write(99,978) ja,xjour,x,
     &      atmcfc11(2),atmcfc12(2),atmcfc11(1),atmcfc12(1)
      endif
 977  format('inter ja xjour   Year:   CFC-11(N): CFC-12(N): CFC-11(S): 
     & CFC-12(S):')
 978  format(i5,1x,f4.0,f9.3,4f11.2)
      return
      end
