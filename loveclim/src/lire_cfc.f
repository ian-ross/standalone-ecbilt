












      subroutine lire_cfc
c
      include 'trace.com'
c
c  Data is given for 66 years, and for the SH (1) and NH (2).
c  Note that year 1989 data is my own extrapolation to ensure last
c    few months of 1988 are forced properly.  STOP PRESS: Now have CFC
c    data going until 1995 based on Elkins et al., 1993 Nature Vol 364
c    pp 780-783.
c  Modified (MFL) June 2011
c  data given for 68 years (1931.5 to 1998.0) and for the SH (1) and
c  NH (2). data from OCMIP2
c  http://ocmip5.ipsl.jussieu.fr/OCMIP/phase2/simulations/CFC/HOWTO-CFC-1.html
c  file cfc1112.atm
c  
c=======================================================================
c    Read in atmospheric CFC histories
c=======================================================================
       nn99=2
       if (nn99.eq.2) write(99,977)
       open(91,status='old',file='cfc1112.atm')
       do i=1,6
         read(91,*)
       enddo
       do 7 nnyear=1931,1998
         iiyear=nnyear-1930
         read(91,979) cfc11(iiyear,2),cfc12(iiyear,2),
     &                cfc11(iiyear,1),cfc12(iiyear,1)
         if (nn99.eq.2) write(99,978) nnyear,
     &                cfc11(iiyear,2),cfc12(iiyear,2),
     &                cfc11(iiyear,1),cfc12(iiyear,1)
7      continue
       close(unit=91)
       if (nn99.eq.2) write(99,*)
977    format('Year:   CFC-11(N): CFC-12(N): CFC-11(S): CFC-12(S):')
978    format(i5, 4f11.2)
979    format(9x,2f8.2,3x,2f8.2)
c
       return
       end
