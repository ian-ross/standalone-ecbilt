c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of atmosphere output files: 
c *** units 90-99 are reserved for initial states
c *** units 1-20 and above 40  are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

 
c *** open data output files

      open(20,file='outputdata/atmos/book'//fini,form='formatted') 

      if (irunatm.eq.1) then 
        open(21, file='outputdata/atmos/atminst'//fini,
     *           form='unformatted')
        open(22, file='outputdata/atmos/atmmmwp'//fini,
     *           form='unformatted')
        open(23, file='outputdata/atmos/atmmmyl'//fini,
     *           form='unformatted')
        open(24, file='outputdata/atmos/atmsmwp'//fini,
     *           form='unformatted')
        open(25, file='outputdata/atmos/atmsmyl'//fini,
     *           form='unformatted')
        open(26,file='outputdata/atmos/result'//fini,
     &                 form='unformatted')
      endif

c *** open error file

      open(100,file='info'//fini,form='formatted')
 
c *** open error file

      open(29,file='error'//fini,form='formatted')
 
c *** open parameter output file

      open(30,file = 'startdata/parlist'//fini,form='formatted')

