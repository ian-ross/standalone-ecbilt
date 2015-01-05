c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of atmosphere output files: 
c *** units 90-99 are reserved for initial states
c *** units 1-20 and above 40  are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

 
c *** open data output files

      open(iuo+20,file='book'//fini,form='formatted') 

      open(iuo+26, file='outputdata/atmos/ocbasin'//fini,
     *           form='formatted')
      open(iuo+27, file='outputdata/atmos/ocheattr'//fini,
     *           form='unformatted')

c *** open error file

      open(iuo+99,file='info'//fini,form='formatted')
 
c *** open error file

      open(iuo+29,file='error'//fini,form='formatted')

      open(iuo+74,file='ipcc'//fini,form='formatted')

 

