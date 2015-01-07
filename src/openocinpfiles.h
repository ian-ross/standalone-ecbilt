c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of ocean files: 
c *** units 90-99 are reserved for initial states
c *** units below 40 are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

      open(unit=40,file='inputdata/ocean/mask.dat')
      open(unit=41,file='inputdata/ocean/ocbas.dat')
      open(unit=42,file='inputdata/ocean/lakemask.dat')
      
