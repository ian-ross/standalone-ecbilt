      open(iuo+46,file='namelist',status='old',form='formatted')

      open(iuo+48,file='inputdata/fractoc.dat',form='formatted')
      open(iuo+50,file='inputdata/darea.dat',form='formatted')

      open(iuo+1,file='inputdata/coef.dat',status='old',form='unformatted')
      open(iuo+7,file='inputdata/gauss.asc',status='old',form='formatted')
      open(iuo+16,file='inputdata/lwrref.dat',form='unformatted')
      open(iuo+17,file='inputdata/lwrcoef.dat',form='unformatted')
      open(iuo+18,file='inputdata/swrref.dat',form='unformatted')
      open(iuo+19,file='inputdata/swrcoef.dat',form='unformatted')

      open(iuo+4,file='inputdata/albsnow.dat')
      open(iuo+10,file='inputdata/labas.dat')
      open(unit=iuo+31,file='inputdata/land_albedo.dat')
      open(unit=iuo+32,file='inputdata/forfr.dat')

      open(unit=iuo+41,file='inputdata/sst_daily.dat',form='unformatted')
      open(unit=iuo+43,file='inputdata/seaice.dat')

      open(iuo+30,file='startdata/parlist'//fini,form='formatted')
      open(iuo+20,file='book'//fini,form='formatted')
      open(iuo+99,file='info'//fini,form='formatted')
      open(iuo+29,file='error'//fini,form='formatted')
      open(iuo+74,file='ipcc'//fini,form='formatted')
