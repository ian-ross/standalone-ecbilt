c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of ocean files:
c *** units 88-99 are reserved for initial states
c *** units below 40 are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

c *** removal of gradsdata files (requirement for power indigo)

      call system('rm -f outputdata/ocean/*'//fini//'*grads')

c *** file 53 only in case of ocean only run

c      open(53,file='outputdata/ocean/ocfor'//fini//'.dat',
c     *      form='unformatted')

      open(54,file='outputdata/ocean/octra'//fini//'.grads',
     *  form='unformatted',access='direct',recl=64*32)

      open(55,file='outputdata/ocean/ocvel'//fini//'.grads',
     *  form='unformatted',access='direct',recl=64*32)

      open(56,file='outputdata/ocean/occomp'//fini//'.grads',
     *  form='unformatted', access='direct',recl=64*32)

      open(57,file='outputdata/ocean/ocro'//fini//'.grads',
     *  form='unformatted',access='direct',recl=64*32)

      open(58,file='outputdata/ocean/oczmtrnsp'//fini//'.grads',
     *  form='unformatted',access='direct',recl=32*12)

      open(59,file='outputdata/ocean/oczmtra'//fini//'.grads',
     *  form='unformatted',access='direct',recl=32*12)

      open(60,file='outputdata/ocean/octottrnsp'//fini//'.grads',
     *  form='unformatted',access='direct',recl=32)

      open(61,file='outputdata/ocean/coup'//fini//'.grads',
     *  form='unformatted',access='direct',recl=64*32)

      open(62,file='outputdata/ocean/lake'//fini//'.grads',
     *  form='unformatted',access='direct',recl=64*32)

      open(65,file='outputdata/ocean/ocmeanov'//fini//'.grads',
     *  form='unformatted',access='direct',recl=4)

      open(71,file='outputdata/ocean/ocmfx'//fini//'.grads',
     *  form='unformatted',access='direct',recl=9)

      open(72,file='outputdata/ocean/lamfx'//fini//'.grads',
     *  form='unformatted',access='direct',recl=8)

      open(73,file='outputdata/ocean/totfx'//fini//'.grads',
     *  form='unformatted',access='direct',recl=12)

      open(76,file='outputdata/ocean/ocmeantra'//fini//'.grads',
     *  form='unformatted',access='direct',recl=12*9)

      open(77,file='outputdata/ocean/lameantra'//fini//'.grads',
     *  form='unformatted',access='direct',recl=6*8)

      open(78,file='outputdata/ocean/coupsn'//fini//'.grads',
     *  form='unformatted',access='direct',recl=64*32)

      open(79,file='outputdata/ocean/ocinst'//fini//'.dat',
     *  form='unformatted')

      open(80,file='outputdata/ocean/ocmm'//fini//'.dat',
     *  form='unformatted')

      open(81,file='outputdata/ocean/ocsm'//fini//'.dat',
     *  form='unformatted')
