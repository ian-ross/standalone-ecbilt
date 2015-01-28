












c
c     ***************************
      function zmaxt(ndeb,nfin,a)
c     ***************************      
c
c		CALCULATE MAXIMUM VALUE OF THE COMPONENTS OF A VECTOR.
      include 'type.com'
c
      dimension a(*)
c
      zmaxt = a(ndeb)
      do 10 ji=ndeb,nfin
        zmaxt = max(zmaxt,a(ji))
10    continue
c
      return
      end
