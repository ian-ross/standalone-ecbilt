












      subroutine raccord(var,spv,krac,ltyp)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  Appele par "class", "?",
c  Mise en place des raccords cycliques et autres raccords
c   pour les "krac" niveaux du tableau "var" .
c tient compte de spv si ltyp > 36 et changement de signe.
c  modif : 25/05/90
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c- dummy variables :
      dimension var(imax,jmax,*)
 
c- variables locales equivalentes :
      dimension ijdl(0:1)
      equivalence (ijdl(0) , ijsdl )
      equivalence (ijdl(1) , ijudl )
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) raccord cyclique .                                              |
c-----------------------------------------------------------------------
 
      if (ltest.ge.1) then
c-
        do 110 k=1,krac
         do 110 j=jcl1,jcl2
          var(ims1-1,j,k) = var(ims2,j,k)
          var(ims2+1,j,k) = var(ims1,j,k)
 110    continue
 
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) raccord des 2 grilles .                                         |
c-----------------------------------------------------------------------
 
      if (ltest.eq.2) then
c-
        llm = mod(ltyp,4) / 3
        llc = ltyp / 12
c-
        if (llc.eq.0) then
          do 210 k=1,krac
           do 210 j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,k) = var(ii,jeq-1,k)
            var(ii,jeq,k) = var(ims1,j,k)
 210      continue
c-
        elseif (llc.eq.1) then
          do 220 k=1,krac
           kk = k + krac
           do 220 j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,k) = var(ii,jeq-1,kk)
            var(ii,jeq,kk) = var(ims1,j,k)
 220      continue
c-
        elseif (llc.eq.2) then
          do 230 k=1,krac
           kk = k + krac
           do 230 j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,kk) = -var(ii,jeq-1,k)
            var(ii,jeq,k) = -var(ims1,j,kk)
 230      continue
c-
        else
          do 240 k=1,krac
           kk = k + krac
           do 240 j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,kk) = -var(ii,jeq-1,k)
            var(ii,jeq,k) = -var(ims1,j,kk)
            if (var(ii,jeq-1,k).eq.spv) var(ims1-1,j,kk) = spv
            if (var(ims1,j,kk) .eq.spv) var(ii,jeq,k) = spv
 240      continue
 
        endif
 
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) raccord pour Bering .                                           |
c-----------------------------------------------------------------------
 
C     if (ltest.eq.3 .and. iberp.ne.ibera) then
      if (ltest.eq.3) then
c-
        ltp = mod(ltyp,4)
c-
c-----position 0 :
        if (ltp.eq.0) then
          do 300 k=1,krac
            var(iberpm,jberp,k) = var(ibera, jberam,k)
            var(iberp, jberp,k) = var(iberam,jberam,k)
            var(iberam,jbera,k) = var(iberp, jberpm,k)
            var(ibera, jbera,k) = var(iberpm,jberpm,k)
 300      continue
c-----position 1 :
        elseif (ltp.eq.1) then
          if (ltyp.lt.12) then
            do 310 k=1,krac
              var(iberp,jberp,k) = var(ibera,jberam,k)
              var(ibera,jbera,k) = var(iberp,jberpm,k)
 310        continue
          elseif (ltyp.lt.36) then
            do 312 k=1,krac
              var(iberp,jberp,k) = -var(ibera,jberam,k)
              var(ibera,jbera,k) = -var(iberp,jberpm,k)
 312        continue
          else
            do 315 k=1,krac
              if (var(iberp,jberpm,k).eq.spv) then
                var(ibera,jbera,k) = spv
              else
                var(ibera,jbera,k) = -var(iberp,jberpm,k)
              endif
              if (var(ibera,jberam,k).eq.spv) then
                var(iberp,jberp,k) = spv
              else
                var(iberp,jberp,k) = -var(ibera,jberam,k)
              endif
 315        continue
          endif
c-----position 2 :
        elseif (ltp.eq.2) then
          if (ltyp.lt.12) then
            do 320 k=1,krac
              var(ibera ,jbera,k) = var(iberpm,jberp,k)
              var(iberam,jbera,k) = var(iberp ,jberp,k)
 320        continue
          elseif (ltyp.lt.36) then
            do 322 k=1,krac
              var(ibera ,jbera,k) = -var(iberpm,jberp,k)
              var(iberam,jbera,k) = -var(iberp ,jberp,k)
 322        continue
          else
            do 325 k=1,krac
              if (var(iberpm,jberp,k).eq.spv) then
                var(ibera, jbera,k) = spv
              else
                var(ibera ,jbera,k) = -var(iberpm,jberp,k)
              endif
              if (var(iberp,jberp,k).eq.spv) then
                var(iberam,jbera,k) = spv
              else
                var(iberam,jbera,k) = -var(iberp ,jberp,k)
              endif
 325        continue
          endif
c-----position 3 :
        elseif (ltyp.lt.12) then
          do 330 k=1,krac
            var(ibera,jbera,k) = var(iberp,jberp,k)
 330      continue
        elseif (ltyp.lt.36) then
          do 332 k=1,krac
            var(ibera,jbera,k) = -var(iberp,jberp,k)
 332      continue
        else
          do 335 k=1,krac
            if (var(iberp,jberp,k).eq.spv) then
              var(ibera,jbera,k) = spv
            else
              var(ibera,jbera,k) = -var(iberp,jberp,k)
            endif
 335      continue
        endif
 
      endif
 
      return
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine raccord -
      end
