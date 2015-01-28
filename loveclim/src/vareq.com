












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "vareq.com" : incorpore par instruction 'include' dans les routines :
c     CLASS, GRADP, STATIS, ncdfout, moyen, local, binout, foroutp,
c      lisstab, (unigl, option)
c  (commun a toutes les routines de sortie (output) de resultats)
c - inclus apres "type.com", "para.com", "bloc.com", complete "varno.com".
c------------------
c Les variables sont rangees par niveau, les uns a la suite des autres,
c   dans un seul tableau 3D : vrl
c   et chaque couche d un tableau 3D est reperee par un indice general (k).
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 30/06/98
 
c-----
c--common regroupant les variables necessaires pour l'ecriture sur fichier :
 
      common / rtowrite /
     & spv(0:nvmax), cmult(0:nvmax), cadit(0:nvmax), cliss(0:nvmax)
 
      common / itowrite /
     & nfrc, nabs(0:nvmax), irn(imax,8), jrn(jmax,8),
     & irl1(0:3), irl2(0:3),jrl1(0:3), jrl2(0:3)
 
      character*10 fmt1
      character*13 titex1, titex2
      character*14 titex3
      character*40 titexp
 
      common / ctowrite /
     & fmt1(-2:nvmax),
     & titex1, titex2, titex3, titexp
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--variables equivalentes :
      dimension t(imax,jmax,kmax), s(imax,jmax,kmax)
      equivalence (scal(1,1,1,1), t(1,1,1)), (scal(1,1,1,2), s(1,1,1))
 
      dimension psh(imax,jmax), hmajc(imax,jmax), egajc(imax,jmax)
      dimension vrl(imax,jmax,krlmin:krlmax)
      dimension vaflux(imax,jmax,kmax), vdflux(imax,jmax,kmax)
      dimension alphxo(imax,jmax,kmax), alphyo(imax,jmax,kmax)
      dimension haterm(imax,jmax,kmax), hdterm(imax,jmax,kmax)
C     dimension vaterm(imax,jmax,kmax), vdterm(imax,jmax,kmax)
      dimension uslpfx(imax,jmax,-1:nsmax), vslpfx(imax,jmax,-1:nsmax)
 
      equivalence (vrl(1,1,krlfw) , fss(1,1,0)     )
      equivalence (vrl(1,1,krlps) , psh(1,1)       )
      equivalence (vrl(1,1,krlet) , eta(1,1)       )
      equivalence (vrl(1,1,krlub) , ub(1,1)        )
      equivalence (vrl(1,1,krlvb) , vb(1,1)        )
      equivalence (vrl(1,1,krlu ) , u(1,1,1)       )
      equivalence (vrl(1,1,krlv ) , v(1,1,1)       )
      equivalence (vrl(1,1,krlt ) , t(1,1,1)       )
      equivalence (vrl(1,1,krls ) , s(1,1,1)       )
      equivalence (vrl(1,1,krlb ) , b(1,1,1)       )
      equivalence (vrl(1,1,krln2) , bvf(1,1,1)     )
      equivalence (vrl(1,1,krlas) , avsdz(1,1,1)   )
      equivalence (vrl(1,1,krlau) , avudz(1,1,1)   )
      equivalence (vrl(1,1,krlw ) , w(1,1,1)       )
      equivalence (vrl(1,1,krltke), q2turb(1,1,1)  )
      equivalence (vrl(1,1,krlajc), fqajc(1,1,1)   )
      equivalence (vrl(1,1,krlvaf), vaflux(1,1,1)  )
      equivalence (vrl(1,1,krlvdf), vdflux(1,1,1)  )
      equivalence (vrl(1,1,krlaxt), alphxo(1,1,1)  )
      equivalence (vrl(1,1,krlayt), alphyo(1,1,1)  )
      equivalence (vrl(1,1,krlhat), haterm(1,1,1)  )
      equivalence (vrl(1,1,krlhdt), hdterm(1,1,1)  )
C     equivalence (vrl(1,1,krlvat), vaterm(1,1,1)  )
C     equivalence (vrl(1,1,krlvdt), vdterm(1,1,1)  )
      equivalence (vrl(1,1,krlhac), hmajc(1,1)     )
      equivalence (vrl(1,1,krlajc), egajc(1,1)     )
      equivalence (vrl(1,1,krlusl), uslpfx(1,1,-1) )
      equivalence (vrl(1,1,krlvsl), vslpfx(1,1,-1) )
 
c--fin du fichier "vareq.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
