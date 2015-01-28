












      subroutine rahmflx(ns,nn99,scaldt)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Calcule (en nitrah ss.iteration) et integre le flux (thermique)
c  issu de la diffusion des anomalies de T. (cf. S.Rahmstorf).
c Rappel Explic. calcule et incorpore apres Diffus.Anom(<-avec ss.iter)
c--
c  modif : 25/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
C     include 'reper.com'
 
c--dummy variables :
      dimension scaldt(imax,jmax,kmax)
 
c--common local (variables a conserver d'un appel a l'autre) :
      common / rhmflx /
     &  xfrz, tfrz, tfrz1, dtfrz, xvois, xxfrz,
     &  cphix(imax,jmax), cphiy(imax,jmax), unsmvs(imax,jmax)
 
c--variables locales :
      dimension tanom(imax,jmax), tanom1(imax,jmax)
      dimension tanom0(imax,jmax), tmfrz(imax,jmax), buffer(nsmax+4)
c-
 
      k = ks2
      if (numit.eq.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Mise en place parametres & coeff. <= only the 1rst time step    |
c-----------------------------------------------------------------------
 
c--Lecture des parametres (additionel) dans le fichier "run.param" :
        open(19, file='run.param', status='old')
        nbline = 10 + 2*nsmax
        do 105 n=1,nbline
          read(19,*)
 105    continue
        read(19,*,err=900) (buffer(n),n=1,nsmax+4), nnloc,
     &             xfrz, tfrz, tfrz0, dtfrz, xvois
        close(19)
        if (nnloc.ne.nitrap) goto 900
        tfrz1 = tfrz0 + dtfrz
c---
c   xfrz = part relative (maximale) de (T*) remplacee par "tfrz"
c   tfrz = temp. de remplacement (a la place de T*)
c   tfrz0,tfrz1 = intervalle de transition (=> decide de remplacer T*)
c   xvois = poids des voisins (centre = 1) dans la moyenne lissee
c   par ex. : xfrz,tfrz,tfrz0,tfrz1,xvois = 1., -2., -2.5, -1.5,  3.
c---
 
        if (nn99.eq.2) then
         write(99,'(2A,I5,1PE14.6)') 'rahmflx : Vers=c(f) ',
     &          '; Rappes apres Ah_Rap ; nitrap,ahrap=', nitrap, ahrap
         write(99,'(A,2F10.6,3F10.5)') ' xfrz,xvois ,tfrz,tfrz0 & 1 :',
     &           xfrz, xvois, tfrz, tfrz0, tfrz1
        endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Definit les constantes une fois pour toutes :
      xxfrz = (xfrz / dtfrz) / (1.0 + xvois)
 
c- precaution : cyclicite (et autre) du tab. ``T.star'' :
      call raccord(scalr(1,1,k,ns), spvr, 1, 8)
 
c--Mise en place des Coeff. pour  Diffus.Anom.Temp. (<-S.Rahmstorf) :
c-  cphix&y <-deja.mult par deltaT/Dz [Rq. Dz etait deja retiree ds run.param]
      unsrap = one / DFLOAT(nitrap)
      ccdif = dts(k) * ahrap * unsdx * unsrap
      do 110 j=js1,js2
       do 110 i=is1(j),1+is2(j)
        ttm1ob = min(tms(i-1,j,k), (scalr(i-1,j,k,1)-spvr),
     &               tms(i,j,k),   (scalr(i,j,k,1)-spvr)   )
        cphix(i,j) = ccdif * smxy(i,j,1) * ttm1ob
 110  continue
      ccdif = dts(k) * ahrap * unsdy * unsrap
      do 130 j=js1,1+js2
       do 130 i=isf1(j),isf2(j)
        ttm2ob = min(tms(i,j-1,k), (scalr(i,j-1,k,1)-spvr),
     &               tms(i,j,k),   (scalr(i,j,k,1)-spvr)   )
        cphiy(i,j) = ccdif * cmxy(i,j,2) * ttm2ob
 130  continue
 
c- Pour le lissage : masque des 4 voisins :
      do 150 j=js1,js2
       do 150 i=is1(j),is2(j)
        sumtm = (tms(i-1,j,k)+tms(i+1,j,k)+tms(i,j-1,k)+tms(i,j+1,k))
        unsmvs(i,j) = (xvois / dtfrz) * tms(i,j,k) / max(epsil,sumtm)
 150  continue
c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Initialisation.                                                 |
c-----------------------------------------------------------------------
 
c- mise en place de "tanom" et du masque "tmfrz" (temp. below freezing point)
      do 210 j=1,jmax
       do 210 i=1,imax
        tanom(i,j) = scalr(i,j,k,ns) - scal(i,j,k,ns)
        tmfrz(i,j) = tms(i,j,k) *
     &           max(zero, min(dtfrz, tfrz1-scal(i,j,k,ns)) )
 210  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Modifie temperature "T*" where T < tfrz ; -> remplacee par (1-X).T* + X.tfrz
c  avec : if T<tfrz, X = xfrz (1 + xvois.S_4(if)/S_4tm)/(1+xvois) ; 0 si non
 
      do 250 j=js1,js2
       do 250 i=is1(j),is2(j)
        ttmfrz = xxfrz * tmfrz(i,j) * ( 1.0 + unsmvs(i,j) *
     &       (tmfrz(i-1,j)+tmfrz(i+1,j)+tmfrz(i,j-1)+tmfrz(i,j+1)) )
        tanom(i,j)  = tanom(i,j)
     &              + ttmfrz * max(zero, tfrz - scalr(i,j,k,ns))
        tanom0(i,j) = tanom(i,j)
 250  continue
c-----
 
C     if (nn99.eq.2 .and. numit.eq.nlast) write(99,'(A,I9,I5,3F12.6)')
C    &   ' numit, nit, tanom =', numit, 0, tanom(icheck,jcheck),
C    &       scalr(icheck,jcheck,k,ns), scal( icheck,jcheck,k,ns)
 
      do 400 nit=1,nitrap
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Calcule et somme (ds tanom) les flux par iteration sur "nit".   |
c-----------------------------------------------------------------------
 
C     call raccord(tanom, zero, 1, 8)
c- "inlining explicite" :
      do 305 j=jcl1,jcl2
        tanom(ims1-1,j) = tanom(ims2,j)
        tanom(ims2+1,j) = tanom(ims1,j)
 305  continue
        tanom(iberpm,jberp) = tanom(ibera, jberam)
        tanom(iberp, jberp) = tanom(iberam,jberam)
        tanom(iberam,jbera) = tanom(iberp, jberpm)
        tanom(ibera, jbera) = tanom(iberpm,jberpm)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Forcage thermique : Diffus.Anom.Temp. (<-S.Rahmstorf) :
c-  cphix & cphiy <-deja.mult par deltaT/Dz et / par nitrap
c----
c- calcul du flux et l'ajoute a "tanom" dans "tanom1":
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,flxrah)
      do 320 j=js1,js2
       do 310 i=is1(j),is2(j)
        flxrah = unsdx*(cphix(i,j)*( tanom(i-1,j) - tanom(i,j) )
     &                 -cphix(i+1,j)*( tanom(i,j) - tanom(i+1,j) ))
     &         + unsdy*(cphiy(i,j)*( tanom(i,j-1) - tanom(i,j) )
     &                 -cphiy(i,j+1)*( tanom(i,j) - tanom(i,j+1) ))
        tanom1(i,j) = tanom(i,j) + smxy(i,j,0) * flxrah
 310   continue
 320  continue
 
c- Mise a jour de "tanom" :
      do 330 j=js1,js2
       do 330 i=is1(j),is2(j)
        tanom(i,j) = tanom1(i,j)
 330  continue
 
C     if (nn99.eq.2 .and. numit.eq.nlast) write(99,'(A,I9,I5,F12.6)')
C    &   ' numit, nit, tanom =', numit, nit, tanom(icheck,jcheck)
c----
 400  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Mise a jour des tableaux fss & scaldt.                          |
c-----------------------------------------------------------------------
 
c  Ajoute la nouvelle evolution de T.surf calculee par nitrap iterations :
      do 450 j=js1,js2
       do 450 i=is1(j),is2(j)
        fflx = (tanom(i,j) - tanom0(i,j)) - rappes(i,j,ns) *
     &        ( tanom(i,j) + (scal(i,j,k,ns) - scalr(i,j,k,ns)) )
        fss(i,j,ns) = fss(i,j,ns) + fflx
        scaldt(i,j,k) = scaldt(i,j,k) - fflx
 450  continue
 
      return
 
 900  continue
c- cas d'erreur :
      write(iuo+66,'(2A,I4)') 'ARRET, rahmflx :',
     &     ' Error in reading "run.param", line', nbline+1
      stop
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine rahmflx -
      end
