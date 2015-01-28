












      subroutine slopez(cslpmx,kslpdw,nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  appelee par "scale" ;
c Detecte les cas de "Down-Sloping" ; Consigne Volume & Niveaux Concernes.
c  modif : 22/08/97
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'comunit.h'
 
c--variables locales equivalentes :
      dimension kfshk(ixjmax)
      dimension cmx1hk(ixjmax), cmx2hk(ixjmax)
      dimension cmy1hk(ixjmax), cmy2hk(ixjmax)
      dimension uhk(ixjmax,kmax), vhk(ixjmax,kmax), bhk(ixjmax,kmax)
      equivalence ( kfshk(1) , kfs(1,1) )
      equivalence ( cmx1hk(1) , cmx(1,1,1) )
      equivalence ( cmx2hk(1) , cmx(1,1,2) )
      equivalence ( cmy1hk(1) , cmy(1,1,1) )
      equivalence ( cmy2hk(1) , cmy(1,1,2) )
      equivalence ( uhk(1,1) , u(1,1,1) )
      equivalence ( vhk(1,1) , v(1,1,1) )
      equivalence ( bhk(1,1) , b(1,1,1) )
 
c--dummy variables :
      dimension kslpdw(nlpmax) 
CKO: , nslp(nlpmax)
CKO   dimension uvcslp(nlpmax), alpslp(nlpmax), cslpmx(kmax)
      dimension cslpmx(kmax)
 
c--variables locales :
      dimension nnkslp(0:kmax,4), uuavr(0:kmax), vvavr(0:kmax)
 
C     if (numit.eq.nstart) then
C       write(iuo+66,*) 'slopez : Uslp = xslop * dz * Db'
C     endif
 
      if (xslop.lt.epsil) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Detecte les cas de DownSloping pour Decentrer (upwing) .        |
c-----------------------------------------------------------------------
 
      nldw = 0
c--Direction X , Dowsloping => Decentrement :
      do 240 nl=1,nxslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        iij = ij + max(0,l)
        u1 = DFLOAT(-l) * (uhk(iij,k) + uhk(iij+imax,k))
        if (u1.gt.zero .and. bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
          kslpdw(nldw) = k
          alpslp(nldw) = 0.25 * cmy1hk(iij) * u1
        endif
 240  continue
      nnx = nldw
 
c--Direction Y , Dowsloping => Decentrement :
      do 260 nl=nxslp+1,nxyslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        ijj = ij + max(0,l)
        v2 = vhk(ijj,k) + vhk(ijj+1,k)
        if (v2*DFLOAT(l).lt.zero .and. bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
          kslpdw(nldw) = k
          alpslp(nldw) = 0.25 * cmx2hk(ijj) * abs(v2)
        endif
 260  continue
      nslpdw = nldw
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Detecte cas de DownSloping pour Decentrer et Permuter (xslop)   |
c-----------------------------------------------------------------------
 
      nldw = 0
c--Direction X , Dowsloping => permutation (xslop) :
      do 340 nl=1,nxslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        iij = ij + max(0,l)
        if ( bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
c--Recheche du 1er Niveau kk / rho_pot(i,j,kk) < rho_pot(i+l,j,k)
          do 330 kk=k-1,kfshk(ij),-1
            if (bhk(ij+l,kk).le.bhk(ij,kk) ) then
              kslpdw(nldw) = kk + 1
              goto 335
            endif
 330      continue
          kslpdw(nldw) = kfshk(ij)
 335      continue
c- calcul de la pente :
C         slopx = ((zw(k) - zw(kfshk(ij)) * smx1hk(iij) * unsdx
c- calcul de la vitesse :
C         uuslp = xslop * (bhk(ij+l,k) - bhk(ij,k))
          uuslp = xslop * dz(k) * (bhk(ij+l,k) - bhk(ij,k))
          uvcslp(nldw) = cmy1hk(iij) *
     &                   min( uuslp, cmx1hk(iij) * cslpmx(k) )
        endif
 340  continue
      nnx = nldw
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--Direction Y , Dowsloping => permutation (xslop) :
      do 360 nl=nxslp+1,nxyslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        ijj = ij + max(0,l)
        if ( bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
c--Recheche du 1er Niveau kk / rho_pot(i,j,kk) < rho_pot(i,j+l,k)
          do 350 kk=k-1,kfshk(ij),-1
            if (bhk(ij+l,kk).le.bhk(ij,kk) ) then
              kslpdw(nldw) = kk + 1
              goto 355
            endif
 350      continue
          kslpdw(nldw) = kfshk(ij)
 355      continue
c- calcul de la pente :
C         slopy = ((zw(k) - zw(kfshk(ij)) * smy1hk(iij) * unsdy
c- calcul de la vitesse :
C         vvslp = xslop * (bhk(ij+l,k) - bhk(ij,k))
          vvslp = xslop * dz(k) * (bhk(ij+l,k) - bhk(ij,k))
          uvcslp(nldw) = cmx2hk(iij) *
     &                   min( vvslp, cmy2hk(iij) * cslpmx(k) )
        endif
 360  continue
      nslpdw = nldw
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Verification a la 1ere itteration :                             |
c-----------------------------------------------------------------------
 
      if (numit.eq.nstart .and. nxyslp.ge.1) then
        nny = nslpdw - nnx
        write(iuo+66,'(2A,F8.3,4I5)') ' slopez(U=xslop*dz*Db) : xslop,',
     &   ' nXslp,nYslp,Nx,Ny=', xslop, nxslp,nxyslp-nxslp, nnx,nny
 
        if (nn99.eq.2) then
c--Evaluation du Nb de "boites" impliquees : (potentiel & effectif)
        do 400 k=0,kmax
         uuavr(k) = 0.
         vvavr(k) = 0.
         do 400 l=1,4
          nnkslp(k,l) = 0
 400    continue
        do 410 nl=1,nxslp
          k = kslp(nl)
          nnkslp(k,1) = nnkslp(k,1) + 1
 410    continue
        do 415 nl=nxslp+1,nxyslp
          k = kslp(nl)
          nnkslp(k,3) = nnkslp(k,3) + 1
 415    continue
c-
        kkx = nnx
        do 420 nldw=1,nnx
          nl = nslp(nldw)
          k = kslp(nl)
          nnkslp(k,2) = nnkslp(k,2) + 1
          kkx = kkx + k - kslpdw(nldw)
          uuavr(k) = uuavr(k) + abs(uvcslp(nldw))
 420    continue
        kky = nny
        do 425 nldw=nnx+1,nslpdw
          nl = nslp(nldw)
          k = kslp(nl)
          nnkslp(k,4) = nnkslp(k,4) + 1
          kky = kky + k - kslpdw(nldw)
          vvavr(k) = vvavr(k) + abs(uvcslp(nldw))
 425    continue
C       write(iuo+66,'(2A,F10.4,4I6)') ' slopez(U=xslop*dz*Db) : ',
C    &      'xslop, Nx,Kx,Ny,Ky=', xslop, nnx, kkx, nny, kky
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Ecriture sur fichier "xyslope" :
        open(70, file='xyslope', status='unknown' )
        write(70,'(A,I11,A,F10.4)') 'Exp '//refexp//', It', numit,
     &    ' , slopez : xslop=', xslop
        write(70,'(A,6I6)') ' nXslp,Nx,Kx, nYslp,Ny,Ky=',
     &    nxslp, nnx, kkx, nxyslp-nxslp, nny, kky
 
        write(70,'(2A)') '  k   Pot.X,Dwn.X   Pot.Y,Dwn.Y ',
     &                   '    U.moy  (m/s)  V.moy'
        do 430 k=kmax,1,-1
          uuavr(k) = uuavr(k) * unsdz(k)
          vvavr(k) = vvavr(k) * unsdz(k)
          uuavr(0) = uuavr(0) + uuavr(k)
          vvavr(0) = vvavr(0) + vvavr(k)
          if (nnkslp(k,2).ge.1) uuavr(k) = uuavr(k)/DFLOAT(nnkslp(k,2))
          if (nnkslp(k,4).ge.1) vvavr(k) = vvavr(k)/DFLOAT(nnkslp(k,4))
          write(70,'(I3,2(2X,2I6),2X,1P2E11.3)')
     &        k, (nnkslp(k,l),l=1,4), uuavr(k), vvavr(k)
         do 430 l=1,4
          nnkslp(0,l) = nnkslp(0,l) + nnkslp(k,l)
 430    continue
          if (nnkslp(0,2).ge.1) uuavr(0) = uuavr(0)/DFLOAT(nnkslp(0,2))
          if (nnkslp(0,4).ge.1) vvavr(0) = vvavr(0)/DFLOAT(nnkslp(0,4))
          write(70,'(A3,2(2X,2I6),2X,1P2E11.3)')
     &     'T: ', (nnkslp(0,l),l=1,4), uuavr(0), vvavr(0)
        write(70,*)
        write(70,'(A)') ' Max V.slp (cslpmx/dz) :'
        write(70,'(1P8E10.3)') (cslpmx(k)*unsdz(k),k=kmax,1,-1)
        write(70,*)
c--
        write(70,'(A)') '  i   j   k   l    b(l)-b     uv.slp  '
     &                //'      u       xslop.eq'
        do 450 nldw=1,nnx
          nl = nslp(nldw)
          ij = ijslp(nl)
          k = kslp(nl)
          l = lslp(nl)
          iij = ij + max(0,l)
          u1 = DFLOAT(-l) * 0.5 * (uhk(iij,k) + uhk(iij+imax,k))
          jj0 = (ij - 1) / imax
          if ( u1.gt.epsil ) then
            cc00 = max(epsil, bhk(ij+l,k)-bhk(ij,k) )
            cc00 = u1 / cc00
            write(70,'(4I4,1P4E11.3)') ij - imax*jj0, 1+jj0, k, l,
     &       bhk(ij+l,k)-bhk(ij,k), uvcslp(nldw)*unsdz(k), u1, cc00
          else
            write(70,'(4I4,1P4E11.3)') ij - imax*jj0, 1+jj0, k, l,
     &       bhk(ij+l,k)-bhk(ij,k), uvcslp(nldw)*unsdz(k)
          endif
 450    continue
        write(70,*)
c--
        write(70,'(A)') '  i   j   k   l    b(l)-b     uv.slp  '
     &                //'      v       xslop.eq'
        do 460 nldw=nnx+1,nslpdw
          nl = nslp(nldw)
          ij = ijslp(nl)
          k = kslp(nl)
          l = lslp(nl)
          ijj = ij + max(0,l)
          v2 = 0.5 * (vhk(ijj,k) + vhk(ijj+1,k))
          v2 = sign(v2, DFLOAT(-l))
          jj0 = (ij - 1) / imax
          if ( v2.gt.epsil ) then
            cc00 = max(epsil, bhk(ij+l,k)-bhk(ij,k) )
            cc00 = v2*dz(k) / cc00
            write(70,'(4I4,1P4E11.3)') ij - imax*jj0, 1+jj0, k, l,
     &       bhk(ij+l,k)-bhk(ij,k), uvcslp(nldw)*unsdz(k), v2, cc00
          else
            write(70,'(4I4,1P4E11.3)') ij - imax*jj0, 1+jj0, k, l,
     &       bhk(ij+l,k)-bhk(ij,k), uvcslp(nldw)*unsdz(k)
          endif
 460    continue
        write(70,*)
c--fin d'ecriture.
        close(70)
        endif
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine slopez -
      end
