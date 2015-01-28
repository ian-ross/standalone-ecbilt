












      subroutine sloptab(uslpfx, vslpfx, valfil, kslpdw, nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Appele par "process".
c Retablit dans 2 tableaux imax X jmax X (2+nsmax)
c    la hauteur de chute et les Transports de masse et de scalaire
c    impliques dans le courant de Dowsloping.
c-----
c  modif : 30/01/98
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c--dummy variables :
      dimension kslpdw(nlpmax)
      dimension uslpfx(imax,jmax,-1:nsmax), vslpfx(imax,jmax,-1:nsmax)
 
c--variables locales equivalentes :
      dimension scalhk(ixjmax,kmax,nsmax)
      equivalence ( scalhk(1,1,1) , scal(1,1,1,1) )
 
c--variables locales :
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation par valfil .
c-----------------------------------------------------------------------
 
c--Debut du remplissage de u/vslpfx :
 
      do 110 n=-1,nsmax
       do 110 j=1,jmax
        do 110 i=1,imax
          uslpfx(i,j,n) = valfil
          vslpfx(i,j,n) = valfil
 110  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Transfert uvslp -> u/v_slpfx(1) .
c-----------------------------------------------------------------------
 
      do 250 n=1,nslpdw
        nl = nslp(n)
        kdw = kslpdw(n)
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        iijj = ij + max(0,l)
        i = 1 + mod(iijj-1,imax)
        j = 1 + (iijj-1) / imax
        dzdw = z(k) - z(kdw)
        uvflow = -l
        uvflow = uvcslp(n) * sign(dx,uvflow)
        if (abs(l).eq.1) then
          uslpfx(i,j,-1) = dzdw
          uslpfx(i,j,0) = uvflow
          do 230 ns=1,nsmax
            uslpfx(i,j,ns) = uvflow
     &                * (scalhk(ij+l,k,ns) - scalhk(ij,k,ns))
 230      continue
        else
          vslpfx(i,j,-1) = dzdw
          vslpfx(i,j,0) = uvflow
          do 240 ns=1,nsmax
            vslpfx(i,j,ns) = uvflow
     &                * (scalhk(ij+l,k,ns) - scalhk(ij,k,ns))
 240      continue
        endif
 250  continue
 
c--Fin du remplissage de u/vslpfx . ---
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine sloptab -
      end
