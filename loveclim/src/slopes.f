












      subroutine slopes(scathk, alphhk, cdtsxz, kslpdw,  ns)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  appelee par "scale" ; pour le scalaire "ns" ,
c Calcule et integre ds scadt le courant de DownSloping .
c  modif : 16/05/96
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c--variables locales equivalentes :
      dimension smxyhk(ixjmax), scalhk(ixjmax,kmax,nsmax)
      equivalence ( smxyhk(1) , smxy(1,1,0) )
      equivalence ( scalhk(1,1,1) , scal(1,1,1,1) )
 
c--dummy variables :
      dimension scathk(ixjmax,kmax)
      dimension alphhk(ixjmax,kmax,2)
CKO   dimension uvcslp(nlpmax), alpslp(nlpmax), cdtsxz(kmax)
      dimension cdtsxz(kmax)
      dimension kslpdw(nlpmax) ! CKO: , nslp(nlpmax)
 
c--variables locales :
 
C     if (numit.eq.nstart) then
C       write(iuo+66,*) 'slopes : ns = ', ns
C     endif
 
c- initialise ds scale : cdtsxz(k) =  dts(k) * unsdx * unsdz(k)
c   ATTENTION :  dx = dy  indispensable !
 
      if (xslop.lt.epsil) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Modification du decentrement (xslop = 0) :
c-----------------------------------------------------------------------
 
      imax1p = imax + 1
      do 550 nldw=1,nslpdw
        ij = ijslp(nslp(nldw))
        kk = kslp(nslp(nldw))
        l  = lslp(nslp(nldw))
c- decentrement = 1 :
        lx = mod((imax1p+l),imax) - 1
        ijc = ij + max(-lx,0) + max(l,lx)
        nnc = 2 - abs(lx)
        alphhk(ijc,kk,nnc) = alpslp(nldw)
 550  continue
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Prise en compte du Flux (xslop > 0 ) :
c-----------------------------------------------------------------------
 
      imax1p = imax + 1
      do 650 nldw=1,nslpdw
        ij = ijslp(nslp(nldw))
        kk = kslp(nslp(nldw))
        l  = lslp(nslp(nldw))
        if (kslpdw(nldw).le.kk) then
c- permutation de la quantite : uvcslp(nl)*dt/dx - de kslpdw a kslp :
          zzslp = smxyhk(ij) * uvcslp(nldw)
          sscal = scalhk(ij+l,kk,ns)
          do 630 k=kslpdw(nldw),kk
            scathk(ij,k) = scathk(ij,k)
     &         + zzslp * cdtsxz(k) * (sscal-scalhk(ij,k,ns))
            sscal = scalhk(ij,k,ns)
 630      continue
          scathk(ij+l,kk) = scathk(ij+l,kk) + smxyhk(ij+l) *
     &      uvcslp(nldw) * cdtsxz(kk) * (sscal-scalhk(ij+l,kk,ns))
        endif
 650  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine slopes -
      end
