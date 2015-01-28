












      subroutine informe(nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c    prepare(titres .en tete...) et ouvre le fichier "evolu".
c    The computation are performed in inforun
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 16/11/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'reper.com'
      include 'ice.com'
      include 'dynami.com'
 
c--variables locales :
      dimension surfs(kmax), surfo(kmax), surfv(kmax)
      character*(nchinf) titinf
      character*8 fmtinf, cc8
      character*30 fmtr, fmtitr
 
c--variables locales a conserver d'un appel a l'autre -> dans "reper.com" .
 
c-----
 
 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation, ouverture et ecriture de l'entete du fichier.   |
c-----------------------------------------------------------------------
 
c--selon ntmoy (<-run.param) , definit "ktsum(nv)" pour "inforun" :
c  ntmoy=0 pas de moy. ; =1 moy. de qq. var. ; =2 moy de ttes les var (ss numit)
c  ktsum(nv) : =1 -> moyenne sur ninfo iter., =0 -> output ponctuel
 
c- initialise ktsum(nv) :
      do 100 nv=1,ninfmx
        ktsum(nv) = 0
 100  continue
 
c--facteur d'echelle pour "D.eta" (unite = 10e-6 m/s) :
        zmdeta = 1.D+6
 
c--lecture des 3 1eres lignes du fichier "run.param" :
        open(10,file='run.param',status='old')
        read(10,'(A6)') refexp
        read(10,*)
        read(10,'(A8,1X,I3)') fmtinf, nfrinf
        close(10)
 
c--definition des titres "titvar":
        nv = 1
        titvar(nv) = 'NoIt'
        nv = nv + 1
        titvar(nv) = 'Tmyr'
        if (ntmoy.ge.1) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'EgAjC'
        ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'VmAjC'
        ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'D.Eta'
        ktsum(nv) = 1
        nv = nv + 1
C       titvar(nv) = 'EgEta'
        titvar(nv) = 'M.Eta'
        if (ntmoy.eq.2) ktsum(nv) = 1
c-------
        nv0 = nv
        nv = nv + 1
        titvar(nv) = 'DrHSF'
        nv = nv + 1
        titvar(nv) = 'InHSF'
        nv = nv + 1
        titvar(nv) = 'BeHSF'
        nv = nv0 + nvhsf + 1
CL30    titvar(nv) = 'FLOS'
CL30    nv = nv + 1
        titvar(nv) = 'DANS'
        nv = nv + 1
        titvar(nv) = 'DANN'
        nv = nv + 1
        titvar(nv) = 'ISCS'
        nv = nv + 1
        titvar(nv) = 'ISCN'
        nv = nv + 1
        titvar(nv) = 'FRAS'
        nv = nv + 1
        titvar(nv) = 'FRAN'
        nv = nv + 1
        titvar(nv) = 'PNWS'
        nv = nv + 1
        titvar(nv) = 'PNWN'
        nv = nv + 1
        titvar(nv) = 'ADGIN'
        nv = nv + 1
        titvar(nv) = 'ADPro'
        nv = nv + 1
        titvar(nv) = 'ADOut'
        nv = nv + 1
        titvar(nv) = 'AABpr'
        nv = nv + 1
        titvar(nv) = 'AABex'
        nv = nv + 1
        titvar(nv) = 'AABat'
        nv = nv + 1
        titvar(nv) = 'Fc30A'
        nv = nv + 1
        titvar(nv) = 'Fs30A'
        nv = nv + 1
        titvar(nv) = 'Fsber'
        if (ntmoy.ge.1) then
          do 105 nn=1+nv0,nv
            ktsum(nn) = 1
 105      continue
        endif
c-------
        nv0 = nv
        nv = nv + 1
Cic0    titvar(nv) = 'T'
        titvar(nv) = 'Tmc'
        scalwr(1) = -273.15
        nv = nv + 1
        titvar(nv) = 'T1mo'
        nv = nv + 1
        titvar(nv) = 'ITmoI'
        nv = nv + 1
        titvar(nv) = 'S_30'
        scalwr(2)  =   -30.
        nv = nv + 1
        titvar(nv) = 'S1mo'
        nv = nv + 1
        titvar(nv) = 'ISmoI'
        nv = nv + 1
        titvar(nv) = 'A'
        nv = nv + 1
        titvar(nv) = 'A1mo'
        nv = nv + 1
        titvar(nv) = 'IAmoI'
        nv = nv + 1
        titvar(nv) = 'B'
        nv = nv + 1
        titvar(nv) = 'B1mo'
        nv = nv + 1
        titvar(nv) = 'IBmoI'
c---
        nv = nv0 + 3*nsmax + 1
        titvar(nv) = 'IwI'
        if (ntmoy.eq.2) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'IuI'
        if (ntmoy.eq.2) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'IvI'
        if (ntmoy.eq.2) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'K.E'
        if (ntmoy.eq.2) ktsum(nv) = 1
C       do 110 k=ks1,ks2
C         nv = nv + 1
C         write(titvar(nv),'(A3,I2)') 'ke.',k
C110    continue
        if (ninfo.ge.0) then
c- Moyenne des valeurs absolues des differences :
          do 111 k=ks1,ks2
            nv = nv + 1
            if(k.lt.10) then
              write(titvar(nv),'(A4,I1)') 'Tmom',k
            else
              write(titvar(nv),'(A3,I2)') 'Tmo',k
            endif
 111      continue
          do 112 k=ks1,ks2
            nv = nv + 1
            if(k.lt.10) then
              write(titvar(nv),'(A4,I1)') 'Smom',k
            else
              write(titvar(nv),'(A3,I2)') 'Smo',k
            endif
 112      continue
        else
            nv = nv + 1
            write(titvar(nv),'(A4,I1)') 'Tmom',1
c- Moyenne des differences :
          do 116 k=1+ks1,ks2
            nv = nv + 1
            if(k.lt.10) then
              write(titvar(nv),'(A2,I1,A2)') 'Tm',k,'mo'
            else
              write(titvar(nv),'(A1,I2,A2)') 'T',k,'mo'
            endif
 116      continue
            nv = nv + 1
            write(titvar(nv),'(A4,I1)') 'Smom',1
          do 117 k=1+ks1,ks2
            nv = nv + 1
            if(k.lt.10) then
              write(titvar(nv),'(A2,I1,A2)') 'Sm',k,'mo'
            else
              write(titvar(nv),'(A1,I2,A2)') 'S',k,'mo'
            endif
 117      continue
        endif
        if (ntmoy.eq.2) then
          do 120 nn=1+nv0,nv
            ktsum(nn) = 1
 120      continue
        endif
c-------
 
        nocean = nv
        nv = nv + 1
        titvar(nv) = 'AIEFN'
        nv = nv + 1
        titvar(nv) = 'AIEFS'
        nv = nv + 1
        titvar(nv) = 'A15N'
        nv = nv + 1
        titvar(nv) = 'A15S'
        nv = nv + 1
        titvar(nv) = 'A85N'
        nv = nv + 1
        titvar(nv) = 'A85S'
        nv = nv + 1
        titvar(nv) = 'ALEN'
        nv = nv + 1
        titvar(nv) = 'ALES'
        nv = nv + 1
        titvar(nv) = 'VOLN'
        nv = nv + 1
        titvar(nv) = 'VOLS'
        nv = nv + 1
        titvar(nv) = 'VONN'
        nv = nv + 1
        titvar(nv) = 'VONS'
        nv = nv + 1
        titvar(nv) = 'ECGN'
        nv = nv + 1
        titvar(nv) = 'ECGS'
        nv = nv + 1
        titvar(nv) = 'FRAG'
        nv = nv + 1
        titvar(nv) = 'SPNG'
        nv = nv + 1
        titvar(nv) = 'BERG'
        nv = nv + 1
        titvar(nv) = 'ThEx'
        nv = nv + 1
        titvar(nv) = 'ISMM'
        nv = nv + 1
        titvar(nv) = 'IcbN'
        nv = nv + 1
        titvar(nv) = 'IcbS'
        if (ntmoy.ge.1) then
          do 125 nn=1+nocean,nv
            ktsum(nn) = 1
 125      continue
        endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        nvinfo = nv
        if(nvinfo.gt.ninfmx) then
          write(iuo+66,*) 'Arret ! Depassement Nombre Max de variables'
     &      //' traitees par la routine "informe"'
          stop
        endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Definition et Ecriture de l'entete :
c- nombre d'enregistrements :
        nninfo = abs(ninfo)
        nn0 = (nstart - 1) / nninfo
        if(nstart.eq.1) nn0 = -1
        nferme = (nstart - 1 + nitrun) / nninfo
        nnntot = nferme - nn0
        nn0 = nninfo * (1 + nn0)
        nferme = nninfo * nferme
 
c- definition des formats :
        write(fmtw,'(A,I3,A,I1,A)')
     &    '(',nfrinf,'(A',nchsep,','//fmtinf//'))'
        write(fmtr,'(A,I3,A,I1,A)')
     &    '(',nfrinf,'(',nchsep,'X,'//fmtinf//'))'
        write(fmtitr,'(A,I3,A,I1,A)') '(',ninfmx,'A',nchinf,')'
 
c--Ouverture du fichier "evolu" :
Cray    irecl = nvinfo*nchinf + nchinf
Cray    open(90,file='evolu',status='unknown',RECL=irecl)
        open(90,file='evolu',status='unknown')
 
c- ecriture de 2 lignes d'entete :
        write(90,1000) fmtr, spvr, nvinfo, nnntot, 0, nfrinf
        xxx = 0.001 * DFLOAT(nninfo)
        xx1 = 0.001 * DFLOAT(nn0)
        write(cc8,'(A,I1)') ' ntmoy=', ntmoy
        write(90,1111) DFLOAT(nchinf), 0., xx1, xxx, 0., 0., 0, cc8
 
c- ecriture de 2 lignes de titre :
        write(90,'(A,I8,A,I8,A,I5)')
     &     'Evolution chronologique - Experience '//refexp
     &   //'   de', nn0, ' a', nferme, ' pas', nninfo
        write(90,fmtitr) (titvar(nv),nv=1,nvinfo)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--preparation de "titvar" pour l'ecriture parmi les valeurs numeriques :
        do 130 nv=2,nvinfo
          titinf = titvar(nv)(:nchinf)
          titvar(nv) = '  '//titinf
 130    continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c initialisation of vwx
        do j=1,57
         do k=1,kmax+1
          do n=0,nbsmax
           vwx(j,k,n)=0.0
          enddo
         enddo
        enddo
        write(99,*) 'vwx',vwx(43,1,3)
c--calcul des differents "volumes" intervenant dans les moyennes :
        ddxy = dx * dy
        vols = 0.
        volo = 0.
        volv = 0.
        volw = 0.
        do 156 k=ks1,ks2
          surfs(k) = 0.
          surfo(k) = 0.
          surfv(k) = 0.
          do 154 j=js1,js2
            do 150 i=is1(j),is2(j)
              surfs(k) = surfs(k) + ctmi(i,j,k,0)
              surfo(k) = surfo(k) +
     &          min( ctmi(i,j,k,0), (scalr(i,j,k,1)-spvr) )
 150        continue
            do 152 i=iu1(j),iu2(j)
              surfv(k) = surfv(k) + ctmi(i,j,k,1)
 152        continue
 154      continue
          vols = vols + surfs(k) * dz(k)
          volo = volo + surfo(k) * dz(k)
          volv = volv + surfv(k) * dz(k)
          if(k.lt.ks2) volw = volw + surfs(k) * dzw(k+1)
 156    continue
        zvols = 0.
        zvolo = 0.
        zvolv = 0.
        zvolw = 0.
        if(vols.gt.epsil) zvols = 1. / vols
        if(volo.gt.epsil) zvolo = 1. / volo
        if(volv.gt.epsil) zvolv = 1. / volv
        if(volw.gt.epsil) zvolw = 1. / volw
        vols = vols * ddxy
        volo = volo * ddxy
        volv = volv * ddxy
        volw = volw * ddxy
        do 158 k=ks1,ks2
          if(surfs(k).gt.epsil) then
            zsurfs(k) = 1. / surfs(k)
          else
            zsurfs(k) = 0.
          endif
          if(surfo(k).gt.epsil) then
            zsurfo(k) = 1. / surfo(k)
          else
            zsurfo(k) = 0.
          endif
          if(surfv(k).gt.epsil) then
            zsurfv(k) = 1. / surfv(k)
          else
            zsurfv(k) = 0.
          endif
          surfs(k)  = surfs(k) * ddxy
          surfo(k)  = surfo(k) * ddxy
          surfv(k)  = surfv(k) * ddxy
 158    continue
        zsurf = zsurfs(ks2)
        zmdeta = zmdeta * zsurf
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--initialisation des variables utilisees uniquement pour "informe":
        do 160 j=1,jmax
         do 160 i=1,imax
          daeta(i,j) = 0.
 160    continue
 
c--Remise a zero de fss, Energ. Aj.Conv & Freq.Aj.Conv (avant 1ere it) :
        do 180 k=1,kmax
         do 180 j=1,jmax
          do 180 i=1,imax
            fqajc(i,j,k) = 0.
 180    continue
        do 185 ns=0,nsmax
         do 185 j=1,jmax
          do 185 i=1,imax
           fss(i,j,ns) = 0.
 185    continue
 
c--Initialisation of the arrays for the accumulation
        do 190 nv=1,nvinfo
          vinfom(nv)=0.
 190    continue
 
c--Definition of the limits for the downsloping budget
      ksud = ks2
      zlim1 = -1000.0
      do kz=ks2,ks1,-1
        if (z(kz).ge.zlim1) ksud = kz
      enddo
      knor = ks2
      zlim1 = -450.0
      do kz=ks2,ks1,-1
        if (z(kz).ge.zlim1) knor = kz
      enddo
      yylats = -55.
      jmsud = nint( (yylats-ylat1)/dlat ) + 1
      yylatn = 70.
      jmnor = nint( (yylatn-ylat1)/dlat ) + 1
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        if (nn99.eq.2) then
c--------
          write(99,*) 'nvinfo, nocean, ntmoy :', nvinfo, nocean, ntmoy
          write(99,*) 'zvols , vols :', zvols, vols
          write(99,*) 'zvolo , volo :', zvolo, volo
          write(99,*) 'zvolv , volv :', zvolv, volv
          write(99,*) 'zvolw , volw :', zvolw, volw
          write(99,*) 'zsurfs, surfs (ks2) :', zsurfs(ks2), surfs(ks2)
          write(99,'(1P6E12.5)') (zsurfs(k),k=ks1,ks2)
          write(99,'(1P6E12.5)') ( surfs(k),k=ks1,ks2)
          write(99,*) 'zsurfo, surfo (ks2) :', zsurfo(ks2), surfo(ks2)
          write(99,'(1P6E12.5)') (zsurfo(k),k=ks1,ks2)
          write(99,'(1P6E12.5)') ( surfo(k),k=ks1,ks2)
          write(99,*) 'zsurfv, surfv (ks2) :', zsurfv(ks2), surfv(ks2)
          write(99,'(1P6E12.5)') (zsurfv(k),k=ks1,ks2)
          write(99,'(1P6E12.5)') ( surfv(k),k=ks1,ks2)
	  write(99,'(A,2I4,2F9.2)')' downsloping jmnor,knor,lat_N,zw :'
     &                   , jmnor, knor, ylat1+(jmnor-1)*dlat, zw(knor)
	  write(99,'(A,2I4,2F9.2)')' downsloping jmsud,ksud,lat_S,zw :'
     &                   , jmsud, ksud, ylat1+(jmsud-1)*dlat, zw(ksud)
          write(99,*)
c--------
          write(99,'(A,I4,2A,I4,A)') 'File "evolu" :',nnntot,' time ',
     &    'records of',nvinfo,' (=nvinfo) var. ; nv,title,ktsum :'
          write(99,'(4(A,I3,1X,A,I2,A1))')
     &      (' no=',n, titvar(n), ktsum(n), ',', n=1,nvinfo)
          write(99,*)
c--------
        endif
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine informe -
      end
