












c
      subroutine acrlbq(kideb,kiut)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
c Cette sous-routine determine les variations de l'epaisseur,
c de la concentration et du contenu thermique de la glace marine
c dues aux processus d'accretion laterale. Elle ajuste egalement
c l'epaisseur de la couche de neige.
c Derniere modification 01/04/94
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'thermo.com'
 
c
      zeps0=1.e-13
      zeps=1.e-20
c
      do 5 ji=kideb,kiut
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1) Accretion laterale (modification de la fraction d'eau libre).    |
c-----------------------------------------------------------------------
c
      zhemis=max(zero,sign(one,albqb(ji)-2.0))
Cchg  hgcri=hgcrit(1+int(zhemis))
      hgcri= hgcolb(ji)
      acri=acrit(1+int(zhemis))
      albqb(ji)=albqb(ji)-2.0*zhemis
      zain=albqb(ji)
      ia=1-max(0,int(sign(1.5*one,zain-1.0+zeps0)))
c     zaf=zain+(qlbqb(ji)-qcmbqb(ji))/(xlg*hgcri)
c     albqb(ji)=max(zaf,acri)
c     zdhb=(albqb(ji)-zaf)*hgcri/(1.0-acri)
c
      xcri      = min(one,(1.0-albqb(ji))/(1.0-acri))
c
Cnld  ziacri    = 1.0
      ziacri    = 1.0
c
c     always use a value of ziacri equal to 1.0
c     users who wish to use the old parameterisation 
c     (with an exponent exld) for the case when
c     the option Cnld is enabled should uncomment the 
c     line below.
c
cCnld  ziacri   = (1.0-xcri**exld)**(1.0/exld)
c
      zaf       = zain+ziacri*(qlbqb(ji)-qcmbqb(ji))/(xlg*hgcri)
      albqb(ji) = max(zaf,acri)
      zdhb      = (albqb(ji)-zaf)*hgcri/(1.0-acri)-
     &            ((1.0-ziacri)/(1.0-albqb(ji)))*
     &            ((qlbqb(ji)-qcmbqb(ji))/xlg)
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2) Ajustement de l'epaisseur de la couche de neige.                 |
c-----------------------------------------------------------------------
c
      hnbqb(ji)=(1.-zain)*hnbqb(ji)/(1.-albqb(ji))
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3) Ajustement de l'epaisseur et du contenu thermique de la glace.   |
c-----------------------------------------------------------------------
c
      zhgold=hgbqb(ji)
      zhg=ia*zhgold+(1-ia)*hgcri
      hgbqb(ji)=((1.-zain)*zhg+(zain-albqb(ji))*hgcri)/(1.-albqb(ji))
      i=1-max(0,int(sign(1.5*one,hgcri-hgbqb(ji))))
      zdh1=max(zero,zhgold-hgbqb(ji)/2.)
      zdh2=max(zero,-zhgold+hgbqb(ji)/2.)
      zdh3=max(zero,hgbqb(ji)-zhgold/2.)
      zdh4=max(zero,-hgbqb(ji)+zhgold/2.)
      zdh5=max(zero,hgcri-hgbqb(ji)/2.)
      tint=i*((zhgold/2.-zdh3)*tbqb(ji,3)+zdh4*tbqb(ji,2))/
     &     max(zeps,hgbqb(ji)-hgcri)+(1-i)*tfub(ji)
      zat1=(1.-zain)*tbqb(ji,2)*ia
      zat2=(1.-zain)*tbqb(ji,3)*ia
      zat3=(1.-zain)*tint*ia
      zat4=(zain-albqb(ji))*tfub(ji)
      zah=(1.-albqb(ji))*hgbqb(ji)/2.
      tbqb(ji,2)=(min(hgbqb(ji),zhgold)/2.*zat1
     &           +(1-i)*(zhgold/2.-zdh1)
     &           *zat2+(i*(hgbqb(ji)/2.-hgcri+zdh5)+
     &           (1-i)*zdh2)*zat3
     &           +min(hgbqb(ji)/2.,hgcri)*zat4)/zah
      tbqb(ji,3)=(i*(hgbqb(ji)/2.-zdh3)*zat1+
     &           (i*zdh3+(1-i)*zdh1)*zat2
     &           +(i*(hgbqb(ji)/2.-zdh5)+(1-i)*
     &           (hgbqb(ji)/2.-zdh1))*zat3
     &           +(i*zdh5+(1-i)*hgbqb(ji)/2.)*zat4)/zah
      qstbqb(ji)=(1.-zain)/(1.-albqb(ji))*qstbqb(ji)
      zhgodp=hgbqb(ji)
      hgbqb(ji)=hgbqb(ji)+zdhb
      tbqb(ji,2)=(zhgodp*tbqb(ji,2)+(hgbqb(ji)-zhgodp)*tbqb(ji,3))
     &           /hgbqb(ji)
      tbqb(ji,3)=((2.*zhgodp-hgbqb(ji))*tbqb(ji,3)+2.*zdhb*tfub(ji))
     s          /hgbqb(ji)
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  4) Calcul des variations du volume et de la masse de glace.         |
c-----------------------------------------------------------------------
c
      dvlbqb(ji)=(1.-albqb(ji))*hgbqb(ji)-(1.-zain)*zhgold
      dmgbqb(ji)=dmgbqb(ji)+rhog*dvlbqb(ji)
c
 5    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine acrlbq -
c
      return
      end
