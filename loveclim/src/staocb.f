












      subroutine staocb(nnt,ccfile)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  redemarrage a partir de l'etat definit par le fichier binaire "rest.om" 
c  et de conditons arbitraires sur la glace.
c  modif : 20/12/93
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
      include 'moment.com'
      include 'comunit.h'
 
c- dummy variables :
      character*(*) ccfile
c- variables locales equivalentes :
      dimension wloc(imax,jmax,kmax)
      equivalence ( w(1,1,1) ,  wloc(1,1,1) )
C     dimension tked(imax,jmax,kmax)
C     equivalence ( tke(1,1,1), tked(1,1,1) )
c- local varaibles :
      character*6 cc6exp

      kmaxp1 = kmax + 1
      if (nnt.ge.2) then
c--initialisation de w(kmax+1) :
        do 130 j=1,jmax
         do 130 i=1,imax
	 w(i,j,kmaxp1) = 0.
130     continue
	endif
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1 ) Bloc d ouverture et de Lecture du fichier binaire "rest.om" .   |
c-----------------------------------------------------------------------
 
      open(unit=60,file=ccfile,status='old',form='UNFORMATTED')
 
c--1.1 lecture de la premiere partie :
c-------------------------------------

      read (60) numit
      read (60) tpstot
      read (60) eta
      read (60) ub
      read (60) vb
      read (60) u
      read (60) v
      read (60) scal
c- option avec ou sans TKE :
C     if (kstart.eq.3) then
C       write(iuo+66,*) 'Pas de lecture de TKE '
C     else
C       read(60) tked
C     endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c correction celcius-kelvin
      do 30 k=1,kmax
       do 20 j=1,jmax
         do 10 i=1,imax
C           scal(i,j,k,1) = scal(i,j,k,1) + 273.15
 10      continue
 20     continue
 30    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (nnt.ge.2) then

c--1.2. lecture de la seconde partie :
c-------------------------------------

        read (60) cc6exp
        read (60) q
        read (60) bvf
        read (60) avsdz
        read (60) avudz
        read (60) wloc
        read (60) egajc
        read (60) hmajc
c--fin de la seconde partie .
        write(iuo+66,'(3A,I11,A)') 'Exp '//refexp//' , File ', ccfile,
     &      ' , Iter', numit, ' : FULL RESULTS read - OK.'
 
      elseif (refexp.eq.'      ') then
        read (60) refexp
        write(iuo+66,'(3A,I11,A)') 'Lecture de ', ccfile, ' terminee (Exp '
     &                    //refexp//', Iter', numit, ' ) .'
      else
        write(iuo+66,'(3A,I11,A)') 'Lecture de ', ccfile, ' terminee (Iter',
     &       numit, ' ) .'
      endif
 
      close(60)
 
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2 ) Traitement du cas d'un changement de stockage .                 |
c-----------------------------------------------------------------------
 
      if(nnt.lt.0) then
        do 300 j=1,jmax
         do 300 i=1,imax
          ub(i,j) = ub(i,j) * hu(i,j)
          vb(i,j) = vb(i,j) * hu(i,j)
 300    continue
      endif

c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3 ) Sea ice initial conditions                                      |
c-----------------------------------------------------------------------
c
      numit=0
C     tpstot=86400*150
      tpstot=0.0
c
c  File of parameters.
c
      open(12,file='inice.param')
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)

      read(12,*) 
      read(12,*) ttest
      read(12,*) 
      read(12,*) hninn
      read(12,*) 
      read(12,*) hginn
      read(12,*) 
      read(12,*) alinn
      read(12,*) 
      read(12,*) hnins
      read(12,*) 
      read(12,*) hgins
      read(12,*) 
      read(12,*) alins
      close(12)

	do 410 j=jeq,js2
	  do 400 i=is1(j),is2(j)
c
c  Northern hemisphere.
c
c  tfu: Melting point of sea water.
c
C
            tfu(i,j)  = abs(273.15-0.0575*scal(i,j,ks2,2)+
     &                  1.710523e-03*sqrt(scal(i,j,ks2,2))**3-
     &	                2.154996e-04*scal(i,j,ks2,2)**2)
c
c  	Criterion for presence (zidto=1) or absence (zidto=0) of ice.
c
	    zidto     = tms(i,j,ks2)*(1.0-max(zero,
     &  	        sign(one,scal(i,j,ks2,1)-tfu(i,j)-ttest)))
c
	    scal(i,j,ks2,1)  = zidto*tfu(i,j)
     &                         +(1.0-zidto)*scal(i,j,ks2,1)
	    scal(i,j,ks2-1,1)  = zidto*tfu(i,j)
     &                           +(1.0-zidto)*scal(i,j,ks2-1,1)
	    hgbq(i,j)   = zidto*hginn
            albq(i,j)   = zidto*alinn+(1.0-zidto)*1.0
            hnbq(i,j)   = zidto*hninn
            ts(i,j)     = tfu(i,j)
            tbq(i,j,1)  = tfu(i,j)
            tbq(i,j,2)  = tfu(i,j)
            tbq(i,j,3)  = tfu(i,j)
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0
            xzo(i,j)    = 0.001d0
	    ug(i,j)     = 0.0
	    vg(i,j)     = 0.0
c  	Moments for advection.
	  do k=1,35
	    vicmom(i,j,k) =0.0
          enddo
c
400       continue
410    continue

c        write(99,*) tfu(10,js2-10),scal(10,js2-10,ks2,1)
c    &               ,scal(10,js2-10,ks2,2)
c           write(99,*) hginn,ttest,hgbq(10,js2-10)

c
c  Southern hemisphere.
c
        do 430 j=js1,jeq-1                                                      
          do 420 i=is1(j),is2(j)
c
            tfu(i,j)  = abs(273.15-0.0575*scal(i,j,ks2,2)+           
     &                  1.710523e-03*sqrt(scal(i,j,ks2,2))**3-
     &                  2.154996e-04*scal(i,j,ks2,2)**2)
	    zidto     = tms(i,j,ks2)*(1.0-max(zero,
     &	                sign(one,scal(i,j,ks2,1)-tfu(i,j)-ttest)))
c
	    scal(i,j,ks2,1) = zidto*tfu(i,j)
     &                        +(1.0-zidto)*scal(i,j,ks2,1)
	    scal(i,j,ks2-1,1) = zidto*tfu(i,j)
     &                          +(1.0-zidto)*scal(i,j,ks2-1,1)
	    scal(i,j,ks2-2,1) = zidto*tfu(i,j)
     &                          +(1.0-zidto)*scal(i,j,ks2-2,1)
	    hgbq(i,j)   = zidto*hgins
            albq(i,j)   = zidto*alins+(1.0-zidto)*1.0
            hnbq(i,j)   = zidto*hnins
            ts(i,j)     = tfu(i,j)
            tbq(i,j,1)  = tfu(i,j)
            tbq(i,j,2)  = tfu(i,j)
            tbq(i,j,3)  = tfu(i,j)
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0
            xzo(i,j)    = 0.001d0
	    ug(i,j)     = 0.0
	    vg(i,j)     = 0.0
c  Moments for advection.
c
            do k=1,35
	      vicmom(i,j,k) =0.0
	    enddo
420	  continue
430      continue
      if (ltest.ge.1) then
c--raccord cyclique pour hgbq,albq,hnbq,ts,tbq,firg,fcsg,fleg,
c                     fsbbq,qstobq,scal:
      call raccord(hgbq(1,1),0.0,1,8)
      call raccord(albq(1,1),0.0,1,8)
      call raccord(hnbq(1,1),0.0,1,8)
      call raccord(ts(1,1),0.0,1,8)
      call raccord(tbq(1,1,1),0.0,nkb0,8)
      call raccord(firg(1,1),0.0,1,8)
      call raccord(fcsg(1,1),0.0,1,8)
      call raccord(fleg(1,1),0.0,1,8)
      call raccord(fsbbq(1,1),0.0,1,8)
      call raccord(qstobq(1,1),0.0,1,8)
      call raccord(scal(1,1,ks2,1),0.0,1,8)
      call raccord(scal(1,1,ks2-1,1),0.0,1,8)
      call raccord(scal(1,1,ks2-2,1),0.0,1,8)
      endif

      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine staocb -
      end
