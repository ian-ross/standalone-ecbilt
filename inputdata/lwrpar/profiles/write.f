      program writes
  
c *** this program produces inputfiles for the ECBilt lwr parameterization:
c *** It reads data from:
c ***   ccECBILT.asc                 monthly climatology of total cloud cover
c ***   ../../berg.dat               orography in [m]
c ***   ../../../ocean/lakemask.dat  land-sea mask
c ***   prof.*                       reference climatology for temperature,
c ***                                moisture, LWR fluxes, GHG concentrations
c ***                                and coefficients used in the
c ***                                parameterisation 
c ***
c *** It outputs data to:
c *** ../lwrref.dat                  monthly reference climatology
c *** ../lwrcoef.dat                 coeficients for lwr parameterisation
c *** ../ecbiltberg.grads            ECBilt topography to inspect in GrADS
c *** ../ecbiltregions.asc           region definition at ECBilt grid: can be
c ***                                plotted in directory ecclio/maps
c *** ../ccECBILT.dat                total cloud cover to inspect in GrADS

      implicit none

      INCLUDE 'gpar.cbl'
      INCLUDE 'greenf0.h'

      integer lev2,i,j,k,l,m,ireg,is,is1,region,irec,ilev,ija,nlon,nlat
      integer ijatm
      parameter (nlat=32,nlon=64,ijatm=nlat*nlon)
      real*4  tncep(19,27,12),qancep(27,12),ghgipcc(19)
      real*4  o3echam4(20,27,4),pecham4(20,27,4),z500ncep(27,12)
      real*4  ccisccp(32,64,12)
      real*4  lwrref(7,27,4,0:1)
      real*4  pncep(17),lwrqts(7,4,27,4,0:1)
      real*4  lwrt(7,18,27,4,0:1),lwrts(7,4,27,4,0:1)
      real*4  lwrghg(7,19,27,4,0:1),lwrqa(7,27,4,0:1)
      integer ir(32,64),ipl(27),irn(32,64,2)
      REAL*4  pt(nlm)
      REAL*4  plr(nlm),tlr(nlm),o3lr(nlm),pwc500,z500,pisccp(27)
      REAL*4  psr,tsr,fur(nfluxl,0:ntype),fdr(nfluxl,0:ntype)
      REAL*4  gfqts(7,0:1,4)
      real*8  oro(nlat,nlon),rg,fracto(nlat,nlon),fractocn(ijatm),epss
      parameter(epss=1e-10)
      integer ioro(nlat,nlon)

      character ls
      INTEGER  flevs(5), levc
      Character*3  flevsc(5)
      integer  irind(27)
      DATA irind /9,9,8,8,7,7,6,6,12,12,11,11,10,10,4,4,3,3,2,2,1,1,
     *            17,16,15,14,13/
      DATA ghgipcc /356.,1714.,311.,268.,503.,82.,20.,5.,122.,0.,0.,
     *              0.,0.,2.,6.,0.,0.,132.,135./
      DATA pncep /10.,20.,30.,50.,70.,100.,150.,200.,250.,300.,400.,
     *500.,600.,700.,850.,925.,1000./
      DATA flevs/350,400,450,500,550/

c      WRITE(*,*)'Which level for ECBILT fluxlevel 3 ?'
c      DO i=1,5
c        WRITE(*,33)'(',i,') ',flevs(i),' hPa'
c      ENDDO
33    FORMAT(A1,I1,A2,I3,A4)
c      READ(*,*)lev2
       lev2=2
      IF(lev2.NE.4)THEN
c       WRITE(*,*)'For which parameter'
c       WRITE(*,34)' (1) both reference fluxes and coefficients at ',
c     &      flevs(lev2),' hPa'
c       WRITE(*,35)'(2) reference fluxes at ',
c     &      flevs(lev2),' hPa, coefficients at 500 hPa'
c       READ(*,*)levc
        levc=2
      ELSE
       levc=1
      ENDIF
34    FORMAT(A47,I3,A4)
35    FORMAT(A25,I3,A29)
      IF(levc.EQ.1)THEN
        levc=lev2
      ELSE
        levc=4
      ENDIF

      DO k=1,17
        pncep(k)=pncep(k)*100.
      ENDDO


      OPEN(1,FILE='../ccECBILT.dat',FORM='UNFORMATTED',ACCESS='DIRECT',
     &  RECL=(4*32*64))
      OPEN(2,FILE='ccECBILT.asc',FORM='FORMATTED')
      IREC=1
      DO M=1,12
        read(2,*)((ccisccp(i,j,m),j=1,64),i=1,32)
        write(1,REC=IREC)((ccisccp(i,j,m),j=1,64),
     &      i=1,32)
        IREC=IREC+1
      ENDDO
      CLOSE(1)
      close(2)

      open(1,file='../../berg.dat',form='unformatted')
      read(1) oro
      close(1)
      
      IREC=1
      open(2,file='../ecbiltberg.grads',form='unformatted',
     &     access='direct',recl=4*32*64)
      write(2,REC=IREC) ((real(oro(i,j)),j=1,64),i=1,32)
      close(2)
      
      rg=1800
      do j=1,64
        do i=1,32
          if (oro(i,j).gt.rg) then
            ioro(i,j)=1
          else
            ioro(i,j)=0
          endif
        enddo
      enddo
      
c *** land/sea/sea-ice fraction
      open(8,file='../../fractoc.dat')
      read(8,*)
      read(8,*) (fractocn(ija),ija=1,ijatm)
      do ija=1,ijatm
           j=int((ija-1)/nlon)+1
           i=ija-(j-1)*nlon
           fracto(j,i)=fractocn(ija)
           if (fracto(j,i).gt.0.990) fracto(j,i)=1.0d0
      enddo
      close(8)
 
      do i=1,32
        do j=1,64
          if (i.gt.29) then
            ir(i,j)=11
          elseif (i.gt.26) then
            ir(i,j)=10
          elseif (i.gt.24) then
            ir(i,j)=9
          elseif (i.gt.21) then
            ir(i,j)=8
          elseif (i.gt.19) then
            ir(i,j)=7
          elseif (i.gt.13) then
            ir(i,j)=6
          elseif (i.gt.11) then
            ir(i,j)=5
          elseif (i.gt.8) then
            ir(i,j)=4
          elseif (i.gt.6) then
            ir(i,j)=3
          elseif (i.gt.3) then
            ir(i,j)=2
          else
            ir(i,j)=1
          endif
        enddo
      enddo
      do i=1,32
        do j=1,64
            irn(i,j,2)=ir(i,j)*2
            irn(i,j,1)=ir(i,j)*2-1
          if (fracto(i,j).lt.epss) then
            ir(i,j)=ir(i,j)*2
          else
            ir(i,j)=ir(i,j)*2-1
          endif
        enddo
      enddo
          
      do i=1,32
        do j=1,64
          if (oro(33-i,j).gt.rg) then
            if (i.lt.7) then
              ir(33-i,j)=27
              irn(33-i,j,2)=27
            endif
            if (i.gt.7.and.i.lt.14) then
              if (j.gt.32) then
                ir(33-i,j)=26
                irn(33-i,j,2)=26
              endif
              if (j.lt.32) then
                ir(33-i,j)=25
                irn(33-i,j,2)=25
              endif
            endif
            if (i.gt.16.and.i.lt.25) then
              ir(33-i,j)=24
              irn(33-i,j,2)=24
            endif
            if (i.gt.25) then
              ir(33-i,j)=23
              irn(33-i,j,2)=23
            endif
          endif
        enddo
      enddo
      
      open(2,file='../ecbiltregions.asc')
      do i=1,32
        write(2,101) (ir(33-i,j),j=1,64)
      enddo
      do i=1,32
        write(2,101) (irn(33-i,j,1),j=1,64)
      enddo
      do i=1,32
        write(2,101) (irn(33-i,j,2),j=1,64)
      enddo
      close(2)
 100  format(64i1)
 101  format(64i3)
      
      

      do ireg=1,27
        if (mod(ireg,2).eq.1) then
          ls='s'
        else
          ls='l'
        endif
        is1=0
        do is=1,12,3
          is1=is1+1
          call readref(irind(ireg),ls,is,lev2,levc,pt,tlr,plr,o3lr,
     &                 psr,tsr,fur,fdr,pwc500,z500,gfqts)
          ipl(ireg)=1
          DOWHILE (pt(ipl(ireg)).LT.psr.AND.ipl(ireg).LT.18)
            ipl(ireg)=ipl(ireg)+1
          ENDDO

          pisccp(ireg)=psr*100.
          do k=1,20
            pecham4(k,ireg,is1)=plr(k)*100.
            o3echam4(k,ireg,is1)=o3lr(k)
          enddo
          qancep(ireg,is)=pwc500*0.001
          do l=0,1
            lwrref(1,ireg,is1,l)=fur(1,l)
            lwrref(2,ireg,is1,l)=fur(2,l)
            lwrref(3,ireg,is1,l)=fur(3,l)
            lwrref(4,ireg,is1,l)=fur(4,l)
            lwrref(5,ireg,is1,l)=fdr(2,l)
            lwrref(6,ireg,is1,l)=fdr(3,l)
            lwrref(7,ireg,is1,l)=fdr(4,l)
            ilev=3
            do k=1,18
              lwrt(1,k,ireg,is1,l)=greenfuat(1,k,l)
              lwrt(2,k,ireg,is1,l)=greenfuat(2,k,l)
              lwrt(3,k,ireg,is1,l)=greenfuat(3,k,l)
              lwrt(4,k,ireg,is1,l)=0.
              lwrt(5,k,ireg,is1,l)=greenfdat(2,k,l)
              lwrt(6,k,ireg,is1,l)=greenfdat(3,k,l)
              lwrt(7,k,ireg,is1,l)=greenfdat(4,k,l)
            enddo
            do k=1,4
              lwrts(1,k,ireg,is1,l)=greenfuts(1,l,k)
              lwrts(2,k,ireg,is1,l)=greenfuts(2,l,k)
              lwrts(3,k,ireg,is1,l)=greenfuts(3,l,k)
              lwrts(4,k,ireg,is1,l)=greenfuts(4,l,k)
              lwrts(5,k,ireg,is1,l)=0.
              lwrts(6,k,ireg,is1,l)=0.
              lwrts(7,k,ireg,is1,l)=greenfdts(l,k)
              lwrqts(1,k,ireg,is1,l)=gfqts(1,l,k)*10.
              lwrqts(2,k,ireg,is1,l)=gfqts(2,l,k)*10.
              lwrqts(3,k,ireg,is1,l)=gfqts(3,l,k)*10.
              lwrqts(4,k,ireg,is1,l)=gfqts(4,l,k)*10.
              lwrqts(5,k,ireg,is1,l)=0.
              lwrqts(6,k,ireg,is1,l)=0.
              lwrqts(7,k,ireg,is1,l)=gfqts(7,l,k)*10.
            enddo
            lwrqa(1,ireg,is1,l)=greenfuaqcol(1,l)*10.
            lwrqa(2,ireg,is1,l)=greenfuaqcol(2,l)*10.
            lwrqa(3,ireg,is1,l)=greenfuaqcol(3,l)*10.
            lwrqa(4,ireg,is1,l)=0.
            lwrqa(5,ireg,is1,l)=greenfdaqcol(2,l)*10.
            lwrqa(6,ireg,is1,l)=greenfdaqcol(3,l)*10.
            lwrqa(7,ireg,is1,l)=greenfdaqcol(4,l)*10.
            lwrghg(1,1,ireg,is1,l)=greenfuco2(1,l,1)
            lwrghg(2,1,ireg,is1,l)=greenfuco2(2,l,1)
            lwrghg(3,1,ireg,is1,l)=greenfuco2(3,l,1)
            lwrghg(4,1,ireg,is1,l)=0.
            lwrghg(5,1,ireg,is1,l)=greenfdco2(2,l,1)
            lwrghg(6,1,ireg,is1,l)=greenfdco2(3,l,1)
            lwrghg(7,1,ireg,is1,l)=greenfdco2(4,l,1)
            lwrghg(1,2,ireg,is1,l)=greenfuch4(1,l,1)
            lwrghg(2,2,ireg,is1,l)=greenfuch4(2,l,1)
            lwrghg(3,2,ireg,is1,l)=greenfdch4(3,l,1)
            lwrghg(4,2,ireg,is1,l)=0.
            lwrghg(5,2,ireg,is1,l)=greenfdch4(2,l,1)
            lwrghg(6,2,ireg,is1,l)=greenfdch4(3,l,1)
            lwrghg(7,2,ireg,is1,l)=greenfdch4(4,l,1)
            lwrghg(1,3,ireg,is1,l)=greenfun2o(1,l,1)
            lwrghg(2,3,ireg,is1,l)=greenfun2o(2,l,1)
            lwrghg(3,3,ireg,is1,l)=greenfun2o(3,l,1)
            lwrghg(4,3,ireg,is1,l)=0.
            lwrghg(5,3,ireg,is1,l)=greenfdn2o(2,l,1)
            lwrghg(6,3,ireg,is1,l)=greenfdn2o(3,l,1)
            lwrghg(7,3,ireg,is1,l)=greenfdn2o(4,l,1)
            do k=4,19
              lwrghg(1,k,ireg,is1,l)=greenfughg(1,k,l)
              lwrghg(2,k,ireg,is1,l)=greenfughg(2,k,l)
              lwrghg(3,k,ireg,is1,l)=greenfughg(3,k,l)
              lwrghg(4,k,ireg,is1,l)=0.
              lwrghg(5,k,ireg,is1,l)=greenfdghg(2,k,l)
              lwrghg(6,k,ireg,is1,l)=greenfdghg(3,k,l)
              lwrghg(7,k,ireg,is1,l)=greenfdghg(4,k,l)
            enddo
          enddo

        enddo
      enddo

      do ireg=1,27
        if (mod(ireg,2).eq.1) then
          ls='s'
        else
          ls='l'
        endif
        do is=1,12
          call readreft(irind(ireg),ls,is,tlr,psr,tsr,pwc500,z500)
          do k=1,18
            tncep(k,ireg,is)=tlr(k)
          enddo
          tncep(19,ireg,is)=tsr
          z500ncep(ireg,is)=z500
          qancep(ireg,is)=pwc500*0.001
        enddo
      enddo


      open(2,file='../lwrref.dat',form='unformatted')

      write(2) irn,ipl,pisccp,pncep,z500ncep
      write(2) tncep,qancep,ghgipcc,ccisccp
      write(2) lwrref
      close(2)

      open(2,file='../lwrcoef.dat',form='unformatted')
      write(2) lwrt,lwrts,lwrqts,lwrqa,lwrghg
      close(2)

      open(2,file='../lwro3echam4.dat',form='unformatted')
      write(2) o3echam4,pecham4
      close(2)
1000  format(x)


      end


      SUBROUTINE ReadRef(lband,ls,monthref,lev,levc,plr,tlr,ptemp,
     &o3lr,psr,tsr,fur,fdr,pwc500,z500,gfqts)
      INCLUDE 'gpar.cbl'
      INCLUDE 'gvar.cbl'
      INCLUDE 'greenf0.h'
c      INCLUDE 'vpar.cbl'
      REAL plr(nlm),pl2r(nlm+1),tlr(nlm),tl2r(nlm+1),clr(nlm),
     &     lwclr(nlm),qlr(nlm),o3lr(nlm),t5r(5),q5r(5),t2mp,q2mp
     &,pwc500,z500
      CHARACTER cg
      INTEGER ntspols,npols,k
      REAL psr,tsr,t2mr,q2mr,dum(4)
      REAL fur(nfluxl,0:ntype),fdr(nfluxl,0:ntype)
      REAL fac
      CHARACTER*3 s,n
      INTEGER lband,monthref,lev,levc,choice,lay,i,j,pi(profnlm),nc
      REAL trr(nlm)
      CHARACTER ls
      REAL tc(profnlm),tsc,qc(profnlm)
      real ptemp(20)
      REAL*4  gfqts(7,0:1,4)
      CHARACTER*2 mon(12)
      CHARACTER*3 num(24)
      DATA mon /'01','02','03','04','05','06','07','08','09',
     &          '10','11','12'/
      DATA num /'75n','90n','60n','75n','45n','60n',
     &          '30n','45n','30s','30n','45s','30s',
     &          '60s','45s','75s','60s','90s','75s',
     &          '15n','30n','15s','15n','30s','15s'/

      cg='y'
      IF(cg.EQ.'y')THEN
        nc=ntype
      ELSE
        nc=0
      ENDIF
C Read reference profile and fluxes
C----------------------------------
c      WRITE(*,*)'Band= ',lband
      IF(lband.LE.12)THEN
        s=num(lband*2-1)
        n=num(lband*2)
        pnamels='prof.'//s//'.'//n//'.'//ls//'.'//mon(monthref)
          OPEN(330,FILE=pnamels//'.ref',
     &            STATUS='OLD',ERR=100)
          OPEN(340,FILE=pnamels//'.reff',
     &            STATUS='OLD',ERR=100)
          open(510,file=pnamels//'.gf.tq'
     &            ,status='old')
          open(511,file=pnamels//'.gf.ts'
     &            ,status='old')
          open(512,file=pnamels//'.gf.q'
     &            ,status='old')
          open(515,file=pnamels//'.gf.o3')
          open(520,file=pnamels//'.gf.ghg'
     &            ,status='old')
        GOTO 200
100     CONTINUE
          WRITE(*,*)'Sorry, this file is not in store'
          WRITE(*,*)'Please look in directory first'
          STOP
200     CONTINUE
      ELSE
        GOTO(201,202,203,204,205)(lband-12)
201       OPEN(330,FILE='prof.greenland.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          OPEN(340,FILE='prof.greenland.'//
     &           mon(monthref)//'.reff',STATUS='UNKNOWN')
          open(510,file='prof.greenland.'//
     &            mon(monthref)//'.gf.tq',status='old')
          open(511,file='prof.greenland.'//
     &            mon(monthref)//'.gf.ts',status='old')
          open(512,file='prof.greenland.'//
     &            mon(monthref)//'.gf.q',status='old')
          open(515,file='prof.greenland.'//
     &            mon(monthref)//'.gf.o3')
          open(520,file='prof.greenland.'//
     &            mon(monthref)//'.gf.ghg',status='old')
          GOTO 206
202       OPEN(330,FILE='prof.rockies.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          OPEN(340,FILE='prof.rockies.'//
     &           mon(monthref)//'.reff',STATUS='UNKNOWN')
          open(510,file='prof.rockies.'//
     &            mon(monthref)//'.gf.tq',status='old')
          open(511,file='prof.rockies.'//
     &            mon(monthref)//'.gf.ts',status='old')
          open(512,file='prof.rockies.'//
     &            mon(monthref)//'.gf.q',status='old')
          open(515,file='prof.rockies.'//
     &            mon(monthref)//'.gf.o3')
          open(520,file='prof.rockies.'//
     &            mon(monthref)//'.gf.ghg',status='old')
          GOTO 206
203       OPEN(330,FILE='prof.himalaya.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          OPEN(340,FILE='prof.himalaya.'//
     &           mon(monthref)//'.reff',STATUS='UNKNOWN')
          open(510,file='prof.himalaya.'//
     &            mon(monthref)//'.gf.tq',status='old')
          open(511,file='prof.himalaya.'//
     &            mon(monthref)//'.gf.ts',status='old')
          open(512,file='prof.himalaya.'//
     &            mon(monthref)//'.gf.q',status='old')
          open(515,file='prof.himalaya.'//
     &            mon(monthref)//'.gf.o3')
          open(520,file='prof.himalaya.'//
     &            mon(monthref)//'.gf.ghg',status='old')
          GOTO 206
204       OPEN(330,FILE='prof.andes.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          OPEN(340,FILE='prof.andes.'//
     &           mon(monthref)//'.reff',STATUS='UNKNOWN')
          open(510,file='prof.andes.'//
     &            mon(monthref)//'.gf.tq',status='old')
          open(511,file='prof.andes.'//
     &            mon(monthref)//'.gf.ts',status='old')
          open(512,file='prof.andes.'//
     &            mon(monthref)//'.gf.q',status='old')
          open(515,file='prof.andes.'//
     &            mon(monthref)//'.gf.o3')
          open(520,file='prof.andes.'//
     &            mon(monthref)//'.gf.ghg',status='old')
          GOTO 206
205       OPEN(330,FILE='prof.antarctica.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          OPEN(340,FILE='prof.antarctica.'//
     &           mon(monthref)//'.reff',STATUS='UNKNOWN')
          open(510,file='prof.antarctica.'//
     &            mon(monthref)//'.gf.tq',status='old')
          open(511,file='prof.antarctica.'//
     &            mon(monthref)//'.gf.ts',status='old')
          open(512,file='prof.antarctica.'//
     &            mon(monthref)//'.gf.q',status='old')
          open(515,file='prof.antarctica.'//
     &            mon(monthref)//'.gf.o3')
          open(520,file='prof.antarctica.'//
     &            mon(monthref)//'.gf.ghg',status='old')
          GOTO 206
206     CONTINUE
      ENDIF
      DO i=1,5
        READ(330,*)
      ENDDO
      DO i=1,profnlm+1
        READ(330,4050)lay,plr(i),tlr(i)
      ENDDO
      DO i=19,20
        plr(i)=-1.
        tlr(i)=-1.
      ENDDO
      READ(330,4060)psr,tsr
      pl2r(nlm+1)=psr
      READ(330,4061)pwc500
      READ(330,4063)z500
      DO i=1,12
        READ(330,*)
      ENDDO
      DO i=1,20
        READ(330,4062)ptemp(i),o3lr(i)
      ENDDO
      CLOSE(330)
      DO k=0,ntype
       READ(340,91)fur(1,k)
       READ(340,*)
       READ(340,*)
       READ(340,*)
       READ(340,91)fur(2,k)
       DO i=1,lev+2
         READ(340,91)fur(3,k)
       ENDDO
       DO i=lev+3,8
         READ(340,91)fur(4,k)
       ENDDO
       READ(340,91)fdr(1,k)
       READ(340,*)
       READ(340,*)
       READ(340,*)
       READ(340,91)fdr(2,k)
       DO i=1,lev+2
         READ(340,91)fdr(3,k)
       ENDDO
       DO i=lev+3,8
         READ(340,91)fdr(4,k)
       ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index k...0 is clear sky, 1 is complete overcast
      ENDDO
      CLOSE(340)

C Read Green's functions
C-----------------------
c atm. temp.:
        do k=0,nc
           READ(510,91)(greenfdat(1,j,k),j=1,18)
	   DO j=1,3*18
  	     READ(510,*)
	   ENDDO
           READ(510,91)(greenfdat(2,j,k),j=1,18)
           DO i=1,levc+2
             READ(510,91)(greenfdat(3,j,k),j=1,18)
           ENDDO
           DO i=levc+3,8
             READ(510,91)(greenfdat(4,j,k),j=1,18)
           ENDDO
           READ(510,91)(greenfuat(1,j,k),j=1,18)
	   DO j=1,3*18
  	     READ(510,*)
	   ENDDO
           READ(510,91)(greenfuat(2,j,k),j=1,18)
           DO i=1,levc+2
             READ(510,91)(greenfuat(3,j,k),j=1,18)
           ENDDO
           DO i=levc+3,8
             READ(510,91)(greenfuat(4,j,k),j=1,18)
           ENDDO
           DO j=1,17
             DO i=1,nfluxl
               IF(greenfdat(i,j+1,k).EQ.-1..AND.
     &            greenfdat(i,j,k).NE.-1.)THEN
                 greenfdat(i,18,k)=greenfdat(i,j,k)
                 greenfdat(i,j,k)=0.
                 greenfuat(i,18,k)=greenfuat(i,j,k)
                 greenfuat(i,j,k)=0.
               ENDIF
             ENDDO
           ENDDO
           DO j=1,17
             DO i=1,nfluxl
               IF(greenfdat(i,j,k).EQ.-1.)greenfdat(i,j,k)=0.
               IF(greenfuat(i,j,k).EQ.-1.)greenfuat(i,j,k)=0.
             ENDDO
           ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index j...atm. temperature layers 10,20,30,40,50,70,100,150,200,250,300,
c           400,500,600,700,850,925,1000 and 2m surface air temp.
c index k...0 is clear sky, 1 is complete overcast
        enddo
c surface temperature:
          do j=1,4
            READ(511,91)(greenfdts(k,j),k=0,nc)
c index j...order of polynomial
c index k...0 is clear sky, 1 is complete overcast
          enddo
          do j=1,4
           READ(511,91)(greenfuts(1,k,j),k=0,nc)
	   DO k=1,3*(nc+1)
	     READ(511,*)
	   ENDDO
           READ(511,91)(greenfuts(2,k,j),k=0,nc)
           DO i=1,levc+2
             READ(511,91)(greenfuts(3,k,j),k=0,nc)
           ENDDO
           DO i=levc+3,8
             READ(511,91)(greenfuts(4,k,j),k=0,nc)
           ENDDO
          enddo
c surface temperature/humidity cross-terms:
          do k=0,nc
            do j=1,4
              gfqts(5,k,j)=0.
              gfqts(6,k,j)=0.
            enddo
            READ(511,91)(gfqts(7,k,j),j=1,4)
c index j...order of polynomial
c index k...0 is clear sky, 1 is complete overcast
            READ(511,91)(gfqts(1,k,j),j=1,4)
	    DO j=1,3*4
	      READ(511,*)
	    ENDDO
            READ(511,91)(gfqts(2,k,j),j=1,4)
            DO i=1,levc+2
              READ(511,91)(gfqts(3,k,j),j=1,4)
            ENDDO
            DO i=levc+3,8
              READ(511,91)(gfqts(4,k,j),j=1,4)
            ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index j...order of polynomial
c index k...0 is clear sky, 1 is complete overcast
          enddo
c READ linear q terms
          do k=0,nc
            READ(512,91)greenfdaqcol(1,k)
	    READ(512,*)
	    READ(512,*)
	    READ(512,*)
            READ(512,91)greenfdaqcol(2,k)
            DO i=1,levc+2
              READ(512,91)greenfdaqcol(3,k)
            ENDDO
            DO i=levc+3,8
              READ(512,91)greenfdaqcol(4,k)
            ENDDO
            READ(512,91)greenfuaqcol(1,k)
	    READ(512,*)
	    READ(512,*)
	    READ(512,*)
            READ(512,91)greenfuaqcol(2,k)
            DO i=1,levc+2
              READ(512,91)greenfuaqcol(3,k)
            ENDDO
            DO i=levc+3,8
              READ(512,91)greenfuaqcol(4,k)
            ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index k...0 is clear sky, 1 is complete overcast
          enddo
c READ NON-linear q terms   (replace values for linear terms)
          do k=0,nc
            READ(512,91)greenfdaqcol(1,k)
	    READ(512,*)
	    READ(512,*)
	    READ(512,*)
            READ(512,91)greenfdaqcol(2,k)
            DO i=1,levc+2
              READ(512,91)greenfdaqcol(3,k)
            ENDDO
            DO i=levc+3,8
              READ(512,91)greenfdaqcol(4,k)
            ENDDO
            READ(512,91)greenfuaqcol(1,k)
	    READ(512,*)
	    READ(512,*)
	    READ(512,*)
            READ(512,91)greenfuaqcol(2,k)
            DO i=1,levc+2
              READ(512,91)greenfuaqcol(3,k)
            ENDDO
            DO i=levc+3,8
              READ(512,91)greenfuaqcol(4,k)
            ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index k...0 is clear sky, 1 is complete overcast
          enddo
      CLOSE(510)
      CLOSE(511)
      CLOSE(512)
c  ozone
c      do k=0,nc
c         READ(515,91)((greenfdo3(i,j,k),j=1,nlm),i=1,nfluxl)
c         READ(515,91)((greenfuo3(i,j,k),j=1,nlm),i=1,nfluxl)
c index i...radiation flux level 0,200,500 mb and surface
c index j...ozone layers (see reference file for pressures)
c index k...0 is clear sky, 1 is complete overcast
c      enddo
      CLOSE(515)
c  Greenhouse Gases
      do k=0,nc
         READ(520,91)(greenfdco2(1,k,j),j=1,1)
	 READ(520,*)
	 READ(520,*)
	 READ(520,*)
         READ(520,91)(greenfdco2(2,k,j),j=1,1)
         DO i=1,levc+2
           READ(520,91)(greenfdco2(3,k,j),j=1,1)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfdco2(4,k,j),j=1,1)
         ENDDO
         READ(520,91)(greenfuco2(1,k,j),j=1,1)
	 READ(520,*)
	 READ(520,*)
	 READ(520,*)
         READ(520,91)(greenfuco2(2,k,j),j=1,1)
         DO i=1,levc+2
           READ(520,91)(greenfuco2(3,k,j),j=1,1)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfuco2(4,k,j),j=1,1)
         ENDDO
         READ(520,91)(greenfdch4(1,k,j),j=1,1)
	 READ(520,*)
	 READ(520,*)
	 READ(520,*)
         READ(520,91)(greenfdch4(2,k,j),j=1,1)
         DO i=1,levc+2
           READ(520,91)(greenfdch4(3,k,j),j=1,1)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfdch4(4,k,j),j=1,1)
         ENDDO
         READ(520,91)(greenfuch4(1,k,j),j=1,1)
	 READ(520,*)
	 READ(520,*)
	 READ(520,*)
         READ(520,91)(greenfuch4(2,k,j),j=1,1)
         DO i=1,levc+2
           READ(520,91)(greenfuch4(3,k,j),j=1,1)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfuch4(4,k,j),j=1,1)
         ENDDO
         READ(520,91)(greenfdn2o(1,k,j),j=1,1)
	 READ(520,*)
	 READ(520,*)
	 READ(520,*)
         READ(520,91)(greenfdn2o(2,k,j),j=1,1)
         DO i=1,levc+2
           READ(520,91)(greenfdn2o(3,k,j),j=1,1)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfdn2o(4,k,j),j=1,1)
         ENDDO
         READ(520,91)(greenfun2o(1,k,j),j=1,1)
	 READ(520,*)
	 READ(520,*)
	 READ(520,*)
         READ(520,91)(greenfun2o(2,k,j),j=1,1)
         DO i=1,levc+2
           READ(520,91)(greenfun2o(3,k,j),j=1,1)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfun2o(4,k,j),j=1,1)
         ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index j...number of functional fit function
c           CO2: CONSTANT*LOG(pCO2(t))/LOG(pCO2(t=0))
c           CH4: CONSTANT*(SQRT(pCH4(t))-SQRT(pCH4(t=0)))
C           N2O: CONSTANT*(SQRT(pN2O(t))-SQRT(pN2O(t=0)))
c index k...0 is clear sky, 1 is complete overcast
         READ(520,91)(greenfdghg(1,j,k),j=4,19)
	 DO i=1,3
	  DO j=4,19
  	   READ(520,*)
	  ENDDO
	 ENDDO
         READ(520,91)(greenfdghg(2,j,k),j=4,19)
         DO i=1,levc+2
           READ(520,91)(greenfdghg(3,j,k),j=4,19)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfdghg(4,j,k),j=4,19)
         ENDDO
         READ(520,91)(greenfughg(1,j,k),j=4,19)
	 DO i=1,3
	  DO j=4,19
  	   READ(520,*)
 	  ENDDO
	 ENDDO
         READ(520,91)(greenfughg(2,j,k),j=4,19)
         DO i=1,levc+2
           READ(520,91)(greenfughg(3,j,k),j=4,19)
         ENDDO
         DO i=levc+3,8
           READ(520,91)(greenfughg(4,j,k),j=4,19)
         ENDDO
c index i...radiation flux level 0,200,500 mb and surface
c index j...number of greenhouse gas for which linear fit was made:
C           1: pCO2-CO2 conc. (ppmv)
C           2: pCH4-methane conc. (ppbv)
C           3: pN2O-nitrous-oxide (ppbv)
C
C           pCFC-concentrations CFCs (pptv) per CFC type (pptv):
C           4: pCFC(1) - CFC-11
C           5: pCFC(2) - CFC-12
C           6: pCFC(3) - CFC-113
C           7: pCFC(4) - CFC-114
C           8: pCFC(5) - CFC-115
C
C           pHCFC-concentrations HCFCs (pptv) per type (pptv):
C           9: pHCFC(1) - HCFC-22
C          10: pHCFC(2) - HCFC-123
C          11: pHCFC(3) - HCFC-124
C          12: pHCFC(4) - HCFC-125
C          13: pHCFC(5) - HCFC-134A
C          14: pHCFC(6) - HCFC-141B
C          15: pHCFC(7) - HCFC-142B
C          16: pHCFC(8) - HCFC-143A
C          17: pHCFC(9) - HCFC-152A
C
C          18: pCTC-concentration Carbon TetraChloride (CCl4) (pptv)
C          19: pMCF-concentration Methyl Chloroform (CH3CCl3) (pptv)
C
c index k...0 is clear sky, 1 is complete overcast
      enddo
      CLOSE(520)

      RETURN
91    format(d20.14)
2040  FORMAT(I4,F9.3,2X,F11.9)
2050  FORMAT(I4,F9.3)
2060  FORMAT(F7.3)
2061  FORMAT(F11.9)
4050  FORMAT(I2,F8.2,1X,F8.3)
4060  FORMAT(2X,F8.2,1X,F8.3)
4061  FORMAT(1X,F8.3)
4062  FORMAT(2X,F8.2,1X,E12.6,1X,E12.6)
4063  FORMAT(4X,F5.0)
5000  FORMAT(I2,F8.2,F10.4,F10.4,F8.2,F8.3)
      END

      SUBROUTINE ReadReft(lband,ls,monthref,tlr,psr,tsr,pwc500,z500)
      INCLUDE 'gpar.cbl'
      INCLUDE 'gvar.cbl'
      INCLUDE 'greenf0.h'
c      INCLUDE 'vpar.cbl'
      REAL plr(nlm),pl2r(nlm+1),tlr(nlm),tl2r(nlm+1),clr(nlm),
     &     lwclr(nlm),qlr(nlm),o3lr(nlm),t5r(5),q5r(5),t2mp,q2mp
     &,pwc500,z500
      CHARACTER cg
      INTEGER ntspols,npols,k
      REAL psr,tsr,t2mr,q2mr,dum(4)
      REAL fur(nfluxl,0:ntype),fdr(nfluxl,0:ntype)
      REAL fac
      CHARACTER*3 s,n
      INTEGER lband,monthref,choice,lay,i,j,pi(profnlm),nc
      REAL trr(nlm)
      CHARACTER ls
      REAL tc(profnlm),tsc,qc(profnlm)
      real ptemp(20)
      CHARACTER*2 mon(12)
      CHARACTER*3 num(24)
      DATA mon /'01','02','03','04','05','06','07','08','09',
     &          '10','11','12'/
      DATA num /'75n','90n','60n','75n','45n','60n',
     &          '30n','45n','30s','30n','45s','30s',
     &          '60s','45s','75s','60s','90s','75s',
     &          '15n','30n','15s','15n','30s','15s'/

      cg='y'
      IF(cg.EQ.'y')THEN
        nc=ntype
      ELSE
        nc=0
      ENDIF
C Read reference profile and fluxes
C----------------------------------
c      WRITE(*,*)'Band= ',lband
      IF(lband.LE.12)THEN
        s=num(lband*2-1)
        n=num(lband*2)
        pnamels='prof.'//s//'.'//n//'.'//ls//'.'//mon(monthref)
          OPEN(330,FILE=pnamels//'.ref',
     &            STATUS='OLD',ERR=100)
        GOTO 200
100     CONTINUE
          WRITE(*,*)'Sorry, this file is not in store'
          WRITE(*,*)'Please look in directory first'
          STOP
200     CONTINUE
      ELSE
        GOTO(201,202,203,204,205)(lband-12)
201       OPEN(330,FILE='prof.greenland.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          GOTO 206
202       OPEN(330,FILE='prof.rockies.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
         GOTO 206
203       OPEN(330,FILE='prof.himalaya.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          GOTO 206
204       OPEN(330,FILE='prof.andes.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          GOTO 206
205       OPEN(330,FILE='prof.antarctica.'//
     &           mon(monthref)//'.ref',STATUS='UNKNOWN')
          GOTO 206
206     CONTINUE
      ENDIF
      DO i=1,5
        READ(330,*)
      ENDDO
      DO i=1,profnlm+1
        READ(330,4050)lay,plr(i),tlr(i)
      ENDDO
      DO i=19,20
        plr(i)=-1.
        tlr(i)=-1.
      ENDDO
      READ(330,4060)psr,tsr
      pl2r(nlm+1)=psr
      READ(330,4061)pwc500
      READ(330,4063)z500
      RETURN
91    format(d20.14)
2040  FORMAT(I4,F9.3,2X,F11.9)
2050  FORMAT(I4,F9.3)
2060  FORMAT(F7.3)
2061  FORMAT(F11.9)
4050  FORMAT(I2,F8.2,1X,F8.3)
4060  FORMAT(2X,F8.2,1X,F8.3)
4061  FORMAT(1X,F8.3)
4062  FORMAT(2X,F8.2,1X,E12.6,1X,E12.6)
4063  FORMAT(4X,F5.0)
5000  FORMAT(I2,F8.2,F10.4,F10.4,F8.2,F8.3)
      END

