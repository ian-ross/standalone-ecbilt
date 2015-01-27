c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comocean.h                                        
c *** Contents: Common and parameter declarations for ocean part of 
c ***           ECbilt 
c *** Nothing should be changed to common posit without changing also
c *** posit in comocouthelp.h
c-----------------------------------------------------------------------
      integer   nm,nl,nz,nsela,nband
  
      parameter (nm=31,nl=63,nz=12)
      parameter (nsela=(nl+1)*(nm+1),nband=nl+nl+3)

      logical refi,mico,ltflux  

      real*8  ronil,grav,omega,pi,ra,rad,tzero

      real*8  fi(0:nm),fiv(0:nm+1)
      real*8  cofi(0:nm),cofiv(0:nm+1),dx(0:nm),dy
      real*8  f(0:nm),fv(0:nm+1),beta(0:nm),betav(0:nm+1)
      real*8  hms,h(nz),hm(0:nz),hbtrop(nz)
      real*8  rkap,rkapv,rkaph,cons,cont 
      real*8  dt,dto,dl,df,rl,rf
      integer ntot,ibstep,ifstep
      real*8  cote1,cosa1,cote2,cosa2,cote3,cosa3,tstar
      real*8  gamma
      real*8  u(0:nl+1,0:nm+1,nz),v(0:nl+1,0:nm+1,nz)
      real*8  um(0:nl+1,0:nm+1),vm(0:nl+1,0:nm+1)
      real*8  ud(0:nl+1,0:nm+1,0:nz+1),vd(0:nl+1,0:nm+1,0:nz+1)
      real*8  w(0:nl,0:nm,0:nz)
      real*8  ro(-1:nl+1,-1:nm+1,nz),rop(-1:nl+1,-1:nm+1,nz)
      real*8  sa(-1:nl+1,-1:nm+1,nz),sap(-1:nl+1,-1:nm+1,nz)
      real*8  te(-1:nl+1,-1:nm+1,nz),tep(-1:nl+1,-1:nm+1,nz)
      real*8  teatm(0:nl,0:nm),toatm(0:nm),saatm(0:nl,0:nm)
      real*8  tau(0:nl+1,0:nm),tauv(0:nl+1,0:nm+1)
      real*8  tav(0:nl+1,0:nm),tavv(0:nl+1,0:nm+1)
      real*8  relafs(-1:nl+1,-1:nm+1),relaft(-1:nl+1,-1:nm+1)
      real*8  relaxt(-1:nl+1,-1:nm+1,nz),relaxs(-1:nl+1,-1:nm+1,nz) 
      real*8  relaxtt(-1:nl+1,-1:nm+1,nz),relaxst(-1:nl+1,-1:nm+1,nz)   
      real*8  relaxtall(-1:nl+1,-1:nm+1),relaxsall(-1:nl+1,-1:nm+1)   
      real*8  hfluxs(-1:nl+1,-1:nm+1,nz),hfluxt(-1:nl+1,-1:nm+1,nz)
      real*8  vfluxs(-1:nl+1,-1:nm+1,nz),vfluxt(-1:nl+1,-1:nm+1,nz)
      real*8  hdifs(-1:nl+1,-1:nm+1,nz),hdift(-1:nl+1,-1:nm+1,nz)
      real*8  vdifs(-1:nl+1,-1:nm+1,nz),vdift(-1:nl+1,-1:nm+1,nz)
      real*8  dcont(-1:nl+1,-1:nm+1,nz),dcons(-1:nl+1,-1:nm+1,nz)
      real*8  teta(-1:nl+1,-1:nm+1,nz-1),accstr 
      real*8  rkiv(0:nl+1,0:nm+1),rkjv(0:nl+1,0:nm+1)
      real*8  rkivs(0:nl+1,0:nm+1),rkjvs(0:nl+1,0:nm+1)

      integer ito(0:nl+1,0:nm+1), iback(nsela), jback(nsela)   
      integer indu(0:nl+1,0:nm+1),isea(-1:nl+1,-1:nm+1) 
      integer lsdistr(0:nl+1,0:nm+1)
      integer ico(0:nl,0:nm,nz-1),isto(0:nl,0:nm,nz-1)  
      integer ndreloct,ndrelocs

      logical sea(0:nl,0:nm)


      common /const/ ronil,grav,omega,pi,ra,rad,tzero,
     *             fv,beta,betav
      common /posit / fi,fiv,dx,dy,cofi,cofiv 
      common /para /  rkap,rkapv,rkaph,cons,cont,accstr 
      common /iparoc / ndreloct,ndrelocs 
      common /discr / dt,dto,dl,df,rl,rf,ntot,ibstep,ifstep
c      common /outpar / nout,mout,mchange,mwork 
c      common /tcount / calctime     
      common /statco/ cote1,cosa1,cote2,cosa2,cote3,cosa3,tstar
c      common /conpar / alocon,alcon   
      common /filt/ gamma
      common /logic/ mico,refi,ltflux

      common /depth/ hms,h,hm,hbtrop
      common /trans/ ito, iback, jback   
      common /masker/ lsdistr,indu,isea,sea,rkiv,rkjv ,rkivs,rkjvs 

      common /veloc / u,v,w
      common /vmod / um,vm,ud,vd
      common /ctracer/ ro,rop,sa,sap,te,tep

      common /atmf / teatm,toatm,saatm,tau,tauv,
     *             tav,tavv
      common /relax / relafs,relaft,relaxt,relaxs, 
     *              relaxtt,relaxst,relaxtall,relaxsall  

      common /con/ ico,isto  
      common /roclin / teta
 
	common /dstcon/ dcont,dcons
      common /advection/ hfluxs,hfluxt,vfluxs,vfluxt
      common /diffusion/ hdifs ,hdift ,vdifs ,vdift
