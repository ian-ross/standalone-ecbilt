c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comice.h                                                    
c *** Contents: Common declarations for ice module of ECbilt
c-----------------------------------------------------------------------
          

      logical ice,icethin,snow
      integer mitetel
      real*8  sigma
      real*8  hsens,elate,qshort,ilw          
      real*8  fs,fb,fa
      real*8  rfluxt,rfluxs
      real*8  tb,tnull,rks,rki,qs,qi,qb,abottom,ho,hsnmax
      real*8  hic(nlat,nlon),hsn(nlat,nlon),brine(nlat,nlon)
      real*8  tmixed,tsice,tatm,hice,hsnow,tflux,sflux
      real*8  dtice,dtis 
      real*8  uvs,alpi,alps,alpv,tster,qgi,cdi
      real*8  tijs(nlat,nlon),Fbtfx(nlat,nlon)


      common /icepara/ sigma
      common /icetest/ mitetel
      common /ijsm/   tijs,Fbtfx
      common /iceatforc/ hsens,elate,qshort,ilw           
      common /icelogic/ ice,icethin,snow
      common /iceflux/ fs,fb,fa
      common /iceflstore/ rfluxt,rfluxs
      common /iceconst/ tb,tnull,rks,rki,qs,qi,qb,abottom,ho,hsnmax
      common /ijsje/ hic,hsn,brine
      common /iceocean/ tmixed,tsice,hice,hsnow,tflux,sflux
      common /positi/ dtice,dtis
      common /iceflxpar/ uvs,alpi,alps,alpv,tster,qgi,cdi
