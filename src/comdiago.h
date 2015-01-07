c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comdiago.h                                                    
c *** Contents: Common declarations for diagnostics of ocean
c-----------------------------------------------------------------------
      integer nto(20),nsa(20),nuo(20),nvo(20),nwo(20),nro(20),nhadvt(20),
     &        nhadvs(20),nvadvt(20),nvadvs(20),nhdift(20),nhdifs(20),
     &        nvdift(20),nvdifs(20),nconvs(20),nconvt(20),nrelaxt(20),
     &        nrelaxs(20),nheex(20),nsaex(20),nbrin(20),ndice(20),
     &        nsds(20),ntice(20),nheico(20)         
      integer ifreq,isnapshot

      integer iss,mlenth,icont,mtype



      common /xwriteo/ nto,nsa,nuo,nvo,nwo,nro,nhadvt,nhadvs,nvadvt,
     &       nvadvs,nhdift,nhdifs,nvdift,nvdifs,nconvs,nconvt,nrelaxt,
     &       nrelaxs,nheex,nsaex,nbrin,ndice,nsds,ntice,nheico        

      common /xwrctlo/ifreq,isnapshot,iss,mlenth,icont,mtype
