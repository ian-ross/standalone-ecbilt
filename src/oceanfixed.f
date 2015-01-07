c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inioceanfixed
c-----------------------------------------------------------------------
c *** initialises and sets parameters of the fixed ocean model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comphys.h'
      include 'comocfix.h'
      include 'comglobal.h'

      integer      i,j,k,id,ia1,ia9,ica,icd,iap
      real*4       sstday4(nlat,nlon)
      integer      ilake(nlat,nlon)
      character*1  ch(nlon),space

c *** reading lakemask and adapt lsmask accordingly
c *** lakes are numbered 2 to 9, 1 corresponds to ocean, 0 to land

      do j=nlat,1,-1
        read (42,310) k,(ilake(j,i),i=1,nlon)
      enddo

      do j=1,nlon
        do i=1,nlat
          if (ilake(i,j).gt.1) then
            if (lsmask(i,j).eq.1) then
              write(29,*) i,j,' already sea'
              call error(10)
            else
              lsmask(i,j)=1
            endif
          endif
        enddo
      enddo

c
c *** read sst on daily basis
c

      do id=1,360
        read(41) (( sstday4(i,j),j=1,nlon),i=1,nlat)
        do i=1,nlat
          do j=1,nlon
            sstday(i,j,id)=sstday4(i,j)
            if (lsmask(i,j).eq.1) then
              if (sstday(i,j,id).ge.310.or.
     *            sstday(i,j,id).lt.tzero-1.) then
                write(100,*) 'sst out of range ',i,j,sstday(i,j,id)
              endif
            endif
          enddo
        enddo
      enddo

c *** read climatological seaice cover 

      ia1=ichar('1') 
      ia9=ichar('9') 
      ica=ichar('A') 
      icd=ichar('D') 
      iap=ichar('.') 

      do i=nlat,1,-1
        read (43,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          lsicebirth(i,j)=ichar(ch(j))
          if (lsicebirth(i,j).eq.iap.and.lsmask(i,j).eq.1) then
            write(29,*) 'birth icemask latlon ',i,j
            call error(17)
          endif
          if (lsicebirth(i,j).ne.iap.and.lsmask(i,j).eq.0) then
            write(29,*) 'birth icemask latlon ',i,j
            call error(17)
          endif
          if (lsicebirth(i,j).ge.ia1.and.lsicebirth(i,j).le.ia9) then
            lsicebirth(i,j)=lsicebirth(i,j)-ia1+1
          else
            if (lsicebirth(i,j).ge.ica.and.lsicebirth(i,j).le.icd) then
              lsicebirth(i,j)=lsicebirth(i,j)-ica+10
            else
              lsicebirth(i,j)=0
            endif
          endif
        enddo
      enddo

      do i=nlat,1,-1
        read (43,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          lsicedeath(i,j)=ichar(ch(j))
          if (lsicedeath(i,j).eq.iap.and.lsmask(i,j).eq.1) then
            write(29,*) 'death icemask latlon ',i,j
            call error(17)
          endif
          if (lsicedeath(i,j).ne.iap.and.lsmask(i,j).eq.0) then
            write(29,*) 'death icemask latlon ',i,j
            call error(17)
          endif
          if (lsicedeath(i,j).ge.ia1.and.lsicedeath(i,j).le.ia9) then
            lsicedeath(i,j)=lsicedeath(i,j)-ia1+1
          else
            if (lsicedeath(i,j).ge.ica.and.lsicedeath(i,j).le.icd) then
              lsicedeath(i,j)=lsicedeath(i,j)-ica+10
            else
              lsicedeath(i,j)=0
            endif
          endif
          if (lsicedeath(i,j).gt.0) then
            if (lsicebirth(i,j).eq.0.or.
     #          lsicebirth(i,j).eq.lsicedeath(i,j)) then
              write(29,*) 'icemask latlon ',i,j
              call error(17)
            endif
          endif
          if (lsicebirth(i,j).gt.0) then
            if (lsicedeath(i,j).eq.0.or.
     #          lsicebirth(i,j).eq.lsicedeath(i,j)) then
              write(29,*) 'icemask latlon ',i,j
              call error(17)
            endif
          endif
        enddo
      enddo


 100  format(i4,65A1)
310   format(i4,i2,90i1)

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine oceanfixed
c-----------------------------------------------------------------------
c *** prescribe sst's and seaice cover over ocean and lakes
c *** calculate seaice temperature assuming all surface fluxes balance
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comocfix.h'
      include 'comphys.h'
      include 'comglobal.h'
      include 'comcoup.h'
      
      integer index,i,j
      real*8  dum,facwin,facsum
    
      index = int((day+0.5*dt)/(iatm*dt)) + 1

c      write(100,*) 'index in oceanfixed ',index

c *** prescription of seaice cover

      facwin=dabs(180.d0-day)/180.d0
      facsum=1d0-facwin

      do i=1,nlat
        do j=1,nlon 
          lseaice(i,j)=0
          if (lsmask(i,j).eq.1) then
            if (lsicebirth(i,j).gt.0) then
              if (lsicebirth(i,j).lt.lsicedeath(i,j)) then
                if (imonth.ge.lsicebirth(i,j).and.
     *              imonth.lt.lsicedeath(i,j)) lseaice(i,j)=1
              endif
              if (lsicebirth(i,j).gt.lsicedeath(i,j)) then
                if (imonth.ge.lsicebirth(i,j).or.
     *              imonth.lt.lsicedeath(i,j)) lseaice(i,j)=1
              endif
            endif
          endif
        enddo
      enddo

c *** calculate seaice temperatures

      call seaicetemp

c *** prescription of sst
 
      do i=1,nlat
        do j=1,nlon 
          if (lsmask(i,j).eq.1.and.lseaice(i,j).eq.0) then
            stmix(i,j)=sstday(i,j,index)
          endif
        enddo
      enddo
     
      return
      end
       
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine seaicetemp
c-----------------------------------------------------------------------
c *** calculates seaice temperatures assuming that all surface fluxes
c *** balance
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'
      include 'comice.h'
      include 'comocfix.h'

      integer  i,j,itetel,mxtetel,il,jl
      real*8   stice,stice1,stice2,tol,fluxsumice,zbrent

      external fluxsumice

c *** tol is the wanted accuracy of the seaice temperature in degrees

      parameter (tol=0.1)

      common /landpoint/il,jl

      mxtetel = 0

      do j=1,nlon
        do i=1,nlat

          if (lseaice(i,j).eq.1) then    

            il=i
            jl=j
            stice=tzero
            stice1=stice - 1.
            stice2=stice 
            
            call zbrac(fluxsumice,stice1,stice2,itetel)
            if (itetel.eq.100) call error(7)
            tijs(i,j)=zbrent(fluxsumice,stice1,stice2,tol,itetel)
            if (itetel.eq.100) call error(8)
            if (itetel.gt.mxtetel) mxtetel=itetel

c *** in case of temperatures above zero, set
c *** surface temperature to meltpoint

            if (tijs(i,j).gt.330.or.tijs(i,j).lt.200) then
              write(100,*) 'ice temperature out of range ',i,j
              write(100,*) lsicebirth(i,j),lsicedeath(i,j),lseaice(i,j)
              write(100,*) uv10(i,j),tempsg(i,j),q10(i,j)
              write(100,*) dlrads(i,j),hesws(i,j)
            endif

            if (tijs(i,j).gt.tzero) then
	      tijs(i,j)=tzero
            endif

          endif
        enddo
      enddo


c      write(100,*) 'in oceanfixed over ice ',mxtetel
      call flush(100)

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      function fluxsumice(stice) 
c-----------------------------------------------------------------------
c *** computes sum of fluxes between the seaice and the atmosphere
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comphys.h'

      integer il,jl
      real*8  fluxsumice,stice,qsatss,qsat,efluxgp,hfluxgp

      common /landpoint/il,jl

c *** sensible heatflux

      hfluxgp=alphad*cdragv(il,jl)*uv10(il,jl)*(stice-tempsg(il,jl))

c *** latent heat flux

      qsatss=qsat(1.d+5,stice)

      if (stice.gt.tzero) then
        efluxgp=alphav*cdragv(il,jl)*uv10(il,jl)*(qsatss-q10(il,jl))
      else
        efluxgp=alphas*cdragv(il,jl)*uv10(il,jl)*(qsatss-q10(il,jl))
      endif      
       
      if (efluxgp.lt.0) efluxgp=0.

c *** sum of all fluxes

      fluxsumice=hesws(il,jl) - hfluxgp + dlrads(il,jl) -
     &            sboltz*stice**4 - efluxgp


      return
      end

