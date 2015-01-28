












c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_oc2co(ist)


c----------------------------------------------------------------------
c *** communicate oceanic data to the coupler
c----------------------------------------------------------------------
      include 'type.com'

      include 'comcouphelp.h'
      include 'comcoup.h'
      include 'comsurf.h'
      include 'comemic.h'

      integer ix,jy,ist

      include 'para.com'
      include 'bloc.com'
      include 'ice.com'

      real*8 zfld(imax,jmax),zalb,zalbp
      real*8 sst(nlat,nlon),albq1(nlat,nlon),ssi(nlat,nlon)
      real*8 seaalb(nlat,nlon),hic(nlat,nlon),hsn(nlat,nlon)
      real*8 tzero


      tzero=273.15d0

c *** SST
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = scal(ix,jy,ks2,1)
        enddo
      enddo
      call oc2at(zfld,sst)
      
c *** ICE fraction
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = 1.0-albq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,albq1)
c *** STI
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = ts(ix,jy)
        enddo
      enddo
      call oc2at(zfld,ssi)
c *** hic
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = hgbq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,hic)
c *** hsn
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = hnbq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,hsn)
      
      call detseaalb(seaalb)

      do ix = 1, nlon
        do jy = 1, nlat
          fractn(jy,ix,noc)=(1.0d0-albq1(jy,ix))*fracto(jy,ix)
          fractn(jy,ix,nse)=albq1(jy,ix)*fracto(jy,ix)
          tsurfn(jy,ix,nse)=min(tzero,ssi(jy,ix))
          tsurfn(jy,ix,noc)=max(tzero-1.8d0,sst(jy,ix))
          call ec_shine(tzero-0.15,tzero-0.25,tsurfn(jy,ix,nse),
     &                       hic(jy,ix),hsn(jy,ix),zalb,zalbp)
          albesn(jy,ix,nse)=(1.0-couptcc(jy,ix))*zalbp +
     &                           couptcc(jy,ix)*zalb
          albesn(jy,ix,noc)=albocef*seaalb(jy,ix)
        enddo
      enddo
      
      
      return

      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine detseaalb(seaalb)
c-----------------------------------------------------------------------
c *** calculates albedos as a function of the time of the year, linearly
c *** interpolating between seasonal mean values
c *** albsea is the albedo of the open sea
c-----------------------------------------------------------------------
  
      implicit none

      include 'comcouphelp.h'
      include 'comemic.h'

      integer i,j,id1,is1,is2
      real*8  albseaz(nlat),albsea(nlat,4),seaalb(nlat,nlon)
      real*8  sfrac
      
      common /albedoclio/albsea

c *** interpolate between seasonal means

      id1=(imonth-1)*30+iday-14
      if (id1.lt.1) id1=id1+360

      is1=(id1+89)/90
      is2=is1+1
      if (is2.eq.5) is2=1

      sfrac=(id1-((is1-1)*90.+1.))/90.
      
      do j=1,nlat
        albseaz(j)=albsea(j,is1)+(albsea(j,is2)-albsea(j,is1))*sfrac
      enddo

      do j=1,nlon
        do i=1,nlat
          seaalb(i,j)=albseaz(i)
        enddo
      enddo
            
      return
      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine initseaalb
      
      implicit none
      
      include 'comcouphelp.h'
      include 'comunit.h'
      
      integer i,is
      real*8  albsea(nlat,4)
      common /albedoclio/albsea

c *** read climatological zonal mean albedos for each season

      open (iuo+49,file='inputdata/albedo.dat')
      
      read(iuo+49,*)
      do i=1,nlat
	read(iuo+49,45)(albsea(i,is),is=1,4)
      enddo
45    format(4(2x,f7.4))
      close(iuo+49)
      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine oc2at(fin,fout)
      
      implicit none
      
      include 'comcouphelp.h'
      include 'comcoup.h'
      
      integer ji,jk,i,j,k
      real*8  fin(ijocn),fout(nlat,nlon),zsum
      
      ji = 0
      do i=1,nlat
        do j=1,nlon
          zsum = 0.
	  ji=ji+1
          do jk = 1, kamax
            zsum = zsum + wo2a(ji,jk) * fin(indo2a(ji,jk))
          enddo 
          fout(i,j) = zsum
	enddo
      enddo
      
      end
      
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_shine(tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp)
c-----------------------------------------------------------------------
c *** This subroutine computes albedo of snow-sea ice following SHINE &
c *** HENDERSSON-SELLERS [1985] 
c-----------------------------------------------------------------------
c
c  tfsn   : melting point temperature of snow (273.15 in CLIO)
c  tfsg   : melting point temperature of ice (273.05 in CLIO)
c  ts     : surface temperature
c  hgbq   : ice thickness
c  hnbq   : snow thickness
c  zalb   : ice/snow albedo for overcast sky
c  zalbp  : ice/snow albedo for clear sky

      include 'comatm.h'
      include 'comphys.h'
      integer ih
      real*8 tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp
!     real*8 albin,albis,albice,alphd,alphdi,alphs,cgren



Cdriess    albin = 0.45
Cdriess    albis = 0.45
Cdriess    albice = 0.45
C     albice = 0.53
C     alphd  = 0.80
C     alphdi = 0.72
C     alphs  = 0.65
Cdriess    alphd  = 0.72
Cdriess    alphdi = 0.64
Cdriess    alphs  = 0.55
!     albin = 0.43
!     albis = 0.43
!     albice = 0.43
!     alphd  = 0.70
!     alphdi = 0.62
!     alphs  = 0.53
!     cgren = 0.04

c  albin: Albedo of melting ice in the arctic.
c  albis: Albedo of melting ice in the antarctic (SHINE 
c         & HENDERSSON-SELLERS, 1985).
c  albice: Albedo of melting ice.
c  alphd : Albedo of snow (thickness > 0.05m)
c  alphdi: Albedo of thick bare ice
c  alphs : Albedo of melting snow
c  cgren: Correction of the snow or ice albedo to take into account
c         effects of cloudiness (GRENFELL & PEROVICH, 1984)
 
      if (hnbq.gt.0.0) then                                       
                                                                        
c ***  Case of ice covered by snow.
                                                                        
        if (ts.lt.tfsn) then                                      
                                                                       
c ***    Freezing snow.
                                                                        
          if (hnbq.gt.0.05) then                                 
            zalbp = alphd                                        
          else                                                      
            if (hgbq.gt.1.5) then                                
              zalbp = alphdi+(hnbq*(alphd-alphdi)/0.05)          
            else if (hgbq.gt.1.0.and.hgbq.le.1.5) then         
                   al = 0.472+2.0*(alphdi-0.472)*(hgbq-1.0)
            else if (hgbq.gt.0.05.and.hgbq.le.1.0) then        
                   al = 0.2467+(0.7049*hgbq)-(0.8608*(hgbq*hgbq))+
     &                 (0.3812*(hgbq*hgbq*hgbq))                     
            else                                                    
              al = 0.1+3.6*hgbq                                  
            endif                                                   
            if (hgbq.le.1.5) zalbp=al+(hnbq*(alphd-al)/0.05)
          endif                                                     
        else                                                        
c                                                                       
c ***    Melting snow.
c                                                                       
          if (hnbq.ge.0.1) then                                  
            zalbp = alphs                                      
          else                                                      
            zalbp = albice+((alphs-albice)/0.1)*hnbq
          endif                                                     
        endif                                                       
      else                                                          
c                                                                       
c *** Case of ice free of snow.
c                                                                       
        if (ts.lt.tfsg) then                                      
c                                                                       
c ***    Freezing ice.
c                                                                       
          if (hgbq.gt.1.5) then                                  
            zalbp = alphdi                                          
          else if (hgbq.gt.1..and.hgbq.le.1.5) then           
	         zalbp = 0.472+2.*(alphdi-0.472)*(hgbq-1.)       
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then          
                 zalbp = 0.2467+                                        
     &                   (0.7049*hgbq)-(0.8608*(hgbq*hgbq))+
     &                   (0.3812*(hgbq*hgbq*hgbq))                  
          else                                                      
            zalbp = 0.1+3.6*hgbq                                
          endif                                                     
        else                                                        
c                                                                       
c *** Melting ice.
c                                                                   
          if (hgbq.gt.1.5) then                                  
            zalbp = albice                                           
          else if (hgbq.gt.1..and.hgbq.le.1.5)  then          
                 zalbp = 0.472+(2.*(albice-0.472)*(hgbq-1.))     
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then          
                 zalbp = 0.2467+0.7049*hgbq                          
     &                  -(0.8608*(hgbq*hgbq))
     &                  +(0.3812*(hgbq*hgbq*hgbq)) 
          else                                                      
            zalbp = 0.1+3.6*hgbq
          endif                                                     
        endif                                                       
      endif                                                         
      zalb=zalbp+cgren                                           
      
      return
      end
      
      
