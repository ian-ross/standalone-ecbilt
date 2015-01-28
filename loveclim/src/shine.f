












c
      subroutine shine(ih,zmue,tfsn,tfsg,ts,hgbq,
     &                 hnbq,zalb,zalcn,zalbp,zalcnp)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
c  Computes albedo of snow-sea ice following SHINE &
c  HENDERSSON-SELLERS [1985].
c
c  albin: Albedo of melting ice in the arctic.
c  albis: Albedo of melting ice in the antarctic (SHINE 
c         & HENDERSSON-SELLERS, 1985).
c  cgren: Correction of the snow or ice albedo to take into account
c         effects of cloudiness (GRENFELL & PEROVICH, 1984)
c
      include 'type.com'
c
C     albin = 0.53
C     albis = 0.53
C     cgren = 0.06
c  albice: Albedo of melting ice.
C     alphd  = 0.80
C     alphdi = 0.72
C     alphs  = 0.65
      albin = 0.45
      albis = 0.45
      albice = 0.45
      alphd  = 0.72
      alphdi = 0.64
      alphs  = 0.55
      cgren = 0.04
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1) Computation of surface albedo.                                   |
c-----------------------------------------------------------------------
c                                                                       
c  zalbp: Albedo for clear sky.
c  zalb:  Albedo for overcast sky.
c                                                                       
c                                                                       
c--1.1 Case of ice or snow.
c--------------------------
c                                                                       
      if (hnbq.gt.0.0) then                                       
c                                                                       
c  a) Case of ice covered by snow.
c                                                                       
        if (ts.lt.tfsn) then                                      
c                                                                       
c     Freezing snow.
c                                                                       
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
c     Melting snow.
c                                                                       
          if (hnbq.ge.0.1) then                                  
            zalbp = 0.65                                           
            zalbp = alphs                                      
          else                                                      
            zalbp = albice+((alphs-albice)/0.1)*hnbq
          endif                                                     
        endif                                                       
      else                                                          
c                                                                       
c  b) Case of ice free of snow.
c                                                                       
        if (ts.lt.tfsg) then                                      
c                                                                       
c     Freezing ice.
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
c     Melting ice.
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
c                                                                       
c--1.2. Case of the ocean.
c-------------------------
c                                                                       
      zalcnp=0.05/(1.1*zmue**1.4+0.15)                
c
c  Parameterization of BRIEGLED AND RAMANATHAN (1982)
c
      zalcn=0.06                                                 
c  see KONDRATYEV (1969) AND PAYNE (1972)
      return                                                            
      end                                                               
