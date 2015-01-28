












      subroutine ocesla
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  Programme simulant un ocean statique agissant uniquement comme un reservoir 
c  de chaleur. Ce programme n'est utilise que pour la mise au point!
c  modif : 01/04/94
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
C     include 'ice.com'
c        write(iuo+66,*)'avant',scal(imax-2,5,ks2,1),scal(imax-2,5,ks2,2)
c      &                    ,phiss(imax-2,5,1),phiss(imax-2,5,2)
      do 30 k=1,nsmax
         do 20 j=js1,js2
            do 10 i=is1(j),is2(j)
               scal(i,j,ks2,k)=scal(i,j,ks2,k) - phiss(i,j,k)
C              scal(i,j,ks2,k)=scal(i,j,ks2,k)
C    &                         -dts(ks2)*unsdz(ks2)*phiss(i,j,k)
 10         continue
 20      continue
 30   continue
c        write(iuo+66,*)'temp apres',scal(imax-2,5,ks2,1),
c    &              scal(imax-2,5,ks2,2)
      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine slaboc -
      end
