














c--parametres lies a la taille du domaine :
c [indispensable pour inclusion des fichiers bloc.com, reper.com, var??.com]
      integer imax,jmax,kmax,nsmax,nbpt,ijkmax,ixjmax
CL15  parameter ( imax = 122 , jmax = 65 , kmax = 15 )
      parameter ( imax = 122 , jmax = 65 , kmax = 20 )
CL30  parameter ( imax = 242 , jmax = 128 , kmax = 30 )

      parameter ( nsmax = 2 )
      parameter ( nbpt=imax*jmax )
      parameter ( ixjmax = imax*jmax , ijkmax = ixjmax*kmax )

c--Nombre maximum de coins (par niveaux et par type) :
      integer ncomax,nlpmax,nrpmax
CL15  parameter ( ncomax = 100 )
      parameter ( ncomax = 100 )
CL30  parameter ( ncomax = 400 )
c--Nombre maximum d'arretes (suivant X, Y, et pour les kmax Niv.)
c      (ordre de grandeur : imax*jmax)
CL15  parameter ( nlpmax = 6000 )
      parameter ( nlpmax = 6000 )
CL30  parameter ( nlpmax = 30000 )
c--Nombre maximum de points de grille avec Rap.Inter. (depend de "forc.corr"):
      parameter ( nrpmax = 10 )

c--determine le type de bassin : 0=ferme , 1=cyclique , 2=grilles raccordees
      integer ltest,jsepar
CL15  parameter ( ltest = 3 , jsepar = 50 )
      parameter ( ltest = 3 , jsepar = 50 )
CL30  parameter ( ltest = 3 , jsepar = 98 )

c--parametres lies a la frequence des donnees pour T , S , tau.x et tau.y
      integer nmois,nseas
      parameter ( nmois = 12 , nseas = 4 )
