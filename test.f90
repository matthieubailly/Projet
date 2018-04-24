
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! Carré original
      DIMENSION IPV(NV),VERTP(NV,2)

! Polygone
      DIMENSION IPV0(NV),VERTP0(NV,2)

! Segment
      DIMENSION X(2),Y(2)

! Maillage carré
      CALL SQUAREMESH(IPV,NTP,NTV,VERTP)

! Calcul de l'aire du polygone
      CALL SURFPOL(IPV,NTV,VERTP,VT)
      WRITE(6,*)'Aire du polygone est:',VT

! Fraction fluide/solide
      F=0.7

! Calcul de la surface de fluide dans le polygone
      V=F*VT

! Composantes du vecteur normal
      XNC=0.7
      YNC=0.7

! Détermination de C:
      DX=DABS(VERTP(1,1)-VERTP(2,1))    !Longueur en x
      DY=DABS(VERTP(2,2)-VERTP(3,2))    !Longueur en y
      CALL CONSTC(C,DX,DY,V,VERTP,XNC,YNC)
      WRITE(6,*)'On a ici C=',C

! Copie du polygone pour définir celui sur lequel on va travailler
      CALL CPPOL(IPV,IPV0,NTP,NTP0,NTV,NTV0,VERTP,VERTP0)

! Configuration du polygone
      CALL NOUVPOL(C,ICONTN,ICONTP,IPV0,NTP0,NTV0,VERTP0,XNC,YNC)

! Composantes du point P pour lequel on va calculer la distance avec l'interface
      XP=0
      YP=0

! Interface défini comme la dernière arête du polygone tronqué 0
      X(1)=VERTP0(NTP0-1,1)
      Y(1)=VERTP0(NTP0-1,2)
      X(2)=VERTP0(NTP0,1)
      Y(2)=VERTP0(NTP0,2)

! Calcul de la distance
      CALL DIST(D,X,Y,XP,YP)
      WRITE(6,*)'Distance de P à l''interface est D=:',D
      END
