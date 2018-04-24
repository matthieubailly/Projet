
!CONSTC - Détermine la constante C (dont la troncature sera fonction)
    !En entrée:
    !    DX   = longueurs des côtés
    !    COOR = coordonnées des sommets du carré
    !    XNC  = composante en x de la normale
    !    YNC  = composante en y de la normale
    !    V    = Volume de fluide
    !  En sortie:
    !    C    = Constante qu'on veut déterminer

      SUBROUTINE  CONSTC(C,DX,DY,V,COOR,XNC,YNC)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    ! Carré original
      DIMENSION COOR(NV,2)
      DOUBLE PRECISION M,M1
      CMIN=1.0D+14
      CMAX=-1.0D+14
      VT=DX*DY
      VBACK=V
      V=V/VT

      DO 10 I=1,4
         CI=-1.0*(COOR(I,1)*XNC+COOR(I,2)*YNC) ! Nombre de points du polygone
         IF(CI.LE.CMIN) THEN
            CMIN=CI
            IMIN=I
         END IF
         IF(CI.GE.CMAX) THEN
            CMAX=CI
            IMAX=I
         END IF
 10   CONTINUE

    ! Si la fraction solide/fluide est supérieure à 0.5, on inverse le problème

      IF((VBACK/VT).LE.0.5) THEN
         CI=CMIN
         I=IMIN
      ELSE
         CI=CMAX
         I=IMAX
         V=1.0-V
      END IF

    ! Normalisation

      SN=DABS(XNC)+DABS(YNC)
      XM=XNC/SN
      YM=YNC/SN
      XMI=XM*DX
      YMI=YM*DY
      SN=DABS(XMI)+DABS(YMI)
      XM=DABS(XMI)/SN
      YM=DABS(YMI)/SN

    ! Conditions aux bords

      M1=DMIN1(XM,YM)
      M=M1
      V1=M/(2.0*(1.0-M))

    ! Solution du problème inversé

      IF(V.GE.0.0.AND.V.LT.V1) THEN
         ALPHA=DSQRT(2.0*M*(1.0-M)*V)
      ELSE
         ALPHA=V*(1.0-M)+M/2.0
      END IF
      IF((VBACK/VT).LE.0.5) THEN
         C=CMIN+ALPHA*DABS(CMAX-CMIN)
      ELSE
         C=CMAX-ALPHA*DABS(CMAX-CMIN)
      END IF
      V=VBACK
      RETURN
      END


!------------------------------------
! CONFIGPOL - CONFIGURATION DU CARRÉ TRONQUÉ (devient alors un polygone)
    !En entrée:
    !     IPV0 = Tableau contenant les indices des sommets du carré original
    !     NTP0 = Indice du sommet
    !     NTV0 = Nombre total de sommets
    !     IA   = 0 si la normale est sortante, 1 sinon
    ! En sortie:
    !     XNCUT = Longueur unitaire des arêtes coupées selon x
    !     YNCUT = Longueur unitaire des arêtes coupées selon y
    !     IPV0  = Tableau contenant les indices des sommets du polygone tronqué
    !     NTP0  = Indice du sommet
    !     NTV0  = Nombre total de sommets
    !     IPIA0 = Indice du sommet de l'arête contenant l'intersection avec IA=0
    !     IPIA1 = Indice du sommet de l'arête contenant l'intersection avec IA=1


      SUBROUTINE CONFIGPOL(IA,IPIA0,IPIA1,IPV0,NTP0,NTV0,COOR0,XNCUT,YNCUT)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    ! Carré original
      DIMENSION IPV0(NV),COOR0(NV,2)

    ! Polygone modifié 1
      DIMENSION IPV1(NV)

    ! Variables
      DIMENSION IA(NV),IPIA0(2),IPIA1(2),XNCUT(2),YNCUT(2)

    ! On détermine les nouvelles arêtes (arêtes originales coupées)
      NTV1=0
      ICUT=0
      DO IV=1,NTV0
         IP=IPV0(IV)
         IV2=IV+1
         IF(IV.EQ.NTV0) IV2=1
         IP2=IPV0(IV2)
         IF(IA(IP).EQ.1) THEN
            NTV1=NTV1+1
            IPV1(NTV1)=IPV0(IV)
         END IF
         IF(IA(IP).NE.IA(IP2)) THEN
            ICUT=ICUT+1
            NTP0=NTP0+1
            NTV1=NTV1+1
            IPV1(NTV1)=NTP0
            IA(NTP0)=0
            IF(IA(IP2).EQ.0) THEN
               IPIA0(ICUT)=IP2
               IPIA1(ICUT)=IP
            ELSE
               IPIA0(ICUT)=IP
               IPIA1(ICUT)=IP2
            END IF
            XV=COOR0(IP2,1)-COOR0(IP,1)
            YV=COOR0(IP2,2)-COOR0(IP,2)
            RMOD=(XV**2.0+YV**2.0)**0.5
            XNCUT(ICUT)=YV/RMOD
            YNCUT(ICUT)=-XV/RMOD
         END IF
      END DO
      NTV0=NTV1
      DO IV=1,NTV1
         IPV0(IV)=IPV1(IV)
      END DO
      RETURN
      END


!-------------------------------------
!NOUVPOL - Sort le nouveau polygone
      !En entrée:
      !    COOR0 = Coordonnées des sommets du carré original
      !    IPV0  = Tableau contenant les indices des sommets du carré original
      !    NTP0  = Indice du sommet
      !    NTV0  = Nombre total de sommets
      !    IA    = 0 si la normale est sortante, 1 sinon
      !En sortie:
      !     IPV0   = Tableau contenant les indices des sommets du polygone tronqué
      !     NTP0   = Indice du sommet
      !     NTV0   = Nombre total de sommets
      !     ICONTN = Nombre de sommets du carré original hors de la zone tronquée
      !     ICONTP = Nombre de sommets du carré original dans la zone tronquée
      !     COOR0  = Coordonnées des sommets du nouveau polygone




      SUBROUTINE NOUVPOL(C,ICONTN,ICONTP,IPV0,NTP0,NTV0,COOR0,XNC,YNC)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    ! Carré original
      DIMENSION IPV0(NV),COOR0(NV,2)

    ! Variables
      DIMENSION IA(NV),IPIA0(2),IPIA1(2),PHIV(NV),XNCUT(2),YNCUT(2)
      ICONTP=0
      ICONTN=0
      TOLP=1.0D-14

    ! Valeurs de IA (0 ou 1)
      DO IV=1,NTV0
         IP=IPV0(IV)
         PHIV(IP)=XNC*COOR0(IP,1)+YNC*COOR0(IP,2)+C
         IF(PHIV(IP).GT.0.0) THEN
            IA(IP)=1
            ICONTP=ICONTP+1
         ELSE
            IA(IP)=0
            ICONTN=ICONTN+1
         END IF
      END DO
      IF(ICONTP.NE.0.AND.ICONTN.NE.0) THEN

    !Construction du nouveau polygone
         CALL CONFIGPOL(IA,IPIA0,IPIA1,IPV0,NTP0,NTV0,COOR0,XNCUT,YNCUT)

    ! Positions des nouveaux sommets
         DO IP=NTP0-1,NTP0
            IP0=IPIA0(IP-NTP0+2)
            IP1=IPIA1(IP-NTP0+2)
            IF(DABS(PHIV(IP1)-PHIV(IP0)).LT.TOLP) THEN
               COOR0(IP,1)=(COOR0(IP0,1)+COOR0(IP1,1))/2.0
               COOR0(IP,2)=(COOR0(IP0,2)+COOR0(IP1,2))/2.0
            ELSE
               COOR0(IP,1)=COOR0(IP0,1)-PHIV(IP0)*(COOR0(IP1,1)-COOR0(IP0,1))/(PHIV(IP1)-PHIV(IP0))
               COOR0(IP,2)=COOR0(IP0,2)-PHIV(IP0)*(COOR0(IP1,2)-COOR0(IP0,2))/(PHIV(IP1)-PHIV(IP0))
            END IF
         END DO
      END IF
      RETURN
      END

!------------------------------------
!SURFPOL - Donne la surface du nouveau polygone
    !En entrée:
    !    COOR  = Coordonnées des sommets du polygone
    !    IPV   = Tableau contenant les indices des sommets du polygone
    !    NTV   = Nombre total de sommets
    !En sortie:
    !    SURF  = Surface du polygone


      SUBROUTINE SURFPOL(IPV,NTV,COOR,SURF)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),COOR(NV,2)

      SUMS=0.0
      IH=INT((ntv-2)/2)
      DO I=2,IH+1
         IP=2*I
         IP1=IP-1
         IP2=IP-2
         XV1=COOR(IPV(IP1),1)-COOR(IPV(1),1)
         YV1=COOR(IPV(IP1),2)-COOR(IPV(1),2)
         XV2=COOR(IPV(IP),1)-COOR(IPV(IP2),1)
         YV2=COOR(IPV(IP),2)-COOR(IPV(IP2),2)
         SUMS=SUMS+XV1*YV2-YV1*XV2
      END DO
      IF(2*(IH+1).LT.ntv) THEN
         XV1=COOR(IPV(ntv),1)-COOR(IPV(1),1)
         YV1=COOR(IPV(ntv),2)-COOR(IPV(1),2)
         XV2=COOR(IPV(1),1)-COOR(IPV(ntv-1),1)
         YV2=COOR(IPV(1),2)-COOR(IPV(ntv-1),2)
         SUMS=SUMS+XV1*YV2-YV1*XV2
      ENDIF
      SURF=SUMS/2.0

      RETURN
      END

!------------------------------------
!CPPOL - Copie le polynome
      !En entrée:
      !    VERTP = Coordonnées des sommets du polygone qu'on veut copier
      !    IPV   = Tableau contenant les indices des sommets du polygone
      !    NTV   = Nombre total de sommets
      !    NTP   = Indice du sommet
      !En sortie:
      !    VERTP1 = Coordonnées des sommets de la copie
      !    IPV1   = Tableau contenant les indices des sommets de la copie
      !    NTV1   = Nombre total de sommets
      !    NTP1   = Indice du sommet


      SUBROUTINE CPPOL(IPV,IPV1,NTP,NTP1,NTV,NTV1,VERTP,VERTP1)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      DIMENSION IPV1(NV),VERTP1(NV,2)
      NTP1=NTP
      NTV1=NTV
      DO IV=1,NTV
         IP=IPV(IV)
         IPV1(IV)=IP
         VERTP1(IP,1)=VERTP(IP,1)
         VERTP1(IP,2)=VERTP(IP,2)
      END DO

      RETURN
      END

!------------------------------------
!RESTPOL - Restaure le polygone précédant la troncature
    !En entrée:
    !    COOR = Coordonnées des sommets du polygone
    !    IPV  = Tableau contenant les indices des sommets du polygone
    !    NTP  = Indice du sommet
    !    NTV  = Nombre total de sommets
    !En sortie:
    !    COOR = Coordonnées des sommets du polygone restauré
    !    IPV  = Tableau contenant les indices des sommets du polygone restauré
    !    NTP  = Indice du sommet
    !    NTV  = Nombre total de sommets


      SUBROUTINE RESTPOL(IPV,NTP,NTV,COOR)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    ! Polygone qu'on veut restaurer
      DIMENSION IPV(NV),COOR(NV,2)

    ! Polygone qu'on veut obtenir
      DIMENSION IPV0(NV),COOR0(NV,2)

    ! Obtention du polygone
      CALL CPPOL(IPV,IPV0,NTP,NTP0,NTV,NTV0,COOR,COOR0)

    ! Les points consécutifs qui ont le même vecteur position sont éliminés
    ! On utilise la tolérance TOLP
      TOLP=1.0D-16
      IVT=0
      DO IV=1,NTV0
         IP=IPV0(IV)
         IV0=IV-1
         IF(IV0.EQ.0) IV0=NTV0
         IP0=IPV0(IV0)
         DMOD=((COOR0(IP,1)-COOR0(IP0,1))**2.0+(COOR0(IP,2)-COOR0(IP0,2))**2.0)**0.5
         IF(DMOD.GT.TOLP) THEN
            IVT=IVT+1
            IPV(IVT)=IVT
            COOR(IVT,1)=COOR0(IP,1)
            COOR(IVT,2)=COOR0(IP,2)
         END IF
      END DO
      NTV=IVT
      NTP=IVT
      RETURN
      END

!------------------------------------
!DIST - Distance entre un point P et un segment
    !En entrée:
    !    X  = Coordonnées en x du segment
    !    Y  = Coordonnées en y du segment
    !    XP = Coordonnées en x du point P
    !    YP = Coordonnées en y du point P
    !En sortie:
    !    D  = Distance du point P au segment


    SUBROUTINE DIST(D,X,Y,XP,YP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    ! Segment
    DIMENSION X(2),Y(2)
    XNT=X(2)-X(1)
    YNT=Y(2)-Y(1)
    VMOD=(XNT**2.0+YNT**2.0)**0.5
    XNT=XNT/VMOD
    YNT=YNT/VMOD
    XN=-YNT
    YN=XNT
    C1=-1.0*(XNT*X(1)+YNT*Y(1))
    C2=1.0*(XNT*X(2)+YNT*Y(2))
    PHI1=XNT*XP+YNT*YP+C1
    PHI2=-XNT*XP-YNT*YP+C2
    IF(PHI1.GE.0.0.AND.PHI2.GE.0.0) THEN
       D=ABS(XN*X(1)+YN*Y(1)-(XN*XP+YN*YP))
    ELSEIF(PHI1.LE.0.0) THEN
       D=((XP-X(1))**2.0+(YP-Y(1))**2.0)**0.5
    ELSE
       D=((XP-X(2))**2.0+(YP-Y(2))**2.0)**0.5
    END IF
    RETURN
    END
