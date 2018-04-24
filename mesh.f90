! Subroutine du maillage carré

      SUBROUTINE SQUAREMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=4
      NTP=4
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=1.0
      VERTP(2,2)=0.0
      VERTP(3,1)=1.0
      VERTP(3,2)=1.0
      VERTP(4,1)=0.0
      VERTP(4,2)=1.0
      RETURN
      END
