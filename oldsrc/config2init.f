CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Config2init converts a configuration to initial configuration.      C
C                                                                     C
C (C) Copyright 1997 by Kenneth Geisshirt (kneth@chem.ruc.dk)         C
C Dept. of Life Sciences and Chemistry, Roskilde University           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PROGRAM Final2Init

      IMPLICIT NONE

      REAL*8  x, y, vx, vy
      REAL*8  density, len2
      INTEGER num, stat, i, label

      PRINT *, 'Number of particles    : '
      READ (*, 1001) num
      PRINT *, 'Density                : '
      READ (*, 1002) density
      len2=0.5*DSQRT(DFLOAT(num)/density)

      CALL MarInit(1, 1)

      OPEN (11, FILE='conf.init', STATUS='NEW')
      OPEN (10, FILE='config', STATUS='OLD')
      DO i=1, num
         READ (10, 1010) x, y, label
         CALL RandGauss(vx, vy)
         x=x/len2
         y=y/len2
         WRITE (11, 1020) x, y, vx, vy, label
      ENDDO
 150  CLOSE (10, STATUS='KEEP')
      CLOSE (11, STATUS='KEEP')

 1001 FORMAT (I15)
 1002 FORMAT (F18.7)
 1010 FORMAT (F18.8, F18.8, I4)
 1020 FORMAT (F12.7, F12.7, F12.7, F12.7, I2)
      END

      INCLUDE 'randgen.f'
