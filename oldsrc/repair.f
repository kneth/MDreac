CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Repair - put broken mdreac simulation together              C
C                                                             C
C (C) Copyright 1997 by Kenneth Geisshirt (kneth@chem.ruc.dk) C
C Depart. of Life Sciences and Chemistry, Roskilde University C
C                                                             C
C Last modified: 9 July 1997                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PROGRAM Repair

      IMPLICIT NONE

      INTEGER sfsites
      PARAMETER (sfsites = 128)
      
      INTEGER  startStep, totalStep, interval
      INTEGER  i, i1, i2, i3, i4, i5, i6, j
      REAL*8   d1, d2, d3, d4, d5, d6
      
      
      PRINT *, 'Broken step '
      READ (*, 3001) startStep
      PRINT *, 'Total step '
      READ (*, 3001) totalStep
      PRINT *, 'Dump interval '
      READ (*, 3001) interval


      OPEN (10, FILE='mdreac.conc',     STATUS='NEW')
      OPEN (21, FILE='mdreac.conc.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.conc.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1020) i1, d1, d2, d3, d4
         WRITE (10, 1020) i, d1, d2, d3, d4
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1020) i1, d1, d2, d4, d4
         WRITE (10, 1020) i, d1, d2, d3, d4
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.sample',     STATUS='NEW')
      OPEN (21, FILE='mdreac.sample.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.sample.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1010) i1, d1, d2
         WRITE (10, 1020) 10*i, d1, d2
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1020) i1, d1, d2
         WRITE (10, 1020) 10*i, d1, d2
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.energy',     STATUS='NEW')
      OPEN (21, FILE='mdreac.energy.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.energy.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1030) i1, d1, d2, d3
         WRITE (10, 1030) 10*i, d1, d2, d3
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1030) i1, d1, d2, d3
         WRITE (10, 1030) 10*i, d1, d2, d3
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.eta',     STATUS='NEW')
      OPEN (21, FILE='mdreac.eta.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.eta.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1070) i1, d1
         WRITE (10, 1070) 10*i, d1
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1070) i1, d1
         WRITE (10, 1070) 10*i, d1
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')

      OPEN (10, FILE='mdreac.diff',     STATUS='NEW')
      OPEN (21, FILE='mdreac.diff.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.diff.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1070) i1, d1
         WRITE (10, 1070) 10*i, d1
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1070) i1, d1
         WRITE (10, 1070) 10*i, d1
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')

      OPEN (10, FILE='mdreac.reac',     STATUS='NEW')
      OPEN (21, FILE='mdreac.reac.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.reac.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1040) i1, i2, i3, i4, i5, i6
         WRITE (10, 1040) 10*i, i2, i3, i4, i5, i6
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1040) i1, i2, i3, i4, i5, i6
         WRITE (10, 1040) 10*i, i2, i3, i4, i5, i6
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.subconc',     STATUS='NEW')
      OPEN (21, FILE='mdreac.subconc.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.subconc.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1020) i1, d1, d2, d3, d4
         WRITE (10, 1020) 10*i, d1, d2, d3, d4
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1020) i1, d1, d2, d3, d4
         WRITE (10, 1020) 10*i, d1, d2, d3, d4
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.enereac',     STATUS='NEW')
      OPEN (21, FILE='mdreac.enereac.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.enereac.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1060) i1, d1, d2, d3, d4, d5, d6
         WRITE (10, 1060) 10*i, d1, d2, d3, d4, d5, d6
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1060) i1, d1, d2, d3, d4, d5, d6
         WRITE (10, 1060) 10*i, d1, d2, d3, d4, d5, d6
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.naborr',     STATUS='NEW')
      OPEN (21, FILE='mdreac.naborr.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.naborr.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1030) i1, d1, d2, d3
         WRITE (10, 1030) 10*i, d1, d2, d3
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1030) i1, d1, d2, d3
         WRITE (10, 1030) 10*i, d1, d2, d3
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.nabor',     STATUS='NEW')
      OPEN (21, FILE='mdreac.nabor.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.nabor.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1030) i1, d1, d2, d3
         WRITE (10, 1030) 10*i, d1, d2, d3
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1030) i1, d1, d2, d3
         WRITE (10, 1030) 10*i, d1, d2, d3
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.order',     STATUS='NEW')
      OPEN (21, FILE='mdreac.order.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.order.2nd', STATUS='OLD')
      DO i=1, startStep
         READ (21, 1030) i1, d1
         WRITE (10, 1030) 10*i, d1
      ENDDO
      DO i=startStep+1, totalStep
         READ (22, 1030) i1, d1
         WRITE (10, 1030) 10*i, d1
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.sf',     STATUS='NEW')
      OPEN (21, FILE='mdreac.sf.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.sf.2nd', STATUS='OLD')
      DO i=25000, startStep, 25000
         READ (21, 3001) i1
         WRITE (10, 3001) i1
         DO j=0, 31
            READ (21, 2010) d1, d2
            WRITE (10, 2010) d1, d2
         ENDDO
      ENDDO
      DO i=startStep+25000, totalStep, 25000
         READ (22, 3001) i1
         WRITE (10, 3001) i1+10*startStep
         DO j=0, 31
            READ (22, 2010) d1, d2
            WRITE (10, 2010) d1, d2
         ENDDO
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


      OPEN (10, FILE='mdreac.cor',     STATUS='NEW')
      OPEN (21, FILE='mdreac.cor.1st', STATUS='OLD')
      OPEN (22, FILE='mdreac.cor.2nd', STATUS='OLD')
      DO i=25000, startStep, 25000
         READ (21, 3001) i1
         WRITE (10, 3001) i1
         DO j=0, 31
            READ (21, 2010) d1, d2
            WRITE (10, 2010) d1, d2
         ENDDO
      ENDDO
      DO i=startStep+25000, totalStep, 25000
         READ (22, 3001) i1
         WRITE (10, 3001) i1+10*startStep
         DO j=0, 31
            READ (22, 2010) d1, d2
            WRITE (10, 2010) d1, d2
         ENDDO
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
      CLOSE (22, STATUS='KEEP')


 1010 FORMAT (I15, F18.7, F18.7)
 1020 FORMAT (I15, F18.7, F18.7, F18.7, F18.7)
 1030 FORMAT (I15, F18.7, F18.7, F18.7)
 1040 FORMAT (I15, I8, I8, I8, I8, I8)
 1050 FORMAT (I15, I8, I8, I8, I8, I8)
 1060 FORMAT (I15, F18.7, F18.7, F18.7, F18.7, F18.7, F18.7)
 1070 FORMAT (I15, F18.7)
 2010 FORMAT (F18.8, F18.8)
 3001 FORMAT (I15)
      END
      



