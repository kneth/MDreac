CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RandGen is a collection of random number generators.                 C
C                                                                      C
C (C) Copyright 1995 by Kenneth Geisshirt (kneth@fatou.ruc.dk)         C
C Dept. of Life Sciences and Chemistry, Roskilde University            C
C                                                                      C
C History.                                                             C
C   9  Oct 1995: MarInit, MarRand, RandUnit                            C
C   11 Oct 1995: RandInt                                               C
C    1 Dec 1995: A few bug fixes.                                      C
C                                                                      C
C References:                                                          C
C [1] G. Marsaglia, Florida State University Report FSU-SCRI-87-50.    C
C [2] The Art of Computer Programming, volume 2. D.E. Knuth. Addison-  C
C     Wesley, 1969.                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MarInit is initialising the random number generator.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE MarInit(ij, kl)

      IMPLICIT NONE
      INCLUDE 'randgen.inc'

      INTEGER ij, kl, s, t
      INTEGER i, j, k, l, ii, jj, m

      IF ((ij .LT. 0) .OR. (ij .GT. 31328) 
     _.OR. (kl .LT. 0) .OR. (kl .GT. 30081)) THEN
         PRINT *, 'The first random number seed must have a value
     _between 0 and 31328.'
         PRINT *, 'The second seed must have a value between 
     _0 and 30081.'
         RETURN
      ELSE
         i=MOD(ij/177, 177)+2
         j=MOD(ij, 177)+2
         k=MOD(kl/169, 178)+1
         l=MOD(kl, 169)
         DO ii=1, 97
            s=0
            t=8388608
            DO jj=1, 24
               m=MOD(MOD(i*j, 179)*k, 179)
               i=j
               j=k
               k=m
               l=MOD(53*l+1, 169)
               IF (MOD(l*m, 64) .GE. 32) THEN
                  s=s+t
               ENDIF
               t=t/2
            ENDDO
            maru(ii)=s
         ENDDO
         marc=362436
         marcd=7654321
         marcm=16777213
         i97=97
         j97=33
      ENDIF

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MarRand is the actual random number generator, i.e. it returns a     C
C pseudo-random integer.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE MarRand(irandy)

      IMPLICIT NONE
      INCLUDE 'randgen.inc'

      INTEGER uni, irandy

      uni=maru(i97)-maru(j97)
      IF (uni .LT. 0) THEN
         uni=uni+16777216
      ENDIF
      maru(i97)=uni
      i97=i97-1
      IF (i97 .EQ. 0) THEN
         i97=97
      ENDIF
      j97=j97-1
      IF (j97 .EQ. 0) THEN
         j97=97
      ENDIF
      marc=marc-marcd
      IF (marc .LT. 0) THEN
         marc=marc+marcm
      ENDIF
      uni=uni-marc
      IF (uni .LT. 0) THEN
         uni=uni+16777216
      ENDIF
      irandy=uni

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RandUnit returns a random number from the unit interval.             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RandUnit(drandy)
      
      IMPLICIT NONE
      INCLUDE 'randgen.inc'
      
      REAL*8 marmax, imarmax
      PARAMETER (marmax=16777216.D0)
      PARAMETER (imarmax=1.D0/marmax)

      REAL*8 drandy
      INTEGER irandy

      CALL MarRand(irandy)
      drandy=FLOAT(irandy)*imarmax
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RandGauss is returning two gaussian distributed numbers (1, 0)       C
C The method implemented is the Box-Muller method.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RandGauss(x, y)

      IMPLICIT NONE
      REAL*8 x, y
      REAL*8 u1, u2, fac, randy, r

      CALL RandUnit(randy)
      u1=2.D0*randy-1.D0
      CALL RandUnit(randy)
      u2=2.D0*randy-1.D0
      r=u1**2+u2**2
 1003 IF ((r .LE. 1.D0) .AND. (r .NE. 0.D0)) GO TO 1002
      CALL RandUnit(randy)
      u1=2.D0*randy-1.D0
      CALL RandUnit(randy)
      u2=2.D0*randy-1.D0
      r=u1**2+u2**2
      GO TO 1003
 1002 fac=DSQRT(-2.D0*DLOG(r)/r)
      x=u1*fac
      y=u2*fac
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RandInt returns a random integer in the interval [1, max].           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RandInt(max, randy)

      IMPLICIT NONE
      INTEGER max, randy
    
      CALL MarRand(randy)
      randy=MOD(randy, max)+1
      END
