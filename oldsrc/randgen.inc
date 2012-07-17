CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The include file for the random number generators.                   C
C                                                                      C
C (C) Copyright 1995 by Kenneth Geisshirt (kneth@fatou.ruc.dk)         C
C Dept. of Life Sciences and Chemistry, Roskilde University            C
C                                                                      C
C History.                                                             C
C   9 Oct 1995: basic structure                                        C
C                                                                      C
C References:                                                          C
C  See randgen.f for details                                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C Variables used by this library
      REAL*8 maru(98), marc, marcd, marcm
      INTEGER i97, j97
      
C Global variables
      COMMON /randgen/ maru, marc, marcd, marcm, i97, j97
