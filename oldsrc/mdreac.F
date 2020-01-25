CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Molecular Dynamics program for reactive Lennard-Jones mixtures.   C
C                                                                   C
C (C) Copyright 1995-2020 by                                        C
C Kenneth Geisshirt <geisshirt@gmail.com>                           C
C                                                                   C
C History:                                                          C
C   5 Oct 1995: main program, ComputeForces, NoseHoover,            C
C               UpdateLinkCell, OneStep, ManyStep, ApplyPerBoundary,C
C               Links.                                              C
C   6 Oct 1995: ReadInit, WriteFinal, Initialise, PrintData.        C
C   7 Oct 1995: correcting.                                         C
C   9 Oct 1995: bug fixing, scaling of initial configuration added, C
C               DeInitialise.                                       C
C  11 Oct 1995: bug fixing.                                         C
C  12 Oct 1995: bug fixing, DoReactions, pressure calculation,      C
C               extended PrintData.                                 C
C  13 Oct 1995: optimisation of the linked cells, tests.            C
C  16 Oct 1995: ComputeConc.                                        C
C  18 Oct 1995: PrintData extended; work on PVM version begun       C
C               InitPVM, changes in main program, SortPart.         C
C  19 Oct 1995: parallellisation of ComputeForces, OneStep.         C
C  20 Oct 1995: PrintData altered so time is the real-world time    C
C               parallellisation of WriteFinal, main program.       C
C  23 Oct 1995: ReadParameters, a bit restructuring so many parame- C
C               ters are read from a file.                          C
C  26 Oct 1995: parallellisation.                                   C
C  27 Oct 1995: bug fixing the PVM version.                         C
C   1 Nov 1995: DumpPos.                                            C
C   2 Nov 1995: DumpPos extended in order to write the gnuplot file C
C               automatically.                                      C
C   3 Nov 1995: DumpPos changed so it write to one file only, code  C
C               for a report file written.                          C
C   6 Nov 1995: bug fixing the PVM version.                         C
C   7 Nov 1995: bug fixing the PVM version, and redesign of the sub-C
C               cell structure for the parallel version.            C
C  15 Nov 1995: printing of real time changed in PrintData, compu-  C
C               tation of the structure factor added.               C
C  16 Nov 1995: correcting structure factor code so it fits this    C
C               MD program.                                         C
C  17 Nov 1995: fixing a bug in DumpPos, parallellisation (SendPos, C
C               RecvPos and changes in ManyStep and DumpPos).       C
C  19 Nov 1995: bug fixing the PVM version.                         C
C  21 Nov 1995: bug fixing the PVM version, and a little            C
C               redesigning (only one temperature, and added as few C
C               constants).                                         C
C  22 Nov 1995: bug fixing the PVM version.                         C
C  24 Nov 1995: bug fixing the PVM version.                         C
C  27 Nov 1995: fixing a bug in the structure factor code, imple-   C
C               mented constant energy (microcanonical) simulation, C
C               bug fixing the PVM version.                         C
C  28 Nov 1995: all COMMON blocks moved to a file called            C
C               mdreac.common, a few small modifications, bug fixingC
C               the PVM version.                                    C
C  29 Nov 1995: bug fixing the PVM version.                         C
C   1 Dec 1995: fixing a few bugs in the structure factor code to-  C
C               gether with Mohamed.                                C
C   4 Dec 1995: fixing a few bugs.                                  C
C  12 Dec 1995: bug fixing the PVM version.                         C
C  13 Dec 1995: bug fixing the PVM version.                         C
C  14 Dec 1995: bug fixing the PVM version.                         C
C  15 Dec 1995: bug fixing the PVM version.                         C
C  18 Dec 1995: bug fixing the PVM version.                         C
C  22 Dec 1995: DoReactionsPVM                                      C
C   4 Jan 1996: extended the program so it also contains a program  C
C               which generates initial configurations, i.e.        C
C               GenInit, redesigned ComputeForces, modified         C
C               WriteFinal a bit.                                   C
C   9 Jan 1996: Optimisation of the code, i.e. reducing the         C
C               communication between processors in the PVM version C
C               and applying misc. tricks. Timing code added.       C
C  10 Jan 1996: code for countering number of reactions occured add.C
C  29 Mar 1996: Fine-tuning the code for IBM RS/6000 (Power2),      C
C               based on advices in [6].                            C
C  23 Apr 1996: Computing energy of reacting particles.             C
C  29 Apr 1996: Code for calculating the diffusion coefficient add. C
C               InitDiffusion, MeanSqrDispl, modification of main,  C
C               ApplyPerBoundary, Initilise, PrintData.             C
C  13 May 1996: DoReactions modified so it dumps the positions if   C
C               a reaction has occured.                             C
C  21 May 1996: Code for counting number of neighbours added.       C
C  30 May 1996: Extended the code above.                            C
C   2 Nov 1996: Extended reaction scheme.                           C
C   3 Dec 1996: Code cleanup.                                       C
C   9 Dec 1996: More code cleanup.                                  C
C  17 Jan 1997: Reactions are placed in a seperate file which makes C
C               it easier to change mechanism. Some code cleanup    C
C               done as well.                                       C
C  31 Jan 1997: The macro FIVE_SPECIES chooses a five-variable      C
C               system.                                             C
C   5 Feb 1997: Fixed a minor bug in GenInit (induced by the exten- C
C               sion from 31 Jan 1997).                             C
C  18 Mar 1997: Printing concentrations in subdomains.              C
C   2 May 1997: Two-thermostat option implemented. Not parallelised.C
C  20 May 1997: NumToStr extended a bit.                            C
C  22 May 1997: NumToStr2 coded.                                    C
C  27 May 1997: Number of reaction parameters increased. Minor bug  C
C               in WriteFinal corrected.                            C
C   1 Jul 1997: Reservoir-model implemented.                        C
C   1 Aug 1997: Verlet list implemented.                            C
C  21 Jun 1999: Dimerisation code added.                            C
C   5 Jul 1999: Cont dimerisation code.                             C
C  30 Dec 2019: Build with a modern Fortran compiler (gcc v8)       C
C                                                                   C
C References:                                                       C
C  [1]  S. Toxvaerd, Molec. Phys., 72(1), 159-168, 1991.            C
C  [2]  J.J. Morales and S. Toxvaerd, Comput. Phys. Commun., 71,    C
C       71-76, 1992.                                                C
C  [3]  D. Brown, J.H.R. Clarke, M. Okuda, and T. Yamazaki, Comput. C
C       Phys. Commun., 74, 67-80, 1993.                             C
C  [4]  PVM: Parallel Virtual Machine. A Users Guide and Tutorial   C
C       for Networked Parallel Computing. A. Geist et al.,          C
C       MIT Press, 1994.                                            C
C  [5]  Numerical Recipes in Fortran, 2nd edition.                  C
C       W.H. Press et al. Cambridge University Press, 1992.         C
C  [6]  High Performance Computing.                                 C
C       K. Dowd. OReilly & Associates, 1993.                        C
C  [7]  Theory of Simple Liquids, 2nd edition.                      C
C       J.P. Hansen and I.R. McDonald. Academic Press, 1986.        C
C  [8]  Algorithms.                                                 C
C       R. Sedgework. Addison-Wesley, 1983.                         C
C  [9]  A.Z. Panagiotopoulos. Molec. Phys. 61(4), 813-826, 1987.    C
C                                                                   C
C Remarks:                                                          C
C  The program simulates a three-component system. The three compo- C
C  nents are called A (label=0), X (label=1) and Y (label=2). The   C
C  reactions implemented are: A+X->2X, X+Y->2Y, Y+A->2A.            C
C  The program comes in two versions: a sequentail and a parallel.  C
C  Which one which is compiled is controlled by the macro PVM_VER,  C
C  i.e. if it is defined then the parallel version is compiled.     C
C  The parallellisation is based on the ideas in [3]. The program-  C
C  ming model is message-parsing based on PVM [4].                  C
C  It is possible to compile the program with two different ways of C
C  implementing the periodic boundary conditions. It is controlled  C
C  by the macro PBC_LIST.                                           C
C  The program for generating initial configurations is compiled if C
C  the macro GEN_INIT is defined.                                   C
C  Timing of the program is enabled by defining the macro TIMING.   C
C  In order to compile the program, it must be preprocessed by cpp. C
C  The macro PHOTO_REACT chooses the the photoinduced scheme, other-C
C  wise it is the simple extended Lotka which is simulated.         C
C  The macro RESERVOIR controls the use of a reservoir coupled to   C
C  reactor. Idea comes from [9].                                    C
C  The macro VERLET controls the use of Verlet list. If the macro   C
C  is undefined, the linked-cell structure is used.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#undef _ONLY_CONC_

#undef  DEBUG
#undef  TIMING
#undef  DO_SORT
#undef  PHOTO_REACT
#undef FIVE_SPECIES
#undef TWO_THERMO
#undef MANY_REAC_PARAMS

#undef RESERVOIR

#ifdef PVM_VER

C Messages parsed during run-time
#define MSG_INIT    1
#define MSG_SORT    2
#define MSG_FPOS    3
#define MSG_FFOR    4
#define MSG_STEP_A  5
#define MSG_PRNT    6
#define MSG_END     7
#define MSG_FINAL   8
#define MSG_STEP_B  9
#define MSG_POSS   10

C Macros used in the PVM version
#ifdef _PLATFORM_SP2_
#  define PACK_TYPE PvmRaw
#else
#  define PACK_TYPE PvmDataDefault
#endif
#define INTEGER_4 4
#define REAL_8  6

C name of executable - may not be too long
#ifdef GEN_INIT
#define EXECNAME "pargeninit"
#else
#define EXECNAME "parmdreac"
#endif
#endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Main program                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PROGRAM mdreac

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
#ifdef PVM_VER
      INCLUDE 'fpvm3.h'
      INTEGER bufid, info
      REAL*4  itime
#endif

      INTEGER i

C The MAIN program

#ifdef DEBUG
      PRINT *, 'Main program has just started'
#endif

#ifdef TIMING
C Start timing
      calctime=etime(timearray)
      commtime=0.0
#endif

#ifdef PVM_VER
      CALL InitPVM
#else
      CALL ReadParameters

#ifdef GEN_INIT
      CALL GenInit
#else
      CALL ReadInit
#endif

      CALL Initialise
      CALL InitDiffusion
#endif

C The main loops
      DO i=1, step
         CALL ManyStep(i)
      ENDDO
      CALL WriteFinal

#ifdef PVM_VER
C Wait until all children have finished
      IF (myparent .EQ. PvmNoParent) THEN
#ifdef TIMING
         itime=etime(timearray)
#endif
         DO i=1, ntask-1
            CALL pvmfrecv(alltids(i), MSG_END, bufid)
            CALL pvmffreebuf(bufid, info)
         ENDDO
#ifdef TIMING
         itime=itime-etime(timearray)
         commtime=commtime+itime
#endif
         CALL DeInitialise
         CALL pvmfexit(info)
      ELSE
         CALL pvmfinitsend(PACK_TYPE, bufid)
         CALL pvmfsend(myparent, MSG_END, info)
         CALL pvmfexit(info)
      ENDIF
#else
      CALL DeInitialise
#endif

#ifdef DEBUG
      PRINT *, 'Main program has finished'
#endif
      END
C end of main program



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NoseHoover is an implementation of the Nose-Hoover thermostat.    C
C The implementation here is based on algorithm 2 in [1].           C
C The constant dt48 is 48*dt; this is a small optimisation trick in C
C order to reduce the number of loops in the ComputeForce routine.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NoseHoover

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i
#ifdef TWO_THERMO
      REAL*8  K(0:2), eta1(0:2), eta2(0:2), sum(0:2)
      INTEGER lab
#else /* TWO_THERMO */
      REAL*8  K, eta1, eta2, sum
#endif /* TWO_THERMO */

#ifdef MICRO_CANON
      eta=0.D0
      eta1=1.D0
      eta2=1.D0
#else /* MICRO_CANON */
#ifdef TWO_THERMO
      DO i=0, 2
         K(i)=0.5D0*dt*eta(i)
         eta1(i)=1.D0-K(i)
         eta2(i)=1.D0/(1.D0+K(i))
      ENDDO
      eta(1)=0.D0
      eta1(1)=1.D0
      eta2(1)=1.D0
#else /* TWO_THERMO */
      K=0.5D0*dt*eta
      eta1=1.D0-K
      eta2=1.D0/(1.D0+K)
#endif /* TWO_THERMO */
#endif /* MICRO_CANON */

C Loop is splitted - SP2 has four-way set ass. cache
#ifdef TWO_THERMO
      DO i=0, 2
         sum(i)=0.D0
      ENDDO
      DO i=1, num
         lab=label(i)
         vx(i)=(vx(i)*eta1(lab)+dt48*fx(i))*eta2(lab)
         vy(i)=(vy(i)*eta1(lab)+dt48*fy(i))*eta2(lab)
      ENDDO
      DO i=1, num
         lab=label(i)
         qx(i)=qx(i)+vx(i)*dt
         qy(i)=qy(i)+vy(i)*dt
         sum(lab)=sum(lab)+vx(i)**2+vy(i)**2
      ENDDO
      Ekin2=sum(0)+sum(1)+sum(2)
      eta(0)=eta(0)+(sum(0)-2.0*DFLOAT(na)*tempA)*invQ*dt
      eta(2)=eta(2)+(sum(2)-2.0*DFLOAT(ny)*tempC)*invQ*dt
#else /* TWO_THERMO */
      sum=0.D0
      DO i=1, num
         vx(i)=(vx(i)*eta1+dt48*fx(i))*eta2
         vy(i)=(vy(i)*eta1+dt48*fy(i))*eta2
      ENDDO
      DO i=1, num
         qx(i)=qx(i)+vx(i)*dt
         qy(i)=qy(i)+vy(i)*dt
         sum=sum+vx(i)**2+vy(i)**2
      ENDDO
#ifdef RESERVOIR
      DO i=1, num
         vxII(i)=(vxII(i)*eta1+dt48*fxII(i))*eta2
         vyII(i)=(vyII(i)*eta1+dt48*fyII(i))*eta2
      ENDDO
      DO i=1, numII
         qxII(i)=qxII(i)+vxII(i)*dt
         qyII(i)=qyII(i)+vyII(i)*dt
         sum=sum+vxII(i)**2+vyII(i)**2
      ENDDO
#endif /* RESERVOIR */
      Ekin2=sum
#ifndef PVM_VER
      eta=eta+(sum-dof*temp)*invQ*dt
#endif /* PVM_VER */
#endif /* TWO_THERMO */

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C UpdateLinkCell finds all particle pairs which may interact in the C
C near future. The implementation is based upon [2].                C
C The routine is exteded in order to handle chemical kinetics.      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE UpdateLinkCell

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, j, icell, cell0, jcell
      INTEGER ik1, ik2, nabor
      REAL*8  r2, rx, ry, rxi, ryi
#ifdef PBC_LIST
      REAL*8  xdiff, ydiff
#endif /* PBC_LIST */
      INTEGER ilabel, jlabel, k

      ik1=0
      ik2=0
      DO icell=1, Mcell2
         i=head(icell)
 2010    IF (i .NE. 0) THEN
            rxi=qx(i)
            ryi=qy(i)
            ilabel=label(i)
            j=list(i)
 2020       IF (j .NE. 0) THEN
               rx=rxi-qx(j)
               ry=ryi-qy(j)
#ifdef PBC_LIST
               IF (rx .GE. Len2) THEN
                  xdiff=-Len
               ELSEIF (rx .LT. Len2m) THEN
                  xdiff=Len
               ELSE
                  xdiff=0.D0
               ENDIF
               rx=rx+xdiff
#ifndef PVM_VER
               IF (ry .GE. Len2) THEN
                  ydiff=-Len
               ELSEIF (ry .LT. Len2m) THEN
                  ydiff=Len
               ELSE
                  ydiff=0.D0
               ENDIF
               ry=ry+ydiff
#endif /* !PVM_VER */
#else /* PBC_LIST */
               rx=rx-Len*DNINT(rx*iLen)
#ifdef PVM_VER
               ry=ry-Len*DNINT(ry*iLen)
#endif /* PVM_VER */
#endif /* PBC_LIST */
               r2=rx**2+ry**2
               IF (r2 .LT. bufzone2) THEN
                  IF (ilabel .EQ. label(j)) THEN
                     ik1=ik1+1
                     ipair1(ik1)=i
                     jpair1(ik1)=j
#ifdef PBC_LIST
                     xper1(ik1)=xdiff
#ifdef PVM_VER
                     yper1(ik1)=0.D0
#else /* PVM_VER */
                     yper1(ik1)=ydiff
#endif /* PVM_VER */
#endif /* PBC_LIST */
                  ELSE
                     ik2=ik2+1
                     ipair2(ik2)=i
                     jpair2(ik2)=j
#ifdef PBC_LIST
                     xper2(ik2)=xdiff
#ifdef PVM_VER
                     yper2(ik2)=0.D0
#else /* PVM_VER */
                     yper2(ik2)=ydiff
#endif /* PVM_VER */
#endif /* PBC_LIST */
                  ENDIF
               ENDIF
               j=list(j)
               GOTO 2020
            ENDIF

C Neighbouring cells
            cell0=4*(icell-1)
            DO nabor=1, 4
               jcell=map(cell0+nabor)
               j=head(jcell)
 2030          IF (j .NE. 0) THEN
                  rx=rxi-qx(j)
                  ry=ryi-qy(j)
#ifdef PBC_LIST
                  IF (rx .GE. Len2) THEN
                     xdiff=-Len
                  ELSEIF (rx .LT. Len2m) THEN
                     xdiff=Len
                  ELSE
                     xdiff=0.D0
                  ENDIF
                  rx=rx+xdiff
#ifndef PVM_VER
                  IF (ry .GE. Len2) THEN
                     ydiff=-Len
                  ELSEIF (ry .LT. Len2m) THEN
                     ydiff=Len
                  ELSE
                     ydiff=0.D0
                  ENDIF
                  ry=ry+ydiff
#endif /* !PVM_VER */
#else /* PBC_LIST */
                  rx=rx-Len*DNINT(rx*iLen)
#ifdef PVM_VER
                  ry=ry-Len*DNINT(ry*iLen)
#endif /* PVM_VER */
#endif /* PBC_LIST */
                  r2=rx**2+ry**2
                  IF (r2 .LT. bufzone2) THEN
                     IF (ilabel .EQ. label(j)) THEN
                        ik1=ik1+1
                        ipair1(ik1)=i
                        jpair1(ik1)=j
#ifdef PBC_LIST
                        xper1(ik1)=xdiff
#ifndef PVM_VER
                        yper1(ik1)=ydiff
#else /* PVM_VER */
                        yper1(ik1)=0.D0
#endif /* PVM_VER */
#endif /* PBC_LIST */
                     ELSE
                        ik2=ik2+1
                        ipair2(ik2)=i
                        jpair2(ik2)=j
#ifdef PBC_LIST
                        xper2(ik2)=xdiff
#ifndef PVM_VER
                        yper2(ik2)=ydiff
#else /* !PVM_VER */
                        yper2(ik2)=0.D0
#endif /* !PVM_VER */
#endif /* PBC_LIST */
                     ENDIF
                  ENDIF
                  j=list(j)
                  GOTO 2030
               ENDIF
            ENDDO
            i=list(i)
            GOTO 2010
         ENDIF
      ENDDO
      numbpair1=ik1
      numbpair2=ik2

#ifndef PVM_VER
      DO i=1, 3
         DO j=1, 3
            nabortable(i, j)=0.D0
         ENDDO
      ENDDO
      DO k=1, numbpair1
         i=ipair1(k)
         j=jpair1(k)
         ilabel=label(i)+1
         jlabel=label(j)+1
#ifdef PBC_LIST
         rx=qx(i)-qx(j)+xper1(k)
         ry=qy(i)-qy(j)+yper1(k)
#else /* PBC_LIST */
         rx=qx(i)-qx(j)
         ry=qy(i)-qy(j)
         rx=rx-Len*DNINT(rx*iLen)
         ry=ry-Len*DNINT(ry*iLen)
#endif /* PBC_LIST */
         r2=rx**2+ry**2
         IF (r2 .LE. R2nabor) THEN
            nabortable(ilabel, jlabel)=nabortable(ilabel, jlabel)+1.0
            nabortable(jlabel, ilabel)=nabortable(jlabel, ilabel)+1.0
         ENDIF
      ENDDO
      DO k=1, numbpair2
         i=ipair2(k)
         j=jpair2(k)
         ilabel=label(i)+1
         jlabel=label(j)+1
#ifdef PBC_LIST
         rx=qx(i)-qx(j)+xper1(k)
         ry=qy(i)-qy(j)+yper1(k)
#else /* PBC_LIST */
         rx=qx(i)-qx(j)
         ry=qy(i)-qy(j)
         rx=rx-Len*DNINT(rx*iLen)
         ry=ry-Len*DNINT(ry*iLen)
#endif /* PBC_LIST */
         r2=rx**2+ry**2
         IF (r2 .LE. R2nabor) THEN
            nabortable(ilabel, jlabel)=nabortable(ilabel, jlabel)+1.0
            nabortable(jlabel, ilabel)=nabortable(jlabel, ilabel)+1.0
         ENDIF
      ENDDO
#endif /* PVM_VER */
      END


#ifdef PVM_VER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CommPart checks which particles should be sent to neighbouring    C
C cpus, and recieves particles from the other cpus.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CommPart

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
      INCLUDE 'fpvm3.h'

      INTEGER i1, i2, i3, j, i, ii
      INTEGER uplist(ntot), downlist(ntot), keeplist(ntot)
      REAL*8  tmpvec1(ntot), tmpvec2(ntot)
      REAL*8  tmpvec3(ntot), tmpvec4(ntot)
      INTEGER iarray(ntot), bufid, info, numb(1)
#ifdef TIMING
      REAL*4  itime
#endif

#ifdef DEBUG
      PRINT *, 'BEGIN of CommPart'
#endif

#ifdef TIMING
      itime=etime(timearray)
#endif
C Check which to be sent
      i1=0
      i2=0
      i3=0
      DO j=1, num
         IF (qy(j) .GE. ymax) THEN
            i1=i1+1
            uplist(i1)=j
         ELSEIF (qy(j) .LT. ymin) THEN
            i2=i2+1
            downlist(i2)=j
         ELSE
            i3=i3+1
            keeplist(i3)=j
         ENDIF
      ENDDO

C Send up
      numb(1)=i1
      DO i=1, i1
         ii=uplist(i)
         IF (qy(ii) .GE. Len2) THEN
            qy(ii)=qy(ii)-Len
         ELSEIF (qy(ii) .LT. Len2m) THEN
            qy(ii)=qy(ii)+Len
         ENDIf
         tmpvec1(i)=qx(ii)
         tmpvec2(i)=qy(ii)
         tmpvec3(i)=vx(ii)
         tmpvec4(i)=vy(ii)
         iarray(i)=label(ii)
      ENDDO
      CALL pvmfinitsend(PACK_TYPE, bufid)
      CALL pvmfpack(INTEGER_4, numb, 1, 1, info)
      CALL pvmfpack(INTEGER_4, iarray, i1, 1, info)
      CALL pvmfpack(REAL_8, tmpvec1, i1, 1, info)
      CALL pvmfpack(REAL_8, tmpvec2, i1, 1, info)
      CALL pvmfpack(REAL_8, tmpvec3, i1, 1, info)
      CALL pvmfpack(REAL_8, tmpvec4, i1, 1, info)
      CALL pvmfsend(CPUup, MSG_SORT, info)
      CALL pvmffreebuf(bufid, info)
#ifdef DEBUG
      PRINT *, mytid, 'sending SORT to ', CPUup
#endif

C Send down
      numb(1)=i2
      DO i=1, i2
         ii=downlist(i)
         IF (qy(ii) .GE. Len2) THEN
            qy(ii)=qy(ii)-Len
         ELSEIF (qy(ii) .LT. Len2m) THEN
            qy(ii)=qy(ii)+Len
         ENDIf
         tmpvec1(i)=qx(ii)
         tmpvec2(i)=qy(ii)
         tmpvec3(i)=vx(ii)
         tmpvec4(i)=vy(ii)
         iarray(i)=label(ii)
      ENDDO
      CALL pvmfinitsend(PACK_TYPE, bufid)
      CALL pvmfpack(INTEGER_4, numb, 1, 1, info)
      CALL pvmfpack(INTEGER_4, iarray, i2, 1, info)
      CALL pvmfpack(REAL_8, tmpvec1, i2, 1, info)
      CALL pvmfpack(REAL_8, tmpvec2, i2, 1, info)
      CALL pvmfpack(REAL_8, tmpvec3, i2, 1, info)
      CALL pvmfpack(REAL_8, tmpvec4, i2, 1, info)
      CALL pvmfsend(CPUdown, MSG_SORT, info)
      CALL pvmffreebuf(bufid, info)
#ifdef DEBUG
      PRINT *, mytid, 'sending SORT to ', CPUdown
#endif

C Keep particles
      i=1
      DO j=1, i3
         ii=keeplist(j)
         qx(i)=qx(ii)
         qy(i)=qy(ii)
         vx(i)=vx(ii)
         vy(i)=vy(ii)
         label(i)=label(ii)
         i=i+1
      ENDDO

C Recieve up
#ifdef DEBUG
      PRINT *, mytid, 'waiting SORT from ', CPUup
#endif
      CALL pvmfrecv(CPUup, MSG_SORT, bufid)
#ifdef DEBUG
      PRINT *, mytid, 'recieved SORT from ', CPUup
#endif
      CALL pvmfunpack(INTEGER_4, numb, 1, 1, info)
      i1=numb(1)
      IF (i1 .NE. 0) THEN
         CALL pvmfunpack(INTEGER_4, label(i3+1), i1, 1, info)
         CALL pvmfunpack(REAL_8, qx(i3+1), i1, 1, info)
         CALL pvmfunpack(REAL_8, qy(i3+1), i1, 1, info)
         CALL pvmfunpack(REAL_8, vx(i3+1), i1, 1, info)
         CALL pvmfunpack(REAL_8, vy(i3+1), i1, 1, info)
      ENDIF
      CALL pvmffreebuf(bufid, info)

C Recieve down
#ifdef DEBUG
      PRINT *, mytid, 'waiting SORT from ', CPUdown
#endif
      CALL pvmfrecv(CPUdown, MSG_SORT, bufid)
#ifdef DEBUG
      PRINT *, mytid, 'recieved SORT from ', CPUdown
#endif
      CALL pvmfunpack(INTEGER_4, numb, 1, 1, info)
      i2=numb(1)
      IF (i2 .NE. 0) THEN
         CALL pvmfunpack(INTEGER_4, label(i3+i1+1), i2, 1, info)
         CALL pvmfunpack(REAL_8, qx(i3+i1+1), i2, 1, info)
         CALL pvmfunpack(REAL_8, qy(i3+i1+1), i2, 1, info)
         CALL pvmfunpack(REAL_8, vx(i3+i1+1), i2, 1, info)
         CALL pvmfunpack(REAL_8, vy(i3+i1+1), i2, 1, info)
      ENDIF
      CALL pvmffreebuf(bufid, info)
      num=i1+i2+i3

#ifdef TIMING
      itime=etime(timearray)-itime
      commtime=commtime+itime
#endif

#ifdef DEBUG
      PRINT *, mytid, 'has now ', num
      PRINT *, 'END of CommPart'
#endif
      END
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Links sorts the particles, i.e. it finds which subcell which      C
C particle belongs to.                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE Links

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, cll
      REAL*8  drandy
#ifdef TIMING
      REAL*4  itime
#endif

#ifdef DEBUG
      PRINT *, 'BEGIN of Links'
#endif

#ifdef PVM_VER
      CALL CommPart
#endif

      DO i=1, Mcell2
         head(i)=0
      ENDDO

      DO i=1, num
         cll=1+INT((qx(i)+Len2)*celli)
     _        +INT((qy(i)+Len2)*celli)*Mcell
         list(i)=head(cll)
         head(cll)=i
      ENDDO

#ifdef PHOTO_REACT
      nr5=0
      DO i=1, num
         IF (label(i) .EQ. 3) THEN
            CALL RandUnit(drandy)
            IF (drandy .LE. prob3) THEN
               label(i)=0
               nr5=nr5+1
            ENDIF
         ENDIF
      ENDDO
#endif

#ifdef PVM_VER
C The particles in the lowest row will always be send
C in ComputeForces
#ifdef TIMING
      itime=etime(timearray)
#endif
      nforce=0
      DO cll=firstcell, firstcell+Mcell-1
         i=head(cll)
 10      IF (i .NE. 0) THEN
            nforce=nforce+1
            forcelist(nforce)=i
            i=list(i)
            GOTO 10
         ENDIF
      ENDDO
#ifdef TIMING
      itime=etime(timearray)-itime
      commtime=commtime+itime
#endif /* TIMING */
#ifdef DEBUG
      PRINT *, mytid, 'presorts particles for ComputeForces:', nforce
#endif /* DEBUG */
#endif /* PVM_VER */

#ifdef DEBUG
      PRINT *, 'END of Links'
#endif /* DEBUG */
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ComputeForces is calculating the force acting on the particles.   C
C Two particles of the same type interact with the LJ potential     C
C while two particles of different types interact with the WCA      C
C potential.                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE ComputeForces(iter)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
#ifdef PVM_VER
      INCLUDE 'fpvm3.h'
#endif

      INTEGER iter

      INTEGER i, j, k
      REAL*8  rxij, ryij, r2, pot, r6, rrr
#ifdef TIMING
      REAL*8  itime
#endif
#ifdef DEBUG
      INTEGER term1, term2
#endif

#ifdef PVM_VER
C Forces and position from neighbouring cpus
      REAL*8  xarray(ntot), yarray(ntot)
      REAL*8  rxi, ryi, fxi, fyi
      INTEGER ilabel, bufid, jp, jc, jx, numb(1)
      INTEGER icell, ipp, i1, info, plabel(ntot)
      INTEGER iarray(ntot)
#endif


#ifdef DEBUG
      term1=0
      term2=0
#endif

C Init. the forces
      DO i=1, ntot
         fx(i)=0.D0
         fy(i)=0.D0
      ENDDO

      Epot=0.D0
      press=0.D0

#ifdef PVM_VER
C sending positions down
#ifdef TIMING
      itime=etime(timearray)
#endif
      DO i=1, nforce
         j=forcelist(i)
         xarray(i)=qx(j)
         yarray(i)=qy(j)
         plabel(i)=label(j)
      ENDDO
      numb(1)=nforce
      CALL pvmfinitsend(PACK_TYPE, bufid)
      CALL pvmfpack(INTEGER_4, numb, 1, 1, info)
      CALL pvmfpack(INTEGER_4, plabel, nforce, 1, info)
      CALL pvmfpack(REAL_8, xarray, nforce, 1, info)
      CALL pvmfpack(REAL_8, yarray, nforce, 1, info)
      CALL pvmfsend(CPUdown, MSG_FPOS, info)
#ifdef DEBUG
      PRINT *, mytid, 'sending FPOS to ', CPUdown
#endif
      CALL pvmffreebuf(bufid, info)

C recieving positions
#ifdef DEBUG
      PRINT *, mytid, 'waiting FPOS from ', CPUup
#endif
      CALL pvmfrecv(CPUup, MSG_FPOS, bufid)
#ifdef DEBUG
      PRINT *, mytid, 'recieved FPOS from ', CPUup
#endif
      CALL pvmfunpack(INTEGER_4, numb, 1, 1, info)
      i1=numb(1)
      IF (i1 .NE. 0) THEN
         CALL pvmfunpack(INTEGER_4, label(num+1), i1, 1, info)
         CALL pvmfunpack(REAL_8, qx(num+1), i1, 1, info)
         CALL pvmfunpack(REAL_8, qy(num+1), i1, 1, info)
      ENDIF
      CALL pvmffreebuf(bufid, info)

#ifdef TIMING
      itime=etime(timearray)-itime
      commtime=commtime+itime
#endif
#endif


#ifdef VERLET
      DO i=1, num-1
         jbeg=list(i)
         jend=list(i+1)-1
         IF (jbeg .LE. jend) THEN
            rxi=qx(i)
            ryi=qy(i)
            fxi=fx(i)
            fyi=fy(i)
            ilab=label(i)
            DO k=jbeg, jend
               j=list(k)
               rxij=rxi-rx(j)
               ryij=ryi-ry(j)
               rxij=rxij-Len*DNINT(rxij*iLen)
               ryij=ryij-Len*DNINT(ryij*iLen)
               r2=rxij**2+ryij**2
               IF (ilab .EQ. label(j)) THEN
                  IF (r2 .LT. r2cut1) THEN
                     r2=1.D0/r2
                     r6=r2*r2*r2
                     pot=r2*r6*(r6-0.5D0)
                     Epot=Epot+r6*(r6-1.D0)
                     press=press+r6*(r6-0.5D0)
                     fxi=fxi+rxij*pot
                     fyi=fyi+ryij*pot
                     fx(j)=fx(j)-rxij*pot
                     fy(j)=fy(j)-ryij*pot
                  ENDIF
               ELSE
                  IF (r2 .LT. r2cut2) THEN
                     r2=1.D0/r2
                     r6=r2*r2*r2
                     pot=r2*r6*(r6-0.5D0)
                     Epot=Epot+r6*(r6-1.D0)
                     press=press+r6*(r6-0.5D0)
                     fxi=fxi+rxij*pot
                     fyi=fyi+ryij*pot
                     fx(j)=fx(j)-rxij*pot
                     fy(j)=fy(j)-ryij*pot
                  ENDIF
               ENDIF
            ENDDO
            fx(i)=fxi
            fy(i)=fyi
         ENDIF
      ENDDO
#else /* not VERLET */
C Loop over even pairs
      DO k=1, numbpair1
         i=ipair1(k)
         j=jpair1(k)
#ifdef PBC_LIST
         rxij=qx(i)-qx(j)+xper1(k)
         ryij=qy(i)-qy(j)+yper1(k)
#else
         rxij=qx(i)-qx(j)
         ryij=qy(i)-qy(j)
         rxij=rxij-Len*DNINT(rxij*iLen)
         ryij=ryij-Len*DNINT(ryij*iLen)
#endif
         r2=rxij**2+ryij**2

C Which potential to use
#ifdef GEN_INIT
         IF (r2 .LT. r2cut1) THEN
            r2=1.D0/r2
            r6=r2*r2*r2
            pot=r2*r6*(r6-0.5D0)
            Epot=Epot+r6*(r6-1.D0)
            press=press+r6*(r6-0.5D0)
            fx(i)=fx(i)+rxij*pot
            fy(i)=fy(i)+ryij*pot
            fx(j)=fx(j)-rxij*pot
            fy(j)=fy(j)-ryij*pot
#ifdef DEBUG
            term1=term1+1
#endif
         ENDIF
#else
         IF (r2 .LT. r2cut1) THEN
            r2=1.D0/r2
            r6=r2*r2*r2
            pot=r2*r6*(r6-0.5D0)
            Epot=Epot+r6*(r6-1.D0)
            press=press+r6*(r6-0.5D0)
            fx(i)=fx(i)+rxij*pot
            fy(i)=fy(i)+ryij*pot
            fx(j)=fx(j)-rxij*pot
            fy(j)=fy(j)-ryij*pot
#ifdef DEBUG
            term1=term1+1
#endif
         ENDIF
#endif
      ENDDO

C Loop over odd pairs
      DO k=1, numbpair2
         i=ipair2(k)
         j=jpair2(k)
#ifdef PBC_LIST
         rxij=qx(i)-qx(j)+xper2(k)
         ryij=qy(i)-qy(j)+yper2(k)
#else
         rxij=qx(i)-qx(j)
         ryij=qy(i)-qy(j)
         rxij=rxij-Len*DNINT(rxij*iLen)
         ryij=ryij-Len*DNINT(ryij*iLen)
#endif
         r2=rxij**2+ryij**2

C Which potential to use
#ifdef GEN_INIT
         IF (r2 .LT. r2cut1) THEN
            r2=1.D0/r2
            r6=r2*r2*r2
            pot=r2*r6*(r6-0.5D0)
            Epot=Epot+r6*(r6-1.D0)
            press=press+r6*(r6-0.5D0)
            fx(i)=fx(i)+rxij*pot
            fy(i)=fy(i)+ryij*pot
            fx(j)=fx(j)-rxij*pot
            fy(j)=fy(j)-ryij*pot
#ifdef DEBUG
            term1=term1+1
#endif
         ENDIF
#else
         IF (r2 .LT. r2cut2) THEN
            r2=1.D0/r2
            r6=r2*r2*r2
            pot=r2*r6*(r6-0.5D0)
            Epot=Epot+r6*(r6-1.D0)
            press=press+r6*(r6-0.5D0)
            fx(i)=fx(i)+rxij*pot
            fy(i)=fy(i)+ryij*pot
            fx(j)=fx(j)-rxij*pot
            fy(j)=fy(j)-ryij*pot
#ifdef DEBUG
            term1=term1+1
#endif
         ENDIF
#endif
      ENDDO

#endif /* not VERLET */


#ifdef PVM_VER
#ifdef TIMING
      itime=etime(timearray)
#endif

#ifdef DEBUG
      PRINT *, 'Parallel Force begins at', mytid
      PRINT *, 'There are', i1, 'particle at', mytid
#endif
C sort particles at neighbouring cpu
      DO i=1, firstcell-1
         head(i)=0
      ENDDO
      DO i=lastcell+1, Mcell2
         head(i)=0
      ENDDO

      DO i=num+1, num+i1
         icell=1+INT((qx(i)+Len2)*celli)
         IF (icell .LE. 0) THEN
            icell=icell+Mcell
         ELSEIF (icell .GT. Mcell) THEN
            icell=icell-Mcell
         ENDIF
         icell=icell+lastcell
         IF (icell .GT. Mcell2) THEN
            icell=icell-Mcell2
         ELSEIF (icell .LE. 0) THEN
            icell=icell+Mcell2
         ENDIF
         list(i)=head(icell)
         head(icell)=i
      ENDDO
#ifdef TIMING
      itime=etime(timearray)-itime
      commtime=commtime+itime
#endif

#ifndef GEN_INIT
      IF (iter .EQ. 1) THEN
         CALL DoReactionsPVM
      ENDIF
#endif

C now compute the forces
      DO icell=lastcell-Mcell+1, lastcell
         ipp=head(icell)
 3020    IF (ipp .NE. 0) THEN
            rxi=qx(ipp)
            ryi=qy(ipp)
            ilabel=label(ipp)
            fxi=fx(ipp)
            fyi=fy(ipp)
            DO jc=1, 3
               jx=map(4*(icell-1)+jc)
               jp=head(jx)
 3030          IF (jp .NE. 0) THEN
                  rxij=rxi-qx(jp)
                  IF (rxij .GE. Len2) THEN
                     rxij=rxij-Len
                  ELSEIF (rxij .LE. Len2m) THEN
                     rxij=rxij+Len
                  ENDIF
                  ryij=ryi-qy(jp)
                  IF (ryij .GE. Len2) THEN
                     ryij=ryij-Len
                  ELSEIF (ryij .LE. Len2m) THEN
                     ryij=ryij+Len
                  ENDIF
                  r2=rxij**2+ryij**2
#ifdef GEN_INIT
                  IF (r2 .LT. r2cut1) THEN
                     r2=1.D0/r2
                     r6=r2*r2*r2
                     pot=r2*r6*(r6-0.5D0)
                     Epot=Epot+r6*(r6-1.D0)
                     press=press+r6*(r6-0.5D0)
                     fxi=fxi+rxij*pot
                     fyi=fyi+ryij*pot
                     fx(jp)=fx(jp)-rxij*pot
                     fy(jp)=fy(jp)-ryij*pot
#ifdef DEBUG
                     term2=term2+1
#endif
                  ENDIF
#else
                  IF (ilabel .EQ. label(jp)) THEN
                     IF (r2 .LT. r2cut1) THEN
                        r2=1.D0/r2
                        r6=r2*r2*r2
                        pot=r2*r6*(r6-0.5D0)
                        Epot=Epot+r6*(r6-1.D0)
                        press=press+r6*(r6-0.5D0)
                        fxi=fxi+rxij*pot
                        fyi=fyi+ryij*pot
                        fx(jp)=fx(jp)-rxij*pot
                        fy(jp)=fy(jp)-ryij*pot
#ifdef DEBUG
                        term2=term2+1
#endif
                     ENDIF
                  ELSE
                     IF (r2 .LT. r2cut2) THEN
                        r2=1.D0/r2
                        r6=r2*r2*r2
                        pot=r2*r6*(r6-0.5D0)
                        Epot=Epot+r6*(r6-1.D0)
                        press=press+r6*(r6-0.5D0)
                        fxi=fxi+rxij*pot
                        fyi=fyi+ryij*pot
                        fx(jp)=fx(jp)-rxij*pot
                        fy(jp)=fy(jp)-ryij*pot
#ifdef DEBUG
                        term2=term2+1
#endif
                     ENDIF
                  ENDIF
#endif
                  jp=list(jp)
                  GOTO 3030
               ENDIF
            ENDDO
            fx(ipp)=fxi
            fy(ipp)=fyi
            ipp=list(ipp)
            GOTO 3020
         ENDIF
      ENDDO

#ifdef DEBUG
      PRINT *, 'Parallel Force ends at', mytid
      PRINT *, 'Seq. terms:', term1, '|| terms:', term2
#endif

C send forces back
#ifdef TIMING
      itime=etime(timearray)
#endif
      CALL pvmfinitsend(PACK_TYPE, bufid)
      CALL pvmfpack(REAL_8, fx(num+1), i1, 1, info)
      CALL pvmfpack(REAL_8, fy(num+1), i1, 1, info)
      IF (iter .EQ. 1) THEN
         CALL pvmfpack(INTEGER_4, label(num+1), i1, 1, info)
      ENDIF
      CALL pvmfsend(CPUup, MSG_FFOR, info)
      CALL pvmffreebuf(bufid, info)
#ifdef DEBUG
      PRINT *, mytid, 'sending FFOR to ', CPUup
#endif

C recieve forces
#ifdef DEBUG
      PRINT *, mytid, 'waiting FFOR from', CPUdown
#endif
      CALL pvmfrecv(CPUdown, MSG_FFOR, bufid)
#ifdef DEBUG
      PRINT *, mytid, 'recieved FFOR from', CPUdown
#endif
      IF (nforce .NE. 0) THEN
         CALL pvmfunpack(REAL_8, xarray, nforce, 1, info)
         CALL pvmfunpack(REAL_8, yarray, nforce, 1, info)
         IF (iter .EQ. 1) THEN
            CALL pvmfunpack(INTEGER_4, iarray, nforce, 1, info)
            DO i=1, nforce
               j=forcelist(i)
               label(j)=iarray(i)
            ENDDO
         ENDIF
         DO i=1, nforce
            j=forcelist(i)
            fx(j)=fx(j)+xarray(i)
            fy(j)=fy(j)+yarray(i)
         ENDDO
      ENDIF
      CALL pvmffreebuf(bufid, info)
#ifdef TIMING
      itime=etime(timearray)-itime
      commtime=commtime+itime
#endif
#endif

      Epot=4.D0*Epot
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C OneStep calls the subroutines to make one time step.              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OneStep(iter)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER iter

#ifdef PVM_VER
      INCLUDE 'fpvm3.h'
      REAL*8  tmpvec(3)
      INTEGER iarray(9), bufid, info
#endif

      INTEGER i
#ifdef TIMING
      REAL*4  itime
#endif

#ifdef DEBUG
      PRINT *, 'BEGIN of OneStep'
#endif

      CALL ApplyPerBoundary
      CALL ComputeForces(iter)
#ifdef RESERVOIR
      CALL ComputeForcesII
#endif /* RESERVOIR */
      CALL NoseHoover

#ifndef MICRO_CANON
#ifdef PVM_VER
      IF (myparent .EQ. PvmNoParent) THEN
         CALL CountSpecies
#ifdef TIMING
         itime=etime(timearray)
#endif
         DO i=1, ntask-1
#ifdef DEBUG
            PRINT *, mytid, 'waiting for MSG_STEP_A from', alltids(i)
#endif
            CALL pvmfrecv(alltids(i), MSG_STEP_A, bufid)
#ifdef DEBUG
            PRINT *, mytid, 'recieved MSG_STEP_A from', alltids(i)
#endif
            CALL pvmfunpack(REAL_8, tmpvec, 3, 1, info)
            CALL pvmfunpack(INTEGER_4, iarray, 9, 1, info)
            CALL pvmffreebuf(bufid, info)
            Ekin2=Ekin2+tmpvec(1)
            press=press+tmpvec(2)
            Epot=Epot+tmpvec(3)
            na=na+iarray(1)
            nx=nx+iarray(2)
            ny=ny+iarray(3)
            nz=nz+iarray(4)
            IF (iter .EQ. 1) THEN
               nr1=nr1+iarray(5)
               nr2=nr2+iarray(6)
               nr3=nr3+iarray(7)
               nr4=nr4+iarray(8)
               nr5=nr5+iarray(9)
            ENDIF
         ENDDO
         eta=eta+(Ekin2-dof*temp)*invQ*dt
         Ekin=0.5D0*Ekin2
         CALL pvmfinitsend(PACK_TYPE, bufid)
         tmpvec(1)=eta
         CALL pvmfpack(REAL_8, tmpvec, 1, 1, info)
         CALL pvmfmcast(ntask-1, alltids, MSG_STEP_B, info)
         CALL pvmffreebuf(bufid, info)
#ifdef TIMING
         itime=etime(timearray)-itime
         commtime=commtime+itime
#endif
      ELSE
         tmpvec(1)=Ekin2
         tmpvec(2)=press
         tmpvec(3)=Epot
         CALL CountSpecies
         iarray(1)=na
         iarray(2)=nx
         iarray(3)=ny
         iarray(4)=nz
         iarray(5)=nr1
         iarray(6)=nr2
         iarray(7)=nr3
         iarray(8)=nr4
         iarray(9)=nr5
         CALL pvmfinitsend(PACK_TYPE, bufid)
         CALL pvmfpack(REAL_8, tmpvec, 3, 1, info)
         CALL pvmfpack(INTEGER_4, iarray, 9, 1, info)
         CALL pvmfsend(myparent, MSG_STEP_A, info)
         CALL pvmffreebuf(bufid, info)
         CALL pvmfrecv(myparent, MSG_STEP_B, bufid)
         CALL pvmfunpack(REAL_8, tmpvec, 1, 1, info)
         eta=tmpvec(1)
         CALL pvmffreebuf(bufid, info)
      ENDIF
#endif
#endif
      Ekin=0.5D0*Ekin2

#ifdef DEBUG
      PRINT *, 'END of OneStep'
#endif
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ManyStep updates the linked-cell structure, makes a number of     C
C time steps, and writes the output to various files.               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ManyStep(start)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
#ifdef PVM_VER
      INCLUDE 'fpvm3.h'
#endif

      INTEGER start, i, time
      REAL*8  sfunc(0:256), corel(0:256), qq, rr

#ifdef DEBUG
      PRINT *, 'BEGIN of ManyStep'
      PRINT *, 'start = ', start
#endif

      CALL ApplyPerBoundary
#ifdef RESERVOIR
      DO i=1, numExchange
         CALL Exchange
      ENDDO
      CALL LinksII
      CALL UpdateLinkCellII
#endif /* RESERVOIR */
#ifdef VERLET
      CALL MakeVerlet
#else /* VERLET */
      CALL Links
      CALL UpdateLinkCell
#endif /* VERLET */
#ifndef GEN_INIT
      CALL DoReactions(start*NLupdate)
#endif
      DO i=1, NLupdate
         CALL OneStep(i)
      ENDDO
      time=start*NLupdate
      CALL PrintData(time)
#ifndef GEN_INIT
      IF (MOD(time, PrnPoss) .EQ. 0) THEN
#ifdef PVM_VER
         IF (myparent .EQ. PvmNoParent) THEN
            CALL RecvPos
            CALL DumpPos(time)
            CALL stfac(sfunc, corel)
            WRITE (18, 1020) time
            WRITE (20, 1020) time
            DO 10 i=0, sfsites/2-1
               qq=DFLOAT(i)*2.D0*PI*iLen
               rr=DFLOAT(i)*Len/DFLOAT(sfsites)
               WRITE (18, 1010) qq, sfunc(i)
               WRITE (20, 1010) rr, corel(i)
 10         CONTINUE
         ELSE
            CALL SendPos
         ENDIF
#else
         CALL DumpPos(time)
         CALL stfac(sfunc, corel)
         WRITE (18, 1020) time
         WRITE (20, 1020) time
         DO 10 i=0, 31
            qq=DFLOAT(i)*2.D0*PI*iLen
            rr=DFLOAT(i)*Len/64.D0
            WRITE (18, 1010) qq, sfunc(i)
            WRITE (20, 1010) rr, corel(i)
 10      CONTINUE
#endif
      ENDIF
#endif

#ifdef DEBUG
      PRINT *, 'END of ManyStep'
#endif

 1010 FORMAT (F18.8, F18.8)
 1020 FORMAT (I15)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ComputeConc is computing the concentration of the three species.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ComputeConc(xa, xx, xy, xz)

      IMPLICIT NONE

      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, number(4)
      REAL*8  xa, xx, xy, xz

      number(1)=0
      number(2)=0
      number(3)=0
      number(4)=0
      DO i=1, num
         number(label(i)+1)=number(label(i)+1)+1
      ENDDO
      xa=DFLOAT(number(1))/DFLOAT(num)
      xx=DFLOAT(number(2))/DFLOAT(num)
      xy=DFLOAT(number(3))/DFLOAT(num)
      xz=DFLOAT(number(4))/DFLOAT(num)

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CountSpecies is finding how many particles there are of each      C
C species.                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CountSpecies

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i

      na=0
      nx=0
      ny=0
      nz=0
      DO i=1, num
         IF (label(i) .EQ. 0) THEN
            na=na+1
         ELSEIF (label(i) .EQ. 1) THEN
            nx=nx+1
         ELSEIF (label(i) .EQ. 2) THEN
            ny=ny+1
         ELSEIF (label(i) .EQ. 3) THEN
            nz=nz+1
         ENDIF
      ENDDO
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PrintData is writing data to the files.                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PrintData(time)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
#ifdef PVM_VER
      INCLUDE 'fpvm3.h'
#endif

      REAL*8  numbA, numbX, numbY, numbZ
      INTEGER time, i
      REAL*8  temperature, pressure, volume

#ifdef PVM_VER
      IF (myparent .EQ. PvmNoParent) THEN
         temperature=2.D0*Ekin/dof
         volume=Len**2
         pressure=rho*temperature+24.D0*rho*press/volume
     _        -3.D0*PI*rho**2*(1.D0/r2cut1)
         Ekin=Ekin/DFLOAT(nall)
         Epot=Epot/DFLOAT(nall)
         WRITE (17, 1030) time, Epot, Ekin, Epot+Ekin
         WRITE (15, 1010) time, temperature, pressure
         numbA=DFLOAT(na)/DFLOAT(nall)
         numbX=DFLOAT(nx)/DFLOAT(nall)
         numbY=DFLOAT(ny)/DFLOAT(nall)
         numbZ=DFLOAT(nz)/DFLOAT(nall)
         WRITE (16, 1020) time, numbA, numbX, numbY, numbZ
         WRITE (21, 1040) time, nr1, nr2, nr3, nr4, nr5
      ENDIF
#else /* NOT PVM_VER */
#ifndef _ONLY_CONC_
      volume=Len**2
      temperature=2.D0*Ekin/dof
      pressure=rho*temperature+24.D0*rho*press/volume
     _-3.D0*PI*rho**2*(1.D0/r2cut1)
      Ekin=Ekin/DFLOAT(nall)
      Epot=Epot/DFLOAT(nall)
      WRITE (17, 1030) time, Epot, Ekin, Epot+Ekin
      WRITE (15, 1010) time, temperature, pressure
#ifdef TWO_THERMO
      WRITE (27, 1030) time, eta(0), eta(1), eta(2)
#else /* TWO_THERMO */
      WRITE (27, 1070) time, eta
#endif /* TWO_THERMO */
#endif /* _ONLY_CONC_ */
      CALL ComputeConc(numbA, numbX, numbY, numbZ)
      CALL DomainConc(time)
      WRITE (16, 1020) time, numbA, numbX, numbY, numbZ
#ifndef _ONLY_CONC_
      WRITE (21, 1040) time, nr1, nr2, nr3, nr4, nr5
      IF (nr1 .NE. 0) THEN
         R1pot=R1pot/DFLOAT(nr1)
         R1kin=R1kin/DFLOAT(nr1)
      ENDIF
      IF (nr2 .NE. 0) THEN
         R2pot=R2pot/DFLOAT(nr2)
         R2kin=R2kin/DFLOAT(nr2)
      ENDIF
      IF (nr3 .NE. 0) THEN
         R3kin=R3kin/DFLOAT(nr3)
         R3pot=R3pot/DFLOAT(nr3)
      ENDIF
      WRITE (22, 1060) time, R1pot, R1kin, R2pot, R2kin, R3pot, R3kin
      CALL MeanSqrDispl
      WRITE (19, 1070) time, meansqr

      CALL CountSpecies
      DO i=1, 3
         nabortable(1, i)=nabortable(1, i)/DFLOAT(na)
         nabortable(2, i)=nabortable(2, i)/DFLOAT(nx)
         nabortable(3, i)=nabortable(3, i)/DFLOAT(ny)
      ENDDO
      WRITE (23, 1030) time, nabortable(1, 1), nabortable(1, 2),
     _     nabortable(1, 3)
      WRITE (24, 1030) time, naborreac(1, 1), naborreac(1, 2),
     _     naborreac(1, 3)
      WRITE (25, 1070) time, meanAA
#endif /* _ONLY_CONC_ */

#endif /* NOT PVM_VER */

 1010 FORMAT (I15, F18.7, F18.7)
 1020 FORMAT (I15, F18.7, F18.7, F18.7, F18.7)
 1030 FORMAT (I15, F18.7, F18.7, F18.7)
 1040 FORMAT (I15, I8, I8, I8, I8, I8)
 1050 FORMAT (I15, F18.7, F18.7)
 1060 FORMAT (I15, F18.7, F18.7, F18.7, F18.7, F18.7, F18.7)
 1070 FORMAT (I15, F18.7)
 1080 FORMAT (I15, I5, I5, I5)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ApplyPerBoundary applies the periodic boundary conditions to all  C
C particles in the system.                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ApplyPerBoundary

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i

      DO i=1, num
         IF (qx(i) .GE. Len2) THEN
            qx(i)=qx(i)-Len
            nboxx(i)=nboxx(i)+1
         ELSEIF (qx(i) .LT. Len2m) THEN
            qx(i)=qx(i)+Len
            nboxx(i)=nboxx(i)-1
         ENDIF
#ifndef PVM_VER
         IF (qy(i) .GE. Len2) THEN
            qy(i)=qy(i)-Len
            nboxy(i)=nboxy(i)+1
         ELSEIF (qy(i) .LT. Len2m) THEN
            qy(i)=qy(i)+Len
            nboxy(i)=nboxy(i)-1
         ENDIF
#endif /* PVM_VER */
      ENDDO

#ifdef RESERVOIR
      DO i=1, numII
         IF (qxII(i) .GE. Len2II) THEN
            qxII(i)=qxII(i)-LenII
         ELSEIF (qxII(i) .LT. Len2mII) THEN
            qxII(i)=qxII(i)+LenII
         ENDIF
         IF (qyII(i) .GE. Len2II) THEN
            qyII(i)=qyII(i)-LenII
         ELSEIF (qyII(i) .LT. Len2mII) THEN
            qyII(i)=qyII(i)+LenII
         ENDIF
      ENDDO
#endif /* RESERVOIR */
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The reactions are include at this point. The source file must at  C
C least has one subroutine called DoReaction.                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "extlotka.reac"


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C WriteFinal writes the final positions and velocities to a file.   C
C The routine also writes other physical quantuties.                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE WriteFinal

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
#ifdef PVM_VER
      INCLUDE 'fpvm3.h'
#endif

      INTEGER i
#ifdef PVM_VER
      INTEGER j, iarray(1), bufid, info
#ifdef TIMING
      REAL*4  itime
#endif
#endif


#ifdef PVM_VER
      IF (myparent .EQ. PvmNoParent) THEN
#ifdef GEN_INIT
         OPEN (10, FILE='conf.init', STATUS='NEW')
         DO i=1, num
            qx(i)=qx(i)/Len2
            qy(i)=qy(i)/Len2
            WRITE (10, 1005) qx(i), qy(i), vx(i), vy(i), label(i)
         ENDDO
#else
         OPEN (10, FILE='conf.final.A', STATUS='NEW')
         OPEN (11, FILE='conf.final.X', STATUS='NEW')
         OPEN (12, FILE='conf.final.Y', STATUS='NEW')
         OPEN (13, FILE='conf.final.Z', STATUS='NEW')
         DO i=1, num
            IF (label(i) .EQ. 0) THEN
               WRITE (10, 1004) qx(i), qy(i), vx(i), vy(i)
            ELSEIF (label(i) .EQ. 1) THEN
               WRITE (11, 1004) qx(i), qy(i), vx(i), vy(i)
            ELSEIF (label(i) .EQ. 2) THEN
               WRITE (12, 1004) qx(i), qy(i), vx(i), vy(i)
            ELSEIF (label(i) .EQ. 3) THEN
               WRITE (13, 1004) qx(i), qy(i), vx(i), vy(i)
            ENDIF
         ENDDO
         DO i=1, ntask-1
#ifdef TIMING
            itime=etime(timearray)
#endif
            CALL pvmfrecv(alltids(i), MSG_FINAL, bufid)
            CALL pvmfunpack(INTEGER_4, iarray, 1, 1, info)
            num=iarray(1)
            IF (num .GT. 0) THEN
               CALL pvmfunpack(REAL_8, qx, num, 1, info)
               CALL pvmfunpack(REAL_8, qy, num, 1, info)
               CALL pvmfunpack(REAL_8, vx, num, 1, info)
               CALL pvmfunpack(REAL_8, vy, num, 1, info)
               CALL pvmfunpack(INTEGER_4, label, num, 1, info)
            ENDIF
            CALL pvmffreebuf(bufid, info)
#ifdef TIMING
            itime=etime(timearray)-itime
            commtime=commtime+itime
#endif
#endif /* GEN_INIT */
#ifdef GEN_INIT
            DO j=1, num
               qx(i)=qx(i)/Len2
               qy(i)=qy(i)/Len2
               WRITE (10, 1005) qx(i), qy(i), vx(i), vy(i), label(i)
            ENDDO
#else /* not GEN_INIT */
             DO j=1, num
               IF (label(j) .EQ. 0) THEN
                  WRITE (10, 1004) qx(j), qy(j), vx(j), vy(j)
               ELSEIF (label(j) .EQ. 1) THEN
                  WRITE (11, 1004) qx(j), qy(j), vx(j), vy(j)
               ELSEIF (label(j) .EQ. 2) THEN
                  WRITE (12, 1004) qx(j), qy(j), vx(j), vy(j)
               ELSEIF (label(j) .EQ. 3) THEN
                  WRITE (13, 1004) qx(j), qy(j), vx(j), vy(j)
               ENDIF
            ENDDO
         ENDDO
         CLOSE (11, STATUS='KEEP')
         CLOSE (12, STATUS='KEEP')
         CLOSE (13, STATUS='KEEP')
#endif /* GEN_INIT */
         CLOSE (10, STATUS='KEEP')
      ELSE
         CALL pvmfinitsend(PACK_TYPE, bufid)
         iarray(1)=num
         CALL pvmfpack(INTEGER_4, iarray, 1, 1, info)
         CALL pvmfpack(REAL_8, qx, num, 1, info)
         CALL pvmfpack(REAL_8, qy, num, 1, info)
         CALL pvmfpack(REAL_8, vx, num, 1, info)
         CALL pvmfpack(REAL_8, vy, num, 1, info)
         CALL pvmfpack(INTEGER_4, label, num, 1, info)
         CALL pvmfsend(myparent, MSG_FINAL, info)
         CALL pvmffreebuf(bufid, info)
      ENDIF
#else /* PVM */
      CALL ApplyPerBoundary
#ifdef GEN_INIT
      OPEN (10, FILE='conf.init', STATUS='NEW')
      DO i=1, num
         qx(i)=qx(i)/Len2
         qy(i)=qy(i)/Len2
         WRITE (10, 1005) qx(i), qy(i), vx(i), vy(i), label(i)
      ENDDO
      CLOSE (10, STATUS='KEEP')
#else
      OPEN (10, FILE='conf.final.A', STATUS='NEW')
      OPEN (11, FILE='conf.final.X', STATUS='NEW')
      OPEN (12, FILE='conf.final.Y', STATUS='NEW')
      OPEN (13, FILE='conf.final.Z', STATUS='NEW')

      DO i=1, num
         IF (label(i) .EQ. 0) THEN
            WRITE (10, 1004) qx(i), qy(i), vx(i), vy(i)
         ELSEIF (label(i) .EQ. 1) THEN
            WRITE (11, 1004) qx(i), qy(i), vx(i), vy(i)
         ELSEIF (label(i) .EQ. 2) THEN
            WRITE (12, 1004) qx(i), qy(i), vx(i), vy(i)
         ELSEIF (label(i) .EQ. 3) THEN
            WRITE (13, 1004) qx(i), qy(i), vx(i), vy(i)
         ENDIF
      ENDDO
      CLOSE (10, STATUS='KEEP')
      CLOSE (11, STATUS='KEEP')
      CLOSE (12, STATUS='KEEP')
      CLOSE (13, STATUS='KEEP')
#endif
#endif

 1004 FORMAT (F12.7, F12.7, F12.7, F12.7)
 1005 FORMAT (F12.7, F12.7, F12.7, F12.7, I2)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ReadInit reads the initial positions and velocities from a file.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ReadInit

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i
#ifdef RESERVOIR
      INTEGER lab
      REAL*8  d1, d2, d3, d4

      num=0
      numII=0
#endif /* RESERVOIR */

      OPEN (10, FILE='conf.init', STATUS='OLD')
      DO i=1, nall
#ifdef RESERVOIR
         READ (10, 1002) d1, d2, d3, d4, lab
         IF (lab .EQ. 5) THEN
            numII=numII+1
            qxII(numII)=d1
            qyII(numII)=d2
            vxII(numII)=d3
            vyII(numII)=d4
            label(numII)=lab
         ELSE
            num=num+1
            qxII(num)=d1
            qyII(num)=d2
            vxII(num)=d3
            vyII(num)=d4
            label(num)=lab
         ENDIF
#else /* RESERVOIR */
         READ (10, 1002) qx(i), qy(i), vx(i), vy(i), label(i)
#endif /* RESERVOIR */
      ENDDO
      CLOSE (10, STATUS='KEEP')

 1002 FORMAT (F12.7, F12.7, F12.7, F12.7, I2)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Initialise is starting up the system, i.e. initialising misc.     C
C variables and data structures. The linked-cell structure is setup C
C this routine. The function Icell is a support function.           C
C Files opened:                                                     C
C  15  - time, temperature, pressure                                C
C  16  - time, x(A), x(X), x(Y)                                     C
C  17  - time, Epot, Ekin, Etot (=Epot+Ekin)                        C
C  18  - the structure factor                                       C
C  19  - mean square displacement                                   C
C  20  - the correlation function                                   C
C  21  - time, #Reaction(1), ...                                    C
C  22  - time, Epot(reac 1), Ekin(reac 1), ...                      C
C  23  - time, number of neighbours for A                           C
C  24  - time, number of neighbours for X                           C
C  25  - time, local order                                          C
C  26  - time, subdomain concentrations                             C
C  27  - time, eta                                                  C
C  28  - time, #in reactor, #in reservoir                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE Initialise

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, j, imap, Icell, ix, iy
      REAL*8  sum, sc

C An inline function
      Icell(ix, iy)=1+MOD(ix+Mcell-1, Mcell)
     _              +MOD(iy+Mcell-1, Mcell)*Mcell


#ifdef RESERVOIR
      dof=2.D0*DFLOAT(nall)-4.D0
#else /* !RESERVOIR */
      dof=2.D0*DFLOAT(nall)-2.D0
#endif /* RESERVOIR */

C Length of the simulation box
      Len=DSQRT(DFLOAT(num)/rho)
      Len2=0.5D0*Len
      Len2m=-0.5D0*Len
      iLen=1.D0/Len

C Link-cell structure
      Mcell=INT(Len/bufzone)
      Mcelld=DINT(Len/bufzone)
      Ncell=Mcell**2
      celli=Mcelld/Len
      Mcell2=Mcell**2

      OPEN (15, FILE='mdreac.sample',  STATUS='NEW')
      OPEN (16, FILE='mdreac.conc',    STATUS='NEW')
      OPEN (17, FILE='mdreac.energy',  STATUS='NEW')
      OPEN (18, FILE='mdreac.sf',      STATUS='NEW')
      OPEN (20, FILE='mdreac.cor',     STATUS='NEW')
      OPEN (21, FILE='mdreac.reac',    STATUS='NEW')
      OPEN (22, FILE='mdreac.enereac', STATUS='NEW')
#ifndef PVM_VER
      OPEN (19, FILE='mdreac.diff',    STATUS='NEW')
      OPEN (23, FILE='mdreac.nabor',   STATUS='NEW')
      OPEN (24, FILE='mdreac.naborr',  STATUS='NEW')
      OPEN (25, FILE='mdreac.order',   STATUS='NEW')
      OPEN (26, FILE='mdreac.subconc', STATUS='NEW')
      OPEN (27, FILE='mdreac.eta',     STATUS='NEW')
#ifdef RESERVOIR
      OPEN (28, FILE='mdreac.extend',  STATUS='NEW')
#endif /* RESERVOIR */
#endif /* PVM_VER */
      DO j=1, Mcell
         DO i=1, Mcell
            imap=4*(Icell(i, j)-1)
            map(imap+1)=Icell(i+1, j-1)
            map(imap+2)=Icell(i+1,   j)
            map(imap+3)=Icell(i+1, j+1)
            map(imap+4)=Icell(i,   j+1)
         ENDDO
      ENDDO

C Scale positions
      DO i=1, num
         qx(i)=qx(i)*Len2
         qy(i)=qy(i)*Len2
      ENDDO
#ifdef RESERVOIR
      DO i=1, numII
         qxII(i)=qxII(i)*Len2II
         qyII(i)=qyII(i)*Len2II
      ENDDO
#endif /* RESERVOIR */

      CALL ApplyPerBoundary

C Scale velocities according to temperature
      sum=0.D0
      DO i=1, num
         sum=sum+vx(i)**2+vy(i)**2
      ENDDO
#ifdef RESERVOIR
      DO i=1, numII
         sum=sum+vxII(i)**2+vyII(i)**2
      ENDDO
#endif /* RESERVOIR */
#ifdef TWO_THERMO
      sc=DSQRT(dof*0.5D0*(tempA+tempC)/sum)
#else /* TWO_THERMO */
      sc=DSQRT(dof*temp/sum)
#endif /* TWO_THERMO */
      DO i=1, num
         vx(i)=vx(i)*sc
         vy(i)=vy(i)*sc
      ENDDO
#ifdef RESERVOIR
      DO i=1, numII
         vxII(i)=vxII(i)*sc
         vyII(i)=vyII(i)*sc
#endif /* RESERVOIR */

      sum=0.D0
      DO i=1, num
         sum=sum+vx(i)
      ENDDO
      sum=sum/DFLOAT(num)
      DO i=1, num
         vx(i)=vx(i)-sum
      ENDDO
      sum=0.D0
      DO i=1, num
         sum=sum+vy(i)
      ENDDO
      sum=sum/DFLOAT(num)
      DO i=1, num
         vy(i)=vy(i)-sum
      ENDDO

#ifdef RESERVOIR
      sum=0.D0
      DO i=1, numII
         sum=sum+vxII(i)
      ENDDO
      sum=sum/DFLOAT(numII)
      DO i=1, numII
         vxII(i)=vxII(i)-sum
      ENDDO
      sum=0.D0
      DO i=1, numII
         sum=sum+vyII(i)
      ENDDO
      sum=sum/DFLOAT(numII)
      DO i=1, numII
         vyII(i)=vyII(i)-sum
      ENDDO
#endif /* RESERVOIR */


C Misc.
#ifdef TWO_THERMO
      eta(0)=0.D0
      eta(1)=0.D0
      eta(2)=0.D0
#else /* TWO_THEOMO */
      eta=0.D0
#endif /* TWO_THEOMO */
      CALL MarInit(1, 1)

C No dimers
      DO i=1, ntot
         dimer(i)=0
      ENDDO

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DeInitialise closes various files, prints statistics, and cleans  C
C up the system.                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DeInitialise

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

C Close files
      CLOSE (15, STATUS='KEEP')
      CLOSE (16, STATUS='KEEP')
      CLOSE (17, STATUS='KEEP')
      CLOSE (18, STATUS='KEEP')
      CLOSE (19, STATUS='KEEP')
      CLOSE (20, STATUS='KEEP')
      CLOSE (21, STATUS='KEEP')
#ifndef PVM_VER
      CLOSE (22, STATUS='KEEP')
      CLOSE (23, STATUS='KEEP')
      CLOSE (24, STATUS='KEEP')
      CLOSE (25, STATUS='KEEP')
      CLOSE (26, STATUS='KEEP')
      CLOSE (27, STATUS='KEEP')
#endif

C Print statistics
#ifdef TIMING
      calctime=etime(timearray)-calctime
#ifdef PVM_VER
      PRINT *, 'Total time         :', DFLOAT(calctime)/100.D0/60.D0
      PRINT *, 'Communication time :', DFLOAT(commtime)/100.D0/60.D0
#else
      PRINT *, 'Total time :', DFLOAT(calctime)/100.D0/60.D0
#endif
#endif
      END


#ifdef PVM_VER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C InitPVM is the initialisation routine for the parallel version of C
C program. This routine calls Initialise.                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE InitPVM

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
      INCLUDE 'fpvm3.h'

      INTEGER i, j, ipp, dummy, iy, ip(MaxTask), ii
      INTEGER info, bufid, ilist(ntot+6, MaxTask)
      REAL*8  flist1(ntot, MaxTask), flist2(ntot, MaxTask)
      REAL*8  flist3(ntot, MaxTask), flist4(ntot, MaxTask)
      REAL*8  flist(10), leny, ileny
      INTEGER iarray(ntot), ix, imap, Icell
      REAL*8  farray(ntot+10)
#ifdef TIMING
      REAL*4  itime
#endif
C An inline function
      Icell(ix, iy)=1+MOD(ix+Mcell-1, Mcell)
     _              +MOD(iy+Mcell-1, Mcell)*Mcell


#ifdef DEBUG
      PRINT *, 'BEGIN of InitPVM'
#endif

      CALL pvmfmytid(mytid)
      CALL pvmfparent(myparent)
      IF (myparent .EQ. PvmNoParent) THEN
C-------------------------------
C Parent process is initialised
C-------------------------------
         CALL ReadParameters
#ifdef GEN_INIT
         CALL GenInit
#else
         CALL ReadInit
#endif
         CALL Initialise
C PBC in the y-direction
         DO i=1, num
            qy(i)=qy(i)-Len*DNINT(qy(i)*iLen)
         ENDDO
#ifdef TIMING
         itime=etime(timearray)
#endif

         nall=na+nx+ny
C Start up slaves
         CALL pvmfspawn(EXECNAME,
     _        0, "*", ntask-1, alltids, dummy)
#ifdef DEBUG
         PRINT *, 'Number of task spawned: ', dummy
#endif
         leny=Len/DFLOAT(ntask)
         ileny=1.D0/leny
         DO i=1, ntask
            ip(i)=0
         ENDDO
         DO i=1, nall
            iy=INT((qy(i)+Len2)*ileny)+1
            ipp=ip(iy)
            flist1(ipp+1, iy)=qx(i)
            flist2(ipp+1, iy)=qy(i)
            flist3(ipp+1, iy)=vx(i)
            flist4(ipp+1, iy)=vy(i)
            ilist(ipp+1, iy)=label(i)
            ip(iy)=ip(iy)+1
         ENDDO
         ymin=Len2m
         ymax=Len2m+leny
         Mcell=INT(leny/bufzone)*ntask
         Mcell2=Mcell**2
         CPUid=0
         celli=DFLOAT(Mcell)/Len

         firstcell=1
         lastcell=Mcell2/ntask
         DO i=1, ip(1)
            qx(i)=flist1(i, 1)
            qy(i)=flist2(i, 1)
            vx(i)=flist3(i, 1)
            vy(i)=flist4(i, 1)
            label(i)=ilist(i, 1)
         ENDDO
         num=ip(1)
         DO ii=2, ntask
            i=ii-1
            ipp=ip(ii)
            flist(1)=leny*DFLOAT(i)-Len2
            flist(2)=leny*DFLOAT(i+1)-Len2
            flist(3)=temp
            flist(4)=r2cut1
            flist(5)=Len
            flist(6)=celli
            flist(7)=r2cut2
            flist(8)=prob1
            flist(9)=prob2
            flist(10)=rreac2
            iarray(1)=ip(ii)
            IF (i .EQ. 1) THEN
               iarray(3)=mytid
               CPUup=alltids(1)
            ELSE
               iarray(3)=alltids(i-1)
            ENDIF
            IF (i .EQ. ntask-1) THEN
               iarray(2)=mytid
               CPUdown=alltids(i)
            ELSE
               iarray(2)=alltids(i+1)
            ENDIF
            iarray(4)=step
            iarray(5)=0
            iarray(6)=Mcell
            iarray(7)=ntask
            iarray(8)=i
            Ncell=Mcell**2/ntask

            CALL pvmfinitsend(PACK_TYPE, bufid)
            CALL pvmfpack(INTEGER_4, iarray, 8, 1, info)
            DO j=1, ipp
               iarray(j)=ilist(j, ii)
            ENDDO
            CALL pvmfpack(INTEGER_4, iarray, ipp, 1, info)
            DO j=1, ipp
               farray(j)=flist1(j, ii)
            ENDDO
            CALL pvmfpack(REAL_8, farray, ipp, 1, info)
            DO j=1, ipp
               farray(j)=flist2(j, ii)
            ENDDO
            CALL pvmfpack(REAL_8, farray, ipp, 1, info)
            DO j=1, ipp
               farray(j)=flist3(j, ii)
            ENDDO
            CALL pvmfpack(REAL_8, farray, ipp, 1, info)
            DO j=1, ipp
               farray(j)=flist4(j, ii)
            ENDDO
            CALL pvmfpack(REAL_8, farray, ipp, 1, info)
            CALL pvmfpack(REAL_8, flist, 10, 1, info)
            CALL pvmfsend(alltids(i), MSG_INIT, info)
            CALL pvmffreebuf(bufid, info)
#ifdef TIMING
            itime=etime(timearray)-itime
            commtime=commtime+itime
#endif
#ifdef DEBUG
            PRINT *, 'Sending INIT to child', alltids(i)
#endif
         ENDDO
      ELSE
C---------------------------------------
C child processes are initialised
C---------------------------------------
#ifdef DEBUG
         PRINT *, 'Process', mytid, 'Waiting for INIT from parent'
#endif
         CALL pvmfrecv(myparent, MSG_INIT, bufid)
#ifdef DEBUG
         PRINT *, 'Process ', mytid, 'Received INIT from parent'
#endif
         CALL pvmfunpack(INTEGER_4, iarray, 8, 1, info)
         num=iarray(1)
         CPUup=iarray(2)
         CPUdown=iarray(3)
         step=iarray(4)
C iarray(5) is not used
         Mcell=iarray(6)
         ntask=iarray(7)
         CPUid=iarray(8)
         Ncell=Mcell**2/ntask
         Mcell2=Mcell**2
         IF (num .GT. 0) THEN
            CALL pvmfunpack(INTEGER_4, label, num, 1, info)
            CALL pvmfunpack(REAL_8, qx, num, 1, info)
            CALL pvmfunpack(REAL_8, qy, num, 1, info)
            CALL pvmfunpack(REAL_8, vx, num, 1, info)
            CALL pvmfunpack(REAL_8, vy, num, 1, info)
         ENDIF
         CALL pvmfunpack(REAL_8, flist, 10, 1, info)
         CALL pvmffreebuf(bufid, info)
         ymin=flist(1)
         ymax=flist(2)
         temp=flist(3)
         r2cut1=flist(4)
         Len=flist(5)
         Len2=0.5D0*Len
         Len2m=-0.5D0*Len
         iLen=1.D0/Len
         celli=flist(6)
         firstcell=CPUid*Mcell2/ntask+1
         lastcell=(CPUid+1)*Mcell2/ntask
         r2cut2=flist(7)
         prob1=flist(8)
         prob2=flist(9)
         rreac2=flist(10)
      ENDIF

#ifdef DEBUG
      PRINT *, 'Paramters for Process ', mytid
      PRINT *, 'CPU number        =', CPUid
      PRINT *, 'Ncell, Mcell      =', Ncell, Mcell
      PRINT *, 'step, temp        =', step, temp
      PRINT *, 'ymin, ymax        =', ymin, ymax
      PRINT *, 'num               =', num
      PRINT *, 'Len, celli        =', Len, celli
      PRINT *, 'Len2, Len2m, iLen =', Len2, Len2m, iLen
      PRINT *, 'first cell        =', firstcell
      PRINT *, 'last cell         =', lastcell
#endif

      DO j=1, Mcell
         DO i=1, Mcell
            imap=4*(Icell(i, j)-1)
            map(imap+1)=Icell(i-1, j+1)
            map(imap+2)=Icell(i,   j+1)
            map(imap+3)=Icell(i+1, j+1)
            map(imap+4)=Icell(i+1, j)
         ENDDO
      ENDDO

      CALL MarInit(1, 1)

#ifdef DEBUG
      PRINT *, 'END of InitPVM'
#endif
      END
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ReadParameters opens and reads the parameters given on disk. The  C
C parameters contain all the information about then simulation.     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ReadParameters

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      CHARACTER*15 dummy

      OPEN (10, FILE='mdreac.param', STATUS='OLD')
C Number of particles
      READ (10, 4020) dummy, na
      READ (10, 4020) dummy, nx
      READ (10, 4020) dummy, ny
      READ (10, 4020) dummy, nz
#ifdef FIVE_SPECIES
      READ (10, 4020) dummy, ns
      num=na+nx+ny+nz+ns
#else
      num=na+nx+ny+nz
#endif
      nall=num
C density and cut-offs
      READ (10, 4010) dummy, rho
      READ (10, 4010) dummy, rcut1
      r2cut1=rcut1**2
      READ (10, 4010) dummy, rcut2
      r2cut2=rcut2**2
C temperatures and number of steps
#ifdef TWO_THERMO
      READ (10, 4010) dummy, tempA
      READ (10, 4010) dummy, tempC
#else /* TWO_THERMO */
      READ (10, 4010) dummy, temp
#endif /* TWO_THERMO */
      READ (10, 4020) dummy, step

      dof=2.D0*DFLOAT(nall)-2.D0
#ifdef TWO_THERMO
      invQ=1.D0/(dof*tau*tau*0.5D0*(tempA+tempC))
#else /* TWO_THERMO */
      invQ=1.D0/(dof*tau*tau*temp)
#endif /* TWO_THERMO */

C reaction parameters
      READ (10, 4010) dummy, prob1
      READ (10, 4010) dummy, prob2
      READ (10, 4010) dummy, prob3
      Ea1=prob1
      Ea2=prob2
      Ea3=prob3
#ifdef MANY_REAC_PARAMS
      READ (10, 4010) dummy, prob4
      READ (10, 4010) dummy, prob5
      READ (10, 4010) dummy, prob6
#endif /* MANY_REAC_PARAMS */
      READ (10, 4010) dummy, rreac
      rreac2=rreac**2

C misc. parameters
c      READ (10, 4010) dummy, Rnabor
c      R2nabor=Rnabor**2
#ifdef PVM_VER
C number of cpus
      READ (10, 4020) dummy, ntask
#endif

#ifdef DEBUG
      PRINT *, 'ReadParameters:'
      PRINT *, 'na, nx, ny      = ', na, nx, ny
      PRINT *, 'rho             = ', rho
      PRINT *, 'step, temp      = ', step, temp
      PRINT *, 'r_c LJ, r_c WCA = ', rcut1, rcut2
      PRINT *, 'P1, P2, Rreac   = ', prob1, prob2, rreac
      PRINT *, 'P3, P4, P5, P6  = ', prob3, prob4, prob5, prob6
      PRINT *, 'dof             = ', dof
#ifdef PVM_VER
      PRINT *, 'ntask           = ', ntask
#endif
#endif

 4010 FORMAT (A15, F12.8)
 4020 FORMAT (A15, I10)

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NumToStr converts an integer to a string, i.e. an array of charac-C
C ters. The input is the parameter i and s is the string while n is C
C the length of string.                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NumToStr(i, s, n)

      IMPLICIT NONE
      INTEGER      i, j, k, n
      CHARACTER*25 s
      CHARACTER*1  swap

      j=i
      n=1
 4010 k=MOD(j, 10)
      s(n:n)=CHAR(k+ICHAR('0'))
      j=j/10
      n=n+1
      IF (j .GE. 10) GOTO 4010
      s(n:n)=CHAR(j+ICHAR('0'))

C The number is backwards
      DO j=1, n/2
         k=n-j+1
         swap(1:1)=s(j:j)
         s(j:j)=s(k:k)
         s(k:k)=swap(1:1)
      ENDDO

C End-of-string for the rest
      DO j=n+1, 25
         s(j:j)=CHAR(0)
      ENDDO
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NumToStr2 converts an integer to a string.                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NumToStr2(hejmeddig, s, col)

      IMPLICIT NONE
      INTEGER  hejmeddig, col
      CHARACTER*(*) s

      WRITE (s(col:), 1000) hejmeddig

 1000 FORMAT (I10)
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DumpPos is writing the positions of the particles to a file. Each C
C species have its own file.                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DumpPos(stepno)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, stepno, n
      CHARACTER*25  filename, str
      CHARACTER*100 prnstr

#ifdef DEBUG
      PRINT *, 'BEGIN of DumpPos'
#endif

      filename(1:7)='mdreac.'
C      CALL NumToStr(stepno, str, n)
C      filename(8:7+n)=str
      CALL NumToStr2(stepno, filename, 8)

      OPEN (10, FILE=filename, STATUS='NEW')

      DO i=1, nall
         WRITE (10, 5010) qx(i), qy(i), label(i)
      ENDDO

      CLOSE (10, STATUS='KEEP')

#ifdef DEBUG
      PRINT *, 'END of DumpPos'
#endif

 5010 FORMAT (F18.8, F18.8, I4)
 5020 FORMAT (A100)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The following subroutine can calculate the structure factor. The C
C routine is written by Mohamed Laradji (mohamed@fki.dtu.dk).      C
C The routine uses a fast fourier transform, which can be found in C
C [5]                                                              C
C Both routines have been modified a bit by Kenneth Geisshirt      C
C (kneth@fatou.ruc.dk).                                            C
C                                                                  C
C In your code, you should incorporate this every som many steps   C
C                                                                  C
C     call stfac(lmn, xlx, yly, sfunc, corel)                      C
C     DO 10 i=0, 31                                                C
C       qq=DFLOAT(i)*pi/xlx                                        C
C       rr=DFLOAT(i)*xlx/64.d0                                     C
C       WRITE (11,9999) qq, sfunc(i)                               C
C       WRITE (12,9999) rr, corel(i)                               C
C 10  CONTINUE                                                     C
C                                                                  C
C In your case, you must have                                      C
C    species=1   for A-component for example                       C
C    species=0   otherwise                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c=============================================================
c Structure factor and correlation function subroutine
c The structure factor is only computed for one of the
c species, namely A.
c=============================================================

      SUBROUTINE stfac(sfunc, corel)

      IMPLICIT DOUBLE PRECISION(a, c-h, o-z)
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      PARAMETER (nn=2, l=sfsites, nsite=l**2)
      REAL*4    data(2*nsite), sf(nsite), cr(nsite),param
      DIMENSION sfunc(0:256), corel(0:256), number(0:256)
      INTEGER   ndim(2), species(ntot)
      REAL*8    factor

 2    xsize=Len/DFLOAT(l)
      param=xsize**2/FLOAT(nsite)
      DO 30 i=1, num
         IF (label(i) .EQ. 0) THEN
            species(i)=1
         ELSE
            species(i)=0
         ENDIF
 30   CONTINUE

      factor=DFLOAT(nsite)/DFLOAT(num)
      DO 10 i=0, l
         sfunc(i)=0.D0
         corel(i)=0.D0
         number(i)=0.D0
 10   CONTINUE

      DO 20 i=1, 2*nsite
         data(i)=0.0
 20   CONTINUE

      ndim(1)=l
      ndim(2)=l
      vol=Len*Len/DFLOAT(nsite)
      DO 90 i=1, num
         ix=INT(((qx(i)+Len2)*iLen)*DFLOAT(l))+1
         iy=INT(((qy(i)+Len2)*iLen)*DFLOAT(l))+1
         j=ix+(iy-1)*l
         ind1=2*j
         ind2=ind1-1
         data(ind1)=0.0
         data(ind2)=data(ind2)+FLOAT(species(i))
 90   CONTINUE
c---------------------------------------------
c First, calculate the structure factor
c---------------------------------------------
      CALL fourn(data, ndim, nn, 1)
      DO 40 i=1, nsite
         ind1=2*i
         ind2=ind1-1
         sf(i)=(data(ind1)**2+data(ind2)**2)*param
 40   CONTINUE
      DO 50 i=1, nsite
         ind1=2*i
         ind2=ind1-1
         data(ind2)=sf(i)
         data(ind1)=0.0
 50   CONTINUE
c---------------------------------------------
c Second, calculate the correlation function
c---------------------------------------------
      CALL fourn(data, ndim, nn, -1)
      DO 60 i=1, nsite
         cr(i)=data(2*i-1)/float(nsite)
 60   CONTINUE

      DO 70 i=1, nsite
         iy=MOD(i-1, l)
         ix=INT((i-1)/l)
         IF (ix .GT. l/2) ix=l-ix
         IF (iy .GT. l/2) iy=l-iy
         kval=INT(DSQRT(DFLOAT(ix**2+iy**2)))
         sfunc(kval)=sfunc(kval)+sf(i)
         corel(kval)=corel(kval)+cr(i)
         number(kval)=number(kval)+1
 70   CONTINUE
      DO 80 i=0, l/2
         IF (number(i) .GT. 0) THEN
            sfunc(i)=sfunc(i)/DFLOAT(number(i))
            corel(i)=corel(i)/(DFLOAT(number(i))*xsize**2)
         ENDIF
 80   CONTINUE
      RETURN
      END


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Fast Fourier transform - from [5].
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE fourn(data, nn, ndim, isign)

      REAL*8 wr, wi, wpr, wpi, wtemp, theta
      DIMENSION nn(ndim), data(*)

      ntot=1
      DO 11 idim=1, ndim
         ntot=ntot*nn(idim)
 11   CONTINUE
      nprev=1
      DO 18 idim=1, ndim
         n=nn(idim)
         nrem=ntot/(n*nprev)
         ip1=2*nprev
         ip2=ip1*n
         ip3=ip2*nrem
         i2rev=1
         DO 14 i2=1, ip2, ip1
            IF (i2 .LT. i2rev) THEN
               DO 13 i1=i2, i2+ip1-2, 2
                  DO 12 i3=i1, ip3, ip2
                     i3rev=i2rev+i3-i2
                     tempr=data(i3)
                     tempi=data(i3+1)
                     data(i3)=data(i3rev)
                     data(i3+1)=data(i3rev+1)
                     data(i3rev)=tempr
                     data(i3rev+1)=tempi
 12               CONTINUE
 13            CONTINUE
            ENDIF
            ibit=ip2/2
 1          IF ((ibit .GE. ip1) .AND. (i2rev .GT. ibit)) THEN
               i2rev=i2rev-ibit
               ibit=ibit/2
               GOTO 1
            ENDIF
            i2rev=i2rev+ibit
 14      CONTINUE
         ifp1=ip1
 2       IF (ifp1 .LT. ip2) THEN
            ifp2=2*ifp1
            theta=isign*6.28318530717959D0/(ifp2/ip1)
            wpr=-2.D0*DSIN(0.5D0*theta)**2
            wpi=DSIN(theta)
            wr=1.D0
            wi=0.D0
            DO 17 i3=1, ifp1, ip1
               DO 16 i1=i3, i3+ip1-2, 2
                  DO 15 i2=i1, ip3, ifp2
                     k1=i2
                     k2=k1+ifp1
                     tempr=SNGL(wr)*data(k2)-SNGL(wi)*data(k2+1)
                     tempi=SNGL(wr)*data(k2+1)+SNGL(wi)*data(k2)
                     data(k2)=data(k1)-tempr
                     data(k2+1)=data(k1+1)-tempi
                     data(k1)=data(k1)+tempr
                     data(k1+1)=data(k1+1)+tempi
 15               CONTINUE
 16            CONTINUE
               wtemp=wr
               wr=wr*wpr-wi*wpi+wr
               wi=wi*wpr+wtemp*wpi+wi
 17         CONTINUE
            ifp1=ifp2
            GOTO 2
         ENDIF
         nprev=n*nprev
 18   CONTINUE
      RETURN
      END

#ifdef PVM_VER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RecvPos receives positions from all the subprocesses, i.e. this   C
C subroutine should only be called by the parent process.           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RecvPos

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
      INCLUDE 'fpvm3.h'

      INTEGER i, info, bufid, iarray(1), j
#ifdef TIMING
      REAL*4  itime
#endif
#ifdef DEBUG
      PRINT *, 'BEGIN RecvPos'
#endif

#ifdef TIMING
      itime=etime(timearray)
#endif
      j=num
      DO i=1, ntask-1
         CALL pvmfrecv(alltids(i), MSG_POSS, bufid)
         CALL pvmfunpack(INTEGER_4, iarray, 1, 1, info)
         IF (iarray(1) .NE. 0) THEN
            CALL pvmfunpack(REAL_8, qx(j+1), iarray(1), 1, info)
            CALL pvmfunpack(REAL_8, qy(j+1), iarray(1), 1, info)
            CALL pvmfunpack(INTEGER_4, label(j+1), iarray(1), 1, info)
            j=j+iarray(1)
         ENDIF
         CALL pvmffreebuf(bufid, info)
      ENDDO

#ifdef TIMING
      itime=etime(timearray)-itime
      commtime=commtime+itime
#endif

#ifdef DEBUG
      PRINT *, 'END RecvPos'
#endif
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SendPos sends the positions to the parent process, i.e. this sub- C
C routine should only be called by a child process.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SendPos

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"
      INCLUDE 'fpvm3.h'

      INTEGER info, bufid, iarray(1)

#ifdef DEBUG
      PRINT *, 'BEGIN SendPos'
#endif

      CALL pvmfinitsend(PACK_TYPE, bufid)
      iarray(1)=num
      CALL pvmfpack(INTEGER_4, iarray, 1, 1, info)
      CALL pvmfpack(REAL_8, qx, num, 1, info)
      CALL pvmfpack(REAL_8, qy, num, 1, info)
      CALL pvmfpack(INTEGER_4, label, num, 1, info)
      CALL pvmfsend(myparent, MSG_POSS, info)
      CALL pvmffreebuf(bufid, info)

#ifdef DEBUG
      PRINT *, 'END SendPos'
#endif
      END

#endif

#ifdef GEN_INIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C GenInit generates initial configurations, i.e. positions, velo-   C
C cities, and labels.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GenInit

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, icell, jcell, nsqr, index
      REAL*8  clen

C---------------------------------------
C Initialisation
C---------------------------------------
      CALL MarInit(1, 1)

C---------------------------------------
C The velocities - gaussian distributed
C---------------------------------------
      DO i=1, num
         CALL RandGauss(vx(i), vy(i))
      ENDDO

C---------------------------------------
C Place the particles lattice-like
C---------------------------------------
      nsqr=INT(DSQRT(DFLOAT(num)))
      clen=2.D0/DFLOAT(nsqr)
      DO i=1, num
         icell=MOD(i-1, nsqr)
         jcell=(i-1)/nsqr
         qx(i)=(0.5D0+DFLOAT(icell))*clen-1.D0
         qy(i)=(0.5D0+DFLOAT(jcell))*clen-1.D0
      ENDDO

C---------------------------------------
C Assign labels
C---------------------------------------
      DO i=1, num
         label(i)=-1
      ENDDO
      DO i=1, na
 5010    CALL RandInt(num, index)
         IF (label(index) .NE. -1) GOTO 5010
         label(index)=0
      ENDDO
      DO i=1, nx
 5020    CALL RandInt(num, index)
         IF (label(index) .NE. -1) GOTO 5020
         label(index)=1
      ENDDO
      DO i=1, ny
 5030    CALL RandInt(num, index)
         IF (label(index) .NE. -1) GOTO 5030
         label(index)=2
      ENDDO
      DO i=1, nz
 5040    CALL RandInt(num, index)
         IF (label(index) .NE. -1) GOTO 5040
         label(index)=3
      ENDDO
#ifdef FIVE_SPECIES
      DO i=1, ns
 5050    CALL RandInt(num, index)
         IF (label(index) .NE. -1) GOTO 5050
         label(index)=5
      ENDDO
#endif
      END

#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C InitDiffusion is initialising the variables which are used in the C
C estimation of the diffusion coefficient. The diffusion coefficientC
C is estimated by mean square displacement [7, pp. 199-206].        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE InitDiffusion

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i

      DO i=1, num
         nboxx(i)=0
         nboxy(i)=0
         R0x(i)=qx(i)
         R0y(i)=qy(i)
      ENDDO
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MeanSqrDispl calculates the mean square displacement.             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE MeanSqrDispl

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i

      meansqr=0.D0
      DO i=1, num
         meansqr=meansqr+(R0x(i)-qx(i)+DFLOAT(nboxx(i))*Len)**2
     _        +(R0y(i)-qy(i)+DFLOAT(nboxy(i))*Len)**2
      ENDDO
      meansqr=meansqr/DFLOAT(num)
      END


#ifndef PVM_VER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SumNabor adds the number of neighbours of a given particle and    C
C calculates the preferred direction of A-A pairs.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SumNabor(k)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER k

      INTEGER cll, i, klabel, nabor1, nabor2, ilabel
      INTEGER cll2, cell0
      INTEGER cllx, clly, cllx2, clly2
      REAL*8  rxk, ryk, r2, rx, ry, rxi, ryi
      REAL*8  ux0, uy0, ux, uy, dot, mean, x0, y0
      INTEGER state, numb
      REAL*8  normx, normy, phi

      rxk=qx(k)
      ryk=qy(k)
      klabel=label(k)+1
      cllx=INT((qx(k)+Len2)*celli)+1
      clly=INT((qy(k)+Len2)*celli)+1

      mean=0.D0
      numb=0
      normx=0.D0
      normy=0.D0
      phi=0.D0

      DO nabor1=-1, 1
         DO nabor2=-1, 1
            cllx2=cllx+nabor1
            clly2=clly+nabor2
            cll=1+MOD(cllx2+Mcell-1, Mcell)
     _           +MOD(clly2+Mcell-1, Mcell)*Mcell
            i=head(cll)
 100        IF (i .NE. 0) THEN
               IF (i .NE. k) THEN
                  rxi=qx(i)
                  ryi=qy(i)
                  ilabel=label(i)+1
                  IF (ilabel .EQ. klabel) THEN
                     rx=rxk-rxi
                     ry=ryk-ryi
                     rx=rx-Len*DNINT(rx*iLen)
                     ry=ry-Len*DNINT(ry*iLen)
                     r2=rx**2+ry**2
                     IF (r2 .LE. r2cut1) THEN
                        phi=phi+DACOS(rx/DSQRT(r2))
                        numb=numb+1
                     ENDIF
                  ENDIF
               ENDIF
               i=list(i)
               GOTO 100
            ENDIF
         ENDDO
      ENDDO

      IF (numb .NE. 0) THEN
         phi=phi/DFLOAT(numb)
      ENDIF


      DO nabor1=-1, 1
         DO nabor2=-1, 1
            cllx2=cllx+nabor1
            clly2=clly+nabor2
            cll=1+MOD(cllx2+Mcell-1, Mcell)
     _           +MOD(clly2+Mcell-1, Mcell)*Mcell
            i=head(cll)
 110        IF (i .NE. 0) THEN
               IF (i .NE. k) THEN
                  rxi=qx(i)
                  ryi=qy(i)
                  ilabel=label(i)+1
                  IF (klabel .EQ. ilabel) THEN
                     rx=rxk-rxi
                     ry=ryk-ryi
                     rx=rx-Len*DNINT(rx*iLen)
                     ry=ry-Len*DNINT(ry*iLen)
                     r2=rx**2+ry**2
                     IF (r2 .LE. r2cut1) THEN
                        mean=mean+(rx*DCOS(phi)+ry*DSIN(phi))/DSQRT(r2)
                     ENDIF
                  ENDIF
               ENDIF
               i=list(i)
               GOTO 110
            ENDIF
         ENDDO
      ENDDO

      IF (numb .NE. 0) THEN
         meanAA=mean/DFLOAT(numb)
         numnabor=numb
      ENDIF
      END
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DomainConc counts the number of particles (A only) is the four    C
C quadrants.                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DomainConc(time)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER time, i, number1(4), number2(4)
      REAL*8  conc(4)

      DO i=1, 4
         number1(i)=0
         number2(i)=0
      ENDDO

      DO i=1, num
         IF (label(i) .EQ. 0) THEN
            IF ((qx(i) .GE. 0.D0) .AND. (qy(i) .GE. 0.D0)) THEN
               number1(1)=number1(1)+1
            ELSEIF ((qx(i) .LT. 0.D0) .AND. (qy(i) .GE. 0.D0)) THEN
               number1(2)=number1(2)+1
            ELSEIF ((qx(i) .LT. 0.D0) .AND. (qy(i) .LT. 0.D0)) THEN
               number1(3)=number1(3)+1
            ELSE
               number1(4)=number1(4)+1
            ENDIF
         ENDIF
      ENDDO

      DO i=1, num
         IF ((qx(i) .GE. 0.D0) .AND. (qy(i) .GE. 0.D0)) THEN
            number2(1)=number2(1)+1
         ELSEIF ((qx(i) .LT. 0.D0) .AND. (qy(i) .GE. 0.D0)) THEN
            number2(2)=number2(2)+1
         ELSEIF ((qx(i) .LT. 0.D0) .AND. (qy(i) .LT. 0.D0)) THEN
            number2(3)=number2(3)+1
         ELSE
            number2(4)=number2(4)+1
         ENDIF
      ENDDO

      DO i=1, 4
         IF (number2(i) .EQ. 0) THEN
            conc(i)=0.D0
         ELSE
            conc(i)=DFLOAT(number1(i))/DFLOAT(number2(i))
         ENDIF
      ENDDO

      WRITE (26, 1000) time, conc(1), conc(2), conc(3), conc(4)

 1000 FORMAT (I15, F18.8, F18.8, F18.8, F18.8)
      END



#ifdef RESERVOIR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The subroutine Exchange makes a Monte Carlo attempt to move par-  C
C ticles from the reservoir to the reactor or vice versa.           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE Exchange

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      REAL*8  drandy, trailX, trailY
      REAL*8  dEI, dEII
      INTEGER indexI, indexII

C Pick particle from reservoir or reactor?
C Less than 0.5 => pick from reservoir
C Greather than 0.5 => pick for reactor
      CALL RandUnit(drandy)
      IF (drandy .GT. 0.5D0) THEN
 1001    CALL RandInt(num, indexI)
         IF (label(indexI) .NE. 5) GOTO 1001
         CALL RandUnit(trailX)
         CALL RandUnit(trailY)
         trailX=Len2*(trailX-0.5D0)
         trailY=Len2*(trailY-0.5D0)
         CALL EnergyItoII(indexI, trailX, trailY, dEI, dEII)
         dE=dEI+dEII
     _        +DFLOAT(num)*temp*DLOG(DFLOAT(num-1)/DFLOAT(num))
     _        +DFLOAT(numII)*temp*DLOG(DFLOAT(numII+1)/DFLOAT(numII))
     _        -2.D0*temp*DLOG(LenII/DFLOAT(numII-1))
     _        +2.D0*temp*DLOG(Len/DFLOAT(num-1))

         CALL RandUnit(drandy)
         IF (dE .LE. 0.D0) THEN
            numII=numII+1
            qxII(numII)=qx(indexI)
            qyII(numII)=qy(indexI)
            vxII(numII)=vx(indexI)
            vyII(numII)=vy(indexI)
            label(num)=0
            DO i=indexI+1, num
               qx(i-1)=qx(i)
               qy(i-1)=qy(i)
               vx(i-1)=vx(i)
               vy(i-1)=vx(i)
               label(i-1)=label(i)
            ENDDO
         ELSEIF (DEXP(-dE/temp) .GT. drandy) THEN
            numII=numII+1
            qxII(numII)=qx(indexI)
            qyII(numII)=qy(indexI)
            vxII(numII)=vx(indexI)
            vyII(numII)=vy(indexI)
            label(num)=0
            DO i=indexI+1, num
               qx(i-1)=qx(i)
               qy(i-1)=qy(i)
               vx(i-1)=vx(i)
               vy(i-1)=vx(i)
               label(i-1)=label(i)
            ENDDO
         ENDIF
      ELSE
         CALL RandInt(numII, indexI)
         CALL RandUnit(trailX)
         CALL RandUnit(trailY)
         trailX=Len2II*(trailX-0.5D0)
         trailY=Len2II*(trailY-0.5D0)
         CALL EnergyIItoI(indexII, trailX, trailY, dEI, dEII)
         dE=dEI+dEII
     _        +DFLOAT(num)*temp*DLOG(DFLOAT(num+1)/DFLOAT(num))
     _        +DFLOAT(numII)*temp*DLOG(DFLOAT(numII-1)/DFLOAT(numII))
     _        +2.D0*temp*DLOG(LenII/DFLOAT(numII-1))
     _        -2.D0*temp*DLOG(Len/DFLOAT(num-1))
         CALL RandUnit(drandy)
         IF (dE .LE. 0.D0) THEN
            num=num+1
            qx(num)=qxII(indexII)
            qy(num)=qyII(indexII)
            vx(num)=vxII(indexII)
            vy(num)=vyII(indexII)
            label(num)=0
            DO i=indexII+1, numII
               qxII(i-1)=qxII(i)
               qyII(i-1)=qyII(i)
               vxII(i-1)=vxII(i)
               vyII(i-1)=vxII(i)
               labelII(i-1)=labelII(i)
            ENDDO
            numII=numII-1
         ELSEIF (DEXP(-dE/temp) .GT. drandy) THEN
            num=num+1
            qx(num)=qxII(indexII)
            qy(num)=qyII(indexII)
            vx(num)=vxII(indexII)
            vy(num)=vyII(indexII)
            label(num)=0
            DO i=indexII+1, numII
               qxII(i-1)=qxII(i)
               qyII(i-1)=qyII(i)
               vxII(i-1)=vxII(i)
               vyII(i-1)=vxII(i)
               labelII(i-1)=labelII(i)
            ENDDO
            numII=numII-1
         ENDIF
      ENDIF
      END


      SUBROUTINE EnergyItoII(index, trailX, trailY, dEI, dEII)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER index
      REAL*8  trailX, trailY, dEI, dEII

      INTEGER i, ilab, c1, c2, cell
      REAL*8  rxi, ryi, rxij, ryij, r2, r6


C Inline function
      Icell(ix, iy)=1+MOD(ix+Mcell-1, Mcell)
     _              +MOD(iy+Mcell-1, Mcell)*Mcell


      dEI=0.D0
      dEII=0.D0

      icellI=Icell(INT((qx(index)+len2)/bufzone),
     _     INT((qy(index)+len2)/bufzone))
      icellII=Icell(INT((trailX+len2II)/bufzone),
     _     INT((trailY+len2II)/bufzone))

C Reactor
      rxi=qx(index)
      ryi=qy(index)
      ilab=label(index)
      DO c1=-1, 1
         DO c2=-1, 1
            cell=Icell(c1+INT((qx(index)+len2)/bufzone),
     _           c2+INT((qy(index)+len2)/bufzone))
            i=head(icell)
 1010       IF (i .EQ. 0) THEN
               IF (i .NE. index) THEN
                  rxij=rxi-qx(i)
                  IF (rxij .GT. len2) THEN
                     rxij=rxij-len2
                  ELSEIF (rxij .LT. len2m) THEN
                     rxij=rxij+len2
                  ENDIF
                  ryij=ryi-qy(i)
                  IF (ryij .GT. len2) THEN
                     ryij=ryij-len2
                  ELSEIF (ryij .LT. len2m) THEN
                     ryij=ryij+len2
                  ENDIF
                  r2=rxij**2+ryij**2
                  IF (ilab .NE. label(i)) THEN
                     IF (r2 .LT. r2cut1) THEN
                        r2=1.D0/r2
                        r6=r2*r2*r2
                        dEI=dEI+r6*(r6-1.D0)
                     ENDIF
                  ELSE
                     IF (r2 .LT. r2cut2) THEN
                        r2=1.D0/r2
                        r6=r2*r2*r2
                        dEI=dEI+r6*(r6-1.D0)
                     ENDIF
                  ENDIF
               ENDIF
               i=list(i)
               GOTO 1010
            ENDIF
         ENDDO
      ENDDO
      dEI=-4.D0*dEI

C Reservoir
      rxi=trailX
      ryi=trailY
      DO c1=-1, 1
         DO c2=-1, 1
            cell=Icell(c1+INT((trailX+len2II)/bufzone),
     _           c2+INT((trailY+len2II)/bufzone))
            i=headII(icell)
 1020       IF (i .EQ. 0) THEN
               IF (i .NE. index) THEN
                  rxij=rxi-qxII(i)
                  IF (rxij .GT. len2II) THEN
                     rxij=rxij-len2II
                  ELSEIF (rxij .LT. len2mII) THEN
                     rxij=rxij+len2II
                  ENDIF
                  ryij=ryi-qyII(i)
                  IF (ryij .GT. len2II) THEN
                     ryij=ryij-len2II
                  ELSEIF (ryij .LT. len2mII) THEN
                     ryij=ryij+len2II
                  ENDIF
                  r2=rxij**2+ryij**2
                  IF (r2 .LT. r2cut2) THEN
                     r2=1.D0/r2
                     r6=r2*r2*r2
                     dEII=dEII+r6*(r6-1.D0)
                  ENDIF
               ENDIF
               i=listII(i)
               GOTO 1020
            ENDIF
         ENDDO
      ENDDO
      dEII=4.D0*dEII
      END


      SUBROUTINE EnergyIItoI(index, trailX, trailY, dEI, dEII)

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER index
      REAL*8  trailX, trailY, dEI, dEII

      INTEGER i, ilab, c1, c2, cell
      REAL*8  rxi, ryi, rxij, ryij, r2, r6


C Inline function
      Icell(ix, iy)=1+MOD(ix+Mcell-1, Mcell)
     _              +MOD(iy+Mcell-1, Mcell)*Mcell


      dEI=0.D0
      dEII=0.D0

      icellI=Icell(INT((qx(index)+len2)/bufzone),
     _     INT((qy(index)+len2)/bufzone))
      icellII=Icell(INT((trailX+len2II)/bufzone),
     _     INT((trailY+len2II)/bufzone))

C Reactor
      rxi=trailX
      ryi=trailY
      ilab=0
      DO c1=-1, 1
         DO c2=-1, 1
            cell=Icell(c1+INT((trailX+len2)/bufzone),
     _           c2+INT((trailY+len2)/bufzone))
            i=head(icell)
 1010       IF (i .EQ. 0) THEN
               IF (i .NE. index) THEN
                  rxij=rxi-qx(i)
                  IF (rxij .GT. len2) THEN
                     rxij=rxij-len2
                  ELSEIF (rxij .LT. len2m) THEN
                     rxij=rxij+len2
                  ENDIF
                  ryij=ryi-qy(i)
                  IF (ryij .GT. len2) THEN
                     ryij=ryij-len2
                  ELSEIF (ryij .LT. len2m) THEN
                     ryij=ryij+len2
                  ENDIF
                  r2=rxij**2+ryij**2
                  IF (ilab .NE. label(i)) THEN
                     IF (r2 .LT. r2cut1) THEN
                        r2=1.D0/r2
                        r6=r2*r2*r2
                        dEI=dEI+r6*(r6-1.D0)
                     ENDIF
                  ELSE
                     IF (r2 .LT. r2cut2) THEN
                        r2=1.D0/r2
                        r6=r2*r2*r2
                        dEI=dEI+r6*(r6-1.D0)
                     ENDIF
                  ENDIF
               ENDIF
               i=list(i)
               GOTO 1010
            ENDIF
         ENDDO
      ENDDO
      dEI=4.D0*dEI

C Reservoir
      rxi=qxII(index)
      ryi=qyII(index)
      DO c1=-1, 1
         DO c2=-1, 1
            cell=Icell(c1+INT((qxII(index)+len2II)/bufzone),
     _           c2+INT((qxII(index)+len2II)/bufzone))
            i=headII(icell)
 1020       IF (i .EQ. 0) THEN
               IF (i .NE. index) THEN
                  rxij=rxi-qxII(i)
                  IF (rxij .GT. len2II) THEN
                     rxij=rxij-len2II
                  ELSEIF (rxij .LT. len2mII) THEN
                     rxij=rxij+len2II
                  ENDIF
                  ryij=ryi-qyII(i)
                  IF (ryij .GT. len2II) THEN
                     ryij=ryij-len2II
                  ELSEIF (ryij .LT. len2mII) THEN
                     ryij=ryij+len2II
                  ENDIF
                  r2=rxij**2+ryij**2
                  IF (r2 .LT. r2cut2) THEN
                     r2=1.D0/r2
                     r6=r2*r2*r2
                     dEII=dEII+r6*(r6-1.D0)
                  ENDIF
               ENDIF
               i=listII(i)
               GOTO 1020
            ENDIF
         ENDDO
      ENDDO
      dEII=-4.D0*dEII
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Same as UpdateLinkCell but for reservoir                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE UpdateLinkCellII

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, j, icell, cell0, jcell
      INTEGER ik1, nabor
      REAL*8  r2, rx, ry, rxi, ryi
      INTEGER k

      ik1=0
      DO icell=1, Mcell2
         i=headII(icell)
 2010    IF (i .NE. 0) THEN
            rxi=qxII(i)
            ryi=qyII(i)
            j=listII(i)
 2020       IF (j .NE. 0) THEN
               rx=rxi-qxII(j)
               ry=ryi-qyII(j)
               IF (rx .GE. Len2) THEN
                  xdiff=-Len
               ELSEIF (rx .LT. Len2m) THEN
                  xdiff=Len
               ELSE
                  xdiff=0.D0
               ENDIF
               rx=rx+xdiff
               IF (ry .GE. Len2) THEN
                  ydiff=-Len
               ELSEIF (ry .LT. Len2m) THEN
                  ydiff=Len
               ELSE
                  ydiff=0.D0
               ENDIF
               ry=ry+ydiff
               r2=rx**2+ry**2
               IF (r2 .LT. bufzone2) THEN
                  ik1=ik1+1
                  ipairII(ik1)=i
                  jpairII(ik1)=j
                  xperII(ik1)=xdiff
                  yperII(ik1)=ydiff
               ENDIF
               j=listII(j)
               GOTO 2020
            ENDIF

C Neighbouring cells
            cell0=4*(icell-1)
            DO nabor=1, 4
               jcell=map(cell0+nabor)
               j=headII(jcell)
 2030          IF (j .NE. 0) THEN
                  rx=rxi-qxII(j)
                  ry=ryi-qyII(j)
                  IF (rx .GE. Len2) THEN
                     xdiff=-Len
                  ELSEIF (rx .LT. Len2m) THEN
                     xdiff=Len
                  ELSE
                     xdiff=0.D0
                  ENDIF
                  rx=rx+xdiff
                  IF (ry .GE. Len2) THEN
                     ydiff=-Len
                  ELSEIF (ry .LT. Len2m) THEN
                     ydiff=Len
                  ELSE
                     ydiff=0.D0
                  ENDIF
                  ry=ry+ydiff
                  r2=rx**2+ry**2
                  IF (r2 .LT. bufzone2) THEN
                     ik1=ik1+1
                     ipairII(ik1)=i
                     jpairII(ik1)=j
                     xperII(ik1)=xdiff
                     yperII(ik1)=ydiff
                  ENDIF
                  j=listII(j)
                  GOTO 2030
               ENDIF
            ENDDO
            i=listII(i)
            GOTO 2010
         ENDIF
      ENDDO
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C LinksII is the same as Links for the reservoir               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE LinksII

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, cll

      DO i=1, Mcell2
         headII(i)=0
      ENDDO

      DO i=1, num
         cll=1+INT((qx(i)+Len2)*celli)
     _        +INT((qy(i)+Len2)*celli)*Mcell
         listII(i)=headII(cll)
         headII(cll)=i
      ENDDO
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ComputeForcesII  - see ComputeForces                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ComputeForcesII

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, j, k
      REAL*8  rxij, ryij, r2, pot, r6

C Init. the forces
      DO i=1, ntot
         fxII(i)=0.D0
         fyII(i)=0.D0
      ENDDO

C Loop over pairs
      DO k=1, numbpairII
         i=ipairII(k)
         j=jpairII(k)
         rxij=qxII(i)-qxII(j)+xperII(k)
         ryij=qyII(i)-qyII(j)+yperII(k)
         r2=rxij**2+ryij**2

C Which potential to use
         IF (r2 .LT. r2cut2) THEN
            r2=1.D0/r2
            r6=r2*r2*r2
            pot=r2*r6*(r6-0.5D0)
            fxII(i)=fxII(i)+rxij*pot
            fyII(i)=fyII(i)+ryij*pot
            fxII(j)=fxII(j)-rxij*pot
            fyII(j)=fyII(j)-ryij*pot
         ENDIF
      ENDDO
      END
#endif /* RESERVOIR */


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MakeVerlet constructs a Verlet list.                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


#ifdef VERLET
      SUBROUTINE MakeVerlet

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i, j
      REAL*8  rxi, ryi, rxij, ryij, r2

      nlist=0
      DO i=1, num-1
         point(i)=nlist+1
         rxi=qx(i)
         ryi=qy(i)
         DO j=i+1, num
            rxij=rxi-qx(j)
            ryij=ryi-qy(j)
            rxij=rxij-Len*DNINT(rxij*iLen)
            ryij=ryij-Len*DNINT(ryij*iLen)
            r2=rxij**2+ryij**2
            IF (r2 .LT. bufzone2) THEN
               nlist=nlist+1
               list(nlist)=j
            ENDIF
         ENDDO
      ENDDO
      point(num)=nlist+1
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Presort determines which particles to be sent to the neighbouring C
C CPU in the ComputeForces routine.                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef PVM_VER
      SUBROUTINE Presort

      IMPLICIT NONE
      INCLUDE 'mdreac.inc'
#include "mdreac.common"

      INTEGER i

      nforce=0
      DO i=1, num
         IF ((qy(i) .GE. ymin) .AND. (qy(i) .GE. ymin+bufzone)) THEN
            nforce=nforce+1
            forcelist(nforce)=i
         ENDIF
      ENDDO
      END

#endif PVM_VER
#endif /* VERLET */


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A number of external routines is included.                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INCLUDE 'randgen.f'
