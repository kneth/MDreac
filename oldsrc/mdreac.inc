CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C (C) Copyright 1995 by Kenneth Geisshirt (kneth@fatou.ruc.dk)      C
C Dept. of Life Sciences and Chemistry, Roskilde University         C
C                                                                   C
C Constants used by mdreac.                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Natural and transcendental constants
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)

C dt = time step
      REAL*8 dt
      PARAMETER (dt = 5.0D-4)

C ntot = maximum number of particles
      INTEGER ntot
      PARAMETER (ntot=64*1024)

C tau = thermostat relax. time, bufzone = radius of buffer zone
C NLupdate = steps between linked-cell updates
C MaxPair  = maximum number of pairs per particle
C PrnPos   = steps between dumping the positions (in units of NLupdate)
C MaxCell  = maximum number of cells
      REAL*8 tau, bufzone
      INTEGER NLupdate, MaxPair, PrnPoss, MaxCell
      PARAMETER (tau = 0.015D0, bufzone = 4.0D0)
      PARAMETER (NLupdate = 10, MaxPair = 250)
      PARAMETER (PrnPoss = 1000*NLupdate)
C      PARAMETER (MaxCell = 48*48)
      PARAMETER (MaxCell = 80*80)
C      PARAMETER (MaxCell = 18*18)

C Derived parameters
      REAL*8    bufzone2, dt48
      PARAMETER (bufzone2=bufzone**2, dt48=48.D0*dt)

C MaxTask = maximal number of tasks
      INTEGER MaxTask
      PARAMETER (Maxtask=4)

C sfsites = number of sites in each direction used in the 
C calculation of the structure factor
      INTEGER sfsites
      PARAMETER (sfsites = 128)

