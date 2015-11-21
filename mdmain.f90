!
! ####### Program for Molecular dynamics of Self-Assembled Monolayer #########
!
!
PROGRAM MDSAM_CELL
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(PT), POINTER:: q(:)
TYPE(CH), POINTER:: Ngroup(:)
TYPE(CL), POINTER:: cell(:,:)
TYPE(TM):: time
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
TYPE(NM):: NS
REAL*8  :: dt, sum
!
! INTERNAL VARIABLES
REAL   :: time0, time1, time2, secnds
INTEGER:: Nloop_max
time1 = secnds(time0)
!
! Open XYZ file for output
OPEN(UNIT=55, file="result.dat")
WRITE(55,200)
200 FORMAT("# time(fs) instant. temp. Kinetic Eenrgy, Potential Energy, Torsional energy, Carbon VDW energy, Total energy")
!
! Parsing input data and allocate pointer variables
CALL Init(NS, q, Ngroup, cell, time, param, sys, crit, dt)
!
! Initialization
time%Ndump = INT((time%tdump+dt*.5)/dt)
time%Nrest = INT((time%trest+dt*.5)/dt)
time%Nsamp = INT((time%tsamp+dt*.5)/dt)
time%Ncell = INT((time%tcell+dt*.5)/dt)
time%Nloop = 0
time%tnow = 0.0
Nloop_max = INT((time%tmax+dt*.5)/dt)
CALL Vinit(NS, Ngroup, q, param, sys, dt)
!CALL PostXYZBIN(NS, q, time, sys, param)
CALL CellAlloc(NS, q, Ngroup, cell, sys)
!
! Stochastic thermostat loop
IF (param%thermo == 'STOCH ')  THEN
   DO WHILE (time%Nloop < Nloop_max)
      time%tnow = time%tnow + dt
      time%Nloop = time%Nloop + 1   
      !
      ! L.1 Cell decomposition and particle allocation
      IF (MOD(time%Nloop, time%Ncell) == 0) THEN
         CALL CellAlloc(NS, q, Ngroup, cell, sys)
      END IF
      !
      ! L.2 Time integration with constraint dynamics     
      CALL VVerletStoch1(NS, q, Ngroup, param, sys, crit, dt)
      CALL Force(NS, q, Ngroup, cell, param, sys)
      CALL Zforce(NS, q, Ngroup, param, sys)    
      CALL VVerletStoch2(NS, q, Ngroup, param, sys, crit, dt)
      !
      IF (MOD(time%Nloop,time%Nrest) == 0) THEN
         CALL RestartBIN(NS, q, Ngroup, time, sys, dt)
      END IF
      IF (MOD(time%Nloop,time%Nsamp) == 0) THEN
         Call Sample(NS, q, param, time, sys, dt)
      END IF
   END DO
!
! Berendsen thermostat loop
ELSE IF (param%thermo == 'BEREND')  THEN
   DO WHILE (time%tnow <= time%tmax)
      time%tnow = time%tnow + dt
      time%Nloop = time%Nloop + 1   
      IF (MOD(time%Nloop, time%Ncell) == 0) THEN
         CALL CellAlloc(NS, q, Ngroup, cell, sys)
      END IF
      CALL VVerletNotemp1(NS, q, Ngroup, param, sys, crit, dt)
      CALL Force(NS, q, Ngroup, cell, param, sys)
      CALL Zforce(NS, q, Ngroup, param, sys)    
      CALL VVerletBerend2(NS, q, Ngroup, param, sys, crit, dt)
      !
      IF (MOD(time%Nloop,time%Nrest) == 0) THEN
         CALL RestartBIN(NS, q, Ngroup, time, sys, dt)
      END IF
      IF (MOD(time%Nloop,time%Nsamp) == 0) THEN
         Call Sample(NS, q, param, time, sys, dt)
      END IF
   END DO
!
! Microcanonical(NVE) ensemble loop
ELSE IF (param%thermo == 'NOTEMP')  THEN
   DO WHILE (time%tnow < time%tmax+dt*0.5)
      time%tnow = time%tnow + dt
      time%Nloop = time%Nloop + 1   
      IF (MOD(time%Nloop, time%Ncell) == 0) THEN
         CALL CellAlloc(NS, q, Ngroup, cell, sys)
      END IF
      CALL VVerletNotemp1(NS, q, Ngroup, param, sys, crit, dt)
      CALL Force(NS, q, Ngroup, cell, param, sys)
      CALL Zforce(NS, q, Ngroup, param, sys)    
      CALL VVerletNotemp2(NS, q, Ngroup, param, sys, crit, dt)
      !
      !IF (MOD(time%Nloop,time%Ndump) == 0) THEN
      !   CALL PostXYZBIN(NS, q, time, sys, param)
      !END IF
      IF (MOD(time%Nloop,time%Nrest) == 0) THEN
         CALL RestartBIN(NS, q, Ngroup, time, sys, dt)
      END IF
      IF (MOD(time%Nloop,time%Nsamp) == 0) THEN
         Call Sample(NS, q, param, time, sys, dt)
      END IF
   END DO
!
!
ELSE
   PRINT *, "Unknown Thermostat"
END IF
!
! Close XYZ file for output
CLOSE(55)
DEALLOCATE(q, Ngroup, cell)
time2 = secnds(time1)
PRINT *, "Wall time is", time2
!
! Main routine ends here. Below are other subroutines
CONTAINS
!
! ####### Initialization and input data parsing routine ######################
! Unit normalization
! Length: 1. means 1Angstrom = 10E-10 m
! Mass: 1. means 1.6605E-27 kg (a.m.u)
! Energy: 1. means 1 eV = 1.6022E-19 J
! Time: 1. means 10.18 fs = 1.018E-14 sec
! 
SUBROUTINE Init(NS, q, Ngroup, cell, time, param, sys, crit, dt)
USE DataStr
USE REALLOC
IMPLICIT NONE
TYPE(PT), POINTER:: q(:)
TYPE(CH), POINTER:: Ngroup(:)
TYPE(CL), POINTER:: cell(:,:)
TYPE(TM):: time
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
TYPE(NM):: NS
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, l, n, length, NV, NH
CHARACTER(len=5):: dummy
REAL :: rand
!
! Initialization
NS%Ncnstr = 0
!
![[[[[[[[[[[[[[[[[ "input.dat" parsing and data arrangement ]]]]]]]]]]]]]]]]]
OPEN (UNIT=11, file="control.dat")
READ(11,*) dummy
!

! Time parameter
! 1.0 = 10.18fs
READ(11,*) dummy
READ(11,*) time%tmax, time%trest, time%tdump, time%tsamp, time%tcell, dt
time%tmax  = time%tmax  / 10.18
time%trest = time%trest / 10.18
time%tdump = time%tdump / 10.18
time%tcell = time%tcell / 10.18
time%tsamp = time%tsamp / 10.18
dt = dt / 10.18
READ(11,*) dummy
READ(11,*) (sys%box(j), j=1,3)
READ(11,*) dummy
READ(11,*) sys%temp
!
! i.1 Number of all particles
READ(11,*) dummy
READ(11,*) NS%Npt
!
!
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) (param%eps(i,j), j=1,Nparam)
   DO j=1, Nparam
      !
      ! Kcal/mol -> eV
      !param%eps(i,j) = param%eps(i,j)/23.060538
      ! K -> eV
      param%eps(i,j) = param%eps(i,j)*8.6173E-5
   END DO
END DO
!
! p.2 Sigma for LJ(12-6) potential
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) (param%sig(i,j), j=1,Nparam)
END DO
!
! p.3 Square of critical radius for LJ(12-6) potential
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) (param%rc2(i,j), j=1,Nparam)
END DO
DO i=1, Nparam
   DO j=1, Nparam
      param%rc2(i,j) = param%rc2(i,j)**2
   END DO
END DO
!
! p.4 Energy cut-off for LJ(12-6) potential
DO i=1,Nparam
   DO j=1, Nparam
      param%ecut(i,j) = param%eps(i,j)*4.* &
           (param%sig(i,j)**6/param%rc2(i,j)**3 - 1.) &
           *param%sig(i,j)**6/param%rc2(i,j)**3
   END DO
END DO
!
! p.5 Mass for each particle kind
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) param%xm(i)
END DO
!
! p.6 Torsion parameter
READ(11,*) dummy
READ(11,*) (param%kt(j), j=1,3)
!
! p.7 Berendsen thermostat/damping constant
READ(11,*) dummy
READ(11,*) param%thermo
READ(11,*) param%tau, param%alpha
!
! p.8 Iteration limit
READ(11,*) dummy
READ(11,*) crit%Nshake, crit%Errshake
!
CLOSE(11)
!
!#################### READ particle position and id data ######################
!
OPEN (UNIT=15, file="sam.bin", FORM="UNFORMATTED", STATUS = "OLD")
ALLOCATE(q(NS%Npt),Ngroup(NS%Npt))
!
NS%Nchain = 0
DO i=1, NS%Npt
   READ(15) q(i)%xx(1), q(i)%xx(2), q(i)%xx(3), q(i)%id
   IF (q(i)%id == 2) THEN
      NS%Nchain = NS%Nchain + 1
      j = 1
      Ngroup(NS%Nchain)%Ncomp(j) = i
   ELSE IF (q(i)%id == 3) THEN
      j = j + 1
      Ngroup(NS%Nchain)%Ncomp(j) = i
      Ngroup(NS%Nchain)%Ng = j
      NS%Ncnstr = NS%Ncnstr + j*2 - 3
   ELSE
      j = j + 1
      Ngroup(NS%Nchain)%Ncomp(j) = i
   END IF
END DO
! Resize q and Ngroup with exact size
Ngroup => REALLOCCH(Ngroup, NS%Nchain)
CLOSE(15)
!
!XXXXXXXXXXXXXXXX Cell size estimation and memory allocation XXXXXXXXXXXXXXXXXX
!
NS%Nx = INT(sys%box(1)/cellunit)
NS%Ny = INT(sys%box(2)/cellunit)
PRINT *, "Number of Cell along x/y dir. is", NS%Nx, NS%Ny
IF ( NS%Nx < 3 .OR. NS%Ny < 3 ) THEN
   PRINT *, "Problem domain is too small for this program"
   STOP
END IF
ALLOCATE(cell(NS%Nx,NS%Ny))
sys%cellx(1) = sys%box(1)/REAL(NS%Nx)
sys%cellx(2) = sys%box(2)/REAL(NS%Ny)
!
![[[[[[[[[[[[[[[[[[[[[[[[[[[ SURROUNDING CELL SEARCH ]]]]]]]]]]]]]]]]]]]]]]]]]
!
DO i = 1, NS%Nx
   DO j = 1, NS%Ny
      l = 1
      !
      ! Top cells x3
      NV = j + 1
      IF (NV > NS%Ny) NV = 1
      DO n = -1, 1
         NH = i+n
         IF (NH > NS%Nx) NH = 1
         IF (NH < 1) NH = NS%Nx
         cell(i,j)%Nx(l) = NH
         cell(i,j)%Ny(l) = NV
         l = l + 1
      END DO
      !
      ! Right cell x1
      NH = i + 1
      IF (NH > NS%Nx) NH = 1
      cell(i,j)%Nx(l) = NH
      cell(i,j)%Ny(l) = j
   END DO
END DO
!
RETURN
END SUBROUTINE Init
!
END PROGRAM MDSAM_CELL
!
! ###################### Veolcity Initialization ##############################
! 
! With given temperature, initial velocity is given to each particle. 
! But distrubtion follows Gaussian(Normalized) random distribution
! 
SUBROUTINE Vinit(NS, Ngroup, q, param, sys, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(CH)::Ngroup(NS%Nchain)
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
REAL*8::dt, xv(NS%Npt,3)
!
INTEGER:: i, j, k, id, Nr
REAL*8:: xm, lambda, T0, kb
!
kb = 8.617343E-5
sys%mv2 = 0.0
IF (param%thermo == 'STOCH ') THEN
   Nr = 0
ELSE 
   Nr = 4
END IF
DO i=1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   DO j=1, 3
      xv(i,j) = fluct(SQRT(sys%temp))
      sys%mv2 = sys%mv2 + xm*xv(i,j)**2 
   END DO
END DO
T0 = sys%mv2/kb/real(3*NS%Npt-NS%Ncnstr - Nr)
lambda = SQRT(sys%temp/T0)
sys%mv2 = 0.0
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   DO j = 1,3
      q(i)%xv(j) = xv(i,j)*lambda
      sys%mv2 = sys%mv2 + xm*q(i)%xv(j)**2
   END DO
END DO
RETURN
END SUBROUTINE Vinit

