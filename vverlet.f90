!
! ####### Verlet time integration routine with Thermal Noise #########
! 
! Velocity verlet with RATTLE - stochastic thermostat implemented
!
SUBROUTINE VVerletStoch1(NS, q, Ngroup, param, sys, crit, dt)
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
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Ng, Na, Nmax, Ncount, Nr
INTEGER, ALLOCATABLE:: Nbond(:,:)
REAL*8:: x, distance, Errmax, gamma, kb, sigma, s(3)
REAL*8, ALLOCATABLE:: xm(:), mu(:), d2(:), qt(:,:), rr(:,:)
!
! Initialization
kb = 8.617343E-5 ! Boltzmann constant in eV/K
Nr = 4           ! 2 for PBC, 2 for vacuum along z-dir.
sigma = 1.       ! Deviation for Gaussian distribution
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
!
! RATTLE 1st stage
DO n=1, NS%Nchain   
   Ng = Ngroup(n)%Ng ! Number of particles in chain
   Na = 2*Ng - 3     ! Number of bonds in chain
   ALLOCATE(xm(Ng), mu(Na), Nbond(Na,2), d2(Na), qt(Ng,3), rr(Na,3))
   d2(1) = 3.2761    ! Square of S-CH2 length
   d2(2) = 6.5280464 ! Square of S-(CH2)-CH2 length
   DO i = 2, Ng-2
      d2(i*2-1) = 2.362369   ! Square of CH2-CH2 length
      d2(i*2) = 6.301887247  ! Square of CH2-(CH2)-CH2 length
   END DO
   d2(Na) = 2.362369
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, Ng-2
      id1 = Ngroup(n)%Ncomp(i)
      id2 = Ngroup(n)%Ncomp(i+1)
      id3 = Ngroup(n)%Ncomp(i+2)
      xm(i) = param%xm(q(id1)%id)
      Nbond(2*i-1,1) = i
      Nbond(2*i-1,2) = i+1
      Nbond(2*i,1) = i
      Nbond(2*i,2) = i+2
      DO k=1,3
!         eta = fluct(sigma)
         qt(i,k) = q(id1)%xv(k)*(1.-0.5*param%alpha*dt/xm(i)) + &
              0.5*dt*q(id1)%ff(k)/xm(i)
         x = q(id2)%xx(k) - q(id1)%xx(k)
         rr(2*i-1,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond vector
         x = q(id3)%xx(k) - q(id1)%xx(k)
         rr(2*i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond vector'
      END DO
      mu(2*i-1) = param%xm(q(id1)%id)*param%xm(q(id2)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id2)%id))
      mu(2*i) = param%xm(q(id1)%id)*param%xm(q(id3)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id3)%id))         
   END DO
   DO i=Ng-1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      xm(i) = param%xm(q(id1)%id)
      DO k=1,3
         qt(i,k) = q(id1)%xv(k)*(1.-0.5*param%alpha*dt/xm(i)) + &
              0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   DO k=1,3
      x = q(Ngroup(n)%Ncomp(Ng))%xx(k) -q(Ngroup(n)%Ncomp(Ng-1))%xx(k) 
      rr(Na,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! bond vector
   END DO
   mu(Na) = xm(Ng)*xm(Ng-1)/(xm(Ng) + xm(Ng-1))
   Nbond(Na,1) = Ng-1
   Nbond(Na,2) = Ng
   !
   ! Iteration loop for gamma
   DO m=1,Nmax
      Ncount = 0
      DO i=1,Na
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            s(k) = rr(i,k) + dt*(qt(id2,k) - qt(id1,k))
         END DO
         distance = s(1)**2 + s(2)**2 + s(3)**2 - d2(i)
         gamma = 0.5*mu(i)*distance/dt/ &
              (s(1)*rr(i,1)+s(2)*rr(i,2)+s(3)*rr(i,3))
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1         
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == Na) GOTO 10
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 1st stage"
   STOP
   !
   ! Position update. Intermediate velocity is available now.
10 DO i=1, Ng
      DO k=1,3
         x = q(Ngroup(n)%Ncomp(i))%xx(k) + dt*qt(i,k)
         q(Ngroup(n)%Ncomp(i))%xx(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         q(Ngroup(n)%Ncomp(i))%xv(k) = qt(i,k)
      END DO
   END DO
   !
   !
   !print *, m, "1st stage iteration"
   DEALLOCATE(xm, mu, Nbond, d2, qt, rr)
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch1
!
!
SUBROUTINE VVerletStoch2(NS, q, Ngroup, param, sys, crit, dt)
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
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Ng, Na, Nmax, Ncount, Nr
INTEGER, ALLOCATABLE:: Nbond(:,:)
REAL*8:: x, distance, Errmax, gamma, kb, T0, eta, beta, sigma, s(3)
REAL*8, ALLOCATABLE:: xm(:), mu(:), d2(:), qt(:,:), rr(:,:)
!
! Initialization
kb = 8.617343E-5 ! Boltzmann constant in eV/K
Nr = 4           ! 2 for PBC, 2 for vacuum along z-dir.
sigma = 1.       ! Deviation for Gaussian distribution
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
beta = SQRT(2.*param%alpha*kb*sys%temp/dt)
!
! RATTLE 1st stage
DO n=1, NS%Nchain   
   Ng = Ngroup(n)%Ng ! Number of particles in chain
   Na = 2*Ng - 3     ! Number of bonds in chain
   ALLOCATE(xm(Ng), mu(Na), Nbond(Na,2), d2(Na), qt(Ng,3), rr(Na,3))
   d2(1) = 3.2761    ! Square of S-CH2 length
   d2(2) = 6.5280464 ! Square of S-(CH2)-CH2 length
   DO i = 2, Ng-2
      d2(i*2-1) = 2.362369   ! Square of CH2-CH2 length
      d2(i*2) = 6.301887247  ! Square of CH2-(CH2)-CH2 length
   END DO
   d2(Na) = 2.362369
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, Ng-2
      id1 = Ngroup(n)%Ncomp(i)
      id2 = Ngroup(n)%Ncomp(i+1)
      id3 = Ngroup(n)%Ncomp(i+2)
      xm(i) = param%xm(q(id1)%id)
      Nbond(2*i-1,1) = i
      Nbond(2*i-1,2) = i+1
      Nbond(2*i,1) = i
      Nbond(2*i,2) = i+2
      mu(2*i-1) = param%xm(q(id1)%id)*param%xm(q(id2)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id2)%id))
      mu(2*i) = param%xm(q(id1)%id)*param%xm(q(id3)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id3)%id))         
   END DO
   DO i=Ng-1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      xm(i) = param%xm(q(id1)%id)
   END DO
   mu(Na) = xm(Ng)*xm(Ng-1)/(xm(Ng) + xm(Ng-1))
   Nbond(Na,1) = Ng-1
   Nbond(Na,2) = Ng
   DO i=1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      DO k=1,3
         eta = fluct(sigma)
         q(id1)%ff(k) = q(id1)%ff(k) + eta*beta
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   !
   ! Iteration loop for gamma
   DO m=1,Nmax
      Ncount = 0
      DO i=1,Na
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            x = q(Ngroup(n)%Ncomp(id2))%xx(k) - q(Ngroup(n)%Ncomp(id1))%xx(k)
            rr(i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond distance
         END DO
         distance = rr(i,1)*(qt(id2,1)-qt(id1,1)) + &
              rr(i,2)*(qt(id2,2)-qt(id1,2)) + rr(i,3)*(qt(id2,3)-qt(id1,3))
         gamma = mu(i)*distance/d2(i)
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == Na) GOTO 20
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 2nd stage"
   STOP
   !
   ! Velocity update (phew! ~~~)
20 DO i=1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      DO k=1,3
         q(id1)%xv(k) = qt(i,k)/(1.+0.5*param%alpha*dt/xm(i))
         sys%mv2 = sys%mv2 + xm(i)*q(id1)%xv(k)**2
      END DO
   END DO
   !
   !
   !print *, m, "iteration number for 2nd stage"
   DEALLOCATE(xm, mu, Nbond, d2, qt, rr)
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch2
!
!
! ####### Verlet time integration routine with Berendsen Thermostat #########
! 
! This routine updates positions of all atoms at every iteration. Using
! the forces, acceleration will be calculated. Then using position verlet
! algorithm, positions are updated. Using intermediate calculations, velocities
! and kinetic energy are available. 
!
! For NVT simulation, temperature should be kept as given one. Using Berendsen
! thermostat, intermediate proportional constant will be estimated and all 
! velocities will be modified with the constant. After modifying velocities
! positions are also modified.
!
! First stage of Verlet scheme comes from no thermostat routine. Then velocity
! update is performed in this routine.
!
SUBROUTINE VVerletBerend2(NS, q, Ngroup, param, sys, crit, dt)
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
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Ng, Na, Nmax, Ncount, Nr
INTEGER, ALLOCATABLE:: Nbond(:,:)
REAL*8:: x, distance, Errmax, gamma,  s(3), lambda, T0, kb
REAL*8, ALLOCATABLE:: xm(:), mu(:), d2(:), qt(:,:), rr(:,:)
!
! 
kb = 8.617343E-5 ! Boltzmann constant in eV/K
Nr = 4           ! 2 for PBC, 2 for vacuum along z-dir.
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
DO n=1, NS%Nchain   
   Ng = Ngroup(n)%Ng ! Number of particles in chain
   Na = 2*Ng - 3     ! Number of bonds in chain
   ALLOCATE(xm(Ng), mu(Na), d2(Na), Nbond(Na,2), qt(Ng,3), rr(Na,3))
   d2(1) = 3.2761    ! Square of S-CH2 length
   d2(2) = 6.5280464 ! Square of S-(CH2)-CH2 length
   DO i = 2, Ng-2
      d2(i*2-1) = 2.362369   ! Square of CH2-CH2 length
      d2(i*2) = 6.301887247  ! Square of CH2-(CH2)-CH2 length
   END DO
   d2(Na) = 2.362369
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, Ng-2
      id1 = Ngroup(n)%Ncomp(i)
      id2 = Ngroup(n)%Ncomp(i+1)
      id3 = Ngroup(n)%Ncomp(i+2)
      xm(i) = param%xm(q(id1)%id)
      Nbond(2*i-1,1) = i
      Nbond(2*i-1,2) = i+1
      Nbond(2*i,1) = i
      Nbond(2*i,2) = i+2
      mu(2*i-1) = param%xm(q(id1)%id)*param%xm(q(id2)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id2)%id))
      mu(2*i) = param%xm(q(id1)%id)*param%xm(q(id3)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id3)%id))         
   END DO
   DO i=Ng-1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      xm(i) = param%xm(q(id1)%id)
   END DO
   mu(Na) = xm(Ng)*xm(Ng-1)/(xm(Ng) + xm(Ng-1))
   Nbond(Na,1) = Ng-1
   Nbond(Na,2) = Ng
   !
   ! Iteration loop for gamma
   !print *, m, "iteration number for 1st stage"
   !
   ! 1st stage of RATTLE is over. Now 2nd stage begins
   DO i=1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      DO k=1,3
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   !goto 20
   DO m=1,Nmax
      Ncount = 0
      DO i=1,Na
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            x = q(Ngroup(n)%Ncomp(id2))%xx(k) - q(Ngroup(n)%Ncomp(id1))%xx(k)
            rr(i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond distance
         END DO
         distance = rr(i,1)*(qt(id2,1)-qt(id1,1)) + &
              rr(i,2)*(qt(id2,2)-qt(id1,2)) + rr(i,3)*(qt(id2,3)-qt(id1,3))
         gamma = distance*mu(i)/d2(i)
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1
         !          print *, i, gamma
         DO k=1,3
             qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
             qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
          END DO
       END DO
       IF (Ncount == Na) GOTO 20
    END DO
    PRINT *, "RATTLE is not converged at", n,"th chain for 2nd stage"
    STOP
    !
    ! Velocity update (phew! ~~~)
20  DO i=1, Ng
       id1 = Ngroup(n)%Ncomp(i)
       DO k=1,3
          q(id1)%xv(k) = qt(i,k)
          sys%mv2 = sys%mv2 + xm(i)*q(id1)%xv(k)**2
       END DO
    END DO
    !
    !print *, m, "iteration number for 2nd stage"
    DEALLOCATE(xm, mu, Nbond, d2, qt, rr)
END DO
!
!
T0 = sys%mv2/kb/REAL(3*NS%Npt-NS%Ncnstr - Nr)
lambda = SQRT(1.+ dt*(sys%temp/T0 - 1.)/param%tau)
sys%mv2 = 0.0
!
! Velocity rescaling
DO i = 1, NS%Npt
   id1 = q(i)%id
   x = param%xm(id1)
   DO j = 1,3
      !
      !
      q(i)%xv(j) = q(i)%xv(j)*lambda
      sys%mv2 = sys%mv2 + x*q(i)%xv(j)**2
   END DO   
ENDDO
RETURN
!
END SUBROUTINE VVerletBerend2
!!
! ####### Pure Verlet time integration routine with No Thermostat #########
! 
! Velocity Verlet algorithm with RATTLE - no thermostat implemented
!
SUBROUTINE VVerletNotemp1(NS, q, Ngroup, param, sys, crit, dt)
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
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Ng, Na, Nmax, Ncount
INTEGER, ALLOCATABLE:: Nbond(:,:)
REAL*8:: x, distance, Errmax, gamma,  s(3)
REAL*8, ALLOCATABLE:: xm(:), mu(:), d2(:), qt(:,:), rr(:,:)
!
! Initialization
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
!
! RATTLE 1st stage
DO n=1, NS%Nchain   
   Ng = Ngroup(n)%Ng ! Number of particles in chain
   Na = 2*Ng - 3     ! Number of bonds in chain
   ALLOCATE(xm(Ng), mu(Na), Nbond(Na,2), d2(Na), qt(Ng,3), rr(Na,3))
   d2(1) = 3.2761    ! Square of S-CH2 length
   d2(2) = 6.5280464 ! Square of S-(CH2)-CH2 length
   DO i = 2, Ng-2
      d2(i*2-1) = 2.362369   ! Square of CH2-CH2 length
      d2(i*2) = 6.301887247  ! Square of CH2-(CH2)-CH2 length
   END DO
   d2(Na) = 2.362369
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, Ng-2
      id1 = Ngroup(n)%Ncomp(i)
      id2 = Ngroup(n)%Ncomp(i+1)
      id3 = Ngroup(n)%Ncomp(i+2)
      xm(i) = param%xm(q(id1)%id)
      Nbond(2*i-1,1) = i
      Nbond(2*i-1,2) = i+1
      Nbond(2*i,1) = i
      Nbond(2*i,2) = i+2
      DO k=1,3
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
         x = q(id2)%xx(k) - q(id1)%xx(k)
         rr(2*i-1,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond vector
         x = q(id3)%xx(k) - q(id1)%xx(k)
         rr(2*i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond vector'
      END DO
      mu(2*i-1) = param%xm(q(id1)%id)*param%xm(q(id2)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id2)%id))
      mu(2*i) = param%xm(q(id1)%id)*param%xm(q(id3)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id3)%id))         
   END DO
   DO i=Ng-1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      xm(i) = param%xm(q(id1)%id)
      DO k=1,3
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   DO k=1,3
      x = q(Ngroup(n)%Ncomp(Ng))%xx(k) - q(Ngroup(n)%Ncomp(Ng-1))%xx(k) 
      rr(Na,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! bond vector
   END DO
   mu(Na) = xm(Ng)*xm(Ng-1)/(xm(Ng) + xm(Ng-1))
   Nbond(Na,1) = Ng-1
   Nbond(Na,2) = Ng
   !
   ! Iteration loop for gamma
   DO m=1,Nmax
      Ncount = 0
      DO i=1,Na
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
           x = rr(i,k) + dt*(qt(id2,k) - qt(id1,k))
           s(k) = x - sys%box(k)*ANINT(x/sys%box(k)) 
         END DO
         distance = s(1)**2 + s(2)**2 + s(3)**2 - d2(i)
         gamma = 0.5*distance*mu(i)/dt/ &
              (s(1)*rr(i,1)+s(2)*rr(i,2)+s(3)*rr(i,3))
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1        
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == Na) GOTO 10
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 1st stage"
   STOP
   !
   ! Position update. Intermediate velocity is available now.
10 DO i=1, Ng
      DO k=1,3
         x = q(Ngroup(n)%Ncomp(i))%xx(k) + dt*qt(i,k)
         q(Ngroup(n)%Ncomp(i))%xx(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         q(Ngroup(n)%Ncomp(i))%xv(k) = qt(i,k)
      END DO
   END DO
   !print *, m, "1st stage iteration" 
   DEALLOCATE(xm, mu, Nbond, d2, qt, rr)
END DO
RETURN
END SUBROUTINE VVerletNotemp1
!
SUBROUTINE VVerletNotemp2(NS, q, Ngroup, param, sys, crit, dt)
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
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Ng, Na, Nmax, Ncount
INTEGER, ALLOCATABLE:: Nbond(:,:)
REAL*8:: x, distance, Errmax, gamma,  ss(3)
REAL*8, ALLOCATABLE:: xm(:), mu(:), d2(:), qt(:,:), rr(:,:)
!
! 
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
DO n=1, NS%Nchain   
   Ng = Ngroup(n)%Ng ! Number of particles in chain
   Na = 2*Ng - 3     ! Number of bonds in chain
   ALLOCATE(xm(Ng), mu(Na), d2(Na), Nbond(Na,2), qt(Ng,3), rr(Na,3))
   d2(1) = 3.2761    ! Square of S-CH2 length
   d2(2) = 6.5280464 ! Square of S-(CH2)-CH2 length
   DO i = 2, Ng-2
      d2(i*2-1) = 2.362369   ! Square of CH2-CH2 length
      d2(i*2) = 6.301887247  ! Square of CH2-(CH2)-CH2 length
   END DO
   d2(Na) = 2.362369
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, Ng-2
      id1 = Ngroup(n)%Ncomp(i)
      id2 = Ngroup(n)%Ncomp(i+1)
      id3 = Ngroup(n)%Ncomp(i+2)
      xm(i) = param%xm(q(id1)%id)
      Nbond(2*i-1,1) = i
      Nbond(2*i-1,2) = i+1
      Nbond(2*i,1) = i
      Nbond(2*i,2) = i+2
      mu(2*i-1) = param%xm(q(id1)%id)*param%xm(q(id2)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id2)%id))
      mu(2*i) = param%xm(q(id1)%id)*param%xm(q(id3)%id)/ &
           (param%xm(q(id1)%id)+param%xm(q(id3)%id))         
   END DO
   DO i=Ng-1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      xm(i) = param%xm(q(id1)%id)
   END DO
   mu(Na) = xm(Ng)*xm(Ng-1)/(xm(Ng) + xm(Ng-1))
   Nbond(Na,1) = Ng-1
   Nbond(Na,2) = Ng
   !
   ! Iteration loop for gamma
   !print *, m, "iteration number for 1st stage"
   !
   ! 1st stage of RATTLE is over. Now 2nd stage begins
   DO i=1, Ng
      id1 = Ngroup(n)%Ncomp(i)
      DO k=1,3
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   !goto 20
   DO m=1,Nmax
      Ncount = 0
      DO i=1,Na
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            x = q(Ngroup(n)%Ncomp(id2))%xx(k) - q(Ngroup(n)%Ncomp(id1))%xx(k)
            rr(i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond distance
         END DO
         distance = rr(i,1)*(qt(id2,1)-qt(id1,1)) + &
              rr(i,2)*(qt(id2,2)-qt(id1,2)) + rr(i,3)*(qt(id2,3)-qt(id1,3))
         gamma = distance*mu(i)/d2(i)
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1
         !          print *, i, gamma
         DO k=1,3
             qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
             qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
          END DO
       END DO
       IF (Ncount == Na) GOTO 20
    END DO
    PRINT *, "RATTLE is not converged at", n,"th chain for 2nd stage"
    STOP
    !
    ! Velocity update (phew! ~~~)
20  DO i=1, Ng
       id1 = Ngroup(n)%Ncomp(i)
       DO k=1,3
          q(id1)%xv(k) = qt(i,k)
          sys%mv2 = sys%mv2 + xm(i)*q(id1)%xv(k)**2
       END DO
    END DO
    !
    !
    !print *, m, "iteration number for 2nd stage"
    DEALLOCATE(xm, mu, Nbond, d2, qt, rr)
END DO
RETURN
!
END SUBROUTINE VVerletNotemp2
!
! ############## Random Gaussian(Normal) Distribution Function ################
! 
! For stochastic thermostat, fluctuation dissipation theorem is implemented.
! Basically, random number generator which follows Gaussian distribution is
! needed to implement this thermal noise force.
! Random number is generated using FORTRAN intrinsic fucntion - RANDOM_SEED
! and RANDOM_NUMBER. But during the implementation, it is found that those
! intrinsic functions may not work well under false seeding - to use this
! routine on new machine or new compiler, please make sure that the 
! distribution follows zero-mean with suitable deviation.
!
! This function provides a random number with zero-mean and deviation x 
! along Gaussian distribution.
! <p> = 0.0
! <p**2> = x**2
! 
FUNCTION fluct(x)
  IMPLICIT NONE
  REAL*8:: fluct, x, r, v1, v2
  REAL*8:: rand1, rand2, ran2
  REAL:: ranf
  !
  ! Initialization
  r=1.
  DO WHILE (r.ge.1.)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.*rand1 - 1.
     v2 = 2.*rand2 - 1.
     r = v1*v1+v2*v2
  END DO
  fluct = v1*SQRT(-2.*log(r)/r)*x
  RETURN
END FUNCTION fluct
!
! ############### personal ANINT function - not used #####################
! 
!FUNCTION BJINT(x)
!  IMPLICIT NONE
!  REAL*8:: x, BJINT
!  IF (x >=0) THEN
!     BJINT = AINT(x+0.5)
!  ELSE
!     BJINT = AINT(x-0.5)
!  END IF
!  !
!  RETURN
!END FUNCTION BJINT
