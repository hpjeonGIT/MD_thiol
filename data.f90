MODULE DataStr
IMPLICIT NONE
INTEGER, PARAMETER:: Nparam = 3, Nlong = 19, Ncell = 20000
REAL*8, PARAMETER :: cellunit = 30.
!
! Particle data type
! xx = current position
! xv = current velocity
! ff = current force
! id = particle type
TYPE PT
   REAL*8:: xx(3), xv(3), ff(3)
   INTEGER:: id
END TYPE PT
!
! Chain data type
! Ng = Number components on each chain
! Ncomp = Particle number along chain
TYPE CH
   INTEGER:: Ng, Ncomp(Nlong)
END TYPE CH
!
! Cell data type
! id = id of each chain in the cell
! Nchain = Number of chains in the cell
! Nx/y/z = Array of surrounding cells
TYPE CL
   INTEGER:: id(Ncell), Nchain, Nx(4), Ny(4)
END TYPE CL
!
! Time data type
! Nloop = Number of loops
! Ndump = Number of dumps
! Nsamp = Number of samplings
! tmax = Maximum simulation time
! tdump = xyz file dump time
! tsamp = sampling time
! tnow = Current time
TYPE TM
   INTEGER:: Nloop, Ndump, Nsamp, Nrest, Ncell
   REAL*8:: tmax, tdump, tsamp, trest, tnow, tcell
END TYPE TM
!
! Parameter data type
! eps = Epsilon for LJ(12-6) potential
! sigma = Sigma for LJ(12-6) potential
! rc2 = Square of critical radius
! ecut = Energy cut-off
! xm = Mass of particle
! kt = Torsion stiffness
! tau = Berendsen thermostat parameter
TYPE PM
   REAL*8:: eps(Nparam,Nparam), sig(Nparam, Nparam), rc2(Nparam, Nparam), &
        ecut(Nparam, Nparam), xm(Nparam), kt(3), tau, alpha
   CHARACTER*6:: thermo
END TYPE PM
!
! System variable data type
! box = Size of simulation cell
! temp = Temperature of (NVT) simulation
! mv2 = Sum of mv^2(=twice of kinetic energy)
! Epot = Potential energy of system
! Etor = Torsional energy of system
! Ecc  = CHx - CHx interaction energy
! cellx = length of cell
TYPE ST
   REAL*8:: box(3), temp, mv2, Epot, Etor, Ecc, cellx(2)
END TYPE ST
!
! Iteration limit data type
! Nshake = Maximum number of SHAKE iterations
! Errshake = Error limit of SHAKE iterations
TYPE LM
   INTEGER:: Nshake
   REAL*8:: Errshake
END TYPE LM
!
! Various number type
! Main objective is transfer specific numbers into print routine
! Ncnstr = Number of all constraints(bond,angle)
! Other sampling variables
! Nx/y = Number of cells along x/y dir.
TYPE NM
   INTEGER:: Npt, Nchain, Ncnstr, Nx, Ny
   REAL*8:: temp, mv2
END TYPE NM
END MODULE DataStr
!
MODULE REALLOC
USE DATASTR
CONTAINS
  !
  ! for array q, resize with n
  FUNCTION REALLOCCH(q,n)
    TYPE(CH), POINTER:: q(:), REALLOCCH(:)
    INTEGER, intent(in):: n
    INTEGER:: ierr
    ALLOCATE(REALLOCCH(1:n), STAT = ierr)
    IF(ierr /=0) STOP "allocate error"
    REALLOCCH(1:n) = q(1:n)
    DEALLOCATE(q)
  END FUNCTION REALLOCCH
END MODULE REALLOC
