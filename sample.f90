!
SUBROUTINE Sample(NS, q, param, time, sys, dt)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(PM):: param
TYPE(TM):: time
TYPE(ST):: sys
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, id, Nr
REAL*8 :: kb, T0
!
! Evaluate the kind of thermostat
IF (param%thermo == 'STOCH ') THEN
   Nr = 0
ELSE
   Nr = 4
END IF
!
!
kb = 8.617343E-5 
!
!
WRITE(55,100) time%tnow*10.18, sys%mv2/kb/real(3*NS%Npt-NS%Ncnstr-Nr), &
     sys%mv2*.5, sys%Epot, sys%Etor, sys%Ecc, sys%mv2*.5 + sys%Etor + sys%Epot
100 FORMAT(7(ES14.6, 1x))
RETURN
!
END SUBROUTINE Sample
!
!##################### Allocating chains into cells ###########################
!
SUBROUTINE CellAlloc(NS, q, Ngroup, cell, sys)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(CL):: cell(NS%Nx,NS%Ny)
TYPE(ST):: sys
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, idx, idy
!
! Initialize the number of atoms in the cell
DO i=1, NS%Nx
   DO j=1, NS%Ny
      cell(i,j)%Nchain = 0
   END DO
END DO
!
! Allocating chains into each cell
DO i = 1, NS%Nchain
   j = INT(Ngroup(i)%Ng/2)
   k = Ngroup(i)%Ncomp(j)
   idx = INT((q(k)%xx(1)+sys%box(1)*.5)/sys%cellx(1)) + 1
   idy = INT((q(k)%xx(2)+sys%box(2)*.5)/sys%cellx(2)) + 1
   IF (idx > NS%Nx) PRINT *, "X-dir allocation error"
   IF (idy > NS%Ny) PRINT *, "Y-dir allocation error"
   !   print *, idx, idy, idz, i
   cell(idx,idy)%Nchain = cell(idx,idy)%Nchain + 1
   IF (cell(idx,idy)%Nchain > Ncell) THEN
      PRINT *, "cell allocation error"
      STOP
   END IF
   cell(idx,idy)%id(cell(idx,idy)%Nchain) = i
END DO
RETURN
END SUBROUTINE CellAlloc
