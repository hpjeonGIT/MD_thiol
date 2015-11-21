!
!
! ####################### Binary restart file generator  ######################
! At this routine, restart file is generated. The format is same to input
! file. To rerun, just change the name "restart.dat" into "input.dat"
!
SUBROUTINE RestartBIN(NS, q, Ngroup, time, sys, dt)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(TM):: time
TYPE(ST):: sys
REAL*8  :: dt
!
!
INTEGER:: i, j, Nfreq
CHARACTER(LEN=10):: NUMBER, FILENAME
DATA NUMBER/'0123456789'/
FILENAME = "REST00.bin"
!
Nfreq = INT(time%Nloop/time%Nrest)
IF(INT(Nfreq/10) < 1 ) THEN
   FILENAME(5:5) = NUMBER(1:1)
   FILENAME(6:6) = NUMBER(Nfreq+1:Nfreq+1)
ELSE
   FILENAME(5:5) = NUMBER(INT(Nfreq/10)+1:INT(Nfreq/10)+1)
   FILENAME(6:6) = NUMBER(Nfreq-INT(Nfreq/10)*10+1:Nfreq-INT(Nfreq/10)*10+1)
END IF
OPEN(UNIT=30,FILE=FILENAME, FORM="UNFORMATTED")
DO i=1,NS%Npt
   WRITE(30) (q(i)%xx(j), j=1,3), q(i)%id
END DO
!
CLOSE(30)
RETURN
END SUBROUTINE RestartBIN
!
! ################### Subroutine for binary data output #######################
!
! Dump intermediate results into binary xyz format.
!
SUBROUTINE PostXYZBIN(NS, q, time, sys, param)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(TM):: time
TYPE(ST):: sys
TYPE(PM):: param
!
! INTERNAL VARIABLES
INTEGER:: i, j, id, Nr
REAL*8:: kb
CHARACTER(len=3):: atom(Nparam), dummy
CHARACTER(len=10):: header1, header2
atom =(/"C  ","S  ", "C  "/)
kb = 8.617343E-5 ! Boltzmann constant in eV/K
header1 = "frame =   "
header2 = "energy =  "
IF (param%thermo == 'STOCH ') THEN
   Nr = 0
ELSE
   Nr = 4
END IF
!
WRITE(25) NS%Npt
WRITE(25) header1, INT(time%Nloop/time%Ndump), header2, sys%Epot+(sys%mv2/2.)
DO i=1, NS%Npt
   id = q(i)%id
   dummy = atom(id)
   WRITE(25) dummy, (q(i)%xx(j), j=1, 3)
END DO
!
!
WRITE(*,40) sys%mv2/REAL(NS%Npt*3- NS%Ncnstr- Nr)/kb, time%Nloop, &
     time%tnow*10.18, time%tmax*10.18
40 FORMAT("Temp.=", F8.1, "K  at ", I8, "th loop with", F8.1,"/",F8.1, "fs")
RETURN
!
END SUBROUTINE PostXYZBIN
