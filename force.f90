!
! ####### Routine for force calculation using Lennard-Jones potential #########
!
! Usually molecular potentials are compsed of bond extension(vibration),
! angle bending, torsion, Van Der Waals interaction, and electrostatics.
! In this program, only torsion and VDW will be considered.
! Assuming alkanechain is rigid - no bond strain, no angle strain - and 
! neglecting long range interaction(electrostatic force), torsion and VDW
! is implemented. 
! 
! For alkanechains, carbon atoms are fully saturated and we can assume that
! there is no electo-charge which will case long-range interactions.
! But if the chain is long, then it is reported that gauche defects will be 
! deteced at high temperature(~300K). Therefore implementing torsion of 
! alkanechains will be efficient way to simulate this gauche defect.
!
! For VDW interactions, conventional Lennard-Jones(12-6) potential is 
! implemented. For saturated carbons(CH2, CH3), conventional potential 
! libraries like MM3, AMBER are available and could be embedded.
!
SUBROUTINE Force(NS, q, Ngroup, cell, param, sys)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(CL):: cell(NS%Nx,NS%Ny)
TYPE(PM):: param
TYPE(ST):: sys
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, ii, jj, kk, l, m, n, nn, id_c, id_a, ix, iy, Nx, Ny
REAL*8:: r2, r2i, r6i, lj, xr(3), epsilon, sigma, ecut, rc2, x, y, z, a,rc, &
     pi, De, Dex, Dey, Re, Rex, Rey, beta, mu, C, rho(NS%Nchain), xp, Eint
!
! I.1 Pi declaration
pi = ATAN(1.0)*4.0
!
! I.2 Initialize potential energy
sys%Epot = 0.0
sys%Ecc  = 0.0
a = 2.884
mu = 5.*a
rc = mu
rho = 0.
!
! I.3 Initialize force
DO i = 1, NS%Npt
   DO j = 1, 3
      q(i)%ff(j) = 0.0
   ENDDO
ENDDO
!goto 30
!
! <<< Lennar-Jones 12-6 potential and force for molecules in the same cell >>>
DO ix = 1, NS%Nx
   DO iy = 1, NS%Ny
      !
      ! Particle interaction within same cell
      DO n = 1, cell(ix,iy)%Nchain-1
         i = cell(ix,iy)%id(n)
         DO ii=1, Ngroup(i)%Ng            
            k = Ngroup(i)%Ncomp(ii)
            id_c = q(k)%id
            DO m = n+1, cell(ix,iy)%Nchain
               j = cell(ix,iy)%id(m)
               DO jj=1, Ngroup(j)%Ng
                  kk = Ngroup(j)%Ncomp(jj)
                  id_a = q(kk)%id
                  rc2 = param%rc2(id_c,id_a)
                  !
                  r2 = 0.0
                  DO l=1,3
                     x = q(k)%xx(l) - q(kk)%xx(l)
                     xr(l) = x - sys%box(l)*ANINT(x/sys%box(l))
                     r2 = r2 + xr(l)**2
                  END DO
                  IF ((id_c == 2) .AND. (id_a == 2)) THEN
                     !
                     ! Thiol potential density
                     IF (r2 < rc**2) THEN
                        lj = (1.-COS(pi*(SQRT(r2)-rc)/mu)**2)
                        rho(i) = rho(i) + lj
                        rho(j) = rho(j) + lj
                     END IF
                  ELSE
                     !
                     ! LJ potential
                     IF (r2 < rc2) THEN
                        epsilon = param%eps(id_c,id_a)
                        sigma = param%sig(id_c,id_a)
                        ecut = param%ecut(id_c,id_a)            
                        r2i= 1./r2
                        r6i = sigma**6*r2i**3
                        lj = 48.*r2i*r6i*(r6i - .5)
                        DO l = 1,3
                           q(k)%ff(l)  =  q(k)%ff(l) + lj*xr(l)*epsilon
                           q(kk)%ff(l) = q(kk)%ff(l) - lj*xr(l)*epsilon
                        END DO
                        Eint = 4.*r6i*(r6i - 1.)*epsilon - ecut
                        sys%Epot = sys%Epot + Eint
                        IF ((id_c /=2).AND.(id_a /=2)) sys%Ecc = sys%Ecc + Eint
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO
      !
      !  For molecules in neighboring cell 
      DO n = 1, cell(ix,iy)%Nchain
         i = cell(ix,iy)%id(n)
         DO ii=1, Ngroup(i)%Ng            
            k = Ngroup(i)%Ncomp(ii)
            id_c = q(k)%id
            DO nn=1, 4
               Nx = cell(ix,iy)%Nx(nn)
               Ny = cell(ix,iy)%Ny(nn)
               DO m = 1, cell(Nx,Ny)%Nchain
                  j = cell(Nx,Ny)%id(m)
                  DO jj=1, Ngroup(j)%Ng
                     kk = Ngroup(j)%Ncomp(jj)
                     id_a = q(kk)%id
                     r2 = 0.0
                     rc2 = param%rc2(id_c,id_a)
                     DO l=1,3
                        x = q(k)%xx(l) - q(kk)%xx(l)
                        xr(l) = x - sys%box(l)*ANINT(x/sys%box(l))
                        r2 = r2 + xr(l)**2
                     END DO
                     IF ((id_c == 2) .AND. (id_a == 2)) THEN
                        !
                        ! Thiol potential density
                        IF (r2 < rc**2) THEN
                           lj = (1.-COS(pi*(SQRT(r2)-rc)/mu)**2)
                           rho(i) = rho(i) + lj
                           rho(j) = rho(j) + lj
                        END IF
                     ELSE
                        !
                        ! LJ potential
                        IF (r2 < rc2) THEN
                           epsilon = param%eps(id_c,id_a)
                           sigma = param%sig(id_c,id_a)
                           ecut = param%ecut(id_c,id_a)            
                           r2i= 1./r2
                           r6i = sigma**6*r2i**3
                           lj = 48.*r2i*r6i*(r6i - .5)
                           DO l = 1,3
                              q(k)%ff(l)  =  q(k)%ff(l) + lj*xr(l)*epsilon
                              q(kk)%ff(l) = q(kk)%ff(l) - lj*xr(l)*epsilon
                           END DO
                           Eint = 4.*r6i*(r6i - 1.)*epsilon - ecut
                           sys%Epot = sys%Epot + Eint
                           IF ((id_c/=2).AND.(id_a/=2)) sys%Ecc=sys%Ecc + Eint
                        ENDIF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
END DO
!
! <<<<<<<<< Lennar-Jones 12-6 potential and force in the same molecule >>>>>>>
DO i = 1, NS%Nchain
   DO j = 1, Ngroup(i)%Ng - 4
      k = Ngroup(i)%Ncomp(j)
      id_c = q(k)%id
      DO jj = j + 4, Ngroup(i)%Ng
         kk = Ngroup(i)%Ncomp(jj)
         id_a = q(kk)%id
         r2 = 0.0
         rc2 = param%rc2(id_c,id_a)
         DO l=1,3
            x = q(k)%xx(l) - q(kk)%xx(l)
            xr(l) = x - sys%box(l)*ANINT(x/sys%box(l))
            r2 = r2 + xr(l)**2
         END DO
         IF (r2 < rc2) THEN
            epsilon = param%eps(id_c,id_a)
            sigma = param%sig(id_c,id_a)
            ecut = param%ecut(id_c,id_a)
            r2i= 1./r2
            r6i = sigma**6*r2i**3
            lj = 48.*r2i*r6i*(r6i - .5)
            DO l = 1,3
               q(k)%ff(l)  =  q(k)%ff(l) + lj*xr(l)*epsilon
               q(kk)%ff(l) = q(kk)%ff(l) - lj*xr(l)*epsilon
            END DO
            Eint = 4.*r6i*(r6i - 1.)*epsilon - ecut
            sys%Epot = sys%Epot + Eint
            IF ( (id_c /=2) .AND. (id_a /=2) ) sys%Ecc = sys%Ecc + Eint
         ENDIF
      END DO
   END DO
END DO
!
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Torsion force >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!30 CALL Torsion2(NS, q, Ngroup, param, sys)
30 CALL Torsion4(NS, q, Ngroup, param, sys)
!
!
! T.6 Normalizing density
rho = rho*0.80229/8.6396
!
! T.7 Effective thiol potential and forces
DO i=1,NS%Nchain
   k = Ngroup(i)%Ncomp(1)
   x = q(k)%xx(1)
   y = q(k)%xx(2)
   z = q(k)%xx(3)
   CALL Thiol_new(rho(i), x, y, De, Dex, Dey, Re, Rex, Rey, beta)
   xp = exp(-beta*(z-Re))
   q(k)%ff(1) = q(k)%ff(1) + Dex*xp*(xp-2.) + De*2.*beta*xp*Rex*(xp-1.)
   q(k)%ff(2) = q(k)%ff(2) + Dey*xp*(xp-2.) + De*2.*beta*xp*Rey*(xp-1.)
   q(k)%ff(3) = q(k)%ff(3) - De*2.*beta*xp*(xp-1.)
   sys%Epot = sys%Epot - De*xp*(xp -2.)
END DO
40 RETURN
END SUBROUTINE Force
!
! ################## Surface potential for CHx group  ######################
!
SUBROUTINE Zforce(NS, q, Ngroup, param, sys)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(CH)::Ngroup(NS%Nchain)
TYPE(PM)::param
TYPE(ST)::sys
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, Ng
!
DO i=1,NS%Nchain
   Ng = Ngroup(i)%Ng
   DO j=3, Ng-1
      k = Ngroup(i)%Ncomp(j)
      q(k)%ff(3) = q(k)%ff(3) + 12.*2.80*8.617343E2/(q(k)%xx(3) - 0.860)**13- &
           3.*1.71*8.617343E-1/(q(k)%xx(3) - 0.860)**4
      sys%Epot = sys%Epot + 2.80*8.617343E2/(q(k)%xx(3) - 0.860)**12 - &
           1.71*8.617343E-1/(q(k)%xx(3) - 0.860)**3
   END DO
   k = Ngroup(i)%Ncomp(Ng)
   q(k)%ff(3) = q(k)%ff(3) + 12.*3.41*8.617343E2/(q(k)%xx(3) - 0.860)**13 - &
        3.*2.08*8.617343E-1/(q(k)%xx(3) - 0.860)**4
   sys%Epot = sys%Epot + 3.41*8.617343E2/(q(k)%xx(3) - 0.860)**12 - &
        2.08*8.617343E-1/(q(k)%xx(3) - 0.860)**3
END DO
!      
RETURN
END SUBROUTINE Zforce
!
!
! ################## TORSION ANGLE ESTIMATION AND FORCE ######################
!
! To simulate gauche defect mechanism, torsional stiffness is implemented.
! In common, k(t-t0)^2/2 is implemented. But this program implemented
! another way to simulate torsion of alkanechains. 
!
! As shown on Figure 2.8, pp.69, Organic Chemistry by Brown and Foote, 
! torsional energy shows simple-symmetric wavy shape. It is minimized
! when the torsion angle is 180 degrees and maximized for 0/360 degrees. 
! And there exist meta-stable points and they work as gauche defects. 
! The curve can be reproduced using combination of trigonometric function, 
! such as cos(x)+cos(3x), and fitting parameter for each constant can be found
! using trial and error. Current energy comes from .03Cos(x)+.07 cos(3*x)
!
! Algorithm is straightforward. For the given chain, go along from bottom atom.
! Then find 3 more atoms - to reproduce torsion, at least 4 atoms are needed -
! and we can build 2 surface using 3lower/3upper atoms. Using cross product
! of each bond vector, surface normal vectors will be available. They provide
! angle between surface using vector product and cosine relation.
! Using the angle, stiffness will be found - torsional energy, and torque.
! Torque will be calculated as force and this should be applied along the
! normal direction of surfaces.
!
!
!############ Algorithm is referred from Allen and Tildesley ##################
!
SUBROUTINE Torsion2(NS, q, Ngroup, param, sys)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, Ngrp, id1, id2, id3, id4
REAL*8:: rx(3,3), x, rx1, rx2, rx3, dot,  &
     C11, C12, C13, C22, C23, C33, D12, D13, D23
!
!
sys%Etor = 0.0
!
!
DO i=1, NS%Nchain
   Ngrp = Ngroup(i)%Ng
   IF (Ngrp > 3) THEN
      id1 = Ngroup(i)%Ncomp(1)
      id2 = Ngroup(i)%Ncomp(2)
      id3 = Ngroup(i)%Ncomp(3)
      DO k=1,3
         x = q(id2)%xx(k) - q(id1)%xx(k)
         rx(2,k) = x - sys%box(k)*ANINT(x/sys%box(k))
         x = q(id3)%xx(k) - q(id2)%xx(k)
         rx(3,k) = x - sys%box(k)*ANINT(x/sys%box(k))
      END DO
      DO j=2, Ngrp-2
         id1 = Ngroup(i)%Ncomp(j-1)
         id2 = Ngroup(i)%Ncomp(j)
         id3 = Ngroup(i)%Ncomp(j+1)
         id4 = Ngroup(i)%Ncomp(j+2)
         DO k=1,3
            rx(1,k) = rx(2,k)
            rx(2,k) = rx(3,k)
            x = q(id4)%xx(k) - q(id3)%xx(k)
            rx(3,k) = x - sys%box(k)*ANINT(x/sys%box(k))
         END DO
         !
         ! COS(Torsion angle)
         dot = - (rx(1,1)*rx(3,1)+rx(1,2)*rx(3,2)+rx(1,3)*rx(3,3)) / &
              SQRT(rx(1,1)**2+rx(1,2)**2+rx(1,3)**2)/ &
              SQRT(rx(3,1)**2+rx(3,2)**2+rx(3,3)**2)
         !
         ! Torsional energy
         sys%Etor = sys%Etor + ((param%kt(1) - 3.*param%kt(2)) &
              + 4.*param%kt(2)*dot*dot)*dot + param%kt(3)
         !
         ! Length of torque arm
         C11 = rx(1,1)**2 + rx(1,2)**2 + rx(1,3)**2
         C22 = rx(2,1)**2 + rx(2,2)**2 + rx(2,3)**2
         C33 = rx(3,1)**2 + rx(3,2)**2 + rx(3,3)**2
         C12 = rx(1,1)*rx(2,1) + rx(1,2)*rx(2,2) + rx(1,3)*rx(2,3)
         C23 = rx(2,1)*rx(3,1) + rx(2,2)*rx(3,2) + rx(2,3)*rx(3,3)
         C13 = rx(1,1)*rx(3,1) + rx(1,2)*rx(3,2) + rx(1,3)*rx(3,3)
         D12 = C11*C22 - C12*C12
         D13 = C11*C33 - C13*C13
         D23 = C22*C33 - C23*C23
         !
         DO k=1,3
           q(id1)%ff(k) = q(id1)%ff(k) + &
                (- param%kt(1) + 3.*param%kt(2) - 12.*param%kt(2)*dot**2)* &
                (-C23*rx(2,k) + C22*rx(3,k) - (C23*C12 - C13*C22)* &
                (-C22*rx(1,k) + C12*rx(2,k))/D12)/SQRT(D23*D12)
           q(id2)%ff(k) = q(id2)%ff(k) + &
                (- param%kt(1) + 3.*param%kt(2) - 12.*param%kt(2)*dot**2)* &
                (-C12*rx(3,k) + C23*rx(2,k) - C23*rx(1,k) - C22*rx(3,k) + &
                2.*C13*rx(2,k) - (C23*C12 - C13*C22)*(C22*rx(1,k) - &
                C11*rx(2,k) - C12*rx(2,k) + C12*rx(1,k))/D12 - (C23*C12 - &
                C13*C22)*(-C33*rx(2,k) + C23*rx(3,k))/D23)/SQRT(D23*D12)
           q(id3)%ff(k) = q(id3)%ff(k) + &
                (- param%kt(1) + 3.*param%kt(2) - 12.*param%kt(2)*dot**2)* &
                (C12*rx(3,k) - C12*rx(2,k) + C23*rx(1,k) + C22*rx(1,k) - &
                2.*C13*rx(2,k) - (C23*C12 - C13*C22)*(C11*rx(2,k) - &
                C12*rx(1,k))/D12 - (C23*C12 - C13*C22)*(C33*rx(2,k) - &
                C22*rx(3,k) - C23*rx(3,k) + C23*rx(2,k))/D23)/SQRT(D23*D12)
           q(id4)%ff(k) = q(id4)%ff(k) + &
                (- param%kt(1) + 3.*param%kt(2) - 12.*param%kt(2)*dot**2)* &
                (C12*rx(2,k) - C22*rx(1,k) - (C23*C12-C13*C22)*(C22*rx(3,k) -&
                C23*rx(2,k))/D23)/SQRT(D23*D12)
         END DO
     END DO
  END IF
END DO
RETURN
END SUBROUTINE Torsion2

!
! Shift Hautmann's potential with 180 angles
SUBROUTINE Torsion4(NS, q, Ngroup, param, sys)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: Ngroup(NS%Nchain)
TYPE(PM):: param
TYPE(ST):: sys
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, Ngrp, id1, id2, id3, id4
REAL*8:: rx(3,3), x, rx1, rx2, rx3, dot,  &
     C11, C12, C13, C22, C23, C33, D12, D13, D23, t(6)
!
!
sys%Etor = 0.0
!
! Torsional parameter from Hautmann and Klein, '89 JCP
t = RESHAPE((/ 1.116, 1.462, -1.578, -0.368, 3.156, -3.788 /), (/6/))
!
! From K -> eV
t = t/11.6045
!
!
DO i=1, NS%Nchain
   Ngrp = Ngroup(i)%Ng
   IF (Ngrp > 3) THEN
      id1 = Ngroup(i)%Ncomp(1)
      id2 = Ngroup(i)%Ncomp(2)
      id3 = Ngroup(i)%Ncomp(3)
      DO k=1,3
         x = q(id2)%xx(k) - q(id1)%xx(k)
         rx(2,k) = x - sys%box(k)*ANINT(x/sys%box(k))
         x = q(id3)%xx(k) - q(id2)%xx(k)
         rx(3,k) = x - sys%box(k)*ANINT(x/sys%box(k))
      END DO
      DO j=2, Ngrp-2
         id1 = Ngroup(i)%Ncomp(j-1)
         id2 = Ngroup(i)%Ncomp(j)
         id3 = Ngroup(i)%Ncomp(j+1)
         id4 = Ngroup(i)%Ncomp(j+2)
         DO k=1,3
            rx(1,k) = rx(2,k)
            rx(2,k) = rx(3,k)
            x = q(id4)%xx(k) - q(id3)%xx(k)
            rx(3,k) = x - sys%box(k)*ANINT(x/sys%box(k))
         END DO
         !
         ! COS(Torsion angle)
         ! Potential is minimum at pi(180 degrees)
         dot = - (rx(1,1)*rx(3,1)+rx(1,2)*rx(3,2)+rx(1,3)*rx(3,3)) / &
              SQRT(rx(1,1)**2+rx(1,2)**2+rx(1,3)**2)/ &
              SQRT(rx(3,1)**2+rx(3,2)**2+rx(3,3)**2)
         !
         ! Torsional energy
         sys%Etor = sys%Etor + t(1) - t(2)*dot + t(3)*dot**2 - &
              t(4)*dot**3 + t(5)*dot**4 - t(6)*dot**5
         !
         ! Length of torque arm
         C11 = rx(1,1)**2 + rx(1,2)**2 + rx(1,3)**2
         C22 = rx(2,1)**2 + rx(2,2)**2 + rx(2,3)**2
         C33 = rx(3,1)**2 + rx(3,2)**2 + rx(3,3)**2
         C12 = rx(1,1)*rx(2,1) + rx(1,2)*rx(2,2) + rx(1,3)*rx(2,3)
         C23 = rx(2,1)*rx(3,1) + rx(2,2)*rx(3,2) + rx(2,3)*rx(3,3)
         C13 = rx(1,1)*rx(3,1) + rx(1,2)*rx(3,2) + rx(1,3)*rx(3,3)
         D12 = C11*C22 - C12*C12
         D13 = C11*C33 - C13*C13
         D23 = C22*C33 - C23*C23
         !
         DO k=1,3
           q(id1)%ff(k) = q(id1)%ff(k) -  &
                ( - t(2) + 2.*t(3)*dot - 3.*t(4)*dot**2 + 4.*t(5)*dot**3 - &
                5.*t(6)*dot**4)* &
                (-C23*rx(2,k) + C22*rx(3,k) - (C23*C12 - C13*C22)* &
                (-C22*rx(1,k) + C12*rx(2,k))/D12)/SQRT(D23*D12)
           q(id2)%ff(k) = q(id2)%ff(k) -  &
                ( - t(2) + 2.*t(3)*dot - 3.*t(4)*dot**2 + 4.*t(5)*dot**3 - &
                5.*t(6)*dot**4)* &
                (-C12*rx(3,k) + C23*rx(2,k) - C23*rx(1,k) - C22*rx(3,k) + &
                2.*C13*rx(2,k) - (C23*C12 - C13*C22)*(C22*rx(1,k) - &
                C11*rx(2,k) - C12*rx(2,k) + C12*rx(1,k))/D12 - (C23*C12 - &
                C13*C22)*(-C33*rx(2,k) + C23*rx(3,k))/D23)/SQRT(D23*D12)
           q(id3)%ff(k) = q(id3)%ff(k) -  &
                ( - t(2) + 2.*t(3)*dot - 3.*t(4)*dot**2 + 4.*t(5)*dot**3 - &
                5.*t(6)*dot**4)* &
                (C12*rx(3,k) - C12*rx(2,k) + C23*rx(1,k) + C22*rx(1,k) - &
                2.*C13*rx(2,k) - (C23*C12 - C13*C22)*(C11*rx(2,k) - &
                C12*rx(1,k))/D12 - (C23*C12 - C13*C22)*(C33*rx(2,k) - &
                C22*rx(3,k) - C23*rx(3,k) + C23*rx(2,k))/D23)/SQRT(D23*D12)
           q(id4)%ff(k) = q(id4)%ff(k)  -  &
                ( - t(2) + 2.*t(3)*dot - 3.*t(4)*dot**2 + 4.*t(5)*dot**3 - &
                5.*t(6)*dot**4)* &
                (C12*rx(2,k) - C22*rx(1,k) - (C23*C12-C13*C22)*(C22*rx(3,k) -&
                C23*rx(2,k))/D23)/SQRT(D23*D12)
         END DO
     END DO
  END IF
END DO
RETURN
END SUBROUTINE Torsion4
