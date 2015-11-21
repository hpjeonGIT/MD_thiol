!
! ############# Effective thiol potential and force routine ################
!
!
SUBROUTINE Thiol(rho, xx, yy, De, Dex, Dey, Re, Rex, Rey, beta)
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
REAL*8:: rho, xx, yy, De, Dex, Dey, Re, Rex, Rey, beta
!
! Internal variables
INTEGER:: i,j
REAL*8:: k0, kvt(3,2), r02, ecub(4,6), rcub(4,5), ccub(4,2), a(6), b(5), c(2),&
     aa, sum(5,3), sqrt3, pi
!
! Initialization
sqrt3 = SQRT(3.)
aa = 2.884
pi = 3.14159265359
k0 = 4.*pi/sqrt3/aa
kvt(1,1) = k0*sqrt3/2.
kvt(1,2) = k0/2.
kvt(2,1) = -kvt(1,1)
kvt(2,2) = k0/2.
kvt(3,1) = 0.
kvt(3,2) = -k0
r02 = .55**2
IF (rho < 1.0) THEN
   !
   ! 3layer of Au(111)
   !ecub = RESHAPE((/ -7.2782, 10.9173, -0.8222, -4.4286, &
   !     -0.5476, 0.8215, -0.0455, -0.2303, &
   !     -0.1444, 0.2166, -0.0158, -0.0693, &
   !     1.9074, -2.8612, 0.1794, 0.8521, &
   !     0.0377, -0.0565, -0.0106, 0.0415, &
   !     0.5287, -0.7930, 0.0643, 0.3000/), (/4,6/))
   !rcub = RESHAPE((/ 0.2376, -0.3564, -0.0421, 2.3931, &
   !     0.0032, -0.0048, 0.0247, 0.0670, &
   !     -0.0009, 0.0014, -0.0115, 0.0104, &
   !     -0.081, 0.0122, -0.0025, -0.0046, &
   !     -0.0079, 0.0118, -0.0042, 0.0028/), (/4,5/))
   !
   ! 4layer of Au(111)
   ecub = RESHAPE((/ 6.6853, -10.0279, 1.1325, -2.1521, &
         0.4260, -0.6390,  0.0908, -0.0605, &
        -0.0375,  0.0562, -0.0129, -0.0293, &
        -1.1700,  1.7549, -0.2214,  0.2727, &
        -0.0155,  0.0232, -0.0205,  0.0314, &
         0.0240, -0.0360, -0.0080,  0.1700 /), (/4,6/))
   rcub = RESHAPE((/ 0.0890, -0.1336, 0.0126, 2.2579, &
        -0.0522,  0.0784,  0.0199, 0.0487, &
         0.0302, -0.0452, -0.0068, 0.0177, &
         0.0054, -0.0081, -0.0105, 0.0084, &
        -0.0022,  0.0032, -0.0105, 0.0123 /), (/4,5/))
   !
   ! commmon
   ccub = RESHAPE((/ 0.5153, -0.7729, 3.9444, -6.8674, &
        -1.3989, 2.0984, -5.1467, 14.1945 /), (/4,2/))
ELSE
   !
   ! 3layer of Au(111)
   !ecub = RESHAPE((/ 10.9324, -43.7145, 53.8096, -22.6392, &
   !     0.8226, -3.2892, 4.0652, -1.6005, &
   !     0.2169, -0.8674, 1.0683, -0.4307, &
   !     -2.8651, 11.4564, -14.1382, 5.6246, &
   !     -0.0566, 0.2263, -0.2935, 0.1358, &
   !     -0.7941, 3.1754, -3.9041, 1.6228/), (/4,6/))
   !rcub = RESHAPE((/ -0.3569, 1.4270, -1.8255, 2.9876, &
   !     -0.0048, 0.0191, 0.0009, 0.0749, &
   !     0.0014, -0.0056, -0.0045, 0.0081, &
   !     0.0122, -0.0489, 0.0587, -0.0250, &
   !     0.0118, -0.0472, 0.0547, -0.0168/), (/4,5/))
   !
   ! 4layer of Au(111)
   ecub = RESHAPE((/ -10.0418, 40.1532, -49.0487, 14.5749, &
        -0.6398, 2.5585, -3.1066, 1.0053, &
         0.0563, -0.2250,  0.2683, -0.1231, &
         1.7574, -7.0270,  8.5605, -2.6547, &
         0.0232, -0.0929,  0.0957, -0.0073, &
        -0.0361,  0.1442, -0.1882,  0.2301 /), (/4,6/))
   rcub = RESHAPE((/ -0.1337, 0.5348, -0.6557, 2.4807, &
         0.0785, -0.3138,  0.4121, -0.0820, &
        -0.0453,  0.1811, -0.2332,  0.0931, &
        -0.0081,  0.0323, -0.0509,  0.0219, &
         0.0032, -0.0129,  0.0057,  0.0069 /), (/4,5/))
   !
   ! common
   ccub = RESHAPE((/ -0.7740, 3.0949, 0.0766, -5.5781, &
        2.1013, -8.4023, 5.3540, 10.6943 /), (/4,2/))
END IF

DO i=1,5
   a(i) = rho**3*ecub(1,i) + rho**2*ecub(2,i) +rho*ecub(3,i) + ecub(4,i)
   b(i) = rho**3*rcub(1,i) + rho**2*rcub(2,i) +rho*rcub(3,i) + rcub(4,i)
END DO
a(6) = rho**3*ecub(1,6) + rho**2*ecub(2,6) +rho*ecub(3,6) + ecub(4,6)
c(1) = rho**3*ccub(1,1) + rho**2*ccub(2,1) +rho*ccub(3,1) + ccub(4,1)
c(2) = rho**3*ccub(1,2) + rho**2*ccub(2,2) +rho*ccub(3,2) + ccub(4,2) 
!
!
sum = 0.0
DO i=1,3
   sum(1,1) = sum(1,1) + cos(kvt(i,1)*xx+kvt(i,2)*yy)
   sum(2,1) = sum(2,1) + cos(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)
   sum(3,1) = sum(3,1) + sin(kvt(i,1)*xx+kvt(i,2)*yy)
   sum(4,1) = sum(4,1) + sin(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)
   sum(5,1) = sum(5,1) + (2.+cos(kvt(i,1)*xx+kvt(i,2)*yy))**.5
   sum(1,2) = sum(1,2) - sin(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,1)
   sum(1,3) = sum(1,3) - sin(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,2)
   sum(3,2) = sum(3,2) + cos(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,1)
   sum(3,3) = sum(3,3) + cos(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,2)
   sum(4,2) = sum(4,2) + 2.*cos(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)*kvt(i,1)
   sum(4,3) = sum(4,3) + 2.*cos(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)*kvt(i,2)
   sum(5,2) = sum(5,2) - sin(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,1)/ &
        (2.+cos(kvt(i,1)*xx+kvt(i,2)*yy))**.5/2.
   sum(5,3) = sum(5,3) - sin(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,2)/ &
        (2.+cos(kvt(i,1)*xx+kvt(i,2)*yy))**.5/2.      
END DO
De = a(1) + a(2)*sum(1,1) + a(3)*sum(3,1) + a(4)*sum(5,1) + a(5)*sum(4,1) + &
     a(6)/(1.+(sum(3,1)-sqrt3*1.5)**2/r02)
Dex = a(2)*sum(1,2) + a(3)*sum(3,2) + a(4)*sum(5,2) + a(5)*sum(4,2) - &
     a(6)*2.*(sum(3,1)-sqrt3*1.5)*sum(3,2)/r02/ &
     (1.+(sum(3,1)-sqrt3*1.5)**2/r02)**2
Dey = a(2)*sum(1,3) + a(3)*sum(3,3) + a(4)*sum(5,3) + a(5)*sum(4,3) - &
     a(6)*2.*(sum(3,1)-sqrt3*1.5)*sum(3,3)/r02/ &
     (1.+(sum(3,1)-sqrt3*1.5)**2/r02)**2
Re = b(1) + b(2)*sum(1,1) + b(3)*sum(2,1) + b(4)*sum(3,1) + b(5)*sum(4,1)
Rex = b(2)*sum(1,2) + b(3)*sum(2,2) + b(4)*sum(3,2) + b(5)*sum(4,2)
Rey = b(2)*sum(1,3) + b(3)*sum(2,3) + b(4)*sum(3,3) + b(5)*sum(4,3)
beta = c(2)*Re/aa/aa + c(1)/aa
RETURN
END SUBROUTINE Thiol
!
SUBROUTINE Thiol_new(rho, xx, yy, De, Dex, Dey, Re, Rex, Rey, beta)
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
REAL*8:: rho, xx, yy, De, Dex, Dey, Re, Rex, Rey, beta
!
! Internal variables
INTEGER:: i,j
REAL*8:: k0, kvt(3,2), r02, ecub(4,6), ccub(4), a(6), b(5), c,&
     aa, coe(5,3), sqrt3, pi
!
! Initialization
sqrt3 = SQRT(3.)
aa = 2.884
pi = 3.14159265359
k0 = 4.*pi/sqrt3/aa
kvt(1,1) = k0*sqrt3/2.
kvt(1,2) = k0/2.
kvt(2,1) = -kvt(1,1)
kvt(2,2) = k0/2.
kvt(3,1) = 0.
kvt(3,2) = -k0
r02 = .55**2
IF (rho < 1.0) THEN
   !
   ! 5layer of Au(111)
   ecub = RESHAPE((/1.76041E-02,-2.64061E-02,3.00502E-01,-1.62700E+00, &
        4.45752E-02,      -6.68628E-02,     -6.54561E-03,    6.89231E-02, &
        -1.85083E-02,      2.77625E-02,     -1.87089E-03,    -4.73690E-04, &
        1.96817E-02,      -2.95225E-02,     7.64664E-03,    -2.38711E-02, &
        -1.63136E-02,      2.44704E-02,     -1.01038E-02,    1.00148E-02, &
        2.32562E-02,      -3.48843E-02,     1.38175E-02,    1.15634E-02  &
        /), (/4,6/))
   ! commmon
   ccub = RESHAPE((/-9.6600E-3, 1.4490E-2, 3.5450E-2, 4.5360E-2/), (/4/))
ELSE
   ecub = RESHAPE((/-2.64425E-02, 1.05734E-01, 1.68362E-01, -1.58295E+00, &
        -6.69550E-02,      2.67728E-01,     -3.41136E-01,    1.80453E-01, &
        2.78008E-02,      -1.11165E-01,     1.37056E-01,    -4.67828E-02, &
        -2.95632E-02,      1.18212E-01,     -1.40088E-01,    2.53738E-02, &
        2.45042E-02,      -9.79829E-02,     1.12349E-01,    -3.08030E-02, &
        -3.49324E-02,      1.39682E-01,     -1.60748E-01,    6.97520E-02  &
        /), (/4,6/))
   !
   ! common
   ccub = RESHAPE((/1.4510E-2, -5.8020E-2, 1.0796E-1, 2.1190E-2/), (/4/))
END IF

DO i=1,5
   a(i) = rho**3*ecub(1,i) + rho**2*ecub(2,i) +rho*ecub(3,i) + ecub(4,i)
END DO
b = RESHAPE((/ 2.26704,    0.0977112,  0.000551831,  -0.00126157,    &
     0.00351557  /), (/5/))
a(6) = rho**3*ecub(1,6) + rho**2*ecub(2,6) +rho*ecub(3,6) + ecub(4,6)
c = rho**3*ccub(1) + rho**2*ccub(2) +rho*ccub(3) + ccub(4)

!
!
coe = 0.0
DO i=1,3
   ! basic
   coe(1,1) = coe(1,1) + cos(kvt(i,1)*xx+kvt(i,2)*yy)
   coe(2,1) = coe(2,1) + sin(kvt(i,1)*xx+kvt(i,2)*yy)
   coe(3,1) = coe(3,1) + cos(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)
   coe(4,1) = coe(4,1) + sin(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)
   coe(5,1) = coe(5,1) + cos(3.*kvt(i,1)*xx+3.*kvt(i,2)*yy)
   ! gradient
   coe(1,2) = coe(1,2) - sin(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,1)
   coe(1,3) = coe(1,3) - sin(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,2)
   coe(2,2) = coe(2,2) + cos(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,1)
   coe(2,3) = coe(2,3) + cos(kvt(i,1)*xx+kvt(i,2)*yy)*kvt(i,2)
   coe(3,2) = coe(3,2) - 2.*sin(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)*kvt(i,1)
   coe(3,3) = coe(3,3) - 2.*sin(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)*kvt(i,2)
   coe(4,2) = coe(4,2) + 2.*cos(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)*kvt(i,1)
   coe(4,3) = coe(4,3) + 2.*cos(2.*kvt(i,1)*xx+2.*kvt(i,2)*yy)*kvt(i,2)
   coe(5,2) = coe(5,2) - 3.*sin(3.*kvt(i,1)*xx+3.*kvt(i,2)*yy)*kvt(i,1)
   coe(5,3) = coe(5,3) - 3.*sin(3.*kvt(i,1)*xx+3.*kvt(i,2)*yy)*kvt(i,2)
END DO
De = a(1) + a(2)*coe(1,1) + a(3)*coe(2,1) + a(4)*coe(3,1) + a(5)*coe(4,1) + &
     a(6)*coe(5,1)
Dex = a(2)*coe(1,2) + a(3)*coe(2,2) + a(4)*coe(3,2) + a(5)*coe(4,2) + &
     a(6)*coe(5,2)
Dey = a(2)*coe(1,3) + a(3)*coe(2,3) + a(4)*coe(3,3) + a(5)*coe(4,3) + &
     a(6)*coe(5,3)
Re = b(1) + b(2)*coe(1,1) + b(3)*coe(2,1) + b(4)*coe(3,1) + b(5)*coe(4,1) 
Rex = b(2)*coe(1,2) + b(3)*coe(2,2) + b(4)*coe(3,2) + b(5)*coe(4,2) 
Rey = b(2)*coe(1,3) + b(3)*coe(2,3) + b(4)*coe(3,3) + b(5)*coe(4,3) 
beta = c*Re
RETURN
END SUBROUTINE Thiol_new
