! =====================================================================================
! LUNA version 1.5
!
! AUTHOR: David Kipping
!         Columbia University
!         Email: dkipping@astro.columbia.edu
!
! CITATION: If using this code, please cite Kipping, D., 2011, 
!           'LUNA: An algorithm for generating dynamic planet-moon transits', 
!           MNRAS, 416, 689
!
! DESCRIPTION: LUNA generates dynamically accurate lightcurves from a 
!              planet-moon pair
!
! CHANGES: * v1.0 is the first new version since v0.1. It is aimed at being a 
!            non-beta code. This entails removing a lot of inefficient code, 
!            comments, unnecessary checks. It also implements several more 
!            efficient versions of the LUNA eqns.
!          * v1.1 corrects some minor errors with gamma/beta
!          * v1.2 corrects psi1 = bb -> ab*DCOS(ib)
!          * v1.3 switches esinw/ecosw to e/w; removes Kstar & error call terms
!          * v1.4 includes some minor speed enhancements
!          * v1.5 cleaned up some redunant parts of the code
!
! =====================================================================================

MODULE lunamod
use mandelmod
implicit none

CONTAINS

! ==============================================================
! === SUBROUTINE: LUNA ===
!
SUBROUTINE luna(t,Pb,T_b,p_in,ab,eb,wb,bb,Ps,phi_s,s_in,as,es,ws,is,Os,Msp_in,&
           u1,u2,fplan,reverse,nz,show,animate,limb,ItotL)

implicit none

 INTEGER :: i, j, nz, error
 INTEGER :: show, animate, limb
 REAL(8) :: tstep
 REAL(8), INTENT(IN) :: Pb, T_b, p_in, ab, eb, wb, bb
 REAL(8), INTENT(IN) :: Ps, phi_s, s_in, as, es, ws, is, Os, Msp_in
 REAL(8), INTENT(IN) :: u1, u2, fplan
 REAL(8) :: p, s, fmoon, Msp
 REAL(8) :: T_s, rb_mid, rs_mid, ib
 REAL(8) :: cosib, sinib, cosis, sinis, coswbOs, sinwbOs
 REAL(8) :: coswb, sinwb, cosws, sinws
 REAL(8) :: pin1, pout1, sin1, sout1, psin1, psout1
 REAL(8) :: temp
 REAL(8), DIMENSION(nz) :: varrhob, varrhos
 REAL(8), DIMENSION(nz), INTENT(IN) :: t
 REAL(8), DIMENSION(nz) :: tb, ts, fb, fs
 REAL(8), DIMENSION(nz) :: cosfb, sinfb, cosfs, sinfs
 REAL(8) :: gamma1, gamma2, beta1, beta2, psi1, epsilon2
 REAL(8), DIMENSION(nz) :: gammy, belly
 REAL(8), DIMENSION(nz) :: S_P, S_S, S_PS
 INTEGER, DIMENSION(nz) :: prinny, subby
 REAL(8), DIMENSION(nz) :: x12, y12, x13, y13, x23, y23
 REAL(8), DIMENSION(nz) :: cost1, sint1, cost2, sint2
 REAL(8), DIMENSION(nz) :: c1, c2, c3
 REAL(8), DIMENSION(nz) :: A_p, A_s
 REAL(8), DIMENSION(nz) :: IpL, IsL              ! Limb darkened output
 REAL(8), DIMENSION(nz), INTENT(OUT) :: ItotL    ! Limb darkened output
 REAL(8), DIMENSION(nz) :: IpN, IsN, ItotN       ! Non-LDed output
 REAL(8), DIMENSION(nz) :: Xp, Yp, Xs, Ys
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0
 REAL(8), PARAMETER :: piI = 0.318309886D0
 REAL(8), PARAMETER :: rad = 0.017453292519943295D0
 REAL(8), PARAMETER :: third = 0.3333333333333333D0
 REAL(8), PARAMETER :: Grv = 6.67428D-11
 LOGICAL :: verb, reverse
 
 verb = .FALSE. ! log file
  IF( verb ) THEN
   open(2,FILE='luna.log',FORM='FORMATTED',ACCESS='APPEND')
   write(2,*) ' - - - - - - - - - - - - - - - - - - - - - - - - -'
  END IF
 error = 0
 p = p_in
 s = s_in
 Msp = Msp_in
 fmoon = halfpi - ws
 IF( fmoon .LT. 0.0D0 ) THEN
  fmoon = fmoon + twopi
 END IF
 coswb = DCOS(wb); sinwb = DSIN(wb)
 cosws = DCOS(ws); sinws = DSIN(ws)
 rb_mid = ab*(1.0D0 - eb**2)/(1.0D0 + eb*sinwb)
 rs_mid = as*(1.0D0 - es**2)/(1.0D0 + es*sinws)
 cosib = bb/rb_mid
 ib = DACOS(cosib)
 sinib = DSIN(ib)
 cosis = DCOS(is); sinis = DSIN(is)
 coswbOs = DCOS(wb+Os); sinwbOs = DSIN(wb+Os)
 pin1   = (1.0D0 - p)
 pout1  = (1.0D0 + p)
 sin1   = (1.0D0 - s)
 sout1  = (1.0D0 + s)
 psin1  = (p - s)
 psout1 = (p + s)
 T_s = 0.5D0*piI*phi_s*Ps
 DO i=1,nz
  tb(i) = t(i) - T_b
  ts(i) = tb(i) - T_s
 END DO
 call kepler(eb,wb,Pb,fplan,nz,tb,fb)
 call kepler(es,ws,Ps,fmoon,nz,ts,fs)
 temp = 1.0D0 - eb**2
 DO i=1,nz
  cosfb(i) = DCOS(fb(i))
  sinfb(i) = DSIN(fb(i))
  cosfs(i) = DCOS(fs(i))
  sinfs(i) = DSIN(fs(i))
  varrhob(i) = temp/(1.0D0 + eb*cosfb(i))
 END DO
 temp = 1.0D0 - es**2
 DO i=1,nz
  varrhos(i) = temp/(1.0D0 + es*cosfs(i))
 END DO
 beta1  = -as*sinis*sinwbOs
 beta2  = as*coswbOs
 gamma1 = as*( cosib*sinis*coswbOs - sinib*cosis )
 gamma2 = as*cosib*sinwbOs
 epsilon2 = ab
 psi1 = ab*cosib
 DO i=1,nz
  belly(i) = beta1*varrhos(i)*( cosws*sinfs(i) + cosfs(i)*sinws ) + &
             beta2*varrhos(i)*( cosfs(i)*cosws - sinfs(i)*sinws )
  gammy(i) = gamma1*varrhos(i)*( cosws*sinfs(i) + cosfs(i)*sinws ) + &
             gamma2*varrhos(i)*( cosfs(i)*cosws - sinfs(i)*sinws )
  S_P(i) = ( epsilon2*varrhob(i)*( cosfb(i)*coswb - sinfb(i)*sinwb ) - &
             Msp*belly(i) )**2 + &
           ( psi1*varrhob(i)*( coswb*sinfb(i) + cosfb(i)*sinwb ) - &
             Msp*gammy(i) )**2
  S_P(i) = DSQRT(S_P(i))
 END DO
 DO i=1,nz
  S_S(i) = ( epsilon2*varrhob(i)*( cosfb(i)*coswb - sinfb(i)*sinwb ) + &
             belly(i) )**2 &
           + ( psi1*varrhob(i)*( coswb*sinfb(i) + cosfb(i)*sinwb ) + &
               gammy(i) )**2
  S_S(i) = DSQRT(S_S(i))
 END DO
 DO i=1,nz
  S_PS(i) = (1.0D0+Msp)**2*(gammy(i)**2+belly(i)**2)
  S_PS(i) = DSQRT(S_PS(i))
 END DO
 DO i=1,nz
  prinny(i) = 0; subby(i) = 0; temp = 0.0D0
  x12(i) = 0.0D0; y12(i) = 0.0D0
  x13(i) = 0.0D0; y13(i) = 0.0D0
  x23(i) = 0.0D0; y23(i) = 0.0D0
  cost1(i) = 0.0D0; sint1(i) = 0.0D0
  cost2(i) = 0.0D0; sint2(i) = 0.0D0
  c1(i) = 0.0D0 ; c2(i) = 0.0D0; c3(i) = 0.0D0
  IF( S_P(i) .GE. pout1 ) THEN
   IF ( S_S(i) .GE. sout1 ) THEN
    prinny(i) = 1
   ELSE IF ( S_S(i) .LE. sin1 ) THEN
    prinny(i) = 7
   ELSE !IF(sin1.LT.S_S(i) .AND. S_S(i).LT.sout1) THEN !<-not needed since only 3 possibilities
    prinny(i) = 4
   END IF
  ELSE IF( S_P(i) .LE. pin1 ) THEN
   IF ( S_S(i) .GE. sout1 ) THEN
    prinny(i) = 19
   ELSE IF ( S_S(i) .LE. sin1 ) THEN
    IF( S_PS(i) .LE. psin1 ) THEN
     prinny(i) = 27
    ELSE IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 25
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 26
    END IF
   ELSE !IF( sin1.LT.S_S(i) .AND. S_S(i).LT.sout1) THEN !<-not needed since only 3 possibilities
    ! Sub-tree: "p_in-s_part" => cases 22,23,24 = case 22 or 23
    ! [case 24 is unphysical]
    IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 22
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 23
    END IF
   END IF
  ELSE !IF(pin1.LT.S_P(i) .AND. S_P(i).LT.pout1) THEN
   IF ( S_S(i) .GE. sout1 ) THEN
    ! Sub-tree: "p_part-s_out" => cases 10,11,12 = case 10
    prinny(i) = 10
   ELSE IF ( S_S(i) .LE. sin1 ) THEN
    IF( S_PS(i) .LE. psin1 ) THEN
     prinny(i) = 18
    ELSE IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 16
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 17
    END IF
   ELSE !IF( sin1.LT.S_S(i) .AND. S_S(i).LT.sout1) THEN !<-not needed since only 3 possibilities
    IF( S_PS(i) .LE. psin1 ) THEN
     prinny(i) = 15
    ELSE IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 13
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 14
     x12(i) = ( 1.0D0 - p**2 + S_P(i)**2 )/( 2.0D0*S_P(i) )
     y12(i) = ( DSQRT(2.0D0*S_P(i)**2*(1.0D0+p**2) &
              - (1.0D0-p**2)**2 - S_P(i)**4) )/( 2.0D0*S_P(i) )
     cost1(i) = ( S_P(i)**2 + S_S(i)**2 - S_PS(i)**2 ) /( 2.0D0*S_P(i)*S_S(i) )
     sint1(i) = DSQRT(1.0D0 - cost1(i)**2)
     IF( ( (x12(i)-S_S(i)*cost1(i))**2 & 
          + (y12(i)+S_S(i)*sint1(i))**2 ) .LT. s**2 ) THEN
      IF( S_P(i) .GT. 1.0D0 ) THEN
       subby(i) = 5  ! 14.3a
      ELSE !IF S_P(i) .LT. 1.0D0 ) THEN
       subby(i) = 6  ! 14.3b
      END IF
     ELSE
      y13(i) = -( DSQRT( 2.0D0*S_S(i)**2*(1.0D0+s**2) &
               - (1.0D0-s**2)**2 - S_S(i)**4 ) )/(2.0D0*S_S(i))
      temp = (1.0D0 - s**2 + S_S(i)**2)/(2.0D0*S_S(i))
      x13(i) = temp*cost1(i) - y13(i)*sint1(i)
      y13(i) = temp*sint1(i) + y13(i)*cost1(i)
      cost2(i) = -(S_P(i)**2 + S_PS(i)**2 - S_S(i)**2)/(2.0D0*S_P(i)*S_PS(i))
      sint2(i) = DSQRT(1.0D0 - cost2(i)**2)
      y23(i) = ( DSQRT(2.0D0*S_PS(i)**2*(p**2+s**2) - (p**2-s**2)**2 & 
               - S_PS(i)**4) )/(2.0D0*S_PS(i))
      temp = (p*p - s*s + S_PS(i)**2)/(2.0D0*S_PS(i))
      x23(i) = temp*cost2(i) - y23(i)*sint2(i) + S_P(i)
      y23(i) = temp*sint2(i) + y23(i)*cost2(i)
      IF( ( (x12(i)-S_S(i)*cost1(i))**2 &
           + (y12(i)-S_S(i)*sint1(i))**2 ) .LT. s*s ) THEN
       IF( (S_S(i)*sint1(i)) .GT. ( y13(i) &
          + ( (y23(i)-y13(i))/(x23(i)-x13(i)) )*( S_S(i)*cost1(i) - x13(i) ) ) ) THEN
        subby(i) = 1
       ELSE
        subby(i) = 2
       END IF
      ELSE
       IF( ((x13(i)-S_P(i))**2 + (y13(i))**2) .LT. p**2 ) THEN
        IF( (S_S(i)-s) .LT. (S_P(i)-p) ) THEN
         subby(i) = 3
        ELSE
         subby(i) = 4
        END IF
       ELSE
        IF( (x23(i)**2 + y23(i)**2) .LT. 1.0D0 ) THEN
         subby(i) = 8
        ELSE
         subby(i) = 7
        END IF
       END IF
      END IF
     END IF ! Condition B
    END IF
   END IF
  END IF
 END DO
 
 DO i=1,nz
  IF( prinny(i) .LE. 9 ) THEN
   A_p(i) = 0.0D0
  ELSE IF( prinny(i) .GE. 19 ) THEN
   A_p(i) = pi*p**2
  ELSE
   call alpha(1.0D0,p,S_P(i),A_p(i))
  END IF
 END DO
 DO i=1,nz
 IF( prinny(i).EQ.1 .OR. prinny(i).EQ.10 .OR. &
 prinny(i).EQ.15 .OR. prinny(i).EQ.18 .OR. &
 prinny(i).EQ.19 .OR. prinny(i).EQ.27 ) THEN
  A_s(i) = 0.0D0
 ELSE IF( prinny(i).EQ.7 .OR. prinny(i).EQ.16 .OR. prinny(i).EQ.25 ) THEN
  A_s(i) = pi*s**2
 ELSE IF( prinny(i).EQ.4 .OR. prinny(i).EQ.13 &
 .OR. prinny(i).EQ.22) THEN
  call alpha(1.0D0,s,S_S(i),A_s(i))
 ELSE IF( prinny(i).EQ.17 .OR. prinny(i).EQ.26 ) THEN
  call alpha(p,s,S_PS(i),A_s(i))
  A_s(i) = pi*s**2 - A_s(i)
 ELSE IF( prinny(i).EQ.23 ) THEN
  call alpha(1.0D0,s,S_S(i),temp)
  call alpha(p,s,S_PS(i),A_s(i))
  A_s(i) = temp - A_s(i)
 ELSE IF( prinny(i).EQ.14 ) THEN
  IF( subby(i).EQ.4 ) THEN
   A_s(i) = 0.0D0
  ELSE IF( subby(i).EQ.3 ) THEN 
   call alpha(p,s,S_PS(i),A_s(i))
   A_s(i) = pi*s**2 - A_s(i)
  ELSE IF( subby(i).EQ.7 ) THEN
   call alpha(1.0D0,s,S_S(i),A_s(i))
  ELSE IF( subby(i).EQ.5 ) THEN
   call alpha(1.0D0,s,S_S(i),temp)
   call alpha(1.0D0,p,S_P(i),A_s(i))
   A_s(i) = temp - A_s(i)
  ELSE IF( subby(i).EQ.6 ) THEN
   call alpha(p,s,S_PS(i),temp)
   call alpha(1.0D0,p,S_P(i),A_s(i))
   A_s(i) = A_s(i) + temp
   call alpha(1.0D0,s,S_S(i),temp)
   A_s(i) = pi*p**2 + temp - A_s(i)
  ELSE IF( subby(i).EQ.8 ) THEN
   call alpha(1.0D0,s,S_S(i),temp)
   call alpha(p,s,S_PS(i),A_s(i))
   A_s(i) = temp - A_s(i)
  ELSE IF( subby(i).EQ.1 .OR. subby(i).EQ.2 ) THEN
   c3(i) = DSQRT( (x13(i)-x23(i))**2 + (y13(i)-y23(i))**2 )
   c1(i) = DSQRT( (x12(i)-x13(i))**2 + (y12(i)-y13(i))**2 )
   c2(i) = DSQRT( (x12(i)-x23(i))**2 + (y12(i)-y23(i))**2 )
   A_s(i) = 0.25D0*DSQRT( (c1(i)+c2(i)+c3(i))*(c2(i)+c3(i)-c1(i))*&
           (c1(i)+c3(i)-c2(i))*(c1(i)+c2(i)-c3(i)) ) + &
           DASIN(0.5D0*c1(i)) + p**2*DASIN(0.5D0*c2(i)/p) + &
           s**2*DASIN(0.5D0*c3(i)/s) - &
           0.25D0*c1(i)*DSQRT(4.0D0-c1(i)**2) - &
           0.25D0*c2(i)*DSQRT(4.0D0*p**2-c2(i)**2)
   IF( subby(i).EQ.1 ) THEN
    A_s(i) = A_s(i) - 0.25D0*c3(i)*DSQRT(4.0D0*s**2-c3(i)**2)
   ELSE
    A_s(i) = A_s(i) + 0.25D0*c3(i)*DSQRT(4.0D0*s**2-c3(i)**2)
   END IF
   call alpha(1.0D0,s,S_S(i),temp)
   A_s(i) = temp - A_s(i)
  ELSE
   !write(2,*) 'ERROR: no sub case identified for point ',i
  END IF
 ! --------------------------------------------
 ELSE
  !write(2,*) 'ERROR: no principal case identified for point ',i
 END IF
 IF( A_s(i) .LT. 0.0D0 ) THEN
  A_s(i) = 0.0D0 ! Override implemented
 END IF
 END DO
 DO i=1,nz
  IsN(i) = 1.0D0 - (A_s(i)*piI)
  ItotN(i) = 1.0D0 - ( (A_p(i)+A_s(i))*piI )
 END DO
 IF( limb .EQ. 1 ) THEN
  call occultquad(S_P,u1,u2,p,IpL,IpN,nz) ! Mandel-Agol call
  call ldsmall(S_S,A_s,0.0D0,u1+2.0D0*u2,0.0D0,-u2,s,IsL,nz)
 ELSE
  DO i=1,nz
   IpN(i) = 1.0D0 - A_p(i)*piI
   IpL(i) = IpN(i)
   IsL(i) = 1.0D0 - A_s(i)*piI
  END DO
 END IF
 IF( reverse ) THEN
   DO i=1,nz
    ItotL(i) = IpL(i) + 1.0D0 - IsL(i)
   END DO
 ELSE
   DO i=1,nz
    ItotL(i) = IpL(i) - 1.0D0 + IsL(i) 
   END DO
 END IF
 IF( verb ) THEN
  close(2)
 END IF

 END SUBROUTINE luna

! ==============================================================
! === SUBROUTINE: KEPLER ===
!
 SUBROUTINE kepler(e,wrad,Pdays,f_ref,n,t,f_)

 implicit none
 INTEGER :: i, j, n!, tt
 REAL(8), INTENT(IN) :: e, wrad, Pdays, f_ref
 REAL(8) :: E_ref, M_ref, ek, Pk, toler
 INTEGER :: ok
 REAL(8), DIMENSION(n), INTENT(IN) :: t
 REAL(8), DIMENSION(n) :: M_, E_, E_p
 REAL(8), DIMENSION(n), INTENT(OUT) :: f_
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0
 REAL(8), PARAMETER :: twopi = 6.2831853071795864769D0
     
 ek = DSQRT((1.0D0+e)/(1.0D0-e))
 Pk = twopi/Pdays
 IF( e .LT. 0.9D-4 ) THEN
   E_ref = f_ref
   M_ref = E_ref
   DO i=1,n
     M_(i) = Pk*t(i) + M_ref
     E_(i) = M_(i)
     f_(i) = E_(i)
   END DO
 ELSE
   E_ref = 2.0D0*DATAN((1.0D0/ek)*DTAN(f_ref*0.5D0))
   IF(E_ref .LT. -halfpi) THEN
     E_ref = E_ref + twopi
   END IF
   M_ref = E_ref - e*DSIN(E_ref)
   toler = 1.1574074074074074D-4*Pk ! Set to 10 second tolerance, divide by 10 for 1 second
   DO i=1,n
     M_(i) = Pk*t(i) + M_ref
   END DO
   E_(:) = M_(:)
   DO i=1,n
     ok = 0
     DO WHILE ( ok .EQ. 0 )
       E_p(i) = E_(i) - ((E_(i) - e*DSIN(E_(i)) &
                - M_(i))/(1.0D0-e*DCOS(E_(i))))
       IF( DABS(E_p(i)-E_(i)) .LT. toler ) THEN
         ok = 1
       END IF
       E_(i) = E_p(i)
     END DO
     f_(i) = 2.0D0*DATAN(ek*DTAN(E_(i)*0.5D0))
   END DO
 END IF

 END SUBROUTINE kepler
! =======================================================

! ==============================================================
! === SUBROUTINE: ALPHA ===
!
 SUBROUTINE alpha(Rb,Rs,S,alph)

 implicit none
 REAL(8), INTENT(IN) :: Rb, Rs, S
 REAL(8) :: kap0, kap1, kap2
 REAL(8), INTENT(OUT) :: alph
     
 kap0 = DACOS( (S**2 + Rs**2 - Rb**2)/(2.0D0*S*Rs) )
 kap1 = DACOS( (S**2 + Rb**2 - Rs**2)/(2.0D0*S*Rb) )
 kap2 = DSQRT(0.25D0*(4.0D0*S**2*Rb**2 - (Rb**2+S**2-Rs**2)**2))
 
 alph = Rs**2*kap0 + Rb**2*kap1 - kap2

 END SUBROUTINE alpha
! =======================================================

! ==============================================================
! === SUBROUTINE: LDSMALL ===
!
 SUBROUTINE ldsmall(S,Ar_occ,c1,c2,c3,c4,r,I0,n)

 implicit none
 REAL(8), DIMENSION(n), INTENT(IN) :: S, Ar_occ
 REAL(8), INTENT(IN) :: c1, c2, c3, c4, r
 INTEGER, INTENT(IN) :: n
 INTEGER :: i
 REAL(8) :: Ftot
 REAL(8), DIMENSION(n) :: am, bm, amR, bmR
 REAL(8), DIMENSION(n) :: Ar_ann, Fl_ann
 REAL(8), DIMENSION(n), INTENT(OUT) :: I0
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: third = 0.33333333333333333333D0
 REAL(8), PARAMETER :: seventh = 0.1428571428571428D0
 
 Ftot = 1.0D0 - 0.2D0*c1 - third*c2 - 3.0D0*seventh*c3 - 0.5D0*c4

 DO i=1,n
 am(i) = (S(i)-r)**2
 bm(i) = (S(i)+r)**2
 amR(i) = DSQRT(DSQRT(1.0D0-am(i)))
 bmR(i) = DSQRT(DSQRT(1.0D0-bm(i)))
  IF( S(i).GT.(1.0D0+r) ) THEN
   Ar_ann(i) = 0.0D0
   Fl_ann(i) = 0.0D0
   I0(i) = 1.0D0
  ELSE IF( S(i).GT.r .AND. S(i).LT.(1.0D0-r) ) THEN 
   Ar_ann(i) = pi*(bm(i) - am(i))
   Fl_ann(i) = (am(i)-bm(i))*(c1+c2+c3+c4-1.0D0) + &
               0.8D0*c1*amr(i)**5 + 2.0D0*third*c2*amr(i)**6 + &
               4.0D0*seventh*c3*amr(i)**7 + 0.5D0*c4*amr(i)**8 - &
               0.8D0*c1*bmr(i)**5 - 2.0D0*third*c2*bmr(i)**6 - &
               4.0D0*seventh*c3*bmr(i)**7 - 0.5D0*c4*bmr(i)**8
    I0(i) = 1.0D0 - (Ar_occ(i)/Ar_ann(i))*(Fl_ann(i)/Ftot)
  ELSE IF( S(i).LT.r ) THEN 
   Ar_ann(i) = pi*bm(i)
   Fl_ann(i) = -bm(i)*(c1+c2+c3+c4-1.0D0) + &
               0.8D0*c1 + 2.0D0*third*c2 + 4.0D0*seventh*c3 + 0.5D0*c4 - &
               0.8D0*c1*bmr(i)**5 - 2.0D0*third*c2*bmr(i)**6 - &
               4.0D0*seventh*c3*bmr(i)**7 - 0.5D0*c4*bmr(i)**8
   I0(i) = 1.0D0 - (Ar_occ(i)/Ar_ann(i))*(Fl_ann(i)/Ftot)
  ELSE !IF( S(i).LT.(1.0D0+r) .AND. S(i).LT.(1.0D0-r) ) THEN
   Ar_ann(i) = pi*(1.0D0-am(i))
   Fl_ann(i) = (am(i)-1.0D0)*(c1+c2+c3+c4-1.0D0) + &
               0.8D0*c1*amr(i)**5 + 2.0D0*third*c2*amr(i)**6 + &
               4.0D0*seventh*c3*amr(i)**7 + 0.5D0*c4*amr(i)**8 
   I0(i) = 1.0D0 - (Ar_occ(i)/Ar_ann(i))*(Fl_ann(i)/Ftot)
  END IF
 END DO

 END SUBROUTINE ldsmall
! =======================================================

END MODULE lunamod