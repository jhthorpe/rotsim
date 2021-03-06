MODULE term
  USE ints

CONTAINS

!------------------------------------------------------------
! term_harm
!       - calculates contribution of term1 (below) to Beff
!         for a given bra,ket, and rotational state
!
!       term1 = 1/2 <i|1/2 Σ_{r,s} μ_{α,β}^{r,s}*q_r*q_s |i>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, rotational axis α
! ket           : 1D int, ket QN's
! mu2           : 4D real*8, second order mu matrix

REAL(KIND=8) FUNCTION term_harm(nvib,a,ket,mu2)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  INTEGER, DIMENSION(0:), INTENT(IN) :: ket
  INTEGER, INTENT(IN) :: nvib,a
  INTEGER, DIMENSION(0:nvib-1) :: bra
  REAL(KIND=8) :: temp
  INTEGER :: j
  temp = 0.0D0
  bra = ket
  DO j=0,nvib-1
    temp = temp + mu2(j,j,a,a)*ints_qq(bra(j),ket(j))
  END DO
  term_harm = 0.25D0*temp

END FUNCTION term_harm

!------------------------------------------------------------
! term_cori
!       - calculates contributions of term2 (below) to Beff
!         of a given bra,ket, and rotational axis
!
!     term2 = 1/I_α^2 * Σ_i!=k <i|π_α|k><k|π_a|i>/(E_i - E_k)
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, rotational axis
! ket           : 1D int, bra quantum numbers
! Be            : 1D real*8, Be's
! phi2          : 1D real*8, quadratic force constants 
! zeta          : 3D real*8, Coriolis zetas (vib,vib,rot)

REAL(KIND=8) FUNCTION term_cori(nvib,a,ket,Be,phi2,zeta)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: zeta
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Be,phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: ket
  INTEGER, INTENT(IN) :: nvib,a
  INTEGER, DIMENSION(0:nvib-1) :: bra
  REAL(KIND=8) :: temp,tol,Ei,Ek,p1,p2
  INTEGER :: i,j
  tol = 1.0D-15
  temp = 0.0D0 

  Ei = 0.0D0
  DO i=0,nvib-1
    Ei = Ei + phi2(i)*(1.0D0*ket(i)+0.5D0) 
  END DO

  !Each vibrational qn can differ by only +/- 1, and two 
  ! different indicies must differ
  !(B/C no diagonal coriolis zetas to give 0,+/-2)
  DO i=0,nvib-2
    DO j=i+1,nvib-1
      IF (ABS(zeta(i,j,a)) .LT. tol) CYCLE
      bra = ket      
 
      !--
      bra(i) = ket(i) - 1
      bra(j) = ket(j) - 1
      IF (bra(i) .GE. 0 .AND. bra(j) .GE. 0) THEN
        Ek = Ei - phi2(i) - phi2(j)
        temp = temp + term_cori_aux(i,j,a,bra,ket,phi2,zeta)/(Ei - Ek)
      END IF

      !+-
      bra(i) = ket(i) + 1
      bra(j) = ket(j) - 1
      IF (bra(j) .GE. 0) THEN
        Ek = Ei + phi2(i) - phi2(j)
        temp = temp + term_cori_aux(i,j,a,bra,ket,phi2,zeta)/(Ei - Ek)
      END IF

      !-+
      bra(i) = ket(i) - 1
      bra(j) = ket(j) + 1
      IF (bra(i) .GE. 0) THEN
        Ek = Ei - phi2(i) + phi2(j)
        temp = temp + term_cori_aux(i,j,a,bra,ket,phi2,zeta)/(Ei - Ek)
      END IF

      !--
      bra(i) = ket(i) + 1
      bra(j) = ket(j) + 1
      Ek = Ei + phi2(i) + phi2(j)
      temp = temp + term_cori_aux(i,j,a,bra,ket,phi2,zeta)/(Ei - Ek)
    END DO
  END DO

  term_cori = 4.0D0*Be(a)*Be(a)*temp

END FUNCTION term_cori

!------------------------------------------------------------
! term_cori_aux
!       - auxilary function that evalutes 
!
!       <ket| π |bra><bra| π |ket>
!
!         for a given bra and ket that differ in the 
!         i'th and j'th indices
!------------------------------------------------------------
! i             : int, first index
! j             : int, second index
! a             : int, rotational axis
! bra           : 1D int, bra quantum numbers
! ket           : 1D int, ket quantum numbers
! phi2          : 1D real*8, quadratic force constants
! zeta          : 3D real*8, coriolis zetas

REAL(KIND=8) FUNCTION term_cori_aux(i,j,a,bra,ket,phi2,zeta)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: zeta
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: i,j,a
  REAL(KIND=8) :: val
  val = 0.0D0
  !ij ij
  val = val + zeta(i,j,a)*zeta(i,j,a)*phi2(j)/phi2(i)*&
              ints_q(ket(i),bra(i))*ints_p(ket(j),bra(j))*&
              ints_q(bra(i),ket(i))*ints_p(bra(j),ket(j))
 ! WRITE(*,*) "ij ij", val
  !ij ji
  val = val + zeta(i,j,a)*zeta(j,i,a)*&
              ints_q(ket(i),bra(i))*ints_p(ket(j),bra(j))*&
              ints_q(bra(j),ket(j))*ints_p(bra(i),ket(i))
 ! WRITE(*,*) "ij ji", val
  !ji ij
  val = val + zeta(j,i,a)*zeta(i,j,a)*&
              ints_q(ket(j),bra(j))*ints_p(ket(i),bra(i))*&
              ints_q(bra(i),ket(i))*ints_p(bra(j),ket(j))
 ! WRITE(*,*) "ji ij", val
  !ji ji
  val = val + zeta(j,i,a)*zeta(j,i,a)*phi2(i)/phi2(j)*&
              ints_q(ket(j),bra(j))*ints_p(ket(i),bra(i))*&
              ints_q(bra(j),ket(j))*ints_p(bra(i),ket(i))
 ! WRITE(*,*) "ji ji", val
  term_cori_aux = -1.0D0*val !- accounts for the two p integrals
END FUNCTION term_cori_aux

!------------------------------------------------------------


!------------------------------------------------------------
! term_cori_old
!       - calculates contributions of term2 (below) to Beff
!         of a given bra,ket, and rotational axis
!
!     term2 = 1/I_α^2 * Σ_i!=k <i|π_α|k><k|π_a|i>/(E_i - E_k)
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, rotational axis
! ket           : 1D int, bra quantum numbers
! Be            : 1D real*8, Be's
! phi2          : 1D real*8, quadratic force constants 
! zeta          : 3D real*8, Coriolis zetas (vib,vib,rot)

REAL(KIND=8) FUNCTION term_cori_old(nvib,a,ket,Be,phi2,zeta)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: zeta
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Be,phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: ket
  INTEGER, INTENT(IN) :: nvib,a
  INTEGER, DIMENSION(0:nvib-1) :: bra
  REAL(KIND=8) :: temp,tol,Ei,Ek,p1,p2
  INTEGER :: i,j
  tol = 1.0D-15
  temp = 0.0D0 

  Ei = 0.0D0
  DO i=0,nvib-1
    Ei = Ei + phi2(i)*(1.0D0*ket(i)+0.5D0) 
  END DO

  !WRITE(*,*) "|",ket(0:nvib-1),">",Ei

  !Each vibrational qn can differ by only +/- 1, and two 
  ! different indicies must differ
  !(B/C no diagonal coriolis zetas to give 0,+/-2)
  DO j=0,nvib-2
    DO i=j+1,nvib-1
      IF (ABS(zeta(i,j,a)) .LT. tol) CYCLE
      bra = ket      

      !-1,-1
      bra(i) = ket(i) - 1
      bra(j) = ket(j) - 1
      IF (bra(i) .GE. 0 .AND. bra(j) .GE. 0) THEN
        Ek = Ei - phi2(i) - phi2(j)
        temp = temp - &
               ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))*&
               ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))/&
               (Ei - Ek)
!        p1 = ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
!                    zeta(j,i,a),phi2(i),phi2(j))
!        p2 =  ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
!                     zeta(j,i,a),phi2(i),phi2(j))
!        WRITE(*,*) "<",bra(0:nvib-1),"|",Ek,p1,p2
      END IF

      !-1,+1
      bra(i) = ket(i) - 1
      bra(j) = ket(j) + 1
      IF (bra(i) .GE. 0) THEN
        Ek = Ei - phi2(i) + phi2(j)
        temp = temp - &
               ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))*&
               ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))/&
               (Ei - Ek)
!        p1 = ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
!                     zeta(j,i,a),phi2(i),phi2(j))
!        p2 =  ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
!                     zeta(j,i,a),phi2(i),phi2(j))
!        WRITE(*,*) "<",bra(0:nvib-1),"|",Ek,p1,p2
      END IF

      !+1,-1
      bra(i) = ket(i) + 1
      bra(j) = ket(j) - 1
      IF (bra(j) .GE. 0) THEN
        Ek = Ei + phi2(i) - phi2(j)
        temp = temp - &
               ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))*&
               ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))/&
               (Ei - Ek)
!        p1 = ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
!                     zeta(j,i,a),phi2(i),phi2(j))
!        p2 =  ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
!                     zeta(j,i,a),phi2(i),phi2(j))
!        WRITE(*,*) "<",bra(0:nvib-1),"|",Ek,p1,p2
      END IF

      !+1,+1
      bra(i) = ket(i) + 1
      bra(j) = ket(j) + 1
      Ek = Ei + phi2(i) + phi2(j)
      temp = temp - &
             ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                     zeta(j,i,a),phi2(i),phi2(j))*&
             ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                     zeta(j,i,a),phi2(i),phi2(j))/&
             (Ei - Ek)
!      p1 = ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
!                   zeta(j,i,a),phi2(i),phi2(j))
!      p2 =  ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
!                   zeta(j,i,a),phi2(i),phi2(j))
!      WRITE(*,*) "<",bra(0:nvib-1),"|",Ek,p1,p2
       
    END DO 
  END DO

  term_cori_old = temp*4.0D0*Be(a)*Be(a)
  
END FUNCTION term_cori_old

!------------------------------------------------------------
! term_anh
!       - evalutes term3 (below) for a given ket and 
!         rotational axis
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, rotational axis α
! ket           : 1D int, ket quantum numbers
! phi2          : 1D real*8, quadratic force constants
! phi3          : 3D real*8, cubic force constants
! mu1           : 3D real*8, first order mu (vib,rot,rot)

REAL(KIND=8) FUNCTION term_anh(nvib,a,ket,phi2,phi3,mu1)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:) :: phi3,mu1
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: ket
  INTEGER, INTENT(IN) :: nvib,a
  INTEGER, DIMENSION(0:nvib-1) :: bra
  REAL(KIND=8) :: temp,tol,Ei,Ek
  INTEGER :: i
  temp = 0.0D0
  tol = 1.0D-15
  
  Ei = 0.0D0
  DO i=0,nvib-1
    Ei = Ei + phi2(i)*(1.0D0*ket(i)+0.5D0) 
  END DO

  !Each vibrational qn can differ by only +/- 1, and one
  ! qn must differ
  DO i=0,nvib-1
    IF (ABS(mu1(i,a,a)) .LT. tol) CYCLE
    bra = ket

    !-1
    bra(i) = ket(i) - 1
    IF (bra(i) .GE. 0) THEN
      Ek = Ei - phi2(i)
      temp = temp + mu1(i,a,a)*ints_q(ket(i),bra(i))*&
             ints_phi3(nvib,i,bra,ket,phi3)/&
             (Ei - Ek) 
    END IF

    !+1
    bra(i) = ket(i) + 1
    Ek = Ei + phi2(i)
    temp = temp + mu1(i,a,a)*ints_q(ket(i),bra(i))*&
           ints_phi3(nvib,i,bra,ket,phi3)/&
           (Ei - Ek) 

  END DO

  term_anh = temp

END FUNCTION term_anh

!------------------------------------------------------------

END MODULE term
