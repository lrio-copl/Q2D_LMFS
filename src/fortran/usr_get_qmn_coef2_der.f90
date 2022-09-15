!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!
!   USING DERIVATIVE
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qmn_coef2_der
SUBROUTINE usr_get_qmn_coef2_der(coefcount, coef, radius, sizemap, mapr, mapt) BIND(C)
   USE qfit
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: coefcount, sizemap
   REAL(KIND=8), INTENT(IN) :: radius
   REAL(KIND=8), DIMENSION(INT(sizemap, KIND=4)), INTENT(IN) :: mapr, mapt
   REAL(KIND=8), DIMENSION(INT(coefcount, KIND=4)), INTENT(OUT) :: coef
   ! REAL(KIND=8), DIMENSION(INT(coefcount, KIND=4)), INTENT(IN) :: calib
   REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: amn, bmn

   CALL qfit_init_from_coef(INT(coefcount, KIND=4) - 2)

   ALLOCATE (amn(m_max + 1, n_max + 1))
   ALLOCATE (bmn(m_max + 1, n_max + 1))

   amn = 0
   bmn = 0

   CALL q_fit_der(mapr, mapt, radius, amn, bmn)
   CALL arr2data(coef(3:), amn, bmn)

   coef(1) = radius
   coef(2) = 1
   coef(3:) = coef(3:)!*calib(3:)
END SUBROUTINE

!DEC$ ATTRIBUTES DLLEXPORT :: cvisthreadsafe
LOGICAL FUNCTION CVISTHREADSAFE() BIND(c)
   CVISTHREADSAFE = .FALSE.
END

!DEC$ ATTRIBUTES DLLEXPORT :: cvumrtype
SUBROUTINE cvumrtype(umrType)
   INTEGER umrType
   CHARACTER*4 ThisType
   ThisType = 'USR '
   READ (ThisType, '(1a4)') umrType
END
