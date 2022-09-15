!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!
!   USING DERIVATIVE
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_remove_curv_der
SUBROUTINE usr_remove_curv_der(coefcount, coef, curv) BIND(C)
   USE qfit
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: coefcount
   REAL(KIND=8), INTENT(OUT) :: curv
   REAL(KIND=8), DIMENSION(INT(coefcount, KIND=4)), INTENT(INOUT) :: coef
   REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: amn, bmn
   ! REAL(KIND=8), DIMENSION(2) :: BACKUPCOEF

   CALL qfit_init_from_coef(INT(coefcount, KIND=4) - 2)

   ALLOCATE (amn(m_max + 1, n_max + 1))
   ALLOCATE (bmn(m_max + 1, n_max + 1))

   amn = 0
   bmn = 0
   CALL data2arr(COEF(3:), &
                 amn, bmn)
   CALL q_fit_remove_curv(amn, bmn, curv)
   curv = curv/(coef(1)**2)
   CALL arr2data(COEF(3:), amn, bmn)
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

