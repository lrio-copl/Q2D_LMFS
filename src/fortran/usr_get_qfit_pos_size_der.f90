!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qfit_pos_size_der
SUBROUTINE usr_get_qfit_pos_size_der(coefcount, sizemap) BIND(C)
   USE qfit
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: coefcount
   REAL(KIND=8), INTENT(OUT) :: sizemap
   CALL qmnp_init_from_coef(INT(coefcount, KIND=4) - 2)
   sizemap = REAL((j_max)*2*(k_max+1), KIND=8)
END SUBROUTINE

!DEC$ ATTRIBUTES DLLEXPORT :: cvisthreadsafe
LOGICAL FUNCTION CVISTHREADSAFE() BIND(c)
   CVISTHREADSAFE = .TRUE.
END

!DEC$ ATTRIBUTES DLLEXPORT :: cvumrtype
SUBROUTINE cvumrtype(umrType)
   INTEGER umrType
   CHARACTER*4 ThisType
   ThisType = 'USR '
   READ (ThisType, '(1a4)') umrType
END
