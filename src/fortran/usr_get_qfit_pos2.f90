!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qfit_pos2
SUBROUTINE usr_get_qfit_pos2(coefcount, input_radius, sizemap, sizeann, xmap, ymap, rmap, tmap, xann, yann) BIND(C)
   USE qfit
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: coefcount
   REAL(KIND=8), INTENT(IN) :: sizemap, sizeann
   REAL(KIND=8), INTENT(IN) :: input_radius
   REAL(KIND=8), DIMENSION(INT(sizemap, KIND=4)), INTENT(OUT) :: xmap, ymap, rmap, tmap
   REAL(KIND=8), DIMENSION(INT(sizeann, KIND=4)), INTENT(OUT) :: xann, yann

   CALL qfit_init_from_coef(INT(coefcount, KIND=4) - 2)

   CALL scan_pos_car(xmap, ymap, xann, yann, input_radius)
   CALL scan_pos_polar(rmap, tmap, input_radius)
   ! rmap = xmap**2 + ymap**2
   ! tmap = ATAN2(ymap, xmap)
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
