!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qfit_pos2_der
SUBROUTINE usr_get_qfit_pos2_der(coefcount, input_radius, sizemap, xmap, ymap, rmap, tmap) BIND(C)
   USE qfit
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: coefcount
   REAL(KIND=8), INTENT(IN) :: sizemap
   REAL(KIND=8), INTENT(IN) :: input_radius
   REAL(KIND=8), DIMENSION(INT(sizemap, KIND=4)), INTENT(OUT) :: xmap, ymap, rmap, tmap

   CALL qfit_init_from_coef(INT(coefcount) - 2)

   CALL scan_pos_car_der(xmap, ymap, input_radius)
   CALL scan_pos_polar_der(rmap, tmap, input_radius)
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
