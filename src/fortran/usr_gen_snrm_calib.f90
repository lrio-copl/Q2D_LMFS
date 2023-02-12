!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_gen_snrm_calib
SUBROUTINE usr_gen_snrm_calib(snrm, snrm_calib_x, snrm_calib_y, snrm_calib_z, mapcount) BIND(c)
   USE math_util
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: mapcount
   REAL(KIND=8), INTENT(IN), DIMENSION(INT(mapcount, KIND=4), 3) :: snrm
   REAL(KIND=8), INTENT(out), DIMENSION(INT(mapcount, KIND=4), 3) :: snrm_calib_x, snrm_calib_y, snrm_calib_z
   INTEGER :: i
   REAL(KIND=8), DIMENSION(3, 3) :: temp_rot

   snrm_calib_z = snrm

   DO i = 1, INT(mapcount, KIND=4)
      temp_rot = get_rot(snrm(i, :), REAL([0, 0, -1], kind=8))
      snrm_calib_x(i, :) = MATMUL(temp_rot, REAL([1, 0, 0], kind=8))
      snrm_calib_y(i, :) = MATMUL(temp_rot, REAL([0, 1, 0], kind=8))
   END DO

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
