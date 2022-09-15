!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_der_rsi
SUBROUTINE usr_get_der_rsi(pos_r, &
                           pos_t, &
                           n_in, &
                           n_out, &
                           maps, &
                           ray_tra_in, &
                           ray_tra_out, &
                           der_r, &
                           der_t, &
                           mapcount) BIND(C)
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: mapcount
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 6), INTENT(IN) :: ray_tra_in, ray_tra_out
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)), INTENT(OUT) :: der_r, der_t
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)), INTENT(IN) :: pos_r, pos_t
   REAL(KIND=8), INTENT(IN) :: n_in, n_out, maps
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)) :: tx1, tx2, ty1, ty2
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)) :: der_x, der_y

   tx1 = ATAN((ray_tra_in(:, 4)/ray_tra_in(:, 6)))
   ty1 = ATAN((ray_tra_in(:, 5)/ray_tra_in(:, 6)))
   tx2 = ATAN((ray_tra_out(:, 4)/ray_tra_out(:, 6)))
   ty2 = ATAN((ray_tra_out(:, 5)/ray_tra_out(:, 6)))


   der_x = (n_in*SIN(tx1) - n_out*SIN(tx2))
   der_y = (n_in*SIN(ty1) - n_out*SIN(ty2))
   IF (n_in /= n_out) THEN
      der_x = der_x/(n_out*COS(tx2) - n_in*COS(tx1))
      der_y = der_y/(n_out*COS(ty2) - n_in*COS(ty1))
   END IF
   der_r = maps*(der_x*COS(pos_t) + der_y*SIN(pos_t))
   der_t = maps*(pos_r*der_y*COS(pos_t) &
                 - pos_r*der_x*SIN(pos_t))
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
