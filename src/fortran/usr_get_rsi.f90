!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_rsi
SUBROUTINE usr_get_rsi( n_in, &
                           n_out, &
                           exp,&
                           ray_tra_in, &
                           ray_tra_out, &
                           map_out, &
                           mapcount) BIND(C)
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: mapcount
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 6), INTENT(IN) :: ray_tra_in, ray_tra_out
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)), INTENT(OUT) :: map_out
   REAL(KIND=8), INTENT(IN) :: n_in, n_out, exp

   ! map_out = exp - (sqrt(ray_tra_out(:,1)**2 + ray_tra_out(:,2)**2)*ray_tra_out(:,6)/(sqrt(1-ray_tra_out(:,6)**2)))
   ! write (*,*) map_out
   write (*,*) ray_tra_out(:,3)
   map_out = 0
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
