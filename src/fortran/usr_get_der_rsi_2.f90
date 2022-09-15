!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_der_rsi_2
SUBROUTINE usr_get_der_rsi_2(pos_r, &
                             pos_t, &
                             n_in, &
                             n_out, &
                             maps, &
                             ray_tra_in, &
                             ray_tra_out, &
                             snrm_calib, &
                             snrm, &
                             der_r, &
                             der_t, &
                             mapcount) BIND(C)
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: mapcount
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 6), INTENT(IN) :: ray_tra_in, ray_tra_out
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 3), INTENT(IN) :: snrm_calib
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 3), INTENT(OUT) :: snrm
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)), INTENT(OUT) :: der_r, der_t
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)), INTENT(IN) :: pos_r, pos_t
   REAL(KIND=8), INTENT(IN) :: n_in, n_out, maps
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)) :: maps_scaled_x, maps_scaled_y
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)) :: tx1, tx2, ty1, ty2, tr1, tt1
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4)) :: der_x_m, der_y_m, der_x, der_y

   REAL(kind=8) :: efl2, alpha, n_in2, n_out2
   REAL(kind=8), DIMENSION(INT(mapcount, kind=4)) :: efl_par, rt, rr, dtx, dty, dtr, dtt
   REAL(kind=8), DIMENSION(INT(mapcount, kind=4), 3) :: obj, out
   REAL(kind=8), DIMENSION(INT(mapcount, kind=4), 2) :: snrm_calib_temp
   LOGICAL, DIMENSION(INT(mapcount, kind=4)) :: filter

   IF (n_in .EQ. n_out) THEN
      n_in2 = n_in
      n_out2 = -n_out
   ELSE
      n_in2 = n_in
      n_out2 = n_out
   END IF
   alpha = SIGN(1.0, n_in2*n_out2)

   rr = SQRT(ray_tra_out(:, 1)**2 + ray_tra_out(:, 2)**2)
   rt = SQRT(1 - ray_tra_in(:, 6)**2)
   efl_par = rr/rt
   filter = (rt <= MINVAL(rt) + 1E-10)
   efl2 = SUM(efl_par*filter)/SUM(filter*1)

   obj = ray_tra_in(:, 4:)
   out(:, :2) = ray_tra_out(:, :2)
   out(:, 3) = efl2
   ! calcul du cosinus directeur (n_out2*dir/norme)
   ! write (*,*) ray_tra_out(:5, 1)
   ! write (*,*)
   ! write (*,*) out(:5, 1)
   ! write (*,*)
   out = n_out2*out/SPREAD(SQRT(SUM(out**2, dim=2)), 2, 3)
   ! write (*,*) ray_tra_in(:5, 4)
   ! write (*,*)
   ! write (*,*) out(:5, 1)
   ! out = out
   ! out =    ! out = out/SPREAD(SQRT(SUM(out**2, dim=2)), 2, 3)
   ! out = n_out2*out

   snrm = obj - ((2*n_out2*obj/n_in2) - out)
   ! snrm = - ((2*n_out2*obj/n_in2) - out)

   ! write (*,*)
   ! write (*,*) snrm(:5, 1)
   ! write (*,*)
   ! write (*,*) "------------------------------------------"


   snrm = snrm/SPREAD(SQRT(SUM(snrm**2, dim=2)), 2, 3)

   IF (ALL(snrm_calib(:, :2) .NE. 0)) THEN
      ! snrm(:, :2) = snrm(:, :2) - snrm_calib(:, :2)
      snrm_calib_temp = snrm_calib(:, 3:1:-2)
      snrm_calib_temp(:, 1) = -snrm_calib_temp(:, 1)
      snrm(:, 1) = SUM(snrm_calib_temp*snrm(:, 1::2), dim=2)
      ! der_x = snrm(:,1)
      snrm_calib_temp = snrm_calib(:, 3:2:-1)
      snrm_calib_temp(:, 1) = -snrm_calib_temp(:, 1)
      snrm(:, 2) = SUM(snrm_calib_temp*snrm(:, 2:), dim=2)
      ! der_y = snrm(:,2)
      snrm(:, 3) = SQRT(1 - SUM(snrm(:, :2)**2, dim=2))
      snrm = snrm/SPREAD(SQRT(SUM(snrm**2, dim=2)), 2, 3)
   ! ELSE
   !    der_x = snrm(:, 1)/snrm(:, 3)
   !    der_y = snrm(:, 2)/snrm(:, 3)
   END IF

   ! write (*,*) "salut"

   der_x = snrm(:, 1)/snrm(:, 3)
   der_y = snrm(:, 2)/snrm(:, 3)

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
