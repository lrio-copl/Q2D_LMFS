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
                             snrm_calib_x, &
                             snrm_calib_y, &
                             snrm_calib_z, &
                             snrm, &
                             der_r, &
                             der_t, &
                             mapcount) BIND(C)
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: mapcount
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 6), INTENT(IN) :: ray_tra_in, ray_tra_out
   REAL(KIND=8), DIMENSION(INT(mapcount, KIND=4), 3), INTENT(IN) :: snrm_calib_x, snrm_calib_y, snrm_calib_z
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
   REAL(kind=8), DIMENSION(3) :: snrm_calib_temp
   REAL(kind=8), DIMENSION(INT(mapcount, kind=4)) :: pos_rr
   REAL(KIND=8) :: curr, prev
   INTEGER(kind=4) :: i
   LOGICAL, DIMENSION(INT(mapcount, kind=4)) :: filter

   IF (n_in .EQ. n_out) THEN
      n_in2 = n_in
      n_out2 = -n_out
   ELSE IF (n_in .EQ. -n_out) THEN
      n_in2 = n_in
      n_out2 = n_out
   ELSE
      n_in2 = n_in
      n_out2 = n_out
      IF (n_in2 < n_out2) THEN
         n_out2 = -n_out2
      END IF
   END IF
   ! n_in2 = 1.5
   ! n_in2 = -1
   ! n_out2 = 1
   alpha = SIGN(1.0, n_in2*n_out2)

   rr = SQRT(ray_tra_out(:, 1)**2 + ray_tra_out(:, 2)**2)
   rt = SIGN(1.0, n_in)*SQRT(1 - ray_tra_in(:, 6)**2)
   efl_par = rr/rt
   filter = (rt <= MINVAL(rt) + 1E-10)
   efl2 = SUM(efl_par*filter)/SUM(filter*1)

   obj = ray_tra_in(:, 4:)
   out(:, :2) = ray_tra_out(:, :2)
   out(:, 3) = efl2
   out = n_out2*out/SPREAD(SQRT(SUM(out**2, dim=2)), 2, 3)
   snrm = obj - ((2*n_out2*obj/n_in2) - out)

   snrm = snrm/SPREAD(SQRT(SUM(snrm**2, dim=2)), 2, 3)

   IF (ALL(SUM(snrm_calib_z(:, :)**2, dim=2) .NE. 0)) THEN
      DO i = 1, INT(mapcount, KIND=4)
         snrm_calib_temp = snrm(i, :)
         snrm(i, 1) = DOT_PRODUCT(snrm_calib_temp, snrm_calib_x(i, :))
         snrm(i, 2) = DOT_PRODUCT(snrm_calib_temp, snrm_calib_y(i, :))
         snrm(i, 3) = DOT_PRODUCT(snrm_calib_temp, snrm_calib_z(i, :))
      END DO
   END IF

   der_x = snrm(:, 1)/snrm(:, 3)
   der_y = snrm(:, 2)/snrm(:, 3)

   curr = pos_r(INT(mapcount, kind=4))
   prev = REAL(0, KIND=8)
   DO i = INT(mapcount, kind=4), 1, -1
      IF (curr .NE. pos_r(i)) THEN
         prev = curr
         curr = pos_r(i)
      END IF
      pos_rr(i) = curr - prev
   END DO

   der_r = maps*(der_x*COS(pos_t) + der_y*SIN(pos_t))
   der_r = der_r - SIN(rt)*SUM(pos_rr*der_r)/SUM(pos_rr*SIN(rt))
   ! dr = dr - test1[:, 0] * np.sum(drr * dr) / np.sum(drr * test1[:, 0])

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
