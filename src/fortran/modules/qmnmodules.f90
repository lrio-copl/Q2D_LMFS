MODULE math_util
   IMPLICIT NONE
   REAL(KIND=8), PARAMETER :: qmnpi = 4.D0*DATAN(1.D0)
CONTAINS

   FUNCTION cross_product(a, b) RESULT(c)
      REAL(KIND=8), INTENT(IN), DIMENSION(3):: a, b
      REAL(KIND=8), DIMENSION(3) :: c

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END FUNCTION

   FUNCTION get_rot(a, b) RESULT(rot)
      REAL(KIND=8), INTENT(IN), DIMENSION(3) :: a, b
      REAL(KIND=8), DIMENSION(3) :: an, bn, cn
      REAL(KIND=8), DIMENSION(3, 3) :: rot, skew_symmetric, identity
      INTEGER(KIND=4) :: i, j

      FORALL (i=1:3, j=1:3) identity(i, j) = REAL((i/j)*(j/i), kind=8)
      an = a/NORM2(a)
      bn = b/NORM2(b)
      cn = CROSS_PRODUCT(an, bn)
      skew_symmetric(1,:) = [REAL(0, KIND=8), -cn(3), cn(2)]

      skew_symmetric(2,:) = [cn(3), REAL(0, KIND=8), -cn(1)]

      skew_symmetric(3,:) = [-cn(2), cn(1), REAL(0, KIND=8)]
      rot = identity + skew_symmetric + MATMUL(skew_symmetric, skew_symmetric)*(1/(1 - DOT_PRODUCT(an, bn)))
   END FUNCTION

   PURE FUNCTION kron_delta(i, j) BIND(c) RESULT(output)

      INTEGER(KIND=4), INTENT(IN), VALUE :: i, j
      INTEGER(KIND=4) :: output

      IF (i == j) THEN
         output = 1
      ELSE
         output = 0
      END IF
      ! output = MERGE(1,0,i==j)
   END FUNCTION

   PURE FUNCTION dct_iv(dctivdata) RESULT(xk)
      REAL(KIND=8), INTENT(IN), VALUE :: dctivdata(:)
      INTEGER(KIND=4) :: i
      REAL(KIND=8) :: nv(SIZE(dctivdata))
      REAL(KIND=8) :: xk(SIZE(dctivdata))

      nv = [(REAL(i + 0.5, kind=8), i=0, SIZE(dctivdata) - 1)]
      xk = 0
      DO i = 1, SIZE(dctivdata)
         xk(i) = SUM(dctivdata*COS((qmnpi*(i - 0.5)/SIZE(dctivdata))*nv))
      END DO
      xk = xk*SQRT(2.0D+00/REAL(SIZE(dctivdata), KIND=8))
   END FUNCTION
END MODULE math_util

MODULE qmn_array_index
   IMPLICIT NONE
   INTEGER(KIND=4) :: max_input
CONTAINS

   PURE FUNCTION num2qmnmax(input_number) RESULT(output) BIND(c)
      INTEGER(KIND=4), INTENT(IN) :: input_number
      ! REAL(KIND=8) :: current_order_float
      INTEGER(KIND=4) :: output
      ! output = CEILING(SQRT(REAL(FLOOR(REAL(MODULO(input_number,2)+(input_number))/2)))+1)
      ! current_order_float = (1.0D+00 + SQRT(REAL(1 + (8 * (input_number)),KIND=8))) / 4.0D+00
      output = FLOOR((1.0D+00 + SQRT(REAL(1 + (8*(input_number - 1)), KIND=8)))/4.0D+00) + 1
   END FUNCTION

   PURE SUBROUTINE num2qmnindex(input_number, table_name, table_index) BIND(c)
      INTEGER(KIND=4), INTENT(IN) :: input_number
      ! CHARACTER(LEN=1), INTENT(OUT) :: table_name
      INTEGER(KIND=4), INTENT(OUT) :: table_name
      INTEGER(KIND=4), DIMENSION(2), INTENT(OUT) :: table_index

      ! INTEGER(KIND=4) :: table_local_index, row_number, column_triangle_top
      ! INTEGER(KIND=4) :: column_triangle, column_triangle_index, row_index
      REAL(KIND=8) :: current_order_float
      INTEGER(KIND=4) :: current_order, current_index, current_index_max

      current_order_float = (1.0D+00 + SQRT(REAL(1 + (8*(input_number)), KIND=8)))/4.0D+00
      current_order = CEILING(current_order_float)
      current_index = (input_number - ((2*(current_order - 1)**2) - current_order + 1)) - 1
      current_index_max = 2*(current_order**2) - current_order - &
                          (2*((current_order - 1)**2) - current_order + 1)

      ! IF (MODULO(current_index, 2)) THEN
      !    table_name = "b"
      ! ELSE
      !    table_name = "a"
      ! END IF
      table_name = MODULO(current_index, 2)

      IF (current_index > 0) THEN
         table_index(1) = MIN(FLOOR(REAL(current_index - 1, KIND=8)/2) + 2, &
                              current_order)
         table_index(2) = MIN(FLOOR(REAL(current_index_max - current_index + 1, KIND=8)/2), &
                              current_order)
      ELSE
         table_index(1) = 1
         table_index(2) = current_order
      END IF
      ! table_index = table_index(ubound(table_index,dim=1)::-1)
   END SUBROUTINE

   SUBROUTINE data2arr_c(DATA, datasize, amnloc, bmnloc, arrlocsize1, arrlocsize2) BIND(c, name="data2arr")
      INTEGER(KIND=4) :: datasize, arrlocsize1, arrlocsize2
      REAL(KIND=8), DIMENSION(datasize), INTENT(IN) :: DATA
      REAL(KIND=8), DIMENSION(arrlocsize1, arrlocsize2), INTENT(OUT) :: amnloc, bmnloc
      CALL data2arr(DATA, amnloc, bmnloc)

   END SUBROUTINE

   SUBROUTINE data2arr(DATA, amnloc, bmnloc)
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: DATA
      REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: amnloc, bmnloc
      INTEGER(KIND=4) :: i, placename
      INTEGER(KIND=4), DIMENSION(2) :: placeind
      ! CHARACTER(LEN=1) :: placename

      DO i = 1, SIZE(DATA)
         CALL num2qmnindex(i, placename, placeind)
         IF (placename) THEN
            bmnloc(placeind(1), placeind(2)) = DATA(i)
         ELSE
            amnloc(placeind(1), placeind(2)) = DATA(i)
         END IF
      END DO
   END SUBROUTINE

   SUBROUTINE arr2data_c(DATA, datasize, amnloc, bmnloc, arrlocsize1, arrlocsize2) BIND(c, name="arr2data")
      INTEGER(KIND=4) :: datasize, arrlocsize1, arrlocsize2
      REAL(KIND=8), DIMENSION(datasize), INTENT(OUT) :: DATA
      REAL(KIND=8), DIMENSION(arrlocsize1, arrlocsize2), INTENT(IN) :: amnloc, bmnloc
      CALL arr2data(DATA, amnloc, bmnloc)

   END SUBROUTINE

   SUBROUTINE arr2data(DATA, amnloc, bmnloc)
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: DATA
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: amnloc, bmnloc
      INTEGER(KIND=4) :: i, PLACENAME
      INTEGER(KIND=4), DIMENSION(2) :: placeind, size_amnloc, size_bmnloc
      ! CHARACTER(LEN=1) :: placename

      size_amnloc = SHAPE(amnloc)
      size_bmnloc = SHAPE(bmnloc)
      DO i = 1, SIZE(DATA)
         CALL num2qmnindex(i, placename, placeind)
         IF ((.NOT. placename) .AND. &
             (placeind(1) <= size_amnloc(1)) .AND. &
             (placeind(2) <= size_amnloc(2))) THEN
            DATA(i) = amnloc(placeind(1), placeind(2))
         ELSE IF ((placename) .AND. &
                  (placeind(1) <= size_bmnloc(1)) .AND. &
                  (placeind(2) <= size_bmnloc(2))) THEN
            DATA(i) = bmnloc(placeind(1), placeind(2))
         ELSE
            DATA(i) = 0
         END IF
      END DO
   END SUBROUTINE

END MODULE qmn_array_index

MODULE asymjacobip
   ! Fortran rewrite of asmjacp from python package Scikit-qfit
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: abc_mn
   INTEGER(KIND=4) :: n_max, m_maxasymjacob
   INTEGER(KIND=4), DIMENSION(2) :: shapeAbc_mn
   PRIVATE :: m_maxasymjacob, shapeAbc_mn

CONTAINS
   SUBROUTINE asymjacobip_init(n_input) BIND(c)

      INTEGER(KIND=4), INTENT(IN), VALUE:: n_input
      REAL(KIND=8), ALLOCATABLE ::  tmpjvec(:)
      INTEGER(KIND=4) :: n0jvec

      n_max = n_input

      IF (n_max < 5) THEN
         WRITE (*, '(a)') ' '
         WRITE (*, '(a)') 'ASYMJACOBIP - Fatal error!'
         WRITE (*, '(a,g14.6)') '  Illegal input value of N_MAX = ', n_max
         WRITE (*, '(a)') '  But N_MAX must be equal or greater than 5.'
         STOP
      END IF
   END SUBROUTINE

   SUBROUTINE build_recursion(m_input) BIND(c)
                !!!!!!!!!!!!!!!!!!!!
      !
      ! Build the recursion coefficients and saves them as a sequence of tuples.
      ! These are the coefficients to build the m type polynominals up to order n.
      !
                !!!!!!!!!!!!!!!!!!!!

      INTEGER(KIND=4), INTENT(IN), VALUE :: m_input
      INTEGER(KIND=4) :: n_min, i
      REAL(KIND=8), ALLOCATABLE :: tmpabc_mn(:, :)
      INTEGER(KIND=4) ::  n0abc_mn(2)
      REAL(KIND=8) :: n, m_real

      m_maxasymjacob = m_input

      ! Make sure that only possible n values are computed
      n_min = 3
      IF (m_maxasymjacob > 1) THEN
         n_min = 1
      END IF

      ! Reallocate memory for abc_mn
      IF (.NOT. ALLOCATED(abc_mn)) THEN
         ALLOCATE (abc_mn(n_max - n_min, 4))
      ELSE
         n0abc_mn = SHAPE(abc_mn)
         IF (.NOT. (n0abc_mn(1) == (n_max - n_min))) THEN
            ALLOCATE (tmpabc_mn(n_max - n_min, 4))
            CALL MOVE_ALLOC(from=tmpabc_mn, to=abc_mn)
         END IF
      END IF
      shapeAbc_mn = SHAPE(abc_mn)

      ! abc_mn = 0
      ! ! Precompute and save polynomial factors
      ! m_real = REAL(m_maxasymjacob, KIND=8)
      ! DO CONCURRENT (i=1:n_max-n_min)
      !         n = REAL(i - 1 + n_min,KIND=8)
      !         abc_mn(i,1) = REAL( n ,KIND=8)
      !         abc_mn(i,2) = ((2*n-1)*&
      !                         (m_real+2*n-2)*&
      !                         (4*n*(m_real+n-2)+(m_real-3)*(2*m_real-1))) /&
      !                         ((4*n*n-1)*(m_real+n-2)*(m_real+2*n-3))
      !         abc_mn(i,3) = (-2*(2*n-1)*&
      !                         (m_real+2*n-1)*&
      !                         (m_real+2*n-2)*(m_real+2*n-3)) /&
      !                        ((4*n*n-1)*(m_real+n-2)*(m_real+2*n-3))
      !         abc_mn(i,4) = (n*(2.0*n-3)*(m_real+2*n-1)*(2*m_real+2*n-3)) /&
      !                        ((4*n*n-1)*(m_real+n-2)*(m_real+2*n-3))
      ! END DO
      CALL build_abc_mn(m_input, abc_mn)
   END SUBROUTINE

   SUBROUTINE build_abc_mn(m_local, abc_mn_loc)
      INTEGER(KIND=4), INTENT(IN) :: m_local
      REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: abc_mn_loc
      INTEGER(KIND=4), DIMENSION(2) :: s_abc_mn_loc
      INTEGER(KIND=4) :: n_min, i
      REAL(KIND=8) :: m_l_r, n_l_r
      s_abc_mn_loc = SHAPE(abc_mn_loc)
      abc_mn_loc = 0
      ! Precompute and save polynomial factors
      n_min = 3
      IF (m_local > 1) THEN
         n_min = 1
      END IF
      m_l_r = REAL(m_local, KIND=8)
      DO i = 1, s_abc_mn_loc(1)
         n_l_r = REAL(i - 1 + n_min, KIND=8)
         abc_mn_loc(i, 1) = REAL(n_l_r, KIND=8)
         abc_mn_loc(i, 2) = ((2*n_l_r - 1)* &
                             (m_l_r + 2*n_l_r - 2)* &
                             (4*n_l_r*(m_l_r + n_l_r - 2) + (m_l_r - 3)*(2*m_l_r - 1)))/ &
                            ((4*n_l_r*n_l_r - 1)*(m_l_r + n_l_r - 2)*(m_l_r + 2*n_l_r - 3))
         abc_mn_loc(i, 3) = (-2*(2*n_l_r - 1)* &
                             (m_l_r + 2*n_l_r - 1)* &
                             (m_l_r + 2*n_l_r - 2)*(m_l_r + 2*n_l_r - 3))/ &
                            ((4*n_l_r*n_l_r - 1)*(m_l_r + n_l_r - 2)*(m_l_r + 2*n_l_r - 3))
         abc_mn_loc(i, 4) = (n_l_r*(2.0*n_l_r - 3)*(m_l_r + 2*n_l_r - 1)*(2*m_l_r + 2*n_l_r - 3))/ &
                            ((4*n_l_r*n_l_r - 1)*(m_l_r + n_l_r - 2)*(m_l_r + 2*n_l_r - 3))
      END DO
   END SUBROUTINE

   SUBROUTINE jmat_x_c(n, xvp, jmat) BIND(C, NAME="jmat_x")
      INTEGER(KIND=4), INTENT(IN) :: n
      REAL(KIND=8), DIMENSION(n), INTENT(IN) :: xvp
      REAL(KIND=8), DIMENSION(n_max + 1, n), INTENT(OUT) :: jmat

      CALL jmat_x(xvp, jmat)
   END SUBROUTINE

   SUBROUTINE jmat_x_raw(abc_mn_loc, m_loc, xv, jmat)
                !!!!!!!!!!!!!!
      !
      ! Builds the asymmetric Jacobi P polynomial as defined in [2] A.1 for
      ! all a vector of x values.
      !
                !!!!!!!!!!!!!!

      REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xv
      REAL(KIND=8), INTENT(OUT), DIMENSION(n_max + 1, SIZE(xv)) :: jmat
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: abc_mn_loc
      INTEGER(KIND=4), INTENT(IN) :: m_loc
      INTEGER(KIND=4), DIMENSION(2) :: s_abc_mn_loc
      REAL(KIND=8) :: vn(SIZE(xv)), vn_(SIZE(xv)), vntmp(SIZE(xv))
      INTEGER(KIND=4) :: i

      s_abc_mn_loc = SHAPE(abc_mn_loc)
      jmat = 0

      IF (m_loc == 1) THEN
         jmat(1, :) = 0.5D+00
         jmat(2, :) = 1.0D+00 - xv/2.0D+00
         jmat(3, :) = (3.0D+00 + xv*(-12.0D+00 + 8.0D+00*xv))/6.0D+00
         jmat(4, :) = (5.0D+00 + xv*(-60.0D+00 + &
                                     (120.0D+00 - 64.0D+00*xv)*xv))/10.0D+00
         vn_ = jmat(3, :)
         vn = jmat(4, :)
      ELSE
         jmat(1, :) = 0.5D+00
         jmat(2, :) = (REAL(m_loc) - 0.5D+00) + (1.0D+00 - REAL(m_loc))*xv
         vn_ = jmat(1, :)
         vn = jmat(2, :)
      END IF

      DO i = 1, s_abc_mn_loc(1)
         vntmp = vn
         vn = ((abc_mn_loc(i, 2) + abc_mn_loc(i, 3)*xv)*vn - abc_mn_loc(i, 4)*vn_)
         vn_ = vntmp
         jmat(INT(abc_mn_loc(i, 1) + 2, KIND=4), :) = vn
      END DO
   END SUBROUTINE

   SUBROUTINE jmat_x(xv, jmat)
      REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xv
      REAL(KIND=8), INTENT(OUT), DIMENSION(n_max + 1, SIZE(xv)) :: jmat
      CALL jmat_x_raw(abc_mn, m_maxasymjacob, xv, jmat)

   END SUBROUTINE

   SUBROUTINE jmat_u_x_c(n, uvp, xvp, jmat) BIND(C, NAME="jmat_u_x")
      INTEGER(KIND=4), INTENT(IN) :: n
      REAL(KIND=8), DIMENSION(n), INTENT(IN) :: xvp, uvp
      REAL(KIND=8), DIMENSION(n_max + 1, n), INTENT(OUT) :: jmat

      CALL jmat_u_x(uvp, xvp, jmat)
   END SUBROUTINE

   SUBROUTINE jmat_u_x_raw(abc_mn_loc, m_loc, uv, xv, jmat)

                !!!!!!!!!!!!!!
      ! Builds the asymmetric Jacobi P polynomial as defined in [2]
      ! A.1 for all of the
      ! x values with the scaling factor of u**m which extends the
      ! range of usable (m,n)
      ! values before an overflow condition occurs.
                !!!!!!!!!!!!!!

      REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xv
      REAL(KIND=8), INTENT(IN), DIMENSION(SIZE(xv)) :: uv
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: abc_mn_loc
      INTEGER(KIND=4), INTENT(IN) :: m_loc
      REAL(KIND=8), INTENT(OUT), DIMENSION(n_max + 1, SIZE(xv)) :: jmat
      REAL(KIND=8), DIMENSION(SIZE(xv)) :: vn, vn_, vntmp, upm
      INTEGER(KIND=4), DIMENSION(2) :: s_abc_mn_loc
      INTEGER(KIND=4) :: i, n_

      s_abc_mn_loc = SHAPE(abc_mn_loc)
      jmat = 0

      IF (m_loc == 1) THEN
         jmat(1, :) = 0.5D+00
         jmat(2, :) = 1.0D+00 - xv/2.0D+00
         jmat(3, :) = (3.0D+00 + xv*(-12.0D+00 + 8.0D+00*xv))/6.0D+00
         jmat(4, :) = (5.0D+00 + xv*(-60.0D+00 + &
                                     (120.0D+00 - 64.0D+00*xv)*xv))/10.0D+00
         vn_ = jmat(3, :)
         vn = jmat(4, :)
      ELSE
         jmat(1, :) = 0.5D+00
         jmat(2, :) = uv*((REAL(m_loc) - 0.5D+00) + (1.0D+00 - REAL(m_loc))*xv)
         vn_ = jmat(1, :)
         vn = jmat(2, :)
      END IF

      DO i = 1, s_abc_mn_loc(1)
         n_ = INT(abc_mn_loc(i, 1) + 1, KIND=4)
         IF (n_ < m_loc) THEN
            vntmp = vn
            vn = uv*(((abc_mn_loc(i, 2) + abc_mn_loc(i, 3)*xv)*vn) - &
                     (uv*abc_mn_loc(i, 4)*vn_))
            vn_ = vntmp
         ELSE IF (n_ == m_loc) THEN
            vntmp = vn
            vn = (((abc_mn_loc(i, 2) + abc_mn_loc(i, 3)*xv)*vn) - &
                  (uv*abc_mn_loc(i, 4)*vn_))
            vn_ = vntmp
         ELSE
            vntmp = vn
            vn = (((abc_mn_loc(i, 2) + abc_mn_loc(i, 3)*xv)*vn) - &
                  (abc_mn_loc(i, 4)*vn_))
            vn_ = vntmp
         END IF
         jmat(INT(abc_mn_loc(i, 1) + 2, KIND=4), :) = vn
      END DO

      upm = uv**(MAX(0, m_loc - 1 - (n_max)))
      DO i = MIN(n_max, m_loc - 1), 0, -1
         jmat(i + 1, :) = jmat(i + 1, :)*upm
         upm = upm*uv
      END DO
   END SUBROUTINE

   SUBROUTINE jmat_u_x(uv, xv, jmat)
      REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xv
      REAL(KIND=8), INTENT(IN), DIMENSION(SIZE(xv)) :: uv
      REAL(KIND=8), INTENT(OUT), DIMENSION(n_max + 1, SIZE(xv)) :: jmat
      CALL jmat_u_x_raw(abc_mn, m_maxasymjacob, uv, xv, jmat)
   END SUBROUTINE

   SUBROUTINE jvec_x_raw(abc_mn_loc, m_loc, x, jvec)

                !!!!!!!!!!!!!!
      !
      ! Builds the asymmetric Jacobi P polynomial as defined in [2] A.1 for
      ! all a vector of x values.
      !
                !!!!!!!!!!!!!!

      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: abc_mn_loc
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(OUT) :: jvec
      INTEGER(KIND=4), INTENT(IN) :: m_loc
      INTEGER(KIND=4), DIMENSION(2) :: s_abc_mn_loc
      REAL(KIND=8) :: vn, vn_, vntmp
      INTEGER(KIND=4) :: i

      s_abc_mn_loc = SHAPE(abc_mn_loc)
      jvec = 0

      IF (m_loc == 1) THEN
         jvec(1) = 0.5
         jvec(2) = 1 - x/2
         jvec(3) = (3 + x*(-12 + 8*x))/6
         jvec(4) = (5 + x*(-60 + (120 - 64*x)*x))/10
         vn_ = jvec(3)
         vn = jvec(4)
      ELSE
         jvec(1) = 0.5
         jvec(2) = (REAL(m_loc) - 0.5) + (1 - REAL(m_loc))*x
         vn_ = jvec(1)
         vn = jvec(2)
      END IF

      DO i = 1, s_abc_mn_loc(1)
         vntmp = vn
         vn = ((abc_mn_loc(i, 2) + abc_mn_loc(i, 3)*x)*vn - abc_mn_loc(i, 4)*vn_)
         vn_ = vntmp
         jvec(INT(abc_mn_loc(i, 1) + 2, KIND=4)) = vn
      END DO
   END SUBROUTINE

   SUBROUTINE jvec_x(x, jvec) BIND(c)
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(OUT) :: jvec
      CALL jvec_x_raw(abc_mn, m_maxasymjacob, x, jvec)
   END SUBROUTINE

   SUBROUTINE jvec_x_b(m_loc, x, jvec)
      INTEGER(KIND=4), INTENT(IN) :: m_loc
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(OUT) :: jvec
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: abc_mn_loc
      INTEGER(KIND=4) :: n_min

      n_min = 3
      IF (m_loc > 1) THEN
         n_min = 1
      END IF

      ALLOCATE (abc_mn_loc(n_max - n_min, 4))
      CALL build_abc_mn(m_loc, abc_mn_loc)
      CALL jvec_x_raw(abc_mn_loc, m_loc, x, jvec)
   END SUBROUTINE

   SUBROUTINE jmat_x_b(m_loc, xv, jmat)
      INTEGER(KIND=4), INTENT(IN) :: m_loc
      REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xv
      REAL(KIND=8), INTENT(OUT), DIMENSION(n_max + 1, SIZE(xv)) :: jmat
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: abc_mn_loc
      INTEGER(KIND=4) :: n_min

      n_min = 3
      IF (m_loc > 1) THEN
         n_min = 1
      END IF

      ALLOCATE (abc_mn_loc(n_max - n_min, 4))
      CALL build_abc_mn(m_loc, abc_mn_loc)
      CALL jmat_x_raw(abc_mn_loc, m_loc, xv, jmat)
   END SUBROUTINE

   SUBROUTINE jmat_u_x_b(m_loc, uv, xv, jmat)
      INTEGER(KIND=4) :: m_loc
      REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xv
      REAL(KIND=8), INTENT(IN), DIMENSION(SIZE(xv)) :: uv
      REAL(KIND=8), INTENT(OUT), DIMENSION(n_max + 1, SIZE(xv)) :: jmat
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: abc_mn_loc
      INTEGER(KIND=4) :: n_min

      n_min = 3
      IF (m_loc > 1) THEN
         n_min = 1
      END IF

      ALLOCATE (abc_mn_loc(n_max - n_min, 4))
      CALL build_abc_mn(m_loc, abc_mn_loc)
      CALL jmat_u_x_raw(abc_mn_loc, m_loc, uv, xv, jmat)
   END SUBROUTINE

   SUBROUTINE jvec_u_x(u, x, jvec) BIND(c)

                !!!!!!!!!!!!!!
      !
      ! Builds the asymmetric Jacobi P polynomial as defined in [2] A.1 for
      ! all a vector of x values.
      !
                !!!!!!!!!!!!!!

      REAL(KIND=8), INTENT(IN) :: x, u
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(OUT) :: jvec
      REAL(KIND=8) :: vn, vn_, vntmp, upm, n_
      INTEGER(KIND=4) :: i

      jvec = 0

      IF (m_maxasymjacob == 1) THEN
         jvec(1) = 0.5
         jvec(2) = 1 - x/2
         jvec(3) = (3 + x*(-12 + 8*x))/6
         jvec(4) = (5 + x*(-60 + (120 - 64*x)*x))/10
         vn_ = jvec(3)
         vn = jvec(4)
      ELSE
         jvec(1) = 0.5
         jvec(2) = u*((REAL(m_maxasymjacob) - 0.5) + (1 - REAL(m_maxasymjacob))*x)
         vn_ = jvec(1)
         vn = jvec(2)
      END IF

      DO i = 1, shapeAbc_mn(1)
         n_ = INT(abc_mn(i, 1) + 1, KIND=4)
         IF (n_ < m_maxasymjacob) THEN
            vntmp = vn
            vn = u*(((abc_mn(i, 2) + abc_mn(i, 3)*x)*vn) - &
                    (u*abc_mn(i, 4)*vn_))
            vn_ = vntmp
         ELSE IF (n_ == m_maxasymjacob) THEN
            vntmp = vn
            vn = (((abc_mn(i, 2) + abc_mn(i, 3)*x)*vn) - &
                  (u*abc_mn(i, 4)*vn_))
            vn_ = vntmp
         ELSE
            vntmp = vn
            vn = (((abc_mn(i, 2) + abc_mn(i, 3)*x)*vn) - &
                  (abc_mn(i, 4)*vn_))
            vn_ = vntmp
         END IF
         jvec(INT(abc_mn(i, 1) + 2, KIND=4)) = vn
      END DO
      upm = u**(MAX(0, m_maxasymjacob - 1 - (n_max)))
      DO i = MIN(n_max, m_maxasymjacob - 1), 0, -1
         jvec(i + 1) = jvec(i + 1)*upm
         upm = upm*u
      END DO
   END SUBROUTINE
END MODULE asymjacobip

MODULE qmnp
   ! Fortran rewrite of qspectre from python package Scikit-qfit
   USE math_util
   USE qmn_array_index
   USE asymjacobip
   ! USE BLAS95
   ! USE F95_PRECISION
   INTEGER(kind=4) :: m_max, k_max, j_max, coefcount_qmnp, n_table, m_table
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: smFn, smGn, smHn
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: bgA, bgB, bgC, bgF, bgG
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: smF, smG
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: smFt, smGt, bgAt, bgBt, bgCt
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: mv0, mv1, mv2
   INTEGER(KIND=4), ALLOCATABLE :: table_name_cache(:), table_index_cache(:, :)

   CHARACTER(10) ::  stringoutint
CONTAINS

   SUBROUTINE qmnp_init(m_input, n_input) BIND(c)

      INTEGER(KIND=4), INTENT(IN), VALUE :: n_input, m_input
      INTEGER(KIND=4) :: n0mv, i
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmpmv0, tmpmv1, tmpmv2

      IF ((.NOT. ((m_max == m_input) .AND. (n_max == (n_input + 5)))) .AND. (error == 0)) THEN
         CALL asymjacobip_init(n_input + INT(5, KIND=4))
         m_max = MAX(m_input, 3)
         k_max = n_max + 2
         j_max = m_max + 1
         n_table = n_max + 6
         m_table = m_max + 1

         IF (.NOT. ALLOCATED(mv0)) THEN
            ALLOCATE (mv0(m_max + 1))
         ELSE
            n0mv = SIZE(mv0)
            IF (.NOT. (n0mv == (m_max + 1))) THEN
               ALLOCATE (tmpmv0(m_max + 1))
               CALL MOVE_ALLOC(FROM=tmpmv0, TO=mv0)
            END IF
         END IF

         IF (.NOT. ALLOCATED(mv1)) THEN
            ALLOCATE (mv1(m_max))
         ELSE
            n0mv = SIZE(mv1)
            IF (.NOT. (n0mv == (m_max))) THEN
               ALLOCATE (tmpmv1(m_max))
               CALL MOVE_ALLOC(FROM=tmpmv1, TO=mv1)
            END IF
         END IF

         IF (.NOT. ALLOCATED(mv2)) THEN
            ALLOCATE (mv2(m_max - 2))
         ELSE
            n0mv = SIZE(mv2)
            IF (.NOT. (n0mv == (m_max - 2))) THEN
               ALLOCATE (tmpmv2(m_max - 2))
               CALL MOVE_ALLOC(FROM=tmpmv2, TO=mv2)
            END IF
         END IF

         mv0 = (/(REAL(i, KIND=8), i=0, m_max)/)
         mv1 = (/(REAL(i, KIND=8), i=1, m_max)/)
         mv2 = (/(REAL(i, KIND=8), i=2, m_max)/)

         ! pre compute tables from [3] A.14, A.15, A.16
         CALL compute_qbfs_tables
         CALL compute_qinv_tables
         CALL compute_freeform_tables
         CALL compute_az_tables
         CALL compute_qmn_array_tables
      END IF
   END SUBROUTINE

   SUBROUTINE compute_qmn_array_tables
      INTEGER(KIND=4) :: n0table_name_cache, i
      INTEGER(KIND=4), DIMENSION(2) :: n0table_index_cache
      INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:)::tmptable_name_cache
      INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:, :)::tmptable_index_cache

      IF (.NOT. ALLOCATED(table_name_cache)) THEN
         ALLOCATE (table_name_cache(coefcount_qmnp))
      ELSE
         n0table_name_cache = SIZE(table_name_cache)
         IF (.NOT. (n0table_name_cache == (coefcount_qmnp))) THEN
            ALLOCATE (tmptable_name_cache(coefcount_qmnp))
            CALL MOVE_ALLOC(FROM=tmptable_name_cache, TO=table_name_cache)
         END IF
      END IF

      IF (.NOT. ALLOCATED(table_index_cache)) THEN
         ALLOCATE (table_index_cache(coefcount_qmnp, INT(2, KIND=4)))
      ELSE
         n0table_index_cache = SHAPE(table_index_cache)
         IF (.NOT. ALL(n0table_index_cache == [coefcount_qmnp, INT(2, KIND=4)])) THEN
            ALLOCATE (tmptable_index_cache(coefcount_qmnp, INT(2, KIND=4)))
            CALL MOVE_ALLOC(FROM=tmptable_index_cache, TO=table_index_cache)
         END IF
      END IF

      DO i = 1, coefcount_qmnp
         CALL num2qmnindex(i, table_name_cache(i), table_index_cache(i, :))
      END DO

   END SUBROUTINE

   SUBROUTINE data2arrfast_c(DATA, datasize, amnloc, bmnloc, arrlocsize1, arrlocsize2) BIND(c, name="data2arrfast")
      INTEGER(KIND=4) :: datasize, arrlocsize1, arrlocsize2
      REAL(KIND=8), DIMENSION(datasize), INTENT(IN) :: DATA
      REAL(KIND=8), DIMENSION(arrlocsize1, arrlocsize2), INTENT(OUT) :: amnloc, bmnloc
      CALL data2arrfast(DATA, amnloc, bmnloc)
   END SUBROUTINE

   SUBROUTINE data2arrfast(DATA, amnloc, bmnloc)
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: DATA
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(OUT) :: amnloc, bmnloc
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1, 2) :: cmnloc

      cmnloc = 0
      DO i = 1, MIN(coefcount_qmnp, SIZE(DATA))
         cmnloc(table_index_cache(i, 1), table_index_cache(i, 2), table_name_cache(i) + 1) = DATA(i)
      END DO
      amnloc = 0
      amnloc = cmnloc(:, :, 1)
      bmnloc = 0
      bmnloc = cmnloc(:, :, 2)
   END SUBROUTINE

   FUNCTION qmnp_get_nm_from_coef(coefcount_qmnp_input) RESULT(ind)
      INTEGER(KIND=4), INTENT(IN) :: coefcount_qmnp_input
      INTEGER(KIND=4) :: maxab
      INTEGER(KIND=4), DIMENSION(2) :: ind
      maxab = num2qmnmax(coefcount_qmnp_input) - 1
      ind = (/MAX(3, maxab), MAX(5, maxab)/)
   END FUNCTION

   SUBROUTINE qmnp_init_from_coef(coefcount_qmnp_input) BIND(C)
      INTEGER(KIND=4), INTENT(IN) :: coefcount_qmnp_input
      INTEGER(KIND=4), DIMENSION(2) :: ind
      IF (.NOT. (coefcount_qmnp_input == coefcount_qmnp)) THEN
         coefcount_qmnp = coefcount_qmnp_input
         ind = qmnp_get_nm_from_coef(coefcount_qmnp_input)
         CALL qmnp_init(ind(1), ind(2))
      END IF
   END SUBROUTINE

   SUBROUTINE compute_az_tables

      ! REAL(KIND=8), DIMENSION(n_table, m_max + 3) :: smFt, smGt, bgAt, bgBt, bgCt
      INTEGER(KIND=4), DIMENSION(2) :: n0tab
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: tmptab
      INTEGER(KIND=4) :: i

      IF (.NOT. ALLOCATED(smFt)) THEN
         ALLOCATE (smFt(n_table, m_max + 3))
      ELSE
         n0tab = SHAPE(smFt)
         IF (.NOT. ALL(n0tab == [n_table, m_max + INT(3, KIND=4)])) THEN
            ALLOCATE (tmptab(n_table, m_max + 3))
            CALL MOVE_ALLOC(FROM=tmptab, TO=smFt)
         END IF
      END IF

      IF (.NOT. ALLOCATED(smGt)) THEN
         ALLOCATE (smGt(n_table, m_max + 3))
      ELSE
         n0tab = SHAPE(smGt)
         IF (.NOT. ALL(n0tab == [n_table, m_max + INT(3, KIND=4)])) THEN
            ALLOCATE (tmptab(n_table, m_max + 3))
            CALL MOVE_ALLOC(FROM=tmptab, TO=smGt)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bgAt)) THEN
         ALLOCATE (bgAt(n_table, m_max + 3))
      ELSE
         n0tab = SHAPE(bgAt)
         IF (.NOT. ALL(n0tab == [n_table, m_max + INT(3, KIND=4)])) THEN
            ALLOCATE (tmptab(n_table, m_max + 3))
            CALL MOVE_ALLOC(FROM=tmptab, TO=bgAt)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bgBt)) THEN
         ALLOCATE (bgBt(n_table, m_max + 3))
      ELSE
         n0tab = SHAPE(bgBt)
         IF (.NOT. ALL(n0tab == [n_table, m_max + INT(3, KIND=4)])) THEN
            ALLOCATE (tmptab(n_table, m_max + 3))
            CALL MOVE_ALLOC(FROM=tmptab, TO=bgBt)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bgCt)) THEN
         ALLOCATE (bgCt(n_table, m_max + 3))
      ELSE
         n0tab = SHAPE(bgCt)
         IF (.NOT. ALL(n0tab == [n_table, m_max + INT(3, KIND=4)])) THEN
            ALLOCATE (tmptab(n_table, m_max + 3))
            CALL MOVE_ALLOC(FROM=tmptab, TO=bgCt)
         END IF
      END IF

      smFt = TRANSPOSE(smF(2:, :))
      smGt = TRANSPOSE(smG(2:, :))
      bgAt = TRANSPOSE(bgA(2:, :))
      bgBt = TRANSPOSE(bgB(2:, :))
      bgCt = TRANSPOSE(bgC(2:, :))

   END SUBROUTINE

   SUBROUTINE compute_qbfs_tables
      INTEGER(KIND=4) :: n0smFn, n0smGn, n0smHn
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmpsmFn, tmpsmGn, tmpsmHn
      INTEGER(KIND=4) :: i

      ! Reallocate memory for smFn, smGn and smHn
      IF (.NOT. ALLOCATED(smFn)) THEN
         ALLOCATE (smFn(n_max + 3))
      ELSE
         n0smFn = SIZE(smFn)
         IF (.NOT. (n0smFn == (n_max + 3))) THEN
            ALLOCATE (tmpsmFn(n_max + 3))
            CALL MOVE_ALLOC(FROM=tmpsmFn, TO=smFn)
         END IF
      END IF

      IF (.NOT. ALLOCATED(smGn)) THEN
         ALLOCATE (smGn(n_max + 2))
      ELSE
         n0smGn = SIZE(smGn)
         IF (.NOT. (n0smGn == (n_max + 2))) THEN
            ALLOCATE (tmpsmGn(n_max + 2))
            CALL MOVE_ALLOC(FROM=tmpsmGn, TO=smGn)
         END IF
      END IF

      IF (.NOT. ALLOCATED(smHn)) THEN
         ALLOCATE (smHn(n_max + 1))
      ELSE
         n0smHn = SIZE(smHn)
         IF (.NOT. (n0smHn == (n_max + 1))) THEN
            ALLOCATE (tmpsmHn(n_max + 1))
            CALL MOVE_ALLOC(FROM=tmpsmHn, TO=smHn)
         END IF
      END IF

      smFn = 0
      smGn = 0
      smHn = 0

      !COMPUTE TABLES
      smFN(1) = 2
      smFN(2) = SQRT(19.0D+00)/2
      smGn(1) = -0.5

      DO i = 3, n_max + 3
         smHn(i - 2) = -(i - 1)*(i - 2)/(2*smFn(i - 2))
         smGn(i - 1) = -(1 + smGn(i - 2)*smHn(i - 2))/(smFn(i - 1))
         smFn(i) = SQRT((i - 1)*(i) + 3 - smGn(i - 1)**2 - smHn(i - 2)**2)
      END DO
   END SUBROUTINE

   SUBROUTINE compute_qinv_tables
                !!!!!!!!!!!!!!
      !
      !     Build the big A, B and C tables as described in [2] (A.3a,b,c,d)
      !
                !!!!!!!!!!!!!

      INTEGER(KIND=4), DIMENSION(2) :: n0bgA, n0bgB, n0bgC
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: tmpbgA, tmpbgB, tmpbgC
      INTEGER(KIND=4) :: i
      REAL(KIND=8), DIMENSION(m_max) :: dv
      ! REAL(KIND=8), DIMENSION(m_max) :: mv1
      REAL(KIND=8) :: d
      REAL(KIND=8) :: ireal

      ! mv1 = (/(REAL(i, KIND=8), i=1, m_max)/)

      ! Reallocate memory for bgA, bgB and bgC
      IF (.NOT. ALLOCATED(bgA)) THEN
         ALLOCATE (bgA(m_table, n_table))
      ELSE
         n0bgA = SHAPE(bgA)
         IF (.NOT. ALL(n0bgA == [m_table, n_table])) THEN
            ALLOCATE (tmpbgA(m_table, n_table))
            CALL MOVE_ALLOC(FROM=tmpbgA, TO=bgA)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bgB)) THEN
         ALLOCATE (bgB(m_table, n_table))
      ELSE
         n0bgB = SHAPE(bgB)
         IF (.NOT. ALL(n0bgB == [m_table, n_table])) THEN
            ALLOCATE (tmpbgB(m_table, n_table))
            CALL MOVE_ALLOC(from=tmpbgB, to=bgB)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bgC)) THEN
         ALLOCATE (bgC(m_table, n_table))
      ELSE
         n0bgC = SHAPE(bgC)
         IF (.NOT. ALL(n0bgC == [m_table, n_table])) THEN
            ALLOCATE (tmpbgC(m_table, n_table))
            CALL MOVE_ALLOC(FROM=tmpbgC, TO=bgC)
         END IF
      END IF

      bgA = 0
      bgB = 0
      bgC = 0

      DO i = 3, n_table
         ireal = REAL(i - 1, KIND=8)
         dv = (4.D0*ireal*ireal - 1.D0)*(mv1 + ireal - 2.D0)* &
              (mv1 + 2.D0*ireal - 3.D0)
         bgA(2:, i) = (2.D0*ireal - 1.D0)*(mv1 + 2.D0*ireal - 2.D0)* &
                      (4.D0*ireal*(mv1 + ireal - 2.D0) + (mv1 - 3.D0)* &
                       (2.D0*mv1 - 1.D0))/dv
         bgB(2:, i) = -2.D0*(2.D0*ireal - 1.D0)*(mv1 + 2.D0*ireal - 1.D0)* &
                      (mv1 + 2.D0*ireal - 2.D0)*(mv1 + 2.D0*ireal - 3.D0)/dv
         bgC(2:, i) = ireal*(2.D0*ireal - 3.D0)*(mv1 + 2.D0*ireal - 1.D0)* &
                      (2.D0*mv1 + 2.D0*ireal - 3.D0)/dv
      END DO

      ! Initialze the special cases using [2] B.7 and B.8

      DO i = 3, m_table
         ireal = REAL(i - 1)
         bgA(i, 1) = 2.D0*ireal - 1.D0
         bgB(i, 1) = 2.D0*(1.D0 - ireal)
         d = 3.D0*(ireal - 1.D0)**2.D0
         bgA(i, 2) = ireal*(4.D0*(ireal - 1.D0) + (ireal - 3.D0)* &
                            (2.D0*ireal - 1.D0))/d
         bgB(i, 2) = -2.D0*(ireal - 1.D0)*ireal*(ireal + 1.D0)/d
         bgC(i, 2) = -(ireal + 1.D0)*(2.D0*ireal - 1.D0)/d
      END DO

      bgA(2, 1) = 2.D0
      bgB(2, 1) = -1.D0
      bgA(2, 2) = -4.D0/3.D0
      bgB(2, 2) = -8.D0/3.D0
      bgC(2, 2) = -11.D0/3.D0
      bgC(2, 3) = 0.D0

   END SUBROUTINE

   PURE FUNCTION gamma_factorial(m, n) RESULT(p)
      INTEGER(KIND=4), INTENT(IN) ::  m, n
      REAL(KIND=8) :: p
      INTEGER(KIND=4) :: i
      IF (m == 0) THEN
         p = REAL(n*(n - 1)*(n - 2), KIND=8)/ &
             REAL(2*(2*n - 1), KIND=8)
      ELSE IF (m == 1) THEN
         p = REAL(n*(n - 1), KIND=8)/4
      ELSE IF (m == 2) THEN
         p = REAL(n*(2*n + 1), KIND=8)/8
      ELSE IF (m == 3) THEN
         p = REAL((2*n + 1)*(2*n + 3), KIND=8)/16
      ELSE
         p = 2**(-REAL(m, KIND=8) - 1)
         DO i = 1, m - 1
            p = p*(2*n + (2*i - 1))
         END DO
         DO i = 1, m - 3
            p = p/(n + i)
         END DO
      END IF
   END FUNCTION

   SUBROUTINE compute_freeform_tables

      INTEGER(KIND=4), DIMENSION(2) :: n0bgF, n0bgG, n0smF, n0smG
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: tmpbgF, tmpbgG, tmpsmF, tmpsmG
      INTEGER(KIND=4) :: i, icalc
      REAL(KIND=8) :: ireal
      REAL(KIND=8), DIMENSION(m_max + 1) ::  fvF, fvG, gv_m, gamma_array
      REAL(KIND=8), DIMENSION(m_max - 1) ::  mv2_sqrd
      REAL(KIND=8) :: facF, facG, g_m

      ! mv0 = (/(REAL(i, KIND=8), i=0, m_max)/)
      ! mv2 = (/(REAL(i, KIND=8), i=2, m_max)/)
      mv2_sqrd = mv2*mv2
      fvF = 1
      fvG = 1
      gv_m = 1
      gamma_array = 0

      ! Reallocate memory for bgF, bgG, smF and smG
      IF (.NOT. ALLOCATED(bgF)) THEN
         ALLOCATE (bgF(m_table, n_table))
      ELSE
         n0bgF = SHAPE(bgF)
         IF (.NOT. ALL(n0bgF == [m_table, n_table])) THEN
            ALLOCATE (tmpbgF(m_table, n_table))
            CALL MOVE_ALLOC(FROM=tmpbgF, TO=bgF)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bgG)) THEN
         ALLOCATE (bgG(m_table, n_table))
      ELSE
         n0bgG = SHAPE(bgG)
         IF (.NOT. ALL(n0bgG == [m_table, n_table])) THEN
            ALLOCATE (tmpbgG(m_table, n_table))
            CALL MOVE_ALLOC(FROM=tmpbgG, TO=bgG)
         END IF
      END IF

      IF (.NOT. ALLOCATED(smF)) THEN
         ALLOCATE (smF(m_table, n_table))
      ELSE
         n0smF = SHAPE(smF)
         IF (.NOT. ALL(n0smF == [m_table, n_table])) THEN
            ALLOCATE (tmpsmF(m_table, n_table))
            CALL MOVE_ALLOC(FROM=tmpsmF, TO=smF)
         END IF
      END IF

      IF (.NOT. ALLOCATED(smG)) THEN
         ALLOCATE (smG(m_table, n_table))
      ELSE
         n0smG = SHAPE(smG)
         IF (.NOT. ALL(n0smG == [m_table, n_table])) THEN
            ALLOCATE (tmpsmG(m_table, n_table))
            CALL MOVE_ALLOC(FROM=tmpsmG, TO=smG)
         END IF
      END IF

      bgF = 0
      bgG = 0
      smF = 0
      smG = 0

      DO i = 2, m_table
         ireal = i - 1
         IF (i == 2) THEN
            facF = 0.25
            facG = 0.25
            g_m = 0.25
         ELSE
            facF = 0.5*((2.0D+00*ireal - 3.0D+00)/(ireal - 1.0D+00))*facF
            facG = 0.5*((2.0D+00*ireal - 1.0D+00)/(ireal - 1.0D+00))*facG
            IF (i > 4) THEN
               g_m = g_m*REAL((2.0D+00*ireal - 3.0D+00)/(2.0D+00*(ireal - 3.0D+00)), kind=8)
            ELSE
               g_m = REAL(3.0D+00/(2**4))
            END IF
            fvF(i) = facF
            fvG(i) = facG
            gv_m(i) = g_m
         END IF
      END DO
      DO i = 1, n_table
         ireal = REAL(i - 1, KIND=8)
         IF (i == 1) THEN
            gamma_array(4) = gamma_factorial(3, 0)
            gamma_array(5:) = gv_m(5:)
         ELSE
            icalc = MAX(1, 6 - i)
            IF (icalc > 1) THEN
               gamma_array(icalc - 1) = gamma_factorial(icalc - INT(2, KIND=4), i - INT(1, KIND=4))
            END IF
            gamma_array(icalc:) = (ireal*(2.D0*mv0(icalc:) + (2.D0*ireal - 3.D0))/((mv0(icalc:) + &
                                                                           (ireal - 3.D0))*(2.D0*ireal - 1.D0)))*gamma_array(icalc:)
         END IF
         IF (i == 1) THEN
            bgF(2, i) = 0.25
            bgG(2, i) = 0.25
            bgF(3:, i) = mv2_sqrd*fvF(3:)
            bgG(3:, i) = fvG(3:)
         ELSE
            icalc = i - 1
            bgF(2, i) = (4.0*((ireal - 1.0)*ireal)**2.0 + 1.0)/(8.0*(2.0*ireal - 1.0)**2.0) + &
                        REAL(kron_delta(i - INT(1, KIND=4), 1), KIND=8)*11.0/32.0
            bgG(2, i) = -(((2.0*ireal*ireal - 1.0)*(ireal*ireal - 1.0))/ &
                          (8.0*(4.0*ireal*ireal - 1.0))) - REAL(kron_delta(i - INT(1, KIND=4), INT(1, KIND=4)), kind=8)/24.0
            bgF(3:, i) = (REAL(2*icalc*(mv2 + (icalc - 2))*(3 - 5*mv2 + 4*icalc*(mv2 + &
                                                        (icalc - 2))) + mv2_sqrd*(3 - mv2 + 4*icalc*(mv2 + (icalc - 2))), kind=8)/ &
                          REAL((2*icalc - 1)*(mv2 + (2*icalc - 3))*(mv2 + (2*icalc - 2))* &
                               (mv2 + (2*icalc - 1)), KIND=8))*gamma_array(3:)
            bgG(3:, i) = -(((2*ireal*(mv2 + (ireal - 1.)) - mv2)*(ireal + 1)* &
                            (2*mv2 + (2*ireal - 1)))/((mv2 + (2*ireal - 2))* &
                                                      (mv2 + (2*ireal - 1))*(mv2 + 2*ireal)*(2*ireal + 1)))*gamma_array(3:)
         END IF
      END DO

      smF(:, 1) = SQRT(bgF(:, 1))
      DO i = 2, n_table
         smG(2:, i - 1) = bgG(2:, i - 1)/smF(2:, i - 1)
         smF(2:, i) = SQRT((bgF(2:, i) - (smG(2:, i - 1))**2))
      END DO
   END SUBROUTINE

   PURE SUBROUTINE radial_sum_int(an, r, radialoutput) BIND(c)
      REAL(KIND=8), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(IN) :: an
      REAL(KIND=8), DIMENSION(2), INTENT(OUT) :: radialoutput
      INTEGER(KIND=4) :: i
      REAL(KIND=8) :: b, b_, tmpb
      REAL(KIND=8) :: alpha, alpha_, tmpalpha
      REAL(KIND=8) :: afp, afp_, tmpafp
      REAL(KIND=8) :: t_4r

      ! if ( size(an) < n_max+1 ) then
      !         write ( *, '(a)' ) ' '
      !         write ( *, '(a)' ) 'RADIAL_SUM - Fatal error!'
      !         write ( *, '(a,g14.6)' ) '  Illegal input value of an = ', an
      !         write ( *, '(a,g14.6,a)' ) '  But an must be greater than n_max+1 = ', n_max+1, '.'
      !         stop
      ! end if

      t_4r = 2.0 - 4.0*r
      b_ = (r*0 + 1)*(an(n_max + 1)/smFn(n_max + 1))
      b = (an(n_max) - smGn(n_max)*b_)/smFn(n_max)
      alpha_ = b_
      alpha = b + t_4r*alpha_
      afp_ = 0*r!arrzero
      afp = -4*alpha_

      DO i = n_max - 1, 1, -1
         tmpb = b
         tmpafp = afp
         tmpalpha = alpha
         b = (an(i) - smGn(i)*b - smHn(i)*b_)/smFn(i)
         alpha = b + (t_4r*alpha) - alpha_
         afp = (t_4r*afp) - afp_ - (4.0D+00*tmpalpha)
         afp_ = tmpafp
         b_ = tmpb
         alpha_ = tmpalpha
      END DO

      radialoutput(1) = 2.0D+00*(alpha + alpha_)
      radialoutput(2) = 2.0D+00*(afp + afp_)
   END SUBROUTINE

   PURE SUBROUTINE roll_u_int(u_cnt, u, upj)
      INTEGER(KIND=4), INTENT(INOUT) :: u_cnt
      REAL(KIND=8), INTENT(INOUT) :: u(:), upj(:)
      INTEGER :: cnd, cnd2

      cnd = (.NOT. ALL(upj == 0))
      cnd2 = (u_cnt > 0)
      u = CSHIFT(u, SHIFT=cnd2, DIM=1)
      u(1) = -1*(cnd2) + u(1)*(1 + cnd2)
      ! u = EOSHIFT(u, SHIFT=cnd2, DIM=1, BOUNDARY=1.0D+00)
      upj = -1*upj*u*cnd + upj*(1 + cnd)
      u_cnt = MAX(u_cnt - 1, 0)
   END SUBROUTINE

   ! PURE SUBROUTINE roll_u_int(u_cnt, u, upj)
   !    INTEGER(KIND=4), INTENT(INOUT) :: u_cnt
   !    REAL(KIND=8), INTENT(INOUT) :: u(:), upj(:)
   !    IF (u_cnt > 0) THEN
   !       u = CSHIFT(u, SHIFT=-1, DIM=1)
   !       u(1) = 1
   !       IF (.NOT. ALL(upj == 0)) THEN
   !          upj = upj*u
   !       END IF
   !       u_cnt = u_cnt - 1
   !    END IF
   ! END SUBROUTINE

   FUNCTION azimuthal_sum_centre(cn) RESULT(k)

      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(IN):: cn
      REAL(KIND=8) :: k
      INTEGER(KIND=4) :: i
      ! REAL(KIND=8), DIMENSION(n_max + 1) :: cnm!, smFt, smGt, bgAt, bgBt, bgCt
      REAL(KIND=8) :: alpha, alpha_, dv, alpha_2, alpha_3

      ! cnm = cn(2, :)
      dv = cn(2, n_max + 1)/smFt(n_max + 1, 1)
      alpha = dv

      dv = (cn(2, n_max) - smGt(n_max, 1)*dv)/smFt(n_max, 1)
      alpha = dv + bgAt(n_max, 1)*alpha
      DO i = n_max - 1, 1, -1
         dv = (cn(2, i) - dv*smGt(i, 1))/smFt(i, 1)
         alpha_3 = alpha_2
         alpha_2 = alpha
         alpha = dv + bgAt(i, 1)*alpha - bgCt(i + 1, 1)*alpha_2
      END DO

      alpha = alpha - 8*alpha_3/10

      k = 0.5*alpha
   END FUNCTION

   SUBROUTINE azimuthal_sum_int_c(cn, r, outarr, outarrsize) BIND(C, name="azimuthal_sum_int")
      INTEGER(KIND=4), INTENT(IN) :: outarrsize
      REAL(KIND=8), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(IN):: cn
      REAL(KIND=8), INTENT(OUT), DIMENSION(2, outarrsize) :: outarr
      CALL azimuthal_sum_int(cn, r, outarr)
   END SUBROUTINE

   SUBROUTINE azimuthal_sum_int(cn, r, out_az)
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(IN):: cn
      REAL(KIND=8), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(2, m_max), INTENT(OUT) :: out_az
      INTEGER(KIND=4) :: uix, i
      REAL(KIND=8), DIMENSION(m_max) :: upj, u, dv
      REAL(KIND=8), DIMENSION(2, m_max) :: alphav, afpv
      ! REAL(KIND=8), DIMENSION(n_table, m_max + 3):: temp2
      REAL(KIND=8) ::alpha2, afp2, alpha3, afp3

      u = SQRT(r)
      upj = u
      uix = m_max

      alphav = 0
      afpv = 0
      dv = 0
      DO i = n_table - 1, n_max + 1, -1
         CALL roll_u_int(uix, u, upj)
      END DO

      DO i = n_max + 1, 1, -1
         CALL roll_u_int(uix, u, upj)
         dv = (cn(2:, i) - dv*smGt(i, :))/smFt(i, :)
         alpha3 = alpha2
         afp3 = afp2
         afp2 = afpv(2, 1)
         alpha2 = alphav(2, 1)
         afpv(2, :) = u*(alphav(1, :)*bgBt(i, :) + &
                         (bgAt(i, :) + (r*bgBt(i, :)))*afpv(1, :) - &
                         ((bgCt(i + 1, :)*u)*afpv(2, :)))
         alphav(2, :) = (dv*upj) + &
                        u*((bgAt(i, :) + (r*bgBt(i, :)))*(alphav(1, :)) - &
                           ((bgCt(i + 1, :)*u)*alphav(2, :)))
         afpv = afpv(2:1:-1, :)
         alphav = alphav(2:1:-1, :)
      END DO

      alphav(1, 1) = alphav(1, 1) - (8.D0*alpha3/10.D0)
      afpv(1, 1) = afpv(1, 1) - (8.D0*afp3/10.D0)

      DO i = 0, uix, 1
         u = CSHIFT(u, SHIFT=-1, DIM=1)
         u(1) = 1
         alphav(1, :) = alphav(1, :)*u
         afpv(1, :) = afpv(1, :)*u
      END DO

      out_az(1, :) = 0.5D0*alphav(1, :)
      out_az(2, :) = 0.5D0*afpv(1, :)
   END SUBROUTINE

   ! SUBROUTINE azimuthal_sum_int(cn, r, out_az, outsize)
   !    REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(IN):: cn
   !    REAL(KIND=8), INTENT(IN) :: r
   !    INTEGER(KIND=4), INTENT(IN) :: outsize
   !    REAL(KIND=8), DIMENSION(2, outsize), INTENT(OUT) :: out_az
   !    INTEGER(KIND=4) :: uix, i, n_loc
   !    REAL(KIND=8), DIMENSION(m_max) :: upj, u , dv
   !    REAL(KIND=8), DIMENSION(m_max, 4) :: alphav, afpv
   !    REAL(KIND=8), DIMENSION(n_table, m_max + 3) :: cnm
   !    REAL(KIND=8), DIMENSION(n_table, m_max+3):: temp2, dvm
   !    REAL(KIND=8), DIMENSION(m_max) :: temp1
   !    ! REAL(KIND=8), DIMENSION(n_table, m_max + 3) :: smFt, smGt, bgAt, bgBt, bgCt

   !    u = SQRT(r)
   !    temp2=bgAt+(r*bgBt)
   !    upj = u
   !    uix = m_max
   !    n_loc = n_table - 1

   !    cnm = 0
   !    DO i = 1, n_max + 1
   !       cnm(i, :) = cn(2:, i)
   !    END DO

   !    dvm(n_loc+1,:) = cnm(n_loc + 1, :)/smFt(n_loc + 1, :)
   !    do i=n_loc,1,-1
   !       dvm(i,:)= (cnm(i, :) - dvm(i+1,:)*smGt(i, :))/smFt(i, :)
   !    end do

   !    afpv(1,:)=0
   !    alphav(1,:)=upj*dvm(n_loc+1,:)

   !    CALL roll_u_int(uix, u, upj)
   !    afpv=cshift(afpv,-1,2)
   !    alphav=cshift(alphav,-1,2)
   !    afpv(1,:) = u*bgBt(n_loc, :)*alphav(2,:)
   !    alphav(1,:) = (upj*dvm(n_loc,:)) + u*(temp2(n_loc,:))*alphav(2,:)
   !    DO i = n_loc - 1, 1, -1
   !       CALL roll_u_int(uix, u, upj)
   !       temp1=u**2*bgCt(i + 1, :)
   !       afpv(:,4) = u*alphav(1,:)*bgBt(i, :) + (temp2(i,:)*u*afpv(1,:)) -&
   !                  (temp1*afpv(2,:))
   !       alphav(:,4) = (dvm(i,:)*upj) + temp2(i, :)* &
   !                     (u*alphav(1,:)) - (temp1*alphav(2,:))
   !       afpv=cshift(afpv,-1,2)
   !       alphav=cshift(alphav,-1,2)
   !    END DO

   !    alphav(1,1) = alphav(1,1) - (8.D0*alphav(1,4)/10.D0)
   !    afpv(1,1) = afpv(1,1) - (8.D0*afpv(1,4)/10.D0)

   !    upj = 0
   !    DO WHILE (uix .GT. 0)
   !       CALL roll_u_int(uix, u, upj)
   !       alphav(1,:) = alphav(1,:)*u
   !       afpv(1,:) = afpv(1,:)*u
   !    END DO

   !    out_az(1, :) = 0.5D0*alphav(1,:)
   !    out_az(2, :) = 0.5D0*afpv(1,:)
   ! END SUBROUTINE

   SUBROUTINE build_map_int2(x, y, curv, conic, radius, coef, output) BIND(c)
      REAL(KIND=8), INTENT(IN) :: x, y, curv, radius, conic
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: coef
      REAL(KIND=8), INTENT(OUT), DIMENSION(3) :: output

      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1) :: anm, bnm
      INTEGER(KIND=4) :: j
      REAL(KIND=8) :: rho, theta
      REAL(KIND=8) :: u, u2, val, radial, azr, asjk, sumjk, psijk, rhojk
      REAL(KIND=8) :: valp, psi, psi2, aztr, aztm_temp, aztr0, azt
      REAL(KIND=8), DIMENSION(2) :: radial_out
      REAL(KIND=8), DIMENSION(2, m_max) :: az_outa, az_outb
      REAL(KIND=8), DIMENSION(m_max + 1) :: mcosf, msinf, sinf, cosf, thetamv

      ! Builds a regular map where rho and theta are the axis values.

      ! Build the radial component first and expand for the theta values
      rho = SQRT(x**2 + y**2)

      u = (rho/radius)
      u2 = u**2

      CALL data2arrfast(coef, &
                        anm, bnm)

      IF (.NOT. ALL(bnm(2:, :) == 0.0D+00)) THEN
         CALL azimuthal_sum_int(bnm, u2, az_outb)
      ELSE
         az_outb = 0.0D+00
      END IF

      theta = ATAN2(y, x)

      IF (rho < 1E-15) THEN
         output(1) = 0.D0
         output(2) = COS(theta)*azimuthal_sum_centre(anm)/(radius)
         output(3) = SIN(theta)*DOT_PRODUCT(mv0(2:), az_outb(1, :))
         RETURN
      END IF
      rhojk = rho
      psijk = 1.0D+00 - (curv**2*rhojk**2)
      IF (psijk < 0.0D+00) THEN
         output = 0.0D+00
         RETURN
      END IF
      psijk = SQRT(psijk)
      IF (u .GE. 1) THEN
         output(1) = (curv*rhojk**2)/(1.0D+00 + psijk)
         output(2) = (curv*x)/psijk
         output(3) = (curv*y)/psijk
         RETURN
      END IF

      IF (.NOT. ALL(anm(2:, :) == 0.0D+00)) THEN
         CALL azimuthal_sum_int(anm, u2, az_outa)
      ELSE
         az_outa = 0.0D+00
      END IF

      IF (.NOT. ALL(anm(1, :) == 0.0D+00)) THEN
         CALL radial_sum_int(anm(1, :), u2, radial_out)
      ELSE
         radial_out = 0.0D+00
      END IF

      val = radial_out(1)*u2*(1.0D+00 - (1 - conic)*u2)

      thetamv = theta*mv0

      sinf = SIN(thetamv)
      cosf = COS(thetamv)

      asjk = DOT_PRODUCT(cosf(2:), az_outa(1, :)) + DOT_PRODUCT(sinf(2:), az_outb(1, :))

      radial = val

      sumjk = (radial + asjk)/psijk
      sumjk = sumjk + ((curv*rhojk**2)/(1.0D+00 + psijk))

      output(1) = sumjk

      !!! DERIVATIVE
      ! Build the derivative in rho and theta maps
      psi = SQRT(1.0D+00 - (curv**2*rho**2))
      psi2 = psi**2
      ! valp =
      ! valp = valp +
      ! valp = valp +

      mcosf = cosf*mv0
      msinf = sinf*mv0

      azr = radial_out(1)*u*(1.0 + psi2 - u2*(1.0 + 3.0*psi2))/(radius*psi*psi2) + &
            (radial_out(2)*2.0*u2*u*(1.0 - u2)/(radius*psi)) + &
            (curv*rho/psi)

      aztr = 2.0D+00*(DOT_PRODUCT(cosf(2:), az_outa(2, :)) + DOT_PRODUCT(sinf(2:), az_outb(2, :)))*u

      IF (u == 0.0D+00) THEN
         aztr = azimuthal_sum_centre(anm)*cosf(1) + azimuthal_sum_centre(bnm)*sinf(1)
      ELSE
         aztr = aztr + (DOT_PRODUCT(mcosf(2:), az_outa(1, :)) + DOT_PRODUCT(msinf(2:), az_outb(1, :)))/u
      END IF

      aztr = aztr/(radius*psijk)
      aztr = aztr + (curv**2*rhojk/psijk**2)*asjk/psijk

      azr = azr + aztr

      azt = (DOT_PRODUCT(mcosf(2:), az_outb(1, :)) - DOT_PRODUCT(msinf(2:), az_outa(1, :)))/psijk

      output(2) = COS(theta)*azr - (1/rho)*SIN(theta)*azt
      output(3) = SIN(theta)*azr + (1/rho)*COS(theta)*azt

   END SUBROUTINE

   SUBROUTINE build_map_int(rho, theta, curv, conic, radius, anm, bnm, output) BIND(c)
      REAL(KIND=8), INTENT(IN) :: rho, theta, curv, radius, conic
      REAL(KIND=8), INTENT(IN), DIMENSION(m_max + 1, n_max + 1) :: anm, bnm
      REAL(KIND=8), INTENT(OUT), DIMENSION(3) :: output

      INTEGER(KIND=4) :: j
      REAL(KIND=8) :: u, u2, val, radial, azr, asjk, sumjk, psijk, rhojk
      REAL(KIND=8) :: valp, psi, psi2, aztr, aztm_temp, aztr0, azt
      REAL(KIND=8), DIMENSION(2) :: radial_out
      REAL(KIND=8), DIMENSION(2, m_max) :: az_outa, az_outb
      REAL(KIND=8), DIMENSION(m_max + 1) :: mcosf, msinf, sinf, cosf, thetamv

      ! Builds a regular map where rho and theta are the axis values.

      ! Build the radial component first and expand for the theta values

      u = (rho/radius)
      u2 = u**2

      IF (.NOT. ALL(bnm(2:, :) == 0.0D+00)) THEN
         CALL azimuthal_sum_int(bnm, u2, az_outb)
      ELSE
         az_outb = 0.0D+00
      END IF

      IF (rho < 1E-15) THEN
         output(1) = 0.D0
         output(2) = azimuthal_sum_centre(anm)/(radius)
         output(3) = DOT_PRODUCT(mv0(2:), az_outb(1, :))
         RETURN
      END IF

      rhojk = rho
      psijk = 1.0D+00 - (curv**2*rhojk**2)

      IF (psijk < 0.0D+00) THEN
         output = 0.0D+00
         RETURN
      END IF

      psijk = SQRT(psijk)

      IF (u .GE. 1) THEN
         output(1) = (curv*rhojk**2)/(1.0D+00 + psijk)
         output(2) = (curv*rhojk)/psijk
         output(3) = 0.0D+00
         RETURN
      END IF

      IF (.NOT. ALL(anm(2:, :) == 0.0D+00)) THEN
         CALL azimuthal_sum_int(anm, u2, az_outa)
      ELSE
         az_outa = 0.0D+00
      END IF

      IF (.NOT. ALL(anm(1, :) == 0.0D+00)) THEN
         CALL radial_sum_int(anm(1, :), u2, radial_out)
      ELSE
         radial_out = 0.0D+00
      END IF

      val = radial_out(1)*u2*(1.0D+00 - (1 - conic)*u2)

      thetamv = theta*mv0

      sinf = SIN(thetamv)
      cosf = COS(thetamv)

      asjk = DOT_PRODUCT(cosf(2:), az_outa(1, :)) + DOT_PRODUCT(sinf(2:), az_outb(1, :))

      radial = val

      sumjk = (radial + asjk)/psijk
      sumjk = sumjk + ((curv*rhojk**2)/(1.0D+00 + psijk))

      output(1) = sumjk

      !!! DERIVATIVE
      ! Build the derivative in rho and theta maps
      psi = SQRT(1.0D+00 - (curv**2*rho**2))
      psi2 = psi**2
      ! valp =
      ! valp = valp +
      ! valp = valp +

      mcosf = cosf*mv0
      msinf = sinf*mv0

      azr = radial_out(1)*u*(1.0 + psi2 - u2*(1.0 + 3.0*psi2))/(radius*psi*psi2) + &
            (radial_out(2)*2.0*u2*u*(1.0 - u2)/(radius*psi)) + &
            (curv*rho/psi)

      aztr = 2.0D+00*(DOT_PRODUCT(cosf(2:), az_outa(2, :)) + DOT_PRODUCT(sinf(2:), az_outb(2, :)))*u

      IF (u == 0.0D+00) THEN
         aztr = azimuthal_sum_centre(anm)*cosf(1) + azimuthal_sum_centre(bnm)*sinf(1)
      ELSE
         aztr = aztr + (DOT_PRODUCT(mcosf(2:), az_outa(1, :)) + DOT_PRODUCT(msinf(2:), az_outb(1, :)))/u
      END IF

      aztr = aztr/(radius*psijk)
      aztr = aztr + (curv**2*rhojk/psijk**2)*asjk/psijk

      azr = azr + aztr

      azt = (DOT_PRODUCT(mcosf(2:), az_outb(1, :)) - DOT_PRODUCT(msinf(2:), az_outa(1, :)))/psijk

      output(2) = azr
      output(3) = azt

   END SUBROUTINE

   SUBROUTINE radial_sum_c(an, rn, r, output) BIND(C, name="radial_sum")
      INTEGER(KIND=4), INTENT(IN) :: rn
      REAL(KIND=8), DIMENSION(rn), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(IN) :: an
      REAL(KIND=8), DIMENSION(2, rn), INTENT(OUT) :: output

      CALL radial_sum(an, r, output)
   END SUBROUTINE

   SUBROUTINE radial_sum(an, r, output)
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(IN) ::  an
      REAL(KIND=8), DIMENSION(2, SIZE(r)), INTENT(OUT) :: output
      ! INTEGER(KIND=4) :: n0qmn_radial
      ! REAL(KIND=8), ALLOCATABLE :: tmpqmn_radial(:)
      REAL(KIND=8), DIMENSION(SIZE(r)) :: b, b_, tmpb, alpha, alpha_, tmpalpha
      REAL(KIND=8), DIMENSION(SIZE(r)) :: afp, afp_, tmpafp, t_4r
      INTEGER(KIND=4) :: i, j

      output = 0

      IF (SIZE(an) < n_max + 1) THEN
         WRITE (*, '(a)') ' '
         WRITE (*, '(a)') 'RADIAL_SUM - Fatal error!'
         WRITE (*, '(a,g14.6)') '  Illegal input value of an = ', an
         WRITE (*, '(a,g14.6,a)') '  But an must be greater than n_max+1 = ', n_max + 1, '.'
         STOP
      END IF

      t_4r = 2.0 - 4.0*r

      b_ = (r*0 + 1)*(an(n_max + 1)/smFn(n_max + 1))
      b = (an(n_max) - smGn(n_max)*b_)/smFn(n_max)
      alpha_ = b_
      alpha = b + t_4r*alpha_
      afp_ = 0*r!arrzero
      afp = -4*alpha_

      DO i = n_max - 1, 1, -1
         tmpb = b
         tmpafp = afp
         tmpalpha = alpha
         DO j = 1, SIZE(r)
            b(j) = (an(i) - smGn(i)*b(j) - smHn(i)*b_(j))/smFn(i)
            alpha(j) = b(j) + (t_4r(j)*alpha(j)) - alpha_(j)
            afp(j) = (t_4r(j)*afp(j)) - afp_(j) - (4.0D+00*tmpalpha(j))
         END DO
         afp_ = tmpafp
         b_ = tmpb
         alpha_ = tmpalpha
      END DO

      output(1, :) = 2.0D+00*(alpha + alpha_)
      output(2, :) = 2.0D+00*(afp + afp_)
   END SUBROUTINE

   SUBROUTINE azimuthal_sum_c(cn, nr, r, output) BIND(C, name="azimuthal_sum")
      INTEGER(KIND=4), INTENT(IN) :: nr
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(IN):: cn
      REAL(KIND=8), DIMENSION(nr), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(2, m_max, nr), INTENT(OUT) :: output
      CALL azimuthal_sum(cn, r, output)
   END SUBROUTINE

   SUBROUTINE azimuthal_sum(cn, r, output)
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(IN):: cn
      REAL(KIND=8), DIMENSION(2, m_max, SIZE(r)), INTENT(OUT) :: output

      INTEGER(KIND=4) :: uix, i, j, n_loc
      REAL(KIND=8), DIMENSION(SIZE(r), m_max) :: upj, u, u2, afp, afp_, alpha_3
      REAL(KIND=8), DIMENSION(SIZE(r), m_max) :: dv_, alpha_, dv, alpha, alpha_2
      REAL(KIND=8), DIMENSION(SIZE(r), m_max) :: afp_2, afp_3
      REAL(KIND=8), DIMENSION(n_table, m_max) :: cnm
      ! REAL(KIND=8), DIMENSION(n_table, m_max) :: smFt, smGt, bgAt, bgBt, bgCt
      REAL(KIND=8), DIMENSION(SIZE(r)) :: sr

      upj = 1
      sr = SQRT(r)

      DO i = 1, m_max
         u(:, i) = sr*upj(1, i)
      END DO

      upj = u*upj
      uix = m_max

      cnm = 0
      DO i = 1, n_max + 1
         cnm(i, :) = cn(2:, i)
      END DO
      ! smFt = TRANSPOSE(smF(2:, :))
      ! smGt = TRANSPOSE(smG(2:, :))
      ! bgAt = TRANSPOSE(bgA(2:, :))
      ! bgBt = TRANSPOSE(bgB(2:, :))
      ! bgCt = TRANSPOSE(bgC(2:, :))

      n_loc = n_table - 1

      DO i = 1, m_max
         dv_(:, i) = cnm(n_loc + 1, i)/smFt(n_loc + 1, i)
      END DO

      dv = 0
      alpha_ = upj*dv_
      afp_ = 0

      IF (n_loc > 0) THEN
         CALL roll_u(uix, u, upj)

         DO i = 1, m_max
            dv(:, i) = (cnm(n_loc, i) - smGt(n_loc, i)*dv_(:, i))/smFt(n_loc, i)
            alpha(:, i) = (upj(:, i)*dv(:, i)) + u(:, i)*(bgAt(n_loc, i) + (r*bgBt(n_loc, i)))*alpha_(:, i)
            afp(:, i) = u(:, i)*bgBt(n_loc, i)*alpha_(:, i)
         END DO

         DO i = n_loc - 1, 1, -1
            CALL roll_u(uix, u, upj)
            u2 = u*u
            alpha_3 = alpha_2
            alpha_2 = alpha_
            alpha_ = alpha
            afp_3 = afp_2
            afp_2 = afp_
            afp_ = afp
            DO j = 1, m_max
               dv(:, j) = (cnm(i, j) - dv(:, j)*smGt(i, j))/smFt(i, j)
               alpha(:, j) = (dv(:, j)*upj(:, j)) + (bgAt(i, j) + (r*bgBt(i, j)))* &
                             (u(:, j)*alpha(:, j)) - (bgCt(i + 1, j)*(u2(:, j)*alpha_2(:, j)))
               afp(:, j) = u(:, j)*alpha_(:, j)*bgBt(i, j) + (bgAt(i, j) + (r*bgBt(i, j)))* &
                           u(:, j)*afp(:, j) - (u2(:, j)*bgCt(i + 1, j)*afp_2(:, j))
            END DO
         END DO
         IF (n_loc > 2) THEN
            alpha(:, 1) = alpha(:, 1) - (8*alpha_3(:, 1)/10)
            afp(:, 1) = afp(:, 1) - (8*afp_3(:, 1)/10)
         END IF
         upj = upj*0
         DO WHILE (uix .GT. 0)
            CALL roll_u(uix, u, upj)
            alpha = alpha*u
            afp = afp*u
         END DO
      ELSE
         alpha = alpha_
         afp = afp_
      END IF

      output(1, :, :) = 0.5*TRANSPOSE(alpha)
      output(2, :, :) = 0.5*TRANSPOSE(afp)
   END SUBROUTINE

   PURE SUBROUTINE roll_u(u_cnt, u, upj)
      INTEGER(KIND=4), INTENT(INOUT) :: u_cnt
      REAL(KIND=8), INTENT(INOUT) :: u(:, :), upj(:, :)
      IF (u_cnt > 0) THEN
         u = CSHIFT(u, SHIFT=-1, DIM=2)
         u(:, 1) = 1
         IF (.NOT. ALL(upj == 0)) THEN
            upj = upj*u
         END IF
         u_cnt = u_cnt - 1
      END IF
   END SUBROUTINE

   SUBROUTINE build_map_c(rhon, rho, thetan, theta, curv, conic, radius, &
                          anmn, anmm, anm, bnmn, bnmm, bnm, output) BIND(C, name="build_map")
      INTEGER(KIND=4), INTENT(IN) ::  anmn, anmm, bnmm, bnmn, rhon, thetan
      REAL(KIND=8), INTENT(IN) :: curv, radius, conic
      REAL(KIND=8), DIMENSION(rhon), INTENT(IN) :: rho
      REAL(KIND=8), DIMENSION(thetan), INTENT(IN) :: theta
      REAL(KIND=8), DIMENSION(anmn, anmm), INTENT(IN) :: anm
      REAL(KIND=8), DIMENSION(bnmn, bnmm), INTENT(IN) :: bnm
      REAL(KIND=8), DIMENSION(3, rhon, thetan), INTENT(OUT) :: output
      CALL build_map(rho, theta, curv, conic, radius, anm, bnm, output)
   END SUBROUTINE

   SUBROUTINE build_map(rho, theta, curv, conic, radius, anm, bnm, output)

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rho, theta
      REAL(KIND=8), INTENT(IN), DIMENSION(m_max + 1, n_max + 1) :: anm, bnm
      REAL(KIND=8), INTENT(IN) :: curv, radius, conic
      REAL(KIND=8), DIMENSION(3, SIZE(rho), SIZE(theta)), INTENT(OUT) :: output

      INTEGER(KIND=4) :: j
      REAL(KIND=8), DIMENSION(2, m_max, SIZE(rho)) :: qmna_output_a, qmna_output_b
      REAL(KIND=8), DIMENSION(2, SIZE(rho)) :: qmnr_output
      REAL(KIND=8), DIMENSION(SIZE(rho)) :: u, u2, val, valp, psi, psi2
      REAL(KIND=8), DIMENSION(SIZE(theta)) :: aztr0
      REAL(KIND=8), DIMENSION(SIZE(theta), SIZE(rho)) :: radial, azr, asjk, sumjk, azt
      REAL(KIND=8), DIMENSION(SIZE(theta), SIZE(rho)) :: psijk, rhojk, aztr, aztm_temp
      REAL(KIND=8), DIMENSION(SIZE(theta), m_max + 1) :: thetamv, sinf, cosf, mcosf, msinf
      ! REAL(KIND=8), DIMENSION(m_max + 1) :: mv0

      ! Builds a regular map where rho and theta are the axis values.

      ! Build the radial component first and expand for the theta values

      u = (rho/radius)
      u2 = u**2
      CALL radial_sum(anm(1, :), u2, qmnr_output)
      val = qmnr_output(1, :)*u2*(1 - (1 + conic)*u2)

      ! The asymmetric terms as [m,k] matrices
      ! mv0 = (/(REAL(j, kind=8), j=0, m_max)/)
      !
      DO j = 1, SIZE(mv0)
         thetamv(:, j) = theta*mv0(j)
      END DO
      !
      sinf = SIN(thetamv)
      cosf = COS(thetamv)

      CALL azimuthal_sum(anm, u2, qmna_output_a)
      CALL azimuthal_sum(bnm, u2, qmna_output_b)

      CALL DGEMM('N','N',SIZE(theta),SIZE(rho),m_max,1.0D0,COSF(:,2:),SIZE(theta),qmna_output_a(1,:,:),m_max,0.0D0,asjk,SIZE(theta))
      CALL DGEMM('N','N',SIZE(theta),SIZE(rho),m_max,1.0D0,SINF(:,2:),SIZE(theta),qmna_output_b(1,:,:),m_max,1.0D0,asjk,SIZE(theta))

      ! asjk = matmul(cosf(:,2:),qmn_aza) + matmul(sinf(:,2:),qmn_azb)
      ! asjk = asjk1+asjk2
      ! write (*,*) asjk(2,4)
      ! asjk = matmul(cosf(:,2:),qmn_aza) + matmul(sinf(:,2:),qmn_azb)
      ! write (*,*) asjk(2,4)

      radial = SPREAD(val, 1, SIZE(theta))
      rhojk = SPREAD(rho, 1, SIZE(theta))

      psijk = SQRT(1 - (curv**2*rhojk**2))
      sumjk = (radial + asjk)/psijk
      sumjk = sumjk + ((curv*rhojk**2)/(1 + psijk))

      output(1, :, :) = TRANSPOSE(sumjk)

                !!! DERIVATIVE
      ! Build the derivative in rho and theta maps
      psi = SQRT(1 - (curv**2*rho**2))
      psi2 = psi**2
      valp = qmnr_output(1, :)*u*(1.0 + psi2 - u2*(1.0 + 3.0*psi2))/(radius*psi*psi2)
      valp = valp + (qmnr_output(2, :)*2.0*u2*u*(1.0 - u2)/(radius*psi))
      valp = valp + (curv*rho/psi)

      mcosf = cosf*SPREAD(mv0, 1, SIZE(theta))
      msinf = sinf*SPREAD(mv0, 1, SIZE(theta))

      CALL DGEMM('N', 'N', SIZE(theta), SIZE(rho), m_max, 1.0D0, cosf(:, 2:), SIZE(theta), &
                 qmna_output_a(2, :, :), m_max, 0.0D0, aztm_temp, SIZE(theta))
      CALL DGEMM('N', 'N', SIZE(theta), SIZE(rho), m_max, 1.0D0, sinf(:, 2:), SIZE(theta), &
                 qmna_output_b(2, :, :), m_max, 1.0D0, aztm_temp, SIZE(theta))
      ! aztm_temp = (matmul(cosf(:, 2:),dqmn_aza) + matmul(sinf(:, 2:),dqmn_azb))
      azr = SPREAD(valp, 1, SIZE(theta))
      aztr = 2*aztm_temp*SPREAD(u, 1, SIZE(theta))

      aztr0 = azimuthal_sum_centre(anm)*COS(theta) + azimuthal_sum_centre(bnm)*SIN(theta)
      CALL DGEMM('N', 'N', SIZE(theta), SIZE(rho), m_max, 1.0D0, mcosf(:, 2:), SIZE(theta), &
                 qmna_output_a(1, :, :), m_max, 0.0D0, aztm_temp, SIZE(theta))
      CALL DGEMM('N', 'N', SIZE(theta), SIZE(rho), m_max, 1.0D0, msinf(:, 2:), SIZE(theta), &
                 qmna_output_b(1, :, :), m_max, 1.0D0, aztm_temp, SIZE(theta))
      ! aztm_temp = matmul(mcosf(:, 2:),qmn_aza) + matmul(msinf(:, 2:),qmn_azb)
      !
      DO j = 1, SIZE(rho)
         IF (u(j) == 0) THEN
            aztr(:, j) = aztr0
         ELSE
            aztr(:, j) = aztr(:, j) + aztm_temp(:, j)/u(j)
         END IF
      END DO
      aztr = aztr/(radius*psijk)
      aztr = aztr + (curv**2*rhojk/psijk**2)*asjk/psijk
      azr = azr + aztr

      CALL DGEMM('N', 'N', SIZE(theta), SIZE(rho), m_max, 1.0D0, -msinf(:, 2:), SIZE(theta), &
                 qmna_output_a(1, :, :), m_max, 0.0D0, azt, SIZE(theta))
      CALL DGEMM('N', 'N', SIZE(theta), SIZE(rho), m_max, 1.0D0, mcosf(:, 2:), SIZE(theta), &
                 qmna_output_b(1, :, :), m_max, 1.0D0, azt, SIZE(theta))
      azt = azt/psijk
      ! azt = (matmul(-msinf(:,2:),qmn_aza)+matmul(mcosf(:,2:),qmn_azb))/psijk

      output(2, :, :) = TRANSPOSE(azr)
      output(3, :, :) = TRANSPOSE(azt)
   END SUBROUTINE
END MODULE qmnp

MODULE qfit
   ! Fortran rewrite of qspectre from python package Scikit-qfit
   USE math_util
   USE qmnp
   USE asymjacobip
   USE MKL_DFTI

   INTEGER(KIND=4), PARAMETER :: bfs_size = 50
   REAL(KIND=8), PARAMETER :: curv_bar = 0.1
   REAL(KIND=8) :: cmn_cur_length
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: phi_kvec, u_vec, u_vec_sqr
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rv_pos_der_norm, thv_pos_der_norm, r_pos_der_norm
   REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: bgK, bgH, smK, smH, smS, smT
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: qfitntable
   REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: amn_cur_ref, bmn_cur_ref
   ! REAL(KIND=8), DIMENSION(2*j_max*SIZE(u_vec)), INTENT(OUT):: rv, thv
   ! REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: arbar, brbar
! , qfitsagpolar
!         REAL(KIND=8) :: bfs_curv
!         REAL(KIND=8) :: qfitradius

CONTAINS

   SUBROUTINE qfit_init(m_input, n_input) BIND(c)
      INTEGER(KIND=4), INTENT(IN), VALUE :: n_input, m_input
      IF ((.NOT. ((m_max == m_input) .AND. (n_max == n_input))) .AND. (error == 0)) THEN
         CALL precompute_factor(m_input, n_input)
         CALL precompute_curvature_impact()
      END IF
   END SUBROUTINE

   SUBROUTINE qfit_init_from_coef(coefcount_qmnp_input)
      INTEGER(KIND=4), INTENT(IN) :: coefcount_qmnp_input
      INTEGER(KIND=4), DIMENSION(2) :: ind
      IF (.NOT. (coefcount_qmnp_input == coefcount_qmnp)) THEN
         ind = qmnp_get_nm_from_coef(coefcount_qmnp_input)
         CALL qfit_init(ind(1), ind(2))
      END IF
   END SUBROUTINE

   SUBROUTINE precompute_factor(m_input, n_input) BIND(c)
      INTEGER(KIND=4), INTENT(IN), VALUE :: m_input, n_input
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmpphi_kvec
      INTEGER(KIND=4) :: n0phi_kvec
      INTEGER(KIND=4) :: i

      CALL qmnp_init(m_input, n_input)

      IF (.NOT. ALLOCATED(phi_kvec)) THEN
         ALLOCATE (phi_kvec(k_max))
      ELSE
         n0phi_kvec = SIZE(phi_kvec)
         IF (.NOT. (n0phi_kvec == (k_max))) THEN
            ALLOCATE (tmpphi_kvec(k_max))
            CALL MOVE_ALLOC(FROM=tmpphi_kvec, TO=phi_kvec)
         END IF
      END IF

      IF (.NOT. ALLOCATED(u_vec)) THEN
         ALLOCATE (u_vec(k_max))
      ELSE
         n0phi_kvec = SIZE(u_vec)
         IF (.NOT. (n0phi_kvec == (k_max))) THEN
            ALLOCATE (tmpphi_kvec(k_max))
            CALL MOVE_ALLOC(FROM=tmpphi_kvec, TO=u_vec)
         END IF
      END IF

      IF (.NOT. ALLOCATED(u_vec_sqr)) THEN
         ALLOCATE (u_vec_sqr(k_max))
      ELSE
         n0phi_kvec = SIZE(u_vec_sqr)
         IF (.NOT. (n0phi_kvec == (k_max))) THEN
            ALLOCATE (tmpphi_kvec(k_max))
            CALL MOVE_ALLOC(FROM=tmpphi_kvec, TO=u_vec_sqr)
         END IF
      END IF
      IF (.NOT. ALLOCATED(r_pos_der_norm)) THEN
         ALLOCATE (r_pos_der_norm(k_max + 1))
      ELSE
         n0phi_kvec = SIZE(r_pos_der_norm)
         IF (.NOT. (n0phi_kvec == (k_max + 1))) THEN
            ALLOCATE (tmpphi_kvec(k_max + 1))
            CALL MOVE_ALLOC(FROM=tmpphi_kvec, TO=r_pos_der_norm)
         END IF
      END IF
      IF (.NOT. ALLOCATED(rv_pos_der_norm)) THEN
         ALLOCATE (rv_pos_der_norm(2*j_max*(k_max + 2)))
      ELSE
         n0phi_kvec = SIZE(rv_pos_der_norm)
         IF (.NOT. (n0phi_kvec == (2*j_max*(k_max + 2)))) THEN
            ALLOCATE (tmpphi_kvec(2*j_max*(k_max + 2)))
            CALL MOVE_ALLOC(FROM=tmpphi_kvec, TO=rv_pos_der_norm)
         END IF
      END IF
      IF (.NOT. ALLOCATED(thv_pos_der_norm)) THEN
         ALLOCATE (thv_pos_der_norm(2*j_max*(k_max + 2)))
      ELSE
         n0phi_kvec = SIZE(thv_pos_der_norm)
         IF (.NOT. (n0phi_kvec == (2*j_max*(k_max + 2)))) THEN
            ALLOCATE (tmpphi_kvec(2*j_max*(k_max + 2)))
            CALL MOVE_ALLOC(FROM=tmpphi_kvec, TO=thv_pos_der_norm)
         END IF
      END IF

      DO i = 1, k_max, 1
         phi_kvec(i) = (2.0D+00*REAL(i) - 1.0D+00)*qmnpi/(4.0D+00*REAL(k_max))
      END DO
      u_vec = SIN(phi_kvec)
      u_vec_sqr = u_vec**2

      CALL scan_pos_radial_der(r_pos_der_norm, 1.00D+00)
      CALL scan_pos_polar_der(rv_pos_der_norm, thv_pos_der_norm, 1.00D+00)

      CALL compute_qfit_tables()
      CALL compute_qfitder_table()

   END SUBROUTINE

   SUBROUTINE precompute_curvature_impact()
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: tmpcmn
      INTEGER(KIND=4), DIMENSION(2) :: n0arr
      REAL(KIND=8), DIMENSION(3, 2*j_max*(k_max + 2), 1) :: output_map
      REAL(KIND=8), DIMENSION(1) :: t_loc
      ! REAL(KIND=8) :: cmn_sum

      IF (.NOT. ALLOCATED(amn_cur_ref)) THEN
         ALLOCATE (amn_cur_ref(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(amn_cur_ref)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmpcmn(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmpcmn, to=amn_cur_ref)
         END IF
      END IF

      IF (.NOT. ALLOCATED(bmn_cur_ref)) THEN
         ALLOCATE (bmn_cur_ref(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(bmn_cur_ref)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmpcmn(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmpcmn, to=bmn_cur_ref)
         END IF
      END IF
      amn_cur_ref = 0
      bmn_cur_ref = 0
      t_loc = 0
      CALL build_map(rv_pos_der_norm, t_loc, &
                     curv_bar, 0.00D+00, 1.00D+00, &
                     amn_cur_ref, bmn_cur_ref, output_map)
      CALL q_fit_der(output_map(2, :, 1), output_map(3, :, 1), 1.00D+00, amn_cur_ref, bmn_cur_ref)
      cmn_cur_length = SUM(amn_cur_ref**2) + SUM(bmn_cur_ref**2)
      ! cmn_cur_length = SUM(amn_cur_ref(1,2:)**2)
      amn_cur_ref = amn_cur_ref
      bmn_cur_ref = bmn_cur_ref
   END SUBROUTINE

   SUBROUTINE compute_qfitder_table()
      INTEGER(KIND=4), DIMENSION(4) :: n0arr
      INTEGER(KIND=4) :: i, j
      REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: tmparr
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1) :: anm, bnm
      REAL(KIND=8), DIMENSION(3, k_max + 1, 1) :: output_map
      REAL(KIND=8), DIMENSION(1) :: t_loc

      ! REALLOCATION OF ARRAYS
      IF (.NOT. ALLOCATED(qfitntable)) THEN
         ALLOCATE (qfitntable(n_max + 1, m_max + 1, k_max + 1, 2))
      ELSE
         n0arr = SHAPE(qfitntable)
         IF (.NOT. (ALL(n0arr == [n_max + 1, m_max + 1, k_max + 1, 2]))) THEN
            ALLOCATE (tmparr(n_max + 1, m_max + 1, k_max + 1, 2))
            CALL MOVE_ALLOC(from=tmparr, to=qfitntable)
         END IF
      END IF

      t_loc = 0

      qfitntable = 0
      DO i = 1, n_max + 1
         DO j = 1, m_max + 1
            anm = 0.0D+00
            anm(j, i) = 1.0D+00
            bnm = 0.0D+00
            output_map = 0.0D+00
            CALL build_map(r_pos_der_norm, t_loc, &
                           0.00D+00, 0.00D+00, 1.00D+0, &
                           anm, bnm, output_map)
            qfitntable(i, j, :, 1) = output_map(2, :, 1)/(REAL(2*(k_max + 1), KIND=8))
     qfitntable(i, j, :, 2) = (REAL(j, KIND=8) - 1.0D+00)*(1.0D+00/r_pos_der_norm)*output_map(1, :, 1)/(REAL(2*(k_max + 1), KIND=8))
         END DO
      END DO

   END SUBROUTINE

   SUBROUTINE compute_qfit_tables()
      INTEGER(KIND=4), DIMENSION(2) :: n0arr
      REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: tmparr
      REAL(KIND=8) :: nfact, ir
      REAL(KIND=8), DIMENSION(m_max) :: v, w
      ! REAL(KIND=8), DIMENSION(m_max) :: mv1
      REAL(KIND=8), DIMENSION(n_max - 1) :: nv, nv2
      REAL(KIND=8), DIMENSION(n_max) :: snv, sn2v
      REAL(KIND=8), DIMENSION(m_max + 1) :: nfv

      ! REALLOCATION OF ARRAYS
      IF (.NOT. ALLOCATED(bgK)) THEN
         ALLOCATE (bgK(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(bgK)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmparr(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmparr, to=bgK)
         END IF
      END IF
      IF (.NOT. ALLOCATED(bgH)) THEN
         ALLOCATE (bgH(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(bgH)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmparr(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmparr, to=bgH)
         END IF
      END IF
      IF (.NOT. ALLOCATED(smK)) THEN
         ALLOCATE (smK(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(smK)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmparr(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmparr, to=smK)
         END IF
      END IF
      IF (.NOT. ALLOCATED(smH)) THEN
         ALLOCATE (smH(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(smH)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmparr(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmparr, to=smH)
         END IF
      END IF
      IF (.NOT. ALLOCATED(smS)) THEN
         ALLOCATE (smS(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(smS)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmparr(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmparr, to=smS)
         END IF
      END IF
      IF (.NOT. ALLOCATED(smT)) THEN
         ALLOCATE (smT(m_max + 1, n_max + 1))
      ELSE
         n0arr = SHAPE(smT)
         IF (.NOT. (ALL(n0arr == [m_max + 1, n_max + 1]))) THEN
            ALLOCATE (tmparr(m_max + 1, n_max + 1))
            CALL MOVE_ALLOC(from=tmparr, to=smT)
         END IF
      END IF

      bgK = 0.0D+00
      bgH = 0.0D+00
      bgK(1, 1) = 3.0D+00/8.0D+00
      bgK(1, 2) = 1.0D+00/24.0D+00
      bgH(1, 1) = 1.0D+00/4.0D+00
      bgH(1, 2) = 19.0D+00/32.0D+00

      ! mv1 = (/(REAL(i, KIND=8), i=1, m_max)/)
      nv = (/(REAL(i, KIND=8), i=2, n_max)/)
      nv2 = nv**2

      bgK(1, 3:) = (nv2 - 1.0D+00)/(32.0D+00*nv2 - 8.0D+00)
      bgH(1, 3:) = (1.0D+00 + 1.0D+00/(1.0D+00 - 2.0D+00*nv)**2)/16.0D+00

      nfv = (/(REAL(i, KIND=8), i=0, m_max + 1)/)

      nfact = 0.5D+00
      DO i = 1, m_max, 1
         nfact = (REAL((2*i) + 1, KIND=8)/REAL((2*i) + 2, KIND=8))*nfact
         nfv(i + 1) = nfact
      END DO

      bgK(2:, 1) = 0.5D+00*nfv(2:)
      bgK(2:, 2) = ((2.0D+00*mv1*(2.0D+00*mv1 + 3.0D+00))/ &
                    (3.0D+00*(mv1 + 3.0D+00)*(mv1 + 2.0D+00))) &
                   *0.5D+00*nfv(2:)
      bgH(2:, 1) = ((mv1 + 1.0D+00)/(2.0D+00*mv1 + 1.0D+00))*0.5D+00*nfv(2:)
      bgH(2:, 2) = ((3.0D+00*mv1 + 2.0D+00)/(mv1 + 2.0D+00))*0.5D+00*nfv(2:)

      v = bgK(2:, 2)
      w = bgH(2:, 2)
      DO i = 2, n_max
         ir = REAL(i, KIND=8)
         bgH(2:, i + 1) = (((mv1 + (2.0D+00*ir - 3.0D+00))* &
                            ((mv1 + (ir - 2.0D+00))*(4.0D+00*ir - 1.0D+00) + 5.0D+00*ir))/ &
                           ((mv1 + (ir - 2.0D+00))*(2.0D+00*ir - 1.0D+00)*(mv1 + 2.0D+00*ir)) &
                           )*v
         v = (((ir + 1.0D+00) &
               *(mv1 + (2.0D+00*ir - 2.0D+00)) &
               *(mv1 + (2.0D+00*ir - 3.0D+00)) &
               *(2.0D+00*mv1 + (2.0D+00*ir + 1.0D+00))) &
              /((2.0D+00*ir + 1.0D+00)* &
                (mv1 + (ir - 2.0D+00))* &
                (mv1 + (2.0D+00*ir + 1.0D+00))* &
                (mv1 + 2.0D+00*ir)))*v
         bgK(2:, i + 1) = v
      END DO

      smK = 0
      smH = 0

      smH(:, 1) = SQRT(bgH(:, 1))
      DO i = 1, n_max
         smK(:, i) = bgK(:, i)/smH(:, i)
         smH(:, i + 1) = SQRT(bgH(:, i + 1) - smK(:, i)**2)
      END DO

      smS = 0
      smT = 0

      snv = (/(REAL(i, KIND=8), i=1, n_max)/)
      sn2v = 2.0D+00*snv
      DO i = 1, m_max
         ir = REAL(i, KIND=8)
         smS(i + 1, 1) = 1.0D+00
         smT(i + 1, 1) = 1.0D+00/ir
         smS(i + 1, 2:) = (snv + (ir - 2.0D+00))/(sn2v + (ir - 2.0D+00))
         smT(i + 1, 2:) = ((1.0D+00 - sn2v)*(snv + 1.0D+00))/((ir + sn2v)*(sn2v + 1.0D+00))
      END DO
      smS(2, 2) = 0.5D+00
      smT(2, 1) = 0.5D+00
   END SUBROUTINE

   SUBROUTINE scan_pos_polar(rv, thv, qfitradius) BIND(c)
      REAL(KIND=8), DIMENSION(2*j_max*SIZE(u_vec)), INTENT(OUT):: rv, thv
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(2*j_max) :: scan_theta
      INTEGER(KIND=4) :: i
      scan_theta = (/(REAL(2*qmnpi*(i - 1), KIND=8)/REAL(2*j_max, KIND=8), i=1, 2*j_max)/)
      DO i = 1, SIZE(u_vec)
         rv(2*j_max*(i - 1) + 1:2*j_max*(i)) = u_vec(i)
      END DO
      rv = qfitradius*rv
      DO i = 1, SIZE(u_vec)
         thv(2*j_max*(i - 1) + 1:2*j_max*(i)) = scan_theta
      END DO
   END SUBROUTINE

   SUBROUTINE scan_pos_polar_der(rv, thv, qfitradius) BIND(c)
      REAL(KIND=8), DIMENSION(2*j_max*(k_max + 1)), INTENT(OUT):: rv, thv
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(2*j_max) :: scan_theta
      INTEGER(KIND=4) :: i
      scan_theta = (/(REAL(2*qmnpi*(i - 1), KIND=8)/REAL(2*j_max, KIND=8), i=1, 2*j_max)/)
      DO i = 1, k_max + 1
         rv(2*j_max*(i - 1) + 1:2*j_max*(i)) = COS((2*(i) - 1)*qmnpi/(4*(k_max + 1)))
      END DO
      rv = qfitradius*rv
      DO i = 1, k_max + 1
         thv(2*j_max*(i - 1) + 1:2*j_max*(i)) = scan_theta
      END DO
   END SUBROUTINE

   SUBROUTINE scan_pos_radial_der(rv, qfitradius)
      REAL(KIND=8), DIMENSION(k_max + 1):: rv
      REAL(KIND=8) :: qfitradius
      INTEGER(KIND=4) :: i
      DO i = 1, k_max + 1
         rv(i) = COS((2*(i) - 1)*qmnpi/(4*(k_max + 1)))*qfitradius
      END DO
   END SUBROUTINE

   SUBROUTINE scan_pos_car(xv, yv, xva, yva, qfitradius) BIND(c)
      REAL(KIND=8), DIMENSION(2*j_max*SIZE(u_vec)), INTENT(OUT):: xv, yv
      REAL(KIND=8), DIMENSION(2*j_max*SIZE(u_vec)):: rv, thv
      REAL(KIND=8), DIMENSION(bfs_size), INTENT(OUT):: xva, yva
      REAL(KIND=8), INTENT(IN) :: qfitradius
      INTEGER(KIND=4) :: i
      REAL(KIND=8), DIMENSION(bfs_size):: theta

      CALL scan_pos_polar(rv, thv, qfitradius)
      xv = rv*COS(thv)
      yv = rv*SIN(thv)

      theta = (/(REAL(2*qmnpi*(i - 1), KIND=8)/REAL(bfs_size, KIND=8), i=1, bfs_size)/)
      xva = qfitradius*COS(theta)
      yva = qfitradius*SIN(theta)
   END SUBROUTINE

   SUBROUTINE scan_pos_car_der(xv, yv, qfitradius) BIND(c)
      REAL(KIND=8), DIMENSION(2*j_max*(k_max + 1)), INTENT(OUT):: xv, yv
      REAL(KIND=8), DIMENSION(2*j_max*(k_max + 1)):: rv, thv
      ! REAL(KIND=8), DIMENSION(bfs_size), INTENT(OUT):: xva, yva
      REAL(KIND=8), INTENT(IN) :: qfitradius
      INTEGER(KIND=4) :: i
      REAL(KIND=8), DIMENSION(bfs_size):: theta

      CALL scan_pos_polar_der(rv, thv, qfitradius)
      xv = rv*COS(thv)
      yv = rv*SIN(thv)

      ! theta = (/(REAL(2*qmnpi*(i - 1), KIND=8)/REAL(bfs_size, KIND=8), i=1, bfs_size)/)
      ! xva = qfitradius*COS(theta)
      ! yva = qfitradius*SIN(theta)
   END SUBROUTINE

   SUBROUTINE get_bfs_curv(datainputrim, qfitradius, bfs_curv) BIND(c)
      REAL(KIND=8), DIMENSION(bfs_size), INTENT(IN) :: datainputrim
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), INTENT(OUT) :: bfs_curv
      REAL(KIND=8) :: sag_rim

      sag_rim = REAL(SUM(datainputrim), KIND=8)/(bfs_size)
      bfs_curv = 2.0D+00*sag_rim/(sag_rim**2 + qfitradius**2)
   END SUBROUTINE

   SUBROUTINE normal_departure(rv, datainput, bfs_curv, intp) BIND(c)
      REAL(KIND=8), DIMENSION(2*j_max*k_max), INTENT(IN) :: rv, datainput
      REAL(KIND=8), INTENT(IN) :: bfs_curv
      REAL(KIND=8), DIMENSION(2*j_max*k_max), INTENT(OUT) :: intp
      REAL(KIND=8), DIMENSION(2*j_max*k_max) :: rv2, fact

      rv2 = rv**2
      fact = SQRT(1.0D+00 - bfs_curv**2*rv2)
      intp = fact*(datainput - bfs_curv*rv2/(1.0D+00 + fact))
   END SUBROUTINE

   SUBROUTINE build_azimuthal_fit(datainputr, datainputt, qfitradius, &
                                  abar, bbar, abat, bbat, bfs_in) BIND(c)
      REAL(KIND=8), DIMENSION((k_max + 1)*2*j_max), INTENT(IN) :: datainputr, datainputt
      REAL(KIND=8), DIMENSION(j_max, k_max + 1), INTENT(OUT):: abar, bbar, abat, bbat
      REAL(KIND=8), DIMENSION((k_max + 1)*2*j_max) :: datainputrdep
      REAL(KIND=8), INTENT(IN) :: bfs_in
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(2*j_max) :: xfft
      INTEGER(KIND=4), DIMENSION(m_max + 1) :: scan_m_0
      COMPLEX(KIND=8), DIMENSION(2*j_max) :: X_out
      COMPLEX(KIND=8), DIMENSION(k_max + 1, 2*j_max) :: X_in
      TYPE(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle
      INTEGER :: status

      ! CALL scan_pos_polar(rv, thv, qfitradius)
      ! CALL normal_departure_der(rv, datainput, bfs_in, intpf)
      ! bfs_in = bfs_curv

      ! intp(2:, :) = TRANSPOSE(RESHAPE(intpf, (/2*j_max, k_max/)))

      ! intp(1, :) = 0

      ! scan_m_0 = (/(i - 1, i=1, m_max + 2)/)
      abar = 0D+00
      bbar = 0D+00
      abat = 0D+00
      bbat = 0D+00
      datainputrdep = datainputr

      status = DftiCreateDescriptor(My_Desc_Handle, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 1, 2*j_max)
      status = DftiSetValue(My_Desc_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      status = DftiCommitDescriptor(My_Desc_Handle)

      X_in = TRANSPOSE(CMPLX(RESHAPE(datainputrdep, (/(2*j_max), k_max + 1/)), KIND=8))
      DO i = 1, k_max + 1
         X_out = 0
         Status = DftiComputeForward(My_Desc_Handle, &
                                     X_in(i, :), &
                                     X_out)
         abar(:m_max + 1, i) = REAL((X_out(:m_max + 1)/j_max), KIND=8)
         bbar(:m_max + 1, i) = -REAL(AIMAG((X_out(:m_max + 1)/j_max)), KIND=8)
      END DO

      X_in = TRANSPOSE(CMPLX(RESHAPE((1/rv_pos_der_norm)*datainputt, (/2*j_max, k_max + 1/)), KIND=8))
      DO i = 1, k_max + 1
         X_out = 0
         Status = DftiComputeForward(My_Desc_Handle, &
                                     X_in(i, :), &
                                     X_out)
         bbat(:m_max + 1, i) = REAL((X_out(:m_max + 1)/j_max), KIND=8)
         abat(:m_max + 1, i) = REAL(AIMAG((X_out(:m_max + 1)/j_max)), KIND=8)
      END DO

      status = DftiFreeDescriptor(My_Desc_Handle)

   END SUBROUTINE

   SUBROUTINE build_abr_bar(datainput, qfitradius, &
                            arbar, brbar, bfs_in) BIND(c)
      REAL(KIND=8), DIMENSION(k_max*2*j_max), INTENT(IN) :: datainput
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(OUT):: arbar, brbar
      REAL(KIND=8), INTENT(IN) :: bfs_in
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(2*j_max) :: scan_theta
      REAL(KIND=8), DIMENSION(2*j_max*k_max) :: rv, thv
      REAL(KIND=8), DIMENSION(k_max*2*j_max) :: intpf
      REAL(KIND=8), DIMENSION(k_max + 1, 2*j_max) :: intp
      REAL(KIND=8), DIMENSION(2*j_max) :: xfft
      INTEGER(KIND=4), DIMENSION(m_max + 1) :: scan_m_0
      REAL(KIND=8), DIMENSION(m_max + 1, k_max + 1) :: abar, bbar
      REAL(KIND=8), DIMENSION(n_max + 1, k_max) :: jmat
      COMPLEX(KIND=8), DIMENSION(2*j_max) :: X_in, X_out
      TYPE(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle
      INTEGER :: status

      CALL scan_pos_polar(rv, thv, qfitradius)
      CALL normal_departure(rv, datainput, bfs_in, intpf)

      intp(2:, :) = TRANSPOSE(RESHAPE(intpf, (/INT(2*j_max, KIND=4), k_max/)))

      intp(1, :) = 0

      scan_m_0 = (/(i - 1, i=1, m_max + 2)/)
      abar = 0D+00
      bbar = 0D+00
      arbar = 0D+00
      brbar = 0D+00

      status = DftiCreateDescriptor(My_Desc_Handle, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 1, 2*j_max)
      status = DftiSetValue(My_Desc_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      status = DftiCommitDescriptor(My_Desc_Handle)

      DO i = 2, k_max + 1
         X_out = 0
         Status = DftiComputeForward(My_Desc_Handle, CMPLX(intp(i, :), KIND=8), X_out)
         abar(:m_max + 1, i) = REAL((X_out(:m_max + 1)/j_max), KIND=8)
         bbar(:m_max + 1, i) = -REAL(AIMAG((X_out(:m_max + 1)/j_max)), KIND=8)
      END DO

      status = DftiFreeDescriptor(My_Desc_Handle)

      DO i = 1, SIZE(scan_m_0)
         CALL jmat_u_x_b(scan_m_0(i) + INT(1, KIND=4), u_vec, u_vec_sqr, jmat)
         CALL DGEMV('N', &
                    SIZE(jmat(:, 1)), &
                    SIZE(jmat(1, :)), &
                    1D+00, &
                    jmat/k_max, &
                    SIZE(jmat(:, 1)), &
                    abar(scan_m_0(i) + 1, 2:), &
                    1, &
                    0D+00, &
                    arbar(scan_m_0(i) + 1, :), &
                    1)
         CALL DGEMV('N', &
                    SIZE(jmat(:, 1)), &
                    SIZE(jmat(1, :)), &
                    1D+00, &
                    jmat/k_max, &
                    SIZE(jmat(:, 1)), &
                    bbar(scan_m_0(i) + 1, 2:), &
                    1, &
                    0D+00, &
                    brbar(scan_m_0(i) + 1, :), &
                    1)
      END DO
   END SUBROUTINE

   SUBROUTINE rbar_to_cbar(rbar, e_bar_0) BIND(c)
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(OUT):: rbar
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(OUT):: e_bar_0
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1) :: sigma_bar, e_bar, d_bar, c_bar
      INTEGER(KIND=4), DIMENSION(m_max + 1) :: scan_m_0
      INTEGER(KIND=4) :: i

      scan_m_0 = (/(i - 1, i=1, m_max + 2)/)
      sigma_bar = 0
      sigma_bar(:, 1) = rbar(:, 1)/smH(1:m_max + 1, 1)
      DO i = 2, n_max + 1
         sigma_bar(:, i) = (rbar(:, i) - smK(1:m_max + 2, i - 1)*sigma_bar(:, i - 1))/ &
                           smH(:m_max + 1, i)
      END DO
      e_bar = 0
      e_bar(:, n_max + 1) = sigma_bar(:, n_max + 1)/smH(1:m_max + 2, n_max + 1)
      DO i = n_max, 1, -1
         e_bar(:, i) = (sigma_bar(:, i) - smK(1:m_max + 2, i)*e_bar(:, i + 1))/ &
                       smh(:m_max + 2, i)
      END DO
      e_bar_0 = e_bar(1, :)

      d_bar = 0
      d_bar(2:, n_max + 1) = e_bar(2:, n_max + 1)/smS(2:m_max + 1, n_max + 1)
      DO i = n_max, 1, -1
         d_bar(2:, i) = (e_bar(2:, i) - smT(2:m_max + 1, i)*d_bar(2:, i + 1))/ &
                        smS(2:m_max + 1, i)
      END DO

      c_bar = 0
      DO i = 1, n_max
         c_bar(2:, i) = (smF(2:m_max + 1, i)*d_bar(2:, i) + &
                         smG(2:m_max + 1, i)*d_bar(2:, i + 1))
      END DO
      c_bar(2:, n_max + 1) = smF(2:m_max + 1, n_max + 1)*d_bar(2:, n_max + 1)
      rbar = c_bar

   END SUBROUTINE

   FUNCTION e_rot_sym_fit(u, e_bar_0) RESULT(uout)
      REAL(KIND=8), INTENT(IN) :: u
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(IN):: e_bar_0
      REAL(KIND=8) :: uout
      REAL(KIND=8), DIMENSION(n_max + 1) :: jvecloc
      CALL jvec_x_b(1, u**2, jvecloc)
      uout = DOT_PRODUCT(jvecloc, e_bar_0)/2
   END FUNCTION

   SUBROUTINE build_avec(e_bar_0, av) BIND(c)
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(IN):: e_bar_0
      REAL(KIND=8), DIMENSION(n_max + 1), INTENT(OUT):: av
      REAL(KIND=8), DIMENSION(n_max + 1, k_max) :: jmatloc
      REAL(KIND=8), DIMENSION(k_max) :: svec, dct
      REAL(KIND=8), DIMENSION(k_max + 1) :: bv
      INTEGER(KIND=4) :: i

      CALL jmat_x_b(1, u_vec_sqr, jmatloc)

      CALL DGEMV('T', &
                 SIZE(jmatloc(:, 1)), &
                 SIZE(jmatloc(1, :)), &
                 0.5D+00, &
                 jmatloc, &
                 SIZE(jmatloc(:, 1)), &
                 e_bar_0, &
                 1, &
                 0.0D+00, &
                 svec, &
                 -1)
      svec = svec - e_rot_sym_fit(0.0D+00, e_bar_0) - &
             (e_rot_sym_fit(1.0D+00, e_bar_0) - &
              e_rot_sym_fit(0.0D+00, e_bar_0))*u_vec(SIZE(u_vec):1:-1)**2
      svec(:k_max + 1) = svec(:k_max + 1)*COS(phi_kvec)/ &
                         (u_vec_sqr*(1 - u_vec_sqr))

      dct = dct_iv(svec)

      bv = 0
      DO i = 1, k_max
         bv(i) = REAL((-1)**(i - 1), KIND=8)*dct(i)
      END DO
      bv = bv*1.0D+00/SQRT(REAL(2*k_max, KIND=8))
      av = 0
      DO i = 1, n_max + 1
         av(i) = bv(i)*smFn(i) + &
                 bv(i + 1)*smGn(i) + &
                 bv(i + 2)*smHn(i)
      END DO
   END SUBROUTINE

   SUBROUTINE q_fit_t(datainput, datainputrim, qfitradius, amn, bmn, bfs_out) BIND(c)
      REAL(KIND=8), DIMENSION(k_max*2*j_max), INTENT(IN) :: datainput
      REAL(KIND=8), DIMENSION(bfs_size), INTENT(IN) :: datainputrim
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(n_max + 1, m_max + 1), INTENT(out):: amn, bmn
      REAL(KIND=8), INTENT(OUT) :: bfs_out
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1) :: arbar, brbar
      REAL(KIND=8), DIMENSION(n_max + 1) :: e_bar_0
      REAL(KIND=8) :: bfs_curv, max_datainput

      max_datainput = MAXVAL(datainput)
      max_datainput = 1.0

      CALL get_bfs_curv(datainputrim, &
                        qfitradius, bfs_curv)
      CALL build_abr_bar(datainput, &
                         qfitradius, arbar, brbar, bfs_curv)
      bfs_out = bfs_curv
      CALL rbar_to_cbar(arbar, e_bar_0)
      CALL build_avec(e_bar_0, arbar(1, :))
      CALL rbar_to_cbar(brbar, e_bar_0)
      amn = arbar
      bmn = brbar
   END SUBROUTINE

   SUBROUTINE q_fit_curv_t(datainput, qfitradius, bfs_curv, amn, bmn) BIND(c)
      REAL(KIND=8), DIMENSION(k_max*2*j_max), INTENT(IN) :: datainput
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(n_max + 1, m_max + 1), INTENT(out):: amn, bmn
      REAL(KIND=8), INTENT(IN) :: bfs_curv
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1) :: arbar, brbar
      REAL(KIND=8), DIMENSION(n_max + 1) :: e_bar_0

      CALL build_abr_bar(datainput, qfitradius, arbar, brbar, bfs_curv)
      bfs_out = bfs_curv
      CALL rbar_to_cbar(arbar, e_bar_0)
      CALL build_avec(e_bar_0, arbar(1, :))
      CALL rbar_to_cbar(brbar, e_bar_0)
      amn = arbar
      bmn = brbar
   END SUBROUTINE

   SUBROUTINE q_fit(datainput, datainputrim, qfitradius, amn, bmn, bfs_out) BIND(c)
      REAL(KIND=8), DIMENSION(k_max*2*j_max), INTENT(IN) :: datainput
      REAL(KIND=8), DIMENSION(bfs_size), INTENT(IN) :: datainputrim
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(n_max + 1, m_max + 1) :: amn2, bmn2
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(out):: amn, bmn
      REAL(KIND=8), INTENT(OUT) :: bfs_out
      REAL(KIND=8) :: bfs_curv

      CALL q_fit_t(datainput, datainputrim, qfitradius, amn2, bmn2, bfs_curv)
      bfs_out = bfs_curv
      amn = TRANSPOSE(amn2)
      bmn = TRANSPOSE(bmn2)
   END SUBROUTINE

   SUBROUTINE q_fit_der(datainputr, datainputt, qfitradius, amn, bmn) BIND(C)
      REAL(KIND=8), DIMENSION((k_max + 1)*2*j_max), INTENT(IN) :: datainputr, datainputt
      REAL(KIND=8), INTENT(IN) :: qfitradius

      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(out) :: amn, bmn
      ! REAL(KIND=8), DIMENSION(m_max + 1):: tmparr

      REAL(KIND=8), DIMENSION(j_max, k_max + 1) :: abar, bbar, abat, bbat
      REAL(KIND=8) :: bfs_in

      amn = 0
      bmn = 0

      CALL build_azimuthal_fit(datainputr, datainputt, qfitradius, &
                               abar, bbar, abat, bbat, bfs_in)

      DO i = 1, n_max + 1
         amn(:, i) = (SUM(qfitntable(i, :, :, 1)*abar, DIM=2) + SUM(qfitntable(i, :, :, 2)*abat, DIM=2))
         bmn(:, i) = (SUM(qfitntable(i, :, :, 1)*bbar, DIM=2) + SUM(qfitntable(i, :, :, 2)*bbat, DIM=2))
      END DO
   END SUBROUTINE

   SUBROUTINE q_fit_remove_curv(amn, bmn, output_curv) BIND(C)
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(inout) :: amn, bmn
      REAL(KIND=8) :: output_curv, curv_fact
      curv_fact = (SUM(amn*amn_cur_ref) + SUM(bmn*bmn_cur_ref))/cmn_cur_length
      ! curv_fact = (SUM(amn(1,2:)*amn_cur_ref(1,2:)))/cmn_cur_length
      output_curv = curv_fact*curv_bar
      amn = amn - (amn_cur_ref*curv_fact)
      bmn = bmn - (bmn_cur_ref*curv_fact)
   END SUBROUTINE

   SUBROUTINE q_fit_curv(datainput, qfitradius, bfs_curv, amn, bmn) BIND(c)
      REAL(KIND=8), DIMENSION(k_max*2*j_max), INTENT(IN) :: datainput
      REAL(KIND=8), INTENT(IN) :: qfitradius
      REAL(KIND=8), DIMENSION(n_max + 1, m_max + 1) :: amn2, bmn2
      REAL(KIND=8), DIMENSION(m_max + 1, n_max + 1), INTENT(out):: amn, bmn
      REAL(KIND=8), INTENT(IN) :: bfs_curv

      CALL q_fit_curv_t(datainput, qfitradius, bfs_curv, amn2, bmn2)
      amn = TRANSPOSE(amn2)
      bmn = TRANSPOSE(bmn2)
   END SUBROUTINE

END MODULE
