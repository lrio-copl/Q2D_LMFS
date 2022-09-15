!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!      ####################################################################
!      ####################################################################
!
!         Copyright (C) 2020 Synopsys, Inc. This Synopsys product and all associated
!         documentation are proprietary to Synopsys, Inc. and may only be used
!         pursuant to the terms and conditions of a written license agreement
!         with Synopsys, Inc. All other use, reproduction, modification, or
!         distribution of the Synopsys product or the associated documentation
!         is strictly prohibited.
!
!      ####################################################################
!      ####################################################################

!DEC$ ATTRIBUTES DLLEXPORT :: uds_qmn2
SUBROUTINE uds_qmn2(x, y, z, curv, f, fx, fy, fz, error, mode, coef, coefCount) BIND(C)

   ! ****************************************************************
   ! * Purpose: Evaluates the function and its derivatives describing
   ! *          the surface in the form F(x,y,z) = 0.  Note that if you
   ! *          have a surface in the form z = f(x,y), then
   ! *          F(x,y,z) = f(x,y)-z = 0;  in this case FX = df/dx,
   ! *          FY = df/dy, FZ = -1.
   ! *
   ! * Parameters:
   ! *   The following is a brief description of the parameters in the
   ! *   call list.  If the parameter is designated as "input", its value
   ! *   is passed to the subroutine by the calling program;  if it is
   ! *   designated as "output", its value is supposed to be calculated
   ! *   or set by this subroutine and passed back to the calling program.
   ! *
   ! *        CURV   - Base curvature of surface (input);  this is
   ! *                 the value of curvature that is entered by you with
   ! *                 your surface data, i.e. CUY for this surface.
   ! *        DATA   - Array containing the special user-defined-surface
   ! *                 parameters (input);  for example, DATA (1)
   ! *                 is the value entered with the command UCO C1 for
   ! *                 the UDS surface,  DATA (2) is entered with the
   ! *                 command UCO C2, etc.
   ! *        X,Y,Z  - Coordinates of the current position of the
   ! *                 ray with respect to the origin of the surface
   ! *                 (input)
   ! *        F      - Value of the function F (output)
   ! *        FX, FY, FZ  - Values of the derivatives of the function F
   ! *                 (output if MODE = 1, otherwise not used)
   ! *        ERROR  - Error code; nonzero if error in surface calculation.
   ! *                 Set to -1 if the limiting value of the real part of
   ! *                 F is negative, else set to 1.  Code V iterates to a
   ! *                 UDS surface by finding the point along the ray where
   ! *                 F changes sign.  When F fails to evaluate (e.g.,
   ! *                 because of division by zero or square root of a
   ! *                 negative number) the program can still converge if
   ! *                 the sign of F can be determined.  (output)
   ! *        MODE   - Flag to the calling program to indicate whether
   ! *                 derivatives are exact functions or finite
   ! *                 differences (output)
   ! *                 MODE = 0 :  CODE V computes derivatives
   ! *                             via finite differences.
   ! *                 MODE = 1 :  CODE V expects this subroutine
   ! *                             to compute the derivatives (user
   ! *                             must program the exact derivative
   ! *                             functions).
   ! *
   ! ******************************************************************

   USE qmnp

   IMPLICIT NONE

   ! INPUT
   ! Note: The following statements defining the type of the variables
   ! in the call list must be left in the subroutine
   INTEGER(kind=4), INTENT(in) :: coefcount
   REAL(kind=8), INTENT(in) :: X, Y, Z, CURV, COEF(coefcount)
   REAL(kind=8), INTENT(inout):: f, fx, fy, fz
   INTEGER(kind=4), INTENT(inout):: error, mode

   !LOCAL
   REAL(KIND=8) :: r2!, t!, psi, arg
   ! REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: amn, bmn
   REAL(KIND=8) :: output(3)

   ! i.e. a formula for the derivatives is programmed.
   mode = 1

   IF ((coef(1) < 1E-13) .OR. (coef(2) == 0) .OR. &
       (ABS(coef(2)/curv) < coef(1)) .OR. (coefcount < 2)) THEN
      error = 1
      RETURN
   END IF
   r2 = x**2 + y**2
   IF ((CURV**2*r2) > 1.0D+00) THEN
      error = INT(SIGN(1.0d+00,(CURV * r2) - z ))
      RETURN 
   END IF

   CALL qmnp_init_from_coef(coefcount - 2)

   ! ALLOCATE (amn(m_max + 1, n_max + 1))
   ! ALLOCATE (bmn(m_max + 1, n_max + 1))

   !*********************************************************************
   !**********                                                ***********
   !**********           COMPUTE THE SAG                      ***********
   !**********                                                ***********
   !*********************************************************************

   ! t = ATAN2(y, x)

   ! amn = 0.D0
   ! bmn = 0.D0

   ! CALL data2arrfast(COEF(3:), &
   !                   amn, bmn)

   ! CALL build_map_int(r, t, curv/coef(2), &
   !                    0.D0, COEF(1), &
   !                    amn, bmn, output)

   ! output = output*coef(2)
   ! F = output(1) - z

   ! FX = COS(t)*output(2)
   ! FY = SIN(t)*output(2)
   ! IF (r > 1E-15) THEN
   !    FX = (FX - (1/r)*SIN(t)*output(3))
   !    FY = (FY + (1/r)*COS(t)*output(3))
   ! END IF
   ! FZ = -1.D0

   ! write (*,*) x, y
   ! write (*,*) "--------"
   CALL build_map_int2(x, y, curv/coef(2), &
                       0.D0, COEF(1), &
                       coef(3:), output)
   output=output*coef(2)
   F = output(1) - z
   FX = output(2)
   FY = output(3)
   FZ = -1.D0
END SUBROUTINE

!DEC$ ATTRIBUTES DLLEXPORT :: cvisthreadsafe
LOGICAL FUNCTION CVISTHREADSAFE() BIND(c)
   CVISTHREADSAFE = .TRUE.
END

!DEC$ ATTRIBUTES DLLEXPORT :: cvumrtype
SUBROUTINE cvumrtype(umrType)
   INTEGER umrType
   CHARACTER*4 ThisType
   ThisType = 'UD1 '
   READ (ThisType, '(1a4)') umrType
END

!DEC$ ATTRIBUTES DLLEXPORT :: cvgetcoefname
SUBROUTINE cvgetcoefname(nCoef, nCoefName) BIND(c)
   USE qmn_array_index
   INTEGER(KIND=4) :: nCoef
   INTEGER(KIND=4), DIMENSION(8) :: nCoefName
   INTEGER(KIND=4), DIMENSION(2) :: table_index
   INTEGER(KIND=4) :: table_indexn
   CHARACTER(LEN=1) :: table_name
   INTEGER(KIND=4) :: table_namen
   CHARACTER(LEN=32) :: sName
   CHARACTER(LEN=32), DIMENSION(2) :: sNames
   sNames = [ &
            'Normalisation radius            ', &
            'Multiplication factor           ' &
            ]

   IF (nCoef .GE. 1 .AND. nCoef .LE. 2) THEN
      sName = sNames(nCoef)
   ELSE
      CALL num2qmnindex(nCoef - 2, table_namen, table_index)
      IF (table_namen) THEN
         table_name = "b"
      ELSE
         table_name = "a"
      END IF
      WRITE (sName, '(a1,a2,i3,a3,i3,a1)') table_name, "_(", table_index(1) - 1, ")_(", table_index(2) - 1, ")"
   END IF

   READ (sName, '(8A4)') nCoefName
END SUBROUTINE

! ! !
! !DEC$ ATTRIBUTES DLLEXPORT :: test_zero_coefs
! SUBROUTINE TEST_ZERO_COEFS
!         real ( kind = 8 ) :: X, Y, Z, CURV, COEF(15)
!         real ( kind = 8 ) :: f, fx, fy, fz
!         integer ( kind = 4 ) :: error, mode
!         character(100) :: message
!
!         ! real (kind = 8) :: x, y, z, curv
!         ! real (kind = 8) :: f, fx, fy, fz
!         ! integer (kind = 4) :: mode, error
!         ! real (kind = 8) :: data
!         ! DOUBLE PRECISION data(81)
!         ! CHARACTER*100 message
!
!         x = 0.0d0
!         y = 0.0d0
!         z = 0.0d0
!         curv = 0.0
!         error = 0
!         mode = 0
!         data = 0
!         coef = 0
!
!         CALL uds_qmn(x, y, z, curv, f, fx, fy, fz, error, mode, coef, size(coef))
!         Write(message,100) error, mode
!         Call CVPUTREC(message)
!         100 Format(1x,'USERSUR test_zero_coefs, ERROR=',i1,' MODE=',i1)
!         Write(message,110) f, fx, fy, fz
!         Call CVPUTREC(message)
!         110 Format(1x,'USERSUR test_zero_coefs, f,fx,fy,fz',4(g15.5,x))
! END subroutine
