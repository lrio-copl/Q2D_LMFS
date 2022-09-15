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

!DEC$ ATTRIBUTES DLLEXPORT :: uds_qmn
SUBROUTINE uds_qmn (x, y, z, curv, f, fx, fy, fz, error, mode, coef, coefCount) BIND(C)

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


        use qmnp
        IMPLICIT NONE



        ! Note: The following statements defining the type of the variables
        ! in the call list must be left in the subroutine
        INTEGER(KIND=4), INTENT(IN) :: coefcount
        REAL(KIND=8), INTENT(IN) :: X, Y, Z, CURV, COEF(coefcount)
        REAL(KIND=8), INTENT(INOUT):: f, fx, fy, fz
        INTEGER(KIND=4), INTENT(INOUT):: error, mode
        REAL(KIND=8), DIMENSION(4, MAX(6,COEFCOUNT-2)+3) :: amn, bmn
        REAL(KIND=8) :: r, t
        REAL(KIND=8) :: finish
        INTEGER(KIND=4) :: coefreal
        REAL(KIND=8) :: coefcount_compute
        REAL(KIND=8) :: output(3)
        REAL(KIND=8) :: output2(2,2,2)
        INTEGER(KIND=4) :: i

        error = 0

        IF ( coefcount < 2 ) THEN
                error = 1
                CALL cverror("Coefficient must be larger than 2")
        ! ELSE IF (coefcount<8) then
        !     coefreal = 6
        ELSE 
            ! coefreal = coefcount-2
            coefreal = MAX(6,coefcount-2)
        END IF

        IF (( .NOT. ( ( m_max==3 ) .AND. ( n_max==coefreal-1 ) ) ) .AND. ( error==0 )) THEN
                call qmnp_init(3,coefreal-1)
        END IF

        !  ****************************************************************

        !     Note: From this point on, you will typically substitute your own
        !           FORTRAN code for the particular surface being programmed
        !  
        !                            Entries in DATA array:

        !                           DATA( )      COEFFICIENT OF
        !                          ------------------------------------------------
        !                                1      Normalization radius 
        !                                2      Conic
        !                                3      Qm #2
        !                                4      Qm #3
        !                                5      Qm #4
        !                                6      Qm #5
        !                                7      Qm #6
        !                                8      Qm #7
        !                                9      Qm #8
        !                                10     Qm #9
        !                                11     Qm #10
        !                                12     Qm #11
        !                                13     Qm #12
        !                                14     Qm #13
        !                                15     Qm #14
        !                                16     Qm #15
        !                                17     Qm #16
        !                                18     Qm #17
        !                                19     Qm #18
        !                                20     Qm #19
        !                                21     Qm #20
        !                                22     Qm #21
        !                                23     Qm #22
        !                                24     Qm #23
        !                                25     Qm #24
        !                                26     Qm #25
        !                                27     Qm #26
        !                                28     Qm #27
        !                                29     Qm #28
        !     ****************************************************************

        ! i.e. a formula for the derivatives is programmed.
        mode = 1

        !*********************************************************************
        !**********                                                ***********
        !**********           COMPUTE THE SAG                      ***********
        !**********                                                ***********
        !*********************************************************************
        r = sqrt(x**2+y**2)
        t = atan2(y,x)

        IF (1/curv < coef(1)) THEN
            error = 1
        END IF

        IF ( ( coef(1) < 1e-13 ) .OR. ( error==1 ) ) THEN
                F = 0
                Fx = 0
                Fy = 0
                Fz = -1
                error = 1
        ELSE
                F = 0
                amn = 0
                bmn = 0
                amn(1,1:coefcount-2) = COEF(3:coefcount)
                IF ( r > COEF(1) ) THEN
                        call build_map_int(1.d0, t, curv, COEF(2), COEF(1), amn, bmn, output)
                ELSE
                        call build_map_int(r, t, curv, COEF(2), COEF(1), amn, bmn, output)
                END IF
                F = output(1) - z
                IF ( r < 1e-15 ) THEN
                    ! IF ( SIZE(drqmn_sum) .NE. 2) THEN
                    CALL build_map([0.0d+00,0.0d+00], [0.0d+00,2.0d+00*Atan(1.0d+00)], curv, COEF(2), COEF(1), amn,bmn,output2)
                    ! END IF
                    Fx = output2(2,1,1)
                    Fy = output2(2,1,2)
                ELSE
                    Fx = cos(t)*output(2) + (1/r)*cos(t)*output(3)
                    fy = sin(t)*output(2) + (1/r)*cos(t)*output(3)
                END IF
                Fz = -1
                error = 0
        END IF
END SUBROUTINE

!DEC$ ATTRIBUTES DLLEXPORT :: cvisthreadsafe
logical function CVISTHREADSAFE() bind(c)
        CVISTHREADSAFE = .TRUE.
END 


!DEC$ ATTRIBUTES DLLEXPORT :: cvumrtype
subroutine cvumrtype(umrType)
    INTEGER umrType
    Character*4 ThisType
    ThisType = 'UD1 '
    READ(ThisType, '(1a4)') umrType
END


!DEC$ ATTRIBUTES DLLEXPORT :: cvgetcoefname
SUBROUTINE cvgetcoefname(nCoef, nCoefName) bind(c)
    INTEGER(kind=4) :: nCoef
    INTEGER(kind=4) :: nCoefName(8)
    CHARACTER(len=32) :: sName
    CHARACTER(len=32) :: sNames(2)
    sNames = [ &
        'Normalisation radius            ',&
        'Multiplication factor           '&
        ]

    IF ( nCoef .GE. 1 .AND. nCoef .LE. 2 ) THEN
        sName = sNames(nCoef)
    ELSE 
        write(sName,'(a7,i3,a1)') "a_(0)_(", nCoef-2, ")"
    ENDIF

    READ(sName,'(8A4)') nCoefName
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
