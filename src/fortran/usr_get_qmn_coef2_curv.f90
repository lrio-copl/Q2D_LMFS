
!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!   
!   USING A KNOW CURVATURE
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qmn_coef2_curv
SUBROUTINE usr_get_qmn_coef2_curv (coefcount, coef, bfs, radius, sizemap, map) BIND(C)
        use qfit
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: coefcount, sizemap
        REAL(KIND=8), INTENT(IN) :: radius
        REAL(KIND=8), DIMENSION(INT(sizemap,KIND=4)), INTENT(IN) :: map
        REAL(KIND=8), INTENT(IN) :: bfs
        REAL(KIND=8), DIMENSION(INT(coefcount,KIND=4)), INTENT(OUT) :: coef
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: amn, bmn
        REAL(KIND=8) :: maxamp

        call qfit_init_from_coef(INT(coefcount,KIND=4)-2)

        ALLOCATE(amn(m_max+1,n_max+1))
        ALLOCATE(bmn(m_max+1,n_max+1))

        maxamp = MAXVAL(map)

        call q_fit_curv_t(map/maxamp, radius, bfs/maxamp, amn, bmn)
        call arr2data(coef(3:), amn, bmn)

        coef(1) = radius
        coef(2) = maxamp
END SUBROUTINE

!DEC$ ATTRIBUTES DLLEXPORT :: cvisthreadsafe
logical FUNCTION CVISTHREADSAFE() BIND(c)
        CVISTHREADSAFE = .TRUE.
END 

!DEC$ ATTRIBUTES DLLEXPORT :: cvumrtype
subroutine cvumrtype(umrType)
    INTEGER umrType
    Character*4 ThisType
    ThisType = 'USR '
    READ(ThisType, '(1a4)') umrType
END
