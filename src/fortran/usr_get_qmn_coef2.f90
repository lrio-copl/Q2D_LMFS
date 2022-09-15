!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qmn_coef2
SUBROUTINE usr_get_qmn_coef2 (coefcount, coef, bfs, radius, sizemap, sizeann, map, ann) BIND(C)
        use qfit
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: coefcount, sizemap, sizeann
        REAL(KIND=8), INTENT(IN) :: radius
        REAL(KIND=8), DIMENSION(INT(sizemap,KIND=4)), INTENT(IN) :: map
        REAL(KIND=8), DIMENSION(INT(sizeann,KIND=4)), INTENT(IN) :: ann
        REAL(KIND=8), INTENT(OUT) :: bfs
        REAL(KIND=8), DIMENSION(INT(coefcount,KIND=4)), INTENT(OUT) :: coef
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: amn, bmn
        REAL(KIND=8) :: maxamp

        call qfit_init_from_coef(INT(coefcount,KIND=4)-2)

        ALLOCATE(amn(m_max+1,n_max+1))
        ALLOCATE(bmn(m_max+1,n_max+1))

        maxamp = MAXVAL(map)

        call q_fit_t(map/maxamp, ann/maxamp, radius, amn, bmn, bfs)
        call arr2data(coef(3:), amn, bmn)

        coef(1) = radius
        coef(2) = maxamp
        bfs = bfs*maxamp
        ! write (*,*) map
        ! write (*,*) ann
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
