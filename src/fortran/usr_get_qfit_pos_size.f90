!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_qfit_pos_size
SUBROUTINE usr_get_qfit_pos_size (coefcount, sizemap, sizeann) BIND(C)
        use qfit
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: coefcount
        REAL(KIND=8), INTENT(OUT) :: sizemap, sizeann
        call qmnp_init_from_coef(INT(coefcount,KIND=4)-2)
        sizeann = REAL(bfs_size, KIND=8)
        sizemap = REAL((m_max+1)*2*(n_max+2), KIND=8)
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
