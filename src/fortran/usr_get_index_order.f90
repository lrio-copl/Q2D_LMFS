!------------------------------------------------------------------------
!   DLLMAIN
!   USEOBJECT:cvputrec
!------------------------------------------------------------------------

!DEC$ ATTRIBUTES DLLEXPORT :: usr_get_index_order
SUBROUTINE usr_get_index_order (order_in, index_out) BIND(C)
        use qmn_array_index
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: order_in
        REAL(KIND=8), INTENT(INOUT) :: index_out
        INTEGER(KIND=4) :: index_local
        INTEGER(KIND=4) :: table_name
        INTEGER(KIND=4), DIMENSION(2) :: table_index
        REAL(KIND=8) :: order_loc

        order_loc = 0
        index_local = 0
        do while ( order_loc < order_in )
            index_local = index_local + 1
            call num2qmnindex(index_local, table_name, table_index)
            if (table_index(1) == 1) then
                order_loc = order_loc + 1
            end if
        end do
        index_out = REAL(index_local, KIND=8)
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
