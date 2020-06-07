module IMSL_IMPL


contains

    subroutine RNSET(ivalue)
        INTEGER :: ivalue
        
        call SRAND(ivalue)
    end subroutine

    subroutine DRNUN(vec_size, random_vector)
        DOUBLE PRECISION :: random_vector(:)
        INTEGER :: vec_size,i

        do i=1, vec_size
            random_vector(i) = rand()
        end do

    end subroutine 

end module
