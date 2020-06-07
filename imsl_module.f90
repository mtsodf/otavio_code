module IMSL_IMPLi


contains

    subroutine RNSET(seed)
        integer seed

        CALL RANDOM_SEED(PUT = seed)
    end subroutine

end module
