      SUBROUTINE magma_ftn_init() BIND(C, NAME='magma_ftn_init')
        use magma
        IMPLICIT NONE

        call magmaf_init()
      end SUBROUTINE magma_ftn_init

      subroutine magma_ftn_finalize() bind(C, name='magma_ftn_finalize')
        use magma
        implicit none

        call magmaf_finalize()
      end subroutine magma_ftn_finalize
