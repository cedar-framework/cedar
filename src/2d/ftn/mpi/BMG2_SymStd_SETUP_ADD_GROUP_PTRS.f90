SUBROUTINE BMG2_SymStd_SETUP_ADD_GROUP_PTRS(&
     & nlines, comms, NOL, ptrs, mpool)&
     & BIND(C, NAME='BMG2_SymStd_SETUP_ADD_GROUP_PTRS')

         USE ModInterface
         use mempool
         IMPLICIT NONE
         include "mpif.h"
! ======================================================
! Intent In Variables
! ======================================================
         integer(c_int), value :: nlines, nol
         type(c_ptr) :: mpool

! ======================================================
! Intent Out Variables
! ======================================================

         integer(c_int) :: ptrs(2, 2:nol)
         integer :: comms(2, 2:nol)

! ======================================================
! Local Variables
! ======================================================

         integer :: kl, gsize, ierr
         integer :: nlcol(2)
         integer(c_int) :: pos(2)

! ======================================================
! Step through levels and record pointers
! ======================================================


         DO kl = NOL, 2, -1
            if (comms(2, kl) .ne. MPI_COMM_NULL) then
               call MPI_Comm_size(comms(2, kl), gsize, ierr)

               nlcol(1) = nlines / 2 + mod(nlines, 2)
               nlcol(2) = nlines / 2

               call cedar_mempool_pos(mpool, nlcol(1) * 8 * gsize, pos(1))
               call cedar_mempool_pos(mpool, nlcol(2) * 8 * gsize, pos(2))

               ptrs(1, kl) = pos(1)
               ptrs(2, kl) = pos(2)
            endif
         END DO

      RETURN
    END SUBROUTINE BMG2_SymStd_SETUP_ADD_GROUP_PTRS
