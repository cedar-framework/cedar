      subroutine BMG2_SymStd_SETUP_comm(COMM, NOL, MyProc, &
           & Min_GSz, mp) &
           & BIND(C, NAME='MPI_BMG2_SymStd_SETUP_comm')

      USE ModInterface
      use message_passing
      IMPLICIT NONE

! ==============================================================
! Include Statements
! ==============================================================

      INCLUDE 'mpif.h'

      INCLUDE 'BMG_parameters_f90.h'

! ==============================================================
! Input Variables
! ==============================================================

      integer(c_int), VALUE :: NOL, MyProc, Min_GSz

      integer :: COMM(2,2:NOL)
      type(c_ptr) :: mp


! ==============================================================
! Local Variables
! ==============================================================

      INTEGER TempProc

      INTEGER INSIZE, GSIZE, HSIZE

      INTEGER ACTIVE, N

      INTEGER IERR

! ==============================================================
! Initialize local variables
! ==============================================================

      call MPI_Comm_size(COMM(1,NOL), INSIZE, IERR)

      N = NOL
      TempProc = MyProc
      ACTIVE = 1

! ==============================================================
! Call comm_coarsen to build communicators
! ==============================================================

      IF (INSIZE .GE. 2*Min_GSz) THEN

10       CONTINUE

         !
         ! If there are any process left in the current
         ! level communicator try to coasren.
         ! This is just a while-loop
         !
         IF (ACTIVE .EQ. 1) THEN

            call comm_coarsen(COMM(1,N), INSIZE, &
     &                        COMM(2,N), COMM(1,N-1),&
     &                        HSIZE, Min_GSz,&
     &                        TempProc, ACTIVE, mp)

            INSIZE = HSIZE
            N = N - 1

            goto 10

         ELSE

            goto 20

         END IF

      END IF

      !
      ! Set COMM(2,1) equal to COMM(1,1).  This facilitates
      ! the loops used in the up and down portions of the
      ! linesolve.
      !

20    CONTINUE
      COMM(2,2) = COMM(1,2)

    END subroutine BMG2_SymStd_SETUP_comm
