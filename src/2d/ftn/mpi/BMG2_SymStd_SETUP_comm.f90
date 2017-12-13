      subroutine BMG2_SymStd_SETUP_comm(COMM, NOL, MyProc, &
           &                                  Min_GSz)&
           BIND(C, NAME='MPI_BMG2_SymStd_SETUP_comm')

      USE ModInterface
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

      call MPI_COMM_SIZE(COMM(1,NOL), INSIZE, IERR)

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
     &                        TempProc, ACTIVE)

            INSIZE = HSIZE
            N = N - 1

            goto 10

         ELSE

            goto 20

         END IF

20    CONTINUE

      END IF

      !
      ! Set COMM(2,1) equal to COMM(1,1).  This facilitates
      ! the loops used in the up and down portions of the
      ! linesolve.
      !

      COMM(2,2) = COMM(1,2)

      END

! ==============================================================
! ==============================================================
! ==============================================================

      SUBROUTINE comm_coarsen(INCOMM, INSIZE, &
     &                        GCOMM, HCOMM,&
     &                        NGPS, Min_GSz, &
     &                        MyProc, ACTIVE)

      IMPLICIT NONE

      INCLUDE 'mpif.h'

      INTEGER INCOMM, INSIZE
      INTEGER GCOMM, HCOMM
      INTEGER NGPS, Min_GSz
      INTEGER ACTIVE

      INTEGER GSIZE, HSIZE

      INTEGER MPI_MY_ID, MyProc
      INTEGER HFLAG, HRANK
      INTEGER MyGroup, MyGroupRank
      INTEGER LBRK, RBRK, REM, CLUMP

      INTEGER EXITFLAG, K, IERR


      !
      ! Compute the number of new groups to form
      !
      NGPS = INSIZE / Min_GSz
      REM = mod(INSIZE, Min_GSz)

      !
      ! This computers CLUMP = ceiling(REM/NGPS)
      !
      IF (real(REM)/NGPS .GT. REM/NGPS) THEN
         CLUMP = REM/NGPS + 1
      ELSe
         CLUMP = REM/NGPS
      END IF

      !
      ! Initialize counters and flags
      !
      LBRK = 0
      RBRK = 0
      EXITFLAG = 0
      K = 1

      !
      ! Determien the size of each group taking
      ! load balancing issues into account
      !
30    CONTINUE

      !
      ! This is just a while-loop
      !
      IF (EXITFLAG .NE. 1) THEN

         IF (REM .GE. CLUMP) THEN
            GSIZE = Min_GSz + CLUMP
            REM = REM - CLUMP
         ELSE
            GSIZE = Min_GSz + REM
            REM = 0
         END IF

         RBRK = RBRK + GSIZE

         !
         ! if local process fits in this
         ! group set group and local
         ! rank flags
         !
         IF (MyProc .LE. RBRK) THEN
            MyGroup = K
            MyGroupRank = MyProc - LBRK - 1
            EXITFLAG = 1
         ELSE
            K = K + 1
         END IF

         LBRK = RBRK

         goto 30

      ELSE

         goto 40

      END IF

40    CONTINUE

      !
      ! Form group communicators
      !
      CALL MPI_COMM_SPLIT(INCOMM, MyGroup, MyGroupRank,&
     &                    GCOMM, IERR)

      !
      ! If I'm the first process in my group designate
      ! me the HEAD of my group.  Also check to see if
      ! I will be part of any other communicators at
      ! lower levels.  If I won't be then set ACTIVE
      ! to 0
      !
      IF (MyGroupRank .EQ. 0) THEN

         HFLAG = 1

         IF (NGPS .GE. 2*Min_GSz) THEN
            ACTIVE = 1
         ELSE
            ACTIVE = 0
         END IF

      ELSE

         HFLAG = MPI_UNDEFINED
         ACTIVE = 0

      END IF

      !
      ! Identify the rank of this process relative
      ! to its new subgroup
      !
      HRANK = MyGroup - 1

      !
      ! Form communicator with all of the HEAD nodes
      ! in it.  Non-head nodes get assigned to a
      ! communicator with type MPI_PROC_NULL
      !
      CALL MPI_COMM_SPLIT(INCOMM, HFLAG, HRANK, &
     &                    HCOMM, IERR)

      !
      ! Update the rank that this process will have
      ! in the new head communicator
      !
      MyProc = MyGroup

      END
