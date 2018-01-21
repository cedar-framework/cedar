      subroutine BMG2_SymStd_SETUP_trisolve(&
     &                NProc, NLines, Min_GSz,&
     &                NOG, NSPACE, NOL)&
     BIND(C, NAME='MPI_BMG2_SymStd_SETUP_trisolve')

      USE ModInterface
      IMPLICIT NONE

! =================================================
! Intent IN Variables
! =================================================

      integer(c_int), value :: NProc, NOG, Min_GSz
      integer(c_int), value :: NLines

! =================================================
! Intent OUT Variables
! =================================================

      integer(c_int) :: NSPACE, NOL

! =================================================
! Local Variables
! =================================================

      INTEGER NGPS, REM, CLUMP, GSIZE, INSIZE
      INTEGER NL, kg

! ==================================================
! Initialize Parameters
! ==================================================

      NL = NLines
      NSPACE = 0

      if (NOG .eq. 1) then
         NOG = NOG + 1
      endif

      DO kg = 2, NOG

         NL = 1 + ( NLines - 1 ) / 2**(kg-2)

         INSIZE = NProc
         NOL = 1

10       CONTINUE

         !
         ! If there are not enough processes in
         ! the current group to form at least
         ! two groups of Min_GSz then after this
         ! level we gather to process 0
         !

         IF (INSIZE .LT. 2*Min_GSz) THEN

            NOL = NOL + 1

            goto 20

         !
         !  If there are enough processes to form
         !  at least two groups of Min_GSz then
         !  we coarsen.
         !

         ELSE IF (INSIZE .GE. 2*Min_GSz) THEN

            !
            ! The number of groups to
            ! divide processes into
            !

            NGPS = INSIZE / Min_GSz

            !
            ! The remaining number of processes
            ! after diving line into groups of
            ! of size Min_GSz
            !
            REM = mod(INSIZE, Min_GSz)

            !
            ! This just computes the ceiling of
            ! REM/NGPS.  CLUMP is the increment
            ! in which the remaining processes
            ! will be distributed
            !
            IF (real(REM)/NGPS .GT. REM/NGPS) THEN
               CLUMP = REM/NGPS + 1
            ELSE
               CLUMP = REM/NGPS
            END IF

            !
            ! Compute the size of the group
            !
            GSIZE = Min_GSz + CLUMP

            !
            ! Reset INSIZE
            !
            INSIZE = NGPS

            !
            ! Increment NSPACE
            !
            NSPACE = NSPACE + (2*GSIZE+2)*NL*4

            !
            !  Increment Number of Levels
            !
            NOL = NOL + 1

            goto 10

         END IF

20       CONTINUE

      END DO

      END
