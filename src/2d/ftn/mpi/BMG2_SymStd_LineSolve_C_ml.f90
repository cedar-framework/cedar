      SUBROUTINE BMG2_SymStd_LineSolve_C_ml( &
     &                D, OD, F,   &
     &                INBUFF, NPts, P,&
     &                MY_ID&
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Solve the local tridiagonal systems.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     P          - Number of processes in row communicator
!
!     Npts       - Number of local unknowns
!
!     INBUFF     - Buffer to read in interface solutions
!
!     DATADIST   - DATADIST(1,I) gives the global index of
!                  the first unknown local to process I.
!                  DATADIST(2,I) gives the global index of
!                  the last unkown local to process I.
!
!     AHAT       - Workspace for linear solve
!     FHAT       - Workspace for linear solve
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================

      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER MY_ID, NPts, P
      REAL*8 D(NPts), OD(NPts), F(0:NPts+1)
      REAL*8 INBUFF(*)

! ----------------------------
!     Local Declarations
!
      INTEGER I, AT, FT, IERR, INFO, rank
      REAL*8  B(NPts-2)

! ======================================================================

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,IERR)

!     write(*,*) 'NPts = ', NPts, ' on proces ', rank

!     IF (rank .EQ. 0) THEN
!        DO I=1,NPts
!           write(*,*) OD(I), D(I)
!        END DO
!     END IF

!     CALL MPI_FINALIZE(IERR)
!     STOP

      !
      ! initialize the solution with the
      ! solution to the interface equations
      !
      IF (My_ID .GT. 0) THEN
         F(0) = INBUFF(1)
      ELSE
         F(0) = 0d0
      ENDIF

      F(1) = INBUFF(2)
      F(NPts) = INBUFF(3)

      IF (MY_ID .LT. P-1) THEN
         F(NPts+1) = INBUFF(4)
      ELSE
         F(NPts+1) = 0d0
      ENDIF

      F(2) = F(2) - OD(2)*F(1)
      F(NPts-1) = F(NPts-1) - OD(NPts)*F(NPts)

      CALL DPTTRS(NPts-2, 1, D(2), OD(3), F(2), NPts-2, INFO)

      IF (INFO .NE. 0) THEN
         write(*,*) 'Uh oh...'
      END IF

!MB      IF (rank .EQ. 0) THEN
!MB         write(*,*) ''

!MB         DO I=1,NPts
!MB            write(*,*) F(I)
!MB         END DO

!        write(*,*) INBUFF(1)
!        write(*,*) INBUFF(2)
!        write(*,*) INBUFF(3)
!        write(*,*) INBUFF(4)

!        write(*,*) F(0)
!        write(*,*) F(1)
!        write(*,*) F(NPts)
!        write(*,*) F(NPts+1)

!MB      END IF


! ======================================================================

      RETURN
      END
