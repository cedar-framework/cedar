      SUBROUTINE BMG2_SymStd_LineSolve_B_ml( &
     &                INBUFF, OUTBUFF, &
     &                AHAT, BHAT, CHAT, FHAT, XHAT, &
     &                Npts, P, DATADIST, MY_ID, INCR, INCR2&
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Solve the interface system for current line.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     Npts - Number of local unknowns
!
!     INBUFF - Buffer to read in interface equation data
!
!     OUTBUFF - Buffer to return interface system solution
!
!     DATADIST - DATADIST(1,I) gives the global index of
!                the first unknown local to process I.
!              - DATADIST(2,I) gives the global index of
!                the last unkown local to process I.
!
!     (A,B,C,F,X)HAT - Workspace for Thomas Algorithm
!
!     AHAT     - Main Diagonal Entries
!     CHAT     - Lower Off Diagonal Entries
!     BHAT     - Upper Off Diagonal Entries
!     FHAT     - Right Hand Side Vector
!     XHAT     - Array to hold solution vector
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
      include 'mpif.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER MY_ID, Npts, P
      INTEGER DATADIST( 1:2 , 0:(P-1) )

      REAL*8  INBUFF(*), OUTBUFF(*), L
      REAL*8  AHAT(2*P+2),  BHAT(2*P+2),  CHAT(2*P+2)
      REAL*8  FHAT(2*P+2),  XHAT(2*P+2)

! ----------------------------
!     Local Declarations
!
      INTEGER I, OFFSET, INCR, INCR2, N, IERR

      INTEGER grank

! ======================================================================

      ! RWORK is workspace with nothing in it on entry
      ! P is the number of processors

      IF (MY_ID .EQ. 0) THEN

! ====================================================================
!     Assemble the tridiagonal system of interface equations
! ====================================================================

         OFFSET = 0
         N = 0
         DO I = 0, P-1

           !
           ! there is always at least one interface equation
           !
           N = N + 1

           CHAT(N) = INBUFF(OFFSET+1)
           AHAT(N) = INBUFF(OFFSET+2)
           BHAT(N) = INBUFF(OFFSET+3)
           FHAT(N) = INBUFF(OFFSET+4)

           !
           ! If there is more than one local unknown on this
           ! process then there are two interface equations
           !
           IF (NPts .GT. 1) THEN

              N = N + 1

              CHAT(N) = INBUFF(OFFSET+5)
              AHAT(N) = INBUFF(OFFSET+6)
              BHAT(N) = INBUFF(OFFSET+7)
              FHAT(N) = INBUFF(OFFSET+8)


           END IF

           OFFSET = OFFSET + INCR
        END DO

        CALL MPI_COMM_RANK(MPI_COMM_WORLD,grank,IERR)

!       IF (grank .EQ. 0) THEN
!          CALL BMG2_SymStd_IF_RES(AHAT,BHAT,CHAT,FHAT,1,1,N,10,
!    &               grank, 0)
!       END IF

!     CALL MPI_FINALIZE(IERR)
!     STOP

!       write(*,*) ''
!       DO I=1,N
!          write(*,*) CHAT(I), AHAT(I), BHAT(I) !, FHAT(I)
!          write(*,*) FHAT(I)
!       END DO

!     CALL MPI_FINALIZE(IERR)
!     STOP


!====================================================================
!      Solve the tridiagonal system using the Thomas Algorithm
!====================================================================

        !
        ! Forward Sweep
        !
        DO I = 2, N

           L = CHAT(I) / AHAT(I-1)
           AHAT(I) = AHAT(I) - BHAT(I-1) * L
           FHAT(I) = FHAT(I) - FHAT(I-1) * L

        END DO

        XHAT(N) = FHAT(N) / AHAT(N)

        !
        ! Backward Sweep
        !
        DO I = N-1, 1, -1
           XHAT(I) = (FHAT(I)-BHAT(I)*XHAT(I+1))/AHAT(I)
        END DO

!====================================================================
!      Copy the solution into OUTBUFF
!====================================================================

        !
        ! Start by copying the solution values corresponding
        ! to local ghost points
        !
        OFFSET = 0
        N = 1

        DO I = 0, P-1

           IF (I .GT. 0) OUTBUFF(OFFSET+1) = XHAT(N-1)

           !
           ! Check
!          IF (I .EQ. 0) OUTBUFF(OFFSET+1) = 0

           IF (NPts .GT. 1) THEN
              N = N + 2
           ELSE
              N = N + 1
           END IF

           IF (I .LT. P-1) OUTBUFF(OFFSET+4) = XHAT(N)

           !
           ! Check
!          IF (I .EQ. P-1) OUTBUFF(OFFSET+4) = 0

           OFFSET = OFFSET + INCR2

        END DO

        !
        ! Then copy solution values corresponding to
        ! local interior points
        !
        OFFSET = 0
        N = 1

        DO I = 0, P-1

           IF (NPts .GT. 1) THEN
              OUTBUFF(OFFSET+2) = XHAT(N)
              OUTBUFF(OFFSET+3) = XHAT(N+1)
              N = N + 2
           ELSE
              OUTBUFF(OFFSET+2) = XHAT(N)
              N = N + 1
           END IF

           OFFSET = OFFSET + INCR2

        END DO

!       write(*,*) ''
!       DO I=0,P-1
!          write(*,*) OUTBUFF(OFFSET+1), OUTBUFF(OFFSET+2),
!    &                OUTBUFF(OFFSET+3), OUTBUFF(OFFSET+4)
!          OFFSET = OFFSET + INCR2
!       END DO

!       write(*,*) OUTBUFF(1), OUTBUFF(2), OUTBUFF(3), OUTBUFF(4)

!       write(*,*) ''
!       DO I=1,N-1
!          write(*,*) XHAT(I)
!       END DO

      END IF

!     CALL MPI_FINALIZE(IERR)
!     STOP


! ======================================================================

      RETURN
      END
