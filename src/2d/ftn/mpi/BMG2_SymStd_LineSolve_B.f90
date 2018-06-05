      SUBROUTINE BMG2_SymStd_LineSolve_B( &
     &                A, B, C, F,  RWORK, RBUFF1, RBUFF2, &
     &                NP, P, DATADIST, MY_ID, INCR, INCR2&
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Solve the interface system for each line.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     P        - number of processes
!     NP       - number of unknowns (or equations) on this process
!     K        - number of tridiagonal systems
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

! ----------------------------
!     Argument Declarations
!
      INTEGER MY_ID, NP, P
      REAL*8 A(NP), B(NP), C(NP), F(NP)
      INTEGER DATADIST( 1:2 , 0:(P-1) )
      REAL*8 RWORK(*), RBUFF1(*), RBUFF2(*), L

! ----------------------------
!     Local Declarations
!
      INTEGER I, OFFSET, INCR, INCR2
      INTEGER AT, BT, CT, FT, XT, N

! ======================================================================

      ! RWORK is workspace with nothing in it on entry
      ! P is the number of processors

      IF (MY_ID.eq.0) THEN

         ! create pointers into the workspace RWORK
         AT = 1
         BT = AT + 2*P +2
         CT = BT + 2*P +2
         FT = CT + 2*P +2
         XT = FT + 2*P +2

         ! assemble the tridiagonal system of
         ! interface equations
         OFFSET = 0
         N = 0
         DO I = 0, P-1

           ! there is always at least one interface equation
           RWORK(CT+N) = RBUFF1(OFFSET+1)
           RWORK(AT+N) = RBUFF1(OFFSET+2)
           RWORK(BT+N) = RBUFF1(OFFSET+3)
           RWORK(FT+N) = RBUFF1(OFFSET+4)
           N = N + 1

           ! figure out if there is more than one interface
           ! euqation, if the end index of interior points
           ! in the direction of the line on processor I is
           ! greater than the start index, there are two interface
           ! equations.
           IF( DATADIST(2,I) .gt. DATADIST(1,I) ) THEN
              RWORK(CT+N) = RBUFF1(OFFSET+5)
              RWORK(AT+N) = RBUFF1(OFFSET+6)
              RWORK(BT+N) = RBUFF1(OFFSET+7)
              RWORK(FT+N) = RBUFF1(OFFSET+8)
              N = N + 1
           END IF

           OFFSET = OFFSET + INCR
        END DO

        ! solve the reduced system on one processor
        ! using the Thomas algorithm
        DO I = 1, N - 1

           L = RWORK(CT+I) / RWORK(AT+I-1)

           RWORK(AT+I) = RWORK(AT+I) -&
     &          RWORK(BT+I-1) * L
           RWORK(FT+I) = RWORK(FT+I) -&
     &          RWORK(FT+I-1) * L
        END DO
        RWORK(XT+N-1) = RWORK(FT+N-1) / RWORK(AT+N-1)
        DO I = N-2, 0, -1
           RWORK(XT+I) = ( RWORK(FT+I) -&
     &          RWORK(BT+I) * RWORK(XT+I+1) ) / RWORK(AT+I)
        END DO

        OFFSET = 0
        N=0

        ! copy the solution into RBUFF2, start with the
        ! solution of the interface equations adjacent
        ! to processor I
        DO I = 0, P-1

           IF (I .GT. 0)    RBUFF2(OFFSET+1) = RWORK(XT+N-1)

           ! account for the case where there is only one interface
           ! equation on processor I
           IF (DATADIST(2,I) .gt. DATADIST(1,I)) THEN
              N = N + 2
           ELSE
              N = N + 1
           END IF
           IF (I .LT. P-1)  RBUFF2(OFFSET+4) = RWORK(XT+N)

           OFFSET = OFFSET + INCR2
        ENDDO


        ! continue with the solution to the interface equation
        ! on processor I
        OFFSET = 0
        N=0
        DO I = 0, P-1

           ! account for the case where there is only
           ! one interface equation on processor I
           IF (DATADIST(2,I) .gt. DATADIST(1,I)) THEN
              RBUFF2(OFFSET+2) = RWORK(XT+N)
              RBUFF2(OFFSET+3) = RWORK(XT+N+1)
              N = N + 2
           ELSE
              RBUFF2(OFFSET+2) = RWORK(XT+N+1)
              N = N + 1
           END IF

           OFFSET = OFFSET + INCR2
        ENDDO

      END IF

! ======================================================================

      RETURN
      END
