      SUBROUTINE BMG2_SymStd_LineSolve_C_ml_eff( &
     &                A, B, C, F, RBUFF, RWORK, &
     &                NP, P, MY_ID&
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
      REAL*8 A(NP), B(NP), C(NP), F(0:NP+1)
      REAL*8 RBUFF(*), RWORK(*), L

! ----------------------------
!     Local Declarations
!
      INTEGER I, AT, FT

! ======================================================================

      ! generate pointers into RWORK
      AT = 1
      FT = AT + NP + 2

      ! initialize the solution with the
      ! solution to the interface equations
      IF (My_ID .GT. 0) THEN
         F(0)    = RBUFF(1)
      ENDIF

      F(1)    = RBUFF(2)
      F(NP)   = RBUFF(3)

      IF (MY_ID .LT. P-1) THEN
         F(NP+1) = RBUFF(4)
      ENDIF

      ! solve the on processor tridiagonal system by
      ! forward and backsubstitution without modifying
      ! A,B,C, and F this requires some workspace (RWORK)
      RWORK(FT+2) = F(2)-C(2)*F(1)
      RWORK(AT+2) = A(2)
      DO I=3,NP-2
         L = -C(I)/RWORK(AT+I-1)
         RWORK(AT+I) = A(I) + L*B(I-1)
         RWORK(FT+I) = F(I) + L*RWORK(FT+I-1)
      ENDDO
      L = -C(NP-1)/RWORK(AT+NP-2)
      RWORK(AT+NP-1) = A(NP-1) + L*B(NP-2)
      RWORK(FT+NP-1) = F(NP-1) - B(NP-1)*F(NP) + L*RWORK(FT+NP-2)

      F(NP-1) = RWORK(FT+NP-1)/RWORK(AT+NP-1)
      DO I=NP-2,2,-1
         F(I) = (RWORK(FT+I)-B(I)*F(I+1))/RWORK(AT+I)
      ENDDO

! ======================================================================

      RETURN
      END
