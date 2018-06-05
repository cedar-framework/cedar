      SUBROUTINE BMG2_SymStd_LineSolve_A_ml_eff( &
     &                D, OD, Q, RWORK, NPts&
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Form the interface system.
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
!     NLINES   - number of lines to be solved on this process row
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

! ---------------------------
!    Includes:
!

! ---------------------------
!    Argument Declarations:
!
      INTEGER  NPts
      REAL*8 D(Npts), OD(NPts+1), Q(NPts)
      REAL*8 RWORK(8)

! ----------------------------
!     Local Declarations
!
      REAL*8  L,A_TMP, F_TMP
      REAL*8  C_TMP, B_TMP
      INTEGER I, J

! ======================================================================

      ! RWORK(OFFSET+1) is the lower off diagonal element
      ! RWORK(OFFSET+2) is the diagonal element
      ! RWORK(OFFSET+3) is the upper off diagonal element
      ! RWORK(OFFSET+4) is the right hand side

      ! JBEG is the index of the first line that is solved
      !      it is computed one level up where the selection
      !      of red or black lines is detemined

      ! NP - THE NUMBER OF UNKNOWNS (OR EQUATIONS) per line ON THIS PROC
      ! NLINES - number of lines to be solved on this process row

      IF (NPts.eq.1) THEN
         ! special case: one local equation on this processor

         RWORK(1) = OD(1)
         RWORK(2) = D(1)
         RWORK(3) = OD(2)
         RWORK(4) = Q(1)

         ! Set flag to tell routine that there was only 1 point
         RWORK(6) = -1.0
      ELSE IF (NPts.eq.2) THEN
         ! special case: two local equations on this processor
         RWORK(1) = OD(1)
         RWORK(2) = D(1)
         RWORK(3) = OD(2)
         RWORK(4) = Q(1)

         RWORK(5) = OD(2)
         RWORK(6) = D(2)
         RWORK(7) = OD(3)
         RWORK(8) = Q(2)
      ELSE
         ! general case: more than two equation on this processor

         ! initialize with coefficients from the second equation
         ! in the line

         C_TMP = OD(2)
         F_TMP = Q(2)
         A_TMP = D(2)

         ! start iterating on the third equation in the line
         ! note that we reference i+1 in the body of the loop
         DO I=2,NPts-1,1
            ! compute the multiplier (lower off diagonal
            ! of euqation i+1 divided by the diagonal of the
            ! previous equation
            L =  -OD(I+1) / A_TMP

            ! compute the non-zero entry in column one of the
            ! matrix local to this processor
            C_TMP = L*C_TMP

            ! compute the right hand side
            F_TMP = Q(I+1) + L*F_TMP

            ! compute the diagonal element
            A_TMP = D(I+1) + L*OD(I+1)
         ENDDO

         ! populate RWORK with the lower interface equation
         RWORK(5) = C_TMP
         RWORK(6) = A_TMP
         RWORK(7) = OD(NPts+1)
         RWORK(8) = F_TMP


         ! initialize with coefficients from the second to last
         ! equation in the line
         A_TMP = D(NPts-1)
         B_TMP = OD(NPts)
         F_TMP = Q(NPts-1)

         ! start iterating on the third to last equation in the line
         ! note that we reference i+1 in the body of the loop
         DO I=NPts-3,0,-1
            ! compute the multiplier (opper off diagonal
            ! of euqation i+1 divided by the diagonal of the
            ! previous equation
            L =  -OD(I+2) / A_TMP

            ! compute the non-zero entry in column NP of the
            ! matrix local to this processor
            B_TMP = L * B_TMP

            ! compute the right hand side
            F_TMP = Q(I+1) + L*F_TMP

            ! compute the diagonal element
            A_TMP = D(I+1) + L*OD(I+2)
         ENDDO

         ! populate RWORK with the upper interface equation
         RWORK(1) = OD(1)
         RWORK(2) = A_TMP
         RWORK(3) = B_TMP
         RWORK(4) = F_TMP
      ENDIF

! ======================================================================

      RETURN
      END
