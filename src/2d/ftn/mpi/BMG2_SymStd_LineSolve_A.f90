      SUBROUTINE BMG2_SymStd_LineSolve_A( &
     &                SOR, Q, II, JJ, JBEG, RWORK, NP, NLINES&
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
      INTEGER  NP, NLINES, JBEG, II, JJ
      REAL*8 SOR(II,JJ,2), Q(II,JJ)
      REAL*8 RWORK(*)

! ----------------------------
!     Local Declarations
!
      REAL*8  L,A_TMP, F_TMP
      REAL*8  C_TMP, B_TMP
      INTEGER I, J, OFFSET

! ======================================================================

      NLINES = 0
      OFFSET = 0

      ! RWORK(OFFSET+1) is the lower off diagonal element
      ! RWORK(OFFSET+2) is the diagonal element
      ! RWORK(OFFSET+3) is the upper off diagonal element
      ! RWORK(OFFSET+4) is the right hand side

      ! JBEG is the index of the first line that is solved
      !      it is computed one level up where the selection
      !      of red or black lines is detemined

      ! NP - THE NUMBER OF UNKNOWNS (OR EQUATIONS) per line ON THIS PROC
      ! NLINES - number of lines to be solved on this process row

      IF (NP.eq.1) THEN

         ! special case: one local equation on this processor
         DO J=JBEG,JJ-1,2 ! loop over all lines

            RWORK(OFFSET+1) = SOR(2,J,2)
            RWORK(OFFSET+2) = SOR(2,J,1)
            RWORK(OFFSET+3) = SOR(3,J,2)
            RWORK(OFFSET+4) = Q(2,J)

            ! count the lines
            NLINES = NLINES+1
            OFFSET = OFFSET+8
         ENDDO

      ELSE IF (NP.eq.2) THEN

         ! special case: two local equations on this processor
         DO J=JBEG,JJ-1,2

            RWORK(OFFSET+1) = SOR(2,J,2)
            RWORK(OFFSET+2) = SOR(2,J,1)
            RWORK(OFFSET+3) = SOR(3,J,2)
            RWORK(OFFSET+4) = Q(2,J)

            RWORK(OFFSET+5) = SOR(3,J,2)
            RWORK(OFFSET+6) = SOR(3,J,1)
            RWORK(OFFSET+7) = SOR(4,J,2)
            RWORK(OFFSET+8) = Q(3,J)

            NLINES = NLINES+1
            OFFSET = OFFSET+8
         ENDDO

      ELSE

         ! general case: more than two equation on this processor
         DO J=JBEG,JJ-1,2

            ! initialize with coefficients from the second equation
            ! in the line

            C_TMP = SOR(3,J,2)
            F_TMP = Q(3,J)
            A_TMP = SOR(3,J,1)

            ! start iterating on the third equation in the line
            ! note that we reference i+1 in the body of the loop
            DO I=3,NP,1
               ! compute the multiplier (lower off diagonal
               ! of euqation i+1 divided by the diagonal of the
               ! previous equation
               L =  -SOR(I+1,J,2) / A_TMP

               ! compute the non-zero entry in column one of the
               ! matrix local to this processor
               C_TMP = L*C_TMP

               ! compute the right hand side
               F_TMP = Q(I+1,J) + L*F_TMP

               ! compute the diagonal element
               A_TMP = SOR(I+1,J,1) + L*SOR(I+1,J,2)
            ENDDO

            ! populate RWORK with the lower interface equation
            RWORK(OFFSET+5) = C_TMP
            RWORK(OFFSET+6) = A_TMP
            RWORK(OFFSET+7) = SOR(NP+2,J,2)
            RWORK(OFFSET+8) = F_TMP


            ! initialize with coefficients from the second to last
            ! equation in the line
            A_TMP = SOR(NP,J,1)
            B_TMP = SOR(NP+1,J,2)
            F_TMP = Q(NP,J)

            ! start iterating on the third to last equation in the line
            ! note that we reference i+1 in the body of the loop
            DO I=NP-2,1,-1
               ! compute the multiplier (opper off diagonal
               ! of euqation i+1 divided by the diagonal of the
               ! previous equation
               L =  -SOR(I+2,J,2) / A_TMP

               ! compute the non-zero entry in column NP of the
               ! matrix local to this processor
               B_TMP = L * B_TMP

               ! compute the right hand side
               F_TMP = Q(I+1,J) + L*F_TMP

               ! compute the diagonal element
               A_TMP = SOR(I+1,J,1) + L*SOR(I+2,J,2)
            ENDDO

            ! populate RWORK with the upper interface equation
            RWORK(OFFSET+1) = SOR(2,J,2)
            RWORK(OFFSET+2) = A_TMP
            RWORK(OFFSET+3) = B_TMP
            RWORK(OFFSET+4) = F_TMP

            ! count the lines
            NLINES = NLINES + 1

            ! increment the offset for the next line
            OFFSET = OFFSET + 8
         ENDDO

      ENDIF

! ======================================================================

      RETURN
      END
