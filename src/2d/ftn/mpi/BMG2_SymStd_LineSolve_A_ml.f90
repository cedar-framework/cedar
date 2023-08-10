      SUBROUTINE BMG2_SymStd_LineSolve_A_ml(&
     &                D, OD, Q, IED, NPts,&
     &                FLG&
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
!     Npts  - Number of local unknowns
!
!     IED   - Buffer to store interface equation nonzeros and rhs
!             IED(OFFSET+1) - the lower off diagonal element
!             IED(OFFSET+2) - the diagonal element
!             IED(OFFSET+3) - the upper off diagonal element
!             IED(OFFSET+4) - the right hand side
!
!     NLINES - Number of lines to be independently relaxed
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

      use iso_c_binding, only: c_bool
      IMPLICIT NONE

! ---------------------------
!    Includes:
!
      INCLUDE 'mpif.h'

! ---------------------------
!    Argument Declarations:
!
      INTEGER  NPts
      REAL*8   D(NPts), OD(NPts+1), Q(NPts)
      REAL*8   IED(8)
      logical(c_bool) :: FLG

! ----------------------------
!     Local Declarations
!
      REAL*8  RHS(NPts-2), Q_TMP
      INTEGER I, J, NLINES, OFFSET
      INTEGER INFO

      INTEGER MyRank, IERR


! ======================================================================

      CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IERR)

! ======================================================================


      NLINES = 0
      OFFSET = 0
      IF (NPts.eq.1) THEN

         ! special case: one local equation on this processor

         IED(1) = OD(1)
         IED(2) = D(1)
         IED(3) = OD(2)
         IED(4) = Q(1)

         ! Set flag to tell routine that there was only 1 point
         IED(6) = -1.0

      ELSE IF (NPts.eq.2) THEN

         ! special case: two local equations on this processor

         IED(1) = OD(1)
         IED(2) = D(1)
         IED(3) = OD(2)
         IED(4) = Q(1)


         IED(5) = OD(2)
         IED(6) = D(2)
         IED(7) = OD(3)
         IED(8) = Q(2)

      ELSE

        ! general case: more than two equation on this processor

!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     This stuff only needs to happen once on each level
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !
        ! Get the Cholesky decomposition of the submatrix
        ! and store it over SOR
        !

        IF (FLG) THEN
           CALL DPTTRF(NPts-2, D(2), OD(3), INFO)
        ENDIF

        !Initialize RHS vector to ZERO
        !
        DO I=1,NPts-2
           RHS(I) = 0D0
        END DO

        !
        ! Set rhs vector to solve for lower interface coefficients
        !
        RHS(1) = -OD(2)

        !
        ! Solve subsystem to get alpha's (for upper interface eqn)
        !
        CALL DPTTRS(NPts-2, 1,  D(2), OD(3), RHS(1), NPts-2, INFO)

        !
        ! Should try to store this RHS vector in SOR somewhere
        ! Then we only have to compute the alpha's and beta's
        ! once.
        !


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !
        ! Compute upper interface eqn and populate IED
        !
        IED(1) = OD(1)
        IED(2) = D(1) + RHS(1)*OD(2)
        IED(3) = RHS(NPts-2)*OD(NPts)

        Q_TMP = Q(1)

        DO I=2,NPts-1
           Q_TMP = Q_TMP + Q(I)*RHS(I-1)
        END DO

        IED(4) = Q_TMP

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     This stuff only needs to happen once on each level
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !
        !Initialize RHS vector to ZERO
        !
        DO I=1,NPts-2
           RHS(I) = 0D0
        END DO

        !
        ! Set rhs vector to solve for lower interface coefficients
        !
        RHS(NPts-2) = -OD(NPts)

        !
        ! Solve subsystem to get betas's
        !
        CALL DPTTRS(NPts-2, 1,  D(2), OD(3), RHS(1), NPts-2, INFO)

        !
        ! Should try to store this RHS vector in SOR somewhere
        ! Then we only have to compute the alpha's and beta's
        ! once.
        !

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !
        ! Compute lower interface eqn and populate IED
        !

        IED(5) = RHS(1)*OD(2)
        IED(6) = D(NPts) + RHS(NPts-2)*OD(NPts)
        IED(7) = OD(NPts+1)

        Q_TMP = Q(NPts)

        DO I=2,NPts-1
           Q_TMP = Q_TMP + Q(I)*RHS(I-1)
        END DO

        IED(8) = Q_TMP

      ENDIF

! ======================================================================

      RETURN
      END
