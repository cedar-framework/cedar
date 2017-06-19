      SUBROUTINE BMG2_SymStd_residual( &
     &                K, KF, NOG, &
     &                SO, QF, Q, RES,&
     &                II, JJ, IFD, NStncl, IRELAX, IRELAX_SYM,&
     &                iWorkMSG, NMSGi, pMSG, &
     &                MSG_Buffer, NMSGr, MPICOMM&
     &                ) BIND(C, NAME='MPI_BMG2_SymStd_residual')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG_SymStd_resl2 calculates the l2 residual on grid K.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     K         index of the current grid
!     KF        index of the finest grid
!
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
!
!     SO        Refer to BMG2_SymStd_SOLVE_boxmg.
!     QF        Refer to BMG2_SymStd_SOLVE_boxmg.
!
!     RES       Refer to BMG2_SymStd_SOLVE_boxmg.
!     IFD       Refer to BMG2_SymStd_SOLVE_boxmg.
!     IRELAX    Refer to BMG2_SymStd_SOLVE_boxmg.
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!     Q         Refer to BMG2_SymStd_SOLVE_boxmg.
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================
      USE ModInterface
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'


      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER(len_t), VALUE :: II, JJ, NMSGi, NMSGr
      INTEGER(C_INT), VALUE :: NOG, NStncl, IFD, IRELAX, IRELAX_SYM, K, KF
      INTEGER :: MPICOMM
      INTEGER(len_t) :: iWorkMSG(NMSGi)
      INTEGER(C_INT) :: pMSG(NBMG_pMSG,NOG)

      REAL(real_t) :: Q(II,JJ), QF(II,JJ),&
           SO(II+1,JJ+1,NStncl), RES(II,JJ), &
           MSG_Buffer(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, J, J1, ptrn, ierror

! ======================================================================

      J1=JJ-1
      I1=II-1

!     ------------------------------------------------------------------

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         !

         DO J=2,J1
            DO I=2,I1
               RES(I,J) = QF(I,J)&
     &                  + SO(I  ,J  ,KW )*Q(I-1,J)&
     &                  + SO(I+1,J  ,KW )*Q(I+1,J)&
     &                  + SO(I  ,J  ,KS )*Q(I  ,J-1)&
     &                  + SO(I  ,J+1,KS )*Q(I  ,J+1)&
     &                  + SO(I  ,J  ,KSW)*Q(I-1,J-1)&
     &                  + SO(I+1,J  ,KNW)*Q(I+1,J-1)&
     &                  + SO(I  ,J+1,KNW)*Q(I-1,J+1)&
     &                  + SO(I+1,J+1,KSW)*Q(I+1,J+1)&
     &                  - SO(I  ,J  ,KO )*Q(I  ,J)
            ENDDO
         ENDDO

         !
      ELSE
         !
         !  5-point stencil
         !
         DO J=2,J1
            DO I=2,I1
               RES(I,J) = QF(I,J)&
     &                  + SO(I  ,J  ,KW)*Q(I-1,J)&
     &                  + SO(I+1,J  ,KW)*Q(I+1,J)&
     &                  + SO(I  ,J  ,KS)*Q(I  ,J-1)&
     &                  + SO(I  ,J+1,KS)*Q(I  ,J+1)&
     &                  - SO(I  ,J  ,KO)*Q(I  ,J)
            ENDDO
         ENDDO
         !
      ENDIF


! ======================================================================

      RETURN
      END
