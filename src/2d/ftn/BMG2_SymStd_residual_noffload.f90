      SUBROUTINE BMG2_SymStd_residual_noffload( &
     &                K, SO, QF, Q, RES, II, JJ,&
     &                KF, IFD, NStncl, IBC,&
     &                IRELAX, IRELAX_SYM, UPDOWN,&
     &                jbeg, jend, ndev&
     &                ) BIND(C,NAME='BMG2_SymStd_residual_noffload')

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
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), intent(in) :: II, JJ
      integer(c_int), value, intent(in) :: ndev
      integer(len_t), intent(in) :: jbeg(ndev), jend(ndev)
      integer(c_int), intent(in) :: NStncl

      integer(c_int), intent(in) :: IBC, IFD, IRELAX, IRELAX_SYM, K, KF, UPDOWN
      real(real_t), intent(in) :: Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl)
      real(real_t), intent(out) :: RES(II,JJ)

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, J, J1
      integer(len_t) :: lbeg, lend
      integer(c_int) ::  di

! ======================================================================

      J1=JJ-1
      I1=II-1
!     ------------------------------------------------------------------

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         !
         !$omp parallel do num_threads(ndev) private(lbeg, lend)
         do di=1,ndev
            lbeg = jbeg(di)
            lend = jend(di)
            !$omp target teams distribute parallel do simd collapse(2) device(di-1)
            DO J=lbeg,lend
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
            !$omp end target teams distribute parallel do
         enddo
         !$omp end parallel do
      ELSE
         !
         !  5-point stencil
         !
         !$omp target teams distribute parallel do simd collapse(2)
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
         !$omp end target teams distribute parallel do
      ENDIF

! ======================================================================

      RETURN
      END
