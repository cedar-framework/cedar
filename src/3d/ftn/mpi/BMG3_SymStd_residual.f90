      SUBROUTINE BMG3_SymStd_residual(&
     &                KG, NOG, IFD,&
     &                Q, QF, SO, RES, II, JJ, KK, NStncl,&
     &                halof&
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_residual')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_residual computes the residual on level kg.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
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
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: II, JJ, KK
      integer(c_int), value :: NOG, KG, IFD, NStncl
      real(real_t) :: q(II,JJ,KK), qf(II,JJ,KK), RES(II,JJ,KK)
      real(real_t) :: so(II+1,JJ+1,KK+1,NStncl)
      type(c_ptr) :: halof
! ----------------------------
!     Local Declarations
!
      INTEGER  i, i1, j, j1, k, k1

! ======================================================================

      i1=II-1
      j1=JJ-1
      k1=KK-1

      IF( KG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)&
     &                       + so(i,j,k,kpw)*q(i-1,j,k)&
     &                       + so(i,j+1,k,kpnw)*q(i-1,j+1,k)&
     &                       + so(i,j+1,k,kps)*q(i,j+1,k)&
     &                       + so(i+1,j+1,k,kpsw) *q(i+1,j+1,k)&
     &                       + so(i+1,j,k,kpw)*q(i+1,j,k)&
     &                       + so(i+1,j,k,kpnw)*q(i+1,j-1,k)&
     &                       + so(i,j,k,kps)*q(i,j-1,k)&
     &                       + so(i,j,k,kpsw)*q(i-1,j-1,k)&
     &                       + so(i,j,k,kb)*q(i,j,k-1)&
     &                       + so(i,j,k,kbw)*q(i-1,j,k-1)&
     &                       + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)&
     &                       + so(i,j+1,k,kbn)*q(i,j+1,k-1)&
     &                       + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)&
     &                       + so(i+1,j,k,kbe)*q(i+1,j,k-1)&
     &                       + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)&
     &                       + so(i,j,k,kbs)*q(i,j-1,k-1)&
     &                       + so(i,j,k,kbsw)*q(i-1,j-1,k-1)&
     &                       + so(i,j,k+1,kb)*q(i,j,k+1)&
     &                       + so(i,j,k+1,kbe)*q(i-1,j,k+1)&
     &                       + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)&
     &                       + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)&
     &                       + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)&
     &                       + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)&
     &                       + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)&
     &                       + so(i,j,k+1,kbn)*q(i,j-1,k+1)&
     &                       + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)&
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO

      ELSE

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)&
     &                       + so(i,j,k,kpw)*q(i-1,j,k)&
     &                       + so(i,j+1,k,kps)*q(i,j+1,k)&
     &                       + so(i+1,j,k,kpw)*q(i+1,j,k)&
     &                       + so(i,j,k,kps)*q(i,j-1,k)&
     &                       + so(i,j,k,kb)*q(i,j,k-1)&
     &                       + so(i,j,k+1,kb)*q(i,j,k+1)&
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      call halo_exchange(KG, RES, halof)

! ======================================================================

      RETURN
      END
