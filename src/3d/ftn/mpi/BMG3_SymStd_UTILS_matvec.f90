      SUBROUTINE BMG3_SymStd_UTILS_matvec( &
     &                       KG, SO, QF, Q, II, JJ, KK,&
     &                       NOG, IFD, NStncl, halof&
     &                       ) BIND(C, NAME='MPI_BMG3_SymStd_UTILS_matvec')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_UTILS_matvec performs a matrix multiplication that is
!     used within the pcg routine called BMG3_SymStd_SOLVE_pcg.
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

      integer(len_t), VALUE :: II, JJ, KK
      integer(c_int), VALUE :: IFD, KG, NOG,&
           NStncl
      real(real_t) :: Q(II,JJ,KK), QF(II,JJ,KK)
      real(real_t) :: SO(II+1,JJ+1,KK+1,NStncl)
      type(c_ptr) :: halof

! ----------------------------
!     Local Declarations
!
      integer i, i1, j, j1, k, k1

! ======================================================================

      i1=II-1
      j1=JJ-1
      k1=KK-1

      IF( KG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  QF(i,j,k) =  SO( i , j , k ,kp  )*Q( i , j , k )&
     &                       - SO( i , j , k ,kpw )*Q(i-1, j , k )&
     &                       - SO( i ,j+1, k ,kpnw)*Q(i-1,j+1, k )&
     &                       - SO( i ,j+1, k ,kps )*Q( i ,j+1, k )&
     &                       - SO(i+1,j+1, k ,kpsw)*Q(i+1,j+1, k )&
     &                       - SO(i+1, j , k ,kpw )*Q(i+1, j , k )&
     &                       - SO(i+1, j , k ,kpnw)*Q(i+1,j-1, k )&
     &                       - SO( i , j , k ,kps )*Q( i ,j-1, k )&
     &                       - SO( i , j , k ,kpsw)*Q(i-1,j-1, k )&
     &                       - SO( i , j , k ,kb  )*Q( i , j ,k-1)&
     &                       - SO( i , j , k ,kbw )*Q(i-1, j ,k-1)&
     &                       - SO( i ,j+1, k ,kbnw)*Q(i-1,j+1,k-1)&
     &                       - SO( i ,j+1, k ,kbn )*Q( i ,j+1,k-1)&
     &                       - SO(i+1,j+1, k ,kbne)*Q(i+1,j+1,k-1)&
     &                       - SO(i+1, j , k ,kbe )*Q(i+1, j ,k-1)&
     &                       - SO(i+1, j , k ,kbse)*Q(i+1,j-1,k-1)&
     &                       - SO( i , j , k ,kbs )*Q( i ,j-1,k-1)&
     &                       - SO( i , j , k ,kbsw)*Q(i-1,j-1,k-1)&
     &                       - SO( i , j ,k+1,kb  )*Q( i , j ,k+1)&
     &                       - SO( i , j ,k+1,kbe )*Q(i-1, j ,k+1)&
     &                       - SO( i ,j+1,k+1,kbse)*Q(i-1,j+1,k+1)&
     &                       - SO( i ,j+1,k+1,kbs )*Q( i ,j+1,k+1)&
     &                       - SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)&
     &                       - SO(i+1, j ,k+1,kbw )*Q(i+1, j ,k+1)&
     &                       - SO(i+1, j ,k+1,kbnw)*Q(i+1,j-1,k+1)&
     &                       - SO( i , j ,k+1,kbn )*Q( i ,j-1,k+1)&
     &                       - SO( i , j ,k+1,kbne)*Q(i-1,j-1,k+1)
               ENDDO
            ENDDO
         ENDDO

      ELSEIF( KG.EQ.NOG .AND. IFD.EQ.BMG_STENCIL_7pt ) THEN

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  QF(i,j,k) =  SO( i , j , k ,kp )*Q( i , j , k )&
     &                       - SO( i , j , k ,kpw)*Q(i-1, j , k )&
     &                       - SO( i ,j+1, k ,kps)*Q( i ,j+1, k )&
     &                       - SO(i+1, j , k ,kpw)*Q(i+1, j , k )&
     &                       - SO( i , j , k ,kps)*Q( i ,j-1, k )&
     &                       - SO( i , j , k ,kb )*Q( i , j ,k-1)&
     &                       - SO( i , j ,k+1,kb )*Q( i , j ,k+1)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

! ======================================================================

      call halo_exchange(KG, Q, halof)

! ======================================================================

      RETURN
      END
