      SUBROUTINE BMG3_SymStd_residual(&
     &                kg, NOG, ifd, q, qf, so, RES, ii, jj, kk, NStncl&
     &                ) BIND(C, NAME='BMG3_SymStd_residual')

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
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: ii, jj, kk
      integer(c_int), value :: ifd, kg, NOG, NStncl
      real(real_t) :: q(ii,jj,kk), qf(ii,jj,kk), RES(ii,jj,kk)
      real(real_t) :: SO(ii,jj,kk,NStncl)

! ----------------------------
!     Local Declarations
!
      integer  i, i1, j, j1, k, k1

! ======================================================================

      i1=ii-1
      j1=jj-1
      k1=kk-1

      IF( kg.lt.NOG .or. IFD.ne.1 ) THEN

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

! ======================================================================

      RETURN
      END
