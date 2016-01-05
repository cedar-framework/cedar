      SUBROUTINE BMG3_SymStd_SETUP_recip( &
     &                SO, SOR, Nx, Ny, Nz, NStncl, NSORv&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_recip')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_recip computes the reciprocal of the central
!     stencil coefficient for use in GS relaxation.
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

! ---------------------------
!    Argument Declarations:
!
      integer(c_int), value :: NSORv, NStncl
      integer(len_t), value :: Nx, Ny, Nz
      real(real_t) :: SO(Nx,Ny,Nz,NStncl), SOR(Nx,Ny,Nz,NSORv)

! --------------------------
!     Local Declarations:
!
      INTEGER   i, j, k

! ======================================================================

      DO k=2, Nz-1
         DO j = 2, Ny-1
            DO i = 2, Nx-1
               SOR(i,j,k,msor)=rONE/SO(i,j,k,kp)
            ENDDO
         ENDDO
      ENDDO

! ======================================================================

      RETURN
      END
