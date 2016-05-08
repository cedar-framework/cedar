      SUBROUTINE BMG2_SymStd_SETUP_recip( &
     &                SO, SOR, Nx, Ny, NStncl, NSOR_v&
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_recip')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_cons_recip.f computes the reciprocal of the central
!     stencil coefficient.
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
      integer(len_t), value :: Nx, Ny
      integer(c_int), value :: NSOR_v, NStncl
      real(real_t) :: SO(Nx,Ny,NStncl), SOR(Nx,Ny,NSOR_v)

! --------------------------
!     Local Declarations:
!
      integer(len_t) :: i, j

! ======================================================================

      DO j = 2, Ny-1
         DO i = 2, Nx-1
            SOR(i,j,msor)=rONE/SO(i,j,ko)
         ENDDO
      ENDDO

! ======================================================================

! =======================

         RETURN
         END