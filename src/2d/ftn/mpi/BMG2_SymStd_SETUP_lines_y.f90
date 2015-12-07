SUBROUTINE BMG2_SymStd_SETUP_lines_y( SO, SOR, Nx, Ny, NStncl)&
     BIND(C, NAME='MPI_BMG2_SymStd_SETUP_lines_y')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   BMG2_SymStd_SETUP_lines_y performs a factorization of the tridiagona
!   matrix that arises in y-line relaxation.  It assumes that the system
!   is diagonally dominant and it works directly with the stencil.
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
      integer(len_t), value :: Nx, Ny
      integer(c_int), value :: NStncl
      real(real_t) :: SO(Nx+1,Ny+1,NStncl), SOR(Ny,Nx,2)

! ----------------------------
!     Local Declarations
!
      INTEGER i, j

! ======================================================================

      DO i=2, Nx-1
         DO j=1, Ny
            SOR(j,i,2) = -SO(i,j,KS)  ! off diagonal
            SOR(j,i,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

! ======================================================================

      RETURN
      END
