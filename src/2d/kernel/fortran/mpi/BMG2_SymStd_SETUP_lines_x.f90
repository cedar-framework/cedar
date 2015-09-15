SUBROUTINE BMG2_SymStd_SETUP_lines_x( SO, SOR, Nx, Ny, NStncl)&
     BIND(C, NAME='MPI_BMG2_SymStd_SETUP_lines_x')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   BMG2_SymStd_SETUP_lines_x performs a factorization of the tridiagona
!   matrix that arises in x-line relaxation.  It assumes that the system
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
      real(real_t) :: SO(Nx+1,Ny+1,NStncl), SOR(Nx,Ny,2)

! ----------------------------
!     Local Declarations
!
      INTEGER i, j

! ======================================================================

      DO j=2, Ny-1
         DO i=1, Nx
            SOR(i,j,2) = -SO(i,j,KW)  ! off diagonal
            SOR(i,j,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

! ======================================================================

      RETURN
      END
