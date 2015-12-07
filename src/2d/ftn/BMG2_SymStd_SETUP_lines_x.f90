      SUBROUTINE BMG2_SymStd_SETUP_lines_x(&
     &                SO, SOR, Nx, Ny, NStncl, JPN&
     &                ) BIND(C,NAME='BMG2_SymStd_SETUP_lines_x')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_lines_x performs a factorization of the tridiago
!     matrix that arises in x-line relaxation.  It assumes that the syst
!     is diagonally dominant and it works directly with the stencil.
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
      INCLUDE 'BMG_parameters_f90.h'
      INCLUDE 'BMG_stencils_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), VALUE :: Nx, Ny
      integer(c_int), VALUE :: NStncl, JPN
      real(real_t) :: SO(Nx,Ny,NStncl), SOR(Nx,Ny,2)

! ----------------------------
!     Local Declarations
!
      INTEGER   i, j
      INTEGER   INFO, IPN, PER_x, Per_xy

! ======================================================================

      IPN = IABS(JPN)

      DO j=2, Ny-1
         DO i=2, Nx-1
            SOR(i,j,2) = -SO(i,j,KW)  ! off diagonal
            SOR(i,j,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

      PER_x = IABS(BMG_BCs_def_per_x)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy) THEN
         DO j = 2, Ny-1
            SOR(2,j,1) = SOR(2,j,1) + SO(2,j,KW)
            SOR(Nx-1,j,1) = SOR(Nx-1,j,1) + SO(Nx,j,KW)
         ENDDO
      ENDIF

!     calculate the L*D*L' factorizations for the lines in x-direction
      DO j=2,Ny-1
         CALL DPTTRF (Nx-2, SOR(2,j,1), SOR(3,j,2), INFO)
      ENDDO


!      DO j=2, Ny-1
!         SOR(2,j,MSOR)=rONE/SO(2,j,KO)
!         DO i=3, Nx-1
!            SOR(i,j,MSOR)=rONE
!     &                   /( SO(i,j,KO) - SOR(i-1,j,MSOR)*SO(i,j,KW)**2
!         ENDDO
!      ENDDO

! ======================================================================

      RETURN
      END
