      SUBROUTINE BMG2_SymStd_SETUP_lines_y(&
     &                SO, SOR, Nx, Ny, NStncl, JPN&
     &                ) BIND(C,NAME='BMG2_SymStd_SETUP_lines_y')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_lines_y.f performs a factorization of the
!     tridiagonal matrix that arises in y-line relaxation.  It assumes
!     that the system is diagonally dominant and it works directly with
!     the stencil.
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
      INTEGER   INFO, IPN, PER_y, PER_xy

! ======================================================================

      IPN = IABS(JPN)

      DO i=2, Nx-1
         DO j=2, Ny-1
            SOR(j,i,2) = -SO(i,j,KS)  ! off diagonal
            SOR(j,i,1) =  SO(i,j,KO)  ! diagonal
         ENDDO
      ENDDO

      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy) THEN
         DO i = 2, Nx-1
            SOR(2,i,1) = SOR(2,i,1) + SO(i,2,KS)
            SOR(Ny-1,i,1) = SOR(Ny-1,i,1) + SO(i,Ny,KS)
         ENDDO
      ENDIF

!     calculate the L*D*L' factorizations for the lines in y-direction
      DO i=2, Nx-1
         CALL DPTTRF (Ny-2, SOR(2,i,1), SOR(3,i,2), INFO)
      ENDDO



!      DO i=2, Nx-1
!         SOR(i,2,MSOS)=rONE/SO(i,2,KO)
!         DO j=3, Ny-1
!            SOR(i,j,MSOS)=rONE
!     &                   /( SO(i,j,KO) - SOR(i,j-1,MSOS)*SO(i,j,KS)**2
!         ENDDO
!      ENDDO

! ======================================================================

      RETURN
      END
