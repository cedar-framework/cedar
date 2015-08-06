      SUBROUTINE BMG2_SymStd_restrict( &
     &                       Q, QC, CI, &
     &                       Nx, Ny, Nxc, Nyc, JPN&
     &                       ) BIND(C, NAME='BMG2_SymStd_restrict')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_restrict computes the restriction of a vector on the
!     fine grid, Q, to a vector on the coarse grid, QC.  The weights
!     involve the transpose of the interpolation operator from the
!     coarse grid to the fine grid.
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
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER(len_t), VALUE ::  Nx, Nxc, Ny, Nyc

      INTEGER(C_INT), VALUE :: JPN
      REAL(real_t) :: CI(Nxc,Nyc,8), Q(Nx,Ny), QC(Nxc,Nyc)

! ----------------------------
!     Local Declarations
!
      INTEGER  i, ic, IPN, j, jc, PER_x, PER_y, PER_xy

! ======================================================================


! ---------------------------------------------
!     Restrict the vector Q -> QC:
! ---------------------------------------------
      IPN = IABS(JPN)
      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     &     THEN

         j=0
         DO jc=2, Nyc-1
            j=j+2
            i=0
            DO ic=2, Nxc-1
               i=i+2
               QC(ic,jc) = CI(ic,jc,LNE)*Q(i-1,j-1)&
     &              + CI(ic,jc,LA)*Q(i,j-1)&
     &              + CI(ic+1,jc,LNW)*Q(i+1,j-1)&
     &              + CI(ic,jc,LR)*Q(i-1,j)&
     &              + Q(i,j)&
     &              + CI(ic+1,jc,LL)*Q(i+1,j)&
     &              + CI(ic,jc+1,LSE)*Q(i-1,j+1)&
     &              + CI(ic,jc+1,LB)*Q(i,j+1)&
     &              + CI(ic+1,jc+1,LSW)*Q(i+1,j+1)
            ENDDO
         ENDDO

      ELSE

         PER_x = IABS(BMG_BCs_def_per_x)
         PER_y = IABS(BMG_BCs_def_per_y)
         PER_xy = IABS(BMG_BCs_def_per_xy)

         IF((IPN.EQ.PER_y.OR.IPN.EQ.PER_xy).AND.Ny/2+1.EQ.Nyc) THEN
            DO I=1,Nx
               Q(I,1) = Q(I,Ny-1)
               Q(I,Ny) = Q(I,2)
            ENDDO
         ENDIF
         IF((IPN.EQ.PER_x.OR.IPN.EQ.PER_xy).AND.Nx/2+1.EQ.Nxc) THEN
            DO J=1,Ny
               Q(1,J) = Q(Nx-1,J)
               Q(Nx,J) = Q(2,J)
            ENDDO
         ENDIF
         j=0
         DO jc=2, Nyc-1
            j=j+2
            i=0
            DO ic=2, Nxc-1
               i=i+2
               QC(ic,jc) = CI(ic,jc,LNE)*Q(i-1,j-1)&
     &              + CI(ic,jc,LA)*Q(i,j-1)&
     &              + CI(ic+1,jc,LNW)*Q(i+1,j-1)&
     &              + CI(ic,jc,LR)*Q(i-1,j)&
     &              + Q(i,j)&
     &              + CI(ic+1,jc,LL)*Q(i+1,j)&
     &              + CI(ic,jc+1,LSE)*Q(i-1,j+1)&
     &              + CI(ic,jc+1,LB)*Q(i,j+1)&
     &              + CI(ic+1,jc+1,LSW)*Q(i+1,j+1)
            ENDDO
         ENDDO

      ENDIF

! ======================================================================

       RETURN
       END
