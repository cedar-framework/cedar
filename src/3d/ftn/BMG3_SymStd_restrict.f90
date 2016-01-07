      SUBROUTINE BMG3_SymStd_restrict(&
           Q, QC, CI, Nx, Ny, Nz, Nxc, Nyc, Nzc, JPN &
           ) BIND(C, NAME='BMG3_SymStd_restrict')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_restrict computes the restriction of a vector on the
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
      integer(c_int), value :: JPN
      integer(len_t), value :: Nx, Nxc, Ny, Nyc, Nz, Nzc
      real(real_t) :: CI(Nxc,Nyc,Nzc,26), Q(Nx,Ny,Nz), QC(Nxc,Nyc,Nzc)

! ----------------------------
!     Local Declarations
!
      INTEGER  i, ic, IPN, j, jc, k, kc
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz
      REAL*8 sum

! ======================================================================

      IPN = IABS(JPN)
      PER_x = BMG_BCs_def_per_x
      PER_y = BMG_BCs_def_per_y
      PER_xy = BMG_BCs_def_per_xy
      PER_z = BMG_BCs_def_per_z
      PER_xz = BMG_BCs_def_per_xz
      PER_yz = BMG_BCs_def_per_yz
      PER_xyz = BMG_BCs_def_per_xyz


      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_XY .OR. IPN.EQ.PER_YZ&
     &     .OR. IPN.EQ.PER_XYZ ) THEN
         DO K = 1,Nz
            DO I = 1,Nx
               Q(I,1,K) = Q(I,Ny-1,K)
               Q(I,Ny,K) = Q(I,2,K)
            ENDDO
         ENDDO
      ENDIF

         IF(IPN.EQ.PER_X .OR. IPN.EQ.PER_XY .OR. IPN.EQ.PER_XZ&
     &        .OR. IPN.EQ.PER_XYZ ) THEN
            DO K = 1,Nz
               DO J = 1,Ny
                  Q(1,J,K) = Q(Nx-1,J,K)
                  Q(Nx,J,K) = Q(2,J,K)
               ENDDO
            ENDDO
         ENDIF

            IF(IPN.EQ.PER_Z .OR. IPN.EQ.PER_XZ .OR. IPN.EQ.PER_YZ&
     &           .OR. IPN.EQ.PER_XYZ ) THEN
               DO J = 1,Ny
                  DO I = 1,Nx
                     Q(I,J,1) = Q(I,J,Nz-1)
                     Q(I,J,Nz) =Q(I,J,2)
                  ENDDO
               ENDDO
            ENDIf

! ---------------------------------------------
!     Restrict the vector Q -> QC:
! ---------------------------------------------


            k=0
            DO kc=2, Nzc-1
               k=k+2
               j=0
               DO jc=2, Nyc-1
                  j=j+2
                  i=0
                  DO ic=2, Nxc-1
                     i=i+2
                     QC(ic,jc,kc) = CI(ic,jc,kc,lxyne)*Q(i-1,j-1,k)&
     &                    + CI(ic,jc,kc,lxya)*Q(i,j-1,k)&
     &                    + CI(ic+1,jc,kc,lxynw)*Q(i+1,j-1,k)&
     &                    + CI(ic,jc,kc,lxyr)*Q(i-1,j,k)&
     &                    + Q(i,j,k)&
     &                    + CI(ic+1,jc,kc,lxyl)*Q(i+1,j,k)&
     &                    + CI(ic,jc+1,kc,lxyse)*Q(i-1,j+1,k)&
     &                    + CI(ic,jc+1,kc,lxyb)*Q(i,j+1,k)&
     &                    + CI(ic+1,jc+1,kc,lxysw)*Q(i+1,j+1,k)&
     &                    + CI(ic,jc,kc,ltne)*Q(i-1,j-1,k-1)&
     &                    + CI(ic,jc,kc,lyznw)*Q(i,j-1,k-1)&
     &                    + CI(ic+1,jc,kc,ltnw)*Q(i+1,j-1,k-1)&
     &                    + CI(ic,jc,kc,lxzne)*Q(i-1,j,k-1)&
     &                    + CI(ic,jc,kc,lxza)*Q(i,j,k-1)&
     &                    + CI(ic+1,jc,kc,lxznw)*Q(i+1,j,k-1)&
     &                    + CI(ic,jc+1,kc,ltse)*Q(i-1,j+1,k-1)&
     &                    + CI(ic,jc+1,kc,lyzne)*Q(i,j+1,k-1)&
     &                    + CI(ic+1,jc+1,kc,ltsw)*Q(i+1,j+1,k-1)&
     &                    + CI(ic,jc,kc+1,lbne)*Q(i-1,j-1,k+1)&
     &                    + CI(ic,jc,kc+1,lyzsw)*Q(i,j-1,k+1)&
     &                    + CI(ic+1,jc,kc+1,lbnw)*Q(i+1,j-1,k+1)&
     &                    + CI(ic,jc,kc+1,lxzse)*Q(i-1,j,k+1)&
     &                    + CI(ic,jc,kc+1,lxzb)*Q(i,j,k+1)&
     &                    + CI(ic+1,jc,kc+1,lxzsw)*Q(i+1,j,k+1)&
     &                    + CI(ic,jc+1,kc+1,lbse)*Q(i-1,j+1,k+1)&
     &                    + CI(ic,jc+1,kc+1,lyzse)*Q(i,j+1,k+1)&
     &                    + CI(ic+1,jc+1,kc+1,lbsw)*Q(i+1,j+1,k+1)
                  enddo

               enddo

            enddo

            sum = rZERO
            DO kc=2, Nzc-1
               DO jc = 2,Nyc-1
                  DO ic = 2,Nxc-1
                     SUM = SUM + QC(ic,jc,kc)
                  ENDDO
               ENDDO
            ENDDO

            !WRITE(*,*) 'sum = ',sum

! ======================================================================

      return
      end
